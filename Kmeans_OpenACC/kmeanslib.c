#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "kmeanslib.h"


/*
 * Read File function
 *
 * mImage->pixels = malloc(mImage->length * sizeof(rgb));
 *
 */
int read_file(char *name, image* mImage)
{
      if((mImage->fp=fopen(name, "r")) == NULL) {
            printf("No fue posible abrir el BMP.\n");
            exit(6);
      }
      if(fseek(mImage->fp, 18, SEEK_SET) != 0) {
            printf("Error fseek.\n");
            exit(6);
      }
      if(fread(&mImage->width, sizeof(mImage->width), 1, mImage->fp) != 1) {
            printf("Error Writing File.\n");
            exit(6);
      }
      if(fseek(mImage->fp, 22, SEEK_SET) != 0) {
            printf("Error fseek.\n");
            exit(6);
      }
      if(fread(&mImage->height, sizeof(mImage->height), 1, mImage->fp) != 1) {
            printf("Error Reading H File.\n");
            exit(6);
      }
 
      mImage->length = (mImage->width * mImage->height);
 
      mImage->pixels  = malloc(mImage->length * sizeof(rgb));

      if(mImage->pixels == NULL){
            printf("Not enough memory.\n");
            return 6;
      }

      if(fseek(mImage->fp, 0, SEEK_SET) != 0){
            printf("Error fseek.\n");
            return 6;
      }

      if(fread(&mImage->header, sizeof(uint8_t), sizeof(mImage->header), mImage->fp) != sizeof(mImage->header)){
            printf("Error reading header.\n");
            return 6;
      }

      int i,c;
 
      for(i = 0; i < mImage->length; i++){
            c = getc(mImage->fp);
            if(c == EOF){
                  printf("Error Reading File.\n");
                  return 5;
            }
            mImage->pixels[i].b = c;
            c = getc(mImage->fp);
            if(c == EOF){
                  printf("Error Reading File.\n");
                  return 5;
            }
            mImage->pixels[i].g = c;
            c = getc(mImage->fp);
            if(c == EOF){
                  printf("Error Reading File.\n");
                  return 5;
            }
            mImage->pixels[i].r = c;
      }
 
      if(fclose(mImage->fp) != 0){
            printf("Error cerrando fichero.\n");
            return 1;
      }
      return (0);
}

int write_file(char *name, image *mImage, cluster*centroids, uint8_t k){
      FILE *fp;
      int i,j;
      uint8_t closest;
      rgb pixel;

      if((fp=fopen(name, "w")) == NULL){
                  printf("Cannot create output file, do you have write permissions for current directory?.\n");
                  return 1;
            }
            if(fwrite(mImage->header, sizeof(uint8_t), sizeof(mImage->header), fp) != sizeof(mImage->header)){
                  printf("Error writing BMP.\n");
                  return 1;
            }
            for(i = 0; i < mImage->width; i++){
                  for(j = 0; j < mImage->height; j++){
                        closest = find_closest_centroid(&mImage->pixels[i * mImage->height + j], centroids, k);
                        pixel.r = centroids[closest].r;
                        pixel.g = centroids[closest].g;
                        pixel.b = centroids[closest].b;
                        
                        if(fwrite(&pixel, sizeof(uint8_t), 3, fp) != 3){
                              printf("Error writing BMP.\n");
                              return 1;
                        }
                  }
            }
            if(fclose(fp) != 0){
                  printf("Error closing file.\n");
                  return 1;
            }
            if((fp=fopen("clusters.txt", "w")) == NULL){
                  printf("Clusters.txt cannot be created. Do you have space | owner in current directory?\n");
                  return 1;
            }
            
            fprintf(fp, "Red\tGreen\tBlue\n");

            for(j = 0; j < k; j++){
                  fprintf(fp, "%d\t%d\t%d\n", centroids[j].r, centroids[j].g, centroids[j].b);
            }

            if(fclose(fp) != 0){
                  printf("Error closing fichero.\n");
                  return 1;
            }
}

/*
 * Compute Checksum
 * Checksum is the total sum of the product (r,g,b) of the centroids final value.
 *
 */
uint32_t getChecksum(cluster* centroids, uint8_t k){
  uint32_t i,j;
  uint32_t sum = 0;
 
  for (i = 0; i < k; i++)
  {
    printf("Centroide %u : R[%u]\tG[%u]\tB[%u]\n", i, centroids[i].r, centroids[i].g, centroids[i].b);
    sum += (centroids[i].r * centroids[i].g * centroids[i].b);
  }
  return sum;
}

/*
 * Return the index number of the closest centroid to a given pixel (p)
 * input @param: p                        --> pixel pointer
 * input @param: centroids          --> cluster array pointer
 * input @param: num_clusters --> amount of centroids in @param: centroids
 *
 * output: uint8_t --> Index of the closest centroid
 *
 */
uint8_t find_closest_centroid(rgb* p, cluster* centroids, uint8_t num_clusters){
      uint32_t min = UINT32_MAX;
      uint32_t dis[num_clusters];
      uint8_t closest = 0, j;
      int16_t diffR, diffG, diffB;  

      for(j = 0; j < num_clusters; j++)
      {
            diffR = centroids[j].r - p->r;
            diffG = centroids[j].g - p->g;
            diffB = centroids[j].b - p->b;
            // No sqrt required.
            dis[j] = diffR*diffR + diffG*diffG + diffB*diffB;
            
            if(dis[j] < min)
            {
                  min = dis[j];
                  closest = j;
            }
      }
      return closest;
}

/*
 * Main function k-means      
 * input @param: K                  --> number of clusters
 * input @param: centroides --> the centroids
 * input @param: num_pixels --> number of total pixels
 * input @param: pixels       --> pinter to array of rgb (pixels)
 */

void kmeans(uint8_t k, cluster* centroides, uint32_t num_pixels, rgb* pixels) {

    uint8_t condition, changed;
    uint32_t i, j, random_num;
    uint32_t* red = malloc(k * sizeof(uint32_t));
    uint32_t* green = malloc(k * sizeof(uint32_t));
    uint32_t* blue = malloc(k * sizeof(uint32_t));
    uint32_t* points = malloc(k * sizeof(uint32_t));
    uint32_t* pixels_aux = malloc(num_pixels * sizeof(uint32_t));

    printf("STEP 1: K = %d\n", k);

    k = MIN(k, num_pixels);
    printf("STEP 2: Init centroids\n");
    for(i = 0; i < k; i++) {
        random_num = rand() % num_pixels;
        centroides[i].r = pixels[random_num].r;
        centroides[i].g = pixels[random_num].g;
        centroides[i].b = pixels[random_num].b;
    }

    printf("STEP 3: Updating centroids\n\n");
    i = 0;
    do {
        #pragma acc data copy(centroides[0:k], pixels[0:num_pixels], pixels_aux[0:num_pixels])
        {
            #pragma acc parallel loop present(centroides[0:k])
            for (j = 0; j < k; j++) {

                centroides[j].media_r = 0;
                centroides[j].media_g = 0;
                centroides[j].media_b = 0;
                centroides[j].num_puntos = 0;
                red[j] = 0;
                green[j] = 0;
                blue[j] = 0;
                points[j] = 0;
            }
            #pragma acc parallel loop present(centroides[0:k], pixels[0:num_pixels], pixels_aux[0:num_pixels])

            for (j = 0; j < num_pixels; j++) {

                int16_t diffR, diffG, diffB;
                uint32_t aux, min = UINT32_MAX, m;
                for (m = 0; m < k; m++) {
                    diffR = centroides[m].r - pixels[j].r;
                    diffG = centroides[m].g - pixels[j].g;
                    diffB = centroides[m].b - pixels[j].b;
                    aux = diffR * diffR + diffG * diffG + diffB * diffB;
                    if (aux < min) {
                        min = aux;
                        pixels_aux[j] = m;
                    }
                }
            }
            }




            for (j = 0; j < num_pixels; j++) {
                uint32_t n_cluster = pixels_aux[j];
                red[n_cluster] += pixels[j].r;
                green[n_cluster] += pixels[j].g;
                blue[n_cluster] += pixels[j].b;
                points[n_cluster] += 1;
            }

            for (j = 0; j < k; j++) {
                centroides[j].media_r = red[j];
                centroides[j].media_g = green[j];
                centroides[j].media_b = blue[j];
                centroides[j].num_puntos = points[j];
            }

            condition = 0;
            for (j = 0; j < k; j++) {
                if (centroides[j].num_puntos == 0) {
                    continue;
                }
                centroides[j].media_r = centroides[j].media_r / centroides[j].num_puntos;
                centroides[j].media_g = centroides[j].media_g / centroides[j].num_puntos;
                centroides[j].media_b = centroides[j].media_b / centroides[j].num_puntos;
                changed = centroides[j].media_r != centroides[j].r || centroides[j].media_g != centroides[j].g || centroides[j].media_b != centroides[j].b;
                condition = condition || changed;
                centroides[j].r = centroides[j].media_r;
                centroides[j].g = centroides[j].media_g;
                centroides[j].b = centroides[j].media_b;
            }
            i++;

    } while (condition);
    printf("Number of K-Means iterations: %d\n\n", i);
    free(red);
    free(green);
    free(blue);
    free(points);
    free(pixels_aux);

}
