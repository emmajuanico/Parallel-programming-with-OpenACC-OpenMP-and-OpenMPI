#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

float stencil ( float v1, float v2, float v3, float v4 )
{
  return (v1 + v2 + v3 + v4) / 4;
}

void laplace_step ( float *in, float *out, int n, int m)
{
  int i, j;
  for ( i=1; i < n-1; i++ )
    for ( j=1; j < m-1; j++ )
      out[i*m+j]= stencil(in[i*m+j+1], in[i*m+j-1], in[(i-1)*m+j], in[(i+1)*m+j]);
}

float laplace_error ( float *old, float *new, int n, int m )
{
  int i, j;
  float error=0.0f;
  for ( i=1; i < n-1; i++ )
    for ( j=1; j < m-1; j++ )
      error = fmaxf( error, sqrtf( fabsf( old[i*m+j] - new[i*m+j] )));
  return error;
}

void laplace_copy ( float *in, float *out, int n, int m )
{
  int i, j;
  for ( i=1; i < n-1; i++ )
    for ( j=1; j < m-1; j++ )
      out[i*m+j]= in[i*m+j];
}


void laplace_init ( float *in, int n, int m )
{
  int i, j;
  const float pi  = 2.0f * asinf(1.0f);
  memset(in, 0, n*m*sizeof(float));
  for (j=0; j<m; j++)  in[    j    ] = 0.f;
  for (j=0; j<m; j++)  in[(n-1)*m+j] = 0.f;
  for (i=0; i<n; i++)  in[   i*m   ] = sinf(pi*i / (n-1));
  for (i=0; i<n; i++)  in[ i*m+m-1 ] = sinf(pi*i / (n-1))*expf(-pi);
}

int main(int argc, char** argv)
{
  int n = 4096, m = 4096;  //mida constant
  const float pi  = 2.0f * asinf(1.0f);
  const float tol = 3.0e-3f;
  int world_size, rank;

  float error= 1.0f, total;

  int i, j, iter_max=1000, iter=0;
  float *A, *Anew;
 
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);
  MPI_Status status;

  // get runtime arguments: n, m and iter_max
  if (argc>1) {  n        = atoi(argv[1]); } //files
  if (argc>2) {  m        = atoi(argv[2]); } //columnes
  if (argc>3) {  iter_max = atoi(argv[3]); } //numero iteracions
  
  

  int num_files = n/world_size; //inicialitzem num_files per a cada proces
  
  //inicialitzem fila inicial i final
 
  int inicial= rank*num_files;  
 
  int final = (rank+1)*num_files - 1;
 
  int init = (rank)? (inicial-1)*m:0;  //inicialitzem posicio inicial

 
 // reservem memoria per a A i Anew
  A    = (float*) malloc( n*m*sizeof(float) );
  Anew = (float*) malloc( n*m*sizeof(float) );
 
  // set boundary conditions
  laplace_init (A, n, m);
 
  if(rank == 0)
    printf("Jacobi relaxation Calculation: %d rows x %d columns mesh,"
         " maximum of %d iterations\n",
         n, m, iter_max );

  double start_time = MPI_Wtime();  // inicialitzem càlcul de temps
  // Main loop: iterate until error <= tol a maximum of iter_max iterations
  
  while ( error > tol && iter < iter_max ) {

  if (rank<world_size-1){  // si el rank és més petit que l'ultim proces
      MPI_Send(&A[m*final], m, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);  //envio fila al proces seguent
      MPI_Recv(&A[m*(final+1)], m, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status); //rebo fila del proces seguent
  }
 
  if (rank>0){ // si el rank es més gran que 0
      MPI_Send(&A[m*inicial], m, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);  //envio fila al proces anterior
      MPI_Recv(&A[m*(inicial-1)], m, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);  //rebo fila del proces anterior
  }
 

    error = 0.0f;  //inicialitzo error
   
// Compute new values using main matrix and writing into auxiliary matrix
  laplace_step(A+init, Anew+init, num_files, m);   
  
  //calculem l'error i fem la reducció
  error = laplace_error(A+init, Anew+init, num_files, m);
  MPI_Allreduce(&error, &total, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  error = total;
 
    // Copy from auxiliary matrix to main matrix
  laplace_copy(Anew+init, A+init, num_files, m);

    // if number of iterations is multiple of 10 then print error on the screen
    iter++;
   
    if (rank == 0 && iter % (iter_max/10) == 0)
       printf("%5d, %0.6f\n", iter, error);
  } // end while

    double end_time = MPI_Wtime();
      if (rank == 0)
    printf("Elapsed time: %f seconds\n", end_time - start_time);
  free(A);
  free(Anew);

  MPI_Finalize();
} 
