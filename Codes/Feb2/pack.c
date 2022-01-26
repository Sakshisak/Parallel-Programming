// Timing codes

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, position=0, count; 
  MPI_Status status;
  double sTime, eTime, time;

  int M = atoi (argv[1]);
  int N = atoi (argv[2]);
  int array2D[M][N];
  int buffer[100];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // initialize data
  for (int i=0; i<M; i++)
   for (int j=0; j<N; j++)
    array2D[i][j] = myrank+i+j;

  sTime = MPI_Wtime();
  if (myrank == 0) {
   // pack the last element of every row (N ints) 
   for (int j=0; j<N; j++) {
     MPI_Pack (&array2D[j][M-1], 1, MPI_INT, buffer, 400, &position, MPI_COMM_WORLD);
     printf ("packed %d %d\n", j, position);
   }
   MPI_Send (buffer, position, MPI_PACKED, 1, 1, MPI_COMM_WORLD);	   
  }
  else {
   // receive N ints
   if (myrank == 1)
    MPI_Recv (buffer, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	   
   // verify
   MPI_Get_count (&status, MPI_INT, &count);
  } 
  eTime = MPI_Wtime();
  time = eTime - sTime;

  printf ("%lf\n", time);

  MPI_Finalize();
  return 0;

}

