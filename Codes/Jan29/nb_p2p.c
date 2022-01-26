
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size;
  double start_time, time, max_time;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Status status[size-1];
  MPI_Request request[size-1];

  int BUFSIZE = atoi(argv[1]);
  int arr[BUFSIZE];

  // send from all ranks to the last rank
  start_time = MPI_Wtime ();
  if (myrank < size-1) 
  {
    MPI_Send(arr, BUFSIZE, MPI_INT, size-1, 99, MPI_COMM_WORLD);
  }
  else 
  {
    int count, recvarr[size][BUFSIZE];
    for (int i=0; i<size-1; i++)
      MPI_Irecv(recvarr[i], BUFSIZE, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request[i]);
    MPI_Waitall (size-1, request, status);
  }
  time = MPI_Wtime () - start_time;

  MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
  if (myrank == size-1) printf ("Max time = %lf\n", max_time);

  MPI_Finalize();
  return 0;
}
