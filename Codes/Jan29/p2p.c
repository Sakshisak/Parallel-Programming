
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>

int main( int argc, char *argv[])
{
  int myrank, size;
  double start_time, time, max_time;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Status status[2];
  MPI_Request request[2];

  int BUFSIZE = 1000;
  int arr[BUFSIZE];

  // send from all ranks to the last rank
  start_time = MPI_Wtime ();
  if (myrank < size-1) 
  {
    MPI_Bsend("good-1", strlen("good-1"), MPI_CHAR, size-1, 1, MPI_COMM_WORLD);
    MPI_Ssend("hello", strlen("hello"), MPI_CHAR, size-1, 2, MPI_COMM_WORLD);
  }
  else 
  {
    int count;
    char recvarr[size][BUFSIZE];
    MPI_Recv(recvarr[0], BUFSIZE, MPI_CHAR, 0, 2, MPI_COMM_WORLD, &status[0]);
    printf("Rank %d of %d received from rank %d with tag %d\n", myrank, size, status[0].MPI_SOURCE, status[0].MPI_TAG);
    MPI_Recv(recvarr[1], BUFSIZE, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status[1]);
    printf("Rank %d of %d received from rank %d with tag %d\n", myrank, size, status[1].MPI_SOURCE, status[1].MPI_TAG);
  }

  time = MPI_Wtime () - start_time;

  MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) printf ("Max time = %lf\n", max_time);

  MPI_Finalize();
  return 0;
}
