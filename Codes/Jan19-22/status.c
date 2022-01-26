
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#define numElements 1000

int main( int argc, char *argv[])
{
  int arr[numElements] = {0};
  int myrank, size;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  if (myrank == 0)    /* code for process zero */
  {
    MPI_Send(arr, 100, MPI_INT, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process one */
  {
    int count, recvarr[numElements];
    MPI_Recv(recvarr, numElements, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_INT, &count);
    printf("Rank %d of %d received %d elements from %d\n", myrank, size, count, status.MPI_SOURCE);

  }

  MPI_Finalize();
  return 0;
}
