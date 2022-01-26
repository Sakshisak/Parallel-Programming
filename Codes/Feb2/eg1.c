
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#define numElements 1000000

int main( int argc, char *argv[])
{
  int arr[numElements];
  int myrank, size;
  MPI_Status status, sstatus;
  MPI_Request request;
  int count, recvarr[numElements];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  MPI_Isend(arr, numElements, MPI_INT, 0, 99, MPI_COMM_SELF, &request);
  MPI_Recv(recvarr, numElements, MPI_INT, 0, 99, MPI_COMM_SELF, &status);
  MPI_Wait (&request, &sstatus);

  MPI_Get_count (&status, MPI_INT, &count);
  printf ("%d ints\n", count);

  MPI_Finalize();
  return 0;
}
