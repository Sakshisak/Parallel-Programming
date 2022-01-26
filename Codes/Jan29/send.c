// Odd ranks receive from even ranks

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int numElements = atoi (argv[1]);
  int loopcount = 1000;
  int arr[numElements];
  int myrank, size, len;
  MPI_Status status;
  int count, recvarr[numElements];
  char hostname[64];

  if (argc>2) loopcount = atoi (argv[2]);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  // get processor name [note: this is not the core ID]
  MPI_Get_processor_name(hostname, &len);

  double time = MPI_Wtime();

  for (int loop=0; loop<loopcount; loop++) {

  if (myrank == 0)    
    MPI_Send(arr, numElements, MPI_BYTE, 1, 99, MPI_COMM_WORLD);
  else  
    MPI_Recv(recvarr, numElements, MPI_BYTE, 0, 99, MPI_COMM_WORLD, &status);

  }

  printf ("Rank %d (%s) %lf\n", myrank, hostname, MPI_Wtime()-time);

  MPI_Finalize();
  return 0;
}
