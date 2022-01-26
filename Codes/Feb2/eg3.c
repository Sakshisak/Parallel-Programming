
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int N = atoi (argv[1]);
  int myrank, count;
  double data[N];
  int recvdata[N]; 
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  if (myrank == 0)       /* code for process 0 */
  {
    MPI_Send(data, N, MPI_DOUBLE, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    MPI_Recv(recvdata, N, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
    MPI_Get_count( &status, MPI_INT, &count );
    printf ("%d\n", count);
  }

  MPI_Finalize();
  return 0;

}

