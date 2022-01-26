
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int N = atoi (argv[1]);
  int myrank, recvCount;
  double data[N]; 
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  //initialize data
  for (int i=0; i<N; i++)
    data[i] = (myrank+i)*1234.0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    MPI_Send(data, N, MPI_INT, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    MPI_Recv(data, N, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_DOUBLE, &recvCount);
  }

  printf ("%d %d %f\n", myrank, recvCount, data[1]);

  MPI_Finalize();
  return 0;

}

