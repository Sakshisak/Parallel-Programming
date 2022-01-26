
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, recvcount;
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Datatype newvtype;
  int N = atoi (argv[1]);
  int column = atoi (argv[2]);
  int count = atoi (argv[3]);
  int blocklen = atoi (argv[4]);
  int stride = atoi (argv[5]);
  int data[N][N], received[N*blocklen]; 

  MPI_Type_vector (count, blocklen, stride, MPI_INT, &newvtype);
  MPI_Type_commit (&newvtype);

  //initialize data
  for (int i=0; i<N; i++)
   for (int j=0; j<N; j++)
    data[i][j]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<N; i++)
     for (int j=0; j<N; j++)
      data[i][j]=column+i+j;
    MPI_Send(&data[0][column], 1, newvtype, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    printf ("\n");
    MPI_Recv(received, count*blocklen, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
    for (int i=0; i<count*blocklen; i++)
       printf ("%d ", received[i]);
    printf ("\n\n");
  }

  MPI_Type_free (&newvtype);

  MPI_Finalize();
  return 0;

}

