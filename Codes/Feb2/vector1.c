
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, recvcount;
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Datatype newvtype1, newvtype2, newvtype3;
  int N = atoi (argv[1]);
  int count = atoi (argv[2]);
  int blocklen = atoi (argv[3]);
  int stride = atoi (argv[4]);
  int numVectors = atoi (argv[5]);
  int data[N];

  MPI_Type_vector (count, blocklen, stride, MPI_INT, &newvtype1);
  MPI_Type_commit (&newvtype1);
  MPI_Type_vector (count, blocklen, stride+1, MPI_INT, &newvtype2);
  MPI_Type_commit (&newvtype2);
  MPI_Type_vector (count, blocklen+1, stride+1, MPI_INT, &newvtype3);
  MPI_Type_commit (&newvtype3);

  //initialize data
  for (int i=0; i<N; i++)
    data[i]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<N; i++)
       data[i]=i;
    MPI_Send(data, numVectors, newvtype1, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    int count3;
    printf("\n");
    MPI_Recv(data, numVectors/2 + 1, newvtype3, 0, 99, MPI_COMM_WORLD, &status);
    for (int i=0; i<N; i++)
       printf ("%d ", data[i]);
    MPI_Get_count (&status, newvtype3 , &count3);
    printf("\nnewvtype3 count = %d \n", count3);
  }

  MPI_Type_free (&newvtype1);
  MPI_Type_free (&newvtype2);

  MPI_Finalize();
  return 0;

}

