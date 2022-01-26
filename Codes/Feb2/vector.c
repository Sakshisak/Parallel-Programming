
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
  int count = atoi (argv[2]);
  int blocklen = atoi (argv[3]);
  int stride = atoi (argv[4]);
  int numVectors = atoi (argv[5]);
  int data[N]; 
  int recvcount, recvcountVector;
  MPI_Type_vector (count, blocklen, stride, MPI_INT, &newvtype);
  MPI_Type_commit (&newvtype);

  //initialize data
  for (int i=0; i<N; i++)
    data[i]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<N; i++)
       data[i]=i;
    MPI_Send(data, numVectors, newvtype, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    printf("\n");
    MPI_Recv(data, numVectors, newvtype, 0, 99, MPI_COMM_WORLD, &status);
    MMPI_Get_count (&status, newvtype, &recvcountVector);
    MPI_Get_count (&status, MPI_INT, &recvcount);
    printf ("Numbers: %d %d\n", recvcountVector, recvcount);


    MPI_Get_elements (&status, newvtype, &recvcountVector);
    MPI_Get_elements (&status, MPI_INT, &recvcount);
    printf ("Numbers: %d %d\n", recvcountVector, recvcount);
    for (int i=0; i<N; i++)
       printf ("%d ", data[i]);
    printf("\nnewvtype count = %d\n",count);
  }

  MPI_Type_free (&newvtype);

  MPI_Finalize();
  return 0;

}

