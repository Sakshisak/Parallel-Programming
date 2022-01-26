
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size;
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Datatype newtype;
  int N = atoi (argv[1]);
  int numElements = atoi (argv[2]);
  int count = 3;
  int blocklengths[] = {1, 2, 3};
  int displacements[] = {0, 3, 6};
  int data[N]; 

  MPI_Type_indexed (count, blocklengths, displacements, MPI_INT, &newtype);
  MPI_Type_commit (&newtype);

  //initialize data
  for (int i=0; i<N; i++)
    data[i]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<N; i++)
       data[i]=i;
    MPI_Send(data, numElements, newtype, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    printf("\n");
    MPI_Recv(data, numElements, newtype, 0, 99, MPI_COMM_WORLD, &status);
    for (int i=0; i<N; i++)
       printf ("%d ", data[i]);
    printf("\n\n");
  }

  MPI_Type_free (&newtype);

  MPI_Finalize();
  return 0;

}

