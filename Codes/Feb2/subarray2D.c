
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, recvcount;
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Datatype newtype;
  int ndims = 2;
  int array_of_sizes[] = {atoi(argv[1]), atoi(argv[2])};
  int array_of_subsizes[] = {2, 2}; 
  int array_of_starts[] = {1, 1}; 
  int data[array_of_sizes[0]][array_of_sizes[1]];

  MPI_Type_create_subarray (ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_C, MPI_INT, &newtype); 
  MPI_Type_commit (&newtype);

  //initialize data
  for (int i=0; i<array_of_sizes[0]; i++)
   for (int j=0; j<array_of_sizes[1]; j++)
    data[i][j]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<array_of_sizes[0]; i++)
     for (int j=0; j<array_of_sizes[1]; j++)
      data[i][j]=i+j;
    MPI_Send(data, 1, newtype, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    printf ("\n");
    MPI_Recv(data, 1, newtype, 0, 99, MPI_COMM_WORLD, &status);
    for (int i=0; i<array_of_sizes[0]; i++) {
     for (int j=0; j<array_of_sizes[1]; j++)
       printf ("%d ", data[i][j]);
     printf ("\n");
    }
    printf ("\n");
  }

  MPI_Type_free (&newtype);

  MPI_Finalize();
  return 0;

}

