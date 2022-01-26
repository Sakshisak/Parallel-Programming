
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
  int ndims = 1;
  int array_of_sizes[] = {atoi (argv[1])};
  int array_of_subsizes[] = {atoi (argv[2])};
  int array_of_starts[] = {atoi (argv[3])};
  int N = array_of_sizes[0];
  int data[N]; 

  MPI_Type_create_subarray (ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_C, MPI_INT, &newtype);
  MPI_Type_commit (&newtype);

  //initialize data
  for (int i=0; i<N; i++)
    data[i]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<N; i++)
       data[i]=i;
    MPI_Send(data, 1, newtype, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    printf("\n");
    MPI_Recv(data, 1, newtype, 0, 99, MPI_COMM_WORLD, &status);
    for (int i=0; i<N; i++)
       printf ("%d ", data[i]);
    printf("\n\n");
  }

  MPI_Type_free (&newtype);

  MPI_Finalize();
  return 0;

}

