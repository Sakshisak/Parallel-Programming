
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
  int N = atoi (argv[1]);
  int data[N][N]; 
  int count = N; 
  int blocklengths[N], displacements[N];

  for (int i=0; i<N; i++)
    blocklengths[i] = i+1, displacements[i] = i*N;
  MPI_Type_indexed (count, blocklengths, displacements, MPI_INT, &newtype);
  MPI_Type_commit (&newtype);

  //initialize data
  for (int i=0; i<N; i++)
   for (int j=0; j<N; j++)
    data[i][j]=0;
  
  if (myrank == 0)    /* code for process 0 */
  {
    for (int i=0; i<N; i++)
     for (int j=0; j<N; j++)
      data[i][j]=i+j;
    MPI_Send(data, 1, newtype, 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank == 1)  /* code for process 1 */
  {
    printf ("\n");
    MPI_Recv(data, 1, newtype, 0, 99, MPI_COMM_WORLD, &status);
    for (int i=0; i<N; i++) {
     for (int j=0; j<N; j++)
       printf ("%d ", data[i][j]);
     printf ("\n");
    }
    printf ("\n");
  }

  MPI_Type_free (&newtype);

  MPI_Finalize();
  return 0;

}

