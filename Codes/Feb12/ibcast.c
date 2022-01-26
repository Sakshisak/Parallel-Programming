
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int count = atoi(argv[1]);
  int myrank, buf[count]; 
  double buf_can_modify[count];
  MPI_Status status;
  MPI_Request request;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  // initialize data
  for (int i=0; i<count; i++)
      buf[i] = myrank + i*i;

  // has to be called by all processes
  MPI_Ibcast(buf, count, MPI_INT, 1, MPI_COMM_WORLD, &request);
  
  // some random expensive computations 
  for (int j=0; j<count; j++)
    for (int i=0; i<count; i++)
      buf_can_modify[i] = i*12022021/(myrank+1)*pow(2,i);

  MPI_Wait (&request, &status);

  MPI_Finalize();
  return 0;

}

