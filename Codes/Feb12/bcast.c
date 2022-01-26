#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int count = atoi(argv[1]);
  int myrank, buf[count]; 
  MPI_Status status;
  MPI_Request request;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  for (int i=0; i<count; i++)
      buf[i] = myrank + i*i;
  if(myrank == 0){
  MPI_Bcast(buf, count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(buf, count-1, MPI_INT, 1, MPI_COMM_WORLD);
  }
  else{
  MPI_Bcast(buf, count, MPI_INT, 1, MPI_COMM_WORLD);
  MPI_Bcast(buf, count-1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;

}