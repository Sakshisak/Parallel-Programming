// Timing codes

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, flag; 
  double sTime, eTime, time;
  MPI_Status status;
  MPI_Request request;

  MPI_Init(&argc, &argv);

  int count = atoi (argv[1]);
  int buf[count];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // initialize data
  for (int i=0; i<count; i++)
   buf[i] = myrank+i;

  sTime = MPI_Wtime();
  if (myrank == 0)
   MPI_Isend (buf, count, MPI_INT, 1, 1, MPI_COMM_WORLD, &request); 
  eTime = MPI_Wtime();
  time = eTime - sTime;

  if (!myrank) {
     MPI_Test (&request, &flag, &status);
     printf ("%lf %d\n", time, flag);
  }

  MPI_Finalize();
  return 0;

}

