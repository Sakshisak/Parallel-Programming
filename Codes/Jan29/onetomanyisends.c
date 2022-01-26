// Timing codes

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size; 
  double sTime, eTime, time, maxTime;

  MPI_Init(&argc, &argv);

  int count = atoi (argv[1]);
  int buf[count];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request request[size];
  MPI_Status status[size];

  // initialize data
  for (int i=0; i<count; i++)
   buf[i] = myrank+i;

  sTime = MPI_Wtime();
  if (myrank == 0) {
   for (int r=1; r<size; r++)
    MPI_Isend (buf, count, MPI_INT, r, r, MPI_COMM_WORLD, &request[r-1]);
   MPI_Waitall (size-1, request, status);
  }
  else
   MPI_Recv (buf, count, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status[0]);	   
  eTime = MPI_Wtime();
  time = eTime - sTime;

  // obtain max time
  MPI_Reduce (&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) printf ("%lf\n", maxTime);

  MPI_Finalize();
  return 0;

}

