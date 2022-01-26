
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int count = atoi(argv[1]);
  int myrank, buf1[count], buf2[count];
  MPI_Request req[2];
  MPI_Status status;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

  // Initialize data
  for (int i=0; i<count; i++)
    buf1[i] = myrank + i*i, buf2[i] = i;

  // Multiple nonblocking collective operations can be outstanding on a single communicator and match in order.

  if (myrank == 0) {
   MPI_Ibcast (buf1, 4, MPI_INT, 0, MPI_COMM_WORLD, &req[0]);
   MPI_Ibcast (buf2, 1, MPI_INT, 0, MPI_COMM_WORLD, &req[1]);
  }
  else if (myrank == 1) {
   MPI_Ibcast (buf1, 4, MPI_INT, 0, MPI_COMM_WORLD, &req[0]);
   MPI_Ibcast (buf2, 1, MPI_INT, 0, MPI_COMM_WORLD, &req[1]);
  }

  MPI_Waitall(2, req, &status);

  MPI_Finalize();
  return 0;

}


