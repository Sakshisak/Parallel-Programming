
#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  
  int myrank, size;
  int ndims=2, dim[2], wrap_around[2], reorder;
  int newrank, newsize, source, dest;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm comm2D;

  dim[0] = 3 /* rows */, dim[1] = 4 /* columns */;
  wrap_around[0] = 0, wrap_around[1] = 0;
  reorder = 0;

  MPI_Cart_create (MPI_COMM_WORLD, ndims, dim, wrap_around, reorder, &comm2D);

  MPI_Comm_rank (comm2D, &newrank);
  MPI_Comm_size (comm2D, &newsize);

  MPI_Cart_shift (comm2D, 0, -1, &source, &dest);
  printf ("Rank %d, new rank %d, source %d dest %d\n", myrank, newrank, source, dest);

  MPI_Comm_free (&comm2D);

  MPI_Finalize();
  return 0;
}
