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
  int rowmax, max, remains[2], comm1D_row, comm1D_col, val;
  dim[0] = 3 /* rows */, dim[1] = 3 /* columns */;
  wrap_around[0] = 0, wrap_around[1] = 0;
  reorder = 0;

  MPI_Cart_create (MPI_COMM_WORLD, ndims, dim, wrap_around, reorder, &comm2D);

  remains[0] = 0, remains[1] = 1;
  MPI_Cart_sub(comm2D, remains, &comm1D_row);
  // printf ("Rank %d, new rank %d, source %d dest %d\n", myrank, newrank, source, dest);
  remains[0] = 1, remains[1] = 0;
  MPI_Cart_sub(comm2D, remains, &comm1D_col);

  MPI_Reduce(&val, &rowmax, 1, MPI_INT, MPI_MAX, 0, comm1D_row);
  MPI_Reduce(&rowmax, &max, 1, MPI_INT, MPI_MAX, 0, comm1D_col);
  printf("max is : %d\n", max);
  MPI_Comm_free (&comm2D);

  MPI_Finalize();
  return 0;
}