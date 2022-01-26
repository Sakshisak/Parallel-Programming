// Parallel sum of array
// 2020-21-II
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

int main(int argc, char *argv[]) 
{
  int numtasks, rank, len, rc, i, sidx;
  double stime, etime, ctime, cotime; 
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;

  int N = atoi(argv[1]);
  double a[N], value, sum;

  // initialize MPI
  MPI_Init (&argc, &argv);

  // get number of tasks
  MPI_Comm_size (MPI_COMM_WORLD, &numtasks);

  // get my rank
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  // initialization
  value = rank; 
  sidx = rank*N/numtasks;
  for (i=sidx; i<sidx+N/numtasks ; i++)
      a[i] = value; 

  // local sum computation 
  sum=0.0;
  stime = MPI_Wtime();
  for (i=sidx; i<sidx+N/numtasks ; i++)
      sum += a[i] * a[i];
  etime = MPI_Wtime();
  ctime = etime - stime;
  
  double globalsum; 

  stime = MPI_Wtime();
  MPI_Reduce (&sum, &globalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  etime = MPI_Wtime();
  cotime = etime - stime;

  if (!rank)
  printf ("%lf %lf %lf\n", globalsum, ctime, cotime);

  // done with MPI
  MPI_Finalize();
}

