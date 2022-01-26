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

  MPI_Get_processor_name (hostname, &len);

  // initialization
  srand(time(NULL));
  value = 2020212000000-rank*rank/N;
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
  
  // receive partial sums at rank 0
  stime = MPI_Wtime();
  if (rank)
  {
    MPI_Send(a+sidx, N/numtasks, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
  }
  else
  {
    int recvarr[numtasks][N/numtasks];
    for (int r=1; r<numtasks; r++)
    {
      MPI_Recv(recvarr[r], N/numtasks, MPI_DOUBLE, r, r, MPI_COMM_WORLD, &status);      
    }
  }
  etime = MPI_Wtime();
  cotime = etime - stime;

  // HOMEWORK
  // Add the partial sums
  

  printf ("%d: Time: %lf %lf\n", rank, ctime, cotime);

  // done with MPI
  MPI_Finalize();
}

