#include "mpi.h"
#include<string.h>
#include<stdlib.h>
#include<stdio.h>

#define numElements 1000
int main( int argc, char *argv[])
{
int arr[numElements] = {0};
int myrank;
int size;
MPI_Status status;
double sTime, eTime, time;

MPI_Init( &argc, &argv );
int count = atoi (argv[1]);
int buf[count];

MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
MPI_Comm_size( MPI_COMM_WORLD, &size);

for(int i=0; i<count; i++){
  buf[i] = myrank+i;
}

sTime = MPI_Wtime();

if (myrank == 0) /* code for process zero */
{
  MPI_Send(buf, count, MPI_INT, 1, 1, MPI_COMM_WORLD);
}
else if (myrank == 1) /* code for process one */
{
  MPI_Recv(buf, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
}

eTime = MPI_Wtime();
printf("%lf\n", time);

MPI_Finalize();
return 0;
}