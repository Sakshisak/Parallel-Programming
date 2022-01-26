#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[]) {

  int rank, numtasks, i;

  // Setup
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  if (argc != 2) {
    printf("usage: %s number_of_elements\n", argv[0]);
    return 1;
  }

  // Allocate message at the senders (may be different sized for different processes)
  int arrSize = atoi(argv[1])*(rank+1);  // synthetically creating different array sizes
  int message[arrSize];
  int countArray[numtasks], displArray[numtasks];

  int displ = 0;                       // note that root process is 0 here

  // this information is needed by the root
  if (!rank) 
   for (i = 0; i < numtasks; i++) {
    countArray[i] = arrSize*(i+1);         // depends on the counts, root may need to get it from the processes
    displArray[i] = displ;
    displ += countArray[i];
    printf ("%d %d %d\n", i, countArray[i], displArray[i]);
   }

  if (!rank) 
    printf ("\n");

  // every process initializes their local array 
  srand(time(NULL));
  for (i = 0; i < arrSize; i++) {
    message[i] = i; // (double)rand() / (double)RAND_MAX;
  }

  int recvMessage[displ];    // significant at the root process
  
  // receive different counts of elements from different processes
  MPI_Gatherv (message, arrSize, MPI_INT, recvMessage, countArray, displArray, MPI_INT, 0, MPI_COMM_WORLD);

  if (!rank) 
   for (i = 0; i < displ; i++) {
     printf ("%d %d\n", i, recvMessage[i]);
   }

  //printf ("\n");

  // finalize
  MPI_Finalize();

  return 0;
}

