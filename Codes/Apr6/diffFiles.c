//From MPI book

#include "mpi.h" 
#include <stdio.h> 
#include <stdlib.h> 

int main(int argc, char *argv[]) {    

	int i, myrank;
    	char filename[128];
    	FILE *myfile;

	int BUFSIZE = atoi(argv[1]);
	int buf[BUFSIZE];

    	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 

    	sprintf(filename, "testfile.%d", myrank);    
	myfile = fopen(filename, "w");    

     	for (i=0; i<BUFSIZE; i++) {  
	   buf[i] = myrank + i;
	   fprintf(myfile, "%d ", buf[i]);
	}   
	fprintf(myfile, "\n");
	fclose(myfile);    

	MPI_Finalize();    
	return 0; 

} 
