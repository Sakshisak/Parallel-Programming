//From MPI book

#include "mpi.h" 
#include <stdio.h> 
#include <string.h> 
#define BUFSIZE 10000 

int main(int argc, char *argv[]) {    

	int i, myrank, buf[BUFSIZE];
    	char filename[128];
    	FILE *myfile;

    	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 

    	strcpy(filename, "testfile");
	myfile = fopen(filename, "w");    

     	for (i=0; i<BUFSIZE; i++) {  
	   buf[i] = myrank + i;
	   fprintf(myfile, "%d\n", buf[i]);
	}   
	fclose(myfile);    

	MPI_Finalize();    
	return 0; 

} 
