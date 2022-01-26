//From MPI book

#include "mpi.h" 
#include <stdio.h> 
#include <string.h> 
#define BUFSIZE 100

int main(int argc, char *argv[]) {    

	int i, myrank, buf[BUFSIZE], rbuf[BUFSIZE];
    	char filename[128];
    	//FILE *myfile;
    	MPI_File myfile;

    	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    	strcpy(filename, "testfile");
	//myfile = fopen(filename, "w");    
	MPI_File_open (MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &myfile);

     	for (i=0; i<BUFSIZE; i++) {  
	   buf[i] = myrank + i;
	}  

	// File write - set process view
	MPI_File_set_view(myfile, myrank * BUFSIZE * sizeof(int), MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write (myfile, buf, BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);

	// File read - set process view
	MPI_File_set_view(myfile, myrank * BUFSIZE * sizeof(int), MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_read (myfile, rbuf, BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);

	MPI_File_close (&myfile);

     	for (i=0; i<BUFSIZE; i++) {  
	  if (buf[i] != rbuf[i]) printf ("%d %d %d\n", i, buf[i], rbuf[i]);
	}

	MPI_Finalize();    
	return 0; 

} 
