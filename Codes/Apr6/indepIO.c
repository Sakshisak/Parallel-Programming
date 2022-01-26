// Set explicit file offset

#include "mpi.h" 
#include <stdio.h> 
#include <string.h> 
#define BUFSIZE 4 

int main(int argc, char *argv[]) {    

  int i, myrank, nprocs, buf[BUFSIZE], rbuf[BUFSIZE], tsize;
  char filename[128];
  MPI_Status status;

  MPI_File fh;  // FILE

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  for (int i=0; i<BUFSIZE ; i++)
    buf[i]=i;

  strcpy(filename, "testfileIO");

  // File open, fh: individual file pointer
  MPI_File_open (MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);	//fopen

  MPI_Offset fo = (MPI_Offset) myrank*BUFSIZE*sizeof(int);

  // File write using explicit offset (independent I/O)
  MPI_File_write_at (fh, fo, buf, BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);	//fwrite

  // File read using explicit offset (independent I/O)
  MPI_File_read_at (fh, fo, rbuf, BUFSIZE, MPI_INT, &status);	//fread

  MPI_File_close (&fh);		//fclose

  for (i=0; i<BUFSIZE ; i++)
    if (buf[i] != rbuf[i]) printf ("Mismatch [%d] %d %d\n", i, buf[i], rbuf[i]);

  MPI_Finalize();    
  return 0; 

} 

