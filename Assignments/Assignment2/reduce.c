#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include "collective.h"

int get_my_node(){
	char host_name[10];
	int name_len, node = 0;
	MPI_Get_processor_name(host_name, &name_len);
	for (int i=5; i<name_len; i++) 
		node = node*10 + (int)(host_name[i] - '0');
	return node;
}

int get_my_group(int node){
	// outliers 
	if (node == 13) return 2;
	else if (node == 31) return 1;
	else if (node == 45) return 4;
	// scan through groups	
	else if (node <= 16) return 1;
	else if (node <= 32) return 2;
	else if (node <= 46) return 3;
	else if (node <= 61) return 4;
	else if (node <= 78) return 5;
	else if (node <= 92) return 6;
	else return -1;
}

void Operation(MPI_Op op, double *itmd_buf, double *sendbuf, int count){
  switch (op)
  {
  case MPI_SUM: for(int i=0; i<count; i++) sendbuf[i] += itmd_buf[i];
                break;
  case MPI_MAX: for(int i=0; i<count; i++) sendbuf[i] = sendbuf[i]>itmd_buf[i]?sendbuf[i]:itmd_buf[i];
                break;
  case MPI_MIN: for(int i=0; i<count; i++) sendbuf[i] = sendbuf[i]<itmd_buf[i]?sendbuf[i]:itmd_buf[i];
                break;
  case MPI_PROD: for(int i=0; i<count; i++) sendbuf[i] *= itmd_buf[i];
                break;
  /*case MPI_LAND: for(int i=0; i<count; i++) sendbuf[i] = sendbuf[i] && itmd_buf[i];
                break;             
  case MPI_BAND: for(int i=0; i<count; i++) sendbuf[i] &= itmd_buf[i];
                break;
  case MPI_LOR: for(int i=0; i<count; i++) sendbuf[i] = sendbuf[i] || itmd_buf[i];
                break;             
  case MPI_BOR: for(int i=0; i<count; i++) sendbuf[i] |= itmd_buf[i];
                break; 
  case MPI_LXOR: for(int i=0; i<count; i++) sendbuf[i] = !sendbuf[i] != !itmd_buf[i];
                break;
  case MPI_BXOR: for(int i=0; i<count; i++) sendbuf[i] = sendbuf[i] ^ itmd_buf[i];
                break;
  case MPI_MAXLOC: for(int i=0; i<count; i++) sendbuf[i]&=itmd_buf[i];
                break;
  case MPI_MINLOC: for(int i=0; i<count; i++) sendbuf[i]&=itmd_buf[i];
                break;*/
  default:
    break;
  }
}

int Binomial_Reduce(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int* ranks, int size, MPI_Comm comm){
    int rank, start=0;
    MPI_Comm_rank(comm, &rank);

    for (int i=0; i<size; i++)
        if(ranks[i] == rank) start = i;

     //typecasting for general datatype

    if(start%2 == 0){   
      MPI_Request request[size];
      MPI_Status status[size];
      int count_r=0, it = 0;

      for(int i=1; i<size; i*=2){
        if(start%(i*2) == 0){
          //printf("i = %d, rank = %d\n",i, rank);
          double itmd_buff[count];
          MPI_Recv(itmd_buff, count, datatype, ranks[(start+i)%size], MPI_ANY_TAG, comm, &status[0]);
          Operation(op, itmd_buff, sendbuf, count);
          if(i == size/2){
            for(int j=0; j<count; j++) recvbuf[j] = itmd_buff[j];
          }
        }
      }
    }

    if(rank != ranks[0]){
      MPI_Request request;
      MPI_Status status;
      for(int i=size/2; i>0 ; i/=2){
        int send_rank;
        if(start%i == 0){
          //send_rank = ranks[start-(int)pow(2,i)];
          MPI_Send(sendbuf, count, datatype, ranks[(start-i+size)%size], rank,comm);
          printf("Send s = %d\n", rank);
          break;
          
        }
      }
    }

    return 0;
}


int MPI_Reduce_optimised(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
    int group, node;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    int allnodes[size], node_leaders[size], node_ranks[size];
    int n_leaders = 0, n_ranks = 0;
    int node_map[100];

    for (int i=0; i<100; i++)
        node_map[i] = 0;
    node = get_my_node();
    MPI_Allgather(&node, 1, MPI_INT, allnodes, 1, MPI_INT, comm);

    for (int i=0; i<size; i++){
        int ptr = (i+root + size)%size;
        if (allnodes[ptr] == node)
            node_ranks[n_ranks++] = ptr;

        if (node_map[allnodes[ptr]] == 0){
            node_map[allnodes[ptr]] = 1;
            node_leaders[n_leaders++] = ptr;
        }
        //if (rank==0) printf("ptr: %d\n", ptr);
    }
/*
    if (node_ranks[0] == rank){
        printf("node leaders: %d\n", rank);
        for (int i=0; i<n_leaders; i++)
                printf("%d ", node_leaders[i]);
        printf("\n");
   }

    printf("node: %d, commsize: %d\n", rank, n_ranks); 
    for(int i=0; i<n_ranks; i++)
        printf("%d ", node_ranks[i]);
    printf("\n");

*/
    double* itmd_buff = (double*)malloc(count*(sizeof(double)));
    Binomial_Reduce(sendbuf, itmd_buff, count, datatype, op, node_ranks, n_ranks, comm);

    if (node_ranks[0] == rank){
        for(int j=0; j<count; j++) sendbuf[j] = itmd_buff[j];
        Binomial_Reduce(sendbuf,recvbuf, count, datatype, op, node_leaders, n_leaders, comm);
    }
    free(itmd_buff);
    return 0;
}

int main(int argc, char** argv){
	int rank, size;
	double stime, etime, ftime;
	double opt_time=0, std_time=0;
	int dsize = atoi(argv[1]);
	double* buff = (double*)malloc(dsize*sizeof(double));
	double* recv = (double*)malloc(dsize*sizeof(double));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	srand(time(0));
	
	for (int i=0; i<dsize; i++)
		buff[i] = rand() % size;

	// int group = get_my_group();
	// printf("rank: %d, group: %d\n", rank, group);
	for (int i=0; i<5; i++){
		stime = MPI_Wtime();
		MPI_Reduce_optimised(buff, recv, dsize, MPI_DOUBLE, MPI_MAX, 2, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		opt_time += etime-stime;

		stime = MPI_Wtime();
		MPI_Reduce(buff, recv, dsize, MPI_DOUBLE, MPI_MAX, 2, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		std_time += etime-stime;
	}

	MPI_Reduce(&opt_time, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)	
		printf("Opt_time = %lf\n", ftime);

	MPI_Reduce(&std_time, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		printf("STD_time = %lf\n", ftime);

	MPI_Finalize();
	return 0;
}

