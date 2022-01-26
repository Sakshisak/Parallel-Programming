#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

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

// Global declaration of communicators. 
MPI_Comm n_comm;        // one node comm per machine
MPI_Comm nl_comm;       // comm for node leaders
MPI_Comm gl_comm;       // comm for group leaders

MPI_Group w_group;
MPI_Group n_group;
MPI_Group nl_group;
MPI_Group gl_group;

int node_ranks[128], node_leaders[100], group_leaders[6];
int n_ranks, n_leaders, g_leaders;

// Number of ranks in n_comm = n_ranks
// Number of ranks in node_leaders = n_leaders
// Number of ranks in group_leaders = g_leaders 
int rank, size;

double get_arrays(int root){
  int allnodes[128];
  int node_map[100] = {0}; 
  int group_map[6] = {0};
  int node = get_my_node();
  int group = get_my_group(node);
  n_ranks = 0; n_leaders = 0; g_leaders = 0;

  double stime = MPI_Wtime();
  // Use allgather to get every rank's node id
  MPI_Allgather(&node, 1, MPI_INT, allnodes, 1, MPI_INT, MPI_COMM_WORLD);
  for (int i=0; i<size; i++){
      int ptr = (i+root + size)%size;
      int grp = get_my_group(allnodes[ptr]);

      if (allnodes[ptr] == node)
          node_ranks[n_ranks++] = ptr;            // Other ranks belonging to the same node 

      if(grp == group){
          if (node_map[allnodes[ptr]] == 0){
              node_map[allnodes[ptr]] = 1;
              node_leaders[n_leaders++] = ptr;    // Other node leaders belonging to the same group
          }
      }
      if(group_map[grp] == 0){
          group_map[grp] = 1;
          group_leaders[g_leaders++] = ptr;       // All the group leaders 
      }
  }
  return MPI_Wtime()-stime;
}

double get_comms(int root, int* is_node_leader, int* is_group_leader){
    double stime = MPI_Wtime();
    *is_node_leader = 0;
    *is_group_leader = 0;

    MPI_Group w_group;
    MPI_Comm_group (MPI_COMM_WORLD, &w_group);

    MPI_Group n_group;
    MPI_Group_incl (w_group, n_ranks, node_ranks, &n_group);
    MPI_Comm_create_group(MPI_COMM_WORLD, n_group, n_ranks, &n_comm);           // Per-node communicator
    
    if (rank == node_ranks[0]){
        *is_node_leader = 1;
        MPI_Group nl_group;
        MPI_Group_incl (w_group, n_leaders, node_leaders, &nl_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, nl_group, n_leaders, &nl_comm);             // Node leader communicator
    }

    if (rank == node_leaders[0]){
        *is_group_leader = 1;
        MPI_Group g_group;
        MPI_Group_incl (w_group, g_leaders, group_leaders, &g_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, g_group, g_leaders, &gl_comm);              // Group leader communicator
    }
    return MPI_Wtime() - stime;
}

// Uses 3 level subcommunicators
int MPI_Reduce_optimised(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int is_node_leader, int is_group_leader){
    double* itmd_buff = (double*)malloc(count*sizeof(double*));
    double* itmd_buff_2 = (double*)malloc(count*sizeof(double*));

    MPI_Reduce(sendbuf, itmd_buff, count, datatype, op, 0, n_comm);
    if (is_node_leader){
        MPI_Reduce((const void *)itmd_buff, itmd_buff_2, count, datatype, op, 0, nl_comm);
    }

    if(is_group_leader){
        MPI_Reduce((const void *)itmd_buff_2, recvbuf, count, datatype, op, 0, gl_comm);
    }

    free(itmd_buff);
    return 0;
}

//function to give reduced value after operation op given sendbuf & itmd_buf 
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
  // other op not applicable for double datatype
  default:
    break;
  }
}

int Binomial_Reduce(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int* ranks, int size, MPI_Comm comm){
    int start = 0;
    
    for (int i=0; i<size; i++)
        if(ranks[i] == rank) start = i;

     //typecasting for general datatype

    if(start%2 == 0){   
      MPI_Request request[size];
      MPI_Status status[size];
      int count_r=0, it = 0;

      for(int i=1; i<size; i*=2){
        if(start%(i*2) == 0){
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
          MPI_Send(sendbuf, count, datatype, ranks[(start-i+size)%size], rank,comm);
          break;
        }
      }
    }

    return 0;
}

// Uses binomial reduce implementation
int MPI_ReduceB_optimised(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
   
    double itmd_buff[count];
    Binomial_Reduce(sendbuf, itmd_buff, count, datatype, op, node_ranks, n_ranks, comm);
    double itmd_buff2[count];
    if (node_ranks[0] == rank){
      Binomial_Reduce(itmd_buff,itmd_buff2, count, datatype, op, node_leaders, n_leaders, comm);
    }
    if(node_leaders[0] == rank){
      Binomial_Reduce(itmd_buff2,recvbuf, count, datatype, op, group_leaders, g_leaders, comm);
    }
    return 0;
}


int main(int argc, char** argv){
  int dsize = atoi(argv[1]);
  int root = 2;
  int is_node_leader=0;
  int is_group_leader=0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double stime, etime, ftime=0, reduce_std, reduce_opt, reduce_bopt;

  double* buff = (double*)malloc(dsize*sizeof(double));
  double* recv = (double*)malloc(dsize*sizeof(double));

  srand(time(0));
	for (int i=0; i<dsize; i++)
		buff[i] = rand() % size;

  double overhead_arrays = get_arrays(root);
  double overhead_comm = get_comms(root, &is_node_leader, &is_group_leader);

  for (int i=0; i<5; i++){
		stime = MPI_Wtime();
		MPI_Reduce(buff, recv, dsize, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		reduce_std += etime-stime;

		stime = MPI_Wtime();
		MPI_Reduce_optimised(buff, recv, dsize, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, is_node_leader, is_group_leader);
		etime = MPI_Wtime();
		reduce_opt += etime-stime;

    stime = MPI_Wtime();
		MPI_ReduceB_optimised(buff, recv, dsize, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		reduce_bopt += etime-stime;
	}

  reduce_bopt+=overhead_arrays;
  reduce_opt+=overhead_arrays+ overhead_comm;

  MPI_Reduce(&reduce_std, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)	
		printf("std_time = %lf\n", ftime);

	MPI_Reduce(&reduce_opt, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		printf("Reduce_opt = %lf\n", ftime);
  
  MPI_Reduce(&reduce_bopt, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		printf("Reduce_binomial_opt = %lf\n", ftime);

	MPI_Finalize();
	return 0;
}
  