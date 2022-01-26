#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

// Global declaration of communicators. 
MPI_Comm n_comm;        // one node comm per machine
MPI_Comm nl_comm;       // comm for node leaders
MPI_Comm gl_comm;       // comm for group leaders

MPI_Group w_group;
MPI_Group n_group;
MPI_Group nl_group;
MPI_Group gl_group;

int node_ranks[128], node_leaders[100], group_leaders[6];

// Number of ranks in n_comm = n_ranks
// Number of ranks in node_leaders = n_leaders
// Number of ranks in group_leaders = g_leaders 
int n_leaders = 0, n_ranks = 0, g_leaders = 0;

int rank, size;
//1 group: 8 nodes
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
	if (node == 13) return 1;
	else if (node == 31) return 0;
	else if (node == 45) return 3;
	// scan through groups	
	else if (node <= 16) return 0;
	else if (node <= 32) return 1;
	else if (node <= 46) return 2;
	else if (node <= 61) return 3;
	else if (node <= 78) return 4;
	else if (node <= 92) return 5;
	else return -1;
}

double get_arrays(int root){
  int allnodes[128];                            // node numbers for all ranks
  int node_map[100] = {0};                      // check if a node leader has been allocated or not
  int group_map[6] = {0};                       // check if a group leader has been allocated or not
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

    if (g_leaders>1 && rank == node_leaders[0]){
        *is_group_leader = 1;
        MPI_Group g_group;
        MPI_Group_incl (w_group, g_leaders, group_leaders, &g_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, g_group, g_leaders, &gl_comm);              // Group leader communicator
    }
    return MPI_Wtime() - stime;
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
  default:
    break;
  }
}

// binomial reduce implementation 
int Binomial_Reduce(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int* ranks, int size, MPI_Comm comm){
    int start = 0;
    for (int i=0; i<size; i++)
        if(ranks[i] == rank) start = i;

    if(start%2 == 0){   
      MPI_Request request[size];
      MPI_Status status[size];
      int count_r=0, it = 0;

      for(int i=1; i<size; i*=2){
        if(start%(i*2) == 0){
          double itmd_buff[count];
          // recieve from children in binomial tree
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
          // pass on to parent in binomial tree
          MPI_Send(sendbuf, count, datatype, ranks[(start-i+size)%size], rank,comm);
          break;
        }
      }
    }

    return 0;
}

// Uses binomial reduce implementataion
int MPI_Reduce_optimised(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
    
    double itmd_buff[count];
    Binomial_Reduce(sendbuf, itmd_buff, count, datatype, op, node_ranks, n_ranks, comm);
    double itmd_buff2[count];
    if (node_ranks[0] == rank){
      if(g_leaders > 1) Binomial_Reduce(itmd_buff,itmd_buff2, count, datatype, op, node_leaders, n_leaders, comm);
      else if (g_leaders == 1) Binomial_Reduce(itmd_buff,recvbuf, count, datatype, op, node_leaders, n_leaders, comm);
    }
    if(g_leaders > 1 && node_leaders[0] == rank){
      Binomial_Reduce(itmd_buff2,recvbuf, count, datatype, op, group_leaders, g_leaders, comm);
    }
    return 0;
}

// Optimization using 3 level subcommunicators
int MPI_Bcast_optimised(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, int is_node_leader, int is_group_leader){
    if (g_leaders > 1 && is_group_leader){
        MPI_Bcast(buffer, count, datatype, 0, gl_comm);
    }
    if (is_node_leader){
        MPI_Bcast(buffer, count, datatype, 0, nl_comm);
    }
    MPI_Bcast(buffer, count, datatype, 0, n_comm);
    return 0;
}

// Optimization using 3 level subcommunicators
int MPI_Gather_optimised(const void *sendbuf, int scount, MPI_Datatype stype, void *recvbuf, int rcount, MPI_Datatype rtype, int root, int* allranks, int is_node_leader, int is_group_leader){
    int stype_size, rtype_size;
    MPI_Type_size(stype, &stype_size);
    MPI_Type_size(rtype, &rtype_size);
    void* itmd_buff = (void*)malloc(rcount*n_ranks*rtype_size);
    void* itmd_buff_2;
    
    MPI_Gather(sendbuf, scount, stype, itmd_buff, scount, stype, 0, n_comm);
    if (is_node_leader){
        itmd_buff_2 = (void*)malloc(scount*n_ranks*n_leaders*stype_size);
        MPI_Gather(itmd_buff, scount*n_ranks, stype, itmd_buff_2, scount*n_ranks, stype, 0, nl_comm);        
    }
    free(itmd_buff);

    if (g_leaders > 1 && is_group_leader){
        itmd_buff = (void*)malloc(scount*size*rtype_size);
        MPI_Gather(itmd_buff_2, scount*n_ranks*n_leaders, stype, itmd_buff, scount*n_ranks*n_leaders, stype, 0, gl_comm);

        if (rank == root){
            for(int i=0; i<size; i++){
                int id = allranks[i];
                memcpy((void*)(recvbuf + id*rcount*rtype_size), (void*)(itmd_buff + i*rcount*rtype_size), rcount*rtype_size);
            }
        }
        free(itmd_buff);
    }
    if(g_leaders == 1 && rank  == root){
      for(int i=0; i<size; i++){
        int id = allranks[i];
        memcpy((void*)(recvbuf + id*rcount*rtype_size), (void*)(itmd_buff_2 + i*rcount*rtype_size), rcount*rtype_size);
      }
    }
    if (is_node_leader) free(itmd_buff_2);
    return 0;
}

int MPI_Alltoallv_optimised( const void* sendbuff,
                             const int* send_counts, 
                             const int* sdispls, 
                             MPI_Datatype stype, 
                             void* recvbuf, 
                             const int* recv_counts, 
                             const int* rdispls, 
                             MPI_Datatype rtype
){
    // every process collects data to be sent 
    int ppn = n_ranks;
    int leader = node_ranks[0];
    int n_id = 0;
    int ssize;
    MPI_Type_size(stype, &ssize);

    int local_scount = 0;
    for(int i=0; i<size; i++){ 
        local_scount += send_counts[i];
    }

    int local_scounts[n_ranks];
    int all_node_scounts[n_ranks*size];
    int node_displs[n_ranks];
    int node_scount = 0;
    int displ=0;
    void* t_recv_buff;

    // how much data to recieve
    MPI_Gather(&local_scount, 1, MPI_INT, local_scounts, 1, MPI_INT, 0, n_comm);
    if (rank == leader){
        for(int i=0; i<n_ranks; i++){
            node_scount+=local_scounts[i];
            node_displs[i] = displ;
            displ+=local_scounts[i];
        }
        t_recv_buff = malloc(node_scount*ssize);
    }   
    // node leader recieves data from all procs in node_ranks
    // all_node_scounts is metadata to rearrange recieved data
    MPI_Gatherv(sendbuff, local_scount, stype, t_recv_buff, local_scounts, node_displs, stype, 0, n_comm);
    MPI_Gather(send_counts, size, MPI_INT, all_node_scounts, size, MPI_INT, 0, n_comm);

    void* t_send_buff;
    int send_sizes[n_ranks];
    int send_displs[n_ranks];
    if (rank == leader){
        t_send_buff = malloc(node_scount*ssize);
        int nl_send_counts[size];                           // analogous to send_counts
        int all_leader_scounts[size];                   // metadata to rearrange second recieve by nodeleaders
        int it = 0;
        //node leader rearranges data first time
        for(int i=0; i<size; i++){
            nl_send_counts[i] = 0;
            for(int j=0; j<n_ranks; j++){
                nl_send_counts[i] += all_node_scounts[j*size + i];
                memcpy((void*)(t_send_buff + it), (void*)(t_recv_buff + node_displs[j]*ssize), all_node_scounts[j*size + i]*ssize);
                it += all_node_scounts[j*size + i]*ssize;
            }
        }

        int nl_sdispls[n_leaders];
        int nl_rdispls[n_leaders];
        int leader_scounts[n_leaders];
        int leader_rcounts[n_leaders];
        int total_recv_count = 0;
        // calculate displacements, send counts
        displ = 0;
        for (int i=0; i<n_leaders; i++){
            nl_sdispls[i] = displ;
            int sum=0;
            for(int j=0; j<n_ranks; j++){
                sum+=nl_send_counts[j];
            }
            displ+=sum;
            leader_scounts[i] = sum;
        }

        free(t_recv_buff);
        // how much data to expect from other node leaders 
        MPI_Alltoall(leader_scounts, 1, MPI_INT, leader_rcounts, 1, MPI_INT, nl_comm);
        displ = 0;
        for(int i=0; i<n_leaders; i++){
            total_recv_count += leader_rcounts[i];
            nl_rdispls[i] = displ;
            displ+=leader_rcounts[i];
        }

        t_recv_buff = malloc(total_recv_count*ssize);
        MPI_Alltoallv(t_send_buff, leader_scounts, nl_sdispls, stype, t_recv_buff, leader_rcounts, nl_rdispls, stype, nl_comm);
        //metadata to rearrange recieved data
        MPI_Alltoall(nl_send_counts, n_ranks, MPI_INT, all_leader_scounts, n_ranks, MPI_INT, nl_comm);
        free(t_send_buff);

        // second rearrangement
        it=0;
        t_send_buff = malloc(total_recv_count*ssize);
        for (int i=0; i<n_ranks; i++){
            for (int j=0; j<n_leaders; j++){
                memcpy((void*)(t_send_buff + it), (void*)(t_recv_buff + nl_rdispls[j]*ssize), all_leader_scounts[j*n_ranks + i]*ssize);
                it+=all_leader_scounts[j*n_ranks + i]*ssize;
            }
        }
        free(t_recv_buff);
        // how much to send to every rank
        int sum = 0;
        for (int i=0; i<n_ranks; i++){
            send_sizes[i] = 0;
            for (int j=0; j<n_leaders; j++){
                send_sizes[i]+=all_leader_scounts[j*n_ranks + i];
            }
            send_displs[i] = sum;
            sum += send_sizes[i];
        }   
    }

    // Scatter
    int total_recv_size=0;
    for(int i=0; i<size; i++) 
        total_recv_size+= recv_counts[i];
    void* itmd_buff = malloc(total_recv_size*ssize);
    MPI_Scatterv(t_send_buff, send_sizes, send_displs, stype, itmd_buff, total_recv_size, stype, 0, n_comm);
    int ptr = 0;
    for(int i=0; i<size; i++){
        memcpy((void*)(recvbuf + rdispls[i]*ssize), (void*)(itmd_buff + ptr), recv_counts[i]*ssize);
        ptr += recv_counts[i]*ssize;
    }
    free(itmd_buff);
    if (rank == leader)
        free(t_send_buff);
    return 0;
}

int main(int argc, char** argv){
    int dsize = atoi(argv[1]);
    int root = 2;                       // root process.
    int is_node_leader=0;
    int is_group_leader=0;

    double stime, etime, ftime=0;
    double overhead_array=0, overhead_comm=0;
    double bcast_opt=0, bcast_std=0;
    double reduce_opt=0, reduce_std=0;
    double gather_opt=0, gather_std=0;
    double alltoall_opt=0, alltoall_std=0;

    // Create buffers 
    double* buff = (double*)malloc(dsize*sizeof(double));
    double* recv = (double*)malloc(dsize*sizeof(double));
	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialise send buffer
	srand(time(0));
	for (int i=0; i<dsize; i++)
		buff[i] = rand() % size;

    // Initialising communicators
    overhead_array = get_arrays(root);
    overhead_comm = get_comms(root, &is_node_leader, &is_group_leader);

    // printf("============== Bcast =============\n");

    for (int i=0; i<5; i++){
        stime = MPI_Wtime();
        if (rank == 2)
            MPI_Bcast(buff, dsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
        else 
            MPI_Bcast(recv, dsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
        etime = MPI_Wtime();
        bcast_std += etime-stime;

        stime = MPI_Wtime();
        if (rank == 2)
            MPI_Bcast_optimised(buff, dsize, MPI_DOUBLE, root, MPI_COMM_WORLD, is_node_leader, is_group_leader);
        else 
            MPI_Bcast_optimised(recv, dsize, MPI_DOUBLE, root, MPI_COMM_WORLD, is_node_leader, is_group_leader);
        etime = MPI_Wtime();
        bcast_opt += etime-stime;
    } 

    // // printf("============== Reduce =============\n");
    for (int i=0; i<5; i++){
		stime = MPI_Wtime();
		MPI_Reduce(buff, recv, dsize, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		reduce_std += etime-stime;

		stime = MPI_Wtime();
		MPI_Reduce_optimised(buff, recv, dsize, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		reduce_opt += etime-stime;
	}

    // printf("============== Gather =============\n");
    double gather_overhead = 0;
    stime = MPI_Wtime();
    int group_ranks[n_ranks*n_leaders];
    int allranks[size];
    double* recvgather = (double*)malloc(dsize*size*sizeof(double));
    if (is_node_leader){
        if(g_leaders > 1) MPI_Gather(node_ranks, n_ranks, MPI_INT, group_ranks, n_ranks, MPI_INT, 0, nl_comm);
        else if(g_leaders == 1) MPI_Gather(node_ranks, n_ranks, MPI_INT, allranks, n_ranks, MPI_INT, 0, nl_comm);
    }
    if (g_leaders > 1 && is_group_leader){
        MPI_Gather(group_ranks, n_ranks*n_leaders, MPI_INT, allranks, n_ranks*n_leaders, MPI_INT, 0, gl_comm);
    }
    gather_overhead = MPI_Wtime() - stime;

    for (int i=0; i<5; i++){
		stime = MPI_Wtime();
		MPI_Gather(buff, dsize, MPI_DOUBLE, recvgather, dsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		gather_std += etime-stime;

    stime = MPI_Wtime();
		MPI_Gather_optimised(buff, dsize, MPI_DOUBLE, recvgather, dsize, MPI_DOUBLE, root, allranks, is_node_leader, is_group_leader);
		etime = MPI_Wtime();
		gather_opt += etime-stime;
    }
    free(recvgather);
    int sendcounts[size], recvcounts[size], senddispls[size], recvdispls[size];
    for(int i=0; i<size; i++){
        sendcounts[i] = dsize/size;
        recvcounts[i] = dsize/size;
        senddispls[i] = (i*dsize)/size;
        recvdispls[i] = (i*dsize)/size;
    }

    // printf("============== Alltoall =============\n");
    for (int i=0; i<5; i++){
		stime = MPI_Wtime();
        MPI_Alltoallv(buff, sendcounts, senddispls, MPI_DOUBLE, recv, recvcounts, recvdispls, MPI_DOUBLE, MPI_COMM_WORLD);
		etime = MPI_Wtime();
		alltoall_std += etime-stime;

		stime = MPI_Wtime();
        MPI_Alltoallv_optimised(buff, sendcounts, senddispls, MPI_DOUBLE, recv, recvcounts, recvdispls, MPI_DOUBLE);
		etime = MPI_Wtime();
		alltoall_opt += etime-stime;
	  }

    if (is_node_leader){
        MPI_Comm_free(&nl_comm);
    }
    if (is_group_leader){
        MPI_Comm_free(&gl_comm);
    }
    MPI_Comm_free(&n_comm);

    bcast_opt += overhead_array + overhead_comm;
    reduce_opt += overhead_array;
    gather_opt += overhead_array + overhead_comm + gather_overhead;
    alltoall_opt += overhead_comm+overhead_array;

    //Get times 
    MPI_Reduce(&bcast_opt, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank==0)	
        printf("%lf\n", ftime);

    MPI_Reduce(&bcast_std, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank==0)
        printf("%lf\n", ftime);

    MPI_Reduce(&reduce_opt, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)	
		printf("%lf\n", ftime);

	MPI_Reduce(&reduce_std, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		printf("%lf\n", ftime);

	MPI_Reduce(&gather_opt, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)	
		printf("%lf\n", ftime);

	MPI_Reduce(&gather_std, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		printf("%lf\n", ftime);

    MPI_Reduce(&alltoall_opt, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)	
		printf("%lf\n", ftime);

	MPI_Reduce(&alltoall_std, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		printf("%lf\n", ftime);

	MPI_Finalize();
	return 0;
}
