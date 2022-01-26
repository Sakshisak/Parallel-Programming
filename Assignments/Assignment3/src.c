#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

int rank, size;
int r;      // ratio for division in 2 nodes

typedef struct TempList{
  char* line;
  int size;
  struct TempList* next;
}TempList;

// helper function to get name of the node.
int get_my_node(){
	char host_name[10];
	int name_len, node = 0;
	MPI_Get_processor_name(host_name, &name_len);
	for (int i=5; i<name_len; i++) 
		node = node*10 + (int)(host_name[i] - '0');
	return node;
}

float Scatter_binomial_one_node(float* data, int M, int N, float* minima){
    int tsize = size;
    int trank = rank;

    int num = tsize + log2(tsize)*tsize/2;
    int den = size + log2(size)*size/2;
    // the remaining values will be processes at root. 
    // [ [####] [#] [##] [#] | (pad) ]
    //     0     1   2    3      0
    int pad = N%den;
    N-=pad;
    if (rank != 0){
        //find size of subtree rooted at this process. 
        while(trank > 0){
            tsize = tsize/2;
            if (trank >= tsize){
                trank = trank - tsize;
            }
        }

        //each process recieves atmost once. Amount of data recieved is M x N*ratio
        MPI_Status stat;
        MPI_Recv(data, M*num*(N/den), MPI_FLOAT, rank-tsize, rank, MPI_COMM_WORLD, &stat);
    }
    
    int num_it = (int)log2(tsize);
    int temp = tsize/2;
    int it = 0;

    MPI_Status stat[num_it];
    MPI_Request req[num_it];
    int offset = num;
    while(temp > 0){
        //the offset from which data is to be sent. [----[###]|     ]       
        int num_t = temp + log2(temp)*temp/2;
        offset = offset - num_t;
        MPI_Isend(&data[(offset*N/den)*M], M*num_t*(N/den), MPI_FLOAT, rank+temp, rank+temp, MPI_COMM_WORLD, &req[it]);
        temp /= 2;
        it++;
    }

    //compute logic
    for (int i=0; i<M; i++){
        float y_min = FLT_MAX;
        for (int j=0; j<tsize*(N/den); j++){
            if (data[j*M + i] < y_min)
                y_min = data[j*M + i];
        }
        minima[i] = y_min;
    }
    //computation for remaining stations(rows) at root
    if (rank == 0){
        for (int i=0; i<M; i++){
           float y_min = minima[i];
            for (int j=N; j<N+pad; j++){
                if (data[j*M + i] < y_min)
                    y_min = data[j*M + i];
            }
            minima[i] = y_min; 
        }
    }

    //compute local minima across all years 
    float g_min = FLT_MAX;
    for (int i=0; i<M; i++){
        g_min = minima[i] < g_min ? minima[i] : g_min;
    }
    MPI_Waitall(num_it, req, stat);
    return g_min;
}

float Scatter_binomial_two_node(float* data, long M, long N, float* minima){
    long tsize = size;
    int trank = rank;
    // calculate weight factors using ratio 'r'
    long num = rank >= size/2 ? 1 : r;
    long den = size > 1 ? (size/2)*(r+1) : 1;
    long pad = N%den;
    N -= pad;

    if (rank != 0){
        while(trank > 0){
            tsize = tsize/2;
            if (trank >= tsize){
                trank = trank - tsize;
            }
        }

        // every process recieves data once. 
        // amount of data recieved in ratio*N/(size/2)
        MPI_Status stat;
        MPI_Recv(data, M*num*tsize*(N/den), MPI_FLOAT, rank-tsize, rank, MPI_COMM_WORLD, &stat);
    }
    
    int num_it = (int)log2(tsize);
    long temp = tsize/2;
    int it = 0;

    MPI_Status stat[num_it];
    MPI_Request req[num_it];
    long offset = num*tsize/2;

    // issue requests to send data further. 
    while(temp > 0){
        long num_t = rank + temp >= size/2 ? temp : temp*r;
        int flag;
        if (it!=0)
            MPI_Test(&req[it-1], &flag, MPI_STATUS_IGNORE);
        MPI_Isend(&data[(offset*N/den)*M], M*num_t*(N/den), MPI_FLOAT, rank+temp, rank+temp, MPI_COMM_WORLD, &req[it]);
        temp /= 2;
        offset /= 2;
        it++;
    }

    // compute logic
    for (int i=0; i<M; i++){
        float y_min = FLT_MAX;
        for (int j=0; j<tsize*(N/den); j++){
            if (data[j*M + i] < y_min)
                y_min = data[j*M + i];
        }
        minima[i] = y_min;
    }

    //computation for remaining stations(rows) at root
    if (rank == 0){
        for (int i=0; i<M; i++){
           float y_min = minima[i];
            for (int j=N; j<N+pad; j++){
                if (data[j*M + i] < y_min)
                    y_min = data[j*M + i];
            }
            minima[i] = y_min; 
        }
    }

    float g_min = FLT_MAX;

    for (int i=0; i<M; i++){
        g_min = minima[i] < g_min ? minima[i] : g_min;
    }

    // fail for all sends to complete
    MPI_Waitall(num_it, req, stat);
    return g_min;
}

int main(int argc, char** argv){
    double stime, etime, ftime = 0;
    float* data;
    long M, N, pad;
    int node1=1, node2=1;
		
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // set value of division ratio.
    // -1 so that denominator is power of 2.
    r = size^2 - 1;

    //file input and processing
    if (rank == 0){
        int bytes = 0, line_ctr = 0;
        size_t buffsize = 4096;
        char* line = malloc(sizeof(char)*buffsize);
        TempList* Head = (TempList*)malloc(sizeof(TempList));
        TempList* ptr = Head;

        FILE* infile = fopen(argv[1], "r");
        bytes = getline(&line, &buffsize, infile);
        // read data into linked list
        while(bytes > 0){
            line_ctr++;	
            TempList* node = (TempList*)malloc(sizeof(TempList));
            node->next = NULL;
            node->line = line;
            node->size = bytes;
            ptr->next = node;
            ptr = node;
            line = malloc(sizeof(char)*buffsize);
            bytes = getline(&line, &buffsize, infile);
        }
        free(line);

        ptr = Head->next;
        const char s[2] = ",";
        char *token;
        token = strtok(ptr->line, s);
        token = strtok(NULL, s);
        token = strtok(NULL, s);
        
        int col = 0;
        while( token != NULL ) {
            col++;
            token = strtok(NULL, s); 
        }
        
        N = line_ctr-1; M = col;

        data = (float*)malloc(sizeof(float)*M*(N));
        
        // populate data array
        int i = 0;
        while(ptr->next){
            token = strtok(ptr->next->line, s);
            token = strtok(NULL, s);
            token = strtok(NULL, s);
            int j = 0;
            while( token != NULL ) {
                data[i*M + j] = atof(token);
                token = strtok(NULL, s); 
                j++;
            }

            TempList* temp = ptr->next;   
            free(ptr->line);
            free(ptr);
            ptr = temp;
            i++;
        }

        free(ptr->line);
        free(ptr);
        free(Head);
        
        fclose(infile);
	}

    // Halt while root process reads data
    MPI_Barrier(MPI_COMM_WORLD);

    // input done, start timing the code.
    stime = MPI_Wtime();

    MPI_Request req;
    MPI_Status stat;

    if(size>1){
      if(rank == size/2){
        node1 = get_my_node();
        MPI_Send(&node1, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
      }
      else if(rank == 0){
        node2 = get_my_node();
        MPI_Recv(&node1, 1, MPI_INT, size/2, size/2, MPI_COMM_WORLD, &stat);
      }
    }

    
    int var[3];
    if (rank == 0){
        var[0] = M; 
        var[1] = N;
        var[2] = (node1 == node2) ?  1 : 2;
    }

    // get values of M and N
    MPI_Bcast(var, 3, MPI_INT, 0, MPI_COMM_WORLD);
    M = var[0];
    N = var[1];


    float minima[M];
    float g_min;

    if (rank!=0){
        // array for root was already allocated.
        data = (float*)malloc(sizeof(float)*M*N);
    }

    if(var[3] == 1)
        g_min = Scatter_binomial_one_node(data, M, N,minima);
    else
        g_min = Scatter_binomial_two_node(data, M, N, minima);

    // Get global minima of all processes. 
    float yearly_min[M];
    float overall_min = FLT_MAX;

    MPI_Ireduce(minima, yearly_min, M, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD, &req);
    
    // get global minimum of temperature across all years
    MPI_Reduce(&g_min, &overall_min, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    free(data);

    MPI_Wait(&req, &stat);
    
    etime = MPI_Wtime() - stime;
    MPI_Reduce(&etime, &ftime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0){
        for (int i=0; i<M; i++){
        printf("%0.2f", yearly_min[i]);
        if (i!=M-1)
            printf(",");
        }
        printf("\n%0.2f\n%lf\n", overall_min, ftime);
    }

    MPI_Finalize();
    return 0;
}





