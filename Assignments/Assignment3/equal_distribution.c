#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

int rank, size;

typedef struct TempList{
  char* line;
  int size;
  struct TempList* next;
}TempList;

float Scatter_binomial(float* data, int M, int N, float* minima){
    int tsize = size;
    int trank = rank;

    if (rank != 0){
        while(trank > 0){
            tsize = tsize/2;
            if (trank >= tsize){
                trank = trank - tsize;
            }
        }
        MPI_Status stat;
        MPI_Recv(data, M*tsize*N/size, MPI_FLOAT, rank-tsize, rank, MPI_COMM_WORLD, &stat);
    }
    
    int num_it = (int)log2(tsize);
    int temp = tsize/2;
    int it = 0;

    MPI_Status stat[num_it];
    MPI_Request req[num_it];
    while(temp > 0){
        MPI_Isend(&data[(temp*N/size)*M], M*temp*N/size, MPI_FLOAT, rank+temp, rank+temp, MPI_COMM_WORLD, &req[it]);
        temp /= 2;
        it++;
    }

    //compute logic
    for (int i=0; i<M; i++){
        float y_min = FLT_MAX;
        for (int j=0; j<N/size; j++){
            if (data[j*M + i] < y_min)
                y_min = data[j*M + i];
        }
        minima[i] = y_min;
    }

    float g_min = FLT_MAX;

    for (int i=0; i<M; i++){
        g_min = minima[i] < g_min ? minima[i] : g_min;
    }
    MPI_Waitall(num_it, req, stat);

    return g_min;
}

int main(int argc, char** argv){
	double stime, etime, ftime = 0;
    float* data;
    int M, N;
		
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    //file input and processing
	if (rank == 0){
        int bytes = 0, line_ctr = 0;
        size_t buffsize = 4096;
        char* line = malloc(sizeof(char)*buffsize);
        TempList* Head = (TempList*)malloc(sizeof(TempList));
        TempList* ptr = Head;

        FILE* infile = fopen(argv[1], "r");
        bytes = getline(&line, &buffsize, infile);
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

        // padding
        int pad = N%size ? size - N%size : 0;
        data = (float*)malloc(sizeof(float)*M*(N+pad));
  
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
        
        // more about padding ...
        for (i=0; i<M; i++){
            for (int j=0; j<pad; j++){
                data[(j+N)*M + i] = FLT_MAX;
            }
        }
        N += pad;
        fclose(infile);
	}

    //Halt while root process reads data
    MPI_Barrier(MPI_COMM_WORLD);

    //input done, start timing the code.
    stime = MPI_Wtime();

    int var[2];
    if (rank == 0){
        var[0] = M; 
        var[1] = N;
    }

    MPI_Bcast(var, 2, MPI_INT, 0, MPI_COMM_WORLD);
    M = var[0];
    N = var[1];

    MPI_Request req;
    MPI_Status stat;

    float minima[M];
    float g_min;

    if (rank!=0){
        data = (float*)malloc(sizeof(float)*M*N);
    }

    g_min = Scatter_binomial(data, M, N, minima);

    // Get global minima of all processes. 
    float yearly_min[M];
    float overall_min = FLT_MAX;

    MPI_Ireduce(minima, yearly_min, M, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD, &req);
    
    MPI_Reduce(&g_min, &overall_min, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

    free(data);

    MPI_Wait(&req, &stat);
    etime = MPI_Wtime() - stime;
    //MPI_Reduce on etime and store in ftime; 
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



