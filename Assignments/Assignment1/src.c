#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

int n, num_time_steps;
MPI_Datatype rowtype, coltype;

void avg_compute(int s, int prow, int pcol, double** arr, double** recv_buff){
    double** buff = (double**)malloc(sizeof(double*)*n);
    for (int i=0; i<n; i++)
        buff[i] = (double*)malloc(sizeof(double)*n);
    for (int i=1; i<n-1; i++){
        for (int j=1; j<n-1; j++){
            buff[i][j] += arr[i-1][j]/4;
            buff[i][j] += arr[i+1][j]/4;
            buff[i][j] += arr[i][j+1]/4;
            buff[i][j] += arr[i][j-1]/4;
        }
    }

    for(int i=1; i<n-1; i++){
        if(prow == 0){
            buff[0][i] = arr[0][i-1]/3 + arr[0][i+1]/3 + arr[1][i]/3;
        }
        else{
            buff[0][i] = recv_buff[2][i]/4 + arr[0][i-1]/4 + arr[0][i+1]/4 + arr[1][i]/4;
        }

        if(prow == s-1){
            buff[n-1][i] = arr[n-1][i-1]/3 + arr[n-1][i+1]/3 + arr[n-2][i]/3;
        }
        else{
            buff[n-1][i] = recv_buff[3][i]/4 + arr[n-1][i-1]/4 + arr[n-1][i+1]/4 + arr[n-2][i]/4;
        }
    }
    
    for(int i=1; i<n-1; i++){
        if(pcol == 0){
            buff[i][0] = arr[i-1][0]/3 + arr[i+1][0]/3 + arr[i][1]/3;
        }
        else{
            buff[i][0] = recv_buff[0][i]/4 + arr[i-1][0]/4 + arr[i+1][0]/4 + arr[i][1]/4;
        }

        if(pcol == s-1){
            buff[i][n-1] = arr[i-1][n-1]/3 + arr[i+1][n-1]/3 + arr[i][n-2]/3;
        }
        else{
            buff[i][n-1] = recv_buff[2][i]/4 + arr[i-1][n-1]/4 + arr[i+1][n-1]/4 + arr[i][n-2]/4;
        }    
    }

    //corners
    if(prow == 0){
      if(pcol == 0){
        buff[0][0] = arr[0][1]/2 + arr[1][0]/2;
        buff[n-1][0] = recv_buff[3][0]/3 + arr[n-2][0]/3 + arr[n-1][1]/3;
        buff[0][n-1] = recv_buff[2][0]/3 + arr[0][n-2]/3 + arr[1][n-1]/3;
        buff[n-1][n-1] = recv_buff[2][n-1]/4 + recv_buff[3][n-1]/4 + arr[n-1][n-2]/4 + arr[n-2][n-1]/4;
      }
      else if(pcol == s-1){
        buff[0][0] = recv_buff[0][0]/3 + arr[0][1]/2 + arr[1][0]/2;
        buff[n-1][0] = recv_buff[3][0]/4 + recv_buff[0][n-1]/4 + arr[n-2][0]/4 + arr[n-1][1]/4;
        buff[0][n-1] = arr[0][n-2]/2 + arr[1][n-1]/2;
        buff[n-1][n-1] = recv_buff[3][n-1]/3 + arr[n-1][n-2]/3 + arr[n-2][n-1]/3;
      }
      else{
        buff[0][0] = recv_buff[0][0]/3 + arr[0][1]/3 + arr[1][0]/3;
        buff[n-1][0] = recv_buff[3][0]/4 + recv_buff[0][n-1]/4 + arr[n-2][0]/4 + arr[n-1][1]/4;
        buff[0][n-1] = recv_buff[2][0]/3 + arr[0][n-2]/3 + arr[1][n-1]/3;
        buff[n-1][n-1] = recv_buff[2][n-1]/4 + recv_buff[3][n-1]/4 + arr[n-1][n-2]/4 + arr[n-2][n-1]/4;
      }
    }
    else if(prow == s-1){
      if(pcol == 0){
        buff[0][0] = recv_buff[1][0]/3 + arr[0][1]/3 + arr[1][0]/3;
        buff[n-1][0] = arr[n-2][0]/2 + arr[n-1][1]/2;
        buff[0][n-1] = recv_buff[2][0]/4 + recv_buff[1][n-1]/4 + arr[0][n-2]/4 + arr[1][n-1]/4;
        buff[n-1][n-1] = recv_buff[2][n-1]/3 + arr[n-1][n-2]/3 + arr[n-2][n-1]/3;
      }
      else if(pcol == s-1){
        buff[0][0] = recv_buff[0][0]/4 + recv_buff[1][0]/4 + arr[0][1]/4 + arr[1][0]/4;
        buff[n-1][0] = recv_buff[0][n-1]/3 + arr[n-2][0]/3 + arr[n-1][1]/3;
        buff[0][n-1] = recv_buff[1][n-1]/3 + arr[0][n-2]/3 + arr[1][n-1]/3;
        buff[n-1][n-1] = arr[n-1][n-2]/2 + arr[n-2][n-1]/2;
      }
      else{
        buff[0][0] = recv_buff[0][0]/4 + recv_buff[1][0]/4 + arr[0][1]/4 + arr[1][0]/4;
        buff[n-1][0] = recv_buff[0][n-1]/3 + arr[n-2][0]/3 + arr[n-1][1]/3;
        buff[0][n-1] = recv_buff[2][0]/4 + recv_buff[1][n-1]/4 + arr[0][n-2]/4 + arr[1][n-1]/4;
        buff[n-1][n-1] = recv_buff[2][n-1]/3 + arr[n-1][n-2]/3 + arr[n-2][n-1]/3;
      }
    }
    else {
      if(pcol == 0){
        buff[0][0] = recv_buff[1][0]/3 + arr[0][1]/3 + arr[1][0]/3;
        buff[n-1][0] = recv_buff[3][0]/3 + arr[n-2][0]/3 + arr[n-1][1]/3;
        buff[0][n-1] = recv_buff[2][0]/4 + recv_buff[1][n-1]/4 + arr[0][n-2]/4 + arr[1][n-1]/4;
        buff[n-1][n-1] = recv_buff[2][n-1]/4 + recv_buff[3][n-1]/4 + arr[n-1][n-2]/4 + arr[n-2][n-1]/4;
      }
      else if(pcol == s-1){
        buff[0][0] = recv_buff[0][0]/4 + recv_buff[1][0]/4 + arr[0][1]/4 + arr[1][0]/4;
        buff[n-1][0] = recv_buff[3][0]/4 + recv_buff[0][n-1]/4 + arr[n-2][0]/4 + arr[n-1][1]/4;
        buff[0][n-1] = recv_buff[1][n-1]/3 + arr[0][n-2]/3 + arr[1][n-1]/3;
        buff[n-1][n-1] = recv_buff[3][n-1]/3 + arr[n-1][n-2]/3 + arr[n-2][n-1]/3;
      }
      else{
        buff[0][0] = recv_buff[0][0]/4 + recv_buff[1][0]/4 + arr[0][1]/4 + arr[1][0]/4;
        buff[n-1][0] = recv_buff[3][0]/4 + recv_buff[0][n-1]/4 + arr[n-2][0]/4 + arr[n-1][1]/4;
        buff[0][n-1] = recv_buff[2][0]/4 + recv_buff[1][n-1]/4 + arr[0][n-2]/4 + arr[1][n-1]/4;
        buff[n-1][n-1] = recv_buff[2][n-1]/4 + recv_buff[3][n-1]/4 + arr[n-1][n-2]/4 + arr[n-2][n-1]/4;
      }
    }

    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            arr[i][j] = buff[i][j];
    //free buffer
    for (int i=0; i<n; i++)
        free(buff[i]);
    free(buff);
}


//halo exchange using only send/recieve
double halo_send(double** arr, int rank, int size){
    //Initialise buffer(s)
    double** recv_buff = (double**)malloc(sizeof(double*)*4);
    for (int i=0; i<4; i++){
        recv_buff[i] = (double*)malloc(sizeof(double)*n);
    }
    int s = sqrt(size);
    int prow = rank/s;
    int pcol = rank%s;
    MPI_Status status;
    double etime = 0;
    //  Add code
    for (int t=0; t<num_time_steps; t++){
        double stime = MPI_Wtime();
        //  left exchange
        if (pcol > 0){
            for (int i=0; i<n; i++)
                MPI_Send(&arr[i][0], 1, MPI_DOUBLE, prow*s + pcol - 1, rank, MPI_COMM_WORLD);
            for (int i=0; i<n; i++)
                MPI_Recv(&recv_buff[0][i], 1, MPI_DOUBLE, prow*s + pcol - 1, prow*s + pcol - 1, MPI_COMM_WORLD, &status);
        }
        // right exchange
        if (pcol < s-1){
            for (int i=0; i<n; i++)
                MPI_Recv(&recv_buff[2][i], 1, MPI_DOUBLE, prow*s + pcol + 1, prow*s + pcol + 1, MPI_COMM_WORLD, &status);
            for (int i=0; i<n; i++)
                MPI_Send(&arr[i][n-1], 1, MPI_DOUBLE, prow*s + pcol + 1, rank, MPI_COMM_WORLD);
        }
        // top exhange
        if (prow > 0){
            for (int i=0; i<n; i++)
                MPI_Send(&arr[0][i], 1, MPI_DOUBLE, (prow-1)*s + pcol, rank, MPI_COMM_WORLD);
            for (int i=0; i<n; i++)
                MPI_Recv(&recv_buff[1][i], 1, MPI_DOUBLE, (prow-1)*s + pcol, (prow-1)*s + pcol, MPI_COMM_WORLD, &status);
        }
        // bottom exchange
        if (prow < s-1){
            for (int i=0; i<n; i++)
                MPI_Recv(&recv_buff[3][i], 1, MPI_DOUBLE, (prow+1)*s + pcol, (prow+1)*s + pcol, MPI_COMM_WORLD, &status);
            for (int i=0; i<n; i++)
                MPI_Send(&arr[n-1][i], 1, MPI_DOUBLE, (prow+1)*s + pcol, rank, MPI_COMM_WORLD);
        }
        etime+=MPI_Wtime() - stime;
        //compute avg
        avg_compute(s, prow, pcol, arr, recv_buff);
    }
    for(int i=0; i<4; i++){
        free(recv_buff[i]);
    }
    free(recv_buff);
    return etime;
}

double halo_packed(double** arr, int rank, int size){
    //halo exchange using MPI_Packed
    //Initialise buffer(s)
    double** recv_buff = (double**)malloc(sizeof(double*)*4);
    double** recvp_buff = (double**)malloc(sizeof(double*)*4);
    double** send_buff = (double**)malloc(sizeof(double*)*4);
    for (int i=0; i<4; i++){
        recv_buff[i] = (double*)malloc(sizeof(double)*n);
        recvp_buff[i] = (double*)malloc(sizeof(double)*n);
        send_buff[i] = (double*)malloc(sizeof(double)*n);
    }

    int s = sqrt(size);
    int prow = rank/s;
    int pcol = rank%s;
    MPI_Status status;
    int position = 0;
    double etime = 0;
    //Add code
    for (int t=0; t<num_time_steps; t++){
        //  left exchange
        double stime = MPI_Wtime();
        if (pcol > 0){
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Pack(&arr[i][0], 1, MPI_DOUBLE, send_buff[0], (n)*8, &position, MPI_COMM_WORLD);
            MPI_Send(send_buff[0], position, MPI_PACKED, prow*s + pcol - 1, rank, MPI_COMM_WORLD);
            MPI_Recv(recvp_buff[0], (n)*8, MPI_PACKED, prow*s + pcol - 1, prow*s + pcol - 1, MPI_COMM_WORLD, &status);
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Unpack (recvp_buff[0], (n)*8, &position, &recv_buff[0][i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        // right exchange
        if (pcol < s-1){
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Pack(&arr[i][n-1], 1, MPI_DOUBLE, send_buff[2], (n)*8, &position, MPI_COMM_WORLD);
            MPI_Recv(recvp_buff[2], (n)*8, MPI_PACKED, prow*s + pcol + 1, prow*s + pcol + 1, MPI_COMM_WORLD, &status);
            MPI_Send(send_buff[2], position, MPI_PACKED, prow*s + pcol + 1, rank, MPI_COMM_WORLD);
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Unpack (recvp_buff[2], (n)*8, &position, &recv_buff[2][i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        // top exhange
        if (prow > 0){
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Pack(&arr[0][i], 1, MPI_DOUBLE, &send_buff[1][0], (n)*8, &position, MPI_COMM_WORLD);
            MPI_Send(send_buff[1], position, MPI_PACKED, (prow-1)*s + pcol, rank, MPI_COMM_WORLD);
            MPI_Recv(recvp_buff[1], (n)*8, MPI_PACKED, (prow-1)*s + pcol, (prow-1)*s + pcol, MPI_COMM_WORLD, &status);
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Unpack (recvp_buff[1], (n)*8, &position, &recv_buff[1][i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        // bottom exchange
        if (prow < s-1){
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Pack(&arr[n-1][i], 1, MPI_DOUBLE, &send_buff[3][0], (n)*8, &position, MPI_COMM_WORLD);
            MPI_Recv(recvp_buff[3], (n)*8, MPI_PACKED, (prow+1)*s + pcol, (prow+1)*s + pcol, MPI_COMM_WORLD, &status);
            MPI_Send(send_buff[3], position, MPI_PACKED, (prow+1)*s + pcol, rank, MPI_COMM_WORLD);
            position = 0;
            for (int i=0; i<n; i++)
              MPI_Unpack (recvp_buff[3], (n)*8, &position, &recv_buff[3][i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }

        //compute avg
        etime += MPI_Wtime() - stime;
        avg_compute(s, prow, pcol, arr, recv_buff);
    }
    for(int i=0; i<4; i++){
        free(recv_buff[i]);
        free(send_buff[i]);
        free(recvp_buff[i]);
    }
    free(recv_buff);
    free(send_buff);
    free(recvp_buff);

    return etime;
}

double halo_type(double** arr, int rank, int size){
    //halo exchange using derived datatypes
    //Initialise buffer(s)
    double** recv_buff = (double**)malloc(sizeof(double*)*4);
    for (int i=0; i<4; i++){
        recv_buff[i] = (double*)malloc(sizeof(double)*n);
    }
    int s = sqrt(size);
    int prow = rank/s;
    int pcol = rank%s;
    MPI_Status status;
    double etime = 0;
    //  Add code
    for (int t=0; t<num_time_steps; t++){
        //  left exchange
        double stime = MPI_Wtime();
        if (pcol > 0){
            MPI_Send(&arr[0][0], 1, coltype, prow*s + pcol - 1, rank, MPI_COMM_WORLD);
            MPI_Recv(recv_buff[0], n, MPI_DOUBLE, prow*s + pcol - 1, prow*s + pcol - 1, MPI_COMM_WORLD, &status);
        }
        // right exchange
        if (pcol < s-1){
            MPI_Recv(recv_buff[2], n, MPI_DOUBLE, prow*s + pcol + 1, prow*s + pcol + 1, MPI_COMM_WORLD, &status);
            MPI_Send(&arr[0][n-1], 1, coltype, prow*s + pcol + 1, rank, MPI_COMM_WORLD);
        }
        // top exhange
        if (prow > 0){
            MPI_Send(arr[0], 1, rowtype, (prow-1)*s + pcol, rank, MPI_COMM_WORLD);
            MPI_Recv(recv_buff[1], n, MPI_DOUBLE, (prow-1)*s + pcol, (prow-1)*s + pcol, MPI_COMM_WORLD, &status);
        }
        // bottom exchange
        if (prow < s-1){
            MPI_Recv(recv_buff[3], n, MPI_DOUBLE, (prow+1)*s + pcol, (prow+1)*s + pcol, MPI_COMM_WORLD, &status);
            MPI_Send(arr[n-1], 1, rowtype, (prow+1)*s + pcol, rank, MPI_COMM_WORLD);
        }
        etime += MPI_Wtime() - stime;
        //compute avg
        avg_compute(s, prow, pcol, arr, recv_buff);
    }
    for(int i=0; i<4; i++){
        free(recv_buff[i]);
    }
    free(recv_buff);
    return etime;
}

int main(int argc, char** argv){
    // TODO: check if number of arguments is correct

    int N = atoi(argv[1]);
    num_time_steps = atoi(argv[2]);
    double exectime, maxtime;
    n = sqrt(N);
    double** data = (double**) malloc(sizeof(double*)*n);
    for(int i=0; i<n; i++){
        data[i] = (double*)malloc(sizeof(double)*n);
    }
    int rank, size;
    MPI_Status status; 

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Type_vector(n, 1, n, MPI_DOUBLE, &coltype);
    MPI_Type_contiguous(n, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&coltype); 
    MPI_Type_commit(&rowtype);

    // TODO: add logic for basic operations
    // Initialisation
    srand(time(0));
    for (int i=0; i<n; i++)
	    for(int j=0; j<n; j++)
	      data[i][j] = rand() % size;

    exectime = halo_send(data, rank, size);
    MPI_Reduce (&exectime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
    if (rank == size-1) 
	    printf ("%lf\n", maxtime);

    exectime = halo_packed(data, rank, size);
    MPI_Reduce (&exectime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
    if (rank == size-1) 
        printf ("%lf\n", maxtime);

    exectime = halo_type(data, rank, size);
    MPI_Reduce (&exectime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
    if (rank == size-1)
        printf ("%lf", maxtime);

    for(int i=0; i<n; i++){
        free(data[i]);
    }
    free(data);

    MPI_Type_free(&coltype);
    MPI_Type_free(&rowtype);
    MPI_Finalize();
    return 0;
}


