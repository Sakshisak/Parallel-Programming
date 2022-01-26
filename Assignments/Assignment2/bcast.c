#include <collective.h>

int Bcast_Binomial(void* buffer, int count, MPI_Datatype datatype, int* ranks, int size, MPI_Comm comm){
    int rank, tsize=size, start=0;
    MPI_Comm_rank(comm, &rank);

    if (rank != ranks[0]){
        MPI_Request request[2];
        MPI_Status status[2];
        MPI_Irecv(buffer, count, datatype, MPI_ANY_SOURCE, rank, comm, &request[0]);
        MPI_Irecv(&tsize, 1, MPI_INT, MPI_ANY_SOURCE, rank, comm, &request[1]);
        for (int i=0; i<size; i++)
            if(ranks[i] == rank) start = i;
        MPI_Waitall(2, request, status);
    }
    tsize = tsize/2;
    MPI_Request req[size];
    MPI_Status stat[size];
    int it=0;
    while(tsize > 0){
        MPI_Isend(buffer, count, datatype, ranks[start + tsize], ranks[start + tsize], comm, &req[it++]);
        MPI_Send(&tsize, 1, MPI_INT, ranks[start + tsize], ranks[start + tsize], comm);
        tsize = tsize/2;
    }
    MPI_Waitall(it, req, stat);
    return 0;
}

int Scatter_Binomial(void* sbuff, int scount, MPI_Datatype stype, void* rbuff, int rcount, int* ranks, int size, MPI_Comm comm){
    int rank, tsize=size, start=0;
    MPI_Comm_rank(comm, &rank);
    double* buffer;
    if (rank != ranks[0]){
        for (int i=0; i<size; i++)
            if(ranks[i] == rank) start = i;
        int id = start;
        while (id > 0){
            tsize = tsize/2;
            if(id >= tsize){
                id = id - tsize;
            }
        }
        buffer = (double*)malloc(sizeof(double)*scount*tsize);
        MPI_Status status;
        MPI_Recv(buffer, scount*tsize, stype, MPI_ANY_SOURCE, rank, comm, &status);
    }
    else
        buffer = sbuff;
    
    for(int i=0; i<rcount; i++){
        ((double*)rbuff)[i] = ((double*)buffer)[i];
    }
    tsize = tsize/2;
    MPI_Request req[size];
    MPI_Status stat[size];
    int it=0;
    while(tsize > 0){
        MPI_Isend(&buffer[tsize*scount], scount*tsize, stype, ranks[start + tsize], ranks[start + tsize], comm, &req[it++]);
        tsize = tsize/2;
    }
    MPI_Waitall(it, req, stat);
    if (rank != ranks[0])
        free(buffer);
    return 0;
}

int Allgather_Ring(void* sbuff, int scount, MPI_Datatype stype, void* rbuff, int rcount, int* ranks, int size, MPI_Comm comm){
    if (size == 1) return 0;
    int rank, start = 0;
    MPI_Comm_rank(comm, &rank);
    for(int i=0; i<size; i++)
        if (ranks[i] == rank) start = i;

    int offset = -1;
    MPI_Status status, stat[size];
    MPI_Request request[size];

    MPI_Isend(sbuff, scount, stype, ranks[(start+1)%size], 0, comm, &request[0]);
    for(int i=0; i<rcount; i++)
        ((double*)rbuff)[start*scount + i] = ((double*)sbuff)[i];

    for (int i=1; i<size-1; i++){
        MPI_Recv(&((double*)rbuff)[((start+offset+size)%size)*scount], scount, stype, ranks[(start-1+size)%size], i-1, comm, &status);
        MPI_Isend(&((double*)rbuff)[((start+offset+size)%size)*scount], scount, stype, ranks[(start+1)%size], i, comm, &request[i]);   
        offset--;
    }
    MPI_Recv(&((double*)rbuff)[((start+offset+size)%size)*scount], scount, stype, ranks[(start-1+size)%size], size-2, comm, &status);
    MPI_Waitall(size-1, request, stat);
    return 0;
}