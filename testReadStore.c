#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "normals.h"
#include "mpiutil.h"
#include "matrixBlockStore.h"

/* Test program 
 * created 23/09/2014
 * author Alex Bombrun
 * 
 * mpicc -Wall -g -o read.mpi testReadStore.c mpiutil.c normals.c matrixBlockStore.c
 * mpirun -n 3 read.mpi
 */
int main(int argc, char **argv) {
  
    FILE* store;
    int p, rank;
    
    int dim[4];
    double * mat;
    
    int t;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("%d/%d: test read store\n",rank,p);
    readStore(&store,rank,"./data/ReducedNormalsTest");
    
    while(readNextBlockDimension(dim,store)!=-1) {
      printf("%d/%d: read %d,%d,%d,%d \n",rank,p,dim[0],dim[1],dim[2],dim[3]);
      mat = malloc((dim[1]-dim[0])*(dim[3]-dim[2]) * sizeof(double));
      readNextBlock(dim[0],dim[1],dim[2],dim[3],mat,store);
      free(mat);
    }
    closeStore(store);
    MPI_Finalize();
    return 0;
}


