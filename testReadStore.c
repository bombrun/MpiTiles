#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "normals.h"
#include "mpiutil.h"
#include "matrixBlockStore.h"

int get_i_grid_ofProcess(int rank, int NG, int MG);
int get_j_grid_ofProcess(int rank, int NG, int MG);
int get_rank_grid(int i_grid, int j_grid, int NG, int MG);
int get_i_grid(int i_block, int j_block, int NG, int MG);
int get_j_grid(int i_block, int j_block, int NG, int MG);
int get_rank_block(int i_block, int j_block, int NG, int MG);


/**
    input c
    assume c = n*n;
    return n
    */
int sqrt(int c);

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
    
    int rankU, rankL;
    
    int NB , MB;
    int N , M;
    int NG, MG; // scalapack grid process   p = NG x MG!
    int n_blocks, n_pTasks;
    int i_block, j_block;
    int i_grid, j_grid;
    int dim[4];
    double * mat;
    
    int t, t_block;
    
    const char* profileG_file_name= "./data/NormalsG/profile.txt";
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("%d/%d: test read store\n",rank,p);
    n_blocks = get_n_blocks(argc,argv,p); // assume that one use the same input number as the one used to generate the matrix
    n_pTasks = get_n_pTasks(p,rank,n_blocks); // the number of tasks per process
  
    
    N = getNumberOfLine(profileG_file_name); // the dimension of the matrix;
    M = N; // square matrix
 
    NB = mpi_get_i1(N,0,n_blocks);
    MB = NB; // square blocks
    
    NG = sqrt(p);
    MG = NG; // square processor grid
   
    
    for(t=0;t<n_pTasks;t++){
      i_block = (rank+t*p); // the diagonal block index depends on the rank index and the task index
      readStore(&store,i_block,"./data/ReducedNormalsTest");
      t_block = 0;
      while(readNextBlockDimension(dim,store)!=-1) {
	printf("%d/%d: read %d,%d,%d,%d \n",rank,p,dim[0],dim[1],dim[2],dim[3]);
	mat = malloc((dim[1]-dim[0])*(dim[3]-dim[2]) * sizeof(double));         
	readNextBlock(dim[0],dim[1],dim[2],dim[3],mat,store);
	
	j_block = mpi_get_diag_block_id(i_block, t_block, n_blocks);
	t_block++;
	
	// the block i_block, i_block should be sent to 1 scalapack grid process
	// the block i_block, j_block should be sent to 2 scalapack grid process (symmetric matrix)
	// problem the last row and last column may have a different dimension than NB x MB
	rankL = get_rank_block(i_block, j_block, NG, MG);
	rankU = get_rank_block(j_block, i_block, NG, MG);
	printf("%d/%d: should send data to %d and %d \n",rank,p,rankL, rankU);
	free(mat);
      }
      closeStore(store);
    }
    MPI_Finalize();
    return 0;
}


int get_i_grid_ofProcess(int rank, int NG, int MG){
    return  rank/MG;  
}


int get_j_grid_ofProcess(int rank, int NG, int MG){
    int i_g= rank/MG;
    return rank-i_g*MG;
}

int get_rank_grid(int i_grid, int j_grid, int NG, int MG){
    return i_grid*MG+j_grid;
}

int get_i_grid(int i_block, int j_block, int NG, int MG){
    return i_block/NG;
}

int get_j_grid(int i_block, int j_block, int NG, int MG){
    return j_block/MG;
}

int get_rank_block(int i_block, int j_block, int NG, int MG){
    int ig = i_block/NG;
    int jg = j_block/MG;
    return get_rank_grid(ig, jg, NG, MG);
}



/**
    input c
    assume c = n*n;
    return n
    */
int sqrt(int c)
{
    int n = 0;
    while( n*n < c)
    {
        n++;
    }
    return n;

}

