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
    double * mat;  // local matrix use for reading
    
    double * matA; // local scalapack matrix
    int NA, MA;
    
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
    
    // one assume full block i.e. N%NB = 0
    // other wise one has to consider the case where the last row and/or last column may have a different dimension than NB x MB
    
    NG = sqrt(p);
    MG = NG; // square processor grid
   
    // one assume that there are more blocks than process i.e (N/NB) > NG
    i_grid = get_i_grid_ofProcess(rank,NG,MG);
    j_grid = get_j_grid_ofProcess(rank,NG,MG);
    int remain = (N/NB)%NG;
    if (i_grid < remain) {
      NA = ((N/NB)/NG + 1) * NB;
    } else {
      NA = (N/NB)/NG * NB;
    }
    printf("%d/%d: size local matrix %d \n",rank,p,NA);
    //matA = malloc(NA*NA * sizeof(double)); // with square grid and square blocks matA should be square
    
    ///
    int MAXT = n_blocks/p+1;
    int MAXB = (n_blocks+1)/2+1;
    int send[p][MAXT][MAXB];
    int rang;
    for(rang=0;rang<p;rang++){ // loop over all process
      for(t=0;t<MAXT;t++){ // loop A over all process tasks
	i_block = (rang+t*p)%p;
	for(t_block=0;t_block<(MAXB);t_block++){ // loop B over all block tasks
	   j_block = mpi_get_diag_block_id(i_block, t_block, n_blocks);
	   send[rang][t][t_block] = get_rank_block(i_block,j_block,NG,MG);
	   if(rank==0) printf("%d/%d: %d should send block (%d,%d) to %d (%d,%d,%d) \n",rank,p,rang, i_block,j_block, send[rang][t][t_block],rang,t,t_block);
	}
      }
    }
   
    // WRONG
    // one should invert send to get recv 
    int recv[p][MAXT][MAXB];
     for(rang=0;rang<p;rang++){ // loop over all process
      for(t=0;t<MAXT;t++){ // loop A over all process tasks
	for(t_block=0;t_block<(MAXB);t_block++){ // loop B over all block tasks
	   recv[send[rang][t][t_block]][t][t_block]=rang;
	}
      }
    }
     
     for(rang=0;rang<p;rang++){ // loop over all process
      for(t=0;t<(MAXT);t++){ // loop A over all process tasks
	for(t_block=0;t_block<(MAXB);t_block++){ // loop B over all block tasks
	     if(rank==0) printf("%d/%d: %d should receive data from %d (%d,%d,%d) \n",rank,p,rang, recv[rang][t][t_block],rang,t,t_block);
	}
      }
    }
    
    
    for(t=0;t<n_pTasks;t++){ // loop A over all process tasks
      i_block = (rank+t*p); // the diagonal block index depends on the rank index and the task index
      readStore(&store,i_block,"./data/ReducedNormalsTest");
      t_block = 0;
      while(readNextBlockDimension(dim,store)!=-1) { // loop B over all block tasks
	
	
	printf("%d/%d: read %d,%d,%d,%d \n",rank,p,dim[0],dim[1],dim[2],dim[3]);
	mat = malloc((dim[1]-dim[0])*(dim[3]-dim[2]) * sizeof(double));         
	readNextBlock(dim[0],dim[1],dim[2],dim[3],mat,store);
	
	j_block = mpi_get_diag_block_id(i_block, t_block, n_blocks);
	
	
	
	// the block i_block, i_block should be sent to 1 scalapack grid process
	// the block i_block, j_block should be sent to 2 scalapack grid process (symmetric matrix)
	rankL = get_rank_block(i_block, j_block, NG, MG);
	rankU = get_rank_block(j_block, i_block, NG, MG);
	printf("%d/%d: should send block (%d,%d) to %d and %d \n",rank,p,i_block,j_block,rankL, rankU);
	
	// TODO there is a bug
	//printf("%d/%d: block %d = %d \n",rank,p,i_block*n_blocks+j_block,send[rank][t][t_block]); // 
	//printf("%d/%d: should receive data from %d (%d,%d,%d) \n",rank,p,recv[rank][t][t_block],rank,t,t_block);
	
	free(mat);
	t_block++;
      }
      closeStore(store);
    }
    
    //free(matA);
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

