#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "mpiutil.h"

/**
 * split ng in p segments according to rang,
 * assume ng>p>rang
 * return the left index of the segment
 */
int mpi_get_i0(int ng, int rang, int p){
    int step = ng/p; // size of the segment
    return rang*step;
}

/**
 * 
 */
int mpi_get_i1(int ng, int rang, int p){
    int step = ng/p;
    if (rang+1==p) {
      return ng;
    } else {
      return (rang+1)*step;
    }
}

int mpi_get_total_blocks(int n){
    return (n*(n+1))/2;
}

int mpi_get_diag_block_id(int i_block,int block_task_id, int n_diag_blocks){
    return (i_block+block_task_id)%n_diag_blocks;
}


int get_n_blocks(int argc, char **argv, int p){
  int n_blocks;
  if (argc == 2) {
      //const char* prog_name = *argv++;
      //const char* nblocks_char = *argv;
      n_blocks= (int) strtol(argv[1], NULL, 10);
       if(n_blocks<p) {
	printf("Configuration error: n_blocks %i is smaller than the number of process %d, set to default (the number of process)\n",n_blocks,p);
	n_blocks = p; 
      }
    }
    else {
      n_blocks = p; // one diagonal block per process
    }
    return n_blocks;
}

int get_n_pTasks(int p, int rank, int n_blocks){
    // diagonal block distribution over all the process
    int remainBlocks = n_blocks%p;
    int n_pTasks = n_blocks/p;
    if (rank<remainBlocks) n_pTasks+=1;
    return n_pTasks;
}

int get_n_blockTasks(int r, int n_blocks){
   // tasks distribution over all the diagonal blocks
    int remainTasks=mpi_get_total_blocks(n_blocks)%n_blocks;
    int nTasks = mpi_get_total_blocks(n_blocks)/n_blocks;
    if (r<remainTasks) nTasks+=1;
    return nTasks;
}