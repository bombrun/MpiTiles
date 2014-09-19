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