#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "normals.h"
#include "mpiutil.h"
#include "matrixBlockStore.h"


#include "mkl_lapacke.h"
#include "mkl.h"

/**
    input c
    assume c = n*n;
    return n
    */
int sqrt(int c);
int saveMatrix(long long int dim, double * mat, const char* fileName);

/* Test program 
 * created 23/09/2014
 * author Alex Bombrun
 * 
 * icc -O1  -o eigen.exe lapackReadStore.c mpiutil.c normals.c matrixBlockStore.c -mkl
 * ./eigen.exe 4 4
 *
 */
int main(int argc, char **argv) {
  
    FILE* store;
    
    int N , M;
    int i, j;
    
    int n_blocks;
    int NB, MB;
    int i_block, j_block;
    int dim[4];
    double * mat;  // local matrix block use for reading
    
    
    long long int sizeA, id;
    double * matA; // lapack matrix
    
    int t, t_block;
    
    const char* profileG_file_name= "./data/NormalsG/profile.txt";
    const char* store_location = "./data/ReducedNormals";
    
     MKL_INT mkl_threads; // on Venus should be set to the number of processor used
     MKL_INT n;
     MKL_INT lda2;
     MKL_INT info;
    
    
     if (argc == 3) {
	//printf("%s %s %s\n", argv[0], argv[1], argv[2]);
	n_blocks= (int) strtol(argv[1], NULL, 10);
	mkl_threads = (int) strtol(argv[2], NULL, 10);
     } else {
	printf("Usage: expect 2 integers \n");
	printf(" first  : the number of diagonal blocks \n");
	printf(" second : the number of threads \n");
	exit( -1);
     }
    
  

    printf("read store\n");
   
    N = getNumberOfLine(profileG_file_name); // the dimension of the matrix;
    M = N; // square matrix
    
    sizeA = N;
    sizeA *= M;
    printf("matrix A : %d x %d \n",N,M);
    matA = malloc(sizeA * sizeof(double));
   
    
    for(i_block=0;i_block<n_blocks;i_block++){
      readStore(&store,i_block,store_location);
      t_block = 0;
      while(readNextBlockDimension(dim,store)!=-1) { // loop B over all block tasks
	j_block = mpi_get_diag_block_id(i_block, t_block, n_blocks);
	mat = malloc((dim[1]-dim[0])*(dim[3]-dim[2]) * sizeof(double));         
	readNextBlock(dim[0],dim[1],dim[2],dim[3],mat,store);
	printf("read block (%d,%d) with global indices (%d,%d,%d,%d) \n",i_block,j_block,dim[0],dim[1],dim[2],dim[3]);
	NB = dim[1]-dim[0];
	MB = dim[3]-dim[2];
	for(i = dim[0];i<dim[1];i++){
	  for(j = dim[2];j<dim[3];j++){
	      id = i;
	      id *= M;
	      id += j; // i*M+j
	      matA[id] = mat[(i-dim[0])*MB+(j-dim[2])];
	  }
	}
	// transpose
	if(j_block != i_block){
	  for(i = dim[0];i<dim[1];i++){
	    for(j = dim[2];j<dim[3];j++){
		id = j;
		id *= M;
		id += i;
		matA[id] = mat[(i-dim[0])*MB+(j-dim[2])];
	    }
	  } 
	}
	
	free(mat);
	t_block++;
      }
      closeStore(store);
    }
    
    
    printf("Start eigenvalues decomposition \n");
   
    MKL_Set_Num_Threads(mkl_threads);
    double w2[N];
    n = N;
    lda2 = N;
    
    info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, matA, lda2, w2 );

    if( info > 0 ) {
	printf( "MKL algorithms failed to compute eigenvalues.\n" );
	exit(-2);
    } else {
       printf("Finished eigenvalues decomposition \n");
      saveMatrix(N,w2,"eigenValues.txt");
      saveMatrix(sizeA,matA,"eigenVectors.txt");
    }
    
    free(matA);
    return 0;
}

int saveMatrix(long long int dim, double * mat, const char* fileName) {
    long i;
    FILE *fp;
    fp = fopen(fileName, "w");
    if (fp == NULL) return -1;
    for(i = 0; i<dim ; i++){
	fprintf(fp,"%f\n",mat[i]);
    }
    fclose(fp);
    return 0;
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

