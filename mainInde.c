#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "normals.h"
#include "mpiutil.h"


/**
 * return the total number of diagonal block
 * if not specify return p
 */
int get_n_blocks(int argc, char **argv, int p);

/**
 *  diagonal block distribution over all the process
 */
int get_n_pTasks(int p, int rank, int n_blocks);
/**
 *  tasks distribution over all the diagonal blocks
 */
int get_n_blockTasks(int i_block, int n_blocks);

/**
 * NOT WORKING
 */
int setMatrixNormal(const char* location,int* profile, double* values);

void reduce(int length, int* profile, double* values, double* matrix, int i0, int i1);

/* Main program 
  a straight forward implementation to read the input matrices and build the reduced normal matrix for the global block
  created 7/09/2014
  author Alex Bombrun
  see TN APB-009
 */
int main(int argc, char **argv) {
  
    // Normal matrix block AB
    const char* profileAB_file_name= "./data/NormalsAB/profile.txt";
    const char* valuesAB_file_name = "./data/NormalsAB/values.txt";
    //const char* AB_location = ".data/NormalsAB";
    
    int* profileAB = NULL;
    int profileAB_length,dimensionAB;
    double* matrixAB = NULL;
    
    // Normal matrix block G
    const char* profileG_file_name= "./data/NormalsG/profile.txt";
    const char* valuesG_file_name = "./data/NormalsG/values.txt";
    //const char* G_location = ".data/NormalsG";
      
    int* profileG = NULL;
    int profileG_length, dimensionG;
    double* matrixG = NULL;
  
    int i, j, t; 	 // loop indices
    int i0, i1;  // main row numbers of matrix (CABG)' to be processed
    int j0, j1;  // secondary row numbers
    int idim,jdim;
    
    int rank, p;  // rank and number of mpi process
    
    int n_blocks; // block size or the number of diagonal blocks
    int i_block, j_block; // the diagonal block index
    
    int n_pTasks, n_blockTasks;
    
    int ierr = 0; // process error
   
    const char* row_prefix = "./data/SparseGtAB/row";
    const char* row_sufix = ".txt";
    
    double* matrixCor = NULL;
    double* matrixCGABi = NULL;
    double* matrixCGABj = NULL;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    n_blocks = get_n_blocks(argc,argv,p);
    n_pTasks = get_n_pTasks(p,rank,n_blocks);
  
    // set up the input matrices
    profileAB_length = getNumberOfLine(profileAB_file_name);
    profileAB = malloc( sizeof(int) * profileAB_length );
    readMatrixInt(profileAB,profileAB_file_name);
    dimensionAB = sumVectorInt(profileAB,profileAB_length);
    matrixAB = calloc(dimensionAB,sizeof(double));
    readMatrixDouble(matrixAB,valuesAB_file_name);
    
    profileG_length = getNumberOfLine(profileG_file_name);
    profileG = malloc( sizeof(int) * profileG_length );
    readMatrixInt(profileG,profileG_file_name);
    dimensionG = sumVectorInt(profileG,profileG_length);
    matrixG = calloc(dimensionG,sizeof(double));
    readMatrixDouble(matrixG,valuesG_file_name);
      
    //profileAB_length = setMatrixNormal(AB_location,profileAB,matrixAB);
    //profileG_length  = setMatrixNormal(G_location,profileG,matrixG);
    
    printf("%d/%d: inde, number of attitude parameters=%d , number of source global parameters=%d , number of blocks=%d\n", rank, p, profileAB_length, profileG_length,n_blocks);

// Memory test
//     i0=mpi_get_i0(profileG_length,0,n_blocks);
//     i1=mpi_get_i1(profileG_length,0,n_blocks);
//     matrixCGABi = calloc(2*(i1-i0)*profileAB_length,sizeof(double));
//    matrixCGABj = calloc(2*(i1-i0)*profileAB_length,sizeof(double));
//     matrixCor = calloc(4*(i1-i0)*(i1-i0),sizeof(double));
     
   for(t=0;t<n_pTasks;t++){
      i_block = (rank+t*p); // the diagonal block index depends on the rank index and the task index
      n_blockTasks = get_n_blockTasks(i_block,n_blocks);
      // diagonal block
      i0=mpi_get_i0(profileG_length,i_block,n_blocks);
      i1=mpi_get_i1(profileG_length,i_block,n_blocks);
      printf("%d/%d: block %d/%d (%d,%d), %d tasks\n",rank,p,i_block,n_blocks,i0,i1,n_blockTasks);
     
      idim = i1-i0;
      jdim = i1-i0;
      printf("%d/%d: process rows from %d to %d\n",rank,p,i0,i1);
      matrixCGABi = calloc((i1-i0)*profileAB_length,sizeof(double));
      reduce(profileAB_length, profileAB, matrixAB, matrixCGABi, i0, i1);
      
      matrixCor = calloc(idim*jdim,sizeof(double));
      setBlockMatrix(matrixCor,i0,i1,i0,i1,matrixG,profileG_length,profileG);   
      dgemmAlex(matrixCGABi,idim,profileAB_length,matrixCGABi,jdim,profileAB_length,matrixCor,idim,jdim);
      saveMatrixBlock(i0,i1,i0,i1,matrixCor,"./data/ReducedBlockMatrixG");
      free(matrixCor);
      printf("%d/%d: finished computing block (%d,%d)x(%d,%d) of the correction\n", rank, p, i0, i1, i0, i1);
      
      for(i=1;i<n_blockTasks;i++){  // off diagonal blocks
	  j_block = mpi_get_diag_block_id(i_block, i, n_blocks);
	  j0 = mpi_get_i0(profileG_length, j_block, n_blocks);
	  j1 = mpi_get_i1(profileG_length, j_block, n_blocks);
	  printf("%d/%d: block %d (%d,%d) linked with block %d (%d,%d) \n",rank, p, i_block, i0, i1, j_block, j0, j1);
	  jdim =j1-j0;
	  matrixCGABj = calloc((j1-j0)*profileAB_length,sizeof(double));
	  reduce(profileAB_length, profileAB, matrixAB, matrixCGABj, j0, j1);
      
	  matrixCor = calloc(idim*jdim,sizeof(double));
	  setBlockMatrix(matrixCor,i0,i1,j0,j1,matrixG,profileG_length,profileG);
	  dgemmAlex(matrixCGABi,idim,profileAB_length,matrixCGABj,jdim,profileAB_length,matrixCor,idim,jdim);
	  saveMatrixBlock(i0,i1,j0,j1,matrixCor,"./data/ReducedBlockMatrixG");
	  printf("%d/%d: block %d (%d,%d) linked with block %d (%d,%d) finished \n",rank, p, i_block, i0, i1, j_block, j0, j1);
	  free(matrixCor);
	  free(matrixCGABj);
	  //MPI_Barrier(MPI_COMM_WORLD); // to prevent memory acess to the same file -> useless does not improve things

      }
      free(matrixCGABi); 
    }
    MPI_Finalize(); // the process are independent no blocking
    return ierr;
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

int setMatrixNormal(const char* location,int* profile, double* values){
  // WOULD BE NICE BUT IT IS NOT WORKING, allocation setting in the same function seems to be tricky in C
    char* profile_file_name=  malloc(sizeof(char)*(strlen(location)+strlen("profile.txt")));
    char* values_file_name =  malloc(sizeof(char)*(strlen(location)+strlen("values.txt")));
    int profile_length, size;
    
    sprintf(profile_file_name,"%s/%s",location,"profile.txt");
    sprintf(values_file_name,"%s/%s",location,"values.txt");
 
    profile_length = getNumberOfLine(profile_file_name);
    profile = calloc( profile_length ,sizeof(int) );
    readMatrixInt(profile,profile_file_name);
    
    size = sumVectorInt(profile,profile_length);
    values = calloc( size, sizeof(double));
    readMatrixDouble(values,values_file_name);
     
    free(profile_file_name);
    free(values_file_name);
    return profile_length;
}


void reduce(int length, int* profile, double* values, double* matrix, int i0, int i1){
      const char* row_prefix = "./data/SparseGtAB/row";
      const char* row_sufix = ".txt";
      int i;
      for(i=i0;i<i1;i++){
	char *row_file_name = malloc(sizeof(char)*(strlen(row_prefix)+strlen(row_sufix)+5));
	sprintf(row_file_name,"%s%d%s",row_prefix,i,row_sufix);
	setRowWithSparseVectorDouble(matrix, i-i0, length ,row_file_name); // MEMORY PROBLEM HERE! : if commented the memory size is stable!!!
	reduceRhs(values, matrix,i-i0,length,profile);
	free(row_file_name);
      } 
}