#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "normals.h"
#include "mpiutil.h"
    
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
    
    int* profileAB = NULL;
    int profileAB_length,dimensionAB;
    double* matrixAB = NULL;
    
    // Normal matrix block G
    const char* profileG_file_name= "./data/NormalsG/profile.txt";
    const char* valuesG_file_name = "./data/NormalsG/values.txt";
      
    int* profileG = NULL;
    int profileG_length, dimensionG;
    double* matrixG = NULL;
  
    int i, j, r; 	 // loop indices
    int i0, i1;  // main row numbers of matrix (CABG)' to be processed
    int j0, j1;  // secondary row numbers
    int idim,jdim;
    
    int rank, p;  // rank and number of mpi process
    
    int n_diag_blocks; // the total number of diagonal blocks
    int i_block, j_block; // the diagonal block index
    
    int ierr = 0; // process error
   
    const char* row_prefix = "./data/SparseGtAB/row";
    const char* row_sufix = ".txt";
    
    double* matrixCor = NULL;
    double* matrixCGABi = NULL;
    double* matrixCGABj = NULL;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n_diag_blocks = p; // one diagonal block per process
    i_block = rank; // the diagonal block index equals rank index
    
    // blocks and tasks distribution over all the process
    int remainTasks=mpi_get_total_blocks(n_diag_blocks)%n_diag_blocks;
    int nTasks[n_diag_blocks]; // number of tasks for each diagonal block : n_total_tasks/p or n_total_tasks/p+1
    int max_tasks_per_diag_block=mpi_get_total_blocks(n_diag_blocks)/n_diag_blocks+1;
    for(r=0;r<n_diag_blocks;r++){
      nTasks[r] = mpi_get_total_blocks(n_diag_blocks)/n_diag_blocks;
      if (r<remainTasks) nTasks[r]+=1;
    }
    printf("%d/%d: inde, started, process diagonal block %d/%d\n",rank,p,i_block,n_diag_blocks);
    printf("%d/%d: process %d tasks\n",rank,p,nTasks[i_block]);
    
    // set up the input matrices
    profileAB_length = getNumberOfLine(profileAB_file_name);
    profileAB = malloc( sizeof(int) * profileAB_length );
    readMatrixInt(profileAB,profileAB_file_name);
    dimensionAB = sumVectorInt(profileAB,profileAB_length);
    allocateMatrixDouble(&matrixAB,dimensionAB);
    readMatrixDouble(matrixAB,valuesAB_file_name);
    
    profileG_length = getNumberOfLine(profileG_file_name);
    allocateMatrixInt(&profileG,profileG_length);
    readMatrixInt(profileG,profileG_file_name);
    dimensionG = sumVectorInt(profileG,profileG_length);
    allocateMatrixDouble(&matrixG,dimensionG);
    readMatrixDouble(matrixG,valuesG_file_name);
    
    printf("%d/%d: number of attitude %d parameters, number of source global %d parameters\n", rank, p, profileAB_length, profileG_length);
        
    // diagonal block
    i0=mpi_get_i0(profileG_length,i_block,n_diag_blocks);
    i1=mpi_get_i1(profileG_length,i_block,n_diag_blocks);
    idim = i1-i0;
    jdim = i1-i0;
    printf("%d/%d: process rows from %d to %d\n",rank,p,i0,i1);
    matrixCGABi = malloc(sizeof(double)*(i1-i0)*profileAB_length);
    for(i=i0;i<i1;i++){
      char *row_file_name = malloc(sizeof(char)*(strlen(row_prefix)+strlen(row_sufix)+5));
      sprintf(row_file_name,"%s%d%s",row_prefix,i,row_sufix);
      setRowWithSparseVectorDouble(matrixCGABi, i-i0, profileAB_length ,row_file_name);
      reduceRhs(matrixAB, matrixCGABi,i-i0,profileAB_length,profileAB);
      free(row_file_name);
    }
    printf("%d/%d: finished computing C-1ABG block (%d,%d) \n",rank,p,i0,i1);
    matrixCor = malloc(sizeof(double)*idim*jdim);
    setBlockMatrix(matrixCor,i0,i1,i0,i1,matrixG,profileG_length,profileG);   
    dgemmAlex(matrixCGABi,idim,profileAB_length,matrixCGABi,jdim,profileAB_length,matrixCor,idim,jdim);
//    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//			  idim, jdim, profileAB_length,
//			  -1.0, matrixCGAB, idim,
//			  matrixCGAB, jdim,
//			  0.0,  matrixCor, idim );
    saveMatrixBlock(i0,i1,i0,i1,matrixCor,"./data/ReducedBlockMatrixG");
    free(matrixCor);
    printf("%d/%d: finished computing block (%d,%d)x(%d,%d) of the correction\n", rank, p, i0, i1, i0, i1);
    
    
    // off diagonal blocks
    for(i=1;i<max_tasks_per_diag_block;i++){
      
      if(i<nTasks[i_block]){
	j_block = mpi_get_diag_block_id(i_block, i, n_diag_blocks);
	j0 = mpi_get_i0(profileG_length, j_block, n_diag_blocks);
	j1 = mpi_get_i1(profileG_length, j_block, n_diag_blocks);
	printf("%d/%d: diagonal block %d (%d,%d) linked with block %d (%d,%d) \n",rank, p, i_block, i0, i1, j_block, j0, j1);
	jdim =j1-j0;
	matrixCGABj = malloc((j1-j0)*profileAB_length*sizeof(double));
	for(j=j0;j<j1;j++){
	  char *row_file_name = malloc(sizeof(char)*(strlen(row_prefix)+strlen(row_sufix)+5));
	  sprintf(row_file_name,"%s%d%s",row_prefix,j,row_sufix);
	  setRowWithSparseVectorDouble(matrixCGABj, j-j0, profileAB_length ,row_file_name);
	  reduceRhs(matrixAB, matrixCGABj,j-j0,profileAB_length,profileAB);
	  free(row_file_name);
	}
	printf("%d/%d: finished computing C-1ABG block (%d,%d) \n",rank,p,j0,j1);
	matrixCor = malloc(sizeof(double)*idim*jdim);
	setBlockMatrix(matrixCor,i0,i1,j0,j1,matrixG,profileG_length,profileG);
	dgemmAlex(matrixCGABi,idim,profileAB_length,matrixCGABj,jdim,profileAB_length,matrixCor,idim,jdim);
	saveMatrixBlock(i0,i1,j0,j1,matrixCor,"./data/ReducedBlockMatrixG");
	printf("%d/%d: finished computing block (%d,%d)x(%d,%d) of the correction\n", rank, p, i0, i1, j0, j1);
	free(matrixCor);
	free(matrixCGABj); 
      }
     
    }
    
    MPI_Finalize(); // the process are independent no blocking
    return ierr;
}


