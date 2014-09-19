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
    
    const char* referenceMatrix_file_name = "./data/reducedNormalMatrix.txt"; 
    double* referenceMatrix = NULL;
    
    int* profileG = NULL;
    int profileG_length, dimensionG;
    double* matrixG = NULL;
    
    int i, r; 	 // loop indices
    int i0, i1;  // main row numbers of matrix (CABG)' to be processed
    int j0, j1;  // secondary row numbers
    int idim,jdim;
    
    int rank, p;  // rank and number of mpi process
    
    int n_diag_blocks; // the total number of diagonal blocks
    int i_diag_block; // the diagonal block index
    
    int ierr = 0; // process error
    
    int test =0;
    
    MPI_Status status,status2; 
    MPI_Request send_request,recv_request;
   
    int  buffsize_in, buffsize_out; // buffer size used for MPI communications
    
    double** matrixCor = NULL;
    double* matrixCGABi = NULL;
    double* matrixCGABj = NULL;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n_diag_blocks = p; // one diagonal block per process
    i_diag_block = rank; // the diagonal block index equals rank index
    printf("%d/%d: started, process diagonal block %d/%d\n",rank,p,i_diag_block,n_diag_blocks);
    
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
    
    referenceMatrix = malloc( sizeof(double) * profileG_length * profileG_length);
    readMatrixDouble(referenceMatrix,referenceMatrix_file_name);
    
    
    i0=mpi_get_i0(profileG_length,rank,p);
    i1=mpi_get_i1(profileG_length,rank,p);
    printf("%d/%d: process rows from %d to %d\n",rank,p,i0,i1);
    
    matrixCGABi = malloc(sizeof(double)*(i1-i0)*profileAB_length);
    const char* row_prefix = "./data/SparseGtAB/row";
    const char* row_sufix = ".txt";
    for(i=i0;i<i1;i++){
      char *row_file_name = malloc(sizeof(char)*(strlen(row_prefix)+strlen(row_sufix)+5));
      sprintf(row_file_name,"%s%d%s",row_prefix,i,row_sufix);
      setRowWithSparseVectorDouble(matrixCGABi, i-i0, profileAB_length ,row_file_name);
      reduceRhs(matrixAB, matrixCGABi,i-i0,profileAB_length,profileAB);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d/%d: finished computing C-1ABG\n",rank,p);
      
    
    int n_tasks=mpi_get_total_blocks(n_diag_blocks)/n_diag_blocks; // n_total_tasks/p or n_total_tasks/p+1
    int remainTasks=mpi_get_total_blocks(n_diag_blocks)%n_diag_blocks;
    
    int nTasks[n_diag_blocks];
    for(r=0;r<p;r++){
      nTasks[r] = n_tasks;
      if (r<remainTasks) nTasks[r]+=1;
    }
    
    int n_com_tasks=mpi_get_total_blocks(n_diag_blocks)/n_diag_blocks+1;
    int recv_tasks[n_diag_blocks][n_com_tasks];
    int send_tasks[n_diag_blocks][n_com_tasks];
     for(r=0;r<n_diag_blocks;r++){
      for(i=0;i<n_com_tasks;i++){
	recv_tasks[r][i]=(r+i)%n_diag_blocks; // rotating permutation, from i=1 the process rank should ask the information at process rank+i mod p
	send_tasks[r][i]=(r+n_diag_blocks-i)%n_diag_blocks;
      }
    }
    printf("%d/%d: process %d tasks over %d communication tasks\n",rank,p,nTasks[rank],n_com_tasks);
    
    matrixCor = malloc(sizeof(double*)*nTasks[i_diag_block]);
    // start by computing the diagonal block 
    // change for dgemm
    
    idim = i1-i0;
    jdim = i1-i0;
    matrixCor[0] = malloc(sizeof(double)*idim*jdim);
    setBlockMatrix(matrixCor[0],i0,i1,i0,i1,matrixG,profileG_length,profileG);
    
//   int * sumProfile = malloc(sizeof(int)*profileG_length);
//    setCumulative(sumProfile,profileG,profileG_length);
//    printf("%d/%d: block (0,0)=%f set to %f\n",rank,p, getMij(matrixG,sumProfile,profileG,i0,i0), matrixCor[0][0]);
 //   printf("%d/%d: block (5,5)=%f set to %f at %d\n",rank,p, getMij(matrixG,sumProfile,profileG,i0+5,i0+5), matrixCor[0][5*jdim+5],5*jdim+5);
    
    
    dgemmAlex(matrixCGABi,idim,profileAB_length,matrixCGABi,jdim,profileAB_length,matrixCor[0],idim,jdim);
//    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//			  idim, jdim, profileAB_length,
//			  -1.0, matrixCGAB, idim,
//			  matrixCGAB, jdim,
//			  0.0,  matrixCor[0], idim );
    saveMatrixBlock(i0,i1,i0,i1,matrixCor[0],"./data/ReducedBlockMatrixG");
    test = compareBlockMatrix(matrixCor[0],i0,i1,i0,i1,referenceMatrix,profileG_length,profileG_length,0.1);
    if(test>0)  printf("%d/%d: ERROR the computed matrix is not equal to the reference matrix \n",rank,p);
    
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d/%d: finished computing diagonal block of the correction\n",i_diag_block,p);
    

    // loop over tasks
    for(i=1;i<n_com_tasks;i++){
      buffsize_out = (i1-i0)*profileAB_length;
      j0 = mpi_get_i0(profileG_length,recv_tasks[i_diag_block][i],p);
      j1 = mpi_get_i1(profileG_length,recv_tasks[i_diag_block][i],p);
      printf("%d/%d: diagonal block %d (%d,%d) linked with block %d (%d,%d) \n",rank, p, i_diag_block, i0, i1, recv_tasks[i_diag_block][i], j0, j1);
      jdim =j1-j0;
      buffsize_in = (j1-j0)*profileAB_length;
           
      matrixCGABj = malloc(buffsize_in*sizeof(double));
      ierr=MPI_Isend(matrixCGABi,buffsize_out,MPI_DOUBLE,send_tasks[rank][i],0,MPI_COMM_WORLD,&send_request); 
      printf("%d/%d: send data to %d\n",rank,p,send_tasks[rank][i]);
      ierr=MPI_Irecv(matrixCGABj,buffsize_in,MPI_DOUBLE,recv_tasks[rank][i],MPI_ANY_TAG,MPI_COMM_WORLD,&recv_request);
      ierr=MPI_Wait(&recv_request,&status);
      ierr=MPI_Wait(&send_request,&status2); 
      printf("%d/%d: received data from %d\n",rank,p,recv_tasks[rank][i]);
      if(i<nTasks[rank]){
	matrixCor[i] = malloc(sizeof(double)*idim*jdim);
	setBlockMatrix(matrixCor[i],i0,i1,j0,j1,matrixG,profileG_length,profileG);
	dgemmAlex(matrixCGABi,idim,profileAB_length,matrixCGABj,jdim,profileAB_length,matrixCor[i],idim,jdim);
	printf("%d/%d: finished computing block %d,%d of the correction\n",rank,p,rank,recv_tasks[rank][i]);
	saveMatrixBlock(i0,i1,j0,j1,matrixCor[i],"./data/ReducedBlockMatrixG");
	
	test = compareBlockMatrix(matrixCor[0],i0,i1,j0,j1,referenceMatrix,profileG_length,profileG_length,0.1);
	if(test>0)  printf("%d/%d: ERROR the computed matrix is not equal to the reference matrix \n",rank,p);
      }
      free(matrixCGABj); 
      MPI_Barrier(MPI_COMM_WORLD); // do not start new tasks while there are running tasks
    }
    
    MPI_Finalize();
    return ierr;
}


