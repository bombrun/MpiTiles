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
    MPI_Init(&argc,&argv);
    int rank, p;
    
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    /* Affecte a rank mon numero de processus */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    const char* profileAB_file_name= "./NormalsAB/profile.txt";
    const char* valuesAB_file_name = "./NormalsAB/values.txt";
    
    const char* profileG_file_name= "./NormalsG/profile.txt";
    const char* valuesG_file_name = "./NormalsG/values.txt";
    
    int i;
    
    int i0,i1; // row numbers of matrix (CABG)' to be processed
    
  
    int* profileAB;
    int profileAB_length,dimensionAB;
    double* matrixAB;
    printf("%d/%d : started\n",rank,p);
    profileAB_length = getNumberOfLine(profileAB_file_name);
    profileAB =(int *)malloc(sizeof(int)*profileAB_length);
    //allocateMatrixInt(&profileAB,profileAB_length);
    readMatrixInt(profileAB,profileAB_file_name);
    dimensionAB = sumVectorInt(profileAB,profileAB_length);
    allocateMatrixDouble(&matrixAB,dimensionAB);
    readMatrixDouble(matrixAB,valuesAB_file_name);
    
   
      int* profileG;
      int profileG_length;
      double* matrixG;
      profileG_length = getNumberOfLine(profileG_file_name);
      allocateMatrixInt(&profileG,profileG_length);
      readMatrixInt(profileG,profileG_file_name);
      allocateMatrixDouble(&matrixG,dimensionAB);
      readMatrixDouble(matrixG,valuesG_file_name);
      printf("%d/%d: number of global %d\n",rank,p,profileG_length);
    
    i0=mpi_get_i0(profileG_length,rank,p);
    i1=mpi_get_i1(profileG_length,rank,p);
    printf("%d/%d: process rows from %d to %d\n",rank,p,i0,i1);
    
    double* matrixCGAB = (double *)malloc(sizeof(double)*(i1-i0)*profileAB_length);
    const char* row_prefix = "./SparseGtAB/row";
    const char* row_sufix = ".txt";
    for(i=i0;i<i1;i++){
      char *row_file_name = (char *)malloc(sizeof(char)*(strlen(row_prefix)+strlen(row_sufix)+5));
      sprintf(row_file_name,"%s%d%s",row_prefix,i,row_sufix);
      setRowWithSparseVectorDouble(matrixCGAB, i-i0, profileAB_length ,row_file_name);
      reduceRhs(matrixAB, matrixCGAB,i-i0,profileAB_length,profileAB);
      //printf("%d/%d : finished %s\n",rank,p,row_file_name);
    }
    //printVectorDouble(matrixCGAB,100);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d/%d: finished computing C-1ABG\n",rank,p);
    
    
    // dispatch blocks of matrixCor between the process
    // each process contains its diagonal block
    // b11	b12	b13	b14
    //		b22	b23	b24
    //			b33	b34
    //				b44
    
    // p1 : b11	b12 b13
    // p2 : b22 b23 b24
    // p3 : b33	b34
    // p4 : b44 b14
    int n_total_tasks = (p*(p+1))/2; //number of tasks to dispatch
    int n_tasks=n_total_tasks/p; // n_total_tasks/p or n_total_tasks/p+1
    int remainTasks=n_total_tasks%p;
    
    int nTasks[p];
    int r;
    int **tasks;
    tasks = (int **)malloc(sizeof(int *)*p);
    for(r=0;r<p;r++){
      nTasks[r] = n_tasks;
      if (r<remainTasks) nTasks[r]+=1;
      tasks[r]=(int *)malloc(sizeof(int)*nTasks[r]);
      for(i=0;i<nTasks[r];i++) tasks[r][i]=(r+i)%p; // rotating permutation, from i=1 the process rank should ask the information at process rank+i mod p
    }
    
    
    int n_com_tasks=n_total_tasks/p+1;
    int recv_tasks[p][n_com_tasks];
    int send_tasks[p][n_com_tasks];
     for(r=0;r<p;r++){
      for(i=0;i<n_com_tasks;i++){
	recv_tasks[r][i]=(r+i)%p; // rotating permutation, from i=1 the process rank should ask the information at process rank+i mod p
	send_tasks[r][i]=(r+p-i)%p;
	//if(rank==0) printf("%d-%d = %d mod %d \n",r,i, (r+p-i)%p,p);
      }
    }
    printf("%d/%d: process %d tasks over %d communication tasks\n",rank,p,nTasks[rank],n_com_tasks);
    
    double** matrixCor = (double **)malloc(sizeof(double*)*nTasks[rank]);
    // start by computing the diagonal block 
    // change for dgemm
    int idim,jdim;
    idim = i1-i0;
    jdim = i1-i0;
    matrixCor[0] = (double *)malloc(sizeof(double)*idim*jdim);
    setBlockMatrix(matrixCor[0],i0,i1,i0,i1,matrixG,profileG_length,profileG);
    
//    dgemm(matrixCGAB,idim,profileAB_length,matrixCGAB,jdim,profileAB_length,matrixCor[0],idim,jdim);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			  idim, jdim, profileAB_length,
			  -1.0, matrixCGAB, idim,
			  matrixCGAB, jdim,
			  0.0,  matrixCor[0], idim );
    saveMatrixBlock(i0,i1,i0,i1,matrixCor[0],"./ReducedBlockMatrixG");
    
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d/%d: finished computing diagonal block of the correction\n",rank,p);
    
    
    // then exchange the information between the process
    // p1 ask p2 
    // p2 ask p3
    // p3 ask p4
    // p4 ask p1
    // compute 
    MPI_Status status,status2; 
    MPI_Request send_request,recv_request;
    int ierr;
    
    // loop over tasks
    for(i=1;i<n_com_tasks;i++){
      /*
      int in;
      int out = rank;
      int buff_s = sizeof(int);
      ierr=MPI_Isend(&out,buff_s,MPI_INTEGER,send_tasks[rank][i],0,MPI_COMM_WORLD,&send_request);   
      ierr=MPI_Irecv(&in,buff_s,MPI_INTEGER,recv_tasks[rank][i],MPI_ANY_TAG,MPI_COMM_WORLD,&recv_request);
      printf("%d/%d: send data to %d and expect data from %d\n",rank,p,send_tasks[rank][i],recv_tasks[rank][i]);
      ierr=MPI_Wait(&recv_request,&status); 
      printf("%d/%d: received %d \n",rank,p,in);
      */
      int buffsize_out = (i1-i0)*profileAB_length;
      int j0 = mpi_get_i0(profileG_length,recv_tasks[rank][i],p);
      int j1 = mpi_get_i1(profileG_length,recv_tasks[rank][i],p);
      printf("%d/%d: i0=%d, i1=%d linked with rank %d j0=%d, j1=%d \n",rank,p,i0,i1,recv_tasks[rank][i],j0,j1);
      jdim =j1-j0;
      int buffsize_in = (j1-j0)*profileAB_length;
      
     
      double* matrixCGABj = (double *)malloc(buffsize_in*sizeof(double));
      ierr=MPI_Isend(matrixCGAB,buffsize_out,MPI_DOUBLE,send_tasks[rank][i],0,MPI_COMM_WORLD,&send_request); 
      printf("%d/%d: send data to %d\n",rank,p,send_tasks[rank][i]);
      ierr=MPI_Irecv(matrixCGABj,buffsize_in,MPI_DOUBLE,recv_tasks[rank][i],MPI_ANY_TAG,MPI_COMM_WORLD,&recv_request);
      ierr=MPI_Wait(&recv_request,&status);
      ierr=MPI_Wait(&send_request,&status2); 
      printf("%d/%d: received data from %d\n",rank,p,recv_tasks[rank][i]);
      if(i<nTasks[rank]){
	matrixCor[i] = (double *)malloc(sizeof(double)*idim*jdim);
	setBlockMatrix(matrixCor[i],i0,i1,j0,j1,matrixG,profileG_length,profileG);
	dgemm(matrixCGAB,idim,profileAB_length,matrixCGABj,jdim,profileAB_length,matrixCor[i],idim,jdim);
	printf("%d/%d: finished computing block %d,%d of the correction\n",rank,p,rank,recv_tasks[rank][i]);
	saveMatrixBlock(i0,i1,j0,j1,matrixCor[i],"./ReducedBlockMatrixG");
	//printVectorDouble(matrixCor[i],10);
	//printMatrixDouble(matrixCor[i],idim,jdim);
	//printf("\n");
      }
      free(matrixCGABj); 
      MPI_Barrier(MPI_COMM_WORLD); // do not start new tasks while there are running tasks
    }
    
     
    // total required memory 2 times ABG + 
    
    
    MPI_Finalize();
    return ierr;
}


