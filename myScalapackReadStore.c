#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "normals.h"
#include "mpiutil.h"
#include "matrixBlockStore.h"

#include "matrixScalapackStore.h"

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);
extern void   Cblacs_gridmap( int* context, int* map, int ld_usermap, int np_row, int np_col);

extern void pdsyev_( char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *info );

/**
    input c
    assume c = n*n;
    return n
    */
int isqrt(int c);
int saveMatrix(long long int dim, double * mat, const char* fileName);

static int max( int a, int b ){
        if (a>b) return(a); else return(b);
}
static int min( int a, int b ){
        if (a<b) return(a); else return(b);
}

/**
 * set the local array store in la with dimenstions mla x nla with the values of matA woth dimenstion m x n
 * mb x nb are the block dimensions
 * myrow, mycol are the current process position in the grid of dimensions mp x np
 * 
 * Compile on linux with scalapack and openmpi 
 * mpicc -O1 -o eigen.exe myScalapackReadStore.c mpiutil.c normals.c matrixBlockStore.c matrixScalapackStore.c -L/opt/scalapack/lib/ -lscalapack -llapack -lrefblas -lgfortran -lm
 * run with
 * mpirun -n 4 eigen.exe 4 4
 * 
 * assume that 
 * mpirun -n 4 bigmatrix.mpi 4
 * was run successfully
 */
void setLocalArray(double *matA, int m, int n, double *la, int mla, int nla,int mb, int nb, int myrow, int mycol, int mp, int np);

//void descinit_(int * idescal, int* m, int *n  , int * mb, int *nb , int *ze, int *ze, int *icon, int * mla, int *ierr);


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
    FILE* scaStore;
    
    int N , M;
    int i, j;
    
    int n_blocks;
    int scalapack_size;
    int NB, MB;
    int i_block, j_block;
    int dim[4];
    double * mat;  // local matrix block use for reading
        
    int t, t_block;
    
    const char* profileG_file_name= "./data/NormalsG/profile.txt";
    const char* store_location = "./data/ReducedNormals";
    const char* scaStore_location ="./data/CholeskyReducedNormals";
    
    int mp;	 // number of rows in the processor grid
    int mla;   // number of rows in the local array
    int mb;    // number of rows in a block
    int np;	 // number of columns in the processor grid
    int nla;   // number of columns in the local array
    int nb;    // number of columns in a block
    
    int mype,npe; // rank and total number of process
    
    int idescal[9]; // matrix descriptors
    double *la; // matrix values: al is the local array
    
    int idesczl[9]; // matrix descriptors
    double *lz; // matrix values: al is the local array
    
    double *w;
   
    int ierr; // error output 
    int mp_ret, np_ret, myrow, mycol; // to store grid info
    
    int zero=0; // value used for the descriptor initialization
    int one=1; // value used for the descriptor initialization
    
    int  m,n; // matrix A dimensions
    double norm, cond;
    double *work = NULL;
    double * work2 = NULL;
    int *iwork = NULL;
    int lwork, liwork;


     float ll,mm,cr,cc;
      int ii,jj,pr,pc,h,g; // ii,jj coordinates of local array element
      int rsrc=0,csrc=0; // assume that 0,0 element should be stored in the 0,0 process
      int n_b = 1;
      int index;
    int icon; // scalapack cblacs context
    char job, jobz, uplo;
    
    double MPIt1, MPIt2, MPIelapsed;
    
    jobz= 'N'; uplo='U';
    Cblacs_pinfo( &mype, &npe );
    
     if (argc == 3) {
	//printf("%s %s %s\n", argv[0], argv[1], argv[2]);
	n_blocks= (int) strtol(argv[1], NULL, 10);
	scalapack_size= (int) strtol(argv[2], NULL, 10);
     } else {
	printf("Usage: expect 2 integers \n");
	printf(" 1 : the number of diagonal blocks \n");
	printf(" 2 : scalapack number to define block size (assume n is divisible by sqrt(p) and that n/sqrt(p) is divisible by this number)\n");
	exit( -1);
     }
    
  

    printf("%d/%d: read store\n",mype,npe);
   
    N = getNumberOfLine(profileG_file_name); // the dimension of the matrix;
    M = N; // square matrix
    
    m=M; //mla*mp;
    n=N; //nla*np;
   
    np = isqrt(npe); // assume that the number of process is a square
    mp = np; // square grid
    
    mla = m/mp; // assume that the matrix dimension if a multiple of the process grid dimension
    nla = n/np;
    
    mb = mla/scalapack_size; // assume that the dimension of the matrix is a multiple of the number of the number of diagonal blocks
    nb = nla/scalapack_size;
    
    // init CBLACS
    Cblacs_get( -1, 0, &icon );
    Cblacs_gridinit( &icon,"c", mp, np ); 
    Cblacs_gridinfo( icon, &mp_ret, &np_ret, &myrow, &mycol);
    
   

    // allocate local matrix
    la=malloc(sizeof(double)*mla*nla);
    printf("%d/%d: full matrix (%d,%d), local matrix (%d,%d), processor grid (%d,%d), block (%d,%d) \n", mype, npe, m, n, mla, nla, np, mp, mb, nb);
    
    for(i_block=0;i_block<n_blocks;i_block++){
      printf("%d/%d: process store block %d \n", mype, npe, i_block);
      readStore(&store,i_block,store_location);
      t_block = 0;
      while(readNextBlockDimension(dim,store)!=-1) { // loop B over all block tasks
	j_block = mpi_get_diag_block_id(i_block, t_block, n_blocks);
	mat = malloc((dim[1]-dim[0])*(dim[3]-dim[2]) * sizeof(double));         
	readNextBlock(dim[0],dim[1],dim[2],dim[3],mat,store);
//	printf("%d/%d: read block (%d,%d) with global indices (%d,%d,%d,%d) \n",mype, npe, i_block,j_block,dim[0],dim[1],dim[2],dim[3]);
	
	NB = dim[1]-dim[0];
	MB = dim[3]-dim[2];
	for(i = dim[0];i<dim[1];i++){
	  for(j = dim[2];j<dim[3];j++){
	      //matA[i*M+j] = mat[(i-dim[0])*MB+(j-dim[2])];
	     // finding out which pe gets this i,j element
              cr = (float)( i/mb );
              h = rsrc+(int)(cr);
              pr = h%np;
              cc = (float)( j/mb );
              g = csrc+(int)(cc);
              pc = g%mp;
	      // check if process should get this element
              if (myrow == pr && mycol==pc){
		  // ii = x + l*mb
		  // jj = y + m*nb
                  ll = (float)( ( i/(np*mb) ) );  // thinks seems to be mixed up does not matter as long as the matrix, the block and the grid is symmetric
                  mm = (float)( ( j/(mp*nb) ) );
                  ii = i%mb + (int)(ll)*mb;
                  jj = j%nb + (int)(mm)*nb;
                  index=jj*mla+ii;   // seems to be the transpose !?
		  //if(index<0) printf("%d/%d: negative index (%d,%d) \n",mype,npe,i,j);
		  //if(index>=mla*nla) printf("%d/%d: too large index (%d,%d) \n",mype,npe,i,j);
                  la[index] = mat[(i-dim[0])*MB+(j-dim[2])];
              }
	  }
	}
	// transpose
	if(j_block != i_block){
	  for(i = dim[0];i<dim[1];i++){
	    for(j = dim[2];j<dim[3];j++){
	      //matA[j*M+i] = mat[(i-dim[0])*MB+(j-dim[2])];
	       // finding out which pe gets this j,i element
              cr = (float)( j/mb );
              h = rsrc+(int)(cr);
              pr = h%np;
              cc = (float)( i/mb );
              g = csrc+(int)(cc);
              pc = g%mp;
	      // check if process should get this element
              if (myrow == pr && mycol==pc){
		  // ii = x + l*mb
		  // jj = y + m*nb
                  ll = (float)( ( j/(np*mb) ) );  // thinks seems to be mixed up does not matter as long as the matrix, the block and the grid is symmetric
                  mm = (float)( ( i/(mp*nb) ) );
                  ii = j%mb + (int)(ll)*mb;
                  jj = i%nb + (int)(mm)*nb;
                  index=jj*mla+ii;   // seems to be the transpose !?
		  //if(index<0) printf("%d/%d: negative index (%d,%d) \n",mype,npe,i,j);
		  //if(index>=mla*nla) printf("%d/%d: too large index (%d,%d) \n",mype,npe,i,j);

                  la[index] = mat[(i-dim[0])*MB+(j-dim[2])];
              }
	    }
	  } 
	}
	
	free(mat);
	t_block++;
      }
      closeStore(store);
    }
    
    
    printf("%d/%d: finished scaterring the matrix \n",mype,npe);
    
    printf("%d/%d: start computing \n",mype,npe);
       // set the matrix descriptor
    ierr=0;
    descinit_(idescal, &m, &n  , &mb, &nb , &zero, &zero, &icon, &mla, &ierr);
    if (mype==0) saveMatrixDescriptor(idescal, scaStore_location);
    
    ierr = 0;
    // compute norm 1 of the reduced normal matrix
    /*
    lwork = 2*mla+2*nla;
    work = malloc(sizeof(double)*lwork);
    job = '1';
    norm = pdlansy_(&job, &uplo, &n, la, &one, &one, idescal, work); 
    printf("%d/%d: norm %f \n",mype,npe,norm);
    free(work);
    */
    
    ierr = 0;
    // compute the cholesky decomposition 
    /*  
    pdpotrf_(&uplo,&n,la,&one,&one,idescal,&ierr);
    printf("%d/%d: finished computing cholesky factor\n",mype,npe);
    openScalapackStore(&scaStore,myrow,mycol,scaStore_location);
    saveLocalMatrix(la,nla,mla,scaStore);
    */
    
    ierr = 0;
    // compute the eigen values
    jobz= 'N'; uplo='U'; // with N z is ignored
    descinit_(idesczl, &m, &n  , &mb, &nb , &zero, &zero, &icon, &mla, &ierr);
    lz = malloc(sizeof(double)*mla*nla);
    w = malloc(sizeof(double)*m);
    lwork = -1;
    work = malloc(sizeof(double)*2);
    pdsyev_( &jobz, &uplo, &n, la, &one, &one, idescal, w, lz, &one, &one, idesczl, work, &lwork, &ierr);   // only compute lwork
    //pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info );
    lwork= (int) work[0];
    free(work);
    work = (double *)calloc(lwork,sizeof(double)) ;
    //MPIt1 = MPI_Wtime();
    pdsyev_( &jobz, &uplo, &n, la, &one, &one, idescal, w, lz, &one, &one, idesczl, work, &lwork, &ierr);   // compute the eigen values
    //MPIt2 = MPI_Wtime();
    //MPIelapsed=MPIt2-MPIt1;
    
    if (mype == 0) {
	saveMatrix(n,w,"eigenvalues.txt");
	//printf("%d/%d: finished job in %8.2fs\n",mype,npe,MPIelapsed); // not working
    }

    ierr = 0;
    // compute the conditioner number assume that the norm and the cholesky decomposition have been computed
    /*
    lwork = 2*mla+3*nla;
    work2 = malloc(sizeof(double)*lwork);
    liwork = 2*mla;
    iwork = malloc(sizeof(int)*liwork);
    pdpocon_("L",&n,la,&one,&one,idescal,&norm,&cond,work2,&lwork,iwork,&liwork,&ierr);
    printf("%d/%d: condition number %f \n",mype,npe,cond);
    */
    
    free(la);
    Cblacs_gridexit(icon);
    Cblacs_exit( 0 );
    return 0;
}

void setLocalArray(double *matA, int m, int n, double *la, int mla, int nla,int mb, int nb, int myrow, int mycol, int np, int mp){
     
    
     float ll,mm,cr,cc;
      int ii,jj,i,j,pr,pc,h,g; // ii,jj coordinates of local array element
      int rsrc=0,csrc=0; // assume that 0,0 element should be stored in the 0,0 process
      int n_b = 1;
      int index;
     
      
      for (i=0;i<m;i++) {
      	for (j=0;j<n;j++) {
	      // finding out which pe gets this i,j element
              cr = (float)( i/mb );
              h = rsrc+(int)(cr);
              pr = h%np;
              cc = (float)( j/mb );
              g = csrc+(int)(cc);
              pc = g%mp;
	      // check if process should get this element
              if (myrow == pr && mycol==pc){
		  // ii = x + l*mb
		  // jj = y + m*nb
                  ll = (float)( ( i/(np*mb) ) );  // thinks seems to be mixed up does not matter as long as the matrix, the block and the grid is symmetric
                  mm = (float)( ( j/(mp*nb) ) );
                  ii = i%mb + (int)(ll)*mb;
                  jj = j%nb + (int)(mm)*nb;
                  index=jj*mla+ii;   // seems to be the transpose !?
                  la[index] = matA[i*n+j];
              }
          }
      }

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
int isqrt(int c)
{
    int n = 0;
    while( n*n < c)
    {
        n++;
    }
    return n;

}

