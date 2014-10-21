#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "normals.h"
#include "mpiutil.h"
#include "matrixBlockStore.h"

#include "matrixScalapackStore.h"


#include <float.h>

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);
extern void   Cblacs_gridmap( int* context, int* map, int ld_usermap, int np_row, int np_col);

extern void pdsyev_( char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *info );

//extern void dgsum2d_(int context, char * scope, char* top, int m, int n, double* A, int lda, int rdest, int cdest);

/**
    input c
    assume c = n*n;
    return n
    */
int isqrt(int c);


int saveMatrix(int dim, double * mat, const char* fileName);

/*
 * solve the reduced normal equation
 * should be called after Cblas_pinfo
 * parameters :
 * 	mpye process 
 * 	npe number of processes
 * 	rhs the value of the right hand side of the equation
 * 	n_blocks the number of diagonal block used to compute the reduced normal matrix
 * 	scalapack_size used to compute the reduced normal matrix
 * 	N the dimension of the reduced normal equation
 * 	scaStore_location the location of the Cholesky factor of the reduced normal matrix
 * 	out will contain the solution
 */
double solveRhs(int mype, int npe, double * rhs, int n_blocks, int scalapack_size, int N, const char* scaStore_location,double **out);

/* 
 * 
 * A program that use precomputed cholesky factor of the reduced normal matrix
 * created 23/09/2014
 * author Alex Bombrun
 * 
 * Compile on linux with scalapack and openmpi 
 * mpicc -O1 -o readChol.mpi myReadCholesky.c mpiutil.c normals.c matrixBlockStore.c matrixScalapackStore.c -L/opt/scalapack/lib/ -lscalapack -llapack -lrefblas -lgfortran -lm
 * run with
 * mpirun -n 4 readChol.mpi 4 4
 * 
 * assume that 
 * mpirun -n 4 cholesky.mpi 4
 * was run successfully
 *
*/
int main(int argc, char **argv) {
    int n_blocks;
    int scalapack_size;
    int N, i;
    double * out;
    double * rhs;
    const char* profileG_file_name= "./data/NormalsG/profile.txt";
    const char* scaStore_location ="./data/CholeskyReducedNormals";
    int mype,npe; // rank and total number of process
    Cblacs_pinfo( &mype, &npe );
    N = getNumberOfLine(profileG_file_name); // the dimension of the matrix;
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
    rhs = malloc(sizeof(double)*N);
    for(i=0;i<N;i++){
      rhs[i]=1.0/N;
    }
    saveMatrix(N,rhs,"rhs.txt");
    out = calloc(sizeof(double),N); // local 
    solveRhs(mype,npe,rhs, n_blocks, scalapack_size, N, scaStore_location,&out);
    if (mype==0) saveMatrix(N,out,"sol.txt");
    Cblacs_exit( 0 );
}

double solveRhs(int mype, int npe, double * rhs, int n_blocks, int scalapack_size, int N, const char* scaStore_location, double ** out){
    FILE* scaStore;
    
    int M;
    
    int mp;	// number of rows in the processor grid
    int mla;    // number of rows in the local array
    int mb;     // number of rows in a block
    int np;	// number of columns in the processor grid
    int nla;    // number of columns in the local array
    int nb;     // number of columns in a block
    
   
    
    int idescal[9]; // matrix descriptors
    double *la; // matrix values: al is the local array
    
    int idescbl[9];
    double *lb;
    double normb;
  
    int ierr; // error output 
    int mp_ret, np_ret, myrow, mycol; // to store grid info
    
    int zero=0; // value used for the descriptor initialization
    int one=1; // value used for the descriptor initialization
    
    int  m,n; // matrix A dimensions
    int icon; // scalapack cblacs context
    char normJob, jobz, uplo, trans, diag;
    
    int rdest,cdest;
    jobz= 'N'; uplo='U';
    
  
    
  
    printf("%d/%d: read store\n",mype,npe);
  
    
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
    
    // read the cholesky precomputed factor
    readScalapackStore(&scaStore,myrow,mycol,scaStore_location);
    readLocalMatrix(la,nla,mla,scaStore);
    fclose(scaStore);
    // testing 
    /*
    double test=0.0;
    for(i=0;i<nla*mla;i++){
	test += la[i]*la[i];
    }
    printf("%d/%d: finished loading the matrix, test=%f \n",mype,npe,test);
    */
    printf("%d/%d: cholesky factor loaded \n",mype,npe);
    // set the matrix descriptors
    ierr=0;
    descinit_(idescal, &m, &n  , &mb, &nb , &zero, &zero, &icon, &mla, &ierr); // processor grip id start at 0
    
    ierr=0;
    descinit_(idescbl, &m, &one  , &mb, &nb , &zero, &zero, &icon, &nla, &ierr); // processor grip id start at 0
    lb = calloc(sizeof(double),mla);
    
    // scatter the right hand side 
    // WRONG
    int i,rsrc,h,pr,ii;
    float cr,ll;
    rsrc=0;
    for(i=0;i<N;i++){
      cr = (float)( i/mb );
      h = rsrc+(int)(cr);
      pr = h%np;
      if(myrow==pr & mycol==0){
	ll = (float)( ( i/(np*mb) ) );
	ii = i%mb + (int)(ll)*mb; 
	lb[ii]=rhs[i];
      }
    }
    printf("%d/%d: right hand side scattered \n",mype,npe);
    normb=0;
    pddot_(&n,&normb,lb,&one,&one,idescbl,&one,lb,&one,&one,idescbl,&one); // norm <- b'b
    if (mype==0)  printf("%d/%d: norm rhs, b\'b=%E \n",mype,npe,normb);  
     
    ierr =0;
    diag = 'N';
   
    pdpotrs_(&uplo, &n , &one , la , &one , &one , idescal , lb , &one , &one , idescbl , &ierr); // b<- A-1 b
    
    normb=0;
    pddot_(&n,&normb,lb,&one,&one,idescbl,&one,lb,&one,&one,idescbl,&one); // norm <- b'b
    if (mype==0)  printf("%d/%d: finish solving, norm rhs b\'b=%E \n",mype,npe,normb);  
    
     // scatter the right hand side
    rsrc=0;
    for(i=0;i<N;i++){
      cr = (float)( i/mb );
      h = rsrc+(int)(cr);
      pr = h%np;
      if(myrow==pr & mycol==0){
	ll = (float)( ( i/(np*mb) ) );
	ii = i%mb + (int)(ll)*mb; 
	(*out)[i]=lb[ii];
      }
    }
    const char * scope = "ALL";
    const char * top = "t";//"1-tree";
    rdest = 0;
    cdest = 0;
    dgsum2d_(&icon,scope,top,&one,&n,(*out),&one,&rdest,&cdest); // combine the local solution and send them to process  (0,0)
    printf("%d/%d: out gather \n",mype,npe);  
    free(la);
    Cblacs_gridexit(icon);
}


int saveMatrix(int dim, double * mat, const char* fileName) {
    long i;
    FILE *fp;
    fp = fopen(fileName, "w");
    if (fp == NULL) return -1;
    for(i = 0; i<dim ; i++){
	fprintf(fp,"%.*e\n",DBL_DIG,mat[i]);
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

