// #define OWNDGEMM
#include <cblas.h>

/* Auxiliary routines prototypes */
extern int getNumberOfLine(const char* fileName);
extern void allocateMatrixDouble(double ** mat,int dimension);
extern void allocateMatrixInt(int ** mat,int dimension);
extern int readMatrixDouble(double * mat,const char* fileName);
extern int readMatrixInt(int * mat,const char* fileName);

/**
 * set the row i0 of the matrix mat with the value of the sparse vector stored in filename
 * mat is stored in row major format,  
 */
extern int setRowWithSparseVectorDouble(double * mat,int row_i, int row_length, const char* fileName);
extern int saveMatrixBlock(int i0, int i1, int j0, int j1, double * mat, const char* outputLocation);
extern void formatMatrix(double * ma,int dim);
extern int printMatrixDouble(double * ma,int idim,int jdim);
extern int printVectorDouble(double * ma,int dim);
extern int printVectorInt(int * ma,int dim);
extern int sumVectorInt(int *ma,int dim);

/**
 * Assume the matrix is a cholesky factor C, compute x, such that C x = vec
 * the value of x is stored in vec
 * */
extern int reduceRhs(double * ma,double * vec,int i, int dim, int * profile);

extern int dgemm(double *matrixA,int n_row_a,int n_col_a, double* matrixB,int n_row_b, int n_col_b, double* matrixC, int n_row_c,int n_col_c);

/**
 * set matOut with the value of matIn
 * assume maout is a sub matrx of matIn in block format i0,i1 x j0,j1
 */
extern void setBlockMatrix(double *matOut, int i0,int i1,int j0, int j1,double * matIn, int dim, int * profile);


extern void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
            const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
            const int K, const double alpha, const double *A,
            const int lda, const double *B, const int ldb,
            const double beta, double *C, const int ldc);