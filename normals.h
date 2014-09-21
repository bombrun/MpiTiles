/* Auxiliary routines prototypes */
extern int getNumberOfLine(const char* fileName);
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

/**
 * a simple dgemm for test
 */
extern int dgemmAlex(double *matrixA,int n_row_a,int n_col_a, double* matrixB,int n_row_b, int n_col_b, double* matrixC, int n_row_c,int n_col_c);

/**
 * set matOut with the value of matIn
 * assume maout is a sub matrx of matIn in block format i0,i1 x j0,j1
 */
extern void setBlockMatrix(double *matOut, int i0,int i1,int j0, int j1,double * matIn, int dim, int * profile);
extern double getMij(double * mat, int * sumProfile, int * profile, int i, int j);
extern void setCumulative(int * sum,int * vec,int dim);
/**
 * compare matOut with the value of matIn
 * assume maout is a sub matrx of matIn in block format i0,i1 x j0,j1
 * assume matIn is in row major format
 * return the number of entries ij where (matOut_ij-matIn_ij)^2 > tol
 */
int compareBlockMatrix(double *matOut, int i0,int i1,int j0, int j1,double * matIn, int idim, int jdim, double tol);