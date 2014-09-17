
/*******************************************************************************
 * a code that compute the eigen value decomposition of a symmetric real matrix stored in an ASCII text file
 * by Alex Bombrun from MKL example
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include "normals.h"

#define OWNDGEMM

int getNumberOfLine(const char* fileName) {
    int i=0;

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        return -1;

    while ((read = getline(&line, &len, fp)) != -1) {
        i++;
    }

    free(line);
    fclose(fp);
    return i;
}

void allocateMatrixDouble(double ** mat,int size) {
    *mat = malloc(size*sizeof(double));
}

void allocateMatrixInt(int ** mat,int size) {
    *mat = malloc(size*sizeof(int));
}

int readMatrixDouble(double * mat, const char* fileName) {
    long i=0;

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        return -1;

    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
        mat[i] = atof(line);
        i++;
    }

    free(line);
    fclose(fp);
    return 0;
}

/**
 * set the row i0 of the matrix mat with the value of the sparse vector stored in filename
 */
int setRowWithSparseVectorDouble(double * mat,int row_i, int row_length, const char* fileName) {
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        return -1;

    int i;
    double val;
    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
      
	i = atoi(strsep(&line, ";"));
	val = atof(strsep(&line, ";"));
	//while ((token = strsep(&line, ";")) != NULL)
	//{
	//  printf("%s\n", token);
	//}
        
        mat[i+row_i*row_length] = val;
        //i++;
    }

    free(line);
    fclose(fp);
    return 0;
}



int readMatrixInt(int * mat, const char* fileName) {
    long i=0;

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        return -1;

    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
        mat[i] = atoi(line);
        i++;
    }

    free(line);
    fclose(fp);
    return 0;
}

int saveMatrixBlock(int i0, int i1, int j0, int j1, double * mat, const char* outputLocation) {
    // check if the output location exists
    struct stat st = {0};
    if (stat(outputLocation, &st) == -1) {
      mkdir(outputLocation, 0700);
    }
    char *file_name = (char *)malloc(sizeof(char)*(strlen(outputLocation)+300));
    sprintf(file_name,"%s/block%d-%d-%d-%d.txt",outputLocation,i0,i1,j0,j1);
    long i;
    FILE *fp;
    fp = fopen(file_name, "w");
    if (fp == NULL) return -1;
    fprintf(fp,"#i0=%d,i1=%d,j0=%d,j1=%d\n",i0,i1,j0,j1);
    int size = (i1-i0)*(j1-j0);
    for(i = 0; i<size ; i++) {
        fprintf(fp,"%f\n",mat[i]);
    }
    fclose(fp);
    return 0;
}

void formatMatrix(double * ma,int dim) {
    int i,j;
    for( i=0; i<dim; i++) {
        for( j=i+1; j<dim; j++) {
            ma[i*dim+j]=0.0;
        }
    }

}


int printMatrixDouble(double * ma,int idim,int jdim) {
    int i,j;
    for( i=0; i<idim; i++) {
        for( j=0; j<jdim; j++) {
            printf( " %e", ma[i*jdim+j] );
        }
        printf( "\n" );
    }
    return 0;
}

int printVectorDouble(double * ma,int dim) {
    int i;
    for( i=0; i<dim; i++) {
        printf( "%i: %e \n",i, ma[i] );
    }
    return 0;
}

int printVectorInt(int * ma,int dim) {
    int i;
    for( i=0; i<dim; i++) {
        printf( "%d", ma[i] );
        printf( "\n" );
    }
    return 0;
}

int sumVectorInt(int *ma,int dim){
   int i,res;
    for( i=0; i<dim; i++) {
        res+= ma[i] ;
    }
    return res;
}

/**
 * return Mij
 */
double get(double * mat, int * sumProfile, int * profile, int i, int j){
  if(i<j) {
    return get(mat,sumProfile,profile,j,i);
  } else {
    int k = i-j;
    if(k>profile[i]) {
      return 0.0;
    } else {
      return mat[sumProfile[i]-profile[i]+k];
    }
  }
  
}

void setCumulative(int * sum,int * vec,int dim){
     int l,total;
     total=0;
      for(l=0;l<dim;l++){
	total+=vec[l];
	sum[l]=total;
      }
}

/**
 * set matOut with the value of matIn
 * assume maout is a sub matrx of matIn in block format i0,i1 x j0,j1
 */
void setBlockMatrix(double *matOut, int i0,int i1,int j0, int j1,double * matIn, int dim, int * profile){
  int * sumProfile = (int*)malloc(sizeof(int)*dim);
  setCumulative(sumProfile,profile,dim);
  int i,j;
  int jdim = j1-j0;
  for(i=i0;i<i1;i++){
    for(j=0;j<j1;j++){
	matOut[(i-i0)*jdim+(j-j0)]=get(matIn,sumProfile,profile,i,j);
    }
  }
}


/**
 * Assume the matrix reduced
 * reduced vector vec
 * */
int reduceRhs(double * ma, double * vec,int i0 ,int dim, int * profile) {
    /*
    double[][] matrix = N.getM();
    	for (int row = 0; row < rhs.length; row++) {
    		int rhsIdx = row;
    		for (int k = 1; k < Math.min(row+1, N.getProfile(row)); k++) {
    			if (rhsIdx-k >= 0) {
    				rhs[rhsIdx] -= matrix[row][k] * rhs[rhsIdx - k];
    			}
    		}
    		if (matrix[row][0] > 0.0) {
    			rhs[rhsIdx] /= matrix[row][0];
    		} else {
    			rhs[row] = 0.0;
    		}
    	}
    	*/
    int i,k,counter;
    double matrix_i_0;
    counter=0;
    for( i=0; i<dim; i++) {
        matrix_i_0 = ma[counter];
        counter++;
        for (k=1; k< profile[i]; k++) {
            vec[i0*dim+i] -= ma[counter] * vec[i0*dim+i-k];
            counter++;
        }
        if (matrix_i_0>0.0) {
            vec[i0*dim+i] /= matrix_i_0;
        } else {
            vec[i0*dim+i] = 0.0;
        }
    }
    return 0;
}

/**
 * assume 
 */
int dgemm(double *matrixA,int n_row_a,int n_col_a, double* matrixB,int n_row_b, int n_col_b, double* matrixC, int n_row_c,int n_col_c){
  int i,j,k;
  double stock;
  for(i=0;i<n_row_a;i++){
     for(j=0;j<n_row_b;j++){
       stock =0;
       for(k=0;k<n_col_a;k++){
	stock += matrixA[i*n_col_a+k]*matrixB[j*n_col_b+k]; //n_col_b == n_col_a
      }
      matrixC[i*n_col_c+j] -= stock;
     }
    }
  return 0;
}


#ifdef OWNDGEMM
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
            const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
            const int K, const double alpha, const double *A,
            const int lda, const double *B, const int ldb,
            const double beta, double *C, const int ldc) {

    int i, j, k;
    double stock;
    
    if ( (TransA == CblasNoTrans) && (TransB == CblasTrans)) {
        for ( i = 0; i < M; i++ ) {
            for ( j = 0; j < N; j++ ) {
	       stock=0;
                for ( k = 0; k < K; k++ )
                {
                    stock += A[ i * K + k ] * B[ j * K + k ];
                }
                C[ i * N + j ]+=alpha* stock;
	    }
	}
    }

}
#endif
