
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "matrixBlockStore.h"

void openStore(FILE** store, int i_block, const char* location){
    // check if the output location exists
    struct stat st = {0};
    if (stat(location, &st) == -1) {
      mkdir(location, 0700);
    }
    char *file_name = malloc(sizeof(char)*(strlen(location)+300));
    sprintf(file_name,"%s/block%d.txt",location,i_block);
    *store = fopen(file_name, "w");
    free(file_name);
}


void readStore(FILE** store, int i_block, const char* location){
    // check if the output location exists
    struct stat st = {0};
    if (stat(location, &st) == -1) {
      mkdir(location, 0700);
    }
    char *file_name = malloc(sizeof(char)*(strlen(location)+300));
    sprintf(file_name,"%s/block%d.txt",location,i_block);
    *store = fopen(file_name, "r");
    free(file_name);
}

void closeStore(FILE* store){
    fclose(store);
}

int saveBlock(int i0, int i1, int j0, int j1, double * mat,FILE* store) {
    long i;
    //if (store == NULL) return -1;
    fprintf(store,"#%d,%d,%d,%d\n",i0,i1,j0,j1);
    int size = (i1-i0)*(j1-j0);
    for(i = 0; i<size ; i++) {
        fprintf(store,"%f\n",mat[i]);
    }
    return 0;
}

/**
 * assume mat allocated
 */
void readNextBlock(int i0, int i1, int j0, int j1, double * mat, FILE* store){
    int i;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    int size = (i1-i0)*(j1-j0);
    for(i = 0; i<size ; i++) {
	if ( (read = getline(&line, &len, store)) != -1) {
	  mat[i] = atof(line);
	}
    }
    free(line);
}  

/**
 * assume dim allocated with length 4
 * read one next of store
 * set dim
 */
int readNextBlockDimension(int* dim, FILE* store){
    char * line= NULL;
    char * ptr = NULL;
    size_t len = 0;
    ssize_t read;
    if ( (read = getline(&line, &len, store)) != -1 ){
	ptr = line;
	strsep(&line, "#");
	dim[0] = atoi(strsep(&line, ","));
	dim[1] = atoi(strsep(&line, ","));
	dim[2] = atoi(strsep(&line, ","));	
	dim[3] = atoi(strsep(&line, ","));
	return 0;
    } else {
      return -1;
    }
    free(ptr);
}