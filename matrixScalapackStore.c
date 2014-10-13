
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "matrixScalapackStore.h"

void openScalapackStore(FILE** store, int myrow, int mycol, const char* location){
    // check if the output location exists
    struct stat st = {0};
    if (stat(location, &st) == -1) {
      mkdir(location, 0700);
    }
    char *file_name = malloc(sizeof(char)*(strlen(location)+300));
    sprintf(file_name,"%s/blockR%dC%d.txt",location,myrow,mycol);
    *store = fopen(file_name, "w");
    free(file_name);
}


void readScalapackStore(FILE** store, int myrow, int mycol, const char* location){
    // check if the output location exists
    struct stat st = {0};
    if (stat(location, &st) == -1) {
      printf("ERROR: the location %s does not exists",location);
      //exit(-1);
    }
    char *file_name = malloc(sizeof(char)*(strlen(location)+300));
    sprintf(file_name,"%s/blockR%dC%d.txt",location,myrow,mycol);
    *store = fopen(file_name, "r");
    free(file_name);
}


int saveLocalMatrix(double* lmat,int nla, int mla, FILE* store) {
    long i;
    for(i = 0; i<(nla*mla) ; i++) {
        fprintf(store,"%f\n",lmat[i]);
    }
    return 0;
}



int readLocalMatrix(double* lmat,int nla, int mla, FILE* store) {
    int i;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int size = nla*mla;
    for(i = 0; i<size ; i++) {
	if ( (read = getline(&line, &len, store)) != -1) {
	  lmat[i] = atof(line);
	}
    }
    free(line);
    return 0;
}


void saveMatrixDescriptor(int * desc, const char* location){
     // check if the output location exists
    struct stat st = {0};
    if (stat(location, &st) == -1) {
      mkdir(location, 0700);
    }
    char *file_name = malloc(sizeof(char)*(strlen(location)+300));
    sprintf(file_name,"%s/descriptor.txt",location);
    FILE* store = fopen(file_name, "w");
    free(file_name);
    int i;
    for(i = 0; i<9 ; i++) {
        fprintf(store,"%d\n",desc[i]);
    }
    fclose(store);
}

void readMatrixDescriptor(int * desc, const char* location){
    char *file_name = malloc(sizeof(char)*(strlen(location)+300));
    sprintf(file_name,"%s/descriptor.txt",location);
    FILE* store = fopen(file_name, "r");
    int i;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int size = 9;
    for(i = 0; i<size ; i++) {
	if ( (read = getline(&line, &len, store)) != -1) {
	  desc[i] = atoi(line);
	}
    }
    free(line);
    fclose(store);
}
