
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

void closeStore(FILE** store){
    fclose(*store);
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