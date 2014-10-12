
/**
 * open with writing autorisation the file corresponding to the processor located at myrow,mycol on the CBLAC processor grid
 * if the location directory does not exist create it 
 */
extern void openScalapackStore(FILE** store, int myrow, int mycol, const char* location);
extern void readScalapackStore(FILE** store, int myrow, int mycol, const char* location);

/**
 * save the local matrix lmat with dimension (nla,mla) corresponding to the process used to open the store
 */
extern int  saveLocalMatrix(double * lmat,int nla, int mla, FILE* store);

/**
 * save the matrix descriptor
 * should be call only once for example by the process of rank 0
 */
extern void saveMatrixDescriptor(int * desc, const char* location);
extern void readMatrixDescriptor(int * desc, const char* location);