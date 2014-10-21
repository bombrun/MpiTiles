// this package implements a simple storage for scalapack distributed matrix

/**
 * open with writing autorisation the file store corresponding to the processor located at myrow,mycol on the CBLAC processor grid
 * if the location directory does not exist create it 
 * each process of the grids should call the functions
 */
extern void openScalapackStore(FILE** store, int myrow, int mycol, const char* location);
/**
 * open withreading autorisation the file corresponding to the processor located at myrow,mycol on the CBLAC processor grid
 * each process of the grids should call the functions
 * assume the grid is the same than the one used to create the store
 */
extern void readScalapackStore(FILE** store, int myrow, int mycol, const char* location);

/**
 * save the local matrix lmat with dimensions (nla,mla) corresponding to the process used to open the store
 * assume the store has been created with openScalapackStore
 * each process of the grids should call the functions
 */
extern int saveLocalMatrix(double * lmat,int nla, int mla, FILE* store);

/**
 * read the local matrix in lmat with dimensions (nla,mla)
 * assume the store has been created with readScalapackStore
 * assume that the file exists and has the correct dimension in the location directory
 * each process of the grids should call the functions
 */
extern int readLocalMatrix(double* lmat,int nla, int mla, FILE* store);

/**
 * save the matrix descriptor
 * should be call only once for example by the process of rank 0
 */
extern void saveMatrixDescriptor(int * desc, const char* location);
extern void readMatrixDescriptor(int * desc, const char* location);