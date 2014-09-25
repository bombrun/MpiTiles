
extern void openStore(FILE** store, int i_block, const char* location);
extern void readStore(FILE** store, int i_block, const char* location);
extern void closeStore(FILE* store);
extern int  saveBlock(int i0, int i1, int j0, int j1, double * mat,FILE* store);
extern int readNextBlockDimension(int* dim, FILE* store);
extern void readNextBlock(int i0, int i1, int j0, int j1, double * mat, FILE* store);