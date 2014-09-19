/**
 * split ng in p segments according to rang,
 * assume ng>p>rang
 * return the left index of the segment
 */
extern int mpi_get_i0(int ng, int rang, int p);
/**
 * split ng in p segments according to rang,
 * assume ng>p>rang
 * return the right index of the segment
 */
extern int mpi_get_i1(int ng, int rang, int p);

/**
 * return the total number of 
 */
extern int mpi_get_total_blocks(int n);



    
    
    // dispatch blocks of matrixCor between the process
    // each process contains its diagonal block
    // b11	b12	b13	b14
    //		b22	b23	b24
    //			b33	b34
    //				b44
    
    // p1 : b11	b12 b13
    // p2 : b22 b23 b24
    // p3 : b33	b34
    // p4 : b44 b14