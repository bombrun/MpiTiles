
/**
 * Tasks managment:
 * split ng in p segments according to rang,
 * assume ng>p>rang
 * return the left index of the segment
 */
extern int mpi_get_i0(int ng, int rang, int p);

/**
 * Tasks managment:
 * split ng in p segments according to rang,
 * assume ng>p>rang
 * return the right index of the segment
 */
extern int mpi_get_i1(int ng, int rang, int p);

/**
 * Tasks managment:
 * return the total number of blocks tasks n(n+1)/2
 */
extern int mpi_get_total_blocks(int n);

/**
 * Tasks managment:
 * return the block id j_block associated to the block i_block and the task_id assuming n diagonal blocks
 * based on a rotating permutation
 */
extern int mpi_get_diag_block_id(int i_diag_block,int diag_taks_id, int n_diag_blocks);


/**
 * Tasks managment:
 * return the total number of diagonal block
 * if not specify return p
 */
int get_n_blocks(int argc, char **argv, int p);

/**
 *  Tasks managment:
 *  diagonal block distribution over all the process
 */
int get_n_pTasks(int p, int rank, int n_blocks);
/**
 *  Tasks managment:
 *  tasks distribution over all the diagonal blocks
 */
int get_n_blockTasks(int i_block, int n_blocks);
    
    
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