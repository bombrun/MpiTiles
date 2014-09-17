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
