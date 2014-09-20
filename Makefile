CC=mpicc
CFLAGS="-Wall"

MLK_CC = icc
MKL_CFLAGS=-Wall -g -O2 -xHost -mkl
MKL_LDFLAGS= -lmpi -mkl
  
debug:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi main.c mpiutil.c normals.c
stable:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi main.c mpiutil.c normals.c
inde:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi mainInde.c mpiutil.c normals.c	
mkl:clean
	$(MKL_CC) $(MKL_CFLAGS) -o bigmatrix.mpi  mainMkl.c mpiutil.c normals.c $(MKL_LDFLAGS)
test:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi test.c mpiutil.c normals.c
clean:
	rm -vfr *~ bigmatrix
