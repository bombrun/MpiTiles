CC=mpicc
CFLAGS="-Wall"

MLKCC=icc
MKLCFLAGS="-Wall -g -O2 -xHost -mkl"
MKLLDFLAGS="-lmpi -mkl"
  
debug:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi main.c mpiutil.c normals.c
	
stable:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi main.c mpiutil.c normals.c
	
inde:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi mainInde.c mpiutil.c normals.c matrixBlockStore.c
	
mkl:clean
	$(MKLCC) $(MKLCFLAGS) -o bigmatrix.mpi  mainMkl.c mpiutil.c normals.c $(MKLLDFLAGS)
	
test:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi test.c mpiutil.c normals.c
	
clean:
	rm -vfr *~ bigmatrix
