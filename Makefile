CC=icc
CFLAGS=-Wall -O2 -xHost -mkl
LDFLAGS=-lmpi -mkl

debug:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi main.c mpiutil.c normals.c matrixBlockStore.c $(LDFLAGS)

stable:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi main.c mpiutil.c normals.c matrixBlockStore.c  $(LDFLAGS)

inde:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi mainInde.c mpiutil.c normals.c matrixBlockStore.c $(LDFLAGS)
	
mkl:clean
	$(CC) $(CFLAGS) -o bigmatrix.mpi  mainMkl.c mpiutil.c normals.c matrixBlockStore.c  $(LDFLAGS)

eigen:
	$(CC) -O2 -o eigen.exe lapackReadStore.c mpiutil.c normals.c matrixBlockStore.c -mkl

test:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi test.c mpiutil.c normals.c matrixBlockStore.c $(LDFLAGS)

clean:
	rm -vfr *~ bigmatrix
