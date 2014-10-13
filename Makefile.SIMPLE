CC=mpicc
CFLAGS="-Wall"

ICC=icc

MKLCFLAGS="-Wall -g -O2 -xHost -mkl"
MKLLDFLAGS="-lmpi -mkl"



  
debug:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi main.c mpiutil.c normals.c matrixBlockStore.c
	
stable:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi main.c mpiutil.c normals.c matrixBlockStore.c
	
inde:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi mainInde.c mpiutil.c normals.c matrixBlockStore.c
	
mkl:clean
	$(ICC) $(MKLCFLAGS) -o bigmatrix.mpi  mainMkl.c mpiutil.c normals.c matrixBlockStore.c $(MKLLDFLAGS)
	
eigen:
	$(ICC) -O2 -o eigen.exe lapackReadStore.c mpiutil.c normals.c matrixBlockStore.c -mkl
	
test:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi test.c mpiutil.c normals.c matrixBlockStore.c
	
clean:
	rm -vfr *~ 
	
cleanAll:
	rm -vrf *~ *.mpi *.exe *.txt
