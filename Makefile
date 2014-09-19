CC=mpicc
CFLAGS="-Wall"

debug:clean
	$(CC) $(CFLAGS) -g -o bigmatrix.mpi main.c mpiutil.c normals.c
stable:clean
	$(CC) $(CFLAGS)    -o bigmatrix.mpi main.c mpiutil.c normals.c
test:clean
	$(CC) $(CFLAGS) -g -o bigmatrixtest.mpi test.c mpiutil.c normals.c
clean:
	rm -vfr *~ bigmatrix
