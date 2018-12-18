all: example
clean:; rm -rf example *.o

example: example.C
	mpicxx -o example example.C -I${TACC_HDF5_INC} -L${TACC_HDF5_LIB} -lhdf5 -lz
