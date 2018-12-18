all: example
clean:; rm -rf example *.o

example: example.f90
	h5fc -o example example.f90 -I${TACC_HDF5_INC} -L${TACC_HDF5_LIB} -lhdf5 -lz
