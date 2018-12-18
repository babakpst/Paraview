all: example2
clean:; rm -rf example2 *.o

example2: example2.f90
	h5fc -o example2 example2.f90 -I${TACC_HDF5_INC} -L${TACC_HDF5_LIB} -lhdf5 -lz
