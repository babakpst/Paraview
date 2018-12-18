all: Myexample1
clean:; rm -rf Myexample1 *.o

Myexample1: Myexample1.f90
	h5fc -o Myexample1 Myexample1.f90 -I${TACC_HDF5_INC} -L${TACC_HDF5_LIB} -lhdf5 -lz
