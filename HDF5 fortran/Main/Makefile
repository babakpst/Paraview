all: example
clean:; rm -rf example *.o

Main: Main.f90
	h5fc -o main Main.f90 Geometry.f90 -I${TACC_HDF5_INC} -L${TACC_HDF5_LIB} -lhdf5 -lz
