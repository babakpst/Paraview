# =================================================================================================================
#
# Developer    : Babak Poursartip 
#                The University of Texas at Austin, TX, USA
# Code         : FEM ?????????????????????
# makefile     : PETSc Library
# Last Update  : 17 July 2012
#
# IMPORTANT NOTICE:
# This make file works only for ".F90". If there exists a file with other extensions, the "OBJSF" MUST be modified.
#
# =================================================================================================================


info ::
	@echo "do:"
	@echo "  make SOURCE"
	@echo "    where SOURCE is one of ${SOURCESF}"

FFLAGS	    = -O3
CPPFLAGS    = 
FPPFLAGS    =
LOCDIR      = /home/babak/Thesis/Program/MAIN

SOURCESF    =   Myexample2.F90

SOURCEH     =

OBJSF       =  $(SOURCESF:%.F90=%.o)

MANSEC      = KSP
#CLEANFILES = rhs.vtk solution.vtk


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

main :  $(OBJSF) chkopts
	@echo Linking...
	-${FLINKER} $(FFLAGS) $(OBJSF) -o main ${PETSC_KSP_LIB}
	@echo Done!

include ${PETSC_DIR}/conf/test



# clean rule for object files
cln ::
	@echo Clean Object files...
	rm -rf *.o main












