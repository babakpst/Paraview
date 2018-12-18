
! This is a toy program for HDF5 in Paraview

     PROGRAM Paraview

     USE HDF5 ! This module contains all necessary modules

     IMPLICIT NONE

     CHARACTER(LEN=8), PARAMETER :: filename = "dsetf.h5" ! File name
     CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name

     INTEGER(HID_T) :: file_id       ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


     INTEGER(HSIZE_T), DIMENSION(2) :: dims = (/4,6/) ! Dataset dimensions
     INTEGER     ::   rank = 2                        ! Dataset rank

     INTEGER     ::   error ! Error flag

     !
     ! Initialize FORTRAN interface.
     !
     CALL h5open_f(error)

     !
     ! Create a new file using default properties.
     !
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

     !
     ! Create the dataspace.
     !
     CALL h5screate_simple_f(rank, dims, dspace_id, error)

     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)

     !
     ! End access to the dataset and release resources used by it.
     !
     CALL h5dclose_f(dset_id, error)

     !
     ! Terminate access to the data space.
     !
     CALL h5sclose_f(dspace_id, error)

     !
     ! Close the file.
     !
     CALL h5fclose_f(file_id, error)

     !
     ! Close FORTRAN interface.
     !
     CALL h5close_f(error)

     END PROGRAM Paraview
