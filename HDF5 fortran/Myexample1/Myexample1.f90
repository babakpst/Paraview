! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the files COPYING and Copyright.html.  COPYING can be found at the root   *
!   of the source code distribution tree; Copyright.html can be found at the  *
!   root level of an installed copy of the electronic HDF5 document set and   *
!   is linked from the top-level documents page.  It can also be found at     *
!   http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
!   access to either file, you may request a copy from help@hdfgroup.org.     *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! The following example shows how to write and read to/from an existing dataset.
! It opens the file created in the previous example, obtains the dataset
! identifier, writes the data to the dataset in the file,
! then reads the dataset  to memory.
!
! This example is used in the HDF5 Tutorial.

PROGRAM H5_RDWT

  USE HDF5 ! This module contains all necessary modules

  IMPLICIT NONE

  CHARACTER(LEN=12), PARAMETER :: filename1 = "dsetf_int.h5" ! File name
  CHARACTER(LEN=13), PARAMETER :: filename2 = "dsetf_real.h5" ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname1 = "dset_int"     ! Dataset name
  CHARACTER(LEN=9), PARAMETER :: dsetname2 = "dset_real"     ! Dataset name

  Integer, parameter :: d1=10,d2=10

  INTEGER(HID_T) :: file_id1       ! File identifier
  INTEGER(HID_T) :: file_id2       ! File identifier
  INTEGER(HID_T) :: dset_id1       ! Dataset identifier
  INTEGER(HID_T) :: dset_id2       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id1     ! Dataspace identifier
  INTEGER(HID_T) :: dspace_id2     ! Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(2) :: dims = (/d1,d2/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  INTEGER     ::   rank = 2                        ! Dataset rank
  INTEGER     ::   error ! Error flag
  INTEGER     ::  i, j


  INTEGER, DIMENSION(d1,d2) :: dset_data_int, data_out_int ! Data buffers
  real, DIMENSION(d1,d2) :: dset_data_real, data_out_real ! Data buffers


  !
  ! Initialize the dset_data array.
  !
  DO i = 1, d1
     DO j = 1, d2
        dset_data_int(i,j) = (i-1)*6 + j
     END DO
  END DO

  DO i = 1, d1
     DO j = 1, d2
        dset_data_real(i,j) = 0.1*i+0.2*j ;
     END DO
  END DO

  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(error)

!  !
!  ! Open an existing file.
!  !
!  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

  !
  ! Create a new file using default properties.
  !
  CALL h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file_id1, error)
  CALL h5fcreate_f(filename2, H5F_ACC_TRUNC_F, file_id2, error)

!  !
!  ! Open an existing dataset.
!  !
!  CALL h5dopen_f(file_id, dsetname, dset_id, error)

  !
  ! Create the dataspace.
  !
  CALL h5screate_simple_f(rank, dims, dspace_id1, error)
  CALL h5screate_simple_f(rank, dims, dspace_id2, error)

  !
  ! Create the dataset with default properties.
  !
  CALL h5dcreate_f(file_id1, dsetname1, H5T_NATIVE_INTEGER, dspace_id1, dset_id1, error)
  CALL h5dcreate_f(file_id2, dsetname2, H5T_NATIVE_REAL, dspace_id2, dset_id2, error)

  !
  ! Write the dataset.
  !
  data_dims(1) = d1 ;
  data_dims(2) = d2 ;
  CALL h5dwrite_f(dset_id1, H5T_NATIVE_INTEGER, dset_data_int, data_dims, error)
  CALL h5dwrite_f(dset_id2, H5T_NATIVE_REAL, dset_data_real, data_dims, error)

  !
  ! Read the dataset.
  !
  CALL h5dread_f(dset_id1, H5T_NATIVE_INTEGER, data_out_int, data_dims, error)
  CALL h5dread_f(dset_id2, H5T_NATIVE_REAL, data_out_real, data_dims, error)

print *, 'int'
print *, data_out_int

print *, 'real'
print *, data_out_real 
  !
  ! Close the dataset.
  !
  CALL h5dclose_f(dset_id1, error)
  CALL h5dclose_f(dset_id2, error)

  !
  ! Close the file.
  !
  CALL h5fclose_f(file_id1, error)
  CALL h5fclose_f(file_id2, error)

  !
  ! Close FORTRAN interface.
  !
  CALL h5close_f(error)

END PROGRAM H5_RDWT

