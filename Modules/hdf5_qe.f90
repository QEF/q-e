! Copyright (C) 2003-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__HDF5)
module hdf5_qe
  !
  USE HDF5
  !USE, intrinsic :: ISO_C_binding
  USE Kinds, ONLY : DP
  !
  implicit none
  !
  ! This module contains some common subroutines used to read and write
  ! in HDF5 format the data produced by Quantum ESPRESSO package.
  !
  ! written by: Nicola Varini 2016
  !  

  TYPE HDF5_type
   INTEGER(HID_T) :: file_id       ! File identifier 
   INTEGER(HID_T) :: dset_id       ! Dataset identifier 
   INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
   INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
   INTEGER(HID_T) :: plist_id      ! Property list identifier 
   INTEGER(HID_T) :: group_id      ! Group identifier 
   CHARACTER(LEN=40) :: dsetname  ! Dataset name
   INTEGER          :: rank
   INTEGER(HSIZE_T), DIMENSION(2) :: counts, counts_g, offset
   INTEGER(HSIZE_T), DIMENSION(1:2) :: size
   INTEGER(HID_T) :: crp_list      ! Dataset creation property identifier 
   INTEGER          :: comm
   INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
   INTEGER(HSIZE_T), DIMENSION(1:2) :: chunk_dim
   character(len=256) filename
  END TYPE HDF5_type

  TYPE(HDF5_type), save :: evc_hdf5, evc_hdf5_write, evq_hdf5_write
  TYPE(HDF5_type), save :: rho_hdf5_write, eig_hdf5_write
  TYPE(HDF5_type), save :: g_hdf5_write, gk_hdf5_write
  
  INTEGER, save ::  off_npw, npw_g, idone_debug


  INTERFACE add_attributes_hdf5
    MODULE PROCEDURE add_attributes_hdf5_i, add_attributes_hdf5_r, &
                     add_attributes_hdf5_c
  END INTERFACE

  INTERFACE read_attributes_hdf5
    MODULE PROCEDURE read_attributes_hdf5_i, read_attributes_hdf5_r
  END INTERFACE


  contains

  subroutine initialize_hdf5()
    implicit none
    integer :: error
    call h5open_f(error)
  end subroutine initialize_hdf5

  subroutine finalize_hdf5(hdf5desc)
    implicit none
    type(HDF5_type), intent(in) :: hdf5desc
    integer :: error
    
    call h5pclose_f(hdf5desc%plist_id,error)
    call h5close_f(error)
  end subroutine finalize_hdf5

  subroutine h5_write_gvecs(hdf5desc, filename, nr1, nr2, nr3, ngm, gamma_only, mill) 
   implicit none 
   type (hdf5_type), intent(inout) :: hdf5desc
   character(len=*), intent (in)   :: filename
   integer, intent(in)             :: nr1, nr2, nr3, ngm
   logical, intent(in)             :: gamma_only
   integer,intent(in)              :: mill (3,ngm)
   continue 
  end subroutine h5_write_gvecs
 
  subroutine setup_file_property_hdf5(hdf5desc ,filename, hyperslab, write, kpoint)
   use parallel_include
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc 
   character(len=*), intent(inout) :: filename
   logical,  intent(in) :: hyperslab, write
   integer, intent(in) ::  kpoint
   integer(HID_T) :: plist_id
   integer :: error, info
   character*12 kstring
   write(kstring,'(I0)') kpoint
   kstring='KPOINT'//kstring


   info = MPI_INFO_NULL



   if(hyperslab .eqv. .true. ) then

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, hdf5desc%plist_id, error) ! Properties for file creation
        CALL h5pset_fapl_mpio_f(hdf5desc%plist_id, hdf5desc%comm, info, error) ! Stores MPI IO communicator information to the file access property list
        if(kpoint.eq.1)then
          CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error, access_prp = hdf5desc%plist_id) ! create the file collectively
        else
          CALL h5fopen_f(filename, H5F_ACC_RDWR_F, hdf5desc%file_id, error, access_prp = hdf5desc%plist_id) ! create the file collectively
        endif
        CALL h5pclose_f(hdf5desc%plist_id, error)
   else

    if(write .eqv. .true.)then
      if(kpoint.eq.1)then
        !CALL h5pcreate_f(H5P_FILE_ACCESS_F, hdf5desc%plist_id, error)
        !CALL h5pset_fapl_mpio_f(hdf5desc%plist_id, hdf5desc%comm, info, error)
        !CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error, access_prp=hdf5desc%plist_id) ! create the file collectively
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error) 
      endif
    else
        CALL h5fopen_f(filename, H5F_ACC_RDWR_F, hdf5desc%file_id, error) 
        IF (error /= 0) CALL errore ('setup_file_property_hdf5','error opening '//filename,error)
    endif
   endif
   
  end subroutine setup_file_property_hdf5


  subroutine define_dataset_hdf5_hyperslab(hdf5desc, kpoint)
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc
   integer,intent(in)             :: kpoint
   integer                        :: error
   character*12                   :: kstring  
   write(kstring,'(I0)') kpoint
   kstring=trim('KPOINT')//kstring
   hdf5desc%dsetname = 'evc'

   CALL h5gcreate_f(hdf5desc%file_id, kstring, hdf5desc%group_id, error)
   CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts_g, hdf5desc%filespace, error) !define HDF5 dataset
   CALL h5dcreate_f(hdf5desc%group_id, hdf5desc%dsetname, H5T_NATIVE_DOUBLE, hdf5desc%filespace, &
                    hdf5desc%dset_id, error)
   CALL h5sclose_f(hdf5desc%filespace, error)

   CALL h5dclose_f(hdf5desc%dset_id, error)
   CALL h5gclose_f(hdf5desc%group_id, error)

   end subroutine define_dataset_hdf5_hyperslab


 
  subroutine  write_data_hdf5(hdf5desc, data,  kpoint)
   USE kinds, ONLY : DP
   USE ISO_C_BINDING

   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc
   complex(kind=dp), target, intent(inout) :: data(:,:)
   integer, intent(in) :: kpoint
   integer :: error, datadim1, datadim2 
   real(kind=dp)   :: tmp
   integer(HID_T)     :: complex_id, double_id
   integer(HSIZE_T)   :: double_size, complex_size
   TYPE(C_PTR)        :: f_ptr
   character*12       :: kstring  
   write(kstring,'(I0)') kpoint
   kstring=trim('KPOINT')//kstring

   
   CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
   CALL h5dopen_f(hdf5desc%group_id, hdf5desc%dsetname, hdf5desc%dset_id, error)
   CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error) 
   CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)

   CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, hdf5desc%offset, hdf5desc%counts, error) ! create hyperslab to read from more than 1 proc

   CALL h5pcreate_f(H5P_DATASET_XFER_F, hdf5desc%plist_id, error)
   CALL h5pset_dxpl_mpio_f(hdf5desc%plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  
   f_ptr = C_LOC(data(1,1))
   CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error,&
                  file_space_id = hdf5desc%filespace, mem_space_id = hdf5desc%memspace, &
                  xfer_prp = hdf5desc%plist_id)

   CALL h5dclose_f(hdf5desc%dset_id, error)
   CALL h5gclose_f(hdf5desc%group_id, error)
  end subroutine write_data_hdf5






  subroutine  read_data_hdf5(hdf5desc, data, kpoint)
   type(HDF5_type), intent(inout) :: hdf5desc
   complex(kind=dp),target, intent(inout) :: data(:,:)
   integer,intent(in) :: kpoint
   integer :: error
   TYPE(C_PTR) :: f_ptr
   character*12       :: kstring  
   write(kstring,'(I0)') kpoint
   kstring=trim('KPOINT')//kstring

   CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
   CALL h5dopen_f(hdf5desc%group_id, hdf5desc%dsetname, hdf5desc%dset_id, error)
   CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)

   CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, hdf5desc%offset, hdf5desc%counts, error)

   CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error) 

   f_ptr = C_LOC(data(1,1))
   CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                  mem_space_id = hdf5desc%memspace, file_space_id = hdf5desc%filespace ,&
                 xfer_prp = hdf5desc%plist_id)

   CALL h5dclose_f(hdf5desc%dset_id, error)
   CALL h5gclose_f(hdf5desc%group_id, error)

  end subroutine read_data_hdf5

  SUBROUTINE prepare_index_hdf5(sendm,recm,globalm,comm,nproc)

   USE parallel_include
   USE mp, ONLY : mp_sum
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: comm, nproc
   INTEGER, INTENT(INOUT) :: sendm, recm, globalm
   INTEGER :: errore

   call mpi_scan(sendm,recm,1,MPI_INTEGER,MPI_SUM,comm,errore)
   recm=recm-sendm
   globalm=sendm
   call mp_sum(globalm,comm)


 END SUBROUTINE prepare_index_hdf5




  subroutine prepare_for_writing_final(hdf5desc,comm,filename_input,kpoint,add_group)
    USE io_files, ONLY : wfc_dir, prefix, tmp_dir
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    character(len=*), intent(in):: filename_input
    integer, intent(in) :: comm
    integer,  intent(in), optional :: kpoint
    logical, intent(in), optional  :: add_group
    character(len=256) filename
    integer :: ik, error
    logical ::   add_group_internal = .true.
    character*12 kstring
    
 
    hdf5desc%comm=comm
    hdf5desc%filename=filename_input
    if ( present (add_group) ) add_group_internal = add_group
    if(present(kpoint)) then
      write(kstring,'(I0)') kpoint
      kstring=trim('KPOINT')//kstring
      IF ( add_group_internal) THEN 
         CALL setup_file_property_hdf5(hdf5desc,hdf5desc%filename ,.false.,.true.,kpoint)
      ELSE 
         CALL setup_file_property_hdf5(hdf5desc,hdf5desc%filename ,.false.,.true.,1)
      END IF 
    
      CALL h5fopen_f(hdf5desc%filename, H5F_ACC_RDWR_F, hdf5desc%file_id, error) ! create the file collectively
      
      CALL h5gcreate_f(hdf5desc%file_id, kstring, hdf5desc%group_id, error)
      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      CALL setup_file_property_hdf5(hdf5desc,hdf5desc%filename ,.false.,.true.,1)
    endif

  end subroutine prepare_for_writing_final



  subroutine prepare_for_reading_final(hdf5desc,comm,filename_input,kpoint)
    USE io_files,             ONLY : wfc_dir, prefix, tmp_dir
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    character(len=*), intent(in):: filename_input
    integer, intent(in) :: comm
    integer, intent(in), optional :: kpoint
    character(len=256) filename
    integer :: ik
 
    hdf5desc%comm=comm
    hdf5desc%rank =1 
    !filename = trim(filename_input) //".wfchdf5"
    filename=filename_input
    if(present(kpoint))then
      CALL setup_file_property_hdf5(hdf5desc,filename ,.false.,.false.,kpoint)
    else
      CALL setup_file_property_hdf5(hdf5desc,filename ,.false.,.false.,1)
    end if

  end subroutine prepare_for_reading_final

  subroutine read_rho(hdf5desc,dsetname,var)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: dsetname
    real(kind=DP), target, intent(in) :: var(:)
    INTEGER(HID_T) :: dspace_id, dset_id, dtype_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    character*12 dset_name
    TYPE(C_PTR) :: f_ptr
    write(dset_name,'(I0)') dsetname
    dset_name='K'//dset_name
    counts=size(var)
    CALL h5dopen_f(hdf5desc%file_id, dset_name, dset_id, error)
    CALL h5dget_type_f(dset_id, dtype_id, error)
    f_ptr = C_LOC(var(1))
    CALL h5dread_f(dset_id, dtype_id, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
  end subroutine read_rho

  subroutine write_rho(hdf5desc,dsetname,var)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: dsetname
    real(kind=DP), target, intent(in) :: var(:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    character*12 dset_name
    TYPE(C_PTR) :: f_ptr
    write(dset_name,'(I0)') dsetname
    dset_name='K'//dset_name
    counts=size(var)
    CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
    CALL h5dcreate_f(hdf5desc%file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
    f_ptr = C_LOC(var(1))
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
  end subroutine write_rho

  subroutine write_eig(hdf5desc,var,kpoint)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: kpoint
    real(kind=DP), target, intent(in) :: var(:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    TYPE(C_PTR) :: f_ptr
    character*12 kstring
    write(kstring,'(I0)') kpoint
    kstring='KPOINT'//kstring

    counts=size(var)
    CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
    CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
    CALL h5dcreate_f(hdf5desc%group_id, kstring, H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
    f_ptr = C_LOC(var(1))
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5gclose_f(hdf5desc%group_id, error)
  end subroutine write_eig

  subroutine read_eig(hdf5desc,var,kpoint)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) ::  kpoint
    real(kind=DP), target, intent(inout) :: var(:)
    INTEGER(HID_T) :: dtype_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    TYPE(C_PTR) :: f_ptr
    character*12 kstring
    character*100 errmsg

    write(kstring,'(I0)') kpoint
    kstring='KPOINT'//kstring

    counts=size(var)
    CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
    CALL h5dopen_f(hdf5desc%group_id, kstring, dset_id, error)
    CALL h5dget_type_f(dset_id, dtype_id, error)
    f_ptr = C_LOC(var(1))
    !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    CALL h5dread_f(dset_id, dtype_id, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(hdf5desc%group_id, error)
  end subroutine read_eig



  subroutine write_evc(hdf5desc,dsetname,var,kpoint)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: dsetname
    integer, intent(in), optional ::  kpoint
    complex(kind=DP), target, intent(in) :: var(:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    character*12 dset_name
    TYPE(C_PTR) :: f_ptr
    character*12 kstring
    write(kstring,'(I0)') kpoint
    kstring='KPOINT'//kstring
    write(dset_name,'(I0)') dsetname
    dset_name='BAND'//dset_name
    counts=size(var)*2  
    CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
    if(present(kpoint))CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
    if(present(kpoint)) then
      CALL h5dcreate_f(hdf5desc%group_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
    else
      CALL h5dcreate_f(hdf5desc%file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
    endif
    f_ptr = C_LOC(var(1))
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    if(present(kpoint))CALL h5gclose_f(hdf5desc%group_id, error)
  end subroutine write_evc

  subroutine read_evc(hdf5desc,dsetname,var,kpoint)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: dsetname, kpoint
    complex(kind=DP), target ,intent(inout) :: var(:)
    INTEGER(HID_T) :: dtype_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    character*12 dset_name
    TYPE(C_PTR) :: f_ptr
    character*12 kstring
    character*100 errmsg

    write(dset_name,'(I0)') dsetname
    write(kstring,'(I0)') kpoint
    kstring='KPOINT'//kstring
    dset_name='BAND'//dset_name

    counts=size(var)*2  
    CALL h5gopen_f(hdf5desc%file_id, kstring, hdf5desc%group_id, error)
    !if(error.ne.0) call errore('error in h5gopen_f',' ',error)
    CALL h5dopen_f(hdf5desc%group_id, dset_name, dset_id, error)
    !if(error.ne.0) call errore('error in h5dopen_f',' ',error)
    !CALL h5dget_type_f(dset_id, dtype_id, error)
    f_ptr = C_LOC(var(1))
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    if(error.ne.0) call errore('error in h5dread_f',' ',error)
    !CALL h5dread_f(dset_id, dtype_id, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(hdf5desc%group_id, error)
     
  end subroutine read_evc

  subroutine write_g(hdf5desc,var,kpoint)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in), optional :: kpoint
    integer, target, intent(in) :: var(:,:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    TYPE(C_PTR) :: f_ptr
    character*12 kstring

    counts=size(var,1)*size(var,2)
    if(present(kpoint))then
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
      CALL h5dcreate_f(hdf5desc%group_id, kstring, H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(var(1,1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)
      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%file_id, 'Miller indexes', H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(var(1,1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)
    endif
  end subroutine write_g



  subroutine write_gkhdf5(hdf5desc,xk,igwk,mill_g,kpoint)
    USE kinds, ONLY : DP
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in), optional :: kpoint
    real(kind=DP), target, intent(in) :: xk(:)
    integer, target, intent(in) :: igwk(:), mill_g(:,:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    TYPE(C_PTR) :: f_ptr
    character*12 kstring

    if(present(kpoint))then
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)

      write(kstring,'(I0)') kpoint
      kstring='xk'//kstring
      counts=size(xk)
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%group_id, kstring, H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(xk(1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)

      write(kstring,'(I0)') kpoint
      kstring='igwk'//kstring
      counts=size(igwk)
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%group_id, kstring, H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(igwk(1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)

      write(kstring,'(I0)') kpoint
      kstring='mill_g'//kstring
      counts=size(mill_g)
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%group_id, kstring, H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(mill_g(1,1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)


      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      counts=size(xk)
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%file_id, 'xk', H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(xk(1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)

      counts=size(igwk)
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%file_id, 'igwk', H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(igwk(1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)

      counts=size(mill_g)
      CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
      CALL h5dcreate_f(hdf5desc%group_id, 'mill_g', H5T_NATIVE_INTEGER, dspace_id, &
                      dset_id, error)
      f_ptr = C_LOC(mill_g(1,1))
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)

    endif
  end subroutine write_gkhdf5


  subroutine initialize_io_hdf5(hdf5desc,comm, data, write,kpoint)
    USE io_files, ONLY : wfc_dir, prefix, tmp_dir
    USE kinds,    ONLY : dp
    USE mp_world, ONLY : nproc
    implicit none
    
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    complex(kind=dp), intent(in) :: data(:,:)
    integer, intent(in) :: comm, kpoint
    logical, intent(in) :: write
    character(len=80) :: filename
    integer :: npwx, nbnd
    npwx=size(data(:,1))
    nbnd=size(data(1,:))
    
    filename=trim(tmp_dir) //TRIM(prefix) //".wfchdf5"
    call initialize_hdf5_array(hdf5desc,comm,npwx,nbnd)
    if(write .eqv. .true.)then
      CALL setup_file_property_hdf5(hdf5desc, filename,.true.,.true.,kpoint)
    else
      CALL setup_file_property_hdf5(hdf5desc, filename,.false.,.false.,kpoint)
    endif
    CALL prepare_index_hdf5(npwx,off_npw,npw_g,hdf5desc%comm,nproc)
    CALL set_index_hdf5(hdf5desc,data,off_npw,npw_g,2)

  end subroutine initialize_io_hdf5

  subroutine initialize_hdf5_array(hdf5desc,comm,n1,n2)
    implicit none
    integer, intent(in) :: n1, n2, comm
    type(HDF5_type), intent(inout) ::  hdf5desc
    hdf5desc%dsetname="evc"
    hdf5desc%comm=comm
    hdf5desc%rank =2 
    hdf5desc%chunk_dim=(/n1,n2/)
    hdf5desc%size(1) = n1*2
    hdf5desc%size(2) = n2
    hdf5desc%offset(1) = 0
    hdf5desc%offset(2) = 0
   
  end subroutine initialize_hdf5_array

  SUBROUTINE set_index_hdf5(hdf5desc, var, offset, nglobal,tsize)
    
    USE kinds, only : DP
    implicit none
    COMPLEX(DP), intent(in) :: var(:,:) 
    type(HDF5_type), intent(inout) :: hdf5desc
    INTEGER, intent(in)  :: offset, nglobal,tsize
    
    hdf5desc%counts(1)   = size(var(:,1))*tsize
    hdf5desc%counts(2)   = size(var(1,:)) 
    hdf5desc%counts_g(1) = nglobal*tsize
    hdf5desc%counts_g(2) = size(var(1,:)) 
    hdf5desc%offset(1)   = offset*tsize
    hdf5desc%offset(2)   = 0
 
  END SUBROUTINE set_index_hdf5


  SUBROUTINE add_attributes_hdf5_i(hdf5desc, attr_data, attr_name, kpoint)
    implicit none
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) ::  attr_data
    integer, intent(in),optional :: kpoint 
    CHARACTER(LEN=*), intent(in) :: attr_name
    character*12 kstring
    integer :: error
    INTEGER     ::   arank = 1                      ! Attribure rank
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
    INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
    INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
    INTEGER(HID_T) :: attr_id      ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  
    data_dims(1) = 1 
    
    if(present(kpoint)) then
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5acreate_f(hdf5desc%group_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(aspace_id, error)
      !
      ! End access to the dataset and release resources used by it.
      !
      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5acreate_f(hdf5desc%file_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(aspace_id, error)
 
    endif
    
  END SUBROUTINE add_attributes_hdf5_i

  SUBROUTINE add_attributes_hdf5_r(hdf5desc, attr_data, attr_name, kpoint)
    implicit none
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in), optional :: kpoint
    real(DP), intent(in) :: attr_data 
    CHARACTER(LEN=*), intent(in) :: attr_name
    character*12 kstring
    integer :: error
    INTEGER     ::   arank = 1                      ! Attribure rank
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
    INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
    INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
    INTEGER(HID_T) :: attr_id      ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  
    data_dims(1) = 1 

    if(present(kpoint)) then 
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
    
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5acreate_f(hdf5desc%group_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(aspace_id, error)
      !
      ! End access to the dataset and release resources used by it.
      !
      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5acreate_f(hdf5desc%file_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(aspace_id, error)
      !
      ! End access to the dataset and release resources used by it.
      !

    endif
    
  END SUBROUTINE add_attributes_hdf5_r

  SUBROUTINE add_attributes_hdf5_c(hdf5desc, attr_data, attr_name, kpoint)
    implicit none
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in), optional :: kpoint
    !LOGICAL, intent(in) :: attr_data 
    CHARACTER(LEN=*), intent(in) :: attr_data
    CHARACTER(LEN=*), intent(in) :: attr_name
    character*100 kstring
    integer :: error
    INTEGER     ::   arank = 1                      ! Attribure rank
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
    INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
    INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
    INTEGER(HID_T) :: attr_id      ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  
    !data_dims(1) = 1 
    data_dims(1) = len(attr_data)
    
    if(present(kpoint)) then
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
      !write(attrdata,'(I0)') attr_data
    
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5acreate_f(hdf5desc%group_id, attr_name, H5T_NATIVE_CHARACTER, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_CHARACTER, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(aspace_id, error)
      !
      ! End access to the dataset and release resources used by it.
      !
      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      !write(attrdata,'(I0)') attr_data
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5acreate_f(hdf5desc%file_id, attr_name, H5T_NATIVE_CHARACTER, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_CHARACTER, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(aspace_id, error)
    endif
    
  END SUBROUTINE add_attributes_hdf5_c


  SUBROUTINE read_attributes_hdf5_i(hdf5desc, attr_data, attr_name, kpoint, debug)
    implicit none
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in), optional  :: kpoint, debug
    integer, intent(out) ::  attr_data 
    CHARACTER(LEN=*), intent(in) :: attr_name
    character*12 kstring
    integer :: error
    INTEGER     ::   arank = 1                      ! Attribure rank
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
    INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
    INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
    INTEGER(HID_T) :: attr_id      ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  
    data_dims(1) = 1 
    attrlen = 1 
    if(present(kpoint))then
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
      CALL h5aopen_name_f(hdf5desc%group_id,attr_name,attr_id,error)
      CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5gclose_f(hdf5desc%group_id, error)
   else
      CALL h5aopen_name_f(hdf5desc%file_id,attr_name,attr_id,error)
      CALL h5aget_type_f(attr_id, atype_id, error)
      CALL h5aread_f(attr_id, atype_id, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
   endif
    
  END SUBROUTINE read_attributes_hdf5_i

  SUBROUTINE read_attributes_hdf5_r(hdf5desc, attr_data, attr_name, kpoint)
    implicit none
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in), optional  :: kpoint
    real(DP), intent(out) ::  attr_data 
    CHARACTER(LEN=*), intent(in) :: attr_name
    character*12 kstring
    integer :: error
    INTEGER     ::   arank = 1                      ! Attribure rank
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
    INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
    INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
    INTEGER(HID_T) :: attr_id      ! Attribute Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  
    data_dims(1) = 1 
    attrlen = 1 
    if(present(kpoint))then
      write(kstring,'(I0)') kpoint
      kstring='KPOINT'//kstring
      CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
      CALL h5aopen_name_f(hdf5desc%group_id,attr_name,attr_id,error)
      CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5gclose_f(hdf5desc%group_id, error)
    else
      CALL h5aopen_name_f(hdf5desc%file_id,attr_name,attr_id,error)
      CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)
      CALL h5aclose_f(attr_id, error)
    endif
    
  END SUBROUTINE read_attributes_hdf5_r

  SUBROUTINE hdf5_close(hdf5desc)
    implicit none
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer :: errore 
    CALL h5fclose_f(hdf5desc%file_id,errore)
  END SUBROUTINE hdf5_close

  SUBROUTINE write_attributes(hdf5desc, ngw, gamma_only, igwx, &
    nbnd, ik, nk, ispin, nspin, scalef)
    implicit none
    INTEGER,            INTENT(IN) :: ik, nk,  ispin, nspin
    REAL(DP),           INTENT(IN) :: scalef    
    LOGICAL,            INTENT(IN) :: gamma_only
    INTEGER,            INTENT(IN) :: nbnd, ngw, igwx
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    integer                        :: gammaonly = 0
    CALL add_attributes_hdf5(hdf5desc,ngw,"ngw",ik)
    IF ( gamma_only) gammaonly = 1 
    CALL add_attributes_hdf5(evc_hdf5_write,gammaonly,"gamma_only",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,igwx,"igwx",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,nbnd,"nbnd",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,ik,"ik",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,nk,"nk",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,ispin,"ispin",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,nspin,"nspin",ik)
    CALL add_attributes_hdf5(evc_hdf5_write,scalef,"scale_factor",ik)


  END SUBROUTINE write_attributes

end module hdf5_qe
#else 
module hdf5_qe
implicit none
integer  :: pippo
end module hdf5_qe
#endif 

