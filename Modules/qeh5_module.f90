!
! Copyright (C) 2016-2017 Quantum ESPRESSO Foundation 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! 
!------------------------------------------------------
MODULE qeh5_base_module
  !---------------------------------------------------
  !
  ! this module contains the basic interface for basic operation for 
  ! serial I/O in HDF5 format. The parallel interface remains in file
  ! hdf5_qe.f90 file. 
  ! 
  ! author N. Varini, P. Delugas
  ! last revision June 2017
  ! 
  ! 
#if defined(__HDF5)
  USE KINDS,   ONLY: DP
  USE hdf5 
  USE ISO_C_BINDING
  IMPLICIT NONE 


  TYPE qeh5_hid
    INTEGER(HID_T)        ::   id
  END TYPE qeh5_hid

  TYPE qeh5_file ! this one is also good for groups 
    INTEGER(HID_T)        ::   id 
    CHARACTER(LEN=256)    ::   name   
  END TYPE qeh5_file
  !
  TYPE qeh5_datatype   
     INTEGER(HID_T)        ::  id 
     INTEGER               ::  rank 
     INTEGER,ALLOCATABLE   ::  dims(:)
  END TYPE qeh5_datatype 
  !
  TYPE qeh5_dataspace    
     INTEGER(HID_T)       ::  id
     INTEGER              ::  rank 
     INTEGER(HSIZE_T),ALLOCATABLE  ::  dims(:),maxdims(:) 
     INTEGER(HSIZE_T),ALLOCATABLE  ::  offset(:), count(:), stride(:), block(:)
  END TYPE qeh5_dataspace  
  !
  TYPE qeh5_dataset
     INTEGER(HID_T)        ::    id 
     CHARACTER(LEN=256)    ::    name 
     TYPE(qeh5_datatype)   ::    datatype
     TYPE(qeh5_dataspace)  ::    filespace
     LOGICAL               ::    memspace_ispresent = .FALSE. 
     TYPE(qeh5_dataspace)  ::    memspace  
  END TYPE qeh5_dataset 


  INTERFACE qeh5_set_space
     MODULE PROCEDURE qeh5_wplan_real, qeh5_wplan_complex, qeh5_wplan_integer   
  END INTERFACE


  INTERFACE qeh5_write_dataset
     MODULE PROCEDURE qeh5_write_real, qeh5_write_real_2, qeh5_write_real_3, &
             qeh5_write_complex, qeh5_write_complex_2,qeh5_write_complex_3, &
             qeh5_write_integer,qeh5_write_integer_2,qeh5_write_integer_3 
  END INTERFACE 


  INTERFACE qeh5_read_dataset
     MODULE PROCEDURE qeh5_read_real, qeh5_read_complex, qeh5_read_integer
  END INTERFACE 


  INTERFACE qeh5_add_attribute
     MODULE PROCEDURE    add_attribute_i, add_array_attribute_i, add_attribute_r, add_array_attribute_r, &
                        add_attribute_string
  END INTERFACE

  INTERFACE qeh5_read_attribute 
     MODULE PROCEDURE  read_real_attribute, read_real_array_attribute, read_string_attribute,&
                     read_integer_attribute, read_integer_array_attribute
  END INTERFACE

  INTERFACE qeh5_to_h5id
     MODULE PROCEDURE get_dataset_hid, get_file_hid
  END INTERFACE

  INTERFACE qeh5_close
     MODULE PROCEDURE qeh5_closefile, close_dataset
  END INTERFACE
  ! 
  INTEGER(HID_T)   :: H5_REALDP_TYPE  
  PRIVATE 
  PUBLIC           :: qeh5_file, qeh5_dataset  
  PUBLIC           :: initialize_hdf5, finalize_hdf5, qeh5_openfile, qeh5_close, qeh5_open_group, &
                    qeh5_open_dataset, qeh5_write_dataset, qeh5_read_dataset, qeh5_set_space ,  & 
                    qeh5_add_attribute, qeh5_read_attribute, qeh5_set_file_hyperslab,           &
                    qeh5_set_memory_hyperslab, qeh5_to_h5id                   
  CONTAINS 
  !
  !----------------------------------------------------------
  SUBROUTINE initialize_hdf5()
    !--------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: ierr 
    CALL h5open_f(ierr)
    H5_REALDP_TYPE = h5kind_to_type( DP, H5_REAL_KIND) 
END SUBROUTINE initialize_hdf5

  !-------------------------------
  SUBROUTINE finalize_hdf5()
    !------------------------------
    IMPLICIT NONE
    INTEGER :: ierr 
    CALL h5close_f(ierr )
END SUBROUTINE finalize_hdf5

  !---------------------------------------------------
  FUNCTION get_dataset_hid( h5_dataset ) RESULT (hid)
    !-------------------------------------------------
    IMPLICIT NONE
    TYPE( qeh5_dataset )               :: h5_dataset
    INTEGER (HID_T)                    :: hid 
    hid  = h5_dataset%id
  END FUNCTION get_dataset_hid

  !---------------------------------------------------
  FUNCTION get_file_hid( h5_file ) RESULT (hid)
    !-------------------------------------------------
    IMPLICIT NONE
    TYPE( qeh5_file )                  :: h5_file
    INTEGER (HID_T)                    :: hid 
    hid  = h5_file%id
  END FUNCTION get_file_hid

  !------------------------------------------------------- 
  SUBROUTINE  qeh5_openfile(h5file, file, action , error)
    !-----------------------------------------------------
    IMPLICIT NONE 
    CHARACTER(LEN=*),INTENT(IN)               :: file, action 
    INTEGER,OPTIONAL,INTENT(OUT)              :: error  
    TYPE(qeh5_file), INTENT(OUT)  :: h5file
    ! 
    INTEGER                       :: ierr, jerr 
    ! 
    
    h5file%name=TRIM(file) 
    IF (PRESENT(error))  THEN
       CALL H5Eset_auto_f( 0, ierr )
    END IF 
    SELECT CASE(TRIM(action)) 
       CASE ('write') 
         CALL H5Fcreate_f( TRIM(file), H5F_ACC_TRUNC_F, h5file%id , ierr )
       CASE ('read' ) 
         CALL H5Fopen_f ( TRIM (file), H5F_ACC_RDONLY_F, h5file%id, ierr )  
       CASE ( 'read-write') 
         CALL H5Fopen_f ( TRIM (file), H5F_ACC_RDWR_F, h5file%id, ierr)  
       CASE default 
         ierr =1 
    END SELECT
    IF ( ierr /=0 ) THEN 
       IF (present (error)) then
          error = ierr
       ELSE 
          CALL H5Eprint_f( jerr ) 
          stop
       END IF
    END IF
    ! //' with action '// trim(action), 1 )  
  END SUBROUTINE  qeh5_openfile  

  !--------------------------------------------
  SUBROUTINE qeh5_closefile ( h5file ) 
    !-----------------------------------------
    IMPLICIT NONE 
    TYPE( qeh5_file ),INTENT(INOUT)             :: h5file
    INTEGER                                     :: h5type, ierr 
    !
    CALL H5Iget_type_f(h5file%id,  h5type, ierr ) 
    IF ( h5type == H5I_FILE_F) THEN 
       CALL H5Fclose_f( h5file%id,  ierr )
    ELSE IF ( h5type == H5I_GROUP_F ) THEN 
       CALL H5Gclose_f( h5file%id,  ierr )
    ELSE 
       ierr = 1 
    END IF 
    h5file%name="" 
  END SUBROUTINE qeh5_closefile
 
  !-------------------------------------------------
  SUBROUTINE set_dataset_name(name, h5_dataset) 
    !-----------------------------------------------
    IMPLICIT NONE
    CHARACTER(LEN=256)      :: name
    TYPE (qeh5_dataset)     :: h5_dataset
    h5_dataset%name = TRIM(name)
  END SUBROUTINE set_dataset_name

  !-------------------------------------------------
  SUBROUTINE close_dataset( h5_dataset )
     !-----------------------------------------------
     IMPLICIT NONE
     TYPE ( qeh5_dataset ),INTENT(INOUT)         :: h5_dataset
     !
     INTEGER                                     :: ierr 
     !
     !
     IF   ( ALLOCATED(h5_dataset%filespace%dims)  )   DEALLOCATE (h5_dataset%filespace%dims) 
     IF   ( ALLOCATED(h5_dataset%filespace%maxdims))  DEALLOCATE (h5_dataset%filespace%maxdims) 
     IF   ( ALLOCATED(h5_dataset%filespace%offset))   DEALLOCATE (h5_dataset%filespace%offset)   
     IF   ( ALLOCATED(h5_dataset%filespace%count) )   DEALLOCATE (h5_dataset%filespace%count)  
     IF   ( ALLOCATED(h5_dataset%filespace%stride))   DEALLOCATE (h5_dataset%filespace%stride) 
     IF   ( ALLOCATED(h5_dataset%filespace%block) )   DEALLOCATE (h5_dataset%filespace%block) 
     h5_dataset%filespace%rank=0
     CALL H5Sclose_f(h5_dataset%filespace%id, ierr )
     h5_dataset%filespace%id = -1_HID_T 
     IF (h5_dataset%memspace_ispresent) THEN  
       IF ( ALLOCATED(h5_dataset%memspace%dims)   ) DEALLOCATE (h5_dataset%memspace%dims)
       IF ( ALLOCATED(h5_dataset%memspace%maxdims)) DEALLOCATE (h5_dataset%memspace%maxdims)
       IF ( ALLOCATED(h5_dataset%memspace%offset) ) DEALLOCATE (h5_dataset%memspace%offset)
       IF ( ALLOCATED(h5_dataset%memspace%count)  ) DEALLOCATE (h5_dataset%memspace%count)
       IF ( ALLOCATED(h5_dataset%memspace%stride) ) DEALLOCATE (h5_dataset%memspace%stride)
       IF ( ALLOCATED(h5_dataset%memspace%block)  ) DEALLOCATE (h5_dataset%memspace%block) 
       h5_dataset%memspace_ispresent = .FALSE.
       CALL H5Sclose_f(h5_dataset%memspace%id, ierr )
       !!print '(A, " memspace closer returned ierr =" ,I4)', TRIM(h5_dataset%name) ,  ierr 
       h5_dataset%memspace%id = -1_HID_T
     END IF

     CALL H5Tclose_f( h5_dataset%datatype%id, ierr ) 
     CALL H5Dclose_f( h5_dataset%id, ierr ) 
     h5_dataset%name = "" 
     h5_dataset%datatype%id = -1_HID_T
  END SUBROUTINE close_dataset 

  !-------------------------------------------------------------------
  SUBROUTINE qeh5_open_group( h5_parent , group_name, h5_group) 
    !-----------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(qeh5_file),INTENT(IN)   :: h5_parent
    CHARACTER(LEN=*),INTENT(IN)  :: group_name
    TYPE(qeh5_file),INTENT(OUT)  :: h5_group
    ! 
    INTEGER                      :: ierr, jerr
    INTEGER(HID_T)               :: loc_id, group_id
    loc_id = h5_parent%id
    CALL h5eset_auto_f(0, jerr )
    CALL h5gopen_f ( loc_id, trim(group_name), group_id, ierr  )
    CALL h5eset_auto_f(1, jerr)
    if (ierr  /= 0 )  CALL h5gcreate_f ( loc_id, trim(group_name), group_id, ierr  ) 
    !
    h5_group%name = trim(group_name) 
    h5_group%id = group_id  
  END SUBROUTINE qeh5_open_group

  !----------------------------------------------------------------------------
  SUBROUTINE qeh5_open_dataset ( h5_parent, h5_dataset, action , name, error ) 
    !--------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(qeh5_file), INTENT(IN)          :: h5_parent
    CHARACTER(LEN=*),INTENT(IN)          :: action
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: name 
    INTEGER,OPTIONAL,INTENT(OUT)          :: error
    TYPE(qeh5_dataset) ,INTENT(INOUT)    :: h5_dataset
    !
    INTEGER                              :: ierr, jerr, rank_ 
    LOGICAL                              :: exists
    IF (PRESENT(name) ) h5_dataset%name = TRIM(name)     
    SELECT CASE ( TRIM(action) ) 
       CASE ('write') 
          CALL H5Lexists_f (h5_parent%id,   TRIM(h5_dataset%name), exists, ierr ) 
          IF ( exists )  CALL H5Ldelete_f ( h5_parent%id, TRIM(h5_dataset%name), ierr )  
          CALL H5Dcreate_f ( h5_parent%id, TRIM(h5_dataset%name) , h5_dataset%datatype%id, &
                          h5_dataset%filespace%id, h5_dataset%id, ierr )
    !
    ! 
      CASE ( 'read', 'get_hid' )
         CALL H5Lexists_f (h5_parent%id,   TRIM(h5_dataset%name), exists, ierr )
         IF ( exists ) THEN 
            CALL H5Dopen_f ( h5_parent%id, TRIM(h5_dataset%name), h5_dataset%id, ierr  )
            CALL H5Dget_space_f(h5_dataset%id, h5_dataset%filespace%id, ierr  )
            CALL H5Sget_simple_extent_ndims_f( h5_dataset%filespace%id, rank_ , ierr )
            ALLOCATE( h5_dataset%filespace%dims( rank_ ) , h5_dataset%filespace%maxdims(rank_) )
            h5_dataset%filespace%rank = rank_
            CALL H5Sget_simple_extent_dims_f( h5_dataset%filespace%id, h5_dataset%filespace%dims, &
                                                                    h5_dataset%filespace%maxdims, ierr )
            CALL H5Dget_type_f(h5_dataset%id, h5_dataset%datatype%id, ierr )
         ELSE 
            ierr = -1
         END IF
      CASE default 
         ierr = -1 
    END SELECT  
    IF ( PRESENT (error) )THEN
         error = ierr
      ELSE
         CALL errore ( 'qeh5_open_datase', 'error opening dataset ' // & 
                               h5_parent%name // '/'//  name // ' with action= ' //TRIM(action), ierr )
    END IF
  END SUBROUTINE qeh5_open_dataset                    
  
  !-------------------------------------------------------------------
  SUBROUTINE prepare_dataspace(  h5_dataset, rank, dimensions, mode )
    !----------------------------------------------------------------- 
    IMPLICIT NONE 
    INTEGER,INTENT(IN)                         :: rank
    INTEGER,INTENT(IN)                         :: dimensions(rank)
    TYPE(qeh5_dataset),INTENT(INOUT)           :: h5_dataset
    CHARACTER(1),OPTIONAL,INTENT(IN)           :: mode
    !
    INTEGER                                    :: ierr 
    INTEGER(HID_T)                             :: dtype_hid
    CHARACTER(1)                               :: what_data_space
    !
    what_data_space = 'f'
    IF (PRESENT ( mode ) ) what_data_space = mode
    SELECT CASE (what_data_space) 
      CASE ('m','M') 
         CALL block ( h5_dataset%memspace )
         h5_dataset%memspace_ispresent = .true.
      CASE DEFAULT 
         CALL block ( h5_dataset%filespace )
    END SELECT
    !
    CONTAINS
    !------------------------------------
       SUBROUTINE block ( dataspace )                                          
          IMPLICIT NONE                                                        
          TYPE(qeh5_dataspace )    :: dataspace                                
          IF (ALLOCATED (dataspace%dims) ) DEALLOCATE ( dataspace%dims )       
          ALLOCATE(dataspace%dims(rank) )                                      
          dataspace%rank = rank                                                
          dataspace%dims(1:rank) = dimensions(1:rank)*1_HSIZE_T                
          CALL H5Screate_simple_f( rank, dataspace%dims, dataspace%id, ierr )  
       END SUBROUTINE block                                                    
    !
  END SUBROUTINE prepare_dataspace
  
  !----------------------------------------------------------------------------
  SUBROUTINE qeh5_wplan_real( h5_dataset, real_data, rank, dimensions, mode ) 
    !--------------------------------------------------------------------------
    IMPLICIT NONE 
    REAL(DP),INTENT(IN)                        :: real_data
    INTEGER,INTENT(IN)                         :: rank
    INTEGER,INTENT(IN)                         :: dimensions(rank)
    TYPE(qeh5_dataset),INTENT(INOUT)           :: h5_dataset
    CHARACTER(1),OPTIONAL,INTENT(IN)           :: mode 
    !
    INTEGER                                    :: ierr 
    !
    CALL H5Tcopy_f  ( H5T_IEEE_F64LE  , h5_dataset%datatype%id , ierr  )
    ! 
    IF ( PRESENT(mode) ) THEN 
       CALL prepare_dataspace( h5_dataset, rank, dimensions, mode )
    ELSE 
       CALL prepare_dataspace( h5_dataset, rank, dimensions )
    END IF
  END SUBROUTINE qeh5_wplan_real 
  
  !-----------------------------------------------------------------------------------
  SUBROUTINE qeh5_wplan_complex  ( h5_dataset, complex_data, rank, dimensions, mode )
    !---------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(DP),INTENT(IN)               :: complex_data
    INTEGER,INTENT(IN)                   :: rank
    INTEGER,INTENT(IN)                   :: dimensions(rank)
    TYPE(qeh5_dataset),INTENT(INOUT)     :: h5_dataset
    CHARACTER(1),OPTIONAL,INTENT(IN)     :: mode 
    ! 
    INTEGER                              :: ierr
    INTEGER,DIMENSION(24)                :: dims
    ! 
    CALL H5Tcopy_f  ( H5T_IEEE_F64LE  , h5_dataset%datatype%id  , ierr ) 
    dims(1:rank) = dimensions(1:rank)
    dims(1) = dims(1)*2 
    IF ( PRESENT( mode ) ) THEN 
       CALL prepare_dataspace( h5_dataset, rank, dims(1:rank), mode )
    ELSE 
       CALL prepare_dataspace( h5_dataset, rank, dims(1:rank) )
    END IF 
  END SUBROUTINE qeh5_wplan_complex 
 
  !------------------------------------------------------------------------------------
  SUBROUTINE qeh5_wplan_integer ( h5_dataset, integer_data, rank, dimensions, mode )  
    !----------------------------------------------------------------------------------
    IMPLICIT NONE 
    INTEGER,TARGET,INTENT(IN)            :: integer_data
    INTEGER,INTENT(IN)                   :: rank
    INTEGER,INTENT(IN)                   :: dimensions(rank)
    TYPE( qeh5_dataset ),INTENT(INOUT)   :: h5_dataset
    CHARACTER(1),OPTIONAL,INTENT(IN)     :: mode
    !  
    INTEGER                              :: ierr 
    !
    CALL H5Tcopy_f  ( H5T_STD_I32LE, h5_dataset%datatype%id , ierr  )
    IF ( PRESENT (mode ) ) THEN 
      CALL prepare_dataspace( h5_dataset, rank, dimensions, mode )
    ELSE 
      CALL prepare_dataspace( h5_dataset, rank, dimensions) 
    END IF
  END SUBROUTINE qeh5_wplan_integer
  
  !----------------------------------------------------------
  SUBROUTINE qeh5_read_real ( real_data, h5_dataset)
    !--------------------------------------------------------
    IMPLICIT NONE
    REAL(DP), TARGET, INTENT(INOUT)                  :: real_data(1) 
    TYPE (qeh5_dataset),INTENT(IN)                   :: h5_dataset
    !
    TYPE(C_PTR)                                      :: ptr
    INTEGER(HID_T)                                   :: mem_hid, file_hid
    INTEGER                                          :: ierr  
    ptr = C_LOC(real_data) 
    mem_hid = H5S_ALL_F
    file_hid = H5S_ALL_F
    IF (ALLOCATED(h5_dataset%filespace%offset)) file_hid = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent ) mem_hid = h5_dataset%memspace%id
    CALL H5Dread_f( h5_dataset%id, h5_realdp_type, ptr, ierr, mem_hid, file_hid, H5P_DEFAULT_F )
    !IF ( ierr /=0)  CALL errore( 'qeh5_read_dataset', 'error reading '//TRIM(h5_descriptor%filename), ierr)
  END SUBROUTINE qeh5_read_real

  !----------------------------------------------------------
  SUBROUTINE qeh5_read_complex ( complex_data, h5_dataset)
    !--------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(DP), TARGET, INTENT(INOUT)               :: complex_data(1)
    TYPE (qeh5_dataset),INTENT(IN)                   :: h5_dataset
    INTEGER(HID_T)                                   :: mem_hid, file_hid
    !
    TYPE(C_PTR)                                      :: ptr
    INTEGER                                          :: ierr 
    ptr = C_LOC(complex_data)
    file_hid = H5S_ALL_F
    mem_hid  = H5S_ALL_F
    IF (ALLOCATED(h5_dataset%filespace%offset)) file_hid = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent ) mem_hid = h5_dataset%memspace%id
    CALL H5Dread_f( h5_dataset%id, H5_REALDP_TYPE, ptr, ierr, mem_hid, file_hid, H5P_DEFAULT_F )
    !IF ( ierr /=0)  CALL errore( 'qeh5_read_dataset', 'error reading '//TRIM(h5_descriptor%filename), ierr)
  END SUBROUTINE qeh5_read_complex          

  !-----------------------------------------------------------
  SUBROUTINE qeh5_read_integer ( integer_data, h5_dataset)
    !---------------------------------------------------------
    IMPLICIT NONE
    INTEGER, TARGET, INTENT(INOUT)                   :: integer_data(1)
    TYPE (qeh5_dataset),INTENT(IN)                   :: h5_dataset
    INTEGER(HID_T)                                   :: mem_hid, file_hid
    !
    TYPE(C_PTR)                                      :: ptr
    INTEGER                                          :: ierr
    LOGICAL                                          :: is_valid 
    ptr = C_LOC(integer_data)
    !
    file_hid = H5S_ALL_F
    mem_hid  = H5S_ALL_F
    IF (ALLOCATED(h5_dataset%filespace%offset)) file_hid = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent ) mem_hid = h5_dataset%memspace%id
    CALL H5Dread_f( h5_dataset%id, H5T_NATIVE_INTEGER , ptr, ierr, mem_hid, file_hid, H5P_DEFAULT_F )
    !IF ( ierr /=0)  CALL errore( 'qeh5_read_dataset', 'error reading '//TRIM(h5_descriptor%filename), ierr)
  END SUBROUTINE  qeh5_read_integer

  !--------------------------------------------------------
  SUBROUTINE qeh5_write_real( real_data, h5_dataset ) 
    !------------------------------------------------------
    IMPLICIT NONE
    REAL(DP), TARGET, INTENT(INOUT)          ::  real_data(1) 
    TYPE(qeh5_dataset),INTENT(IN)            ::  h5_dataset  
    ! 
    TYPE(C_PTR)                              ::  buf
    INTEGER                                  ::  ierr, jerr 
    INTEGER(HID_T)                           ::  memspace_, filespace_
    ! 
    buf = C_LOC(real_data)  
    filespace_ = H5S_ALL_F
    memspace_  = H5S_ALL_F
    IF ( ALLOCATED (h5_dataset%filespace%offset) ) filespace_ = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id 
    CALL H5Dwrite_f ( h5_dataset%id, H5_REALDP_TYPE, buf  , ierr, memspace_,&
                      filespace_, H5P_DEFAULT_F )
  END SUBROUTINE qeh5_write_real  

  !-------------------------------------------------------
  SUBROUTINE qeh5_write_real_2( real_data, h5_dataset ) 
    !-----------------------------------------------------
    IMPLICIT NONE
    REAL(DP), TARGET, INTENT(INOUT)          ::  real_data(1,1) 
    TYPE(qeh5_dataset),INTENT(IN)            ::  h5_dataset  
    ! 
    TYPE(C_PTR)                              ::  buf
    INTEGER                                  ::  ierr, jerr 
    INTEGER(HID_T)                           ::  memspace_, filespace_
    ! 
    buf = C_LOC(real_data)  
    filespace_ = H5S_ALL_F
    memspace_  = H5S_ALL_F
    IF( ALLOCATED (h5_dataset%filespace%offset)) filespace_ = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id 
    CALL H5Dwrite_f ( h5_dataset%id, H5_REALDP_TYPE, buf  , ierr, memspace_,&
                      filespace_, H5P_DEFAULT_F )
  END SUBROUTINE qeh5_write_real_2  
  

  SUBROUTINE qeh5_write_real_3( real_data, h5_dataset ) 
    IMPLICIT NONE
    REAL(DP), TARGET, INTENT(INOUT)          ::  real_data(1,1,1) 
    TYPE(qeh5_dataset),INTENT(IN)            ::  h5_dataset  
    ! 
    TYPE(C_PTR)                              ::  buf
    INTEGER                                  ::  ierr, jerr 
    INTEGER(HID_T)                           ::  memspace_, filespace_
    ! 
    buf = C_LOC(real_data)  
    filespace_ = H5S_ALL_F
    memspace_  = H5S_ALL_F
    IF( ALLOCATED (h5_dataset%filespace%offset)) filespace_ = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id 
    CALL H5Dwrite_f ( h5_dataset%id, H5_REALDP_TYPE, buf  , ierr, memspace_,&
                      filespace_, H5P_DEFAULT_F )
  END SUBROUTINE qeh5_write_real_3  

  !--------------------------------------------------------------
  SUBROUTINE qeh5_write_complex( complex_data, h5_dataset )
   !-------------------------------------------------------------
   IMPLICIT NONE
   COMPLEX(DP), TARGET, INTENT(IN)          ::  complex_data(1)
   TYPE(qeh5_dataset),INTENT(IN)            ::  h5_dataset
   ! 
   TYPE(C_PTR)                              ::  buf
   INTEGER                                  ::  ierr, jerr
   INTEGER(HID_T)                           ::  memspace_, filespace_
   ! 
   buf = C_LOC(complex_data)
   memspace_ =  H5S_ALL_F
   filespace_ = H5S_ALL_F
   IF ( ALLOCATED (h5_dataset%filespace%offset) ) filespace_ = h5_dataset%filespace%id
   IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id
   CALL H5Dwrite_f ( h5_dataset%id, H5_REALDP_TYPE, buf  , ierr, memspace_,&
                     filespace_, H5P_DEFAULT_F )
   
  END SUBROUTINE qeh5_write_complex

  !------------------------------------------------------------
  SUBROUTINE qeh5_write_complex_2( complex_data, h5_dataset )
   !-----------------------------------------------------------
   IMPLICIT NONE
   COMPLEX(DP), TARGET, INTENT(INOUT)       ::  complex_data(1,1)
   TYPE(qeh5_dataset),INTENT(IN)            ::  h5_dataset
   ! 
   TYPE(C_PTR)                              ::  buf
   INTEGER                                  ::  ierr, jerr
   INTEGER(HID_T)                           ::  memspace_, filespace_
   ! 
   buf = C_LOC(complex_data)
   memspace_ =  H5S_ALL_F
   filespace_ = H5S_ALL_F
   IF (ALLOCATED (h5_dataset%filespace%offset)) filespace_ = h5_dataset%filespace%id
   IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id
   CALL H5Dwrite_f ( h5_dataset%id, H5_REALDP_TYPE, buf  , ierr, memspace_,&
                     filespace_ , H5P_DEFAULT_F )
  END SUBROUTINE qeh5_write_complex_2


  SUBROUTINE qeh5_write_complex_3( complex_data, h5_dataset )
    IMPLICIT NONE
    COMPLEX(DP), TARGET, INTENT(INOUT)       ::  complex_data(1,1,1)
    TYPE(qeh5_dataset),INTENT(IN)            ::  h5_dataset
    ! 
    TYPE(C_PTR)                              ::  buf
    INTEGER                                  ::  ierr, jerr
    INTEGER(HID_T)                           ::  memspace_, filespace_
    ! 
    buf = C_LOC(complex_data)
    memspace_ =  H5S_ALL_F
    filespace_ = H5S_ALL_F
    IF (ALLOCATED (h5_dataset%filespace%offset) ) filespace_ = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id
   
    CALL H5Dwrite_f ( h5_dataset%id, h5_realdp_type, buf  , ierr, memspace_,&
                     filespace_ , H5P_DEFAULT_F  ) 
  END SUBROUTINE qeh5_write_complex_3

  !----------------------------------------------------------
  SUBROUTINE qeh5_write_integer( integer_data, h5_dataset )
   !---------------------------------------------------------
   IMPLICIT NONE
   INTEGER, TARGET   , INTENT(INOUT)        ::  integer_data(1)
   TYPE(qeh5_dataset), INTENT(IN)           ::  h5_dataset
   ! 
   TYPE(C_PTR)                              ::  buf
   INTEGER                                  ::  ierr, jerr
   INTEGER(HID_T)                           ::  memspace_, filespace_
   ! 
   buf = C_LOC(integer_data)
   memspace_ =  H5S_ALL_F
   filespace_ = H5S_ALL_F
   IF (ALLOCATED (h5_dataset%filespace%offset)) filespace_ = h5_dataset%filespace%id
   IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id
   CALL H5Dwrite_f ( h5_dataset%id, H5T_NATIVE_INTEGER , buf  , ierr, memspace_,&
                     filespace_ , H5P_DEFAULT_F   )
   
  END SUBROUTINE qeh5_write_integer

  !-------------------------------------------------------------
  SUBROUTINE qeh5_write_integer_2( integer_data, h5_dataset )
   !------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, TARGET   , INTENT(INOUT)        ::  integer_data(1,1)
   TYPE(qeh5_dataset), INTENT(IN)           ::  h5_dataset
   ! 
   TYPE(C_PTR)                              ::  buf
   INTEGER                                  ::  ierr, jerr
   INTEGER(HID_T)                           ::  memspace_, filespace_
   ! 
   buf = C_LOC(integer_data)
   memspace_ =  H5S_ALL_F
   filespace_ = H5S_ALL_F
   IF (ALLOCATED (h5_dataset%filespace%offset)) filespace_ = h5_dataset%filespace%id
   IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id
   CALL H5Dwrite_f ( h5_dataset%id, H5T_NATIVE_INTEGER , buf  , ierr, memspace_,&
                     filespace_ , H5P_DEFAULT_F   )
  END SUBROUTINE qeh5_write_integer_2

  !--------------------------------------------------------------
  SUBROUTINE qeh5_write_integer_3( integer_data, h5_dataset )
    !------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, TARGET   , INTENT(INOUT)        ::  integer_data(1,1,1)
    TYPE(qeh5_dataset), INTENT(IN)           ::  h5_dataset
   ! 
    TYPE(C_PTR)                              ::  buf
    INTEGER                                  ::  ierr, jerr
    INTEGER(HID_T)                           ::  memspace_, filespace_
   ! 
    buf = C_LOC(integer_data)
    memspace_ =  H5S_ALL_F
    filespace_ = H5S_ALL_F
    IF (ALLOCATED (h5_dataset%filespace%offset)) filespace_ = h5_dataset%filespace%id
    IF ( h5_dataset%memspace_ispresent)  memspace_ = h5_dataset%memspace%id
    CALL H5Dwrite_f ( h5_dataset%id, H5T_NATIVE_INTEGER , buf  , ierr, memspace_,&
                     filespace_ , H5P_DEFAULT_F   )
  END SUBROUTINE qeh5_write_integer_3
!

  !-------------------------------------------------------
  SUBROUTINE add_attribute_string ( h5_hid, name, text) 
   !------------------------------------------------------
   IMPLICIT NONE
   INTEGER(HID_T),INTENT(IN)    :: h5_hid
   CHARACTER(LEN=*)             :: name, text 
   TARGET                       :: text
   ! 
   INTEGER(HSIZE_T)             :: text_length
   INTEGER(HID_T)               :: sp_id, attr_id, string_type
   LOGICAL                      :: exists
   INTEGER                      :: ierr
   TYPE(C_PTR)                  :: buf
   ! 
   text_length = LEN(TRIM(text))*1_HSIZE_T
   buf = c_loc(text)
   CALL H5Screate_f( H5S_SCALAR_F, sp_id, ierr )
   CALL H5Tcopy_f(H5T_FORTRAN_S1, string_type, ierr )
   CALL H5Tset_size_f( string_type,text_length, ierr )
   CALL H5Aexists_by_name_f(h5_hid, '.', trim(name), exists, ierr )
   if (exists ) call H5Adelete_by_name_f( h5_hid, '.', trim(name), ierr )
   !
   call H5Acreate_f( h5_hid, trim(name),  string_type, sp_id, attr_id, ierr)
   call H5Awrite_f (attr_id, string_type, buf, ierr )
   
   CALL H5Sclose_f(sp_id,   ierr )
   CALL H5Aclose_f(attr_id, ierr )
   !
  END SUBROUTINE add_attribute_string

  !---------------------------------------------------------------------
  SUBROUTINE add_array_attribute_i ( h5_hid, name, value, rank, dims )
     !--------------------------------------------------------------------
     IMPLICIT NONE 
     INTEGER(HID_T)  ,INTENT(IN)  :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)  ::   name
     INTEGER,TARGET, INTENT(IN)   :: value(:)
     INTEGER,INTENT(IN)           :: rank, dims(:)
     CALL add_attribute_i ( h5_hid, name, value(1), rank, dims) 
  END SUBROUTINE add_array_attribute_i

  !---------------------------------------------------------------
  SUBROUTINE add_attribute_i ( h5_hid, name, value, rank, dims )
     !------------------------------------------------------------
     IMPLICIT NONE
     INTEGER(HID_T)  ,INTENT(IN)  :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)  ::   name
     INTEGER,TARGET, INTENT(IN)   :: value
     INTEGER,OPTIONAL, INTENT(IN)           :: rank, dims(:)
     
     INTEGER(HID_T)               ::   loc_id, sp_id, attr_id, data_type, mem_type
     INTEGER(HSIZE_T),ALLOCATABLE ::   hdims(:)
     INTEGER                      ::   ierr
     LOGICAL                      ::   exists
     TYPE(C_PTR)                  ::   buf 
     INTEGER                      ::   i
     !
     buf = C_LOC(value) 
     
     IF (PRESENT(rank)) THEN
        ALLOCATE(hdims(rank)) 
        DO i =1, rank 
           hdims(i) =dims(i)*1_HSIZE_T 
        END DO
        CALL H5Tarray_create_f( H5T_STD_I32LE , rank, hdims, data_type, ierr )
        CALL H5Tarray_create_f( H5T_NATIVE_INTEGER, rank, hdims, mem_type, ierr  ) 
     ELSE 
        CALL H5Tcopy_f( H5T_STD_I32LE,data_type, ierr)  
        CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr ) 
     END IF
     
     loc_id = h5_hid
  
     CALL H5Screate_f      ( H5S_SCALAR_F, sp_id, ierr )
     CALL H5Aexists_by_name_f(loc_id, '.', trim(name), exists, ierr )
     if (exists ) call H5Adelete_by_name_f( loc_id, '.', trim(name), ierr )
     !
     CALL H5Acreate_f( loc_id, TRIM(name), data_type, sp_id, attr_id, ierr )
     CALL H5Awrite_f (attr_id, mem_type, buf , ierr )
     
     CALL H5Tclose_f(data_type, ierr) 
     CALL H5Tclose_f(mem_type, ierr)
     CALL H5Sclose_f(sp_id, ierr )
     CALL H5Aclose_f(attr_id, ierr )
  END SUBROUTINE add_attribute_i

  !-------------------------------------------------------------------
  SUBROUTINE add_array_attribute_r ( h5_hid, name, value, rank, dims)
     !---------------------------------------------------------------
     IMPLICIT NONE
     INTEGER(HID_T),INTENT(IN)    :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)  :: name
     REAL(8),TARGET, INTENT(IN)   :: value(:)
     INTEGER,INTENT(IN)           :: rank, dims(:)
     CALL add_attribute_r( h5_hid, name, value(1), rank, dims )
  END SUBROUTINE add_array_attribute_r

  !---------------------------------------------------------------
  SUBROUTINE add_attribute_r ( h5_hid, name, value, rank, dims )
     !-----------------------------------------------------------
     IMPLICIT NONE
     INTEGER(HID_T),INTENT(IN)    :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)  :: name
     REAL(8),TARGET, INTENT(IN)   :: value
     INTEGER,OPTIONAL,INTENT(IN)           :: rank, dims(:)
     !
     INTEGER(HID_T)               ::   loc_id, sp_id, attr_id, data_type, real8_type, mem_type
     LOGICAL                      ::   exists
     INTEGER(HSIZE_T),ALLOCATABLE ::   hdims(:)
     INTEGER                      ::   ierr
     TYPE(C_PTR)                  ::   buf 
     INTEGER                      ::   i
     !
     loc_id = h5_hid
     buf = C_LOC(value) 
     IF (PRESENT(rank) ) THEN
        ALLOCATE(hdims(rank)) 
        DO i =1, rank 
           hdims(i) =dims(i)*1_HSIZE_T 
        END DO
        CALL H5Tarray_create_f(H5T_IEEE_F64LE, rank, hdims, data_type, ierr)  
        CALL H5Tarray_create_f(H5_REALDP_TYPE, rank, hdims, mem_type, ierr ) 
     ELSE 
        CALL H5Tcopy_f(H5T_IEEE_F64LE, data_type, ierr)
        CALL H5Tcopy_f(H5_REALDP_TYPE, mem_type, ierr)  
     END IF 
     CALL H5Screate_f(  H5S_SCALAR_F, sp_id, ierr )
     CALL H5Aexists_by_name_f(loc_id, '.', trim(name), exists, ierr )
     if (exists ) call H5Adelete_by_name_f( loc_id, '.', trim(name), ierr )
     CALL H5Acreate_f( loc_id, TRIM(name), data_type, sp_id, attr_id, ierr )
     !
     !  
     CALL H5Awrite_f (attr_id, mem_type, buf , ierr )
     !
     CALL H5Tclose_f(mem_type,  ierr)
     CALL H5Tclose_f(data_type, ierr)
     CALL H5Sclose_f(sp_id,     ierr)
     CALL H5Aclose_f(attr_id,   ierr)
  END SUBROUTINE add_attribute_r

  !-------------------------------------------------------------------------------------
  SUBROUTINE read_integer_array_attribute ( h5_hid, attribute, value, rank, dimensions) 
     !----------------------------------------------------------------------------------
     IMPLICIT NONE 
     INTEGER(HID_T) , INTENT(IN)                  :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)                  :: attribute
     INTEGER,TARGET,INTENT(OUT)                   :: value(:) 
     INTEGER, TARGET,INTENT(IN)                   :: rank, dimensions(:)
     !
     !
     CALL read_integer_attribute( h5_hid, attribute, value(1), rank, dimensions) 
  END SUBROUTINE read_integer_array_attribute     

  !-------------------------------------------------------------------------------
  SUBROUTINE read_integer_attribute ( h5_hid, attribute, value, rank, dimensions) 
     !----------------------------------------------------------------------------
     IMPLICIT NONE 
     INTEGER(HID_T) , INTENT(IN)                  :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)                  :: attribute
     INTEGER,TARGET,INTENT(OUT)                   :: value 
     INTEGER,OPTIONAL, INTENT(IN)                 :: rank, dimensions(:)
     ! 
     INTEGER(HID_T)                               :: loc_id, attr_id, data_type, mem_type
     INTEGER(HSIZE_T),ALLOCATABLE                 :: hdims(:)
     INTEGER                                      :: i, ierr
     TYPE(C_PTR)                                  :: buf
     !
     loc_id = h5_hid
     buf = C_LOC(value) 
     IF (PRESENT (rank ) ) THEN 
        ALLOCATE(hdims(rank)) 
        DO i =1, rank 
          hdims(i) =dimensions(i)*1_HSIZE_T 
        END DO
        CALL H5Tarray_create_f(H5T_NATIVE_INTEGER, rank, hdims, mem_type, ierr  )
     ELSE 
        CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr )
     END IF 
     !
     ! 
     CALL H5Aopen_by_name_f( loc_id, '.', TRIM(attribute), attr_id, ierr) 
     CALL H5Aread_f( attr_id, mem_type, buf, ierr )
     !
     CALL H5Tclose_f( mem_type, ierr) 
     CALL H5Aclose_f( attr_id, ierr) 
     !
  END SUBROUTINE read_integer_attribute         

  !------------------------------------------------------------------------------------
  SUBROUTINE read_real_array_attribute ( h5_hid, attribute, value, rank, dimensions )
     !---------------------------------------------------------------------------------
     IMPLICIT NONE 
     INTEGER(HID_T), INTENT(IN)                   :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)                  :: attribute
     REAL(8),TARGET,INTENT(OUT)                   :: value(:) 
     INTEGER, INTENT(IN)                          :: rank, dimensions(:)
     ! 
     CALL read_real_attribute ( h5_hid, attribute, value(1), rank, dimensions)
  END SUBROUTINE read_real_array_attribute

  !------------------------------------------------------------------------------  
  SUBROUTINE read_real_attribute ( h5_hid, attribute, value, rank, dimensions) 
     !--------------------------------------------------------------------------
     IMPLICIT NONE 
     INTEGER(HID_T), INTENT(IN)                   :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)                  :: attribute
     REAL(8),TARGET,INTENT(OUT)                   :: value 
     INTEGER,OPTIONAL, INTENT(IN)                 :: rank, dimensions(:)
     ! 
     INTEGER(HID_T)                               :: loc_id, attr_id, data_type, mem_type
     INTEGER(HSIZE_T),ALLOCATABLE                 :: hdims(:)
     INTEGER                                      :: i, ierr
     TYPE(C_PTR)                                  :: buf
     
     !
     loc_id = h5_hid
     buf = C_LOC(value) 
     IF (PRESENT(rank) )  THEN
        ALLOCATE(hdims(rank)) 
        DO i =1, rank 
           hdims(i) =dimensions(i)*1_HSIZE_T 
        END DO
        CALL H5Tarray_create_f( H5_REALDP_TYPE, rank, hdims, mem_type, ierr  ) 
     ELSE 
        CALL H5Tcopy_f(H5_REALDP_TYPE , mem_type, ierr) 
     END IF
     !
     ! 
     CALL H5Aopen_by_name_f( loc_id, '.', TRIM(attribute), attr_id, ierr) 
     CALL H5Aread_f( attr_id, mem_type, buf, ierr )
     !
     CALL H5Tclose_f( mem_type, ierr) 
     CALL H5Aclose_f( attr_id, ierr) 
     !
  END SUBROUTINE read_real_attribute      
  
  !-----------------------------------------------------------------
  SUBROUTINE read_string_attribute(h5_hid, attribute, text, maxlen) 
     !--------------------------------------------------------------
     IMPLICIT NONE 
     INTEGER(HID_T), INTENT(IN)                   :: h5_hid
     CHARACTER(LEN=*),INTENT(IN)                  :: attribute
     INTEGER,INTENT(IN)                           :: maxlen
     CHARACTER(LEN=*),INTENT(OUT)                 :: text 
     ! 
     INTEGER(HID_T)                               :: loc_id, attr_id, data_type
     INTEGER(HSIZE_T)                             :: text_length
     INTEGER                                      :: i, ierr
     TYPE(C_PTR)                                  :: buf
     CHARACTER,TARGET,ALLOCATABLE                 :: chars(:) 
     !
     !
     ! 
     loc_id = h5_hid
     text=""
     ALLOCATE(chars(maxlen))
     buf = C_LOC(chars) 
     CALL H5Aopen_by_name_f( loc_id, '.', trim(attribute), attr_id, ierr) 
     CALL H5Aget_type_f( attr_id, data_type, ierr) 
     CALL H5Tget_size_f( data_type, text_length, ierr) 
     IF ( text_length > maxlen*1_HSIZE_T) &
                        CALL infomsg (TRIM (attribute) // ' text too long will be truncated on reading') 
     CALL H5Aread_f( attr_id, data_type, buf, ierr )
     DO I = 1, maxlen
        IF (i > text_length) EXIT 
        text(i:i) = chars(i)
     END DO
     !
     DEALLOCATE(chars) 
     buf = C_NULL_PTR
     CALL H5Tclose_f( data_type, ierr) 
     CALL H5Aclose_f( attr_id, ierr) 
     !
  END SUBROUTINE read_string_attribute 
  
  !--------------------------------------------------------------------
  SUBROUTINE set_hyperslab (dataspace, offset, count, stride, block )
     !-----------------------------------------------------------------
     IMPLICIT NONE 
     TYPE(qeh5_dataspace),INTENT(INOUT)            :: dataspace
     INTEGER,DIMENSION(:),INTENT(IN)               :: offset, count
     INTEGER,DIMENSION(:),OPTIONAL, INTENT(IN)     :: stride, block
     ! 
     INTEGER                             :: rank, ierr 
     rank = dataspace%rank
     !
     IF (ALLOCATED(dataspace%offset) ) DEALLOCATE (dataspace%offset) 
     IF (ALLOCATED(dataspace%count ) ) DEALLOCATE (dataspace%count) 
     IF (ALLOCATED(dataspace%stride) ) DEALLOCATE (dataspace%stride) 
     IF (ALLOCATED(dataspace%block ) ) DEALLOCATE (dataspace%block) 
     ALLOCATE ( dataspace%offset(rank), dataspace%count(rank), dataspace%stride(rank), dataspace%block(rank) )
     !
     dataspace%offset(1:rank) = offset(1:rank) * 1_HSIZE_T
     dataspace%count (1:rank) = count (1:rank) * 1_HSIZE_T   
     IF (PRESENT(stride) )  THEN
        dataspace%stride(1:rank) = stride(1:rank) * 1_HSIZE_T
     ELSE 
        dataspace%stride(1:rank)  =  1_HSIZE_T
     END IF 
     IF (PRESENT( block ) ) THEN 
        dataspace%block (1:rank) = block (1:rank) * 1_HSIZE_T
     ELSE 
        dataspace%block (1:rank)  =  1_HSIZE_T
     END IF
     CALL H5Sselect_hyperslab_f( dataspace%id,  H5S_SELECT_SET_F, dataspace%offset, dataspace%count, &
                                 ierr, dataspace%stride, dataspace%block )    
  END SUBROUTINE set_hyperslab
   
  !------------------------------------------------------------------------------- 
  SUBROUTINE qeh5_set_memory_hyperslab (dataset, offset, count, stride, block ) 
     !----------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(qeh5_dataset),INTENT(INOUT)              :: dataset
     INTEGER,DIMENSION(:),INTENT(IN)               :: offset, count
     INTEGER,DIMENSION(:),OPTIONAL, INTENT(IN)     :: stride, block
     CALL set_hyperslab (dataset%memspace, offset, count, stride, block )
     END SUBROUTINE qeh5_set_memory_hyperslab
     ! 
     SUBROUTINE qeh5_set_file_hyperslab (dataset, offset, count, stride, block )
     IMPLICIT NONE
     TYPE(qeh5_dataset),INTENT(INOUT)            :: dataset
     INTEGER,DIMENSION(:),INTENT(IN)               :: offset, count
     INTEGER,DIMENSION(:),OPTIONAL, INTENT(IN)     :: stride, block
     CALL set_hyperslab (dataset%filespace, offset, count, stride, block )
  END SUBROUTINE qeh5_set_file_hyperslab
#endif
END MODULE qeh5_base_module
 
