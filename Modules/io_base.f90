!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains some common subroutines used to read and write
  ! ... data produced by the Quantum ESPRESSO package
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(iotk_attlenx)  :: attr
  !
  PRIVATE
  PUBLIC :: write_wfc, read_wfc
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, ik, nk, ispin, nspin, wfc, ngw,   &
                          gamma_only, nbnd, igl, ngwl, filename, scalef, &
                          ionode_in_group, root_in_group, intra_group_comm)
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE mp,         ONLY : mp_size, mp_rank, mp_max
      !
#if defined(__HDF5)
      USE hdf5_qe,    ONLY : prepare_for_writing_final, add_attributes_hdf5, &
           write_evc, h5fclose_f, hdf5_type              
      USE HDF5
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wfc(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=*),   INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      LOGICAL,            INTENT(IN) :: ionode_in_group
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm
      !
      INTEGER                  :: j, ierr
      INTEGER                  :: igwx, npwx, npol
      INTEGER                  :: me_in_group, nproc_in_group, my_group
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
#if defined(__HDF5) 
      TYPE (hdf5_type),ALLOCATABLE    :: h5_write_desc
      ! 
      IF ( ionode_in_group ) ALLOCATE (h5_write_desc) 
#endif 
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      npol = 1
      IF ( nspin == 4 ) npol = 2
      npwx = SIZE( wfc, 1 ) / npol
      ALLOCATE( wtmp( MAX( npol*igwx, 1 ) ) )
      !
      wtmp = 0.0_DP
      !
      IF ( ionode_in_group ) THEN
#if defined  __HDF5
         CALL prepare_for_writing_final ( h5_write_desc, 0, &
              TRIM(filename)//'.hdf5',ik, ADD_GROUP = .false.)
         CALL add_attributes_hdf5(h5_write_desc, ngw,"ngw",ik)
         CALL add_attributes_hdf5(h5_write_desc, gamma_only,"gamma_only",ik)
         CALL add_attributes_hdf5(h5_write_desc, igwx,"igwx",ik)
         CALL add_attributes_hdf5(h5_write_desc, nbnd,"nbnd",ik)
         CALL add_attributes_hdf5(h5_write_desc, ik,"ik",ik)
         CALL add_attributes_hdf5(h5_write_desc, nk,"nk",ik)
         CALL add_attributes_hdf5(h5_write_desc, ispin,"ispin",ik)
         CALL add_attributes_hdf5(h5_write_desc, nspin,"nspin",ik)
         CALL add_attributes_hdf5(h5_write_desc, scalef,"scale_factor",ik)
         !
#else
         !
         CALL iotk_open_write( iuni, FILE = TRIM(filename)//'.dat', &
              ROOT="WFC", BINARY = .TRUE. )
         !
         CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "igwx",         igwx )
         CALL iotk_write_attr( attr, "gamma_only",   gamma_only )
         CALL iotk_write_attr( attr, "nbnd",         nbnd )
         CALL iotk_write_attr( attr, "ik",           ik )
         CALL iotk_write_attr( attr, "nk",           nk )
         CALL iotk_write_attr( attr, "ispin",        ispin )
         CALL iotk_write_attr( attr, "nspin",        nspin )
         CALL iotk_write_attr( attr, "scale_factor", scalef )
         !
         CALL iotk_write_empty( iuni, "INFO", attr )
#endif
         !
      END IF
      !
      DO j = 1, nbnd
         !
         IF ( npol == 2 ) THEN
            !
            ! Quick-and-dirty noncolinear case - mergewf should be modified
            !
            CALL mergewf( wfc(1:npwx,       j), wtmp(1:igwx),       ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            CALL mergewf( wfc(npwx+1:2*npwx,j), wtmp(igwx+1:2*igwx), ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
         ELSE
            !
            CALL mergewf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                 nproc_in_group, root_in_group, intra_group_comm )
            !
         END IF
         !
         IF ( ionode_in_group ) &
#if defined(__HDF5)
            CALL write_evc(h5_write_desc, j, wtmp(1:npol*igwx), ik) 
#else
            CALL iotk_write_dat &
              ( iuni, "evc" // iotk_index( j ), wtmp(1:npol*igwx))
#endif
         !
      END DO
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL h5fclose_f(h5_write_desc%file_id, ierr)
         DEALLOCATE ( h5_write_desc)  
#else 
         CALL iotk_close_write( iuni )
#endif
     END IF
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, ik, nk, ispin, nspin, wfc, ngw, nbnd, &
                         igl, ngwl, filename, scalef, &
                         ionode_in_group, root_in_group, intra_group_comm, &
                         ierr )
      ! if ierr is present, return 0 if everything is ok, /= 0 if not
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf
      USE mp,        ONLY : mp_bcast, mp_size, mp_rank, mp_max
      !
#if defined  __HDF5
      USE hdf5_qe,  ONLY  : prepare_for_reading_final, read_attributes_hdf5, read_evc, &
                            h5fclose_f, hdf5_type
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      COMPLEX(DP),        INTENT(OUT)   :: wfc(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      CHARACTER(LEN=*),   INTENT(IN)    :: filename
      REAL(DP),           INTENT(OUT)   :: scalef
      LOGICAL,            INTENT(IN)    :: ionode_in_group
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      COMPLEX(DP), ALLOCATABLE          :: wtmp(:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx, igwx_, npwx, npol, ik_, nk_
      INTEGER                           :: me_in_group, nproc_in_group
#if defined(__HDF5)
      TYPE (hdf5_type),ALLOCATABLE      :: h5_read_desc
      ! 
      if (ionode_in_group ) ALLOCATE ( h5_read_desc) 
#endif  
      !
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      !
#if !defined __HDF5
      IF ( ionode_in_group ) CALL iotk_open_read( iuni, &
           FILE = TRIM(filename)//'.dat', BINARY = .TRUE., IERR = ierr_ )
      CALL mp_bcast( ierr_, root_in_group, intra_group_comm )
      IF ( PRESENT(ierr) ) THEN
         ierr = ierr_
         IF ( ierr /= 0 ) RETURN
      ELSE
         CALL errore( 'read_wfc ', &
              'cannot open restart file for reading', ierr_ )
      END IF
#endif
      IF ( ionode_in_group ) THEN
          !
#if defined  __HDF5
         IF ( PRESENT (ierr )) THEN 
            CALL prepare_for_reading_final(h5_read_desc, 0, &
               TRIM(filename)//'.hdf5',KPOINT = ik, IERR = ierr )
            IF (ierr /= 0 ) RETURN 
         ELSE
           CALL prepare_for_reading_final(h5_read_desc, 0, &
               TRIM(filename)//'.hdf5',KPOINT = ik)
         END IF 
         CALL read_attributes_hdf5(h5_read_desc, ngw,"ngw",ik)
         CALL read_attributes_hdf5(h5_read_desc, nbnd,"nbnd",ik)
         CALL read_attributes_hdf5(h5_read_desc, ik_,"ik",ik)
         CALL read_attributes_hdf5(h5_read_desc, nk_,"ik",ik)
         CALL read_attributes_hdf5(h5_read_desc, ispin,"ispin",ik)
         CALL read_attributes_hdf5(h5_read_desc, nspin,"nspin",ik)
         CALL read_attributes_hdf5(h5_read_desc, igwx_,"igwx",ik)
         CALL read_attributes_hdf5(h5_read_desc, scalef,"scale_factor",ik)
#else
         CALL iotk_scan_empty( iuni, "INFO", attr )
          !
         CALL iotk_scan_attr( attr, "ngw",          ngw )
          CALL iotk_scan_attr( attr, "nbnd",         nbnd )
          CALL iotk_scan_attr( attr, "ik",           ik_ )
          CALL iotk_scan_attr( attr, "nk",           nk_ )
          CALL iotk_scan_attr( attr, "ispin",        ispin )
          CALL iotk_scan_attr( attr, "nspin",        nspin )
          CALL iotk_scan_attr( attr, "igwx",         igwx_ )
          CALL iotk_scan_attr( attr, "scale_factor", scalef )
          !
#endif
      END IF
      !
      CALL mp_bcast( ngw,    root_in_group, intra_group_comm )
      CALL mp_bcast( nbnd,   root_in_group, intra_group_comm )
      CALL mp_bcast( ik_,    root_in_group, intra_group_comm )
      CALL mp_bcast( nk_,    root_in_group, intra_group_comm )
      CALL mp_bcast( ispin,  root_in_group, intra_group_comm )
      CALL mp_bcast( nspin,  root_in_group, intra_group_comm )
      CALL mp_bcast( igwx_,  root_in_group, intra_group_comm )
      CALL mp_bcast( scalef, root_in_group, intra_group_comm )
      !
      npol = 1
      IF ( nspin == 4 ) npol = 2
      ALLOCATE( wtmp( npol*MAX( igwx_, igwx ) ) )
      npwx = SIZE( wfc, 1 ) / npol
      !
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wfc, 2 ) ) THEN
            !
            IF ( ionode_in_group ) THEN 
#if defined __HDF5
               CALL read_evc(h5_read_desc, j, wtmp(1:npol*igwx_),ik)
#else
               CALL iotk_scan_dat( iuni, &
                                   "evc" // iotk_index(j), wtmp(1:npol*igwx_) )
 
#endif
               !
               IF ( igwx > igwx_ ) wtmp((npol*igwx_+1):npol*igwx) = 0.0_DP
               !
            END IF
            !
            IF ( npol == 2 ) THEN
               CALL splitwf( wfc(1:npwx,       j), wtmp(1:igwx_       ),   &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
               CALL splitwf( wfc(npwx+1:2*npwx,j), wtmp(igwx_+1:2*igwx_),  &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
            ELSE
               CALL splitwf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                    nproc_in_group, root_in_group, intra_group_comm )
            END IF
            !
         END IF
         !
      END DO
      !
#if !defined (__HDF5)
      IF ( ionode_in_group ) CALL iotk_close_read( iuni )
#else
      IF ( ionode_in_group ) THEN 
         CALL h5fclose_f(h5_read_desc%file_id, ierr_)
         DEALLOCATE (h5_read_desc)
      END IF
#endif
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !        
  END MODULE io_base
