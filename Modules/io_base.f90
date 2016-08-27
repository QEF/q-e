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
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, &
                        current_fmt_version => qexml_version
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
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
    SUBROUTINE write_wfc( iuni, ik, nk, ispin, nspin, wf0, ngw,   &
                          gamma_only, nbnd, igl, ngwl, filename, scalef, &
                          ionode, root_in_group, intra_group_comm)
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE mp,         ONLY : mp_get, mp_size, mp_rank, mp_sum, mp_max
      !
#if defined  __HDF5
      USE hdf5_qe,    ONLY : prepare_for_writing_final, add_attributes_hdf5, write_evc, h5fclose_f, evc_hdf5_write              
      USE mp_global,    ONLY : inter_pool_comm, world_comm
      USE HDF5
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wf0(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=256), INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      LOGICAL,            INTENT(IN) :: ionode
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm
      !
      INTEGER                  :: j, ierr
      INTEGER                  :: igwx
      INTEGER                  :: me_in_group, nproc_in_group, my_group
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
#if defined __HDF5
      CHARACTER(LEN=256) :: filename_hdf5
      INTEGER            :: gammaonly
#endif

      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
      !
      wtmp = 0.0_DP
      !
      IF ( ionode ) THEN
#if defined  __HDF5
      filename_hdf5=trim(tmp_dir) //"evc.hdf5"
      CALL prepare_for_writing_final(evc_hdf5_write,inter_pool_comm,filename_hdf5,ik)
      CALL add_attributes_hdf5(evc_hdf5_write,ngw,"ngw",ik)
      gammaonly = 0 
      IF (gamma_only) gammaonly = 1 
      CALL add_attributes_hdf5(evc_hdf5_write,gammaonly,"gamma_only",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,igwx,"igwx",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,nbnd,"nbnd",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,ik,"ik",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,nk,"nk",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,ispin,"ispin",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,nspin,"nspin",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,scalef,"scale_factor",ik)
         !
#else

         !
         CALL iotk_open_write( iuni, FILE = TRIM( filename ), ROOT="WFC", BINARY = .TRUE. )
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
         CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_in_group, &
              nproc_in_group, root_in_group, intra_group_comm )
         !
#ifdef __HDF5
         IF (ionode ) CALL write_evc(evc_hdf5_write,j, wtmp(1:igwx), ik) 
#else
         IF ( ionode ) &
            CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
#endif
         !
      END DO
#if defined __HDF5
      IF (ionode ) CALL h5fclose_f(evc_hdf5_write%file_id, ierr) 
#else 
      IF ( ionode ) CALL iotk_close_write( iuni )
#endif
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, ik, nk, ispin, &
                         nspin, wf, ngw, nbnd, igl, ngwl, filename, scalef, &
                         ionode, root_in_group, intra_group_comm)
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf
      USE mp,        ONLY : mp_put, mp_size, mp_rank, mp_sum, mp_max
      !
#if defined  __HDF5
      USE hdf5_qe
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      COMPLEX(DP),        INTENT(OUT)   :: wf(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      CHARACTER(LEN=256), INTENT(IN)    :: filename
      REAL(DP),           INTENT(OUT)   :: scalef
      LOGICAL,            INTENT(IN)    :: ionode
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      !
      INTEGER                  :: j
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER                  :: ierr
      INTEGER                  :: igwx, igwx_, ik_, nk_
      INTEGER                  :: me_in_group, nproc_in_group
      CHARACTER(LEN=256)       :: filename_hdf5
      !
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      !
      ierr = 0
      !
#if !defined __HDF5
      IF ( ionode ) CALL iotk_open_read( iuni, FILE = filename, &
                              BINARY = .TRUE., IERR = ierr )
      !
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      CALL errore( 'read_wfc ', &
                   'cannot open restart file for reading', ierr )
      !
#endif
      IF ( ionode ) THEN
          !
#if defined  __HDF5
          filename_hdf5=filename
          CALL prepare_for_reading_final(evc_hdf5_write,evc_hdf5_write%comm,filename_hdf5,ik)
          CALL read_attributes_hdf5(evc_hdf5_write,ngw,"ngw",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,nbnd,"nbnd",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,ik_,"ik",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,nk_,"ik",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,ispin,"ispin",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,nspin,"nspin",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,igwx_,"igwx",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,scalef,"scale_factor",ik)
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
      ALLOCATE( wtmp( MAX( igwx_, igwx ) ) )
      !
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wf, 2 ) ) THEN
            !
            IF ( ionode ) THEN 
#if defined __HDF5
             CALL read_evc(evc_hdf5_write,j,wtmp(1:igwx_),ik)
             !  CALL iotk_scan_dat( iuni, &
             !                      "evc" // iotk_index( j ), wtmp(1:igwx_) )
#else
               CALL iotk_scan_dat( iuni, &
                                   "evc" // iotk_index( j ), wtmp(1:igwx_) )
 
#endif

               !
               IF ( igwx > igwx_ ) wtmp((igwx_+1):igwx) = 0.0_DP
               !
            END IF
            !
            CALL splitwf( wf(:,j), wtmp, ngwl, igl, &
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
         END IF
         !
      END DO
      !
#if !defined __HDF5
      IF ( ionode ) CALL iotk_close_read( iuni )
#endif
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !        
  END MODULE io_base
