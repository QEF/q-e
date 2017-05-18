!
! Copyright (C) 2016-2017 Quantum ESPRESSO Foundation 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_base
  !----------------------------------------------------------------------------
  !
  ! ... subroutines used to read and write binary data produced by QE
  ! ... Author: Paolo Giannozzi, based on previous work by Carlo Cavazzoni
  !
  USE kinds,     ONLY : dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: write_wfc, read_wfc, write_rhog, read_rhog
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, filename, ik, nk, ispin, nspin, wfc, ngw,   &
                          gamma_only, nbnd, igl, ngwl, mill_k, scalef, &
                          ionode_in_group, root_in_group, intra_group_comm)
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf, mergekg
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
      CHARACTER(LEN=*),   INTENT(IN) :: filename
      INTEGER,            INTENT(IN) :: ik, nk, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wfc(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      INTEGER,            INTENT(IN) :: mill_k(:,:)
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      LOGICAL,            INTENT(IN) :: ionode_in_group
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm
      !
      INTEGER                  :: j, ierr
      INTEGER                  :: igwx, npwx, npol
      INTEGER                  :: me_in_group, nproc_in_group, my_group
      INTEGER, ALLOCATABLE     :: itmp(:,:)
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
         OPEN ( UNIT = iuni, FILE = TRIM(filename)//'.dat', &
              FORM='unformatted', STATUS = 'unknown' )
         !
         WRITE(iuni) ik, nk, ispin, nspin
         WRITE(iuni) gamma_only, scalef
         WRITE(iuni) ngw, igwx, nbnd
         !
#endif
         !
      END IF
      !
      ALLOCATE( itmp( 3, MAX (igwx,1) ) )
      itmp (:,:) = 0
      CALL mergekg( mill_k, itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      IF ( ionode_in_group ) WRITE(iuni) itmp(1:3,1:igwx)
      DEALLOCATE( itmp )
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
            WRITE(iuni) wtmp(1:npol*igwx)
#endif
         !
      END DO
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL h5fclose_f(h5_write_desc%file_id, ierr)
         DEALLOCATE ( h5_write_desc)  
#else 
         CLOSE (UNIT = iuni, STATUS = 'keep' )
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
    SUBROUTINE read_wfc( iuni, filename, ik, nk, ispin, nspin, wfc, ngw, nbnd, &
                         igl, ngwl, mill_k, scalef, &
                         ionode_in_group, root_in_group, intra_group_comm, &
                         ierr )
      ! if ierr is present, return 0 if everything is ok, /= 0 if not
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf, splitkg
      USE mp,        ONLY : mp_bcast, mp_size, mp_rank, mp_max
      !
#if defined  __HDF5
      USE hdf5_qe,  ONLY  : prepare_for_reading_final, read_attributes_hdf5, read_evc, &
                            h5fclose_f, hdf5_type
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      CHARACTER(LEN=*),   INTENT(IN)    :: filename
      COMPLEX(DP),        INTENT(OUT)   :: wfc(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      REAL(DP),           INTENT(OUT)   :: scalef
      INTEGER,            INTENT(OUT)   :: mill_k(:,:)
      LOGICAL,            INTENT(IN)    :: ionode_in_group
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      COMPLEX(DP), ALLOCATABLE          :: wtmp(:)
      INTEGER, ALLOCATABLE              :: itmp(:,:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx, igwx_, npwx, npol, ik_, nk_
      INTEGER                           :: me_in_group, nproc_in_group
      LOGICAL                           :: gamma_only_
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
      IF ( ionode_in_group ) OPEN ( UNIT = iuni, FILE=TRIM(filename)//'.dat', &
           FORM='unformatted', STATUS = 'old', IOSTAT = ierr_)
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
         IF ( PRESENT (ierr) ) THEN 
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
         READ (iuni) ik_, nk_, ispin, nspin
         READ(iuni) gamma_only_, scalef
         READ (iuni) ngw, igwx_, nbnd
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
      ALLOCATE( itmp( 3,MAX( igwx_, igwx ) ) )
      IF ( ionode_in_group ) THEN 
         READ(iuni) itmp(1:3,1:igwx_)
         IF ( igwx > igwx_ ) itmp(1:3,igwx_+1:igwx) = 0
      ELSE
         itmp (:,:) = 0
      END IF
      CALL splitkg( mill_k(:,:), itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      DEALLOCATE (itmp)
      !
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wfc, 2 ) ) THEN
            !
            IF ( ionode_in_group ) THEN 
#if defined __HDF5
               CALL read_evc(h5_read_desc, j, wtmp(1:npol*igwx_),ik)
#else
               READ (iuni) wtmp(1:npol*igwx_) 
#endif
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
      IF ( ionode_in_group ) CLOSE ( UNIT = iuni, STATUS = 'keep' )
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
    !------------------------------------------------------------------------
    SUBROUTINE write_rhog ( dirname, b1, b2, b3, gamma_only, mill, ig_l2g, &
         rho )
      !------------------------------------------------------------------------
      !! Write rho(G) in reciprocal space and related information to file
      !! 'charge-density.*' (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Quick-and-dirty version, allocates a large array on all mpi processes
      !
      USE mp,                   ONLY : mp_sum, mp_bcast
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE io_global,            ONLY : ionode, ionode_id
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !! directory name where file is written - must end by '/'
      REAL(dp),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !!  b1, b2, b3 are the three primitive vectors in a.u.
      INTEGER,          INTENT(IN) :: mill(:,:)
      !! Miller indices for local G-vectors
      !! G = mill(1)*b1 + mill(2)*b2 + mill(3)*b3
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      LOGICAL,          INTENT(IN) :: gamma_only
      !! if true, only the upper half of G-vectors (z >=0) is present
      COMPLEX(dp),      INTENT(IN) :: rho(:,:)
      !! rho(G) on this processor
      !
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      !! Global rho(G) for all processors
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !! Global Miller indices for all processors
      INTEGER                  :: ngm, nspin, ngm_g
      INTEGER                  :: iun, ns, ig, ierr
      CHARACTER(LEN=320)       :: filename
      !
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (mill, 2) .OR. ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('write_rhog', 'inconsistent input dimensions', 1)
      nspin= SIZE (rho, 2)
      iun  = 4
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      filename = TRIM( dirname ) // 'charge-density.dat'
      ierr = 0
      IF ( ionode ) OPEN ( UNIT = iun, FILE = TRIM( filename ), &
                FORM = 'unformatted', STATUS = 'unknown', iostat = ierr )
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error opening file ' &
           & // TRIM( filename ), 1 )
      IF ( ionode ) THEN
          WRITE (iun, iostat=ierr) gamma_only, ngm_g, nspin
          WRITE (iun, iostat=ierr) b1, b2, b3
      END IF
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 1 )
      !
      ! ... collect all G-vectors across processors within the pools
      !
      ALLOCATE( mill_g( 3, ngm_g ) )
      mill_g = 0
      DO ig = 1, ngm
         !
         ! ... ig is the local index for local G-vectors
         ! ... ig_l2g(ig) is the position of G-vector ig in the global list
         !
         mill_g(1,ig_l2g(ig)) = mill(1,ig)
         mill_g(2,ig_l2g(ig)) = mill(2,ig)
         mill_g(3,ig_l2g(ig)) = mill(3,ig)
         !
      END DO
      !
      CALL mp_sum( mill_g, intra_bgrp_comm )
      !
      ! ... write G-vectors
      !
      IF ( ionode ) WRITE (iun, iostat=ierr) mill_g(1:3,1:ngm_g)
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... deallocate to save memory
      !
      DEALLOCATE( mill_g )
      !
      ! ... now collect all G-vector components of the charge density
      ! ... (one spin at the time to save memory) using the same logic
      !
      ALLOCATE( rho_g(ngm_g ) )
      !
      DO ns = 1, nspin
         !
         rho_g = 0
         DO ig = 1, ngm
            rho_g(ig_l2g(ig)) = rho(ig,ns)
         END DO
         !
         CALL mp_sum( rho_g, intra_bgrp_comm )
         !
         IF ( ionode ) WRITE (iun,iostat=ierr) rho_g(1:ngm_g)
         CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
         IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
              & // TRIM( filename ), 2+ns )
         !
      END DO
      !
      IF (ionode) CLOSE (UNIT = iun, status ='keep' )
      !
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE write_rhog
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rhog ( dirname, ig_l2g, nspin, rho )
      !------------------------------------------------------------------------
      !! Read rho(G) in reciprocal space from file  'charge-density.*' 
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Quick-and-dirty version, allocates a large array on all mpi processes
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE io_global,            ONLY : ionode, ionode_id
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !! directory name where file is read - must end by '/'
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      INTEGER,          INTENT(IN) :: nspin
      !! read up to npsin components
      COMPLEX(dp),  INTENT(INOUT) :: rho(:,:)
      !
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      REAL(dp)                 :: b1(3), b2(3), b3(3)
      INTEGER                  :: ngm, nspin_, ngm_g
      INTEGER                  :: iun, mill_dum, ns, ig, ierr
      LOGICAL                  :: gamma_only
      CHARACTER(LEN=320)       :: filename
       !
      iun  = 4
      !
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('read_rhog', 'inconsistent input dimensions', 1)
      !
      filename = TRIM( dirname ) // 'charge-density.dat'
      ierr = 0
      IF ( ionode ) OPEN ( UNIT = iun, FILE = TRIM( filename ), &
                FORM = 'unformatted', STATUS = 'old', iostat = ierr )
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error opening file ' & !
           & // TRIM( filename ), 1 )
      IF ( ionode ) THEN
         READ (iun, iostat=ierr) gamma_only, ngm_g, nspin_
         READ (iun, iostat=ierr) b1, b2, b3
      END IF
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), 1 )
      CALL mp_bcast( ngm_g, ionode_id, intra_bgrp_comm )
      CALL mp_bcast( nspin_, ionode_id, intra_bgrp_comm )
      IF ( nspin > nspin_ ) &
         CALL errore('read_rhog', 'not enough spin components found', 1)
      !
      IF ( ngm_g < MAXVAL (ig_l2g(:)) ) &
           CALL errore('read_rhog', 'some G-vectors are missing', 1)
      !
      ! ... skip record containing G-vector indices
      !
      IF ( ionode ) READ (iun, iostat = ierr) mill_dum
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... now read, broadcast and re-order G-vector components
      ! ... of the charge density (one spin at the time to save memory)
      !
      ALLOCATE( rho_g(ngm_g) )
      !
      DO ns = 1, nspin
         !
         IF ( ionode ) READ (iun, iostat=ierr) rho_g(1:ngm_g)
         CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
         IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
              & // TRIM( filename ), 2+ns )
         CALL mp_bcast( rho_g, ionode_id, intra_bgrp_comm )
         !
         DO ig = 1, ngm
            rho(ig,ns) = rho_g(ig_l2g(ig))
         END DO
         !
      END DO
      !
      IF (ionode) CLOSE (UNIT = iun, status ='keep' )
      !
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE read_rhog
    !
  END MODULE io_base
