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
    SUBROUTINE write_wfc( iuni, filename, ik, xk, ispin, nspin, wfc, ngw,   &
                          gamma_only, nbnd, igl, ngwl, b1,b2,b3, mill_k,    &
                          scalef, ionode_in_group, root_in_group, intra_group_comm)
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
      INTEGER,            INTENT(IN) :: ik, ispin, nspin
      REAL(DP),           INTENT(IN) :: xk(:)
      COMPLEX(DP),        INTENT(IN) :: wfc(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      INTEGER,            INTENT(IN) :: mill_k(:,:)
      REAL(DP),           INTENT(IN) :: b1(3), b2(3), b3(3)    
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
      !
      IF ( ionode_in_group ) THEN
#if defined  __HDF5
         CALL prepare_for_writing_final ( h5_write_desc, 0, &
              TRIM(filename)//'.hdf5',ik, ADD_GROUP = .false.)
         CALL add_attributes_hdf5(h5_write_desc, ik,"ik",ik)
         CALL add_attributes_hdf5(h5_write_desc, xk,"xk",xk)
         CALL add_attributes_hdf5(h5_write_desc, ispin,"ispin",ik)
         CALL add_attributes_hdf5(h5_write_desc, gamma_only,"gamma_only",ik)
         CALL add_attributes_hdf5(h5_write_desc, scalef,"scale_factor",ik)
         CALL add_attributes_hdf5(h5_write_desc, ngw,"ngw",ik)
         CALL add_attributes_hdf5(h5_write_desc, igwx,"igwx",ik)
         CALL add_attributes_hdf5(h5_write_desc, npol,"npol",ik)
         CALL add_attributes_hdf5(h5_write_desc, nbnd,"nbnd",ik)
#else
         OPEN ( UNIT = iuni, FILE = TRIM(filename)//'.dat', &
              FORM='unformatted', STATUS = 'unknown' )
         WRITE(iuni) ik, xk, ispin, gamma_only, scalef
         WRITE(iuni) ngw, igwx, npol, nbnd
#endif
         !
      END IF
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( itmp( 3, MAX (igwx,1) ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( itmp( 3, 1 ) )
      ENDIF
      itmp (:,:) = 0
      CALL mergekg( mill_k, itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL errore('write_wfc', 'hdf5 not yet ready',1)
#else
         WRITE(iuni) b1, b2, b3
         WRITE(iuni) itmp(1:3,1:igwx)
#endif
      END IF
      DEALLOCATE( itmp )
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( wtmp( MAX( npol*igwx, 1 ) ) )
      ELSE
         ALLOCATE( wtmp( 1 ) )
      ENDIF
      wtmp = 0.0_DP
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
         IF ( ionode_in_group ) THEN
#if defined(__HDF5)
            CALL write_evc(h5_write_desc, j, wtmp(1:npol*igwx), ik) 
#else
            WRITE(iuni) wtmp(1:npol*igwx)
#endif
         END IF
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
    SUBROUTINE read_wfc( iuni, filename, ik, xk, ispin, npol, wfc, ngw, &
                         gamma_only, nbnd, igl, ngwl, b1, b2, b3, mill_k,&
                         scalef, ionode_in_group, root_in_group, &
                         intra_group_comm, ierr )
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
      INTEGER,            INTENT(IN)    :: ik
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, npol
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      REAL(DP),           INTENT(OUT)   :: scalef
      REAL(DP),           INTENT(OUT)   :: xk(3)
      REAL(DP),           INTENT(OUT)   :: b1(3), b2(3), b3(3)
      INTEGER,            INTENT(OUT)   :: mill_k(:,:)
      LOGICAL,            INTENT(OUT)   :: gamma_only
      LOGICAL,            INTENT(IN)    :: ionode_in_group
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      COMPLEX(DP), ALLOCATABLE          :: wtmp(:)
      INTEGER, ALLOCATABLE              :: itmp(:,:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx, igwx_, npwx, ik_
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
         CALL read_attributes_hdf5(h5_read_desc, ik_,"ik",ik)
         CALL read_attributes_hdf5(h5_read_desc, xk,"xk",ik)
         CALL read_attributes_hdf5(h5_read_desc, ispin,"ispin",ik)
         CALL read_attributes_hdf5(h5_read_desc, gamma_only,"gamma_only",ik)
         CALL read_attributes_hdf5(h5_read_desc, scalef,"scale_factor",ik)
         CALL read_attributes_hdf5(h5_read_desc, ngw,"ngw",ik)
         CALL read_attributes_hdf5(h5_read_desc, nbnd,"nbnd",ik)
         CALL read_attributes_hdf5(h5_read_desc, npol,"npol",ik)
         CALL read_attributes_hdf5(h5_read_desc, igwx_,"igwx",ik)
#else
         READ (iuni) ik_, xk, ispin, gamma_only, scalef
         READ (iuni) ngw, igwx_, npol, nbnd
#endif
      END IF
      !
      CALL mp_bcast( ik_,    root_in_group, intra_group_comm )
      CALL mp_bcast( xk,     root_in_group, intra_group_comm )
      CALL mp_bcast( ispin,  root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only, root_in_group, intra_group_comm )
      CALL mp_bcast( scalef, root_in_group, intra_group_comm )
      CALL mp_bcast( ngw,    root_in_group, intra_group_comm )
      CALL mp_bcast( igwx_,  root_in_group, intra_group_comm )
      CALL mp_bcast( npol,   root_in_group, intra_group_comm )
      CALL mp_bcast( nbnd,   root_in_group, intra_group_comm )
      !
      npwx = SIZE( wfc, 1 ) / npol
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( itmp( 3,MAX( igwx_, igwx ) ) )
#if defined(__HDF5)
         CALL errore('read_wfc', 'hdf5 not yet ready',1)
#else
         READ (iuni) b1, b2, b3
         READ (iuni) itmp(1:3,1:igwx_)
#endif
         IF ( igwx > igwx_ ) itmp(1:3,igwx_+1:igwx) = 0
      ELSE
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      CALL splitkg( mill_k(:,:), itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      DEALLOCATE (itmp)
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( wtmp( npol*MAX( igwx_, igwx ) ) )
      ELSE
         ALLOCATE( wtmp(1) )
      ENDIF
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
      IF ( ionode_in_group ) THEN
#if defined (__HDF5)
         CALL h5fclose_f(h5_read_desc%file_id, ierr_)
         DEALLOCATE (h5_read_desc)
#else
         CLOSE ( UNIT = iuni, STATUS = 'keep' )
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
      !! Processor "ionode" collects data from band group, writes to file
      !
      USE mp,                   ONLY : mp_sum, mp_bcast
      USE mp_bands,             ONLY : me_bgrp, root_bgrp, nproc_bgrp, &
                                       intra_bgrp_comm
      USE mp_wave,              ONLY : mergewf, mergekg
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
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      !! Local rho(G), with LSDA workaround
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      !! Global rho(G) collected on root proc
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !! Global Miller indices collected on root proc
      INTEGER                  :: ngm, nspin, ngm_g, igwx
      INTEGER                  :: iun, ns, ig, ierr
      CHARACTER(LEN=320)       :: filename
      !
#if defined __HDF5
      CALL errore('write_rhog', 'hdf5 not yet ready',1)
#endif
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
      ! ... collect all G-vectors across processors within the band group
      !
      IF ( me_bgrp == root_bgrp ) THEN
         ALLOCATE( mill_g( 3, ngm_g ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( mill_g( 3, 1 ) )
      END IF
      !
      ! ... mergekg collects distributed array mill(1:3,ig) where ig is the
      ! ... local index, into array mill_g(1:3,ig_g), where ig_g=ig_l2g(ig)
      ! ... is the global index. mill_g is collected on root_bgrp only
      !
      CALL mergekg( mill, mill_g, ngm, ig_l2g, me_bgrp, &
           nproc_bgrp, root_bgrp, intra_bgrp_comm )
      !
      ! ... write G-vectors
      !
      IF ( ionode ) THEN
#if defined(__HDF5)
         CALL errore('write_rhog', 'hdf5 not yet ready',1)
#else
         WRITE (iun, iostat=ierr) mill_g(1:3,1:ngm_g)
#endif
      END IF
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
      IF ( me_bgrp == root_bgrp ) THEN
         ALLOCATE( rho_g( ngm_g ) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
      !
      DO ns = 1, nspin
         !
         ! Workaround for LSDA, while waiting for much-needed harmonization:
         ! we have rhoup and rhodw, we write rhotot=up+dw and rhodif=up-dw
         ! 
         IF ( ns == 1 .AND. nspin == 2 ) THEN
            DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns) + rho(ig,ns+1)
            END DO
         ELSE IF ( ns == 2 .AND. nspin == 2 ) THEN
            DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns-1) - rho(ig,ns)
            END DO
        ELSE
            DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns)
            END DO
         END IF
         !
         rho_g = 0
         CALL mergewf( rhoaux, rho_g, ngm, ig_l2g, me_bgrp, &
              nproc_bgrp, root_bgrp, intra_bgrp_comm )
         !
         IF ( ionode ) THEN
#if defined(__HDF5)
            CALL errore('write_rhog', 'hdf5 not yet ready',2)
#else
            WRITE (iun, iostat=ierr) rho_g(1:ngm_g)
#endif
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
         IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
              & // TRIM( filename ), 2+ns )
         !
      END DO
      !
      IF (ionode) CLOSE (UNIT = iun, status ='keep' )
      !
      DEALLOCATE( rhoaux )
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
      !! All pools read the file: 1 proc reads, broadcasts to other procs
      !! Works only if there is a single band group per pool
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_wave,              ONLY : splitwf
      USE mp_pools,             ONLY : me_pool, root_pool, nproc_pool, &
                                       intra_pool_comm
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !! directory name where file is read - must end by '/'
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      INTEGER,          INTENT(IN) :: nspin
      !! read up to nspin components
      COMPLEX(dp),  INTENT(INOUT) :: rho(:,:)
      !
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      COMPLEX(dp)              :: rhoup, rhodw
      REAL(dp)                 :: b1(3), b2(3), b3(3)
      INTEGER                  :: ngm, nspin_, ngm_g, isup, isdw
      INTEGER                  :: iun, mill_dum, ns, ig, ierr
      LOGICAL                  :: ionode_k, gamma_only
      CHARACTER(LEN=320)       :: filename
      !
#if defined __HDF5
      CALL errore('read_rhog', 'hdf5 not yet ready',1)
#endif
      !
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('read_rhog', 'inconsistent input dimensions', 1)
      !
      iun  = 4
      filename = TRIM( dirname ) // 'charge-density.dat'
      ierr = 0
      !
      ! ... the root processor of each pool reads
      !
      ionode_k = (me_pool == root_pool)
      !
      IF ( ionode_k ) THEN
         OPEN ( UNIT = iun, FILE = TRIM( filename ), &
              FORM = 'unformatted', STATUS = 'old', iostat = ierr )
         IF ( ierr /= 0 ) THEN
            ierr = 1
            GO TO 10
         END IF
         READ (iun, iostat=ierr) gamma_only, ngm_g, nspin_
         IF ( ierr /= 0 ) THEN
            ierr = 2
            GO TO 10
         END IF
         READ (iun, iostat=ierr) b1, b2, b3
         IF ( ierr /= 0 ) ierr = 3
10       CONTINUE
      END IF
      !
      CALL mp_bcast( ierr, root_pool, intra_pool_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), ierr )
      CALL mp_bcast( ngm_g, root_pool, intra_pool_comm )
      CALL mp_bcast( nspin_, root_pool, intra_pool_comm )
      !
      IF ( nspin > nspin_ ) &
         CALL infomsg('read_rhog', 'some spin components not found')
      IF ( ngm_g < MAXVAL (ig_l2g(:)) ) &
           CALL errore('read_rhog', 'some G-vectors are missing', 1)
      !
      ! ... skip record containing G-vector indices
      !
      IF ( ionode_k ) THEN
#if defined(__HDF5)
         CALL errore('write_rhog', 'hdf5 not yet ready',2)
#else
         READ (iun, iostat=ierr) mill_dum
#endif
      END IF
      CALL mp_bcast( ierr, root_pool, intra_pool_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... now read, broadcast and re-order G-vector components
      ! ... of the charge density (one spin at the time to save memory)
      !
      IF ( ionode_k ) THEN
         ALLOCATE( rho_g( ngm_g ) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
      !
      DO ns = 1, nspin
         !
         IF ( ionode_k ) THEN
#if defined(__HDF5)
            CALL errore('write_rhog', 'hdf5 not yet ready',2)
#else
            READ (iun, iostat=ierr) rho_g(1:ngm_g)
#endif
         END IF
         CALL mp_bcast( ierr, root_pool, intra_pool_comm )
         IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
              & // TRIM( filename ), 2+ns )
         !
         CALL splitwf( rhoaux, rho_g, ngm, ig_l2g, me_pool, &
              nproc_pool, root_pool, intra_pool_comm )
         DO ig = 1, ngm
            rho(ig,ns) = rhoaux(ig)
         END DO
         !
         ! Workaround for LSDA, while waiting for much-needed harmonization:
         ! if file contains rhotot=up+dw and rhodif=up-dw (nspin_=2), and
         ! if we want rhoup and rho down (nspin=2), convert 
         ! 
         IF ( nspin_ == 2 .AND. nspin == 2 .AND. ns == 2 ) THEN
            DO ig = 1, ngm
               rhoup = (rho(ig,ns-1) + rhoaux(ig)) / 2.0_dp
               rhodw = (rho(ig,ns-1) - rhoaux(ig)) / 2.0_dp
               rho(ig,ns-1)= rhoup
               rho(ig,ns  )= rhodw
            END DO
         END IF
      END DO
      !
      IF ( ionode_k ) CLOSE (UNIT = iun, status ='keep' )
      !
      DEALLOCATE( rhoaux )
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE read_rhog
    !
  END MODULE io_base
