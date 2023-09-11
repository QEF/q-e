!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_kcw
  !----------------------------------------------------------------------------
  !
  !! this module contains some common subroutines used to read and write
  !! the data produced by KCW. 
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, iunpun
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  PRIVATE
  !
  LOGICAL,       SAVE      :: rho_binary = .TRUE.
  !
  PUBLIC :: rho_binary
  PUBLIC :: write_rhowann, read_rhowann, write_mlwf, read_mlwf
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rhowann( file_base, rho, fft_desc, ionode, inter_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho in real-space, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE mp,        ONLY : mp_get, mp_sum, mp_rank, mp_size
#if defined(__HDF5)
      USE qeh5_base_module,  ONLY  : qeh5_file, qeh5_dataset, qeh5_openfile, qeh5_open_dataset, &
                             qeh5_add_attribute, qeh5_write_dataset, qeh5_close, qeh5_set_space, &
                             qeh5_set_file_hyperslab              
#endif
      USE fft_types
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN) :: file_base
      COMPLEX(DP),       INTENT(IN) :: rho(:)
      TYPE(fft_type_descriptor),INTENT(IN) :: fft_desc
      INTEGER,           INTENT(IN) :: inter_group_comm
      LOGICAL,           INTENT(IN) :: ionode
      !
      INTEGER               :: nr1,nr2,nr3, nr1x, nr2x,nr3x
      INTEGER               :: rhounit, ierr, i, j, jj, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: rho_file
      CHARACTER(LEN=256)    :: rho_file_hdf5
      CHARACTER(LEN=10)     :: rho_extension
      COMPLEX(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: my_group_id, me_group, me_group2, me_group3, &
                                            nproc_group, nproc_group2, nproc_group3, &
                                            io_group_id, io_group2, io_group3
      INTEGER,  EXTERNAL    :: find_free_unit
      !
#if defined(__HDF5) 
      TYPE (qeh5_file)         :: h5file
      TYPE (qeh5_dataset)      :: rhowann_dset
      !  
#endif 
      !
      my_group_id = mp_rank( inter_group_comm )

      me_group = fft_desc%mype ; me_group2 = fft_desc%mype2 ; me_group3 = fft_desc%mype3
      nproc_group = fft_desc%nproc ; nproc_group2 = fft_desc%nproc2 ; nproc_group3 = fft_desc%nproc3
      !
      nr1  = fft_desc%nr1  ; nr2  = fft_desc%nr2  ; nr3  = fft_desc%nr3
      nr1x = fft_desc%nr1x ; nr2x = fft_desc%nr2x ; nr3x = fft_desc%nr3x
      !
      rho_extension = '.dat'
      IF ( .NOT. rho_binary ) rho_extension = '.xml'
      !
      rho_file = TRIM( file_base ) // TRIM( rho_extension )
      rhounit = find_free_unit ()
      !
      IF ( ionode ) THEN 
#if defined(__HDF5)
         CALL qeh5_openfile(h5file, TRIM(file_base)//'.hdf5',action = 'write') 
         CALL qeh5_add_attribute( h5file%id, "nr1", nr1 )
         CALL qeh5_add_attribute( h5file%id, "nr2", nr2 )
         CALL qeh5_add_attribute( h5file%id, "nr3", nr3 )
#else
         IF (rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='unformatted')
         IF (.NOT. rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='formatted')
         CALL errore( 'write_rhowann', 'cannot open ' // TRIM( rho_file ) // ' file for writing', ierr )
#endif
      END IF 
      !
#if !defined(__HDF5)
      IF ( ionode ) THEN
         !
         IF (rho_binary) THEN 
           WRITE(rhounit) nr1, nr2, nr3
         ELSE
           WRITE(rhounit, '("mesh ", 3i10)') nr1, nr2, nr3
         ENDIF
         !
      END IF
#endif
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      rho_plane = (0.D0, 0.D0)
      ALLOCATE( kowner( nr3 ) )
      !
#if defined(__HDF5)
      IF ( ionode ) THEN 
         CALL qeh5_set_space ( rhowann_dset, rho_plane(1), 2, [nr1*nr2, nr3], MODE = 'f')
         CALL qeh5_set_space ( rhowann_dset, rho_plane(1), 1, [nr1*nr2], MODE = 'm')
         CALL qeh5_open_dataset (h5file, rhowann_dset, ACTION = 'write', NAME = 'rhowann' )
      END IF
#endif 
      !
      ! ... find the index of the group (pool) that will write rho
      !
      io_group_id = 0
      !
      IF ( ionode ) io_group_id = my_group_id
      !
      CALL mp_sum( io_group_id, fft_desc%comm )
      CALL mp_sum( io_group_id, inter_group_comm ) ! io_group_id is the (pool) group that contains the ionode
      !
      ! ... find the index of the ionode within Y and Z  groups
      !
      io_group2 = 0 ; IF ( ionode ) io_group2 = me_group2
      CALL mp_sum( io_group2, fft_desc%comm )  ! io_group2 is the group index of the ionode in the Y group (nproc2)
      io_group3 = 0 ; IF ( ionode ) io_group3 = me_group3
      CALL mp_sum( io_group3, fft_desc%comm )  ! io_group3 is the group index of the ionode in the Z group (nproc3)
      !
      ! ... find out the owner of each "z" plane
      !
      DO ip = 1, nproc_group3
         !
         kowner( (fft_desc%i0r3p(ip)+1):(fft_desc%i0r3p(ip)+fft_desc%nr3p(ip)) ) = ip - 1
         !
      END DO
      !
      ldr = nr1x*fft_desc%my_nr2p
      !
      IF ( ( my_group_id == io_group_id ) ) THEN ! only the group of ionode collects and writes the data
         !
         DO k = 1, nr3
            !
            !  Only one subgroup write the charge density
            ! 
            rho_plane = (0.D0, 0.D0)
            IF( ( kowner(k) == me_group3 ) ) THEN
               !
               kk = k - fft_desc%my_i0r3p
               ! 
               DO jj = 1, fft_desc%my_nr2p
                  !
                  j = jj + fft_desc%my_i0r2p
                  DO i = 1, nr1
                     !
                     rho_plane(i+(j-1)*nr1) = rho(i+(jj-1)*nr1x+(kk-1)*ldr)
                     !
                  END DO
                  !
               END DO
               call mp_sum(rho_plane, fft_desc%comm2 ) ! collect the data over the Y group (nproc2)
               !
            END IF
            !
            ! if this processor is in the same comm3 group as ionode (me_group2==io_group2) 
            IF ( kowner(k) /= io_group3 .and. me_group2==io_group2) & 
               CALL mp_get( rho_plane, rho_plane, me_group3, io_group3, kowner(k), k, fft_desc%comm3 )
            !
            IF ( ionode ) THEN
#if defined(__HDF5)
              CALL qeh5_set_file_hyperslab ( rhowann_dset,  OFFSET = [0,k-1], COUNT = [2*nr1*nr2,1] ) 
              CALL qeh5_write_dataset ( rho_plane, rhowann_dset)   
#else
              !
              IF (rho_binary) THEN 
                 WRITE(rhounit) k
                 WRITE(rhounit) rho_plane
              ELSE
                 WRITE(rhounit, '("z= ", i5)') k
                 WRITE(rhounit, '(E23.16)') rho_plane
              ENDIF
              !
#endif
            ENDIF
            !
         END DO
         !
      END IF
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
#if defined(__HDF5)
         CALL qeh5_close (rhowann_dset) 
         CALL qeh5_close (h5file)   
#else
         CLOSE (rhounit, STATUS='keep')
#endif       
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rhowann
    !
   !------------------------------------------------------------------------
   SUBROUTINE read_rhowann( rho_file_base, fft_desc, rho )
     !------------------------------------------------------------------------
     !
     ! ... Reads charge density rho, one plane at a time, to avoid 
     ! ... collecting the entire charge density on a single processor
     !
     USE io_global, ONLY : ionode, ionode_id
     USE mp_images, ONLY : intra_image_comm
     USE mp,        ONLY : mp_put, mp_sum, mp_rank, mp_size
#if defined(__HDF5)
      USE  qeh5_base_module
#endif
     USE fft_types
     USE io_files,  ONLY : check_file_exist
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=*),  INTENT(IN)  :: rho_file_base
     TYPE(fft_type_descriptor),INTENT(IN) :: fft_desc
     COMPLEX(DP),          INTENT(OUT) :: rho(:)
     !
     INTEGER               :: rhounit, ierr, i, j, jj, k, kk, ldr, ip, k_
     INTEGER               :: nr( 3 ), nr1, nr2, nr3, nr1x, nr2x, nr3x
     INTEGER               :: me_group, me_group2, me_group3, &
                              nproc_group, nproc_group2, nproc_group3
     CHARACTER(LEN=256)    :: rho_file
     CHARACTER(LEN=256)    :: rho_file_hdf5
     COMPLEX(DP), ALLOCATABLE :: rho_plane(:)
     INTEGER,  ALLOCATABLE :: kowner(:)
     LOGICAL               :: exst
     INTEGER,  EXTERNAL    :: find_free_unit
     CHARACTER(LEN=256) :: string
     CHARACTER(LEN=10)     :: rho_extension
#if defined(__HDF5)
     INTEGER             ::   nr1_, nr2_, nr3_
     TYPE (qeh5_file)    ::   h5file
     TYPE (qeh5_dataset) ::   rhowann_dset
#endif  
     !
     me_group = fft_desc%mype ; me_group2 = fft_desc%mype2 ; me_group3 = fft_desc%mype3
     nproc_group = fft_desc%nproc ; nproc_group2 = fft_desc%nproc2 ; nproc_group3 = fft_desc%nproc3
     !
     nr1  = fft_desc%nr1  ; nr2  = fft_desc%nr2  ; nr3  = fft_desc%nr3
     nr1x = fft_desc%nr1x ; nr2x = fft_desc%nr2x ; nr3x = fft_desc%nr3x
     !
#if defined(__HDF5)
     rho_file_hdf5 = TRIM( rho_file_base ) // '.hdf5'
     exst = check_file_exist(TRIM(rho_file_hdf5))
     IF ( .NOT. exst ) CALL errore ('read_rhowann', 'searching for '// TRIM(rho_file_hdf5),10)
#else 
     rho_extension = '.dat'
     IF ( .NOT. rho_binary ) rho_extension = '.xml'
     rhounit = find_free_unit ( )
     rho_file = TRIM( rho_file_base ) // TRIM( rho_extension )
     exst = check_file_exist( TRIM(rho_file) ) 
     !
     IF ( .NOT. exst ) CALL errore('read_rhowann', 'searching for '//TRIM(rho_file), 10)
#endif
     !
     IF ( ionode ) THEN
#if defined(__HDF5)
        !ALLOCATE ( h5desc)
        !CALL prepare_for_reading_final(h5desc, 0 ,rho_file_hdf5)
        !CALL read_attributes_hdf5(h5desc, nr1_,"nr1")
        !CALL read_attributes_hdf5(h5desc, nr2_,"nr2")
        !CALL read_attributes_hdf5(h5desc, nr3_,"nr3")
        !nr = [nr1_,nr2_,nr3_]
        CALL qeh5_openfile( h5file, TRIM(rho_file_hdf5), ACTION = 'read', ERROR = ierr)
        CALL errore( 'read_rhowann', 'cannot open ' // TRIM( rho_file_hdf5 ) // ' file for reading', ierr )
        CALL qeh5_read_attribute (h5file%id, "nr1", nr1_)
        CALL qeh5_read_attribute (h5file%id, "nr2", nr2_)
        CALL qeh5_read_attribute (h5file%id, "nr3", nr3_)
        nr = [nr1_,nr2_,nr3_]
#else
        IF (rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='unformatted', STATUS='old')
        IF (.NOT. rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='formatted', STATUS='old')
        CALL errore( 'read_rhowann', 'cannot open ' // TRIM( rho_file ) // ' file for reading', ierr )
        IF (rho_binary) READ(rhounit) nr(1), nr(2), nr(3)
        IF (.NOT. rho_binary) READ(rhounit,*) string,nr(1), nr(2), nr(3)
#endif
     !
     END IF
     !
     CALL mp_bcast( nr, ionode_id, intra_image_comm )
     !
     IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
        CALL errore( 'read_rhowann', 'dimensions do not match', 1 )
     !
     !
     ALLOCATE( rho_plane( nr1*nr2 ) )
     ALLOCATE( kowner( nr3 ) )
     !
#if defined(__HDF5) 
     IF (ionode ) THEN
       CALL qeh5_open_dataset( h5file, rhowann_dset, ACTION = 'read', NAME = 'rhowann')
       CALL qeh5_set_space ( rhowann_dset, rho_plane(1), RANK = 1, DIMENSIONS = [nr1*nr2], MODE = 'm') 
     ENDIF
#endif
     !
     DO ip = 1, nproc_group3
        !
        kowner( (fft_desc%i0r3p(ip)+1):(fft_desc%i0r3p(ip)+fft_desc%nr3p(ip)) ) = ip - 1
        !
     END DO
     !
     ldr = nr1x*fft_desc%my_nr2p
     !
     ! ... explicit initialization to zero is needed because the physical
     ! ... dimensions of rho may exceed the true size of the FFT grid 
     !
     rho(:) = (0.0_DP, 0.0_DP)
     !
     DO k = 1, nr3
        !
        ! ... only ionode reads the charge planes
        !
        IF ( ionode ) THEN
#if defined(__HDF5)
           !CALL  read_rho_hdf5(h5desc , k,rho_plane)
           CALL qeh5_set_file_hyperslab (rhowann_dset, OFFSET = [0,k-1], COUNT = [2*nr1*nr2,1] )
           CALL qeh5_read_dataset (rho_plane, rhowann_dset )
#else
           IF (rho_binary) THEN 
              READ(rhounit) k_
              READ(rhounit) rho_plane
           ELSE
              READ(rhounit,*) string, k_
              READ(rhounit, '(E23.16)') rho_plane
           ENDIF
#endif
        ENDIF
        !
        ! ... planes are sent to the destination processor (all processors in this image)
        !
        CALL mp_bcast( rho_plane, ionode_id, intra_image_comm )
        !
        IF( kowner(k) == me_group3 ) THEN
           !
           kk = k - fft_desc%my_i0r3p
           DO jj = 1, fft_desc%my_nr2p
              j = jj + fft_desc%my_i0r2p
              DO i = 1, nr1
                 rho(i+(jj-1)*nr1x+(kk-1)*ldr) = rho_plane(i+(j-1)*nr1)
              END DO
           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( rho_plane )
     DEALLOCATE( kowner )
     !
     IF ( ionode ) THEN
        !
#if defined(__HDF5)
        CALL qeh5_close(rhowann_dset) 
        CALL qeh5_close(h5file)
        !CALL h5fclose_f(h5desc%file_id,ierr)
        !DEALLOCATE ( h5desc)
#else
        CLOSE (rhounit)
#endif    
     END IF
     !
     RETURN
     !
   END SUBROUTINE read_rhowann
   !
   !
   ! NsC: Adapted from pw_write_binaries inside PW/scr/pw_restart_new.f90
   SUBROUTINE write_mlwf( )
     !------------------------------------------------------------------------
     !
     USE mp,                   ONLY : mp_sum, mp_max
     USE io_base,              ONLY : write_wfc
     USE cell_base,            ONLY : tpiba, bg
     USE control_flags,        ONLY : gamma_only
     USE gvect,                ONLY : ig_l2g
     USE buffers,              ONLY : get_buffer
     USE klist,                ONLY : nks, nkstot, xk, ngk, igk_k
     USE gvect,                ONLY : mill
     USE wvfct,                ONLY : npwx, nbnd
     USE lsda_mod,             ONLY : nspin, isk, lsda, nspin
     USE mp_pools,             ONLY : intra_pool_comm
     USE mp_bands,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, &
                                      root_bgrp_id, my_bgrp_id
     USE clib_wrappers,        ONLY : f_mkdir_safe
     !
     USE control_kcw,          ONLY : tmp_dir_kcw
     USE control_kcw,          ONLY : num_wann, spin_component, evc0, iuwfc_wann
     !
     IMPLICIT NONE
     !
     INTEGER               :: ios, ig, ispin
     INTEGER               :: ik, ik_g, ike, iks, npw_g
     INTEGER, EXTERNAL     :: global_kpoint_index
     INTEGER,  ALLOCATABLE :: ngk_g(:), mill_k(:,:)
     INTEGER,  ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
     CHARACTER(LEN=256)    :: dirname
     CHARACTER(LEN=320)    :: filename
     CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
     !
     INTEGER :: ik_eff
     INTEGER :: lrwfc 
     !
     dirname = TRIM (tmp_dir_kcw) 
     !
     ! ... check that restart_dir exists on all processors that write
     ! ... wavefunctions; create one if restart_dir is not found. This
     ! ... is needed for k-point parallelization, in the case of non-parallel
     ! ... scratch file systems, that are not visible to all processors
     !
     IF ( my_bgrp_id == root_bgrp_id .AND. me_bgrp == root_bgrp ) THEN
        ios = f_mkdir_safe( TRIM(dirname) )
     END IF
     !
     ! ... write wavefunctions and k+G vectors
     !
     ! ... These are the global index of the first and last k-points in this pool
     iks = global_kpoint_index (nkstot, 1)
     ike = iks + nks - 1
     !
     ! ... ngk_g: global number of k+G vectors
     !
     ALLOCATE( ngk_g( nks ) )
     ngk_g(1:nks) = ngk(1:nks)
     CALL mp_sum( ngk_g, intra_bgrp_comm)
     !
     ! ... The igk_l2g array yields the correspondence between the
     ! ... local k+G index and the global G index
     !
     ALLOCATE ( igk_l2g( npwx ) )
     !
     ! ... the igk_l2g_kdip local-to-global map yields the correspondence
     ! ... between the global order of k+G and the local index for k+G.
     !
     ALLOCATE ( igk_l2g_kdip( npwx ) )
     !
     ALLOCATE ( mill_k( 3, npwx ) )
     !
     k_points_loop: DO ik = 1, nks
        !
        ! I deal with one component at the time in KCW
        IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
        !
        ! ik_g is the index of k-point ik in the global list
        !
        ik_g = ik + iks - 1
        !
        ! ... Compute the igk_l2g array from previously computed arrays
        ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
        !
        igk_l2g = 0
        DO ig = 1, ngk (ik)
           igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
        END DO
        !
        ! ... npw_g is the maximum G vector index among all processors
        !
        npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
        CALL mp_max( npw_g, intra_pool_comm )
        !
        igk_l2g_kdip = 0
        CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik), igk_l2g, &
                             igk_l2g_kdip )
        !
        ! ... mill_k(:,i) contains Miller indices for (k+G)_i
        !
        DO ig = 1, ngk (ik)
           mill_k(:,ig) = mill(:,igk_k(ig,ik))
        END DO
        !
        ! ... read wavefunctions - do not read if already in memory (nsk==1)
        !
        lrwfc = num_wann * npwx 
        ik_eff = ik-(spin_component-1)*nkstot/nspin
        CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
        !
        IF ( nspin == 2 ) THEN
           !
           ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
           !
           ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
           ispin = isk(ik)
           filename = TRIM(dirname) // 'wfc_wann_' // updw(ispin) // '_ik' //&
                & TRIM(int_to_char(ik_g))
           !
        ELSE
           !
           ispin = 1
           filename = TRIM(dirname) // 'wfc_wann_ik' // TRIM(int_to_char(ik_g))
           !
        ENDIF
        !
        ! ... Only the first band group of each pool writes
        ! ... No warranty it works for more than one band group
        !
        IF ( my_bgrp_id == root_bgrp_id ) CALL write_wfc( iunpun, &
             filename, root_bgrp, intra_bgrp_comm, ik_g, tpiba*xk(:,ik), &
             ispin, nspin, evc0, npw_g, gamma_only, num_wann, &
             igk_l2g_kdip(:), ngk(ik), tpiba*bg(:,1), tpiba*bg(:,2), &
             tpiba*bg(:,3), mill_k, 1.D0 )
        !
     END DO k_points_loop
     !
     DEALLOCATE ( mill_k )
     DEALLOCATE ( igk_l2g_kdip )
     DEALLOCATE ( igk_l2g )
     DEALLOCATE ( ngk_g )
     !
     RETURN
     !
   END SUBROUTINE write_mlwf
   !
   !
   SUBROUTINE gk_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip, igwk )
     !-----------------------------------------------------------------------
     !
     ! ... This subroutine maps local G+k index to the global G vector index
     ! ... the mapping is used to collect wavefunctions subsets distributed
     ! ... across processors.
     ! ... This map is used to obtained the G+k grids related to each kpt
     !
     USE mp_bands,             ONLY : intra_bgrp_comm
     USE mp,                   ONLY : mp_sum
     !
     IMPLICIT NONE
     !
     ! ... Here the dummy variables
     !
     INTEGER, INTENT(IN)  :: npw_g, ngk_g, ngk
     INTEGER, INTENT(IN)  :: igk_l2g(ngk)
     INTEGER, INTENT(OUT) :: igk_l2g_kdip(ngk)
     INTEGER, OPTIONAL, INTENT(OUT) :: igwk(ngk_g)
     !
     INTEGER, ALLOCATABLE :: igwk_(:), itmp(:), igwk_lup(:)
     INTEGER              :: ig, ig_, ngg
     !
     !
     ALLOCATE( itmp( npw_g ) )
     ALLOCATE( igwk_( ngk_g ) )
     !
     itmp(:)  = 0
     igwk_(:) = 0
     !
     DO ig = 1, ngk
        itmp(igk_l2g(ig)) = igk_l2g(ig)
     END DO
     !
     CALL mp_sum( itmp, intra_bgrp_comm )
     !
     ngg = 0
     DO ig = 1, npw_g
        !
        IF ( itmp(ig) == ig ) THEN
           !
           ngg = ngg + 1
           igwk_(ngg) = ig
           !
        END IF
        !
     END DO
     !
     IF ( ngg /= ngk_g ) &
        CALL errore( 'gk_l2gmap_kdip', 'unexpected dimension in ngg', 1 )
     !
     IF ( PRESENT( igwk ) ) THEN
        !
        igwk(1:ngk_g) = igwk_(1:ngk_g)
        !
     END IF
     !
     ALLOCATE( igwk_lup( npw_g ) )
     !
!$omp parallel private(ig_, ig)
!$omp workshare
     igwk_lup = 0
!$omp end workshare
!$omp do
     DO ig_ = 1, ngk_g
        igwk_lup(igwk_(ig_)) = ig_
     END DO
!$omp end do
!$omp do
     DO ig = 1, ngk
        igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
     END DO
!$omp end do
!$omp end parallel
     !
     DEALLOCATE( igwk_lup )
     !
     DEALLOCATE( itmp, igwk_ )
     !
     RETURN
     !
   END SUBROUTINE gk_l2gmap_kdip
   !
   !
   ! NsC: Adapted from read_collected_wfc inside PW/scr/pw_restart_new.f90
   SUBROUTINE read_mlwf ( dirname, ik, evc )
     !------------------------------------------------------------------------
     !
     ! ... reads from directory "dirname" (new file format) for k-point "ik"
     ! ... wavefunctions from collected format into distributed array "evc"
     !
     USE control_flags,        ONLY : gamma_only
     USE klist,                ONLY : nkstot, nks, ngk, igk_k
     USE wvfct,                ONLY : npwx, nbnd
     USE gvect,                ONLY : ig_l2g
     USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm
     USE mp_pools,             ONLY : me_pool, root_pool, &
                                      intra_pool_comm
     USE mp,                   ONLY : mp_sum, mp_max
     USE io_base,              ONLY : read_wfc
     USE lsda_mod,             ONLY : nspin, isk, nspin
     USE control_kcw,          ONLY : num_wann
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=*), INTENT(IN) :: dirname
     INTEGER, INTENT(IN) :: ik
     COMPLEX(dp), INTENT(OUT) :: evc(:,:)
     !
     CHARACTER(LEN=320)   :: filename, msg
     INTEGER              :: ik_g, ig
     INTEGER              :: npol_, nbnd_
     INTEGER              :: ike, iks, ngk_g, npw_g, ispin
     INTEGER, EXTERNAL    :: global_kpoint_index
     INTEGER, ALLOCATABLE :: mill_k(:,:)
     INTEGER, ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
     LOGICAL              :: ionode_k
     REAL(DP)             :: scalef, xk_(3), b1(3), b2(3), b3(3)
     CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
     !
     ! ... the root processor of each pool reads
     !
     ionode_k = (me_pool == root_pool)
     !
     iks = global_kpoint_index (nkstot, 1)
     ike = iks + nks - 1
     !
     ! ik_g: index of k-point ik in the global list
     !
     ik_g = ik + iks - 1
     !
     ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
     !
     ALLOCATE ( igk_l2g_kdip( npwx ) )
     !
     ! ... The igk_l2g array yields the correspondence between the
     ! ... local k+G index and the global G index - requires arrays
     ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
     !
     ALLOCATE ( igk_l2g( npwx ) )
     igk_l2g = 0
     DO ig = 1, ngk(ik)
        igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
     END DO
     !
     ! ... npw_g: the maximum G vector index among all processors
     ! ... ngk_g: global number of k+G vectors for all k points
     !
     npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
     CALL mp_max( npw_g, intra_pool_comm )
     ngk_g = ngk(ik)
     CALL mp_sum( ngk_g, intra_bgrp_comm)
     !
     ! ... now compute the igk_l2g_kdip local-to-global map
     !
     igk_l2g_kdip = 0
     CALL gk_l2gmap_kdip( npw_g, ngk_g, ngk(ik), igk_l2g, &
          igk_l2g_kdip )
     DEALLOCATE ( igk_l2g )
     !
     IF ( nspin == 2 ) THEN
        !
        ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
        !
        ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
        ispin = isk(ik)
        filename = TRIM(dirname) // 'wfc_wann_' // updw(ispin) // '_ik' // &
             & TRIM(int_to_char(ik_g))
        !
     ELSE
        !
        filename = TRIM(dirname) // 'wfc_wann_ik' // TRIM(int_to_char(ik_g))
        !
     ENDIF
     !
     ! ... Miller indices are read from file (but not used)
     !
     ALLOCATE( mill_k ( 3,npwx ) )
     !
     evc=(0.0_DP, 0.0_DP)
     !
     CALL read_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
          ik_g, xk_, ispin, npol_, evc, npw_g, gamma_only, nbnd_, &
          igk_l2g_kdip(:), ngk(ik), b1, b2, b3, mill_k, scalef )
     !
     DEALLOCATE ( mill_k )
     DEALLOCATE ( igk_l2g_kdip )
     !
     ! ... here one should check for consistency between what is read
     ! ... and what is expected
     !
     IF ( nbnd_ < num_wann ) THEN
        WRITE (msg,'("The number of bands for this run is",I6,", but only",&
             & I6," bands were read from file")')  nbnd, nbnd_  
        CALL errore ('pw_restart - read_collected_wfc', msg, 1 )
     END IF
     !
     RETURN
     !
   END SUBROUTINE read_mlwf

END MODULE io_kcw




