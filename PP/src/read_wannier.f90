!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!-----------------------------------------------------------------------
MODULE read_wannier
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_wannier_chk
  !
  ! Remember to broadcast possible additions
  !
  INTEGER, PUBLIC :: num_bands          ! num of PC bands for wannierization
  INTEGER, PUBLIC :: num_wann           ! num of PC Wannier functions
  INTEGER, PUBLIC :: num_kpts           ! num of k-points
  INTEGER, PUBLIC :: kgrid(3)           ! MP grid
  LOGICAL, ALLOCATABLE, PUBLIC :: excluded_band(:)
  LOGICAL, ALLOCATABLE, PUBLIC :: lwindow(:,:)  ! disentanglement parameters
  INTEGER, ALLOCATABLE, PUBLIC :: ndimwin(:)    ! disentanglement parameters
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: u_mat(:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: u_mat_opt(:,:,:)
  LOGICAL, PUBLIC :: have_disentangled
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE read_wannier_chk( seedname )
    !-------------------------------------------------------------------
    !
    ! ...  parser for the Wannier90 chk file
    !
    USE io_global,           ONLY : ionode, ionode_id
    USE mp_global,           ONLY : intra_image_comm
    USE mp,                  ONLY : mp_bcast
    USE cell_base,           ONLY : bg
    USE klist,               ONLY : nkstot
    USE wvfct,               ONLY : nbnd
    USE lsda_mod,            ONLY : nspin
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(IN) :: seedname
    !
    CHARACTER(LEN=33) :: header
    INTEGER :: i, j, nkp, nn
    INTEGER :: chk_unit=124
    INTEGER :: num_exclude_bands, nntot
    INTEGER, ALLOCATABLE :: exclude_bands(:)
    LOGICAL :: exst
    LOGICAL :: checkpoint
    REAL(DP) :: at_(3,3), bg_(3,3)
    REAL(DP), ALLOCATABLE :: kpt_latt(:,:)
    REAL(DP), ALLOCATABLE :: centers(:,:)
    REAL(DP), ALLOCATABLE :: spreads(:)
    REAL(DP) :: omega_invariant
    COMPLEX(DP), ALLOCATABLE :: m_mat(:,:,:,:)
    !
    !
    !
    IF ( ionode ) THEN
      !
      INQUIRE( FILE=trim(seedname)//'.chk', EXIST=exst )
      IF ( .not. exst ) CALL errore( 'read_wannier_chk', 'chk file not found', 1 )
      !
      OPEN( UNIT=chk_unit, FILE=trim(seedname)//'.chk', STATUS='old', FORM='unformatted' )
      !
      READ( chk_unit ) header                   ! date and time
      READ( chk_unit ) num_bands                ! number of bands
      READ( chk_unit ) num_exclude_bands        ! number of excluded bands
      !
      IF ( num_exclude_bands .lt. 0 .or. num_exclude_bands .ne. (nbnd - num_bands) ) &
        CALL  errore( 'read_wannier_chk', 'Invalid value for num_exclude_bands', &
                      num_exclude_bands )
      !
    ENDIF
    !
    CALL mp_bcast( num_bands, ionode_id, intra_image_comm )
    CALL mp_bcast( num_exclude_bands, ionode_id, intra_image_comm )
    !
    ALLOCATE( exclude_bands(num_exclude_bands) )
    ALLOCATE( excluded_band(nbnd) )
    !
    IF ( ionode ) THEN 
      !
      READ( chk_unit ) ( exclude_bands(i), i=1,num_exclude_bands )   ! list of excluded bands
      excluded_band(:) = .false.
      DO i = 1, num_exclude_bands
        excluded_band( exclude_bands(i) ) = .true.
      ENDDO
      !
      READ( chk_unit ) (( at_(i,j), i=1,3 ), j=1,3 )                 ! prim real latt vectors
      READ( chk_unit ) (( bg_(i,j), i=1,3 ), j=1,3 )                 ! prim recip latt vectors
      READ( chk_unit ) num_kpts                                      ! num of k-points
      !
      IF ( nspin == 2 ) THEN
        IF ( num_kpts .ne. nkstot/2 ) &
        CALL errore( 'read_wannier_chk', 'Invalid value for num_kpts', num_kpts )
      ELSE
        IF ( num_kpts .ne. nkstot ) &
        CALL errore( 'read_wannier_chk', 'Invalid value for num_kpts', num_kpts )
      ENDIF
      !
      READ( chk_unit ) ( kgrid(i), i=1,3 )                           ! MP grid
      !
    ENDIF
    !
    CALL mp_bcast( excluded_band, ionode_id, intra_image_comm )
    CALL mp_bcast( num_kpts, ionode_id, intra_image_comm )
    CALL mp_bcast( kgrid, ionode_id, intra_image_comm )
    !
    ALLOCATE( kpt_latt(3,num_kpts) )
    !
    IF ( ionode ) THEN
      !
      READ( chk_unit ) (( kpt_latt(i,nkp), i=1,3 ), nkp=1,num_kpts )
      READ( chk_unit ) nntot                                         ! nntot
      READ( chk_unit ) num_wann                                      ! num of WFs
      READ( chk_unit ) checkpoint                                    ! checkpoint
      READ( chk_unit ) have_disentangled     ! .true. if disentanglement has been performed
      !
    ENDIF
    !
    CALL mp_bcast( num_wann, ionode_id, intra_image_comm )
    CALL mp_bcast( have_disentangled, ionode_id, intra_image_comm )
    !
    IF ( have_disentangled ) THEN
      !
      ALLOCATE( lwindow(num_bands,num_kpts) )
      ALLOCATE( ndimwin(num_kpts) )
      ALLOCATE( u_mat_opt(num_bands,num_wann,num_kpts) )
      !
      IF ( ionode ) THEN
        !
        READ( chk_unit ) omega_invariant              
        READ( chk_unit ) (( lwindow(i,nkp), i=1,num_bands ), nkp=1,num_kpts )
        READ( chk_unit ) ( ndimwin(nkp), nkp=1,num_kpts )
        READ( chk_unit ) ((( u_mat_opt(i,j,nkp), i=1,num_bands ), &   ! optimal U-matrix
                                                  j=1,num_wann ), &
                                                  nkp=1,num_kpts )
        !
      ENDIF
      !
      CALL mp_bcast( lwindow, ionode_id, intra_image_comm )
      CALL mp_bcast( ndimwin, ionode_id, intra_image_comm )
      CALL mp_bcast( u_mat_opt, ionode_id, intra_image_comm )
      !
    ELSE
      !
      IF ( num_wann .ne. num_bands ) &
        CALL errore( 'read_wannier_chk', 'mismatch between num_bands and num_wann', &
                                                  num_bands-num_wann )
      !
    ENDIF
    !
    ALLOCATE( u_mat(num_wann,num_wann,num_kpts) )
    !
    IF ( ionode ) THEN
      !
      READ ( chk_unit ) ((( u_mat(i,j,nkp), i=1,num_wann ), &         ! U-matrix
                                            j=1,num_wann ), &
                                            nkp=1,num_kpts )
      !
      CALL check_u_unitary          ! checks u_mat is unitary
      !
      ALLOCATE( m_mat(num_wann,num_wann,nntot,num_kpts) )  
      ALLOCATE( centers(3,num_wann) )                                ! Wannier centers
      ALLOCATE( spreads(num_wann) )                                  ! Wannier spreads
      !
      READ( chk_unit ) (((( m_mat(i,j,nn,nkp), i=1,num_wann ), &
                                               j=1,num_wann ), &
                                               nn=1,nntot ), &
                                               nkp=1,num_kpts )
      !
      READ( chk_unit ) (( centers(i,j), i=1,3 ), j=1,num_wann )
      READ( chk_unit ) ( spreads(i), i=1,num_wann )
      !
      CLOSE( chk_unit )
      !
    ENDIF
    !
    CALL mp_bcast( u_mat, ionode_id, intra_image_comm )
    !
    !
  END SUBROUTINE read_wannier_chk
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE check_u_unitary
    !-------------------------------------------------------------------
    !
    ! ...  check that u_mat is unitary
    !
    USE constants,           ONLY : eps8
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: nkp, i, j
    COMPLEX(DP) :: uu_prod(num_wann,num_wann)
    !
    !  
    DO nkp = 1, num_kpts
      !
      uu_prod = MATMUL( u_mat(:,:,nkp), CONJG(TRANSPOSE( u_mat(:,:,nkp) )) )
      !
      DO i = 1, num_wann
        DO j = 1, num_wann
          !
          IF ( i == j ) THEN
            !
            IF ( ( ABS(DBLE( uu_prod(i,j) - 1 )) .ge. eps8 ) .or. &
                ( ABS(AIMAG( uu_prod(i,j) )) .ge. eps8 ) ) &
              CALL errore( 'read_wannier_chk', 'u_mat is not unitary', nkp )
            !
          ELSE
            !
            IF ( ( ABS(DBLE( uu_prod(i,j) )) .ge. eps8 ) .or. &
                ( ABS(AIMAG( uu_prod(i,j) )) .ge. eps8 ) ) &
              CALL errore( 'read_wannier_chk', 'u_mat is not unitary', nkp )
            !
          ENDIF
          !
        ENDDO
      ENDDO
      !
    ENDDO
    !
    !
  END SUBROUTINE check_u_unitary
  !
  !
END MODULE read_wannier
