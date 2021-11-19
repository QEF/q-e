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
!---------------------------------------------------------------------
MODULE cp_files
  !-------------------------------------------------------------------
  !
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: write_wannier_cp
  !
  CONTAINS
  !
  !-------------------------------------------------------------------
  SUBROUTINE write_wannier_cp( iun, nword, nwann, ks_only, typ )
    !-----------------------------------------------------------------
    !
    ! ...  This routine takes the Wannier functions in input and 
    ! ...  writes them into a file, readable by the CP-Koopmans code.
    !
    ! ...  split_evc_file = .false. : only one file called 'evcw.dat'
    ! ...                             containing the 2 spin channels
    ! ...                             is produced in output
    !
    ! ...  split_evc_file = .false. : two files called 'evcw1.dat' and 
    ! ...                             'evcw2.dat' containing respectively
    ! ...                             the spin up and spin down components
    ! ...                             are produced 
    !
    ! ...  For ks_only = .true. the KS states are directly written into the 
    ! ...  files 'evc_occupied.dat', 'evc0_empty1.dat' and 'evc0_empty2.dat'.
    !
    ! ...  NB: for the moment the spin down component is a copy of the spin up!
    !
    USE kinds,               ONLY : DP
    USE klist,               ONLY : nelec
    USE io_global,           ONLY : ionode, ionode_id
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE mp_wave,             ONLY : mergewf
    USE mp_world,            ONLY : mpime, nproc
    USE mp,                  ONLY : mp_sum
    USE noncollin_module,    ONLY : npol
    USE buffers,             ONLY : get_buffer
    USE wannier,             ONLY : split_evc_file
    USE read_wannier,        ONLY : num_kpts
    USE fft_supercell,       ONLY : npwxcp, ig_l2g_cp
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iun                   ! unit to the WFs buffer
    INTEGER, INTENT(IN) :: nword                 ! record length WF file
    INTEGER, INTENT(IN) :: nwann                 ! num of (primitive cell) WFs
    LOGICAL, INTENT(IN) :: ks_only
    CHARACTER(LEN=3), OPTIONAL, INTENT(IN) :: typ
    !
    CHARACTER(LEN=20) :: filename
    INTEGER :: io_level = 1
    INTEGER :: cp_unit = 125
    INTEGER :: npw_g                             ! global number of PWs
    INTEGER :: ir, ibnd, ibnd_, ispin, ipw
    INTEGER :: nwannx
    COMPLEX(DP), ALLOCATABLE :: evc(:,:)
    COMPLEX(DP), ALLOCATABLE :: evc_g(:)
    !
    !
    ALLOCATE( evc(npwxcp*npol,nwann) )
    !
    IF ( ks_only ) THEN
      !
      IF ( .not. PRESENT(typ) ) &
        CALL errore( 'write_wannier_cp', 'ks_only=.true. needs typ', 1 )
      !
      IF ( trim(typ) == 'occ' ) THEN
        !
        nwannx = nelec / 2 * num_kpts
        split_evc_file = .false.
        !
      ELSEIF ( trim(typ) == 'emp' ) THEN
        !
        nwannx = ( nwann - nelec / 2 ) * num_kpts
        split_evc_file = .true.
        !
      ELSE
        !
        CALL errore( 'write_wannier_cp', 'Wrong value for typ', 1 )
        !
      ENDIF
      !
    ELSE
      !
      nwannx = nwann * num_kpts
      !
    ENDIF
    !
    npw_g = npwxcp
    CALL mp_sum( npw_g, intra_bgrp_comm )
    ALLOCATE( evc_g(npw_g) )
    !
    IF ( .not. split_evc_file .and. ionode ) THEN
      !
      IF ( ks_only ) THEN
        filename = 'evc_occupied.dat'
      ELSE
        filename = 'evcw.dat'
      ENDIF
      !
      OPEN( UNIT=cp_unit, FILE=trim(filename), STATUS='unknown', FORM='unformatted' )
      WRITE( cp_unit ) npw_g, nwannx*2
      !
    ENDIF
    !
    ! ... here we gather the wfc from all the processes
    ! ... and we write it to file (nspin=2 in CP-Koopmans)
    !
    DO ispin = 1, 2
      !
      IF ( split_evc_file .and. ionode ) THEN
        !
        IF ( ks_only ) THEN
          WRITE( filename, 100 ) ispin
        ELSE
          WRITE( filename, 101 ) ispin
        ENDIF
        !
        OPEN( UNIT=cp_unit, FILE=trim(filename), STATUS='unknown', FORM='unformatted' )
        WRITE( cp_unit ) npw_g, nwannx
        !
      ENDIF
      !
      DO ir = 1, num_kpts
        !
        CALL get_buffer( evc, nword, iun, ir )
        !
        DO ibnd = 1, nwannx/num_kpts
          !
          IF ( ks_only .and. typ == 'emp' ) THEN
            ibnd_ = ibnd + nelec / 2
          ELSE
            ibnd_ = ibnd
          ENDIF
          !
          evc_g(:) = ( 0.D0, 0.D0 )
          CALL mergewf( evc(:,ibnd_), evc_g, npwxcp, ig_l2g_cp, mpime, &
                        nproc, ionode_id, intra_bgrp_comm )
          !
          IF ( ionode ) THEN
            !
            WRITE( cp_unit ) ( evc_g(ipw), ipw=1,npw_g )
            !
          ENDIF
          !
        ENDDO
        !
      ENDDO
      !
      IF ( split_evc_file .and. ionode ) CLOSE ( cp_unit )
      !
    ENDDO
    !
    IF ( .not. split_evc_file .and. ionode ) CLOSE( cp_unit )
    !
    !
100 FORMAT( 'evc0_empty', I1, '.dat' )
101 FORMAT( 'evcw', I1, '.dat' )
    !
    !
  END SUBROUTINE write_wannier_cp
  !
  !
END MODULE cp_files
