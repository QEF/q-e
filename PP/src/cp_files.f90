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
  SUBROUTINE write_wannier_cp( iun, nword, norb, ks_only, typ )
    !-----------------------------------------------------------------
    !
    ! ...  This routine takes the Wannier/KS functions in input and 
    ! ...  writes them into file(2), readable by the CP-Koopmans code.
    ! ...  The output files are structured as below:
    ! ...
    ! ...  - ks_only=.false. and nspin=2
    ! ...    the Wannier functions of the calculated spin component
    ! ...    are written to a file called evcw.dat
    ! ...
    ! ...  - ks_only=.false. and nspin=1
    ! ...    the Wannier functions of the only spin component present
    ! ...    are written twice into two distinct files called evcw1.dat
    ! ...    and evcw2.dat
    ! ...
    ! ...  - ks_only=.true. and nspin=2
    ! ...    the KS wave functions of each spin component are written
    ! ...    to two files: evc_occupied1.dat and evc_occupied2.dat (for
    ! ...    occupied states), evc0_empty1.dat and evc0_empty2.dat (for 
    ! ...    empty states)
    ! ...
    ! ...  - ks_only=.true. and nspin=1
    ! ...    the KS wave functions of the spin component are written
    ! ...    twice into two distinct files: evc_occupied1.dat and 
    ! ...    evc_occupied2.dat (for occupied states), evc0_empty1.dat and
    ! ...    evc0_empty2.dat (for empty states)
    !
    USE kinds,               ONLY : DP
    USE klist,               ONLY : nelec
    USE io_global,           ONLY : ionode, ionode_id
    USE lsda_mod,            ONLY : nspin
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE mp_wave,             ONLY : mergewf
    USE mp_world,            ONLY : mpime, nproc
    USE mp,                  ONLY : mp_sum
    USE noncollin_module,    ONLY : npol
    USE buffers,             ONLY : get_buffer
    USE read_wannier,        ONLY : num_kpts
    USE fft_supercell,       ONLY : npwxcp, ig_l2g_cp
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iun                     ! unit to the WFs buffer
    INTEGER, INTENT(IN) :: nword                   ! record length WF file
    INTEGER, INTENT(IN) :: norb                    ! num of (primitive cell) WFs
    LOGICAL, INTENT(IN) :: ks_only
    CHARACTER(LEN=3), OPTIONAL, INTENT(IN) :: typ  ! required when ks_only=.true.
    !
    CHARACTER(LEN=20) :: filename
    INTEGER :: io_level = 1
    INTEGER :: cp_unit = 125
    INTEGER :: num_files
    INTEGER :: npw_g                             ! global number of PWs
    INTEGER :: ir, ir_, ibnd, ibnd_, ifile, ipw
    INTEGER :: norbx, nrtot
    COMPLEX(DP), ALLOCATABLE :: evc(:,:)
    COMPLEX(DP), ALLOCATABLE :: evc_g(:)
    !
    !
    ALLOCATE( evc(npwxcp*npol,norb) )
    !
    IF ( ks_only ) THEN
      !
      IF ( .not. PRESENT(typ) ) &
        CALL errore( 'write_wannier_cp', 'ks_only=.true. needs typ', 1 )
      !
      nrtot = num_kpts / nspin
      !
      IF ( typ == 'occ' ) THEN
        !
        norbx = nelec / 2 * nrtot
        !
      ELSEIF ( typ == 'emp' ) THEN
        !
        norbx = ( norb - nelec / 2 ) * nrtot
        !
      ELSE
        !
        CALL errore( 'write_wannier_cp', 'Wrong value for typ', 1 )
        !
      ENDIF
      !
    ELSE
      !
      nrtot = num_kpts
      norbx = norb * nrtot
      !
    ENDIF
    !
    npw_g = npwxcp
    CALL mp_sum( npw_g, intra_bgrp_comm )
    ALLOCATE( evc_g(npw_g) )
    !
    ! ... the only case where we want to write one file only is when
    ! ... we are calculating Wannier functions in the spin polarized
    ! ... case, since each spin channel is calculated independently.
    ! ... In all the other cases we want two files:
    ! ...
    ! ...   - evcw1.dat and evcw2.dat, which are one the copy of the
    ! ...     other, when calculating Wannier functions (ks_only=.false.)
    ! ...
    ! ...   - evc_XXX1.dat and evc_XXX2.dat (where XXX is 'occupied' 
    ! ...     or empty) when calculating KS states (ks_only=.true.)
    ! ...     where evc_XXX2.dat is whether a copy of evc_XXX1.dat
    ! ...     or a file written independently depending on the value
    ! ...     of nspin
    !
    IF ( nspin == 1 .or. ks_only ) THEN
      num_files = 2
    ELSE
      num_files = 1
    ENDIF
    !
    ! ... here we gather the wfc from all the processes
    ! ... and we write it to file(s) (in CP Koopmans nspin=2 always!).
    !
    DO ifile = 1, num_files
      !
      IF ( ks_only ) THEN
        !
        IF ( typ == 'occ' ) THEN
          WRITE( filename, 101 ) ifile
        ELSE
          WRITE( filename, 102 ) ifile
        ENDIF
        !
      ELSE
        !
        IF ( num_files == 1 ) THEN
          filename = 'evcw.dat'
        ELSE
          WRITE( filename, 100 ) ifile
        ENDIF
        !
      ENDIF
      !
      IF ( ionode ) THEN
        OPEN( UNIT=cp_unit, FILE=filename, STATUS='unknown', FORM='unformatted' )
        WRITE( cp_unit ) npw_g, norbx
      ENDIF
      !
      !
      DO ir = 1, nrtot
        !
        IF ( ks_only .and. nspin == 2 ) THEN
          ir_ = ir + (ifile-1)*nrtot
        ELSE
          ir_ = ir
        ENDIF
        !
        CALL get_buffer( evc, nword, iun, ir_ )
        !
        DO ibnd = 1, norbx/nrtot
          !
          ibnd_ = ibnd
          IF ( ks_only ) THEN
            IF ( typ == 'emp' ) ibnd_ = ibnd + nelec / 2
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
      IF ( ionode ) CLOSE ( cp_unit )
      !
    ENDDO
    !
    !
100 FORMAT( 'evcw', I1, '.dat' )
101 FORMAT( 'evc_occupied', I1, '.dat' )
102 FORMAT( 'evc0_empty', I1, '.dat' )
    !
    !
  END SUBROUTINE write_wannier_cp
  !
  !
END MODULE cp_files
