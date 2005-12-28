!
! Copyright (C) 2002-2005 quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE io_base
!=----------------------------------------------------------------------------=!

! Written by Carlo Cavazzoni - December 2002

! This module contains low level subroutine to write data to the restart file.
! The layout of the restart file is as follow:
!
! -----------------------------------------------------------------------------!
! Data Sections            !  Write Subroutine          Read Subroutine        !
! -----------------------------------------------------------------------------!
! HEADER                   !  write_restart_header      read_restart_header    !
! XDIM                     !  write_restart_xdim        read_restart_xdim      !
! CELL                     !  write_restart_cell        read_restart_cell      !
! IONS                     !  write_restart_ions        read_restart_ions      !
! LDA+U                    !  write_restart_ldaU        read_restart_ldaU
! SYMMETRIES               !  write_restart_symmetry    read_restart_symmetry  !
! do i = 1, ntyp           !                                                   !
!   PSEUDOPOTENTIALS( i )  !  write_restart_pseudo      read_restart_pseudo    !
! end do                   !                                                   !
! do j = 1, nspin          !                                                   !
!   do i = 1, nk           !                                                   !
!     OCCUPATIONS( i, j )  !  write_restart_electrons   read_restart_electrons !
!   end do                 !                                                   !
! end do                   !                                                   !
! G-VECTORS                !  write_restart_gvec        read_restart_gvec      !
! do i = 1, nk             !                                                   !
!   (G+k)-VECTORS( i )     !  write_restart_gkvec       read_restart_gkvec     !
! end do                   !                                                   !
! TETRAHEDRA               !  write_restart_tetra       read_restart_tetra     !
! DENSITY_AND_POTENTIAL    !  write_restart_charge      read_restart_charge    !
! do i = 1, nk             !                                                   !
!   WAVEFUNCTIONS( i )     !  write_restart_wfc         read_restart_wfc       !
! end do                   !                                                   !
! -----------------------------------------------------------------------------!
!
! All Data Sections are independent from each other, arrays are always stored
! togeter with their dimensions. Write and Read dummy subroutines should be
! used to skip Data Sections both writing and reading. The dummy subroutines
! have the same name as the effective subroutines (overloading), and are
! accessed through f90 interfaces.
!
! The following arguments are common to all subroutines:
! - iuni   - integer - the I/O fortran unit associated with the restart file
! - twrite - logical - true, write effective data; false, write dummy data
!
! All Data Sections have a well defined number of records, and the first
! record always contains the following _NON_ _DUMMY_ information:
! "twrite" "file_version" 
!
!

  USE io_global,  ONLY : stdout
  USE kinds
  USE parameters, ONLY: nsx

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: file_version = 299
  INTEGER :: restart_module_verbosity = 0

  INTERFACE write_restart_header
    MODULE PROCEDURE write_restart_header1, write_restart_header2
  END INTERFACE
  INTERFACE read_restart_header
    MODULE PROCEDURE read_restart_header1, read_restart_header2
  END INTERFACE

  INTERFACE write_restart_ions
    MODULE PROCEDURE write_restart_ions1, write_restart_ions2
  END INTERFACE
  INTERFACE read_restart_ions
    MODULE PROCEDURE read_restart_ions1, read_restart_ions2
  END INTERFACE

  INTERFACE write_restart_cell
    MODULE PROCEDURE write_restart_cell1, write_restart_cell2
  END INTERFACE
  INTERFACE read_restart_cell
    MODULE PROCEDURE read_restart_cell1, read_restart_cell2
  END INTERFACE

  INTERFACE write_restart_electrons
    MODULE PROCEDURE write_restart_electrons1, write_restart_electrons2
  END INTERFACE
  INTERFACE read_restart_electrons
    MODULE PROCEDURE read_restart_electrons1, read_restart_electrons2
  END INTERFACE

  INTERFACE write_restart_ldaU
    MODULE PROCEDURE write_restart_ldaU1, write_restart_ldaU2
  END INTERFACE
  INTERFACE read_restart_ldaU
    MODULE PROCEDURE read_restart_ldaU1, read_restart_ldaU2
  END INTERFACE

  INTERFACE write_restart_symmetry
    MODULE PROCEDURE write_restart_symmetry1, write_restart_symmetry2
  END INTERFACE
  INTERFACE read_restart_symmetry
    MODULE PROCEDURE read_restart_symmetry1, read_restart_symmetry2
  END INTERFACE

  INTERFACE write_restart_pseudo
    MODULE PROCEDURE write_restart_pseudo1, write_restart_pseudo2, write_restart_pseudo3
  END INTERFACE
  INTERFACE read_restart_pseudo
    MODULE PROCEDURE read_restart_pseudo1, read_restart_pseudo2, read_restart_pseudo3
  END INTERFACE

  INTERFACE write_restart_xdim
    MODULE PROCEDURE write_restart_xdim1, write_restart_xdim2
  END INTERFACE
  INTERFACE read_restart_xdim
    MODULE PROCEDURE read_restart_xdim1, read_restart_xdim2
  END INTERFACE

  INTERFACE write_restart_gvec
    MODULE PROCEDURE write_restart_gvec1, write_restart_gvec2
  END INTERFACE
  INTERFACE read_restart_gvec
    MODULE PROCEDURE read_restart_gvec1, read_restart_gvec2
  END INTERFACE

  INTERFACE write_restart_gkvec
    MODULE PROCEDURE write_restart_gkvec1, write_restart_gkvec2
  END INTERFACE
  INTERFACE read_restart_gkvec
    MODULE PROCEDURE read_restart_gkvec1, read_restart_gkvec2
  END INTERFACE

  INTERFACE write_restart_tetra
    MODULE PROCEDURE write_restart_tetra1, write_restart_tetra2
  END INTERFACE
  INTERFACE read_restart_tetra
    MODULE PROCEDURE read_restart_tetra1, read_restart_tetra2
  END INTERFACE

  INTERFACE write_restart_charge
    MODULE PROCEDURE write_restart_charge1, write_restart_charge2
  END INTERFACE
  INTERFACE read_restart_charge
    MODULE PROCEDURE read_restart_charge1, read_restart_charge2
  END INTERFACE

  INTERFACE write_restart_wfc
    MODULE PROCEDURE write_restart_wfc1, write_restart_wfc2
  END INTERFACE
  INTERFACE read_restart_wfc
    MODULE PROCEDURE read_restart_wfc1, read_restart_wfc2
  END INTERFACE


!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk dimensions, and status variables
!
! .. Where:
! iuni          = Restart file I/O fortran unit
! nfi           = Step counter
! trutim        = true time (in a.u.) since last 'from_scratch'
! nr1, nr2, nr3 = dims of the real space grid
! ng            = number of reciprocal vectors
! nk            = number of k points
! ngwk(.)       = number of wave functions G-vectors for k points
! nspin         = number of spin
! nbnd          = number of electronic bands
! nel           = total number of electrons
! nelu          = number of electrons with spin up
! neld          = number of electrons with spin down
! nat           = number of atoms
! ntyp          = number of atomic species
! na(.)         = number of atom for each species
! acc(.)        = accomulators
! nacc          = number of accumulators
! ecutwfc       = wave function cutoff
! ecutrho       = charge density cutoff cutoff

    SUBROUTINE write_restart_header1(iuni, &
      nfi, trutim, nr1, nr2, nr3, nr1s, nr2s, nr3s, ng_g, nk_g, &
      ngwk_g, nspin, nbnd, nel, nelu, neld, nat, ntyp, na, acc, nacc, &
      ecutwfc, ecutrho, alat, ekinc, kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, &
      ngauss, lgauss, ntetra, ltetra, natomwfc, gcutm, gcuts, dual, doublegrid, &
      modenum, lforce, lstres, title, crystal, tmp_dir, tupf, gamma_only, &
      noncolin, lspinorb, lda_plus_u, &
      tfixed_occ, tefield, dipfield, edir, emaxpos, eopreg, eamp, twfcollect )
!
      USE io_global, ONLY: ionode
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni  
      INTEGER, INTENT(IN) :: nfi   
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      INTEGER, INTENT(IN) :: nr1s, nr2s, nr3s
      INTEGER, INTENT(IN) :: ng_g
      REAL(DP), INTENT(IN) :: trutim  ! true time since last 'from_scratch'
      REAL(DP), INTENT(IN) :: ecutwfc, ecutrho  ! wfc and density cutoff
      REAL(DP), INTENT(IN) :: nel
      INTEGER, INTENT(IN) :: nk_g  ! global number of k points
      INTEGER, INTENT(IN) :: nspin, nbnd, nelu, neld, nat, ntyp
      INTEGER, INTENT(IN) :: ngwk_g(:)
      INTEGER, INTENT(IN) :: na(:)
      REAL(DP), INTENT(IN) :: acc(:)
      INTEGER, INTENT(IN) :: nacc
      REAL(DP), INTENT(IN) :: alat
      REAL(DP), INTENT(IN) :: ekinc

      INTEGER, INTENT(IN) ::  kunit, k1, k2, k3, nk1, nk2, nk3
      REAL(DP), INTENT(IN) :: dgauss
      INTEGER, INTENT(IN) :: ngauss
      LOGICAL, INTENT(IN) :: lgauss
      INTEGER, INTENT(IN) :: ntetra
      LOGICAL, INTENT(IN) :: ltetra

      INTEGER, INTENT(IN) :: natomwfc
      LOGICAL, INTENT(IN) :: doublegrid
      REAL(DP), INTENT(IN) :: gcutm, gcuts, dual
      INTEGER, INTENT(IN) :: modenum

      LOGICAL, INTENT(IN) :: lstres
      LOGICAL, INTENT(IN) :: lforce
      CHARACTER(LEN=*), INTENT(IN) :: title
      CHARACTER(LEN=*), INTENT(IN) :: crystal
      CHARACTER(LEN=*), INTENT(IN) :: tmp_dir

      !  tupf is .TRUE. if pseudo in restart file are saved in upf format
      LOGICAL, INTENT(IN) :: tupf        

      !  gamma_only is .TRUE. if calculation is at gamma (G-vecs span only half space)
      LOGICAL, INTENT(IN) :: gamma_only  
      LOGICAL, INTENT(IN) :: lda_plus_u
      LOGICAL, INTENT(IN) :: noncolin
      LOGICAL, INTENT(IN) :: lspinorb

      LOGICAL, INTENT(IN) :: tfixed_occ
      LOGICAL, INTENT(IN) :: tefield
      LOGICAL, INTENT(IN) :: dipfield
      INTEGER, INTENT(IN) :: edir
      REAL(DP), INTENT(IN) :: emaxpos
      REAL(DP), INTENT(IN) :: eopreg
      REAL(DP), INTENT(IN) :: eamp
      LOGICAL, INTENT(IN) :: twfcollect

      INTEGER :: i, j
      CHARACTER(LEN=80)  :: t_ , c_ 
      CHARACTER(LEN=256) :: tmp_dir_
      CHARACTER(LEN=30)  :: sub_name = ' write_restart_header '
      CHARACTER(LEN=20)  :: section_name = 'header'
      LOGICAL :: twrite = .TRUE.

      t_ = title
      c_ = crystal
      tmp_dir_ = tmp_dir

      IF( ntyp > SIZE( na ) ) &
        CALL errore(sub_name, ' wrong size: na ', 1 )
      IF( nk_g > SIZE( ngwk_g ) ) &
        CALL errore(sub_name, ' wrong size: ngwk_g ', 3 )
      IF( nacc > SIZE( acc ) ) &
        CALL errore(sub_name, ' wrong size: acc ', 4 )

      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) nfi, nr1, nr2, nr3, nr1s, nr2s, nr3s, ng_g, nk_g, &
          nspin, nbnd, nel, nelu, neld, nat, ntyp, nacc, trutim, ecutwfc, ecutrho, alat, ekinc,   &
          kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra,  &
          natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, lforce, tupf, &
          gamma_only, noncolin, lspinorb, lda_plus_u, & 
          tfixed_occ, tefield, dipfield, edir, emaxpos, eopreg, eamp, twfcollect
        WRITE(iuni) (na(i),i=1,ntyp), (ngwk_g(i),i=1,nk_g), (acc(i),i=1,nacc)
        WRITE(iuni) t_ , c_ , tmp_dir_ 
      END IF

      RETURN
    END SUBROUTINE write_restart_header1


    SUBROUTINE write_restart_header2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'header'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_header2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine read from disk dimensions, and status variables
!
    SUBROUTINE read_restart_header1(iuni, nfi, trutim, nr1, nr2, nr3, &
      nr1s, nr2s, nr3s, ng_g, nk_g, ngwk_g, nspin, nbnd, nel, nelu, neld, &
      nat, ntyp, na, acc, nacc, ecutwfc, ecutrho, alat, ekinc, kunit, &
      k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra, &
      natomwfc, gcutm, gcuts, dual, doublegrid, modenum, &
      lforce, lstres, title, crystal, tmp_dir, tupf, gamma_only, noncolin, &
      lspinorb, lda_plus_u,&
      tfixed_occ, tefield, dipfield, edir, emaxpos, eopreg, eamp, twfcollect )

!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(OUT) :: nfi
      INTEGER, INTENT(OUT) :: nr1, nr2, nr3, ng_g
      INTEGER, INTENT(OUT) :: nr1s, nr2s, nr3s
      REAL(DP), INTENT(OUT) :: trutim
      REAL(DP), INTENT(OUT) :: ecutwfc, ecutrho
      REAL(DP), INTENT(OUT) :: nel
      INTEGER, INTENT(OUT) :: nk_g, nspin, nbnd, nelu, neld, nat, ntyp
      INTEGER, INTENT(OUT) :: ngwk_g(:)
      INTEGER, INTENT(OUT) :: na(:)
      REAL(DP), INTENT(OUT) :: acc(:)
      INTEGER, INTENT(OUT) :: nacc
      REAL(DP), INTENT(OUT) :: alat
      REAL(DP), INTENT(OUT) :: ekinc
      INTEGER, INTENT(OUT) ::  kunit, k1, k2, k3, nk1, nk2, nk3
      REAL(DP), INTENT(OUT) :: dgauss
      INTEGER, INTENT(OUT) :: ngauss
      LOGICAL, INTENT(OUT) :: lgauss
      INTEGER, INTENT(OUT) :: ntetra
      LOGICAL, INTENT(OUT) :: ltetra

      INTEGER, INTENT(OUT) :: natomwfc
      LOGICAL, INTENT(OUT) :: doublegrid
      REAL(DP), INTENT(OUT) :: gcutm, gcuts, dual
      INTEGER, INTENT(OUT) :: modenum

      LOGICAL, INTENT(OUT) :: lstres
      LOGICAL, INTENT(OUT) :: lforce
      CHARACTER(LEN=*), INTENT(OUT) :: title
      CHARACTER(LEN=*), INTENT(OUT) :: crystal
      CHARACTER(LEN=*), INTENT(OUT) :: tmp_dir
      LOGICAL, INTENT(OUT) :: tupf
      LOGICAL, INTENT(OUT) :: gamma_only
      LOGICAL, INTENT(OUT) :: lda_plus_u
      LOGICAL, INTENT(OUT) :: noncolin
      LOGICAL, INTENT(OUT) :: lspinorb

      LOGICAL, INTENT(OUT) :: tfixed_occ
      LOGICAL, INTENT(OUT) :: tefield
      LOGICAL, INTENT(OUT) :: dipfield
      INTEGER, INTENT(OUT) :: edir
      REAL(DP), INTENT(OUT) :: emaxpos
      REAL(DP), INTENT(OUT) :: eopreg
      REAL(DP), INTENT(OUT) :: eamp
      LOGICAL, INTENT(OUT) :: twfcollect

      CHARACTER(LEN=80)  :: t_, c_
      CHARACTER(LEN=256) :: tmp_dir_
!
      INTEGER :: i, j, ierr
      INTEGER :: idum = 0
      LOGICAL :: twrite_
      CHARACTER(LEN=30) :: sub_name = ' read_restart_header '
      CHARACTER(LEN=20) :: section_name = 'header'
      CHARACTER(LEN=20) :: section_name_

!
! ... Subroutine Body
!
      CALL data_section_head( iuni, section_name_, twrite_, ierr )
!
      IF( .NOT. twrite_ ) &
        CALL errore(sub_name,' Data Section not present in restart file ', 1)
!
        IF( ionode ) THEN
          READ(iuni) nfi, nr1, nr2, nr3, nr1s, nr2s, nr3s, ng_g, nk_g, &
            nspin, nbnd, nel, nelu, neld, nat, ntyp, nacc, trutim, ecutwfc, ecutrho, &
            alat, ekinc, kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, &
            ntetra, ltetra, natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, &
            lforce, tupf, gamma_only, noncolin, lspinorb, lda_plus_u,&
            tfixed_occ, tefield, dipfield, edir, emaxpos, eopreg, eamp, twfcollect
        END IF
!
        CALL mp_bcast( nfi, ionode_id )
        CALL mp_bcast( nr1, ionode_id )
        CALL mp_bcast( nr2, ionode_id )
        CALL mp_bcast( nr3, ionode_id )
        CALL mp_bcast( nr1s, ionode_id )
        CALL mp_bcast( nr2s, ionode_id )
        CALL mp_bcast( nr3s, ionode_id )
        CALL mp_bcast( ng_g, ionode_id )
        CALL mp_bcast( nk_g, ionode_id )
        CALL mp_bcast( nspin, ionode_id )
        CALL mp_bcast( nbnd, ionode_id )
        CALL mp_bcast( nel, ionode_id )
        CALL mp_bcast( nelu, ionode_id )
        CALL mp_bcast( neld, ionode_id )
        CALL mp_bcast( nat, ionode_id )
        CALL mp_bcast( ntyp, ionode_id )
        CALL mp_bcast( nacc, ionode_id )
        CALL mp_bcast( trutim, ionode_id )
        CALL mp_bcast( ecutwfc, ionode_id )
        CALL mp_bcast( ecutrho, ionode_id )
        CALL mp_bcast( alat, ionode_id )
        CALL mp_bcast( ekinc, ionode_id )
        CALL mp_bcast( kunit, ionode_id )
        CALL mp_bcast( k1, ionode_id )
        CALL mp_bcast( k2, ionode_id )
        CALL mp_bcast( k3, ionode_id )
        CALL mp_bcast( nk1, ionode_id )
        CALL mp_bcast( nk2, ionode_id )
        CALL mp_bcast( nk3, ionode_id )
        CALL mp_bcast( dgauss, ionode_id )
        CALL mp_bcast( ngauss, ionode_id )
        CALL mp_bcast( lgauss, ionode_id )
        CALL mp_bcast( ntetra, ionode_id )
        CALL mp_bcast( ltetra, ionode_id )
        CALL mp_bcast( natomwfc, ionode_id )
        CALL mp_bcast( gcutm, ionode_id )
        CALL mp_bcast( gcuts, ionode_id )
        CALL mp_bcast( dual, ionode_id )
        CALL mp_bcast( doublegrid, ionode_id )
        CALL mp_bcast( modenum, ionode_id )
        CALL mp_bcast( lstres, ionode_id )
        CALL mp_bcast( lforce, ionode_id )
        CALL mp_bcast( tupf, ionode_id )
        CALL mp_bcast( gamma_only, ionode_id )
        CALL mp_bcast( noncolin, ionode_id )
        CALL mp_bcast( lspinorb, ionode_id )
        CALL mp_bcast( lda_plus_u, ionode_id )

        CALL mp_bcast( tfixed_occ, ionode_id ) 
        CALL mp_bcast( tefield, ionode_id ) 
        CALL mp_bcast( dipfield, ionode_id ) 
        CALL mp_bcast( edir, ionode_id ) 
        CALL mp_bcast( emaxpos, ionode_id ) 
        CALL mp_bcast( eopreg, ionode_id ) 
        CALL mp_bcast( eamp, ionode_id ) 
        CALL mp_bcast( twfcollect, ionode_id )
!
        IF( ntyp > SIZE( na ) ) &
          CALL errore(sub_name,' too many types ', ntyp )
        IF( nk_g > SIZE( ngwk_g ) ) &
          CALL errore(sub_name,' too many k points ', nk_g )
        IF( nacc > SIZE( acc ) ) &
          CALL errore(sub_name,' too many accumulators ', nacc )
!
        IF( ionode ) THEN
          READ(iuni) (na(i),i=1,ntyp), (ngwk_g(i),i=1,nk_g), (acc(i),i=1,nacc)
        END IF

        CALL mp_bcast( na, ionode_id )
        CALL mp_bcast( ngwk_g, ionode_id )
        CALL mp_bcast( acc, ionode_id )

        IF( ionode ) THEN
          READ(iuni) t_, c_, tmp_dir_ 
        END IF

        CALL mp_bcast( t_ , ionode_id )
        CALL mp_bcast( c_ , ionode_id )
        CALL mp_bcast( tmp_dir_ , ionode_id )

        title   = t_
        crystal = c_
        tmp_dir = tmp_dir_

      RETURN
    END SUBROUTINE read_restart_header1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_header2(iuni)
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum, ierr
      CHARACTER(LEN=30) :: sub_name = ' read_restart_header '
      CHARACTER(LEN=20) :: section_name = 'header'
      CHARACTER(LEN=20) :: section_name_
!
      CALL data_section_head( iuni, section_name_, twrite_, ierr )
      IF( ierr == 1 ) &
        CALL errore( sub_name, ' Wrong Data Section, '//section_name_//' instead of '//section_name, 1)
      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF

      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout, fmt = " (3X,'W: read_restart_header, header not read from restart ' ) " )
!
      RETURN
    END SUBROUTINE read_restart_header2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_xdim1(iuni, &
      npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: npwx, nbndx
      INTEGER, INTENT(IN) :: nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
      CHARACTER(LEN=20) :: section_name = 'xdim'
      LOGICAL :: twrite = .TRUE.
!
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
      END IF
!
      RETURN
    END SUBROUTINE write_restart_xdim1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_xdim2(iuni)
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'xdim'
!
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
      END IF 
!
      RETURN 
    END SUBROUTINE write_restart_xdim2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_xdim1(iuni, &
      npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs )

!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(OUT) :: npwx, nbndx
      INTEGER, INTENT(OUT) :: nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
      LOGICAL :: twrite_ 
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'xdim'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!
      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_xdim ',' Data Section not present in restart file ', 1)

        IF( ionode ) THEN
          READ(iuni) npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
        END IF
        CALL mp_bcast( npwx, ionode_id )
        CALL mp_bcast( nbndx, ionode_id )
        CALL mp_bcast( nrx1, ionode_id )
        CALL mp_bcast( nrx2, ionode_id )
        CALL mp_bcast( nrx3, ionode_id )
        CALL mp_bcast( nrxx, ionode_id )
        CALL mp_bcast( nrx1s, ionode_id )
        CALL mp_bcast( nrx2s, ionode_id )
        CALL mp_bcast( nrx3s, ionode_id )
        CALL mp_bcast( nrxxs, ionode_id )

      RETURN
    END SUBROUTINE read_restart_xdim1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_xdim2(iuni)

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum, ierr
      CHARACTER(LEN=20) :: section_name = 'xdim'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
      END IF

      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_xdim, xdim not read from restart ' )")

      RETURN
    END SUBROUTINE read_restart_xdim2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variable related to the lda+U calculation
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_ldaU1(iuni, &
         ntyp, Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha)
!
      USE io_global, ONLY: ionode
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ntyp
      INTEGER, INTENT(IN) :: Hubbard_lmax
      INTEGER, INTENT(IN) :: Hubbard_l(:)
      REAL(DP), INTENT(IN) :: Hubbard_U(:)
      REAL(DP), INTENT(IN) :: Hubbard_alpha(:)

      INTEGER :: i
      CHARACTER(LEN=30) :: sub_name = ' write_restart_ldaU '
      CHARACTER(LEN=20) :: section_name = 'ldaU'

      LOGICAL :: twrite = .TRUE.

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) ntyp
          WRITE(iuni) Hubbard_lmax
          WRITE(iuni) (Hubbard_l(i),i=1,ntyp)
          WRITE(iuni) (Hubbard_U(i),i=1,ntyp)
          WRITE(iuni) (Hubbard_alpha(i),i=1,ntyp)
       ENDIF

      RETURN
    END SUBROUTINE write_restart_ldau1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_ldaU2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'ldaU'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_ldau2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_ldaU1(iuni, &
         ntyp, Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha)
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(OUT) :: ntyp
      INTEGER, INTENT(OUT) :: Hubbard_lmax
      INTEGER, INTENT(OUT) :: Hubbard_l(:)
      REAL(DP), INTENT(OUT) :: Hubbard_U(:)
      REAL(DP), INTENT(OUT) :: Hubbard_alpha(:)

      INTEGER :: i
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_ldaU '
      CHARACTER(LEN=20) :: section_name = 'ldaU'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_ldaU ',' Data Section not present in restart file ', 1)
!
! some checks should be added here
!

        IF( ionode ) THEN
          READ(iuni) ntyp
        END IF
        CALL mp_bcast(ntyp, ionode_id)

        IF( ionode ) THEN
          READ(iuni) Hubbard_lmax
          READ(iuni) (Hubbard_l(i),i=1,ntyp)
          READ(iuni) (Hubbard_U(i),i=1,ntyp)
          READ(iuni) (Hubbard_alpha(i),i=1,ntyp)
       ENDIF
        CALL mp_bcast(Hubbard_lmax, ionode_id)
        CALL mp_bcast(Hubbard_l, ionode_id)
        CALL mp_bcast(Hubbard_U, ionode_id)
        CALL mp_bcast(Hubbard_alpha, ionode_id)

      RETURN
    END SUBROUTINE read_restart_ldau1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_ldaU2(iuni)

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'ldaU'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_ldaU, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_ldau2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_symmetry1(iuni, &
      symm_type, sname, s, irt, nat, ftau, nsym, invsym, noinv )

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      CHARACTER(LEN=9), INTENT(IN) :: symm_type
      INTEGER, INTENT(IN) :: s(:,:,:)
      CHARACTER(LEN=45), INTENT(IN) :: sname(:)
      INTEGER, INTENT(IN) :: ftau(:,:)
      INTEGER, INTENT(IN) :: irt(:,:)
      INTEGER, INTENT(IN) :: nsym
      INTEGER, INTENT(IN) :: nat
      LOGICAL, INTENT(IN) :: invsym
      LOGICAL, INTENT(IN) :: noinv
      INTEGER :: i,j
      CHARACTER(LEN=30) :: sub_name = ' write_restart_symmetry '
      CHARACTER(LEN=20) :: section_name = 'symmetry'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

        IF( SIZE(s,3) < nsym ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( ( SIZE(irt,1) < nsym ) .OR. ( SIZE(irt,2) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE(ftau,2) < nsym ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( SIZE(sname) < nsym ) &
          CALL errore( sub_name, ' wrong size ', 4 )

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) symm_type, nsym, invsym, noinv, nat
          !
          ! We write all possible sym.ops., and not only true symmetries,
          ! because in the third-order code D3 we may need all of them
          !
          WRITE(iuni) (s(:,:,i),i=1,48), ((irt(i,j),i=1,48),j=1,nat),  &
            (ftau(:,i),i=1,48), (sname(i),i=1,48)
        END IF

      RETURN
    END SUBROUTINE write_restart_symmetry1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_symmetry2( iuni )
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'symmetry'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_symmetry2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_symmetry1(iuni, &
      symm_type, sname, s, irt, nat, ftau, nsym, invsym, noinv )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      CHARACTER(LEN=9),  INTENT(OUT) :: symm_type
      CHARACTER(LEN=45), INTENT(OUT) :: sname(:)
      INTEGER, INTENT(OUT) :: s(:,:,:)
      INTEGER, INTENT(OUT) :: irt(:,:)
      INTEGER, INTENT(OUT) :: ftau(:,:)
      INTEGER, INTENT(OUT) :: nsym
      INTEGER, INTENT(OUT) :: nat
      LOGICAL, INTENT(OUT) :: invsym
      LOGICAL, INTENT(OUT) :: noinv

      LOGICAL :: twrite_
      INTEGER :: i, j
      INTEGER :: idum, ierr
      CHARACTER(LEN=30) :: sub_name = ' read_restart_symmetry '
      CHARACTER(LEN=20) :: section_name = 'symmetry'
      CHARACTER(LEN=20) :: section_name_ 
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_symmetry ',' symmetries not present in restart file ', 1)

        IF( ionode ) THEN
          READ(iuni) symm_type, nsym, invsym, noinv, nat
        END IF
        CALL mp_bcast( symm_type, ionode_id )
        CALL mp_bcast( nsym, ionode_id )
        CALL mp_bcast( invsym, ionode_id )
        CALL mp_bcast( noinv, ionode_id )
        CALL mp_bcast( nat, ionode_id )
!
! ...   Check dummy variables
        IF( SIZE(s,3) < nsym ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( SIZE(ftau,2) < nsym ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE(sname) < nsym ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( ( SIZE(irt,1) < nsym ) .OR. ( SIZE(irt,2) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 5 )
!
        IF( ionode ) THEN
          READ(iuni) (s(:,:,i),i=1,48), ((irt(i,j),i=1,48),j=1,nat),  &
            (ftau(:,i),i=1,48), (sname(i),i=1,48)
        END IF
        CALL mp_bcast( s, ionode_id )
        CALL mp_bcast( sname(:), ionode_id )
        CALL mp_bcast( ftau, ionode_id )
        CALL mp_bcast( irt, ionode_id )
!
      RETURN
    END SUBROUTINE read_restart_symmetry1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_symmetry2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum, ierr
      CHARACTER(LEN=20) :: section_name = 'symmetry'
      CHARACTER(LEN=20) :: section_name_ 
      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )
      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_symmetry, symmetries not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_symmetry2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_pseudo1(iuni, &
      zmesh, xmin, dx, r, rab, vloc_at, chi, oc, rho_at, &
      rho_atc, mesh, msh, nchi, lchi, jchi, numeric, cc, alpc, zp, aps, alps, zv, nlc, &
      nnl, lmax, lloc, dion, betar, qqq, qfunc, qfcoef, rinner, nh, nbeta, &
      kkbeta, nqf, nqlc, ifqopt, lll, jjj, iver, tvanp, okvan, newpseudo, iexch, icorr, &
      igcx, igcc, lsda, a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd)
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      REAL(DP), INTENT(IN) :: zmesh, xmin, dx
      REAL(DP), INTENT(IN) :: r(:), rab(:), vloc_at(:), chi(:,:)
      REAL(DP), INTENT(IN) :: oc(:), rho_at(:), rho_atc(:)
      INTEGER, INTENT(IN) :: mesh, msh, nchi, lchi(:)
      REAL(DP), INTENT(IN) :: jchi(:)
      LOGICAL, INTENT(IN) :: numeric
      REAL(DP), INTENT(IN) :: cc(2), alpc(2), zp, aps(6,0:3), alps(3,0:3), zv
      INTEGER, INTENT(IN) :: nlc, nnl, lmax, lloc
      REAL(DP), INTENT(IN) :: dion(:,:), betar(:,:), qqq(:,:), qfunc(:,:,:)
      REAL(DP), INTENT(IN) :: qfcoef(:,:,:,:), rinner(:)
      INTEGER, INTENT(IN) :: nh, nbeta, kkbeta, nqf, nqlc, ifqopt, lll(:), iver(3)
      REAL(DP), INTENT(IN) :: jjj(:)
      LOGICAL, INTENT(IN) :: tvanp, okvan, newpseudo
      INTEGER, INTENT(IN) :: iexch, icorr, igcx, igcc
      LOGICAL, INTENT(IN) :: lsda
      REAL(DP), INTENT(IN) :: a_nlcc, b_nlcc, alpha_nlcc
      LOGICAL, INTENT(IN) :: nlcc
      CHARACTER(LEN=2), INTENT(IN) :: psd
!
      INTEGER :: mesh_, lloc_, nchi_, nbeta_, nqf_, nqlc_
      CHARACTER(LEN=30) :: sub_name = ' write_restart_pseudo '
      CHARACTER(LEN=20) :: section_name = 'pseudo'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

! ...   Check dummy variables
        mesh_  = MAX( mesh, 1 )
        lloc_  = MAX( lloc, 1 )
        nchi_  = MAX( nchi, 1 )
        nbeta_ = MAX( nbeta, 1 )
        nqf_   = MAX( nqf, 1 )
        nqlc_  = MAX( nqlc, 1 )

        IF( SIZE(r) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( SIZE(rab) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE(vloc_at) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( ( SIZE(chi,1) < mesh_ ) .OR. ( SIZE(chi,2) < nchi_ ) ) &
          CALL errore( sub_name, ' wrong size ', 4 )
        IF( SIZE(oc) < nchi_ ) &
          CALL errore( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_at) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 6 )
        IF( SIZE(rho_atc) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 7 )
        IF( SIZE(lchi) < nchi_ ) &
          CALL errore( sub_name, ' wrong size ', 8 )
        IF( ( SIZE(dion,1) < nbeta_ ) .OR. ( SIZE(dion,2) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(betar,1) < mesh_ ) .OR. ( SIZE(betar,2) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(qqq,1) < nbeta_  ) .OR. ( SIZE(qqq,2) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( ( SIZE(qfunc,1) < mesh_ ) .OR. ( SIZE(qfunc,2) < nbeta_ ) .OR. &
            ( SIZE(qfunc,3) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 11 )
        IF( ( SIZE(qfcoef,1) < nqf_ ) .OR. ( SIZE(qfcoef,2) < nqlc_ ) .OR. &
            ( SIZE(qfcoef,3) < nbeta_ ) .OR. ( SIZE(qfcoef,4) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 12 )
        IF( SIZE(rinner) < nqlc_ ) &
          CALL errore( sub_name, ' wrong size ', 13 )
        IF( SIZE(lll) < nbeta_ ) &
          CALL errore( sub_name, ' wrong size ', 14 )

        IF( ionode ) THEN

          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) zmesh, xmin, dx, mesh, msh, nchi, numeric, zp, zv, &
            nlc, nnl, lmax, lloc,  nh, nbeta, kkbeta, nqf, nqlc, ifqopt, &
            tvanp, okvan, newpseudo, iexch, icorr, igcx, igcc, lsda,     &
            a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd
          WRITE(iuni) r( 1:mesh_ ), rab( 1:mesh_ ), &
            vloc_at( 1:mesh_ ), chi( 1:mesh_, 1:nchi_ ), &
            oc( 1:nchi_ ), rho_at( 1:mesh_ ), rho_atc( 1:mesh_ ), &
            lchi( 1:nchi_ ), jchi(1:nchi_)
          WRITE(iuni) cc(1:2), alpc(1:2), aps(1:6,0:3), alps(1:3,0:3)
          WRITE(iuni) dion( 1:nbeta_, 1:nbeta_ ), betar( 1:mesh_, 1:nbeta_ ), &
            qqq( 1:nbeta_, 1:nbeta_ ), qfunc( 1:mesh_, 1:nbeta_, 1:nbeta_ ), &
            qfcoef( 1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_ ), &
            rinner( 1:nqlc_ ), lll( 1:nbeta_ ), jjj(1:nbeta_), iver(1:3)

        END IF

      RETURN
    END SUBROUTINE write_restart_pseudo1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_pseudo3(iuni, &
      generated, date_author, comment, psd, typ, tvanp, nlcc, dft, zp, etotps, &
      ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, els, lchi, jchi, &
      oc, r, rab, &
      rho_atc, vloc, lll, jjj, kkbeta, beta, nd, dion, nqf, nqlc, rinner, qqq, &
      qfunc, qfcoef, chi, rho_at )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
!
      CHARACTER(LEN=80):: generated   ! 
      CHARACTER(LEN=80):: date_author ! Misc info
      CHARACTER(LEN=80):: comment     !
      CHARACTER(LEN=2) :: psd       ! Element label
      CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
      LOGICAL :: tvanp             ! .true. if Ultrasoft
      LOGICAL :: nlcc               ! Non linear core corrections
      CHARACTER(LEN=20) :: dft      ! Exch-Corr type
      REAL(DP) :: zp               ! z valence
      REAL(DP) :: etotps           ! total energy
      REAL(DP) :: ecutwfc          ! suggested cut-off for wfc
      REAL(DP) :: ecutrho          ! suggested cut-off for rho
      INTEGER :: nv                 ! UPF file version number
      INTEGER :: lmax               ! maximum angular momentum component
      INTEGER :: mesh               ! number of point in the radial mesh
      INTEGER :: nwfc               ! number of wavefunctions
      INTEGER :: nbeta              ! number of projectors
      CHARACTER(LEN=2) :: els(:)  ! els(nwfc)
      INTEGER :: lchi(:)   ! lchi(nwfc)
      REAL(DP) :: jchi(:)   ! jchi(nwfc)
      REAL(DP) :: oc(:)   ! oc(nwfc)
      REAL(DP) :: r(:)    ! r(mesh)
      REAL(DP) :: rab(:)  ! rab(mesh)
      REAL(DP) :: rho_atc(:) ! rho_atc(mesh)
      REAL(DP) :: vloc(:)    ! vloc(mesh)
      INTEGER :: lll(:)       ! lll(nbeta)
      REAL(DP) :: jjj(:)     ! jjj(nbeta)
      INTEGER :: kkbeta(:)    ! kkbeta(nbeta)
      REAL(DP) :: beta(:,:)  ! beta(mesh,nbeta)
      INTEGER :: nd
      REAL(DP) :: dion(:,:)  ! dion(nbeta,nbeta)
      INTEGER :: nqf
      INTEGER :: nqlc
      REAL(DP) :: rinner(:)  ! rinner(0:2*lmax)
      REAL(DP) :: qqq(:,:)   ! qqq(nbeta,nbeta)
      REAL(DP) :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
      REAL(DP) :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
      REAL(DP) :: chi(:,:) !  chi(mesh,nwfc)
      REAL(DP) :: rho_at(:) !  rho_at(mesh)
!
      INTEGER :: idum = 0
      CHARACTER(LEN=30) :: sub_name = ' write_restart_pseudo '
      CHARACTER(LEN=20) :: section_name = 'pseudo'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

! ...   Check dummy variables
        IF( SIZE(els) < nwfc ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( SIZE(lchi) < nwfc ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE(oc) < nwfc ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( SIZE(r) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 4 )
        IF( SIZE(rab) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_atc) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 6 )
        IF(  SIZE(vloc) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 7 )
        IF( SIZE(lll) < nbeta ) &
          CALL errore( sub_name, ' wrong size ', 8 )
        IF( SIZE(kkbeta) < nbeta ) &
          CALL errore( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(beta,1) < mesh ) .OR. ( SIZE(beta,2) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(dion,1) < nbeta ) .OR. ( SIZE(dion,2) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 11 )
        IF( SIZE(rinner) < nqlc ) &
          CALL errore( sub_name, ' wrong size ', 12 )
        IF( ( SIZE(qqq,1) < nbeta ) .OR. ( SIZE(qqq,2) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 13 )
        IF( ( SIZE(qfunc,1) < mesh ) .OR. ( SIZE(qfunc,2) < nbeta ) .OR. &
            ( SIZE(qfunc,3) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 14 )
        IF( ( SIZE(qfcoef,1) < nqf ) .OR. ( SIZE(qfcoef,2) < nqlc ) .OR. &
            ( SIZE(qfcoef,3) < nbeta ) .OR. ( SIZE(qfcoef,4) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 15 )
        IF( ( SIZE(chi,1) < mesh ) .OR. ( SIZE(chi,2) < nwfc ) ) &
          CALL errore( sub_name, ' wrong size ', 16 )
        IF( SIZE(rho_at) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 17 )

        IF( ionode ) THEN

          WRITE(iuni) twrite, file_version, section_name
!
          WRITE(iuni) generated, date_author, comment, psd, typ, tvanp, nlcc, dft, &
           zp, etotps, ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, nd, nqf, nqlc
!           
          WRITE(iuni) els(1:nwfc), lchi(1:nwfc), jchi(1:nwfc), oc(1:nwfc), &
            r(1:mesh), rab(1:mesh), &
            rho_atc(1:mesh), vloc(1:mesh), lll(1:nbeta), jjj(1:nbeta), &
            kkbeta(1:nbeta), &
            beta(1:mesh,1:nbeta), &
            dion(1:nbeta,1:nbeta), rinner(1:nqlc), qqq(1:nbeta,1:nbeta), &
            qfunc(1:mesh, 1:nbeta, 1:nbeta), qfcoef(1:nqf, 1:nqlc, 1:nbeta, 1:nbeta), &
            chi(1:mesh, 1:nwfc), rho_at(1:mesh) 

          WRITE(iuni) idum
          WRITE(iuni) idum

        END IF

      RETURN
    END SUBROUTINE write_restart_pseudo3

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_pseudo2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'pseudo'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_pseudo2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_pseudo1(iuni, &
      zmesh, xmin, dx, r, rab, vloc_at, chi, oc, rho_at, &
      rho_atc, mesh, msh, nchi, lchi, jchi, numeric, cc, alpc, zp, aps, alps, zv, nlc, &
      nnl, lmax, lloc, dion, betar, qqq, qfunc, qfcoef, rinner, nh, nbeta, &
      kkbeta, nqf, nqlc, ifqopt, lll, jjj, iver, tvanp, okvan, newpseudo, iexch, icorr, &
      igcx, igcc, lsda, a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd )

!

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      REAL(DP), INTENT(OUT) :: zmesh, xmin, dx
      REAL(DP), INTENT(OUT) :: r(:), rab(:), vloc_at(:), chi(:,:)
      REAL(DP), INTENT(OUT) :: oc(:), rho_at(:), rho_atc(:)
      INTEGER, INTENT(OUT) :: mesh, msh, nchi, lchi(:)
      REAL(DP), INTENT(OUT) :: jchi(:)
      LOGICAL, INTENT(OUT) :: numeric
      REAL(DP), INTENT(OUT) :: cc(2), alpc(2), zp, aps(6,0:3), alps(3,0:3), zv
      INTEGER, INTENT(OUT) :: nlc, nnl, lmax, lloc
      REAL(DP), INTENT(OUT) :: dion(:,:), betar(:,:), qqq(:,:), qfunc(:,:,:)
      REAL(DP), INTENT(OUT) :: qfcoef(:,:,:,:), rinner(:)
      INTEGER, INTENT(OUT) :: nh, nbeta, kkbeta, nqf, nqlc, ifqopt, &
                              lll(:), iver(:)
      REAL(DP), INTENT(OUT) :: jjj(:)
      LOGICAL, INTENT(OUT) :: tvanp, okvan, newpseudo
      INTEGER, INTENT(OUT) :: iexch, icorr, igcx, igcc
      LOGICAL, INTENT(OUT) :: lsda
      REAL(DP), INTENT(OUT) :: a_nlcc, b_nlcc, alpha_nlcc
      LOGICAL, INTENT(OUT) :: nlcc
      CHARACTER(LEN=2), INTENT(OUT) :: psd
     
!
      LOGICAL :: twrite_
!
      INTEGER :: idum, ierr
      CHARACTER(LEN=30) :: sub_name = ' read_restart_pseudo '
      CHARACTER(LEN=20) :: section_name = 'pseudo'
      CHARACTER(LEN=20) :: section_name_

      INTEGER :: mesh_, lloc_, nchi_, nbeta_, nqf_, nqlc_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_pseudo ',' pseudo not present in restart file ', 1)


        IF( ionode ) THEN
          READ(iuni) zmesh, xmin, dx, mesh, msh, nchi, numeric, zp, zv, nlc, nnl, lmax, &
            lloc, nh, nbeta, kkbeta, nqf, nqlc, ifqopt, tvanp, okvan, newpseudo, &
            iexch, icorr, igcx, igcc, lsda, a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd
        END IF

        CALL mp_bcast( zmesh, ionode_id )
        CALL mp_bcast( xmin, ionode_id )
        CALL mp_bcast( dx, ionode_id )
        CALL mp_bcast( mesh, ionode_id )
        CALL mp_bcast( msh, ionode_id )
        CALL mp_bcast( nchi, ionode_id )
        CALL mp_bcast( numeric, ionode_id )
        CALL mp_bcast( zp, ionode_id )
        CALL mp_bcast( zv, ionode_id )
        CALL mp_bcast( nlc, ionode_id )
        CALL mp_bcast( nnl, ionode_id )
        CALL mp_bcast( lmax, ionode_id )
        CALL mp_bcast( lloc, ionode_id )
        CALL mp_bcast( nh, ionode_id )
        CALL mp_bcast( nbeta, ionode_id )
        CALL mp_bcast( kkbeta, ionode_id )
        CALL mp_bcast( nqf, ionode_id )
        CALL mp_bcast( nqlc, ionode_id )
        CALL mp_bcast( ifqopt, ionode_id )
        CALL mp_bcast( tvanp, ionode_id )
        CALL mp_bcast( okvan, ionode_id )
        CALL mp_bcast( newpseudo, ionode_id )
        CALL mp_bcast( iexch, ionode_id )
        CALL mp_bcast( icorr, ionode_id )
        CALL mp_bcast( igcx, ionode_id )
        CALL mp_bcast( igcc, ionode_id )
        CALL mp_bcast( lsda, ionode_id )
        CALL mp_bcast( a_nlcc, ionode_id )
        CALL mp_bcast( b_nlcc, ionode_id )
        CALL mp_bcast( alpha_nlcc, ionode_id )
        CALL mp_bcast( nlcc, ionode_id )
        CALL mp_bcast( psd, ionode_id )


        mesh_  = MAX( mesh, 1 )
        lloc_  = MAX( lloc, 1 )
        nchi_  = MAX( nchi, 1 )
        nbeta_ = MAX( nbeta, 1 )
        nqf_   = MAX( nqf, 1 )
        nqlc_  = MAX( nqlc, 1 )

        IF( mesh < 0 ) &
          CALL errore( sub_name, ' wrong value ', 1 )
        IF( lloc < 0 ) &
          CALL errore( sub_name, ' wrong value ', 2 )
        IF( nchi < 0 ) &
          CALL errore( sub_name, ' wrong value ', 3 )
        IF( nbeta < 0 ) &
          CALL errore( sub_name, ' wrong value ', 4 )
        IF( nqf < 0 ) &
          CALL errore( sub_name, ' wrong value ', 5 )
        IF( nqlc < 0 ) &
          CALL errore( sub_name, ' wrong value ', 6 )

! ...   Check dummy variables
        IF( SIZE(r) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( SIZE(rab) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE(vloc_at) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( ( SIZE(chi,1) < mesh_ ) .OR. ( SIZE(chi,2) < nchi_ ) ) &
          CALL errore( sub_name, ' wrong size ', 4 )
        IF( SIZE(oc) < nchi_ ) &
          CALL errore( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_at) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 6 )
        IF( SIZE(rho_atc) < mesh_ ) &
          CALL errore( sub_name, ' wrong size ', 7 )
        IF( SIZE(lchi) < nchi_ ) &
          CALL errore( sub_name, ' wrong size ', 8 )
        IF( ( SIZE(dion,1) < nbeta_ ) .OR. ( SIZE(dion,2) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(betar,1) < mesh_ ) .OR. ( SIZE(betar,2) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(qqq,1) < nbeta_ ) .OR. ( SIZE(qqq,2) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( ( SIZE(qfunc,1) < mesh_ ) .OR. ( SIZE(qfunc,2) < nbeta_ ) .OR. &
            ( SIZE(qfunc,3) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 11 )
        IF( ( SIZE(qfcoef,1) < nqf_ ) .OR. ( SIZE(qfcoef,2) < nqlc_ ) .OR. &
            ( SIZE(qfcoef,3) < nbeta_ ) .OR. ( SIZE(qfcoef,4) < nbeta_ ) ) &
          CALL errore( sub_name, ' wrong size ', 12 )
        IF( SIZE(rinner) < nqlc_ ) &
          CALL errore( sub_name, ' wrong size ', 13 )
        IF( SIZE(lll) < nbeta_ ) &
          CALL errore( sub_name, ' wrong size ', 14 )

        IF( ionode ) THEN
          READ(iuni) r(1:mesh_), rab(1:mesh_), vloc_at(1:mesh_), chi(1:mesh_,1:nchi_), &
            oc(1:nchi_), rho_at(1:mesh_), rho_atc(1:mesh_), lchi(1:nchi_), &
            jchi(1:nchi_)
          READ(iuni) cc(1:2), alpc(1:2), aps(1:6,0:3), alps(1:3,0:3)
          READ(iuni) dion(1:nbeta_,1:nbeta_), betar(1:mesh_,1:nbeta_), qqq(1:nbeta_,1:nbeta_), &
            qfunc(1:mesh_, 1:nbeta_, 1:nbeta_), qfcoef(1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_), &
            rinner(1:nqlc_), lll(1:nbeta_), jjj(1:nbeta_), iver(1:3)
        END IF

        CALL mp_bcast( r, ionode_id )
        CALL mp_bcast( rab, ionode_id )
        CALL mp_bcast( vloc_at, ionode_id )
        CALL mp_bcast( chi, ionode_id )
        CALL mp_bcast( oc, ionode_id )
        CALL mp_bcast( rho_at, ionode_id )
        CALL mp_bcast( rho_atc, ionode_id )
        CALL mp_bcast( lchi, ionode_id )
        CALL mp_bcast( jchi, ionode_id )
        CALL mp_bcast( cc, ionode_id )
        CALL mp_bcast( alpc, ionode_id )
        CALL mp_bcast( aps, ionode_id )
        CALL mp_bcast( alps, ionode_id )
        CALL mp_bcast( dion, ionode_id )
        CALL mp_bcast( betar, ionode_id )
        CALL mp_bcast( qqq, ionode_id )
        CALL mp_bcast( qfunc, ionode_id )
        CALL mp_bcast( qfcoef, ionode_id )
        CALL mp_bcast( rinner, ionode_id )
        CALL mp_bcast( lll, ionode_id )
        CALL mp_bcast( jjj, ionode_id )
        CALL mp_bcast( iver, ionode_id )

      RETURN
    END SUBROUTINE read_restart_pseudo1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_pseudo3(iuni, &
      generated, date_author, comment, psd, typ, tvanp, nlcc, dft, zp, etotps, &
      ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, els, lchi, jchi, &
      oc, r, rab, rho_atc, vloc, lll, jjj, kkbeta, beta, nd, dion, nqf, &
      nqlc, rinner, qqq, qfunc, qfcoef, chi, rho_at )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
!
      CHARACTER(LEN=80):: generated   ! 
      CHARACTER(LEN=80):: date_author ! Misc info
      CHARACTER(LEN=80):: comment     !
      CHARACTER(LEN=2) :: psd       ! Element label
      CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
      LOGICAL :: tvanp             ! .true. if Ultrasoft
      LOGICAL :: nlcc               ! Non linear core corrections
      CHARACTER(LEN=20) :: dft      ! Exch-Corr type
      REAL(DP) :: zp               ! z valence
      REAL(DP) :: etotps           ! total energy
      REAL(DP) :: ecutwfc          ! suggested cut-off for wfc
      REAL(DP) :: ecutrho          ! suggested cut-off for rho
      INTEGER :: nv                 ! UPF file version number
      INTEGER :: lmax               ! maximum angular momentum component
      INTEGER :: mesh               ! number of point in the radial mesh
      INTEGER :: nwfc               ! number of wavefunctions
      INTEGER :: nbeta              ! number of projectors
      CHARACTER(LEN=2) :: els(:)  ! els(nwfc)
      INTEGER :: lchi(:)   ! lchi(nwfc)
      REAL(DP) :: jchi(:)   ! jchi(nwfc)
      REAL(DP) :: oc(:)   ! oc(nwfc)
      REAL(DP) :: r(:)    ! r(mesh)
      REAL(DP) :: rab(:)  ! rab(mesh)
      REAL(DP) :: rho_atc(:) ! rho_atc(mesh)
      REAL(DP) :: vloc(:)    ! vloc(mesh)
      INTEGER :: lll(:)       ! lll(nbeta)
      REAL(DP) :: jjj(:)    ! jjj(nbeta)
      INTEGER :: kkbeta(:)    ! kkbeta(nbeta)
      REAL(DP) :: beta(:,:)  ! beta(mesh,nbeta)
      INTEGER :: nd
      REAL(DP) :: dion(:,:)  ! dion(nbeta,nbeta)
      INTEGER :: nqf
      INTEGER :: nqlc
      REAL(DP) :: rinner(:)  ! rinner(0:2*lmax)
      REAL(DP) :: qqq(:,:)   ! qqq(nbeta,nbeta)
      REAL(DP) :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
      REAL(DP) :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
      REAL(DP) :: chi(:,:) !  chi(mesh,nwfc)
      REAL(DP) :: rho_at(:) !  rho_at(mesh)
!
!
!
      LOGICAL :: twrite_
      INTEGER :: idum, ierr
      CHARACTER(LEN=30) :: sub_name = ' read_restart_pseudo '
      CHARACTER(LEN=20) :: section_name = 'pseudo'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(sub_name, ' pseudo not present in restart file ', 1)


        IF( ionode ) THEN
          READ(iuni) generated, date_author, comment, psd, typ, tvanp, nlcc, dft, &
           zp, etotps, ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, nd, nqf, nqlc
        END IF
        CALL mp_bcast( generated, ionode_id )
        CALL mp_bcast( date_author, ionode_id )
        CALL mp_bcast( comment, ionode_id )
        CALL mp_bcast( psd, ionode_id )
        CALL mp_bcast( typ, ionode_id )
        CALL mp_bcast( tvanp, ionode_id )
        CALL mp_bcast( nlcc, ionode_id )
        CALL mp_bcast( dft, ionode_id )
        CALL mp_bcast( zp, ionode_id )
        CALL mp_bcast( etotps, ionode_id )
        CALL mp_bcast( ecutwfc, ionode_id )
        CALL mp_bcast( ecutrho, ionode_id )
        CALL mp_bcast( nv, ionode_id )
        CALL mp_bcast( lmax, ionode_id )
        CALL mp_bcast( mesh, ionode_id )
        CALL mp_bcast( nwfc, ionode_id )
        CALL mp_bcast( nbeta, ionode_id )
        CALL mp_bcast( nd, ionode_id )
        CALL mp_bcast( nqf, ionode_id )
        CALL mp_bcast( nqlc, ionode_id )

        IF( mesh < 0 ) &
          CALL errore( sub_name, ' wrong value ', 1 )
        IF( nwfc < 1 ) &
          CALL errore( sub_name, ' wrong value ', 2 )
        IF( nbeta < 1 ) &
          CALL errore( sub_name, ' wrong value ', 3 )
        IF( nqf < 1 ) &
          CALL errore( sub_name, ' wrong value ', 4 )
        IF( nqlc < 1 ) &
          CALL errore( sub_name, ' wrong value ', 5 )


! ...   Check dummy variables
        IF( SIZE(els) < nwfc ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( SIZE(lchi) < nwfc ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE(oc) < nwfc ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( SIZE(r) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 4 )
        IF( SIZE(rab) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_atc) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 6 )
        IF(  SIZE(vloc) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 7 )
        IF( SIZE(lll) < nbeta ) &
          CALL errore( sub_name, ' wrong size ', 8 )
        IF( SIZE(kkbeta) < nbeta ) &
          CALL errore( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(beta,1) < mesh ) .OR. ( SIZE(beta,2) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(dion,1) < nbeta ) .OR. ( SIZE(dion,2) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 11 )
        IF( SIZE(rinner) < nqlc ) &
          CALL errore( sub_name, ' wrong size ', 12 )
        IF( ( SIZE(qqq,1) < nbeta ) .OR. ( SIZE(qqq,2) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 13 )
        IF( ( SIZE(qfunc,1) < mesh ) .OR. ( SIZE(qfunc,2) < nbeta ) .OR. &
            ( SIZE(qfunc,3) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 14 )
        IF( ( SIZE(qfcoef,1) < nqf ) .OR. ( SIZE(qfcoef,2) < nqlc ) .OR. &
            ( SIZE(qfcoef,3) < nbeta ) .OR. ( SIZE(qfcoef,4) < nbeta ) ) &
          CALL errore( sub_name, ' wrong size ', 15 )
        IF( ( SIZE(chi,1) < mesh ) .OR. ( SIZE(chi,2) < nwfc ) ) &
          CALL errore( sub_name, ' wrong size ', 16 )
        IF( SIZE(rho_at) < mesh ) &
          CALL errore( sub_name, ' wrong size ', 17 )

        IF( ionode ) THEN
!           
          READ(iuni) els(1:nwfc), lchi(1:nwfc), jchi(1:nwfc), oc(1:nwfc), &
            r(1:mesh), rab(1:mesh), &
            rho_atc(1:mesh), vloc(1:mesh), lll(1:nbeta), jjj(1:nbeta), &
            kkbeta(1:nbeta), &
            beta(1:mesh,1:nbeta), &
            dion(1:nbeta,1:nbeta), rinner(1:nqlc), qqq(1:nbeta,1:nbeta), &
            qfunc(1:mesh, 1:nbeta, 1:nbeta), &
            qfcoef(1:nqf, 1:nqlc, 1:nbeta, 1:nbeta), &
            chi(1:mesh, 1:nwfc), rho_at(1:mesh) 

          READ(iuni) idum
          READ(iuni) idum

        END IF

        CALL mp_bcast( els(1:nwfc), ionode_id ) 
        CALL mp_bcast( lchi(1:nwfc), ionode_id ) 
        CALL mp_bcast( jchi(1:nwfc), ionode_id ) 
        CALL mp_bcast( oc(1:nwfc), ionode_id ) 
        CALL mp_bcast( r(1:mesh), ionode_id ) 
        CALL mp_bcast( rab(1:mesh), ionode_id ) 
        CALL mp_bcast( rho_atc(1:mesh), ionode_id ) 
        CALL mp_bcast( vloc(1:mesh), ionode_id ) 
        CALL mp_bcast( lll(1:nbeta), ionode_id ) 
        CALL mp_bcast( jjj(1:nbeta), ionode_id ) 
        CALL mp_bcast( kkbeta(1:nbeta), ionode_id ) 
        CALL mp_bcast( beta(1:mesh,1:nbeta), ionode_id ) 
        CALL mp_bcast( dion(1:nbeta,1:nbeta), ionode_id ) 
        CALL mp_bcast( rinner(1:nqlc), ionode_id ) 
        CALL mp_bcast( qqq(1:nbeta,1:nbeta), ionode_id ) 
        CALL mp_bcast( qfunc(1:mesh, 1:nbeta, 1:nbeta), ionode_id )
        CALL mp_bcast( qfcoef(1:nqf, 1:nqlc, 1:nbeta, 1:nbeta), ionode_id ) 
        CALL mp_bcast( chi(1:mesh, 1:nwfc), ionode_id ) 
        CALL mp_bcast( rho_at(1:mesh), ionode_id )
!
      RETURN
    END SUBROUTINE read_restart_pseudo3

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_pseudo2( iuni )

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'pseudo'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum 
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_pseudo, pseudo not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_pseudo2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!
! ..  This subroutine write to disk variables related to the reciprocal 
! ..  space mesh
! .. Where:
! iuni        = Restart file I/O fortran unit
! ng          = number of g vectors
! bi1, bi2, bi3  = initial reciprocal space base vectors (used to determine ng)
! b1, b2, b3     = actual reciprocal space base vectors, to be used to determine
!                  the square modulus of G-vectors 
! mill        = miller index of the G-vectors 
!               Gx(i) = mill(1,i)*b1(1)+mill(2,i)*b2(1)+mill(3,i)*b3(1)
!               Gy(i) = mill(1,i)*b1(2)+mill(2,i)*b2(2)+mill(3,i)*b3(2)
!               Gz(i) = mill(1,i)*b1(3)+mill(2,i)*b2(3)+mill(3,i)*b3(3)
!
    SUBROUTINE write_restart_gvec1(iuni, &
      ng, bi1, bi2, bi3, b1, b2, b3, tmill, mill )

      USE io_global, ONLY: ionode

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ng
      REAL(DP), INTENT(IN) :: bi1(3), bi2(3), bi3(3)
      REAL(DP), INTENT(IN) :: b1(3), b2(3), b3(3)
      INTEGER, INTENT(IN) :: mill(:,:)
      LOGICAL, INTENT(IN) :: tmill
      INTEGER :: idum = 0
      INTEGER :: i, j
      CHARACTER(LEN=30) :: sub_name = ' write_restart_gvec '
      CHARACTER(LEN=20) :: section_name = 'gvec'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!
        IF( tmill ) THEN
          IF( ( SIZE(mill,1) < 3 ) .OR. ( SIZE(mill,2) < ng) ) &
            CALL errore( sub_name, ' wrong size ', 1 )
        END IF

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) ng, tmill
          WRITE(iuni) bi1, bi2, bi3, b1, b2, b3
          IF( tmill ) THEN
            WRITE(iuni) ((mill(i,j),i=1,3),j=1,ng)
          ELSE
            WRITE(iuni) idum
          END IF
        END IF

      RETURN
    END SUBROUTINE write_restart_gvec1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_gvec2( iuni )
 
      USE io_global, ONLY: ionode

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'gvec'

      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum 
      END IF

      RETURN
    END SUBROUTINE write_restart_gvec2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!


    SUBROUTINE read_restart_gvec1(iuni, &
      ng, bi1, bi2, bi3, b1, b2, b3, tmill, mill )

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

! .. Subroutine output:
!    if tmill is true "mill" array is read from file
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tmill
      INTEGER,   INTENT(OUT) :: ng
      REAL(DP), INTENT(OUT) :: b1(3), b2(3), b3(3)
      REAL(DP), INTENT(OUT) :: bi1(3), bi2(3), bi3(3)
      INTEGER,   INTENT(OUT) :: mill(:,:)

      INTEGER :: i, j
      LOGICAL :: twrite_, tmill_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_gvec '
      CHARACTER(LEN=20) :: section_name = 'gvec'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(sub_name, ' Data Section not present in restart file ', 1)

        IF( ionode ) THEN
          READ(iuni) ng, tmill_
          READ(iuni) bi1, bi2, bi3, b1, b2, b3
        END IF
        CALL mp_bcast( ng, ionode_id )
        CALL mp_bcast( tmill_, ionode_id )
        CALL mp_bcast( bi1, ionode_id )
        CALL mp_bcast( bi2, ionode_id )
        CALL mp_bcast( bi3, ionode_id )
        CALL mp_bcast( b1, ionode_id )
        CALL mp_bcast( b2, ionode_id )
        CALL mp_bcast( b3, ionode_id )

        IF( tmill .AND. .NOT. tmill_ ) &
          CALL errore(sub_name, ' mill indexes not present in restart file ', 1)

        IF( tmill ) THEN

          IF( ( SIZE( mill, 2) < ng ) .OR. ( SIZE( mill, 1 ) < 3 ) ) &
            CALL errore(sub_name, ' mill array too small ', 1)

          IF( ionode ) THEN
            READ(iuni) ( ( mill( i, j ), i = 1, 3 ), j = 1, ng )
          END IF
          CALL mp_bcast( mill, ionode_id )

        ELSE

          IF( ionode ) THEN
            READ(iuni) idum
          END IF

        END IF
  
      RETURN
    END SUBROUTINE read_restart_gvec1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_gvec2( iuni )

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'gvec'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_gvec, data not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_gvec2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variables related to a single k point
! .. Where:
! iuni    = Restart file I/O fortran unit
! ngwk    = number of wavefunctions G vectors, for this k point 
! igk(.)  = for each G+k, igk is the index of the G vectors
! xk(.)   = k point coordinate
! wk      = k point weight

    SUBROUTINE write_restart_gkvec1(iuni, ik, nk, ngwk, xk, wk, isk)

      USE io_global, ONLY: ionode

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      ! ... INTEGER, INTENT(IN) :: igk(:)
      INTEGER, INTENT(IN) :: ik, nk, ngwk, isk
      REAL(DP), INTENT(IN) :: xk(3)
      REAL(DP), INTENT(IN) :: wk
      INTEGER :: i, idum = 0
      CHARACTER(LEN=20) :: section_name = 'gkvec'

      INTEGER :: ised( 4 )

      LOGICAL :: twrite = .TRUE.

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) ik, nk, ngwk, ised(1:4), isk
          WRITE(iuni) (xk(i),i=1,3), wk
          WRITE(iuni) idum ! (igk(i),i=1,ngwk)
        END IF

      RETURN
    END SUBROUTINE write_restart_gkvec1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_gkvec2(iuni)
      USE io_global, ONLY: ionode
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'gkvec'
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) idum
          WRITE(iuni) idum
          WRITE(iuni) idum
        END IF
      RETURN
    END SUBROUTINE write_restart_gkvec2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_gkvec1(iuni, &
      ik, nk, ngwk, xk, wk, isk )

!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      ! ... INTEGER, INTENT(INOUT) :: igk(:)
      INTEGER,   INTENT(OUT) :: ngwk, ik, nk, isk
      REAL(DP), INTENT(OUT) :: xk(3)
      REAL(DP), INTENT(OUT) :: wk

      INTEGER :: ised( 4 )

      INTEGER :: i

      INTEGER :: idum, nigk
      LOGICAL :: twrite_
      INTEGER :: ierr
      CHARACTER(LEN=20) :: section_name = 'gkvec'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_gkvec ',' Data Section not present in restart file ', 1)

        IF( ionode ) THEN
          READ(iuni) ik, nk, ngwk, ised( 1 : 4 ), isk
          READ(iuni) ( xk ( i ), i = 1, 3 ), wk
        END IF
        CALL mp_bcast( ngwk, ionode_id )
        CALL mp_bcast( ised, ionode_id )
        CALL mp_bcast( ik, ionode_id )
        CALL mp_bcast( isk, ionode_id )
        CALL mp_bcast( nk, ionode_id )
        CALL mp_bcast( xk, ionode_id )
        CALL mp_bcast( wk, ionode_id )

        IF( ionode ) THEN
          READ(iuni) idum ! .. (igk_(i),i=1,ngwk_)
        END IF
        ! .. CALL mp_bcast( igk_, ionode_id )

      RETURN
    END SUBROUTINE read_restart_gkvec1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_gkvec2(iuni)

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'gkvec'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF

      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_gkvec, xdim not read from restart ' )")

      RETURN
    END SUBROUTINE read_restart_gkvec2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variable related to the simulation cell 
! .. Where:
! iuni    = Restart file I/O fortran unit
! ibrav   = index of the bravais lattice
! celldm  = starting values used to generate the crystal
! ht0     = cell parameters at simulation time "t"
! htm     = cell parameters at simulation time "t-dt"
! htm2    = cell parameters at simulation time "t-2*dt"
! xnosp   = nose thermostat variable at simulation time "t+dt"
! xnos0   = nose thermostat variable at simulation time "t"
! xnosm   = nose thermostat variable at simulation time "t-dt"
! xnosm2  = nose thermostat variable at simulation time "t-2*dt"
!
    SUBROUTINE write_restart_cell1( iuni, &
      ibrav, celldm, ht0, htm, htm2, htvel, xnosp, xnos0, xnosm, xnosm2)
      USE io_global, ONLY: ionode
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ibrav
      REAL(DP), INTENT(IN) :: celldm(6)
      REAL(DP), INTENT(IN) :: ht0(3,3)
      REAL(DP), INTENT(IN) :: htm(3,3)
      REAL(DP), INTENT(IN) :: htm2(3,3)
      REAL(DP), INTENT(IN) :: htvel(3,3)
      REAL(DP), INTENT(IN) :: xnosp(3,3)
      REAL(DP), INTENT(IN) :: xnos0(3,3)
      REAL(DP), INTENT(IN) :: xnosm(3,3)
      REAL(DP), INTENT(IN) :: xnosm2(3,3)
      INTEGER :: i
      CHARACTER(LEN=20) :: section_name = 'cell'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) ibrav, (celldm(i), i=1,6)
          WRITE(iuni) ht0, htm, htm2, htvel
          WRITE(iuni) xnosp, xnos0, xnosm, xnosm2
        END IF

      RETURN
    END SUBROUTINE write_restart_cell1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_cell2(iuni) 
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'cell'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_cell2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_cell1( iuni, &
      ibrav, celldm, ht0, htm, htm2, htvel, xnosp, xnos0, xnosm, xnosm2)

!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      INTEGER,   INTENT(OUT) :: ibrav
      REAL(DP), INTENT(OUT) :: celldm(6)
      REAL(DP), INTENT(OUT) :: ht0(3,3)
      REAL(DP), INTENT(OUT) :: htm(3,3)
      REAL(DP), INTENT(OUT) :: htm2(3,3)
      REAL(DP), INTENT(OUT) :: htvel(3,3)
      REAL(DP), INTENT(OUT) :: xnosp(3,3)
      REAL(DP), INTENT(OUT) :: xnos0(3,3)
      REAL(DP), INTENT(OUT) :: xnosm(3,3)
      REAL(DP), INTENT(OUT) :: xnosm2(3,3)

      INTEGER :: i
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'cell'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_cell ', ' Section not present in restart file ', 1)

        IF( ionode ) THEN
          READ(iuni) ibrav, ( celldm(i), i = 1, 6 )
          READ(iuni) ht0, htm, htm2, htvel
          READ(iuni) xnosp, xnos0, xnosm, xnosm2
        END IF
        CALL mp_bcast( ibrav, ionode_id )
        CALL mp_bcast( celldm, ionode_id )
        CALL mp_bcast( ht0, ionode_id )
        CALL mp_bcast( htm, ionode_id )
        CALL mp_bcast( htm2, ionode_id )
        CALL mp_bcast( htvel, ionode_id )
        CALL mp_bcast( xnosp, ionode_id )
        CALL mp_bcast( xnos0, ionode_id )
        CALL mp_bcast( xnosm, ionode_id )
        CALL mp_bcast( xnosm2, ionode_id )

      RETURN
    END SUBROUTINE read_restart_cell1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_cell2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'cell'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_cell, xdim not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_cell2




!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variable related to the ion types 
! ..  positions, and velocities
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_ions1(iuni, &
      label, tscal, stau0, svel0, staum, svelm, taui, fion, &
      cdmi, nat, ntyp, ityp, na, mass, xnosp, xnos0, xnosm, xnosm2)
!
      USE io_global, ONLY: ionode
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tscal
      CHARACTER(LEN=*), INTENT(IN) :: label(:)
      REAL(DP), INTENT(IN) :: stau0(:,:)
      REAL(DP), INTENT(IN) :: svel0(:,:)
      REAL(DP), INTENT(IN) :: staum(:,:)
      REAL(DP), INTENT(IN) :: svelm(:,:)
      REAL(DP), INTENT(IN) :: taui(:,:)
      REAL(DP), INTENT(IN) :: fion(:,:)
      REAL(DP), INTENT(IN) :: cdmi(:)
      INTEGER, INTENT(IN) :: nat
      INTEGER, INTENT(IN) :: ntyp
      INTEGER, INTENT(IN) :: ityp(:)
      INTEGER, INTENT(IN) :: na(:)
      REAL(DP), INTENT(IN) :: mass(:)
      REAL(DP), INTENT(IN) :: xnosp
      REAL(DP), INTENT(IN) :: xnos0
      REAL(DP), INTENT(IN) :: xnosm
      REAL(DP), INTENT(IN) :: xnosm2

      INTEGER :: i,j
      CHARACTER(LEN=4) :: label_(ntyp)
      CHARACTER(LEN=30) :: sub_name = ' write_restart_ions '
      CHARACTER(LEN=20) :: section_name = 'ions'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

        IF( SIZE( label ) < ntyp ) &
          CALL errore( sub_name, ' wrong size ', 1 )
        IF( SIZE( ityp ) < nat ) &
          CALL errore( sub_name, ' wrong size ', 2 )
        IF( SIZE( na ) < ntyp ) &
          CALL errore( sub_name, ' wrong size ', 3 )
        IF( SIZE( mass ) < ntyp ) &
          CALL errore( sub_name, ' wrong size ', 4 )
        IF( ( SIZE( stau0, 1 ) < 3 ) .OR. ( SIZE( stau0, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 5 )
        IF( ( SIZE( svel0, 1 ) < 3 ) .OR. ( SIZE( svel0, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 6 )
        IF( ( SIZE( staum, 1 ) < 3 ) .OR. ( SIZE( staum, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 7 )
        IF( ( SIZE( svelm, 1 ) < 3 ) .OR. ( SIZE( svelm, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 8 )
        IF( ( SIZE( taui, 1 ) < 3 ) .OR. ( SIZE( taui, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 9 )
        IF( ( SIZE( fion, 1 ) < 3 ) .OR. ( SIZE( fion, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 10 )
        IF( SIZE( cdmi ) < 3 ) &
          CALL errore( sub_name, ' wrong size ', 11 )

        label_ = label(1:ntyp)

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version, section_name
          WRITE(iuni) nat, ntyp, tscal
          WRITE(iuni) (ityp(i),i=1,nat), (na(i),i=1,ntyp), (label_(i),i=1,ntyp)
          WRITE(iuni) (mass(i),i=1,ntyp)
          WRITE(iuni) ((stau0(i,j),i=1,3),j=1,nat)
          WRITE(iuni) ((svel0(i,j),i=1,3),j=1,nat)
          WRITE(iuni) ((staum(i,j),i=1,3),j=1,nat)
          WRITE(iuni) ((svelm(i,j),i=1,3),j=1,nat)
          WRITE(iuni) ((taui(i,j),i=1,3),j=1,nat)
          WRITE(iuni) ((fion(i,j),i=1,3),j=1,nat)
          WRITE(iuni) (cdmi(i),i=1,3)
          WRITE(iuni) xnosp, xnos0, xnosm, xnosm2
        END IF

      RETURN
    END SUBROUTINE write_restart_ions1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_ions2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'ions'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_ions2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_ions1(iuni, &
      label, tscal, stau0, svel0, staum, svelm, taui, fion, &
      cdmi, nat, ntyp, ityp, na, mass, xnosp, xnos0, xnosm, xnosm2)
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(OUT) :: tscal
      CHARACTER(LEN=*), INTENT(OUT) :: label(:)
      REAL(DP), INTENT(OUT) :: stau0(:,:)
      REAL(DP), INTENT(OUT) :: svel0(:,:)
      REAL(DP), INTENT(OUT) :: staum(:,:)
      REAL(DP), INTENT(OUT) :: svelm(:,:)
      REAL(DP), INTENT(OUT) :: taui(:,:)
      REAL(DP), INTENT(OUT) :: fion(:,:)
      REAL(DP), INTENT(OUT) :: cdmi(:)
      INTEGER, INTENT(OUT) :: nat
      INTEGER, INTENT(OUT) :: ntyp
      INTEGER, INTENT(OUT) :: ityp(:)
      INTEGER, INTENT(OUT) :: na(:)
      REAL(DP), INTENT(OUT) :: mass(:)
      REAL(DP), INTENT(OUT) :: xnosp
      REAL(DP), INTENT(OUT) :: xnos0
      REAL(DP), INTENT(OUT) :: xnosm
      REAL(DP), INTENT(OUT) :: xnosm2

      INTEGER   :: i, j
      CHARACTER(LEN=4) :: label_( nsx )

      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_ions '
      CHARACTER(LEN=20) :: section_name = 'ions'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_ions ',' Data Section not present in restart file ', 1)

        IF( ionode ) THEN
          READ(iuni) nat, ntyp, tscal
        END IF
        CALL mp_bcast(nat, ionode_id)
        CALL mp_bcast(ntyp, ionode_id)
        CALL mp_bcast(tscal, ionode_id)

        IF( ntyp > SIZE( na ) ) &
            CALL errore( ' read_restart_ions ', ' too many types ', ntyp )
        IF( nat > SIZE( ityp ) ) &
          CALL errore( ' read_restart_ions ', ' too many atoms ', nat )
        IF( ( SIZE( label ) < ntyp ) .OR. ( SIZE( label_ ) < ntyp ) ) &
          CALL errore( sub_name, ' wrong size for label ', 1 )
        IF( SIZE( mass ) < ntyp ) &
          CALL errore( sub_name, ' wrong size for mass ', 4 )

        IF( ionode ) THEN
          READ(iuni) ( ityp(i), i = 1, nat  ), ( na(i), i = 1, ntyp ), ( label_(i), i = 1, ntyp )
          READ(iuni) ( mass(i), i = 1, ntyp )
        END IF
        CALL mp_bcast( ityp  , ionode_id )
        CALL mp_bcast( na    , ionode_id )
        CALL mp_bcast( label_ , ionode_id )
        CALL mp_bcast( mass  , ionode_id )

        label( 1 : ntyp ) = label_( 1 : ntyp )

        IF( ( SIZE( stau0, 1 ) < 3 ) .OR. ( SIZE( stau0, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 5 )
        IF( ( SIZE( svel0, 1 ) < 3 ) .OR. ( SIZE( svel0, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 6 )
        IF( ( SIZE( staum, 1 ) < 3 ) .OR. ( SIZE( staum, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 7 )
        IF( ( SIZE( svelm, 1 ) < 3 ) .OR. ( SIZE( svelm, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 8 )
        IF( ( SIZE( taui, 1 ) < 3 ) .OR. ( SIZE( taui, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 9 )
        IF( ( SIZE( fion, 1 ) < 3 ) .OR. ( SIZE( fion, 2 ) < nat ) ) &
          CALL errore( sub_name, ' wrong size ', 10 )
        IF( SIZE( cdmi ) < 3 ) &
          CALL errore( sub_name, ' wrong size ', 11 )

        IF( ionode ) THEN
          READ(iuni) ( ( stau0(i,j), i = 1, 3 ), j = 1, nat )
          READ(iuni) ( ( svel0(i,j), i = 1, 3 ), j = 1, nat )
          READ(iuni) ( ( staum(i,j), i = 1, 3 ), j = 1, nat )
          READ(iuni) ( ( svelm(i,j), i = 1, 3 ), j = 1, nat )
          READ(iuni) ( ( taui(i,j), i = 1, 3 ), j = 1, nat )
          READ(iuni) ( ( fion(i,j), i = 1, 3 ), j = 1, nat )
        END IF
        CALL mp_bcast(stau0, ionode_id)
        CALL mp_bcast(svel0, ionode_id)
        CALL mp_bcast(staum, ionode_id)
        CALL mp_bcast(svelm, ionode_id)
        CALL mp_bcast(taui, ionode_id)
        CALL mp_bcast(fion, ionode_id)

        IF( ionode ) THEN
          READ(iuni) ( cdmi(i), i = 1, 3 )
          READ(iuni) xnosp, xnos0, xnosm, xnosm2
        END IF
        CALL mp_bcast(cdmi, ionode_id)
        CALL mp_bcast(xnosp, ionode_id)
        CALL mp_bcast(xnos0, ionode_id)
        CALL mp_bcast(xnosm, ionode_id)
        CALL mp_bcast(xnosm2, ionode_id)
  
      RETURN
    END SUBROUTINE read_restart_ions1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_ions2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'ions'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_ions, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_ions2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variable related to electronic band
! ..  structure (NOT the wavefunctions)
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_electrons1( iuni, &
      occ, occm, tocc, lambda, lambdam, ldim, tlam, nbnd, ispin, nspin, ik, nk, nel, nelu, &
      neld, xnosp, xnos0, xnosm, xnosm2, ef, teig, eig, weig)
!
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      REAL(DP), INTENT(IN) :: occ(:)
      REAL(DP), INTENT(IN) :: occm(:)
      REAL(DP), INTENT(IN) :: lambda(:,:)
      REAL(DP), INTENT(IN) :: lambdam(:,:)
      REAL(DP), INTENT(IN) :: eig(:)
      REAL(DP), INTENT(IN) :: weig(:)
      LOGICAL, INTENT(IN) :: tocc, tlam, teig
      INTEGER, INTENT(IN) :: nbnd, ldim
      INTEGER, INTENT(IN) :: ispin
      INTEGER, INTENT(IN) :: nspin
      INTEGER, INTENT(IN) :: ik
      INTEGER, INTENT(IN) :: nk
      REAL(DP), INTENT(IN) :: nel
      INTEGER, INTENT(IN) :: nelu
      INTEGER, INTENT(IN) :: neld
      REAL(DP), INTENT(IN) :: xnosp
      REAL(DP), INTENT(IN) :: xnos0
      REAL(DP), INTENT(IN) :: xnosm
      REAL(DP), INTENT(IN) :: xnosm2
      REAL(DP), INTENT(IN) :: ef

      INTEGER :: i, l, idum = 0
      CHARACTER(LEN=30) :: sub_name = ' write_restart_electrons '
      CHARACTER(LEN=20) :: section_name = 'electrons'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

        IF( ionode ) WRITE(iuni) twrite, file_version, section_name
        IF( ionode ) WRITE(iuni) nbnd, ispin, nspin, ik, nk, nel, nelu, neld, ldim
        IF( ionode ) WRITE(iuni) tocc

        IF( tocc ) THEN
         
          IF( SIZE( occ ) < nbnd ) &
            CALL errore(sub_name, ' wrong size ', 1 )
          IF( SIZE( occm ) < nbnd ) &
            CALL errore(sub_name, ' wrong size ', 2 )
          
          IF( ionode ) WRITE(iuni) (occ(i),i=1,nbnd)
          IF( ionode ) WRITE(iuni) (occm(i),i=1,nbnd)

        ELSE

          IF( ionode ) WRITE(iuni) idum
          IF( ionode ) WRITE(iuni) idum

        END IF

        IF( ionode ) WRITE(iuni) tlam

        IF( tlam ) THEN

          IF( ( SIZE( lambda, 1 ) < ldim ) .OR. ( SIZE( lambda, 2 ) < ldim ) ) &
            CALL errore(sub_name, ' wrong size ', 3 )
          IF( ( SIZE( lambdam, 1 ) < ldim ) .OR. ( SIZE( lambdam, 2 ) < ldim ) ) &
            CALL errore(sub_name, ' wrong size ', 4 )

          IF( ionode ) WRITE(iuni) ((lambda(l,i),l=1,ldim),i=1,ldim)
          IF( ionode ) WRITE(iuni) ((lambdam(l,i),l=1,ldim),i=1,ldim)

        ELSE

          IF( ionode ) WRITE(iuni) idum
          IF( ionode ) WRITE(iuni) idum

        END IF

        IF( ionode ) WRITE(iuni) xnosp, xnos0, xnosm, xnosm2
        IF( ionode ) WRITE(iuni) ef

        IF( ionode ) WRITE(iuni) teig

        IF( teig ) THEN

          IF( SIZE( eig ) < nbnd ) &
            CALL errore(sub_name, ' wrong size ', 5 )
          IF( SIZE( weig ) < nbnd ) &
            CALL errore(sub_name, ' wrong size ', 6 )

          IF( ionode ) WRITE(iuni) (eig(i),i=1,nbnd)
          IF( ionode ) WRITE(iuni) (weig(i),i=1,nbnd)

        ELSE

          IF( ionode ) WRITE(iuni) idum
          IF( ionode ) WRITE(iuni) idum

        END IF

      RETURN
    END SUBROUTINE write_restart_electrons1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

   SUBROUTINE write_restart_electrons2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'electrons'
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_electrons2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_electrons1( iuni, &
      occ, occm, tocc, lambda, lambdam, ldim, tlam, nbnd, ispin, nspin, ik, nk, nel, nelu, &
      neld, xnosp, xnos0, xnosm, xnosm2, ef, teig, eig, weig)

! .. Subroutine output:
!    if tocc is true then "occ" and "occm" are overwritten
!    if tlam is true then "lambda" and "lambdam" are overwritten
!    if teig is true then "eig" and "weig" are overwritten
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      REAL(DP), INTENT(OUT) :: occ(:)
      REAL(DP), INTENT(OUT) :: occm(:)
      REAL(DP), INTENT(OUT) :: eig(:)
      REAL(DP), INTENT(OUT) :: weig(:)
      REAL(DP), INTENT(OUT) :: lambda(:,:)
      REAL(DP), INTENT(OUT) :: lambdam(:,:)
      INTEGER, INTENT(OUT) :: ldim
      INTEGER, INTENT(OUT) :: nbnd
      INTEGER, INTENT(OUT) :: ispin
      INTEGER, INTENT(OUT) :: nspin
      INTEGER, INTENT(OUT) :: ik
      INTEGER, INTENT(OUT) :: nk
      REAL(DP), INTENT(OUT) :: nel
      INTEGER, INTENT(OUT) :: nelu
      INTEGER, INTENT(OUT) :: neld
      REAL(DP), INTENT(OUT) :: xnosp
      REAL(DP), INTENT(OUT) :: xnos0
      REAL(DP), INTENT(OUT) :: xnosm
      REAL(DP), INTENT(OUT) :: xnosm2
      LOGICAL, INTENT(IN) :: tocc, tlam, teig
      REAL(DP), INTENT(OUT) :: ef

      INTEGER :: i, j, k, l
      LOGICAL :: tocc_, tlam_, teig_
      INTEGER :: idum

      LOGICAL :: twrite_
      INTEGER :: ierr
      CHARACTER(LEN=30) :: sub_name = ' read_restart_electrons '
      CHARACTER(LEN=20) :: section_name = 'electrons'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!
      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_electrons ',' Data Section not present in restart file ', 1)

        IF( ionode ) READ(iuni) nbnd, ispin, nspin, ik, nk, nel, nelu, neld, ldim

        CALL mp_bcast( nbnd, ionode_id)
        CALL mp_bcast( ispin, ionode_id)
        CALL mp_bcast( nspin, ionode_id)
        CALL mp_bcast( ik, ionode_id)
        CALL mp_bcast( nk, ionode_id)
        CALL mp_bcast( nel, ionode_id)
        CALL mp_bcast( nelu, ionode_id)
        CALL mp_bcast( neld, ionode_id)
        CALL mp_bcast( ldim, ionode_id)

!
! ..    Manage occ and occm

        IF( ionode ) READ(iuni) tocc_
        CALL mp_bcast( tocc_, ionode_id)

        IF( tocc .AND. .NOT. tocc_ ) &
          CALL errore( ' read_restart_electrons ',' occupation number not present in restart ', 1)

        IF( tocc ) THEN

          IF( nbnd >  SIZE( occ ) ) &
            CALL errore( ' read_restart_electrons ',' wrong dimensions for occ ', 1)
          IF( nbnd >  SIZE( occm ) ) &
            CALL errore( ' read_restart_electrons ',' wrong dimensions for occm ', 1)

          IF( ionode ) READ(iuni) ( occ(i), i = 1, nbnd )
          CALL mp_bcast( occ, ionode_id )

          IF( ionode ) READ(iuni) ( occm(i), i = 1, nbnd )
          CALL mp_bcast( occm, ionode_id )

        ELSE

          IF( ionode ) READ(iuni) idum
          IF( ionode ) READ(iuni) idum

        END IF

! ..    Manage lambda and lambdam

        IF( ionode ) READ(iuni) tlam_
        CALL mp_bcast( tlam_, ionode_id)

        IF( tlam .AND. .NOT. tlam_ ) &
          CALL errore( ' read_restart_electrons ',' lambda matrix not present in restart ', 1)

        IF( tlam ) THEN

          IF( ldim > SIZE( lambda, 1 ) .OR. ldim > SIZE( lambda, 2 ) ) &
            CALL errore( ' read_restart_electrons ',' wrong dimensions for lambda ', 1)
          IF( ldim > SIZE( lambdam, 1 ) .OR. ldim > SIZE( lambdam, 2 ) ) &
            CALL errore( ' read_restart_electrons ',' wrong dimensions for lambdam ', 1)

          IF( ionode ) READ(iuni) ( ( lambda(l,i), l = 1, ldim ), i = 1, ldim )
          CALL mp_bcast( lambda, ionode_id)

          IF( ionode ) READ(iuni) ( ( lambdam(l,i), l = 1, ldim ), i = 1, ldim )
          CALL mp_bcast( lambdam, ionode_id)

        ELSE

          IF( ionode ) READ(iuni) idum
          IF( ionode ) READ(iuni) idum

        END IF

        IF( ionode ) READ(iuni) xnosp, xnos0, xnosm, xnosm2
        IF( ionode ) READ(iuni) ef
        CALL mp_bcast( xnosp, ionode_id)
        CALL mp_bcast( xnos0, ionode_id)
        CALL mp_bcast( xnosm, ionode_id)
        CALL mp_bcast( xnosm2, ionode_id)
        CALL mp_bcast( ef, ionode_id )

        IF( ionode ) READ(iuni) teig_
        CALL mp_bcast( teig_, ionode_id)

        IF( teig .AND. .NOT. teig_ ) &
          CALL errore( ' read_restart_electrons ',' occupation number not present in restart ', 1)

        IF( teig ) THEN

          IF( nbnd > SIZE( eig ) ) &
            CALL errore( ' read_restart_electrons ',' wrong dimensions for eig ', 1)
          IF( nbnd > SIZE( weig ) ) &
            CALL errore( ' read_restart_electrons ',' wrong dimensions for weig ', 1)

          IF( ionode ) READ(iuni) ( eig(i), i = 1, nbnd )
          CALL mp_bcast( eig, ionode_id )

          IF( ionode ) READ(iuni) ( weig(i), i = 1, nbnd )
          CALL mp_bcast( weig, ionode_id )

        ELSE

          IF( ionode ) READ(iuni) idum
          IF( ionode ) READ(iuni) idum

        END IF
  
      RETURN
    END SUBROUTINE read_restart_electrons1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_electrons2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum
      CHARACTER(LEN=20) :: section_name = 'electrons'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_electrons, Data Sections not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_electrons2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write wavefunctions to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_wfc1(iuni, &
      ik, nk, kunit, ispin, nspin, scal, wf0, t0, wfm, tm, ngw, nbnd, igl, ngwl )
!
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_get, mp_bcast, mp_max
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, my_image_id
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP), INTENT(IN) :: wf0(:,:)
      COMPLEX(DP), INTENT(IN) :: wfm(:,:)
      INTEGER, INTENT(IN) :: ngw   ! 
      INTEGER, INTENT(IN) :: nbnd
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)
      REAL(DP), INTENT(IN) :: scal
      LOGICAL, INTENT(IN) :: t0, tm

      INTEGER :: i, j, ierr, idum = 0
      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx
      INTEGER :: npool, ipmask( nproc ), ipsour
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)

      CHARACTER(LEN=20) :: section_name = 'wfc'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

        IF( ionode ) WRITE(iuni) twrite, file_version, section_name

        ! set working variables for k point index (ikt) and k points number (nkt)
        ikt = ik
        nkt = nk

        !  find out the number of pools
        npool = nproc / nproc_pool 

        !  find out number of k points blocks
        nkbl = nkt / kunit  

        !  k points per pool
        nkl = kunit * ( nkbl / npool )

        !  find out the reminder
        nkr = ( nkt - nkl * npool ) / kunit

        !  Assign the reminder to the first nkr pools
        IF( my_pool_id < nkr ) nkl = nkl + kunit

        !  find out the index of the first k point in this pool
        iks = nkl * my_pool_id + 1
        IF( my_pool_id >= nkr ) iks = iks + nkr * kunit
      
        !  find out the index of the last k point in this pool
        ike = iks + nkl - 1

        ipmask = 0
        ipsour = ionode_id

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool > 1 ) THEN
          IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
          END IF
          CALL mp_sum( ipmask )
          DO i = 1, nproc
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
          END DO
        END IF

        igwx = 0
        ierr = 0
        IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
          IF( ngwl > SIZE( igl ) ) THEN
            ierr = 1
          ELSE
            igwx = MAXVAL( igl(1:ngwl) )
          END IF
        END IF

        ! get the maximum index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm ) 

        ! now notify all procs if an error has been found 
        !
        CALL mp_max( ierr ) 

        IF( ierr > 0 ) &
          CALL errore(' write_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipsour /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipsour, 1 )
        END IF

        IF( ionode ) WRITE(iuni) ngw, nbnd, ik, nk, kunit, ispin, nspin, scal
        IF( ionode ) WRITE(iuni) igwx

        ! write(200+mpime+ik*10,*) mpime, nproc, root, me_pool, my_pool_id, nproc_pool, intra_pool_comm, root_pool, npool
        ! write(200+mpime+ik*10,*) ngwl, nkbl, kunit, iks, ike, ngw, nbnd, ik, nk, kunit, ispin, nspin, scal, igwx, ierr
        ! close(200+mpime+ik*10)

        ALLOCATE( wtmp( MAX(igwx,1) ) )
        wtmp = 0.0d0

        IF( ionode ) WRITE(iuni) t0

        DO j = 1, nbnd
          IF( t0 ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
                CALL mergewf(wf0(:,j), wtmp, ngwl, igl, me_pool, nproc_pool, root_pool, intra_pool_comm)
              END IF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j )
              END IF
            ELSE
              CALL mergewf(wf0(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id)
            END IF
            IF( ionode ) WRITE(iuni) ( wtmp(i), i=1,igwx )
          ELSE
            IF( ionode ) WRITE(iuni) j
          END IF
        END DO

        IF( ionode ) WRITE(iuni) tm

        DO j = 1, nbnd
          IF( tm ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
                CALL mergewf(wfm(:,j), wtmp, ngwl, igl, me_pool, nproc_pool, root_pool, intra_pool_comm)
              END IF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j )
              END IF
            ELSE
              CALL mergewf(wfm(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id)
            END IF
            IF( ionode ) WRITE(iuni) (wtmp(i),i=1,igwx)
          ELSE
            IF( ionode ) WRITE(iuni) j
          END IF
        END DO

        DEALLOCATE( wtmp )

      RETURN
    END SUBROUTINE write_restart_wfc1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_wfc2(iuni, nbnd)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni, nbnd
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum, i
      CHARACTER(LEN=20) :: section_name = 'wfc'
      idum = nbnd
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum, idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        DO i = 1, nbnd
          WRITE(iuni) idum
        END DO
        WRITE(iuni) idum
        DO i = 1, nbnd
          WRITE(iuni) idum
        END DO
      END IF
      RETURN
    END SUBROUTINE write_restart_wfc2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_wfc1(iuni, &
      ik, nk, kunit, ispin, nspin, scal, wf0, t0, wfm, tm, ngw, nbnd, igl, ngwl )
!
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_put, mp_bcast, mp_max, mp_get
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, my_image_id
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      COMPLEX(DP), INTENT(INOUT) :: wf0(:,:)
      COMPLEX(DP), INTENT(INOUT) :: wfm(:,:)
      INTEGER, INTENT(IN) :: ik, nk, kunit
      INTEGER, INTENT(OUT) :: ngw, nbnd, ispin, nspin
      REAL(DP), INTENT(OUT) :: scal
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)
      LOGICAL, INTENT(INOUT) :: t0, tm

      INTEGER :: i, j, idum
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)
      LOGICAL :: t0_, tm_
      LOGICAL :: twrite_
      INTEGER :: ierr

      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx, igwx_
      INTEGER :: ik_, nk_, kunit_
      INTEGER :: npool, ipmask( nproc ), ipdest
      CHARACTER(LEN=20) :: section_name = 'wfc'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!
      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_wfc ',' Data Section not present in restart file ', 1)

        IF( ionode ) READ(iuni) ngw, nbnd, ik_, nk_, kunit_, ispin, nspin, scal
        IF( ionode ) READ(iuni) igwx_
        CALL mp_bcast( ngw, ionode_id )
        CALL mp_bcast( nbnd, ionode_id )
        CALL mp_bcast( ik_, ionode_id )
        CALL mp_bcast( nk_, ionode_id )
        CALL mp_bcast( kunit_, ionode_id )
        CALL mp_bcast( ispin, ionode_id )
        CALL mp_bcast( nspin, ionode_id )
        CALL mp_bcast( scal, ionode_id )
        CALL mp_bcast( igwx_, ionode_id )

        ! set working variables for k point index (ikt) and k points number (nkt)
        ikt = ik
        nkt = nk

        !  find out the number of pools
        npool = nproc / nproc_pool

        !  find out number of k points blocks (each block contains kunit k points)
        nkbl = nkt / kunit

        !  k points per pool
        nkl = kunit * ( nkbl / npool )

        !  find out the reminder
        nkr = ( nkt - nkl * npool ) / kunit

        !  Assign the reminder to the first nkr pools
        IF( my_pool_id < nkr ) nkl = nkl + kunit

        !  find out the index of the first k point in this pool
        iks = nkl * my_pool_id + 1
        IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

        !  find out the index of the last k point in this pool
        ike = iks + nkl - 1

        ipmask = 0
        ipdest = ionode_id

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool > 1 ) THEN
          IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
          END IF
          CALL mp_sum( ipmask )
          DO i = 1, nproc
            IF( ipmask(i) == 1 ) ipdest = ( i - 1 )
          END DO
        END IF

        igwx = 0
        ierr = 0
        IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
          IF( ngwl > SIZE( igl ) ) THEN
            ierr = 1
          ELSE
            igwx = MAXVAL( igl(1:ngwl) )
          END IF
        END IF

        ! get the maximum index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm ) 

        ! now notify all procs if an error has been found 
        !
        CALL mp_max( ierr ) 

        IF( ierr > 0 ) &
          CALL errore(' read_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipdest /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipdest, 1 )
        END IF

!
! ...   Here read wave function at time t
!  
        IF( ionode ) READ(iuni) t0_
        CALL mp_bcast( t0_, ionode_id )

        IF( .NOT. ( t0_ .AND. t0 ) .AND. ( restart_module_verbosity > 1000 ) ) &
          WRITE( stdout,fmt="(3X,'W: read_restart_wfc, wf0 not read from restart ' )")

        ! ... WRITE( stdout,*) ' #### ', igwx_, igwx, ngwl, iks, ikt, ike ! DEBUG

        DO j = 1, nbnd

          IF( t0_ .AND. t0 ) THEN

            ALLOCATE( wtmp( MAX(igwx_, igwx) ) )

            IF( ionode ) READ(iuni) ( wtmp(i), i=1,igwx_ )
            IF( igwx > igwx_ ) wtmp( (igwx_ + 1) : igwx ) = 0.0d0

            IF( npool > 1 ) THEN
              IF( ipdest /= ionode_id ) THEN
                CALL mp_put( wtmp, wtmp, mpime, ionode_id, ipdest, j )
              END IF
              IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
                CALL splitwf(wf0(:,j), wtmp, ngwl, igl, me_pool, nproc_pool, root_pool, intra_pool_comm)
              END IF
            ELSE
              CALL splitwf(wf0(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id)
            END IF

            DEALLOCATE( wtmp )

          ELSE

            IF( ionode ) READ(iuni) idum

          END IF

        END DO

!
! ...   Here read wave function at time t-dt
!
        IF( ionode ) READ(iuni) tm_
        CALL mp_bcast( tm_, ionode_id )

        IF( .NOT. ( tm_ .AND. tm ) .AND. ( restart_module_verbosity > 1000 ) ) &
          WRITE( stdout,fmt="(3X,'W: read_restart_wfc, wfm not read from restart ' )")

        DO j = 1, nbnd

          IF( tm_ .AND. tm ) THEN

            ALLOCATE( wtmp( MAX(igwx_, igwx) ) )

            IF( ionode ) READ(iuni) ( wtmp(i), i=1,igwx_ )
            IF( igwx > igwx_ ) wtmp( (igwx_ + 1) : igwx ) = 0.0d0

            IF( npool > 1 ) THEN
              IF( ipdest /= ionode_id ) THEN
                CALL mp_put( wtmp, wtmp, mpime, ionode_id, ipdest, j )
              END IF
              IF( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
                CALL splitwf(wfm(:,j), wtmp, ngwl, igl, me_pool, nproc_pool, root_pool, intra_pool_comm)
              END IF
            ELSE
              CALL splitwf(wfm(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id)
            END IF

            DEALLOCATE( wtmp )

          ELSE

            IF( ionode ) READ(iuni) idum

          END IF

        END DO

! ...   this is to inform the calling subroutine on what has been read
!
        t0   = t0_
        tm   = tm_

      RETURN
    END SUBROUTINE read_restart_wfc1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_wfc2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum, i, nbnd_
      INTEGER :: ierr
      CHARACTER(LEN=20) :: section_name = 'wfc'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum, nbnd_
        READ(iuni) idum
        READ(iuni) idum   ! t0
        DO i = 1, nbnd_
          READ(iuni) idum
        END DO
        READ(iuni) idum   ! t1
        DO i = 1, nbnd_
          READ(iuni) idum
        END DO
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_wfc, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_wfc2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write potential and charge density to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_charge1(iuni, &
      rhog, tr, vg, tv, ng, ispin, nspin, igl, ngl)

      USE mp_wave
      USE mp_global, ONLY: mpime, nproc, root
      USE io_global, ONLY: ionode, ionode_id
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      COMPLEX(DP), INTENT(IN) :: rhog(:)
      COMPLEX(DP), INTENT(IN) :: vg(:)
      INTEGER, INTENT(IN) :: ispin, nspin, ng, ngl, igl(:)
      LOGICAL, INTENT(IN) :: tr, tv
      INTEGER :: i, is
      COMPLEX(DP), ALLOCATABLE :: vtmp(:)
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'charge'

      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!

        ALLOCATE( vtmp (ng) )
        IF( ionode ) WRITE(iuni) twrite, file_version, section_name
        IF( ionode ) WRITE(iuni) ng, ispin, nspin
        IF( ionode ) WRITE(iuni) tr
        IF( tr ) THEN
          CALL mergewf(rhog(:), vtmp, ngl, igl, mpime, nproc, ionode_id)
          IF( ionode ) WRITE(iuni) (vtmp(i),i=1,ng)
        ELSE
          IF( ionode ) WRITE(iuni) idum
        END IF
        IF( ionode ) WRITE(iuni) tv
        IF( tv ) THEN
          CALL mergewf(vg(:), vtmp, ngl, igl, mpime, nproc, ionode_id)
          IF( ionode ) WRITE(iuni) (vtmp(i),i=1,ng)
        ELSE
          IF( ionode ) WRITE(iuni) idum
        END IF
        DEALLOCATE( vtmp )

      RETURN
    END SUBROUTINE write_restart_charge1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_charge2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'charge'
      idum    = 0
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version, section_name
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE write_restart_charge2


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_charge1(iuni, &
      rhog, tr, vg, tv, ng, ispin, nspin, igl, ngl)

      USE mp_wave
      USE mp_global, ONLY: mpime, nproc, root
      USE io_global, ONLY: ionode, ionode_id
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      COMPLEX(DP), INTENT(OUT) :: rhog(:)
      COMPLEX(DP), INTENT(OUT) :: vg(:)
      INTEGER, INTENT(IN) :: ngl, igl(:)
      INTEGER, INTENT(OUT) :: ispin, nspin, ng
      LOGICAL, INTENT(INOUT) :: tr, tv
      INTEGER :: i, j, k, is
      LOGICAL :: tr_, tv_
      COMPLEX(DP), ALLOCATABLE :: vtmp(:)

      LOGICAL :: twrite_
      INTEGER :: idum, ierr
      CHARACTER(LEN=20) :: section_name = 'charge'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore(' read_restart_charge ',' Data Section not present in restart file ', 1)


        IF( ionode ) READ(iuni) ng, ispin, nspin
        CALL mp_bcast( ng, ionode_id )
        CALL mp_bcast( ispin, ionode_id )
        CALL mp_bcast( nspin, ionode_id )

        IF( ionode ) READ(iuni) tr_
        CALL mp_bcast( tr_, ionode_id )

        IF( tr .AND. .NOT. tr_ ) &
          CALL errore(' read_restart_charge ',' rho not present in restart ', 1)  

        IF( tr_ ) THEN
          ALLOCATE( vtmp( ng ) )
          IF( ionode ) READ(iuni) (vtmp(i),i=1,ng)
          CALL splitwf(rhog(:), vtmp, ngl, igl, mpime, nproc, ionode_id)
          DEALLOCATE( vtmp )
        ELSE
          IF( ionode ) READ(iuni) idum
        END IF

        IF( ionode ) READ(iuni) tv_
        CALL mp_bcast( tv_, ionode_id )

        IF( tv .AND. .NOT. tv_ ) &
          CALL errore(' read_restart_charge ',' V not present in restart ', 1)  

        IF( tv_ ) THEN
          ALLOCATE( vtmp( ng ) )
          IF( ionode ) READ(iuni) (vtmp(i),i=1,ng)
          CALL splitwf(vg(:), vtmp, ngl, igl, mpime, nproc, ionode_id)
          DEALLOCATE( vtmp )
        ELSE
          IF( ionode ) READ(iuni) idum
        END IF

      RETURN
    END SUBROUTINE read_restart_charge1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_charge2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum, i, nspin_
      INTEGER :: ierr
      CHARACTER(LEN=20) :: section_name = 'charge'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_charge, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_charge2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! .. This subroutine write information about tetrahedra to the disk
! ..   Where:
!      iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_tetra1( iuni, ltetra, ntetra, tetra )

      USE mp_wave
      USE io_global, ONLY: ionode, ionode_id

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: ltetra
      INTEGER, INTENT(IN) :: ntetra, tetra(:,:)
      INTEGER :: i, j, idum(4,1)
      CHARACTER(LEN=20) :: section_name = 'tetra'
      LOGICAL :: twrite = .TRUE.
!
! ... Subroutine Body
!
      IF( ionode ) WRITE(iuni) twrite, file_version, section_name
      IF( ionode ) WRITE(iuni) ltetra, ntetra
      IF( ltetra ) THEN
        IF( ionode ) WRITE(iuni) ( ( tetra(i,j), i = 1, 4 ) , j = 1, ntetra )
      ELSE
        idum = 0
        IF( ionode ) WRITE(iuni) ( ( idum(i,j), i = 1, 4 ) , j = 1, 1 ) 
      END IF 

      RETURN
    END SUBROUTINE write_restart_tetra1

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

! .. This subroutine write information about tetrahedra to the disk
! ..   Where:
!      iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_tetra2( iuni )

      USE io_global, ONLY: ionode, ionode_id

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      INTEGER :: idum = 0
      CHARACTER(LEN=20) :: section_name = 'tetra'
      LOGICAL :: twrite = .FALSE.
!
! ... Subroutine Body
!
      IF( ionode ) WRITE(iuni) twrite, file_version, section_name
      IF( ionode ) WRITE(iuni) idum
      IF( ionode ) WRITE(iuni) idum

      RETURN
    END SUBROUTINE write_restart_tetra2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!


    SUBROUTINE read_restart_tetra1( iuni, ltetra, ntetra, tetra )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(OUT) :: ltetra
      INTEGER, INTENT(OUT) :: ntetra, tetra(:,:)

      INTEGER :: i, j
      LOGICAL :: twrite_
      INTEGER :: ierr
      INTEGER :: idum(4,1)
      CHARACTER(LEN=30) :: sub_name = ' read_restart_tetra1 '
      CHARACTER(LEN=20) :: section_name = 'tetra'
      CHARACTER(LEN=20) :: section_name_
!
! ... Subroutine Body
!
      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( .NOT. twrite_ ) &
        CALL errore( sub_name , ' Data Section not present in restart file ', 1)

      IF( ionode ) THEN
        READ(iuni) ltetra, ntetra
      END IF
      CALL mp_bcast(ltetra, ionode_id)
      CALL mp_bcast(ntetra, ionode_id)

      IF( ltetra ) THEN
        IF( ionode ) READ(iuni) ( ( tetra(i,j), i = 1, 4 ) , j = 1, ntetra )
        CALL mp_bcast( tetra, ionode_id )
      ELSE
        IF( ionode ) READ(iuni) ( ( idum(i,j), i = 1, 4 ) , j = 1, 1 ) 
      END IF 

      RETURN
    END SUBROUTINE read_restart_tetra1


!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!


    SUBROUTINE read_restart_tetra2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum
      INTEGER :: ierr
      CHARACTER(LEN=20) :: section_name = 'tetra'
      CHARACTER(LEN=20) :: section_name_

      CALL data_section_head( iuni, section_name_ , twrite_ , ierr )

      IF( ionode ) THEN
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE( stdout,fmt="(3X,'W: read_restart_tetra, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE read_restart_tetra2

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE data_section_head( iuni, section_name, twrite, ierr )

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=20), INTENT(OUT) :: section_name
      LOGICAL, INTENT(OUT) :: twrite

      INTEGER :: file_version_
!
      ierr = 0
      IF( ionode ) THEN
        READ(iuni) twrite, file_version_ , section_name
        IF( file_version_ /= file_version ) ierr = 2
      END IF
!
      CALL mp_bcast( ierr, ionode_id )
      IF( ierr == 2 ) &
        CALL errore( ' data_section_head ', ' Restart file versions do not match ', 1)

      CALL mp_bcast( twrite, ionode_id )
      CALL mp_bcast( section_name, ionode_id )

      RETURN
    END SUBROUTINE data_section_head

!=----------------------------------------------------------------------------=!
  END MODULE io_base
!=----------------------------------------------------------------------------=!
