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
! - tread  - logical - true, read effective data; false, read dummy data
! - tovrw  - logical - true, overvrite all input variables with the read data
!
! All Data Sections have a well defined number of records, and the first
! record always contains the following _NON_ _DUMMY_ informations:
! "twrite" "file_version" 
!
!


  USE kinds
  USE parameters
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: file_version = 200
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
! nbeg          = run flags
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
! na(.)         = numbero of atom for each species
! acc(.)        = accomulators
! nacc          = number of accomulators
! ecutwfc       = wave function cutoff
! ecutrho       = charge density cutoff cutoff

    SUBROUTINE write_restart_header1(iuni, twrite, &
      nfi, trutim, nbeg, nr1, nr2, nr3, nr1s, nr2s, nr3s, ng_l, ng_g, nk_l, nk_g, &
      ngwk_l, ngwk_g, nspin, nbnd, nel, nelu, neld, nat, ntyp, na, acc, nacc, &
      ecutwfc, ecutrho, alat, ekinc, kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, &
      ngauss, lgauss, ntetra, ltetra, natomwfc, gcutm, gcuts, dual, doublegrid, &
      modenum, lforce, lstres, title, crystal, tmp_dir, tupf, gamma_only )
!
      USE io_global, ONLY: ionode
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni  
      LOGICAL, INTENT(IN) :: twrite
      INTEGER, INTENT(IN) :: nfi   
      INTEGER, INTENT(IN) :: nbeg  
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      INTEGER, INTENT(IN) :: nr1s, nr2s, nr3s
      INTEGER, INTENT(IN) :: ng_l, ng_g
      REAL(dbl), INTENT(IN) :: trutim  ! true time since last 'from_scratch'
      REAL(dbl), INTENT(IN) :: ecutwfc, ecutrho  ! wfc and density cutoff
      REAL(dbl), INTENT(IN) :: nel
      INTEGER, INTENT(IN) :: nk_l  ! local number of k points (k points in the pool)
      INTEGER, INTENT(IN) :: nk_g  ! global number of k points
      INTEGER, INTENT(IN) :: nspin, nbnd, nelu, neld, nat, ntyp
      INTEGER, INTENT(IN) :: ngwk_l(:), ngwk_g(:)
      INTEGER, INTENT(IN) :: na(:)
      REAL(dbl), INTENT(IN) :: acc(:)
      INTEGER, INTENT(IN) :: nacc
      REAL(dbl), INTENT(IN) :: alat
      REAL(dbl), INTENT(IN) :: ekinc

      INTEGER, INTENT(IN) ::  kunit, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl), INTENT(IN) :: dgauss
      INTEGER, INTENT(IN) :: ngauss
      LOGICAL, INTENT(IN) :: lgauss
      INTEGER, INTENT(IN) :: ntetra
      LOGICAL, INTENT(IN) :: ltetra

      INTEGER, INTENT(IN) :: natomwfc
      LOGICAL, INTENT(IN) :: doublegrid
      REAL(dbl), INTENT(IN) :: gcutm, gcuts, dual
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

      INTEGER :: i
      INTEGER :: idum
      CHARACTER(LEN=80) :: t_, c_, tmp_dir_
      CHARACTER(LEN=30) :: sub_name = ' write_restart_header '

      t_ = title
      c_ = crystal
      tmp_dir_ = tmp_dir

      IF( ntyp > SIZE( na ) ) &
        CALL error(sub_name, ' wrong size ', 1 )
      IF( nk_g > SIZE( ngwk_l ) ) &
        CALL error(sub_name, ' wrong size ', 2 )
      IF( nk_g > SIZE( ngwk_g ) ) &
        CALL error(sub_name, ' wrong size ', 3 )
      IF( nacc > SIZE( acc ) ) &
        CALL error(sub_name, ' wrong size ', 4 )

      IF( ionode ) THEN
        IF( twrite ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) nfi, nbeg, nr1, nr2, nr3, nr1s, nr2s, nr3s, ng_l, ng_g, nk_l, nk_g, &
            nspin, nbnd, &
            nel, nelu, neld, nat, ntyp, nacc, trutim, ecutwfc, ecutrho, alat, ekinc,   &
            kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra,  &
            natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, lforce, tupf, gamma_only 
          WRITE(iuni) (na(i),i=1,ntyp), (ngwk_l(i),i=1,nk_g), (ngwk_g(i),i=1,nk_g), (acc(i),i=1,nacc)
          WRITE(iuni) t_, c_, tmp_dir_
        ELSE
          CALL write_restart_header2(iuni)
        END IF
      END IF

      ! CALL mp_bcast(ios, root, group)
      ! IF( ios /= 0 ) CALL error(' writefile ',' writing integral_time ', ios)

      RETURN
    END SUBROUTINE


    SUBROUTINE write_restart_header2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

! ..  This subroutine read from disk dimensions, and status variables
!
    SUBROUTINE read_restart_header1(iuni, tovrw, tread, nfi, trutim, nbeg, nr1, nr2, nr3, &
      nr1s, nr2s, nr3s, ng_l, ng_g, nk_l, nk_g, ngwk_l, ngwk_g, nspin, nbnd, nel, nelu, neld, &
      nat, ntyp, na, acc, nacc, ecutwfc, ecutrho, alat, ekinc, kunit, &
      k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra, &
      natomwfc, gcutm, gcuts, dual, doublegrid, modenum, &
      lforce, lstres, title, crystal, tmp_dir, tupf, gamma_only )

! .. Subroutine output:
!    if tread is true then on output the following variables are overwritten
!      acc, trutim, nfi, ekinc
!    if tovrw is true all variables are overvritten with values read from file
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      INTEGER, INTENT(INOUT) :: nfi
      INTEGER, INTENT(INOUT) :: nbeg, nr1, nr2, nr3, ng_l, ng_g
      INTEGER, INTENT(INOUT) :: nr1s, nr2s, nr3s
      REAL(dbl), INTENT(OUT) :: trutim
      REAL(dbl), INTENT(INOUT) :: ecutwfc, ecutrho
      REAL(dbl), INTENT(INOUT) :: nel
      INTEGER, INTENT(INOUT) :: nk_l, nk_g, nspin, nbnd, nelu, neld, nat, ntyp
      INTEGER, INTENT(INOUT) :: ngwk_l(:), ngwk_g(:)
      INTEGER, INTENT(INOUT) :: na(:)
      REAL(dbl), INTENT(INOUT) :: acc(:)
      INTEGER, INTENT(INOUT) :: nacc
      REAL(dbl), INTENT(INOUT) :: alat
      REAL(dbl), INTENT(INOUT) :: ekinc
      INTEGER, INTENT(INOUT) ::  kunit, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl), INTENT(INOUT) :: dgauss
      INTEGER, INTENT(INOUT) :: ngauss
      LOGICAL, INTENT(INOUT) :: lgauss
      INTEGER, INTENT(INOUT) :: ntetra
      LOGICAL, INTENT(INOUT) :: ltetra

      INTEGER, INTENT(INOUT) :: natomwfc
      LOGICAL, INTENT(INOUT) :: doublegrid
      REAL(dbl), INTENT(INOUT) :: gcutm, gcuts, dual
      INTEGER, INTENT(INOUT) :: modenum

      LOGICAL, INTENT(INOUT) :: lstres
      LOGICAL, INTENT(INOUT) :: lforce
      CHARACTER(LEN=*), INTENT(INOUT) :: title
      CHARACTER(LEN=*), INTENT(INOUT) :: crystal
      CHARACTER(LEN=*), INTENT(INOUT) :: tmp_dir
      LOGICAL, INTENT(INOUT) :: tupf
      LOGICAL, INTENT(INOUT) :: gamma_only

      INTEGER :: nfi_, nbeg_, nr1_, nr2_, nr3_, ngl_, ngg_, nkl_, nkg_, nspin_
      INTEGER :: nr1s_, nr2s_, nr3s_
      INTEGER :: nbnd_, nelu_, neld_, nat_, ntyp_, nacc_
      INTEGER :: file_version_
      INTEGER :: na_(nsx)
      INTEGER :: ngwkl_(npk)
      INTEGER :: ngwkg_(npk)
      REAL(dbl) :: acc_(nacx)
      REAL(dbl) :: trutim_, ecutwfc_, ecutrho_, alat_, ekinc_
      INTEGER ::  kunit_, k1_, k2_, k3_, nk1_, nk2_, nk3_
      REAL(dbl) :: dgauss_, nel_
      INTEGER :: ngauss_
      LOGICAL :: lgauss_
      INTEGER :: ntetra_
      LOGICAL :: ltetra_
      INTEGER :: natomwfc_ 
      REAL(dbl) :: gcutm_, gcuts_, dual_
      LOGICAL :: doublegrid_
      INTEGER :: modenum_ 
      LOGICAL :: lstres_, lforce_
      CHARACTER(LEN=80) :: t_, c_, tmp_dir_
      LOGICAL :: tupf_, gamma_only_
!
      INTEGER :: i
      INTEGER :: idum
      LOGICAL :: twrite_
      CHARACTER(LEN=30) :: sub_name = ' read_restart_header '
!
! ... Subroutine Body
!
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version_, ionode_id )
!
      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(sub_name,' Data Section not present in restart file ', 1)
!
      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) nfi_, nbeg_, nr1_, nr2_, nr3_, nr1s_, nr2s_, nr3s_, ngl_, ngg_, nkl_, nkg_, &
            nspin_, nbnd_, nel_, nelu_, neld_, nat_, ntyp_, nacc_, trutim_, ecutwfc_, ecutrho_, &
            alat_, ekinc_, kunit_, k1_, k2_, k3_, nk1_, nk2_, nk3_, dgauss_, ngauss_, lgauss_, &
            ntetra_, ltetra_, natomwfc_, gcutm_, gcuts_, dual_, doublegrid_, modenum_, lstres_, &
            lforce_, tupf_, gamma_only_
        END IF
!
        CALL mp_bcast( nfi_, ionode_id )
        CALL mp_bcast( nbeg_, ionode_id )
        CALL mp_bcast( nr1_, ionode_id )
        CALL mp_bcast( nr2_, ionode_id )
        CALL mp_bcast( nr3_, ionode_id )
        CALL mp_bcast( nr1s_, ionode_id )
        CALL mp_bcast( nr2s_, ionode_id )
        CALL mp_bcast( nr3s_, ionode_id )
        CALL mp_bcast( ngl_, ionode_id )
        CALL mp_bcast( ngg_, ionode_id )
        CALL mp_bcast( nkl_, ionode_id )
        CALL mp_bcast( nkg_, ionode_id )
        CALL mp_bcast( nspin_, ionode_id )
        CALL mp_bcast( nbnd_, ionode_id )
        CALL mp_bcast( nel_, ionode_id )
        CALL mp_bcast( nelu_, ionode_id )
        CALL mp_bcast( neld_, ionode_id )
        CALL mp_bcast( nat_, ionode_id )
        CALL mp_bcast( ntyp_, ionode_id )
        CALL mp_bcast( nacc_, ionode_id )
        CALL mp_bcast( trutim_, ionode_id )
        CALL mp_bcast( ecutwfc_, ionode_id )
        CALL mp_bcast( ecutrho_, ionode_id )
        CALL mp_bcast( alat_, ionode_id )
        CALL mp_bcast( ekinc_, ionode_id )
        CALL mp_bcast( kunit_, ionode_id )
        CALL mp_bcast( k1_, ionode_id )
        CALL mp_bcast( k2_, ionode_id )
        CALL mp_bcast( k3_, ionode_id )
        CALL mp_bcast( nk1_, ionode_id )
        CALL mp_bcast( nk2_, ionode_id )
        CALL mp_bcast( nk3_, ionode_id )
        CALL mp_bcast( dgauss_, ionode_id )
        CALL mp_bcast( ngauss_, ionode_id )
        CALL mp_bcast( lgauss_, ionode_id )
        CALL mp_bcast( ntetra_, ionode_id )
        CALL mp_bcast( ltetra_, ionode_id )
        CALL mp_bcast( natomwfc_, ionode_id )
        CALL mp_bcast( gcutm_, ionode_id )
        CALL mp_bcast( gcuts_, ionode_id )
        CALL mp_bcast( dual_, ionode_id )
        CALL mp_bcast( doublegrid_, ionode_id )
        CALL mp_bcast( modenum_, ionode_id )
        CALL mp_bcast( lstres_, ionode_id )
        CALL mp_bcast( lforce_, ionode_id )
        CALL mp_bcast( tupf_, ionode_id )
        CALL mp_bcast( gamma_only_, ionode_id )
!
        IF( ntyp_ > SIZE( na_ ) ) &
          CALL error(sub_name,' too many types ', ntyp_ )
        IF( ( nkg_ > SIZE( ngwkl_ ) ) .OR. ( nkg_ > SIZE( ngwkg_ ) ) ) &
          CALL error(sub_name,' too many k points ', nkg_ )
        IF( nacc_ > SIZE( acc_ ) ) &
          CALL error(sub_name,' too many accumulators ', nacc_ )
!
        IF( ntyp_ > SIZE( na ) ) &
          CALL error(sub_name,' wrong size for na ', ntyp_ )
        IF( ( nkg_ > SIZE( ngwk_l ) ) .OR. ( nkg_ > SIZE( ngwk_g ) ) ) &
          CALL error(sub_name,' wrong size for ngwk_l or ngwk_g ', nkg_ )
        IF( nacc_ > SIZE( acc ) ) &
          CALL error(sub_name,' wrong size for acc ', nacc_ )
!
        IF( .NOT. tovrw ) THEN
          IF( nat_ /= nat ) &
            CALL error(sub_name,' atoms numbers differ ',nat_)
          IF( ntyp_ /= ntyp ) &
            CALL error(sub_name,' types numbers differ ',ntyp_)
          IF( nkg_ /= nk_g ) &
            CALL error(sub_name,' k points numbers differ ',nkg_)
          IF( nspin_ /= nspin ) &
            CALL error(sub_name,' nspin differs ',nspin_)
          IF( nacc_ > nacc ) &
            CALL error(sub_name,' too many accumulators ',nacc_)
        END IF

        IF( ionode ) THEN
          READ(iuni) (na_(i),i=1,ntyp_), (ngwkl_(i),i=1,nkg_), (ngwkg_(i),i=1,nkg_), (acc_(i),i=1,nacc_)
        END IF

        CALL mp_bcast( na_, ionode_id )
        CALL mp_bcast( ngwkl_, ionode_id )
        CALL mp_bcast( ngwkg_, ionode_id )
        CALL mp_bcast( acc_, ionode_id )

        IF( .NOT. tovrw ) THEN
          DO i = 1, MIN( ntyp, ntyp_ )
            IF( na_(i) /= na(i) ) &
              CALL error(' read_restart_header ',' inconsistent number of atoms ',na_(i))
          END DO
        END IF

        IF( ionode ) THEN
          READ(iuni) t_, c_, tmp_dir_
        END IF

        CALL mp_bcast( t_, ionode_id )
        CALL mp_bcast( c_, ionode_id )
        CALL mp_bcast( tmp_dir_, ionode_id )

        !  Variables acc, trutim, nfi, ekinc, tupf are ALWAYS overwritten on output 
        acc(1:nacc_) = acc_(1:nacc_)
        trutim       = trutim_
        nfi          = nfi_
        ekinc        = ekinc_
        tupf         = tupf_
        gamma_only   = gamma_only_
 
        IF( tovrw ) THEN
          nfi = nfi_
          nbeg = nbeg_
          nr1 =  nr1_
          nr2 =  nr2_
          nr3 =  nr3_
          nr1s =  nr1s_
          nr2s =  nr2s_
          nr3s =  nr3s_
          ng_l = ngl_
          ng_g = ngg_
          nk_l = nkl_
          nk_g = nkg_
          nspin = nspin_
          nbnd = nbnd_
          nel = nel_
          nelu = nelu_
          neld = neld_
          nat = nat_
          ntyp = ntyp_
          nacc = nacc_
          trutim = trutim_
          ecutwfc = ecutwfc_
          ecutrho = ecutrho_
          alat = alat_
          ekinc = ekinc_
          kunit =  kunit_
          k1 =  k1_
          k2 = k2_
          k3 = k3_
          nk1 = nk1_
          nk2 = nk2_
          nk3 = nk3_
          dgauss = dgauss_
          ngauss = ngauss_
          lgauss = lgauss_
          ntetra = ntetra_
          ltetra =ltetra_
          natomwfc =  natomwfc_
          gcutm = gcutm_
          gcuts = gcuts_
          dual = dual_
          doublegrid = doublegrid_
          modenum = modenum_
          lstres = lstres_
          lforce = lforce_
          na(1:ntyp_) = na_(1:ntyp_)
          ngwk_l(1:nkg_) = ngwkl_(1:nkg_)
          ngwk_g(1:nkg_) = ngwkg_(1:nkg_)
          acc(1:nacc_) = acc_(1:nacc_)
          title = t_
          crystal = c_
          tmp_dir = tmp_dir_
        END IF

      ELSE

        IF( ionode ) THEN
          WRITE(iuni) idum
          WRITE(iuni) idum
          WRITE(iuni) idum
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_header, data not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

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
      INTEGER :: file_version_
      INTEGER :: idum
!
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
!
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_header, xdim not read from restart ' )")
!
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_xdim1(iuni, twrite, &
      npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER :: idum
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      INTEGER, INTENT(IN) :: npwx, nbndx
      INTEGER, INTENT(IN) :: nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
!
      IF( twrite ) THEN
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
        END IF
      ELSE
        CALL write_restart_xdim2(iuni)
      END IF
!
      RETURN
    END SUBROUTINE

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
      INTEGER :: idum
!
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
        WRITE(iuni) idum
      END IF 
!
      RETURN 
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_xdim1(iuni, tovrw, tread, &
      npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs )

! .. Subroutine output:
!    if tread is true then variables are read from file but the values are
!      not copied on output variables
!    if tovrw is true all variables are overvritten with values read from file
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      INTEGER, INTENT(INOUT) :: npwx, nbndx
      INTEGER, INTENT(INOUT) :: nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs
      INTEGER :: npwx_, nbndx_
      INTEGER :: nrx1_, nrx2_, nrx3_, nrxx_, nrx1s_, nrx2s_, nrx3s_, nrxxs_
      LOGICAL :: twrite_ 
      INTEGER :: file_version_
      INTEGER :: idum
!
! ... Subroutine Body
!
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version_, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_xdim ',' Data Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) npwx_, nbndx_, nrx1_, nrx2_, nrx3_, nrxx_, nrx1s_, nrx2s_, nrx3s_, nrxxs_
        END IF
        CALL mp_bcast( npwx_, ionode_id )
        CALL mp_bcast( nbndx_, ionode_id )
        CALL mp_bcast( nrx1_, ionode_id )
        CALL mp_bcast( nrx2_, ionode_id )
        CALL mp_bcast( nrx3_, ionode_id )
        CALL mp_bcast( nrxx_, ionode_id )
        CALL mp_bcast( nrx1s_, ionode_id )
        CALL mp_bcast( nrx2s_, ionode_id )
        CALL mp_bcast( nrx3s_, ionode_id )
        CALL mp_bcast( nrxxs_, ionode_id )

        IF( tovrw ) THEN
          npwx = npwx_
          nbndx = nbndx_
          nrx1 = nrx1_
          nrx2 =  nrx2_
          nrx3 =  nrx3_
          nrxx =  nrxx_
          nrx1s =  nrx1s_
          nrx2s =  nrx2s_
          nrx3s =  nrx3s_
          nrxxs =  nrxxs_
        END IF

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_xdim, xdim not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_xdim2(iuni)

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
      END IF

      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_xdim, xdim not read from restart ' )")

      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_symmetry1(iuni, twrite, &
      symm_type, sname, s, irt, nat, ftau, nsym, invsym, noinv )

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      CHARACTER(LEN=9), INTENT(IN) :: symm_type
      INTEGER, INTENT(IN) :: s(:,:,:)
      CHARACTER(LEN=45), INTENT(IN) :: sname(:)
      INTEGER, INTENT(IN) :: ftau(:,:)
      INTEGER, INTENT(IN) :: irt(:,:)
      INTEGER, INTENT(IN) :: nsym
      INTEGER, INTENT(IN) :: nat
      LOGICAL, INTENT(IN) :: invsym
      LOGICAL, INTENT(IN) :: noinv
      INTEGER :: idum
      INTEGER :: i,j
      CHARACTER(LEN=30) :: sub_name = ' write_restart_symmetry '
!
! ... Subroutine Body
!

      IF( twrite ) THEN

        IF( SIZE(s,3) < nsym ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( ( SIZE(irt,1) < nsym ) .OR. ( SIZE(irt,2) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( SIZE(ftau,2) < nsym ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( SIZE(sname) < nsym ) &
          CALL error( sub_name, ' wrong size ', 4 )

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) symm_type, nsym, invsym, noinv, nat
          WRITE(iuni) (s(:,:,i),i=1,nsym), ((irt(i,j),i=1,nsym),j=1,nat),  &
            (ftau(:,i),i=1,nsym), (sname(i),i=1,nsym)
        END IF

      ELSE

        CALL write_restart_symmetry2(iuni)

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_symmetry2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_symmetry1(iuni, tovrw, tread, &
      symm_type, sname, s, irt, nat, ftau, nsym, invsym, noinv )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast, mp_end

! .. Subroutine output:
!    if tread is true then variables are read from file but the values are
!      not copied on output variables
!    if tovrw is true all variables are overvritten with values read from file
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      CHARACTER(LEN=9), INTENT(INOUT) :: symm_type
      CHARACTER(LEN=45), INTENT(INOUT) :: sname(:)
      INTEGER, INTENT(INOUT) :: s(:,:,:)
      INTEGER, INTENT(INOUT) :: irt(:,:)
      INTEGER, INTENT(INOUT) :: ftau(:,:)
      INTEGER, INTENT(INOUT) :: nsym
      INTEGER, INTENT(INOUT) :: nat
      LOGICAL, INTENT(INOUT) :: invsym
      LOGICAL, INTENT(INOUT) :: noinv

      LOGICAL :: noinv_, invsym_, twrite_
      INTEGER :: nsym_, nat_, s_(3,3,48), ftau_(3,48)
      CHARACTER(LEN=9) :: symm_type_
      CHARACTER(LEN=45) :: sname_(48)
      INTEGER :: i, j, file_version_
      INTEGER :: idum
      INTEGER, ALLOCATABLE :: irt_(:,:)
      CHARACTER(LEN=30) :: sub_name = ' read_restart_symmetry '
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_symmetry ',' symmetries not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) symm_type_, nsym_, invsym_, noinv_, nat_
        END IF
        CALL mp_bcast( symm_type_, ionode_id )
        CALL mp_bcast( nsym_, ionode_id )
        CALL mp_bcast( invsym_, ionode_id )
        CALL mp_bcast( noinv_, ionode_id )
        CALL mp_bcast( nat_, ionode_id )
!
! ...   Check temporary local variables
        IF( SIZE(s_,3) < nsym_ ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE(ftau_,2) < nsym_ ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( SIZE(sname_) < nsym_ ) &
          CALL error( sub_name, ' wrong size ', 3 )
!
! ...   Check dummy variables
        IF( SIZE(s,3) < nsym_ ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( ( SIZE(irt,1) < nsym_ ) .OR. ( SIZE(irt,2) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( SIZE(ftau,2) < nsym_ ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF( SIZE(sname) < nsym_ ) &
          CALL error( sub_name, ' wrong size ', 7 )
!
        ALLOCATE( irt_(nsym_, nat_) )
!
        IF( ionode ) THEN
          READ(iuni) (s_(:,:,i),i=1,nsym_), ((irt_(i,j),i=1,nsym_),j=1,nat_),  &
            (ftau_(:,i),i=1,nsym_), (sname_(i),i=1,nsym_)
        END IF
        CALL mp_bcast( s_, ionode_id )
        CALL mp_bcast( sname_(:), ionode_id )
        CALL mp_bcast( ftau_, ionode_id )
        CALL mp_bcast( irt_, ionode_id )
!
        IF( tovrw ) THEN
          symm_type = symm_type_
          nsym      = nsym_
          invsym    = invsym_
          noinv     = noinv_
          nat       = nat_
          s(:,:,1:nsym_) = s_(:,:,1:nsym_)
          ftau(:,1:nsym_) = ftau_(:,1:nsym_)
          irt(1:nsym_, 1:nat_) = irt_(1:nsym_, 1:nat_)
          sname(1:nsym_) = sname_(1:nsym_)
        END IF
!
        DEALLOCATE( irt_ )

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum
          READ(iuni) idum
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_symmetry, symmetries not read from restart ' )")

      END IF
        
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_symmetry2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_symmetry, symmetries not read from restart ' )")
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_pseudo1(iuni, twrite, &
      zmesh, xmin, dx, r, rab, vnl, chi, oc, rho_at, &
      rho_atc, mesh, msh, nchi, lchi, numeric, cc, alpc, zp, aps, alps, zv, nlc, &
      nnl, lmax, lloc, bhstype, dion, betar, qqq, qfunc, qfcoef, rinner, nh, nbeta, &
      kkbeta, nqf, nqlc, ifqopt, lll, iver, tvanp, okvan, newpseudo, iexch, icorr, &
      igcx, igcc, lsda, a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd)
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      REAL(dbl), INTENT(IN) :: zmesh, xmin, dx
      REAL(dbl), INTENT(IN) :: r(0:), rab(0:), vnl(0:,0:), chi(0:,:)
      REAL(dbl), INTENT(IN) :: oc(:), rho_at(0:), rho_atc(0:)
      INTEGER, INTENT(IN) :: mesh, msh, nchi, lchi(:)
      LOGICAL, INTENT(IN) :: numeric
      REAL(dbl), INTENT(IN) :: cc(2), alpc(2), zp, aps(6,0:3), alps(3,0:3), zv
      INTEGER, INTENT(IN) :: nlc, nnl, lmax, lloc
      LOGICAL, INTENT(IN) :: bhstype
      REAL(dbl), INTENT(IN) :: dion(:,:), betar(0:,:), qqq(:,:), qfunc(0:,:,:)
      REAL(dbl), INTENT(IN) :: qfcoef(:,:,:,:), rinner(:)
      INTEGER, INTENT(IN) :: nh, nbeta, kkbeta, nqf, nqlc, ifqopt, lll(:), iver(3)
      LOGICAL, INTENT(IN) :: tvanp, okvan, newpseudo
      INTEGER, INTENT(IN) :: iexch, icorr, igcx, igcc
      LOGICAL, INTENT(IN) :: lsda
      REAL(dbl), INTENT(IN) :: a_nlcc, b_nlcc, alpha_nlcc
      LOGICAL, INTENT(IN) :: nlcc
      CHARACTER(LEN=2), INTENT(IN) :: psd
!
      INTEGER :: idum
      INTEGER :: mesh_, lloc_, nchi_, nbeta_, nqf_, nqlc_
      CHARACTER(LEN=30) :: sub_name = ' write_restart_pseudo '
!
! ... Subroutine Body
!
      IF( twrite ) THEN

! ...   Check dummy variables
        mesh_  = MAX( mesh, 1 )
        lloc_  = MAX( lloc, 1 )
        nchi_  = MAX( nchi, 1 )
        nbeta_ = MAX( nbeta, 1 )
        nqf_   = MAX( nqf, 1 )
        nqlc_  = MAX( nqlc, 1 )

        IF( SIZE(r) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE(rab) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( ( SIZE(vnl,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(vnl,2) < ( lloc_ + 1 ) ) ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( ( SIZE(chi,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(chi,2) < nchi_ ) ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( SIZE(oc) < nchi_ ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_at) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF( SIZE(rho_atc) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 7 )
        IF( SIZE(lchi) < nchi_ ) &
          CALL error( sub_name, ' wrong size ', 8 )
        IF( ( SIZE(dion,1) < nbeta_ ) .OR. ( SIZE(dion,2) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(betar,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(betar,2) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(qqq,1) < nbeta_  ) .OR. ( SIZE(qqq,2) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( ( SIZE(qfunc,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(qfunc,2) < nbeta_ ) .OR. &
            ( SIZE(qfunc,3) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 11 )
        IF( ( SIZE(qfcoef,1) < nqf_ ) .OR. ( SIZE(qfcoef,2) < nqlc_ ) .OR. &
            ( SIZE(qfcoef,3) < nbeta_ ) .OR. ( SIZE(qfcoef,4) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 12 )
        IF( SIZE(rinner) < nqlc_ ) &
          CALL error( sub_name, ' wrong size ', 13 )
        IF( SIZE(lll) < nbeta_ ) &
          CALL error( sub_name, ' wrong size ', 14 )

        IF( ionode ) THEN

          WRITE(iuni) twrite, file_version
          WRITE(iuni) zmesh, xmin, dx, mesh, msh, nchi, numeric, zp, zv, nlc, nnl, lmax, lloc, &
            bhstype, nh, nbeta, kkbeta, nqf, nqlc, ifqopt, tvanp, okvan, newpseudo, &
            iexch, icorr, igcx, igcc, lsda, a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd
          WRITE(iuni) r( 0:mesh_ ), rab( 0:mesh_ ), &
            vnl( 0:mesh_, 0:lloc_ ), chi( 0:mesh_, 1:nchi_ ), &
            oc( 1:nchi_ ), rho_at( 0:mesh_ ), rho_atc( 0:mesh_ ), lchi( 1:nchi_ )
          WRITE(iuni) cc(1:2), alpc(1:2), aps(1:6,0:3), alps(1:3,0:3)
          WRITE(iuni) dion( 1:nbeta_, 1:nbeta_ ), betar( 0:mesh_, 1:nbeta_ ), &
            qqq( 1:nbeta_, 1:nbeta_ ), qfunc( 0:mesh_, 1:nbeta_, 1:nbeta_ ), &
            qfcoef( 1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_ ), &
            rinner( 1:nqlc_ ), lll( 1:nbeta_ ), iver(1:3)

        END IF

      ELSE

        CALL write_restart_pseudo2(iuni)

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_pseudo3(iuni, twrite, &
      generated, date_author, comment, psd, typ, tvanp, nlcc, dft, zp, etotps, &
      ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, els, lchi, oc, r, rab, &
      rho_atc, vloc, lll, kkbeta, beta, nd, dion, nqf, nqlc, rinner, qqq, &
      qfunc, qfcoef, chi, rho_at )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
!
      CHARACTER(LEN=80):: generated   ! 
      CHARACTER(LEN=80):: date_author ! Misc info
      CHARACTER(LEN=80):: comment     !
      CHARACTER(LEN=2) :: psd       ! Element label
      CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
      LOGICAL :: tvanp             ! .true. if Ultrasoft
      LOGICAL :: nlcc               ! Non linear core corrections
      CHARACTER(LEN=20) :: dft      ! Exch-Corr type
      REAL(dbl) :: zp               ! z valence
      REAL(dbl) :: etotps           ! total energy
      REAL(dbl) :: ecutwfc          ! suggested cut-off for wfc
      REAL(dbl) :: ecutrho          ! suggested cut-off for rho
      INTEGER :: nv                 ! UPF file version number
      INTEGER :: lmax               ! maximum angular momentum component
      INTEGER :: mesh               ! number of point in the radial mesh
      INTEGER :: nwfc               ! number of wavefunctions
      INTEGER :: nbeta              ! number of projectors
      CHARACTER(LEN=2) :: els(:)  ! els(nwfc)
      INTEGER :: lchi(:)   ! lchi(nwfc)
      REAL(dbl) :: oc(:)   ! oc(nwfc)
      REAL(dbl) :: r(:)    ! r(mesh)
      REAL(dbl) :: rab(:)  ! rab(mesh)
      REAL(dbl) :: rho_atc(:) ! rho_atc(mesh)
      REAL(dbl) :: vloc(:)    ! vloc(mesh)
      INTEGER :: lll(:)       ! lll(nbeta)
      INTEGER :: kkbeta(:)    ! kkbeta(nbeta)
      REAL(dbl) :: beta(:,:)  ! beta(mesh,nbeta)
      INTEGER :: nd
      REAL(dbl) :: dion(:,:)  ! dion(nbeta,nbeta)
      INTEGER :: nqf
      INTEGER :: nqlc
      REAL(dbl) :: rinner(:)  ! rinner(0:2*lmax)
      REAL(dbl) :: qqq(:,:)   ! qqq(nbeta,nbeta)
      REAL(dbl) :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
      REAL(dbl) :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
      REAL(dbl) :: chi(:,:) !  chi(mesh,nwfc)
      REAL(dbl) :: rho_at(:) !  rho_at(mesh)
!
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' write_restart_pseudo '
!
! ... Subroutine Body
!
      IF( twrite ) THEN

! ...   Check dummy variables
        IF( SIZE(els) < nwfc ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE(lchi) < nwfc ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( SIZE(oc) < nwfc ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( SIZE(r) < ( mesh + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( SIZE(rab) < ( mesh + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_atc) < ( mesh + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF(  SIZE(vloc) < ( mesh + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 7 )
        IF( SIZE(lll) < nbeta ) &
          CALL error( sub_name, ' wrong size ', 8 )
        IF( SIZE(kkbeta) < nbeta ) &
          CALL error( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(beta,1) < ( mesh + 1 ) ) .OR. ( SIZE(beta,2) < nbeta ) ) &
          CALL error( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(dion,1) < nbeta ) .OR. ( SIZE(dion,2) < nbeta ) ) &
          CALL error( sub_name, ' wrong size ', 11 )
        IF( SIZE(rinner) < nqlc ) &
          CALL error( sub_name, ' wrong size ', 12 )
        IF( ( SIZE(qqq,1) < nbeta ) .OR. ( SIZE(qqq,2) < nbeta ) ) &
          CALL error( sub_name, ' wrong size ', 13 )
        IF( ( SIZE(qfunc,1) < ( mesh + 1 ) ) .OR. ( SIZE(qfunc,2) < nbeta ) .OR. &
            ( SIZE(qfunc,3) < nbeta ) ) &
          CALL error( sub_name, ' wrong size ', 14 )
        IF( ( SIZE(qfcoef,1) < nqf ) .OR. ( SIZE(qfcoef,2) < nqlc ) .OR. &
            ( SIZE(qfcoef,3) < nbeta ) .OR. ( SIZE(qfcoef,4) < nbeta ) ) &
          CALL error( sub_name, ' wrong size ', 15 )
        IF( ( SIZE(chi,1) < ( mesh + 1 ) ) .OR. ( SIZE(chi,2) < nwfc ) ) &
          CALL error( sub_name, ' wrong size ', 16 )
        IF( SIZE(rho_at) < ( mesh + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 17 )

        IF( ionode ) THEN

          WRITE(iuni) twrite, file_version
!
          WRITE(iuni) generated, date_author, comment, psd, typ, tvanp, nlcc, dft, &
           zp, etotps, ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, nd, nqf, nqlc
!           
          WRITE(iuni) els(1:nwfc), lchi(nwfc), oc(nwfc), r(0:mesh), rab(0:mesh), &
            rho_atc(0:mesh), vloc(0:mesh), lll(1:nbeta), kkbeta(1:nbeta), beta(0:mesh,1:nbeta), &
            dion(1:nbeta,1:nbeta), rinner(1:nqlc), qqq(1:nbeta,1:nbeta), &
            qfunc(0:mesh, 1:nbeta, 1:nbeta), qfcoef(1:nqf, 1:nqlc, 1:nbeta, 1:nbeta), &
            chi(0:mesh, nwfc), rho_at(0:mesh) 

          WRITE(iuni) idum
          WRITE(iuni) idum

        END IF

      ELSE

        CALL write_restart_pseudo2(iuni)

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_pseudo2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_pseudo1(iuni, tovrw, tread, &
      zmesh, xmin, dx, r, rab, vnl, chi, oc, rho_at, &
      rho_atc, mesh, msh, nchi, lchi, numeric, cc, alpc, zp, aps, alps, zv, nlc, &
      nnl, lmax, lloc, bhstype, dion, betar, qqq, qfunc, qfcoef, rinner, nh, nbeta, &
      kkbeta, nqf, nqlc, ifqopt, lll, iver, tvanp, okvan, newpseudo, iexch, icorr, &
      igcx, igcc, lsda, a_nlcc, b_nlcc, alpha_nlcc, nlcc, psd )

! .. Subroutine output:
!    if tread is true then variables are read from file but the values are
!      not copied on output variables
!    if tovrw is true all variables are overvritten with values read from file
!

      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      REAL(dbl), INTENT(INOUT) :: zmesh, xmin, dx
      REAL(dbl), INTENT(INOUT) :: r(0:), rab(0:), vnl(0:,0:), chi(0:,:)
      REAL(dbl), INTENT(INOUT) :: oc(:), rho_at(0:), rho_atc(0:)
      INTEGER, INTENT(INOUT) :: mesh, msh, nchi, lchi(:)
      LOGICAL, INTENT(INOUT) :: numeric
      REAL(dbl), INTENT(INOUT) :: cc(2), alpc(2), zp, aps(6,0:3), alps(3,0:3), zv
      INTEGER, INTENT(INOUT) :: nlc, nnl, lmax, lloc
      LOGICAL, INTENT(INOUT) :: bhstype
      REAL(dbl), INTENT(INOUT) :: dion(:,:), betar(0:,:), qqq(:,:), qfunc(0:,:,:)
      REAL(dbl), INTENT(INOUT) :: qfcoef(:,:,:,:), rinner(:)
      INTEGER, INTENT(INOUT) :: nh, nbeta, kkbeta, nqf, nqlc, ifqopt, lll(:), iver(:)
      LOGICAL, INTENT(INOUT) :: tvanp, okvan, newpseudo
      INTEGER, INTENT(INOUT) :: iexch, icorr, igcx, igcc
      LOGICAL, INTENT(INOUT) :: lsda
      REAL(dbl), INTENT(INOUT) :: a_nlcc, b_nlcc, alpha_nlcc
      LOGICAL, INTENT(INOUT) :: nlcc
      CHARACTER(LEN=2), INTENT(INOUT) :: psd
     
!
      REAL(dbl) :: zmesh_, xmin_, dx_
      REAL(dbl), ALLOCATABLE :: r_(:), rab_(:), vnl_(:,:), chi_(:,:)
      REAL(dbl), ALLOCATABLE :: oc_(:), rho_at_(:), rho_atc_(:)
      INTEGER, ALLOCATABLE :: lchi_(:)
      INTEGER :: mesh_, msh_, nchi_
      LOGICAL :: numeric_
      REAL(dbl) :: cc_(2), alpc_(2), zp_, aps_(6,0:3), alps_(3,0:3), zv_
      INTEGER :: nlc_, nnl_, lmax_, lloc_
      LOGICAL :: bhstype_
      REAL(dbl), ALLOCATABLE :: dion_(:,:), betar_(:,:), qqq_(:,:)
      REAL(dbl), ALLOCATABLE :: qfunc_(:,:,:)
      REAL(dbl), ALLOCATABLE :: qfcoef_(:,:,:,:), rinner_(:)
      INTEGER, ALLOCATABLE :: lll_(:)
      INTEGER :: nh_, nbeta_, kkbeta_, nqf_, nqlc_, ifqopt_, iver_(3)
      LOGICAL :: tvanp_, okvan_, newpseudo_
      INTEGER :: iexch_, icorr_, igcx_, igcc_
      LOGICAL :: lsda_
      REAL(dbl) :: a_nlcc_, b_nlcc_, alpha_nlcc_
      LOGICAL :: nlcc_
      CHARACTER(LEN=2) :: psd_
      LOGICAL :: twrite_
!
      LOGICAL :: file_version_
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_pseudo '
      INTEGER :: mesh__, lloc__, nchi__, nbeta__, nqf__, nqlc__
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_pseudo ',' pseudo not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) zmesh_, xmin_, dx_, mesh_, msh_, nchi_, numeric_, zp_, zv_, nlc_, nnl_, lmax_, &
            lloc_, bhstype_, nh_, nbeta_, kkbeta_, nqf_, nqlc_, ifqopt_, tvanp_, okvan_, newpseudo_, &
            iexch_, icorr_, igcx_, igcc_, lsda_, a_nlcc_, b_nlcc_, alpha_nlcc_, nlcc_, psd_
        END IF

        CALL mp_bcast( zmesh_, ionode_id )
        CALL mp_bcast( xmin_, ionode_id )
        CALL mp_bcast( dx_, ionode_id )
        CALL mp_bcast( mesh_, ionode_id )
        CALL mp_bcast( msh_, ionode_id )
        CALL mp_bcast( nchi_, ionode_id )
        CALL mp_bcast( numeric_, ionode_id )
        CALL mp_bcast( zp_, ionode_id )
        CALL mp_bcast( zv_, ionode_id )
        CALL mp_bcast( nlc_, ionode_id )
        CALL mp_bcast( nnl_, ionode_id )
        CALL mp_bcast( lmax_, ionode_id )
        CALL mp_bcast( lloc_, ionode_id )
        CALL mp_bcast( bhstype_, ionode_id )
        CALL mp_bcast( nh_, ionode_id )
        CALL mp_bcast( nbeta_, ionode_id )
        CALL mp_bcast( kkbeta_, ionode_id )
        CALL mp_bcast( nqf_, ionode_id )
        CALL mp_bcast( nqlc_, ionode_id )
        CALL mp_bcast( ifqopt_, ionode_id )
        CALL mp_bcast( tvanp_, ionode_id )
        CALL mp_bcast( okvan_, ionode_id )
        CALL mp_bcast( newpseudo_, ionode_id )
        CALL mp_bcast( iexch_, ionode_id )
        CALL mp_bcast( icorr_, ionode_id )
        CALL mp_bcast( igcx_, ionode_id )
        CALL mp_bcast( igcc_, ionode_id )
        CALL mp_bcast( lsda_, ionode_id )
        CALL mp_bcast( a_nlcc_, ionode_id )
        CALL mp_bcast( b_nlcc_, ionode_id )
        CALL mp_bcast( alpha_nlcc_, ionode_id )
        CALL mp_bcast( nlcc_, ionode_id )
        CALL mp_bcast( psd_, ionode_id )

        mesh__  = MAX( mesh_, 1 )
        lloc__  = MAX( lloc_, 1 )
        nchi__  = MAX( nchi_, 1 )
        nbeta__ = MAX( nbeta_, 1 )
        nqf__   = MAX( nqf_, 1 )
        nqlc__  = MAX( nqlc_, 1 )

! ...   Check dummy variables
        IF( SIZE(r) < ( mesh__ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE(rab) < ( mesh__ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( ( SIZE(vnl,1) < ( mesh__ + 1 ) ) .OR. ( SIZE(vnl,2) < ( lloc__ + 1 ) ) ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( ( SIZE(chi,1) < ( mesh__ + 1 ) ) .OR. ( SIZE(chi,2) < nchi__ ) ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( SIZE(oc) < nchi__ ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_at) < ( mesh__ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF( SIZE(rho_atc) < ( mesh__ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 7 )
        IF( SIZE(lchi) < nchi__ ) &
          CALL error( sub_name, ' wrong size ', 8 )
        IF( ( SIZE(dion,1) < nbeta__ ) .OR. ( SIZE(dion,2) < nbeta__ ) ) &
          CALL error( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(betar,1) < ( mesh__ + 1 ) ) .OR. ( SIZE(betar,2) < nbeta__ ) ) &
          CALL error( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(qqq,1) < nbeta__ ) .OR. ( SIZE(qqq,2) < nbeta__ ) ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( ( SIZE(qfunc,1) < ( mesh__ + 1 ) ) .OR. ( SIZE(qfunc,2) < nbeta__ ) .OR. &
            ( SIZE(qfunc,3) < nbeta__ ) ) &
          CALL error( sub_name, ' wrong size ', 11 )
        IF( ( SIZE(qfcoef,1) < nqf__ ) .OR. ( SIZE(qfcoef,2) < nqlc__ ) .OR. &
            ( SIZE(qfcoef,3) < nbeta__ ) .OR. ( SIZE(qfcoef,4) < nbeta__ ) ) &
          CALL error( sub_name, ' wrong size ', 12 )
        IF( SIZE(rinner) < nqlc__ ) &
          CALL error( sub_name, ' wrong size ', 13 )
        IF( SIZE(lll) < nbeta__ ) &
          CALL error( sub_name, ' wrong size ', 14 )

        IF( mesh__ < 0 ) &
          CALL error( sub_name, ' wrong value ', 1 )
        IF( lloc__ < 0 ) &
          CALL error( sub_name, ' wrong value ', 2 )
        IF( nchi__ < 1 ) &
          CALL error( sub_name, ' wrong value ', 3 )
        IF( nbeta__ < 1 ) &
          CALL error( sub_name, ' wrong value ', 4 )
        IF( nqf__ < 1 ) &
          CALL error( sub_name, ' wrong value ', 5 )
        IF( nqlc__ < 1 ) &
          CALL error( sub_name, ' wrong value ', 6 )

        ALLOCATE ( r_(0:mesh__), rab_(0:mesh__), vnl_(0:mesh__,0:lloc__), chi_(0:mesh__,nchi__) )
        ALLOCATE ( oc_(nchi__), rho_at_(0:mesh__), rho_atc_(0:mesh__) )
        ALLOCATE ( lchi_(nchi__) )
        ALLOCATE ( dion_(nbeta__,nbeta__), betar_(0:mesh__,nbeta__), qqq_(nbeta__,nbeta__) )
        ALLOCATE ( qfunc_(0:mesh__,nbeta__,nbeta__) )
        ALLOCATE ( qfcoef_(nqf__,nqlc__,nbeta__,nbeta__), rinner_(nqlc__) )
        ALLOCATE ( lll_(nbeta__) )

        IF( ionode ) THEN
          READ(iuni) r_(0:mesh__), rab_(0:mesh__), vnl_(0:mesh__,0:lloc__), chi_(0:mesh__,1:nchi__), &
            oc_(1:nchi__), rho_at_(0:mesh__), rho_atc_(0:mesh__), lchi_(1:nchi__)
          READ(iuni) cc_(1:2), alpc_(1:2), aps_(1:6,0:3), alps_(1:3,0:3)
          READ(iuni) dion_(1:nbeta__,1:nbeta__), betar_(0:mesh__,1:nbeta__), qqq_(1:nbeta__,1:nbeta__), &
            qfunc_(0:mesh_, 1:nbeta__, 1:nbeta__), qfcoef_(1:nqf__, 1:nqlc__, 1:nbeta__, 1:nbeta__), &
            rinner_(1:nqlc__), lll_(1:nbeta__), iver_(1:3)
        END IF

        CALL mp_bcast( r_, ionode_id )
        CALL mp_bcast( rab_, ionode_id )
        CALL mp_bcast( vnl_, ionode_id )
        CALL mp_bcast( chi_, ionode_id )
        CALL mp_bcast( oc_, ionode_id )
        CALL mp_bcast( rho_at_, ionode_id )
        CALL mp_bcast( rho_atc_, ionode_id )
        CALL mp_bcast( lchi_, ionode_id )
        CALL mp_bcast( cc_, ionode_id )
        CALL mp_bcast( alpc_, ionode_id )
        CALL mp_bcast( aps_, ionode_id )
        CALL mp_bcast( alps_, ionode_id )
        CALL mp_bcast( dion_, ionode_id )
        CALL mp_bcast( betar_, ionode_id )
        CALL mp_bcast( qqq_, ionode_id )
        CALL mp_bcast( qfunc_, ionode_id )
        CALL mp_bcast( qfcoef_, ionode_id )
        CALL mp_bcast( rinner_, ionode_id )
        CALL mp_bcast( lll_, ionode_id )
        CALL mp_bcast( iver_, ionode_id )

        IF( tovrw ) THEN
          zmesh     = zmesh_
          xmin      = xmin_
          dx        = dx_
          mesh      = mesh_
          msh       = msh_
          nchi      = nchi_
          numeric   = numeric_
          zp        = zp_
          zv        = zv_
          nlc       = nlc_
          nnl       = nnl_
          psd       = psd_
          lmax      = lmax_
          lloc      = lloc_
          bhstype   = bhstype_
          nh        = nh_
          nbeta     = nbeta_
          kkbeta    = kkbeta_
          nqf       = nqf_
          nqlc      = nqlc_
          ifqopt    = ifqopt_
          tvanp     = tvanp_
          okvan     = okvan_
          newpseudo = newpseudo_
          iexch     = iexch_
          icorr     = icorr_
          igcx      = igcx_
          igcc      = igcc_
          lsda      = lsda_
          a_nlcc    = a_nlcc_
          b_nlcc    = b_nlcc_
          alpha_nlcc = alpha_nlcc_
          nlcc      = nlcc_
          r(0:mesh__) = r_(0:mesh__)
          rab(0:mesh__) = rab_(0:mesh__)
          vnl(0:mesh__,0:lloc__) = vnl_(0:mesh__,0:lloc__)
          chi(0:mesh__,1:nchi__) = chi_(0:mesh__,1:nchi__)
          oc(1:nchi__) = oc_(1:nchi__)
          rho_at(0:mesh__) = rho_at_(0:mesh__)
          rho_atc(0:mesh__) = rho_atc_(0:mesh__)
          lchi(1:nchi__) = lchi_(1:nchi__)
          cc(1:2) = cc_(1:2)
          alpc(1:2) = alpc_(1:2)
          aps(1:6,0:3) = aps_(1:6,0:3)
          alps(1:3,0:3) = alps_(1:3,0:3)
          dion(1:nbeta__,1:nbeta__) = dion_(1:nbeta__,1:nbeta__)
          betar(0:mesh__,1:nbeta__) = betar_(0:mesh__,1:nbeta__)
          qqq(1:nbeta__,1:nbeta__) = qqq_(1:nbeta__,1:nbeta__)
          qfunc(0:mesh__, 1:nbeta__, 1:nbeta__) = qfunc_(0:mesh__, 1:nbeta__, 1:nbeta__)
          qfcoef(1:nqf__, 1:nqlc__, 1:nbeta__, 1:nbeta__) = qfcoef_(1:nqf__, 1:nqlc__, 1:nbeta__, 1:nbeta__)
          rinner(1:nqlc__) = rinner_(1:nqlc__)
          lll(1:nbeta__) = lll_(1:nbeta__)
          iver(1:3) = iver_(1:3)
        END IF

        DEALLOCATE ( r_, rab_, vnl_, chi_ )
        DEALLOCATE ( oc_, rho_at_, rho_atc_ )
        DEALLOCATE ( lchi_ )
        DEALLOCATE ( dion_, betar_, qqq_ )
        DEALLOCATE ( qfunc_ )
        DEALLOCATE ( qfcoef_, rinner_ )
        DEALLOCATE ( lll_ )

      ELSE

        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum

        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_pseudo, pseudo not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_pseudo3(iuni, tovrw, tread, &
      generated, date_author, comment, psd, typ, tvanp, nlcc, dft, zp, etotps, &
      ecutwfc, ecutrho, nv, lmax, mesh, nwfc, nbeta, els, lchi, oc, r, rab, &
      rho_atc, vloc, lll, kkbeta, beta, nd, dion, nqf, nqlc, rinner, qqq, &
      qfunc, qfcoef, chi, rho_at )
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw, tread
!
      CHARACTER(LEN=80):: generated   ! 
      CHARACTER(LEN=80):: date_author ! Misc info
      CHARACTER(LEN=80):: comment     !
      CHARACTER(LEN=2) :: psd       ! Element label
      CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
      LOGICAL :: tvanp             ! .true. if Ultrasoft
      LOGICAL :: nlcc               ! Non linear core corrections
      CHARACTER(LEN=20) :: dft      ! Exch-Corr type
      REAL(dbl) :: zp               ! z valence
      REAL(dbl) :: etotps           ! total energy
      REAL(dbl) :: ecutwfc          ! suggested cut-off for wfc
      REAL(dbl) :: ecutrho          ! suggested cut-off for rho
      INTEGER :: nv                 ! UPF file version number
      INTEGER :: lmax               ! maximum angular momentum component
      INTEGER :: mesh               ! number of point in the radial mesh
      INTEGER :: nwfc               ! number of wavefunctions
      INTEGER :: nbeta              ! number of projectors
      CHARACTER(LEN=2) :: els(:)  ! els(nwfc)
      INTEGER :: lchi(:)   ! lchi(nwfc)
      REAL(dbl) :: oc(:)   ! oc(nwfc)
      REAL(dbl) :: r(0:)    ! r(mesh)
      REAL(dbl) :: rab(0:)  ! rab(mesh)
      REAL(dbl) :: rho_atc(0:) ! rho_atc(mesh)
      REAL(dbl) :: vloc(0:)    ! vloc(mesh)
      INTEGER :: lll(:)       ! lll(nbeta)
      INTEGER :: kkbeta(:)    ! kkbeta(nbeta)
      REAL(dbl) :: beta(0:,:)  ! beta(mesh,nbeta)
      INTEGER :: nd
      REAL(dbl) :: dion(:,:)  ! dion(nbeta,nbeta)
      INTEGER :: nqf
      INTEGER :: nqlc
      REAL(dbl) :: rinner(:)  ! rinner(0:2*lmax)
      REAL(dbl) :: qqq(:,:)   ! qqq(nbeta,nbeta)
      REAL(dbl) :: qfunc(0:,:,:) ! qfunc(mesh,nbeta,nbeta)
      REAL(dbl) :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
      REAL(dbl) :: chi(0:,:) !  chi(mesh,nwfc)
      REAL(dbl) :: rho_at(0:) !  rho_at(mesh)
!
!
      CHARACTER(LEN=80):: generated_   ! 
      CHARACTER(LEN=80):: date_author_ ! Misc info
      CHARACTER(LEN=80):: comment_     !
      CHARACTER(LEN=2) :: psd_       ! Element label
      CHARACTER(LEN=20) :: typ_      ! Pseudo type ( NC or US )
      LOGICAL  :: tvanp_             ! .true. if Ultrasoft
      LOGICAL :: nlcc_               ! Non linear core corrections
      CHARACTER(LEN=20) :: dft_      ! Exch-Corr type
      REAL(dbl) :: zp_               ! z valence
      REAL(dbl) :: etotps_           ! total energy
      REAL(dbl) :: ecutwfc_          ! suggested cut-off for wfc
      REAL(dbl) :: ecutrho_          ! suggested cut-off for rho
      INTEGER :: nv_                 ! UPF file version number
      INTEGER :: lmax_               ! maximum angular momentum component
      INTEGER :: mesh_               ! number of point in the radial mesh
      INTEGER :: nwfc_               ! number of wavefunctions
      INTEGER :: nbeta_              ! number of projectors
      CHARACTER(LEN=2), ALLOCATABLE :: els_(:)  ! els(nwfc)
      INTEGER, ALLOCATABLE :: lchi_(:)   ! lchi(nwfc)
      REAL(dbl), ALLOCATABLE :: oc_(:)   ! oc(nwfc)
      REAL(dbl), ALLOCATABLE :: r_(:)    ! r(mesh)
      REAL(dbl), ALLOCATABLE :: rab_(:)  ! rab(mesh)
      REAL(dbl), ALLOCATABLE :: rho_atc_(:) ! rho_atc(mesh)
      REAL(dbl), ALLOCATABLE :: vloc_(:)    ! vloc(mesh)
      INTEGER, ALLOCATABLE :: lll_(:)       ! lll(nbeta)
      INTEGER, ALLOCATABLE :: kkbeta_(:)    ! kkbeta(nbeta)
      REAL(dbl), ALLOCATABLE :: beta_(:,:)  ! beta(mesh,nbeta)
      INTEGER :: nd_
      REAL(dbl), ALLOCATABLE :: dion_(:,:)  ! dion(nbeta,nbeta)
      INTEGER :: nqf_
      INTEGER :: nqlc_
      REAL(dbl), ALLOCATABLE :: rinner_(:)  ! rinner(0:2*lmax)
      REAL(dbl), ALLOCATABLE :: qqq_(:,:)   ! qqq(nbeta,nbeta)
      REAL(dbl), ALLOCATABLE :: qfunc_(:,:,:) ! qfunc(mesh,nbeta,nbeta)
      REAL(dbl), ALLOCATABLE :: qfcoef_(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
      REAL(dbl), ALLOCATABLE :: chi_(:,:) !  chi(mesh,nwfc)
      REAL(dbl), ALLOCATABLE :: rho_at_(:) !  rho_at(mesh)
!
!
      LOGICAL :: file_version_, twrite_
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_pseudo '
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(sub_name, ' pseudo not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) generated_, date_author_, comment_, psd_, typ_, tvanp_, nlcc_, dft_, &
           zp_, etotps_, ecutwfc_, ecutrho_, nv_, lmax_, mesh_, nwfc_, nbeta_, nd_, nqf_, nqlc_
        END IF
        CALL mp_bcast( generated_, ionode_id )
        CALL mp_bcast( date_author_, ionode_id )
        CALL mp_bcast( comment_, ionode_id )
        CALL mp_bcast( psd_, ionode_id )
        CALL mp_bcast( typ_, ionode_id )
        CALL mp_bcast( tvanp_, ionode_id )
        CALL mp_bcast( nlcc_, ionode_id )
        CALL mp_bcast( dft_, ionode_id )
        CALL mp_bcast( zp_, ionode_id )
        CALL mp_bcast( etotps_, ionode_id )
        CALL mp_bcast( ecutwfc_, ionode_id )
        CALL mp_bcast( ecutrho_, ionode_id )
        CALL mp_bcast( nv_, ionode_id )
        CALL mp_bcast( lmax_, ionode_id )
        CALL mp_bcast( mesh_, ionode_id )
        CALL mp_bcast( nwfc_, ionode_id )
        CALL mp_bcast( nbeta_, ionode_id )
        CALL mp_bcast( nd_, ionode_id )
        CALL mp_bcast( nqf_, ionode_id )
        CALL mp_bcast( nqlc_, ionode_id )

        IF( mesh_ < 0 ) &
          CALL error( sub_name, ' wrong value ', 1 )
        IF( nwfc_ < 1 ) &
          CALL error( sub_name, ' wrong value ', 2 )
        IF( nbeta_ < 1 ) &
          CALL error( sub_name, ' wrong value ', 3 )
        IF( nqf_ < 1 ) &
          CALL error( sub_name, ' wrong value ', 4 )
        IF( nqlc_ < 1 ) &
          CALL error( sub_name, ' wrong value ', 5 )


! ...   Check dummy variables
        IF( SIZE(els) < nwfc_ ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE(lchi) < nwfc_ ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( SIZE(oc) < nwfc_ ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( SIZE(r) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( SIZE(rab) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( SIZE(rho_atc) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF(  SIZE(vloc) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 7 )
        IF( SIZE(lll) < nbeta_ ) &
          CALL error( sub_name, ' wrong size ', 8 )
        IF( SIZE(kkbeta) < nbeta_ ) &
          CALL error( sub_name, ' wrong size ', 9 )
        IF( ( SIZE(beta,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(beta,2) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 10 )
        IF( ( SIZE(dion,1) < nbeta_ ) .OR. ( SIZE(dion,2) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 11 )
        IF( SIZE(rinner) < nqlc_ ) &
          CALL error( sub_name, ' wrong size ', 12 )
        IF( ( SIZE(qqq,1) < nbeta_ ) .OR. ( SIZE(qqq,2) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 13 )
        IF( ( SIZE(qfunc,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(qfunc,2) < nbeta_ ) .OR. &
            ( SIZE(qfunc,3) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 14 )
        IF( ( SIZE(qfcoef,1) < nqf_ ) .OR. ( SIZE(qfcoef,2) < nqlc_ ) .OR. &
            ( SIZE(qfcoef,3) < nbeta_ ) .OR. ( SIZE(qfcoef,4) < nbeta_ ) ) &
          CALL error( sub_name, ' wrong size ', 15 )
        IF( ( SIZE(chi,1) < ( mesh_ + 1 ) ) .OR. ( SIZE(chi,2) < nwfc_ ) ) &
          CALL error( sub_name, ' wrong size ', 16 )
        IF( SIZE(rho_at) < ( mesh_ + 1 ) ) &
          CALL error( sub_name, ' wrong size ', 17 )

        ALLOCATE( els_( nwfc_ ),  lchi_(nwfc_), oc_(nwfc_), r_(0:mesh_), rab_(0:mesh_), &
            rho_atc_(0:mesh_), vloc_(0:mesh_), lll_(1:nbeta_), kkbeta_(1:nbeta_), &
            beta_(0:mesh_,1:nbeta_), &
            dion_(1:nbeta_,1:nbeta_), rinner_(1:nqlc_), qqq_(1:nbeta_,1:nbeta_), &
            qfunc_(0:mesh_, 1:nbeta_, 1:nbeta_), qfcoef_(1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_), &
            chi_(0:mesh_, nwfc_), rho_at_(0:mesh_) &
        )
 

        IF( ionode ) THEN
!           
          READ(iuni) els_(1:nwfc_), lchi_(nwfc_), oc_(nwfc_), r_(0:mesh_), rab_(0:mesh_), &
            rho_atc_(0:mesh_), vloc_(0:mesh_), lll_(1:nbeta_), kkbeta_(1:nbeta_), &
            beta_(0:mesh_,1:nbeta_), &
            dion_(1:nbeta_,1:nbeta_), rinner_(1:nqlc_), qqq_(1:nbeta_,1:nbeta_), &
            qfunc_(0:mesh_, 1:nbeta_, 1:nbeta_), qfcoef_(1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_), &
            chi_(0:mesh_, nwfc_), rho_at_(0:mesh_) 

          READ(iuni) idum
          READ(iuni) idum

        END IF

        CALL mp_bcast( els_(1:nwfc_), ionode_id ) 
        CALL mp_bcast( lchi_(nwfc_), ionode_id ) 
        CALL mp_bcast( oc_(nwfc_), ionode_id ) 
        CALL mp_bcast( r_(0:mesh_), ionode_id ) 
        CALL mp_bcast( rab_(0:mesh_), ionode_id ) 
        CALL mp_bcast( rho_atc_(0:mesh_), ionode_id ) 
        CALL mp_bcast( vloc_(0:mesh_), ionode_id ) 
        CALL mp_bcast( lll_(1:nbeta_), ionode_id ) 
        CALL mp_bcast( kkbeta_(1:nbeta_), ionode_id ) 
        CALL mp_bcast( beta_(0:mesh_,1:nbeta_), ionode_id ) 
        CALL mp_bcast( dion_(1:nbeta_,1:nbeta_), ionode_id ) 
        CALL mp_bcast( rinner_(1:nqlc_), ionode_id ) 
        CALL mp_bcast( qqq_(1:nbeta_,1:nbeta_), ionode_id ) 
        CALL mp_bcast( qfunc_(0:mesh_, 1:nbeta_, 1:nbeta_), ionode_id )
        CALL mp_bcast( qfcoef_(1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_), ionode_id ) 
        CALL mp_bcast( chi_(0:mesh_, nwfc_), ionode_id ) 
        CALL mp_bcast( rho_at_(0:mesh_), ionode_id )
!
!
        generated   = generated_
        date_author = date_author_
        comment     = comment_
        psd         = psd_
        typ         = typ_
        tvanp       = tvanp_
        nlcc        = nlcc_
        dft         = dft_
        zp          = zp_
        etotps      = etotps_
        ecutwfc     = ecutwfc_
        ecutrho     = ecutrho_
        nv          = nv_
        lmax        = lmax_
        mesh        = mesh_
        nwfc        = nwfc_
        nbeta       = nbeta_
        nd          = nd_
        nqf         = nqf_
        nqlc        = nqlc_

        IF( tovrw ) THEN

          els(1:nwfc_) = els_(1:nwfc_)
          lchi(nwfc_)  = lchi_(nwfc_)
          oc(nwfc_)    = oc_(nwfc_)
          r(0:mesh_)   = r_(0:mesh_)
          rab(0:mesh_) = rab_(0:mesh_)
          rho_atc(0:mesh_) = rho_atc_(0:mesh_)
          vloc(0:mesh_) = vloc_(0:mesh_)
          lll(1:nbeta_) = lll_(1:nbeta_)
          kkbeta(1:nbeta_) = kkbeta_(1:nbeta_)
          beta(0:mesh_,1:nbeta_) = beta_(0:mesh_,1:nbeta_)
          dion(1:nbeta_,1:nbeta_) = dion_(1:nbeta_,1:nbeta_)
          rinner(1:nqlc_) = rinner_(1:nqlc_)
          qqq(1:nbeta_,1:nbeta_) = qqq_(1:nbeta_,1:nbeta_)
          qfunc(0:mesh_, 1:nbeta_, 1:nbeta_) = qfunc_(0:mesh_, 1:nbeta_, 1:nbeta_)
          qfcoef(1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_) = qfcoef_(1:nqf_, 1:nqlc_, 1:nbeta_, 1:nbeta_)
          chi(0:mesh_, nwfc_) = chi_(0:mesh_, nwfc_)
          rho_at(0:mesh_) = rho_at_(0:mesh_)

        END IF

        DEALLOCATE( els_,  lchi_, oc_, r_, rab_, rho_atc_, vloc_, lll_, kkbeta_, &
            beta_, dion_, rinner_, qqq_, qfunc_, qfcoef_, chi_, rho_at_ )

      ELSE

        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum

        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_pseudo, pseudo not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_pseudo2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum 
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_pseudo, pseudo not read from restart ' )")
      RETURN
    END SUBROUTINE


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
    SUBROUTINE write_restart_gvec1(iuni, twrite, &
      ng, bi1, bi2, bi3, b1, b2, b3, tmill, mill )

      USE io_global, ONLY: ionode

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ng
      REAL(dbl), INTENT(IN) :: bi1(3), bi2(3), bi3(3)
      REAL(dbl), INTENT(IN) :: b1(3), b2(3), b3(3)
      INTEGER, INTENT(IN) :: mill(:,:)
      LOGICAL, INTENT(IN) :: twrite, tmill
      INTEGER :: idum
      INTEGER :: i, j
      CHARACTER(LEN=30) :: sub_name = ' write_restart_gvec '
!
! ... Subroutine Body
!
      IF( twrite ) THEN

        IF( tmill ) THEN
          IF( ( SIZE(mill,1) < 3 ) .OR. ( SIZE(mill,2) < ng) ) &
            CALL error( sub_name, ' wrong size ', 1 )
        END IF

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) ng, tmill
          WRITE(iuni) bi1, bi2, bi3, b1, b2, b3
          IF( tmill ) THEN
            WRITE(iuni) ((mill(i,j),i=1,3),j=1,ng)
          ELSE
            WRITE(iuni) idum
          END IF
        END IF
      ELSE
        CALL write_restart_gvec2(iuni)
      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_gvec2(iuni)
      USE io_global, ONLY: ionode
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum 
      END IF
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!


    SUBROUTINE read_restart_gvec1(iuni, tovrw, tread, &
      ng, bi1, bi2, bi3, b1, b2, b3, tmill, mill )
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast

! .. Subroutine output:
!    if tread is true then variables are read from file but values are
!      not copied on output variables
!    if tmill is true "mill" array is read from file
!    if tovrw is true all variables other than "mill" are overvritten with values read from file
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      LOGICAL, INTENT(INOUT) :: tmill
      INTEGER, INTENT(INOUT) :: ng
      REAL(dbl), INTENT(INOUT) :: b1(3), b2(3), b3(3)
      REAL(dbl), INTENT(INOUT) :: bi1(3), bi2(3), bi3(3)
      INTEGER, INTENT(INOUT) :: mill(:,:)

      INTEGER :: i, j
      INTEGER :: ng_
      REAL(dbl) :: b1_(3), b2_(3), b3_(3)
      REAL(dbl) :: bi1_(3), bi2_(3), bi3_(3)
      INTEGER, ALLOCATABLE :: mill_(:,:)
      LOGICAL :: twrite_, tmill_
      INTEGER :: file_version_
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_gvec '
!
! ... Subroutine Body
!


      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version_, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(sub_name, ' Data Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) ng_, tmill_
          READ(iuni) bi1_, bi2_, bi3_, b1_, b2_, b3_
        END IF
        CALL mp_bcast( ng_, ionode_id )
        CALL mp_bcast( tmill_, ionode_id )
        CALL mp_bcast( bi1_, ionode_id )
        CALL mp_bcast( bi2_, ionode_id )
        CALL mp_bcast( bi3_, ionode_id )
        CALL mp_bcast( b1_, ionode_id )
        CALL mp_bcast( b2_, ionode_id )
        CALL mp_bcast( b3_, ionode_id )

        IF( tmill .AND. .NOT. tmill_ ) &
          CALL error(sub_name, ' mill indexes not present in restart file ', 1)

        IF( tmill ) THEN

          IF( ( SIZE( mill, 2) < ng_ ) .OR. ( SIZE( mill, 1 ) < 3 ) ) &
            CALL error(sub_name, ' mill array too small ', 1)

          IF( ionode ) THEN
            READ(iuni) ((mill(i,j),i=1,3),j=1,ng_)
          END IF
          CALL mp_bcast( mill, ionode_id )

        ELSE

          IF( ionode ) THEN
            READ(iuni) idum
          END IF

        END IF
  
        IF( .NOT. tovrw ) THEN

          IF( ng_ /= ng ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, ng changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") ng_, ng
          END IF
          IF( ANY( bi1_ /= bi1 ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, bi1 changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") bi1_, bi1
          END IF
          IF( ANY( bi2_ /= bi2 ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, bi2 changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") bi2_, bi2
          END IF
          IF( ANY( bi3_ /= bi3 ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, bi3 changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") bi3_, bi3
          END IF
          IF( ANY( b1_ /= b1 ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, b1 changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") b1_, b1
          END IF
          IF( ANY( b2_ /= b2 ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, b2 changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") b2_, b2
          END IF
          IF( ANY( b3_ /= b3 ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gvec, b3 changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") b3_, b3
          END IF
        END IF

        IF( tovrw ) THEN
          ng  = ng_
          bi1 = bi1_
          bi2 = bi2_
          bi3 = bi3_
          b1  = b1_
          b2  = b2_
          b3  = b3_
        END IF

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum
          READ(iuni) idum
          READ(iuni) idum
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_gvec, data not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_gvec2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_gvec, data not read from restart ' )")
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variables related to a single k point
! .. Where:
! iuni    = Restart file I/O fortran unit
! ngwk    = number of wavefunctions G vectors, for this k point 
! igk(.)  = for each G+k, igk is the index of the G vectors
! xk(.)   = k point coordinate
! wk      = k point weight

    SUBROUTINE write_restart_gkvec1(iuni, twrite, &
      ik, nk, ngwk, xk, wk, tetra, isk)

      USE io_global, ONLY: ionode

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      ! ... INTEGER, INTENT(IN) :: igk(:)
      INTEGER, INTENT(IN) :: ik, nk, ngwk, tetra(4), isk
      REAL(dbl), INTENT(IN) :: xk(3)
      REAL(dbl), INTENT(IN) :: wk
      INTEGER :: idum, i

      IF( twrite ) THEN
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) ik, nk, ngwk, tetra(1:4), isk
          WRITE(iuni) (xk(i),i=1,3), wk
          WRITE(iuni) idum ! (igk(i),i=1,ngwk)
        END IF
      ELSE
        CALL write_restart_gkvec2(iuni)
      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_gkvec2(iuni)
      USE io_global, ONLY: ionode
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( twrite ) THEN
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) idum
          WRITE(iuni) idum
          WRITE(iuni) idum
        END IF
      END IF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_gkvec1(iuni, tovrw, tread, &
      ik, nk, ngwk, xk, wk, tetra, isk )

! .. Subroutine output:
!    if tread is true then variables are read from file but values are
!      not copied on output variables
!    if tovrw is true all variables are overvritten with values read from file
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      ! ... INTEGER, INTENT(INOUT) :: igk(:)
      INTEGER, INTENT(INOUT) :: ngwk, ik, nk, tetra(4), isk
      REAL(dbl), INTENT(INOUT) :: xk(3)
      REAL(dbl), INTENT(INOUT) :: wk

      INTEGER :: i
      INTEGER :: ngwk_, ik_, nk_, tetra_(4), isk_
      REAL(dbl) :: xk_(3)
      REAL(dbl) :: wk_
      INTEGER, ALLOCATABLE :: igk_(:)

      INTEGER :: idum, nigk
      LOGICAL :: twrite_
      INTEGER :: file_version_
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version_, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_gkvec ',' Data Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) ik_, nk_, ngwk_, tetra_(1:4), isk_
          READ(iuni) (xk_(i),i=1,3), wk_
        END IF
        CALL mp_bcast( ngwk_, ionode_id )
        CALL mp_bcast( tetra_, ionode_id )
        CALL mp_bcast( ik_, ionode_id )
        CALL mp_bcast( isk_, ionode_id )
        CALL mp_bcast( nk_, ionode_id )
        CALL mp_bcast( xk_, ionode_id )
        CALL mp_bcast( wk_, ionode_id )

        ! .. ALLOCATE( igk_( ngwk_ ) )

        IF( ionode ) THEN
          READ(iuni) idum ! .. (igk_(i),i=1,ngwk_)
        END IF
        ! .. CALL mp_bcast( igk_, ionode_id )
 
        IF( .NOT. tovrw ) THEN
          IF( ngwk_ /= ngwk ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gkvec, ngwk changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") ngwk_, ngwk
          END IF
          IF( ANY( xk_ /= xk ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gkvec, xk changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") xk_, xk
          END IF
          IF( wk_ /= wk ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_gkvec, wk changed' )")
            WRITE(6,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") wk_, wk
          END IF
        END IF

        IF( tovrw ) THEN
          ngwk = ngwk_
          ik   = ik_
          nk   = nk_
          tetra = tetra_
          isk  = isk_
          xk   = xk_
          wk   = wk_
          ! ... nigk = MIN( SIZE(igk), ngwk_ )
          ! ... igk(1:nigk) = igk_(1:nigk)
        END IF

        ! ... DEALLOCATE( igk_ )

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum
          READ(iuni) idum
          READ(iuni) idum
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_gkvec, xdim not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_gkvec2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_gkvec, xdim not read from restart ' )")
      RETURN
    END SUBROUTINE


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
    SUBROUTINE write_restart_cell1( iuni, twrite, &
      ibrav, celldm, ht0, htm, htm2, htvel, xnosp, xnos0, xnosm, xnosm2)
      USE io_global, ONLY: ionode
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      INTEGER, INTENT(IN) :: ibrav
      REAL(dbl), INTENT(IN) :: celldm(6)
      REAL(dbl), INTENT(IN) :: ht0(3,3)
      REAL(dbl), INTENT(IN) :: htm(3,3)
      REAL(dbl), INTENT(IN) :: htm2(3,3)
      REAL(dbl), INTENT(OUT) :: htvel(3,3)
      REAL(dbl), INTENT(OUT) :: xnosp(3,3)
      REAL(dbl), INTENT(OUT) :: xnos0(3,3)
      REAL(dbl), INTENT(OUT) :: xnosm(3,3)
      REAL(dbl), INTENT(OUT) :: xnosm2(3,3)
      INTEGER :: i
!
! ... Subroutine Body
!
      IF( twrite ) THEN
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
          WRITE(iuni) ibrav, (celldm(i), i=1,6)
          WRITE(iuni) ht0, htm, htm2, htvel
          WRITE(iuni) xnosp, xnos0, xnosm, xnosm2
        END IF
      ELSE
        CALL  write_restart_cell2( iuni )
      END IF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_cell2(iuni) 
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_cell1( iuni, tovrw, tread, &
      ibrav, celldm, ht0, htm, htm2, htvel, xnosp, xnos0, xnosm, xnosm2)

! .. Subroutine output:
!    if tread is true then the following variables are overwritten on output
!      ht0, htm, htm2, htvel, xnosp, xnos0, xnosm, xnosm2
!    if tovrw is true all variables are overvritten with values read from file
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      INTEGER, INTENT(INOUT) :: ibrav
      REAL(dbl), INTENT(INOUT) :: celldm(6)
      REAL(dbl), INTENT(INOUT) :: ht0(3,3)
      REAL(dbl), INTENT(INOUT) :: htm(3,3)
      REAL(dbl), INTENT(INOUT) :: htm2(3,3)
      REAL(dbl), INTENT(INOUT) :: htvel(3,3)
      REAL(dbl), INTENT(INOUT) :: xnosp(3,3)
      REAL(dbl), INTENT(INOUT) :: xnos0(3,3)
      REAL(dbl), INTENT(INOUT) :: xnosm(3,3)
      REAL(dbl), INTENT(INOUT) :: xnosm2(3,3)

      INTEGER :: i
      INTEGER :: ibrav_
      REAL(dbl) :: ht0_(3,3), htm_(3,3), htm2_(3,3), htvel_(3,3)
      REAL(dbl) :: xnosp_(3,3), xnos0_(3,3), xnosm_(3,3), xnosm2_(3,3)
      REAL(dbl) :: celldm_(6)
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version, ionode_id )
      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_cell ', ' Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) ibrav_, (celldm_(i), i=1,6)
          READ(iuni) ht0_, htm_, htm2_, htvel_
          READ(iuni) xnosp_, xnos0_, xnosm_, xnosm2_
        END IF
        CALL mp_bcast( ibrav_, ionode_id )
        CALL mp_bcast( celldm_, ionode_id )
        CALL mp_bcast( ht0_, ionode_id )
        CALL mp_bcast( htm_, ionode_id )
        CALL mp_bcast( htm2_, ionode_id )
        CALL mp_bcast( htvel_, ionode_id )
        CALL mp_bcast( xnosp_, ionode_id )
        CALL mp_bcast( xnos0_, ionode_id )
        CALL mp_bcast( xnosm_, ionode_id )
        CALL mp_bcast( xnosm2_, ionode_id )

        IF( .NOT. tovrw ) THEN
          IF( ibrav_ /= ibrav ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_cell, ibrav changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") ibrav_, ibrav
          END IF
          IF( ANY( celldm_ /= celldm ) ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_cell, celldm changed' )")
            WRITE(6,fmt="(3X,'W: old = ',6F14.8 )") celldm_
            WRITE(6,fmt="(3X,'W: new = ',6F14.8 )") celldm
          END IF
        END IF

        ht0 = ht0_
        htm = htm_
        htm2 = htm2_
        htvel = htvel_
        xnosp = xnosp_
        xnos0 = xnos0_
        xnosm = xnosm_
        xnosm2 = xnosm2_

        IF( tovrw ) THEN
          ibrav = ibrav_
          celldm = celldm_
        END IF

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum
          READ(iuni) idum
          READ(iuni) idum
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_cell, xdim not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_cell2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_cell, xdim not read from restart ' )")
      RETURN
    END SUBROUTINE




!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variable related to the ion types 
! ..  positions, and velocities
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_ions1(iuni, twrite, &
      label, tscal, stau0, svel0, staum, svelm, taui, fion, &
      cdmi, nat, ntyp, ityp, na, mass, xnosp, xnos0, xnosm, xnosm2)
!
      USE io_global, ONLY: ionode
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      LOGICAL, INTENT(IN) :: tscal
      CHARACTER(LEN=*), INTENT(IN) :: label(:)
      REAL(dbl), INTENT(IN) :: stau0(:,:)
      REAL(dbl), INTENT(IN) :: svel0(:,:)
      REAL(dbl), INTENT(IN) :: staum(:,:)
      REAL(dbl), INTENT(IN) :: svelm(:,:)
      REAL(dbl), INTENT(IN) :: taui(:,:)
      REAL(dbl), INTENT(IN) :: fion(:,:)
      REAL(dbl), INTENT(IN) :: cdmi(:)
      INTEGER, INTENT(IN) :: nat
      INTEGER, INTENT(IN) :: ntyp
      INTEGER, INTENT(IN) :: ityp(:)
      INTEGER, INTENT(IN) :: na(:)
      REAL(dbl), INTENT(IN) :: mass(:)
      REAL(dbl), INTENT(IN) :: xnosp
      REAL(dbl), INTENT(IN) :: xnos0
      REAL(dbl), INTENT(IN) :: xnosm
      REAL(dbl), INTENT(IN) :: xnosm2

      INTEGER :: i,j
      CHARACTER(LEN=4) :: label_(ntyp)
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' write_restart_ions '
!
! ... Subroutine Body
!

      IF( twrite ) THEN

        IF( SIZE( label ) < ntyp ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE( ityp ) < nat ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( SIZE( na ) < ntyp ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( SIZE( mass ) < ntyp ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( ( SIZE( stau0, 1 ) < 3 ) .OR. ( SIZE( stau0, 2 ) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( ( SIZE( svel0, 1 ) < 3 ) .OR. ( SIZE( svel0, 2 ) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF( ( SIZE( staum, 1 ) < 3 ) .OR. ( SIZE( staum, 2 ) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 7 )
        IF( ( SIZE( svelm, 1 ) < 3 ) .OR. ( SIZE( svelm, 2 ) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 8 )
        IF( ( SIZE( taui, 1 ) < 3 ) .OR. ( SIZE( taui, 2 ) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 9 )
        IF( ( SIZE( fion, 1 ) < 3 ) .OR. ( SIZE( fion, 2 ) < nat ) ) &
          CALL error( sub_name, ' wrong size ', 10 )
        IF( SIZE( cdmi ) < 3 ) &
          CALL error( sub_name, ' wrong size ', 11 )

        label_ = label(1:ntyp)

        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
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

      ELSE

        CALL write_restart_ions2(iuni)

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_ions2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
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
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_ions1(iuni, tovrw, tread, &
      label, tscal, stau0, svel0, staum, svelm, taui, fion, &
      cdmi, nat, ntyp, ityp, na, mass, xnosp, xnos0, xnosm, xnosm2)

! .. Subroutine output:
!    if tread is true then the following variables are overwritten on output
!      tscal, stau0, svel0, staum, svelm, taui, fion, cdmi, xnosp, xnos0, xnosm, xnosm2
!    if tovrw is true all variables are overvritten with values read from file
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      LOGICAL, INTENT(OUT) :: tscal
      CHARACTER(LEN=*), INTENT(INOUT) :: label(:)
      REAL(dbl), INTENT(OUT) :: stau0(:,:)
      REAL(dbl), INTENT(OUT) :: svel0(:,:)
      REAL(dbl), INTENT(OUT) :: staum(:,:)
      REAL(dbl), INTENT(OUT) :: svelm(:,:)
      REAL(dbl), INTENT(OUT) :: taui(:,:)
      REAL(dbl), INTENT(OUT) :: fion(:,:)
      REAL(dbl), INTENT(OUT) :: cdmi(:)
      INTEGER, INTENT(INOUT) :: nat
      INTEGER, INTENT(INOUT) :: ntyp
      INTEGER, INTENT(INOUT) :: ityp(:)
      INTEGER, INTENT(INOUT) :: na(:)
      REAL(dbl), INTENT(INOUT) :: mass(:)
      REAL(dbl), INTENT(OUT) :: xnosp
      REAL(dbl), INTENT(OUT) :: xnos0
      REAL(dbl), INTENT(OUT) :: xnosm
      REAL(dbl), INTENT(OUT) :: xnosm2

      INTEGER   :: i, j
      INTEGER   :: nat_, ntyp_, ntau
      INTEGER   :: na_(nsx)
      INTEGER   :: ityp_(natx)
      REAL(dbl) :: mass_(nsx)
      REAL(dbl) :: cdmi_(3)
      REAL(dbl) :: xnosp_, xnos0_, xnosm_, xnosm2_
      CHARACTER(LEN=4) :: label_(nsx)

      REAL(dbl), ALLOCATABLE :: stmp_(:,:)
      LOGICAL :: twrite_, tscal_
      INTEGER :: file_version_
      INTEGER :: idum
      CHARACTER(LEN=30) :: sub_name = ' read_restart_ions '
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_ions ',' Data Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) THEN
          READ(iuni) nat_, ntyp_, tscal_
        END IF
        CALL mp_bcast(nat_, ionode_id)
        CALL mp_bcast(ntyp_, ionode_id)
        CALL mp_bcast(tscal_, ionode_id)
        tscal = tscal_

        IF( ntyp_ > SIZE(na_) ) &
            CALL error(' read_restart_ions ',' too many types ',ntyp_)
        IF( nat_ > SIZE(ityp_) ) &
          CALL error(' read_restart_ions ',' too many atoms ',nat_)

        IF( .NOT. tovrw ) THEN
          IF( nat_ /= nat ) &
            CALL error(' read_restart_ions ',' atoms numbers differ ',nat_)
          IF( ntyp_ /= ntyp ) &
            CALL error(' read_restart_ions ',' types numbers differ ',ntyp_)
        END IF
       
        IF( ionode ) THEN
          READ(iuni) (ityp_(i),i=1,nat_), (na_(i),i=1,ntyp_), (label_(i),i=1,ntyp_)
          READ(iuni) (mass_(i),i=1,ntyp_)
        END IF
        CALL mp_bcast(ityp_, ionode_id)
        CALL mp_bcast(na_, ionode_id)
        CALL mp_bcast(label_, ionode_id)
        CALL mp_bcast(mass_, ionode_id)
 
        IF( .NOT. tovrw ) THEN
          DO i = 1, ntyp_
            IF( na_(i) /= na(i) ) THEN
              WRITE(6,fmt="(3X,'W: read_restart_ions, Number of atoms chenged' )")
              WRITE(6,fmt="(3X,'W: is = ',I3,' old = ',I10,' new = ',I10 )") i, na_(i), na(i)
            END IF
            IF( mass_(i) /= mass(i) ) THEN
              WRITE(6,fmt="(3X,'W: read_restart_ions, Atomic mass chenged' )")
              WRITE(6,fmt="(3X,'W: is = ',I3,' old = ',F18.8,' new = ',F18.8 )") i, mass_(i), mass(i)
            END IF
          END DO
        END IF

        IF( SIZE( label ) < ntyp_ ) &
          CALL error( sub_name, ' wrong size ', 1 )
        IF( SIZE( ityp ) < nat_ ) &
          CALL error( sub_name, ' wrong size ', 2 )
        IF( SIZE( na ) < ntyp_ ) &
          CALL error( sub_name, ' wrong size ', 3 )
        IF( SIZE( mass ) < ntyp_ ) &
          CALL error( sub_name, ' wrong size ', 4 )
        IF( ( SIZE( stau0, 1 ) < 3 ) .OR. ( SIZE( stau0, 2 ) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 5 )
        IF( ( SIZE( svel0, 1 ) < 3 ) .OR. ( SIZE( svel0, 2 ) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 6 )
        IF( ( SIZE( staum, 1 ) < 3 ) .OR. ( SIZE( staum, 2 ) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 7 )
        IF( ( SIZE( svelm, 1 ) < 3 ) .OR. ( SIZE( svelm, 2 ) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 8 )
        IF( ( SIZE( taui, 1 ) < 3 ) .OR. ( SIZE( taui, 2 ) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 9 )
        IF( ( SIZE( fion, 1 ) < 3 ) .OR. ( SIZE( fion, 2 ) < nat_ ) ) &
          CALL error( sub_name, ' wrong size ', 10 )
        IF( SIZE( cdmi ) < 3 ) &
          CALL error( sub_name, ' wrong size ', 11 )

        ALLOCATE( stmp_( 3, nat_ ) )
        IF( ionode ) THEN
          READ(iuni) ((stmp_(i,j),i=1,3),j=1,nat_)
          stau0 = stmp_
          READ(iuni) ((stmp_(i,j),i=1,3),j=1,nat_)
          svel0 = stmp_
          READ(iuni) ((stmp_(i,j),i=1,3),j=1,nat_)
          staum = stmp_
          READ(iuni) ((stmp_(i,j),i=1,3),j=1,nat_)
          svelm = stmp_
          READ(iuni) ((stmp_(i,j),i=1,3),j=1,nat_)
          taui = stmp_
          READ(iuni) ((stmp_(i,j),i=1,3),j=1,nat_)
          fion = stmp_
        END IF
        CALL mp_bcast(stau0, ionode_id)
        CALL mp_bcast(svel0, ionode_id)
        CALL mp_bcast(staum, ionode_id)
        CALL mp_bcast(svelm, ionode_id)
        CALL mp_bcast(taui, ionode_id)
        CALL mp_bcast(fion, ionode_id)
        DEALLOCATE( stmp_ )

        IF( ionode ) THEN
          READ(iuni) (cdmi_(i),i=1,3)
          READ(iuni) xnosp_, xnos0_, xnosm_, xnosm2_
        END IF
        CALL mp_bcast(cdmi_, ionode_id)
        CALL mp_bcast(xnosp_, ionode_id)
        CALL mp_bcast(xnos0_, ionode_id)
        CALL mp_bcast(xnosm_, ionode_id)
        CALL mp_bcast(xnosm2_, ionode_id)
  
        cdmi = cdmi_
        xnosp = xnosp_
        xnos0 = xnos0_
        xnosm = xnosm_
        xnosm2 = xnosm2_

        IF( tovrw ) THEN
          nat = nat_
          ntyp = ntyp_
          mass(1:ntyp) = mass_
          na(1:ntyp) = na_
          ityp(1:nat) = ityp_(1:nat)
          label(1:ntyp) = label_
        END IF

      ELSE

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
          WRITE(6,fmt="(3X,'W: read_restart_ions, Data Section not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_ions2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
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
        WRITE(6,fmt="(3X,'W: read_restart_ions, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

! ..  This subroutine write to disk variable related to electronic band
! ..  structure (NOT the wavefunctions)
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_electrons1( iuni, twrite, &
      occ, occm, tocc, lambda, lambdam, ldim, tlam, nbnd, ispin, nspin, ik, nk, nel, nelu, &
      neld, xnosp, xnos0, xnosm, xnosm2, ef, teig, eig, weig)
!
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      REAL(dbl), INTENT(IN) :: occ(:)
      REAL(dbl), INTENT(IN) :: occm(:)
      REAL(dbl), INTENT(IN) :: lambda(:,:)
      REAL(dbl), INTENT(IN) :: lambdam(:,:)
      REAL(dbl), INTENT(IN) :: eig(:)
      REAL(dbl), INTENT(IN) :: weig(:)
      LOGICAL, INTENT(IN) :: tocc, tlam, teig
      INTEGER, INTENT(IN) :: nbnd, ldim
      INTEGER, INTENT(IN) :: ispin
      INTEGER, INTENT(IN) :: nspin
      INTEGER, INTENT(IN) :: ik
      INTEGER, INTENT(IN) :: nk
      REAL(dbl), INTENT(IN) :: nel
      INTEGER, INTENT(IN) :: nelu
      INTEGER, INTENT(IN) :: neld
      REAL(dbl), INTENT(IN) :: xnosp
      REAL(dbl), INTENT(IN) :: xnos0
      REAL(dbl), INTENT(IN) :: xnosm
      REAL(dbl), INTENT(IN) :: xnosm2
      REAL(dbl), INTENT(IN) :: ef

      INTEGER :: i, l, idum
      CHARACTER(LEN=30) :: sub_name = ' write_restart_electrons '
!
! ... Subroutine Body
!

      IF( twrite ) THEN

        IF( ionode ) WRITE(iuni) twrite, file_version
        IF( ionode ) WRITE(iuni) nbnd, ispin, nspin, ik, nk, nel, nelu, neld, ldim

        IF( ionode ) WRITE(iuni) tocc

        IF( tocc ) THEN
         
          IF( SIZE( occ ) < nbnd ) &
            CALL error(sub_name, ' wrong size ', 1 )
          IF( SIZE( occm ) < nbnd ) &
            CALL error(sub_name, ' wrong size ', 2 )
          
          IF( ionode ) WRITE(iuni) (occ(i),i=1,nbnd)
          IF( ionode ) WRITE(iuni) (occm(i),i=1,nbnd)

        ELSE

          IF( ionode ) WRITE(iuni) idum
          IF( ionode ) WRITE(iuni) idum

        END IF

        IF( ionode ) WRITE(iuni) tlam

        IF( tlam ) THEN

          IF( ( SIZE( lambda, 1 ) < ldim ) .OR. ( SIZE( lambda, 2 ) < ldim ) ) &
            CALL error(sub_name, ' wrong size ', 3 )
          IF( ( SIZE( lambdam, 1 ) < ldim ) .OR. ( SIZE( lambdam, 2 ) < ldim ) ) &
            CALL error(sub_name, ' wrong size ', 4 )

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
            CALL error(sub_name, ' wrong size ', 5 )
          IF( SIZE( weig ) < nbnd ) &
            CALL error(sub_name, ' wrong size ', 6 )

          IF( ionode ) WRITE(iuni) (eig(i),i=1,nbnd)
          IF( ionode ) WRITE(iuni) (weig(i),i=1,nbnd)

        ELSE

          IF( ionode ) WRITE(iuni) idum
          IF( ionode ) WRITE(iuni) idum

        END IF

      ELSE

        CALL write_restart_electrons2( iuni )

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

   SUBROUTINE write_restart_electrons2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
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
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_electrons1( iuni, tovrw, tread, &
      occ, occm, tocc, lambda, lambdam, ldim, tlam, nbnd, ispin, nspin, ik, nk, nel, nelu, &
      neld, xnosp, xnos0, xnosm, xnosm2, ef, teig, eig, weig)

! .. Subroutine output:
!    if tread is true then variables xnosp, xnos0, xnosm, xnosm2 are overwritten with 
!      values read from file
!    if tocc is true then "occ" and "occm" are overwritten
!    if tlam is true then "lambda" and "lambdam" are overwritten
!    if teig is true then "eig" and "weig" are overwritten
!    if tovrw is true all variables but occ, occm, lambda, lambdam, eig, weig,
!      are overvritten with values read from file
!
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      REAL(dbl), INTENT(OUT) :: occ(:)
      REAL(dbl), INTENT(OUT) :: occm(:)
      REAL(dbl), INTENT(OUT) :: eig(:)
      REAL(dbl), INTENT(OUT) :: weig(:)
      REAL(dbl), INTENT(OUT) :: lambda(:,:)
      REAL(dbl), INTENT(OUT) :: lambdam(:,:)
      INTEGER, INTENT(INOUT) :: ldim
      INTEGER, INTENT(INOUT) :: nbnd
      INTEGER, INTENT(INOUT) :: ispin
      INTEGER, INTENT(INOUT) :: nspin
      INTEGER, INTENT(INOUT) :: ik
      INTEGER, INTENT(INOUT) :: nk
      REAL(dbl), INTENT(INOUT) :: nel
      INTEGER, INTENT(INOUT) :: nelu
      INTEGER, INTENT(INOUT) :: neld
      REAL(dbl), INTENT(OUT) :: xnosp
      REAL(dbl), INTENT(OUT) :: xnos0
      REAL(dbl), INTENT(OUT) :: xnosm
      REAL(dbl), INTENT(OUT) :: xnosm2
      LOGICAL, INTENT(IN) :: tocc, tlam, teig
      REAL(dbl), INTENT(INOUT) :: ef

      INTEGER :: i,j,k,l, nsiz
      LOGICAL :: tocc_, tlam_, teig_
      REAL(dbl) :: nel_
      INTEGER :: nbnd_, ispin_, nspin_, ik_, nk_, nelu_, neld_, ldim_
      REAL(dbl) :: xnosp_, xnos0_, xnosm_, xnosm2_, ef_
      REAL(dbl), ALLOCATABLE :: otmp_(:)
      REAL(dbl), ALLOCATABLE :: ltmp_(:,:)
      INTEGER :: idum

      LOGICAL :: twrite_
      INTEGER :: file_version_
      CHARACTER(LEN=30) :: sub_name = ' read_restart_electrons '
!
! ... Subroutine Body
!
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version, ionode_id )
      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_electrons ',' Data Section not present in restart file ', 1)

      IF( tread ) THEN
        IF( ionode ) READ(iuni) nbnd_, ispin_, nspin_, ik_, nk_, nel_, nelu_, neld_, ldim_
        CALL mp_bcast( nbnd_, ionode_id)
        CALL mp_bcast( ispin_, ionode_id)
        CALL mp_bcast( nspin_, ionode_id)
        CALL mp_bcast( ik_, ionode_id)
        CALL mp_bcast( nk_, ionode_id)
        CALL mp_bcast( nel_, ionode_id)
        CALL mp_bcast( nelu_, ionode_id)
        CALL mp_bcast( neld_, ionode_id)
        CALL mp_bcast( ldim_, ionode_id)

        IF( .NOT. tovrw ) THEN
          IF( nbnd_ /= nbnd ) &
            CALL error( ' read_restart_electrons ',' number of bands differ ', 1)
          IF( nspin_ /= nspin ) &
            CALL error( ' read_restart_electrons ',' number of spin differ ', 1)
          IF( nk_ /= nk ) &
            CALL error( ' read_restart_electrons ',' number of k points differ ', 1)
        END IF
!
! ..    Manage occ and occm

        IF( ionode ) READ(iuni) tocc_
        CALL mp_bcast( tocc_, ionode_id)

        IF( tocc .AND. .NOT. tocc_ ) &
          CALL error( ' read_restart_electrons ',' occupation number not present in restart ', 1)

        IF( tocc ) THEN

          nsiz = MIN( nbnd_, SIZE( occ ) )
          ALLOCATE( otmp_( nbnd_) )

          IF( ionode ) READ(iuni) (otmp_(i),i=1,nbnd_)
          CALL mp_bcast( otmp_, ionode_id)

          occ(1:nsiz) = otmp_(1:nsiz)

          IF( ionode ) READ(iuni) (otmp_(i),i=1,nbnd_)
          CALL mp_bcast( otmp_, ionode_id)

          occm(1:nsiz) = otmp_(1:nsiz)

          DEALLOCATE( otmp_ )

        ELSE

          IF( ionode ) READ(iuni) idum
          IF( ionode ) READ(iuni) idum

        END IF

! ..    Manage lambda and lambdam

        IF( ionode ) READ(iuni) tlam_
        CALL mp_bcast( tlam_, ionode_id)

        IF( tlam .AND. .NOT. tlam_ ) &
          CALL error( ' read_restart_electrons ',' lambda matrix not present in restart ', 1)

        IF( tlam ) THEN

          nsiz = MIN( ldim_, SIZE( lambda, 1 ) )
          ALLOCATE( ltmp_( ldim_, ldim_ ) )

          IF( ionode ) READ(iuni) ((ltmp_(l,i),l=1,ldim_),i=1,ldim_)
          CALL mp_bcast( ltmp_, ionode_id)

          lambda(1:nsiz, 1:nsiz) = ltmp_(1:nsiz, 1:nsiz)

          IF( ionode ) READ(iuni) ((ltmp_(l,i),l=1,ldim_),i=1,ldim_)
          CALL mp_bcast( ltmp_, ionode_id)

          lambdam(1:nsiz, 1:nsiz) = ltmp_(1:nsiz, 1:nsiz)

          DEALLOCATE( ltmp_ )

        ELSE

          IF( ionode ) READ(iuni) idum
          IF( ionode ) READ(iuni) idum

        END IF


        IF( ionode ) READ(iuni) xnosp_, xnos0_, xnosm_, xnosm2_
        CALL mp_bcast( xnosp_, ionode_id)
        CALL mp_bcast( xnos0_, ionode_id)
        CALL mp_bcast( xnosm_, ionode_id)
        CALL mp_bcast( xnosm2_, ionode_id)
  
        IF( ionode ) READ(iuni) ef_
        CALL mp_bcast( ef_, ionode_id)

        xnosp = xnosp_
        xnos0 = xnos0_
        xnosm = xnosm_
        xnosm2 = xnosm2_

        IF( ionode ) READ(iuni) teig_
        CALL mp_bcast( teig_, ionode_id)

        IF( teig .AND. .NOT. teig_ ) &
          CALL error( ' read_restart_electrons ',' occupation number not present in restart ', 1)

        IF( teig ) THEN

          nsiz = MIN( nbnd_, SIZE( eig ) )
          ALLOCATE( otmp_( nbnd_ ) )

          IF( ionode ) READ(iuni) (otmp_(i),i=1,nbnd_)
          CALL mp_bcast( otmp_, ionode_id)

          eig(1:nsiz) = otmp_(1:nsiz)

          IF( ionode ) READ(iuni) (otmp_(i),i=1,nbnd_)
          CALL mp_bcast( otmp_, ionode_id)

          weig(1:nsiz) = otmp_(1:nsiz)

          DEALLOCATE( otmp_ )

        ELSE

          IF( ionode ) READ(iuni) idum
          IF( ionode ) READ(iuni) idum

        END IF
  
        IF( tovrw ) THEN
          nbnd = nbnd_
          ispin = ispin_
          nspin = nspin_
          ik    = ik_
          nk    = nk_ 
          nel   = nel_
          nelu  = nelu_
          neld  = neld_
          ldim  = ldim_
          ef    = ef_
        END IF

      ELSE

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
          WRITE(6,fmt="(3X,'W: read_restart_electrons, Data Sections not read from restart ' )")

      ENDIF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_electrons2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
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
        WRITE(6,fmt="(3X,'W: read_restart_electrons, Data Sections not read from restart ' )")
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

! ..  This subroutine write wavefunctions to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_wfc1(iuni, twrite, &
      ik, nk, kunit, ispin, nspin, scal, wf0, t0, wfm, tm, ngw, nbnd, igl, ngwl )
!
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_get, mp_bcast, mp_max
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      INTEGER, INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(dbl), INTENT(IN) :: wf0(:,:)
      COMPLEX(dbl), INTENT(IN) :: wfm(:,:)
      INTEGER, INTENT(IN) :: ngw   ! 
      INTEGER, INTENT(IN) :: nbnd
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)
      REAL(dbl), INTENT(IN) :: scal
      LOGICAL, INTENT(IN) :: t0, tm

      INTEGER :: i, j, idum, ierr
      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx
      INTEGER :: npool, ipmask( nproc ), ipsour
      COMPLEX(dbl), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)
!
! ... Subroutine Body
!

      IF( twrite ) THEN

        IF( ionode ) WRITE(iuni) twrite, file_version


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
          CALL error(' write_restart_wfc ',' wrong size ngl ', ierr )

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
            IF( ionode ) WRITE(iuni) idum
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
            IF( ionode ) WRITE(iuni) idum
          END IF
        END DO

        DEALLOCATE( wtmp )

      ELSE

        CALL write_restart_wfc2(iuni, nbnd)

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_wfc2(iuni, nbnd)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni, nbnd
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum, i
      idum = nbnd
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
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
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_wfc1(iuni, tovrw, tread, &
      ik, nk, kunit, ispin, nspin, scal, wf0, t0, wfm, tm, ngw, nbnd, igl, tigl, ngwl )
!
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_put, mp_bcast, mp_max, mp_get
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      COMPLEX(dbl), INTENT(INOUT) :: wf0(:,:)
      COMPLEX(dbl), INTENT(INOUT) :: wfm(:,:)
      INTEGER, INTENT(INOUT) :: ngw, nbnd, ik, nk, kunit, ispin, nspin
      REAL(dbl), INTENT(INOUT) :: scal
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)
      LOGICAL, INTENT(INOUT) :: t0, tm, tigl

! ... If tigl is true then use the global to local mapping stored in the file

      INTEGER :: i, j, idum, ierr
      COMPLEX(dbl), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)
      LOGICAL :: t0_, tm_
      INTEGER :: ngw_, nbnd_, ik_, nk_, ispin_, nspin_, kunit_
      LOGICAL :: twrite_
      INTEGER :: file_version_

      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx, igwx_
      INTEGER :: npool, ipmask( nproc ), ipdest
      REAL(dbl) :: scal_
!
! ... Subroutine Body
!
      IF( ionode ) READ(iuni) twrite_, file_version_
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version, ionode_id )

      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_wfc ',' Data Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) READ(iuni) ngw_, nbnd_, ik_, nk_, kunit_, ispin_, nspin_, scal_
        IF( ionode ) READ(iuni) igwx_
        CALL mp_bcast( ngw_, ionode_id )
        CALL mp_bcast( nbnd_, ionode_id )
        CALL mp_bcast( ik_, ionode_id )
        CALL mp_bcast( nk_, ionode_id )
        CALL mp_bcast( kunit_, ionode_id )
        CALL mp_bcast( ispin_, ionode_id )
        CALL mp_bcast( nspin_, ionode_id )
        CALL mp_bcast( scal_, ionode_id )
        CALL mp_bcast( igwx_, ionode_id )

        IF( tovrw ) THEN

! ...     Use values read from the file
          ngw   = ngw_
          nbnd  = nbnd_
          ik    = ik_
          nk    = nk_
          kunit = kunit_
          ispin = ispin_
          nspin = nspin_ 

        ELSE

          ! IF( ngw_ /= ngw ) THEN
          !   WRITE(6,fmt="(3X,'W: read_restart_wfc, ngw changed' )")
          !   WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") ngw_, ngw
          ! END IF
          IF( nbnd_ /= nbnd ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_wfc, nbnd changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") nbnd_, nbnd
          END IF
          IF( nspin_ /= nspin ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_wfc, nspin changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") nspin_, nspin
          END IF
          IF( nk_ /= nk ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_wfc, nk changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") nk_, nk
          END IF
          IF( kunit_ /= kunit ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_wfc, kunit changed' )")
            WRITE(6,fmt="(3X,'W: old = ',I10,' new = ',I10 )") kunit_, kunit
          END IF

        END IF

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
          CALL error(' read_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipdest /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipdest, 1 )
        END IF

!
! ...   Here read wave function at time t
!  
        IF( ionode ) READ(iuni) t0_
        CALL mp_bcast( t0_, ionode_id )

        IF( .NOT. ( t0_ .AND. t0 ) .AND. ( restart_module_verbosity > 1000 ) ) &
          WRITE(6,fmt="(3X,'W: read_restart_wfc, wf0 not read from restart ' )")

        ! ... WRITE(6,*) ' #### ', igwx_, igwx, ngwl, iks, ikt, ike ! DEBUG

        DO j = 1, nbnd_

          IF( t0_ .AND. t0 ) THEN

            ALLOCATE( wtmp( MAX(igwx_, igwx) ) )

            IF( ionode ) READ(iuni) ( wtmp(i), i=1,igwx_ )
            IF( igwx > igwx_ ) wtmp( (igwx_ + 1) : igwx ) = 0.0d0

            IF( j <=  nbnd ) THEN

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

            END IF

            DEALLOCATE( wtmp )

          ELSE

            IF( ionode ) READ(iuni) idum

          END IF

        END DO

        IF( ( nbnd_ > nbnd ) .AND. ( restart_module_verbosity > 100 ) ) THEN
          WRITE(6,fmt="(3X,'W: read_restart_wfc, not all wf0 read from restart ' )")
        END IF

!
! ...   Here read wave function at time t-dt
!
        IF( ionode ) READ(iuni) tm_
        CALL mp_bcast( tm_, ionode_id )

        IF( .NOT. ( tm_ .AND. tm ) .AND. ( restart_module_verbosity > 1000 ) ) &
          WRITE(6,fmt="(3X,'W: read_restart_wfc, wfm not read from restart ' )")

        DO j = 1, nbnd_

          IF( tm_ .AND. tm ) THEN

            ALLOCATE( wtmp( MAX(igwx_, igwx) ) )

            IF( ionode ) READ(iuni) ( wtmp(i), i=1,igwx_ )
            IF( igwx > igwx_ ) wtmp( (igwx_ + 1) : igwx ) = 0.0d0

            IF( j <=  nbnd ) THEN

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

            END IF

            DEALLOCATE( wtmp )

          ELSE

            IF( ionode ) READ(iuni) idum

          END IF

        END DO

        IF( ( nbnd_ > nbnd ) .AND. ( restart_module_verbosity > 100 ) ) THEN
          WRITE(6,fmt="(3X,'W: read_restart_wfc, not all wfm read from restart ' )")
        END IF

! ...   this is to inform the calling subroutine on what has been read
!
        t0   = t0_
        tm   = tm_
        scal = scal_

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum, nbnd_
          READ(iuni) idum
          READ(iuni) idum    ! t0
          DO i = 1, nbnd_
            READ(iuni) idum
          END DO
          READ(iuni) idum    ! t1
          DO i = 1, nbnd_
            READ(iuni) idum
          END DO
        END IF
        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_wfc, Data Section not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_wfc2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum, i, nbnd_
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
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
        WRITE(6,fmt="(3X,'W: read_restart_wfc, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

! ..  This subroutine write potential and charge density to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_charge1(iuni, twrite, &
      rhog, tr, vg, tv, ng, ispin, nspin, igl, ngl)

      USE mp_wave
      USE mp_global, ONLY: mpime, nproc, root
      USE io_global, ONLY: ionode, ionode_id
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      COMPLEX(dbl), INTENT(IN) :: rhog(:)
      COMPLEX(dbl), INTENT(IN) :: vg(:)
      INTEGER, INTENT(IN) :: ispin, nspin, ng, ngl, igl(:)
      LOGICAL, INTENT(IN) :: tr, tv
      INTEGER :: i, is
      COMPLEX(dbl), ALLOCATABLE :: vtmp(:)
      INTEGER :: idum
!
! ... Subroutine Body
!

      IF ( twrite ) THEN

        ALLOCATE( vtmp (ng) )
        IF( ionode ) WRITE(iuni) twrite, file_version
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

      ELSE
        
        CALL write_restart_charge2( iuni )

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_charge2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      idum    = 0
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
        WRITE(iuni) idum
      END IF
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_charge1(iuni, tovrw, tread, &
      rhog, tr, vg, tv, ng, ispin, nspin, igl, ngl)

      USE mp_wave
      USE mp_global, ONLY: mpime, nproc, root
      USE io_global, ONLY: ionode, ionode_id
      USE mp, ONLY: mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      COMPLEX(dbl), INTENT(OUT) :: rhog(:)
      COMPLEX(dbl), INTENT(OUT) :: vg(:)
      INTEGER, INTENT(INOUT) :: ispin, nspin, ng, ngl, igl(:)
      LOGICAL, INTENT(INOUT) :: tr, tv
      INTEGER :: i, j, k, is, ng_, nspin_, ispin_
      LOGICAL :: tr_, tv_
      COMPLEX(dbl), ALLOCATABLE :: vtmp(:)

      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
!
! ... Subroutine Body
!

      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version, ionode_id )
      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_charge ',' Data Section not present in restart file ', 1)

      IF( tread ) THEN

        IF( ionode ) READ(iuni) ng_, ispin_, nspin_
        CALL mp_bcast( ng_, ionode_id )
        CALL mp_bcast( ispin_, ionode_id )
        CALL mp_bcast( nspin_, ionode_id )

        IF( .NOT. tovrw ) THEN
          IF( nspin /= nspin_ ) &
            CALL error(' read_restart_charge ',' wrong spin ', 1)  
          IF( ng_ /= ng ) THEN
            WRITE(6,fmt="(3X,'W: read_restart_charge, number of G-vecs changed ' )")
          END IF
        END IF

        IF( ionode ) READ(iuni) tr_
        CALL mp_bcast( tr_, ionode_id )

        IF( tr .AND. .NOT. tr_ ) &
          CALL error(' read_restart_charge ',' rho not present in restart ', 1)  

        IF( tr_ ) THEN
          ALLOCATE( vtmp( MAX(ng_, ng) ) )
          IF( ionode ) THEN
            READ(iuni) (vtmp(i),i=1,ng_)
            IF( ng_ < ng ) vtmp( ng_+1 : ng ) = 0.0d0
          END IF
          CALL splitwf(rhog(:), vtmp, ngl, igl, mpime, nproc, ionode_id)
          DEALLOCATE( vtmp )
        ELSE
          IF( ionode ) READ(iuni) idum
        END IF

        IF( ionode ) READ(iuni) tv_
        CALL mp_bcast( tv_, ionode_id )

        IF( tv .AND. .NOT. tv_ ) &
          CALL error(' read_restart_charge ',' V not present in restart ', 1)  

        IF( tv_ ) THEN
          ALLOCATE( vtmp( MAX(ng_, ng) ) )
          IF( ionode ) THEN
            READ(iuni) (vtmp(i),i=1,ng_)
            IF( ng_ < ng ) vtmp( ng_+1 : ng ) = 0.0d0
          END IF
          CALL splitwf(vg(:), vtmp, ngl, igl, mpime, nproc, ionode_id)
          DEALLOCATE( vtmp )
        ELSE
          IF( ionode ) READ(iuni) idum
        END IF
  

        IF( tovrw ) THEN
          ng    = ng_
          ispin = ispin_
          nspin = nspin_
        END IF

      ELSE

        IF( ionode ) THEN
          READ(iuni) idum
          READ(iuni) idum
          READ(iuni) idum
          READ(iuni) idum
          READ(iuni) idum
        END IF

        IF( restart_module_verbosity > 1000 ) &
          WRITE(6,fmt="(3X,'W: read_restart_charge, Data Section not read from restart ' )")

      END IF

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_charge2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum, i, nspin_
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
        READ(iuni) idum
      END IF
      IF( restart_module_verbosity > 1000 ) &
        WRITE(6,fmt="(3X,'W: read_restart_charge, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_template1(iuni, twrite)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: twrite
      INTEGER :: idum
      IF( twrite ) THEN
        IF( ionode ) THEN
          WRITE(iuni) twrite, file_version
        END IF
      ELSE
        CALL write_restart_template2(iuni)
      END IF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_template2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum
      IF( ionode ) THEN
        WRITE(iuni) twrite, file_version
      END IF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_template1(iuni, tovrw, tread)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL, INTENT(IN) :: tovrw
      LOGICAL, INTENT(IN) :: tread
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      CALL mp_bcast( twrite_, ionode_id )
      CALL mp_bcast( file_version, ionode_id )
      IF( tread .AND. .NOT. twrite_ ) &
        CALL error(' read_restart_template ',' Data Section not present in restart file ', 1)
      IF( tread ) THEN
        IF( ionode ) THEN
        END IF
        IF( tovrw ) THEN
        END IF
      ELSE
        IF( ionode ) THEN
        END IF
        WRITE(6,fmt="(3X,'W: read_restart_template, Data Section not read from restart ' )")
      ENDIF
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE read_restart_template2(iuni)
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: group
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni
      LOGICAL :: twrite_
      INTEGER :: file_version_
      INTEGER :: idum
      IF( ionode ) THEN
        READ(iuni) twrite_, file_version_
      END IF
      WRITE(6,fmt="(3X,'W: read_restart_template, Data Section not read from restart ' )")
      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!
  END MODULE
!=----------------------------------------------------------------------------=!
