!
! Copyright (C) 2002-2003 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  extracted from module "pseudo_types" of FPMD
!

      MODULE pseudo_types

!  this module contains the definitions of several TYPE structures,
!  together with their allocation/deallocation routines

        USE kinds, ONLY: dbl
        USE parameters, ONLY: mmaxx, cp_lmax

        IMPLICIT NONE
        SAVE

!  BEGIN manual
!  TYPE DEFINITIONS

        TYPE pseudo_upf
          CHARACTER(LEN=80):: generated   ! 
          CHARACTER(LEN=80):: date_author ! Misc info
          CHARACTER(LEN=80):: comment     !
          CHARACTER(LEN=2) :: psd       ! Element label
          CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
          LOGICAL  :: tvanp             ! .true. if Ultrasoft
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
          CHARACTER(LEN=2), POINTER :: els(:)  ! els(nwfc)
          INTEGER, POINTER :: lchi(:)   ! lchi(nwfc)
          REAL(dbl), POINTER :: oc(:)   ! oc(nwfc)
          REAL(dbl), POINTER :: r(:)    ! r(mesh)
          REAL(dbl), POINTER :: rab(:)  ! rab(mesh)
          REAL(dbl), POINTER :: rho_atc(:) ! rho_atc(mesh)
          REAL(dbl), POINTER :: vloc(:)    ! vloc(mesh)
          INTEGER, POINTER :: lll(:)       ! lll(nbeta)
          INTEGER, POINTER :: kkbeta(:)    ! kkbeta(nbeta)
          REAL(dbl), POINTER :: beta(:,:)  ! beta(mesh,nbeta)
          INTEGER :: nd
          REAL(dbl), POINTER :: dion(:,:)  ! dion(nbeta,nbeta)
          INTEGER :: nqf
          INTEGER :: nqlc
          REAL(dbl), POINTER :: rinner(:)  ! rinner(0:2*lmax)
          REAL(dbl), POINTER :: qqq(:,:)   ! qqq(nbeta,nbeta)
          REAL(dbl), POINTER :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
          REAL(dbl), POINTER :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
          REAL(dbl), POINTER :: chi(:,:) !  chi(mesh,nwfc)
          REAL(dbl), POINTER :: rho_at(:) !  rho_at(mesh)
        END TYPE


        TYPE pseudo_ncpp
          CHARACTER(LEN=4) :: psd         ! Element label
          CHARACTER(LEN=20) :: pottyp     ! Potential type
          LOGICAL :: tmix
          LOGICAL :: tnlcc
          INTEGER :: igau
          INTEGER :: lloc
          INTEGER :: lnl 
          INTEGER :: indl(mmaxx)
          INTEGER :: nchan
          INTEGER :: mesh
          REAL(dbl) ::  zv
          REAL(dbl) ::  raggio
          REAL(dbl) ::  dx            ! r(i) = cost * EXP( xmin + dx * (i-1) )
          REAL(dbl) ::  rab(mmaxx)
          REAL(dbl) ::  rw(mmaxx)
          REAL(dbl) ::  vnl(mmaxx, cp_lmax)
          REAL(dbl) ::  vloc(mmaxx)
          REAL(dbl) ::  vrps(mmaxx, cp_lmax)
          REAL(dbl) ::  wgv(cp_lmax)
          REAL(dbl) ::  rc(2)
          REAL(dbl) ::  wrc(2)
          REAL(dbl) ::  rcl(3,3)
          REAL(dbl) ::  al(3,3)
          REAL(dbl) ::  bl(3,3)
          INTEGER :: nrps                     ! number of atomic wave function
          INTEGER :: lrps(cp_lmax)            ! angular momentum
          REAL(dbl) :: oc(cp_lmax)            ! occupation for each rps
          REAL(dbl) :: rps(mmaxx, cp_lmax)  ! atomic pseudo wave function
          REAL(dbl) :: rhoc(mmaxx)          ! core charge
        END TYPE pseudo_ncpp

!  ----------------------------------------------
!  END manual

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines

!  ----------------------------------------------

        SUBROUTINE nullify_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          NULLIFY( upf%els, upf%lchi, upf%oc )
          NULLIFY( upf%r, upf%rab )  
          NULLIFY( upf%rho_atc, upf%vloc )  
          NULLIFY( upf%lll, upf%kkbeta, upf%beta, upf%dion )  
          NULLIFY( upf%rinner, upf%qqq, upf%qfunc, upf%qfcoef )  
          NULLIFY( upf%chi )  
          NULLIFY( upf%rho_at )  
          RETURN
        END SUBROUTINE

        SUBROUTINE deallocate_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          IF( ASSOCIATED( upf%els ) ) DEALLOCATE( upf%els )
          IF( ASSOCIATED( upf%lchi ) ) DEALLOCATE( upf%lchi )
          IF( ASSOCIATED( upf%oc ) ) DEALLOCATE( upf%oc )
          IF( ASSOCIATED( upf%r ) ) DEALLOCATE( upf%rab )
          IF( ASSOCIATED( upf%rho_atc ) ) DEALLOCATE( upf%rho_atc )
          IF( ASSOCIATED( upf%vloc ) ) DEALLOCATE( upf%vloc )
          IF( ASSOCIATED( upf%lll ) ) DEALLOCATE( upf%lll )
          IF( ASSOCIATED( upf%kkbeta ) ) DEALLOCATE( upf%kkbeta )
          IF( ASSOCIATED( upf%beta ) ) DEALLOCATE( upf%beta )
          IF( ASSOCIATED( upf%dion ) ) DEALLOCATE( upf%dion )
          IF( ASSOCIATED( upf%rinner ) ) DEALLOCATE( upf%rinner )
          IF( ASSOCIATED( upf%qqq ) ) DEALLOCATE( upf%qqq )
          IF( ASSOCIATED( upf%qfunc ) ) DEALLOCATE( upf%qfunc )
          IF( ASSOCIATED( upf%qfcoef ) ) DEALLOCATE( upf%qfcoef )
          IF( ASSOCIATED( upf%chi ) ) DEALLOCATE( upf%chi )
          IF( ASSOCIATED( upf%rho_at ) ) DEALLOCATE( upf%rho_at )
          RETURN
        END SUBROUTINE

      END MODULE pseudo_types

