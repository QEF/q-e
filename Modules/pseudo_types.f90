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

        USE kinds, ONLY: DP
        USE parameters, ONLY: ndmx, lmaxx

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
          REAL(DP) :: zp               ! z valence
          REAL(DP) :: etotps           ! total energy
          REAL(DP) :: ecutwfc          ! suggested cut-off for wfc
          REAL(DP) :: ecutrho          ! suggested cut-off for rho

          LOGICAL :: has_so             ! if .true. includes spin-orbit
          REAL(DP) :: xmin             ! the minimum x of the linear mesh
          REAL(DP) :: rmax             ! the maximum radius of the mesh
          REAL(DP) :: zmesh            ! the nuclear charge used for mesh
          REAL(DP) :: dx               ! the deltax of the linear mesh
          INTEGER, POINTER :: nn(:)      ! nn(nwfc)
          REAL(DP), POINTER :: rcut(:)  ! cut-off radius(nbeta)
          REAL(DP), POINTER :: rcutus(:)! cut-off ultrasoft radius (nbeta)
          REAL(DP), POINTER :: epseu(:) ! energy (nwfc)
          REAL(DP), POINTER :: jchi(:)  ! jchi(nwfc)
          REAL(DP), POINTER :: jjj(:)   ! jjj(nbeta)

          INTEGER :: nv                 ! UPF file version number
          INTEGER :: lmax               ! maximum angular momentum component
          INTEGER :: mesh               ! number of point in the radial mesh
          INTEGER :: nwfc               ! number of wavefunctions
          INTEGER :: nbeta              ! number of projectors
          CHARACTER(LEN=2), POINTER :: els(:)  ! els(nwfc)
          CHARACTER(LEN=2), POINTER :: els_beta(:)  ! els(nbeta)
          INTEGER, POINTER :: lchi(:)   ! lchi(nwfc)
          REAL(DP), POINTER :: oc(:)   ! oc(nwfc)
          REAL(DP), POINTER :: r(:)    ! r(mesh)
          REAL(DP), POINTER :: rab(:)  ! rab(mesh)
          REAL(DP), POINTER :: rho_atc(:) ! rho_atc(mesh)
          REAL(DP), POINTER :: vloc(:)    ! vloc(mesh)
          INTEGER, POINTER :: lll(:)       ! lll(nbeta)
          INTEGER, POINTER :: kkbeta(:)    ! kkbeta(nbeta)
          REAL(DP), POINTER :: beta(:,:)  ! beta(mesh,nbeta)
          INTEGER :: nd
          REAL(DP), POINTER :: dion(:,:)  ! dion(nbeta,nbeta)
          INTEGER :: nqf
          INTEGER :: nqlc
          REAL(DP), POINTER :: rinner(:)  ! rinner(0:2*lmax)
          REAL(DP), POINTER :: qqq(:,:)   ! qqq(nbeta,nbeta)
          REAL(DP), POINTER :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
          REAL(DP), POINTER :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
          REAL(DP), POINTER :: chi(:,:) !  chi(mesh,nwfc)
          REAL(DP), POINTER :: rho_at(:) !  rho_at(mesh)
        END TYPE


        TYPE pseudo_ncpp
          CHARACTER(LEN=4) :: psd         ! Element label
          CHARACTER(LEN=20) :: pottyp     ! Potential type
          LOGICAL :: tmix
          LOGICAL :: tnlcc
          INTEGER :: igau
          INTEGER :: lloc
          INTEGER :: nbeta 
          INTEGER :: lll(lmaxx+1)
          INTEGER :: nchan
          INTEGER :: mesh
          REAL(DP) ::  zv
          REAL(DP) ::  dx            ! r(i) = cost * EXP( xmin + dx * (i-1) )
          REAL(DP) ::  rab(ndmx)
          REAL(DP) ::  rw(ndmx)
          REAL(DP) ::  vnl(ndmx, lmaxx+1)
          REAL(DP) ::  vloc(ndmx)
          REAL(DP) ::  vrps(ndmx, lmaxx+1)
          REAL(DP) ::  wgv(lmaxx+1)
          REAL(DP) ::  rc(2)
          REAL(DP) ::  wrc(2)
          REAL(DP) ::  rcl(3,3)
          REAL(DP) ::  al(3,3)
          REAL(DP) ::  bl(3,3)
          INTEGER :: nrps                     ! number of atomic wave function
          INTEGER :: lrps(lmaxx+1)            ! angular momentum
          REAL(DP) :: oc(lmaxx+1)            ! occupation for each rps
          REAL(DP) :: rps(ndmx, lmaxx+1)  ! atomic pseudo wave function
          REAL(DP) :: rhoc(ndmx)          ! core charge
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
          NULLIFY( upf%els, upf%lchi, upf%jchi, upf%oc )
          NULLIFY( upf%r, upf%rab )  
          NULLIFY( upf%rho_atc, upf%vloc )  
          NULLIFY( upf%nn, upf%rcut)
          NULLIFY( upf%els_beta)
          NULLIFY( upf%rcutus, upf%epseu)
          NULLIFY( upf%lll, upf%jjj, upf%kkbeta, upf%beta, upf%dion )  
          NULLIFY( upf%rinner, upf%qqq, upf%qfunc, upf%qfcoef )  
          NULLIFY( upf%chi )  
          NULLIFY( upf%rho_at )  
          RETURN
        END SUBROUTINE nullify_pseudo_upf

        SUBROUTINE deallocate_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          IF( ASSOCIATED( upf%els ) ) DEALLOCATE( upf%els )
          IF( ASSOCIATED( upf%lchi ) ) DEALLOCATE( upf%lchi )
          IF( ASSOCIATED( upf%jchi ) ) DEALLOCATE( upf%jchi )
          IF( ASSOCIATED( upf%oc ) ) DEALLOCATE( upf%oc )
          IF( ASSOCIATED( upf%r ) ) DEALLOCATE( upf%r )
          IF( ASSOCIATED( upf%rab ) ) DEALLOCATE( upf%rab )
          IF( ASSOCIATED( upf%nn ) ) DEALLOCATE( upf%nn )
          IF( ASSOCIATED( upf%els_beta ) ) DEALLOCATE( upf%els_beta )
          IF( ASSOCIATED( upf%rcut ) ) DEALLOCATE( upf%rcut )
          IF( ASSOCIATED( upf%rcutus ) ) DEALLOCATE( upf%rcutus )
          IF( ASSOCIATED( upf%epseu ) ) DEALLOCATE( upf%epseu )
          IF( ASSOCIATED( upf%rho_atc ) ) DEALLOCATE( upf%rho_atc )
          IF( ASSOCIATED( upf%vloc ) ) DEALLOCATE( upf%vloc )
          IF( ASSOCIATED( upf%lll ) ) DEALLOCATE( upf%lll )
          IF( ASSOCIATED( upf%jjj ) ) DEALLOCATE( upf%jjj )
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
        END SUBROUTINE deallocate_pseudo_upf

      END MODULE pseudo_types

