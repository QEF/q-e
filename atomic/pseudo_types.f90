!
! Copyright (C) 2002-2003 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  extracted from module "pseudo_types" of FPMD
!

      MODULE pseudo_types_ld1

!  this module contains the definitions of several TYPE structures,
!  together with their allocation/deallocation routines

        USE kinds, only: dbl

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
          LOGICAL ::  nlcc              ! Non linear core corrections
          CHARACTER(LEN=20) :: dft      ! Exch-Corr type
          REAL(dbl) :: zp               ! z valence
          REAL(dbl) :: etotps           ! total energy
          REAL(dbl) :: ecutwfc          ! suggested cut-off for wfc
          REAL(dbl) :: ecutrho          ! suggested cut-off for rho
          INTEGER :: nv                 ! UPF file version number
          INTEGER :: mesh               ! number of point in the radial mesh
          REAL(dbl), POINTER :: r(:)    ! r(mesh)
          REAL(dbl), POINTER :: rab(:)  ! rab(mesh)
          REAL(dbl) :: xmin             ! the minimum x of the linear mesh
          REAL(dbl) :: rmax             ! the maximum radius of the mesh
          REAL(dbl) :: zmesh            ! the nuclear charge used for mesh
          REAL(dbl) :: dx               ! the deltax of the linear mesh
          INTEGER :: lmax               ! maximum angular momentum component
          INTEGER :: nwfc               ! number of wavefunctions
          CHARACTER(LEN=2), POINTER :: els(:)  ! els(nwfc)
          INTEGER, POINTER :: nn(:)      ! n of the wavefunction     (nwfc)
          INTEGER, POINTER :: lchi(:)    ! l of the wavefunction     (nwfc)
          REAL(dbl), POINTER :: jchi(:)  ! j total angular momentum  (nwfc)
          REAL(dbl), POINTER :: oc(:)    ! the occupation            (nwfc)
          REAL(dbl), POINTER :: epseu(:) ! energy of wavefunction    (nwfc)
          REAL(dbl), POINTER :: rcut(:)  ! cut-off radius            (nwfc)
          REAL(dbl), POINTER :: rcutus(:)! cut-off ultrasoft radius  (nwfc)
          REAL(dbl), POINTER :: chi(:,:) ! the wavefunctions chi(mesh,nwfc)            
          INTEGER :: nbeta              ! number of projectors
          INTEGER, POINTER :: kkbeta(:) ! kkbeta(nbeta)
          INTEGER, POINTER :: lll(:)    ! lll(nbeta)
          REAL(dbl), POINTER :: jjj(:)  ! total angular momentum j (nbeta)
          REAL(dbl), POINTER :: beta(:,:)! beta(mesh,nbeta)
          REAL(dbl), POINTER :: dion(:,:)! dion(nbeta,nbeta)
          INTEGER :: nd
          INTEGER :: nqf
          INTEGER :: nqlc
          REAL(dbl), POINTER :: rinner(:)  ! rinner(0:2*lmax)
          REAL(dbl), POINTER :: qqq(:,:)   ! qqq(nbeta,nbeta)
          REAL(dbl), POINTER :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
          REAL(dbl), POINTER :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
          REAL(dbl), POINTER :: vloc(:)    ! vloc(mesh)
          REAL(dbl), POINTER :: rho_at(:)  ! rho_at(mesh)
          REAL(dbl), POINTER :: rho_atc(:) ! rho_atc(mesh)
        END TYPE

!  ----------------------------------------------
!  END manual

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines

!  ----------------------------------------------

        SUBROUTINE nullify_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          NULLIFY( upf%r, upf%rab )  
          NULLIFY( upf%els, upf%nn, upf%lchi, upf%jchi, upf%oc )
          NULLIFY( upf%epseu, upf%rcut, upf%rcutus )
          NULLIFY( upf%chi )
          NULLIFY( upf%kkbeta, upf%lll, upf%jjj, upf%beta, upf%dion )
          NULLIFY( upf%rinner, upf%qqq, upf%qfunc, upf%qfcoef )
          NULLIFY( upf%vloc )
          NULLIFY( upf%rho_at )
          NULLIFY( upf%rho_atc )
          RETURN
        END SUBROUTINE

        SUBROUTINE deallocate_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          IF( ASSOCIATED( upf%r ) )      DEALLOCATE( upf%r )
          IF( ASSOCIATED( upf%rab ) )    DEALLOCATE( upf%rab )
          IF( ASSOCIATED( upf%els ) )    DEALLOCATE( upf%els )
          IF( ASSOCIATED( upf%nn ) )     DEALLOCATE( upf%nn )
          IF( ASSOCIATED( upf%lchi ) )   DEALLOCATE( upf%lchi )
          IF( ASSOCIATED( upf%jchi ) )   DEALLOCATE( upf%jchi )
          IF( ASSOCIATED( upf%oc ) )     DEALLOCATE( upf%oc )
          IF( ASSOCIATED( upf%epseu ) )  DEALLOCATE( upf%epseu )
          IF( ASSOCIATED( upf%rcut ) )   DEALLOCATE( upf%rcut )
          IF( ASSOCIATED( upf%rcutus ) ) DEALLOCATE( upf%rcutus )
          IF( ASSOCIATED( upf%chi ) )    DEALLOCATE( upf%chi )
          IF( ASSOCIATED( upf%kkbeta ) ) DEALLOCATE( upf%kkbeta )
          IF( ASSOCIATED( upf%lll ) )    DEALLOCATE( upf%lll )
          IF( ASSOCIATED( upf%jjj ) )    DEALLOCATE( upf%jjj )
          IF( ASSOCIATED( upf%beta ) )   DEALLOCATE( upf%beta )
          IF( ASSOCIATED( upf%dion ) )   DEALLOCATE( upf%dion )
          IF( ASSOCIATED( upf%rinner ) ) DEALLOCATE( upf%rinner )
          IF( ASSOCIATED( upf%qqq ) )    DEALLOCATE( upf%qqq )
          IF( ASSOCIATED( upf%qfunc ) )  DEALLOCATE( upf%qfunc )
          IF( ASSOCIATED( upf%qfcoef ) ) DEALLOCATE( upf%qfcoef )
          IF( ASSOCIATED( upf%vloc ) )   DEALLOCATE( upf%vloc )
          IF( ASSOCIATED( upf%rho_at ) ) DEALLOCATE( upf%rho_at )
          IF( ASSOCIATED( upf%rho_atc ) ) DEALLOCATE( upf%rho_atc )
          RETURN
        END SUBROUTINE

      END MODULE pseudo_types_ld1
