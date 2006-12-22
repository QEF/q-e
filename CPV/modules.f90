!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module bhs
  !     analytical BHS pseudopotential parameters
  use parameters, only: nsx
  implicit none
  save
  real(8) :: rc1(nsx), rc2(nsx), wrc1(nsx), wrc2(nsx), &
       rcl(3,nsx,3), al(3,nsx,3), bl(3,nsx,3)
  integer :: lloc(nsx)
end module bhs

!     f    = occupation numbers
!     qbac = background neutralizing charge
!     nspin = number of spins (1=no spin, 2=LSDA)
!     nel(nspin) = number of electrons (up, down)
!     nupdwn= number of states with spin up (1) and down (2)
!     iupdwn=      first state with spin (1) and down (2)
!     n     = total number of electronic states
!     nx    = if n is even, nx=n ; if it is odd, nx=n+1
!             nx is used only to dimension arrays
!     ispin = spin of each state
!

!     tpiba   = 2*pi/alat
!     tpiba2  = (2*pi/alat)**2
!     ng      = number of G vectors for density and potential
!     ngl     = number of shells of G

!     G-vector quantities for the thick grid - see also doc in ggen 
!     g       = G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
!     gl      = shells of G^2           ( "   "   "    "      "      )
!     gx      = G-vectors               ( "   "   "  tpiba =(2pi/a)  )
!
!     g2_g    = all G^2 in increasing order, replicated on all procs
!     mill_g  = miller index of G vecs (increasing order), replicated on all procs
!     mill_l  = miller index of G vecs local to the processors
!     ig_l2g  = "l2g" means local to global, this array convert a local
!               G-vector index into the global index, in other words
!               the index of the G-v. in the overall array of G-vectors
!     bi?     = base vector used to generate the reciprocal space
!
!     np      = fft index for G>
!     nm      = fft index for G<
!     mill_l  = G components in crystal axis
!


!
!  lqmax:  maximum angular momentum of Q (Vanderbilt augmentation charges)
!  nqfx :  maximum number of coefficients in Q smoothing
!  nbrx :  maximum number of distinct radial beta functions
!  ndmx:  maximum number of points in the radial grid
! 

!  nbeta    number of beta functions (sum over all l)
!  kkbeta   last radial mesh point used to describe functions
!                 which vanish outside core
!  nqf      coefficients in Q smoothing
!  nqlc     angular momenta present in Q smoothing
!  lll      lll(j) is l quantum number of j'th beta function
!  lmaxq      highest angular momentum that is present in Q functions
!  lmaxkb   highest angular momentum that is present in beta functions
!  dion     bare pseudopotential D_{\mu,\nu} parameters
!              (ionic and screening parts subtracted out)
!  betar    the beta function on a r grid (actually, r*beta)
!  qqq      Q_ij matrix
!  qfunc    Q_ij(r) function (for r>rinner)
!  rinner   radius at which to cut off partial core or Q_ij
!
!  qfcoef   coefficients to pseudize qfunc for different total
!              angular momentum (for r<rinner)
!  vloc_at  local potential for each atom


module local_pseudo
  implicit none
  save
  !
  !    rhops = ionic pseudocharges (for Ewald term)
  !    vps   = local pseudopotential in G space for each species
  !
  real(8), allocatable:: rhops(:,:), vps(:,:)
  !
  !    drhops = derivative of rhops respect to G^2
  !    dvps   = derivative of vps respect to G^2
  !
  real(8),allocatable:: dvps(:,:), drhops(:,:)
  !
contains
  !
  subroutine allocate_local_pseudo( ng, nsp )
      integer, intent(in) :: ng, nsp
      call deallocate_local_pseudo()
      ALLOCATE( rhops( ng, nsp ) )
      ALLOCATE( vps( ng, nsp ) )
      ALLOCATE( drhops( ng, nsp ) )
      ALLOCATE( dvps( ng, nsp ) )
  end subroutine
  !
  subroutine deallocate_local_pseudo
      IF( ALLOCATED( rhops ) ) DEALLOCATE( rhops )
      IF( ALLOCATED( vps ) ) DEALLOCATE( vps )
      IF( ALLOCATED( dvps ) ) DEALLOCATE( dvps )
      IF( ALLOCATED( drhops ) ) DEALLOCATE( drhops )
  end subroutine
  !
end module local_pseudo

module qgb_mod
  implicit none
  save
  complex(8), allocatable :: qgb(:,:,:)
contains
  subroutine deallocate_qgb_mod
      IF( ALLOCATED( qgb ) ) DEALLOCATE( qgb )
  end subroutine deallocate_qgb_mod
end module qgb_mod

module qradb_mod
  implicit none
  save
  real(8), allocatable:: qradb(:,:,:,:)
contains
  subroutine deallocate_qradb_mod
      IF( ALLOCATED( qradb ) ) DEALLOCATE( qradb )
  end subroutine deallocate_qradb_mod
end module qradb_mod

! Variable cell
module derho
  use kinds, only: DP
  implicit none
  save
  complex(DP),allocatable:: drhog(:,:,:,:)
  real(DP),allocatable::    drhor(:,:,:,:)
contains
  subroutine deallocate_derho
      IF( ALLOCATED( drhog ) ) DEALLOCATE( drhog )
      IF( ALLOCATED( drhor ) ) DEALLOCATE( drhor )
  end subroutine deallocate_derho
end module derho

MODULE metagga  !metagga
  !the variables needed for meta-GGA
  REAL(8), ALLOCATABLE :: &
       kedtaus(:,:), &! KineticEnergyDensity in real space,smooth grid
       kedtaur(:,:), &! real space, density grid
       crosstaus(:,:,:), &!used by stress tensor,in smooth grid
       dkedtaus(:,:,:,:)  !derivative of kedtau wrt h on smooth grid
  COMPLEX(8) , ALLOCATABLE :: &
       kedtaug(:,:),    & !KineticEnergyDensity in G space
       gradwfc(:,:)    !used by stress tensor
contains
  subroutine deallocate_metagga
      IF( ALLOCATED(crosstaus))DEALLOCATE(crosstaus)
      IF( ALLOCATED(dkedtaus)) DEALLOCATE(dkedtaus)
      IF( ALLOCATED(gradwfc))  DEALLOCATE(gradwfc)
  end subroutine deallocate_metagga
END MODULE metagga  !end metagga

MODULE dener
  USE kinds, ONLY: DP
  IMPLICIT NONE
  SAVE
  REAL(DP) :: dekin(3,3)
  REAL(DP) :: dh(3,3)
  REAL(DP) :: dps(3,3)
  REAL(DP) :: denl(3,3)
  REAL(DP) :: dxc(3,3)
  REAL(DP) :: dsr(3,3)
  REAL(DP) :: detot(3,3)
  REAL(DP) :: dekin6(6)
  REAL(DP) :: dh6(6)
  REAL(DP) :: dps6(6)
  REAL(DP) :: denl6(6)
  REAL(DP) :: dxc6(6)
  REAL(DP) :: dsr6(6)
  REAL(DP) :: detot6(6)
END MODULE dener

module dqgb_mod
  implicit none
  save
  complex(8),allocatable:: dqgb(:,:,:,:,:)
contains
  subroutine deallocate_dqgb_mod
      IF( ALLOCATED( dqgb ) ) DEALLOCATE( dqgb )
  end subroutine deallocate_dqgb_mod
end module dqgb_mod

MODULE cdvan
  USE kinds, ONLY: DP
  IMPLICIT NONE
  SAVE
  REAL(DP), ALLOCATABLE :: dbeta(:,:,:,:,:), dbec(:,:,:,:), &
                             drhovan(:,:,:,:,:)
CONTAINS
  SUBROUTINE deallocate_cdvan
      IF( ALLOCATED( dbeta ) ) DEALLOCATE( dbeta )
      IF( ALLOCATED( dbec ) ) DEALLOCATE( dbec )
      IF( ALLOCATED( drhovan ) ) DEALLOCATE( drhovan )
  END SUBROUTINE deallocate_cdvan
END MODULE cdvan


MODULE ncpp
  !
  ! norm-conserving pseudo-potentials, Kleinman-Bylander factors 
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  SAVE
  REAL(DP), ALLOCATABLE :: wsg(:,:)     ! inverse of Kleinman-Bylander
                                         !   denominators
                                         ! <Y phi | V | phi Y>**(-1)
                                         !   first index: orbital
                                         !   second index: atomic species
  REAL(DP), ALLOCATABLE :: wnl(:,:,:,:) ! Kleinman-Bylander products
                                         ! <Y phi | V | exp(i(k+G) dot r)>
                                         !   first index: G vector
                                         !   second index: orbital
                                         !   third index: atomic species
                                         !   fourth index: k point
CONTAINS

  SUBROUTINE allocate_ncpp( nsp, ngw, nbetax, nhm, nk )
    INTEGER, INTENT(IN) :: nsp, nbetax, nhm, ngw, nk
    INTEGER :: ierr

    ALLOCATE( wnl( ngw, nbetax, nsp, nk ), STAT=ierr)
    IF( ierr /= 0 ) CALL errore(' allocate_ncpp ', ' allocating wnl ', ierr )
    ALLOCATE( wsg( nhm, nsp ), STAT=ierr)
    IF( ierr /= 0 ) CALL errore(' allocate_ncpp ', ' allocating wsg ', ierr )
    RETURN
  END SUBROUTINE allocate_ncpp

  SUBROUTINE deallocate_ncpp
    IF( ALLOCATED( wsg ) ) DEALLOCATE( wsg )
    IF( ALLOCATED( wnl ) ) DEALLOCATE( wnl )
    RETURN
  END SUBROUTINE deallocate_ncpp

END MODULE ncpp

module cvan

  ! this file contains common subroutines and modules between
  ! CP and FPMD

  !     ionic pseudo-potential variables
  use parameters, only: nsx
  implicit none
  save
  integer nvb, ish(nsx)
  !     nvb    = number of species with Vanderbilt PPs
  !     ish(is)= used for indexing the nonlocal projectors betae
  !              with contiguous indices inl=ish(is)+(iv-1)*na(is)+1
  !              where "is" is the species and iv=1,nh(is)
  !
  !     indlm: indlm(ind,is)=Y_lm for projector ind
  integer, allocatable:: indlm(:,:)
contains

  subroutine allocate_cvan( nind, ns )
    integer, intent(in) :: nind, ns
    allocate( indlm( nind, ns ) )
  end subroutine allocate_cvan

  subroutine deallocate_cvan( )
    if( allocated(indlm) ) deallocate( indlm )
  end subroutine deallocate_cvan

end module cvan


MODULE stress_param

   USE kinds, ONLY : DP

   IMPLICIT NONE
   SAVE

   INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
   INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)

   REAL(DP),  DIMENSION(3,3), PARAMETER :: delta = reshape &
         ( (/ 1.0_DP, 0.0_DP, 0.0_DP, &
              0.0_DP, 1.0_DP, 0.0_DP, &
              0.0_DP, 0.0_DP, 1.0_DP  &
            /), (/ 3, 3 /) )

   ! ...  dalbe(:) = delta(alpha(:),beta(:))
   !
   REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)

END MODULE



MODULE core
   !
   USE kinds
   ! 
   IMPLICIT NONE
   SAVE
   !     nlcc_any = 0 no core correction on any atom
   !     rhocb  = core charge in G space (box grid)
   !     rhoc   = core charge in real space  (dense grid)
   !     rhocg  = core charge in G space  (dense grid)
   !     drhocg = derivative of core charge in G space (used for stress)
   !
   LOGICAL :: nlcc_any
   REAL(DP), ALLOCATABLE:: rhocb(:,:)
   REAL(DP), ALLOCATABLE:: rhoc(:)
   REAL(DP), ALLOCATABLE:: rhocg(:,:)
   REAL(DP), ALLOCATABLE:: drhocg(:,:)
   !
CONTAINS
   !
   SUBROUTINE allocate_core( nnrx, ngm, ngb, nsp ) 
     INTEGER, INTENT(IN) :: nnrx, ngm, ngb, nsp
     IF ( nlcc_any ) THEN    
        !
        ALLOCATE( rhoc( nnrx ) )
        ALLOCATE( rhocb( ngb, nsp ) )
        ALLOCATE( rhocg( ngm, nsp ) )
        ALLOCATE( drhocg( ngm, nsp ) )
        !
     ELSE
        !
        ! ... dummy allocation required because this array appears in the
        ! ... list of arguments of some routines
        !
        ALLOCATE( rhoc( 1 ) )
        !
     END IF
   END SUBROUTINE allocate_core
   !
   SUBROUTINE deallocate_core()
      IF( ALLOCATED( rhocb  ) ) DEALLOCATE( rhocb )
      IF( ALLOCATED( rhoc   ) ) DEALLOCATE( rhoc  )
      IF( ALLOCATED( rhocg  ) ) DEALLOCATE( rhocg  )
      IF( ALLOCATED( drhocg ) ) DEALLOCATE( drhocg )
   END SUBROUTINE deallocate_core
   !
END MODULE core
!@@@@
module elct
  USE kinds
  use electrons_base, only: nspin, nel, nupdwn, iupdwn
  use electrons_base, only: n => nbnd, nx => nbndx, ispin, f
  implicit none
  save
  !     f    = occupation numbers
  !     qbac = background neutralizing charge
  real(DP) qbac
  !     nspin = number of spins (1=no spin, 2=LSDA)
  !     nel(nspin) = number of electrons (up, down)
  !     nupdwn= number of states with spin up (1) and down (2)
  !     iupdwn=      first state with spin (1) and down (2)
  !     n     = total number of electronic states
  !     nx    = if n is even, nx=n ; if it is odd, nx=n+1
  !            nx is used only to dimension arrays
  !     ispin = spin of each state
!  integer, allocatable:: ispin(:)
  !
contains
                                                                                                 
  subroutine deallocate_elct()
      IF( ALLOCATED( f ) ) DEALLOCATE( f )
      IF( ALLOCATED( ispin ) ) DEALLOCATE( ispin )
      return
  end subroutine
  !
end module elct
!
module ldaU
  use parameters, only: nsx
  USE kinds
  implicit none
  complex(DP), allocatable :: atomwfc(:,:)
  complex(DP), allocatable :: swfcatom(:,:)
  real(DP) :: Hubbard_U(nsx), Hubbard_lambda(nsx,2), ns0(nsx,2),   &
     & Hubbard_alpha(nsx)
  real(DP) :: e_hubbard = 0.d0, e_lambda = 0.d0
  real(DP), allocatable :: ns(:,:,:,:)
  integer :: Hubbard_l(nsx), Hubbard_lmax=0, n_atomic_wfc
  logical lda_plus_u
  COMPLEX(DP), allocatable::  vupsi(:,:) !@@@@
contains
  !
  subroutine deallocate_lda_plus_u()
     !
     IF( ALLOCATED( atomwfc ) ) DEALLOCATE( atomwfc )
     IF( ALLOCATED( swfcatom ) ) DEALLOCATE( swfcatom )
     IF( ALLOCATED( ns ) ) DEALLOCATE( ns )
     IF( ALLOCATED( vupsi ) ) DEALLOCATE( vupsi )
     !
  end subroutine
  !
end module ldaU
!
! Occupation constraint ...to be implemented...
!
module step_constraint
  use parameters, only: natx_ => natx
  USE kinds
  implicit none
  real(DP) :: E_con
  real(DP) :: A_con(natx_,2), sigma_con(natx_), alpha_con(natx_)
  logical :: step_con
  ! complex(DP), allocatable:: vpsi_con(:,:)
  complex(DP) :: vpsi_con(1,1)
end module step_constraint
!
!@@@@

