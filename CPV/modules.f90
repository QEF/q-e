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
  real(kind=8) :: rc1(nsx), rc2(nsx), wrc1(nsx), wrc2(nsx), &
       rcl(3,nsx,3), al(3,nsx,3), bl(3,nsx,3)
  integer :: lloc(nsx)
end module bhs

module core
  implicit none
  save
  !     nlcc_any = 0 no core correction on any atom
  !     rhocb  = core charge in G space (box grid)
  !     rhoc   = core charge in real space  (dense grid)
  !     rhocg  = core charge in G space  (dense grid)
  !     drhocg = derivative of core charge in G space (used for stress) 
  logical :: nlcc_any
  real(kind=8), allocatable:: rhocb(:,:)
  real(kind=8), allocatable:: rhoc(:)
  real(kind=8), allocatable:: rhocg(:,:)
  real(kind=8), allocatable:: drhocg(:,:)
contains
  subroutine deallocate_core()
      IF( ALLOCATED( rhocb  ) ) DEALLOCATE( rhocb )
      IF( ALLOCATED( rhoc   ) ) DEALLOCATE( rhoc  )
      IF( ALLOCATED( rhocg  ) ) DEALLOCATE( rhocg  )
      IF( ALLOCATED( drhocg ) ) DEALLOCATE( drhocg )
  end subroutine deallocate_core
end module core

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
  real(kind=8), allocatable:: rhops(:,:), vps(:,:)
  !
  !    drhops = derivative of rhops respect to G^2
  !    dvps   = derivative of vps respect to G^2
  !
  real(kind=8),allocatable:: dvps(:,:), drhops(:,:)
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
  complex(kind=8), allocatable :: qgb(:,:,:)
contains
  subroutine deallocate_qgb_mod
      IF( ALLOCATED( qgb ) ) DEALLOCATE( qgb )
  end subroutine deallocate_qgb_mod
end module qgb_mod

module qradb_mod
  implicit none
  save
  real(kind=8), allocatable:: qradb(:,:,:,:,:)
contains
  subroutine deallocate_qradb_mod
      IF( ALLOCATED( qradb ) ) DEALLOCATE( qradb )
  end subroutine deallocate_qradb_mod
end module qradb_mod

! Variable cell
module derho
  implicit none
  save
  complex(kind=8),allocatable:: drhog(:,:,:,:)
  real(kind=8),allocatable::     drhor(:,:,:,:)
contains
  subroutine deallocate_derho
      IF( ALLOCATED( drhog ) ) DEALLOCATE( drhog )
      IF( ALLOCATED( drhor ) ) DEALLOCATE( drhor )
  end subroutine deallocate_derho
end module derho

module dener
  implicit none
  save
  real(kind=8) detot(3,3), dekin(3,3), dh(3,3), dps(3,3), &
  &       denl(3,3), dxc(3,3), dsr(3,3)
end module dener

module dqgb_mod
  implicit none
  save
  complex(kind=8),allocatable:: dqgb(:,:,:,:,:)
contains
  subroutine deallocate_dqgb_mod
      IF( ALLOCATED( dqgb ) ) DEALLOCATE( dqgb )
  end subroutine deallocate_dqgb_mod
end module dqgb_mod

module cdvan
  implicit none
  save
  real(kind=8),allocatable:: dbeta(:,:,:,:,:), dbec(:,:,:,:), &
                             drhovan(:,:,:,:,:)
contains
  subroutine deallocate_cdvan
      IF( ALLOCATED( dbeta ) ) DEALLOCATE( dbeta )
      IF( ALLOCATED( dbec ) ) DEALLOCATE( dbec )
      IF( ALLOCATED( drhovan ) ) DEALLOCATE( drhovan )
  end subroutine deallocate_cdvan
end module cdvan
