!
! Copyright (C) 2002 CP90 group
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
  !     rhocb = core charge in G space (box grid)
  logical :: nlcc_any
  real(kind=8), allocatable:: rhocb(:,:)
contains
  subroutine deallocate_core()
      IF( ALLOCATED( rhocb ) ) DEALLOCATE( rhocb )
  end subroutine
end module core

module elct
  use electrons_base, only: nspin, nel, nupdwn, iupdwn
  use electrons_base, only: n => nbnd, nx => nbndx
  use electrons_base, only: f, qbac, ispin => fspin
  use electrons_base, only: deallocate_elct
  implicit none
  save
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
end module elct

module gvec

  use cell_base, only: tpiba, tpiba2
  use reciprocal_vectors, only: &
        gl, g, gx, g2_g, mill_g, mill_l, ig_l2g, igl, bi1, bi2, bi3
  use reciprocal_vectors, only: deallocate_recvecs
  use recvecs_indexes, only: np, nm, deallocate_recvecs_indexes
  use gvecp, only: &
        ng => ngm, &
        ngl => ngml, &
        ng_g => ngmt

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
  implicit none
  save

contains
  subroutine deallocate_gvec
      CALL deallocate_recvecs( )
      CALL deallocate_recvecs_indexes( )
  end subroutine
end module gvec


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

module pseu
  implicit none
  save
  !    rhops = ionic pseudocharges (for Ewald term)
  !    vps   = local pseudopotential in G space for each species
  real(kind=8), allocatable:: rhops(:,:), vps(:,:)
contains
  subroutine deallocate_pseu
      IF( ALLOCATED( rhops ) ) DEALLOCATE( rhops )
      IF( ALLOCATED( vps ) ) DEALLOCATE( vps )
  end subroutine
end module pseu

module qgb_mod
  implicit none
  save
  complex(kind=8), allocatable :: qgb(:,:,:)
contains
  subroutine deallocate_qgb_mod
      IF( ALLOCATED( qgb ) ) DEALLOCATE( qgb )
  end subroutine
end module qgb_mod

module qradb_mod
  implicit none
  save
  real(kind=8), allocatable:: qradb(:,:,:,:,:)
contains
  subroutine deallocate_qradb_mod
      IF( ALLOCATED( qradb ) ) DEALLOCATE( qradb )
  end subroutine
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
  end subroutine
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
  end subroutine
end module dqgb_mod

module dpseu
  implicit none
  save
  real(kind=8),allocatable:: dvps(:,:), drhops(:,:)
contains
  subroutine deallocate_dpseu
      IF( ALLOCATED( dvps ) ) DEALLOCATE( dvps )
      IF( ALLOCATED( drhops ) ) DEALLOCATE( drhops )
  end subroutine
end module dpseu

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
  end subroutine
end module cdvan
