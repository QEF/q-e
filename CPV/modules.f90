!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!     n     = total number of electronic states
!     nx    = if n is even, nx=n ; if it is odd, nx=n+1
!             nx is used only to dimension arrays

!     tpiba   = 2*pi/alat
!     tpiba2  = (2*pi/alat)**2
!     ng      = number of G vectors for density and potential
!     ngl     = number of shells of G

!     G-vector quantities for the thick grid - see also doc in ggen 
!     g2      = G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
!     gl      = shells of G^2           ( "   "   "    "      "      )
!     g       = G-vectors               ( "   "   "  tpiba =(2pi/a)  )
!     g2_g    = all G^2 in increasing order, replicated on all procs
!
!     mill    = miller index of G vecs (local to the processors)
!     mill_g  = miller index of G vecs (increasing order), replicated on all procs
!     ig_l2g  = "l2g" means local to global, this array convert a local
!               G-vector index into the global index, in other words
!               the index of the G-v. in the overall array of G-vectors
!     bi?     = base vector used to generate the reciprocal space
!
!     np      = fft index for G>
!     nm      = fft index for G<
!


!
!  lqmax:  maximum angular momentum of Q (Vanderbilt augmentation charges)
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
  use kinds, only: DP
  implicit none
  save
  !
  !    rhops = ionic pseudocharges (for Ewald term)
  !    vps   = local pseudopotential in G space for each species
  !
  real(DP), allocatable:: rhops(:,:), vps(:,:)
  !
  !    drhops = derivative of rhops respect to G^2
  !    dvps   = derivative of vps respect to G^2
  !
  real(DP),allocatable:: dvps(:,:), drhops(:,:)
  !
  !    vps0  = correction factors needed to align V(0) to the "traditional"
  !            value used by other plane-wave codes - one per species
  !
  real(DP),allocatable:: vps0(:)
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
      ALLOCATE( vps0( nsp ) )
  end subroutine
  !
  subroutine deallocate_local_pseudo
      IF( ALLOCATED( vps0 ) ) DEALLOCATE( vps0 )
      IF( ALLOCATED( dvps ) ) DEALLOCATE( dvps )
      IF( ALLOCATED( drhops ) ) DEALLOCATE( drhops )
      IF( ALLOCATED( vps ) ) DEALLOCATE( vps )
      IF( ALLOCATED( rhops ) ) DEALLOCATE( rhops )
  end subroutine
  !
end module local_pseudo

module qgb_mod
  USE kinds, ONLY: DP
  implicit none
  save
  complex(DP), allocatable :: qgb(:,:,:)
  complex(DP), allocatable :: dqgb(:,:,:,:,:)
contains
  subroutine deallocate_qgb_mod
      IF( ALLOCATED( qgb ) ) DEALLOCATE( qgb )
      IF( ALLOCATED( dqgb ) ) DEALLOCATE( dqgb )
  end subroutine deallocate_qgb_mod
end module qgb_mod


MODULE metagga  !metagga
  USE kinds, ONLY: DP
  implicit none
  !the variables needed for meta-GGA
  REAL(DP), ALLOCATABLE :: &
       kedtaus(:,:), &! KineticEnergyDensity in real space,smooth grid
       kedtaur(:,:), &! real space, density grid
       crosstaus(:,:,:), &!used by stress tensor,in smooth grid
       dkedtaus(:,:,:,:)  !derivative of kedtau wrt h on smooth grid
  COMPLEX(DP) , ALLOCATABLE :: &
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
   SUBROUTINE allocate_core( nrxx, ngm, ngb, nsp ) 
     INTEGER, INTENT(IN) :: nrxx, ngm, ngb, nsp
     IF ( nlcc_any ) THEN    
        !
        ALLOCATE( rhoc( nrxx ) )
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
!
