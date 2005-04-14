!
! Copyright (C) 2002-2004 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

MODULE cp_main_variables

  USE kinds, ONLY: dbl
  USE parameters, ONLY: natx, nsx, nacx
  
  IMPLICIT NONE
  SAVE

!
! structure factors e^{-ig*R}
!
      complex(dbl), allocatable:: ei1(:,:),  ei2(:,:),  ei3(:,:)
      complex(dbl), allocatable:: eigr(:,:)
!
! structure factors (summed over atoms of the same kind)
!
      complex(dbl), allocatable:: sfac(:,:)
!
! indexes, positions, and structure factors for the box grid
!
      real(dbl)  :: taub(3,natx)
      complex(dbl), allocatable:: eigrb(:,:)
      integer,      allocatable:: irb(:,:)
! 
! charge densities and potentials
!     rhog  = charge density in g space
!     rhor  = charge density in r space (dense grid)
!     rhos  = charge density in r space (smooth grid)
!     rhoc  = core charge density in real space (dense grid)
!
      complex(dbl), allocatable:: rhog(:,:)
      real(dbl),    allocatable:: rhor(:,:), rhos(:,:), rhoc(:)
!
! nonlocal projectors:
!     bec   = scalar product of projectors and wave functions
!     betae = nonlocal projectors in g space = beta x e^(-ig.R) 
!     becdr = <betae|g|psi> used in force calculation
!     rhovan= \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
!     deeq  = \int V_eff(r) q_lm(r) dr
!
      real(dbl), allocatable:: bec(:,:), becdr(:,:,:)
      real(dbl), allocatable:: bephi(:,:), becp(:,:)
!
!  mass preconditioning
!
      real(dbl), allocatable:: ema0bg(:)
!
!  constraints (lambda at t, lambdam at t-dt, lambdap at t+dt)
!
      real(dbl), allocatable:: lambda(:,:), lambdam(:,:), lambdap(:,:)
!
      real(dbl) acc(nacx)


CONTAINS


   SUBROUTINE allocate_mainvar &
      ( ngw, ngb, ngs, ng, nr1, nr2, nr3, nnr, nnrsx, nat, nax, nsp, nspin, n, nx, nhsa, &
        nlcc_any, smd )
      INTEGER, INTENT(IN) :: ngw, ngb, ngs, ng, nr1, nr2, nr3, nnr, nnrsx, nat, nax, nsp, nspin, n, nx, nhsa
      LOGICAL, INTENT(IN) :: nlcc_any
      LOGICAL, OPTIONAL, INTENT(IN) :: smd
      !
      !     allocation of all arrays not already allocated in init and nlinit
      !
      LOGICAL :: nosmd
      
      nosmd = .TRUE.
      IF( PRESENT( smd ) ) THEN
        IF( smd ) nosmd = .FALSE.
      END IF

      allocate(eigr(ngw,nat))
      allocate(eigrb(ngb,nat))
      allocate(irb(3,nat))
      allocate(sfac(ngs,nsp))
      allocate(ei1(-nr1:nr1,nat))
      allocate(ei2(-nr2:nr2,nat))
      allocate(ei3(-nr3:nr3,nat))
      allocate(rhor(nnr,nspin))
      allocate(rhos(nnrsx,nspin))
      allocate(rhog(ng,nspin))
      if ( nlcc_any ) allocate(rhoc(nnr))
      allocate(ema0bg(ngw))
      IF( nosmd ) THEN
        allocate(lambda(nx,nx))
        allocate(lambdam(nx,nx))
        allocate(lambdap(nx,nx))
      END IF
      allocate(becdr(nhsa,n,3))
      IF( nosmd ) THEN
        allocate(bec  (nhsa,n))
      END IF
      allocate(bephi(nhsa,n))
      allocate(becp (nhsa,n))
      RETURN
   END SUBROUTINE


   SUBROUTINE deallocate_mainvar( )

      IF( ALLOCATED( ei1 ) ) DEALLOCATE( ei1 )
      IF( ALLOCATED( ei2 ) ) DEALLOCATE( ei2 )
      IF( ALLOCATED( ei3 ) ) DEALLOCATE( ei3 )
      IF( ALLOCATED( eigr ) ) DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) ) DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) ) DEALLOCATE( eigrb )
      IF( ALLOCATED( irb ) ) DEALLOCATE( irb )
      IF( ALLOCATED( rhor ) ) DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) ) DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) ) DEALLOCATE( rhog )
      IF( ALLOCATED( rhoc ) ) DEALLOCATE( rhoc )
      IF( ALLOCATED( bec ) ) DEALLOCATE( bec )
      IF( ALLOCATED( becdr ) ) DEALLOCATE( becdr )
      IF( ALLOCATED( bephi ) ) DEALLOCATE( bephi )
      IF( ALLOCATED( becp ) ) DEALLOCATE( becp )
      IF( ALLOCATED( ema0bg ) ) DEALLOCATE( ema0bg )
      IF( ALLOCATED( lambda ) ) DEALLOCATE( lambda )
      IF( ALLOCATED( lambdam ) ) DEALLOCATE( lambdam )
      IF( ALLOCATED( lambdap ) ) DEALLOCATE( lambdap )

      RETURN
   END SUBROUTINE

END MODULE
