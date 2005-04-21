!
! Copyright (C) 2002-2004 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE cp_main_variables
  !----------------------------------------------------------------------------
  !
  USE kinds,      ONLY : dbl
  USE parameters, ONLY : natx, nsx, nacx
  !
  IMPLICIT NONE
  SAVE
  !
  ! ... structure factors e^{-ig*R}
  !
  COMPLEX(dbl), ALLOCATABLE:: ei1(:,:), ei2(:,:), ei3(:,:)
  COMPLEX(dbl), ALLOCATABLE:: eigr(:,:)
  !
  ! ... structure factors (summed over atoms of the same kind)
  !
  COMPLEX(dbl), ALLOCATABLE:: sfac(:,:)
  !
  ! ... indexes, positions, and structure factors for the box grid
  !
  REAL(dbl)                :: taub(3,natx)
  COMPLEX(dbl), ALLOCATABLE:: eigrb(:,:)
  INTEGER,      ALLOCATABLE:: irb(:,:)
  ! 
  ! ... charge densities and potentials
  ! ...    rhog  = charge density in g space
  ! ...    rhor  = charge density in r space (dense grid)
  ! ...    rhos  = charge density in r space (smooth grid)
  ! ...    rhoc  = core charge density in real space (dense grid)
  !
  COMPLEX(dbl), ALLOCATABLE:: rhog(:,:)
  REAL(dbl),    ALLOCATABLE:: rhor(:,:), rhos(:,:), rhoc(:)
  !
  ! ... nonlocal projectors:
  ! ...    bec   = scalar product of projectors and wave functions
  ! ...    betae = nonlocal projectors in g space = beta x e^(-ig.R) 
  ! ...    becdr = <betae|g|psi> used in force calculation
  ! ...    rhovan= \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  ! ...    deeq  = \int V_eff(r) q_lm(r) dr
  !
  REAL(dbl), ALLOCATABLE:: bec(:,:), becdr(:,:,:)
  REAL(dbl), ALLOCATABLE:: bephi(:,:), becp(:,:)
  !
  ! ... mass preconditioning
  !
  REAL(dbl), ALLOCATABLE:: ema0bg(:)
  !
  ! ... constraints (lambda at t, lambdam at t-dt, lambdap at t+dt)
  !
  REAL(dbl), ALLOCATABLE:: lambda(:,:), lambdam(:,:), lambdap(:,:)
  !
  REAL(dbl) :: acc(nacx)
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_mainvar( ngw, ngb, ngs, ng, nr1, nr2, nr3, nnr, &
                                 nnrsx, nat, nax, nsp, nspin, n, nx, nhsa, &
                                 nlcc_any, smd )
      !------------------------------------------------------------------------
      !
      INTEGER,           INTENT(IN) :: ngw, ngb, ngs, ng, nr1, nr2, nr3, &
                                       nnr, nnrsx, nat, nax, nsp, nspin, &
                                       n, nx, nhsa
      LOGICAL,           INTENT(IN) :: nlcc_any
      LOGICAL, OPTIONAL, INTENT(IN) :: smd
      LOGICAL                       :: nosmd
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      nosmd = .TRUE.
      !
      IF ( PRESENT( smd ) ) THEN
         !
         IF( smd ) nosmd = .FALSE.
         !
      END IF
      !
      ALLOCATE( eigr(  ngw, nat ) )
      ALLOCATE( eigrb( ngb, nat ) )
      ALLOCATE( irb( 3, nat ) )
      ALLOCATE( sfac( ngs, nsp ) )
      ALLOCATE( ei1( -nr1:nr1, nat ) )
      ALLOCATE( ei2( -nr2:nr2, nat ) )
      ALLOCATE( ei3( -nr3:nr3, nat ) )
      ALLOCATE( rhor( nnr,   nspin ) )
      ALLOCATE( rhos( nnrsx, nspin ) )
      ALLOCATE( rhog( ng,    nspin ) )
      !
      IF ( nlcc_any ) THEN
         !
         ALLOCATE( rhoc( nnr ) )
         !
      ELSE
         !
         ! ... dummy allocation required because this array appears in the
         ! ... list of arguments of some routines
         !
         ALLOCATE( rhoc( 1 ) )
         !
      END IF
      !
      ALLOCATE( ema0bg( ngw ) )
      !
      IF ( nosmd ) THEN
         !
         ALLOCATE( lambda(  nx, nx ) )
         ALLOCATE( lambdam( nx, nx ) )
         ALLOCATE( lambdap( nx, nx ) )
         !
      END IF
      !
      ALLOCATE( becdr( nhsa, n, 3 ) )
      !
      IF ( nosmd ) ALLOCATE( bec( nhsa, n ) )
      !
      ALLOCATE( bephi( nhsa, n ) )
      ALLOCATE( becp(  nhsa, n ) )
      !
      RETURN
      !
    END SUBROUTINE allocate_mainvar
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_mainvar()
      !------------------------------------------------------------------------
      !
      IF( ALLOCATED( ei1 ) )     DEALLOCATE( ei1 )
      IF( ALLOCATED( ei2 ) )     DEALLOCATE( ei2 )
      IF( ALLOCATED( ei3 ) )     DEALLOCATE( ei3 )
      IF( ALLOCATED( eigr ) )    DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) )    DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) )   DEALLOCATE( eigrb )
      IF( ALLOCATED( irb ) )     DEALLOCATE( irb )
      IF( ALLOCATED( rhor ) )    DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) )    DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) )    DEALLOCATE( rhog )
      IF( ALLOCATED( rhoc ) )    DEALLOCATE( rhoc )
      IF( ALLOCATED( bec ) )     DEALLOCATE( bec )
      IF( ALLOCATED( becdr ) )   DEALLOCATE( becdr )
      IF( ALLOCATED( bephi ) )   DEALLOCATE( bephi )
      IF( ALLOCATED( becp ) )    DEALLOCATE( becp )
      IF( ALLOCATED( ema0bg ) )  DEALLOCATE( ema0bg )
      IF( ALLOCATED( lambda ) )  DEALLOCATE( lambda )
      IF( ALLOCATED( lambdam ) ) DEALLOCATE( lambdam )
      IF( ALLOCATED( lambdap ) ) DEALLOCATE( lambdap )
      !
      RETURN
      !
    END SUBROUTINE deallocate_mainvar
    !
END MODULE cp_main_variables
