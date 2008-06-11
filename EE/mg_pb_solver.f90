!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!--------------------------------------------------------------------
      SUBROUTINE mg_pb_solver( rhotot, pot, eps, ekap2, ixmax,         &
                               iymax, izmax, deltax, deltay, deltaz )
!--------------------------------------------------------------------
      !
      ! ... Solves the Poisson equation
      ! ... nabla ^ 2 pot = - 4 * pi * e ^ 2 * rho
      ! ... on a three-dimensional mesh using an 
      ! ... multigrid Poisson-Boltzmann solver
      ! ...
      ! ... dielectric constant to be added
      !
      USE kinds,              ONLY: DP     
      USE constants,          ONLY: fpi, e2
      USE ee_mod,             ONLY: whichbc,itmax,nlev,errtol
      !
      IMPLICIT NONE
      !
      INTEGER                     :: ixmax
      INTEGER                     :: iymax
      INTEGER                     :: izmax
      !
      REAL( DP )                  :: pot( ixmax * iymax * izmax )
      REAL( DP )                  :: rhotot( ixmax * iymax * izmax )
      REAL( DP )                  :: eps( ixmax * iymax * izmax )
      REAL( DP )                  :: ekap2( ixmax * iymax * izmax )
      REAL( DP )                  :: deltax
      REAL( DP )                  :: deltay
      REAL( DP )                  :: deltaz
      !
      REAL( DP ), ALLOCATABLE     :: source( : )   
      REAL( DP ), ALLOCATABLE     :: epsmg( :, :, :, : )
      !
      INTEGER                     :: ix, iy, iz
      INTEGER                     :: ir, jr
      !
      INTEGER, EXTERNAL       :: compindex
      !

      ALLOCATE( source( ixmax * iymax * izmax ) )
      ALLOCATE( epsmg( ixmax, iymax, izmax, 3 ) )
      !
      DO ix = 1, ixmax
        DO iy = 1, iymax
          DO iz = 1, izmax
            ir = compindex( ix, iy, iz, ixmax, iymax, izmax )
            jr = compindex( ix + 1, iy, iz, ixmax, iymax, izmax )
            epsmg( ix, iy, iz, 1 ) = 0.5D0 * ( eps( ir ) + eps( jr ) )
            jr = compindex( ix, iy + 1, iz, ixmax, iymax, izmax )
            epsmg( ix, iy, iz, 2 ) = 0.5D0 * ( eps( ir ) + eps( jr ) )
            jr = compindex( ix, iy, iz + 1, ixmax, iymax, izmax )
            epsmg( ix, iy, iz, 3 ) = 0.5D0 * ( eps( ir ) + eps( jr ) )
          END DO
        END DO
      END DO
      !
      CALL start_clock( 'mg_pb' )
      !
      source( : ) = fpi * e2 * rhotot( : )
      !

      CALL mgmain( source, pot, epsmg, ekap2, ixmax, iymax, izmax,     &
                   deltax, deltay, deltaz,whichbc,nlev,errtol,itmax )

      !
      CALL stop_clock( 'mg_pb' )
      !
      DEALLOCATE( source )
      DEALLOCATE( epsmg )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE mg_pb_solver
!--------------------------------------------------------------------
