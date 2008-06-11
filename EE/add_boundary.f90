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
      SUBROUTINE add_boundary( rho, pot, ixmax, iymax, izmax,          &
                               deltax, deltay, deltaz )
!--------------------------------------------------------------------
      !
      ! ... Adds the Dirichlet boundary conditions corresponding to
      ! ... the charge density rho
      !
      USE kinds,             ONLY : DP
      USE constants,         ONLY : fpi, e2
      USE cell_base,         ONLY : alat, omega, at
      USE mp_global,         ONLY : me_pool, intra_pool_comm
      USE mp,                ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER                    :: ixmax
      INTEGER                    :: iymax
      INTEGER                    :: izmax
      REAL( DP )                 :: deltax
      REAL( DP )                 :: deltay
      REAL( DP )                 :: deltaz
      REAL( DP )                 :: pot( ixmax * iymax * izmax )
      REAL( DP )                 :: rho( ixmax * iymax * izmax )
      !
      INTEGER                    :: i, ii, ifirst, ilast
      INTEGER                    :: ix
      INTEGER                    :: iy
      INTEGER                    :: iz 
      INTEGER                    :: j
      INTEGER                    :: jx
      INTEGER                    :: jy
      INTEGER                    :: jz
      INTEGER                    :: na
      INTEGER                    :: nt
      !
      REAL( DP )                 :: dist
      REAL( DP )                 :: fact
      REAL( DP )                 :: tx
      REAL( DP )                 :: ty
      REAL( DP )                 :: tz
      REAL( DP )                 :: vg
      !
      REAL( DP ), PARAMETER      :: vanishing_density = 1.D-7
      REAL( DP ), PARAMETER      :: vanishing_dist = 1.D-5
!! DCC 
      REAL( DP ), ALLOCATABLE    :: aux (:)
      INTEGER, EXTERNAL           :: compindex
      !
      CALL start_clock( 'boundary' )
      !

      fact = deltax * deltay * deltaz * e2
      !
      !
!     DO iz = 1, izmax
!     DO iy = 1, iymax
!     DO ix = 1, ixmax
!       i = compindex( ix, iy, iz, ixmax, iymax, izmax )
!
! this should be equivalent
!
      allocate (aux(ixmax*iymax*izmax))
      aux(:) = 0.d0
      call divide (ixmax*iymax*izmax, ifirst, ilast)
      DO i =ifirst, ilast 
         ii = i - 1
         iz = ii / (ixmax*iymax) + 1
         ii = ii - (ixmax*iymax) * (iz-1)
         iy = ii / ixmax + 1
         ii = ii - ixmax * (iy-1)
         ix = ii + 1

        IF( ABS( rho( i ) ) > vanishing_density ) THEN
          DO jz = 1, izmax
          DO jy = 1, iymax
          DO jx = 1, ixmax
            IF( ( jx == 1 ) .OR. ( jx == ixmax ) .OR.                  &
                ( jy == 1 ) .OR. ( jy == iymax ) .OR.                  &
                ( jz == 1 ) .OR. ( jz == izmax ) ) THEN
              j = compindex( jx, jy, jz, ixmax, iymax, izmax )
              dist =  SQRT( deltax ** 2 * DBLE( ix - jx ) ** 2         &
                          + deltay ** 2 * DBLE( iy - jy ) ** 2         &
                          + deltaz ** 2 * DBLE( iz - jz ) ** 2 ) 
                IF( dist > vanishing_dist )                            &
                  aux( j ) = aux( j ) + rho( i ) * fact / dist 
            END IF
          END DO
          END DO
          END DO
        END IF
!      END DO
!      END DO
      END DO

#ifdef __PARA
      call mp_sum(  aux, intra_pool_comm )
#endif
      pot(:) = pot(:) + aux(:)
      deallocate (aux)
#ifdef DEBUG
      do i=1, ixmax*iymax*izmax
         write (80+me_pool, *)  rho(i), pot(i)
      end do
      write (80+me_pool,*) "END "
#endif
      !
      CALL stop_clock( 'boundary' )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE add_boundary
!--------------------------------------------------------------------

