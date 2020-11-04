!
! Copyright (C) 2003-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------------
MODULE path_reparametrisation
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the reparametrisation of the path in the string method
  !
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  !
  USE kinds,     ONLY : DP
  USE path_io_units_module,  ONLY : iunpath
  USE io_global, ONLY : meta_ionode, meta_ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  !
  USE basic_algebra_routines
  !
  PRIVATE
  !
  PUBLIC :: reparametrise, spline_interpolation
  !
  INTERFACE spline_interpolation
     !
     MODULE PROCEDURE spline_interpolation_1D, spline_interpolation_2D
     !
  END INTERFACE
  !
  CONTAINS
    !
    ! ... reparametrisation routines in real space
    !
    !------------------------------------------------------------------------
    SUBROUTINE reparametrise()
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : pos
      USE path_variables, ONLY : nim => num_of_images
      USE path_variables, ONLY : climbing
      !
      IMPLICIT NONE
      !
      INTEGER  :: i, ni, nf
      !
      !
      IF ( meta_ionode ) THEN
         !
         IF ( ANY( climbing(:) ) ) THEN
            !
            ni = 1
            !
            DO i = 2, nim
               !
               IF ( .NOT. climbing(i) ) CYCLE
               !
               nf = i
               !
               CALL spline_interpolation( pos, ni, nf )
               !
               ni = nf
               !
            END DO
            !
            nf = nim
            !
            CALL spline_interpolation( pos, ni, nf )
            !
         ELSE
            !
            ni = 1
            nf = nim
            !
            CALL spline_interpolation( pos, ni, nf )
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( pos, meta_ionode_id, world_comm )
      !
      RETURN
      !
    END SUBROUTINE reparametrise
    !
    !--------------------------------------------------------------------
    SUBROUTINE spline_interpolation_1D( vec, ni, nf, nim )
      !--------------------------------------------------------------------
      !
      USE splinelib, ONLY : dosplineint
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(INOUT)        :: vec(:)
      INTEGER,  INTENT(IN)           :: ni, nf
      INTEGER,  INTENT(IN), OPTIONAL :: nim
      !
      INTEGER                :: i, j
      INTEGER                :: nio, nfo
      REAL(DP)               :: delta, length
      REAL(DP), ALLOCATABLE  :: new_vec(:)
      REAL(DP), ALLOCATABLE  :: old_mesh(:), new_mesh(:)
      !
      !
      IF ( PRESENT( nim ) ) THEN
         !
         nio = 1
         nfo = nim
         !
      ELSE
         !
         nio = ni
         nfo = nf
         !
      END IF
      !
      ! ... cubic spline interpolation
      !
      ALLOCATE( new_vec( ni:nf ) )
      !
      ALLOCATE( old_mesh( nio:nfo ) )
      ALLOCATE( new_mesh( ni:nf ) )
      !
      old_mesh(:) = 0.0_DP
      new_mesh(:) = 0.0_DP
      !
      DO i = nio, nfo - 1
         !
         old_mesh(i+1) = old_mesh(i) + ABS( vec(i+1) - vec(i) )
         !
      END DO
      !
      length = old_mesh(nfo)
      !
      delta = length / DBLE( nf - ni )
      !
      DO j = 0, nf - ni
         !
         new_mesh(j+ni) = DBLE(j) * delta
         !
      END DO
      !
      old_mesh(:) = old_mesh(:) / length
      new_mesh(:) = new_mesh(:) / length
      !
      CALL dosplineint( old_mesh(:), vec(nio:nfo), new_mesh(:), new_vec(:) )
      !
      vec(ni:nf) = new_vec(:)
      !
      DEALLOCATE( new_vec, old_mesh, new_mesh )
      !
      RETURN
      !
    END SUBROUTINE spline_interpolation_1D
    !
    !--------------------------------------------------------------------
    SUBROUTINE spline_interpolation_2D( vec, ni, nf, nim )
      !--------------------------------------------------------------------
      !
      USE splinelib, ONLY : dosplineint
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(INOUT)        :: vec(:,:)
      INTEGER,  INTENT(IN)           :: ni, nf
      INTEGER,  INTENT(IN), OPTIONAL :: nim
      !
      INTEGER                :: i, j
      INTEGER                :: nio, nfo
      INTEGER                :: dim1
      REAL(DP)               :: delta, length
      REAL(DP), ALLOCATABLE  :: new_vec(:,:)
      REAL(DP), ALLOCATABLE  :: old_mesh(:), new_mesh(:)
      !
      !
      dim1 = SIZE( vec, 1 )
      !
      IF ( PRESENT( nim ) ) THEN
         !
         nio = 1
         nfo = nim
         !
      ELSE
         !
         nio = ni
         nfo = nf
         !
      END IF
      !
      ! ... cubic spline interpolation
      !
      ALLOCATE( new_vec( dim1, ni:nf ) )
      !
      ALLOCATE( old_mesh( nio:nfo ) )
      ALLOCATE( new_mesh( ni:nf ) )
      !
      old_mesh(:) = 0.0_DP
      new_mesh(:) = 0.0_DP
      !
      DO i = nio, nfo - 1
         !
         old_mesh(i+1) = old_mesh(i) + norm( vec(:,i+1) - vec(:,i) )
         !
      END DO
      !
      length = old_mesh(nfo)
      !
      delta = length / DBLE( nf - ni )
      !
      DO j = 0, nf - ni
         !
         new_mesh(j+ni) = DBLE(j) * delta
         !
      END DO
      !
      old_mesh(:) = old_mesh(:) / length
      new_mesh(:) = new_mesh(:) / length
      !
      CALL dosplineint( old_mesh(:), vec(:,nio:nfo), new_mesh(:), new_vec(:,:) )
      !
      vec(:,ni:nf) = new_vec(:,:)
      !
      DEALLOCATE( new_vec, old_mesh, new_mesh )
      !
      RETURN
      !
    END SUBROUTINE spline_interpolation_2D
    !
    !--------------------------------------------------------------------
    SUBROUTINE cubic_interpolation( ni, nf )
      !--------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim1, pos
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ni, nf
      !
      INTEGER               :: i, j
      REAL(DP)              :: r, delta, x
      REAL(DP), ALLOCATABLE :: a(:,:), b(:,:), c(:,:), d(:,:), t(:,:), s(:)
      !
      ALLOCATE( a( dim1, ni:nf-1 ) )
      ALLOCATE( b( dim1, ni:nf-1 ) )
      ALLOCATE( c( dim1, ni:nf-1 ) )
      ALLOCATE( d( dim1, ni:nf-1 ) )
      ALLOCATE( t( dim1, ni:nf ) )
      ALLOCATE( s( ni:nf ) )
      !
      t(:,ni) = pos(:,ni+1) - pos(:,ni)
      t(:,nf) = pos(:,nf) - pos(:,nf-1)
      !
      DO i = ni+1, nf - 1
         !
         t(:,i) = ( pos(:,i+1) - pos(:,i-1) ) / 2.0_DP
         !
      END DO
      !
      s(ni) = 0.0_DP
      !
      DO i = ni, nf - 1
         !
         r = norm( pos(:,i+1) - pos(:,i) )
         !
         s(i+1) = s(i) + r
         !
         ! ... cubic interpolation
         !
         a(:,i) = 2.0_DP *( pos(:,i) - pos(:,i+1) ) / r**3 + &
                           ( t(:,i) + t(:,i+1) ) / r**2
         !
         b(:,i) = 3.0_DP *( pos(:,i+1) - pos(:,i) ) / r**2 - &
                      ( 2.0_DP*t(:,i) + t(:,i+1) ) / r
         !
         c(:,i) = t(:,i)
         !
         d(:,i) = pos(:,i)
         !
      END DO
      !
      i = ni
      !
      delta = s(nf) / DBLE( nf - ni )
      !
      DO j = ni, nf
         !
         r = DBLE( j - ni ) * delta
         !
         IF ( r >= s(i+1) .AND. i < nf - 1 ) i = i + 1
         !
         x = r - s(i)
         !
         pos(:,j) = a(:,i)*x**3 + b(:,i)*x**2 + c(:,i)*x + d(:,i)
         !
      END DO
      !
      DEALLOCATE( a, b, c, d, t, s )
      !
      RETURN
      !
    END SUBROUTINE cubic_interpolation
    !
END MODULE path_reparametrisation
