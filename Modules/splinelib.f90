!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE splinelib
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: dosplineint, spline, splint
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE spline( xdata, ydata, startu, startd, d2y )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL (DP), INTENT(IN)  :: xdata(:), ydata(:), startu, startd 
      REAL (DP), INTENT(OUT) :: d2y(:)
      INTEGER                     :: i, k, old_num_of_images
      REAL (DP)              :: p, qn, sig, un
      REAL (DP), ALLOCATABLE :: u(:)
      !
      !
      old_num_of_images = SIZE( ydata )
      !
      allocate(u(old_num_of_images))
      u(1)   = startu
      d2y(1) = startd
      !
      DO  i = 2, ( old_num_of_images - 1 ) 
         !
         sig    = ( xdata(i) - xdata(i - 1) ) / ( xdata(i + 1) - xdata(i - 1) ) 
         p      = sig * d2y(i - 1) + 2.D0 
         d2y(i) = ( sig - 1.D0 ) / p 
         u(i)   = ( 6.D0 * ( (ydata(i + 1) - ydata(i) ) / &
                  ( xdata(i + 1) - xdata(i) ) - ( ydata(i) - ydata(i - 1) ) / &
                  ( xdata(i) - xdata(i - 1) ) ) / &
                  ( xdata(i + 1) - xdata(i - 1) ) - sig * u(i - 1) ) / p 
         !       
      END DO
      !
      d2y(old_num_of_images) = 0  
      !
      DO  k = ( old_num_of_images - 1 ), 1, -1 
         !
         d2y(k) = d2y(k) * d2y(k + 1) + u(k) 
         !
      END DO
      !
      DEALLOCATE (u)
      !
    END SUBROUTINE spline
    !
    !
    !------------------------------------------------------------------------
    FUNCTION splint( xdata, ydata, d2y, x )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL (DP), INTENT(IN)  :: xdata(:), ydata(:), d2y(:)
      REAL (DP), INTENT(IN)  :: x
      REAL (DP)              :: splint
      INTEGER                     :: k, khi, klo, dim
      REAL (DP)              :: a, b, h
      !
      !
      dim = SIZE( xdata )
      klo = 1
      khi = dim
      !
      klo = MAX( MIN( locate( xdata , x ) , ( dim - 1 ) ) , 1 )
      !
      khi = klo + 1
      !
      h = xdata(khi) - xdata(klo)
      !
      a = ( xdata(khi) - x ) / h
      b = ( x - xdata(klo) ) / h
      !
      splint = a * ydata(klo) + b * ydata(khi) + &
               ( ( a**3 - a ) * d2y(klo) + ( b**3 - b ) * d2y(khi) ) * &
               ( h**2 ) / 6.D0
      !
      CONTAINS
         !
         !-------------------------------------------------------------------
         FUNCTION locate( xx , x )
           !-------------------------------------------------------------------
           !
           IMPLICIT NONE
           !
           REAL (DP), INTENT(IN)  :: xx(:)
           REAL (DP), INTENT(IN)  :: x
           INTEGER                     :: locate
           INTEGER                     :: n, jl, jm, ju
           LOGICAL                     :: ascnd
           !
           !
           n     = SIZE( xx )
           ascnd = ( xx(n) >= xx(1) )
           jl    = 0
           ju    = n + 1
           !
           main_loop: DO
              !
              IF ( ( ju - jl ) <= 1 ) EXIT main_loop
              ! 
              jm = ( ju + jl ) / 2
              !
              IF ( ascnd .EQV. ( x >= xx(jm) ) ) THEN
                 !
                 jl = jm
                 !
              ELSE
                 !
                 ju = jm
                 !
              END IF
              !
           END DO main_loop
           !
           IF ( x == xx(1) ) THEN
              !
              locate = 1
              !
           ELSE IF ( x == xx(n) ) THEN
              !
              locate = n - 1
              !
           ELSE 
              !
              locate = jl
              !
           END IF
           !
         END FUNCTION locate      
         !
    END FUNCTION splint
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE dosplineint( old_mesh, old_vect, new_mesh, new_vect )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL (DP), INTENT(IN)   :: old_mesh(:), new_mesh(:)
      REAL (DP), INTENT(IN)   :: old_vect(:,:)
      REAL (DP), INTENT(OUT)  :: new_vect(:,:)
      REAL (DP), ALLOCATABLE  :: d2y(:)
      INTEGER                      :: dim, i, j
      INTEGER                      :: old_num_of_images, new_num_of_images
      !
      !
      dim = SIZE( old_vect , 1 )
      !
      IF( dim /= SIZE( new_vect , 1 ) ) THEN
         !
         WRITE(*,'("ERROR in dosplineint: ", &
                  &"dimensions of old_vect and new_vect")')
         WRITE(*,'("                      do not match")')
         STOP
         !
      END IF  
      !
      old_num_of_images = SIZE( old_vect , 2 )
      new_num_of_images = SIZE( new_vect , 2 )
      !
      IF ( old_num_of_images /= SIZE( old_mesh , 1 ) ) THEN
         !
         WRITE(*,'("ERROR in dosplineint: ", &
                  &"dimensions of old_mesh and old_vect")')
         WRITE(*,'("                      do not match")')
         STOP
         !
      ELSE IF( new_num_of_images /= SIZE( new_mesh , 1 ) ) THEN
         !
         WRITE(*,'("ERROR in dosplineint: ", &
                  &"dimensions of new_mesh and new_vect")')
         WRITE(*,'("                      do not match")')
         STOP
         !
      END IF
      !
#if defined ( _DEBUG ) || ( _DEBUG_SPLINELIB )
      !
      PRINT *, dim            
      PRINT *, old_num_of_images, new_num_of_images
      !
#endif      
      !
      ALLOCATE( d2y( old_num_of_images ) )
      !
      DO i = 1, dim
         ! 
         d2y = 0
         !
         CALL spline( old_mesh , old_vect(i,:), 0.d0, 0.d0, d2y  ) 
         !
         DO j = 1, new_num_of_images
            !
            new_vect(i,j) = splint( old_mesh, old_vect(i,:), d2y, new_mesh(j) )
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( d2y )
      !
    END SUBROUTINE dosplineint
    !
END MODULE splinelib
