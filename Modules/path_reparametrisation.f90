!
! Copyright (C) 2003-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!---------------------------------------------------------------------------
MODULE path_reparametrisation
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the reparametrisation of the path in the smd method
  !
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
  !
  USE io_files,  ONLY : iunpath
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi
  !
  USE basic_algebra_routines
  ! moved here to work around a mips compiler bug
  USE input_parameters, ONLY : num_of_images_inp => num_of_images
  !
  CONTAINS
    !
    ! ... module procedures    
    !
    ! ... in direct space
    !
    !------------------------------------------------------------------------
    SUBROUTINE spline_reparametrisation()
      !------------------------------------------------------------------------
      !
      USE splinelib,        ONLY : dosplineint
      USE path_variables,   ONLY : pos, num_of_images, dim
      USE path_variables,   ONLY : istep_path, path_thr, use_multistep, &
                                   err_max, frozen, vel
      USE io_global,        ONLY : meta_ionode
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER  :: i, j
      REAL(DP) :: delta_R
      REAL(DP) :: change_image_thr
      REAL(DP) :: multistep_coeff = 2.D0
      INTEGER  :: new_num_of_images
      INTEGER  :: N_in, N_fin
      INTEGER  :: init_num_of_images = 3
      !
      REAL(DP), ALLOCATABLE  :: new_pos(:,:)
      REAL(DP), ALLOCATABLE  :: old_mesh(:), new_mesh(:)
      !
      !
      new_num_of_images = num_of_images
      !
      IF ( use_multistep ) THEN
         !
         change_image_thr = multistep_coeff * &
                            DBLE( num_of_images_inp - num_of_images ) * path_thr
         !
         IF ( istep_path == 0 ) THEN
            !
            ! ... initialisation
            !
            IF ( meta_ionode ) &
               WRITE( UNIT = iunpath, &
                    & FMT = '(5X,"initial number of images = ",I3,/)' ) &
                   init_num_of_images
            !
            new_num_of_images = init_num_of_images
            !
            CALL redispose_last_image( new_num_of_images )
            !
         ELSE IF ( err_max < change_image_thr ) THEN
            !
            new_num_of_images = MIN( num_of_images_inp, num_of_images + 2 )
            !
            IF ( new_num_of_images > num_of_images ) THEN
               !
               IF ( meta_ionode ) &
                  WRITE( UNIT = iunpath, &
                       & FMT = '(5X,"new number of images = ",I3,/)' ) &
                      new_num_of_images
               !
               CALL redispose_last_image( new_num_of_images )
               !
               N_in  = 2
               N_fin = new_num_of_images - 1
               !
               vel(:,N_in:N_fin) = 0.D0
               !
               frozen(N_in:N_fin) = .FALSE.
               !
            END IF
            !
         END IF
         !
      END IF
      !
      ! ... cubic spline reparametrisation
      !
      ALLOCATE( new_pos( dim, new_num_of_images ) )
      !
      ALLOCATE( old_mesh( num_of_images ) )
      ALLOCATE( new_mesh( new_num_of_images ) )
      !
      old_mesh(1) = 0.D0
      !
      DO i = 1, ( num_of_images - 1 )
         !
         old_mesh(i+1) = old_mesh(i) + norm( pos(:,i+1) - pos(:,i) )
         !
      END DO
      !
      delta_R = old_mesh(num_of_images) / DBLE( new_num_of_images - 1 )
      ! 
      DO j = 0, ( new_num_of_images - 1 )
         !
         new_mesh(j+1) = DBLE(j) * delta_R 
         !
      END DO
      !
      new_mesh = new_mesh / old_mesh(num_of_images)  
      old_mesh = old_mesh / old_mesh(num_of_images)
      !
      CALL dosplineint( old_mesh , pos(:,1:num_of_images) , new_mesh , new_pos )
      !
      num_of_images = new_num_of_images
      !
      pos(:,1:num_of_images) = new_pos(:,1:num_of_images)
      !
      DEALLOCATE( new_pos, old_mesh, new_mesh )
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE redispose_last_image( n )
          !--------------------------------------------------------------------
          !
          USE path_variables, ONLY : error, grad, pes, grad_pes
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: n
          !
          !
          pes(n)        = pes(num_of_images)
          grad_pes(:,n) = grad_pes(:,num_of_images)
          error(n)      = error(num_of_images)
          vel(:,n)      = vel(:,num_of_images)
          grad(:,n)     = grad(:,num_of_images)
          !
          RETURN
          !
        END SUBROUTINE redispose_last_image
      !
    END SUBROUTINE spline_reparametrisation
    !
    ! ... in reciprocal space
    !
    !-----------------------------------------------------------------------
    SUBROUTINE update_num_of_images()
      !-----------------------------------------------------------------------
      !
      USE path_variables,   ONLY : istep_path, num_of_images, &
                                   path_thr, Nft, ft_coeff, pos, pes, &
                                   use_multistep, grad_pes, err_max,  &
                                   frozen, vel
      USE io_global,        ONLY : meta_ionode
      !
      IMPLICIT NONE
      !
      REAL(DP) :: change_image_thr
      REAL(DP) :: multistep_coeff = 2.D0
      INTEGER  :: new_num_of_images
      LOGICAL  :: images_updated
      INTEGER  :: N_in, N_fin
      INTEGER  :: init_num_of_images = 3
      !
      !
      IF ( .NOT. use_multistep ) RETURN
      !
      images_updated = .FALSE.
      !
      change_image_thr = multistep_coeff * &
                         DBLE( num_of_images_inp - num_of_images ) * path_thr
      !
      IF ( istep_path == 0 ) THEN
         !
         ! ... initialisation
         !
         IF ( meta_ionode ) &
            WRITE( UNIT = iunpath, &
                 & FMT = '(5X,"initial number of images = ",I3,/)' ) &
                init_num_of_images
         !
         CALL redispose_last_image( init_num_of_images )
         !
         num_of_images = init_num_of_images
         !
         images_updated = .TRUE.
         !
      ELSE IF ( err_max < change_image_thr ) THEN
         !
         new_num_of_images = MIN( num_of_images_inp, num_of_images + 2 )
         !
         IF ( new_num_of_images > num_of_images ) THEN
            !
            IF ( meta_ionode ) &
               WRITE( UNIT = iunpath, &
                    & FMT = '(5X,"new number of images = ",I3,/)' ) &
                   new_num_of_images
            !
            CALL redispose_last_image( new_num_of_images )
            !
            N_in  = 2
            N_fin = new_num_of_images - 1
            !
            vel(:,N_in:N_fin) = 0.D0
            !
            frozen(N_in:N_fin) = .FALSE.
            !
            images_updated = .TRUE.
            !
            num_of_images = new_num_of_images
            !
         END IF
         !
      END IF
      !
      IF ( images_updated ) THEN
         !
         ! ... reciprocal space dimensions updated
         !
         Nft = ( num_of_images - 1 )
         !
         ft_coeff = 2.D0 / DBLE( Nft )
         !
      END IF
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE redispose_last_image( n )
          !--------------------------------------------------------------------
          !
          USE path_variables, ONLY : error, grad
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: n
          !
          !
          pos(:,n)      = pos(:,num_of_images)
          pes(n)        = pes(num_of_images)
          grad_pes(:,n) = grad_pes(:,num_of_images)
          error(n)      = error(num_of_images)
          vel(:,n)      = vel(:,num_of_images)
          grad(:,n)     = grad(:,num_of_images)
          !
          RETURN
          !
        END SUBROUTINE redispose_last_image
        !
    END SUBROUTINE update_num_of_images
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_path_length()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim, path_length, pos, ft_pos, Nft, &
                                 Nft_smooth ,num_of_images, num_of_modes
      !
      IMPLICIT NONE
      !
      REAL(DP), ALLOCATABLE :: r_h(:), r_n(:), delta_pos(:)
      REAL(DP)              :: x, delta_x
      INTEGER               :: i, j, n
      !
      !
      ALLOCATE( r_h(       dim ) )
      ALLOCATE( r_n(       dim ) )
      ALLOCATE( delta_pos( dim ) )
      !
      delta_pos(:) = ( pos(:,num_of_images) - pos(:,1) )
      !
      path_length = 0.D0
      !
      r_h(:) = pos(:,1)
      !
      delta_x = 1.D0 / DBLE( Nft_smooth * Nft )
      !
      DO i = 1, Nft
         !
         DO j = 1, Nft_smooth
            !
            x = delta_x * DBLE( Nft_smooth * ( i - 1 ) + j )
            !
            r_n(:) = pos(:,1) + x * delta_pos(:)
            !
            DO n = 1, num_of_modes
               !
               r_n(:) = r_n(:) + ft_pos(:,n) * SIN( DBLE( n ) * pi * x )
               !
            END DO
            !
            path_length = path_length + norm( r_n - r_h )
            !
            r_h(:) = r_n(:)
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( r_h )
      DEALLOCATE( r_n )
      DEALLOCATE( delta_pos )
      !
      RETURN
      !
    END SUBROUTINE compute_path_length
    !
    !-----------------------------------------------------------------------
    SUBROUTINE to_real_space()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_modes, num_of_images, dim, &
                                 pos, ft_pos, Nft, Nft_smooth, path_length
      !
      IMPLICIT NONE
      !
      REAL(DP), ALLOCATABLE :: r_h(:), r_n(:), delta_pos(:)
      REAL(DP)              :: x, delta_x, s, s_image, pi_n
      INTEGER               :: i, j, n, image
      !
      !
      ALLOCATE( r_h(       dim ) )
      ALLOCATE( r_n(       dim ) )
      ALLOCATE( delta_pos( dim ) )
      !
      delta_pos(:) = ( pos(:,num_of_images) - pos(:,1) )
      !
      s = 0.D0
      !
      image = 1
      !
      s_image = path_length / DBLE( Nft )
      !
      r_h(:) = pos(:,1)
      !
      delta_x = 1.D0 / DBLE( Nft_smooth * Nft )
      !
      DO i = 1, Nft
         !
         DO j = 1, Nft_smooth
            !
            x = delta_x * DBLE( Nft_smooth * ( i - 1 ) + j )
            !
            r_n(:) = pos(:,1) + x * delta_pos(:)
            !
            DO n = 1, num_of_modes
               !
               r_n(:) = r_n(:) + ft_pos(:,n) * SIN( DBLE( n ) * pi * x )
               !
            END DO
            !
            s = s + norm( r_n - r_h )
            !
            IF ( s >= s_image ) THEN
               !
               image = image + 1
               !
               pos(:,image) = r_n(:)
               !
               s_image = s_image + path_length / DBLE( Nft )
               !
            END IF
            !
            r_h(:) = r_n(:)
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( r_h )
      DEALLOCATE( r_n )
      DEALLOCATE( delta_pos )
      !
      RETURN
      !
    END SUBROUTINE to_real_space
    !
    !------------------------------------------------------------------------
    SUBROUTINE to_reciprocal_space( f, ft_f )
      !------------------------------------------------------------------------
      !
      ! ... the function that is fourier transformed is defined here :
      !
      ! ... f(j), j = 1,...,N   is the function in real space (it starts from 1)
      !
      ! ... where :  f(1) /= f(N)
      !
      ! ... f_star(j) = f(j+1) - f(0) + j/(N-1)*( f(N) - f(0) )
      !
      ! ... with index running in  j = 0,..., N-1  and :
      !
      ! ... f_star(0) = f_star(N-1) = 0   so that f_star(j) is stroed between
      ! ...                               0  and  Nft - 1 ( = N - 2 )
      !
      USE path_variables, ONLY : nim => num_of_images
      USE path_variables, ONLY : dim, num_of_modes, Nft, Nft_smooth, ft_coeff
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)  :: f(:,:)
      REAL(DP), INTENT(OUT) :: ft_f(:,:)
      !
      INTEGER               :: j, n
      REAL(DP), ALLOCATABLE :: f_star(:,:)
      REAL(DP)              :: x, coeff, inv_Nft
      !
      !
      ALLOCATE( f_star( dim, 0:( Nft - 1 ) ) )
      !
      inv_Nft = 1.D0 / DBLE( Nft )
      !
      DO j = 0, ( Nft - 1 )
         !
         x = DBLE( j ) * inv_Nft
         !
         f_star(:,j) = f(:,j+1) - f(:,1) - x * ( f(:,nim) - f(:,1) )
         !
      END DO
      !
      ! ... fourier components of f_star are computed
      !
      ft_f = 0.D0
      !
      DO n = 1, num_of_modes
         !
         coeff = DBLE( n ) * pi * inv_Nft
         !
         DO j = 0, ( Nft - 1 )
            !
            x = DBLE( j )
            !
            ft_f(:,n) = ft_f(:,n) + f_star(:,j) * SIN( coeff * x )
            !
         END DO
         !
      END DO
      !
      ! ... normalisation
      !
      ft_f = ft_f * ft_coeff
      !
      DEALLOCATE( f_star )
      !
      RETURN
      !
    END SUBROUTINE to_reciprocal_space
    !
END MODULE path_reparametrisation
