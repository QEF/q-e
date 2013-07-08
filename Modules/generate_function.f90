!
! Copyright (C) 2006-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! Module to generate functions on the real space dense grid
! Written by Oliviero Andreussi
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
MODULE generate_function
!=----------------------------------------------------------------------=!

  USE kinds, ONLY: DP

  IMPLICIT NONE
 
CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE generate_gaussian( nnr, charge, spread, pos, rho )
!----------------------------------------------------------------------
      !
      USE kinds,            ONLY : DP
      USE constants,        ONLY : sqrtpi
      USE cell_base,        ONLY : at, bg, alat
      USE fft_base,         ONLY : dfftp
      USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)       :: nnr
      REAL( DP ), INTENT(IN)    :: charge, spread
      REAL( DP ), INTENT(IN)    :: pos( 3 )
      REAL( DP ), INTENT(INOUT) :: rho( nnr )
      !
      ! ... Local variables
      !
      INTEGER                   :: i, j, k, ir, ir_end, ip
      INTEGER                   :: index0
      !
      REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
      REAL( DP )                :: scale, spr2, dist
      REAL( DP )                :: r( 3 ), s( 3 )
      REAL( DP ), ALLOCATABLE   :: rholocal ( : )
      !
      inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
      inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
      inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
      !
      index0 = 0
      !
#if defined (__MPI)
      DO i = 1, me_bgrp
        index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
#endif  
      !
#if defined (__MPI)
      ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
      ir_end = nnr
#endif  
      !
      scale = charge / ( sqrtpi * spread )**3
      spr2 = ( spread / alat )**2
      ALLOCATE( rholocal( nnr ) )
      rholocal = 0.D0
      !
      DO ir = 1, ir_end
         !
         ! ... three dimensional indexes
         !
         i = index0 + ir - 1
         k = i / (dfftp%nr1x*dfftp%nr2x)
         i = i - (dfftp%nr1x*dfftp%nr2x)*k
         j = i / dfftp%nr1x
         i = i - dfftp%nr1x*j
         !
         DO ip = 1, 3
            r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                    DBLE( j )*inv_nr2*at(ip,2) + &
                    DBLE( k )*inv_nr3*at(ip,3)
         END DO
         !
         r(:) = pos(:) - r(:) 
         !
         ! ... minimum image convention
         !
         s(:) = MATMUL( r(:), bg(:,:) )
         s(:) = s(:) - ANINT(s(:))
         r(:) = MATMUL( at(:,:), s(:) )
         !
         dist = SUM( r * r ) 
         !
         rholocal( ir ) = scale * EXP(-dist/spr2) 
         !      
      END DO
      !
      rho = rho + rholocal
      DEALLOCATE( rholocal )
      !
      RETURN
      !
!----------------------------------------------------------------------
      END SUBROUTINE generate_gaussian
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE generate_gradgaussian( nnr, charge, spread, pos, gradrho )
!----------------------------------------------------------------------
      !
      USE kinds,            ONLY: DP
      USE constants,        ONLY: sqrtpi
      USE cell_base,        ONLY : at, bg, alat
      USE fft_base,         ONLY : dfftp
      USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)       :: nnr
      REAL( DP ), INTENT(IN)    :: charge, spread
      REAL( DP ), INTENT(IN)    :: pos( 3 )
      REAL( DP ), INTENT(INOUT) :: gradrho( 3, nnr )
      !
      ! ... Local variables
      !
      INTEGER                   :: i, j, k, ir, ir_end, ip
      INTEGER                   :: index0
      !
      REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
      REAL( DP )                :: scale, spr2, dist
      REAL( DP )                :: r( 3 ), s( 3 )
      REAL( DP ), ALLOCATABLE   :: gradrholocal ( :, : )
      !
      inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
      inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
      inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
      !
      index0 = 0
      !
#if defined (__MPI)
      DO i = 1, me_bgrp
        index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
#endif
      !
#if defined (__MPI)
      ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
      ir_end = nnr
#endif  
      !
      scale = 2.d0 * charge / sqrtpi**3 / spread**5
      spr2 = ( spread / alat )**2
      ALLOCATE( gradrholocal( 3, nnr ) )
      gradrholocal = 0.D0
      !
      DO ir = 1, ir_end
         !
         ! ... three dimensional indexes
         !
         i = index0 + ir - 1
         k = i / (dfftp%nr1x*dfftp%nr2x)
         i = i - (dfftp%nr1x*dfftp%nr2x)*k
         j = i / dfftp%nr1x
         i = i - dfftp%nr1x*j
         !
         DO ip = 1, 3
            r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                    DBLE( j )*inv_nr2*at(ip,2) + &
                    DBLE( k )*inv_nr3*at(ip,3)
         END DO
         !
         r(:) = pos(:) - r(:) 
         !
         ! ... minimum image convention
         !
         s(:) = MATMUL( r(:), bg(:,:) )
         s(:) = s(:) - ANINT(s(:))
         r(:) = MATMUL( at(:,:), s(:) )
         !
         dist = SUM( r * r ) 
         !
         gradrholocal( :, ir ) = scale * EXP(-dist/spr2) * r(:) * alat
         !      
      END DO
      !
      gradrho = gradrho + gradrholocal
      DEALLOCATE( gradrholocal )
      !
      RETURN
      !
!----------------------------------------------------------------------
      END SUBROUTINE generate_gradgaussian
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE generate_exponential( nnr, spread, pos, rho )
!----------------------------------------------------------------------
      !
      USE kinds,            ONLY : DP
      USE cell_base,        ONLY : at, bg, alat
      USE fft_base,         ONLY : dfftp
      USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)       :: nnr
      REAL( DP ), INTENT(IN)    :: spread
      REAL( DP ), INTENT(IN)    :: pos( 3 )
      REAL( DP ), INTENT(INOUT) :: rho( nnr )
      !
      ! ... Local variables
      !
      INTEGER                   :: i, j, k, ir, ir_end, ip
      INTEGER                   :: index0
      !
      REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
      REAL( DP )                :: dist, arg
      REAL( DP )                :: r( 3 ), s( 3 )
      REAL( DP ), ALLOCATABLE   :: rholocal ( : )
      REAL( DP ), PARAMETER     :: exp_arg_limit = 25.D0
      !
      inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
      inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
      inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
      !
      index0 = 0
      !
#if defined (__MPI)
      DO i = 1, me_bgrp
        index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
#endif
      !
#if defined (__MPI)
      ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
      ir_end = nnr
#endif
      !
      ALLOCATE( rholocal( nnr ) )
      rholocal = 0.D0
      !
      DO ir = 1, ir_end
         !
         ! ... three dimensional indexes
         !
         i = index0 + ir - 1
         k = i / (dfftp%nr1x*dfftp%nr2x)
         i = i - (dfftp%nr1x*dfftp%nr2x)*k
         j = i / dfftp%nr1x
         i = i - dfftp%nr1x*j
         r = 0.D0
         !
         DO ip = 1, 3
            r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                    DBLE( j )*inv_nr2*at(ip,2) + &
                    DBLE( k )*inv_nr3*at(ip,3)
         END DO
         !
         r(:) = pos(:) - r(:) 
         !
         ! ... minimum image convention
         !
         s(:) = MATMUL( r(:), bg(:,:) )
         s(:) = s(:) - ANINT(s(:))
         r(:) = MATMUL( at(:,:), s(:) )
         !
         dist = SQRT(SUM( r * r )) * alat
         arg = dist - spread
         !
         IF( ABS( arg ) .LT. exp_arg_limit ) THEN
           rholocal( ir ) = EXP( - arg ) 
         ELSE 
           rholocal( ir ) = 0.D0
         END IF
         !      
      END DO
      !
      rho = rho + rholocal
      DEALLOCATE( rholocal )
      !
      RETURN
      !
!----------------------------------------------------------------------
      END SUBROUTINE generate_exponential
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE generate_gradexponential( nnr, spread, pos, gradrho )
!----------------------------------------------------------------------
      !
      USE kinds,            ONLY : DP
      USE cell_base,        ONLY : at, bg, alat
      USE fft_base,         ONLY : dfftp
      USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)       :: nnr
      REAL( DP ), INTENT(IN)    :: spread
      REAL( DP ), INTENT(IN)    :: pos( 3 )
      REAL( DP ), INTENT(INOUT) :: gradrho( 3, nnr )
      !
      ! ... Local variables
      !
      INTEGER                   :: i, j, k, ir, ir_end, ip
      INTEGER                   :: index0
      !
      REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
      REAL( DP )                :: dist, arg
      REAL( DP )                :: r( 3 ), s( 3 )
      REAL( DP ), ALLOCATABLE   :: gradrholocal ( :, : )
      REAL( DP ), PARAMETER     :: exp_arg_limit = 25.D0
      !
      inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
      inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
      inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
      !
      index0 = 0
      !
#if defined (__MPI)
      DO i = 1, me_bgrp
        index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
#endif
      !
#if defined (__MPI)
      ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
      ir_end = nnr
#endif
      !
      ALLOCATE( gradrholocal( 3, nnr ) )
      gradrholocal = 0.D0
      !
      DO ir = 1, ir_end
         !
         ! ... three dimensional indexes
         !
         i = index0 + ir - 1
         k = i / (dfftp%nr1x*dfftp%nr2x)
         i = i - (dfftp%nr1x*dfftp%nr2x)*k
         j = i / dfftp%nr1x
         i = i - dfftp%nr1x*j
         !
         DO ip = 1, 3
            r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                    DBLE( j )*inv_nr2*at(ip,2) + &
                    DBLE( k )*inv_nr3*at(ip,3)
         END DO
         !
         r(:) = pos(:) - r(:) 
         !
         ! ... minimum image convention
         !
         s(:) = MATMUL( r(:), bg(:,:) )
         s(:) = s(:) - ANINT(s(:))
         r(:) = MATMUL( at(:,:), s(:) )
         !
         dist = SQRT(SUM( r * r )) * alat 
         arg = dist - spread
         IF ( dist .GT. 1.D-6 .AND. ABS( arg ) .LT. exp_arg_limit ) THEN
           gradrholocal( :, ir ) = r(:) * alat / dist * EXP( - arg ) 
         ELSE
           gradrholocal( :, ir ) = 0.D0
         ENDIF
         !      
      END DO
      !
      gradrho = gradrho + gradrholocal
      DEALLOCATE( gradrholocal )
      !
      RETURN
      !
!----------------------------------------------------------------------
      END SUBROUTINE generate_gradexponential
!----------------------------------------------------------------------
!----------------------------------------------------------------------
   SUBROUTINE generate_axis( nnr, icor, pos, axis )
!----------------------------------------------------------------------
   USE kinds,            ONLY : DP
   USE cell_base,        ONLY : at, bg, alat
   USE fft_base,         ONLY : dfftp
   USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
  !
  INTEGER, INTENT(IN) :: nnr
  INTEGER, INTENT(IN) :: icor
  REAL(DP), INTENT(IN) :: pos(3)
  REAL(DP), INTENT(OUT) :: axis( dfftp%nnr )
  !
  INTEGER  :: i, j, k, ir, ir_end, ip, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  REAL(DP) :: r(3), s(3)
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
  index0 = 0
  !
#if defined (__MPI)
  DO i = 1, me_bgrp
    index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
#endif
  !
#if defined (__MPI)
  ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  ir_end = nnr
#endif
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
     END DO
     !
     r(:) = r(:) - pos(:)  
     !
     ! ... minimum image convention
     !
     CALL cryst_to_cart( 1, r, bg, -1 )
     !
     r(:) = r(:) - ANINT( r(:) )
     !
     CALL cryst_to_cart( 1, r, at, 1 )
     !
     axis(ir) = r(icor)
     !
  END DO
  !
  axis = axis * alat
  !
  RETURN
  !
!----------------------------------------------------------------------
  END SUBROUTINE generate_axis 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
   SUBROUTINE generate_distance( nnr, pos, distance )
!----------------------------------------------------------------------
   USE kinds,            ONLY : DP
   USE cell_base,        ONLY : at, bg, alat
   USE fft_base,         ONLY : dfftp
   USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
  !
  INTEGER, INTENT(IN) :: nnr
  REAL(DP), INTENT(IN) :: pos(3)
  REAL(DP), INTENT(OUT) :: distance( 3, dfftp%nnr )
  !
  INTEGER  :: i, j, k, ir, ir_end, ip, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  REAL(DP) :: r(3), s(3)
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
  index0 = 0
  !
#if defined (__MPI)
  DO i = 1, me_bgrp
    index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
#endif
  !
#if defined (__MPI)
  ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  ir_end = nnr
#endif
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
     END DO
     !
     r(:) = r(:) - pos(:)  
     !
     ! ... minimum image convention
     !
     CALL cryst_to_cart( 1, r, bg, -1 )
     !
     r(:) = r(:) - ANINT( r(:) )
     !
     CALL cryst_to_cart( 1, r, at, 1 )
     !
     distance(:,ir) = r(:)
     !
  END DO
  !
  distance = distance * alat
  !
  RETURN
  !
!----------------------------------------------------------------------
  END SUBROUTINE generate_distance
!----------------------------------------------------------------------
!=----------------------------------------------------------------------=!
END MODULE generate_function
!=----------------------------------------------------------------------=!


