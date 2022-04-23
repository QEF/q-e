!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!-----------------------------------------------------------------------
SUBROUTINE calculate_phase (gvect, phase)
  !-----------------------------------------------------------------------
  !
  !! This routine computes phase(r) = exp(-i gvect*r)
  !! gvect is a G vector from input. In output phase(r)
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE fft_support,          ONLY : good_fft_dimension
  USE cell_base,            ONLY : at
  USE constants,            ONLY : tpi
  !
  IMPLICIT NONE
  !
  !REAL(DP), INTENT(IN) :: gvect(3)
  !
  REAL(DP) :: gvect(3)
  ! ... the G vector in input 
  !
  COMPLEX(DP), INTENT(OUT) :: phase(dffts%nnr)
  ! ... e^{i G \dot r } in OUTPUT
  !
  REAL(DP),    ALLOCATABLE :: r(:,:)
  ! ... position 
  !
  INTEGER  :: i, j, k, ip, ir, ir_end, idx, j0, k0
  ! ... counters 
  !
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3, dot_prod
  ! ... invers of the real grid dimensions
  ! 
  ALLOCATE(r(dffts%nnr,3))
  r(:,:) = 0.d0
  !
  ! The reciprocal vector in crystall coordinate
  !CALL cryst_to_cart(1, gvect, at, -1)
  !
  !gvect (1) = gvect(1)*bg(1,1) + gvect(2)*bg(1,2) + gvect(3)*bg(1,3)
  !gvect (2) = gvect(1)*bg(2,1) + gvect(2)*bg(2,2) + gvect(3)*bg(2,3)
  !gvect (3) = gvect(1)*bg(3,1) + gvect(2)*bg(3,2) + gvect(3)*bg(3,3)
  !WRITE(*,'("NICOLA gvect=", 3f8.4)') gvect
  !WRITE(*,'(3f8.4)') sum(at(:,1)*bg(:,1)), sum(at(:,1)*bg(:,2)), sum(at(:,1)*bg(:,3))
  !WRITE(*,'(3f8.4)') sum(at(:,2)*bg(:,1)), sum(at(:,2)*bg(:,2)), sum(at(:,2)*bg(:,3))
  !WRITE(*,'(3f8.4)') sum(at(:,3)*bg(:,1)), sum(at(:,3)*bg(:,2)), sum(at(:,3)*bg(:,3))
  !
  ! Calculate r 
  !
  inv_nr1 = 1.D0 / DBLE( dffts%nr1 )
  inv_nr2 = 1.D0 / DBLE( dffts%nr2 )
  inv_nr3 = 1.D0 / DBLE( dffts%nr3 )
  !
#if defined (__MPI)
  j0 = dffts%my_i0r2p ; k0 = dffts%my_i0r3p
  IF(dffts%nr1x == 0 ) dffts%nr1x = good_fft_dimension( dffts%nr1 )
  ir_end = MIN(dffts%nnr,dffts%nr1x*dffts%my_nr2p*dffts%my_nr3p)
#else
  j0 = 0 ; k0 = 0
  ir_end = dffts%nnr
#endif
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     idx = ir -1
     k   = idx / (dffts%nr1x*dffts%my_nr2p)
     idx = idx - (dffts%nr1x*dffts%my_nr2p)*k
     k   = k + k0
     IF ( k .GE. dffts%nr3 ) CYCLE
     j   = idx / dffts%nr1x
     idx = idx - dffts%nr1x * j
     j   = j + j0
     IF ( j .GE. dffts%nr2 ) CYCLE
     i   = idx
     IF ( i .GE. dffts%nr1 ) CYCLE
     !
     DO ip = 1, 3
        r(ir,ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                   DBLE( j )*inv_nr2*at(ip,2) + &
                   DBLE( k )*inv_nr3*at(ip,3)
     ENDDO
     !
     dot_prod = tpi*sum(r(ir,:)*gvect(:))
     phase(ir) = CMPLX( COS(dot_prod), -1.d0*SIN(dot_prod),KIND=dp)

     !
  ENDDO
  DEALLOCATE(r)
  !
END subroutine
