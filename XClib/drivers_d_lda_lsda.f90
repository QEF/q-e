!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dmxc_lda_l( length, rho_in, dmuxc )
  !---------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density.
  !
  USE dft_par_mod
  USE exch_lda_l,   ONLY: slater_l
  USE kind_l,       ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! the charge density ( positive )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: dmuxc
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: ex, vx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: arho, rhoaux, dr
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: ec, vc
  !
  REAL(DP) :: rho, rs, ex_s, vx_s
  REAL(DP) :: dpz_l
  INTEGER  :: iflg, ir, i1, i2, f1, f2
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP,        &
                         pi34 = 0.75_DP/3.141592653589793_DP,   &
                         third = 1.0_DP/3.0_DP, rho_trash = 0.5_DP, &
                         rs_trash = 1.0_DP
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  dmuxc = 0.0_DP
  !
  ! ... first case: analytical derivatives available
  !
  IF (iexch == 1 .AND. icorr == 1) THEN
  !
!$omp parallel if(ntids==1)
!$omp do private( rs, rho, ex_s, vx_s , iflg)
     DO ir = 1, length
        !
        rho = rho_in(ir)
        IF ( rho < -small ) rho = -rho_in(ir)
        !
        IF ( rho > small ) THEN
           rs = (pi34 / rho)**third
        ELSE
           dmuxc(ir) = 0.0_DP
           CYCLE
        ENDIF
        !
        CALL slater_l( rs, ex_s, vx_s )
        dmuxc(ir) = vx_s / (3.0_DP * rho)
        !
        iflg = 2
        IF (rs < 1.0_DP) iflg = 1
        dmuxc(ir) = dmuxc(ir) + dpz_l( rs, iflg )
        dmuxc(ir) = dmuxc(ir) * SIGN(1.0_DP,rho_in(ir))
        !
     ENDDO
!$omp end do
!$omp end parallel
     !
  ELSE
     !
     ! ... second case: numerical derivatives
     !
     ALLOCATE( ex(2*length), vx(2*length)  )
     ALLOCATE( ec(2*length), vc(2*length)  )
     ALLOCATE( arho(length), dr(length), rhoaux(2*length) )
     !
     i1 = 1         ;  f1 = length             !two blocks:  [ rho+dr ]
     i2 = length+1  ;  f2 = 2*length           !             [ rho-dr ]              
     !
     arho = ABS(rho_in)
     dr = 0.0_DP
     WHERE ( arho > small ) dr = MIN( 1.E-6_DP, 1.E-4_DP * rho_in )
     !
     rhoaux(i1:f1) = arho+dr
     rhoaux(i2:f2) = arho-dr
     !
     CALL xc_lda_l( length*2, rhoaux, ex, ec, vx, vc )
     !
     WHERE ( arho < small ) dr = 1.0_DP ! ... to avoid NaN in the next operation
     !
     dmuxc(:) = (vx(i1:f1) + vc(i1:f1) - vx(i2:f2) - vc(i2:f2)) / &
                (2.0_DP * dr(:))
     !
     DEALLOCATE( ex, vx  )
     DEALLOCATE( ec, vc  )
     DEALLOCATE( dr, rhoaux )
     !
     WHERE ( arho < small ) dmuxc = 0.0_DP
     ! however a higher threshold is already present in xc_lda()
     dmuxc(:) = dmuxc(:) * SIGN(1.0_DP,rho_in(:))
     !
     DEALLOCATE( arho )
     !
  ENDIF
  !
  ! bring to rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  RETURN
  !
END SUBROUTINE dmxc_lda_l
!
!
!-----------------------------------------------------------------------
FUNCTION dpz_l( rs, iflg )
  !-----------------------------------------------------------------------
  !!  Derivative of the correlation potential with respect to local density
  !!  Perdew and Zunger parameterization of the Ceperley-Alder functional.
  !
  USE kind_l,      ONLY: DP
  USE constants_l, ONLY: pi, fpi
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  INTEGER,  INTENT(IN) :: iflg
  REAL(DP) :: dpz_l
  !
  !  ... local variables
  !  a,b,c,d,gc,b1,b2 are the parameters defining the functional
  !
  REAL(DP), PARAMETER :: a = 0.0311d0, b = -0.048d0, c = 0.0020d0, &
       d = -0.0116d0, gc = -0.1423d0, b1 = 1.0529d0, b2 = 0.3334d0,&
       a1 = 7.0d0 * b1 / 6.d0, a2 = 4.d0 * b2 / 3.d0
  REAL(DP) :: x, den, dmx, dmrs
  !
  IF (iflg == 1) THEN
     dmrs = a / rs + 2.d0 / 3.d0 * c * (LOG(rs) + 1.d0) + &
          (2.d0 * d-c) / 3.d0
  ELSE
     x = SQRT(rs)
     den = 1.d0 + x * (b1 + x * b2)
     dmx = gc * ( (a1 + 2.d0 * a2 * x) * den - 2.d0 * (b1 + 2.d0 * &
           b2 * x) * (1.d0 + x * (a1 + x * a2) ) ) / den**3
     dmrs = 0.5d0 * dmx / x
  ENDIF
  !
  dpz_l = - fpi * rs**4.d0 / 9.d0 * dmrs
  !
  RETURN
  !
END FUNCTION dpz_l
!
