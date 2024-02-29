!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE stres_ewa( alat, nat, ntyp, ityp, zv, at, bg, tau, &
                      omega, g, gg, ngm, gstart, gamma_only,  &
                      gcutm, sigmaewa )
  !---------------------------------------------------------------------
  !! Ewald contribution. Both real- and reciprocal-space terms are 
  !! present.
  !
  USE kinds
  USE constants,    ONLY : tpi, e2, eps6
  USE mp_bands,     ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp,           ONLY : mp_sum
  USE Coul_cut_2D,  ONLY : do_cutoff_2D, cutoff_stres_sigmaewa
  !
  IMPLICIT NONE
  !
  INTEGER :: nat
  !! input: number of atoms in the unit cell
  INTEGER :: ntyp
  !! input: number of different types of atoms
  INTEGER :: ityp(nat)
  !! input: the type of each atom
  INTEGER :: ngm
  !! input: number of plane waves for G sum
  INTEGER :: gstart
  !! input: first nonzero g vector
  LOGICAL, INTENT(IN) :: gamma_only
  !! gamma point only
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! input: the positions of the atoms in the cell
  REAL(DP), INTENT(IN) :: g(3,ngm)
  !! input: the coordinates of G vectors
  REAL(DP), INTENT(IN) :: gg(ngm)
  !! input: the square moduli of G vectors
  REAL(DP), INTENT(IN) :: zv(ntyp)
  !! input: the charge of each type of atoms
  REAL(DP), INTENT(IN) :: at(3,3)
  !! input: the direct lattice vectors
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! input: the reciprocal lattice vectors
  REAL(DP), INTENT(IN) :: omega
  !! input: the volume of the unit cell
  REAL(DP), INTENT(IN) :: alat
  !! input: measure of length
  REAL(DP), INTENT(IN) :: gcutm
  !! input: cut-off of g vectors
  REAL(DP), INTENT(OUT) :: sigmaewa(3,3)
  !! output: the ewald stress
  !
  ! ... local variables
  !
  INTEGER, PARAMETER :: mxr = 50
  ! the maximum number of R vectors included in r sum
  integer :: ng,  nr, na, nb, l, m, nrm
  ! counter over reciprocal G vectors
  ! counter over direct vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atoms
  ! number of R vectors included in r sum
  INTEGER :: na_s, na_e, mykey
  !
  REAL(DP) :: charge, arg, tpiba2, dtau(3), alpha, r(3,mxr),     &
              r2(mxr), rmax, rr, upperbound, fact, fac, g2, g2a, &
              sdewald, sewald
  ! total ionic charge in the cell
  ! the argument of the phase
  ! length in reciprocal space
  ! the difference tau_s - tau_s'
  ! alpha term in ewald sum
  ! input of the rgen routine ( not used here )
  ! the square modulus of R_j-tau_s-tau_s'
  ! the maximum radius to consider real space sum
  ! buffer variable
  ! used to optimize alpha
  ! auxiliary variables
  ! diagonal term
  ! nondiagonal term
  COMPLEX(DP) :: rhostar
  REAL(DP) :: sigma11, sigma21, sigma22, sigma31, sigma32, sigma33
  !
  !$acc data present( g, gg )
  !
  tpiba2 = (tpi / alat)**2
  sigmaewa(:,:) = 0.d0
  charge = 0.d0
  !
  DO na = 1, nat
     charge = charge + zv(ityp(na))
  ENDDO
  !
  ! ... choose alpha in order to have convergence in the sum over G
  !     upperbound is a safe upper bound for the error ON THE ENERGY
  !
  alpha = 2.9d0
12 alpha = alpha - 0.1d0
  !
  IF (alpha==0.0) CALL errore( 'stres_ew', 'optimal alpha not found', 1 )
  upperbound = e2 * charge**2 * SQRT(2 * alpha / tpi) * &
               ERFC( SQRT(tpiba2 * gcutm / 4.0d0 / alpha) )
  !
  IF (upperbound > 1d-7) GOTO 12
  !
  ! ... Determine if this processor contains G=0 and set the constant term
  !     sdewald is the diagonal term
  IF (gstart == 2) THEN
     sdewald = tpi * e2 / 4.d0 / alpha * (charge / omega)**2
  ELSE
     sdewald = 0.d0
  ENDIF
  !
  IF (gamma_only) THEN
    fact = 2.d0
  ELSE
    fact = 1.d0
  ENDIF
  !
  ! ... G-space sum here below
  !
  IF (do_cutoff_2D) THEN 
     !
     CALL cutoff_stres_sigmaewa( gamma_only, alpha, sdewald, sigmaewa )
     !
  ELSE
     !
     sigma11 = 0._DP ; sigma21 = 0._DP ; sigma22 = 0._DP
     sigma31 = 0._DP ; sigma32 = 0._DP ; sigma33 = 0._DP
     !
#if !defined(_OPENACC)
!$omp parallel do default(none) shared(gstart, ngm, g, gg, tpiba2, alpha, tau,&
!$omp   nat, zv, ityp, omega, fact) private(g2,g2a, rhostar, na, arg,&
!$omp   sewald) reduction(+:sdewald,sigma11,sigma21,sigma22,sigma31,&
!$omp   sigma32,sigma33)
#else
!$acc parallel loop copyin(tau,zv,ityp) reduction(+:sigma11,sigma21,sigma22,&
!$acc                                               sigma31,sigma32,sigma33,sdewald)
#endif
     DO ng = gstart, ngm
        g2 = gg(ng) * tpiba2
        g2a = g2 / 4._DP / alpha
        rhostar = (0._DP,0._DP)
        DO na = 1, nat
           arg = (g(1,ng) * tau(1,na) + g(2,ng) * tau(2,na) + &
                  g(3,ng) * tau(3,na) ) * tpi
           rhostar = rhostar + CMPLX(zv(ityp(na)), KIND=DP) * CMPLX(COS(arg), SIN(arg), KIND=DP)
        ENDDO
        rhostar = rhostar / CMPLX(omega, KIND=DP)
        sewald = fact * tpi * e2 * EXP(-g2a) / g2 * ABS(rhostar)**2
        sdewald = sdewald - sewald
        !
        sigma11 = sigma11 + sewald * tpiba2 * 2._DP * &
                              g(1,ng) * g(1,ng) / g2 * (g2a + 1)
        sigma21 = sigma21 + sewald * tpiba2 * 2._DP * &
                              g(2,ng) * g(1,ng) / g2 * (g2a + 1)
        sigma22 = sigma22 + sewald * tpiba2 * 2._DP * &
                              g(2,ng) * g(2,ng) / g2 * (g2a + 1)
        sigma31 = sigma31 + sewald * tpiba2 * 2._DP * &
                              g(3,ng) * g(1,ng) / g2 * (g2a + 1)
        sigma32 = sigma32 + sewald * tpiba2 * 2._DP * &
                              g(3,ng) * g(2,ng) / g2 * (g2a + 1)
        sigma33 = sigma33 + sewald * tpiba2 * 2._DP * &
                              g(3,ng) * g(3,ng) / g2 * (g2a + 1)
     ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
     !
     sigmaewa(1,1) = sigmaewa(1,1) + sigma11
     sigmaewa(2,1) = sigmaewa(2,1) + sigma21
     sigmaewa(2,2) = sigmaewa(2,2) + sigma22
     sigmaewa(3,1) = sigmaewa(3,1) + sigma31
     sigmaewa(3,2) = sigmaewa(3,2) + sigma32
     sigmaewa(3,3) = sigmaewa(3,3) + sigma33
     !
  ENDIF
  !
  !$acc end data
  !
  DO l = 1, 3
     sigmaewa(l,l) = sigmaewa(l,l) + sdewald
  ENDDO
  !
  ! ... R-space sum here (see ewald.f90 for details on parallelization)
  !
  CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
  !
  IF ( mykey == 0 ) THEN
     rmax = 4.0d0 / SQRT(alpha) / alat
     !
     ! ... with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
     !
!$omp parallel do default(none) shared(na_s,na_e,nat,tau,rmax,at,bg,alat,ityp,alpha,omega,zv)&
!$omp                          &private(nb,dtau,r,r2,nrm,nr,rr,fac,l,m) reduction(+:sigmaewa)
     DO na = na_s, na_e
        DO nb = 1, nat
           dtau(:) = tau(:,na) - tau(:,nb)
           !
           ! ... generates nearest-neighbors shells r(i)=R(i)-dtau(i)
           !
           CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
           !
           DO nr = 1, nrm
              rr = SQRT(r2 (nr) ) * alat
              fac = - e2 / 2.0_DP/ omega * alat**2 * zv(ityp(na)) * &
                    zv(ityp(nb)) / rr**3 * (erfc(SQRT(alpha) * rr) + &
                    rr * SQRT(8.0_DP * alpha / tpi) * EXP( - alpha * rr**2) )
              DO l = 1, 3
                 DO m = 1, l
                    sigmaewa(l,m) = sigmaewa(l,m) + fac * r(l,nr) * r(m,nr)
                 ENDDO
              ENDDO
           ENDDO
           !
        ENDDO
     ENDDO
!$omp end parallel do
  ENDIF
  !
  DO l = 1, 3
     DO m = 1, l - 1
        sigmaewa(m,l) = sigmaewa(l,m)
     ENDDO
  ENDDO
  !
  DO l = 1, 3
     DO m = 1, 3
        sigmaewa(l,m) = - sigmaewa(l,m)
     ENDDO
  ENDDO
  !
  CALL mp_sum( sigmaewa, intra_bgrp_comm )
  !
  RETURN
  !
END SUBROUTINE stres_ewa

