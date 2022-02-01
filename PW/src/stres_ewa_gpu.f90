!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE stres_ewa_gpu( alat, nat, ntyp, ityp, zv, at, bg, tau,    &
                          omega, g_d, gg_d, ngm, gstart, gamma_only, &
                          gcutm, sigmaewa )
  !---------------------------------------------------------------------
  !! Ewald contribution. Both real- and reciprocal-space terms are 
  !! present.
  !
  USE kinds
  USE constants,           ONLY : tpi, e2, eps6
  USE mp_bands,            ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp,                  ONLY : mp_sum
  USE Coul_cut_2D,         ONLY : do_cutoff_2D, cutoff_stres_sigmaewa_gpu
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
  REAL(DP), INTENT(IN) :: g_d(3,ngm)
  !! input: the coordinates of G vectors
  REAL(DP), INTENT(IN) :: gg_d(ngm)
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
  ! output: the ewald stress
  !
  ! ... local variables
  !
  INTEGER, PARAMETER :: mxr = 50
  ! the maximum number of R vectors included in r sum
  INTEGER :: ng,  nr, na, nb, l, m, nrm
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
  !
  INTEGER :: ierr(2)
  REAL(DP) :: sigma11, sigma21, sigma22, sigma31, sigma32, sigma33
  COMPLEX(DP) :: rhostar
  !
  INTEGER,  ALLOCATABLE :: ityp_d(:)
  REAL(DP), ALLOCATABLE :: zv_d(:), tau_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: g_d, gg_d, zv_d, ityp_d, tau_d
#endif  
  !
  tpiba2 = (tpi / alat)**2
  sigmaewa(:,:) = 0._DP
  charge = 0._DP
  !
  ALLOCATE( zv_d(ntyp), tau_d(3,nat) )
  zv_d  = zv
  tau_d = tau
  ALLOCATE( ityp_d(nat) )
  ityp_d = ityp
  !
  DO na = 1, nat
     charge = charge + zv(ityp(na))
  ENDDO
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error ON THE ENERGY
  !
  alpha = 2.9_DP
12 alpha = alpha - 0.1_DP
  !
  IF (alpha==0.0) CALL errore( 'stres_ew', 'optimal alpha not found', 1 )
  upperbound = e2 * charge**2 * SQRT(2 * alpha / tpi) * &
               erfc ( SQRT(tpiba2 * gcutm / 4._DP / alpha) )
  !
  IF (upperbound > 1d-7) GOTO 12
  !
  ! G-space sum here
  !
  ! Determine if this processor contains G=0 and set the constant term
  !
  IF (gstart == 2) THEN
     sdewald = tpi * e2 / 4._DP / alpha * (charge / omega)**2
  ELSE
     sdewald = 0._DP
  ENDIF
  !
  ! sdewald is the diagonal term
  IF ( gamma_only ) THEN
     fact = 2._DP
  ELSE
     fact = 1._DP
  ENDIF
  !
  IF ( do_cutoff_2D ) THEN 
     !
     CALL cutoff_stres_sigmaewa_gpu( alpha, sdewald, sigmaewa )
     !
  ELSE
     !
     sigma11 = 0._DP ; sigma21 = 0._DP ; sigma22 = 0._DP
     sigma31 = 0._DP ; sigma32 = 0._DP ; sigma33 = 0._DP
     !
     !$cuf kernel do (1) <<<*,*>>>
     DO ng = gstart, ngm
        g2 = gg_d(ng) * tpiba2
        g2a = g2 / 4._DP / alpha
        rhostar = (0._DP,0._DP)
        DO na = 1, nat
           arg = (g_d(1,ng) * tau_d(1,na) + g_d(2,ng) * tau_d(2,na) + &
                  g_d(3,ng) * tau_d(3,na) ) * tpi
           rhostar = rhostar + CMPLX(zv_d(ityp_d(na))) * CMPLX(COS(arg), SIN(arg), KIND=DP)
        ENDDO
        rhostar = rhostar / CMPLX(omega)
        sewald = fact * tpi * e2 * EXP(-g2a) / g2 * ABS(rhostar)**2
        sdewald = sdewald - sewald
        !
        sigma11 = sigma11 + sewald * tpiba2 * 2._DP * &
                              g_d(1,ng) * g_d(1,ng) / g2 * (g2a + 1)
        sigma21 = sigma21 + sewald * tpiba2 * 2._DP * &
                              g_d(2,ng) * g_d(1,ng) / g2 * (g2a + 1)
        sigma22 = sigma22 + sewald * tpiba2 * 2._DP * &
                              g_d(2,ng) * g_d(2,ng) / g2 * (g2a + 1)
        sigma31 = sigma31 + sewald * tpiba2 * 2._DP * &
                              g_d(3,ng) * g_d(1,ng) / g2 * (g2a + 1)
        sigma32 = sigma32 + sewald * tpiba2 * 2._DP * &
                              g_d(3,ng) * g_d(2,ng) / g2 * (g2a + 1)
        sigma33 = sigma33 + sewald * tpiba2 * 2._DP * &
                              g_d(3,ng) * g_d(3,ng) / g2 * (g2a + 1)
        !
     ENDDO
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
  DO l = 1, 3
     sigmaewa(l,l) = sigmaewa(l,l) + sdewald
  ENDDO
  !
  ! R-space sum here (see ewald.f90 for details on parallelization)
  !
  CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
  !
  IF ( mykey == 0 ) THEN
     rmax = 4.0d0 / SQRT(alpha) / alat
     !
     ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
     !
     DO na = na_s, na_e
        DO nb = 1, nat
           dtau(:) = tau(:,na) - tau(:,nb)
           !
           !     generates nearest-neighbors shells r(i)=R(i)-dtau(i)
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
  DEALLOCATE( zv_d, tau_d )
  DEALLOCATE( ityp_d )
  !
  CALL mp_sum( sigmaewa, intra_bgrp_comm )
  !
  RETURN
  !
END SUBROUTINE stres_ewa_gpu

