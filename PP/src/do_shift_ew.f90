!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE do_shift_ew (alat, nat, ntyp, ityp, zv, delta_zv, at, bg, tau, &
     omega, g, gg, ngm, gcutm, gstart, gamma_only, shift_ion)
  !-----------------------------------------------------------------------
  !
  ! Calculates Ewald energy with both G- and R-space terms.
  ! Determines optimal alpha. Should hopefully work for any structure.
  !
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi, e2
  USE mp_global, ONLY : intra_pool_comm
  USE mp, ONLY : mp_sum
  IMPLICIT NONE
  !
  !   first the dummy variables
  !

  INTEGER :: nat, ntyp, ityp (nat), ngm, gstart
  ! input: number of atoms in the unit cell
  ! input: number of different types of atoms
  ! input: the type of each atom
  ! input: number of plane waves for G sum
  ! input: first non-zero G vector

  LOGICAL :: gamma_only

  real(DP) :: tau (3, nat), g (3, ngm), gg (ngm), zv (ntyp), &
       at (3, 3), bg (3, 3), omega, alat, gcutm, delta_zv(ntyp),  &
       shift_ion(nat)
  ! input: the positions of the atoms in the cell
  ! input: the coordinates of G vectors
  ! input: the square moduli of G vectors
  ! input: the charge of each type of atoms
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! input: the volume of the unit cell
  ! input: lattice parameter
  ! input: cut-off of g vectors
  real(DP) :: ewald
  ! output: the ewald energy
  !
  !    here the local variables
  !
  INTEGER, PARAMETER :: mxr = 50
  ! the maximum number of R vectors included in r
  INTEGER :: ng, nr, na, nb, nt, nrm, ipol
  ! counter over reciprocal G vectors
  ! counter over direct vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atomic types
  ! number of R vectors included in r sum
  ! counter on polarization

  real(DP) :: charge, tpiba2, ewaldg, ewaldr, dtau (3), alpha, &
       r (3, mxr), r2 (mxr), rmax, rr, upperbound, fact, arg
  ! total ionic charge in the cell
  ! length in reciprocal space
  ! ewald energy computed in reciprocal space
  ! ewald energy computed in real space
  ! the difference tau_s - tau_s'
  ! alpha term in ewald sum
  ! input of the rgen routine ( not used here )
  ! the square modulus of R_j-tau_s-tau_s'
  ! the maximum radius to consider real space sum
  ! buffer variable
  ! used to optimize alpha
  COMPLEX(DP), ALLOCATABLE :: rhon(:)
  real(DP), EXTERNAL :: qe_erfc

  ALLOCATE (rhon(ngm))

  shift_ion(:) = 0.d0

  tpiba2 = (tpi / alat) **2
  charge = 0.d0
  DO na = 1, nat
     charge = charge+zv (ityp (na) )
  ENDDO
  alpha = 2.9d0
100 alpha = alpha - 0.1d0
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  !
  IF (alpha<=0.d0) CALL errore ('do_shift_ew', 'optimal alpha not found', 1)
  upperbound = 2.d0 * charge**2 * sqrt (2.d0 * alpha / tpi) * qe_erfc ( &
       sqrt (tpiba2 * gcutm / 4.d0 / alpha) )
  IF (upperbound>1.0d-7) GOTO 100
  !
  ! G-space sum here.
  ! Determine if this processor contains G=0 and set the constant term
  !
  IF (gstart==2) THEN
     DO na =1,nat
        shift_ion(na) = - charge * delta_zv(ityp(na)) /alpha/ 4.0d0
     ENDDO
  ENDIF
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  DO ng = gstart, ngm
     rhon(ng) = (0.d0, 0.d0)
     DO na =1, nat
        arg = (g (1, ng) * tau (1, na) + &
               g (2, ng) * tau (2, na) + &
               g (3, ng) * tau (3, na) ) * tpi
        rhon(ng) = rhon(ng) + zv (ityp(na)) * cmplx(cos (arg), -sin (arg),kind=DP)
     ENDDO
  ENDDO
  DO na=1,nat
     DO ng=gstart, ngm
        arg = (g (1, ng) * tau (1, na) + g (2, ng) * tau (2, na) &
             + g (3, ng) * tau (3, na) ) * tpi
        shift_ion(na) = shift_ion(na) + fact * delta_zv(ityp(na)) * &
                        conjg(rhon(ng)) * cmplx(cos (arg), -sin (arg),kind=DP) * &
                        exp ( -gg(ng)*tpiba2/alpha/4.d0) / gg(ng)/tpiba2
     ENDDO
  ENDDO
  shift_ion(:) = 2.d0 * tpi / omega * shift_ion(:)
  !
  !  Here add the other constant term
  !
  IF (gstart==2) THEN
     DO na = 1, nat
        shift_ion(na) = shift_ion(na) - &
                        zv (ityp (na) ) * delta_zv(ityp(na)) *  &
                        sqrt (8.d0/tpi*alpha)
     ENDDO
  ENDIF
  !
  ! R-space sum here (only for the processor that contains G=0)
  !
  IF (gstart==2) THEN
     rmax = 4.d0 / sqrt (alpha) / alat
     !
     ! with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
     !
     DO na = 1, nat
        DO nb = 1, nat
           DO ipol = 1, 3
              dtau (ipol) = tau (ipol, na) - tau (ipol, nb)
           ENDDO
           !
           ! generates nearest-neighbors shells
           !
           CALL rgen (dtau, rmax, mxr, at, bg, r, r2, nrm)
           !
           ! and sum to the real space part
           !
           DO nr = 1, nrm
              rr = sqrt (r2 (nr) ) * alat
              shift_ion(na) = shift_ion(na) + &
                       delta_zv(ityp(na)) * zv (ityp (nb) ) *  &
                       qe_erfc ( sqrt (alpha) * rr) / rr
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  shift_ion(:) = e2 * shift_ion(:)

  CALL mp_sum ( shift_ion, intra_pool_comm )

  DEALLOCATE (rhon)
  RETURN
END SUBROUTINE do_shift_ew

