!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine d2ionq (nat, ntyp, ityp, zv, tau, alat, omega, q, at, &
     bg, g, gg, ngm, gcutm, nmodes, u, dyn)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the contribution of the ions to the
  !    dynamical matrix. Both the real and reciprocal space terms
  !    are included.
  !
  !    The original routine was from C. Bungaro.
  !    Revised 16 oct. 1995 by Andrea Dal Corso.
  !    April 1997: parallel stuff added (SdG)
  !
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  !  first the dummy variables
  !

  integer :: nat, ntyp, ngm, ityp (nat), nmodes
  ! input: the number of atoms
  ! input: the number of types of atoms
  ! input: the number of G vectors
  ! input: the type of each atom
  ! input: the number of modes

  real(kind=DP) :: tau (3, nat), g (3, ngm), gg (ngm), zv (ntyp), &
       at (3, 3), bg (3, 3), omega, alat, gcutm, q (3)
  ! input: the positions of the atoms
  ! input: the coordinates of g vectors
  ! input: the modulus of g vectors
  ! input: the charge of each type
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! input: the volume of the unit cell
  ! input: the length scale
  ! input: cut-off of g vectors
  ! input: the q vector
  complex(kind=DP) :: dyn (3 * nat, nmodes), u (3 * nat, nmodes)
  ! output: the ionic part of the dyn. mat
  ! input: the pattern of the modes
  !
  !   three parameters
  !

  integer :: mxr
  ! the maximum number of r shells

  real(kind=DP) :: e2, tpi
  ! the electron charge
  ! twp times pi

  parameter (mxr = 100, e2 = 2.d0, tpi = 2.d0 * 3.14159265358979d0)

  integer :: nu_i, nu_j, na, nb, nta, ntb, ng, nrm, nr, icart, &
       jcart, na_icart, na_jcart, nb_icart, nb_jcart
  ! counter on the modes
  ! counter on the modes
  ! counter on atoms
  ! counter on atoms
  ! the type of atom na
  ! the type of atom nb
  ! counter on G vectors
  ! the real number of R shells
  ! counter on r shells
  ! counter on cartesian coordinates
  ! counter on cartesian coordinates
  ! counter on the modes
  ! counter on the modes
  ! counter on the modes
  ! counter on the modes
  real(kind=DP) :: arg, fpi, argq, tpiba2, alpha, r (3, mxr), r2 (mxr), &
       dtau (3), erfc, rmax, rr, upperbound, charge, fac, df, d2f, ar, &
       gtq2, gt2, facq, qrg
  ! the argument of the phase
  ! 4 times pi
  ! the argument of the phase
  ! two pi / alat
  ! the parameter alpha used to converge the summ
  ! the R_i-\tau_s-\tau_s' parameter
  ! the square of this parameter
  ! delta tau
  ! the erfc function
  ! the maximum r on real space summations
  ! the rmax
  ! an estimate of the error in the real term
  ! the total charge
  ! the phase factor
  ! auxiliary variable
  ! auxiliary variable
  !
  ! the modulus of |q+G|
  ! the modulus of G
  ! a phase factor
  !


  complex(kind=DP) :: dy1 (3 * nat, nmodes), dy2 (3 * nat, nmodes), &
       dy3 (3 * nat, nmodes), facg, fnat, work
  ! the first term dynamical matrix
  ! the second dynamical matrix
  ! the third term of the dynamical matrix
  ! a phase with g
  ! a phase with the atoms
  ! working space

  call start_clock ('d2ionq')

  tpiba2 = (tpi / alat) **2

  fpi = 2.d0 * tpi
  charge = 0.d0
  do na = 1, nat
     charge = charge+zv (ityp (na) )
  enddo
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is an upper bound for the error in the sum over G
  ! estimated for the energy (empirical trust!)
  !
  alpha = 2.9d0

11 alpha = alpha - 0.1d0

  if (alpha.eq.0.d0) call errore ('d2ion', 'optimal alpha not found', &
       1)

  upperbound = 2.d0 * charge**2 * sqrt (2.d0 * alpha / tpi) * erfc ( &
       sqrt (tpiba2 * gcutm / 4.d0 / alpha) )

  if (upperbound.gt.1.d-9) goto 11

  WRITE( stdout, '(/5x,"Alpha used in Ewald sum = ",f8.4)') alpha
  call setv (2 * 3 * nat * nmodes, 0.d0, dy1, 1)
  call setv (2 * 3 * nat * nmodes, 0.d0, dy2, 1)
  call setv (2 * 3 * nat * nmodes, 0.d0, dy3, 1)
  !
  ! G-space sums here
  !
  do ng = 1, ngm
     !
     !     The reciprocal space term has two contributions. The first one
     !
     gtq2 = ( (g (1, ng) + q (1) ) **2 + (g (2, ng) + q (2) ) **2 + &
          (g (3, ng) + q (3) ) **2) * tpiba2
     if (abs (gtq2) .gt.1.d-8) then
        facq = - e2 * fpi * tpiba2 / omega * exp ( - gtq2 / alpha / &
             4.d0) / gtq2
     else
        facq = 0.d0

     endif
     do na = 1, nat
        nta = ityp (na)
        do nb = 1, nat
           ntb = ityp (nb)
           argq = tpi * ( (g (1, ng) + q (1) ) * (tau (1, na) - tau (1, nb) ) &
                + (g (2, ng) + q (2) ) * (tau (2, na) - tau (2, nb) ) + (g (3, ng) &
                + q (3) ) * (tau (3, na) - tau (3, nb) ) )
           facg = facq * zv (nta) * zv (ntb) * DCMPLX (cos (argq), sin (argq) &
                )
           do icart = 1, 3
              nu_i = 3 * (na - 1) + icart
              do jcart = 1, 3
                 nu_j = 3 * (nb - 1) + jcart
                 dy1 (nu_i, nu_j) = dy1 (nu_i, nu_j) + facg * (q (icart) + g ( &
                      icart, ng) ) * (q (jcart) + g (jcart, ng) )
              enddo
           enddo
        enddo
     enddo
     !
     !   the second term
     !

     gt2 = gg (ng) * tpiba2
     if (abs (gt2) .gt.1.d-8) then
        fac = - e2 * fpi * tpiba2 / omega * exp ( - gt2 / alpha / 4.d0) &
             / gt2
     else
        fac = 0.d0

     endif
     do na = 1, nat
        nta = ityp (na)
        fnat = (0.d0, 0.d0)
        do nb = 1, nat
           ntb = ityp (nb)
           arg = tpi * ( (g (1, ng) ) * (tau (1, na) - tau (1, nb) ) + &
                (g (2, ng) ) * (tau (2, na) - tau (2, nb) ) + (g (3, ng) ) &
                * (tau (3, na) - tau (3, nb) ) )
           facg = fac * zv (nta) * zv (ntb) * DCMPLX (cos (arg), 0.d0)
           fnat = fnat + facg

        enddo
        do icart = 1, 3
           nu_i = 3 * (na - 1) + icart
           do jcart = 1, 3
              nu_j = 3 * (na - 1) + jcart
              dy2 (nu_i, nu_j) = dy2 (nu_i, nu_j) + fnat * g (icart, ng) &
                   * g (jcart, ng)
           enddo
        enddo
     enddo

  enddo
  do nu_i = 1, nmodes
     do nu_j = 1, nmodes
        dy3 (nu_i, nu_j) = dy3 (nu_i, nu_j) + dy1 (nu_i, nu_j) - dy2 ( &
             nu_i, nu_j)
     enddo
  enddo
  !
  !  Then there is also a part in real space which is computed here.
#ifdef __PARA
  !   ... only by the node that contains  G=0
  !
  if (gg (1) .gt.1.d-8) goto 100
#endif
  !
  rmax = 5.d0 / sqrt (alpha) / alat
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
  !
  do na = 1, nat
     nta = ityp (na)
     do nb = 1, nat
        ntb = ityp (nb)
        do icart = 1, 3
           dtau (icart) = tau (icart, na) - tau (icart, nb)
        enddo
        !
        ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
        !
        call rgen (dtau, rmax, mxr, at, bg, r, r2, nrm)
        do nr = 1, nrm
           rr = sqrt (r2 (nr) ) * alat
           ar = sqrt (alpha) * rr
           qrg = tpi * (q (1) * (r (1, nr) + dtau (1) ) + q (2) * (r (2, nr) &
                + dtau (2) ) + q (3) * (r (3, nr) + dtau (3) ) )
           d2f = (3.d0 * erfc (ar) + sqrt (8.d0 / tpi) * ar * (3.d0 + 2.d0 * &
                ar**2) * exp ( - ar**2) ) / rr**5
           df = ( - erfc (ar) - sqrt (8.d0 / tpi) * ar * exp ( - ar**2) ) &
                / rr**3
           do icart = 1, 3
              na_icart = 3 * (na - 1) + icart
              nb_icart = 3 * (nb - 1) + icart
              do jcart = 1, 3
                 nb_jcart = 3 * (nb - 1) + jcart
                 na_jcart = 3 * (na - 1) + jcart
                 dy3 (na_icart, nb_jcart) = dy3 (na_icart, nb_jcart) + e2 * zv ( &
                      nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * (d2f * alat * &
                      r (icart, nr) * alat * r (jcart, nr) )
                 dy3 (na_icart, na_jcart) = dy3 (na_icart, na_jcart) - e2 * zv ( &
                      nta) * zv (ntb) * (d2f * alat * r (icart, nr) * alat * r (jcart, &
                      nr) )
              enddo
              dy3 (na_icart, nb_icart) = dy3 (na_icart, nb_icart) + e2 * zv ( &
                   nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * df
              dy3 (na_icart, na_icart) = dy3 (na_icart, na_icart) - e2 * zv ( &
                   nta) * zv (ntb) * df
           enddo
        enddo
     enddo

  enddo
#ifdef __PARA
100 continue

  call reduce (18 * nat * nat, dy3)
#endif
  !
  !   The dynamical matrix was computed in cartesian axis and now we put
  !   it on the basis of the modes
  !
  do nu_i = 1, nmodes
     do nu_j = 1, nmodes
        work = (0.d0, 0.d0)
        do nb_jcart = 1, 3 * nat
           do na_icart = 1, 3 * nat
              work = work + conjg (u (na_icart, nu_i) ) * dy3 (na_icart, &
                   nb_jcart) * u (nb_jcart, nu_j)
           enddo
        enddo
        dyn (nu_i, nu_j) = dyn (nu_i, nu_j) - work
     enddo

  enddo
  call stop_clock ('d2ionq')
  return
end subroutine d2ionq
