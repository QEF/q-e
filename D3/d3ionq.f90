!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine d3ionq (nat, ntyp, ityp, zv, tau, alat, omega, q, at, &
     bg, g, gg, ngm, gcutm, nmodes, u, ug0, npert_i, npert_f, q0mode, &
     d3dyn)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the contribution of the ions to the
  !    third order derivative of the total energy.
  !    Both the real and reciprocal space terms
  !    are included.
  !
#include "machine.h"
  implicit none
  !
  !  first the I/O variables
  !

  integer :: nat, ntyp, ngm, ityp (nat), nmodes, npert_i, npert_f
  ! input: the number of atoms
  ! input: the number of types of atoms
  ! input: the number of G vectors
  ! input: the type of each atom
  ! input: the number of modes

  real (8) :: tau (3, nat), g (3, ngm), gg (ngm), zv (ntyp), &
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

  complex (8) :: d3dyn (3 * nat, nmodes, 3 * nat), &
       u (3 * nat, nmodes), ug0 (3 * nat, nmodes)
  ! output: derivative of the dyn. matrix
  ! input: the pattern of the modes
  ! input: the pattern of the modes (q=0)
  logical :: q0mode (300)
  ! input: if .true. this mode is to be co
  !
  !   three parameters
  !

  integer :: mxr
  ! the maximum number of r shells

  real (8) :: e2, tpi
  ! the electron charge
  ! twp times pi

  parameter (mxr = 100, e2 = 2.d0, tpi = 2.d0 * 3.14159265358979d0)

  integer :: nu_i, nu_j, nu_k, na, nb, nta, ntb, ng, nrm, nr, icart, &
       jcart, kcart, na_icart, nb_jcart, nc_kcart
  ! counter on the modes
  ! counter on the modes
  ! counter on the modes
  ! counter on atoms
  ! counter on atoms
  ! the type of atom na
  ! the type of atom nb
  ! counter on G vectors
  ! the real number of R shells
  ! counter on r shells
  ! counters on cartesian coordinates

  real (8) :: arg, fpi, argq, tpiba2, tpiba3, alpha, erfc, &
       upperbound, charge, fac, gtq2, gt2, facq, d2f, d3f, rmax, r (3, &
       mxr), r2 (mxr), dtau (3), rr, ar, qrg
  ! the argument of the phase
  ! 4 times pi
  ! the argument of the phase
  ! two pi / alat
  ! two pi / alat
  ! the parameter alpha used to converge the summ
  ! the erfc function
  ! an estimate of the error in the real term
  ! the total charge
  ! the phase factor
  ! the modulus of |q+G|
  ! the modulus of G
  ! a phase factor
  ! derivatives for the real-space term
  ! cutoff for real space term
  ! the R_i-\tau_s-\tau_s' parameter
  ! the square of this parameter

  complex (8), allocatable :: d3dy1 (:,:,:), d3dy2 (:,:,:), d3dy3 (:,:,:)
  ! first term dynamical matrix
  ! second term dynamical matrix
  ! third term dynamical matrix
  complex (8) :: facg, fnat,  work
  ! a phase with g
  ! a phase with the atoms
  ! working space

  allocate  (d3dy1( 3 * nat, nmodes, 3 * nat))    
  allocate  (d3dy2( 3 * nat, nmodes, 3 * nat))    
  allocate  (d3dy3( 3 * nat, nmodes, 3 * nat))    
  tpiba2 = (tpi / alat) **2
  tpiba3 = (tpi / alat) **3

  fpi = 2.d0 * tpi
  charge = 0.d0
  do na = 1, nat
     charge = charge+zv (ityp (na) )
  enddo
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is an upper bound for the error in the sum over G
  ! estimated for the energy (empirical trust!)
  ! (PG: appropriate for c60 = 0.2)
  !
  alpha = 2.9d0

11 alpha = alpha - 0.1d0

  if (alpha.eq.0.d0) call error ('d3ion', 'optimal alpha not found', &
       1)

  upperbound = 2.d0 * charge**2 * sqrt (2.d0 * alpha / tpi) * erfc ( &
       sqrt (tpiba2 * gcutm / 4.d0 / alpha) )

  if (upperbound.gt.1.d-9) goto 11

  write (6, '(/5x,"Alpha used in Ewald sum = ",f8.4)') alpha
  call setv (2 * 9 * nat * nat * nmodes, 0.d0, d3dy1, 1)
  call setv (2 * 9 * nat * nat * nmodes, 0.d0, d3dy2, 1)

  call setv (2 * 9 * nat * nat * nmodes, 0.d0, d3dy3, 1)
  do ng = 1, ngm

     gt2 = gg (ng) * tpiba2
     if (abs (gt2) .gt.1.d-8) then
        fac = e2 * fpi * tpiba3 / omega * exp ( - gt2 / alpha / 4.d0) &
             / gt2
     else
        fac = 0.d0

     endif
     do nu_i = npert_i, npert_f
        na = (nu_i - 1) / 3 + 1
        icart = nu_i - 3 * (na - 1)

        nta = ityp (na)
        if (nu_i.eq.npert_i.or.icart.eq.1) then
           fnat = (0.d0, 0.d0)
           do nb = 1, nat
              ntb = ityp (nb)
              arg = tpi * ( (g (1, ng) ) * (tau (1, na) - tau (1, nb) ) &
                   + (g (2, ng) ) * (tau (2, na) - tau (2, nb) ) + (g (3, ng) ) &
                   * (tau (3, na) - tau (3, nb) ) )
              facg = fac * zv (nta) * zv (ntb) * cmplx (sin (arg), 0.d0)
              fnat = fnat + facg
           enddo

        endif
        do jcart = 1, 3
           nu_j = 3 * (na - 1) + jcart
           do kcart = 1, 3
              nu_k = 3 * (na - 1) + kcart
              d3dy1 (nu_i, nu_j, nu_k) = d3dy1 (nu_i, nu_j, nu_k) + fnat * g ( &
                   icart, ng) * g (jcart, ng) * g (kcart, ng)
           enddo
        enddo

     enddo
     do nu_i = npert_i, npert_f
        na = (nu_i - 1) / 3 + 1
        icart = nu_i - 3 * (na - 1)

        nta = ityp (na)
        do nb = 1, nat
           ntb = ityp (nb)
           arg = tpi * ( (g (1, ng) ) * (tau (1, na) - tau (1, nb) ) + &
                (g (2, ng) ) * (tau (2, na) - tau (2, nb) ) + (g (3, ng) ) &
                * (tau (3, na) - tau (3, nb) ) )
           fnat = fac * zv (nta) * zv (ntb) * cmplx (sin (arg), 0.d0)
           do jcart = 1, 3
              nu_j = 3 * (nb - 1) + jcart
              do kcart = 1, 3
                 nu_k = 3 * (nb - 1) + kcart
                 d3dy1 (nu_i, nu_j, nu_k) = d3dy1 (nu_i, nu_j, nu_k) + fnat * g ( &
                      icart, ng) * g (jcart, ng) * g (kcart, ng)
              enddo
           enddo
        enddo

     enddo
     gtq2 = ( (g (1, ng) + q (1) ) **2 + (g (2, ng) + q (2) ) **2 + &
          (g (3, ng) + q (3) ) **2) * tpiba2
     if (abs (gtq2) .gt.1.d-8) then
        facq = e2 * fpi * tpiba3 / omega * exp ( - gtq2 / alpha / 4.d0) &
             / gtq2
     else
        facq = 0.d0

     endif
     do nu_i = npert_i, npert_f
        na = (nu_i - 1) / 3 + 1
        icart = nu_i - 3 * (na - 1)

        nta = ityp (na)
        do nb = 1, nat
           ntb = ityp (nb)
           argq = tpi * ( (g (1, ng) + q (1) ) * (tau (1, nb) - tau (1, na) ) &
                + (g (2, ng) + q (2) ) * (tau (2, nb) - tau (2, na) ) + (g (3, ng) &
                + q (3) ) * (tau (3, nb) - tau (3, na) ) )
           facg = facq * zv (nta) * zv (ntb) * cmplx ( - sin (argq), - cos ( &
                argq) )
           do jcart = 1, 3
              nu_j = 3 * (nb - 1) + jcart
              do kcart = 1, 3
                 nu_k = 3 * (nb - 1) + kcart
                 d3dy2 (nu_j, nu_i, nu_k) = d3dy2 (nu_j, nu_i, nu_k) + facg * &
                      (q (icart) + g (icart, ng) ) * (q (jcart) + g (jcart, ng) ) &
                      * (q (kcart) + g (kcart, ng) )
                 d3dy3 (nu_j, nu_k, nu_i) = d3dy3 (nu_j, nu_k, nu_i) - conjg (facg) &
                      * (q (icart) + g (icart, ng) ) * (q (jcart) + g (jcart, ng) ) &
                      * (q (kcart) + g (kcart, ng) )
              enddo
           enddo
        enddo

     enddo
  enddo
  call DAXPY (2 * 27 * nat * nat * nat, 1.d0, d3dy2, 1, d3dy1, 1)
  call DAXPY (2 * 27 * nat * nat * nat, - 1.d0, d3dy3, 1, d3dy1, 1)
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! uncomment next line if you want to skip real space term for testing pu
  !      if (.true.) goto 100
  !
  !  Then there is also a part in real space which is computed here.
#ifdef PARA
  !   ... only by the node that contains  G=0
  !
  if (gg (1) .gt.1.d-8) goto 100
#endif
  call setv (2 * 27 * nat * nat * nat, 0.d0, d3dy2, 1)
  !
  rmax = 5.d0 / sqrt (alpha) / alat
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
  !
  do nu_i = npert_i, npert_f
     na = (nu_i - 1) / 3 + 1
     icart = nu_i - (na - 1) * 3
     nta = ityp (na)
     do nb = 1, nat
        ntb = ityp (nb)
        do jcart = 1, 3
           dtau (jcart) = tau (jcart, na) - tau (jcart, nb)
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
           d2f = (3.d0 * erfc (ar) + sqrt (8.d0 / tpi) * exp ( - ar**2) &
                * ar * (3.d0 + 2.d0 * ar**2) ) / rr**5

           d3f = ( - 15.d0 * erfc (ar) - sqrt (8.d0 / tpi) * exp ( - ar**2) &
                * ar * (15.d0 + 10.d0 * ar**2 + 4.d0 * ar**4) ) / rr**7
           do jcart = 1, 3
              nu_j = (nb - 1) * 3 + jcart
              !
              ! nc = nb case
              !
              do kcart = 1, 3
                 nu_k = (nb - 1) * 3 + kcart
                 d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) - e2 * zv ( &
                      nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * d3f * alat**3 * &
                      r (icart, nr) * r (jcart, nr) * r (kcart, nr)
                 if (icart.eq.jcart) d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, &
                      nu_k) - e2 * zv (nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) &
                      * d2f * alat * r (kcart, nr)
              enddo
              nu_k = (nb - 1) * 3 + icart
              d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) - e2 * zv ( &
                   nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * d2f * alat * r &
                   (jcart, nr)
              nu_k = (nb - 1) * 3 + jcart
              d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) - e2 * zv ( &
                   nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * d2f * alat * r &
                   (icart, nr)
              !
              ! nc = na case
              !
              do kcart = 1, 3
                 nu_k = (na - 1) * 3 + kcart
                 d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) + e2 * zv ( &
                      nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * d3f * alat**3 * &
                      r (icart, nr) * r (jcart, nr) * r (kcart, nr)
                 if (icart.eq.jcart) d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, &
                      nu_k) + e2 * zv (nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) &
                      * d2f * alat * r (kcart, nr)
              enddo
              nu_k = (na - 1) * 3 + icart
              d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) + e2 * zv ( &
                   nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * d2f * alat * r &
                   (jcart, nr)
              nu_k = (na - 1) * 3 + jcart
              d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) + e2 * zv ( &
                   nta) * zv (ntb) * DCMPLX (cos (qrg), sin (qrg) ) * d2f * alat * r &
                   (icart, nr)
              !
              ! na = nb case (NB: role of nu_k and nu_j are interchanged)
              !
              do kcart = 1, 3
                 nu_k = (na - 1) * 3 + kcart
                 d3dy2 (nu_i, nu_k, nu_j) = d3dy2 (nu_i, nu_k, nu_j) + e2 * zv ( &
                      nta) * zv (ntb) * d3f * alat**3 * r (icart, nr) * r (jcart, nr) &
                      * r (kcart, nr)
                 if (icart.eq.jcart) d3dy2 (nu_i, nu_k, nu_j) = d3dy2 (nu_i, nu_k, &
                      nu_j) + e2 * zv (nta) * zv (ntb) * d2f * alat * r (kcart, nr)
              enddo
              nu_k = (na - 1) * 3 + icart
              d3dy2 (nu_i, nu_k, nu_j) = d3dy2 (nu_i, nu_k, nu_j) + e2 * zv ( &
                   nta) * zv (ntb) * d2f * alat * r (jcart, nr)
              nu_k = (na - 1) * 3 + jcart
              d3dy2 (nu_i, nu_k, nu_j) = d3dy2 (nu_i, nu_k, nu_j) + e2 * zv ( &
                   nta) * zv (ntb) * d2f * alat * r (icart, nr)
              !
              ! case na=nb=nc
              !
              nu_j = (na - 1) * 3 + jcart
              do kcart = 1, 3
                 nu_k = (na - 1) * 3 + kcart
                 d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) - e2 * zv ( &
                      nta) * zv (ntb) * d3f * alat**3 * r (icart, nr) * r (jcart, nr) &
                      * r (kcart, nr)
                 if (icart.eq.jcart) d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, &
                      nu_k) - e2 * zv (nta) * zv (ntb) * d2f * alat * r (kcart, nr)
              enddo
              nu_k = (na - 1) * 3 + icart
              d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) - e2 * zv ( &
                   nta) * zv (ntb) * d2f * alat * r (jcart, nr)
              nu_k = (na - 1) * 3 + jcart
              d3dy2 (nu_i, nu_j, nu_k) = d3dy2 (nu_i, nu_j, nu_k) - e2 * zv ( &
                   nta) * zv (ntb) * d2f * alat * r (icart, nr)
           enddo
        enddo
     enddo
  enddo

  call DAXPY (2 * 27 * nat * nat * nat, 1.d0, d3dy2, 1, d3dy1, 1)
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

100 continue
#ifdef PARA
  call reduce (2 * 27 * nat * nat * nat, d3dy1)
  call poolreduce (2 * 27 * nat * nat * nat, d3dy1)
#endif
  !
  !   The dynamical matrix was computed in cartesian axis and now we put
  !   it on the basis of the modes; d3dy2 used as working array
  !
  call setv (2 * 27 * nat * nat * nat, 0.d0, d3dy2, 1)
  do nu_k = npert_i, npert_f
     if (q0mode (nu_k) ) then
        do nu_i = 1, 3 * nat
           do nu_j = 1, 3 * nat
              work = (0.d0, 0.d0)
              do nc_kcart = 1, 3 * nat
                 do na_icart = 1, 3 * nat
                    do nb_jcart = 1, 3 * nat
                       work = work + ug0 (nc_kcart, nu_k) * conjg (u (na_icart, nu_i) &
                            ) * d3dy1 (nc_kcart, na_icart, nb_jcart) * u (nb_jcart, nu_j)
                    enddo
                 enddo
              enddo
              d3dy2 (nu_k, nu_i, nu_j) = work
           enddo
        enddo
     endif
  enddo
#ifdef PARA
  call poolreduce (2 * 27 * nat * nat * nat, d3dy2)
#endif


  call DAXPY (2 * 27 * nat * nat * nat, 1.d0, d3dy2, 1, d3dyn, 1)
  !      call writed3dyn_5(d3dy2,'d3ewald',-1)
  deallocate (d3dy3)
  deallocate (d3dy2)
  deallocate (d3dy1)
  return
end subroutine d3ionq
