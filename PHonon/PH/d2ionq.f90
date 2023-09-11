!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d2ionq( nat, ntyp, ityp, zv, tau, alat, omega, q, at, &
   !-----------------------------------------------------------------------
                   bg, g, gg, ngm, gcutm, nmodes, u, dyn )
  !! This routine computes the contribution of the ions to the
  !! dynamical matrix. Both the real and reciprocal space terms
  !! are included.
  !
  !! The original routine was from C. Bungaro.
  !! Revised 16 oct. 1995 by Andrea Dal Corso.
  !! April 1997: parallel stuff added (SdG).
  !! 2021-2022: Refactored by P. Bonfa with parallelization over atoms and openMP parallelization
  !! 2022: restrict outer loops to selected atoms if nat_todo specified in input (P. Delugas)
  !
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE constants,      ONLY : e2, tpi, fpi
  USE mp_bands,       ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp,             ONLY : mp_sum
  USE Coul_cut_2D,    ONLY : do_cutoff_2D, cutoff_2D 
  USE Coul_cut_2D_ph, ONLY : cutoff_2D_qg
  USE partial,        ONLY : atomo, nat_todo, nat_todo_input, set_local_atomo
  USE symm_base,      ONLY : irt 
  USE lr_symm_base,   ONLY : nsymq
  !
  implicit none
  !
  integer :: nat
  !! input: the number of atoms
  integer :: ntyp
  !! input: the number of types of atoms
  integer :: ngm
  !! input: the number of G vectors
  integer :: ityp(nat)
  !! input: the type of each atom
  integer :: nmodes
  !! input: the number of modes
  real(DP) :: tau(3,nat)
  !! input: the positions of the atoms
  real(DP) :: g(3,ngm)
  !! input: the coordinates of g vectors
  real(DP) :: gg(ngm)
  !! input: the modulus of g vectors
  real(DP) :: zv(ntyp)
  !! input: the charge of each type
  real(DP) :: at(3,3)
  !! input: the direct lattice vectors
  real(DP) :: bg(3,3)
  !! input: the reciprocal lattice vectors
  real(DP) :: omega
  !! input: the volume of the unit cell
  real(DP) :: alat
  !! input: the length scale
  real(DP) :: gcutm
  !! input: cut-off of g vectors
  real(DP) :: q(3)
  !! input: the q vector
  complex(DP) :: dyn(3*nat,nmodes)
  !! output: the ionic part of the dynamical matrix
  complex(DP) :: u(3*nat,nmodes)
  !! input: the pattern of the modes
  !
  ! ... local variables
  !
  integer, parameter :: mxr = 100
  ! the maximum number of r shells
  integer ::  nat_l
  ! total number of effectively displaced atoms 
  integer, allocatable :: atomo_l(:) 
  ! list of atoms effectively displaced and work space 
  integer :: na, nb, nta, ntb, ng, nrm, nr, icart, &
       jcart, na_icart, na_jcart, nb_icart, nb_jcart, na_l, nb_l 
  ! counters
  real(DP) :: arg, argq, tpiba2, alpha, r (3, mxr), r2 (mxr), &
       dtau (3), rmax, rr, upperbound, charge, df, d2f, ar, &
       gtq2, gt2, qrg, zvab
  real(DP), allocatable :: fac(:), facq(:)
  ! auxiliary variables

  complex(DP) ::  dy1(3,3), dy2 (3,3), facg, fnat
  complex(DP), allocatable :: dy3 (: , :)
  ! work spaces, factors
  real(DP), external :: qe_erfc
  integer            :: na_s, na_e, mykey
  integer            :: na_i, nb_i
  integer            :: allocstat

  call start_clock ('d2ionq')
  !
  if (nat_todo_input > 0 ) then 
    call set_local_atomo(nat, nat_todo, atomo, nsymq, irt, nat_l, atomo_l) 
  else 
    nat_l = nat 
  end if 
  !  
  allocate(dy3(3*nat, nmodes) , stat = allocstat)
  call errore ("d2ionq:", "failed allocation of workspace", allocstat)
  !
  tpiba2 = (tpi / alat) **2
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

  if (alpha == 0.d0) call errore ('d2ionq', 'optimal alpha not found',1)

  upperbound = 2.d0 * charge**2 * sqrt (2.d0 * alpha / tpi) * &
       erfc ( sqrt (tpiba2 * gcutm / 4.d0 / alpha) )

  if (upperbound > 1.d-9) goto 11

  WRITE( stdout, '(/5x,"Alpha used in Ewald sum = ",f8.4)') alpha
  dy3 (:,:) = (0.d0, 0.d0)
  !
  ! Prepare coefficients
  !
  allocate(facq(ngm))
  allocate(fac(ngm))
!$omp parallel do private(gtq2, gt2) shared(facq, fac) if(ngm > 0)
  do ng = 1, ngm
     !
     !
     gtq2 = dot_product(g(:, ng) + q, g(:, ng) + q) * tpiba2
     if (abs (gtq2) > 1.d-8) then
        facq(ng) = - e2*fpi * tpiba2 / omega * exp ( - gtq2 / alpha / 4.d0) / gtq2
        if (do_cutoff_2D) facq(ng) = facq(ng)*cutoff_2D_qg(ng)
     else
        facq(ng) = 0.d0
     endif
     !
     !
     gt2 = gg (ng) * tpiba2
     if (abs (gt2) > 1.d-8) then
        fac(ng) = - e2 * fpi * tpiba2 / omega * exp ( - gt2 / alpha / 4.d0) / gt2
        if (do_cutoff_2D) fac(ng) = fac(ng)*cutoff_2D(ng)
     else
        fac(ng) = 0.d0
     endif
  end do
!$omp end parallel do
  !
  ! G-space sums here
  !
!$omp parallel do &
!$omp private(na, na_l, na_i, nta, dy2, nb, nb_i, ntb, zvab, dtau, argq, dy1) &
!$omp private(ng, arg, fnat, facg, icart, jcart) shared(dy3,atomo_l) &
!$omp schedule(static) if(ngm > 0)
  do na_l = 1, nat_l
     if (nat_l < nat) then 
       na = atomo_l(na_l) 
     else 
       na = na_l 
     end if 
     na_i = 3 * (na - 1)
     nta = ityp (na)
     !
     dy2 (:,:) = (0.d0, 0.d0)
     do nb = 1, nat
        nb_i = 3 * (nb - 1)
        ntb = ityp (nb)
        zvab = zv (nta)  * zv (ntb)
        dtau = tau (:, na) - tau (:, nb)

        argq = tpi * dot_product( q , dtau )

        dy1 (:,:) = (0.d0, 0.d0)

        do ng = 1, ngm
           !
           !if (fac(ng) == 0.d0  .and. facq(ng) == 0.d0) cycle
           !
           arg = tpi * dot_product( g (:, ng) , dtau )
           !
           fnat = fac(ng) * zvab * CMPLX(cos (arg), 0.d0, kind=DP)
           facg = facq(ng) * zvab * CMPLX(cos (arg + argq), sin (arg + argq), kind=DP)
           do icart = 1, 3
              do jcart = 1, 3
                 dy1 (icart, jcart) = dy1 (icart, jcart) + facg * (q (icart) + &
                   g (icart, ng) ) * (q (jcart) + g (jcart, ng) )
                 dy2(icart , jcart) = dy2(icart,jcart) +  fnat * g(icart, ng) * g(jcart, ng)
              enddo
           enddo
        enddo
        dy3( na_i + 1: na_i+3,  nb_i + 1:nb_i+3) = dy3( na_i+1: na_i+3,  nb_i+1:nb_i+3) + dy1(1:3,1:3)
     enddo
     dy3( na_i + 1: na_i+3,  na_i + 1:na_i+3) = dy3( na_i+1: na_i+3,  na_i+1:na_i+3) - dy2(1:3,1:3)
  enddo
!$omp end parallel do
  deallocate(facq, fac)
  !
  CALL block_distribute(nat_l, me_bgrp, nproc_bgrp, na_s, na_e, mykey ) 
  !  Then there is also a part in real space which is computed here.
  !   ... only by the nodes for which mykey == 0 see PW/src/ewald.f90 for more.
  !
  if ( mykey>0 ) goto 100
  !
  rmax = 5.d0 / sqrt (alpha) / alat
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
  !
  
  do na_l  = na_s, na_e 
     do nb = 1, nat
        if (nat_l < nat ) then 
           na = atomo_l(na_l)
        else 
           na = na_l
        end if 
        nta = ityp (na)
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
           qrg = tpi * (q (1) * (r (1, nr) + dtau (1) ) + &
                        q (2) * (r (2, nr) + dtau (2) ) + &
                        q (3) * (r (3, nr) + dtau (3) ) )
           d2f = (3.d0 * erfc (ar) + sqrt (8.d0 / tpi) * ar * &
                    (3.d0 + 2.d0 * ar**2) * exp ( - ar**2) ) / rr**5
           df = ( - erfc (ar) - sqrt (8.d0 / tpi) * ar * exp ( - ar**2) ) &
                / rr**3
           do icart = 1, 3
              na_icart = 3 * (na - 1) + icart
              nb_icart = 3 * (nb - 1) + icart
              do jcart = 1, 3
                 nb_jcart = 3 * (nb - 1) + jcart
                 na_jcart = 3 * (na - 1) + jcart
                 dy3 (na_icart, nb_jcart) = dy3 (na_icart, nb_jcart) + &
                      e2 * zv (nta) * zv (ntb) * CMPLX(cos (qrg), sin (qrg),kind=DP)&
                      * (d2f * alat * r (icart, nr) * alat * r (jcart, nr) )
                 dy3 (na_icart, na_jcart) = dy3 (na_icart, na_jcart) - &
                      e2 * zv (nta) * zv (ntb) * (d2f * alat * r (icart, nr) *&
                      alat * r (jcart, nr) )
              enddo
              dy3 (na_icart, nb_icart) = dy3 (na_icart, nb_icart) + e2 * &
                   zv (nta) * zv (ntb) * CMPLX(cos (qrg), sin (qrg) ,kind=DP) * df
              dy3 (na_icart, na_icart) = dy3 (na_icart, na_icart) - e2 * &
                   zv (nta) * zv (ntb) * df
           enddo
        enddo 
     enddo
  enddo
100 continue
  call mp_sum ( dy3, intra_bgrp_comm )
  !
  !   The dynamical matrix was computed in Cartesian axis and now we put
  !   it on the basis of the modes
  !
  dy3 = -dy3
  !
  CALL rotate_pattern_add(nat, u, dyn, dy3)
  !
  deallocate(dy3)
  call stop_clock ('d2ionq')
  return
end subroutine d2ionq
