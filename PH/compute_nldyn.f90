!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_nldyn (wdyn, wgg, becq, alpq)
  !-----------------------------------------------------------------------
  !
  !
  !  This routine compute the term of the dynamical matrix due to
  !  the orthogonality constraint. Only the part which is due to
  !  the nonlocal terms is computed here
  !
#include "machine.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE kinds, only : DP
  USE uspp_param, ONLY: nh
  use phcom
  implicit none

  complex(kind=DP) :: becq (nkb, nbnd, nksq), alpq (nkb, nbnd, 3, nksq), &
       wdyn (3 * nat, 3 * nat)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k}
  ! output: the term of the dynamical matrix

  real(kind=DP) :: wgg (nbnd, nbnd, nksq)
  ! input: the weights

  complex(kind=DP) :: ps, aux1 (nbnd), aux2 (nbnd)
  complex(kind=DP), allocatable ::  ps1 (:,:),&
       ps2 (:,:,:), ps3 (:,:), ps4 (:,:,:)
  ! work space
  complex(kind=DP) ::  dynwrk (3 * nat, 3 * nat)
  ! auxiliary dynamical matrix

  integer :: ik, ikk, ikq, ibnd, jbnd, ijkb0, ijkb0b, ih, jh, ikb, &
       jkb, ipol, jpol, startb, lastb, na, nb, nt, ntb, nu_i, nu_j, &
       na_icart, na_jcart, mu, nu
  ! counters

  allocate (ps1 (  nkb , nbnd))    
  allocate (ps2 (  nkb , nbnd , 3))    
  allocate (ps3 (  nkb , nbnd))    
  allocate (ps4 (  nkb , nbnd , 3))    

  dynwrk (:,:) = (0.d0, 0.d0)
  call divide (nbnd, startb, lastb)
  do ik = 1, nksq
     if (lgamma) then
        ikk = ik
        ikq = ik
     else
        ikk = 2 * ik - 1
        ikq = ikk + 1
     endif

     if (lsda) current_spin = isk (ikk)
     ps1 (:,:) = (0.d0, 0.d0)
     ps2 (:,:,:) = (0.d0, 0.d0)
     ps3 (:,:) = (0.d0, 0.d0)
     ps4 (:,:,:) = (0.d0, 0.d0)
     !
     !   Here we prepare the two terms
     !
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) == nt) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd
                       ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                            (deeq (ih, jh, na, current_spin) - &
                            et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                            becp1 (jkb, ibnd, ik)
                       ps3 (ikb, ibnd) = ps3 (ikb, ibnd) - &
                            qq (ih, jh, nt) * becq (jkb, ibnd, ik)
                       do ipol = 1, 3
                          ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                               (deeq (ih, jh,na, current_spin) - &
                               et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                               alphap (jkb, ibnd, ipol, ik) + &
                               int1 (ih, jh, ipol, na, current_spin) * &
                               becp1 (jkb, ibnd, ik)
                          ps4 (ikb, ibnd, ipol) = ps4 (ikb, ibnd, ipol) - &
                               qq (ih, jh, nt) * alpq (jkb, ibnd, ipol, ik)
                       enddo
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     enddo
     !
     !     Here starts the loop on the atoms (rows)
     !
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ipol = 1, 3
                 mu = 3 * (na - 1) + ipol
                 do ibnd = 1, nbnd_occ (ikk)
                    aux1 (:) = (0.d0, 0.d0)
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       do jbnd = startb, lastb
                          aux1 (jbnd) = aux1 (jbnd) + &
                               conjg(alpq(ikb,jbnd,ipol,ik)) * ps1(ikb,ibnd)+&
                               conjg(becq(ikb,jbnd,ik)) * ps2(ikb,ibnd,ipol)
                       enddo
                    enddo
                    ijkb0b = 0
                    do ntb = 1, ntyp
                       do nb = 1, nat
                          if (ityp (nb) == ntb) then
                             do ih = 1, nh (ntb)
                                ikb = ijkb0b + ih
                                ps = (0.d0, 0.d0)
                                do jh = 1, nh (ntb)
                                   jkb = ijkb0b + jh
                                   ps = ps + int2 (ih, jh, ipol, na, nb) * &
                                             becp1 (jkb, ibnd,ik)
                                enddo
                                do jbnd = startb, lastb
                                   aux1 (jbnd) = aux1 (jbnd) + &
                                        ps * conjg (becq (ikb, jbnd, ik))
                                enddo
                             enddo
                             ijkb0b = ijkb0b + nh (ntb)
                          endif
                       enddo
                    enddo
                    !
                    !     here starts the second loop on the atoms
                    !
                    ijkb0b = 0
                    do ntb = 1, ntyp
                       do nb = 1, nat
                          if (ityp (nb) == ntb) then
                             do jpol = 1, 3
                                nu = 3 * (nb - 1) + jpol
                                aux2 (:) = (0.d0, 0.d0)
                                do ih = 1, nh (ntb)
                                   ikb = ijkb0b + ih
                                   do jbnd = startb, lastb
                                      aux2 (jbnd) = aux2 (jbnd) + &
                                           wgg (ibnd, jbnd, ik) * &
                                           (conjg(alphap(ikb,ibnd,jpol,ik)) * &
                                            ps3 (ikb, jbnd) + &
                                            conjg(becp1 (ikb, ibnd, ik) ) * &
                                            ps4 (ikb, jbnd, jpol) )
                                   enddo
                                enddo
                                do jbnd = startb, lastb
                                   dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0*wk(ikk) * aux2(jbnd) * aux1(jbnd)
                                enddo
                             enddo
                             ijkb0b = ijkb0b + nh (ntb)
                          endif
                       enddo
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     enddo
  enddo
#ifdef __PARA
  call reduce (2 * 3 * nat * 3 * nat, dynwrk)
#endif
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        ps = (0.0d0, 0.0d0)
        do na_jcart = 1, 3 * nat
           do na_icart = 1, 3 * nat
              ps = ps + conjg (u (na_icart, nu_i) ) * dynwrk (na_icart, &
                   na_jcart) * u (na_jcart, nu_j)
           enddo
        enddo
        wdyn (nu_i, nu_j) = wdyn (nu_i, nu_j) + ps
     enddo
  enddo
  !      call tra_write_matrix('nldyn wdyn',wdyn,u,nat)
  !      call stop_ph(.true.)
  deallocate (ps4)
  deallocate (ps3)
  deallocate (ps2)
  deallocate (ps1)
  return
end subroutine compute_nldyn
