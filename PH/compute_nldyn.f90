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
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  USE uspp_param, ONLY: nh
  use phcom
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  implicit none

  complex(DP) :: becq (nkb, npol, nbnd, nksq), alpq(nkb, npol, nbnd, 3, nksq), &
       wdyn (3 * nat, 3 * nat)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k}
  ! output: the term of the dynamical matrix

  real(DP) :: wgg (nbnd, nbnd, nksq)
  ! input: the weights

  complex(DP) :: ps, aux1 (nbnd), aux2 (nbnd)
  complex(DP), allocatable ::  ps1 (:,:), ps2 (:,:,:), ps3 (:,:), ps4 (:,:,:)
  complex(DP), allocatable ::  ps1_nc(:,:,:), ps2_nc(:,:,:,:), &
                               ps3_nc (:,:,:), ps4_nc (:,:,:,:)
  ! work space
  complex(DP) ::  dynwrk (3 * nat, 3 * nat), ps_nc(2)
  ! auxiliary dynamical matrix

  integer :: ik, ikk, ikq, ibnd, jbnd, ijkb0, ijkb0b, ih, jh, ikb, &
       jkb, ipol, jpol, startb, lastb, na, nb, nt, ntb, nu_i, nu_j, &
       na_icart, na_jcart, mu, nu
  ! counters

  IF (noncolin) THEN
     allocate (ps1_nc (  nkb, npol, nbnd))    
     allocate (ps2_nc (  nkb, npol, nbnd , 3))    
     allocate (ps3_nc (  nkb, npol, nbnd))    
     allocate (ps4_nc (  nkb, npol, nbnd , 3))    
  ELSE
     allocate (ps1 (  nkb, nbnd))    
     allocate (ps2 (  nkb, nbnd , 3))    
     allocate (ps3 (  nkb, nbnd))    
     allocate (ps4 (  nkb, nbnd , 3))    
  END IF

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
     IF (noncolin) THEN
        ps1_nc = (0.d0, 0.d0)
        ps2_nc = (0.d0, 0.d0)
        ps3_nc = (0.d0, 0.d0)
        ps4_nc = (0.d0, 0.d0)
     ELSE
        ps1 = (0.d0, 0.d0)
        ps2 = (0.d0, 0.d0)
        ps3 = (0.d0, 0.d0)
        ps4 = (0.d0, 0.d0)
     END IF
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
                       IF (noncolin) THEN
                          IF (lspinorb) THEN
                             ps1_nc (ikb, 1, ibnd) = ps1_nc (ikb, 1, ibnd) +  &
                                (deeq_nc (ih, jh, na, 1) -                    &
                                   et (ibnd, ikk) * qq_so (ih, jh, 1, nt) ) * &
                                 becp1_nc (jkb, 1, ibnd, ik) +                &
                                (deeq_nc (ih, jh, na, 2) -                    &
                                   et (ibnd, ikk) * qq_so (ih, jh, 2, nt) ) * &
                                 becp1_nc (jkb, 2, ibnd, ik) 
                             ps1_nc (ikb, 2, ibnd) = ps1_nc (ikb, 2, ibnd) + &
                                (deeq_nc (ih, jh, na, 3) -                  &
                                   et (ibnd, ikk) * qq_so (ih, jh, 3, nt) ) * &
                                 becp1_nc (jkb, 1, ibnd, ik) +              &
                                (deeq_nc (ih, jh, na, 4) -                  &
                                   et (ibnd, ikk) * qq_so (ih, jh, 4, nt) ) * &
                                  becp1_nc (jkb, 2, ibnd, ik) 
                             ps3_nc (ikb, 1, ibnd) = ps3_nc (ikb, 1, ibnd) - &
                                qq_so(ih,jh,1,nt)*becq(jkb,1,ibnd,ik) -      &
                                qq_so(ih,jh,2,nt)*becq(jkb,2,ibnd,ik)
                             ps3_nc (ikb, 2, ibnd) = ps3_nc (ikb, 2, ibnd) - &
                                qq_so(ih,jh,3,nt)*becq(jkb,1,ibnd,ik) -      &
                                qq_so(ih,jh,4,nt)*becq(jkb,2,ibnd,ik)
                          ELSE
                             ps1_nc (ikb, 1, ibnd) = ps1_nc (ikb, 1, ibnd) + &
                                 (deeq_nc (ih, jh, na, 1) - &
                                 et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                 becp1_nc (jkb, 1, ibnd, ik) +  &
                                 deeq_nc (ih, jh, na, 2) *  &
                                 becp1_nc (jkb, 2, ibnd, ik) 
                             ps1_nc (ikb, 2, ibnd) = ps1_nc (ikb, 2, ibnd) + &
                                  (deeq_nc (ih, jh, na, 4) - &
                                  et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                  becp1_nc (jkb, 2, ibnd, ik) +  &
                                  deeq_nc (ih, jh, na, 3) *      &
                                  becp1_nc (jkb, 1, ibnd, ik) 
                             ps3_nc (ikb, 1, ibnd) = ps3_nc (ikb, 1, ibnd) - &
                                  qq (ih, jh, nt) * becq (jkb, 1, ibnd, ik)
                             ps3_nc (ikb, 2, ibnd) = ps3_nc (ikb, 2, ibnd) - &
                                  qq (ih, jh, nt) * becq (jkb, 2, ibnd, ik)
                          END IF
                       ELSE
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                            (deeq (ih, jh, na, current_spin) - &
                            et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                            becp1 (jkb, ibnd, ik)
                          ps3 (ikb, ibnd) = ps3 (ikb, ibnd) - &
                            qq (ih, jh, nt) * becq (jkb, 1, ibnd, ik)
                       END IF
                       do ipol = 1, 3
                          IF (noncolin) THEN
                             IF (lspinorb) THEN
                                ps2_nc(ikb,1,ibnd,ipol) =                   &
                                    ps2_nc(ikb,1,ibnd,ipol) +               &
                                   ( deeq_nc(ih, jh, na, 1) -               &
                                     et(ibnd,ikk)*qq_so(ih,jh,1,nt) ) *     &
                                     alphap_nc (jkb, 1, ibnd, ipol, ik) +   &
                                     (deeq_nc(ih, jh, na, 2) -              &
                                      et(ibnd,ikk)*qq_so(ih,jh,2,nt) ) *    &
                                     alphap_nc (jkb, 2, ibnd, ipol, ik) +   &
                                     int1_nc(ih, jh, ipol, na, 1) *         &
                                     becp1_nc (jkb, 1, ibnd, ik)  +         &
                                     int1_nc(ih, jh, ipol, na, 2) *         &
                                     becp1_nc (jkb, 2, ibnd, ik)
                                ps2_nc(ikb,2,ibnd,ipol)= &
                                     ps2_nc(ikb,2,ibnd,ipol)+ &
                                    (deeq_nc(ih, jh, na, 3) -              &
                                        et(ibnd,ikk)*qq_so(ih,jh,3,nt))*   &
                                    alphap_nc (jkb, 1, ibnd, ipol, ik) +   &
                                    (deeq_nc(ih, jh, na, 4) -              &
                                        et(ibnd,ikk)*qq_so(ih,jh,4,nt))*   &
                                    alphap_nc (jkb, 2, ibnd, ipol, ik) +   &
                                    int1_nc(ih, jh, ipol, na, 3) *         &
                                    becp1_nc (jkb, 1, ibnd, ik)  +         &
                                    int1_nc(ih, jh, ipol, na, 4) *         &
                                    becp1_nc (jkb, 2, ibnd, ik)
                                ps4_nc(ikb,1,ibnd,ipol) =               &
                                  ps4_nc(ikb,1,ibnd,ipol)-              &
                                   qq_so(ih,jh,1,nt)*alpq(jkb,1,ibnd,ipol,ik)-&
                                   qq_so(ih,jh,2,nt)*alpq(jkb,2,ibnd,ipol,ik)
                                ps4_nc(ikb,2,ibnd,ipol) =       &
                                   ps4_nc(ikb,2,ibnd,ipol)-     &
                                   qq_so(ih,jh,3,nt)*alpq(jkb,1,ibnd,ipol,ik)-&
                                   qq_so(ih,jh,4,nt)*alpq(jkb,2,ibnd,ipol,ik)
                             ELSE
                                ps2_nc(ikb,1,ibnd,ipol)= &
                                    ps2_nc(ikb,1,ibnd,ipol)+ &
                                    (deeq_nc(ih, jh, na, 1) -         &
                                     et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                     alphap_nc (jkb, 1, ibnd, ipol, ik) +   &
                                     deeq_nc(ih, jh, na, 2) *               &
                                     alphap_nc (jkb, 2, ibnd, ipol, ik) +   &
                                     int1_nc(ih, jh, ipol, na, 1) *           &
                                     becp1_nc (jkb, 1, ibnd, ik) +            &
                                     int1_nc(ih, jh, ipol, na, 2) *           &
                                     becp1_nc (jkb, 2, ibnd, ik)
                                ps2_nc(ikb,2,ibnd,ipol)= &
                                     ps2_nc(ikb,2,ibnd,ipol)+ &
                                    (deeq_nc(ih, jh, na, 4) -              &
                                    et (ibnd, ikk) * qq (ih, jh, nt) ) *   &
                                    alphap_nc (jkb, 2, ibnd, ipol, ik) +   &
                                    deeq_nc(ih, jh, na, 3) *               &
                                    alphap_nc (jkb, 1, ibnd, ipol, ik) +   &
                                    int1_nc(ih, jh, ipol, na, 3) *         &
                                    becp1_nc (jkb, 1, ibnd, ik)   +        &
                                    int1_nc(ih, jh, ipol, na, 4) *         &
                                    becp1_nc (jkb, 2, ibnd, ik)
                                ps4_nc(ikb,1,ibnd,ipol) =                  &
                                     ps4_nc(ikb,1,ibnd,ipol)-              &
                                     qq(ih,jh,nt)*alpq(jkb,1,ibnd,ipol,ik)
                                ps4_nc(ikb,2,ibnd,ipol) =                  &
                                     ps4_nc(ikb,2,ibnd,ipol)-              &
                                     qq(ih,jh,nt)*alpq(jkb,2,ibnd,ipol,ik)
                             END IF
                          ELSE
                             ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                               (deeq (ih, jh,na, current_spin) -             &
                               et (ibnd, ikk) * qq (ih, jh, nt) ) *          &
                               alphap (jkb, ibnd, ipol, ik) +                &
                               int1 (ih, jh, ipol, na, current_spin) *       &
                               becp1 (jkb, ibnd, ik)
                             ps4 (ikb, ibnd, ipol) = ps4 (ikb, ibnd, ipol) - &
                               qq (ih, jh, nt) * alpq (jkb, 1, ibnd, ipol, ik)
                          END IF
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
                          IF (noncolin) THEN
                             aux1 (jbnd) = aux1 (jbnd) + &
                            CONJG(alpq(ikb,1,jbnd,ipol,ik))*ps1_nc(ikb,1,ibnd)+&
                            CONJG(becq(ikb,1,jbnd,ik))*ps2_nc(ikb,1,ibnd,ipol)+&
                            CONJG(alpq(ikb,2,jbnd,ipol,ik))*ps1_nc(ikb,2,ibnd)+&
                            CONJG(becq(ikb,2,jbnd,ik))*ps2_nc(ikb,2,ibnd,ipol)
                          ELSE
                             aux1 (jbnd) = aux1 (jbnd) + &
                               CONJG(alpq(ikb,1,jbnd,ipol,ik))*ps1(ikb,ibnd)+&
                               CONJG(becq(ikb,1,jbnd,ik))*ps2(ikb,ibnd,ipol)
                          END IF
                       enddo
                    enddo
                    ijkb0b = 0
                    do ntb = 1, ntyp
                       do nb = 1, nat
                          if (ityp (nb) == ntb) then
                             do ih = 1, nh (ntb)
                                ikb = ijkb0b + ih
                                ps_nc =(0.d0,0.d0)
                                ps = (0.d0, 0.d0)
                                do jh = 1, nh (ntb)
                                   jkb = ijkb0b + jh
                                   IF (noncolin) THEN
                                      IF (lspinorb) THEN
                                         ps_nc(1) = ps_nc(1) + &
                                               int2_so(ih,jh,ipol,na,nb,1)*&
                                               becp1_nc(jkb,1,ibnd,ik) +   &
                                               int2_so(ih,jh,ipol,na,nb,2)*&
                                               becp1_nc(jkb,2,ibnd,ik)
                                         ps_nc(2)=ps_nc(2) + &
                                               int2_so(ih,jh,ipol,na,nb,3)*&
                                               becp1_nc(jkb,1,ibnd,ik) +   &
                                               int2_so(ih,jh,ipol,na,nb,4)*&
                                               becp1_nc(jkb,2,ibnd,ik)
                                      ELSE
                                         ps_nc(1) = ps_nc(1) + &
                                               int2(ih,jh,ipol,na,nb)*&
                                               becp1_nc(jkb,1,ibnd,ik)
                                         ps_nc(2) = ps_nc(2) + &
                                               int2(ih,jh,ipol,na,nb)*&
                                               becp1_nc(jkb,2,ibnd,ik)
                                      ENDIF
                                   ELSE
                                      ps = ps + int2 (ih, jh, ipol, na, nb) * &
                                             becp1 (jkb, ibnd,ik)
                                   END IF
                                enddo
                                do jbnd = startb, lastb
                                   IF (noncolin) THEN
                                      aux1(jbnd) = aux1 (jbnd) + &
                                        ps_nc(1)*CONJG(becq(ikb,1,jbnd,ik))+&
                                        ps_nc(2)*CONJG(becq(ikb,2,jbnd,ik))
                                   ELSE
                                      aux1(jbnd) = aux1 (jbnd) + &
                                        ps * CONJG(becq(ikb,1,jbnd,ik))
                                   END IF
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
                                      IF (noncolin) THEN
                                         aux2 (jbnd) = aux2 (jbnd) + &
                                           wgg(ibnd, jbnd, ik) * &
                                        (CONJG(alphap_nc(ikb,1,ibnd,jpol,ik))*&
                                            ps3_nc (ikb, 1, jbnd) + &
                                            CONJG(becp1_nc (ikb,1,ibnd, ik)) * &
                                            ps4_nc (ikb, 1, jbnd, jpol) +  &
                                         CONJG(alphap_nc(ikb,2,ibnd,jpol,ik))* &
                                            ps3_nc (ikb,2,jbnd) + &
                                            CONJG(becp1_nc (ikb,2,ibnd,ik)) * &
                                            ps4_nc (ikb, 2, jbnd, jpol) )
                                      ELSE
                                         aux2 (jbnd) = aux2 (jbnd) + &
                                           wgg (ibnd, jbnd, ik) * &
                                           (CONJG(alphap(ikb,ibnd,jpol,ik)) * &
                                            ps3 (ikb, jbnd) + &
                                            CONJG(becp1 (ikb, ibnd, ik) ) * &
                                            ps4 (ikb, jbnd, jpol) )
                                      END IF
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
  call mp_sum ( dynwrk, intra_pool_comm )
#endif
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        ps = (0.0d0, 0.0d0)
        do na_jcart = 1, 3 * nat
           do na_icart = 1, 3 * nat
              ps = ps + CONJG(u (na_icart, nu_i) ) * dynwrk (na_icart, &
                   na_jcart) * u (na_jcart, nu_j)
           enddo
        enddo
        wdyn (nu_i, nu_j) = wdyn (nu_i, nu_j) + ps
     enddo
  enddo
  !      call tra_write_matrix('nldyn wdyn',wdyn,u,nat)
  !      call stop_ph(.true.)
  IF (noncolin) THEN
     deallocate (ps4_nc)
     deallocate (ps3_nc)
     deallocate (ps2_nc)
     deallocate (ps1_nc)
  ELSE
     deallocate (ps4)
     deallocate (ps3)
     deallocate (ps2)
     deallocate (ps1)
  END IF
  return
end subroutine compute_nldyn
