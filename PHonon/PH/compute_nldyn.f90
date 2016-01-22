!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
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
  !  This routine computes the term of the dynamical matrix due to
  !  the orthogonality constraint. Only the part which is due to
  !  the nonlocal terms is computed here
  !
  USE kinds,     ONLY : DP
  USE klist,     ONLY : wk
  USE lsda_mod,  ONLY : lsda, current_spin, isk, nspin
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp,      ONLY : nkb, qq, qq_so
  USE uspp_param,ONLY : nh, nhm
  USE spin_orb,  ONLY : lspinorb
  USE wvfct,     ONLY : nbnd, et

  USE modes,     ONLY : u
  USE phus,      ONLY : alphap, int1, int2, &
                        int2_so, int1_nc
  USE control_ph, ONLY : rec_code_read

  USE lrus,       ONLY : becp1
  USE qpoint,     ONLY : nksq, ikks, ikqs
  USE control_lr, ONLY : nbnd_occ

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE becmod,    ONLY : bec_type

  implicit none

  type (bec_type) :: becq (nksq),   & ! input: the becp with psi_{k+q}
                     alpq(3, nksq)
  complex(DP) :: wdyn (3 * nat, 3 * nat)
  ! input: the alphap with psi_{k}
  ! output: the term of the dynamical matrix

  real(DP) :: wgg (nbnd, nbnd, nksq)
  ! input: the weights

  complex(DP) :: ps, aux1 (nbnd), aux2 (nbnd)
  complex(DP), allocatable ::  ps1 (:,:), ps2 (:,:,:), ps3 (:,:), ps4 (:,:,:)
  complex(DP), allocatable ::  ps1_nc(:,:,:), ps2_nc(:,:,:,:), &
                               ps3_nc (:,:,:), ps4_nc (:,:,:,:), &
                               deff_nc(:,:,:,:)
  real(DP), allocatable :: deff(:,:,:)
  ! work space
  complex(DP) ::  dynwrk (3 * nat, 3 * nat), ps_nc(2)
  ! auxiliary dynamical matrix

  integer :: ik, ikk, ikq, ibnd, jbnd, ijkb0, ijkb0b, ih, jh, ikb, &
       jkb, ipol, jpol, startb, lastb, na, nb, nt, ntb, nu_i, nu_j, &
       na_icart, na_jcart, mu, nu, is, js, ijs
  ! counters

  IF (rec_code_read >=-20) return

  IF (noncolin) THEN
     allocate (ps1_nc (  nkb, npol, nbnd))
     allocate (ps2_nc (  nkb, npol, nbnd , 3))
     allocate (ps3_nc (  nkb, npol, nbnd))
     allocate (ps4_nc (  nkb, npol, nbnd , 3))
     allocate (deff_nc (  nhm, nhm, nat, nspin))
  ELSE
     allocate (ps1 (  nkb, nbnd))
     allocate (ps2 (  nkb, nbnd , 3))
     allocate (ps3 (  nkb, nbnd))
     allocate (ps4 (  nkb, nbnd , 3))
     allocate (deff ( nhm, nhm, nat ))
  END IF

  dynwrk (:,:) = (0.d0, 0.d0)
  call divide (intra_bgrp_comm, nbnd, startb, lastb)
  do ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)

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
     do ibnd = 1, nbnd
        IF (noncolin) THEN
           CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
        ELSE
           CALL compute_deff(deff,et(ibnd,ikk))
        ENDIF
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              if (ityp (na) == nt) then
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       IF (noncolin) THEN
                          ijs=0
                          DO is=1,npol
                             DO js=1,npol
                                ijs=ijs+1
                                ps1_nc (ikb, is, ibnd) =      &
                                   ps1_nc (ikb, is, ibnd) +  &
                                   deff_nc(ih,jh,na,ijs)*    &
                                   becp1(ik)%nc (jkb, js, ibnd)
                             END DO
                          END DO
                          IF (lspinorb) THEN
                             ijs=0
                             DO is=1,npol
                                DO js=1,npol
                                   ijs=ijs+1
                                   ps3_nc (ikb, is, ibnd) = &
                                       ps3_nc (ikb, is, ibnd) - &
                                      qq_so(ih,jh,ijs,nt)*becq(ik)%nc(jkb,js,ibnd)
                                END DO
                             END DO
                          ELSE
                             DO is=1,npol
                                ps3_nc(ikb,is,ibnd)=ps3_nc(ikb,is,ibnd) - &
                                  qq (ih, jh, nt) * becq(ik)%nc (jkb, is, ibnd)
                             ENDDO
                          END IF
                       ELSE
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                            deff(ih,jh,na) *                  &
                            becp1(ik)%k (jkb, ibnd)
                          ps3 (ikb, ibnd) = ps3 (ikb, ibnd) - &
                            qq (ih, jh, nt) * becq(ik)%k (jkb, ibnd)
                       END IF
                       do ipol = 1, 3
                          IF (noncolin) THEN
                             ijs=0
                             DO is=1,npol
                                DO js=1,npol
                                   ijs=ijs+1
                                   ps2_nc(ikb,is,ibnd,ipol) =               &
                                       ps2_nc(ikb,is,ibnd,ipol) +           &
                                       deff_nc(ih,jh,na,ijs) *              &
                                       alphap(ipol,ik)%nc(jkb, js, ibnd)+ &
                                       int1_nc(ih, jh, ipol, na, ijs) *     &
                                        becp1(ik)%nc (jkb, js, ibnd)
                                END DO
                             END DO
                             IF (lspinorb) THEN
                                ijs=0
                                DO is=1,npol
                                   DO js=1,npol
                                      ijs=ijs+1
                                      ps4_nc(ikb,is,ibnd,ipol) =          &
                                             ps4_nc(ikb,is,ibnd,ipol)-    &
                                             qq_so(ih,jh,ijs,nt) *        &
                                             alpq(ipol,ik)%nc(jkb,js,ibnd)
                                   END DO
                                END DO
                             ELSE
                                DO is=1,npol
                                   ps4_nc(ikb,is,ibnd,ipol) =                  &
                                     ps4_nc(ikb,is,ibnd,ipol)-              &
                                     qq(ih,jh,nt)*alpq(ipol,ik)%nc(jkb,is,ibnd)
                                END DO
                             END IF
                          ELSE
                             ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                                deff (ih, jh, na) *                          &
                                   alphap(ipol,ik)%k(jkb, ibnd) +            &
                               int1 (ih, jh, ipol, na, current_spin) *       &
                               becp1(ik)%k (jkb, ibnd)
                             ps4 (ikb, ibnd, ipol) = ps4 (ikb, ibnd, ipol) - &
                               qq (ih, jh, nt) * alpq(ipol,ik)%k (jkb,ibnd)
                          END IF
                       enddo  ! ipol
                    enddo
                 enddo
                 ijkb0 = ijkb0 + nh (nt)
              endif
           enddo
        enddo
     END DO
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
                            CONJG(alpq(ipol,ik)%nc(ikb,1,jbnd))*ps1_nc(ikb,1,ibnd)+&
                            CONJG(becq(ik)%nc(ikb,1,jbnd))*ps2_nc(ikb,1,ibnd,ipol)+&
                            CONJG(alpq(ipol,ik)%nc(ikb,2,jbnd))*ps1_nc(ikb,2,ibnd)+&
                            CONJG(becq(ik)%nc(ikb,2,jbnd))*ps2_nc(ikb,2,ibnd,ipol)
                          ELSE
                             aux1 (jbnd) = aux1 (jbnd) + &
                               CONJG(alpq(ipol,ik)%k(ikb,jbnd))*ps1(ikb,ibnd)+&
                               CONJG(becq(ik)%k(ikb,jbnd))*ps2(ikb,ibnd,ipol)
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
                                         ijs=0
                                         DO is=1,npol
                                            DO js=1,npol
                                               ijs=ijs+1
                                               ps_nc(is) = ps_nc(is) + &
                                               int2_so(ih,jh,ipol,na,nb,ijs)*&
                                               becp1(ik)%nc(jkb,js,ibnd)
                                            END DO
                                         END DO
                                      ELSE
                                         DO is=1,npol
                                            ps_nc(is) = ps_nc(is) + &
                                               int2(ih,jh,ipol,na,nb)*&
                                               becp1(ik)%nc(jkb,is,ibnd)
                                         END DO
                                      ENDIF
                                   ELSE
                                      ps = ps + int2 (ih, jh, ipol, na, nb) * &
                                             becp1(ik)%k (jkb, ibnd)
                                   END IF
                                enddo
                                do jbnd = startb, lastb
                                   IF (noncolin) THEN
                                      aux1(jbnd) = aux1 (jbnd) + &
                                        ps_nc(1)*CONJG(becq(ik)%nc(ikb,1,jbnd))+&
                                        ps_nc(2)*CONJG(becq(ik)%nc(ikb,2,jbnd))
                                   ELSE
                                      aux1(jbnd) = aux1 (jbnd) + &
                                        ps * CONJG(becq(ik)%k(ikb,jbnd))
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
                                        (CONJG(alphap(jpol,ik)%nc(ikb,1,ibnd))*&
                                            ps3_nc (ikb, 1, jbnd) + &
                                            CONJG(becp1(ik)%nc (ikb,1,ibnd))* &
                                            ps4_nc (ikb, 1, jbnd, jpol) +  &
                                         CONJG(alphap(jpol,ik)%nc(ikb,2,ibnd))*&
                                            ps3_nc (ikb,2,jbnd) + &
                                            CONJG(becp1(ik)%nc (ikb,2,ibnd)) * &
                                            ps4_nc (ikb, 2, jbnd, jpol) )
                                      ELSE
                                         aux2 (jbnd) = aux2 (jbnd) + &
                                           wgg (ibnd, jbnd, ik) * &
                                           (CONJG(alphap(jpol,ik)%k(ikb,ibnd))*&
                                            ps3 (ikb, jbnd) + &
                                            CONJG(becp1(ik)%k (ikb, ibnd) ) * &
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

  call mp_sum ( dynwrk, intra_bgrp_comm )

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
     deallocate (deff_nc)
  ELSE
     deallocate (ps4)
     deallocate (ps3)
     deallocate (ps2)
     deallocate (ps1)
     deallocate (deff)
  END IF
  return
end subroutine compute_nldyn
