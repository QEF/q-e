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
  !! This routine computes the term of the dynamical matrix due to
  !! the orthogonality constraint. Only the part which is due to
  !! the nonlocal terms is computed here.
  !
  USE kinds,            ONLY : DP
  USE klist,            ONLY : wk
  USE lsda_mod,         ONLY : lsda, current_spin, isk, nspin
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module, ONLY : noncolin, npol, lspinorb
  USE uspp,             ONLY : nkb, qq_nt, qq_so, ofsbeta
  USE uspp_param,       ONLY : nh, nhm
  USE wvfct,            ONLY : nbnd, et

  USE modes,            ONLY : u
  USE phus,             ONLY : alphap, int1, int2, &
                               int2_so, int1_nc
  USE control_ph,       ONLY : rec_code_read

  USE lrus,             ONLY : becp1
  USE qpoint,           ONLY : nksq, ikks, ikqs
  USE control_lr,       ONLY : nbnd_occ

  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE becmod,           ONLY : bec_type

  implicit none

  type (bec_type) :: becq(nksq)
  !! input: the becp with \(\text{psi}_{k+q}\)
  type (bec_type) :: alpq(3,nksq)
  !! input: the alphap with \(\text{psi}_{k}\)
  complex(DP) :: wdyn(3*nat,3*nat)
  !! output: the term of the dynamical matrix
  real(DP) :: wgg(nbnd,nbnd,nksq)
  !! input: the weights
  
  ! ... local variables
  
  complex(DP) :: ps 
  complex(DP),parameter  :: ONE = cmplx(1._DP, 0._DP, kind=DP)
  complex(DP), allocatable ::  aux1 (:), aux2 (:)
  complex(DP), allocatable ::  ps1 (:,:), ps2 (:,:,:), ps3 (:,:), ps4 (:,:,:)
  complex(DP), allocatable ::  ps1_nc(:,:,:), ps2_nc(:,:,:,:), &
                               ps3_nc (:,:,:), ps4_nc (:,:,:,:), &
                               deff_nc(:,:,:,:), auxdyn(:,:), psv(:), psv_nc(:,:)
  real(DP), allocatable :: deff(:,:,:)
  ! work space
  complex(DP) ::  dynwrk (3 * nat, 3 * nat), ps_nc(2)
  ! auxiliary dynamical matrix

  integer :: ik, ikk, ikq, ibnd, jbnd, ijkb0, ijkbnh, ijkb0b, ijkbnhb, ih, jh, ikb, &
       jkb, ipol, jpol, startb, lastb, na, nb, nt, ntb, nu_i, nu_j, &
       na_icart, na_jcart, mu, nu, is, js, ijs, jbnd_loc, ijkb1b, ijkbnb
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
  allocate (psv_nc(npol,nhm), psv(nhm)) 
  dynwrk (:,:) = (0.d0, 0.d0)
  call divide (intra_bgrp_comm, nbnd, startb, lastb)
  allocate(aux1(lastb - startb + 1), aux2(lastb - startb + 1)) 
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
     call start_clock('c_nldynb1')
     do ibnd = 1, nbnd
        IF (noncolin) THEN
           CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
        ELSE
           CALL compute_deff(deff,et(ibnd,ikk))
        ENDIF
        call start_clock('c_nldyn_loop1')
        do nt = 1, ntyp
           do na = 1, nat
              if (ityp (na) == nt) then
                 ijkb0 = ofsbeta(na) 
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
                                  qq_nt (ih, jh, nt) * becq(ik)%nc (jkb, is, ibnd)
                             ENDDO
                          END IF
                       ELSE
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                            deff(ih,jh,na) *                  &
                            becp1(ik)%k (jkb, ibnd)
                          ps3 (ikb, ibnd) = ps3 (ikb, ibnd) - &
                            qq_nt (ih, jh, nt) * becq(ik)%k (jkb, ibnd)
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
                                     qq_nt(ih,jh,nt)*alpq(ipol,ik)%nc(jkb,is,ibnd)
                                END DO
                             END IF
                          ELSE
                             ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                                deff (ih, jh, na) *                          &
                                   alphap(ipol,ik)%k(jkb, ibnd) +            &
                               int1 (ih, jh, ipol, na, current_spin) *       &
                               becp1(ik)%k (jkb, ibnd)
                             ps4 (ikb, ibnd, ipol) = ps4 (ikb, ibnd, ipol) - &
                               qq_nt (ih, jh, nt) * alpq(ipol,ik)%k (jkb,ibnd)
                          END IF
                       enddo  ! ipol
                    enddo
                 enddo
              endif
           enddo
        enddo
        call stop_clock('c_nldyn_loop1')
     END DO
     call stop_clock('c_nldynb1') 
     !
     !     Here starts the loop on the atoms (rows)
     ! 
     call start_clock('c_nldynb2')
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              ijkb0 = ofsbeta(na)
              ijkbnh = ijkb0 + nh(nt)
              do ipol = 1, 3
                 mu = 3 * (na - 1) + ipol
                 do ibnd = 1, nbnd_occ (ikk)
                     aux1 (:) = (0.d0, 0.d0)
                     do jbnd = startb, lastb
                        jbnd_loc = jbnd - startb + 1 
                        IF (noncolin) THEN
                           aux1 (jbnd_loc) = aux1 (jbnd_loc) + &
                           dot_product(alpq(ipol,ik)%nc(ijkb0 + 1 :ijkbnh,1,jbnd),ps1_nc(ijkb0 + 1 :ijkbnh,1, ibnd))+& 
                           dot_product(becq(ik)%nc(ijkb0 + 1 :ijkbnh,1,jbnd),ps2_nc(ijkb0 + 1 :ijkbnh,1,ibnd,ipol))+& 
                           dot_product(alpq(ipol,ik)%nc(ijkb0 + 1 :ijkbnh,2,jbnd),ps1_nc(ijkb0 + 1 :ijkbnh,2, ibnd))+&
                           dot_product(becq(ik)%nc(ijkb0 + 1 :ijkbnh,1,jbnd),ps2_nc(ijkb0 + 1:ijkbnh,1,ibnd,ipol))
                        ELSE
                           aux1 (jbnd_loc) = aux1 (jbnd_loc) + &
                           dot_product(alpq(ipol,ik)%k(ijkb0 + 1:ijkbnh,jbnd),ps1(ijkb0 + 1:ijkbnh,ibnd))+&
                           dot_product(becq(ik)%k(ijkb0 + 1:ijkbnh,jbnd),ps2(ijkb0 + 1:ijkbnh,ibnd,ipol))
                        END IF
                     enddo
                     call start_clock('nldyn_b21')
                     do ntb = 1, ntyp
                        do nb = 1, nat
                          if (ityp (nb) == ntb) then
                             psv_nc(:, :) =(0.d0,0.d0)
                             psv(:) = (0.d0, 0.d0)
                             ijkb0b = ofsbeta(nb)
                             IF (noncolin) THEN
                                 IF (lspinorb) THEN
                                    ijs=0
                                    DO is=1,npol
                                       DO js=1,npol
                                          ijs=ijs+1
                                          call zgemv('N', nh(ntb), nh(ntb), ONE, int2_so(1,1,ipol,na,nb,ijs), nhm, & 
                                                      becp1(ik)%nc(ijkb0b +1, js, ibnd), 1, ONE, psv_nc(is, 1), npol)  
                                       END DO
                                    END DO
                                 ELSE
                                    DO is=1,npol
                                       call zgemv('N', nh(ntb), nh(ntb), ONE, int2(1,1,ipol,na,nb), nhm, & 
                                                       becp1(ik)%nc(ijkb0b + 1, js, ibnd), 1, & 
                                                       ONE, psv_nc(is,1), npol) 
                                    END DO
                                 ENDIF
                             ELSE
                                 call zgemv('N', nh(ntb), nh(ntb), ONE, int2(1,1,ipol, na,nb), nhm, & 
                                                 becp1(ik)%k(ijkb0b + 1, ibnd),1 , ONE,  psv(1), 1) 
                             END IF
                             do jbnd = startb, lastb
                                 jbnd_loc = jbnd - startb + 1 
                                 IF (noncolin) THEN
                                    aux1(jbnd_loc) = aux1 (jbnd_loc) + &
                                    dot_product(becq(ik)%nc(ijkb0b + 1: ijkb0b + nh(ntb), 1, jbnd), & 
                                                    psv_nc(1,1:nh(ntb))) + & 
                                    dot_product(becq(ik)%nc(ijkb0b + 1: ijkb0b + nh(ntb), 2, jbnd), & 
                                                    psv_nc(2,1:nh(ntb)))   
                                 ELSE
                                    aux1(jbnd_loc) = aux1 (jbnd_loc) + &
                                    dot_product(becq(ik)%k(ijkb0b + 1: ijkb0b + nh(ntb),jbnd), & 
                                                   psv(1:nh(ntb)))
                                 END IF
                              enddo
                          endif
                       enddo
                    enddo
                    call stop_clock('nldyn_b21')
                    !
                    !     here starts the second loop on the atoms
                    !
                    call start_clock("nldyn_b22")
                    do ntb = 1, ntyp
                       do nb = 1, nat
                          if (ityp (nb) == ntb) then
                             ijkb0b = ofsbeta(nb)
                             ijkb1b = ijkb0b + 1 
                             ijkbnb = ijkb0b + nh(ntb)
                             do jpol = 1, 3
                                nu = 3 * (nb - 1) + jpol
                                aux2 (:) = (0.d0, 0.d0)
                                do jbnd = startb, lastb
                                    jbnd_loc = jbnd - startb + 1 
                                    IF (noncolin) THEN 
                                       aux2 (jbnd_loc) = aux2 (jbnd_loc) + &
                                       wgg(ibnd, jbnd, ik) * ( &
                                       dot_product(alphap(jpol,ik)%nc(ijkb1b:ijkbnb, 1 ,ibnd), &  
                                                       ps3_nc(ijkb1b:ijkbnb,1, jbnd)) + & 
                                       dot_product(becp1(ik)%nc(ijkb1b: ijkbnb, 1, ibnd), &
                                                        ps4_nc(ijkb1b:ijkbnb, 1, jbnd, jpol)) + & 
                                       dot_product(alphap(jpol,ik)%nc(ijkb1b:ijkbnb, 2, ibnd), & 
                                                       ps3_nc(ijkb1b:ijkbnb, 2, jbnd)) + &  
                                       dot_product(becp1(ik)%nc(ijkb1b:ijkbnb, 2, ibnd), & 
                                                      ps4_nc(ijkb1b:ijkbnb, 2, jbnd, jpol)))  
                                      ELSE
                                         aux2 (jbnd_loc) = aux2 (jbnd_loc) + &
                                           wgg (ibnd, jbnd, ik) * ( & 
                                           dot_product(alphap(jpol,ik)%k(ijkb1b:ijkbnb, ibnd), & 
                                                       ps3(ijkb1b:ijkbnb, jbnd)) + & 
                                           dot_product(becp1(ik)%k(ijkb1b:ijkbnb, ibnd), & 
                                                      ps4(ijkb1b:ijkbnb, jbnd, jpol)) & 
                                           )
                                      END IF
                                   enddo
                                do jbnd = startb, lastb
                                   jbnd_loc = jbnd - startb + 1 
                                   dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0*wk(ikk) * aux2(jbnd_loc) * aux1(jbnd_loc)
                                enddo
                             enddo
                          endif
                       enddo
                    enddo
                    call stop_clock('nldyn_b22')
                 enddo
              enddo
           endif
        enddo
     enddo
  enddo
  call stop_clock('c_nldynb2')
  call mp_sum ( dynwrk, intra_bgrp_comm )
  allocate(auxdyn(3*nat,3*nat))
  call zgemm('C', 'N', 3*nat, 3*nat, 3*nat, cmplx(1.d0,0.d0,kind=dp),  u,     3*nat, dynwrk, 3*nat,& 
                                            cmplx(0.d0,0.d0,kind=dp), auxdyn, 3*nat) 
  call zgemm('N','N', 3*nat, 3* nat, 3*nat, cmplx(1.d0,0.d0,kind=dp), auxdyn, 3*nat, u     , 3*nat,& 
                                            cmplx(1.d0,0.d0,kind=dp),  wdyn , 3*nat)

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
