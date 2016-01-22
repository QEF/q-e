!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drhodvnl (ik, ikk, nper, nu_i0, wdyn, dbecq, dalpq)
  !-----------------------------------------------------------------------
  !
  !  This subroutine computes the electronic term 2 <dpsi|dv-e ds|psi> of 
  !  the dynamical matrix. It can be used both for KB and for US 
  !  pseudopotentials. All the nonlocal (and overlap matrix) terms
  !  are computed here. The contribution of the local potential is not
  !  computed here. This routine must be called for each k point and
  !  accumulates in wdyn the contribution of each k point.
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp,      ONLY : okvan, nkb
  USE uspp_param,ONLY : nh, nhm
  USE becmod,    ONLY : bec_type
  USE wvfct,     ONLY : nbnd, et
  USE klist,     ONLY : wk
  USE lsda_mod,  ONLY : current_spin, nspin
  USE spin_orb,  ONLY : lspinorb
  USE phus,      ONLY : int1, int1_nc, int2, int2_so, alphap

  USE lrus,      ONLY : becp1

  USE mp_bands, ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  implicit none
  integer :: ik, ikk, nper, nu_i0
  ! input: the current k point
  ! input: the number of perturbations
  ! input: the initial mode

  TYPE(bec_type) :: dbecq(nper), dalpq(3,nper)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k}
  complex(DP) :: wdyn (3 * nat, 3 * nat)
  ! output: the term of the dynamical matryx

  complex(DP) :: ps, ps_nc(npol), dynwrk (3 * nat, 3 * nat)
  ! dynamical matrix
  complex(DP) , allocatable :: ps1 (:,:), ps2 (:,:,:)
  complex(DP) , allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:), &
                               deff_nc(:,:,:,:)
  real(DP), allocatable :: deff(:,:,:)

  integer :: ibnd, ijkb0, ijkb0b, ih, jh, ikb, jkb, ipol, &
       startb, lastb, iper, na, nb, nt, ntb, mu, nu, is, js, ijs
  ! counters

  IF (noncolin) THEN
     allocate (ps1_nc ( nkb, npol, nbnd))
     allocate (ps2_nc ( nkb, npol, nbnd, 3))
     allocate (deff_nc ( nhm, nhm, nat, nspin ))
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
  ELSE
     allocate (ps1 (  nkb , nbnd))
     allocate (ps2 (  nkb , nbnd , 3))
     allocate (deff ( nhm, nhm, nat ))
     ps1 = (0.d0, 0.d0)
     ps2 = (0.d0, 0.d0)
  END IF

  dynwrk (:, :) = (0.d0, 0.d0)

  call divide (intra_bgrp_comm, nbnd, startb, lastb)
  !
  !   Here we prepare the two terms
  !
  do ibnd = startb, lastb
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
                       DO is=1, npol
                          DO js=1, npol
                             ijs=ijs+1
                             ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd) + &
                               deff_nc(ih,jh,na,ijs) * becp1(ik)%nc(jkb,js,ibnd)
                          END DO
                       END DO
                    ELSE
                       ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                         deff(ih,jh,na)*becp1(ik)%k(jkb,ibnd)
                    END IF
                    do ipol = 1, 3
                       IF (noncolin) THEN
                          ijs=0
                          DO is=1, npol
                             DO js=1, npol
                                ijs=ijs+1
                                ps2_nc(ikb,is,ibnd,ipol) =         &
                                       ps2_nc(ikb,is,ibnd,ipol)+   &
                                       deff_nc(ih,jh,na,ijs)   *   &
                                       alphap(ipol,ik)%nc(jkb,js,ibnd)
                             END DO
                          END DO
                       ELSE
                          ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                          deff(ih,jh,na) * alphap(ipol,ik)%k(jkb,ibnd)
                       END IF

                       IF (okvan) THEN
                          IF (noncolin) THEN
                             ijs=0
                             DO is=1, npol
                                DO js=1, npol
                                   ijs=ijs+1
                                   ps2_nc (ikb, is, ibnd, ipol) =        &
                                        ps2_nc (ikb, is, ibnd, ipol)   + &
                                        int1_nc(ih, jh, ipol, na, ijs) * &
                                        becp1(ik)%nc (jkb, js, ibnd)
                                END DO
                             END DO
                          ELSE
                             ps2 (ikb, ibnd, ipol) = &
                                   ps2 (ikb, ibnd, ipol) + &
                                   int1 (ih, jh, ipol, na, current_spin) * &
                                   becp1(ik)%k (jkb, ibnd)
                          END IF
                       END IF
                    enddo  ! ipol
                 enddo  ! jh
              enddo  ! ih
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo  ! na
     enddo  ! nt
  enddo  ! nbnd
  !
  !     Here starts the loop on the atoms (rows)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) == nt) then
           do ipol = 1, 3
              mu = 3 * (na - 1) + ipol
              do ibnd = startb, lastb
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    do iper = 1, nper
                       nu = nu_i0 + iper
                       IF (noncolin) THEN
                          DO is=1, npol
                             dynwrk (nu, mu) = dynwrk (nu, mu) +2.d0*wk(ikk)* &
                                (ps2_nc(ikb,is,ibnd,ipol)*       &
                                CONJG(dbecq(iper)%nc(ikb,is,ibnd))+  &
                                ps1_nc(ikb,is,ibnd)*CONJG(       &
                                dalpq(ipol,iper)%nc(ikb,is,ibnd)) )
                          END DO
                       ELSE
                          dynwrk (nu, mu) = dynwrk (nu, mu) + &
                            2.d0 * wk (ikk) * (ps2 (ikb, ibnd, ipol) * &
                            CONJG(dbecq(iper)%k(ikb, ibnd) ) + &
                            ps1(ikb,ibnd) * CONJG(dalpq(ipol,iper)%k(ikb,ibnd)))
                       END IF
                    enddo
                 enddo
                 if (okvan) then
                    ijkb0b = 0
                    do ntb = 1, ntyp
                       do nb = 1, nat
                          if (ityp (nb) == ntb) then
                             do ih = 1, nh (ntb)
                                ikb = ijkb0b + ih
                                IF (noncolin) THEN
                                   ps_nc = (0.d0, 0.d0)
                                ELSE
                                   ps = (0.d0, 0.d0)
                                END IF
                                do jh = 1, nh (ntb)
                                   jkb = ijkb0b + jh
                                   IF (noncolin) THEN
                                      IF (lspinorb) THEN
                                         ijs=0
                                         DO is=1, npol
                                            DO js=1, npol
                                               ijs=ijs+1
                                               ps_nc(is)=ps_nc(is)+ &
                                                int2_so(ih,jh,ipol,na,nb,ijs)*&
                                              becp1(ik)%nc(jkb, js, ibnd)
                                            END DO
                                         END DO
                                      ELSE
                                         DO is=1, npol
                                            ps_nc(is)=ps_nc(is)+ &
                                               int2(ih,jh,ipol,na,nb)*&
                                              becp1(ik)%nc(jkb, is, ibnd)
                                         END DO
                                      END IF
                                   ELSE
                                      ps = ps + int2 (ih, jh, ipol, na, nb) * &
                                        becp1(ik)%k (jkb, ibnd)
                                   ENDIF
                                enddo
                                do iper = 1, nper
                                   nu = nu_i0 + iper
                                   IF (noncolin) THEN
                                      DO is=1, npol
                                         dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0 * wk (ikk) * ps_nc(is) * &
                                        CONJG(dbecq(iper)%nc(ikb, is, ibnd))
                                      END DO
                                   ELSE
                                      dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0 * wk (ikk) * ps * &
                                        CONJG(dbecq(iper)%k(ikb,ibnd) )
                                   END IF
                                enddo
                             enddo
                             ijkb0b = ijkb0b + nh (ntb)
                          endif
                       enddo
                    enddo
                 endif
              enddo
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  call mp_sum ( dynwrk, intra_bgrp_comm )
  wdyn (:,:) = wdyn (:,:) + dynwrk (:,:)

  IF (noncolin) THEN
     deallocate (ps2_nc)
     deallocate (ps1_nc)
     deallocate (deff_nc)
  ELSE
     deallocate (ps2)
     deallocate (deff)
  END IF
  return
end subroutine drhodvnl
