!
! Copyright (C) 2001 PWSCF group
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
  !  This routine compute the term of the dynamical matrix due to
  !  the orthogonality constraint. Only the part which is due to
  !  the nonlocal terms is computed here
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ntyp => nsp, ityp 
  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  USE uspp_param, only: nh
  use phcom
  implicit none
  integer :: ik, ikk, nper, nu_i0
  ! input: the current k point
  ! input: the number of perturbations
  ! input: the initial mode

  complex(DP) :: dbecq(nkb,npol,nbnd,nper), dalpq(nkb,npol,nbnd,3,nper),&
          wdyn (3 * nat, 3 * nat)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k}
  ! output: the term of the dynamical matryx

  complex(DP) :: ps, ps_nc(npol), dynwrk (3 * nat, 3 * nat)
  ! dynamical matrix
  complex(DP) , allocatable :: ps1 (:,:), ps2 (:,:,:)
  complex(DP) , allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)

  integer :: ibnd, ijkb0, ijkb0b, ih, jh, ikb, jkb, ipol, jpol, &
       startb, lastb, iper, na, nb, nt, ntb, mu, nu
  ! counters

  IF (noncolin) THEN
     allocate (ps1_nc ( nkb, npol, nbnd))
     allocate (ps2_nc ( nkb, npol, nbnd, 3))    
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
  ELSE
     allocate (ps1 (  nkb , nbnd))
     allocate (ps2 (  nkb , nbnd , 3))    
     ps1 = (0.d0, 0.d0)
     ps2 = (0.d0, 0.d0)
  END IF

  dynwrk (:, :) = (0.d0, 0.d0)

  call divide (nbnd, startb, lastb)
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
                 do ibnd = startb, lastb
                    IF (noncolin) THEN 
                       IF (lspinorb) THEN
                          ps1_nc (ikb, 1, ibnd) = ps1_nc (ikb, 1, ibnd) + &
                            (deeq_nc(ih,jh,na,1) -                        &
                                   et(ibnd,ikk)*qq_so(ih,jh,1,nt))*       &
                              becp1_nc(jkb,1,ibnd,ik)                     &
                            + (deeq_nc(ih,jh,na,2) -                      &
                                   et(ibnd,ikk)*qq_so(ih,jh,2,nt)) *      &
                              becp1_nc(jkb,2,ibnd,ik)
                          ps1_nc(ikb,2,ibnd)=ps1_nc(ikb,2,ibnd) + &
                             (deeq_nc(ih,jh,na,3) -                       &
                                    et(ibnd,ikk)*qq_so(ih,jh,3,nt))       &
                            *  becp1_nc(jkb,1,ibnd,ik)                    &
                            + (deeq_nc(ih,jh,na,4) -                      &
                                  et(ibnd,ikk)*qq_so(ih,jh,4,nt)) *       &
                                 becp1_nc(jkb,2,ibnd,ik)
                       ELSE
                          ps1_nc (ikb, 1, ibnd) = ps1_nc (ikb, 1, ibnd) + &
                         (deeq_nc(ih,jh,na,1) - &
                          et (ibnd, ikk)*qq(ih,jh,nt))*becp1_nc(jkb,1,ibnd,ik)&
                          +deeq_nc(ih,jh,na,2)*becp1_nc(jkb,2,ibnd,ik)
                          ps1_nc (ikb, 2, ibnd) = ps1_nc (ikb, 2, ibnd) + &
                         ( deeq_nc(ih,jh,na,4) - &
                           et(ibnd,ikk)*qq(ih,jh,nt))*becp1_nc(jkb,2,ibnd,ik)&
                         + deeq_nc(ih,jh,na,3)*becp1_nc(jkb,1,ibnd,ik)
                       END IF
                    ELSE
                       ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                           (deeq(ih, jh, na,current_spin) - &
                           et(ibnd,ikk)*qq(ih,jh,nt))*becp1(jkb,ibnd,ik)
                    END IF
                    do ipol = 1, 3
                       IF (noncolin) THEN
                          IF (lspinorb) THEN
                             ps2_nc(ikb,1,ibnd,ipol)=ps2_nc(ikb,1,ibnd,ipol)+ &
                                (deeq_nc(ih, jh, na, 1) -                   &
                                  et (ibnd, ikk)*qq_so(ih, jh, 1, nt) )*    &
                                  alphap_nc(jkb, 1, ibnd, ipol, ik) +       &
                                (deeq_nc(ih, jh, na, 2) -                   &
                                  et (ibnd, ikk)*qq_so(ih, jh, 2, nt) )*    &
                                  alphap_nc(jkb, 2, ibnd, ipol, ik)
                             ps2_nc(ikb,2,ibnd,ipol)=ps2_nc(ikb,2,ibnd,ipol)+ &
                                (deeq_nc(ih, jh, na, 3) -                   &
                                   et(ibnd,ikk)*qq_so(ih,jh,3,nt))*         &
                                   alphap_nc(jkb, 1, ibnd, ipol, ik) +      &
                                (deeq_nc(ih,jh,na,4) -                      &
                                   et (ibnd,ikk)*qq_so(ih,jh,4,nt) )*      & 
                                   alphap_nc(jkb, 2, ibnd, ipol, ik)
                          ELSE
                             ps2_nc(ikb,1,ibnd,ipol)=ps2_nc(ikb,1,ibnd,ipol)+ &
                                (deeq_nc(ih, jh, na, 1) -                   &
                                  et (ibnd, ikk)*qq(ih, jh, nt) )*    &
                                  alphap_nc(jkb, 1, ibnd, ipol, ik) +       &
                                 deeq_nc(ih, jh, na, 2) *                   &
                                  alphap_nc(jkb, 2, ibnd, ipol, ik)
                             ps2_nc(ikb,2,ibnd,ipol)=ps2_nc(ikb,2,ibnd,ipol)+ &
                                 deeq_nc(ih, jh, na, 3) *                   &
                                   alphap_nc(jkb, 1, ibnd, ipol, ik) +      &
                                (deeq_nc(ih,jh,na,4) -                      &
                                   et (ibnd,ikk)*qq(ih,jh,nt) )*            & 
                                   alphap_nc(jkb, 2, ibnd, ipol, ik)
                          ENDIF
                       ELSE
                          ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                            (deeq (ih, jh,na, current_spin) - &
                             et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                            alphap (jkb, ibnd, ipol, ik)
                       END IF

                       IF (okvan) THEN
                          IF (noncolin) THEN
                             ps2_nc (ikb, 1, ibnd, ipol) = &
                                 ps2_nc (ikb, 1, ibnd, ipol)  + &
                                 int1_nc(ih, jh, ipol, na, 1) * &
                                 becp1_nc (jkb, 1, ibnd, ik)  + &
                                 int1_nc(ih, jh, ipol, na, 2) * &
                                 becp1_nc (jkb, 2, ibnd, ik)
                             ps2_nc (ikb, 2, ibnd, ipol) = &
                                 ps2_nc (ikb, 2, ibnd, ipol)  + &
                                 int1_nc(ih, jh, ipol, na, 3) * &
                                 becp1_nc (jkb, 1, ibnd, ik)  + &
                                 int1_nc(ih, jh, ipol, na, 4) * &
                                 becp1_nc (jkb, 2, ibnd, ik)
                          ELSE
                             ps2 (ikb, ibnd, ipol) = &
                                   ps2 (ikb, ibnd, ipol) + &
                                   int1 (ih, jh, ipol, na, current_spin) * &
                                   becp1 (jkb, ibnd, ik)
                          END IF
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
        if (ityp (na) == nt) then
           do ipol = 1, 3
              mu = 3 * (na - 1) + ipol
              do ibnd = startb, lastb
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    do iper = 1, nper
                       nu = nu_i0 + iper
                       IF (noncolin) THEN
                          dynwrk (nu, mu) = dynwrk (nu, mu) +2.d0*wk(ikk)* &
                            (ps2_nc(ikb,1,ibnd,ipol)*       &
                            CONJG(dbecq(ikb,1,ibnd,iper))+  &
                            ps1_nc(ikb,1,ibnd)*CONJG(       &
                            dalpq(ikb,1,ibnd,ipol,iper)) +  &
                            ps2_nc(ikb,2,ibnd,ipol)*        &
                            CONJG(dbecq(ikb,2,ibnd,iper))+  &
                            ps1_nc(ikb,2,ibnd)*CONJG(       &
                            dalpq(ikb,2,ibnd,ipol,iper)) )
                       ELSE
                          dynwrk (nu, mu) = dynwrk (nu, mu) + &
                            2.d0 * wk (ikk) * (ps2 (ikb, ibnd, ipol) * &
                            CONJG(dbecq (ikb, 1, ibnd, iper) ) + &
                            ps1(ikb,ibnd) * CONJG(dalpq(ikb,1,ibnd,ipol,iper)) )
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
                                         ps_nc(1)=ps_nc(1)+ &
                                           int2_so(ih,jh,ipol,na,nb,1)*&
                                              becp1_nc(jkb, 1, ibnd,ik) + &
                                           int2_so(ih,jh,ipol,na,nb,2)*&
                                              becp1_nc(jkb, 2, ibnd,ik)
                                         ps_nc(2)=ps_nc(2)+ &
                                           int2_so(ih,jh,ipol,na,nb,3)*&
                                              becp1_nc(jkb, 1, ibnd,ik) + &
                                           int2_so(ih,jh,ipol,na,nb,4)*&
                                              becp1_nc(jkb, 2, ibnd,ik)
                                      ELSE
                                         ps_nc(1)=ps_nc(1)+ &
                                           int2(ih,jh,ipol,na,nb)*&
                                              becp1_nc(jkb, 1, ibnd,ik) 
                                         ps_nc(2)=ps_nc(2)+ &
                                           int2(ih,jh,ipol,na,nb)*&
                                              becp1_nc(jkb, 2, ibnd,ik) 
                                      END IF
                                   ELSE 
                                      ps = ps + int2 (ih, jh, ipol, na, nb) * &
                                        becp1 (jkb, ibnd,ik)
                                   ENDIF
                                enddo
                                do iper = 1, nper
                                   nu = nu_i0 + iper
                                   IF (noncolin) THEN
                                      dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0 * wk (ikk) * (ps_nc(1) * &
                                        CONJG(dbecq (ikb, 1, ibnd, iper)) + &
                                         ps_nc(2) * &
                                        CONJG(dbecq (ikb, 2, ibnd, iper) ) )
                                   ELSE
                                      dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0 * wk (ikk) * ps * &
                                        CONJG(dbecq (ikb, 1, ibnd, iper) )
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
#ifdef __PARA
  call reduce (2 * 3 * nat * 3 * nat, dynwrk)
#endif
  wdyn (:,:) = wdyn (:,:) + dynwrk (:,:)

  IF (noncolin) THEN
     deallocate (ps2_nc)
     deallocate (ps1_nc)
  ELSE
     deallocate (ps2)
     deallocate (ps1)
  END IF
  return
end subroutine drhodvnl
