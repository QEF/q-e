!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us_only (ik, mode, uact)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector uact.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  USE uspp_param, ONLY: nh
  use phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, mode
  ! input: the k point
  ! input: the actual perturbation
  complex(DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! the point k+q
  ! counter on G vectors
  ! auxiliary counter on G vectors
  ! counter on atomic types
  ! counter on bands
  ! auxiliary variable for counting
  ! counter on becp functions
  ! counter on becp functions
  ! counter on n index
  ! counter on m index
  ! counter on polarizations

  real(DP), parameter :: eps = 1.d-12

  complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:)
  complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  ! work space

  logical :: ok

  call start_clock ('dvqpsi_us_on')
  if (noncolin) then
     allocate (ps1_nc(nkb , npol, nbnd))    
     allocate (ps2_nc(nkb , npol, nbnd , 3))    
  else
     allocate (ps1 ( nkb , nbnd))    
     allocate (ps2 ( nkb , nbnd , 3))    
  end if
  allocate (aux ( npwx))    
  if (lgamma) then
     ikk = ik
     ikq = ik
  else
     ikk = 2 * ik - 1
     ikq = ikk + 1
  endif
  if (lsda) current_spin = isk (ikk)
  !
  !   we first compute the coefficients of the vectors
  !
  if (noncolin) then
     ps1_nc(:,:,:)   = (0.d0, 0.d0)
     ps2_nc(:,:,:,:) = (0.d0, 0.d0)
  else
     ps1(:,:)   = (0.d0, 0.d0)
     ps2(:,:,:) = (0.d0, 0.d0)
  end if
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           mu = 3 * (na - 1)
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ipol = 1, 3
                    do ibnd = 1, nbnd
                       if ( abs (uact (mu + 1) ) + &
                            abs (uact (mu + 2) ) + &
                            abs (uact (mu + 3) ) > eps) then
                          IF (noncolin) THEN
                             IF (lspinorb) THEN
                                ps1_nc(ikb,1,ibnd)=ps1_nc(ikb,1,ibnd)+    &
                                    ((deeq_nc(ih,jh,na,1) -               &
                                       et(ibnd,ikk)*qq_so(ih,jh,1,nt))*   &
                                       alphap_nc(jkb,1,ibnd, ipol, ik) +  &
                                    (deeq_nc(ih, jh, na, 2)-              &
                                       et(ibnd,ikk)*qq_so(ih,jh,2,nt))*   &
                                       alphap_nc(jkb,2,ibnd,ipol,ik))*    &
                                       uact(mu + ipol) 
                                ps1_nc(ikb,2,ibnd)=ps1_nc(ikb,2,ibnd)  +  &
                                   ((deeq_nc(ih,jh,na,3) -                &
                                       et(ibnd,ikk)*qq_so(ih,jh,3,nt)) *  &
                                       alphap_nc(jkb,1,ibnd,ipol,ik)   +  & 
                                    (deeq_nc(ih,jh,na,4) -                &
                                       et(ibnd,ikk)*qq_so(ih,jh,4,nt)) *  &
                                       alphap_nc(jkb,2,ibnd,ipol,ik))  *  &
                                       uact(mu + ipol) 
                                ps2_nc(ikb,1,ibnd,ipol)=                  &
                                                 ps2_nc(ikb,1,ibnd,ipol)+ &
                                   ((deeq_nc(ih,jh,na,1) -                &
                                       et(ibnd,ikk)*qq_so(ih,jh,1,nt))  * &
                                       becp1_nc(jkb,1,ibnd,ik) +          &
                                    (deeq_nc(ih,jh,na,2) -                &
                                       et(ibnd,ikk)*qq_so(ih,jh,2,nt))  * &
                                       becp1_nc(jkb,2,ibnd,ik) )*         &
                                      (0.d0,-1.d0)*uact(mu+ipol)*tpiba
                                ps2_nc(ikb,2,ibnd,ipol) =                 &
                                                 ps2_nc(ikb,2,ibnd,ipol)+ &
                                   ((deeq_nc(ih,jh,na,3) -                &
                                       et(ibnd,ikk)*qq_so(ih,jh,3,nt))*   &
                                       becp1_nc(jkb,1,ibnd,ik)        +   &
                                    (deeq_nc(ih,jh,na,4) -                &
                                       et(ibnd,ikk)*qq_so(ih,jh,4,nt))*   &
                                       becp1_nc(jkb,2,ibnd,ik) )*         &
                                      (0.d0,-1.d0)*uact(mu+ipol)*tpiba
                             ELSE
                                ps1_nc(ikb,1,ibnd) = ps1_nc(ikb,1,ibnd) +   &
                                      ((deeq_nc(ih,jh,na,1)-                &
                                                et(ibnd,ikk)*qq(ih,jh,nt))* &
                                        alphap_nc(jkb,1,ibnd,ipol,ik) +     &
                                        deeq_nc(ih, jh, na, 2)*             &
                                        alphap_nc(jkb,2,ibnd,ipol,ik) ) *   &
                                        uact(mu + ipol) 
                                ps1_nc(ikb,2,ibnd)=ps1_nc(ikb, 2, ibnd) +   &
                                      (deeq_nc(ih,jh,na,3)*                 &
                                             alphap_nc(jkb,1,ibnd,ipol,ik)+ & 
                                      (deeq_nc(ih,jh,na,4) -                &
                                             et(ibnd,ikk)*qq(ih,jh,nt) )*   &
                                       alphap_nc(jkb,2,ibnd,ipol,ik))*      &
                                                           uact(mu+ipol) 
                                ps2_nc(ikb,1,ibnd,ipol)=                    &
                                             ps2_nc(ikb,1,ibnd,ipol)+       &
                                    ((deeq_nc(ih,jh,na,1) -                 &
                                             et(ibnd,ikk)*qq(ih,jh,nt))*    &
                                             becp1_nc(jkb,1,ibnd,ik) +      &
                                      deeq_nc(ih,jh,na,2)*                  &
                                            becp1_nc(jkb,2,ibnd,ik) )*      &
                                      (0.d0,-1.d0)*uact(mu+ipol)*tpiba
                                ps2_nc(ikb,2,ibnd,ipol) =                   &
                                             ps2_nc(ikb,2,ibnd,ipol)+       &
                                    ((deeq_nc(ih,jh,na,4) -                 &
                                             et(ibnd,ikk)*qq(ih,jh,nt))*    &
                                             becp1_nc(jkb,2,ibnd,ik) +      &
                                      deeq_nc(ih,jh,na,3)*                  &
                                             becp1_nc(jkb,1,ibnd,ik) ) *    &
                                      (0.d0,-1.d0)*uact(mu+ipol)*tpiba
                             ENDIF
                          ELSE
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                               (deeq (ih, jh, na, current_spin) - &
                                et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                alphap(jkb, ibnd, ipol, ik) * uact (mu + ipol)
                             ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) +&
                                 (deeq (ih,jh, na, current_spin) - &
                                  et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                  (0.d0, -1.d0) * becp1 (jkb, ibnd, ik) * &
                                  uact (mu + ipol) * tpiba
                          ENDIF
                          if (okvan) then
                             IF (noncolin) THEN
                                ps1_nc(ikb,1,ibnd)=ps1_nc(ikb,1,ibnd) +  &
                                   (int1_nc(ih,jh,ipol,na,1) *            &
                                       becp1_nc(jkb,1,ibnd,ik) +          &
                                    int1_nc(ih,jh,ipol,na,2) *            &
                                       becp1_nc(jkb,2,ibnd,ik))*uact(mu+ipol)
                                 ps1_nc(ikb,2,ibnd)=ps1_nc(ikb,2,ibnd) +  &
                                   (int1_nc(ih,jh,ipol,na,3) *            &
                                       becp1_nc(jkb,1,ibnd,ik) +          &
                                    int1_nc(ih,jh,ipol,na,4) *            &
                                       becp1_nc(jkb,2,ibnd,ik))*uact(mu+ipol)
                             ELSE
                                ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int1 (ih, jh, ipol,na, current_spin) * &
                                  becp1 (jkb, ibnd, ik) ) * uact (mu +ipol)
                             END IF
                          endif
                       endif
                       if (okvan) then
                          do nb = 1, nat
                             nu = 3 * (nb - 1)
                             IF (noncolin) THEN
                                IF (lspinorb) THEN
                                   ps1_nc(ikb,1,ibnd)=ps1_nc(ikb,1,ibnd) + &
                                   (int2_so(ih,jh,ipol,nb,na,1) *          &
                                         becp1_nc(jkb,1,ibnd,ik)   +       &
                                    int2_so(ih,jh,ipol,nb,na,2) *          &
                                         becp1_nc(jkb,2,ibnd,ik))*uact(nu+ipol)
                                   ps1_nc(ikb,2,ibnd)=ps1_nc(ikb,2,ibnd) + &
                                   (int2_so(ih,jh,ipol,nb,na,3) *          &
                                         becp1_nc(jkb,1,ibnd,ik)   +       &
                                    int2_so(ih,jh,ipol,nb,na,4) *          &
                                         becp1_nc(jkb,2,ibnd,ik))*uact(nu+ipol)
                                ELSE
                                   ps1_nc(ikb,1,ibnd) = ps1_nc(ikb,1,ibnd) + &
                                   (int2(ih,jh,ipol,nb,na) *                 &
                                         becp1_nc(jkb,1,ibnd,ik) )*uact(nu+ipol)
                                   ps1_nc(ikb,2,ibnd) = ps1_nc(ikb,2,ibnd) + &
                                   (int2(ih,jh,ipol,nb,na) *                 &
                                         becp1_nc(jkb,2,ibnd,ik) )*uact(nu+ipol)
                                END IF
                             ELSE
                                ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                    (int2 (ih, jh, ipol, nb, na) * &
                                     becp1 (jkb, ibnd, ik) ) * uact (nu + ipol)
                             END IF
                          enddo
                       endif
                    enddo
                 enddo
              enddo
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  if (nkb.gt.0) then
     if (noncolin) then
        call ZGEMM ('N', 'N', npwq, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     else
        call ZGEMM ('N', 'N', npwq, nbnd, nkb, &
         (1.d0, 0.d0) , vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     end if
  end if
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        IF (noncolin) THEN
           do ibnd = 1, nbnd
              ok = ok.or.(abs (ps2_nc (ikb, 1, ibnd, ipol) ).gt.eps).or. &
                         (abs (ps2_nc (ikb, 2, ibnd, ipol) ).gt.eps)
           end do
        ELSE
           do ibnd = 1, nbnd
              ok = ok.or. (abs (ps2 (ikb, ibnd, ipol) ) .gt.eps)
           enddo
        ENDIF
        if (ok) then
           do ig = 1, npwq
              igg = igkq (ig)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, ikq) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              IF (noncolin) THEN
                 call ZAXPY(npwq,ps2_nc(ikb,1,ibnd,ipol),aux,1,dvpsi(1,ibnd),1)
                 call ZAXPY(npwq,ps2_nc(ikb,2,ibnd,ipol),aux,1, &
                                                         dvpsi(1+npwx,ibnd),1)
              ELSE
                 call ZAXPY (npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1)
              END IF
           enddo
        endif
     enddo

  enddo
  deallocate (aux)
  IF (noncolin) THEN
     deallocate (ps2_nc)
     deallocate (ps1_nc)
  ELSE
     deallocate (ps2)
     deallocate (ps1)
  END IF

  call stop_clock ('dvqpsi_us_on')
  return
end subroutine dvqpsi_us_only
