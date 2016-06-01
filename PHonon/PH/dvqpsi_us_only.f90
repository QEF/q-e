!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us_only (ik, uact)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector uact.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  ! This routine implements Eq. B29 of PRB 64, 235118 (2001).
  ! Only the contribution of the nonlocal potential is calculated here.
  !
  !
  USE kinds, only : DP
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : g
  USE klist,     ONLY : xk, ngk, igk_k
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,  ONLY : lsda, current_spin, isk, nspin
  USE spin_orb,  ONLY : lspinorb
  USE wvfct,     ONLY : nbnd, npwx, et
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp, ONLY: okvan, nkb, vkb
  USE uspp_param, ONLY: nh, nhm
  USE phus,      ONLY : int1, int1_nc, int2, int2_so, alphap

  USE lrus,       ONLY : becp1
  USE qpoint,     ONLY : ikks, ikqs
  USE eqv,        ONLY : dvpsi
  USE control_lr, ONLY : lgamma

  implicit none
  !
  !   The dummy variables
  !

  integer :: ik
  ! input: the k point
  complex(DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol, is, js, ijs, npwq
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

  complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:), deff_nc(:,:,:,:)
  real(DP), allocatable :: deff(:,:,:)
  complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  ! work space

  logical :: ok

  call start_clock ('dvqpsi_us_on')
  if (noncolin) then
     allocate (ps1_nc(nkb , npol, nbnd))
     allocate (ps2_nc(nkb , npol, nbnd , 3))
     allocate (deff_nc(nhm, nhm, nat, nspin))
  else
     allocate (ps1 ( nkb , nbnd))
     allocate (ps2 ( nkb , nbnd , 3))
     allocate (deff(nhm, nhm, nat))
  end if
  allocate (aux ( npwx))
  ikk = ikks(ik)
  ikq = ikqs(ik)
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
  do ibnd = 1, nbnd
     IF (noncolin) THEN
        CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
     ELSE
        CALL compute_deff(deff,et(ibnd,ikk))
     ENDIF
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
                       if ( abs (uact (mu + 1) ) + &
                            abs (uact (mu + 2) ) + &
                            abs (uact (mu + 3) ) > eps) then
                          IF (noncolin) THEN
                             ijs=0
                             DO is=1,npol
                                DO js=1,npol
                                   ijs=ijs+1
                                   ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd) +  &
                                      deff_nc(ih,jh,na,ijs) * &
                                      alphap(ipol, ik)%nc(jkb,js,ibnd)* &
                                       uact(mu + ipol)
                                   ps2_nc(ikb,is,ibnd,ipol)=               &
                                          ps2_nc(ikb,is,ibnd,ipol)+        &
                                          deff_nc(ih,jh,na,ijs) *          &
                                          becp1(ik)%nc(jkb,js,ibnd) *      &
                                          (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                                END DO
                             END DO
                          ELSE
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +      &
                                        deff(ih, jh, na) *            &
                                alphap(ipol, ik)%k(jkb, ibnd) * uact (mu + ipol)
                             ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) +&
                                  deff(ih,jh,na)*becp1(ik)%k (jkb, ibnd) *  &
                                  (0.0_DP,-1.0_DP) * uact (mu + ipol) * tpiba
                          ENDIF
                          IF (okvan) THEN
                             IF (noncolin) THEN
                                ijs=0
                                DO is=1,npol
                                   DO js=1,npol
                                      ijs=ijs+1
                                      ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd)+ &
                                         int1_nc(ih,jh,ipol,na,ijs) *     &
                                         becp1(ik)%nc(jkb,js,ibnd)*uact(mu+ipol)
                                   END DO
                                END DO
                             ELSE
                                ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int1 (ih, jh, ipol,na, current_spin) * &
                                  becp1(ik)%k (jkb, ibnd) ) * uact (mu +ipol)
                             END IF
                          END IF
                       END IF  ! uact>0
                       if (okvan) then
                          do nb = 1, nat
                             nu = 3 * (nb - 1)
                             IF (noncolin) THEN
                                IF (lspinorb) THEN
                                   ijs=0
                                   DO is=1,npol
                                      DO js=1,npol
                                         ijs=ijs+1
                                         ps1_nc(ikb,is,ibnd)= &
                                                   ps1_nc(ikb,is,ibnd)+ &
                                         int2_so(ih,jh,ipol,nb,na,ijs)* &
                                          becp1(ik)%nc(jkb,js,ibnd)*uact(nu+ipol)
                                      END DO
                                   END DO
                                ELSE
                                   DO is=1,npol
                                      ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd)+ &
                                         int2(ih,jh,ipol,nb,na) * &
                                         becp1(ik)%nc(jkb,is,ibnd)*uact(nu+ipol)
                                   END DO
                                END IF
                             ELSE
                                ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                    (int2 (ih, jh, ipol, nb, na) * &
                                     becp1(ik)%k (jkb, ibnd) ) * uact (nu + ipol)
                             END IF
                          enddo
                       endif  ! okvan
                    enddo ! ipol
                 enddo ! jh
              enddo ! ih
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo  ! na
     enddo ! nt
  enddo ! nbnd
  !
  !      This term is proportional to beta(k+q+G)
  ! 
  npwq = ngk(ikq)
  if (nkb.gt.0) then
     if (noncolin) then
        call zgemm ('N', 'N', npwq, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     else
        call zgemm ('N', 'N', npwq, nbnd, nkb, &
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
              igg = igk_k (ig,ikq)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, ikq) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              IF (noncolin) THEN
                 call zaxpy(npwq,ps2_nc(ikb,1,ibnd,ipol),aux,1,dvpsi(1,ibnd),1)
                 call zaxpy(npwq,ps2_nc(ikb,2,ibnd,ipol),aux,1, &
                                                         dvpsi(1+npwx,ibnd),1)
              ELSE
                 call zaxpy (npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1)
              END IF
           enddo
        endif
     enddo

  enddo
  deallocate (aux)
  IF (noncolin) THEN
     deallocate (ps2_nc)
     deallocate (ps1_nc)
     deallocate (deff_nc)
  ELSE
     deallocate (ps2)
     deallocate (ps1)
     deallocate (deff)
  END IF

  call stop_clock ('dvqpsi_us_on')
  return
end subroutine dvqpsi_us_only
