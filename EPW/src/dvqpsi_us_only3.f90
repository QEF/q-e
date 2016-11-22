  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! 
  ! Copyright (C) 2001-2008 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! adapted from PH/dvqpsi_us_only (QE)
  !
  !----------------------------------------------------------------------
  subroutine dvqpsi_us_only3 (ik, uact,xxk)
  !----------------------------------------------------------------------
  !!
  !! This routine calculates dV_bare/dtau * psi for one perturbation
  !! with a given q. The displacements are described by a vector uact.
  !! The result is stored in dvpsi. The routine is called for each k point
  !! and for each pattern u. It computes simultaneously all the bands.
  !! This routine implements Eq. B29 of PRB 64, 235118 (2001).
  !! Only the contribution of the nonlocal potential is calculated here.
  !!
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : tpiba
  USE gvect,      ONLY : g
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,   ONLY : lsda, current_spin, isk, nspin
  USE spin_orb,   ONLY : lspinorb
  USE wvfct,      ONLY : nbnd, npwx, et
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp,       ONLY : okvan, nkb, vkb
  USE uspp_param, ONLY : nh, nhm
  USE qpoint,     ONLY : npwq
  USE phus,       ONLY : int1, int1_nc, int2, int2_so, alphap
  USE lrus,       ONLY : becp1
  USE eqv,        ONLY : dvpsi
  USE elph2,      ONLY : igkq, lower_band, upper_band

  implicit none
  !
  INTEGER, INTENT(in) :: ik
  !! Input: the k point
  REAL(kind=DP), INTENT(in) :: xxk(3) 
  !! input: the k point (cartesian coordinates)
  COMPLEX(kind=DP), INTENT(in) :: uact (3 * nat)
  !! input: the pattern of displacements
  !
  !   And the local variables
  !
  INTEGER :: na
  !! Counter on atoms
  INTEGER :: nb
  !! Counter on atoms
  INTEGER :: mu
  !! Counter on modes
  INTEGER :: nu
  !! Counter on modes
  INTEGER :: ig
  !! Counter on G vectors
  INTEGER :: igg
  !! Auxiliary counter on G vectors
  INTEGER :: nt
  !! Counter on atomic types
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: ijkb0
  !! Auxiliary variable for counting
  INTEGER :: ikb
  !! Counter on becp functions
  INTEGER :: jkb
  !! Counter on becp functions
  INTEGER :: ipol
  !! Counter on polarizations
  INTEGER :: ih
  !! Counter on nh
  INTEGER :: jh
  !! Counter on nh
  INTEGER :: is
  !! Counter on polarization
  INTEGER :: js
  !! Counter on polarization
  INTEGER ::  ijs
  !! Counter on combined is and js polarization

  REAL(kind=DP), parameter :: eps = 1.d-12

  complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:), deff_nc(:,:,:,:)
  real(DP), allocatable :: deff(:,:,:)
  complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  ! work space

  logical :: ok

  call start_clock ('dvqpsi_us_on')
  if (noncolin) then
     allocate (ps1_nc(nkb , npol, lower_band: upper_band))
     allocate (ps2_nc(nkb , npol, lower_band: upper_band , 3))
     allocate (deff_nc(nhm, nhm, nat, nspin))
  else
     allocate (ps1 ( nkb , lower_band: upper_band))
     allocate (ps2 ( nkb , lower_band: upper_band , 3))
     allocate (deff(nhm, nhm, nat))
  end if
  allocate (aux ( npwx))
  if (lsda) current_spin = isk (ik)
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
  do ibnd = lower_band, upper_band
     IF (noncolin) THEN
        CALL compute_deff_nc(deff_nc,et(ibnd,ik))
     ELSE
        CALL compute_deff(deff,et(ibnd,ik))
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
  if (nkb.gt.0) then
     if (noncolin) then
        call zgemm ('N', 'N', npwq, (upper_band-lower_band+1)*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     else
        call zgemm ('N', 'N', npwq, (upper_band-lower_band+1), nkb, &
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
           do ibnd = lower_band, upper_band
              ok = ok.or.(abs (ps2_nc (ikb, 1, ibnd, ipol) ).gt.eps).or. &
                         (abs (ps2_nc (ikb, 2, ibnd, ipol) ).gt.eps)
           end do
        ELSE
           do ibnd = lower_band, upper_band
              ok = ok.or. (abs (ps2 (ikb, ibnd, ipol) ) .gt.eps)
           enddo
        ENDIF
        if (ok) then
           do ig = 1, npwq
              igg = igkq (ig)
              !aux (ig) =  vkb(ig, ikb) * (xk(ipol,ikq) + g(ipol, igg) )
              aux (ig) =  vkb(ig, ikb) * (xxk(ipol) + g(ipol, igg) )
           enddo
           do ibnd = lower_band, upper_band
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
end subroutine dvqpsi_us_only3
