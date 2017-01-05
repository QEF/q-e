!
! Copyright (C) 2005-2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE local_dos_mag(spin_component, kpoint, kband, raux)
  !----------------------------------------------------------------------------
  !
  ! ... compute the contribution of band "kband" at k-point "kpoint"
  ! ... to the noncolinear magnetization for the given "spin_component"
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, g
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : nls, doublegrid
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE scf,                  ONLY : rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE wvfct,                ONLY : nbnd, npwx
  USE becmod,               ONLY : calbec
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER, INTENT (IN) :: spin_component, kpoint, kband
  REAL(DP), INTENT(OUT):: raux(dfftp%nnr)

  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np, npw
  ! counters on beta functions, atoms, pseudopotentials
  INTEGER :: ir, is, ig, ibnd, ik
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands
  ! counter on k points
  !
  REAL(DP) :: w1
  ! weights
  COMPLEX(DP), ALLOCATABLE :: becp_nc(:,:,:)
  ! contains <beta|psi>
  !
  COMPLEX(DP), ALLOCATABLE :: be1(:,:), be2(:,:)
  !
  INTEGER :: ipol, kh, kkb, is1, is2

  becsum(:,:,:) = 0.D0
  rho%of_r(:,:) = 0.D0
  w1=1.D0/omega

  ALLOCATE( becp_nc( nkb, npol, nbnd ) )
  IF (lspinorb) ALLOCATE(be1(nhm,2), be2(nhm,2))
  !
  ! ... here we compute, for the specified k-point and band,
  ! ... the magnetization for the specified spin component
  ! ... Following code is a stripped-down version of "sum_band",
  ! ... without summation over k-points and bands, without symmetrization
  !
  ik = kpoint
  !
     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     IF (nkb > 0) CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     CALL calbec ( npw, vkb, evc, becp_nc)
     !
     ibnd = kband
        !
        psic_nc = (0.D0,0.D0)
        DO ig = 1, npw
           psic_nc(nls(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
           psic_nc(nls(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
        ENDDO
        DO ipol=1,npol
           CALL invfft ('Wave', psic_nc(:,ipol), dffts)
        ENDDO
        IF (spin_component==1) THEN
           DO ir = 1,dffts%nnr
              rho%of_r(ir,2) = rho%of_r(ir,2) + 2.D0*w1* &
                   (dble(psic_nc(ir,1))* dble(psic_nc(ir,2)) + &
                   aimag(psic_nc(ir,1))*aimag(psic_nc(ir,2)))
           ENDDO
        ENDIF
        IF (spin_component==2) THEN
           DO ir = 1,dffts%nnr
              rho%of_r(ir,3) = rho%of_r(ir,3) + 2.D0*w1* &
                   (dble(psic_nc(ir,1))*aimag(psic_nc(ir,2)) - &
                   dble(psic_nc(ir,2))*aimag(psic_nc(ir,1)))
           ENDDO
        ENDIF
        IF (spin_component==3) THEN
           DO ir = 1,dffts%nnr
              rho%of_r(ir,4) = rho%of_r(ir,4) + w1* &
                   (dble(psic_nc(ir,1))**2+aimag(psic_nc(ir,1))**2 &
                   -dble(psic_nc(ir,2))**2-aimag(psic_nc(ir,2))**2)
           ENDDO
        ENDIF
        !
        ! ... spin-orbit contribution
        !
        ijkb0 = 0
        DO np = 1, ntyp
           !
           IF ( upf(np)%tvanp ) THEN
              !
              DO na = 1, nat
                 !
                 IF (ityp(na)==np) THEN
                    !
                    IF (upf(np)%has_so) THEN
                       be1=(0.d0,0.d0)
                       be2=(0.d0,0.d0)
                       DO ih = 1, nh(np)
                          ikb = ijkb0 + ih
                          DO kh = 1, nh(np)
                             IF ((nhtol(kh,np)==nhtol(ih,np)).and. &
                                  (nhtoj(kh,np)==nhtoj(ih,np)).and.     &
                                  (indv(kh,np)==indv(ih,np))) THEN
                                kkb=ijkb0 + kh
                                DO is1=1,2
                                   DO is2=1,2
                                      be1(ih,is1)=be1(ih,is1)+  &
                                           fcoef(ih,kh,is1,is2,np)*  &
                                           becp_nc(kkb,is2,ibnd)
                                      be2(ih,is1)=be2(ih,is1)+ &
                                           fcoef(kh,ih,is2,is1,np)* &
                                           conjg(becp_nc(kkb,is2,ibnd))
                                   ENDDO
                                ENDDO
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDIF
                    ijh = 1
                    !
                    DO ih = 1, nh(np)
                       !
                       ikb = ijkb0 + ih
                       !
                       IF (upf(np)%has_so) THEN
                          IF (spin_component==1) &
                               becsum(ijh,na,2)=becsum(ijh,na,2)+ &
                               (be1(ih,2)*be2(ih,1)+ be1(ih,1)*be2(ih,2))
                          IF (spin_component==2) &
                               becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                               (0.d0,-1.d0)*      &
                               (be1(ih,2)*be2(ih,1)-be1(ih,1)*be2(ih,2))
                          IF (spin_component==3) &
                               becsum(ijh,na,4)=becsum(ijh,na,4)+ &
                               (be1(ih,1)*be2(ih,1)-be1(ih,2)*be2(ih,2))
                       ELSE
                          IF (spin_component==1) &
                               becsum(ijh,na,2)=becsum(ijh,na,2)  &
                               + (conjg(becp_nc(ikb,2,ibnd))   &
                               *becp_nc(ikb,1,ibnd)   &
                               +     conjg(becp_nc(ikb,1,ibnd))   &
                               *becp_nc(ikb,2,ibnd) )
                          IF (spin_component==2) &
                               becsum(ijh,na,3)=becsum(ijh,na,3)+2.d0   &
                               *aimag(conjg(becp_nc(ikb,1,ibnd))* &
                               becp_nc(ikb,2,ibnd) )
                          IF (spin_component==3) &
                               becsum(ijh,na,4) = becsum(ijh,na,4)    &
                               + ( conjg(becp_nc(ikb,1,ibnd)) &
                               *becp_nc(ikb,1,ibnd)  &
                               -      conjg(becp_nc(ikb,2,ibnd)) &
                               *becp_nc(ikb,2,ibnd) )
                       ENDIF
                       !
                       ijh = ijh + 1
                       !
                       DO jh = ( ih + 1 ), nh(np)
                          !
                          jkb = ijkb0 + jh
                          !
                          IF (upf(np)%has_so) THEN
                             IF (spin_component==1) &
                                  becsum(ijh,na,2)=becsum(ijh,na,2)+( &
                                  (be1(jh,2)*be2(ih,1)+be1(jh,1)*be2(ih,2))+&
                                  (be1(ih,2)*be2(jh,1)+be1(ih,1)*be2(jh,2)))
                             IF (spin_component==2) &
                                  becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                  (0.d0,-1.d0)*((be1(jh,2)*&
                                  be2(ih,1)-be1(jh,1)*be2(ih,2))+ &
                                  (be1(ih,2)*be2(jh,1)-be1(ih,1)*be2(jh,2)))
                             IF (spin_component==3) &
                                  becsum(ijh,na,4)=becsum(ijh,na,4)+ &
                                  ((be1(jh,1)*be2(ih,1)- &
                                  be1(jh,2)*be2(ih,2))+  &
                                  (be1(ih,1)*be2(jh,1)-  &
                                  be1(ih,2)*be2(jh,2)) )
                          ELSE
                             IF (spin_component==1) &
                                  becsum(ijh,na,2)=becsum(ijh,na,2)+ 2.d0* &
                                  dble(conjg(becp_nc(ikb,2,ibnd))* &
                                  becp_nc(jkb,1,ibnd) + &
                                  conjg(becp_nc(ikb,1,ibnd))* &
                                  becp_nc(jkb,2,ibnd) )
                             IF (spin_component==2) &
                                  becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                  2.d0*aimag(conjg(becp_nc(ikb,1,ibnd))* &
                                  becp_nc(jkb,2,ibnd) + &
                                  conjg(becp_nc(ikb,1,ibnd))* &
                                  becp_nc(jkb,2,ibnd) )
                             IF (spin_component==3) &
                                  becsum(ijh,na,4)=becsum(ijh,na,4)+ 2.d0* &
                                  dble(conjg(becp_nc(ikb,1,ibnd))* &
                                  becp_nc(jkb,1,ibnd) - &
                                  conjg(becp_nc(ikb,2,ibnd))* &
                                  becp_nc(jkb,2,ibnd) )
                          ENDIF
                          !
                          ijh = ijh + 1
                          !
                       ENDDO
                       !
                    ENDDO
                    !
                    ijkb0 = ijkb0 + nh(np)
                    !
                 ENDIF
                 !
              ENDDO
              !
           ELSE
              !
              DO na = 1, nat
                 !
                 IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
        !
     !
  !
  IF ( doublegrid ) THEN
    is=spin_component+1
    CALL interpolate( rho%of_r(1,is), rho%of_r(1,is), 1 )
  ENDIF
  !
  ! ... Here we add the Ultrasoft contribution to the charge and magnetization
  !
  IF ( okvan ) CALL addusdens(rho%of_r(:,:))

  DO ir=1,dfftp%nnr
     raux(ir)=rho%of_r(ir,spin_component+1)
  ENDDO
  !
  IF (lspinorb) DEALLOCATE(be1, be2)
  DEALLOCATE( becp_nc )
  RETURN
  !
END SUBROUTINE local_dos_mag
