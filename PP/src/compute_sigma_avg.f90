!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_sigma_avg(sigma_avg,becp_nc,ik,lsigma)
  !
  ! This subroutine calculates the average value of the spin on
  ! the spinor wavefunctions.
  !
  USE kinds,                ONLY : DP
  USE noncollin_module,     ONLY : noncolin, npol
  USE cell_base,            ONLY : alat, at, tpiba, omega
  USE spin_orb,             ONLY : fcoef
  USE uspp,                 ONLY : nkb,qq,vkb,nhtol,nhtoj,nhtolm,indv
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wvfct,                ONLY : nbnd, npwx
  USE wavefunctions_module, ONLY : evc, psic_nc
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : g,gg
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE scf,                  ONLY : rho
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE mp_global,            ONLY : me_pool, intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft


  IMPLICIT NONE

  LOGICAL :: lsigma(4)
  ! if true the expectation value in this direction is calculated
  COMPLEX(DP) :: becp_nc(nkb,npol,nbnd)
  !
  REAL(kind=DP) :: sigma_avg(4,nbnd)
  INTEGER, INTENT(in) :: ik

  INTEGER :: ibnd, ig, ir, ijkb0, na, np, ih, ikb, jh
  INTEGER :: ipol, kh, kkb, is1, is2, npw, npwi, npwf
  INTEGER :: li, mi, lj, mj, mi1, i, j, k, ijk
  REAL(DP) :: magtot1(4), magtot2(4)
  REAL(DP) :: x0, y0, dx, dy, r_cut, r_aux, xx, yy
  COMPLEX(DP), ALLOCATABLE :: be1(:,:), qq_lz(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dfx(:), dfy(:)

  COMPLEX(DP) :: c_aux, zdotc

  IF (.not.(lsigma(1).or.lsigma(2).or.lsigma(3).or.lsigma(4))) RETURN

  ALLOCATE(be1(nhm,2))
  ALLOCATE(dfx(dffts%nnr), dfy(dffts%nnr))
  ALLOCATE(qq_lz(nhm,nhm,ntyp))

  sigma_avg=0.d0

  r_cut = 7.d0
  x0 = 0.5d0*at(1,1)*alat
  y0 = 0.5d0*at(2,2)*alat
  dx = at(1,1)*alat/dffts%nr1
  dy = at(2,2)*alat/dffts%nr2

  qq_lz = 0.d0

  DO np=1, ntyp
     DO ih = 1, nh (np)
        li = nhtol(ih,np)
        mi = nhtolm(ih,np) - li**2
        IF (mi==2) THEN
           mi1 = 3
           c_aux = -(0.d0,1.d0)
        ELSEIF (mi==3) THEN
           mi1 = 2
           c_aux = (0.d0,1.d0)
        ELSEIF (mi==4) THEN
           mi1 = 5
           c_aux = -(0.d0,2.d0)
        ELSEIF (mi==5) THEN
           mi1 = 4
           c_aux = (0.d0,2.d0)
        ENDIF
        DO jh = ih+1, nh (np)
           lj = nhtol(jh,np)
           mj = nhtolm(jh,np) - lj**2
           IF (lj==li.and.mj==mi1) THEN
              IF (mj>mi) THEN
                 r_aux = qq(ih,jh-1,np)
              ELSE
                 r_aux = qq(ih,jh+1,np)
              ENDIF
              qq_lz(ih,jh,np) = c_aux * r_aux
           ENDIF
        ENDDO
     ENDDO

     DO ih = 1, nh (np)
        DO jh = 1, ih-1
           qq_lz(ih,jh,np) = conjg(qq_lz(jh,ih,np))
        ENDDO
     ENDDO
  ENDDO

  npw = ngk(ik)
  DO ibnd = 1, nbnd
     rho%of_r = 0.d0
     magtot1 = 0.d0
     magtot2 = 0.d0

     !--  Pseudo part
     psic_nc = (0.D0,0.D0)
     DO ig = 1, npw
        psic_nc(nls(igk_k(ig,ik)), 1)=evc(ig     ,ibnd)
        psic_nc(nls(igk_k(ig,ik)), 2)=evc(ig+npwx,ibnd)
     ENDDO
     DO ipol=1,npol
        CALL invfft ('Wave', psic_nc(:,ipol), dffts)
     ENDDO
     !
     ! Calculate the three components of the magnetization
     ! (stored in rho%of_r(ir,2-4) )
     !
     IF (lsigma(1)) THEN
        DO ir = 1,dffts%nnr
           rho%of_r(ir,2) = rho%of_r(ir,2) + 2.D0* &
                (REAL(psic_nc(ir,1))*REAL(psic_nc(ir,2)) + &
                aimag(psic_nc(ir,1))*aimag(psic_nc(ir,2)))
        ENDDO
        IF (doublegrid) CALL interpolate( rho%of_r(1,2), rho%of_r(1,2), 1 )
     ENDIF
     IF (lsigma(2)) THEN
        DO ir = 1,dffts%nnr
           rho%of_r(ir,3) = rho%of_r(ir,3) + 2.D0* &
                (REAL(psic_nc(ir,1))*aimag(psic_nc(ir,2)) - &
                REAL(psic_nc(ir,2))*aimag(psic_nc(ir,1)))
        ENDDO
        IF (doublegrid) CALL interpolate( rho%of_r(1,3), rho%of_r(1,3), 1 )
     ENDIF
     IF (lsigma(3)) THEN
        DO ir = 1,dffts%nnr
           rho%of_r(ir,4) = rho%of_r(ir,4) + &
                (REAL(psic_nc(ir,1))**2+aimag(psic_nc(ir,1))**2 &
                -REAL(psic_nc(ir,2))**2-aimag(psic_nc(ir,2))**2)
        ENDDO
        IF (doublegrid) CALL interpolate( rho%of_r(1,4), rho%of_r(1,4), 1 )
     ENDIF

     IF (lsigma(4)) THEN
        !-- Calculate pseudo part of L_z
        DO ipol = 1, npol
           dfx = 0.d0
           dfy = 0.d0
           npwi=(ipol-1)*npwx+1
           npwf=(ipol-1)*npwx+npw
           dfx(nls(igk_k(1:npw,ik))) = (xk(1,ik)+g(1,igk_k(1:npw,ik)))*tpiba* &
                (0.d0,1.d0)*evc(npwi:npwf,ibnd)
           dfy(nls(igk_k(1:npw,ik))) = (xk(2,ik)+g(2,igk_k(1:npw,ik)))*tpiba* &
                (0.d0,1.d0)*evc(npwi:npwf,ibnd)
           CALL invfft ('Wave', dfx, dffts)
           CALL invfft ('Wave', dfy, dffts)
           DO i = 1, dffts%nr1
              xx = (i-1)*dx - x0
              DO j = 1, dffts%nr2
                 yy = (j-1)*dy - y0
                 r_aux = DSQRT (xx**2 + yy**2)
                 IF (r_aux<=r_cut) THEN
                    DO k = 1, dffts%npp(me_pool+1)
                       ijk = i + (j-1)*dffts%nr1x + (k-1)*dffts%nr1x*dffts%nr2x
                       dfx(ijk) = xx * dfy(ijk) - yy * dfx(ijk)
                    ENDDO
                 ELSE
                    DO k = 1, dffts%npp(me_pool+1)
                       ijk = i + (j-1)*dffts%nr1x + (k-1)*dffts%nr1x*dffts%nr2x
                       dfx (ijk) = 0.d0
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
           c_aux = zdotc(dffts%nnr, psic_nc(1,ipol), 1, dfx, 1)
           magtot1(4) = magtot1(4) + aimag(c_aux)
        ENDDO
        CALL mp_sum( magtot1(4), intra_bgrp_comm )
        magtot1(4) = magtot1(4)/(dffts%nr1*dffts%nr2*dffts%nr3)
     ENDIF

     DO ipol=1,3
        IF (lsigma(ipol)) THEN
           DO ir = 1,dfftp%nnr
              magtot1(ipol) = magtot1(ipol) + rho%of_r(ir,ipol+1)
           ENDDO
           CALL mp_sum( magtot1(ipol), intra_bgrp_comm )
           magtot1(ipol) = magtot1(ipol) / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
        ENDIF
     ENDDO

     !-- Augmentation part

     ijkb0 = 0
     !
     DO np = 1, ntyp
        !
        IF ( upf(np)%tvanp ) THEN
           !
           DO na = 1, nat
              !
              IF (ityp(na)==np) THEN
                 !
                 be1 = 0.d0
                 DO ih = 1, nh(np)
                    ikb = ijkb0 + ih
                    IF (upf(np)%has_so) THEN
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
                                ENDDO
                             ENDDO
                          ENDIF
                       ENDDO
                    ELSE
                       DO is1=1,2
                          be1(ih,is1) = becp_nc(ikb,is1,ibnd)
                       ENDDO
                    ENDIF
                 ENDDO
                 IF (lsigma(1)) THEN
                    DO ih = 1, nh(np)
                       magtot2(1)=magtot2(1)+ 2.d0*qq(ih,ih,np)  &
                            * REAL( be1(ih,2)*conjg(be1(ih,1)) )
                       DO jh = ih + 1, nh(np)
                          magtot2(1)=magtot2(1)+2.d0*qq(ih,jh,np) &
                               * REAL( be1(jh,2)*conjg(be1(ih,1))+ &
                               be1(jh,1)*conjg(be1(ih,2)) )

                       ENDDO
                    ENDDO
                 ENDIF
                 IF (lsigma(2)) THEN
                    DO ih = 1, nh(np)
                       magtot2(2)=magtot2(2)+ 2.d0*qq(ih,ih,np)*aimag   &
                            ( be1(ih,2)*conjg(be1(ih,1)) )
                       DO jh = ih + 1, nh(np)
                          magtot2(2)=magtot2(2) + 2.d0*qq(ih,jh,np)*aimag &
                               (  be1(jh,2) * conjg(be1(ih,1)) &
                               - be1(jh,1) * conjg(be1(ih,2)) )
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (lsigma(3)) THEN
                    DO ih = 1, nh(np)
                       magtot2(3) = magtot2(3) + qq(ih,ih,np)*              &
                            ( abs(be1(ih,1))**2 - abs(be1(ih,2))**2 )
                       DO jh = ih + 1, nh(np)
                          magtot2(3) = magtot2(3) + 2.d0*qq(ih,jh,np) &
                               * REAL( be1(jh,1)*conjg(be1(ih,1)) &
                               -be1(jh,2)*conjg(be1(ih,2)) )
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (lsigma(4)) THEN
                    DO ih = 1, nh(np)
                       DO jh = ih + 1, nh(np)
                          magtot2(4)= magtot2(4)+2.d0*REAL(qq_lz(ih,jh,np)*  &
                               ( conjg(be1(ih,1))*be1(jh,1) +           &
                               conjg(be1(ih,2))*be1(jh,2) ) )
                       ENDDO
                    ENDDO
                 ENDIF
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

     DO ipol=1,3
        IF (lsigma(ipol)) &
             sigma_avg(ipol,ibnd) = 0.5d0 * ( magtot1(ipol) + magtot2(ipol) )
     ENDDO
     IF (lsigma(4)) &
          sigma_avg(4,ibnd) =  magtot1(4) + magtot2(4) + sigma_avg(3,ibnd)

  ENDDO

  DEALLOCATE(be1)
  DEALLOCATE(dfx,dfy)
  DEALLOCATE(qq_lz)

  RETURN
END SUBROUTINE compute_sigma_avg
