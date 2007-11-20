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
USE spin_orb,             ONLY : so, fcoef
USE uspp,                 ONLY : nkb,qq,vkb,nhtol,nhtoj,nhtolm,indv
USE uspp_param,           ONLY : upf, nh, nhm
USE wvfct,                ONLY : nbnd, npwx, npw, igk 
USE wavefunctions_module, ONLY : evc, psic_nc
USE klist,                ONLY : nks, xk
USE gvect,                ONLY : g,gg,nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx
USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                  nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
USE scf,                  ONLY : rho
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE mp_global,            ONLY : me_pool
USE pffts,                ONLY : npps


IMPLICIT NONE

LOGICAL :: lsigma(4)
! if true the expectation value in this direction is calculated
COMPLEX(DP) :: becp_nc(nkb,npol,nbnd)
! 
REAL(KIND=DP) :: sigma_avg(4,nbnd)
INTEGER :: ik    

INTEGER :: ibnd, ig, ir, ijkb0, na, np, ih, ikb, ijh, jh, jkb    
INTEGER :: ipol, kh, kkb, is1, is2, npwi, npwf
INTEGER :: li, mi, lj, mj, mi1, i, j, k, ijk
REAL(DP) :: magtot1(4), magtot2(4)
REAL(DP) :: x0, y0, dx, dy, r_cut, r_aux, xx, yy
COMPLEX(DP), ALLOCATABLE :: be1(:,:), qq_lz(:,:,:)
COMPLEX(DP), ALLOCATABLE :: dfx(:), dfy(:)

COMPLEX(DP) :: c_aux, ZDOTC

IF (.NOT.(lsigma(1).OR.lsigma(2).OR.lsigma(3).OR.lsigma(4))) RETURN

ALLOCATE(be1(nhm,2))
ALLOCATE(dfx(nrxxs), dfy(nrxxs))
ALLOCATE(qq_lz(nhm,nhm,ntyp))

sigma_avg=0.d0

r_cut = 7.d0
x0 = 0.5d0*at(1,1)*alat
y0 = 0.5d0*at(2,2)*alat
dx = at(1,1)*alat/nr1s
dy = at(2,2)*alat/nr2s

qq_lz = 0.d0

DO np=1, ntyp
   DO ih = 1, nh (np)
      li = nhtol(ih,np)
      mi = nhtolm(ih,np) - li**2
      IF (mi.EQ.2) THEN
         mi1 = 3
         c_aux = -(0.d0,1.d0)
      ELSE IF (mi.EQ.3) THEN
         mi1 = 2
         c_aux = (0.d0,1.d0)
      ELSE IF (mi.EQ.4) THEN
         mi1 = 5
         c_aux = -(0.d0,2.d0)
      ELSE IF (mi.EQ.5) THEN
         mi1 = 4
         c_aux = (0.d0,2.d0)
      END IF
      DO jh = ih+1, nh (np)
         lj = nhtol(jh,np)
         mj = nhtolm(jh,np) - lj**2
         IF (lj.EQ.li.AND.mj.EQ.mi1) THEN
            IF (mj.GT.mi) THEN
               r_aux = qq(ih,jh-1,np)
            ELSE
               r_aux = qq(ih,jh+1,np)
            END IF
            qq_lz(ih,jh,np) = c_aux * r_aux
         END IF
      END DO
   END DO

   DO ih = 1, nh (np)
      DO jh = 1, ih-1
         qq_lz(ih,jh,np) = CONJG(qq_lz(jh,ih,np))
      END DO
   END DO
END DO

DO ibnd = 1, nbnd
   rho%of_r = 0.d0
   magtot1 = 0.d0
   magtot2 = 0.d0 

!--  Pseudo part
   psic_nc = (0.D0,0.D0)
   DO ig = 1, npw
      psic_nc(nls(igk(ig)), 1)=evc(ig     ,ibnd)
      psic_nc(nls(igk(ig)), 2)=evc(ig+npwx,ibnd)
   END DO
   DO ipol=1,npol
      call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
   ENDDO
   !
   ! Calculate the three components of the magnetization 
   ! (stored in rho%of_r(ir,2-4) )
   !
   IF (lsigma(1)) THEN
      DO ir = 1,nrxxs
         rho%of_r(ir,2) = rho%of_r(ir,2) + 2.D0* &
                   (REAL(psic_nc(ir,1))*REAL(psic_nc(ir,2)) + &
                   AIMAG(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)))
      END DO
      IF (doublegrid) CALL interpolate( rho%of_r(1,2), rho%of_r(1,2), 1 )
   END IF
   IF (lsigma(2)) THEN
      DO ir = 1,nrxxs
         rho%of_r(ir,3) = rho%of_r(ir,3) + 2.D0* &
                   (REAL(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)) - &
                    REAL(psic_nc(ir,2))*AIMAG(psic_nc(ir,1)))
      END DO
      IF (doublegrid) CALL interpolate( rho%of_r(1,3), rho%of_r(1,3), 1 )
   END IF
   IF (lsigma(3)) THEN
      DO ir = 1,nrxxs
         rho%of_r(ir,4) = rho%of_r(ir,4) + &
                   (REAL(psic_nc(ir,1))**2+AIMAG(psic_nc(ir,1))**2 &
                   -REAL(psic_nc(ir,2))**2-AIMAG(psic_nc(ir,2))**2)
      END DO
      IF (doublegrid) CALL interpolate( rho%of_r(1,4), rho%of_r(1,4), 1 )
   END IF

   IF (lsigma(4)) THEN
      !-- Calculate pseudo part of L_z
      DO ipol = 1, npol
         dfx = 0.d0
         dfy = 0.d0
         npwi=(ipol-1)*npwx+1
         npwf=(ipol-1)*npwx+npw
         dfx(nls(igk(1:npw))) = (xk(1,ik)+g(1,igk(1:npw)))*tpiba* &
                          (0.d0,1.d0)*evc(npwi:npwf,ibnd)
         dfy(nls(igk(1:npw))) = (xk(2,ik)+g(2,igk(1:npw)))*tpiba* &
                          (0.d0,1.d0)*evc(npwi:npwf,ibnd)
         CALL cft3s( dfx, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
         CALL cft3s( dfy, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
         DO i = 1, nr1s
            xx = (i-1)*dx - x0
            DO j = 1, nr2s
               yy = (j-1)*dy - y0
               r_aux = DSQRT (xx**2 + yy**2)
               IF (r_aux.LE.r_cut) THEN
                  DO k = 1, npps(me_pool+1)
                     ijk = i + (j-1) * nrx1s + (k-1) * nrx1s * nrx2s
                     dfx(ijk) = xx * dfy(ijk) - yy * dfx(ijk)
                  END DO
               ELSE
                  DO k = 1, npps(me_pool+1)
                     ijk = i + (j-1) * nrx1s + (k-1) * nrx1s * nrx2s
                     dfx (ijk) = 0.d0
                  END DO
               END IF
            END DO
         END DO
         c_aux = ZDOTC(nrxxs, psic_nc(1,ipol), 1, dfx, 1)
         magtot1(4) = magtot1(4) + AIMAG(c_aux)
      END DO
      CALL reduce( 1, magtot1(4) )
      magtot1(4) = magtot1(4)/(nr1s*nr2s*nr3s)
   END IF

   DO ipol=1,3
      IF (lsigma(ipol)) THEN
         DO ir = 1,nrxx
            magtot1(ipol) = magtot1(ipol) + rho%of_r(ir,ipol+1)
         END DO
         CALL reduce( 1, magtot1(ipol) )
         magtot1(ipol) = magtot1(ipol) / ( nr1 * nr2 * nr3 )
      END IF
   END DO

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
                  IF (so(np)) THEN
                     DO kh = 1, nh(np)
                        IF ((nhtol(kh,np)==nhtol(ih,np)).AND. &
                            (nhtoj(kh,np)==nhtoj(ih,np)).AND.     &
                            (indv(kh,np)==indv(ih,np))) THEN
                           kkb=ijkb0 + kh
                           DO is1=1,2
                              DO is2=1,2
                                 be1(ih,is1)=be1(ih,is1)+  &
                                    fcoef(ih,kh,is1,is2,np)*  &
                                    becp_nc(kkb,is2,ibnd)
                              END DO
                           END DO
                        END IF
                     END DO
                  ELSE
                     DO is1=1,2 
                        be1(ih,is1) = becp_nc(ikb,is1,ibnd)
                     END DO
                  END IF
               END DO
               IF (lsigma(1)) THEN
                  DO ih = 1, nh(np)
                     magtot2(1)=magtot2(1)+ 2.d0*qq(ih,ih,np)  &
                              * REAL( be1(ih,2)*CONJG(be1(ih,1)) )
                     DO jh = ih + 1, nh(np)   
                        magtot2(1)=magtot2(1)+2.d0*qq(ih,jh,np) &
                              * REAL( be1(jh,2)*CONJG(be1(ih,1))+ &
                                   be1(jh,1)*CONJG(be1(ih,2)) )

                     ENDDO
                  ENDDO
               ENDIF
               IF (lsigma(2)) THEN
                  DO ih = 1, nh(np)
                     magtot2(2)=magtot2(2)+ 2.d0*qq(ih,ih,np)*AIMAG   &
                               ( be1(ih,2)*CONJG(be1(ih,1)) )
                     DO jh = ih + 1, nh(np)   
                        magtot2(2)=magtot2(2) + 2.d0*qq(ih,jh,np)*AIMAG &
                               (  be1(jh,2) * CONJG(be1(ih,1)) &
                                - be1(jh,1) * CONJG(be1(ih,2)) )
                     END DO
                  END DO
               END IF
               IF (lsigma(3)) THEN
                  DO ih = 1, nh(np)
                     magtot2(3) = magtot2(3) + qq(ih,ih,np)*              &
                            ( ABS(be1(ih,1))**2 - ABS(be1(ih,2))**2 )
                     DO jh = ih + 1, nh(np)   
                        magtot2(3) = magtot2(3) + 2.d0*qq(ih,jh,np) &
                                   * REAL( be1(jh,1)*CONJG(be1(ih,1)) &
                                        -be1(jh,2)*CONJG(be1(ih,2)) ) 
                     END DO 
                  END DO
               END IF
               IF (lsigma(4)) THEN
                  DO ih = 1, nh(np)
                     DO jh = ih + 1, nh(np)   
                        magtot2(4)= magtot2(4)+2.d0*REAL(qq_lz(ih,jh,np)*  &
                                 ( CONJG(be1(ih,1))*be1(jh,1) +           &
                                   CONJG(be1(ih,2))*be1(jh,2) ) )
                     END DO
                  END DO
               END IF
               !
               ijkb0 = ijkb0 + nh(np)
               !
            END IF
            !
         END DO
         !
      ELSE
         !
         DO na = 1, nat
            !
            IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
            !
         END DO
         !
      END IF
      !
   END DO

   DO ipol=1,3
      IF (lsigma(ipol)) &
         sigma_avg(ipol,ibnd) = 0.5d0 * ( magtot1(ipol) + magtot2(ipol) )
   END DO
      IF (lsigma(4)) &
         sigma_avg(4,ibnd) =  magtot1(4) + magtot2(4) + sigma_avg(3,ibnd)

END DO

DEALLOCATE(be1)
DEALLOCATE(dfx,dfy)
DEALLOCATE(qq_lz)

RETURN
END SUBROUTINE compute_sigma_avg
