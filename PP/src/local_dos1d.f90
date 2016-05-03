
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE local_dos1d (ik, kband, plan)
  !--------------------------------------------------------------------
  !
  !     calculates |psi|^2 for band kband at point ik
  !
  USE kinds,     ONLY: dp
  USE cell_base, ONLY: omega
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp
  USE fft_base,  ONLY: dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvecs,   ONLY : nls, doublegrid
  USE lsda_mod, ONLY: current_spin
  USE uspp, ONLY: becsum, indv, nhtol, nhtoj
  USE uspp_param, ONLY: upf, nh, nhm
  USE wvfct, ONLY: npwx, wg
  USE klist, ONLY: ngk, igk_k
  USE noncollin_module, ONLY: noncolin, npol
  USE spin_orb, ONLY: lspinorb, fcoef
  USE wavefunctions_module,  ONLY: evc, psic, psic_nc
  USE becmod, ONLY: bec_type, becp
  IMPLICIT NONE
  !
  ! input variables
  !
  INTEGER :: ik, kband
  ! input: the k point
  ! input: the band

  real(DP) :: plan (dfftp%nr3)
  ! output: the planar average of this state
  !
  !    Additional local variables for Ultrasoft PP's
  !

  INTEGER :: npw, ikb, jkb, ijkb0, ih, jh, na, ijh, ipol, np
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for ijkb0
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on atoms
  ! counter on composite beta functions
  ! the pseudopotential
  !
  !    And here the local variables
  !
  INTEGER :: ir, ig, ibnd, is1, is2, kkb, kh
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands

  real(DP) :: w1
  ! the weight of one k point
  real(DP), ALLOCATABLE :: aux (:)
  ! auxiliary for rho

  COMPLEX(DP), ALLOCATABLE :: prho (:), be1(:,:), be2(:,:)
  ! complex charge for fft

  ALLOCATE (prho(dfftp%nnr))
  ALLOCATE (aux(dfftp%nnr))
  IF (lspinorb) THEN
     ALLOCATE(be1(nhm,2))
     ALLOCATE(be2(nhm,2))
  ENDIF

  aux(:) = 0.d0
  becsum(:,:,:) = 0.d0

  npw = ngk(ik)
  wg (kband, ik) = 1.d0
  !
  !
  !     First compute the square modulus of the state kband,ik on the smooth
  !     mesh
  !
  IF (noncolin) THEN
     psic_nc = (0.d0,0.d0)
     DO ig = 1, npw
        psic_nc (nls (igk_k (ig,ik) ), 1 ) = evc (ig     , kband)
        psic_nc (nls (igk_k (ig,ik) ), 2 ) = evc (ig+npwx, kband)
     ENDDO
     DO ipol=1,npol
        CALL invfft ('Wave', psic_nc(:,ipol), dffts)
     ENDDO

     w1 = wg (kband, ik) / omega
     DO ipol=1,npol
        DO ir = 1, dffts%nnr
           aux(ir) = aux(ir) + w1 * ( dble(psic_nc(ir,ipol))**2 + &
                                     aimag(psic_nc(ir,ipol))**2 )
        ENDDO
     ENDDO
  ELSE
     psic(1:dffts%nnr) = (0.d0,0.d0)
     DO ig = 1, npw
        psic (nls (igk_k (ig,ik) ) ) = evc (ig, kband)
     ENDDO
     CALL invfft ('Wave', psic, dffts)

     w1 = wg (kband, ik) / omega
     DO ir = 1, dffts%nnr
        aux(ir) = aux(ir) + w1 * (dble(psic(ir))**2 + aimag(psic(ir))**2)
     ENDDO
  ENDIF

  !
  !    If we have a US pseudopotential we compute here the becsum term
  !
  ibnd = kband

  w1 = wg (ibnd, ik)
  ijkb0 = 0
  DO np = 1, ntyp
     IF (upf(np)%tvanp) THEN
        DO na = 1, nat
           IF (ityp (na) == np) THEN
              IF (noncolin) THEN
                 IF (upf(np)%has_so) THEN
                    be1=(0.d0,0.d0)
                    be2=(0.d0,0.d0)
                    DO ih = 1, nh(np)
                       ikb = ijkb0 + ih
                       DO kh = 1, nh(np)
                          IF ((nhtol(kh,np)==nhtol(ih,np)).and. &
                              (nhtoj(kh,np)==nhtoj(ih,np)).and. &
                              (indv(kh,np)==indv(ih,np))) THEN
                             kkb=ijkb0 + kh
                             DO is1=1,2
                                DO is2=1,2
                                   be1(ih,is1)=be1(ih,is1)+ &
                                        fcoef(ih,kh,is1,is2,np)* &
                                        becp%nc(kkb,is2,ibnd)
                                   be2(ih,is1)=be2(ih,is1)+ &
                                        fcoef(kh,ih,is2,is1,np)* &
                                     conjg(becp%nc(kkb,is2,ibnd))
                                ENDDO
                             ENDDO
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDIF
              ENDIF
              ijh = 1
              DO ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 IF (noncolin) THEN
                    IF (upf(np)%has_so) THEN
                        becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*    &
                            (be1(ih,1)*be2(ih,1)+be1(ih,2)*be2(ih,2))
                    ELSE
                       DO ipol=1,npol
                          becsum(ijh,na,current_spin) = &
                             becsum(ijh,na,current_spin) + w1 * &
                               dble( conjg(becp%nc(ikb,ipol,ibnd)) * &
                                           becp%nc(ikb,ipol,ibnd) )
                       ENDDO
                    ENDIF
                 ELSE
                    becsum(ijh,na,current_spin) = &
                        becsum(ijh,na,current_spin) + w1 * &
                        dble( conjg(becp%k(ikb,ibnd)) * becp%k(ikb,ibnd) )
                 ENDIF
                 ijh = ijh + 1
                 DO jh = ih + 1, nh (np)
                    jkb = ijkb0 + jh
                    IF (noncolin) THEN
                       IF (upf(np)%has_so) THEN
                          becsum(ijh,na,1)=becsum(ijh,na,1) &
                              + w1*((be1(jh,1)*be2(ih,1)+   &
                                     be1(jh,2)*be2(ih,2))+  &
                                    (be1(ih,1)*be2(jh,1)+   &
                                     be1(ih,2)*be2(jh,2)) )
                       ELSE
                          DO ipol=1,npol
                             becsum(ijh,na,current_spin) = &
                                becsum(ijh,na,current_spin) + w1 * 2.d0 * &
                                dble( conjg(becp%nc(ikb,ipol,ibnd))  &
                                        * becp%nc(jkb,ipol,ibnd) )
                          ENDDO
                       ENDIF
                    ELSE
                       becsum(ijh,na,current_spin) = &
                           becsum(ijh,na,current_spin) + w1 * 2.d0 * &
                           dble( conjg(becp%k(ikb,ibnd)) * becp%k(jkb,ibnd) )
                    ENDIF
                    ijh = ijh + 1
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (np)
           ENDIF
        ENDDO
     ELSE
        DO na = 1, nat
           IF (ityp (na) ==np) ijkb0 = ijkb0 + nh (np)
        ENDDO
     ENDIF
  ENDDO
  !
  !    Interpolate on the thick mesh and pass to reciprocal space
  !
  IF (doublegrid) THEN
     CALL interpolate (aux, aux, 1)
  ENDIF
  DO ir = 1, dfftp%nnr
     prho (ir) = cmplx(aux (ir), 0.d0,kind=DP)
  ENDDO
  CALL fwfft ('Dense', prho, dfftp)
  !
  !    Here we add the US contribution to the charge for the atoms which n
  !    it. Or compute the planar average in the NC case.
  !
  CALL addusdens1d (plan, prho)
  !
  DEALLOCATE (aux)
  DEALLOCATE (prho)
  IF (lspinorb) THEN
     DEALLOCATE(be1)
     DEALLOCATE(be2)
  ENDIF
  !
  RETURN
END SUBROUTINE local_dos1d
