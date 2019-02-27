!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
SUBROUTINE dvqhub_barepsi_us2 (ik, dvqhbar, dvqhbar_orth, dvqhbar_orth_lm)
  !----------------------------------------------------------------------------
  !
  ! DFPT+U: This routine calculate several terms entering the 
  !         Hubbard dynamical matrix calculated in dynmat_hub_scf.f90   
  !         These terms are in the cartesian coordinates.
  ! 
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode
  USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : U_projection, Hubbard_l, is_hubbard, Hubbard_J0, offsetU, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, swfcatomk, swfcatomkpq, dwfcatomkpq,  &
                            sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
                            proj1, proj2, effU 
  USE wvfct,         ONLY : npwx, nbnd
  USE uspp,          ONLY : vkb, nkb, okvan 
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE control_lr,    ONLY : lgamma, ofsbeta
  USE units_lr,      ONLY : iuatwfc, iuatswfc
  USE uspp_param,    ONLY : nh
  USE lsda_mod,      ONLY : lsda, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE eqv,           ONLY : dvpsi
  USE scf,           ONLY : rho
  USE mp_bands,      ONLY : intra_bgrp_comm       
  USE mp,            ONLY : mp_sum  
  USE buffers,       ONLY : get_buffer
  !  
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  ! the k point under consideration
  COMPLEX(DP), INTENT(INOUT) :: dvqhbar(npwx,nbnd,3,nat),        & 
                                dvqhbar_orth(npwx,nbnd,3,nat),   &
                                dvqhbar_orth_lm(npwx,nbnd,3,nat) 
  !
  ! Local variables
  !
  COMPLEX(DP), ALLOCATABLE :: dqsphi(:,:), dmqsphi(:,:), dwfcatom_(:), dvqi(:,:), &
                              dvqi_orth(:,:), dvqi_orth_lm(:,:), aux1(:), aux2(:)
  INTEGER :: i, j, k, icart, na, nt, l, ih, n, mu, ig, npw, npwq, &
             ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, op_spin, ikk, ikq, ibeta
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock( 'dvqhub_barepsi_us2' )
  !
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  ALLOCATE (dqsphi(npwx,nwfcU))
  ALLOCATE (dmqsphi(npwx,nwfcU))
  ALLOCATE (dwfcatom_(npwx))
  ALLOCATE (dvqi(npwx,nbnd))
  ALLOCATE (dvqi_orth(npwx,nbnd))
  ALLOCATE (dvqi_orth_lm(npwx,nbnd))
  ALLOCATE (aux1(npwx))
  ALLOCATE (aux2(npwx))
  !
  proj1 = (0.d0, 0.d0)
  proj2 = (0.d0, 0.d0)
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  IF (lsda) THEN
     current_spin = isk(ikk)
     IF (current_spin==1) THEN
        op_spin = 2
     ELSE
        op_spin = 1
     END IF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Compute beta functions at k and k+q     
  !
  CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb)
  IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(:,ikq), vkbkpq)
  !
  ! Calculate the derivatives of beta functions 
  ! d^{icart}beta at k and k+q for all the bands and for 
  ! the 3 cartesian directions
  !
  DO icart = 1, 3
     DO na = 1, nat 
        nt = ityp(na)
        DO ih = 1, nh(nt)
           !
           ibeta = ofsbeta(na) + ih
           !
           CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                      vkb(:,ibeta), dvkb(:,ibeta,icart))
           IF (.NOT.lgamma) &
           CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, &
                      vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
           !
        ENDDO
     ENDDO
  ENDDO
  !
  ! Read \phi at k and k+q from file (unit iuatwfc)
  ! 
  CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
  IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
  !
  ! Read S*\phi at k and k+q from file (unit iuatswfc)
  !
  CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  !
  dvqhbar         = (0.d0, 0.d0) 
  dvqhbar_orth    = (0.d0, 0.d0)
  dvqhbar_orth_lm = (0.d0, 0.d0)
  !   
  DO na = 1, nat
     !   
     DO icart = 1, 3 
        !   
        dqsphi  = (0.d0, 0.d0)    
        dmqsphi = (0.d0, 0.d0)
        !   
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           IF (is_hubbard(nt)) THEN
              !   
              DO m = 1, 2*Hubbard_l(nt)+1
                 !   
                 ihubst = offsetU(nah) + m   ! I m index
                 !
                 IF (nah==na) THEN
                    !
                    ! Calculate | d_icart\phi_(k,I,m)) >
                    !
                    CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                               wfcatomk(:,ihubst), dwfcatom_) 
                    !   
                    ! Calculate  | S d_^(I,icart)\phi_(k,I,m) >
                    !
                    CALL swfc (npw, 1, vkb, dwfcatom_, sdwfcatomk(:,ihubst))
                    !   
                    IF (.NOT.lgamma) THEN
                       !
                       ! Calculate |d_icart\phi_(k+q,I,m)) >
                       !
                       CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, &
                                  wfcatomkpq(:,ihubst), dwfcatom_) 
                       !
                       ! Calculate | S d_^(I,icart)\phi_(k+q,I,m) >
                       !                       
                       CALL swfc (npwq, 1, vkbkpq, dwfcatom_, sdwfcatomkpq(:,ihubst))
                       !
                     ENDIF
                     !  
                 ENDIF
                 !
                 ! Calculate |\Delta_q(S_k \phi_(k,I,m)) >  
                 ! and |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                 !   
                 CALL delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq, &
                                  sdwfcatomk, sdwfcatomkpq, vkb, vkbkpq, dvkb(:,:,icart), &
                                  dvkbkpq(:,:,icart), dqsphi, dmqsphi, 1) 
                 !   
                 ! Calculate: 
                 ! proj1 (ihubst, ibnd) = < S_{k}\phi_(k,I,m) | psi(inbd,k) >
                 ! proj2 (ihubst, ibnd) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(inbd,k) > 
                 !
                 DO ibnd = 1, nbnd
                    proj1(ibnd,ihubst) = ZDOTC (npw, swfcatomk(:,ihubst), 1, evc(:,ibnd), 1)
                    proj2(ibnd,ihubst) = ZDOTC (npw, dmqsphi(:,ihubst), 1, evc(:,ibnd), 1)
                 ENDDO
                 !
              ENDDO ! m
              !    
           ENDIF
           !
        ENDDO ! nah
        !   
        CALL mp_sum(proj1, intra_bgrp_comm) 
        CALL mp_sum(proj2, intra_bgrp_comm)
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           ! For Hubbard_U - Hubbard_J0
           ! 
           IF (is_hubbard(nt)) THEN
              !
              dvqi         = (0.d0, 0.d0)
              dvqi_orth    = (0.d0, 0.d0)
              dvqi_orth_lm = (0.d0, 0.d0)
              !   
              DO ibnd = 1, nbnd
                 !      
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    !   
                    ihubst1 = offsetU(nah) + m1 
                    !   
                    DO ig = 1, npwq
                       !  
                       ! Note the factor 2, we are considering a diagonal term
                       ! 
                       aux1(ig) = 2.d0 * ( dqsphi(ig,ihubst1) * proj1(ibnd,ihubst1) + &
                                           swfcatomkpq(ig,ihubst1)* proj2(ibnd,ihubst1) )  
                       !                             
                       dvqi(ig,ibnd) = dvqi(ig,ibnd)+ 0.5d0 * aux1(ig) 
                       !                             
                    ENDDO 
                    !
                    ! USPP case
                    !
                    IF (okvan) THEN
                       !   
                       DO ig = 1, npwq
                          !                                
                          aux1(ig) = dqsphi(ig,ihubst1) * proj1(ibnd,ihubst1) + & 
                                     swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst1) 
                          !                             
                          dvqi_orth(ig,ibnd) = dvqi_orth(ig,ibnd) + 0.5d0 * aux1(ig)  
                          !   
                          dvqi_orth_lm(ig,ibnd) = dvqi_orth_lm(ig,ibnd) + 0.5d0 * aux1(ig) 
                          !                             
                       ENDDO 
                       !   
                    ENDIF
                    !   
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       !   
                       ihubst2 = offsetU(nah) + m2
                       !   
                       DO ig = 1, npwq          
                          !      
                          aux1(ig) = rho%ns(m1,m2,current_spin,nah) * &
                                     ( dqsphi(ig,ihubst1) * proj1(ibnd,ihubst2)      + &
                                       dqsphi(ig,ihubst2) * proj1(ibnd,ihubst1)      + &
                                       swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst2) + &
                                       swfcatomkpq(ig,ihubst2) * proj2(ibnd,ihubst1) )
                          !   
                          dvqi(ig,ibnd) = dvqi(ig,ibnd) - aux1(ig)   
                          !   
                       ENDDO
                       !
                       ! USPP case
                       !
                       IF (okvan) THEN
                          !     
                          DO ig = 1, npwq   
                             !   
                             aux1(ig) = rho%ns(m1,m2,current_spin,nah) * &
                                        ( dqsphi(ig,ihubst2) * proj1(ibnd,ihubst1) + &
                                          swfcatomkpq(ig,ihubst2) * proj2(ibnd,ihubst1) )
                             !   
                             dvqi_orth(ig,ibnd) = dvqi_orth(ig,ibnd) - aux1(ig) 
                             !   
                             aux1(ig) = rho%ns(m1,m2,current_spin,nah) * &
                                        ( dqsphi(ig,ihubst1) * proj1(ibnd,ihubst2) + & 
                                          swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst2) )  
                             !
                             ! The conjg will be taken in dynmat_hub_scf
                             !   
                             dvqi_orth_lm(ig,ibnd) = dvqi_orth_lm(ig,ibnd) - aux1(ig) 
                             !   
                          ENDDO
                          !
                       ENDIF 
                       !
                    ENDDO ! m2
                    !
                 ENDDO ! m1
                 !
              ENDDO ! ibnd
              !   
              ! effU = Hubbard_U - Hubbard_J0
              ! 
              dvqi = dvqi * effU(nt)
              !   
              DO ig = 1, npwq
                 dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
              ENDDO
              !
              ! USPP case
              !   
              IF (okvan) THEN
                 !   
                 dvqi_orth    = dvqi_orth * effU(nt)
                 dvqi_orth_lm = dvqi_orth_lm * effU(nt)
                 !   
                 DO ig = 1, npwq
                    dvqhbar_orth(ig,:,icart,na)    = dvqhbar_orth(ig,:,icart,na)    + &
                                                     dvqi_orth(ig,:)
                    dvqhbar_orth_lm(ig,:,icart,na) = dvqhbar_orth_lm(ig,:,icart,na) + &
                                                     dvqi_orth_lm(ig,:)
                 ENDDO
                 !   
              ENDIF 
              !   
           ENDIF
           ! 
           ! For Hubbard_J0 
           !
           IF (Hubbard_J0(nt).NE.0.d0) THEN
              !      
              dvqi         = (0.d0, 0.d0)
              dvqi_orth    = (0.d0, 0.d0)
              dvqi_orth_lm = (0.d0, 0.d0)
              !   
              DO ibnd = 1, nbnd
                 !      
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    !   
                    ihubst1 = offsetU(nah) + m1 
                    !
                    ! No diagonal term for J0
                    !   
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       !   
                       ihubst2 = offsetU(nah) + m2
                       !   
                       DO ig = 1, npwq
                          !          
                          aux1(ig) = rho%ns(m1,m2,op_spin,nah) * &
                                     ( dqsphi(ig,ihubst1) * proj1(ibnd,ihubst2)      + &
                                       dqsphi(ig,ihubst2) * proj1(ibnd,ihubst1)      + &
                                       swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst2) + &
                                       swfcatomkpq(ig,ihubst2) * proj2(ibnd,ihubst1) )
                          !  
                          ! Note the sign change w.r.t. the case above
                          ! 
                          dvqi(ig,ibnd) = dvqi(ig,ibnd) + aux1(ig)
                          !   
                       ENDDO
                       !  
                       ! USPP case
                       ! 
                       IF (okvan) THEN
                          !   
                          DO ig = 1, npwq   
                             !    
                             aux1(ig) = rho%ns(m1,m2,op_spin,nah) * &
                                        ( dqsphi(ig,ihubst2) * proj1(ibnd,ihubst1) + &
                                          swfcatomkpq(ig,ihubst2) * proj2(ibnd,ihubst1) )
                             !   
                             dvqi_orth(ig,ibnd) = dvqi_orth(ig,ibnd) + aux1(ig) ! sign change
                             !   
                             aux1(ig) = rho%ns(m1,m2,op_spin,nah) * &
                                        ( dqsphi(ig,ihubst1) * proj1(ibnd,ihubst2) + &   
                                          swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst2) )
                             !          
                             ! The conjg will be taken in dynmat_hub_scf
                             !   
                             dvqi_orth_lm(ig,ibnd) = dvqi_orth_lm(ig,ibnd) + aux1(ig) ! sign change
                             !   
                          ENDDO
                          !   
                       ENDIF
                       !   
                    ENDDO ! m2
                    !   
                 ENDDO ! m1
                 !   
              ENDDO ! ibnd
              !   
              dvqi = dvqi * Hubbard_J0(nt)
              !   
              DO ig = 1, npwq
                 dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
              ENDDO
              !
              ! USPP case   
              !   
              IF (okvan) THEN
                 !   
                 dvqi_orth = dvqi_orth * Hubbard_J0(nt)
                 dvqi_orth_lm = dvqi_orth_lm * Hubbard_J0(nt)
                 !   
                 DO ig = 1, npwq
                    dvqhbar_orth(ig,:,icart,na)    = dvqhbar_orth(ig,:,icart,na)    + &
                                                     dvqi_orth(ig,:)
                    dvqhbar_orth_lm(ig,:,icart,na) = dvqhbar_orth_lm(ig,:,icart,na) + &
                                                     dvqi_orth_lm(ig,:)
                 ENDDO
                 !   
              ENDIF
              !   
           ENDIF
           !
        ENDDO ! nah
        !   
     ENDDO ! icart
     !
  ENDDO ! na
  !   
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dwfcatom_)
  DEALLOCATE (dvqi)
  DEALLOCATE (dvqi_orth)
  DEALLOCATE (dvqi_orth_lm)
  DEALLOCATE (aux1)
  DEALLOCATE (aux2)
  !
  CALL stop_clock( 'dvqhub_barepsi_us2' )
  !   
  RETURN
  !   
END SUBROUTINE dvqhub_barepsi_us2
!----------------------------------------------------------------------------------   
