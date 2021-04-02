!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dvqhub_barepsi_us (ik, uact)
  !-----------------------------------------------------------------------
  !
  ! DFPT+U 
  ! This routines calculates the BARE derivative of the Hubbard potential times psi.
  ! is  = current_spin
  ! isi = opposite of the current_spin 
  !
  ! |Delta V_BARE_(k+q,is) psi(ibnd,k,is)> = 
  !     + \sum_(I,m1,m2) Hubbard_U(I) * [0.5\delta_(m1,m2) - ns(m1,m2,is,I)] *
  !                         { |dqsphi(imode,I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)> + 
  !                           |S\phi(I,k+q,m1)><dmqsphi(imode,I,k,m2)|psi(ibnd,k,is)> } 
  !     - \sum_(I,m1,m2) Hubbard_U(I) * dnsbare(m1,m2,is,I,imode) *
  !                           |S\phi(I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)>
  !
  ! Addition of the J0 terms:
  !
  !     + \sum_(I,m1,m2) Hubbard_J0(I) * ns(m1,m2,isi,I) *
  !                         { |dqsphi(imode,I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)> + 
  !                           |S\phi(I,k+q,m1)><dmqsphi(imode,I,k,m2)|psi(ibnd,k,is)> } 
  !     + \sum_(I,m1,m2) Hubbard_J0(I) * dnsbare(m1,m2,isi,I,imode) *
  !                           |S\phi(I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)>
  !
  ! Important: in this routine vkb is a beta function at k+q, and vkb_ is beta at k.
  ! This is done so because vkb is calculated at k+q in solve_linter (i.e. before calling
  ! this routine), so we keep the same attribution here.
  ! 
  ! Written by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !                 
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode
  USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : U_projection, Hubbard_l, is_hubbard, Hubbard_J0, offsetU, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, swfcatomk, swfcatomkpq, dwfcatomkpq, &
                            sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
                            proj1, proj2, dnsbare, effU 
  USE wvfct,         ONLY : npwx, nbnd
  USE uspp,          ONLY : vkb,nkb 
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
  COMPLEX(DP), INTENT(IN) :: uact(3*nat)
  ! the pattern of displacements
  !
  ! Local variables
  !
  INTEGER :: i, j, k, icart, counter, na, nt, l, ih, n, mu, ig, &
             ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, op_spin, &
             ikk, ikq, npw, npwq, ibeta
  COMPLEX(DP), ALLOCATABLE :: aux1(:), aux2(:), aux3(:), aux4(:), aux5(:), &
                              dqsphi(:,:), dmqsphi(:,:), dvqi(:,:), dvqhbar(:,:,:,:), &
                              vkb_(:,:), dwfcatom_(:)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !  
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  ALLOCATE (aux1(npwx))
  ALLOCATE (aux2(npwx))
  ALLOCATE (aux3(npwx))
  ALLOCATE (aux4(npwx)) 
  ALLOCATE (aux5(npwx))
  ALLOCATE (dqsphi(npwx,nwfcU))
  ALLOCATE (dmqsphi(npwx,nwfcU))
  ALLOCATE (dvqi(npwx,nbnd))
  ALLOCATE (dvqhbar(npwx,nbnd,3,nat))
  ALLOCATE (vkb_(npwx,nkb))
  ALLOCATE (dwfcatom_(npwx))
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
     ENDIF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Compute the beta function at k and put the result in vkb_
  !
  IF (.NOT.lgamma) THEN   
     CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb_)
  ELSE
     vkb_ = vkb
  ENDIF
  !
  ! The beta function at k+q. Let us put it in the proper array vkbkpq,
  ! because in the following of this routine the array vkb will be 
  ! used as a workspace in the routine swfc. 
  !
  vkbkpq = vkb
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
                      vkb_(:,ibeta), dvkb(:,ibeta,icart))
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
  dvqhbar = (0.d0, 0.d0)
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
           ! For Hubbard_U - Hubbard_J0
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
                    ! Apply the S matrix: | S d_^(I,icart)\phi_(k,I,m) >
                    !
                    CALL swfc (npw, 1, vkb_, dwfcatom_, sdwfcatomk(:,ihubst))
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
                 CALL delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq,  &
                                  sdwfcatomk, sdwfcatomkpq, vkb_, vkbkpq, dvkb(:,:,icart), &
                                  dvkbkpq(:,:,icart), dqsphi, dmqsphi, 1)  
                 !
                 ! Calculate:
                 ! proj1 (ihubst,ibnd) = < S_{k}\phi_(k,I,m)| psi(ibnd,k) >
                 ! proj2 (ihubst,ibnd) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k) > 
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
        CALL mp_sum (proj1, intra_bgrp_comm) 
        CALL mp_sum (proj2, intra_bgrp_comm)
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           dvqi = (0.d0, 0.d0) 
           !
           IF (is_hubbard(nt)) THEN
              !
              DO ibnd = 1, nbnd
                 !
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    !
                    ihubst1 = offsetU(nah) + m1 
                    !
                    DO ig = 1, npwq
                       !
                       aux1(ig) = dqsphi(ig,ihubst1) * proj1(ibnd,ihubst1) 
                       !
                       aux3(ig) = swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst1)  
                       !
                       dvqi(ig,ibnd) = dvqi(ig,ibnd) + 0.5d0 * (aux1(ig)+aux3(ig))
                       !
                    ENDDO
                    !
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       !
                       ihubst2 = offsetU(nah) + m2
                       !
                       DO ig = 1, npwq 
                          !                         
                          aux2(ig) = dqsphi(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          aux4(ig) = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                     * proj2(ibnd, ihubst2)
                          aux5(ig) = swfcatomkpq(ig,ihubst1) &
                                     * dnsbare(m1,m2,current_spin,nah,icart,na) &
                                     * proj1(ibnd, ihubst2)
                          !
                          dvqi(ig, ibnd) = dvqi(ig, ibnd) - aux2(ig) - aux4(ig) - aux5(ig) 
                          ! 
                       ENDDO
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
           ENDIF
           !
           ! Hubbard_J0 \= 0 case 
           !
           dvqi = (0.d0, 0.d0) 
           ! 
           IF (Hubbard_J0(nt).NE.0.d0) THEN
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
                          aux2(ig) = dqsphi(ig, ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          aux4(ig) = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                     * proj2(ibnd, ihubst2)
                          aux5(ig) = swfcatomkpq(ig,ihubst1) &
                                     * dnsbare (m1,m2,op_spin,nah,icart,na) &
                                     * proj1 (ibnd, ihubst2)
                          !
                          ! Note the sign change w.r.t. the case above
                          !
                          dvqi(ig, ibnd) = dvqi(ig, ibnd) + aux2(ig) + aux4(ig) + aux5(ig)
                          ! 
                       ENDDO
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
           ENDIF
           !
        ENDDO ! nah
        !
     ENDDO ! icart
     ! 
  ENDDO ! na
  ! 
  ! Compute the displacement along the pattern uact (for each band).
  ! The result is stored in aux1.
  !
  DO ibnd = 1, nbnd    
     !   
     aux1 = (0.d0, 0.d0)
     !
     DO na = 1, nat
        mu = 3 * (na - 1)
        ! Here is the basis transformation from cartesian to pattern uact
        DO icart = 1, 3 
           DO ig = 1, npwq
              aux1(ig) = aux1(ig) + dvqhbar(ig,ibnd,icart,na) * uact(mu+icart) 
           ENDDO
        ENDDO
     ENDDO
     !
     ! Add the result to dvpsi
     !
     DO ig = 1, npwq
        dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + aux1(ig)
     ENDDO
     !
  ENDDO 
  !
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)  
  DEALLOCATE (aux1) 
  DEALLOCATE (aux2)
  DEALLOCATE (aux3) 
  DEALLOCATE (aux4)
  DEALLOCATE (aux5)
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dvqi)
  DEALLOCATE (dvqhbar)
  DEALLOCATE (vkb_)
  DEALLOCATE (dwfcatom_) 
  !
  RETURN
  !
END SUBROUTINE dvqhub_barepsi_us
