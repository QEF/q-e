!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------------
MODULE doubleprojqq_module
!
CONTAINS
!
!----------------------------------------------------------------------------------
SUBROUTINE doubleprojqq (na, vec1, vec2, vec3, vec4, npw1, npw2, dpqq) 
   !--------------------------------------------------------------------------------
   !
   ! This routine calculates for all ibnd: 
   ! dpqq(ibnd) = \sum{l1 l2} < vec1(ibnd)  | vec2(na,l1) > * qq(na, l1 ,l2) * &
   !                          < vec3(na,l2) | vec4 >
   ! 
   USE kinds,       ONLY : DP
   USE uspp_param,  ONLY : nh
   USE ions_base,   ONLY : ityp
   USE uspp,        ONLY : qq_nt
   USE wvfct,       ONLY : npwx, nbnd
   USE mp_pools,    ONLY : intra_pool_comm
   USE mp,          ONLY : mp_sum
   USE control_lr,  ONLY : ofsbeta
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: na
   ! index of the displaced atom
   COMPLEX(DP), INTENT(IN) :: vec1(:,:), & ! (npwx,nkb)
                              vec2(:,:), & ! (npwx,nkb)
                              vec3(:,:), & ! (npwx,nkb)
                              vec4(:)      ! (npwx)
   INTEGER, INTENT(IN) :: npw1, npw2
   COMPLEX(DP), INTENT(OUT) :: dpqq(:)     ! (nbnd)
   !
   ! Local variables
   !
   INTEGER :: nt, l1, l2, ibeta1, ibeta2, ibnd
   COMPLEX(DP) :: projauxvec4
   COMPLEX(DP), ALLOCATABLE :: aux1(:), projvec1vec2(:)
   COMPLEX(DP), EXTERNAL :: ZDOTC
   ! 
   CALL start_clock ( 'doubleprojqq' )
   !
   ALLOCATE (aux1(npwx))
   ALLOCATE (projvec1vec2(nbnd))
   !
   dpqq = (0.d0, 0.d0)
   !
   nt = ityp(na)  
   !
   DO l1 = 1, nh(nt)
      ! 
      ibeta1 = ofsbeta(na) + l1
      !     
      ! Calculate: projvec1vec2(ibnd) = < vec1(ibnd) | vec2 > for each l1 
      !
      DO ibnd = 1, nbnd
         projvec1vec2(ibnd) = ZDOTC (npw1, vec1(:,ibnd), 1, vec2(:,ibeta1), 1)
      ENDDO
      !
#if defined(__MPI)
      CALL mp_sum(projvec1vec2, intra_pool_comm)
#endif
      !
      aux1 =  (0.d0, 0.d0)
      !
      ! aux1 = \sum_l2 qq_nt(l1,l2,nt) * |vec3_(na,l2)>
      !
      DO l2 = 1, nh(nt)
         ibeta2 = ofsbeta(na) + l2
         aux1(:) = aux1(:) + qq_nt(l1,l2,nt) * vec3(:,ibeta2)    
      ENDDO
      !
      ! Calculate projauxvec4 = < aux1 | vec4 >
      !
      projauxvec4 = ZDOTC (npw2, aux1, 1, vec4, 1)
      !
#if defined(__MPI)
      CALL mp_sum(projauxvec4, intra_pool_comm)
#endif
      !
      ! Summing on l1 for each band ibnd
      !
      dpqq(:) = dpqq(:) + projvec1vec2(:) * projauxvec4
      !
   ENDDO
   !
   DEALLOCATE (aux1)
   DEALLOCATE (projvec1vec2)
   !
   CALL stop_clock ( 'doubleprojqq' )
   !
   RETURN
   !
END SUBROUTINE doubleprojqq
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE doubleprojqq2 (na, proj, vec3, vec4, npw2, dpqq) 
   !
   ! This routine calculates for all ibnd:
   ! dpqq(ibnd) = \sum{l1 l2} proj(ibnd,na,l1) * qq_nt(na, l1 ,l2) * &
   !                          < vec3 (na,l2) | vec4 >
   !   
   USE kinds,      ONLY : DP
   USE uspp_param, ONLY : nh
   USE ions_base,  ONLY : ityp
   USE uspp,       ONLY : qq_nt
   USE wvfct,      ONLY : npwx, nbnd
   USE mp_pools,   ONLY : intra_pool_comm
   USE mp,         ONLY : mp_sum
   USE control_lr, ONLY : ofsbeta
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: na
   ! index of the displaced atom
   COMPLEX(DP), INTENT(IN) :: proj(:,:), & ! (nbnd,nkb)
                              vec3(:,:), & ! (npwx,nkb)
                              vec4 (:)     ! (npwx)
   INTEGER, INTENT (IN) :: npw2
   COMPLEX(DP), INTENT(OUT) :: dpqq(:)     ! (nbnd)
   !
   ! Local variables
   !
   INTEGER :: nt, l1, l2, ibeta1, ibeta2, ibnd
   COMPLEX(DP), ALLOCATABLE :: aux1(:)
   COMPLEX(DP) ::  projauxvec4
   COMPLEX(DP), EXTERNAL :: ZDOTC
   !
   CALL start_clock ( 'doubleprojqq2' )
   ! 
   ALLOCATE (aux1(npwx))
   !
   dpqq = (0.d0, 0.d0)
   !
   nt = ityp(na)  
   !
   DO l1 = 1, nh(nt)
      ! 
      ibeta1 = ofsbeta(na) + l1
      !
      aux1 = (0.d0, 0.d0)
      !
      DO l2 = 1, nh(nt)
         ibeta2 = ofsbeta(na) + l2
         aux1(:) = aux1(:) + qq_nt(l1,l2,nt) * vec3(:,ibeta2)    
      ENDDO
      !
      ! Calculate projauxvec4 = < aux1 | vec4 >
      !
      projauxvec4 = ZDOTC (npw2, aux1, 1, vec4, 1)
      !
#if defined(__MPI)
      CALL mp_sum(projauxvec4, intra_pool_comm)
#endif
      !
      ! Summing over l1 for each band ibnd
      !
      dpqq(:) = dpqq(:) + proj(:,ibeta1) * projauxvec4 
      !
   ENDDO
   !
   DEALLOCATE (aux1)
   !
   CALL stop_clock ( 'doubleprojqq2' )
   !
   RETURN
   !
END SUBROUTINE doubleprojqq2
!-------------------------------------------------------------- 

END MODULE doubleprojqq_module
!--------------------------------------------------------------


!--------------------------------------------------------
MODULE term_one_1_module
!--------------------------------------------------------  
  USE mp_pools,   ONLY : intra_pool_comm
  USE mp,         ONLY:  mp_sum  
!  
CONTAINS
!
!--------------------------------------------------------  
SUBROUTINE term_one_1 (ik, icart, jcart, evc_, &
                       wfcatom_, proj_, vkb_, resone_1)
    !----------------------------------------------------
    !
    USE kinds,      ONLY : DP
    USE wvfct,      ONLY : npwx, nbnd, wg 
    USE uspp,       ONLY : vkb, nkb 
    USE klist,      ONLY : ngk, igk_k
    USE qpoint,     ONLY : ikks
    USE doubleprojqq_module  
    !
    IMPLICIT NONE
    !
    INTEGER,     INTENT(IN) :: ik, icart, jcart
    COMPLEX(DP), INTENT(IN) :: evc_(:,:),   & ! (npwx,nbnd)
                               wfcatom_(:), & ! (npwx)
                               proj_(:),    & ! (nbnd)
                               vkb_(:,:)      ! (npwx,nkb)
    COMPLEX(DP), INTENT(INOUT) :: resone_1
    !
    ! Local variables
    !
    INTEGER     :: npw, ikk, ibnd
    COMPLEX(DP), ALLOCATABLE :: d2wfcatomk(:), sd2wfcatomk(:), projd2(:)
    COMPLEX(DP), EXTERNAL :: ZDOTC
    !
    ALLOCATE(d2wfcatomk(npwx))
    ALLOCATE(sd2wfcatomk(npwx))
    ALLOCATE(projd2(nbnd))
    !
    resone_1 = (0.d0, 0.d0) 
    !
    ikk = ikks(ik)    
    npw = ngk(ikk)
    !
    ! Calculate the 2nd derivative of the atomic orbitals:
    ! | d2_^(I,icart,jcart) \phi_(k,I,m) >
    !
    CALL d2wfc (npw, igk_k(1,ikk), ikk, icart, jcart, &
                wfcatom_, d2wfcatomk) 
    ! 
    ! Apply the S operator to the result above:
    ! | S d2_^(I,icart,jcart) \phi_(k,I,m) >
    !
    CALL swfc (npw, 1, vkb_, d2wfcatomk, sd2wfcatomk)  
    !    
    ! Calculate projd2(ibnd) = < psi(ibnd,k) | S d2_^(I,icart,jcart) \phi_(k,I,m) > 
    ! at ihubst1 (i.e. I m) 
    !
    DO ibnd = 1, nbnd
       projd2(ibnd) = ZDOTC (npw, evc_(:,ibnd), 1, sd2wfcatomk, 1)
    ENDDO
    !
    CALL mp_sum(projd2, intra_pool_comm)
    !
    DO ibnd = 1, nbnd
       resone_1 = resone_1 + wg(ibnd,ikk) * projd2(ibnd) * proj_(ibnd)
    ENDDO
    !
    DEALLOCATE(d2wfcatomk)
    DEALLOCATE(sd2wfcatomk)
    DEALLOCATE(projd2)   
    !
    RETURN
    !
END SUBROUTINE term_one_1
!------------------------------------------------------------------------------
  
END MODULE term_one_1_module
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
MODULE term_three_1_module
!  
CONTAINS
!
!-------------------------------------------------------------------------------  
SUBROUTINE  term_three_1 (ik, icart, jcart, ihubst1, ihubst2, &
                          projdphi, resthree_1)
    !---------------------------------------------------------------------------
    !
    USE kinds,      ONLY : DP
    USE wvfct,      ONLY : wg, nbnd
    USE qpoint,     ONLY : ikks
    !  
    IMPLICIT NONE
    !
    INTEGER,     INTENT(IN) :: ik, icart, jcart, ihubst1, ihubst2 
    COMPLEX(DP), INTENT(IN) :: projdphi(:,:,:) ! (nbnd,nwfcU,3)
    COMPLEX(DP), INTENT(INOUT) :: resthree_1
    !
    ! Local variables
    !
    INTEGER :: ibnd, ikk
    !
    resthree_1 = (0.d0, 0.d0)
    !
    ikk = ikks(ik)
    !
    DO ibnd = 1, nbnd
       resthree_1 = resthree_1 + wg(ibnd,ikk) * projdphi(ibnd,ihubst1,icart) * &
                                 CONJG(projdphi (ibnd,ihubst2, jcart))
    ENDDO
    ! 
    RETURN
    !
END SUBROUTINE term_three_1
!----------------------------------------------------------------------------------
  
END MODULE term_three_1_module
!----------------------------------------------------------------------------------

!-----------------------------------------------------------------------
MODULE term_one_module
!-----------------------------------------------------------------------
  USE doubleprojqq_module
!  
CONTAINS
!
!------------------------------------------------------------------------  
SUBROUTINE term_one (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2, &
                     evc_, wfcatomk, swfcatomk, vkb_, vkbkpq_, dvkb_,    &
                     dvkbkpq_, dwfcatomkpq_, res_one) 
    !--------------------------------------------------------------------
    ! 
    USE kinds,      ONLY : DP
    USE uspp,       ONLY : nkb, okvan
    USE wvfct,      ONLY : npwx, nbnd, wg 
    USE uspp_param, ONLY : nh
    USE ions_base,  ONLY : ityp
    USE control_lr, ONLY : ofsbeta
    USE ldaU_ph,    ONLY : proj1, projpb, projpdb
    USE klist,      ONLY : ngk, igk_k
    USE qpoint,     ONLY : ikks, ikqs
    USE doubleprojqq_module
    USE term_one_1_module
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik, icart, jcart, na, nap, nah, ihubst1, ihubst2
    COMPLEX(DP), INTENT(IN) :: evc_(:,:),        & ! (npwx,nbnd)
                               wfcatomk(:,:),    & ! (npwx,nwfcU)
                               swfcatomk(:,:),   & ! (npwx,nwfcU)
                               vkb_(:,:),        & ! (npwx,nkb)
                               vkbkpq_(:,:),     & ! (npwx,nkb)
                               dvkb_(:,:,:),     & ! (npwx,nkb,3)
                               dvkbkpq_(:,:,:),  & ! (npwx,nkb,3)
                               dwfcatomkpq_(:,:,:) ! (npwx,nwfcU,3)
    COMPLEX(DP), INTENT(INOUT) :: res_one
    !
    ! Local variables
    !
    INTEGER     :: npw, npwq, ikk, ikq, ibnd, nt, l1, l2, l, ibeta
    COMPLEX(DP) :: resone_1, resone_2,resone_3, resone_4, resone_5, &
                   resone_6_9
    COMPLEX(DP), ALLOCATABLE :: dpqq(:), dpqq1(:), dpqq2(:), &
                                dpqq3(:), dpqq4(:), d2vkb(:,:)
    COMPLEX(DP), EXTERNAL :: ZDOTC  
    !
    res_one = 0.d0
    !
    ALLOCATE(dpqq(nbnd))
    ALLOCATE(dpqq1(nbnd))
    ALLOCATE(dpqq2(nbnd))
    ALLOCATE(dpqq3(nbnd))
    ALLOCATE(dpqq4(nbnd))
    ALLOCATE(d2vkb(npwx,nkb))
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    npw = ngk(ikk)
    npwq= ngk(ikq)
    !
    IF ((na==nap) .AND. (nah==na)) THEN    
       ! term_one_1 contains a delta_na_nap 
       !
       ! Calculate term_one_1
       !
       CALL term_one_1 (ik, icart, jcart, evc_, wfcatomk(:,ihubst1), &
                        proj1(:,ihubst2), vkb_, resone_1)
       !
       res_one = res_one + resone_1
       !
    ENDIF
    !
    ! USPP case
    !
    IF (okvan) THEN  
       !       
       IF (nah==nap) THEN
          ! 
          ! Calculate term_one_2
          !
          resone_2 = (0.d0, 0.d0) 

          CALL doubleprojqq2 (na, projpdb(:,:,icart), vkbkpq_, &
                              dwfcatomkpq_(:,ihubst1,jcart), npwq, dpqq)
          !          
          DO ibnd = 1, nbnd
             resone_2 = resone_2 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + CONJG(resone_2)
          !
          ! Calculate term_one_3
          ! 
          resone_3 = (0.d0, 0.d0) 
          !
          CALL doubleprojqq2 (na, projpb, dvkbkpq_(:,:,icart), &
                              dwfcatomkpq_(:,ihubst1,jcart), npwq, dpqq)
          !          
          DO ibnd = 1, nbnd
             resone_3 = resone_3 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + CONJG(resone_3)
          !
       ENDIF 
       !       
       IF (nah==na) THEN
          !
          ! Calculate term_one_4
          !
          resone_4 = (0.d0, 0.d0) 
          !
          CALL doubleprojqq2 (nap, projpdb(:,:,jcart), vkbkpq_, &
                              dwfcatomkpq_(:,ihubst1,icart), npwq, dpqq)
          !
          DO ibnd = 1, nbnd
             resone_4 = resone_4 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + resone_4
          !
          ! Calculate term_one_5
          ! 
          resone_5 = (0.d0, 0.d0)
          !
          CALL doubleprojqq2 (nap, projpb, dvkbkpq_(:,:,jcart), &
                              dwfcatomkpq_(:,ihubst1,icart), npwq, dpqq)
          !          
          DO ibnd = 1, nbnd
             resone_5 = resone_5 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + resone_5
          !
       ENDIF  
       !
       IF (na==nap) THEN
          !   
          resone_6_9 = (0.d0, 0.d0)
          !
          d2vkb = (0.d0, 0.d0)
          !
          nt = ityp(na)
          !
          DO l = 1, nh(nt)
             !
             ibeta = ofsbeta(na) + l
             !     
             ! Calculate the 2nd derivative of the beta functions for 
             ! all l states of atom na 
             ! 
             CALL d2wfc (npw, igk_k(1,ikk), ikk, icart, jcart, &
                         vkb_(:,ibeta), d2vkb(:,ibeta)) 
             !
             ! d2vkb is always the 2nd derivative at icart and jcart, 
             ! displacing the atom na and looking at the beta of atom j=na_
             !
          ENDDO
          !         
          ! doubleprojqq, unlike doubleprojqq2, calculates the first proj1 inside
          ! 
          CALL doubleprojqq  (na, evc_, d2vkb, vkb_, wfcatomk(:,ihubst1), &
                              npw, npw,dpqq1)     
          !
          CALL doubleprojqq2 (na, projpb, d2vkb, wfcatomk(:,ihubst1),     &
                              npw,dpqq2)
          !
          CALL doubleprojqq2 (na, projpdb(:,:,icart), dvkb_(:,:,jcart),  &
                              wfcatomk(:,ihubst1), npw, dpqq3)
          !
          CALL doubleprojqq2 (na, projpdb(:,:,jcart), dvkb_(:,:,icart),  &
                              wfcatomk(:,ihubst1), npw,dpqq4)
          !
          DO ibnd = 1, nbnd
             resone_6_9 = resone_6_9 + ( dpqq1(ibnd) + dpqq2(ibnd) +   &
                                         dpqq3(ibnd) + dpqq4(ibnd) ) * &
                                       proj1(ibnd,ihubst2) * wg(ibnd,ikk) 
          ENDDO
          ! 
          res_one = res_one + resone_6_9
          !
       ENDIF
       ! 
    ENDIF 
    !
    DEALLOCATE(dpqq)
    DEALLOCATE(dpqq1)
    DEALLOCATE(dpqq2)
    DEALLOCATE(dpqq3)
    DEALLOCATE(dpqq4)
    DEALLOCATE(d2vkb) 
    !
    RETURN
    !
END SUBROUTINE term_one
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
SUBROUTINE term_one_diag (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2, &
                          evc_, wfcatomk, swfcatomk, vkb_, vkbkpq_, dvkb_,   &
                          dvkbkpq_, dwfcatomkpq_, res_one) 
    !------------------------------------------------------------------------
    ! 
    USE kinds,      ONLY : DP
    USE uspp,       ONLY : nkb, okvan
    USE wvfct,      ONLY : npwx, nbnd, wg 
    USE uspp_param, ONLY : nh
    USE ions_base,  ONLY : ityp
    USE control_lr, ONLY : ofsbeta
    USE ldaU_ph,    ONLY : proj1, projpb, projpdb
    USE klist,      ONLY : ngk, igk_k
    USE qpoint,     ONLY : ikks, ikqs
    USE doubleprojqq_module
    USE term_one_1_module
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik, icart, jcart, na, nap, nah, ihubst1, ihubst2
    COMPLEX(DP), INTENT(IN) :: evc_(:,:),        & ! (npwx,nbnd)
                               wfcatomk(:,:),    & ! (npwx,nwfcU)
                               swfcatomk(:,:),   & ! (npwx,nwfcU)
                               vkb_(:,:),        & ! (npwx,nkb)
                               vkbkpq_(:,:),     & ! (npwx,nkb)
                               dvkb_(:,:,:),     & ! (npwx,nkb,3)
                               dvkbkpq_(:,:,:),  & ! (npwx,nkb,3)
                               dwfcatomkpq_(:,:,:) ! (npwx,nwfcU,3)
    COMPLEX(DP), INTENT(INOUT) :: res_one
    !
    ! Local variables
    !
    INTEGER :: npw, npwq, ikk, ikq, ibnd, nt, l1, l2, l, ibeta
    COMPLEX(DP), ALLOCATABLE :: dpqq(:), dpqq1(:), dpqq2(:), dpqq3(:), dpqq4(:), d2vkb(:,:)
    COMPLEX(DP) :: resone_1, resone_2,resone_3, resone_4, resone_5, &
                   resone_6_9
    COMPLEX(DP), EXTERNAL :: ZDOTC  
    !
    res_one = 0.d0
    !
    ALLOCATE(dpqq(nbnd))
    ALLOCATE(dpqq1(nbnd))
    ALLOCATE(dpqq2(nbnd))
    ALLOCATE(dpqq3(nbnd))
    ALLOCATE(dpqq4(nbnd))
    ALLOCATE(d2vkb(npwx,nkb))
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    npw = ngk(ikk)
    npwq= ngk(ikq)
    !
    ! In the diagonal approximation J=HUBBARD_I, 
    ! all terms are such that na=nap=nah
    !    
    IF ((na==nap) .AND. (nah==na)) THEN    
       ! term_one_1 contains a delta_na_nap 
       ! 
       ! Calculate term_one_1
       !
       CALL term_one_1 (ik, icart, jcart, evc_, wfcatomk(:,ihubst1), &
                        proj1(:,ihubst2), vkb_, resone_1)
       !
       res_one = res_one + resone_1
       !
       ! USPP case
       !
       IF (okvan) THEN  
          !          
          ! Calculate term_one_2
          !
          resone_2 = (0.d0, 0.d0) 
          !
          CALL doubleprojqq2 (na, projpdb(:,:,icart), vkbkpq_, &
                              dwfcatomkpq_(:,ihubst1,jcart), npwq, dpqq)
          !          
          DO ibnd = 1, nbnd
             resone_2 = resone_2 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + CONJG(resone_2)
          ! 
          ! Calculate term_one_3
          ! 
          resone_3 = (0.d0, 0.d0) 
          !
          CALL doubleprojqq2 (na, projpb, dvkbkpq_(:,:,icart), & 
                              dwfcatomkpq_(:,ihubst1,jcart), npwq, dpqq)
          !          
          DO ibnd = 1, nbnd
             resone_3 = resone_3 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + CONJG(resone_3)
          !
          ! Calculate term_one_4
          !
          resone_4 = (0.d0, 0.d0) 
          !
          CALL doubleprojqq2 (nap, projpdb(:,:,jcart), vkbkpq_, &
                              dwfcatomkpq_(:,ihubst1,icart), npwq, dpqq)
          !
          DO ibnd = 1, nbnd
             resone_4 = resone_4 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + resone_4
          ! 
          ! Calculate term_one_5
          !
          resone_5 = (0.d0, 0.d0)
          !
          CALL doubleprojqq2 (nap, projpb, dvkbkpq_(:,:,jcart), &
                              dwfcatomkpq_(:,ihubst1,icart), npwq, dpqq)
          !          
          DO ibnd = 1, nbnd
             resone_5 = resone_5 + wg(ibnd,ikk) * dpqq(ibnd) * proj1(ibnd,ihubst2)
          ENDDO
          !
          res_one = res_one + resone_5
          !
          resone_6_9 = (0.d0, 0.d0)
          !
          d2vkb = (0.d0, 0.d0)
          !
          nt = ityp(na)
          !
          DO l = 1,  nh(nt)
             !
             ibeta = ofsbeta(na) + l
             !     
             ! Calculate the 2nd derivative of the beta functions 
             ! for all l states of atom na 
             ! 
             CALL d2wfc (npw, igk_k(1,ikk), ikk, icart, jcart, &
                         vkb_(:,ibeta), d2vkb(:,ibeta)) 
             ! 
             ! d2vkb is always the 2nd derivative at icart and jcart, 
             ! displacing the atom na and looking at the beta of atom j=na_
             !
          ENDDO
          !          
          ! doubleprojqq, unlike doubleprojqq2, calculates the first proj1 inside
          !
          CALL doubleprojqq  (na, evc_, d2vkb, vkb_, wfcatomk(:,ihubst1), &
                              npw, npw, dpqq1)     
          !
          CALL doubleprojqq2 (na, projpb, d2vkb, wfcatomk(:,ihubst1),     &
                              npw, dpqq2)
          !
          CALL doubleprojqq2 (na, projpdb(:,:,icart), dvkb_(:,:,jcart),  & 
                              wfcatomk(:,ihubst1), npw, dpqq3)
          !
          CALL doubleprojqq2 (na, projpdb(:,:,jcart), dvkb_(:,:,icart),  &
                              wfcatomk(:,ihubst1), npw, dpqq4)
          !
          DO ibnd = 1, nbnd
             resone_6_9 = resone_6_9 + ( dpqq1(ibnd) + dpqq2(ibnd) +   &
                                         dpqq3(ibnd) + dpqq4(ibnd) ) * &
                                       proj1(ibnd,ihubst2) * wg(ibnd,ikk) 
          ENDDO
          !
          res_one = res_one + resone_6_9
          !
       ENDIF
       !
    ENDIF
    ! 
    DEALLOCATE(dpqq)
    DEALLOCATE(dpqq1)
    DEALLOCATE(dpqq2)
    DEALLOCATE(dpqq3)
    DEALLOCATE(dpqq4)
    DEALLOCATE(d2vkb)   
    !
    RETURN
    !
END SUBROUTINE term_one_diag
!------------------------------------------------------------------------  
  
END MODULE term_one_module
!------------------------------------------------------------------------

!-------------------------------------------------------------------------
MODULE term_three_module
!-------------------------------------------------------------------------
  USE mp_pools,   ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum 
!
CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE term_three (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2, &
                       evc_, wfcatomk, dwfcatomk, vkb_, dvkb_, wfcatomkpq, &
                       vkbkpq_, dvkbkpq_, res_three) 
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE uspp,       ONLY : vkb, nkb, okvan
  USE wvfct,      ONLY : npwx, nbnd, wg 
  USE ldaU,       ONLY : nwfcU
  USE ldaU_ph,    ONLY : projpb, projpdb
  USE klist,      ONLY : ngk, igk_k
  USE qpoint,     ONLY : ikks, ikqs
  USE doubleprojqq_module
  USE term_three_1_module
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik, icart, jcart, na, nap, nah, ihubst1, ihubst2
  COMPLEX(DP), INTENT(IN) :: evc_(:,:),        & ! (npwx,nbnd)
                             wfcatomk(:,:),    & ! (npwx,nwfcU)
                             dwfcatomk(:,:,:), & ! (npwx,nwfcU,3)
                             vkb_(:,:),        & ! (npwx,nkb)
                             dvkb_(:,:,:),     & ! (npwx,nkb,3)
                             wfcatomkpq(:,:),  & ! (npwx,nwfcU)
                             vkbkpq_(:,:),     & ! (npwx,nkb)
                             dvkbkpq_(:,:,:)     ! (npwx,nkb,3)
  COMPLEX(DP), INTENT(INOUT) :: res_three
  !
  ! Local variables
  !
  INTEGER :: npw, npwq, ikk, ikq, icar, ibnd
  COMPLEX(DP) :: resthree_1, resthree_2, resthree_3, resthree_4
  COMPLEX(DP), ALLOCATABLE :: sdwfcatomk(:,:,:), projdphi(:,:,:), &
                              dpqq1(:), dpqq2(:), dpqq3(:),     &
                              dpqq4(:), aux(:), aux2(:)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  ALLOCATE (sdwfcatomk(npwx,nwfcU,3)) 
  ALLOCATE (projdphi(nbnd,nwfcU,3)) 
  ALLOCATE (dpqq1(nbnd)) 
  ALLOCATE (dpqq2(nbnd)) 
  ALLOCATE (dpqq3(nbnd)) 
  ALLOCATE (dpqq4(nbnd)) 
  ALLOCATE (aux(nbnd)) 
  ALLOCATE (aux2(nbnd)) 
  !
  res_three=0.d0
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  DO icar = 1, 3
     IF ((icar==icart) .OR. (icar==jcart)) THEN 
        ! we want only icart jcart
        ! 
        ! Calculate | S d^{icar} \phi_(k,I,m) >
        !
        CALL swfc (npw, 1, vkb_, dwfcatomk(:,ihubst1,icar), sdwfcatomk(:,ihubst1,icar))
        !
        ! Calculate | S d^{icar} \phi_(k,I,m') >  
        !
        CALL swfc (npw, 1, vkb_, dwfcatomk(:,ihubst2,icar), sdwfcatomk(:,ihubst2,icar))
        !
        ! Calculate projdphi(ibnd) = < psi(inbd,k) | S d_^(I,icart) \phi_(k,I,m) > 
        ! at ihubst1 (i.e. I m). 
        ! 
        DO ibnd = 1, nbnd
           projdphi(ibnd, ihubst1, icar) = &
                ZDOTC (npw, evc_(:,ibnd), 1, sdwfcatomk(:,ihubst1,icar), 1)
           projdphi(ibnd, ihubst2, icar) = &
                ZDOTC (npw, evc_(:,ibnd), 1, sdwfcatomk(:,ihubst2,icar), 1)
        ENDDO
        ! 
     ENDIF 
  ENDDO  
  !
  CALL mp_sum(projdphi, intra_pool_comm)
  ! 
  ! Calculate term_three_1
  !
  IF ((na==nap) .AND. (nah==na)) THEN
     !
     resthree_1 = (0.d0, 0.d0)
     !
     CALL term_three_1 (ik, icart, jcart, ihubst1, ihubst2, projdphi, resthree_1)
     !
     res_three = res_three + resthree_1
     !
  ENDIF
  !
  ! USPP case
  !
  IF (okvan) THEN
     !     
     ! Calculate term_three_2
     !
     resthree_2 = (0.d0, 0.d0)  
     !
     CALL doubleprojqq2 (na, projpdb(:,:,icart), vkbkpq_, wfcatomkpq(:,ihubst1), &
                         npwq, dpqq1)
     !
     CALL doubleprojqq2 (na, projpb, dvkbkpq_(:,:,icart), wfcatomkpq(:,ihubst1), &
                         npwq, dpqq2)
     ! 
     aux = dpqq1 + dpqq2
     !
     CALL doubleprojqq2 (nap, projpdb(:,:,jcart), vkbkpq_, wfcatomkpq(:,ihubst2), &
                         npwq, dpqq3)
     !
     CALL doubleprojqq2 (nap, projpb, dvkbkpq_(:,:,jcart), wfcatomkpq(:,ihubst2), &
                         npwq, dpqq4)
     !
     aux2 = dpqq3 + dpqq4
     !
     DO ibnd = 1, nbnd
        resthree_2 =  resthree_2 + wg(ibnd,ikk) * CONJG(aux(ibnd)) * aux2(ibnd)
     ENDDO
     !
     res_three = res_three + resthree_2
     !
     ! Calculate term_three_3
     !
     IF (nah == na) THEN
        !
        resthree_3 = (0.d0, 0.d0)
        !
        ! Calculate \sum {l1 l2} [ < psi | \beta(k,na_,l1) > qq(na_, l1 ,l2) * &
        !                          < \dbeta^jcar(k+q,na_,l2) | phi_(k+q,nah,m) > ] 
        !
        CALL doubleprojqq2 (nap, projpb, dvkbkpq_(:,:,jcart), wfcatomkpq(:,ihubst2), &
                            npwq, dpqq1)
        !
        ! Calculate \sum {l1 l2} [ < psi| \dbeta^jcar_(k,na_,l1)> qq(na_, l1 ,l2) * &
        !                          < \beta_(k+q,na_,l2) | phi_(k+q,nah,m) > ]
        !
        CALL doubleprojqq2 (nap, projpdb(:,:,jcart), vkbkpq_, wfcatomkpq(:,ihubst2), &
                            npwq, dpqq2)
        !
        DO ibnd = 1, nbnd
           resthree_3 = resthree_3 + wg(ibnd,ikk) * (dpqq1(ibnd)+dpqq2(ibnd)) * &
                                     CONJG(projdphi(ibnd,ihubst1,icart))
        ENDDO
        ! 
        res_three = res_three + resthree_3
        !
     ENDIF
     ! 
     ! Calculate term_three_4 
     !
     IF (nah == nap) THEN
        !   
        resthree_4 = (0.d0, 0.d0)     
        !
        CALL doubleprojqq2 (na, projpb, dvkbkpq_(:,:,icart), wfcatomkpq(:,ihubst1), &
                            npwq, dpqq2)
        !
        CALL doubleprojqq2 (na, projpdb(:,:,icart), vkbkpq_, wfcatomkpq(:,ihubst1), &
                            npwq, dpqq1)
        !
        DO ibnd = 1, nbnd
           resthree_4 = resthree_4  +  wg(ibnd,ikk) * (CONJG(dpqq1(ibnd))+CONJG(dpqq2(ibnd))) * &
                                       projdphi(ibnd,ihubst2,jcart)
        ENDDO
        ! 
        res_three = res_three + resthree_4
        !
     ENDIF
     !
  ENDIF 
  !
  DEALLOCATE (sdwfcatomk)
  DEALLOCATE (projdphi) 
  DEALLOCATE (dpqq1)            
  DEALLOCATE (dpqq2) 
  DEALLOCATE (dpqq3) 
  DEALLOCATE (dpqq4) 
  DEALLOCATE (aux)  
  DEALLOCATE (aux2)
  !
  RETURN
  !
END SUBROUTINE term_three
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
SUBROUTINE term_three_diag (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2, &
                            evc_, wfcatomk, dwfcatomk, vkb_, dvkb_, wfcatomkpq, &
                            vkbkpq_, dvkbkpq_, res_three) 
  !------------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE uspp,       ONLY : vkb, nkb, okvan
  USE wvfct,      ONLY : npwx, nbnd, wg 
  USE ldaU,       ONLY : nwfcU
  USE ldaU_ph,    ONLY : projpb, projpdb
  USE klist,      ONLY : ngk, igk_k
  USE qpoint,     ONLY : ikks, ikqs
  USE doubleprojqq_module
  USE term_three_1_module

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik, icart, jcart, na, nap, nah, ihubst1, ihubst2
  COMPLEX(DP), INTENT(IN) :: evc_(:,:),        & ! (npwx,nbnd)
                             wfcatomk(:,:),    & ! (npwx,nwfcU)
                             dwfcatomk(:,:,:), & ! (npwx,nwfcU,3)
                             vkb_(:,:),        & ! (npwx,nkb)
                             dvkb_(:,:,:),     & ! (npwx,nkb,3)
                             wfcatomkpq(:,:),  & ! (npwx,nwfcU)
                             vkbkpq_(:,:),     & ! (npwx,nkb)
                             dvkbkpq_(:,:,:)     ! (npwx,nkb,3)
  COMPLEX(DP), INTENT(INOUT) :: res_three
  !
  ! Local variables
  !
  INTEGER     :: ikk, ikq, npw, npwq, icar, ibnd
  COMPLEX(DP) :: resthree_1, resthree_2, resthree_3, resthree_4
  COMPLEX(DP), ALLOCATABLE :: sdwfcatomk(:,:,:), projdphi(:,:,:),     &
                              dpqq1(:), dpqq2(:), dpqq3(:), dpqq4(:), &
                              aux(:), aux2(:)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  ALLOCATE (sdwfcatomk(npwx,nwfcU,3))
  ALLOCATE (projdphi(nbnd,nwfcU,3)) 
  ALLOCATE (dpqq1(nbnd))
  ALLOCATE (dpqq2(nbnd))
  ALLOCATE (dpqq3(nbnd))
  ALLOCATE (dpqq4(nbnd))
  ALLOCATE (aux(nbnd))
  ALLOCATE (aux2(nbnd))
  !
  res_three = 0.d0
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  DO icar = 1, 3
     IF ((icar == icart) .OR. (icar == jcart)) THEN 
        ! we want only icart jcart
        !
        ! Calculate | S d^{icar} \phi_(k,I,m) >
        !
        CALL swfc (npw, 1, vkb_, dwfcatomk(:,ihubst1,icar), sdwfcatomk(:,ihubst1,icar))
        !
        ! Calculate | S d^{icar} \phi_(k,I,m') >
        !  
        CALL swfc (npw, 1, vkb_, dwfcatomk(:,ihubst2,icar), sdwfcatomk(:,ihubst2,icar))
        !
        ! Calculate projdphi(ibnd) = < \psi(inbd,k) | S d_^(I,icart) \phi_(k,I,m) > 
        ! at ihubst1 (i.e. I m). 
        !
        DO ibnd = 1, nbnd
           projdphi(ibnd, ihubst1, icar) = &
                & ZDOTC (npw, evc_(:,ibnd), 1, sdwfcatomk(:,ihubst1,icar), 1)
           projdphi(ibnd, ihubst2, icar) = &
                & ZDOTC (npw, evc_(:,ibnd), 1, sdwfcatomk(:,ihubst2,icar), 1)
        ENDDO
        !
     ENDIF 
  ENDDO  
  !
  CALL mp_sum(projdphi, intra_pool_comm)
  !
  IF ((na==nap) .AND. (nah==na)) THEN  
     !
     ! Calculate term_three_1
     !
     resthree_1 = (0.d0, 0.d0)
     !
     CALL term_three_1 (ik, icart, jcart, ihubst1, ihubst2, projdphi, resthree_1)
     !
     res_three = res_three + resthree_1
     !
     ! USPP case
     !
     IF (okvan) THEN
        !     
        ! Calculate term_three_2
        !
        resthree_2 = (0.d0, 0.d0)  
        !
        CALL doubleprojqq2 (na, projpdb(:,:,icart), vkbkpq_, wfcatomkpq(:,ihubst1), &
                            npwq, dpqq1)
        !
        CALL doubleprojqq2 (na, projpb, dvkbkpq_(:,:,icart), wfcatomkpq(:,ihubst1), &
                            npwq, dpqq2)
        !
        aux = dpqq1 + dpqq2
        !
        CALL doubleprojqq2 (nap, projpdb(:,:,jcart), vkbkpq_, wfcatomkpq(:,ihubst2), &
                            npwq, dpqq3)
        !
        CALL doubleprojqq2 (nap, projpb, dvkbkpq_(:,:,jcart), wfcatomkpq(:,ihubst2), &
                            npwq, dpqq4)
        !
        aux2 = dpqq3 + dpqq4
        !
        DO ibnd = 1, nbnd
           resthree_2 = resthree_2 + wg(ibnd,ikk) * CONJG(aux(ibnd)) * aux2(ibnd)
        ENDDO
        !
        res_three = res_three + resthree_2
        ! 
        ! Calculate term_three_3  
        !
        resthree_3 = (0.d0, 0.d0)
        !
        ! Calculate \sum {l1 l2} [ < psi | \beta(k,na_,l1) > qq_nt(na_, l1 ,l2) * &
        !                          < \dbeta^jcar(k+q ,na_,l2) | phi_(k+q,nah,m) > ] 
        !
        CALL doubleprojqq2 (nap, projpb, dvkbkpq_(:,:,jcart), wfcatomkpq(:,ihubst2), &
                            npwq, dpqq1)
        !
        ! Calculate \sum {l1 l2} [ < psi | \dbeta^jcar_(k,na_,l1)> qq_nt(na_, l1 ,l2) * &
        !                          < \beta_(k+q,na_,l2) | phi_(k+q,nah,m) > ]
        !
        CALL doubleprojqq2 (nap, projpdb(:,:,jcart), vkbkpq_, wfcatomkpq(:,ihubst2), &
                            npwq, dpqq2)
        !
        DO ibnd = 1, nbnd
           resthree_3 = resthree_3 + wg(ibnd,ikk) * (dpqq1(ibnd)+dpqq2(ibnd)) * &
                                     CONJG(projdphi(ibnd,ihubst1,icart))
        ENDDO
        !
        res_three = res_three + resthree_3
        !
        ! Calculate term_three_4
        !
        resthree_4 = (0.d0, 0.d0)     
        !
        CALL doubleprojqq2 (na, projpb, dvkbkpq_(:,:,icart), wfcatomkpq(:,ihubst1), &
                            npwq, dpqq2)
        !
        CALL doubleprojqq2 (na, projpdb(:,:,icart), vkbkpq_, wfcatomkpq(:,ihubst1), &
                            npwq, dpqq1)
        !
        DO ibnd = 1, nbnd
           resthree_4 = resthree_4 + wg(ibnd,ikk) * (CONJG(dpqq1(ibnd))+CONJG(dpqq2(ibnd))) * &
                                     projdphi (ibnd,ihubst2,jcart)
        ENDDO
        !
        res_three = res_three + resthree_4
        !
     ENDIF 
     ! 
  ENDIF
  ! 
  DEALLOCATE (sdwfcatomk)
  DEALLOCATE (projdphi)
  DEALLOCATE (dpqq1)
  DEALLOCATE (dpqq2)
  DEALLOCATE (dpqq3)
  DEALLOCATE (dpqq4)
  DEALLOCATE (aux)
  DEALLOCATE (aux2) 
  ! 
  RETURN
  !
END SUBROUTINE term_three_diag
!--------------------------------------------------------

END MODULE term_three_module
!--------------------------------------------------------

!------------------------------------------------------------
MODULE d2nsq_bare_module
!------------------------------------------------------------
!
CONTAINS
!  
!------------------------------------------------------------
SUBROUTINE d2nsq_bare_k (ik, icart, jcart, na, nap, nah, &
                         & ihubst1, ihubst2, d2ns_bare_k) 
    !---------------------------------------------------------
    ! 
    ! DFPT+U: This routine calculates the second bare derivative 
    !         of the occupation matrix ns. ns is derived 
    !         two times w.r.t. the atomic positions, using the 
    !         unperturbed wfc's
    !
    ! Written  by A. Floris
    ! Modified by I. Timrov (01.10.2018)
    ! 
    USE kinds,         ONLY : DP
    USE units_lr,      ONLY : iuwfc, lrwfc
    USE ions_base,     ONLY : nat, ityp, ntyp => nsp
    USE klist,         ONLY : xk, ngk, igk_k
    USE ldaU_ph,       ONLY : wfcatomk, swfcatomk, wfcatomkpq, dwfcatomk, dwfcatomkpq, &
                              dvkb, vkbkpq, dvkbkpq, proj1, d2ns_type
    USE wvfct,         ONLY : npwx, nbnd, wg
    USE uspp,          ONLY : vkb, nkb
    USE qpoint,        ONLY : nksq, ikks, ikqs
    USE control_lr,    ONLY : lgamma
    USE uspp_param,    ONLY : nh
    USE lsda_mod,      ONLY : lsda, isk, nspin
    USE io_global,     ONLY : stdout
    USE wavefunctions, ONLY : evc
    USE term_one_module
    USE term_three_module
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik, icart, jcart, na, nap, nah, ihubst1, ihubst2 
    COMPLEX(DP), INTENT(INOUT) :: d2ns_bare_k
    ! k point index 
    ! cartesian component
    ! cartesian component
    ! displaced atom index
    ! displaced atom index
    ! hubbard atom index
    ! atomic state
    ! atomic state
    ! second bare derivative of the occupation matrix
    !
    ! Local variables
    !
    INTEGER  :: ikk, ikq, npw, npwq, icar, &
                nt, ic, nti, ina, ih, ibeta, ibnd 
    COMPLEX(DP) :: res_one, res_two, res_three, res_four 
    COMPLEX(DP), EXTERNAL :: ZDOTC 
    !
    CALL start_clock( 'd2nsq_bare_k' )
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    npw = ngk(ikk)
    npwq= ngk(ikq)
    !
    d2ns_bare_k = (0.d0, 0.d0)
    !
    ! Calculate the derivatives of atomic functions at k and k+q
    !
    DO icar = 1, 3
       !
       IF ( (icar == icart) .OR. (icar == jcart)) THEN 
          ! we want only icart jcart
          !
          IF ((nah == na) .OR. (nah == nap)) THEN 
             ! we want only dphi at na or nap
             !
             ! Calculate |d_icart\phi_(k,I,m))>        
             !              
             CALL dwfc (npw, igk_k(1,ikk), ikk, icar, &
                        wfcatomk(:,ihubst1), dwfcatomk(:,ihubst1,icar)) 
             !
             ! Calculate |d_icart\phi_(k,I,m')>        
             !
             CALL dwfc (npw, igk_k(1,ikk), ikk, icar, &
                        wfcatomk(:,ihubst2), dwfcatomk(:,ihubst2,icar)) 
             !
             IF (.NOT.lgamma) THEN
                !
                ! Calculate |d_icart\phi_(k+q,I,m))>
                !
                CALL dwfc (npwq, igk_k(1,ikq), ikq, icar, &
                           wfcatomkpq(:,ihubst1), dwfcatomkpq(:,ihubst1,icar)) 
                !                 
                ! Calculate |d_icart\phi_(k+q,I,m'))>
                !
                CALL dwfc (npwq, igk_k(1,ikq), ikq, icar, &
                           wfcatomkpq(:,ihubst2), dwfcatomkpq(:,ihubst2,icar)) 
                !
             ENDIF
             !
          ENDIF 
          !
       ENDIF
       !
    ENDDO 
    !
    !-------------------- term1 ------------------------------------------
    !
    CALL term_one (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2, &
                   evc, wfcatomk, swfcatomk, vkb, vkbkpq, dvkb,       &
                   dvkbkpq, dwfcatomkpq, res_one) 
    ! 
    d2ns_bare_k = d2ns_bare_k + res_one
    ! 
    IF (d2ns_type == 'fmmp') THEN 
       !
       ! fmmp approximation, we just have a factor of 2 
       ! for both term_one and term_three,
       ! without recalculating term_two and term_four. 
       ! 
       d2ns_bare_k = d2ns_bare_k + res_one
       !
    ELSE
       !
       ! term2 = term1 with the exchange of indices m <=> m' 
       !
       IF (ihubst1==ihubst2) THEN
          d2ns_bare_k  = d2ns_bare_k + res_one
       ELSE
          CALL term_one (ik, icart, jcart, na, nap, nah, ihubst2, ihubst1, &
                         evc, wfcatomk, swfcatomk, vkb, vkbkpq, dvkb,        &
                         dvkbkpq, dwfcatomkpq, res_two)
          d2ns_bare_k = d2ns_bare_k + res_two    
       ENDIF
       ! 
    ENDIF
    !
    !-------------------- term3 ------------------------------------------
    !
    CALL term_three (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2,      &
                     evc, wfcatomk, dwfcatomk, vkb, dvkb, wfcatomkpq, vkbkpq, &
                     dvkbkpq, res_three) 
    ! 
    d2ns_bare_k = d2ns_bare_k + res_three
    !
    IF  (d2ns_type == 'fmmp') THEN
       ! 
       ! fmmp approximation, we just have a factor of 2 
       ! for both term_one and term_three, 
       ! without recalculating term_two and term_four. 
       ! 
       d2ns_bare_k = d2ns_bare_k + res_three
       !
    ELSE
       !
       ! term4 = term3 with the exchange of indices m <=> m' 
       !
       IF (ihubst1==ihubst2) THEN
          d2ns_bare_k = d2ns_bare_k + res_three
       ELSE
          CALL term_three (ik, icart, jcart, na, nap, nah, ihubst2, ihubst1,      &
                           evc, wfcatomk, dwfcatomk, vkb, dvkb, wfcatomkpq, vkbkpq, &
                           dvkbkpq, res_four)
          d2ns_bare_k = d2ns_bare_k + res_four
       ENDIF
       ! 
    ENDIF
    ! 
    CALL stop_clock( 'd2nsq_bare_k' )
    !
    RETURN
    ! 
END SUBROUTINE d2nsq_bare_k
!------------------------------------------------------------------
  
!------------------------------------------------------------------
SUBROUTINE d2nsq_bare_k_diag (ik, icart, jcart, na, nap, nah, &
                              & ihubst1, ihubst2, d2ns_bare_k) 
    !--------------------------------------------------------------
    ! 
    ! DFPT+U: This routine calculates the second bare derivative 
    !         of the occupation matrix ns. ns is derived    
    !         two times w.r.t. the atomic positions, using the 
    !         unperturbed wfc's. 
    !         This routines does an approximate calculation.
    !
    !  d2ns_type='diag': if okvan=.true. the d2ns_bare matrix 
    !                    is calculated retaining only 
    !                    the <\beta_J|\phi_I> products on the 
    !                    same atomic site, i.e. for J==I. 
    !                    WARNING: Check against 'full'   
    !
    !  d2ns_type='dmmp': same as 'diag', but also assuming a m <=> m' 
    !                    symmetry in the various contributions of 
    !                    the d2ns_bare_k matrix. 
    !                    WARNING: Check against 'full'  
    !
    ! Written  by A. Floris
    ! Modified by I. Timrov (01.10.2018)
    ! 
    USE kinds,           ONLY : DP
    USE units_lr,        ONLY : iuwfc, lrwfc
    USE ions_base,       ONLY : nat, ityp, ntyp => nsp
    USE klist,           ONLY : xk, ngk, igk_k
    USE ldaU_ph,         ONLY : wfcatomk, swfcatomk, wfcatomkpq, dwfcatomk, dwfcatomkpq, &
                                dvkb, vkbkpq, dvkbkpq, proj1, d2ns_type
    USE wvfct,           ONLY : npwx, nbnd, wg
    USE uspp,            ONLY : vkb, nkb
    USE qpoint,          ONLY : nksq, ikks, ikqs
    USE control_lr,      ONLY : lgamma
    USE uspp_param,      ONLY : nh
    USE lsda_mod,        ONLY : lsda, isk, nspin
    USE io_global,       ONLY : stdout
    USE wavefunctions,   ONLY : evc
    USE term_one_module
    USE term_three_module
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik, icart, jcart, na, nap, nah, ihubst1, ihubst2 
    COMPLEX(DP), INTENT(INOUT) :: d2ns_bare_k 
    ! k point index
    ! cartesian component
    ! cartesian component
    ! displaced atom index
    ! displaced atom index
    ! hubbard atom index
    ! atomic state
    ! atomic state
    ! second bare derivative of the occupation matrix
    !
    ! Local variables
    !
    INTEGER :: ikk, ikq, npw, npwq, icar, nt, ic, nti, ina, ih, ibeta, ibnd 
    COMPLEX(DP) :: res_one,res_two,res_three,res_four 
    COMPLEX(DP), EXTERNAL :: ZDOTC 
    !
    CALL start_clock( 'd2nsq_bare_k_diag' )
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    npw = ngk(ikk)
    npwq= ngk(ikq)
    !
    d2ns_bare_k = (0.d0, 0.d0)
    !
    ! Calculate the derivatives of atomic functions at k and k+q
    !
    DO icar = 1, 3
       !
       IF ((icar==icart) .OR. (icar==jcart)) THEN 
          ! Wwe want only icart jcart
          ! In the diagonal approximation J=HUBBARD_I
          !
          IF ((na == nah).AND. (nap == nah)) THEN 
             ! we want only NA=NAP=NAH
             ! we want only dphi at na or nap
             !
             ! Calculate |d_icart\phi_(k,I,m))>        
             !              
             CALL dwfc (npw, igk_k(1,ikk), ikk, icar, &
                        wfcatomk(:,ihubst1), dwfcatomk(:,ihubst1,icar)) 
             !
             ! Calculate |d_icart\phi_(k,I,m')>        
             !
             CALL dwfc (npw, igk_k(1,ikk), ikk, icar, &
                        wfcatomk(:,ihubst2), dwfcatomk(:,ihubst2,icar)) 
             ! 
             IF (.NOT.lgamma) THEN
                !
                ! Calculate |d_icart\phi_(k+q,I,m))>
                !
                CALL dwfc (npwq, igk_k(1,ikq), ikq, icar, &
                           wfcatomkpq(:,ihubst1), dwfcatomkpq(:,ihubst1,icar)) 
                !                 
                ! calculate |d_icart\fi_(k+q,I,m')) >
                !
                CALL dwfc (npwq, igk_k(1,ikq), ikq, icar, &
                           wfcatomkpq(:,ihubst2), dwfcatomkpq(:,ihubst2,icar)) 
                !
             ENDIF
             !       
          ENDIF
          !          
       ENDIF
       !
    ENDDO
    !
    !-------------------- term1 ---------------------------------------------
    !
    CALL term_one_diag (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2, &
                        evc, wfcatomk, swfcatomk, vkb, vkbkpq, dvkb,        &
                        dvkbkpq, dwfcatomkpq, res_one) 
    !
    d2ns_bare_k = d2ns_bare_k + res_one
    !
    IF (d2ns_type == 'dmmp') THEN 
       ! 
       ! dmmp approximation, we just have a factor of 2 
       ! for both term_one and term_three,
       ! without recalculating term_two and term_four. 
       !
       d2ns_bare_k = d2ns_bare_k + res_one
       !
    ELSE
       !
       ! term2 = term1 with the exchange of indices m <=> m' 
       !
       IF (ihubst1==ihubst2) THEN
          d2ns_bare_k = d2ns_bare_k + res_one
       ELSE
          CALL term_one_diag (ik, icart, jcart, na, nap, nah, ihubst2, ihubst1, &
                              evc, wfcatomk, swfcatomk, vkb, vkbkpq, dvkb,        &
                              dvkbkpq, dwfcatomkpq, res_two)
          d2ns_bare_k  = d2ns_bare_k + res_two    
       ENDIF
       ! 
    ENDIF 
    !
    !-------------------- term3 -------------------------------------------------
    !
    CALL term_three_diag (ik, icart, jcart, na, nap, nah, ihubst1, ihubst2,      &
                          evc, wfcatomk, dwfcatomk, vkb, dvkb, wfcatomkpq, vkbkpq, &
                          dvkbkpq, res_three) 
    !
    d2ns_bare_k = d2ns_bare_k + res_three
    !
    ! term4 = term3 with the exchange of indices m <=> m' 
    !
    IF (d2ns_type == 'dmmp') THEN 
       ! 
       ! dmmp approximation, we just have a factor of 2 
       ! for both term_one and term_three,
       ! without recalculating term_two and term_four. 
       !
       d2ns_bare_k = d2ns_bare_k + res_three
       !
    ELSE
       ! 
       IF (ihubst1==ihubst2) THEN
          d2ns_bare_k = d2ns_bare_k + res_three
       ELSE
          CALL term_three_diag (ik, icart, jcart, na, nap, nah, ihubst2, ihubst1,      &
                                evc, wfcatomk, dwfcatomk, vkb, dvkb, wfcatomkpq, vkbkpq, &
                                dvkbkpq, res_four)
          d2ns_bare_k = d2ns_bare_k + res_four
       ENDIF
       ! 
    ENDIF
    !
    CALL stop_clock( 'd2nsq_bare_k_diag' )
    ! 
    RETURN
    ! 
END SUBROUTINE d2nsq_bare_k_diag
!-----------------------------------------------------------------------
  
END MODULE d2nsq_bare_module
!-----------------------------------------------------------------------
