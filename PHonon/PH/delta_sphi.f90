!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------------------
SUBROUTINE delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk_, wfcatomkpq_,   & 
                       sdwfcatomk_, sdwfcatomkpq_, vkb_, vkbkpq_, dvkb_, dvkbkpq_, &
                       dqsphi, dmqsphi, iflag)  
  !---------------------------------------------------------------------------------
  !
  ! DFPT+U: This routine calculates a vector at k :
  !
  ! |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) > = S_{k} | \d^{icart}phi_(k,nah,m) > + 
  !       \sum{l1,l2} [ | \dbeta^{icar_}_(k,na_,l1) > * qq_nt(na_, l1 ,l2) * 
  !                           < \beta_(k+q ,na_,l2) | phi_(k+q,nah,m)> + 
  !                      | \beta_(k,na_,l1)> * qq_nt(na_, l1 ,l2) * 
  !                           < \dbeta^{icar_}_(k+q,na_,l2) | phi_(k+q,nah,m) > ]
  !
  ! and also a vector at k+q :
  !
  ! |\Delta_q(S_{k} \phi_(k,I,m)) > = S_{k+q}| \d^{icart}phi_(k+q,nah,m) > + 
  !       \sum{l1,l2} [ | \dbeta^{icar_}_(k+q,na_,l1) > * qq_nt(na_, l1 ,l2) * 
  !                           < \beta_(k ,na_,l2) | phi_(k,nah,m)> + 
  !                     | \beta_(k+q,na_,l1)> * qq_nt(na_, l1 ,l2) * 
  !                           < \dbeta^{icar_}_(k,na_,l2) | phi_(k ,nah,m) > ]
  !
  !  iflag = 1 : calculate |\Delta_q(S_{k} \phi_(k,I,m)) >  AND 
  !                        |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) > 
  !  iflag = 0 : calculate ONLY |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) > 
  !
  ! Written by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  ! 
  USE kinds,      ONLY : DP
  USE uspp_param, ONLY : nh, nhm
  USE ions_base,  ONLY : nat, ityp
  USE uspp,       ONLY : nkb, qq_nt, okvan
  USE ldaU,       ONLY : nwfcU
  USE wvfct,      ONLY : npwx
  USE mp_pools,   ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum 
  USE klist,      ONLY : ngk
  USE io_global,  ONLY : stdout
  USE control_lr, ONLY : ofsbeta 
  !  
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ikk, ikq, na, icart, nah, ihubst 
  ! index of k point
  ! index of k+q point
  ! index of displaced atom
  ! index of cartesian direction of displacement
  ! nah identifies the Hubbard atom I 
  ! index of (I,m) atomic function to which the Dq is applied 
  !
  COMPLEX(DP), INTENT(IN) :: wfcatomk_(npwx,nwfcU),     & 
                             sdwfcatomk_(npwx,nwfcU),   &
                             wfcatomkpq_(npwx,nwfcU),   & 
                             sdwfcatomkpq_(npwx,nwfcU), &
                             vkb_(npwx,nkb),            &
                             vkbkpq_(npwx,nkb),         &
                             dvkb_(npwx,nkb),           &
                             dvkbkpq_(npwx,nkb)
  COMPLEX(DP), INTENT(INOUT) :: dqsphi(npwx,nwfcU), &
                                dmqsphi(npwx,nwfcU)
  INTEGER, INTENT(IN) :: iflag
  !
  ! Local variables
  !
  INTEGER :: nt, ih, m3, m4, ig, l, npw, npwq
  COMPLEX(DP), ALLOCATABLE :: sc1(:), sc2(:), aux1(:), aux2(:)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  ! 
  CALL start_clock( 'delta_sphi' )
  !
  ALLOCATE (sc1(nhm))
  ALLOCATE (sc2(nhm))
  ALLOCATE (aux1(npwx))
  ALLOCATE (aux2(npwx))
  !
  npw  = ngk(ikk)
  npwq = ngk(ikq)
  !
  nt = ityp(na)  
  !  
  ! Calculation of |\Delta_q(S_{k} \phi_(k,I,m)) > 
  ! 
  IF (iflag == 1) THEN
     ! 
     aux1 = (0.d0, 0.d0)
     aux2 = (0.d0, 0.d0)
     !
     ! USPP case 
     !
     IF ( okvan ) THEN
        !
        ! Scalar products in the m3 m4 sum  
        !
        DO ih = 1, nh(nt)
           sc1(ih) = ZDOTC (npw, vkb_(:,ih+ofsbeta(na)),  1, wfcatomk_(:,ihubst), 1)
           sc2(ih) = ZDOTC (npw, dvkb_(:,ih+ofsbeta(na)), 1, wfcatomk_(:,ihubst), 1)
        ENDDO
        !
        CALL mp_sum(sc1, intra_pool_comm)
        CALL mp_sum(sc2, intra_pool_comm)
        ! 
     ENDIF
     ! 
     ! Add to Dq the term |S_{k+q} d_^(na,icart)\phi_(k+q,I,m) > * dkroneker Ina
     ! 
     IF (nah==na) THEN
        DO ig = 1, npwq
           dqsphi(ig,ihubst) = dqsphi(ig,ihubst) + sdwfcatomkpq_(ig,ihubst)
        ENDDO
     ENDIF
     !
     ! USPP case 
     !
     IF ( okvan ) THEN
        DO m3 = 1, nh(nt)
           DO m4 = 1, nh(nt)  
              DO ig = 1, npwq
                 aux1(ig) = dvkbkpq_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc1(m4)
                 aux2(ig) =  vkbkpq_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc2(m4)
                 dqsphi(ig,ihubst) = dqsphi(ig,ihubst) + aux1(ig) + aux2(ig)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
  ENDIF
  !
  ! Calculation of |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
  !
  IF (iflag == 0 .OR. iflag == 1) THEN
     !  
     aux1 = (0.d0, 0.d0)
     aux2 = (0.d0, 0.d0)  
     !
     ! USPP case
     ! 
     IF ( okvan ) THEN  
        !
        ! Scalar products in the m3 m4 sum
        !
        DO ih = 1, nh(nt)
           sc1(ih) = ZDOTC (npwq, vkbkpq_(:,ih+ofsbeta(na)),  1, wfcatomkpq_(:,ihubst), 1)
           sc2(ih) = ZDOTC (npwq, dvkbkpq_(:,ih+ofsbeta(na)), 1, wfcatomkpq_(:,ihubst), 1)
        ENDDO
        !
        CALL mp_sum(sc1, intra_pool_comm)
        CALL mp_sum(sc2, intra_pool_comm)
        ! 
     ENDIF
     !
     ! Add to D-q the term |S_{k} d_^(na,icart)\phi_(k,I,m) > * dkroneker Ina
     !
     IF (nah==na) THEN
        DO ig = 1, npw
           dmqsphi(ig,ihubst) = dmqsphi(ig,ihubst) + sdwfcatomk_(ig,ihubst)
        ENDDO
     ENDIF
     !
     ! USPP case
     !
     IF ( okvan ) THEN
        DO m3 = 1, nh(nt)
           DO m4 = 1, nh(nt)
              DO ig = 1, npw  
                 aux1(ig) = dvkb_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc1(m4)
                 aux2(ig) =  vkb_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc2(m4)
                 dmqsphi(ig,ihubst) = dmqsphi(ig,ihubst) + aux1(ig) + aux2(ig)
                 !
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
  ELSE 
     CALL errore ("delta_sphi"," wrong iflag", 1)
  ENDIF
  !  
  DEALLOCATE (sc1)
  DEALLOCATE (sc2)
  DEALLOCATE (aux1)
  DEALLOCATE (aux2)  
  !
  CALL stop_clock( 'delta_sphi' )
  !
  RETURN
  !
END SUBROUTINE delta_sphi
