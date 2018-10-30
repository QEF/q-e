!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE dpsi_orth (ik, wgg, dpsi_orth_cart) 
  !-----------------------------------------------------------------------
  ! 
  ! DFPT+U: This routine calculates for USPP, for each k point, due to
  ! the othogonality contraints, the vector at k+q :
  !
  ! |dpsi_orth_cart(na,icart,ibnd,ispin,ig)> = 
  !     \sum_{n'} wgg(ibnd,n',k) * |psi(n',k+q,ispin)> 
  !           * \sum_{l1,l2} [ <psi(n',k+q,ispin)| d_{na icart}beta(k+q,L,l1)> 
  !                             * qq_nt(L,l1,l2) * <beta(k,L,l2)| psi(n,k,ispin)> 
  !                          + <psi(n',k+q,ispin)| beta(k+q,L,l1)>  
  !                             * qq_nt(L,l1,l2) * <d_{na,icart}beta(k,L,l2)| psi(n,k,ispin)> ]                                      !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,         ONLY : DP
  USE units_lr,      ONLY : iuwfc, lrwfc
  USE ions_base,     ONLY : nat
  USE klist,         ONLY : ngk, igk_k
  USE wvfct,         ONLY : npwx, nbnd
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE wavefunctions, ONLY : evc
  USE eqv,           ONLY : evq
  USE uspp,          ONLY : vkb
  USE ldaU_ph,       ONLY : vkbkpq, dvkb, dvkbkpq
  USE buffers,       ONLY : get_buffer
  USE doubleprojqq_module  
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: ik   
  REAL(DP),    INTENT(IN)    :: wgg(nbnd,nbnd,nksq)
  COMPLEX(DP), INTENT(INOUT) :: dpsi_orth_cart(npwx,nbnd,3,nat)
  !
  ! Local variables
  !
  COMPLEX(DP), ALLOCATABLE :: dpqq(:), dpqq1(:), sum_dpqq(:,:,:,:)
  INTEGER :: ikk, ikq, npw, npwq, i, j, k, icart, jcart, &
             nt, na, nap, l, n, ibnd, jbnd, is
  ! 
  CALL start_clock( 'dpsi_orth' )
  !
  ALLOCATE (dpqq(nbnd))
  ALLOCATE (dpqq1(nbnd))
  ALLOCATE (sum_dpqq(nat,3,nbnd,nbnd))
  !  
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  ! Read the unperturbed KS wave functions at k+q from file
  ! 
  IF (nksq.GT.1) CALL get_buffer (evq, lrwfc, iuwfc, ikq)
  ! 
  sum_dpqq = (0.d0, 0.d0) 
  !
  ! Calculate for all bands jbnd (in evq) and for the band ibnd:
  ! dpqq  = \sum_{l1,l2} <psi(n',k+q,ispin) | d_{nap,jcart}beta(k+q,nap,l1)> * 
  !                      qq_nt(nap,l1,l2) * <beta(k,nap,l2) | psi(n,k)>
  ! dpqq1 = \sum_{l1,l2} <psi(n',k+q,ispin) | beta(k+q,nap,l1)> *
  !                      qq_nt(nap,l1,l2) * <d_{nap,jcart}beta(k,nap,l2) | psi(n,k)>
  ! Note that the derivatives of beta functions are calculated in dvhub_barepsi_us2.
  ! 
  DO nap = 1, nat
     DO jcart = 1, 3 
        DO ibnd = 1, nbnd
           !      
           CALL doubleprojqq (nap, evq, dvkbkpq(:,:,jcart), vkb, evc(:,ibnd), &
                              npwq, npw, dpqq)
           !
           CALL doubleprojqq (nap, evq, vkbkpq, dvkb(:,:,jcart), evc(:,ibnd), &
                              npwq, npw, dpqq1)
           ! 
           sum_dpqq(nap,jcart,ibnd,:) = dpqq + dpqq1
           !
        ENDDO
     ENDDO
  ENDDO
  ! 
  dpsi_orth_cart = (0.d0, 0.d0)
  !
  ! Finally, calculate dpsi_orth_cart
  !
  DO nap = 1, nat
     DO jcart = 1, 3 
        DO ibnd = 1, nbnd
           DO jbnd = 1, nbnd
              dpsi_orth_cart(:,ibnd,jcart,nap) = dpsi_orth_cart(:,ibnd,jcart,nap) + &
                     wgg(ibnd,jbnd,ik) * sum_dpqq(nap,jcart,ibnd,jbnd) * evq(:,jbnd)
           ENDDO 
        ENDDO 
     ENDDO 
  ENDDO 
  !
  DEALLOCATE (dpqq)
  DEALLOCATE (dpqq1)
  DEALLOCATE (sum_dpqq)
  !  
  CALL stop_clock( 'dpsi_orth' )
  !  
  RETURN
  !
END SUBROUTINE dpsi_orth

