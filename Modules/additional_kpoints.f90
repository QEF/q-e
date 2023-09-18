!
! Copyright (C) 2020-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Lorenzo Paulatto <paulatz@gmail.com>, July 2020
!
MODULE additional_kpoints
  USE kinds, ONLY : DP
  USE parameters,   ONLY : npk
  IMPLICIT NONE
  REAL(DP),ALLOCATABLE :: xk_add(:,:) !, wk_add(:)
  CHARACTER(len=80) :: k_points_add = 'bogus'
  INTEGER :: nkstot_add=0


  CONTAINS
  !
  SUBROUTINE bcast_additional_kpoints
    USE mp,             ONLY : mp_bcast
    USE io_global,      ONLY : ionode_id
    USE mp_images,      ONLY: intra_image_comm
    !
    IMPLICIT NONE
    CALL mp_bcast(nkstot_add, ionode_id, intra_image_comm)
    IF(nkstot_add==0) RETURN
    CALL mp_bcast(xk_add,     ionode_id, intra_image_comm)

  END SUBROUTINE
  !
  SUBROUTINE add_additional_kpoints(nkstot, xk, wk)
     USE input_parameters, ONLY : nqx1, nqx2, nqx3
     USE cell_base,        ONLY : bg
     USE io_global,        ONLY : stdout
     IMPLICIT NONE
     INTEGER,INTENT(inout)  :: nkstot
     REAL(DP),INTENT(inout) :: xk(3,npk), wk(npk)
     !
     REAL(DP),ALLOCATABLE :: xk_old(:,:), wk_old(:)
     INTEGER :: nkstot_old
     INTEGER :: nk1_old, nk2_old, nk3_old
     INTEGER :: k1_old,  k2_old,  k3_old
     INTEGER :: nqtot, i,j,k, iq, jq
     REAL(DP) :: xq(3), rq(3)
     REAL(DP) :: almost_zero = 1.d-8
     LOGICAL, EXTERNAL  :: matches
     !
!     IF(.not.allocated(xk) .or. .not.allocated(wk))&
!       CALL errore("add_kpoints", "K-points not ready yet",1)
     CALL bcast_additional_kpoints()
     IF(nkstot_add==0) RETURN

     IF(matches("crystal",k_points_add))THEN
         CALL cryst_to_cart(nkstot_add,xk_add,bg,+1)
     ENDIF
     !
     ! Back-up existing points
     nkstot_old = nkstot
     ALLOCATE(xk_old(3,nkstot_old))
     ALLOCATE(wk_old(nkstot_old))
     xk_old = xk(1:3, 1:nkstot)
     wk_old = wk(1:nkstot)
     nkstot = 0
     !
     ! Simple case: EXX not used or used with self-exchange only: 
     IF( nqx1<=1 .and. nqx2<=1 .and. nqx3<=1 ) THEN
       nkstot = nkstot_old + nkstot_add
       IF(nkstot>npk) CALL errore("add_kpoint", "Maximum number of k-points exceeded: increase npk in pwcom", 1)
       xk(:,1:nkstot_old) = xk_old
       xk(:,nkstot_old+1:nkstot_old+nkstot_add) = xk_add
       wk(1:nkstot_old) = wk_old
       wk(nkstot_old+1:nkstot_old+nkstot_add) = almost_zero
       nqtot=1
     ELSE
     ! Difficult case: EXX with a finite grid of q-points. Ideally, we would want to use
     ! The grid from module EXX, but it may not have been initialized at this points.
     ! Furthermore, the q-point grid is obtained by opening the k-points one, so this would
     ! be a dog wagging its own tails
       nqtot = nqx1*nqx2*nqx3
       nkstot = nkstot_old + nkstot_add *nqtot
       IF(nkstot>npk) CALL errore("add_kpoint", "Maximum number of k-points exceeded: increase npk in pwcom", 1)
       xk(:,1:nkstot_old) = xk_old
       wk(1:nkstot_old) = wk_old
       
       rq = (/nqx1,nqx2,nqx3/)
       rq = 1._dp / rq
       iq = nkstot_old
       ! 
       DO  i = 0,nqx1-1
       DO  j = 0,nqx2-1
       DO  k = 0,nqx3-1
         xq = rq*(/i,j,k/)
         CALL cryst_to_cart(1,xq,bg,+1)
         DO jq = 1, nkstot_add
            iq = iq + 1
           xk(:,iq) = xk_add(:,jq)  + xq
           wk(iq)=almost_zero
         ENDDO
       ENDDO
       ENDDO
       ENDDO
     ENDIF
  
     ! renormalize 
     wk(1:nkstot) = wk(1:nkstot)/SUM(wk(1:nkstot))
     
     WRITE(stdout,"(5x,a)")    " --- Additional k-points: --- "
     WRITE(stdout,"(5x,a,i6,a)")    "A request of ",nkstot_add," k-points with zero weight added to list"
     IF(nqtot>1) WRITE(stdout,"(5x,a,i6,a)") "Furthermore,  ",nkstot_add*nqtot, " k-points where added to perform EXX calculation"
!     WRITE(stdout,"(5x,a)")         "They can be extracted with bands.x using:"
     WRITE(stdout,"(5x,a,i6,a,i6)") "first=",nkstot_old+1,",  last=",nkstot_old+nkstot_add
     WRITE(stdout,*)

  END SUBROUTINE
  !
END MODULE
