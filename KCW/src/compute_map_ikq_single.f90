!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!#define DEBUG
!
!-----------------------------------------------------------------------
SUBROUTINE compute_map_ikq_single (iq)
  !-----------------------------------------------------------------------
  !
  !! This routine compute the map to go from (iq,ik) --> ikq_1BZ
  !! In general k+q lies outside the 1BZ, but if k and q belongs to the 1BZ 
  !! there exist a G vector such that k + q = k'+ G with k' in the 1BZ
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, nkstot
  USE io_global,            ONLY : ionode, ionode_id
  USE control_kcw,          ONLY : map_ikq, shift_1bz, kcw_iverbosity
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_image_comm 
  USE io_global,            ONLY : stdout
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : at
  !
#ifdef DEBUG
  USE mp_world,             ONLY : mpime
#endif
  !
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN) :: iq 
  ! The index of the q point
  !
  INTEGER :: ik, ikq
  ! Counter for the k and k+q points in the BZ
  !
  REAL(DP) :: xkq(3), gvect(3), gvect_(3)
  ! the k+q coordinate and the G vector that shift it into the 1 BZ
  !
#ifdef DEBUG
  REAL(DP) :: xk_(3), xq_(3), xkq_(3)
#endif
  !
  CALL start_clock ('map')
  IF ( ALLOCATED (map_ikq) ) DEALLOCATE (map_ikq)
  IF ( ALLOCATED (shift_1bz) ) DEALLOCATE (shift_1bz)
  !
  ALLOCATE ( map_ikq(nkstot) )
  ALLOCATE ( shift_1bz(3,nkstot) )
  map_ikq (:) = 0
  shift_1bz(:,:) = 0.D0 
  !
  IF (ionode) THEN 
    !
    DO ik = 1, nkstot/nspin
      !
      xkq(:) = xk(:,ik)+xk(:,iq)
      ! the k+q coordinate (cart)
      !
#ifdef DEBUG
      xq_(:) = xk(:,iq)
      CALL cryst_to_cart(1, xq_, at, -1)
      xk_(:) = xk(:,ik)
      CALL cryst_to_cart(1, xk_, at, -1)
      WRITE(mpime+100,'("iq, xq", i5, 3f8.4)') iq, xq_(:)
      WRITE(mpime+100,'("ik, xk", i5, 3f8.4)') ik, xk_(:)
#endif
      gvect(:) = 0.D0
      ! the index of the corresponding k point in the IBZ
      !
      !CALL find_index_1bz(xkq, gvect, ikq)
      CALL find_index_1bz_iterate(xkq, gvect, ikq)
      !
      ! ... Store the result in a global variable 
      map_ikq(ik) = ikq
      shift_1bz(:,ik) = gvect(:)
      gvect_(:) = gvect(:)
      CALL cryst_to_cart(1, gvect_, at, -1)
#ifdef DEBUG
      xkq_(:) = xkq(:)
      CALL cryst_to_cart(1, xkq_, at, -1)
      WRITE(mpime+100,'("The map (iq,ik) --> ikq", 2i, 8x, i)') iq,ik, ikq
      WRITE(mpime+100,'("xkq, shift" , 5x, 3f8.4, 5x, 3f8.4)') xkq_(:), gvect_(:)
      WRITE(mpime+100,'("ikq, xkq_1BZ", i5, 3f8.4,/)') ikq, xq_(:) + xk_(:) -gvect_(:)
#endif
      !
      IF (kcw_iverbosity .gt. 1) WRITE(stdout,'(8X, "The map (iq,ik) --> ip + G", 5x, & 
                                                &  " (", 2i3, "  ) " , 5x, i3 , 8x, "+", 3f8.4, " [Cryst]" )') iq, ik, ikq, gvect_
    ENDDO
    !
  ENDIF 
  !
  CALL mp_bcast ( map_ikq, ionode_id, intra_image_comm )
  CALL mp_bcast ( shift_1bz, ionode_id, intra_image_comm )
  !
  IF (kcw_iverbosity .gt. 1) WRITE(stdout, *) ""
  WRITE( stdout, '(8X,"INFO: Map k+q -> p in 1BZ DONE  ",/)') 
  CALL stop_clock ('map')
  !
END subroutine



!-----------------------------------------------------------------------
SUBROUTINE find_index_1bz(xkq, gvect, ikq)
  !-----------------------------------------------------------------------
  !
  ! Obsolote. Does not work properly if k points have negative cooerdinates
  ! Use find_index_1bz_iterate instead. More roboust, but slower.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg
  USE klist,                ONLY : xk, nkstot
  USE lsda_mod,             ONLY : nspin
  !USE mp_world,             ONLY : mpime
  
  !USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xkq(3)
  REAL(DP), INTENT(OUT) :: gvect(3)
  REAL(DP) :: xkq_(3), dist, xk_(3)
  INTEGER :: ikq, ik, cntr
  LOGICAL :: found(nkstot)
  !
  !
  found(:) = .false.
  xkq_(:) = xkq(:)
  ! The k+q in crystal coordinates 
  CALL cryst_to_cart(1, xkq_, at, -1)
  ! The integer part of the k+q (says if is inside the 1BZ (=0) or not (/=0))
  !gvect(:) = REAL( INT( xkq_(:) + 1e-12) ) ! the small addend to prevent 0.9999999999 -> 0 and not 1
  gvect(:) = REAL( INT( xkq_(:) ) )
  !WRITE(100,*) 'xkq', xkq_
  !WRITE(100,*) 'shift', gvect
  ! The k+q shifted inside the 1BZ. 
  ! NB: Not sure if this works when xkq is negative (not sure either if it can be negative, tough) 
  xkq_(:) = xkq_(:) - gvect(:)
  !WRITE(100,*) 'xkq_1BZ', xkq_
  !
  cntr = 0
  DO ik = 1, nkstot/nspin
    !
    !
    ! The k point in crystal coordinates
    xk_(:)  = xk(:,ik)
    CALL cryst_to_cart(1, xk_, at, -1)  
    !
    dist = sum( (xk_(:)-xkq_(:))**2 )
    dist = SQRT(dist)
    !
    IF (dist .lt. 1.d-6) THEN
      !
      ikq=ik
      found(ik) = .true.
      cntr = cntr + 1
      !WRITE(100,'("Match Found for xkq=", 3f8.4, 2x, "ikq =", i5)') xkq_, ikq
      !
    ENDIF
    !
    !WRITE(100,'("ik, xk, xkq, dist", i10, 3(3f12.8,3x), f12.8)') ik, xk_, xkq_, gvect, dist
  ENDDO
  ! In cart coordinates
  CALL cryst_to_cart(1, xkq_, bg, 1)
  CALL cryst_to_cart(1, gvect, bg, 1)
  !
  !WRITE(100,*) (found(ik), ik = 1, nkstot), cntr
  IF ( cntr .gt. 1 ) CALL errore('find_index_1bz','More than 1 match Found! ',cntr)
  IF ( cntr .eq. 0 ) CALL errore('find_index_1bz','No match Found! ',1)
  !
END subroutine


!-----------------------------------------------------------------------
SUBROUTINE find_index_1bz_iterate(xkq, gvect, ikq)
  !-----------------------------------------------------------------------
  !
  ! Find  g such that the k+q = p+g with k and p belonging 
  ! to the original mesh of k points. Take 3x3x3=27 g' vect around
  ! g=0 and loop over p until p=k+q-g'. (k+q from input) 
  ! 
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg
  USE klist,                ONLY : xk, nkstot
  USE lsda_mod,             ONLY : nspin
  !USE mp_world,             ONLY : mpime
  
  !USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xkq(3)
  REAL(DP), INTENT(OUT) :: gvect(3)
  REAL(DP) :: xkq_(3), dist, xk_(3)
  INTEGER :: ikq, ik, cntr
  LOGICAL :: found(nkstot)
  INTEGER :: g1, g2 ,g3, count_g
  !
  !
  count_g = 0
  !
  DO g1=-2,2
     DO g2=-2,2
       DO g3=-2,2
         !
         xkq_(:) = xkq(:) ! the k+q from input
         ! The k+q in crystal coordinates 
         CALL cryst_to_cart(1, xkq_, at, -1)
         found = .false. 
         gvect = (/REAL(g1), REAL(g2), REAL(g3)/)
         !WRITE(100,*) "gvect ", gvect
         !
         xkq_(:) = xkq_(:) - gvect(:)
         !WRITE(100,*) 'xkq_trial', xkq_
         !
         cntr = 0
         DO ik = 1, nkstot/nspin
           ! 
           ! The k point in crystal coordinates
           xk_(:)  = xk(:,ik)
           CALL cryst_to_cart(1, xk_, at, -1)
           !
           dist = sum( (xk_(:)-xkq_(:))**2 )
           dist = SQRT(dist)
           !
           IF (dist .lt. 1.d-6) THEN
             !
             ikq=ik
             found(ik) = .true.
             cntr = cntr + 1
             !WRITE(100,'("Match Found for xkq=", 3f8.4, 2x, "ikq =", i5)') xkq_, ikq
             !
           ENDIF
           !
           !WRITE(100,'("ik, xk, xkq, dist", i10, 3(3f12.8,3x), f12.8)') ik, xk_, xkq_,gvect, dist
         ENDDO
         !
         !WRITE(100,*) (found(ik), ik = 1, nkstot), cntr
         IF ( cntr .gt. 1 ) CALL errore('find_index_1bz','More than 1 match Found!',cntr)
         IF ( cntr .eq. 1 ) THEN 
           !WRITE(100,'("Match Found for xkq=", 3f8.4, 2x, "ikq =", i5, "G =", 3f8.4)') xkq_, ikq, gvect
           GOTO 100
         ELSE  
            count_g = count_g + 1
            !write(100, '(" Next iG =", 3x, 3i)') count_g+1 
         ENDIF
         !
       ENDDO
    ENDDO
  ENDDO
  ! if you arrive here it means none of the G was sucesfull. No match found ,
  ! call error
  CALL errore('find_index_1bz','No match Found! ',1)
  !
100 continue
  ! In cart coordinates
  CALL cryst_to_cart(1, xkq_, bg, 1)
  CALL cryst_to_cart(1, gvect, bg, 1)
  !
  RETURN
  !
END subroutine
