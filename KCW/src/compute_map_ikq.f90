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
SUBROUTINE compute_map_ikq ()
  !-----------------------------------------------------------------------
  !
  !! This routine compute the map to go from (iq,ik) --> ikq_1BZ
  !! In general k+q lies outside the 1BZ, but if k and q belongs to the 1BZ 
  !! there exist a G vector such that k + q = k'+ G with k' in the 1BZ
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, nkstot
  USE io_global,            ONLY : ionode, ionode_id
  USE klist,                ONLY : nkstot
  USE control_kcw,          ONLY : map_ikq, shift_1bz
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_image_comm
  !
#ifdef DEBUG
  USE mp_world,             ONLY : mpime
#endif
  !
  IMPLICIT NONE
  ! 
  INTEGER :: iq, ik, ikq
  ! Counters for the k/q points in the BZ
  !
  REAL(DP) :: xkq(3), gvect(3)
  ! the k+q coordinate and the G vector that shift it into the 1 BZ
  !
  IF ( ALLOCATED (map_ikq) ) DEALLOCATE (map_ikq)
  IF ( ALLOCATED (shift_1bz) ) DEALLOCATE (shift_1bz)
  !
  ALLOCATE ( map_ikq(nkstot, nkstot) )
  ALLOCATE ( shift_1bz(3, nkstot, nkstot) )
  map_ikq (:,:) = 0
  shift_1bz(:,:,:) = 0.D0 
  !
  ! CHECK: does xk (for the ionode) have all the kpoint coordinates when npool>1 ? 
  IF (ionode) THEN 
    !
    DO iq = 1, nkstot
      !
      DO ik = 1, nkstot
        !
        ! the k+q coordinate (cart)
        xkq(:) = xk(:,ik)+xk(:,iq)
#ifdef DEBUG
        WRITE(mpime+100,'("iq, xq", i5, 3f8.4)') iq, xk(:,iq)
        WRITE(mpime+100,'("ik, xk", i5, 3f8.4)') ik, xk(:,ik)
#endif
        ! the index of the corresponding k point in the IBZ
        gvect(:) = 0.D0
        CALL find_index_1bz(xkq, gvect, ikq)
        map_ikq(iq,ik) = ikq
        shift_1bz(:,iq,ik) = gvect(:)
#ifdef DEBUG
        WRITE(mpime+100,'("The map (iq,ik) --> ikq", 2i, 8x, i)') iq,ik, ikq
        WRITE(mpime+100,'("xkq, shift" , 5x, 3f8.4, 5x, 3f8.4)') xkq(:), gvect(:)
        WRITE(mpime+100,'("ikq, xkq_1BZ", i5, 3f8.4,/)') ikq, xkq(:)-gvect(:)
#endif
        !
      ENDDO
      !
    ENDDO
    !
  ENDIF 
  CALL mp_bcast ( map_ikq, ionode_id, intra_image_comm )
  !
END subroutine



!-----------------------------------------------------------------------
SUBROUTINE find_index_1bz(xkq, gvect, ikq)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg
  USE klist,                ONLY : xk, nkstot
  !USE mp_world,             ONLY : mpime
  
  !USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xkq(3)
  REAL(DP), INTENT(OUT) :: gvect(3)
  REAL(DP) :: xkq_(3), dist
  INTEGER :: ikq, ik, cntr
  LOGICAL :: found(nkstot)
  !
  !
  found(:) = .false.
  xkq_(:)=xkq(:)
  ! The k+q in crystal coordinates 
  !
  CALL cryst_to_cart(1, xkq_, at, -1)
  ! The integer part of the k+q (says if is inside the 1BZ (=0) or not (/=0))
  !
  gvect(:) = REAL( INT( xkq_(:) ) )
  ! The k+q shifted inside the 1BZ
  !
  xkq_(:) = xkq_(:) - int(xkq_(:) )
  ! In cart coordinates
  !
  CALL cryst_to_cart(1, xkq_, bg, 1)
  CALL cryst_to_cart(1, gvect, bg, 1)
  !
  cntr = 0
  DO ik = 1, nkstot
    !
    dist = sum( (xk(:,ik)-xkq_(:))**2 )
    dist = SQRT(dist)
    IF (dist .lt. 1.d-8) THEN
      !
      ikq=ik
      found(ik) = .true.
      cntr = cntr + 1
      !WRITE(mpime+100,'("Match Found for xkq=", 3f8.4, 2x, "ikq =", i5)') xkq_, ikq
      !
    ENDIF
    !
  ENDDO
  !
  !WRITE(*,*) (found(ik), ik = 1, nkstot), cntr
  IF ( cntr .gt. 1 ) CALL errore('find_indez_1BZ','More than 1 match Found! ',1)
  IF ( cntr == 0 )   CALL errore('find_indez_1BZ','No match Found! ',1)
  !
END subroutine

