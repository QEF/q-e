!
! Copyright (C) 2001 - 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  !-----------------------------------------------------------------------
subroutine set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the possible vectors G associated
  ! to the symmetries of the small group of q: Sq -> q + G
  ! Furthermore if minus_q and irotmq are set if finds the G for Sq -> -q+G.
  !
  !  The dummy variables
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : bg, at
  USE control_lr, ONLY : lgamma
  implicit none

  REAL(DP), PARAMETER :: accep=1.e-5_dp

  real(DP), INTENT(IN) :: xq (3)
  ! input: the q point 
  real(DP), INTENT(OUT) ::gi (3, 48), gimq (3)
  ! output: the G associated to a symmetry:[S(irotq)*q - q]
  ! output: the G associated to:  [S(irotmq)*q + q]

  LOGICAL, INTENT(IN) :: minus_q
  ! input: .t. if there is sym.ops. such that Sq=-q+G 
  INTEGER, INTENT(IN) :: s (3, 3, 48), nsymq, nsym
  ! input: the symmetry matrices
  ! input: dimension of the small group of q

  INTEGER, INTENT(OUT) :: irotmq
  ! input: op. symmetry: s_irotmq(q)=-q+G

  real(DP) :: wrk (3), aq (3), raq (3), zero (3)
  ! additional space to compute gi and gimq
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: isym, ipol, jpol
  ! counter on symmetry operations
  ! counter on polarizations
  ! counter on polarizations

  logical :: eqvect
  ! logical function, check if two vectors are equal
  !
  !  Set to zero some variables and transform xq to the crystal basis
  !
  zero = 0.d0
  gi   = 0.d0
  gimq = 0.d0
  irotmq = 0
  IF (lgamma) THEN
     irotmq=1
     RETURN
  ENDIF
  aq = xq
  call cryst_to_cart (1, aq, at, - 1)
  !
  !   test all symmetries to see if the operation S sends q in q+G ...
  !
  do isym = 1, nsymq
     raq = 0.d0
     do ipol = 1, 3
        do jpol = 1, 3
           raq (ipol) = raq (ipol) + DBLE (s (ipol, jpol, isym) ) * &
                aq (jpol)
        enddo
     enddo
     if (.NOT. eqvect (raq, aq, zero, accep) ) CALL errore('set_giq',&
                            'problems with the input group',1)
     do ipol = 1, 3
        wrk (ipol) = raq (ipol) - aq (ipol)
     enddo
     call cryst_to_cart (1, wrk, bg, 1)
     gi (:, isym) = wrk (:)
     IF (irotmq == 0) THEN
        raq=-raq
        IF (eqvect (raq, aq, zero, accep)) THEN
           irotmq=isym
           wrk = aq - raq 
           call cryst_to_cart (1, wrk, bg, 1)
           gimq = wrk 
        ENDIF
     ENDIF
  enddo
  !
  !   ... and in -q+G
  !
  if (minus_q.and.irotmq==0) then
     do isym = nsymq+1,nsym
        raq = 0.d0
        do ipol = 1, 3
           do jpol = 1, 3
              raq (ipol) = raq (ipol) + DBLE (s (ipol, jpol, isym) ) * &
                   aq (jpol)
           enddo
        enddo
        raq=-raq
        if (eqvect (raq, aq, zero, accep) ) then
           wrk = aq - raq 
           call cryst_to_cart (1, wrk, bg, 1)
           gimq (:) = wrk (:)
           irotmq=isym
        endif
        if (irotmq /= 0 ) exit
     enddo
  endif
  IF (minus_q.AND. irotmq == 0 ) &
     CALL errore('set_giq','problem with minus_q',1)
  !
  return
end subroutine set_giq

