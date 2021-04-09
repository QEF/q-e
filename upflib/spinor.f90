!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
FUNCTION spinor( l, j, m, spin )
  !------------------------------------------------------------------------
  !! This function calculates the numerical coefficient of a spinor
  !! with orbital angular momentum l, total angular momentum j, 
  !! projection along z of the total angular momentum m+-1/2. Spin selects
  !! the up (spin=1) or down (spin=2) coefficient.
  !
  USE upf_kinds, only: dp
  !
  IMPLICIT NONE
  !
  REAL(DP) :: spinor
  !! the spinor coefficient
  INTEGER :: l
  !! orbital angular momentum
  INTEGER :: m
  !! projection of the total angular momentum+-1/2
  INTEGER :: spin
  !! 1 or 2 select the component
  REAL(DP) :: j
  !! total angular momentum
  !
  ! ... local variables
  !
  REAL(DP) :: denom     ! denominator
  !
  IF ( spin/=1 .AND. spin/=2 ) CALL upf_error( 'spinor', 'spin direction unknown', 1 )
  IF ( m<-l-1  .OR.  m>l )     CALL upf_error( 'spinor', 'm not allowed', 1 )
  !
  denom = 1.d0 / (2.d0*l+1.d0)
  !
  IF ( ABS(j-l-0.5d0) < 1.d-8 ) THEN
     IF (spin==1) spinor = SQRT((l+m+1.d0)*denom)
     IF (spin==2) spinor = SQRT((l-m)*denom)
  ELSEIF (ABS(j-l+0.5d0) < 1.d-8) THEN
     IF (m < -l+1) THEN
        spinor=0.d0
     ELSE
        IF (spin == 1) spinor = SQRT((l-m+1.d0)*denom)
        IF (spin == 2) spinor = -SQRT((l+m)*denom)
     ENDIF
  ELSE
     CALL upf_error( 'spinor', 'j and l not compatible', 1 )
  ENDIF
  !
  RETURN
  !
END FUNCTION spinor
