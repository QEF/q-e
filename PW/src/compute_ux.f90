!
! Copyright (C) 2007-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
   SUBROUTINE compute_ux(m_loc,ux,nat)
!
!   This subroutine determines the direction of a fixed quantization axis
!   from the starting magnetization.
!
   USE kinds, ONLY : dp
   USE constants, ONLY: pi, eps12
   USE io_global,  ONLY :  stdout
   USE noncollin_module, ONLY : lsign

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nat    ! number of atoms
   REAL(DP), INTENT(OUT) :: ux(3)     ! fixed direction to calculate signs
   REAL(DP), INTENT(IN) :: m_loc(3,nat)     ! local moments
         
   REAL(DP) :: amag, uxmod ! modulus of the magnetization and of ux 

   INTEGER :: na          ! counter on atoms
   INTEGER :: starting_na ! auxiliary variable
   LOGICAL :: is_parallel ! external function true if two vectors are parallel
!
!  Do not use the sign feature in the general case
!
   lsign=.FALSE.
   ux=0.0_DP

   starting_na=0
   DO na=1,nat
      amag=m_loc(1,na)**2+m_loc(2,na)**2+m_loc(3,na)**2
      IF (amag > eps12) THEN
         ux(:)=m_loc(:,na)
         starting_na=na
         lsign=.TRUE.
         GOTO 20
      ENDIF
   ENDDO
20 CONTINUE
!
!  The sign feature is used only when all initial magnetizations are parallel
!  to a fixed direction that is taken as the quantization axis.
!
   DO na=starting_na+1, nat
      lsign=lsign.AND.is_parallel(ux,m_loc(:,na))
   ENDDO

   IF (lsign) THEN
      uxmod=ux(1)**2+ux(2)**2+ux(3)**2
      IF (uxmod<eps12) CALL errore('compute_ux','strange uxmod',1) 
      ux=ux/SQRT(uxmod)
      WRITE( stdout,'(/,5x,"Fixed quantization axis for GGA: ", 3f12.6)') &
                          ux(1), ux(2), ux(3)
   ENDIF
   RETURN
   END SUBROUTINE compute_ux
