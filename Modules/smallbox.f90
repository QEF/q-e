!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE small_box
!------------------------------------------------------------------------------!

      !  This module contains the basis vector of the small sub-cell (small box)
      !  used for charge augmentation process

      USE kinds, ONLY : dbl
!
      IMPLICIT NONE
      SAVE

        !  a1, a2 and a3 are the simulation cell base vector as calculated from celldm

      REAL(dbl) :: a1b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
      REAL(dbl) :: a2b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
      REAL(dbl) :: a3b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)

      REAL(dbl) :: b1b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
      REAL(dbl) :: b2b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
      REAL(dbl) :: b3b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)

      REAL(dbl) :: ainvb(3,3) = 0.0d0

      REAl(dbl) :: omegab = 0.0d0  !  volume of the small boxes 

      REAL(dbl) :: tpibab = 0.0d0

      REAL(dbl) :: alatb  = 0.0d0

!------------------------------------------------------------------------------!
   CONTAINS
!------------------------------------------------------------------------------!
!

     SUBROUTINE small_box_set( alat, omega, a1, a2, a3, rat1, rat2, rat3 )
       USE constants, ONLY: pi
       USE io_global, ONLY: stdout
       IMPLICIT NONE
       REAL(dbl), INTENT(IN) :: alat, omega, a1(3), a2(3), a3(3), rat1, rat2, rat3
       INTEGER :: i

       alatb  = alat * rat1
       IF( alatb <= 0.0d0 ) CALL errore(' small_box_set ', ' alatb <= 0 ', 1 )
       tpibab = 2.d0 * pi / alatb
       a1b = a1 * rat1
       a2b = a2 * rat2
       a3b = a3 * rat3
       omegab = omega * rat1 * rat2 * rat3
!
       CALL recips( a1b, a2b, a3b, b1b, b2b, b3b )
       b1b = b1b * alatb
       b2b = b2b * alatb
       b3b = b3b * alatb

       WRITE( stdout,*)
       WRITE( stdout,220)
220    format(' unit vectors of box grid cell',/,                        &
     &       ' in real space:',25x,'in reciprocal space:')
       WRITE( stdout,'(3f10.4,10x,3f10.4)') a1b, b1b
       WRITE( stdout,'(3f10.4,10x,3f10.4)') a2b, b2b
       WRITE( stdout,'(3f10.4,10x,3f10.4)') a3b, b3b

       ainvb(1,:) = b1b(:) / alatb
       ainvb(2,:) = b2b(:) / alatb
       ainvb(3,:) = b3b(:) / alatb

       RETURN
     END SUBROUTINE
!
!------------------------------------------------------------------------------!
   END MODULE small_box
!------------------------------------------------------------------------------!
