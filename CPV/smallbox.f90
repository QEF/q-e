!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE small_box
!------------------------------------------------------------------------------!

      !  This module contains the basis vector of the small sub-cell 
      !  (small box) used for charge augmentation process

      USE kinds, ONLY : DP
!
      IMPLICIT NONE
      SAVE
      !  atb: lattice crystal vectors in "alatb" units (equivalent to "at")
      REAL(DP) :: atb(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
      !  bgb: reciprocal lattice vectors, in 2pi/alatb units (equiv to "bg")
      REAL(DP) :: bgb(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
      !  omegab: volume of the small boxes 
      REAl(DP) :: omegab = 0.0_DP
      !  alatb: lattice parameter of the boxes (the equivalent of "alat")
      REAL(DP) :: alatb  = 0.0_DP
      !  tpibab: 2pi/alatb
      REAL(DP) :: tpibab = 0.0_DP

!------------------------------------------------------------------------------!
   CONTAINS
!------------------------------------------------------------------------------!
!

     SUBROUTINE small_box_set( alat, omega, at, rat1, rat2, rat3, tprint )
       USE constants, ONLY: pi
       USE io_global, ONLY: stdout, ionode
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: alat, omega, at(3,3), rat1, rat2, rat3
       LOGICAL,  INTENT(IN) :: tprint

       alatb  = alat * rat1
       IF( alatb <= 0.0_DP ) CALL errore(' small_box_set ', ' alatb <= 0 ', 1 )
       tpibab = 2.0_DP * pi / alatb
       atb(:,1) = at(:,1)*alat * rat1 / alatb
       atb(:,2) = at(:,2)*alat * rat2 / alatb
       atb(:,3) = at(:,3)*alat * rat3 / alatb 
       omegab = omega * rat1 * rat2 * rat3
!
       CALL recips( atb(1,1), atb(1,2), atb(1,3), bgb(1,1), bgb(1,2), bgb(1,3) )

       IF( tprint .AND. ionode ) THEN
          WRITE( stdout,*)
          WRITE( stdout,220)
220       format( 3X, 'unit vectors of box grid cell',/,                        &
     &            3X, 'in real space:',25x,'in reciprocal space:')
          WRITE( stdout,'(3X,3f10.4,10x,3f10.4)') atb(:,1)*alatb, bgb(:,1)
          WRITE( stdout,'(3X,3f10.4,10x,3f10.4)') atb(:,2)*alatb, bgb(:,2)
          WRITE( stdout,'(3X,3f10.4,10x,3f10.4)') atb(:,3)*alatb, bgb(:,3)
       END IF

       RETURN
     END SUBROUTINE small_box_set
!
!------------------------------------------------------------------------------!
   END MODULE small_box
!------------------------------------------------------------------------------!
