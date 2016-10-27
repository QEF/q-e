!
! Copyright (C) 2014 Federico Zadra
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE wy_pos
USE kinds,  ONLY : DP
IMPLICIT NONE

SAVE
PRIVATE

PUBLIC wypos

CONTAINS
   SUBROUTINE wypos(tau,wp,inp,space_group_number,uniqueb,&
                                      rhombohedral,origin_choice)

   !-----------------------------------------------------------
   ! Convert atomic positions given in Wyckoff convention:
   ! multiplicity-letter + parameter(s), to crystal positions
   !   wp = Wyckoff label (e.g. 8c)
   !   inp(3) = parameter(s) (if needed)
   !-----------------------------------------------------------

      REAL(DP), DIMENSION(3), INTENT(OUT) :: tau
      REAL(DP), INTENT(IN) :: inp(3)
      CHARACTER(LEN=*), INTENT (IN) :: wp

      INTEGER, INTENT(IN) :: space_group_number      
      LOGICAL, INTENT(IN) :: uniqueb, rhombohedral
      INTEGER, INTENT(IN) :: origin_choice
      
      tau=1.d5 

      SELECT CASE (space_group_number)
         CASE (2) !P-1
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ENDIF

         CASE (3) !P2
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ENDIF
               
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (5) !C2
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ENDIF
               
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (6) !Pm
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=inp(2)
               ENDIF
               
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.5_DP
               ENDIF
            ENDIF
         
         CASE (8) !Cm
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF
               
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ENDIF
            ENDIF

         CASE (10) !P2/m
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2k') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2l') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2m') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='2n') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=inp(2)
               ENDIF
              
   
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2i') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2j') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2k') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2l') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2n') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.5_DP
               ENDIF
            ENDIF
         
         CASE (11) !P2(1)/m
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=inp(2)               
               ENDIF
              
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.25_DP
               ENDIF   
            ENDIF

         CASE (12) !C2/m
           IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF
              

            ELSEIF (.NOT.uniqueb) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ENDIF   
            ENDIF

         CASE (13) !P2/c
            IF (uniqueb) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2f') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
              

            ELSEIF (.NOT.uniqueb) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
              
            ENDIF

         CASE (14) !-P2(1)/c
            IF (uniqueb) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF
              

            ELSEIF (.NOT.uniqueb) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ENDIF
              
            ENDIF

         CASE (15) !C2/c
            IF (uniqueb) THEN

              
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25
               ENDIF
              

            ELSEIF (.NOT.uniqueb) THEN

              
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
              
            ENDIF

         CASE (16) !P222
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2o') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2q') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2r') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2s') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2t') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           

         CASE (17) !P222(1)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF

         CASE (18) !P2(1)2(1)2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (20) !C222(1)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF

         CASE (21) !C222
           
           
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF
     
         
         CASE (22) !F222   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
           

         CASE (23) !I222
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           
         CASE (24) !I2(1)2(1)2(1)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF

         CASE (25) !Pmm2
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (26) !Pmc2(1)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (27) !Pcc2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (28) !Pma2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (30) !Pca2(1)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (31) !Pmn2(1)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (32) !Pba2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (34) !Pnn2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF 
            
         CASE (35) !Cmm2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (36) !Cmc2(1)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (37) !Ccc2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF

         CASE (38) !Amm2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (39) !Aem2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF

         CASE (40) !Ama2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (41) !Aea2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (42) !Fmm2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

         CASE (43) !Fdd2
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (44) !Imm2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (45) !Iba2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (46) !Ima2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (47) !Pmmm
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2o') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2q') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2r') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2s') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2t') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4u') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4v') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4w') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4x') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4y') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4z') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

         CASE (48) !Pnnn
            IF (origin_choice==1) THEN 
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.75_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF              
            ENDIF
         
         CASE (49) !Pccm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4p') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           
         CASE (50) !Pban
            IF (origin_choice==1) THEN           
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (51) !Pmma
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (52) !Pnna
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (53) !Pmna
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           

         CASE (54) !Pcca
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           
         CASE (55) !Pbam
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

         CASE (56) !Pccn
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.75_DP
               tau(3)=inp(1)
            ENDIF

         CASE (57) !Pbcm
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
           
         CASE (58) !Pnnm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF


      
         CASE (59) !Pmmn
            IF (origin_choice==1) THEN           
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=inp(2)
               ENDIF
            ENDIF
        
         CASE (60) !Pbcn
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF
           

         CASE (61) !Pbca
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
           

         CASE (62) !Pnma
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF
           

         CASE (63) !Cmcm
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
           
         CASE (64) !Cmce
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (65) !Cmmm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
           

         CASE (66) !Cccm     
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.75_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF 

         CASE (67) !Cmma
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF
           
         CASE (68) !Ccce
            IF (origin_choice==1) THEN   
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (69) !Fmmm
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='16j') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16o') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           
         CASE (70) !Fddd
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
              

            ELSEIF (origin_choice==2) THEN  
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (71) !Immm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           
         CASE (72) !Ibam
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           
         CASE (73) !Ibca
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF
           

         CASE (74) !Imma
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF

         CASE (75) !P4
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (77) !P4(2)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (79) !I4(2)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (80) !I4(1)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (81) !P-4
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           
         CASE (82) !I-4
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           
         CASE (83) !P4/m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

         CASE (84)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           
         CASE (85)
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (86)
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (87) !I4/m
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           

         CASE (88) !I4(1)/a
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN 
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (89) !P422
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ENDIF
           
         CASE (90) !P42(1)2        
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ENDIF
           
         CASE (91) !P4(1)22        
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.375_DP
            ENDIF

         CASE (92) !P4(1)2(1)2        
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ENDIF

         CASE (93) !P4(2)22
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.75_DP
            ENDIF
           
         CASE (94) !P4(2)2(1)2 
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ENDIF

         CASE (95) !P4(3)22 
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.625_DP
            ENDIF

         CASE (96) !P4(2)2(1)2 
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ENDIF

         CASE (97) !I422
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (98) !I4(1)22
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=-inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.125_DP
            ENDIF

         CASE (99) !P4mm
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF           

         CASE (100) !P4bm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF

         CASE (101) !P4(2)cm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (102) !P4(2)nm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (103) !P4cc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (104) !P4nc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (105) !P4(2)mc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF

         CASE (106) !P4(2)bc
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (107) !I4mm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

         CASE (108) !I4cm
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF

         CASE (109) !I4(1)md
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (110) !I4(1)cd
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (111) !P-42m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (112) !P-42c
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           

         CASE (113) !P-42(1)m
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
           
         CASE (114) !P-42(1)c           
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           

         CASE (115) !P-4m2
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF
           

         CASE (116) !P4c2   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           

         CASE (117) !P-4b2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.5_DP
            ENDIF
           

         CASE (118) !P-4n2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
           
         CASE (119) !I-4m2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
           

         CASE (120) !I-4c2
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ENDIF
           
         CASE (121) !I-42m
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (122) !I-42d
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.125_DP
            ENDIF
           
          
         CASE (123) !P4/mmm
           
           
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8r') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8s') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8t') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF
           

         CASE (124) !P4/mmc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP 
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           

         CASE (125) !P/nbm
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

         CASE (126)
            IF (origin_choice==1) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

         CASE (127) !P4/mbm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
           

         CASE (128) !P4/mnc         
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           
         CASE (129) !P4/nmm
            IF (origin_choice==1) THEN   
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF
              

            ELSEIF (origin_choice==2) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
      
         CASE (130) !P4/ncc
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

         CASE (131) !P4(2)/mmc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
           

         CASE (132) !P4(2)mcm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (133) !P4(2)/nbc
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.0_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.00_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

         CASE (134) !P4(2)/nnm
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
              
            ENDIF

         CASE (135) !P3(2)/mbc
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

         CASE (136) !P4(2)/mnm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (137) !P4(2)/nmc
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF     
            ENDIF

         CASE (138) !P4(2)/ncm
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF    
            ENDIF

         CASE (139) !I4/mmm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16m') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (140) !I4/mcm
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
           
         CASE (141) !I4(1)/amd
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.25_DP
                  tau(3)=0.875_DP
               ELSEIF (TRIM(wp)=='16h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
              
            ENDIF

         CASE (142) !I4(1)/acd
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.25_DP
                  tau(3)=0.125_DP
               ENDIF
              
            ENDIF

         CASE(143) !P3
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (146) !R3
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
            ELSE !If HEXAGONAL
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE(147) !P-3
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
           

         CASE (148) !R-3
           IF (rhombohedral) THEN     
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF
            ENDIF

         CASE (149) !P312
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ENDIF
           
         CASE (150) !P321
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF

         CASE (151) !P3(1)12
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=1.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=5.0_DP/6.0_DP
            ENDIF

         CASE (152) !P3(1)21
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=1.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=5.0_DP/6.0_DP
            ENDIF

         CASE (153) !P3(2)12
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=2.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=1.0_DP/6.0_DP
            ENDIF

         CASE (154) !3(2)21
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=2.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=1.0_DP/6.0_DP
            ENDIF

         CASE (155) !R32
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF      
            ENDIF

         CASE (156) !P-3m1
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (157) !P31m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

         CASE (158) !P3c1
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (159) !P31c
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF        

         CASE (160) !R3m
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9b') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF      
            ENDIF

         CASE (161) !R3c
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF      
            ENDIF

         CASE (162) !P-31m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

         CASE (163) !P-31c
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.25_DP
            ENDIF

         CASE (164) !P-3m1
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (165) !P-3c1
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF

         CASE (166) !R-3m
           IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN               
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='18h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

         CASE (167) !R-3c
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6e') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)+0.5_DP
                  tau(3)=0.25_DP
               ENDIF
              
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='18d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

         CASE (168) !P6
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (171) !P6/m
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (172) !P6(4)
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

         CASE (173) !P6(3)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (174) !P-6
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

         CASE (175) !P6/m
           
           
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

         CASE (176) !P6(3)/m
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF

         CASE (177) !P622
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ENDIF
           
         CASE (178) !P6(1)22
            IF (TRIM(wp)=='6a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ENDIF

         CASE (179) !P6(5)22
            IF (TRIM(wp)=='6a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.75_DP
            ENDIF

         CASE (180) !P6(2)22
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ENDIF

         CASE (181) !P6(4)22
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ENDIF

         CASE (182) !P6(3)22
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ENDIF

         CASE (183) !P6mm
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (184) !P6cc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

         CASE (185) !P6(3)cm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

         CASE (186) !P6(3)mc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (187) !P-6m2
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6n') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (188) !P-6c2
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF

         CASE (189) !P-62m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

         CASE (190) !P-62c
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
           

         CASE (191) !P6/mmm
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12o') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
           
         CASE (192) !P6/mcc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP       
            ENDIF
           
         CASE (193) !P6(3)/mcm
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF


         CASE (194) !P6(3)mmc
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (195) !P23
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ENDIF

         CASE (196) !F23
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (197) !I23
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ENDIF

         CASE (198) !P2(1)3
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF

         CASE (199) !I2(1)3
            IF (TRIM(wp)=='8a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF

         CASE (200) !Pm-3
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (201) !Pn-3
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

         CASE (202) !Fm-3
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (203) !Fd-3
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ENDIF
            ENDIF

         CASE (204) ! Im-3
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (205) !Pa-3
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF
           
         CASE (206) !Ia-3
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (207) !P432
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF
           

         CASE (208) !P4(2)32
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='12l') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.5_DP
            ENDIF
           
         CASE (209) !F432
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (210) !F4(1)32
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='16d') THEN
               tau(1)=0.625_DP
               tau(2)=0.625_DP
               tau(3)=0.625_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
           
         CASE (211) !I432
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ENDIF
           
         CASE (212) !P4(3)32
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.625_DP
               tau(2)=0.625_DP
               tau(3)=0.625_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
           
         CASE (213) !P4(1)32
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.375_DP
               tau(2)=0.375_DP
               tau(3)=0.375_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.875_DP
               tau(2)=0.875_DP
               tau(3)=0.875_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.25_DP
            ENDIF

         CASE (214) !I4(I)32
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.875_DP
               tau(2)=0.875_DP
               tau(3)=0.875_DP
            ELSEIF (TRIM(wp)=='12c') THEN
               tau(1)=0.125_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.625_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.25_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
           

         CASE (215) !P-43m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (216) !F-43m
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (217) !I-43m
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (218) !P-43n
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
           
         CASE (219) !F-43c
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (220) !I-43d
            IF (TRIM(wp)=='12a') THEN
               tau(1)=0.375_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12b') THEN
               tau(1)=0.875_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF
           
         CASE (221) !Pm-3m
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='24l') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='24m') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           

         CASE (222) !Pn-3n
            IF (origin_choice==1) THEN             
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='12e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='24g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='24g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
            ENDIF

         CASE (223) !Pm-3n
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24j') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='24k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (224) !Pn-3m
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='24i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24j') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='24i') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ELSEIF (TRIM(wp)=='24k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

         CASE (225) !Fm-3m
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='96k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (226) !Fm-3c
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='64g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96h') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
           
         CASE (227) !Fd-3m
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='96h') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.25_DP
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.375_DP
                  tau(2)=0.375_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='96h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
            ENDIF

         CASE (228) !Fd-3c
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='16a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='32b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='32c') THEN
                  tau(1)=0.375_DP
                  tau(2)=0.375_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='48d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='64e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='96f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.25_DP
               ENDIF
              
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='16a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='32b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='32c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='48d') THEN
                  tau(1)=0.875_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='64e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='96f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
            ENDIF

         CASE (229) !Im-3m
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='48j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='48k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

         CASE (230) !Ia-3d
            IF (TRIM(wp)=='16a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16b') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.125_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.375_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF

          CASE DEFAULT
            CALL errore('wypos','group not recognized',1)
          END SELECT

         IF (tau(1)==1.d5.OR.tau(2)==1.d5.OR.tau(3)==1.d5) THEN
            IF (inp(1)==1.d5.OR.inp(2)==1.d5.OR.inp(3)==1.d5) THEN
               CALL errore('wypos','wyckoff position not found',1)
            ELSE
               CALL infomsg('wypos','wyckoff position not found, assuming x y z')
               tau(:)=inp(:)
            END IF
         END IF

      END SUBROUTINE wypos
END MODULE
