MODULE wy_pos
USE kinds,  ONLY : DP
IMPLICIT NONE

SAVE
PRIVATE

PUBLIC wypos

CONTAINS
   SUBROUTINE wypos(tau,wp,space_group_number,uniqueb,rhombhoedral,origin_choice)

   !-----------------------------------------------------------
   !Convert the Wycoff positions(multiplicity-letter) in crystall
   !positions, only for well defined posistion: where the
   !multiplicity and the Wyckoff letter assigne 3 coordinates.
   !-----------------------------------------------------------

      REAL(DP), DIMENSION(3), INTENT(INOUT) :: tau
      CHARACTER(LEN=*), INTENT (IN) :: wp

      INTEGER, INTENT(IN) :: space_group_number      
      LOGICAL, INTENT(IN) :: uniqueb, rhombhoedral
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
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
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
               ENDIF
              

            ELSEIF (origin_choice==2) THEN

              
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
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
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
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
            ENDIF
           

         CASE (73) !Ibca
           
           
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
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
            ENDIF
           

         CASE (82) !I-4
           
           
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
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
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
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
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
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
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
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
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
               ENDIF
              

            ELSEIF (origin_choice==2) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ENDIF
           

         CASE (90) !P42(1)2
           
           
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
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
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
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
            ENDIF
           

         CASE (111) !P-42m
           
           
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
            ENDIF
           

         CASE (115)
           
           
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
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
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
               tau(3)=0.75_DP
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
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
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
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
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
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
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
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
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
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
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
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
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
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
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
               ENDIF
              

            ELSEIF (origin_choice==2) THEN
 
               
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
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
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
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
           IF (rhombhoedral) THEN

              
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ENDIF
              

            ELSEIF (.NOT.rhombhoedral) THEN
 
               
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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
               tau(1)=0.3333333333333333333333_DP
               tau(2)=0.6666666666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.3333333333333333333333_DP
               tau(2)=0.6666666666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.6666666666666666666667_DP
               tau(2)=0.3333333333333333333333_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.6666666666666666666667_DP
               tau(2)=0.3333333333333333333333_DP
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
            ENDIF
           

         CASE (155) !R32
            IF (rhombhoedral) THEN

              
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ENDIF
              

            ELSEIF (.NOT.rhombhoedral) THEN
 
               
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
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
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.66666666666666667_DP
               tau(2)=0.33333333333333333_DP
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
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
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
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ENDIF
           

         CASE (166) !R-3m
           IF (rhombhoedral) THEN

              
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ENDIF
              

            ELSEIF (.NOT.rhombhoedral) THEN
 
               
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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

         CASE (167) !R-3c
            IF (rhombhoedral) THEN

              
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF
              

            ELSEIF (.NOT.rhombhoedral) THEN
 
               
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF
              
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
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.66666666666666667_DP
               tau(2)=0.33333333333333333_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.66666666666666667_DP
               tau(2)=0.33333333333333333_DP
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
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
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
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.66666666666666667_DP
               tau(2)=0.33333333333333333_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
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
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.33333333333333333_DP
               tau(2)=0.66666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.75_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.6666666666666667_DP
               tau(2)=0.3333333333333333_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.6666666666666667_DP
               tau(2)=0.3333333333333333_DP
               tau(3)=0.5_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.6666666666666667_DP
               tau(2)=0.3333333333333333_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.6666666666666667_DP
               tau(2)=0.3333333333333333_DP
               tau(3)=0.25_DP
            ENDIF
           

         CASE (189) !P-62m
           
           
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.6666666666666667_DP
               tau(2)=0.3333333333333333_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
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
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.3333333333333333_DP
               tau(2)=0.6666666666666667_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
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
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
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
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
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
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
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
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
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
            ENDIF
          CASE DEFAULT
             CALL errore('wypos','group not recognized',1)          
         END SELECT

         IF (tau(1)==1.d5.OR.tau(2)==1.d5.OR.tau(3)==1.d5) &
            CALL errore('wypos','position not available',1)

         RETURN
      END SUBROUTINE wypos
END MODULE wy_pos
