SUBROUTINE kp_strings &
           ( nppstr,gdir,nrot,s,bg,npk,k1,k2,k3,nk1,nk2,nk3,nks,xk,wk)

#include "machine.h"

!  --- Usage of modules ---
   USE parameters, ONLY: dp

!  --- No implicit definitions ---
   IMPLICIT NONE

!  --- Input arguments ---
   INTEGER , INTENT(IN) :: k1
   INTEGER , INTENT(IN) :: k2
   INTEGER , INTENT(IN) :: k3
   INTEGER , INTENT(IN) :: nk1
   INTEGER , INTENT(IN) :: nk2
   INTEGER , INTENT(IN) :: nk3
   INTEGER , INTENT(IN) :: nppstr
   INTEGER , INTENT(IN) :: npk
   INTEGER , INTENT(IN) :: nrot
   INTEGER , INTENT(IN) :: gdir
   INTEGER , INTENT(IN) :: s(3,3,48)
   REAL(dp) , INTENT(IN) :: bg(3,3)

!  --- Output arguments ---
   INTEGER , INTENT(OUT) :: nks
   REAL(dp), INTENT(OUT) :: xk(3,npk)
   REAL(dp), INTENT(OUT) :: wk(npk)

!  --- Internal definitions ---
   INTEGER :: ipar
   INTEGER :: iort
   INTEGER :: kindex
   INTEGER :: nkpol
   REAL(dp) :: dk
   REAL(dp) :: xk0(3,npk)
   REAL(dp) :: wk0(npk)

!  --- Generate a k-point grid in the two dimensions other than gdir ---
   IF (gdir == 1) THEN
      CALL kpoint_grid(nrot,s,bg,npk,k1,k2,k3,1,nk2,nk3,nks,xk0,wk0) 
   ELSE IF (gdir == 2) THEN
      CALL kpoint_grid(nrot,s,bg,npk,k1,k2,k3,nk1,1,nk3,nks,xk0,wk0) 
   ELSE IF (gdir == 3) THEN
      CALL kpoint_grid(nrot,s,bg,npk,k1,k2,k3,nk1,nk2,1,nks,xk0,wk0) 
   ELSE
      CALL errore('kp_strings','gdir different from 1, 2, or 3',1)
   END IF

!  --- Generate a string of k-points for every k-point in the 2D grid ---
   kindex=0
   dk=1.0_dp/REAL(nppstr-1,dp)
   DO iort=1,nks
      DO ipar=1,nppstr
         kindex=kindex+1
         xk(1,kindex)=xk0(1,iort)
         xk(2,kindex)=xk0(2,iort)
         xk(3,kindex)=xk0(3,iort)
         xk(gdir,kindex)=REAL(ipar-1,dp)*dk
         wk(kindex)=wk0(iort)/REAL(nppstr,dp)
      END DO
   END DO
   nks=nks*nppstr

END SUBROUTINE kp_strings
