!
! Copyright (C) 2004 Vanderbilt's group at Rutgers University, NJ
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE kp_strings ( nppstr, gdir, nrot, s, bg, npk, &
                        k1,k2,k3, nk1,nk2,nk3, nks, xk, wk )


!  --- Usage of modules ---
   USE kinds, ONLY: dp

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
!  time reversal and no magnetic symmetries assumed
   INTEGER :: t_rev(48) = 0
   LOGICAL :: time_reversal = .true., skip_equivalence=.FALSE.
   REAL(dp) :: dk(3)
   REAL(dp) :: xk0(3,npk)
   REAL(dp) :: wk0(npk)

!  --- Generate a k-point grid in the two dimensions other than gdir ---
   IF (gdir == 1) THEN
      CALL kpoint_grid (nrot, time_reversal, skip_equivalence, s, t_rev, bg, &
                        npk, k1,k2,k3, 1,nk2,nk3, nks, xk0, wk0 ) 

   ELSE IF (gdir == 2) THEN
      CALL kpoint_grid (nrot, time_reversal, skip_equivalence, s, t_rev, bg, &
                        npk, k1,k2,k3, nk1,1,nk3, nks, xk0, wk0 ) 
   ELSE IF (gdir == 3) THEN
      CALL kpoint_grid (nrot, time_reversal, skip_equivalence, s, t_rev, bg, &
                        npk, k1,k2,k3, nk1,nk2,1, nks, xk0, wk0 ) 
   ELSE
      CALL errore('kp_strings','gdir different from 1, 2, or 3',1)
   END IF

!  --- Generate a string of k-points for every k-point in the 2D grid ---
   kindex=0
   dk(1)=bg(1,gdir)/REAL(nppstr-1,dp)
   dk(2)=bg(2,gdir)/REAL(nppstr-1,dp)
   dk(3)=bg(3,gdir)/REAL(nppstr-1,dp)
   DO iort=1,nks
      DO ipar=1,nppstr
         kindex=kindex+1
         xk(1,kindex)=xk0(1,iort)+REAL(ipar-1,dp)*dk(1)
         xk(2,kindex)=xk0(2,iort)+REAL(ipar-1,dp)*dk(2)
         xk(3,kindex)=xk0(3,iort)+REAL(ipar-1,dp)*dk(3)
         wk(kindex)=wk0(iort)/REAL(nppstr,dp)
      END DO
   END DO
   nks=nks*nppstr

END SUBROUTINE kp_strings
