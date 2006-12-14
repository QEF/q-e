!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Tue Nov 30 10:59:55 MET 1999
!  ----------------------------------------------
!  BEGIN manual

       SUBROUTINE entropy(f,temp,nx,ent)

!  this routine computes the entropic contribution due to the finite
!  temperature assigned to electrons when computing occupation numbers
!  ----------------------------------------------
!  END manual

       USE kinds
       IMPLICIT NONE

! ...  declare subroutine arguments
       INTEGER nx
       REAL(DP)  f(nx),temp,ent

! ...  declare other variables
       INTEGER i
       REAL(DP) fm
       REAL(DP), PARAMETER :: eps = 1.0d-10

!  end of declarations
!  ----------------------------------------------

       ent=0.d0
       DO i=1,nx
         fm=0.5d0*f(i)
         ent = ent+ fm*log(eps+fm)+(1.d0-fm)*log(eps+1.d0-fm)
       END DO
       ent=-2.d0*temp*ent

       RETURN
       END SUBROUTINE entropy

       subroutine entropy_s(f,temp,nx,ent)
       use kinds
       implicit none
       integer nx
       integer i
       real(DP)  f(nx),temp,ent, fm,eps
       parameter(eps=1.d-10)

       ent=0.d0
       do i=1,nx
         fm=0.5d0*f(i)
         ent = ent+ fm*log(eps+fm)+(1.d0-fm)*log(eps+1.d0-fm)
       enddo
       ent=-2.d0*temp*ent

       return
       end subroutine entropy_s

