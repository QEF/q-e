!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine init_ns  
   !-----------------------------------------------------------------------
   !
   ! This routine computes the starting ns (for lda+U calculation) filling
   ! up the d states (the only interested by the on-site potential for the
   ! moment) according to the Hund's rule (valid for the isolated atoms on
   ! which starting potential is built), and to the starting_magnetization:
   ! majority spin levels are populated first, then the remaining electrons
   ! are equally distributed among the minority spin states
   !
#include "machine.h"

   use pwcom  
   implicit none

   real(kind=DP) :: totoc  

   integer :: na, nt, is, m1, majs, mins, i  
   logical :: nm        ! true if the atom is non magnetic

   ns(:,:,:,:) = 0.d0

   do na = 1, nat  
      nt = ityp (na)  
      if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
         call tabd (nt, totoc)  
         nm=.false.
         if (starting_magnetization (nt) .gt.0.d0) then  
            majs = 1  
            mins = 2  
         elseif (starting_magnetization (nt) .lt.0.d0) then  
            majs = 2  
            mins = 1  
         else  
            nm = .true.  
         endif  
         if (.not.nm) then  
            if (totoc.gt.5.d0) then  
               do m1 = 1, 5  
                  ns (na, majs, m1, m1) = 1.d0  
                  ns (na, mins, m1, m1) = (totoc - 5.d0) / 5.d0  
               enddo  
            else  
               do m1 = 1, 5  
                  ns (na, majs, m1, m1) = totoc / 5.d0  
               enddo  
            endif  
         else  
            do m1 = 1, 5  
               ns (na, 1, m1, m1) = totoc / 10.d0  
               ns (na, 2, m1, m1) = ns (na, 1, m1, m1)  
            enddo  
         endif  
      endif  
   enddo  
   return  
end subroutine init_ns

