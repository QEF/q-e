!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine init_run  
  !-----------------------------------------------------------------------
  !
  use pwcom  
  implicit none

  call start_clock ('init_run')  
  write (6, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")')  

  gamma_only =.false.
  write (6, '(5x,"Complex Hamiltonian")')  
  write (6, 9010) ntypx, npsx, lmaxx, npk  
  write (6, 9200) nbrx, lqmax, nqfm  

  call iosys
  call setup  
  call summary  
  call allocate_nlpot  
  call allocate_locpot  
  call allocate_wfc  

  call openfil  

  call hinit0  
  call potinit  

  call newd  

  call wfcinit  
  call stop_clock ('init_run')  
  write (6, * )  
  call show_memory ()  

  return  
9010 format (/5x,'current dimensions of program pwscf are:'/5x, &
       &'ntypx=',i5,' npsx =',i5,' lmax =',i5,' npk =',i5)

9200 format( 5x,'nbrx =',i6,' lqmax =',i6,' nqfm =',i6 )  
end subroutine init_run

