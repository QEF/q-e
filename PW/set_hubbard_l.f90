!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
integer function set_hubbard_l(psd) result (hubbard_l)
!
implicit none
character(len=2) :: psd
!
! TRANSITION METALS
!
if (psd.eq.'V'  .or. psd.eq.'Cr' .or. psd .eq.'Mn' .or. psd.eq.'Fe' .or. &
    psd.eq.'Co' .or. psd.eq.'Ni' .or. psd .eq.'Cu' ) then  
    hubbard_l = 2  
!
! RARE EARTHS
!
elseif (psd .eq.'Ce') then  
   hubbard_l =  3
!
! OTHER ELEMENTS
!
elseif (psd .eq.'C' .or. psd .eq. 'O') then  
   hubbard_l =  1
elseif (psd .eq.'H') then  
   hubbard_l =  0
else  
   hubbard_l = -1
   call errore ('set_hubbard_l','pseudopotential not yet inserted', 1)
endif  
return  

end function set_Hubbard_l
