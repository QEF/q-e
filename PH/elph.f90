!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module el_phon
use parameters, only : DP
complex(kind=DP), pointer ::  el_ph_mat(:,:,:,:) !nbnd,nbnd,nks,3*nat)
end module el_phon
