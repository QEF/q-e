!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module funct
character (len=20) :: dft
                         ! exch-corr name (in human-readable format)

integer :: iexch, icorr, igcx, igcc
                         ! type of exchange
                         ! type of correlation
                         ! type of gradient correction on exchange
                         ! type of gradient correction on correlations
end module funct
