! dftd3 program for computing the dispersion energy and forces from cart
! and atomic numbers as described in
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
!
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456
! (for BJ-damping)
!
! Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germany
!
! Repackaging of the original code without any change in the functionality:
!
! Copyright (C) 2016, Bálint Aradi
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

module dftd3_sizes
  implicit none

  integer, parameter :: max_elem = 94

  ! maximum coordination number references per element
  integer, parameter :: maxc = 5

  integer, parameter :: nlines = 32385

  integer, parameter :: npars = 161925

end module dftd3_sizes
