  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  program kgen
  !----------------------------------------------------------------------
  ! generate a uniform mesh of kpoints in the first BZ [0,1)
  !
  !----------------------------------------------------------------------
  use kinds,         only : DP
  implicit none
  integer, parameter :: n = 100
  integer :: i, j, k
  real(kind=DP) :: x,y,z,w
  !
  WRITE(*,*) n*n*n
  w = 2.d0/float(n*n*n)
  DO k=1,n
  DO i=1,n
  DO j=1,n
     x = (float(i)-1.d0)/float(n)
     y = (float(j)-1.d0)/float(n)
     z = (float(k)-1.d0)/float(n)
     WRITE(*,'4f20.15') x,y,z,w
  ENDDO
  ENDDO
  ENDDO
  !
  stop
  !
  !----------------------------------------------------------------------
  end program kgen
  !----------------------------------------------------------------------

