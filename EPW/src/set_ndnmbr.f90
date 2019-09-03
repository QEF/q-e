  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------
  SUBROUTINE set_ndnmbr(pool, proc, procp, npool, ndlab)
  !----------------------------------------------------------------
  !!
  !!  create ndlab label from pool and proc numbers
  !!
  !!  The rule for deciding the node number is based on 
  !!  the restriction set in startup.f90 that every pool 
  !!  has the same number of procs.
  !!
  !----------------------------------------------------------------
  IMPLICIT NONE
  ! 
  CHARACTER(LEN = 3), INTENT(out) :: ndlab 
  !! Label
  INTEGER, INTENT(in) :: pool
  !! Pool = 1,..., npool
  INTEGER, INTENT(in) :: proc
  !! Processor = 0,..., nproc_pool-1
  INTEGER, INTENT(in) :: procp
  !! 
  INTEGER, INTENT(in) :: npool
  !! 
  ! Local variables
  INTEGER :: node
  !! Number of nodes 
  INTEGER :: nprocs
  !! 
  !
  nprocs = npool * procp
  !
  node = (pool - 1) * procp + proc + 1
  !
  ndlab = '   '
  IF (nprocs < 10) THEN
    WRITE(ndlab(1:1), '(i1)') node
  ELSEIF (nprocs < 100 ) then
    IF (node < 10) THEN
      WRITE(ndlab(1:1), '(i1)') node
    ELSE
      WRITE(ndlab(1:2), '(i2)') node
    ENDIF
  ELSEIF (nprocs < 100) THEN
    IF (node < 10) THEN
      WRITE(ndlab(1:1), '(i1)') node
    ELSEIF (node < 100) THEN
      WRITE(ndlab(1:2), '(i2)') node
    ELSE
      WRITE(ndlab(1:3), '(i3)') node
    ENDIF
  ELSE
    IF (node < 10) THEN
      WRITE(ndlab(1:1), '(i1)') node
    ELSEIF (node < 100) THEN
      WRITE(ndlab(1:2), '(i2)') node
    ELSEIF (node < 1000) THEN
      WRITE(ndlab(1:3), '(i3)') node
    ELSE
      WRITE(ndlab, '(i4)') node
    ENDIF
  ENDIF
  !
  !----------------------------------------------------------------
  END SUBROUTINE set_ndnmbr
  !----------------------------------------------------------------

