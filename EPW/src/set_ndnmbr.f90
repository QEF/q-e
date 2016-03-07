  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------
  subroutine set_ndnmbr ( pool, proc, procp, npool, ndlab)
  !----------------------------------------------------------------
  !
  !  create ndlab label from pool and proc numbers
  !
  !  The rule for deciding the node number is based on 
  !  the restriction set in startup.f90 that every pool 
  !  has the same number of procs.
  !
  !----------------------------------------------------------------
  implicit none
  character(len=3) :: ndlab
  integer :: pool, proc, procp, node, nprocs, npool 
  ! pool = 1,...,npool
  ! proc = 0,...,nproc_pool-1
  !
  nprocs = npool * procp
  !
  node = (pool-1)*procp + proc + 1
  !
  ndlab = '   '
  IF ( nprocs < 10 ) then
    WRITE (ndlab(1:1),'(i1)') node
  ELSEif ( nprocs < 100 ) then
! SP: Seems like QE change its convention to
    IF ( node < 10 ) then
!       ndlab = '0'
!       WRITE (ndlab(2:2),'(i1)') node
       WRITE (ndlab(1:1),'(i1)') node
    ELSE
       WRITE (ndlab(1:2),'(i2)') node
    ENDIF
  ELSE
    IF ( node < 10 ) then
 !      ndlab = '00'
!       WRITE (ndlab(3:3),'(i1)') node
       WRITE (ndlab(1:1),'(i1)') node
    ELSEif ( node < 100 ) then
 !      ndlab = '0'
!       WRITE (ndlab(2:3),'(i2)') node
       WRITE (ndlab(1:2),'(i2)') node
    ELSE
       WRITE (ndlab,'(i3)') node
    ENDIF
  ENDIF
  !
  end subroutine set_ndnmbr

