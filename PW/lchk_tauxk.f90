!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
logical function lchk_tauxk (nvec, vec, trmat)  
  !-----------------------------------------------------------------------
  !
  !     This routine tests that the atomic coordinates and/or k-points
  !     are different (except for a lattice translation).
  !
  !     Last revision  15 June 1995  by A. Di Pomponio
  !
#include "machine.h"
  !
  use parameters
  use allocate 
  implicit none  
  !
  !     first the dummy variables
  !
  integer :: nvec  
  ! input: number of vectors (atom. pos. or k-p
  real(kind=DP) :: vec (3, nvec), trmat (3, 3)  
  ! input: cryst./cart. coord. of the vectors
  !               (atom. pos. or k-points)
  ! input: transf. matrix
  !        ( = bg , basis of the real-space lat
  !                 for atoms or
  !          = at , basis of the rec.-space lat
  !                 for k-points
  !
  !    here the local variables
  !
  integer :: nv1, nv2, kpol  
  ! first counter on vectors
  ! second counter on vectors
  ! counter on polarizations
  real(kind=DP), pointer :: vaux (:,:)
  real(kind=DP) :: vdf (3), accep  
  ! auxiliary vectors (atom. coord. or k-points
  !                    in cryst. units)
  ! auxiliary vector
  ! acceptance parameter  
  !
  !   Here, set the value of the acceptance parameter
  !
  parameter (accep = 1.0d-5)  
  !
  !
  !   allocate work space
  !
  call mallocate(vaux, 3 , nvec)  
  !
  !   The vectors are in cart. coordinates; they are transformed
  !   in crystall. units
  !
  do nv1 = 1, nvec  
     do kpol = 1, 3  
        vaux (kpol, nv1) = trmat (1, kpol) * vec (1, nv1) + trmat (2, &
             kpol) * vec (2, nv1) + trmat (3, kpol) * vec (3, nv1)
     enddo
  enddo
  !
  !   Test that all the atomic coordinates and/or k-points are different
  !   (if true, set the logical function  lchk_tauxk  as  .true.)
  !
  do nv1 = 1, nvec - 1  
     do nv2 = nv1 + 1, nvec  
        do kpol = 1, 3  
           vdf (kpol) = vaux (kpol, nv2) - vaux (kpol, nv1)  
           vdf (kpol) = abs (vdf (kpol) - float (nint (vdf (kpol) ) ) )  
        enddo
        if ( (vdf (1) .lt.accep) .and. (vdf (2) .lt.accep) .and. (vdf (3) &
             .lt.accep) ) then
           lchk_tauxk = .false.  
           goto 10  
        endif
     enddo
  enddo
  !
  lchk_tauxk = .true.  
  !
  !   deallocate work space
  !
10 call mfree(vaux)  
  return  
end function lchk_tauxk

