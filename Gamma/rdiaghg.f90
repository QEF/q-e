!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine rdiaghg (n, m, h, s, ldh, e, v)  
  !     =================
  !
  !   calculates eigenvalues and eigenvectors of the generalized problem
  !   Hv=eSv, with H symmetric matrix, S overlap matrix .
  !   On output both matrix are unchanged
  !   Uses LAPACK routines
  !-----------------------------------------------------------------------
#include "machine.h"
  use parameters
  use allocate 
#ifdef PARA
  use para
#endif
  implicit none  

  !  on INPUT
  integer :: n, m, ldh  
  ! dimension of the matrix to be diagonalized
  ! number of eigenstates to be calculated
  ! leading dimension of h, as declared in the calling pgm unit
  real(kind=DP) :: h (ldh, n)  
  ! matrix to be diagonalized
  real(kind=DP) :: s (ldh, n)  
  ! overlap matrix

  !  on OUTPUT
  real(kind=DP) :: e (n)  
  ! eigenvalues
  real(kind=DP) :: v (ldh, m)  
  ! eigenvectors (column-wise)

  !  LOCAL variables
  integer :: lwork, nb, ILAENV, mm, info  
  ! ILAENV returns optimal block size "nb"
  ! mm = number of calculated eigenvectors
  external ILAENV

  integer, pointer :: iwork (:), ifail (:)  
  real(kind=DP), pointer :: sdum (:,:), hdum (:,:),  work (:) 
  logical :: all_eigenvalues  


  call start_clock ('cdiaghg')  
  all_eigenvalues = m.eq.n  
  !
  !     check for optimal block size
  !
  nb = ILAENV (1, 'DSYTRD', 'U', n, - 1, - 1, - 1)  
  if (nb.lt.1.or.nb.ge.n) then  
     lwork = 8*n
  else  
     lwork = (nb + 3) * n  
  endif
  !
  ! allocate workspace
  !
  call mallocate(work, lwork)  
  call mallocate(sdum, ldh, n)  
  if (.not.all_eigenvalues) then  
     call mallocate(hdum, ldh, n)  
     call mallocate(iwork,  5 * n)  
     call mallocate(ifail,  n)  
  endif
  !
  ! input s and (see below) h are copied so that they are not destroyed
  !
  call DCOPY (ldh * n, s, 1, sdum, 1)  
#ifdef PARA
  !
  ! only the first processor diagonalize the matrix
  !
  if (me.eq.1) then  
#endif
#ifdef HAS_DSYGVX
     if (all_eigenvalues) then  
#endif
        !
        ! calculate all eigenvalues
        !
        call DCOPY (ldh * n, h, 1, v, 1)  
#ifdef AIX
        !
        ! there is a name conflict between essl and lapack ...
        !
        call DSYGV (1, v, ldh, sdum, ldh, e, v, ldh, n, work, lwork)
        info = 0
#else
        call DSYGV (1, 'V', 'U', n, v, ldh, sdum, ldh, e, work, &
             lwork, info)
#endif
#ifdef HAS_DSYGVX
     else  
        !
        ! calculate only m lowest eigenvalues
        !
        call DCOPY (ldh * n, h, 1, hdum, 1)  
        call DSYGVX (1, 'V', 'I', 'U', n, hdum, ldh, sdum, ldh, &
             0.0D0, 0.0D0, 1, m, 0.d0, mm, e, v, ldh, work, lwork, &
             iwork, ifail, info)
     endif
#endif
     call error ('rdiaghg', 'info =/= 0', abs (info) )  
#ifdef PARA
  endif
  !
  ! broadcast eigenvectors and eigenvalues to all other processors
  !
  call broadcast (n, e)  
  call broadcast (ldh * m, v)  
#endif
  !
  ! deallocate workspace
  !
  if (.not.all_eigenvalues) then  
     call mfree (ifail)  
     call mfree (iwork)  
     call mfree (hdum)  
  endif
  call mfree (sdum)  
  call mfree (work)  
  !
  call stop_clock ('cdiaghg')  
  !
  return  
end subroutine rdiaghg

