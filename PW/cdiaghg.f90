!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine cdiaghg (n, m, h, s, ldh, e, v)
  !     =================
  !
  !   calculates eigenvalues and eigenvectors of the generalized problem
  !   Hv=eSv, with H hermitean matrix, S overlap matrix .
  !   On output both matrix are unchanged
  !
  !   LAPACK version - may use both ZHEGV and ZHEGVX
  !   ZHEGVX should be faster but it is not available on many machines
  !   define HAS_ZHEGVX in the preprocessing options to use ZHEGVX
  !-----------------------------------------------------------------------
#include "machine.h"
  use parameters
#ifdef __PARA
  use para
#endif
  implicit none

  !  on INPUT
  integer :: n, m, ldh
  ! dimension of the matrix to be diagonalized
  ! number of eigenstates to be calculate
  ! leading dimension of h, as declared in the calling pgm unit
  complex(kind=DP) :: h (ldh, n)
  ! matrix to be diagonalized
  complex(kind=DP) :: s (ldh, n)
  ! overlap matrix

  !  on OUTPUT
  real(kind=DP) :: e (n)
  ! eigenvalues
  complex(kind=DP) :: v (ldh, m)
  ! eigenvectors (column-wise)

  !  LOCAL variables
  integer :: lwork, nb, ILAENV, mm, info
  ! ILAENV returns optimal block size "nb"
  ! mm = number of calculated eigenvectors
  external ZCOPY, ZHEGV, ILAENV, error

  integer, allocatable :: iwork (:), ifail (:)
  real(kind=DP), allocatable :: rwork (:)
  complex(kind=DP), allocatable :: sdum (:,:), hdum (:,:),  work (:)
  logical :: all_eigenvalues


  call start_clock ('cdiaghg')
#ifdef __PARA
#ifdef __T3E
  !
  !   NB: 150 has been determined empirically on the T3E as the point
  !       where it is convenient to use a parallel routines.
  !
  if (npool.eq.1.and.n.gt.150) then
     call scala_cdiaghg (n, h, ldh, s, ldh, e, v, ldh)
     call stop_clock ('cdiaghg')
     return
  endif
#endif
#endif
#ifdef HAS_ZHEGVX
  all_eigenvalues = m.eq.n
#else
  all_eigenvalues = .true.
#endif
  !
  !     check for optimal block size
  !
  nb = ILAENV (1, 'ZHETRD', 'U', n, - 1, - 1, - 1)
  if (nb.lt.1) nb = max (1, n)
  if (nb.eq.1.or.nb.ge.n) then
     lwork = 2 * n - 1
  else
     lwork = (nb + 1) * n
  endif
  !
  ! allocate workspace
  !
  allocate (work( lwork))    
  allocate (sdum( ldh, n))    
  if (all_eigenvalues) then
     allocate (rwork( (3 * n - 2) ))    
  else
     allocate (rwork( 7 * n))    
     allocate (hdum( ldh, n))    
     allocate (iwork(  5 * n))    
     allocate (ifail(  n))    
  endif
  !
  ! input s and (see below) h are copied so that they are not destroyed
  !
  call ZCOPY (ldh * n, s, 1, sdum, 1)
#ifdef __PARA
  !
  ! only the first processor diagonalize the matrix
  !
  if (me.eq.1) then
#endif
     if (all_eigenvalues) then
        !
        ! calculate all eigenvalues
        !
        call ZCOPY (ldh * n, h, 1, v, 1)
        call ZHEGV (1, 'V', 'U', n, v, ldh, sdum, ldh, e, work, &
             lwork, rwork, info)
     else
#ifdef HAS_ZHEGVX
        !
        ! calculate only m lowest eigenvalues
        !
        call ZCOPY (ldh * n, h, 1, hdum, 1)
        call ZHEGVX (1, 'V', 'I', 'U', n, hdum, ldh, sdum, ldh, &
             0.0D0, 0.0D0, 1, m, 0.d0, mm, e, v, ldh, work, lwork, rwork, &
             iwork, ifail, info)
#endif
     endif
     call errore ('cdiaghg', 'info =/= 0', abs (info) )
#ifdef __PARA
  endif
  !
  ! broadcast the eigenvectors and the eigenvalues
  !
  call broadcast (n, e)
  call broadcast (2 * ldh * m, v)
#endif
  !
  ! deallocate workspace
  !
  if (.not.all_eigenvalues) then
     deallocate (ifail)
     deallocate (iwork)
     deallocate (hdum)
  endif
  deallocate (sdum)
  deallocate (rwork)
  deallocate (work)
  !
  call stop_clock ('cdiaghg')
  !
  return
end subroutine cdiaghg

