!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#ifdef FFTW
subroutine bidon_sun  
  stop 'cft_sun'  
end subroutine bidon_sun
#else
#ifdef SUN
!----------------------------------------------------------------------

subroutine cft_1 (f, m, n, nx, sgn, fout)  
  !     ===============
  !     driver routine for m 1d complex fft's of lenght n
  !     nx is the actual dimension of f (may differ from n)
  !     SUN (using sunperf library)
  !----------------------------------------------------------------------
#include "machine.h"
  use parameters, only : DP
  implicit none  

  integer :: m, n, nx, sgn  
  complex (kind=DP) :: f (nx * m), fout (nx * m)  
!
! Local variables
!
  integer :: on (2), naux1, isign, itype, i  
  parameter (naux1 = 20000)  
  real (kind=DP) :: aux1 (naux1, 2)  
  external zffti, zfftb, zfftf, zdscal
  data on / 0, 0 /  


  save on, aux1  
  isign = sign (1, sgn)  
  itype = abs (sgn)  

  if (itype.le.0.or.itype.gt.2) call error ('cft_1', 'wrong call', 1)

  if (n.ne.on (itype) ) then  
     call zffti (n, aux1 (1, itype) )  
     on (itype) = n  
  endif

  if (isign.eq. 1) then
    do i = 1, m  
       call zfftb ( n, f (1 + (i - 1) * nx), aux1 ( 1, itype) )  
    enddo
  else 
    do i = 1, m  
       call zfftf ( n, f (1 + (i - 1) * nx), aux1 ( 1, itype) )  
    enddo
    call zdscal ( nx*m, 1.d0/n, f, 1)  
  endif
  !
  !
  call zcopy ( nx*m, f, 1, fout, 1)  

  return  
end subroutine cft_1
!
!----------------------------------------------------------------------

subroutine cft_2 (f, mplane, n1, n2, nx1, nx2, sgn)  
  !     ===============
  !     driver routine for mplane 2d complex fft's of lenghts n1 and n2
  !     nx1=n1+1 is allowed (in order to avoid memory conflicts)
  !     for compatibility: nx2=n2, nx2 is not used
  !     SUN (using sunperf lib)
  !
  !----------------------------------------------------------------------
  !
#include "machine.h"
use parameters, only : DP
  implicit none  
  integer :: n1, n2, mplane, nx1, nx2, sgn  
  complex (kind=DP) :: f (nx1 * nx2 * mplane)  
!
! Local variables
!
  integer :: isign, itype, on1 (2), on2 (2), m, i, k, istrt, naux1, &
       naux2
  parameter (naux1 = 20000, naux2 = 10000)  
 real (kind=DP) :: aux1 (naux1, 2, 2), fj (naux2)  
  !
  !
  external zffti, zfftb, zfftf, zdscal
  data on1 / 0, 0 /, on2 / 0, 0 /  
  save on1, on2, aux1  
  !
  !
  isign = sign (1, sgn)  
  if (isign.ne. - 1.and.isign.ne.1) call error ('cft_2', 'wrong call', 1)
  itype = abs (sgn)  

  if (itype.le.0.or.itype.gt.2) call error ('cft_2', 'wrong call', &
       2)

  if (n2.ne.nx2) call error ('cft_2', 'no longer implemented', 1)  

  if (n1.ne.on1 (itype) ) then  
     call zffti (n1, aux1 (1, 1, itype) )  
     on1 (itype) = n1  
  endif

  if (n2.ne.on2 (itype) ) then  
     call zffti (n2, aux1 (1, 2, itype) )  
     on2 (itype) = n2  
  endif

  !  i - direction ...
  m = n2 * mplane  
  if ( isign.eq. 1 ) then
    do i = 1, m  
       call zfftb ( n1, f (1 + (i - 1) * nx1), aux1 (1, 1, itype) )
    enddo
  else
    do i = 1, m  
       call zfftf ( n1, f (1 + (i - 1) * nx1), aux1 ( 1, 1, itype) )  
    enddo
  endif  

  ! ... j-direction ...
  m = n1  
  if( isign.eq. 1) then
    do k = 1, mplane  
       istrt = 1 + (k - 1) * nx1 * n2  
       do i = 1, m  
          call zcopy (n2, f (istrt + i - 1), nx1, fj, 1)  
          call zfftb (n2, fj, aux1 (1, 2, itype) )  
          call zcopy (n2, fj, 1, f (istrt + i - 1), nx1)  
       enddo
    enddo
  else
    do k = 1, mplane  
       istrt = 1 + (k - 1) * nx1 * n2  
       do i = 1, m  
          call zcopy (n2, f (istrt + i - 1), nx1, fj, 1)  
          call zfftf (n2, fj, aux1 (1, 2, itype) )  
          call zcopy (n2, fj, 1, f (istrt + i - 1), nx1)  
       enddo
    enddo
    call zdscal ( nx1 * n2 * mplane, 1d0/(n1 * n2), f, 1)
  endif
   
  return  
end subroutine cft_2
#else
subroutine bidon_sun  
  stop 'cft_sun'  
end subroutine bidon_sun
#endif
#endif
