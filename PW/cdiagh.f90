!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cdiagh (n, h, ldh, e, v)  
!-----------------------------------------------------------------------
!
!   calculates all the eigenvalues and eigenvectors of a complex
!   hermitean matrix H . On output, the matrix is unchanged
!
#include "machine.h"
  use parameters
  use allocate 
#ifdef PARA
  use para  
#endif
  implicit none
! on INPUT
  integer :: n, &                 ! dimension of the matrix to be diagonalized
             ldh                  ! leading dimension of h, as declared in 
                                  ! the calling pgm unit
  complex(kind=DP) :: &
           h (ldh, n)             ! matrix to be diagonalized
! on OUTPUT
  real(kind=DP) :: e (n)          ! eigenvalues
  complex(kind=DP) :: v (ldh, n)  ! eigenvectors (column-wise)
#ifdef AIX
! LOCAL variables (ESSL version)
  integer :: naux, i, j, ij  
  external ZHPEV  
  complex(kind=DP), pointer :: hp (:), aux (:)  

  call start_clock ('cdiagh')  
  naux = 4 * n  

  call mallocate(hp,  n * (n + 1) / 2)  
  call mallocate(aux, naux)  

! copy to upper triangular packed matrix
  ij = 0  
  do j = 1, n  
     do i = 1, j  
        ij = ij + 1  
        hp (ij) = h (i, j)  
     enddo
  enddo
#ifdef PARA
!cc        call pcdiagh( n, h, ldh, e, v, ij )
!cc        goto 10
!
!  only the first processor diagonalize the matrix
!
  if (me.eq.1) then  
#endif
     call ZHPEV (21, hp, e, v, ldh, n, aux, naux)  
#ifdef PARA
  endif
  call broadcast (n, e)  
  call broadcast (2 * ldh * n, v)  
#endif
10 call mfree (aux)  
  call mfree (hp)  
#else
#if defined(CRAYY)
! LOCAL variables (Cray Eispack/Scilib version)

  integer :: i, j, k, info  
!
  real(kind=DP) :: ar (ldh, n), ai (ldh, n), zr (ldh, n), zi (ldh, n)  
  ! real and imaginary part of  h(ldh,n) and of  v(ldh,n)
  ! (used as auxiliary arrays)
  real(kind=DP) :: rwork (2, ldh), work (ldh)  
  !
  call start_clock ('cdiagh')  
  do i = 1, n  
     do j = 1, ldh  
        ar (j, i) = DREAL (h (j, i) )  
        ai (j, i) = DIMAG (h (j, i) )  
     enddo
  enddo
  call ch (ldh, n, ar, ai, e, 1, zr, zi, work, work, rwork, info)  
  call error ('cdiagh', 'info =/= 0', abs (info) )  
  do i = 1, n  
     do j = 1, ldh  
        v (j, i) = DCMPLX (zr (j, i), zi (j, i) )  
     enddo
  enddo
#else
! LOCAL variables (LAPACK version)
  integer :: lwork, ILAENV, nb, info  
! ILAENV returns optimal block size "nb"
  real(kind=DP), pointer :: rwork ( : )  
  complex(kind=DP), pointer:: work(:)
!
  call start_clock ('cdiagh')  
!
!     check for the block size
!
  nb = ILAENV (1, 'ZHETRD', 'U', n, - 1, - 1, - 1)  
  if (nb.lt.1) nb = max (1, n)  
  if (nb.eq.1.or.nb.ge.n) then  
     lwork = 2 * n - 1  
  else  
     lwork = (nb + 1) * n  

  endif
#ifdef PARA
!
!  if scalapack library is present and we have just one pool
!  and the matrix is larger than 130 we use the scalapack driver
!
#ifdef T3D
  if (npool.eq.1.and.n.gt.130) then  
     call scala_cdiag (n, h, ldh, e, v, ldh)  
     goto 10  

  endif
#endif
!
!  else only the first processor diagonalize the matrix
!

  if (me.eq.1) then  
#endif
!
! allocate workspace
!
     call ZCOPY (n * ldh, h, 1, v, 1)  
     call mallocate(work, lwork)  
     call mallocate(rwork, (3 * n - 2) )  
     call ZHEEV ('V', 'U', n, v, ldh, e, work, lwork, rwork, info)  
     call error ('cdiagh', 'info =/= 0', abs (info) )  
! deallocate workspace
     call mfree (rwork)  
     call mfree (work)  
#ifdef PARA
  endif
  call broadcast (n, e)  
  call broadcast (2 * ldh * n, v)  
#endif

10 continue  
#endif
#endif

  call stop_clock ('cdiagh')  
  return  
end subroutine cdiagh

#ifdef MKL
! ILAENV is missing in the Intel Mathematical Kernel Library (mkl)
integer function ILAENV ()
  ILAENV=64
end function ILAENV
#endif
