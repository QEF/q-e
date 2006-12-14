!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE io_global,      ONLY : stdout
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE becmod, ONLY: rbecp
  use cgcom
  !
  implicit none
  !
  integer :: ipol, nrec, i, ibnd, jbnd, info, iter, kpoint
  real(DP), allocatable ::diag(:)
  complex(DP), allocatable :: gr(:,:), h(:,:), work(:,:)
  real(DP), allocatable :: overlap(:,:)
  logical :: orthonormal, precondition,startwith0,here
  character(len=7) :: fildwf, filbar
  external A_h
  !
  call start_clock('solve_e')
  !
  allocate ( rbecp( nkb,nbnd) )
  allocate ( diag( npwx) )
  allocate ( overlap( nbnd, nbnd) )
  allocate ( work( npwx, nbnd) )
  allocate ( gr  ( npwx, nbnd) )
  allocate ( h   ( npwx, nbnd) )
  !
  kpoint = 1
  do i = 1,npw
     g2kin(i) = ( (xk(1,kpoint)+g(1,igk(i)))**2 +                   &
                  (xk(2,kpoint)+g(2,igk(i)))**2 +                   &
                  (xk(3,kpoint)+g(3,igk(i)))**2 ) * tpiba2
  end do
  !
  orthonormal = .false.
  precondition= .true.
 !
  if (precondition) then
     do i = 1,npw
        diag(i) = 1.0d0/max(1.d0,g2kin(i))
     end do
     call zvscal(npw,npwx,nbnd,diag,evc,work)
     call pw_gemm ('Y',nbnd, nbnd, npw, work, npwx, evc, npwx, overlap, nbnd)
     call DPOTRF('U',nbnd,overlap,nbnd,info)
     if (info.ne.0) call errore('solve_e','cannot factorize',info)
  end if
  !
  WRITE( stdout,'(/" ***  Starting Conjugate Gradient minimization",   &
       &            9x,"***")')
  nrec=0
  !
  do ipol = 1,3
     !  read |b> = dV/dtau*psi
     iubar=ipol
     write(filbar,'("filbar",i1)') ipol
     call seqopn (iubar,filbar,'unformatted',here)
     if (.not.here) call errore('solve_e','file '//filbar//          &
          &        'mysteriously vanished',ipol)
     read (iubar) dvpsi
     close(unit=iubar,status='keep')
     !
     iudwf=10+ipol
     write(fildwf,'("fildwx",i1)') ipol
     call  seqopn (iudwf,fildwf,'unformatted',here)
!!!         if (.not.here) then
     !  calculate Delta*psi  (if not already done)
     dpsi(:,:) = (0.d0, 0.d0)
     startwith0= .true.
!!!         else
     !  otherwise restart from Delta*psi that is found on file
!!!            read(iudwf) dpsi
!!!         end if
     call cgsolve (A_h,npw,evc,npwx,nbnd,overlap,nbnd, &
                   orthonormal,precondition,diag,      &
                   startwith0,et(1,kpoint),dvpsi,gr,h, &
                   dvpsi,work,niter_ph,tr2_ph,iter,dpsi)
     !  write Delta*psi for an electric field
     rewind (iudwf)
     write (iudwf) dpsi
     close(unit=iudwf)
     !
     WRITE( stdout,'(" ***  pol. # ",i3," : ",i3," iterations")')  &
          &              ipol, iter
  end do
  !
  deallocate(h)
  deallocate(gr)
  deallocate(overlap)
  deallocate(work)
  deallocate(diag)
  deallocate(rbecp)
  !
  call stop_clock('solve_e')
  !
  return
end subroutine solve_e
