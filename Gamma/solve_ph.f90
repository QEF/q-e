!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_ph
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE io_global,      ONLY : stdout
  use pwcom
  USE wavefunctions,  ONLY : evc
  USE rbecmod,        ONLY : becp, becp_
  use cgcom
#ifdef __PARA
  use para
#endif
  integer :: nu, i, ibnd, jbnd, info, iter, mode_done, kpoint
  real(kind=DP), allocatable ::  diag(:)
  complex(kind=DP), allocatable :: gr(:,:), h(:,:), work(:,:)
  real(kind=DP), allocatable :: overlap(:,:)
  logical :: orthonormal, precondition, startwith0
  external A_h
  !
  call start_clock('solve_ph')
  !
  allocate ( becp( nkb,nbnd), becp_(nkb,nbnd) )
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
        diag(i) = 1.0/max(1.d0,g2kin(i))
     end do
     call zvscal(npw,npwx,nbnd,diag,evc,work)
     call pw_gemm ('Y',nbnd, nbnd, npw, work, npwx, evc, npwx, overlap, nbnd)
     call DPOTRF('U',nbnd,overlap,nbnd,info)
     if (info.ne.0) call errore('solve_ph','cannot factorize',info)
  end if
  !
  WRITE( stdout,'(/" ***  Starting Conjugate Gradient minimization",   &
       &            9x,"***")')
  !
  !  check if a restart file exists
  !
  open (unit=iunres,file='restartph',form='formatted',status='unknown')
  read (iunres,*,err=1,end=1) mode_done
  read (iunres,*,err=1,end=1) dyn
  close(unit=iunres)
  print '("  Phonon: modes up to mode ",i3," are done")',       &
       &     mode_done
  go to 2
1 close(unit=iunres)
  !  restart failed or file not found
  call  dynmat_init
  mode_done=0
2 continue
  !
  do nu = 1, nmodes
     if ( has_equivalent((nu-1)/3+1).eq.1) then
        ! calculate only independent modes
        WRITE( stdout,'(" ***  mode # ",i3," : using symmetry")') nu
        goto 10
     end if
     if ( nu.le.mode_done) then
        ! do not recalculate modes already done
        WRITE( stdout,'(" ***  mode # ",i3," : using previous run")') nu
        goto 10
     end if
     if ( asr .and. (nu-1)/3+1.eq.nasr ) then
        ! impose ASR on last atom instead of calculating mode
        WRITE( stdout,'(" ***  mode # ",i3," : using asr")') nu
        goto 10
     end if
     ! calculate |b> = dV/dtau*psi
     call dvpsi_kb(kpoint,nu)
     ! initialize delta psi
     startwith0=.true.
     dpsi(:,:) = (0.d0, 0.d0)
     ! solve the linear system
     ! NB: dvpsi is used also as work space and is destroyed by cgsolve
     call cgsolve (A_h,npw,evc,npwx,nbnd,overlap,nbnd, &
                   orthonormal,precondition,diag,      &
                   startwith0,et(1,kpoint),dvpsi,gr,h, &
                   dvpsi,work,niter_ph,tr2_ph,iter,dpsi)
     ! < DeltaPsi | DeltaV | Psi > contribution to the dynamical matrix
     call drhodv(nu)
     ! save partial result
#ifdef __PARA
     if (me.eq.1) then
#endif
        open (unit=iunres,file='restartph',form='formatted',status='unknown')
        write(iunres,*) nu
        write(iunres,*) dyn
        close(unit=iunres)
#ifdef __PARA
     end if
#endif
     WRITE( stdout,'(" ***  mode # ",i3," : ",i3," iterations")')  &
          &          nu, iter
10   continue
  end do
  !
  deallocate(h)
  deallocate(gr)
  deallocate(overlap)
  deallocate(work)
  deallocate(diag)
  deallocate(becp, becp_)
  !
  call stop_clock('solve_ph')
  !
  return
end subroutine solve_ph

subroutine set_asr(nat,nasr,dyn)
  !
  ! Impose Acoustic Sum Rule on the dynamical matrix
  ! We assume that (3*nat-1) columns have been calculated
  ! and that the missing column corresponds to atom nasr
  !
  implicit none
  integer nat, nasr
  real(kind=8) :: dyn(3*nat,3*nat)
  !
  integer na, nb, i,j
  real(kind=8) :: sum

  if (nasr.le.0 .or. nasr.gt.nat) return
  do j=1,3
     do i=1,3
        do nb=1,nat
           sum=0.d0
           do na=1,nat
              if (na.ne.nasr) sum = sum + dyn(3*(na-1)+i,3*(nb-1)+j)
           end do
           dyn(3*(nasr-1)+i,3*(nb-1)+j)= -sum
        end do
     end do
  end do

  return
end subroutine set_asr
