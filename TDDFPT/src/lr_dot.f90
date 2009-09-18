!-----------------------------------------------------------------------
!
! OBM
!  060608 gamma_only correction
!  160608 reduce --> mp_sum
function lr_dot(x,y)
  !---------------------------------------------------------------------
  ! Brent Walker, ICTP, 2004
  !---------------------------------------------------------------------
  ! ... wrapper for PWSCF linear response inner product routines 
  ! ... sums over the bands
  ! ... call for each k-point with arguments:
  ! ... call lr_dot(npw_k(ik),evc1(1,1,ik,1),1,evc1(1,1,ik,2),1)
  !---------------------------------------------------------------------
#include "f_defs.h"
  !
  use io_global,            only : stdout
  use kinds,                only : dp
  use klist,                only : nks
  !use lr_variables,         only : npw_k
  use realus,               only : npw_k
  use lsda_mod,             only : nspin
  use wvfct,                only : npwx,nbnd,wg
  use control_flags,        only : gamma_only
  use gvect,                only : gstart
  use mp,                   only : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE lr_variables,   ONLY : lr_verbosity
   USE io_global,      ONLY : stdout
 !
  implicit none
  !
  complex(kind=dp) :: x(npwx,nbnd,nks),y(npwx,nbnd,nks)
  complex(kind=dp) :: lr_dot
  complex(kind=dp) :: temp_k
  real(kind=dp) :: temp_gamma
  real(kind=dp) :: degspin
  integer :: ibnd
  integer :: ik
  real(kind=dp), external    :: DDOT
  complex(kind=dp), external :: ZDOTC
  !
  If (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_dot>")')
  endif
  call start_clock ('lr_dot')
  ! 
  lr_dot=(0.0d0,0.0d0)
  !
  temp_gamma=0.0d0
  temp_k=(0.0d0,0.0d0)
  !
  degspin=2.0d0
  if(nspin==2) degspin=1.0d0
  !
  if(gamma_only) then
     call lr_dot_gamma()
     lr_dot=cmplx(temp_gamma,0.0d0,dp)
  else
     call lr_dot_k()
     lr_dot=temp_k
  endif
  !
  lr_dot=lr_dot/degspin
  !
  if (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dot>")') 
  call stop_clock ('lr_dot')
  !
  return
  !
contains
  !
  subroutine lr_dot_gamma
    !
    do ibnd=1,nbnd
       !
       temp_gamma = temp_gamma + 2.D0*wg(ibnd,1)*DDOT(2*npw_k(1),x(:,ibnd,1),1,y(:,ibnd,1),1)
       if (gstart==2) temp_gamma = temp_gamma - wg(ibnd,1)*dble(x(1,ibnd,1))*dble(y(1,ibnd,1))
       !
    enddo
    !
#ifdef __PARA
    !call reduce(1,temp_gamma)
    call mp_sum(temp_gamma, intra_pool_comm)
#endif
    !
    return
  end subroutine lr_dot_gamma
  !
  subroutine lr_dot_k
    !
    do ik=1,nks
       do ibnd=1,nbnd
          !
          temp_k=temp_k+wg(ibnd,ik)*ZDOTC(npw_k(ik),x(1,ibnd,ik),1,y(1,ibnd,ik),1)
          !
       enddo
    enddo
#ifdef __PARA
    !call poolreduce(2,temp_k)
    call mp_sum(temp_k,inter_pool_comm)
    !call reduce(2,temp_k)
    call mp_sum(temp_k, intra_pool_comm)
#endif
    !
    return
  end subroutine lr_dot_k
  !
end function lr_dot
!-----------------------------------------------------------------------
