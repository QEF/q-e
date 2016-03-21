subroutine sdescent(i_state,vstatesd,vstate_rsd,cstate,cstate_r,fc,en)

use exciton
use bse_basic_structures
USE fft_custom_gwl
USE io_global, ONLY : stdout,ionode
USE wvfct,    ONLY : npw
use bse_wannier, ONLY:num_nbndv,eps,lambda,eps_eig
USE mp,          ONLY :mp_barrier
USE mp_world,             ONLY : world_comm
USE constants,        ONLY: RYTOEV



implicit none

type(exc) :: a_in
type(exc) :: a_out
type(v_state) :: vstatesd
type(v_state_r) :: vstate_rsd
type(c_state) :: cstate
type(c_state_r) :: cstate_r
type(fft_cus) :: fc

real(kind=DP), intent(out) :: en
real(kind=DP) :: eigout
real(kind=DP) :: eig
real(kind=DP) ::delta,delta_eig,hsquare

integer :: it,is,i,i_state

call start_clock('sdescent')

!create a random excitonic wavefunction vector a_exc
!and normalize it
call initialize_exc(a_in)
a_in%label=50
a_in%npw=npw
a_in%numb_v=num_nbndv(1)
allocate(a_in%a(a_in%npw,a_in%numb_v))

call random_exc(a_in)
!project into the conduction manifold

do is = 1,vstatesd%nspin
   call pc_operator_exc(a_in,vstatesd,is)
enddo

!project out all the previous found state
call pout_operator_exc(a_in,i_state)

call normalize_exc(a_in)
CALL mp_barrier(world_comm)

call initialize_exc(a_out)
a_out%label=1
a_out%npw=npw
a_out%numb_v=num_nbndv(1)
allocate(a_out%a(a_out%npw,a_out%numb_v))

eig=0.d0
eigout=100.d0
delta_eig=100.d0
delta=100.d0

if(ionode) write(stdout,*) 'Steepest descent started.'
if(ionode) write(stdout,*) 'lambda=',lambda
if(ionode) write(stdout,*) 'eps',eps

it=0
!do while(((abs(delta))>=eps)) 
do while(((abs(delta))>=eps).or.(abs(delta_eig)>eps_eig)) 
!   write(*,*) 's descent, iteration, delta=',it,delta
! |a_out>=H|a_in>
   call exc_h_a(a_in,a_out,vstatesd,vstate_rsd,cstate,cstate_r,fc) 
   call mp_barrier(world_comm)
!   call normalize_exc(a_out)
   
!  eigout= <psi_(it)|H|psi_(it)>/<psi_(it)|psi_(it)>  
   call sproduct_exc(a_out,a_in,eigout)
   eigout=eigout*RYTOEV
   write(*,*) 'sd. eig# =',i_state, 'it=', it, 'E(eV)=', eigout 
  
!  check how good it is as an eigenstate 
   call sproduct_exc(a_out,a_out,hsquare)
   hsquare=hsquare*RYTOEV*RYTOEV

   delta_eig=hsquare-eigout**2
 
   eigout=eigout*RYTOEV
! compute |psi_(it+1)>
   a_out%a(1:a_out%npw,1:a_out%numb_v)=(1.d0+lambda*eigout)*a_in%a(1:a_in%npw,1:a_in%numb_v)&
                          -lambda*a_out%a(1:a_out%npw,1:a_out%numb_v)

   !project into the conduction manifold

   do is = 1,vstatesd%nspin
      call pc_operator_exc(a_out,vstatesd,is)
   enddo

   !project out all the previous found state
   call pout_operator_exc(a_out,i_state)
   call normalize_exc(a_out)

   a_in%a(1:a_out%npw,1:a_out%numb_v)= a_out%a(1:a_out%npw,1:a_out%numb_v)

   call mp_barrier(world_comm)

   it=it+1
   delta=eig-eigout
   eig=eigout
enddo

bse_spectrum(i_state)%a(1:bse_spectrum(i_state)%npw,1:bse_spectrum(i_state)%numb_v)=&
                         a_out%a(1:a_out%npw,1:a_out%numb_v)


bse_spectrum(i_state)%e=eigout
en=eigout

!if(ionode) write(stdout,*) 'Lowest eigenvalue=',eig


!free memory
call free_memory_exc_a(a_out)
call free_memory_exc_a(a_in)

call stop_clock('sdescent')
return
end subroutine sdescent
