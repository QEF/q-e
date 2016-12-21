subroutine conjgrad(i_state,vstatecg,vstate_rcg,cstate,wcstate,fc,en)
use exciton 
use bse_basic_structures 
USE fft_custom_gwl 
USE io_global, ONLY : stdout,ionode
USE wvfct,    ONLY : npw
use bse_wannier, ONLY:num_nbndv,eps,maxit,eps_eig,lm_delta
USE mp,          ONLY :mp_barrier
USE mp_world,             ONLY : world_comm
USE constants,        ONLY: RYTOEV

implicit none

type(exc) :: a_exc
type(exc) :: x
type(exc) :: g
type(exc) :: h

type(v_state) :: vstatecg
type(v_state_r) :: vstate_rcg
type(c_state) :: cstate
type(c_state) :: wcstate
type(fft_cus) :: fc

!real(kind=DP) :: ha,hb,hc

real(kind=dp), intent(out) :: en

real(kind=DP) :: fp,fp_out,dfp,cg_lambda,dfg,step_rel
real(kind=DP) :: eigout,gg,dgg,gam,delta,hsquare,delta_eig,modh
integer :: i, is, i_state

logical :: restart
logical :: cgstatus

call start_clock('conjgrad')

restart=.false.
cgstatus=.false.

if(ionode) write(stdout,*) 'Conjugate gradient started.'
if(ionode) write(stdout,*) 'Looking for eigenvalue number:',i_state
if(ionode) write(stdout,*) 'eps=',eps
if(ionode) write(stdout,*) 'eps_eig=',eps_eig

!create initial random excitonic wavefunction vector a_exc 
call initialize_exc(a_exc)
a_exc%label=1
a_exc%npw=npw
a_exc%numb_v=num_nbndv(1)
allocate(a_exc%a(a_exc%npw,a_exc%numb_v))

call random_exc(a_exc)

!project into the conduction manifold
do is = 1,vstatecg%nspin
   call pc_operator_exc(a_exc,vstatecg,is)
enddo

!project out all the previous found state
call pout_operator_exc(a_exc,i_state)
 
!and normalize it
call normalize_exc(a_exc)

!initialize the other (internal) vectors
!within the cc iteration x is the gradient vector 
call initialize_exc(x)
x%label=1
x%npw=npw
x%numb_v=num_nbndv(1)
allocate(x%a(x%npw,x%numb_v))

!within the cg iteration g is the vector that stores (minus) the previous-step's gradient  
call initialize_exc(g)
g%label=1
g%npw=npw
g%numb_v=num_nbndv(1)
allocate(g%a(g%npw,g%numb_v))

!within the cg gradient iteration h is storing the search direction 
!check if really needed
call initialize_exc(h)
h%label=1
h%npw=npw
h%numb_v=num_nbndv(1)
allocate(h%a(h%npw,h%numb_v))


!Compute gradient at the initial guess position  
call exc_h_a(a_exc,x,vstatecg,vstate_rcg,cstate,wcstate,fc) 

!Compute function value fp at the initial guess position a_exc
call sproduct_exc(a_exc,x,fp)
!call sproduct_exc(x,x,hsquare)

x%a(1:x%npw,1:x%numb_v)=x%a(1:x%npw,1:x%numb_v)&
                          -fp*a_exc%a(1:a_exc%npw,1:a_exc%numb_v)

!Project the gradient into the conduction states manifold
!Remove any component along the eigenstates already found
do is = 1,vstatecg%nspin
   call pc_operator_exc(x,vstatecg,is)
enddo
call pout_operator_exc(x,i_state)

!Compute dfp, i.e. the magnitude of the gradient at the initial guess position a_exc
call sproduct_exc(x,x,dfp)
dfp=sqrt(dfp)
dfg=dfp
modh=dfp

!Set initial g and h values 
g%a(1:g%npw,1:g%numb_v)=-x%a(1:x%npw,1:x%numb_v)
h%a(1:h%npw,1:h%numb_v)=g%a(1:g%npw,1:g%numb_v)

!Start CG iterations
step_rel=lm_delta
delta=100
delta_eig=dfp-fp**2
i=1
do while ((i<=maxit).and.((delta>=eps).or.(delta_eig>eps_eig)))
   
   if(i>1)delta=abs(eigout-fp*RYTOEV)
!  note that delta_eig is calculated for each step within linmin
   
   eigout=fp*RYTOEV
   if(ionode) write(stdout,*) 'CG: eig#',i_state,'it=', i, 'Eig (eV)=',eigout 
!   if(ionode) write(stdout,*) 'CG: eig#',i_state,'it=', i, 'delta_eig=',delta_eig

!   if(ionode) write(stdout,*) 'CG: it=', i, 'Delta(eV)=',delta

   call linmin(i_state,a_exc,fp,dfp,h,modh,x,vstatecg,vstate_rcg,cstate,wcstate,fc,delta_eig,restart,i,cgstatus,step_rel)

   if(cgstatus) then   

      gg=dfg**2.d0 
      if (gg==0.d0) exit

      call sproduct_exc(g,x,dgg)
      dgg=dfp**2+dgg

      gam=max(dgg/gg,0.d0)

      g%a(1:g%npw,1:g%numb_v)=-x%a(1:x%npw,1:x%numb_v)
      dfg=dfp

      h%a(1:h%npw,1:h%numb_v)= g%a(1:g%npw,1:g%numb_v)+gam*h%a(1:h%npw,1:h%numb_v)
!      x%a(1:x%npw,1:x%numb_v)=h%a(1:h%npw,1:h%numb_v)
     
      call sproduct_exc(h,h,modh)
      modh=sqrt(modh)

   else 
   ! every 20 iterations restart the CG or if parabolic fit didn't work
   
 
      h%a(1:h%npw,1:h%numb_v)=-x%a(1:x%npw,1:x%numb_v)
      g%a(1:g%npw,1:g%numb_v)=-x%a(1:x%npw,1:x%numb_v)
      dfg=dfp
      modh=dfp


   endif
   
   i=i+1

   call mp_barrier(world_comm)
enddo

bse_spectrum(i_state)%a(1:bse_spectrum(i_state)%npw,1:bse_spectrum(i_state)%numb_v)=&
                         a_exc%a(1:a_exc%npw,1:a_exc%numb_v)


bse_spectrum(i_state)%e=fp*RYTOEV

en=fp


if(i==maxit)  then
   if(ionode) write(stdout,*) 'WARNING Conjugate gradient: Max iteration reached'
   if(ionode) write(stdout,*) 'Please increase the max iteration number or decrease accuracy' 
endif



!free memory
call free_memory_exc_a(a_exc)
call free_memory_exc_a(h)
call free_memory_exc_a(g)
call free_memory_exc_a(x)


call stop_clock('conjgrad')

if(ionode) write(stdout,*) 'Conjugate gradient ended.'
return
end subroutine
