subroutine linmin(i_state,a_exc,fpa,dfa,h,modh,x,vstatecg,vstate_rcg,cstate,wcstate,fc,delta_eig,restart,iter,cgstatus,step_rel)
! This subroutine approximates the minimum of the functional
! f=<a_exc|Hexc|a_exc>-lambda*(<a_exc|a_exc>-1)
! along a given direction |h> with the vertex of a 
! parabola passing through the points (a_exc,fpa), and (b_exc=a_exc+delta*grad f,fpb) and having
! slope -dfa at the a_exc point.
! it updates the vector |a_exc> setting it at the 
! found minimum position, and it updates the value of the functional fpa, and returns its 
! gradient, |x>, at the new position  

USE exciton 
USE fft_custom_gwl 
USE io_global, ONLY : stdout,ionode
USE wvfct,    ONLY : npw
use bse_wannier, ONLY:num_nbndv,lm_delta,eps,maxit
USE bse_basic_structures  
!USE mp,          ONLY :mp_barrier
!USE constants,        ONLY: RYTOEV

implicit none
type(exc), intent(inout) :: a_exc
type(exc), intent(out) :: x
type(exc), intent(in) :: h
real(kind=DP), intent(inout) :: dfa,fpa
real(kind=DP), intent(inout) :: modh
real(kind=DP), intent(inout):: step_rel
logical, intent(inout) ::cgstatus

type(v_state) :: vstatecg
type(v_state_r) :: vstate_rcg
type(c_state) :: cstate
type(c_state) :: wcstate
type(fft_cus) :: fc

real(kind=DP) :: fp_out, fpb, fpv,delta_new,b,c,a,delta_eig,hsquare,fpc

type(exc) :: b_exc
type(exc) :: c_exc
type(exc) :: v_exc
type(exc) :: xb
type(exc) :: xv

integer :: is,j,i,i_state,iter
real(kind=DP):: step,fptmp
real(kind=DP) ::e_expected,pb,pc
logical ::restart,fit_3p
integer :: fit_ok

real(kind=DP) :: deri,sca

call start_clock('linmin')
fit_3p=.true.


!initialize internal vectors
!second point for the qudratic fit
call initialize_exc(b_exc)
b_exc%label=1
b_exc%npw=npw
b_exc%numb_v=num_nbndv(1)
allocate(b_exc%a(b_exc%npw,b_exc%numb_v))

if(fit_3p) then
!this variable is needed if we want to find the 
!parabola equations given three points
   call initialize_exc(c_exc)
   c_exc%label=1
   c_exc%npw=npw
   c_exc%numb_v=num_nbndv(1)
   allocate(c_exc%a(c_exc%npw,c_exc%numb_v))
endif

!gradient at the second point for the qudratic fit
call initialize_exc(xb)
xb%label=1
xb%npw=npw
xb%numb_v=num_nbndv(1)
allocate(xb%a(xb%npw,xb%numb_v))

!position of the vertex of the qudratic fit
call initialize_exc(v_exc)
v_exc%label=1
v_exc%npw=npw
v_exc%numb_v=num_nbndv(1)
allocate(v_exc%a(v_exc%npw,v_exc%numb_v))

!gradient at the position of the vertex of the qudratic fit
call initialize_exc(xv)
xv%label=1
xv%npw=npw
xv%numb_v=num_nbndv(1)
allocate(xv%a(xv%npw,xv%numb_v))

!for debug purpose plot the energies along the search line
!if(.not.cgstatus) then
!   do j=0,50
!    step=step_rel/25
!      b_exc%a(1:b_exc%npw,1:b_exc%numb_v)=a_exc%a(1:a_exc%npw,1:a_exc%numb_v)+&
!                             step*dble(j)*h%a(1:h%npw,1:h%numb_v) 
!      do is = 1,vstatecg%nspin
!         call pc_operator_exc(b_exc,vstatecg,is)
!      enddo
!      call pout_operator_exc(b_exc,i_state)
!
!      call normalize_exc(b_exc)
!
!!compute the function (and the gradient) at the new b_exc position
!      call exc_h_a(b_exc,xb,vstatecg,vstate_rcg,fc) 
!      call sproduct_exc(b_exc,xb,fpb)
!
!      write(stdout,*) 'en_prof', iter, dble(j), fpb
!   enddo
!endif



!take one step along the search direction
b_exc%a(1:b_exc%npw,1:b_exc%numb_v)=a_exc%a(1:a_exc%npw,1:a_exc%numb_v)+step_rel*h%a(1:h%npw,1:h%numb_v) 


!project into the conduction manifold,remove any component along previous
!eigenstates, and normalize b_exc
do is = 1,vstatecg%nspin
   call pc_operator_exc(b_exc,vstatecg,is)
enddo

call pout_operator_exc(b_exc,i_state)!DEBUG

call normalize_exc(b_exc)

!compute the function at the new b_exc position
call exc_h_a(b_exc,xb,vstatecg,vstate_rcg,cstate,wcstate,fc) 
call sproduct_exc(b_exc,xb,fpb)


call exc_h_a(a_exc,xb,vstatecg,vstate_rcg,cstate,wcstate,fc)
call sproduct_exc(a_exc,xb,sca)
xb%a(1:xb%npw,1:xb%numb_v)=xb%a(1:xb%npw,1:xb%numb_v)&
                          -sca*a_exc%a(1:b_exc%npw,1:b_exc%numb_v)
call sproduct_exc(xb,h,deri)
deri=deri*2.d0

!compute the quadratic fit
if(.not.fit_3p) then 
   b=deri!-dfa DEBUG
   c=fpa
   pb=step_rel
   a= (fpb-b*pb-c)/((pb)**2.d0)

endif

if(fit_3p) then
! find the third point C and compute there the function fpc
   c_exc%a(1:c_exc%npw,1:c_exc%numb_v)=a_exc%a(1:a_exc%npw,1:a_exc%numb_v)+0.5d0*step_rel*h%a(1:h%npw,1:h%numb_v) 
   do is = 1,vstatecg%nspin
      call pc_operator_exc(c_exc,vstatecg,is)
   enddo

   call pout_operator_exc(c_exc,i_state)

   call normalize_exc(c_exc)

!  xv is used here a tmp variable
   call exc_h_a(c_exc,xv,vstatecg,vstate_rcg,cstate,wcstate,fc) 
   call sproduct_exc(c_exc,xv,fpc)

   pb=step_rel
   pc=0.5d0*step_rel
   
   c=fpa 
   a=(fpb*pc-fpc*pb-fpa*(pc-pb))/(pc*pb*(pb-pc))
   b=(-fpb*pc*pc+fpc*pb*pb+fpa*(pc*pc-pb*pb))/(pc*pb*(pb-pc))

endif


delta_new=-b/(2.d0*a)

!compute the coordinate of the parabola vertex
v_exc%a(:,:)=a_exc%a(:,:)+delta_new*h%a(:,:) 

!write(stdout,*) 'delta_new',iter, delta_new

!compute the expected value for the minimum from the parabolic fit
e_expected=-b**2/(4*a)+c

!project into the conduction manifold, remove any component along previous
!eigenstates and normalize v_exc
do is = 1,vstatecg%nspin
   call pc_operator_exc(v_exc,vstatecg,is)
enddo

call pout_operator_exc(v_exc,i_state)
call normalize_exc(v_exc)

!compute the function (and the gradient) at the new v_exc position
call exc_h_a(v_exc,xv,vstatecg,vstate_rcg,cstate,wcstate,fc) 
call sproduct_exc(v_exc,xv,fpv)
xv%a(1:xv%npw,1:xv%numb_v)=xv%a(1:xv%npw,1:xv%numb_v)&
                          -fpv*v_exc%a(1:v_exc%npw,1:v_exc%numb_v)

!write(stdout,*) 'min', iter, e_expected
!write(stdout,*) 'fpa', iter, fpa
!write(stdout,*) 'fpb', iter, fpb
!write(stdout,*) 'fpv', iter, fpv
!write(stdout,*) 'FIT PARAMETER A', iter, a
!write(stdout,*) 'FIT PARAMETER B', iter, b
!write(stdout,*) 'FIT PARAMETER C', iter, c

!find where our functional is minimum among the last three points

if ((fpv<=fpa).and.(fpv<=fpb)) then
!     the quadratic fit found a good position update the information
   fit_ok=1
   cgstatus=.true.
   a_exc%a(:,:)=v_exc%a(:,:)
!   write(stdout,*) 'Quadratic fit: V new point, fpv=',fpv, 'fpa=',fpa  
   fpa=fpv      
   x%a(1:x%npw,1:x%numb_v)=xv%a(1:xv%npw,1:xv%numb_v)
else
!take a steepest descent step and restart the cg
   write(stdout,*) 'WARNING:I restart the cg'
   cgstatus=.false.
   fit_ok=0
   j=1
   fptmp=fpa+0.1d0
   do while (fptmp>=fpa)
      step=(0.3d0)*(0.1d0**j)
      b_exc%a(:,:)=a_exc%a(:,:)+step_rel*step*h%a(:,:) 

!  project into the conduction manifold remove any component along previous
!  eigenstates, and normalize b_exc
      do is = 1,vstatecg%nspin
         call pc_operator_exc(b_exc,vstatecg,is)
      enddo
      call pout_operator_exc(b_exc,i_state)
      call normalize_exc(b_exc)

!  compute the function (and the gradient) at the new b_exc position
      call exc_h_a(b_exc,xb,vstatecg,vstate_rcg,cstate,wcstate,fc) 
      call sproduct_exc(b_exc,xb,fptmp)
      xb%a(1:xb%npw,1:xb%numb_v)=xb%a(1:xb%npw,1:xb%numb_v)&
                          -fptmp*b_exc%a(1:b_exc%npw,1:b_exc%numb_v)
  
      if(fptmp>=fpa) then! look in the opposite direction
         b_exc%a(:,:)=a_exc%a(:,:)-step_rel*step*h%a(:,:)
!        project into the conduction manifold remove any component along previous
!        eigenstates, and normalize b_exc
         do is = 1,vstatecg%nspin
            call pc_operator_exc(b_exc,vstatecg,is)
         enddo
         call pout_operator_exc(b_exc,i_state)
         call normalize_exc(b_exc)

!  compute the function (and the gradient) at the new b_exc position
         call exc_h_a(b_exc,xb,vstatecg,vstate_rcg,cstate,wcstate,fc) 
         call sproduct_exc(b_exc,xb,fptmp)
         xb%a(1:xb%npw,1:xb%numb_v)=xb%a(1:xb%npw,1:xb%numb_v)&
                          -fptmp*b_exc%a(1:b_exc%npw,1:b_exc%numb_v)
    
      endif
      j=j+1
   enddo 
   
   a_exc%a(:,:)=b_exc%a(:,:)
   fpa=fptmp      
   x%a(1:x%npw,1:x%numb_v)=xb%a(1:xb%npw,1:xb%numb_v)
endif   


!project into the conduction manifold
do is = 1,vstatecg%nspin
   call pc_operator_exc(x,vstatecg,is)
enddo
call pout_operator_exc(x,i_state)
!compute dfa, i.e. the magnitude of the gradient at the position a_exc
call sproduct_exc(x,x,dfa)
dfa=sqrt(dfa)


!compute delta_eig (using as a tmp variable b)
call exc_h_a(a_exc,b_exc,vstatecg,vstate_rcg,cstate,wcstate,fc) 

call sproduct_exc(b_exc,b_exc,hsquare)

delta_eig=hsquare-fpa**2

! set the magnitude of the next step 
if(cgstatus) then 
   step_rel=2.d0*delta_new 
else
   step_rel=(0.3d0)*(0.1d0**(j-1))
endif

!write(stdout,*) 'fit_ok', iter, fit_ok
!free memory
call free_memory_exc_a(b_exc)
if(fit_3p) call free_memory_exc_a(c_exc)
call free_memory_exc_a(v_exc)
call free_memory_exc_a(xb)
call free_memory_exc_a(xv)

call stop_clock('linmin')
return
end subroutine 
