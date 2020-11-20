subroutine exc_h_a(a_in,a_out,vstate,vstate_r,cstate,wcstate,fc) 
!this subroutine applies the excitonic hamiltonian exc_h_a on a given vector
!(a_in) and returns the transformed vector (a_out)
!if l_gtrick==.true. uses vstate only and does not use vstate_r

use exciton
use bse_basic_structures

use bse_wannier, ONLY:num_nbndv,l_truncated_coulomb, &
                      truncation_radius,l_fullbse,l_tdhf,l_lf,l_rpa,&
                      l_contraction
use pwcom
USE wvfct,     ONLY : npwx
!use io_files,  ONLY : find_free_unit,diropn
use io_files,  ONLY : diropn
USE io_global, ONLY : stdout,ionode
USE fft_custom_gwl
USE mp,          ONLY :mp_barrier
USE mp_world,             ONLY : world_comm
USE contract_w

implicit none
INTEGER, EXTERNAL :: find_free_unit

type(exc) :: a_in
type(exc) :: a_out ! to be initialized and allocated outside this subroutine
type(v_state) :: vstate
type(v_state_r) :: vstate_r
type(c_state) :: cstate
type(c_state) :: wcstate
type(fft_cus) :: fc

type(exc) :: a_excdiago,a_exchange 
type(exc):: a_dirv
type(exc):: a_dirw
type(exc):: a_rot


integer iuv
logical exst
logical debug

!variables introduced for debug purposes
real(kind=DP) prod_test1 
real(kind=DP) prod_test2 

call start_clock('exc_h_a')
debug=.false.

  if(debug) then 
     write(stdout,*) 'Starting exc_h_a subroutine'
  endif

  call mp_barrier(world_comm)

! initialize and nullify all needed excitonic vectors
  call initialize_exc(a_excdiago)
  a_excdiago%label=2
  a_excdiago%npw=npw
  a_excdiago%numb_v=num_nbndv(1)
  allocate(a_excdiago%a(a_excdiago%npw,a_excdiago%numb_v))
  a_excdiago%a(1:a_excdiago%npw,1:a_excdiago%numb_v)=dcmplx(0.d0,0.d0)

  call initialize_exc(a_exchange)
  a_exchange%label=3
  a_exchange%npw=npw
  a_exchange%numb_v=num_nbndv(1)
  allocate(a_exchange%a(a_exchange%npw,a_exchange%numb_v))
  a_exchange%a(1:a_exchange%npw,1:a_exchange%numb_v)=dcmplx(0.d0,0.d0)

  call initialize_exc(a_dirv)
  a_dirv%label=4
  a_dirv%npw=npw
  a_dirv%numb_v=num_nbndv(1)
  allocate(a_dirv%a(a_dirv%npw,a_dirv%numb_v))
  a_dirv%a(1:a_dirv%npw,1:a_dirv%numb_v)=dcmplx(0.d0,0.d0)

  call initialize_exc(a_dirw)
  a_dirw%label=6
  a_dirw%npw=npw
  a_dirw%numb_v=num_nbndv(1)
  allocate(a_dirw%a(a_dirw%npw,a_dirw%numb_v))
  a_dirw%a(1:a_dirw%npw,1:a_dirw%numb_v)=dcmplx(0.d0,0.d0)

! apply the diagonal part of the excitonic Hamiltonian to a copy of a_exc
  a_excdiago%a(1:a_excdiago%npw,1:a_excdiago%numb_v)=a_in%a(1:a_in%npw,1:a_in%numb_v)
  call diago_exc(a_excdiago,vstate,cstate,wcstate)

  if(debug) then 
     write(stdout,*) 'Diagonal part computed'
  endif
  call mp_barrier(world_comm)
  call initialize_exc(a_rot)

! apply the exchange term of the Hamiltonian
  if(.not.l_rpa) then
!     if(.not.l_truncated_coulomb) then
!        iuv = find_free_unit()
!        CALL diropn( iuv, 'vgq', npwx, exst )
!        CALL davcio(vg_q,npwx,iuv,1,-1)
!        close(iuv)
!     endif


     if(debug) then 
        write(stdout,*) 'vg_q read'
     call mp_barrier(world_comm)
     endif

     call exchange_exc(a_in,vstate,vstate_r,fc,a_exchange)

     if(debug) then 
        write(stdout,*) 'Exchange part computed'
     endif

     call mp_barrier(world_comm)

!    apply the direct term of the Hamiltonian (v part)
     
     if(.not.l_lf) then
   
        call initialize_exc(a_rot)
        a_rot%label=5
        a_rot%npw=npw
        a_rot%numb_v=num_nbndv(1)
        allocate(a_rot%a(a_rot%npw,a_rot%numb_v))

!       first rotate the excitonic wave function wave vector to use the wannier
!       wavefunctions

        if(debug) write(stdout,*)  'DEBUG1'
        Call urot_a(a_in,a_rot,0)
        if(debug) write(stdout,*)  'DEBUG2'
        if(.not.l_contraction) then
           call direct_v_exc(a_rot,fc,a_dirv)
        else
           call contract_v_apply(a_rot,fc,a_dirv)
        endif
        if(debug) write(stdout,*)  'DEBUG3'
        call pc_operator_exc(a_dirv,vstate,1)
        if(debug) write(stdout,*)  'DEBUG4'
!        and rotate back
        call urot_a(a_dirv,a_rot,1)
        a_dirv%a(1:a_dirv%npw,1:a_dirv%numb_v)=a_rot%a(1:a_rot%npw,1:a_rot%numb_v)

        if(debug) write(stdout,*)  'DEBUG5'
        call mp_barrier(world_comm)
!
        if(debug) then 
            write(stdout,*) 'After direct_v_exc'
        endif

! apply the direct term of the excitonic Hamiltonian (Wc part)

  
        if(.not.l_tdhf) then
           call urot_a(a_in,a_rot,0)

           if(debug) then 
              write(stdout,*) 'Before direct_W_exc'
           endif

           !if(.true.) then !DEBUG
           if(.not.l_contraction) then
              call direct_w_exc(a_rot,fc,a_dirw)
           else
              call contract_w_apply(a_rot,fc,a_dirw)
           endif
           call pc_operator_exc(a_dirw,vstate,1)

!           and rotate back
           call urot_a(a_dirw,a_rot,1)
           a_dirw%a(1:a_dirw%npw,1:a_dirw%numb_v)=a_rot%a(1:a_rot%npw,1:a_rot%numb_v)

           call mp_barrier(world_comm)
           if(debug) then 
              write(stdout,*) 'After direct_W_exc'
           endif




        endif ! .not.l_tdhf
     endif ! .not.l_lf
  endif  ! .not. l_rpa  

! sum up all the terms
! only valid for spin-singlet class of solutions (non-spin polarized case)
  a_out%a(1:a_out%npw,1:a_out%numb_v)=a_excdiago%a(1:a_out%npw,1:a_out%numb_v)&
                                      -a_dirv%a(1:a_out%npw,1:a_out%numb_v)&
                                     -a_dirw%a(1:a_out%npw,1:a_out%numb_v)&
                                     +2.d0*a_exchange%a(1:a_out%npw,1:a_out%numb_v)

! free memory
  call free_memory_exc_a(a_rot)
  call free_memory_exc_a(a_excdiago)
  call free_memory_exc_a(a_exchange)
  call free_memory_exc_a(a_dirv)
  call free_memory_exc_a(a_dirw)



  call stop_clock('exc_h_a')

if(debug) then
  call print_clock('exc_h_a')
  call print_clock('diago_exc')
  call print_clock('exchange_exc')
  call print_clock('direct_v_exc')
  call print_clock('direct_w_contract')
  call print_clock('wdirect_fftback')
  call print_clock('contract_w_dgemv')
  call print_clock('direct_v_contract')
  call print_clock('d_v_fft')
endif

 return
end subroutine exc_h_a





