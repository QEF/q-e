program bse_punch

use io_global, ONLY : stdout, ionode, ionode_id
use io_files,  ONLY : psfile, pseudo_dir,diropn
use io_files,  ONLY : prefix,tmp_dir,iunwfc
use mp_world, ONLY : mpime
use mp_pools, ONLY : kunit
USE wvfct,     ONLY : nbnd, et, npwx
USE gvecw,              ONLY : ecutwfc
USE gvecs,              ONLY : doublegrid
use pwcom
USE wavefunctions_module, ONLY : evc
use mp, ONLY: mp_bcast
USE mp_world, ONLY : world_comm
USE fft_base,             ONLY : dfftp
use scf, only : vrs, vltot, v, kedtau
USE fft_custom_gwl
use bse_basic_structures
use exciton
USE constants,   ONLY: RYTOEV
USE mp,          ONLY: mp_barrier
USE qpe_exc,         ONLY: qpc
use bse_wannier, ONLY: num_nbndv,&
           l_truncated_coulomb,&
           truncation_radius, &
           numw_prod,&
           dual_bse,&
           l_verbose, &
           lambda,eps,&
           l_cgrad,maxit,cg_nreset,lm_delta,n_eig,eps_eig, scissor,&
           l_plotexc,plotn_min,plotn_max,r_hole,l_plotaverage,&
           l_tspace,l_finite,r_pola,&
           spectra_e_min,spectra_e_max,spectra_nstep,spectra_broad,&
           l_restart,n_eig_start, nit_lcz,l_lanczos, l_restart_lcz, nlcz_restart,&
           l_tdhf,l_fullbse,l_lf,l_rpa, l_contraction, l_gtrick, qpe_imin, qpe_imax,&
           l_scissor,l_dielectric
implicit none
INTEGER, EXTERNAL :: find_free_unit

integer :: i, kunittmp, ios, is
CHARACTER(LEN=256), EXTERNAL :: trimcheck
CHARACTER(LEN=256) :: outdir
character(len=200) :: pp_file
logical ::  uspp_spsi, ascii, single_file, raw


type(v_state) :: vstate
type(v_state_r) :: vstate_r
type(c_state)   :: cstate
type(c_state)   :: wcstate
type(exc) :: a_exc
type(exc) :: b_exc
!type(exc) :: a_excdiago,a_exchange 
!type(exc):: a_dirv
!type(exc):: a_dirw
!type(exc):: a_rot
type(fft_cus) :: fc


logical exst
integer iuv

logical :: debug

real(kind=DP) :: sdeig



NAMELIST /inputbse/ prefix,num_nbndv,dual_bse,outdir,l_truncated_coulomb,&
                    truncation_radius, numw_prod, l_verbose,lambda,eps,&
                    l_cgrad,maxit,cg_nreset,lm_delta,n_eig,eps_eig,&
                    scissor,l_plotexc,plotn_min,plotn_max,r_hole,&
                    l_plotaverage,l_tspace,l_finite,r_pola,&
                    spectra_e_min,spectra_e_max,spectra_nstep,spectra_broad,&
                    l_restart,n_eig_start, nit_lcz,l_lanczos, l_restart_lcz, nlcz_restart,&
                    l_fullbse,l_tdhf,l_lf,l_rpa,l_contraction,l_gtrick, qpe_imin, qpe_imax,&
                    l_scissor,l_dielectric 

debug=.false.
            
call start_bse( )
call start_clock('bse_main')


!
!   set default values for variables in namelist
!
prefix='export'
 CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )	
IF ( TRIM( outdir ) == ' ' ) outdir = './'
  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.



num_nbndv(1:2) = 1
l_truncated_coulomb = .false.
truncation_radius = 10.d0
numw_prod=1
dual_bse=1.d0
l_verbose=.false.
lambda=0.001d0
eps=0.0001d0
maxit=400
cg_nreset=25
l_cgrad=.true.
lm_delta=0.3
n_eig=1
eps_eig=0.0000000001
scissor=0.d0
l_plotexc=.false.
plotn_min=0 
plotn_max=0 
r_hole(1:3)=0.d0
l_plotaverage=.false.
l_tspace=.false.
l_finite=.false.
r_pola(1:3)=1.d0
spectra_e_min=0.d0
spectra_e_max=10.d0
spectra_nstep=100
spectra_broad=0.01d0
l_restart=0
n_eig_start=0
nit_lcz=100
l_lanczos=.false.
l_restart_lcz=.false.
nlcz_restart=1
l_fullbse=.true.
l_tdhf=.false.
l_lf=.false.
l_rpa=.false.
l_contraction=.true.
l_gtrick=.true.
qpe_imin=1
qpe_imax=1
l_scissor=.true.
l_dielectric=.false.
!
!    Reading input file
!
IF ( ionode ) THEN
      !
      CALL input_from_file ( )
      !
      READ(5,inputbse,IOSTAT=ios)
      !
!      call read_namelists( 'PW4GWW' )
      !
      IF (ios /= 0) CALL errore ('pw4gww', 'reading inputbse namelist', ABS(ios) )
      scissor=scissor/RYTOEV
ENDIF


!-------------------------------------------------------------------------
! ... Broadcasting variables
!------------------------------------------------------------------------

   

  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id , world_comm)
  CALL mp_bcast( num_nbndv,     ionode_id , world_comm)
  CALL mp_bcast(l_truncated_coulomb, ionode_id, world_comm)
  CALL mp_bcast(truncation_radius, ionode_id, world_comm)
  call mp_bcast(numw_prod, ionode_id, world_comm)
  CALL mp_bcast(dual_bse, ionode_id, world_comm)
  CALL mp_bcast( pp_file, ionode_id , world_comm)
  CALL mp_bcast( uspp_spsi, ionode_id , world_comm)
  CALL mp_bcast( ascii, ionode_id , world_comm)
  CALL mp_bcast( single_file, ionode_id , world_comm)
  CALL mp_bcast( raw, ionode_id , world_comm)
  CALL mp_bcast( pseudo_dir, ionode_id , world_comm)
  CALL mp_bcast( psfile, ionode_id , world_comm)
  CALL mp_bcast( lambda, ionode_id , world_comm)
  CALL mp_bcast( eps, ionode_id , world_comm)
  CALL mp_bcast( maxit, ionode_id , world_comm)
  CALL mp_bcast( cg_nreset, ionode_id , world_comm)
  CALL mp_bcast( l_cgrad, ionode_id , world_comm)
  CALL mp_bcast( lm_delta, ionode_id , world_comm)
  CALL mp_bcast( n_eig, ionode_id , world_comm)
  CALL mp_bcast( eps_eig, ionode_id , world_comm)
  CALL mp_bcast( scissor, ionode_id , world_comm)
  CALL mp_bcast( l_plotexc, ionode_id , world_comm)
  CALL mp_bcast( plotn_min, ionode_id , world_comm) 
  CALL mp_bcast( plotn_max, ionode_id , world_comm) 
  CALL mp_bcast( r_hole, ionode_id, world_comm ) 
  CALL mp_bcast( l_plotaverage, ionode_id, world_comm ) 
  CALL mp_bcast( l_tspace, ionode_id, world_comm ) 
  CALL mp_bcast( l_finite, ionode_id, world_comm ) 
  CALL mp_bcast( r_pola, ionode_id, world_comm ) 
  CALL mp_bcast( spectra_e_min, ionode_id, world_comm ) 
  CALL mp_bcast( spectra_e_max, ionode_id, world_comm ) 
  CALL mp_bcast( spectra_nstep, ionode_id, world_comm ) 
  CALL mp_bcast( spectra_broad, ionode_id, world_comm ) 
  CALL mp_bcast( l_restart, ionode_id, world_comm ) 
  CALL mp_bcast( n_eig_start, ionode_id, world_comm ) 
  CALL mp_bcast( nit_lcz, ionode_id, world_comm ) 
  CALL mp_bcast( l_lanczos, ionode_id, world_comm ) 
  CALL mp_bcast( l_restart_lcz, ionode_id, world_comm ) 
  CALL mp_bcast( nlcz_restart, ionode_id, world_comm ) 
  CALL mp_bcast( l_fullbse, ionode_id, world_comm ) 
  CALL mp_bcast( l_tdhf, ionode_id, world_comm ) 
  CALL mp_bcast( l_lf, ionode_id, world_comm ) 
  CALL mp_bcast( l_rpa, ionode_id, world_comm ) 
  CALL mp_bcast( l_contraction, ionode_id, world_comm)
  CALL mp_bcast( l_gtrick, ionode_id, world_comm)
  CALL mp_bcast( l_scissor, ionode_id, world_comm)
  CALL mp_bcast( qpe_imin, ionode_id, world_comm)
  CALL mp_bcast( qpe_imax, ionode_id, world_comm)
  CALL mp_bcast( l_dielectric, ionode_id, world_comm)

  call read_file 
! after read_file everything is known

#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  call openfil_bse

  call read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)
  call summary()  
  call print_bseinfo()

  CALL hinit0()
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )

  if(l_verbose) write(stdout,*) 'To check, we print the KS eigenvalues:'
  FLUSH( stdout )
  !
  CALL print_ks_energies()

! inizialize dual grid once for all
  fc%dual_t=dual_bse
  fc%ecutt=ecutwfc
  call initialize_fft_custom(fc)

! read ks wavefunction, allocate and fill up the v_state object and c_state object
  if (allocated(evc)) deallocate (evc)
  call initialize_v_state(vstate)
  call make_v_state(num_nbndv,vstate)
  
  call initialize_c_state(cstate)
  call make_c_state(num_nbndv,cstate)

  call initialize_c_state(wcstate)
  call make_c_state(num_nbndv,wcstate)

! FFT the valence states vector into r-space (using dual_bse)
  call initialize_v_state_r(vstate_r)
  call v_wfng_to_wfnr(vstate,fc,vstate_r)
 
  if(.not.l_gtrick) then
     call v_wfng_to_wfnr(vstate,fc,vstate_r)
  endif

! if debug mode check polarizability basis orthonormality
  if(debug) then
     call check_basis(numw_prod,npw) 
  endif

! allocate once for all vg_q
  allocate(vg_q(npwx))
  if(.not.l_truncated_coulomb) then
     iuv = find_free_unit()
     CALL diropn( iuv, 'vgq', npwx, exst )
     CALL davcio(vg_q,npwx,iuv,1,-1)
     close(iuv)
  endif
     

  if(debug) write(*,*) 'vgq allocated'

! QP corrections read and used to prepare wcstate
  if(.not.l_scissor) then
     allocate(qpc(qpe_imax))
     call qpcorrections(wcstate)
  endif

  if(l_tspace) then
!    solve the BSE in transition space
     call tspace_diago(vstate,vstate_r,fc)
  else
!    solve the BSE with cg, or steepest descent procedure, compute the optical spectrum,
!    and the excitonic wfns 
     if(l_lanczos) then
        if(debug) write(*,*) 'Solve using Lanczos'
        if(l_gtrick) call v_wfng_to_wfnr(vstate,fc,vstate_r)
        call lanczos(vstate,vstate_r,cstate,wcstate,fc)
     else
        if(l_gtrick) call v_wfng_to_wfnr(vstate,fc,vstate_r)!still to be implemented
        call find_eig(vstate,vstate_r,cstate,wcstate,fc)
     endif
     call mp_barrier(world_comm)
  endif 

! 


! free memory
  call free_v_state_r(vstate_r)
  call free_v_state(vstate)
  call free_c_state(cstate)
  call free_c_state(wcstate)
  if(.not.l_scissor) deallocate(qpc)

  write(stdout,*) 'BSE COMPLETED'
  call stop_clock('bse_main')
  call print_clock('bse_main')
  call print_clock('fft')
  call print_clock('ffts')
  call print_clock('fftw')
  call print_clock('cft3t')
  call print_clock('davcio')
  call print_clock('make_v_state')
  call print_clock('make_c_state')
  call print_clock('v_wfng_to_wfnr')
  call print_clock('c_wfng_to_wfnr')
  call print_clock('c_times_exc')
  call print_clock('pc_operator_exc')
  call print_clock('sproduct_exc')
  call print_clock('normalize_exc')
  call print_clock('pout_operator_exc')
  call print_clock('fft_a_exc')
  call print_clock('fftback_a_exc')
  call print_clock('urot_a')
  CALL print_clock('cgsolve')
  call print_clock('conjgrad')
  call print_clock('linmin')
  call print_clock('diago_exc')
  call print_clock('direct_v_exc')
  call print_clock('direct_w_exc')
  call print_clock('direct_w_dgemv')
  call print_clock('dgemv1')
  call print_clock('dgemv2')
  call print_clock('dgemv3')
  call print_clock('dgemv4')
  call print_clock('direct_w_cft3t')
  call print_clock('wdirect_fftback')
  call print_clock('exchange_exc')
  call print_clock('direct_w_contract')
  call print_clock('direct_v_contract')
  call print_clock('dvpsi_e')
  call print_clock('exc_h_a')
  call print_clock('find_eig')
  call print_clock('h_h')
  call print_clock('lanczos')
  call print_clock('lanczos_iterations')
  call print_clock('lanczos_cf')
  call print_clock('plot_excwfn')
  call print_clock('print_spectrum')
  call print_clock('read_export')
  call print_clock('rotate_wannier_gamma_bse')
  call print_clock('sdescent')
  call print_clock('build_spectrum')
  call print_clock('absorption')
  call print_clock('amplitude_finite')
  call print_clock('amplitude')
  call print_clock('tspace_diago')
  call print_clock('build_exch')
  call print_clock('read_wannier_matrix')

  FLUSH( stdout )

  call stop_pp

  stop
end program bse_punch





