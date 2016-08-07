!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutine build a set of fake_conduction states
!then it find the optimal basis set for representing 
!the products with valence wannier functions


MODULE fake_cond_mod
 USE kinds, ONLY : DP
 IMPLICIT NONE
 SAVE

 INTEGER :: fcw_number!number of "producs of fake conduction with valence wannier" states for O matrix method
 COMPLEX(kind=DP), ALLOCATABLE, DIMENSION(:,:) :: fcw_state! fcw states for O matrix method
 REAL(kind=DP), ALLOCATABLE, DIMENSION(:,:) :: fcw_mat! "fcw matrix 


CONTAINS
  subroutine fake_conduction_wannier( cutoff, s_cutoff,ks_wfcs ,l_frac, ks_wfcs_diag,l_cond)
 !IT WORKS ONLY FOR NORMCONSERVING PSEUDOPOTENTIALS
  !the valence states in G space must be in evc
  ! Gamma point version

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et, wg
   USE gvecw,    ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_pools, ONLY : intra_pool_comm
   USE mp_world, ONLY: world_comm, mpime, nproc
   USE gvecs,              ONLY : nls, nlsm, doublegrid

   USE kinds, ONLY : DP
   USE io_files, ONLY : prefix, tmp_dir, diropn
   USE g_psi_mod,            ONLY : h_diag, s_diag
   USE noncollin_module,     ONLY : noncolin, npol
   USE becmod,           ONLY : becp
   USE uspp,                 ONLY : vkb, nkb, okvan
   USE klist,                ONLY : xk,igk_k
   USE fft_custom_gwl
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp
   USE lsda_mod,             ONLY : nspin
   

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

!  INTEGER,INTENT(out) :: fcw_number!number of "fake conduction" states for O matrix method
!  COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: fcw_state! "fake conduction" states for O matrix method
!  REAL(kind=DP),    POINTER, DIMENSION(:,:) :: fcw_mat! "fake conduction" matrix

  REAL(kind=DP), INTENT(in) :: cutoff!cutoff for planewaves
  REAL(kind=DP), INTENT(in) :: s_cutoff!cutoff for orthonormalization
  COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npwx,nbnd,nspin)!Kohn-Sham or Wannier wavefunctios
  LOGICAL, INTENT(in) :: l_frac!if true consider fractional occupancies
  COMPLEX(kind=DP), INTENT(in) :: ks_wfcs_diag(npwx,nbnd,nspin)!Kohn-Sham wavefunctios 
  LOGICAL, INTENT(in) :: l_cond!if true consider also conduction states for the construction of the polarizability basis


  COMPLEX(kind=DP), ALLOCATABLE :: state_fc(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: state_g(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: fcw_state_old(:,:) 
  COMPLEX(kind=DP), ALLOCATABLE :: h_state_fc(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:),evc_t(:,:,:),state_fc_t(:,:,:),state_g_t(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: fcw_state_n(:,:)
  REAL(kind=DP), ALLOCATABLE :: wv_real(:),state_real(:),wv_real_all(:,:),state_real_tmp(:)
  REAL(kind=DP), ALLOCATABLE :: state_real_tmp2(:),state_real2(:)
  REAL(kind=DP), ALLOCATABLE :: omat(:,:)
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:),tmp_mat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: omat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: hmat(:,:)
  REAL(kind=DP), ALLOCATABLE :: e_fake(:), vec_fake(:,:)
  REAL(kind=DP), ALLOCATABLE :: gap(:)
  REAL(kind=DP), ALLOCATABLE :: hmat_i(:,:),hmat_o(:,:), omat_i(:,:)
  REAL(kind=DP), ALLOCATABLE :: ovec(:)
  REAL(kind=DP), ALLOCATABLE :: g2kint(:)
  INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
  INTEGER, ALLOCATABLE :: isuppz(:)
  INTEGER, ALLOCATABLE :: iclustr(:)
  INTEGER, ALLOCATABLE :: igkt(:)

  REAL(kind=DP):: sca1,sca2

  LOGICAL :: l_test=.false.!if true test the completness of the basis
  LOGICAL :: l_dsyevr=.true.!if true uses dsyevr instead of dsyev
  LOGICAL :: l_diago_cg=.true.!if true uses diago_cg instead of  dsyevr ATTENZIONE
  LOGICAL :: exst
  LOGICAL :: l_dsygvx=.false.!if .true. uses serial dsygvx instead of parallel diago_cg_g
  LOGICAL :: l_gramsc=.true.!if true orthonormalization through gram-schimdt
  LOGICAL :: l_diago_para=.true.!if true uses parallel diago_cg
  LOGICAL :: l_fft_custom=.false.

  INTEGER :: ig,ip, ii, iv, jj, iw, ir, is
  INTEGER :: num_fc!number of fake conduction states
  INTEGER :: lwork,info,liwork
  INTEGER :: n_out
  INTEGER :: fcw_number_old
  INTEGER :: l_blk,nbegin,nend
  INTEGER :: max_state
  INTEGER :: iunfcw
  INTEGER :: nsize
  INTEGER :: nbegin_loc,nend_loc,nsize_loc
  INTEGER :: n_found_state
!variables for scalapack
  INTEGER :: num_fc_r,num_fc_c,num_fc_dimr,num_fc_dimc
  INTEGER :: m,nz,icrow,iccol,iproc,ilrow,ilcol
  INTEGER :: desc_a(9),desc_b(9),desc_c(9)
  INTEGER :: n_computed
  INTEGER :: num_built!number of states already built
  INTEGER :: num_out
  INTEGER ::  kilobytes
  INTEGER ::  kb_old, kb_new

  INTEGER, EXTERNAL :: indxg2p,indxg2l

  TYPE(optimal_options) :: options



  INTEGER :: bufferx,fcw_numberx,fcw_number_oldx, fcw_numberx_tmp


  LOGICAL :: l_restart0!if true restart is enabled
  INTEGER :: iunrestart0, iv_start,iunfsr
  REAL(kind=DP), ALLOCATABLE :: state_fc_r(:,:,:)
  INTEGER :: num_nbndv_max, num_fc_spin
  INTEGER :: num_fc_eff(2),num_fc_eff_max
 
  LOGICAL :: l_do_optimal
  INTEGER :: iun_oap

  TYPE(fft_cus) :: fc

!determine bufferx,fcw_numberx
  bufferx=num_nbndv(1)*300/4
  bufferx=max(1000,bufferx)
!ONLY FOR PROJECT ON JADE
  bufferx=5000
  fcw_numberx=bufferx
  fcw_number_oldx=bufferx
  fcw_numberx_tmp=bufferx

!generate fake conduction states
!!determine number of states

!generate custom in grid in case can be equal to norm-conserving grid

  fc%ecutt=ecutwfc
  fc%dual_t=dual_pb

  write(stdout,*) 'Call initialize_fft_custom'
  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory0', kilobytes
  FLUSH(stdout)

  call initialize_fft_custom(fc)

  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory0.0', kilobytes
  FLUSH(stdout)

  ! this is for compatibility

  allocate( igkt( fc%npwt ) )
  do ig=1,fc%npwt
     igkt(ig)=ig
  enddo

  !allocate( evc_g( fc%ngmt_g ) )

  !plane waves basis set

  !state_fc are first obtained on the ordering of the normconserving grid

  g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
       ( g(2,igk_k(1:npw,1)) )**2 + &
       ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2
  
  num_fc=0
  do ig=1,npw
     if(g2kin(ig) <= cutoff) num_fc=num_fc+1
  enddo
  call mp_sum(num_fc,world_comm)
  num_fc=(num_fc-1)*2+1

  if(.not.l_cond) then
     if(.not.l_frac) then
        num_fc_eff(1:2)=num_fc
        num_fc_eff_max=num_fc
     else
!NOT_TO_BE_INCLUDED_START
        num_fc_eff(1:2)=num_fc+num_nbndv(1:2)-num_nbndv_min(1:2)
        num_fc_eff_max=max(num_fc_eff(1),num_fc_eff(2))
!NOT_TO_BE_INCLUDED_END
     endif
  else
      if(.not.l_frac) then
         num_fc_eff(1:2)=num_fc+num_nbnds-num_nbndv(1:2)
         num_fc_eff_max=num_fc+num_nbnds-min(num_nbndv(1),num_nbndv(2))
     else
!NOT_TO_BE_INCLUDED_START
        num_fc_eff(1:2)=num_fc+num_nbndv(1:2)-num_nbndv_min(1:2)+num_nbnds-num_nbndv(1:2)
        num_fc_eff_max=max(num_fc_eff(1),num_fc_eff(2))
!NOT_TO_BE_INCLUDED_END
     endif
  endif
  allocate( state_fc( npw, num_fc_eff_max, nspin ) )

  state_fc(:,:,:)=(0.d0,0.d0)
  write(stdout,*) "Number of projected orthonormalized plane waves:", num_fc
  CALL memstat( kilobytes )
  if(l_verbose)  write(stdout,*) 'memory0.1', kilobytes, ' new kb = ', &
       &(SIZE( state_fc )*16 + SIZE( igkt )*4)/1024
  FLUSH(stdout)
  
  ii=0
  do ip=0,nproc-1
     if(mpime==ip) then
        do ig=gstart,npw
           if(g2kin(ig) <= cutoff) then
              ii=ii+1
              state_fc(ig,ii,1)=cmplx(dsqrt(0.5d0),0.d0)
              ii=ii+1
              state_fc(ig,ii,1)=cmplx(0.d0,dsqrt(0.5d0))
           endif
        enddo
        if(gstart==2) then
           ii=ii+1
           state_fc(1,ii,1)=(1.d0,0.d0)
        endif
     else
        ii=0
     endif
     call mp_sum(ii,world_comm)
  enddo

  if(ii/=num_fc) then 
     write(stdout,*) 'ERRORE FAKE CONDUCTION',ii
     FLUSH(stdout)
     stop
    return
  endif
  if(l_verbose)  write(stdout,*) 'FAKE1'
  FLUSH(stdout)
  
  if(nspin==2) state_fc(:,1:num_fc,2)=state_fc(:,1:num_fc,1)
  do is=1,nspin

!!project out of valence space
     do ii=1,num_fc
        evc(1:npw,1:num_nbndv(is))=ks_wfcs(1:npw,1:num_nbndv(is),is)!for calling pc_operator
        call pc_operator(state_fc(:,ii,is),is,l_cond)
      enddo
  enddo


!!add partially occupied states
  if(l_frac) then
!NOT_TO_BE_INCLUDED_START
     do is=1,nspin
        do ii=num_nbndv_min(is)+1,num_nbndv(is)
           state_fc(1:npw,num_fc+ii-num_nbndv_min(is),is)=ks_wfcs_diag(1:npw,ii,is)
        enddo
     enddo
!NOT_TO_BE_INCLUDED_END
  endif


!!add conduction states if required
  if(l_cond) then
     if(.not.l_frac) then
        do is=1,nspin
           do ii=num_nbndv(is)+1,num_nbnds
              state_fc(1:npw,num_fc+ii-num_nbndv(is),is)=ks_wfcs_diag(1:npw,ii,is)
           enddo
        enddo
     else
!NOT_TO_BE_INCLUDED_START
        do is=1,nspin
           do ii=num_nbndv(is)+1,num_nbnds
              state_fc(1:npw,num_fc+num_nbndv(is)-num_nbndv_min(is)+ii-num_nbndv(is),is)=ks_wfcs_diag(1:npw,ii,is)
           enddo
        enddo
!NOT_TO_BE_INCLUDED_END
     endif
  endif
!orthonormalize fake_conduction states 

 

     
!for the moment finds all the first fcw_fast_n eigenstates
     
     if(l_verbose)  write(stdout,*) 'CASE ORTHONORMALIZATION ONLY'
     FLUSH(stdout)
     
!if required orthonormalize the projected plane_waves or read from disk

     l_do_optimal=.false.
     inquire(file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', exist = exst)
     if(.not. exst) then
        l_do_optimal=.true.
     else
        iunrestart0 =  find_free_unit()
        open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='old')
        read(iunrestart0,*) iv_start
        close(iunrestart0)
        if(iv_start<1 ) l_do_optimal=.true.
     endif
           


     if(l_do_optimal) then
        if(l_verbose) write(stdout,*) 'Call optimal driver'
        FLUSH(stdout)
        options%l_complete=.true.
        options%idiago=0
        do is=1,nspin
           call optimal_driver(num_fc_eff(is),state_fc(1,1,is),npw,options,num_out, info)
        enddo

!read orthonormalized projected plane-waves from disk

     endif
        CALL memstat( kilobytes )
        if(l_verbose)  write(stdout,*) 'memory0.3', kilobytes
        FLUSH(stdout)
        
       
     
  !now state_fc are put on the ordering of the redueced grid, if required
  allocate(state_fc_t(fc%npwt,num_fc_eff_max,nspin))
 
 if(l_do_optimal) then
    if(fc%dual_t==4.d0) then
        do is=1,nspin
           state_fc_t(1:fc%npwt,1:num_fc_eff(is),is)=state_fc(1:fc%npwt,1:num_fc_eff(is),is)
        enddo
     else
        do is=1,nspin
           call reorderwfp_col (num_fc_eff(is),npw,fc%npwt,state_fc(1,1,is),state_fc_t(1,1,is), npw,fc%npwt, &
                & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
        enddo
     endif
    

     iun_oap = find_free_unit()
     CALL diropn( iun_oap, 'oap', fc%npwt*2, exst )
     do ii=1,num_fc_eff(1)
        CALL davcio( state_fc_t(:,ii,1), 2*fc%npwt, iun_oap, ii, 1 )
     enddo
     close(iun_oap)

  else
     if(l_verbose) write(stdout,*) 'Read OAP from disk'
     FLUSH(stdout)
     iun_oap = find_free_unit()
     CALL diropn( iun_oap, 'oap', fc%npwt*2, exst )
     do ii=1,num_fc_eff(1)
        CALL davcio( state_fc_t(:,ii,1), 2*fc%npwt, iun_oap, ii, -1 )
     enddo
     close(iun_oap)


  endif
 
  deallocate(state_fc)
  if(l_iter_algorithm) then
      allocate(state_fc_r(fc%nrxxt,num_fc_eff_max,nspin))
      do is=1,nspin
         do ii=1,num_fc_eff(is),2
            psic(:)=(0.d0,0.d0)
            if(ii==num_fc_eff(is)) then
               psic(fc%nlt(1:fc%npwt))  = state_fc_t(1:fc%npwt,ii,is)
               psic(fc%nltm(1:fc%npwt)) = CONJG( state_fc_t(1:fc%npwt,ii,is) )
            else
               psic(fc%nlt(1:fc%npwt))=state_fc_t(1:fc%npwt,ii,is)+(0.d0,1.d0)*state_fc_t(1:fc%npwt,ii+1,is)
               psic(fc%nltm(1:fc%npwt)) = CONJG( state_fc_t(1:fc%npwt,ii,is) )+(0.d0,1.d0)*CONJG( state_fc_t(1:fc%npwt,ii+1,is) )
            endif
            

            CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
           

            state_fc_r(1:fc%nrxxt,ii,is)= DBLE(psic(1:fc%nrxxt))
            if(ii/=num_fc_eff(is)) state_fc_r(1:fc%nrxxt,ii+1,is)= DIMAG(psic(1:fc%nrxxt))
         enddo
      enddo
      deallocate(state_fc_t) 
   endif

  CALL memstat( kilobytes )
  if(l_verbose)  write(stdout,*) 'memory0.4', kilobytes, ' new kb = ', (SIZE( state_fc_t ))/64
  FLUSH(stdout)

!set maximum number of valence states for both spin channels
  if(nspin==1) then
     num_nbndv_max=num_nbndv(1)
  else
     num_nbndv_max=max(num_nbndv(1),num_nbndv(2))
  endif

!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbndv_max,nspin))
  if(fc%dual_t==4.d0) then
     evc_t(1:fc%npwt,1:num_nbndv_max,1:nspin)=ks_wfcs(1:fc%npwt,1:num_nbndv_max,1:nspin)
  else
     do is=1,nspin
        call reorderwfp_col(num_nbndv(is),npw,fc%npwt,ks_wfcs(1,1,is),evc_t(1,1,is), npw,fc%npwt, &
             & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
     !   do iv=1,num_nbndv(is)
     !      call mergewf(ks_wfcs(:,iv,is),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
     !      call splitwf(evc_t(:,iv,is),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     !   enddo
     enddo
  endif

  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory0.5', kilobytes, ' new kb = ', (SIZE( evc_t ))/64
  FLUSH(stdout)

!cycle on v 
!! product in real space with wannier
!! orthonormalize and take N most important
!! gram-schmidt like 

!calculate D matrix 

!  l_blk= (num_fc)/nproc
!  if(l_blk*nproc < (num_fc)) l_blk = l_blk+1
!  nbegin=mpime*l_blk+1
!  nend=nbegin+l_blk-1
!  if(nend > num_fc) nend=num_fc
!  nsize=nend-nbegin+1

!check for restart
 if(ionode) then

     inquire(file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', exist = exst)
     if(.not. exst) then
        iv_start=1
     else
        iunrestart0 =  find_free_unit()
        open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='old')
        read(iunrestart0,*) iv_start
        read(iunrestart0,*) fcw_number
        read(iunrestart0,*) fcw_numberx
        close(iunrestart0)
        if(iv_start<1 ) then
           iv_start=1
        else
           iv_start=iv_start+1
        endif
     endif
  endif
  call mp_bcast(iv_start,ionode_id,world_comm)

  if(iv_start/=1) then
     call mp_bcast(fcw_number,ionode_id,world_comm)
     call mp_bcast(fcw_numberx,ionode_id,world_comm)
     fcw_number_oldx=fcw_numberx
     fcw_numberx_tmp=fcw_numberx
     fcw_number_old=fcw_number
     allocate(fcw_state(fc%npwt,fcw_numberx))
     allocate(fcw_state_old(fc%npwt,fcw_numberx))
 
 !read them from file
     iunfsr = find_free_unit()
 
     CALL diropn( iunfsr, 'fsr', fc%npwt*2, exst )
     do ii=1,fcw_number
        CALL davcio( fcw_state(1,ii), 2*fc%npwt, iunfsr, ii, -1 )
        fcw_state_old(:,ii)=fcw_state(:,ii)
     enddo
 
     close(iunfsr)

  endif

  allocate (wv_real(fc%nrxxt),state_real(fc%nrxxt),state_real2(fc%nrxxt),state_g(fc%npwt,num_fc_eff_max*nspin))


  
  FIRST_LOOP: do iv=iv_start,num_nbndv_max
     call mp_barrier( world_comm )
     write(stdout,*) 'FK state:', iv,fc%nrxxt,fc%npwt,num_fc
     CALL memstat( kilobytes )
     if(l_verbose)  write(stdout,*) 'memory1', kilobytes
     FLUSH(stdout)
     

     call mp_barrier( world_comm )
     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory2', kilobytes, ' new kb = ', &
                     (SIZE(wv_real)+SIZE(state_real)+SIZE(state_real2)+SIZE(state_g)*2)/128
     FLUSH(stdout)

     num_fc_spin=0
     do is=1,nspin

        if(iv<= num_nbndv(is)) then
           psic(:)=(0.d0,0.d0)
           psic(fc%nlt(igkt(1:fc%npwt)))  = evc_t(1:fc%npwt,iv,is)
           psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( evc_t(1:fc%npwt,iv,is) )

           call mp_barrier( world_comm )
           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory3', kilobytes
           FLUSH(stdout)

           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
           call mp_barrier( world_comm )
           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory4', kilobytes
           FLUSH(stdout)

     
           wv_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))
           if(l_verbose) then
              if(fc%gstart_t==2) write(stdout,*) 'FAKE modulus valence', iv, evc_t(1,iv,is)
           endif
!loop on fake conduction states
           if(l_verbose) write(stdout,*) 'Start FFTs part'
           FLUSH(stdout)
           do ii=1,num_fc_eff(is),2
!fourier transform each state to real space
              if(.not.l_iter_algorithm) then
                 psic(:)=(0.d0,0.d0)
                 if(ii==num_fc_eff(is)) then
                    psic(fc%nlt(igkt(1:fc%npwt)))  = state_fc_t(1:fc%npwt,ii,is)
                    psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( state_fc_t(1:fc%npwt,ii,is) )
                 else
                    psic(fc%nlt(igkt(1:fc%npwt)))=state_fc_t(1:fc%npwt,ii,is)+(0.d0,1.d0)*state_fc_t(1:fc%npwt,ii+1,is)
                    psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( state_fc_t(1:fc%npwt,ii,is) )+&
                         &(0.d0,1.d0)*CONJG( state_fc_t(1:fc%npwt,ii+1,is) )
                 endif
                 CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
                 state_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))
                 if(ii<num_fc_eff(is)) state_real2(1:fc%nrxxt)=dimag(psic(1:fc%nrxxt))
!form product in real space
                 state_real(1:fc%nrxxt)=state_real(1:fc%nrxxt)*wv_real(1:fc%nrxxt)
                 if(ii<num_fc_eff(is)) state_real2(1:fc%nrxxt)=state_real2(1:fc%nrxxt)*wv_real(1:fc%nrxxt)
!back to G space
                 if(ii==num_fc_eff(is)) then
                    psic(1:fc%nrxxt)=dcmplx(state_real(1:fc%nrxxt),0.d0)
                 else
                    psic(1:fc%nrxxt)=dcmplx(state_real(1:fc%nrxxt),state_real2(1:fc%nrxxt))
                 endif
              else
                 state_real(1:fc%nrxxt)=state_fc_r(1:fc%nrxxt,ii,is)*wv_real(1:fc%nrxxt)
                 if(ii==num_fc_eff(is)) then
                    psic(1:fc%nrxxt)=dcmplx(state_real(1:fc%nrxxt),0.d0)
                 else
                    state_real2(1:fc%nrxxt)=state_fc_r(1:fc%nrxxt,ii+1,is)*wv_real(1:fc%nrxxt)
                    psic(1:fc%nrxxt)=dcmplx(state_real(1:fc%nrxxt),state_real2(1:fc%nrxxt))
                 endif
              endif
              CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
              if(ii==num_fc_eff(is)) then
                 state_g(1:fc%npwt, ii+num_fc_spin) = psic(fc%nlt(igkt(1:fc%npwt)))
                 if(fc%gstart_t==2) state_g(1,ii+num_fc_spin)=(0.d0,0.d0)
              else
                 state_g(1:fc%npwt, ii+num_fc_spin)= 0.5d0*(psic(fc%nlt(igkt(1:fc%npwt)))+conjg( psic(fc%nltm(igkt(1:fc%npwt)))))
                 state_g(1:fc%npwt, ii+1+num_fc_spin)= (0.d0,-0.5d0)*(psic(fc%nlt(igkt(1:fc%npwt))) - &
                      &conjg(psic(fc%nltm(igkt(1:fc%npwt)))))
                 if(fc%gstart_t==2) state_g(1,ii+num_fc_spin)=(0.d0,0.d0)
                 if(fc%gstart_t==2) state_g(1,ii+1+num_fc_spin)=(0.d0,0.d0)
              endif
       
       
           enddo
           num_fc_spin=num_fc_spin+num_fc_eff(is)
        endif
     enddo!on spin

     if(l_verbose) write(stdout,*) 'End FFTs part'
     call mp_barrier( world_comm )
     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory5', kilobytes
     FLUSH(stdout)

!if not first valence wannier project the products out of partial orthonormal basis

     if(l_verbose) write(stdout,*) 'Start Projection part'
     FLUSH(stdout)

     if(iv/=1) then
        ! build overlap matrix
        
        if(iv==2 .or. (iv_start/=1.and. iv==iv_start)) then
           allocate(tmp_mat2(fcw_numberx_tmp,num_fc_eff_max*nspin))
        else
           if(fcw_number>fcw_numberx_tmp) then
              deallocate(tmp_mat2)
              fcw_numberx_tmp=fcw_numberx_tmp+bufferx
              allocate(tmp_mat2(fcw_numberx_tmp,num_fc_eff_max*nspin))
              if(l_verbose) write(stdout,*) 'Updated dimension of tmp_mat2', fcw_numberx_tmp
           endif
        endif
        
        call dgemm('T','N',fcw_number,num_fc_spin,2*fc%npwt,2.d0,fcw_state,2*fc%npwt,&
             &state_g,2*fc%npwt,0.d0,tmp_mat2,fcw_numberx_tmp)
        if(fc%gstart_t==2) then
           do ii=1,num_fc_spin
              do jj=1,fcw_number
                 tmp_mat2(jj,ii)=tmp_mat2(jj,ii)-dble(conjg(fcw_state(1,jj))*state_g(1,ii))
              enddo
           enddo
        endif
        do ii=1,num_fc_spin
           call mp_sum(tmp_mat2(1:fcw_number,ii),world_comm)
        enddo
        !call mp_sum(tmp_mat2)

        call dgemm('N','N',2*fc%npwt, num_fc_spin,fcw_number,-1.d0,fcw_state,2*fc%npwt,tmp_mat2,&
        &fcw_numberx_tmp,1.d0,state_g,2*fc%npwt)
        
        if(iv==num_nbndv_max) deallocate(tmp_mat2)
     endif

     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory6', kilobytes
     if(l_verbose) write(stdout,*) 'End Projection part'
     FLUSH(stdout)

!calculate overlap matrix

     if(l_verbose) write(stdout,*) 'FK2'!ATTENZIONE
     FLUSH(stdout)

     max_state=max(300,num_fc/20)
     if(max_state > num_fc) max_state=num_fc/2

     l_blk= (num_fc_spin)/nproc
     if(l_blk*nproc < (num_fc_spin)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > num_fc_spin) nend=num_fc
     nsize=nend-nbegin+1



     if( .not. l_iter_algorithm ) then

        if(.not.l_diago_para) then
           allocate(omat(num_fc_spin,num_fc_spin))
           if(l_dsyevr) then
              allocate(omat2(num_fc_spin,num_fc_spin))
           else if(l_diago_cg) then
              allocate(omat2(num_fc_spin,max_state))
           endif
        else
           allocate(omat(num_fc_spin,l_blk))
           allocate(omat2(num_fc_spin,max_state))
        endif

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.1', kilobytes, ' new kb = ', (SIZE(omat)+SIZE(omat2))/128
        FLUSH(stdout)

        if(.not.l_diago_para) then
           call dgemm('T','N',num_fc_spin,num_fc_spin,2*fc%npwt,2.d0,state_g,2*fc%npwt,state_g,2*fc%npwt,0.d0,omat,num_fc_spin)
           if(fc%gstart_t==2) then
              do ii=1,num_fc_spin
                 do jj=1,num_fc_spin
                    omat(jj,ii)=omat(jj,ii)-dble(conjg(state_g(1,jj))*state_g(1,ii))
                 enddo
              enddo
           endif
           call mp_sum(omat,world_comm)
        else
           allocate(tmp_mat(num_fc_spin,l_blk))
           do ip=0,nproc-1
              nbegin_loc=ip*l_blk+1
              nend_loc=nbegin_loc+l_blk-1
              if(nend_loc > num_fc_spin) nend_loc=num_fc_spin
              nsize_loc=nend_loc-nbegin_loc+1
              if(nsize_loc >0) then
               
                 call dgemm('T','N',num_fc_spin,nsize_loc,2*fc%npwt,2.d0,state_g,2*fc%npwt,&
       &state_g(1,nbegin_loc),2*fc%npwt,0.d0,tmp_mat,num_fc_spin)
                 if(fc%gstart_t==2) then
                    do ii=nbegin_loc,nend_loc
                       do jj=1,num_fc_spin
                          tmp_mat(jj,ii-nbegin_loc+1)=tmp_mat(jj,ii-nbegin_loc+1)-dble(conjg(state_g(1,jj))*state_g(1,ii))
                       enddo
                    enddo
                 endif
                 do ii=1,nsize_loc
                    call mp_sum(tmp_mat(:,ii),world_comm)
                 enddo
                 if(ip==mpime) then
                    omat(:,1:nsize_loc)=tmp_mat(:,1:nsize_loc)
                 endif
              endif
              
           enddo
           deallocate(tmp_mat)
        endif

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.2', kilobytes
        FLUSH(stdout)

!solve eigenvalues problem
        allocate(eigen(num_fc_spin))
        if(.not.l_diago_cg) then
           if(ionode) then
              if(.not.l_dsyevr) then 
                 allocate(work(1))
                 call  DSYEV( 'V', 'U', num_fc_spin, omat, num_fc_spin, eigen, work, -1, info )
                 lwork=work(1)
                 deallocate(work)
                 allocate(work(lwork))
                 call  DSYEV( 'V', 'U', num_fc_spin, omat, num_fc_spin, eigen, work, lwork, info )
                 deallocate(work)
                 if(info/=0) then
                    write(stdout,*) 'ROUTINE fake_conduction_wannier, INFO:', info
                    stop
                 endif
              else
                 allocate(isuppz(2*num_fc_spin))
                 allocate(work(1),iwork(1))
                 call DSYEVR('V','V','U',num_fc_spin,omat,num_fc_spin,s_cutoff,1.d6,1,1,0.001d0*s_cutoff,&
  &n_out,eigen,omat2,num_fc_spin,isuppz,work,-1,iwork,-1,info)
                 lwork=work(1)
                 liwork=iwork(1)
                 deallocate(work,iwork)
                 allocate(work(lwork))
                 allocate(iwork(liwork))
                 call DSYEVR('V','V','U',num_fc_spin,omat,num_fc_spin,s_cutoff,10000d0,1,1,&
        &0.001d0*s_cutoff,n_out,eigen,omat2,num_fc_spin,isuppz,work,lwork,iwork,liwork,info)
                 if(info/=0) then
                    write(stdout,*) 'ROUTINE fake_conduction_wannier, INFO:', info
                    stop
                 endif
                 deallocate(work,iwork)
                 deallocate(isuppz)
              endif
           else
              omat(:,:)=0.d0
              if(l_dsyevr) then
                 omat2(:,:)=0.d0
                 n_out=0
              endif
              eigen(:)=0.d0
           endif
           if(l_verbose) write(stdout,*) 'FK3'!ATTENZIONE
           FLUSH(stdout)
         
           if(l_dsyevr) then
              call mp_sum(n_out,world_comm)
              do iw=1,n_out
                 call mp_sum(omat2(:,iw),world_comm)
              enddo
           else
              call mp_sum(omat,world_comm)
           endif
      !call mp_bcast(eigen(:), ionode_id,world_comm)
           call mp_sum(eigen(:),world_comm)
        else
           if(l_verbose) write(stdout,*) 'Before diago_cg',max_state
           FLUSH(stdout)
           if(l_diago_para) then
              call diago_cg(num_fc_spin,omat,1000,max_state,eigen,omat2,s_cutoff,0.000001d0*s_cutoff,n_out,.true.)
           else
              call diago_cg(num_fc_spin,omat,1000,max_state,eigen,omat2,s_cutoff,0.000001d0*s_cutoff,n_out,.false.)
           endif
           if(l_verbose) write(stdout,*) 'After diago_cg'
           FLUSH(stdout)
         
        endif

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.3', kilobytes, ' new kb = ', SIZE(eigen)/128
        FLUSH(stdout)

        !if first valence state construct first basis

        if(iv==1) then

           !construct orthonormal basis set
           !  state_out(:,:)=(0.d0,0.d0)

           if(.not.(l_dsyevr.or.l_diago_cg)) then
              n_out=0
              do ii=1,num_fc_spin
                 if(l_verbose) write(stdout,*) 'FK eigen:', eigen(ii)
                 if(eigen(ii) >= s_cutoff) then
                    n_out=n_out+1
                 endif
              enddo
           else
              do ii=1,n_out
                 if(l_verbose) write(stdout,*) 'FK eigen:', eigen(ii)
              enddo
           endif

           if(l_verbose) write(stdout,*) 'FK orthonormal states:', n_out, num_fc_spin
           FLUSH(stdout)

           
           if(.not.(l_dsyevr.or.l_diago_cg)) then
              do ii=num_fc_spin-n_out+1,num_fc_spin
                 omat(1:num_fc_spin,ii)=omat(1:num_fc_spin,ii)/dsqrt(eigen(ii))
              enddo
           else
              do ii=1,n_out
                 omat2(1:num_fc_spin,ii)=omat2(1:num_fc_spin,ii)/dsqrt(abs(eigen(ii)))
              enddo
           endif
           
           fcw_number=n_out
           if(fcw_number<fcw_numberx) then
              allocate(fcw_state(fc%npwt,fcw_numberx))
           else
              fcw_numberx=fcw_numberx+bufferx
              allocate(fcw_state(fc%npwt,fcw_numberx))
           endif

           
           if(.not.(l_dsyevr.or.l_diago_cg)) then

              call dgemm('N','N',2*fc%npwt,n_out,num_fc_spin,1.d0,state_g,2*fc%npwt,&
       &omat(1,num_fc_spin-n_out+1),num_fc_spin,0.d0,fcw_state,2*fc%npwt)
           else

              call dgemm('N','N',2*fc%npwt,n_out,num_fc_spin,1.d0,state_g,2*fc%npwt,omat2(1,1),num_fc_spin,0.d0,fcw_state,2*fc%npwt)
           endif
           if(l_verbose) write(stdout,*) 'FK4'!ATTENZIONE
           FLUSH(stdout)

           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory6.3.0', kilobytes, ' new kb = ', SIZE( fcw_state ) / 64
           FLUSH(stdout)

 !write restart on file
           iunfsr = find_free_unit()
           CALL diropn( iunfsr, 'fsr', fc%npwt*2, exst )
           do ii=1,fcw_number
              CALL davcio( fcw_state(1,ii), 2*fc%npwt, iunfsr, ii, 1 )
           enddo
           close(iunfsr)

        else


           if(.not.(l_dsyevr.or.l_diago_cg)) then
              n_out=0
              do ii=1,num_fc_spin
                 if(eigen(ii) >= s_cutoff) then
                    n_out=n_out+1
                 endif
              enddo
           endif

           if(l_verbose) write(stdout,*) 'FK orthonormal states:', n_out
           FLUSH(stdout)

         
           if(.not.(l_dsyevr.or.l_diago_cg)) then
              do ii=num_fc_spin-n_out+1,num_fc_spin
                 omat(1:num_fc_spin,ii)=omat(1:num_fc_spin,ii)/dsqrt(eigen(ii))
              enddo
           else
              do ii=1,n_out
                 omat2(1:num_fc_spin,ii)=omat2(1:num_fc_spin,ii)/dsqrt(abs(eigen(ii)))
              enddo
           endif

           fcw_number=n_out+fcw_number_old

           kb_old = SIZE( fcw_state ) 

           if(fcw_number>fcw_numberx) then
              fcw_numberx=fcw_numberx+bufferx
              deallocate(fcw_state)
              allocate(fcw_state(fc%npwt,fcw_numberx))
           endif

           kb_new = SIZE( fcw_state )

           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory6.3.1', kilobytes, ' new kb = ', ( kb_new - kb_old )/64
           FLUSH(stdout)

           if(.not.(l_dsyevr.or.l_diago_cg)) then
              call dgemm('N','N',2*fc%npwt,n_out,num_fc_spin,1.d0,state_g,2*fc%npwt,&
                    &omat(1,num_fc_spin-n_out+1),num_fc_spin,0.d0,fcw_state,2*fc%npwt)
           else
              call dgemm('N','N',2*fc%npwt,n_out,num_fc_spin,1.d0,state_g,2*fc%npwt,omat2(1,1),num_fc_spin,0.d0,fcw_state,2*fc%npwt)
           endif
           fcw_state(:,n_out+1:fcw_number)=fcw_state_old(:,1:fcw_number_old)

           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory6.3.2', kilobytes
           FLUSH(stdout)

 !write restart on file
           iunfsr = find_free_unit()
           CALL diropn( iunfsr, 'fsr', fc%npwt*2, exst )
           do ii=1,n_out
              CALL davcio( fcw_state(1,ii), 2*fc%npwt, iunfsr, fcw_number_old+ii, 1 )
           enddo
           close(iunfsr)
           
        endif


        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.4', kilobytes
        FLUSH(stdout)


        ! if iv is not the last save a copy of the basis

        if(iv/=num_nbndv_max) then
           if(l_verbose) write(stdout,*) 'FK5'!ATTENZIONE
           FLUSH(stdout)
           
           kb_old = 0
           if( iv/=1 .and. allocated( fcw_state_old ) ) then 
              kb_old = kb_old + SIZE( fcw_state_old )

              if(fcw_number > fcw_number_oldx)then
                 fcw_number_oldx=fcw_number_oldx+bufferx
                 deallocate(fcw_state_old)
                 allocate(fcw_state_old(fc%npwt,fcw_number_oldx))
              endif
           else
              allocate(fcw_state_old(fc%npwt,fcw_number_oldx))
           endif

           fcw_state_old(1:fc%npwt,1:fcw_number_oldx)=fcw_state(1:fc%npwt,1:fcw_number_oldx)
           kb_new = SIZE(fcw_state_old)

           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory6.4.1', kilobytes, ' new kb = ', (kb_new-kb_old)/64
           FLUSH(stdout)

           fcw_number_old=fcw_number

        endif

        if(l_verbose) write(stdout,*) 'FK6'!ATTENZIONE
        FLUSH(stdout)

        kb_old = SIZE( omat ) + SIZE( eigen )
        deallocate( omat, eigen )
        if( allocated( omat2 ) ) then
           kb_old = kb_old + SIZE( omat2 )
           deallocate(omat2)
        endif
        kb_old = kb_old + SIZE( wv_real ) + SIZE( state_real ) + SIZE( state_real2 ) + 2*SIZE( state_g )
        !deallocate( wv_real, state_real, state_real2, state_g )
        
        if(l_verbose) write(stdout,*) 'FK7'!ATTENZIONE
        FLUSH(stdout)

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.5', kilobytes, ' old kb = ', kb_old / 128
        FLUSH(stdout)
    
     else  ! -------------------------iter algorithm 

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.6', kilobytes
        FLUSH(stdout)

        !uses iterative algorithm
        !allocate max number of new states

        !deallocate(wv_real,state_real,state_real2)
        !gram shimdt like
        allocate(ovec(num_fc_spin))
        num_built=0


        do ii=1,num_fc_spin
           
           if(num_built>0) then
              call dgemm('T','N',num_built,1,2*fc%npwt,2.d0,state_g,2*fc%npwt,state_g(1,ii),2*fc%npwt,0.d0,ovec,num_fc_spin)
              if(fc%gstart_t==2) then
                 do jj=1,num_built
                    ovec(jj)=ovec(jj) -dble(conjg(state_g(1,jj))*state_g(1,ii))
                 enddo
              endif
              call mp_sum(ovec(1:num_built),world_comm)
              call dgemm('T','N',1,1,2*fc%npwt,2.d0,state_g(1,ii),2*fc%npwt,state_g(1,ii),2*fc%npwt,0.d0,sca2,1)
              if(fc%gstart_t==2) sca2=sca2-dble(conjg(state_g(1,ii))*state_g(1,ii))
              call mp_sum(sca2,world_comm)
              sca1=0.d0
              do jj=1,num_built
                 sca1=sca1+ovec(jj)**2.d0
              enddo
              if(abs(sca2-sca1) >= s_cutoff) then
                 if(num_built+1 /= ii) state_g(1:fc%npwt,num_built+1)=state_g(1:fc%npwt,ii)
                 call dgemm('N','N',2*fc%npwt,1,num_built,-1.d0,state_g,2*fc%npwt,&
                      &ovec,num_fc_spin,1.d0,state_g(1,num_built+1),2*fc%npwt)
                 num_built=num_built+1
                 call dgemm('T','N',1,1,2*fc%npwt,2.d0,state_g(1,num_built),&
                      &2*fc%npwt,state_g(1,num_built),2*fc%npwt,0.d0,ovec(num_built),1)
                 if(fc%gstart_t==2) ovec(num_built)=ovec(num_built)-dble(conjg(state_g(1,num_built))*state_g(1,num_built))
                 call mp_sum(ovec(num_built),world_comm)
                 ovec(num_built)=1.d0/dsqrt(ovec(num_built))
                 state_g(1:fc%npwt,num_built)=state_g(1:fc%npwt,num_built)*ovec(num_built)

              endif
           else
              call dgemm('T','N',1,1,2*fc%npwt,2.d0,state_g(1,ii),&
                   &2*fc%npwt,state_g(1,ii),2*fc%npwt,0.d0,ovec(1),1)
              if(fc%gstart_t==2) ovec(1)=ovec(1)-dble(conjg(state_g(1,ii))*state_g(1,ii))
              call mp_sum(ovec(1),world_comm)
              if(ovec(1) >= s_cutoff) then
                 num_built=1
                 ovec(num_built)=1.d0/dsqrt(ovec(num_built))
                 state_g(1:fc%npwt,num_built)=state_g(1:fc%npwt,ii)*ovec(num_built)                
              endif
           endif
        enddo

        deallocate(ovec)   
        
        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory6.7', kilobytes
        write(stdout,*) 'FK GS', num_built
        FLUSH(stdout)


        !opportune basis
        
        if( iv == 1 ) then
           fcw_number=num_built
           allocate( fcw_state( fc%npwt, fcw_numberx ) )
           allocate( fcw_state_old( fc%npwt, fcw_number_oldx ) )
           if(num_built>0) then
              fcw_state(1:fc%npwt,1:num_built) = state_g(1:fc%npwt,1:num_built)
              iunfsr = find_free_unit()
              CALL diropn( iunfsr, 'fsr', 2*fc%npwt, exst )
              do ii=1,num_built
                 CALL davcio( fcw_state(1,ii), 2*fc%npwt, iunfsr, ii, 1 )
              enddo
              close(iunfsr)
           endif
        else
           if(fcw_number+num_built>fcw_number_oldx) then
              fcw_number_oldx=fcw_number_oldx+bufferx
              deallocate(fcw_state_old)
              allocate( fcw_state_old( fc%npwt, fcw_number_oldx ) )
           endif
           
           fcw_state_old(1:fc%npwt,1:fcw_number) = fcw_state(1:fc%npwt,1:fcw_number)
           if(fcw_number+num_built>fcw_numberx) then
              fcw_numberx=fcw_numberx+bufferx
              deallocate(fcw_state)
              allocate(fcw_state(fc%npwt,fcw_numberx))
           endif
            
           fcw_state(1:fc%npwt,1:fcw_number)=fcw_state_old(1:fc%npwt,1:fcw_number)
           if(num_built> 0) then
              fcw_state(1:fc%npwt,fcw_number+1:fcw_number+num_built)=state_g(1:fc%npwt,1:num_built)
              CALL diropn( iunfsr, 'fsr', 2*fc%npwt, exst )
              do ii=1,num_built
                 CALL davcio( fcw_state(1,ii+fcw_number), 2*fc%npwt, iunfsr, ii+fcw_number, 1 )
              enddo
              close(iunfsr)
           endif
           fcw_number=fcw_number+num_built

        endif
     end if

        

     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory6.8', kilobytes
     FLUSH(stdout)
     iunrestart0 =  find_free_unit()
     open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='unknown')
     write(iunrestart0,*) iv
     write(iunrestart0,*) fcw_number
     write(iunrestart0,*) fcw_numberx
     close(iunrestart0)

  end do FIRST_LOOP
 
  deallocate( wv_real, state_real, state_real2, state_g )
 
  if(l_verbose) write(stdout,*) 'FK8'
  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory7', kilobytes
  FLUSH(stdout)
  
  if(num_nbndv(1)/=1 ) deallocate(fcw_state_old)
 
! calculate D matrix distributed among processors

  if(fcw_number < nproc) then
     write(stdout,*) 'too many processors'
     stop
  endif

  l_blk= (fcw_number)/nproc
  if(l_blk*nproc < (fcw_number)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > fcw_number) nend=fcw_number
  nsize=nend-nbegin+1

  allocate(fcw_mat(fcw_number,l_blk))
  fcw_mat(:,:)=0.d0

  if(l_verbose) write(stdout,*) 'FK9'
  FLUSH(stdout)

  allocate(wv_real_all(fc%nrxxt,num_nbndv_max))

  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory8', kilobytes
  
  allocate(state_real(fc%nrxxt),state_real_tmp(fc%nrxxt),state_real_tmp2(fc%nrxxt))
  
  allocate(state_g(fc%npwt,num_nbndv_max))
  
  allocate(tmp_mat(fcw_number,num_nbndv_max))

  write(stdout,*) 'Calculate FK matrix'
  FLUSH(stdout)

  do is=1,nspin
     do iv=1,num_nbndv(is)
      
     
        psic(:)=(0.d0,0.d0)
        psic(fc%nlt(igkt(1:fc%npwt)))  = evc_t(1:fc%npwt,iv,is)
        psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( evc_t(1:fc%npwt,iv,is) )
           
        CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
        wv_real_all(1:fc%nrxxt,iv)= DBLE(psic(1:fc%nrxxt))
!check for modulus
        sca1=0.d0
        do ir=1,fc%nrxxt
           sca1=sca1+wv_real_all(ir,iv)**2.d0
        enddo
        call mp_sum(sca1,world_comm)
        if(l_verbose) write(stdout,*) 'Modulus:',fc%nrxxt,fc%nr1t*fc%nr2t*fc%nr3t, sca1/(dble(fc%nr1t*fc%nr2t*fc%nr3t))
        
  
       
     enddo
   
!loop on fake conduction states
 


     call mp_barrier( world_comm )
        
     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory9', kilobytes
  
     do ii=1,num_fc_eff(is)
        if(.not.l_iter_algorithm) then
           psic(:)=(0.d0,0.d0)
           psic(fc%nlt(igkt(1:fc%npwt)))  = state_fc_t(1:fc%npwt,ii,is)
           psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( state_fc_t(1:fc%npwt,ii,is) )
           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
           state_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))
        endif

   
     
        do iv=1,num_nbndv(is),2
!form product in real space
           if(.not.l_iter_algorithm) then
              state_real_tmp(:)=state_real(:)*wv_real_all(:,iv)
              if(iv < num_nbndv(is)) then
                 state_real_tmp2(:)=state_real(:)*wv_real_all(:,iv+1)
              else
                 state_real_tmp2(:)=0.d0
              endif
           else
              state_real_tmp(1:fc%nrxxt)=state_fc_r(1:fc%nrxxt,ii,is)*wv_real_all(1:fc%nrxxt,iv)
              if(iv < num_nbndv(is)) then
                 state_real_tmp2(1:fc%nrxxt)=state_fc_r(1:fc%nrxxt,ii,is)*wv_real_all(1:fc%nrxxt,iv+1)
              else
                 state_real_tmp2(1:fc%nrxxt)=0.d0
              endif
           endif
!back to G space


           psic(1:fc%nrxxt)=dcmplx(state_real_tmp(1:fc%nrxxt),state_real_tmp2(1:fc%nrxxt))
           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
              

           if(iv < num_nbndv(is)) then
              state_g(1:fc%npwt, iv)= 0.5d0*(psic(fc%nlt(igkt(1:fc%npwt)))+conjg( psic(fc%nltm(igkt(1:fc%npwt)))))
              state_g(1:fc%npwt, iv+1)= (0.d0,-0.5d0)*(psic(fc%nlt(igkt(1:fc%npwt))) - conjg(psic(fc%nltm(igkt(1:fc%npwt)))))
           else
              state_g(1:fc%npwt, iv) = psic(fc%nlt(igkt(1:fc%npwt)))
           endif
           
        enddo
    



!do products with fcw states


        call dgemm('T','N',fcw_number,num_nbndv(is),2*fc%npwt,2.d0,fcw_state,2*fc%npwt,state_g,2*fc%npwt,0.d0,tmp_mat,fcw_number)
        if(fc%gstart_t==2) then
           do iv=1,num_nbndv(is)
              do jj=1,fcw_number
                 tmp_mat(jj,iv)=tmp_mat(jj,iv)-dble(conjg(fcw_state(1,jj))*state_g(1,iv))
              enddo
           enddo
        endif
 
        call mp_sum(tmp_mat,world_comm)

        if(l_frac) then
           if(ii<=num_fc) then
              do iv=1,num_nbndv(is)
                 if(nspin==1) then
                    sca1=dsqrt(abs(wg(iv,is))/2.d0)
                 else
                    sca1=dsqrt(abs(wg(iv,is)))
                 endif
                 tmp_mat(1:fcw_number,iv)=tmp_mat(1:fcw_number,iv)*sca1
              enddo
           else
              do iv=1,num_nbndv(is)
                 if(nspin==1) then
                    sca1=dsqrt(abs(wg(ii-num_fc+num_nbndv_min(is),is)-wg(iv,is))/2.d0)
                 else
                    sca1=dsqrt(abs(wg(ii-num_fc+num_nbndv_min(is),is)-wg(iv,is)))
                 endif
                 tmp_mat(1:fcw_number,iv)=tmp_mat(1:fcw_number,iv)*sca1
              enddo
           endif
        endif

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory10', kilobytes
        if(l_verbose) write(stdout,*) 'TOTAL NUMBER OF FCW STATES:', fcw_number,ii,dfftp%nnr,fc%nrxxt,wg(1,is)
        FLUSH(stdout)
        call mp_barrier( world_comm )
     
!calculate contribution to D matrix
        if(nsize>0) then
           call dgemm('N','T',fcw_number,nend-nbegin+1,num_nbndv(is),1.d0,tmp_mat,fcw_number,&
                &tmp_mat(nbegin:nend,1:num_nbndv(is)),nend-nbegin+1,1.d0,fcw_mat,fcw_number)
           endif
           if(l_test) then
              do iv=1,num_nbndv(is)
                 sca1=0.d0
                 do jj=1,fcw_number
                    sca1=sca1+tmp_mat(jj,iv)**2.d0
                 enddo
                 sca2=0.d0
                 do ig=1,fc%npwt
                    sca2=sca2+2.d0*dble(conjg(state_g(ig,iv))*state_g(ig,iv))
                 enddo
                 if(fc%gstart_t==2) sca2=sca2-dble(state_g(1,iv))**2.d0
                 call mp_sum(sca2,world_comm)
                 write(stdout,*) 'Projection', ii,iv,sca1/sca2
              enddo
           endif
        enddo
     enddo!on spin

!if required put fcw_state on normconserving ordering

     allocate(fcw_state_n(npw,fcw_number))

     if(fc%dual_t==4.d0) then
        fcw_state_n(1:fc%npwt,1:fcw_number)=fcw_state(1:fc%npwt,1:fcw_number)
     else
        call reorderwfp_col(fcw_number,fc%npwt,npw,fcw_state(1,1),fcw_state_n(1,1),fc%npwt,npw, &
             & fc%ig_l2gt,ig_l2g,fc%ngmt_g,mpime, nproc,intra_pool_comm )
   !     do ii=1,fcw_number
   !          call mergewf(fcw_state(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
   !         call splitwf(fcw_state_n(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
   !     enddo
     endif

   CALL memstat( kilobytes )
   if(l_verbose) write(stdout,*) 'memory11', kilobytes
 

  !save on file

  iunfcw = find_free_unit()
  CALL diropn( iunfcw, 'fcw', npw*2, exst )
  do ii=1,fcw_number
     CALL davcio( fcw_state_n(1,ii), 2*npw, iunfcw, ii, 1 )
  enddo
  close(iunfcw)

  !write number of states

  if(ionode) then
     open(unit=iunfcw,file=trim(tmp_dir)//trim(prefix)//'.nfcws',status='unknown')
     write(iunfcw,*) fcw_number
     close(iunfcw)
  endif

  CALL diropn( iunfcw, 'fmat',fcw_number, exst )
  do ii=1,nsize
     CALL davcio( fcw_mat(1,ii), fcw_number, iunfcw, ii, 1 )
  enddo
  close(iunfcw)

  if(l_verbose) write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  call deallocate_fft_custom(fc)



  iunrestart0 =  find_free_unit()
  open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='unknown')
  write(iunrestart0,*) -1
  write(iunrestart0,*) fcw_number
  write(iunrestart0,*) fcw_numberx
  close(iunrestart0)


  deallocate(wv_real_all)

  
  
  if( allocated( state_fc_t ) ) deallocate( state_fc_t )
  deallocate(state_real,state_g,state_real_tmp,state_real_tmp2)
  deallocate(tmp_mat)
  if(allocated(e_fake)) deallocate(e_fake)
  deallocate(fcw_state_n)
  deallocate(evc_t)


  if( allocated( state_fc ) ) deallocate( state_fc )
  if( allocated( state_g ) ) deallocate( state_g )
  if( allocated( fcw_state_old ) ) deallocate( fcw_state_old )
  if( allocated( h_state_fc ) ) deallocate( h_state_fc ) 
  if( allocated( evc_g ) ) deallocate( evc_g )
  if( allocated( evc_t ) ) deallocate( evc_t )
  if( allocated( state_fc_t ) ) deallocate( state_fc_t )
  if( allocated( state_g_t ) ) deallocate( state_g_t ) 
  if( allocated( fcw_state_n ) ) deallocate( fcw_state_n )
  if( allocated( wv_real ) ) deallocate( wv_real )
  if( allocated( state_real ) ) deallocate( state_real )
  if( allocated( wv_real_all ) ) deallocate( wv_real_all )
  if( allocated( state_real_tmp ) ) deallocate( state_real_tmp )
  if( allocated( state_real_tmp2 ) ) deallocate( state_real_tmp2 )
  if( allocated( state_real2 ) ) deallocate( state_real2 )
  if( allocated( omat ) ) deallocate( omat )
  if( allocated( eigen ) ) deallocate( eigen )
  if( allocated( work ) ) deallocate( work )
  if( allocated( tmp_mat ) ) deallocate( tmp_mat )
  if( allocated( omat2 ) ) deallocate( omat2 )
  if( allocated( hmat ) ) deallocate( hmat )
  if( allocated( e_fake ) ) deallocate( e_fake )
  if( allocated( vec_fake ) ) deallocate( vec_fake )
  if( allocated( gap ) ) deallocate( gap )
  if( allocated( hmat_i ) ) deallocate( hmat_i )
  if( allocated( hmat_o ) ) deallocate( hmat_o )
  if( allocated( omat_i ) ) deallocate( omat_i )
  if( allocated( ovec ) ) deallocate( ovec )
  if( allocated( g2kint ) ) deallocate( g2kint )
  !
  if( allocated( iwork ) )   deallocate( iwork )
  if( allocated( ifail ) )   deallocate( ifail )
  if( allocated( isuppz ) )  deallocate( isuppz )
  if( allocated( iclustr ) ) deallocate( iclustr )
  if( allocated( igkt ) )    deallocate( igkt )
  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory12', kilobytes
  if(l_verbose) write(stdout,*) 'memory fcw_state = ',  SIZE( fcw_state ) / 64 , ' kb'
  if(l_verbose) write(stdout,*) 'memory fcw_mat   = ',  SIZE( fcw_mat   ) / 64 , ' kb'
  FLUSH(stdout)

  return

end subroutine fake_conduction_wannier







subroutine fake_conduction_wannier_real( cutoff, s_cutoff )

 !IT WORKS ONLY FOR NORMCONSERVING PSEUDOPOTENTIALS
  !the valence states in G space must be in evc
  ! Gamma point version
 !real space version

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
   USE gvecw,    ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm, mpime, nproc
   USE mp_pools, ONLY : intra_pool_comm
   USE gvecs,              ONLY : nls, nlsm,  doublegrid
   USE kinds, ONLY : DP
   USE io_files, ONLY : prefix, tmp_dir, diropn
   USE g_psi_mod,            ONLY : h_diag, s_diag
   USE noncollin_module,     ONLY : noncolin, npol
   USE becmod,           ONLY : becp
   USE uspp,                 ONLY : vkb, nkb, okvan
   USE klist,                ONLY : xk,igk_k
   USE fft_custom_gwl
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

!  INTEGER,INTENT(out) :: fcw_number!number of "fake conduction" states for O matrix method
!  COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: fcw_state! "fake conduction" states for O matrix method
!  REAL(kind=DP),    POINTER, DIMENSION(:,:) :: fcw_mat! "fake conduction" matrix

  REAL(kind=DP), INTENT(in) :: cutoff!cutoff for planewaves
  REAL(kind=DP), INTENT(in) :: s_cutoff!cutoff for orthonormalization
!NOT_TO_BE_INCLUDED_START

  COMPLEX(kind=DP), ALLOCATABLE :: state_fc(:,:)
  !COMPLEX(kind=DP), ALLOCATABLE :: state_g(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: fcw_state_old(:,:) 
  COMPLEX(kind=DP), ALLOCATABLE :: h_state_fc(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:),evc_t(:,:),state_fc_t(:,:),state_g_t(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: fcw_state_n(:,:)
  REAL(kind=DP), ALLOCATABLE :: wv_real(:),state_real(:),state_real_tmp(:)
  REAL(kind=DP), ALLOCATABLE :: state_real_tmp2(:),state_real2(:)
  REAL(kind=DP), ALLOCATABLE :: omat(:,:)
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:),tmp_mat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: omat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: hmat(:,:)
  REAL(kind=DP), ALLOCATABLE :: e_fake(:), vec_fake(:,:)
  REAL(kind=DP), ALLOCATABLE :: gap(:)
  REAL(kind=DP), ALLOCATABLE :: hmat_i(:,:),hmat_o(:,:), omat_i(:,:)
  REAL(kind=DP), ALLOCATABLE :: ovec(:)
  REAL(kind=DP), ALLOCATABLE :: g2kint(:)
  INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
  INTEGER, ALLOCATABLE :: isuppz(:)
  INTEGER, ALLOCATABLE :: iclustr(:)
  INTEGER, ALLOCATABLE :: igkt(:)

  REAL(kind=DP):: sca1,sca2

  LOGICAL :: l_test=.false.!if true test the completness of the basis
  LOGICAL :: l_dsyevr=.true.!if true uses dsyevr instead of dsyev
  LOGICAL :: l_diago_cg=.true.!if true uses diago_cg instead of  dsyevr ATTENZIONE
  LOGICAL :: exst
  LOGICAL :: l_dsygvx=.false.!if .true. uses serial dsygvx instead of parallel diago_cg_g
  LOGICAL :: l_gramsc=.true.!if true orthonormalization through gram-schimdt
  LOGICAL :: l_diago_para=.true.!if true uses parallel diago_cg
  LOGICAL :: l_fft_custom=.false.

  INTEGER :: ig,ip, ii, iv, jj, iw, ir
  INTEGER :: num_fc!number of fake conduction states
  INTEGER :: lwork,info,liwork
  INTEGER :: n_out
  INTEGER :: fcw_number_old
  INTEGER :: l_blk,nbegin,nend
  INTEGER :: max_state
  INTEGER :: iunfcw
  INTEGER :: nsize
  INTEGER :: nbegin_loc,nend_loc,nsize_loc
  INTEGER :: n_found_state
!variables for scalapack
  INTEGER :: num_fc_r,num_fc_c,num_fc_dimr,num_fc_dimc
  INTEGER :: m,nz,icrow,iccol,iproc,ilrow,ilcol
  INTEGER :: desc_a(9),desc_b(9),desc_c(9)
  INTEGER :: n_computed
  INTEGER :: num_built!number of states already built
  INTEGER :: num_out
  INTEGER ::  kilobytes
  INTEGER ::  kb_old, kb_new

  INTEGER, EXTERNAL :: indxg2p,indxg2l

  TYPE(optimal_options) :: options



  INTEGER :: bufferx,fcw_numberx,fcw_number_oldx, fcw_numberx_tmp


  LOGICAL :: l_restart0!if true restart is enabled
  INTEGER :: iunrestart0, iv_start,iunfsr
  REAL(kind=DP), ALLOCATABLE :: state_fc_r(:,:)

  INTEGER, ALLOCATABLE :: g_to_loc(:) !global to local correspondance
  INTEGER, ALLOCATABLE :: loc_to_g(:)
  INTEGER :: n_loc!number of local r points
  REAL(kind=DP), ALLOCATABLE:: state_r(:,:),wv_real_loc(:),state_loc(:)
  REAL(kind=DP), ALLOCATABLE :: fcw_state_r(:,:), fcw_state_r_loc(:,:)
  INTEGER :: nmod
  REAL(kind=DP), ALLOCATABLE :: tmp_vec(:)


  TYPE(fft_cus) :: fc

!determine bufferx,fcw_numberx
  bufferx=num_nbndv(1)*300/4
  bufferx=max(1000,bufferx)
!ONLY FOR PROJECT ON JADE
  bufferx=5000
  fcw_numberx=bufferx
  fcw_number_oldx=bufferx
  fcw_numberx_tmp=bufferx

!generate fake conduction states
!!determine number of states

!generate custom in grid in case can be equal to norm-conserving grid

  fc%ecutt=ecutwfc
  fc%dual_t=dual_pb

  write(stdout,*) 'Call initialize_fft_custom'
  CALL memstat( kilobytes )
  write(stdout,*) 'memory0', kilobytes
  FLUSH(stdout)

  call initialize_fft_custom(fc)

  CALL memstat( kilobytes )
  write(stdout,*) 'memory0.0', kilobytes
  FLUSH(stdout)

  ! this is for compatibility

  allocate( igkt( fc%npwt ) )
  do ig=1,fc%npwt
     igkt(ig)=ig
  enddo

  allocate( evc_g( fc%ngmt_g ) )

  !plane waves basis set

  !state_fc are first obtained on the ordering of the normconserving grid

  g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
       ( g(2,igk_k(1:npw,1)) )**2 + &
       ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2
  
  num_fc=0
  do ig=1,npw
     if(g2kin(ig) <= cutoff) num_fc=num_fc+1
  enddo
  call mp_sum(num_fc,world_comm)
  num_fc=(num_fc-1)*2+1

  allocate( state_fc( npw, num_fc ) )

  state_fc(:,:)=(0.d0,0.d0)
  write(stdout,*) "Number of fake conduction states:", num_fc
  CALL memstat( kilobytes )
  write(stdout,*) 'memory0.1', kilobytes, ' new kb = ', (SIZE( state_fc )*16 + SIZE( evc_g )*16 + SIZE( igkt )*4)/1024
  FLUSH(stdout)
  
  ii=0
  do ip=0,nproc-1
     if(mpime==ip) then
        do ig=gstart,npw
           if(g2kin(ig) <= cutoff) then
              ii=ii+1
              state_fc(ig,ii)=cmplx(dsqrt(0.5d0),0.d0)
              ii=ii+1
              state_fc(ig,ii)=cmplx(0.d0,dsqrt(0.5d0))
           endif
        enddo
        if(gstart==2) then
           ii=ii+1
           state_fc(1,ii)=(1.d0,0.d0)
        endif
     else
        ii=0
     endif
     call mp_sum(ii,world_comm)
  enddo

  if(ii/=num_fc) then 
     write(stdout,*) 'ERRORE FAKE CONDUCTION',ii
     FLUSH(stdout)
     stop
    return
  endif
  write(stdout,*) 'FAKE1'
  FLUSH(stdout)

!!project out of valence space
  do ii=1,num_fc
     call pc_operator(state_fc(:,ii),1,.false.)!ATTENZIONE spin not implemented yet here
    ! if(gstart==2) write(stdout,*) 'FAKE modulus', ii, state_fc(1,ii)
  enddo



!orthonormalize fake_conduction states 

 

     
!for the moment finds all the first fcw_fast_n eigenstates
     
     write(stdout,*) 'CASE ORTHONORMALIZATION ONLY'
     FLUSH(stdout)
     
     options%l_complete=.true.
     options%idiago=0
     call optimal_driver(num_fc,state_fc,npw,options,num_out, info)
      
  CALL memstat( kilobytes )
  write(stdout,*) 'memory0.3', kilobytes
  FLUSH(stdout)
        
       
     
  !now state_fc are put on the ordering of the redueced grid, if required
  allocate(state_fc_t(fc%npwt,num_fc))
  if(fc%dual_t==4.d0) then
     state_fc_t(:,:)=state_fc(:,:)
  else
     do ii=1,num_fc
        call mergewf(state_fc(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
        call splitwf(state_fc_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     enddo
  endif

  allocate(state_fc_r(fc%nrxxt,num_fc))
  do ii=1,num_fc,2
     psic(:)=(0.d0,0.d0)
     if(ii==num_fc) then
        psic(fc%nlt(1:fc%npwt))  = state_fc_t(1:fc%npwt,ii)
        psic(fc%nltm(1:fc%npwt)) = CONJG( state_fc_t(1:fc%npwt,ii) )
     else
        psic(fc%nlt(1:fc%npwt))=state_fc_t(1:fc%npwt,ii)+(0.d0,1.d0)*state_fc_t(1:fc%npwt,ii+1)
        psic(fc%nltm(1:fc%npwt)) = CONJG( state_fc_t(1:fc%npwt,ii) )+(0.d0,1.d0)*CONJG( state_fc_t(1:fc%npwt,ii+1) )
     endif
     CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
     state_fc_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
     if(ii/=num_fc) state_fc_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
  enddo
  

  CALL memstat( kilobytes )
  write(stdout,*) 'memory0.4', kilobytes, ' new kb = ', (SIZE( state_fc_t ))/64
  FLUSH(stdout)

!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbndv(1)))
  if(fc%dual_t==4.d0) then
     evc_t(:,1:num_nbndv(1))=evc(:,1:num_nbndv(1))
  else
     do iv=1,num_nbndv(1)
        call mergewf(evc(:,iv),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
        call splitwf(evc_t(:,iv),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     enddo
  endif

  CALL memstat( kilobytes )
  write(stdout,*) 'memory0.5', kilobytes, ' new kb = ', (SIZE( evc_t ))/64
  FLUSH(stdout)

!cycle on v 
!! product in real space with wannier
!! orthonormalize and take N most important
!! gram-schmidt like 

!calculate D matrix 

  l_blk= (num_fc)/nproc
  if(l_blk*nproc < (num_fc)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > num_fc) nend=num_fc
  nsize=nend-nbegin+1

!check for restart
  allocate(fcw_state_r(fc%nrxxt,fcw_numberx))

  

  allocate (wv_real(fc%nrxxt))!,state_g(fc%npwt,num_fc))
  allocate(wv_real_loc(fc%nrxxt))
!  allocate(state_r(fc%nrxxt,num_fc))
  allocate(g_to_loc(fc%nrxxt),loc_to_g(fc%nrxxt))

  fcw_number=0
  FIRST_LOOP: do iv=1,num_nbndv(1)
     call mp_barrier( world_comm )
     write(stdout,*) 'FK state:', iv,fc%nrxxt,fc%npwt,num_fc
     CALL memstat( kilobytes )
     write(stdout,*) 'memory1', kilobytes
     FLUSH(stdout)


     

  !   allocate (wv_real(fc%nrxxt),state_real(fc%nrxxt),state_real2(fc%nrxxt),state_g(fc%npwt,num_fc))
  !   if(l_iter_algorithm) allocate (state_g_r(fc%nrxxt,num_fc))
     call mp_barrier( world_comm )
     CALL memstat( kilobytes )
     write(stdout,*) 'memory2', kilobytes, ' new kb = ', &
                     (SIZE(wv_real))/128
     FLUSH(stdout)

     psic(:)=(0.d0,0.d0)
     psic(fc%nlt(igkt(1:fc%npwt)))  = evc_t(1:fc%npwt,iv)
     psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( evc_t(1:fc%npwt,iv) )

     call mp_barrier( world_comm )
     CALL memstat( kilobytes )
     write(stdout,*) 'memory3', kilobytes
     FLUSH(stdout)

    CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
    call mp_barrier( world_comm )
     CALL memstat( kilobytes )
     write(stdout,*) 'memory4', kilobytes
     FLUSH(stdout)

     
     wv_real(:)= DBLE(psic(:))
     if(fc%gstart_t==2) write(stdout,*) 'FAKE modulus valence', iv, evc_t(1,iv)
!loop on fake conduction states

     
!find global to local correspodance
     n_loc=0
     do ir=1,fc%nrxxt
        if(wv_real(ir)**2.d0 >= wannier_thres) then
           n_loc=n_loc+1
           g_to_loc(ir)=n_loc
           loc_to_g(n_loc)=ir
        else
           g_to_loc(ir)=0
        endif
     enddo

     write(stdout,*) 'Start products',n_loc,fc%nrxxt,loc_to_g(n_loc)
     FLUSH(stdout)
    
     do ir=1,n_loc
        wv_real_loc(ir)=wv_real(loc_to_g(ir))
     enddo

     if(n_loc>= 1) then
        allocate(state_loc(n_loc))
     else
        allocate(state_loc(1))
     endif


     

     write(stdout,*) 'End products part'
     call mp_barrier( world_comm )
     CALL memstat( kilobytes )
     write(stdout,*) 'memory5', kilobytes
     FLUSH(stdout)

!if not first valence wannier project the products out of partial orthonormal basis


!loop on fake conduction states
     allocate(tmp_vec(fcw_numberx))
     if(n_loc >=1 ) then
        allocate(fcw_state_r_loc(n_loc,fcw_numberx))
     else
        allocate(fcw_state_r_loc(1,fcw_numberx))
     endif
!$OMP PARALLEL SHARED(fcw_number,n_loc,fcw_state_r,fcw_state_r_loc,loc_to_g) PRIVATE(ii,ir)
!$OMP DO
     do ii=1,fcw_number
        do ir=1,n_loc
           fcw_state_r_loc(ir,ii)=fcw_state_r(loc_to_g(ir),ii)
        enddo
     enddo
!$OMP END DO 
!$OMP END PARALLEL
 
     

     do ii=1,num_fc

        
 !global to local trasform
        do ir=1,n_loc
           state_loc(ir)=state_fc_r(loc_to_g(ir),ii)
        enddo
        do ir=1,n_loc
           state_loc(ir)=state_loc(ir)*wv_real_loc(ir)
        enddo

     
    
        if(iv==1 .and. ii==1) then
!put the first one as it is
!calculate modulus
           
           if(n_loc >=1) then
              call dgemm('T','N',1,1,n_loc,1.d0,state_loc,n_loc,state_loc,n_loc,0.d0,sca2,1)
           else
              sca2=0.d0
           endif
           call mp_sum(sca2,world_comm)
           sca2=sca2/dble(fc%nr1t*fc%nr2t*fc%nr3t)
           sca2=1.d0/dsqrt(sca2)
           if(n_loc >= 1) fcw_state_r_loc(1:n_loc,1)=state_loc(1:n_loc)*sca2
           fcw_state_r(:,1)=0.d0
!$OMP PARALLEL SHARED(fcw_state_r,fcw_state_r_loc,loc_to_g,n_loc) PRIVATE(ir)
!$OMP DO
           do ir=1,n_loc
              fcw_state_r(loc_to_g(ir),1)=fcw_state_r_loc(ir,1)
           enddo
!$OMP END DO 
!$OMP END PARALLEL
           fcw_number=1
        else
         
            if(n_loc >=1) then
               call dgemm('T','N',fcw_number,1,n_loc,1.d0,fcw_state_r_loc,n_loc,state_loc,n_loc,0.d0,tmp_vec,fcw_numberx)
            else
               tmp_vec(1:fcw_number)=0.d0
            endif
           call mp_sum(tmp_vec,world_comm)
           tmp_vec(:)=tmp_vec(:)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
           if(n_loc >=1) then
              call dgemm('T','N',1,1,n_loc,1.d0,state_loc,n_loc,state_loc,n_loc,0.d0,sca2,1)
           else
              sca2=0.d0
           endif
           call mp_sum(sca2,world_comm)
           sca2=sca2/dble(fc%nr1t*fc%nr2t*fc%nr3t)
           sca1=0.d0
          
     
           do jj=1,fcw_number
              sca1=sca1+tmp_vec(jj)**2.d0
           enddo
           if(abs(sca2-sca1) >= s_cutoff) then
              fcw_state_r(:,fcw_number+1)=0.d0
!$OMP PARALLEL SHARED(fcw_state_r,state_loc,loc_to_g,fcw_number,n_loc) PRIVATE(ir)
!$OMP DO
              do ir=1,n_loc
                 fcw_state_r(loc_to_g(ir),fcw_number+1)=state_loc(ir)
              enddo
!$OMP END DO 
!$OMP END PARALLEL
              call dgemm('N','N',fc%nrxxt, 1,fcw_number,-1.d0,fcw_state_r,fc%nrxxt,tmp_vec,&
     &fcw_numberx,1.d0,fcw_state_r(1,fcw_number+1),fc%nrxxt)
              sca1=1.d0/(dsqrt(abs(sca2-sca1)))
              fcw_state_r(:,fcw_number+1)= fcw_state_r(:,fcw_number+1)*sca1
!$OMP PARALLEL SHARED(fcw_state_r,fcw_state_r_loc,loc_to_g,fcw_number,n_loc) PRIVATE(ir)
!$OMP DO
              do ir=1,n_loc
                 fcw_state_r_loc(ir,fcw_number+1)=fcw_state_r(loc_to_g(ir),fcw_number+1)
              enddo
!$OMP END DO 
!$OMP END PARALLEL
              fcw_number=fcw_number+1
           endif
        endif


     enddo

     deallocate(tmp_vec)
     deallocate(fcw_state_r_loc)
     deallocate(state_loc)

     FLUSH(stdout)

     write(stdout,*) 'memory6.8', kilobytes
     FLUSH(stdout)
  end do FIRST_LOOP
 

 
  write(stdout,*) 'FK8'
  CALL memstat( kilobytes )
  write(stdout,*) 'memory7', kilobytes
  FLUSH(stdout)
  

 
! calculate D matrix distributed among processors

  if(fcw_number < nproc) then
     write(stdout,*) 'too many processors'
     stop
  endif

  l_blk= (fcw_number)/nproc
  if(l_blk*nproc < (fcw_number)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > fcw_number) nend=fcw_number
  nsize=nend-nbegin+1

  allocate(fcw_mat(fcw_number,l_blk))
  fcw_mat(:,:)=0.d0

  write(stdout,*) 'FK9'
  FLUSH(stdout)

 

  CALL memstat( kilobytes )
  write(stdout,*) 'memory8', kilobytes
  
  allocate(tmp_mat(fcw_number,100))
  do iv=1,num_nbndv(1)
      
     
     psic(:)=(0.d0,0.d0)
     psic(fc%nlt(1:fc%npwt))  = evc_t(1:fc%npwt,iv)
     psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )
     
     CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
     wv_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))

     
!find global to local correspodance
     n_loc=0
     do ir=1,fc%nrxxt
        if(wv_real(ir)**2.d0 >= wannier_thres) then
           n_loc=n_loc+1
           g_to_loc(ir)=n_loc
           loc_to_g(n_loc)=ir
        else
           g_to_loc(ir)=0
        endif
     enddo
     
        
     do ir=1,n_loc
        wv_real_loc(ir)=wv_real(loc_to_g(ir))
     enddo

!put fcw_state_r and state_fc_r on local grid

     if(n_loc >= 1) then
        allocate(state_r(n_loc,num_fc))
        allocate(fcw_state_r_loc(n_loc, fcw_number))
     else
         allocate(state_r(n_loc,num_fc))
        allocate(fcw_state_r_loc(n_loc, fcw_number))
     endif

     do ii=1,num_fc
 !global to local trasform
        do ir=1,n_loc
           state_r(ir,ii)=state_fc_r(loc_to_g(ir),ii)*wv_real_loc(ir)
        enddo
     enddo

     do ii=1,fcw_number
        do ir=1,n_loc
           fcw_state_r_loc(ir,ii)=fcw_state_r(loc_to_g(ir),ii)
        enddo
     enddo
     
     do ii=1,num_fc,100
        nmod=min(ii+100-1,num_fc)-ii+1
        if(n_loc >= 1 ) then
           call dgemm('T','N',fcw_number,nmod,n_loc,1.d0,fcw_state_r_loc,n_loc,state_r(1,ii),n_loc,0.d0,tmp_mat,fcw_number)
        else
           tmp_mat(:,:)=0.d0
        endif
        call mp_sum(tmp_mat,world_comm)
        tmp_mat(:,:)=tmp_mat(:,:)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
     
        CALL memstat( kilobytes )
        write(stdout,*) 'memory10', kilobytes
        write(stdout,*) 'TOTAL NUMBER OF FCW STATES:', fcw_number,ii,dfftp%nnr,fc%nrxxt,n_loc
        FLUSH(stdout)
        if(nsize>0) then
           call dgemm('N','T',fcw_number,nend-nbegin+1,nmod,1.d0,tmp_mat,fcw_number,&
     &tmp_mat(nbegin:nend,1:nmod),nend-nbegin+1,1.d0,fcw_mat,fcw_number)
        endif
     enddo
     deallocate(state_r)
     deallocate(fcw_state_r_loc)
  enddo
        

!trasform fcw_state_r to fcw_state

  allocate(fcw_state(fc%npwt,fcw_number))
  do ii=1,fcw_number,2
     if(ii==fcw_number) then
        psic(1:fc%nrxxt)=dcmplx(fcw_state_r(1:fc%nrxxt,ii),0.d0)
     else
        psic(1:fc%nrxxt)=dcmplx(fcw_state_r(1:fc%nrxxt,ii),fcw_state_r(1:fc%nrxxt,ii+1))
     endif
     CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
     if(ii==fcw_number) then
        fcw_state(1:fc%npwt, ii) = psic(fc%nlt(1:fc%npwt))
        if(fc%gstart_t==2) fcw_state(1,ii)=(0.d0,0.d0)
     else
        fcw_state(1:fc%npwt, ii)= 0.5d0*(psic(fc%nlt(igkt(1:fc%npwt)))+conjg( psic(fc%nltm(igkt(1:fc%npwt)))))
        fcw_state(1:fc%npwt, ii+1)= (0.d0,-0.5d0)*(psic(fc%nlt(igkt(1:fc%npwt))) - conjg(psic(fc%nltm(igkt(1:fc%npwt)))))
        if(fc%gstart_t==2) fcw_state(1,ii)=(0.d0,0.d0)
        if(fc%gstart_t==2) fcw_state(1,ii+1)=(0.d0,0.d0)
     endif
  enddo



    

  
  

  !write(stdout,*) 'Att0'
  !   FLUSH(stdout)
  call mp_barrier( world_comm )
  
  CALL memstat( kilobytes )
  write(stdout,*) 'memory9', kilobytes
  FLUSH(stdout)

!if required put fcw_state on normconserving ordering

  deallocate(fcw_state_r)
  allocate(fcw_state_n(npw,fcw_number))

  if(fc%dual_t==4.d0) then
     fcw_state_n(:,:)=fcw_state(:,:)
  else
     do ii=1,fcw_number
         call mergewf(fcw_state(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
         call splitwf(fcw_state_n(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
      enddo
   endif

   CALL memstat( kilobytes )
   write(stdout,*) 'memory11', kilobytes
 

  !save on file

  iunfcw = find_free_unit()
  CALL diropn( iunfcw, 'fcw', npw*2, exst )
  do ii=1,fcw_number
     CALL davcio( fcw_state_n(1,ii), 2*npw, iunfcw, ii, 1 )
  enddo
  close(iunfcw)

  !write number of states

  if(ionode) then
     open(unit=iunfcw,file=trim(tmp_dir)//trim(prefix)//'.nfcws',status='unknown')
     write(iunfcw,*) fcw_number
     close(iunfcw)
  endif

  CALL diropn( iunfcw, 'fmat',fcw_number, exst )
  do ii=1,nsize
     CALL davcio( fcw_mat(1,ii), fcw_number, iunfcw, ii, 1 )
  enddo
  close(iunfcw)

  write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  call deallocate_fft_custom(fc)


 
  
  deallocate(g_to_loc,loc_to_g)
  deallocate(wv_real_loc)
  
  if( allocated( state_fc_t ) ) deallocate( state_fc_t )
  deallocate(tmp_mat)
  if(allocated(e_fake)) deallocate(e_fake)
  deallocate(fcw_state_n)
  deallocate(evc_g,evc_t)


  if( allocated( state_fc ) ) deallocate( state_fc )
  if( allocated( fcw_state_old ) ) deallocate( fcw_state_old )
  if( allocated( h_state_fc ) ) deallocate( h_state_fc ) 
  if( allocated( evc_g ) ) deallocate( evc_g )
  if( allocated( evc_t ) ) deallocate( evc_t )
  if( allocated( state_fc_t ) ) deallocate( state_fc_t )
  if( allocated( state_g_t ) ) deallocate( state_g_t ) 
  if( allocated( fcw_state_n ) ) deallocate( fcw_state_n )
  if( allocated( wv_real ) ) deallocate( wv_real )
  if( allocated( state_real ) ) deallocate( state_real )
  if( allocated( state_real_tmp ) ) deallocate( state_real_tmp )
  if( allocated( state_real_tmp2 ) ) deallocate( state_real_tmp2 )
  if( allocated( state_real2 ) ) deallocate( state_real2 )
  if( allocated( omat ) ) deallocate( omat )
  if( allocated( eigen ) ) deallocate( eigen )
  if( allocated( work ) ) deallocate( work )
  if( allocated( tmp_mat ) ) deallocate( tmp_mat )
  if( allocated( omat2 ) ) deallocate( omat2 )
  if( allocated( hmat ) ) deallocate( hmat )
  if( allocated( e_fake ) ) deallocate( e_fake )
  if( allocated( vec_fake ) ) deallocate( vec_fake )
  if( allocated( gap ) ) deallocate( gap )
  if( allocated( hmat_i ) ) deallocate( hmat_i )
  if( allocated( hmat_o ) ) deallocate( hmat_o )
  if( allocated( omat_i ) ) deallocate( omat_i )
  if( allocated( ovec ) ) deallocate( ovec )
  if( allocated( g2kint ) ) deallocate( g2kint )
  !
  if( allocated( iwork ) )   deallocate( iwork )
  if( allocated( ifail ) )   deallocate( ifail )
  if( allocated( isuppz ) )  deallocate( isuppz )
  if( allocated( iclustr ) ) deallocate( iclustr )
  if( allocated( igkt ) )    deallocate( igkt )

  CALL memstat( kilobytes )
  write(stdout,*) 'memory12', kilobytes
  write(stdout,*) 'memory fcw_state = ',  SIZE( fcw_state ) / 64 , ' kb'
  write(stdout,*) 'memory fcw_mat   = ',  SIZE( fcw_mat   ) / 64 , ' kb'
  FLUSH(stdout)

  return
!NOT_TO_BE_INCLUDED_END
end subroutine fake_conduction_wannier_real


  subroutine fake_conduction_real( cutoff, s_cutoff,ks_wfcs ,l_frac, ks_wfcs_diag,l_cond)

 !IT WORKS ONLY FOR NORMCONSERVING PSEUDOPOTENTIALS
  !the valence states in G space must be in evc
  ! Gamma point version

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et, wg
   USE gvecw,    ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm, mpime, nproc
   USE mp_pools, ONLY : intra_pool_comm
   USE gvecs,              ONLY : nls, nlsm, doublegrid

   USE kinds, ONLY : DP
   USE io_files, ONLY : prefix, tmp_dir, diropn
   USE g_psi_mod,            ONLY : h_diag, s_diag
   USE noncollin_module,     ONLY : noncolin, npol
   USE becmod,           ONLY : becp
   USE uspp,                 ONLY : vkb, nkb, okvan
   USE klist,                ONLY : xk,igk_k
   USE fft_custom_gwl
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp
   USE lsda_mod,             ONLY : nspin
   USE mp_wave_parallel


  implicit none

  INTEGER, EXTERNAL :: find_free_unit

!  INTEGER,INTENT(out) :: fcw_number!number of "fake conduction" states for O matrix method
!  COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: fcw_state! "fake conduction" states for O matrix method
!  REAL(kind=DP),    POINTER, DIMENSION(:,:) :: fcw_mat! "fake conduction" matrix

  REAL(kind=DP), INTENT(in) :: cutoff!cutoff for planewaves
  REAL(kind=DP), INTENT(in) :: s_cutoff!cutoff for orthonormalization
  COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npwx,nbnd,nspin)!Kohn-Sham or Wannier wavefunctios
  LOGICAL, INTENT(in) :: l_frac!if true consider fractional occupancies
  COMPLEX(kind=DP), INTENT(in) :: ks_wfcs_diag(npwx,nbnd,nspin)!Kohn-Sham wavefunctios 
  LOGICAL, INTENT(in) :: l_cond!if true consider also conduction states for the construction of the polarizability basis
!NOT_TO_BE_INCLUDED_START

  COMPLEX(kind=DP), ALLOCATABLE :: state_fc(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: state_g(:,:)
  REAL(kind=DP), ALLOCATABLE :: fcw_state_r(:,:)
  REAL(kind=DP), ALLOCATABLE :: fcw_state_old_r(:,:) 
  COMPLEX(kind=DP), ALLOCATABLE :: h_state_fc(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:),evc_t(:,:,:),state_fc_t(:,:,:),state_g_t(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: fcw_state_n(:,:)
  REAL(kind=DP), ALLOCATABLE :: wv_real(:),state_real(:),wv_real_all(:,:),state_real_tmp(:)
  REAL(kind=DP), ALLOCATABLE :: state_real_tmp2(:),state_real2(:)
  REAL(kind=DP), ALLOCATABLE :: omat(:,:)
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:),tmp_mat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: omat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: hmat(:,:)
  REAL(kind=DP), ALLOCATABLE :: e_fake(:), vec_fake(:,:)
  REAL(kind=DP), ALLOCATABLE :: gap(:)
  REAL(kind=DP), ALLOCATABLE :: hmat_i(:,:),hmat_o(:,:), omat_i(:,:)
  REAL(kind=DP), ALLOCATABLE :: ovec(:)
  REAL(kind=DP), ALLOCATABLE :: g2kint(:)
  INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
  INTEGER, ALLOCATABLE :: isuppz(:)
  INTEGER, ALLOCATABLE :: iclustr(:)
  INTEGER, ALLOCATABLE :: igkt(:)

  REAL(kind=DP):: sca1,sca2

  LOGICAL :: l_test=.false.!if true test the completness of the basis
  LOGICAL :: l_dsyevr=.true.!if true uses dsyevr instead of dsyev
  LOGICAL :: l_diago_cg=.true.!if true uses diago_cg instead of  dsyevr ATTENZIONE
  LOGICAL :: exst
  LOGICAL :: l_dsygvx=.false.!if .true. uses serial dsygvx instead of parallel diago_cg_g
  LOGICAL :: l_gramsc=.true.!if true orthonormalization through gram-schimdt
  LOGICAL :: l_diago_para=.true.!if true uses parallel diago_cg
  LOGICAL :: l_fft_custom=.false.

  INTEGER :: ig,ip, ii, iv, jj, iw, ir, is
  INTEGER :: num_fc!number of fake conduction states
  INTEGER :: lwork,info,liwork
  INTEGER :: n_out
  INTEGER :: fcw_number_old
  INTEGER :: l_blk,nbegin,nend
  INTEGER :: max_state
  INTEGER :: iunfcw
  INTEGER :: nsize
  INTEGER :: nbegin_loc,nend_loc,nsize_loc
  INTEGER :: n_found_state
!variables for scalapack
  INTEGER :: num_fc_r,num_fc_c,num_fc_dimr,num_fc_dimc
  INTEGER :: m,nz,icrow,iccol,iproc,ilrow,ilcol
  INTEGER :: desc_a(9),desc_b(9),desc_c(9)
  INTEGER :: n_computed
  INTEGER :: num_built!number of states already built
  INTEGER :: num_out
  INTEGER ::  kilobytes
  INTEGER ::  kb_old, kb_new

  INTEGER, EXTERNAL :: indxg2p,indxg2l

  TYPE(optimal_options) :: options



  INTEGER :: bufferx,fcw_numberx,fcw_number_oldx, fcw_numberx_tmp


  LOGICAL :: l_restart0!if true restart is enabled
  INTEGER :: iunrestart0, iv_start,iunfsr
  REAL(kind=DP), ALLOCATABLE :: state_fc_r(:,:,:)
  INTEGER :: num_nbndv_max, num_fc_spin
  INTEGER :: num_fc_eff(2),num_fc_eff_max
 
  LOGICAL :: l_do_optimal
  INTEGER :: iun_oap

  COMPLEX(kind=DP), ALLOCATABLE :: cbuf(:,:)


  TYPE(fft_cus) :: fc

!determine bufferx,fcw_numberx
  bufferx=num_nbndv(1)*300/4
  bufferx=max(1000,bufferx)
!ONLY FOR PROJECT ON JADE
  bufferx=5000
  fcw_numberx=bufferx
  fcw_number_oldx=bufferx
  fcw_numberx_tmp=bufferx

!generate fake conduction states
!!determine number of states

!generate custom in grid in case can be equal to norm-conserving grid

  fc%ecutt=ecutwfc
  fc%dual_t=dual_pb

  write(stdout,*) 'Call initialize_fft_custom'
  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory0', kilobytes
  FLUSH(stdout)

  call initialize_fft_custom(fc)

  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory0.0', kilobytes
  FLUSH(stdout)

  ! this is for compatibility

  allocate( igkt( fc%npwt ) )
  do ig=1,fc%npwt
     igkt(ig)=ig
  enddo

  allocate( evc_g( fc%ngmt_g ) )

  !plane waves basis set

  !state_fc are first obtained on the ordering of the normconserving grid

  g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
       ( g(2,igk_k(1:npw,1)) )**2 + &
       ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2
  
  num_fc=0
  do ig=1,npw
     if(g2kin(ig) <= cutoff) num_fc=num_fc+1
  enddo
  call start_clock('mpsum')
  call mp_sum(num_fc,world_comm)
  call stop_clock('mpsum')
  num_fc=(num_fc-1)*2+1

  if(.not.l_cond) then
     if(.not.l_frac) then
        num_fc_eff(1:2)=num_fc
        num_fc_eff_max=num_fc
     else
        num_fc_eff(1:2)=num_fc+num_nbndv(1:2)-num_nbndv_min(1:2)
        num_fc_eff_max=max(num_fc_eff(1),num_fc_eff(2))
     endif
  else
      if(.not.l_frac) then
         num_fc_eff(1:2)=num_fc+num_nbnds-num_nbndv(1:2)
         num_fc_eff_max=num_fc+num_nbnds-min(num_nbndv(1),num_nbndv(2))
     else
        num_fc_eff(1:2)=num_fc+num_nbndv(1:2)-num_nbndv_min(1:2)+num_nbnds-num_nbndv(1:2)
        num_fc_eff_max=max(num_fc_eff(1),num_fc_eff(2))
     endif
  endif
  allocate( state_fc( npw, num_fc_eff_max, nspin ) )

  state_fc(:,:,:)=(0.d0,0.d0)
  write(stdout,*) "Number of projected orthonormalized plane waves:", num_fc
  CALL memstat( kilobytes )
  if(l_verbose)  write(stdout,*) 'memory0.1', kilobytes, ' new kb = ', &
       &(SIZE( state_fc )*16 + SIZE( evc_g )*16 + SIZE( igkt )*4)/1024
  FLUSH(stdout)
  
  ii=0
  do ip=0,nproc-1
     if(mpime==ip) then
        do ig=gstart,npw
           if(g2kin(ig) <= cutoff) then
              ii=ii+1
              state_fc(ig,ii,1)=cmplx(dsqrt(0.5d0),0.d0)
              ii=ii+1
              state_fc(ig,ii,1)=cmplx(0.d0,dsqrt(0.5d0))
           endif
        enddo
        if(gstart==2) then
           ii=ii+1
           state_fc(1,ii,1)=(1.d0,0.d0)
        endif
     else
        ii=0
     endif
     call start_clock('mpsum')
     call mp_sum(ii,world_comm)
     call stop_clock('mpsum')
  enddo

  if(ii/=num_fc) then 
     write(stdout,*) 'ERRORE FAKE CONDUCTION',ii
     FLUSH(stdout)
     stop
    return
  endif
  if(l_verbose)  write(stdout,*) 'FAKE1'
  FLUSH(stdout)
  
  if(nspin==2) state_fc(:,1:num_fc,2)=state_fc(:,1:num_fc,1)
  do is=1,nspin

!!project out of valence space
     do ii=1,num_fc
        evc(1:npw,1:num_nbndv(is))=ks_wfcs(1:npw,1:num_nbndv(is),is)!for calling pc_operator
        call pc_operator(state_fc(:,ii,is),is,l_cond)
      enddo
  enddo


!!add partially occupied states
  if(l_frac) then
     do is=1,nspin
        do ii=num_nbndv_min(is)+1,num_nbndv(is)
           state_fc(1:npw,num_fc+ii-num_nbndv_min(is),is)=ks_wfcs_diag(1:npw,ii,is)
        enddo
     enddo
  endif


!!add conduction states if required
  if(l_cond) then
     if(.not.l_frac) then
        do is=1,nspin
           do ii=num_nbndv(is)+1,num_nbnds
              state_fc(1:npw,num_fc+ii-num_nbndv(is),is)=ks_wfcs_diag(1:npw,ii,is)
           enddo
        enddo
     else
        do is=1,nspin
           do ii=num_nbndv(is)+1,num_nbnds
              state_fc(1:npw,num_fc+num_nbndv(is)-num_nbndv_min(is)+ii-num_nbndv(is),is)=ks_wfcs_diag(1:npw,ii,is)
           enddo
        enddo
     endif
  endif
!orthonormalize fake_conduction states 

 

     
!for the moment finds all the first fcw_fast_n eigenstates
     
     if(l_verbose)  write(stdout,*) 'CASE ORTHONORMALIZATION ONLY'
     FLUSH(stdout)
     
!if required orthonormalize the projected plane_waves or read from disk

     l_do_optimal=.false.
     inquire(file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', exist = exst)
     if(.not. exst) then
        l_do_optimal=.true.
     else
        iunrestart0 =  find_free_unit()
        open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='old')
        read(iunrestart0,*) iv_start
        close(iunrestart0)
        if(iv_start<1 ) l_do_optimal=.true.
     endif
           


     if(l_do_optimal) then
        if(l_verbose) write(stdout,*) 'Call optimal driver'
        FLUSH(stdout)
        options%l_complete=.true.
        options%idiago=0
        call start_clock('fc_optimal')
        do is=1,nspin
           call optimal_driver(num_fc_eff(is),state_fc(1,1,is),npw,options,num_out, info)
        enddo
        call stop_clock('fc_optimal')  

!read orthonormalized projected plane-waves from disk

     endif
        CALL memstat( kilobytes )
        if(l_verbose)  write(stdout,*) 'memory0.3', kilobytes
        FLUSH(stdout)
        
       
     
  !now state_fc are put on the ordering of the redueced grid, if required
  allocate(state_fc_t(fc%npwt,num_fc_eff_max,nspin))
 
  if(l_do_optimal) then
     if(fc%dual_t==4.d0) then
        do is=1,nspin
           state_fc_t(1:fc%npwt,1:num_fc_eff(is),is)=state_fc(1:fc%npwt,1:num_fc_eff(is),is)
        enddo
     else
        call start_clock('fc_merge')
        do is=1,nspin
           call reorderwfp (num_fc_eff(is),npw, fc%npwt,state_fc(:,:,is),state_fc_t(:,:,is), &
              &npw,fc%npwt, ig_l2g,fc%ig_l2gt, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

     
        !   do ii=1,num_fc_eff(is)
        !      call mergewf(state_fc(:,ii,is),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
        !      call splitwf(state_fc_t(:,ii,is),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
        !   enddo
        enddo
        call stop_clock('fc_merge')
     endif
     iun_oap = find_free_unit()
     CALL diropn( iun_oap, 'oap', fc%npwt*2, exst )
     do ii=1,num_fc_eff(1)
        CALL davcio( state_fc_t(:,ii,1), 2*fc%npwt, iun_oap, ii, 1 )
     enddo
     close(iun_oap)

  else
     if(l_verbose) write(stdout,*) 'Read OAP from disk'
     FLUSH(stdout)
     iun_oap = find_free_unit()
     CALL diropn( iun_oap, 'oap', fc%npwt*2, exst )
     do ii=1,num_fc_eff(1)
        CALL davcio( state_fc_t(:,ii,1), 2*fc%npwt, iun_oap, ii, -1 )
     enddo
     close(iun_oap)


  endif
  deallocate(state_fc)
 
  allocate(state_fc_r(fc%nrxxt,num_fc_eff_max,nspin))
  do is=1,nspin
     do ii=1,num_fc_eff(is),2
        psic(:)=(0.d0,0.d0)
        if(ii==num_fc_eff(is)) then
           psic(fc%nlt(1:fc%npwt))  = state_fc_t(1:fc%npwt,ii,is)
           psic(fc%nltm(1:fc%npwt)) = CONJG( state_fc_t(1:fc%npwt,ii,is) )
        else
           psic(fc%nlt(1:fc%npwt))=state_fc_t(1:fc%npwt,ii,is)+(0.d0,1.d0)*state_fc_t(1:fc%npwt,ii+1,is)
           psic(fc%nltm(1:fc%npwt)) = CONJG( state_fc_t(1:fc%npwt,ii,is) )+(0.d0,1.d0)*CONJG( state_fc_t(1:fc%npwt,ii+1,is) )
        endif
        CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
        state_fc_r(1:fc%nrxxt,ii,is)= DBLE(psic(1:fc%nrxxt))
        if(ii/=num_fc_eff(is)) state_fc_r(1:fc%nrxxt,ii+1,is)= DIMAG(psic(1:fc%nrxxt))
     enddo
  enddo
  deallocate(state_fc_t) 
   

  CALL memstat( kilobytes )
  if(l_verbose)  write(stdout,*) 'memory0.4', kilobytes, ' new kb = ', (SIZE( state_fc_t ))/64
  FLUSH(stdout)

!set maximum number of valence states for both spin channels
  if(nspin==1) then
     num_nbndv_max=num_nbndv(1)
  else
     num_nbndv_max=max(num_nbndv(1),num_nbndv(2))
  endif

!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbndv_max,nspin))
  if(fc%dual_t==4.d0) then
     evc_t(1:fc%npwt,1:num_nbndv_max,1:nspin)=ks_wfcs(1:fc%npwt,1:num_nbndv_max,1:nspin)
  else
     call start_clock('fc_merge')
     do is=1,nspin
        call reorderwfp (num_nbndv(is),npw, fc%npwt,ks_wfcs(:,:,is),evc_t(:,:,is), &
             &npw,fc%npwt, ig_l2g,fc%ig_l2gt, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

      !  do iv=1,num_nbndv(is)
      !     call mergewf(ks_wfcs(:,iv,is),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
      !     call splitwf(evc_t(:,iv,is),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
      !  enddo
     enddo
     call stop_clock('fc_merge')
  endif

  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory0.5', kilobytes, ' new kb = ', (SIZE( evc_t ))/64
  FLUSH(stdout)

!cycle on v 
!! product in real space with wannier
!! orthonormalize and take N most important
!! gram-schmidt like 

!calculate D matrix 

!  l_blk= (num_fc)/nproc
!  if(l_blk*nproc < (num_fc)) l_blk = l_blk+1
!  nbegin=mpime*l_blk+1
!  nend=nbegin+l_blk-1
!  if(nend > num_fc) nend=num_fc
!  nsize=nend-nbegin+1

!check for restart
 if(ionode) then

     inquire(file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', exist = exst)
     if(.not. exst) then
        iv_start=1
     else
        iunrestart0 =  find_free_unit()
        open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='old')
        read(iunrestart0,*) iv_start
        read(iunrestart0,*) fcw_number
        read(iunrestart0,*) fcw_numberx
        close(iunrestart0)
        if(iv_start<1 ) then
           iv_start=1
        else
           iv_start=iv_start+1
        endif
     endif
  endif
  call mp_bcast(iv_start,ionode_id,world_comm)

  if(iv_start/=1) then
     call mp_bcast(fcw_number,ionode_id,world_comm)
     call mp_bcast(fcw_numberx,ionode_id,world_comm)
     fcw_number_oldx=fcw_numberx
     fcw_numberx_tmp=fcw_numberx
     fcw_number_old=fcw_number
     allocate(fcw_state_r(fc%nrxxt,fcw_numberx))
     allocate(fcw_state_old_r(fc%nrxxt,fcw_numberx))
 
 !read them from file
     iunfsr = find_free_unit()
 
     CALL diropn( iunfsr, 'fsr', fc%nrxxt, exst )
     do ii=1,fcw_number
        CALL davcio( fcw_state_r(1,ii), fc%nrxxt, iunfsr, ii, -1 )
        fcw_state_old_r(1:fc%nrxxt,ii)=fcw_state_r(1:fc%nrxxt,ii)
     enddo
 
     close(iunfsr)

  endif

  allocate (wv_real(fc%nrxxt),state_real(fc%nrxxt),state_real2(fc%nrxxt),state_g(fc%nrxxt,num_fc_eff_max*nspin))


  
  FIRST_LOOP: do iv=iv_start,num_nbndv_max
     call start_clock('fc_loop')
     write(stdout,*) 'FK state:', iv,fc%nrxxt,fc%npwt,num_fc
     CALL memstat( kilobytes )
     if(l_verbose)  write(stdout,*) 'memory1', kilobytes
     FLUSH(stdout)
     

     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory2', kilobytes, ' new kb = ', &
                     (SIZE(wv_real)+SIZE(state_real)+SIZE(state_real2)+SIZE(state_g))/128
     FLUSH(stdout)

     num_fc_spin=0
     do is=1,nspin

        if(iv<= num_nbndv(is)) then
           psic(:)=(0.d0,0.d0)
           psic(fc%nlt(igkt(1:fc%npwt)))  = evc_t(1:fc%npwt,iv,is)
           psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( evc_t(1:fc%npwt,iv,is) )

           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory3', kilobytes
           FLUSH(stdout)

           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
           CALL memstat( kilobytes )
           if(l_verbose) write(stdout,*) 'memory4', kilobytes
           FLUSH(stdout)

     
           wv_real(:)= DBLE(psic(:))
           if(l_verbose) then
              if(fc%gstart_t==2) write(stdout,*) 'FAKE modulus valence', iv, evc_t(1,iv,is)
           endif
!loop on fake conduction states

           FLUSH(stdout)
           do ii=1,num_fc_eff(is)
              state_g(1:fc%nrxxt, ii+num_fc_spin)=state_fc_r(1:fc%nrxxt,ii,is)*wv_real(1:fc%nrxxt)
           enddo

           num_fc_spin=num_fc_spin+num_fc_eff(is)
        endif
     enddo!on spin


     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory5', kilobytes
     FLUSH(stdout)

!if not first valence wannier project the products out of partial orthonormal basis

     if(l_verbose) write(stdout,*) 'Start Projection part'
     FLUSH(stdout)

     if(iv/=1) then
        ! build overlap matrix
        
        if(iv==2 .or. (iv_start/=1.and. iv==iv_start)) then
           allocate(tmp_mat2(fcw_numberx_tmp,num_fc_eff_max*nspin))
        else
           if(fcw_number>fcw_numberx_tmp) then
              deallocate(tmp_mat2)
              fcw_numberx_tmp=fcw_numberx_tmp+bufferx
              allocate(tmp_mat2(fcw_numberx_tmp,num_fc_eff_max*nspin))
              if(l_verbose) write(stdout,*) 'Updated dimension of tmp_mat2', fcw_numberx_tmp
           endif
        endif
        
        call start_clock('fc_dgemm')
        call dgemm('T','N',fcw_number,num_fc_spin,fc%nrxxt,1.d0,fcw_state_r,fc%nrxxt,&
             &state_g,fc%nrxxt,0.d0,tmp_mat2,fcw_numberx_tmp)
        call stop_clock('fc_dgemm')
        do ii=1,num_fc_spin
           call start_clock('mpsum')
           call mp_sum(tmp_mat2(1:fcw_number,ii),world_comm)
           call stop_clock('mpsum')
           tmp_mat2(1:fcw_number,ii)=tmp_mat2(1:fcw_number,ii)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
        enddo
        !call mp_sum(tmp_mat2,world_comm)
        
        call start_clock('fc_dgemm')
        call dgemm('N','N',fc%nrxxt, num_fc_spin,fcw_number,-1.d0,fcw_state_r,fc%nrxxt,tmp_mat2,&
        &fcw_numberx_tmp,1.d0,state_g,fc%nrxxt)
        call stop_clock('fc_dgemm')
        if(iv==num_nbndv_max) deallocate(tmp_mat2)
     endif

     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory6', kilobytes
     if(l_verbose) write(stdout,*) 'End Projection part'
     FLUSH(stdout)

!calculate overlap matrix

     if(l_verbose) write(stdout,*) 'FK2'!ATTENZIONE
     FLUSH(stdout)

     max_state=max(300,num_fc/20)
     if(max_state > num_fc) max_state=num_fc/2

     l_blk= (num_fc_spin)/nproc
     if(l_blk*nproc < (num_fc_spin)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > num_fc_spin) nend=num_fc
     nsize=nend-nbegin+1



     

     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory6.6', kilobytes
     FLUSH(stdout)

        !uses iterative algorithm
        !allocate max number of new states

        !deallocate(wv_real,state_real,state_real2)
        !gram shimdt like
     allocate(ovec(num_fc_spin))
     num_built=0


     do ii=1,num_fc_spin
           
        if(num_built>0) then
           call start_clock('fc_dgemm')
           call dgemm('T','N',num_built,1,fc%nrxxt,1.d0,state_g,fc%nrxxt,state_g(1,ii),fc%nrxxt,0.d0,ovec,num_fc_spin)
           call stop_clock('fc_dgemm')
           call start_clock('mpsum')
           call mp_sum(ovec(1:num_built),world_comm)
           call stop_clock('mpsum')
           ovec(1:num_built)=ovec(1:num_built)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
           call start_clock('fc_dgemm')
           call dgemm('T','N',1,1,fc%nrxxt,1.d0,state_g(1,ii),fc%nrxxt,state_g(1,ii),fc%nrxxt,0.d0,sca2,1)
           call stop_clock('fc_dgemm')
           call start_clock('mpsum')
           call mp_sum(sca2,world_comm)
           call stop_clock('mpsum')
           sca2=sca2/dble(fc%nr1t*fc%nr2t*fc%nr3t)
           sca1=0.d0
           do jj=1,num_built
              sca1=sca1+ovec(jj)**2.d0
           enddo
           if(abs(sca2-sca1) >= s_cutoff) then
              if(num_built+1 /= ii) state_g(1:fc%nrxxt,num_built+1)=state_g(1:fc%nrxxt,ii)
              call start_clock('fc_dgemm')
              call dgemm('N','N',fc%nrxxt,1,num_built,-1.d0,state_g,fc%nrxxt,&
                   &ovec,num_fc_spin,1.d0,state_g(1,num_built+1),fc%nrxxt)
              call stop_clock('fc_dgemm')
              num_built=num_built+1
              call start_clock('fc_dgemm')
              call dgemm('T','N',1,1,fc%nrxxt,1.d0,state_g(1,num_built),&
                   &fc%nrxxt,state_g(1,num_built),fc%nrxxt,0.d0,ovec(num_built),1)
              call stop_clock('fc_dgemm')
              call start_clock('mpsum')
              call mp_sum(ovec(num_built),world_comm)
              call stop_clock('mpsum')
              ovec(num_built)=ovec(num_built)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
              ovec(num_built)=1.d0/dsqrt(ovec(num_built))
              state_g(1:fc%nrxxt,num_built)=state_g(1:fc%nrxxt,num_built)*ovec(num_built)
              
           endif
        else
           call start_clock('fc_dgemm')
           call dgemm('T','N',1,1,fc%nrxxt,1.d0,state_g(1,ii),&
                &fc%nrxxt,state_g(1,ii),fc%nrxxt,0.d0,ovec(1),1)
           call stop_clock('fc_dgemm')
           call start_clock('mpsum')
           call mp_sum(ovec(1),world_comm)
           call stop_clock('mpsum')
           ovec(1)=ovec(1)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
           if(ovec(1) >= s_cutoff) then
              num_built=1
              ovec(num_built)=1.d0/dsqrt(ovec(num_built))
              state_g(1:fc%nrxxt,num_built)=state_g(1:fc%nrxxt,ii)*ovec(num_built)                
           endif
        endif
     enddo

     deallocate(ovec)   
        
     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory6.7', kilobytes
     write(stdout,*) 'FK GS', num_built
     FLUSH(stdout)


        !opportune basis
        
     if( iv == 1 ) then
        fcw_number=num_built
        allocate( fcw_state_r( fc%nrxxt, fcw_numberx ) )
        allocate( fcw_state_old_r( fc%nrxxt, fcw_number_oldx ) )
        if(num_built>0) then
           fcw_state_r(1:fc%nrxxt,1:num_built) = state_g(1:fc%nrxxt,1:num_built)
           iunfsr = find_free_unit()
           CALL diropn( iunfsr, 'fsr', fc%nrxxt, exst )
           do ii=1,num_built
              CALL davcio( fcw_state_r(1,ii), fc%nrxxt, iunfsr, ii, 1 )
           enddo
           close(iunfsr)
        endif
     else
        if(fcw_number+num_built>fcw_number_oldx) then
           fcw_number_oldx=fcw_number_oldx+bufferx
           deallocate(fcw_state_old_r)
           allocate( fcw_state_old_r( fc%nrxxt, fcw_number_oldx ) )
        endif
           
        fcw_state_old_r(1:fc%nrxxt,1:fcw_number) = fcw_state_r(1:fc%nrxxt,1:fcw_number)
        if(fcw_number+num_built>fcw_numberx) then
           fcw_numberx=fcw_numberx+bufferx
           deallocate(fcw_state_r)
           allocate(fcw_state_r(fc%nrxxt,fcw_numberx))
        endif
        
        fcw_state_r(1:fc%nrxxt,1:fcw_number)=fcw_state_old_r(1:fc%nrxxt,1:fcw_number)
        if(num_built> 0) then
           fcw_state_r(1:fc%nrxxt,fcw_number+1:fcw_number+num_built)=state_g(1:fc%nrxxt,1:num_built)
           CALL diropn( iunfsr, 'fsr', fc%nrxxt, exst )
           do ii=1,num_built
              CALL davcio( fcw_state_r(1,ii+fcw_number), fc%nrxxt, iunfsr, ii+fcw_number, 1 )
           enddo
           close(iunfsr)
        endif
        fcw_number=fcw_number+num_built
        
     endif
 

        

     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory6.8', kilobytes
     FLUSH(stdout)
     iunrestart0 =  find_free_unit()
     open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='unknown')
     write(iunrestart0,*) iv
     write(iunrestart0,*) fcw_number
     write(iunrestart0,*) fcw_numberx
     close(iunrestart0)
     call stop_clock('fc_loop')
  end do FIRST_LOOP
 
  deallocate( wv_real, state_real, state_real2, state_g )
 
  if(l_verbose) write(stdout,*) 'FK8'
  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory7', kilobytes
  FLUSH(stdout)
  
  if(num_nbndv(1)/=1 ) deallocate(fcw_state_old_r)
 
! calculate D matrix distributed among processors

  if(fcw_number < nproc) then
     write(stdout,*) 'too many processors'
     stop
  endif

  l_blk= (fcw_number)/nproc
  if(l_blk*nproc < (fcw_number)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > fcw_number) nend=fcw_number
  nsize=nend-nbegin+1

  allocate(fcw_mat(fcw_number,l_blk))
  fcw_mat(:,:)=0.d0

  if(l_verbose) write(stdout,*) 'FK9'
  FLUSH(stdout)

  allocate(wv_real_all(fc%nrxxt,num_nbndv_max))

  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory8', kilobytes
  
  allocate(state_real(fc%nrxxt),state_real_tmp(fc%nrxxt),state_real_tmp2(fc%nrxxt))
  
  allocate(state_g(fc%nrxxt,num_nbndv_max))
  
  allocate(tmp_mat(fcw_number,num_nbndv_max))

  write(stdout,*) 'Calculate FK matrix'
  FLUSH(stdout)

  do is=1,nspin
     do iv=1,num_nbndv(is)
      
     
        psic(:)=(0.d0,0.d0)
        psic(fc%nlt(igkt(1:fc%npwt)))  = evc_t(1:fc%npwt,iv,is)
        psic(fc%nltm(igkt(1:fc%npwt))) = CONJG( evc_t(1:fc%npwt,iv,is) )
           
        CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
        wv_real_all(1:fc%nrxxt,iv)= DBLE(psic(1:fc%nrxxt))
!check for modulus
        sca1=0.d0
        do ir=1,fc%nrxxt
           sca1=sca1+wv_real_all(ir,iv)**2.d0
        enddo
        call start_clock('mpsum')
        call mp_sum(sca1,world_comm)
        call stop_clock('mpsum')
        if(l_verbose) write(stdout,*) 'Modulus:',fc%nrxxt,fc%nr1t*fc%nr2t*fc%nr3t, sca1/(dble(fc%nr1t*fc%nr2t*fc%nr3t))
        
  
       
     enddo
   
!loop on fake conduction states
 


           
     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory9', kilobytes
  
     do ii=1,num_fc_eff(is)
        do iv=1,num_nbndv(is)
           state_g(1:fc%nrxxt,iv)=state_fc_r(1:fc%nrxxt,ii,is)*wv_real_all(1:fc%nrxxt,iv)
        enddo


!do products with fcw states

        call start_clock('fc_dgemm')
        call dgemm('T','N',fcw_number,num_nbndv(is),fc%nrxxt,1.d0,fcw_state_r,fc%nrxxt,state_g,fc%nrxxt,0.d0,tmp_mat,fcw_number)
        call stop_clock('fc_dgemm')
        call start_clock('mpsum')
        call mp_sum(tmp_mat,world_comm)
        call stop_clock('mpsum')
        tmp_mat=tmp_mat/dble(fc%nr1t*fc%nr2t*fc%nr3t)

        if(l_frac) then
           if(ii<=num_fc) then
              do iv=1,num_nbndv(is)
                 if(nspin==1) then
                    sca1=dsqrt(abs(wg(iv,is))/2.d0)
                 else
                    sca1=dsqrt(abs(wg(iv,is)))
                 endif
                 tmp_mat(1:fcw_number,iv)=tmp_mat(1:fcw_number,iv)*sca1
              enddo
           else
              do iv=1,num_nbndv(is)
                 if(nspin==1) then
                    sca1=dsqrt(abs(wg(ii-num_fc+num_nbndv_min(is),is)-wg(iv,is))/2.d0)
                 else
                    sca1=dsqrt(abs(wg(ii-num_fc+num_nbndv_min(is),is)-wg(iv,is)))
                 endif
                 tmp_mat(1:fcw_number,iv)=tmp_mat(1:fcw_number,iv)*sca1
              enddo
           endif
        endif

        CALL memstat( kilobytes )
        if(l_verbose) write(stdout,*) 'memory10', kilobytes
        if(l_verbose) write(stdout,*) 'TOTAL NUMBER OF FCW STATES:', fcw_number,ii,dfftp%nnr,fc%nrxxt,wg(1,is)
        FLUSH(stdout)
            
!calculate contribution to D matrix
        if(nsize>0) then
           call start_clock('fc_dgemm')
           call dgemm('N','T',fcw_number,nend-nbegin+1,num_nbndv(is),1.d0,tmp_mat,fcw_number,&
                &tmp_mat(nbegin:nend,1:num_nbndv(is)),nend-nbegin+1,1.d0,fcw_mat,fcw_number)
           call stop_clock('fc_dgemm')
        endif
       
        enddo
     enddo!on spin

!if required put fcw_state on normconserving ordering

     allocate(fcw_state_n(npw,fcw_number))
     allocate(fcw_state(fc%npwt,fcw_number))!ATTENZIONE the use of memory could be reduced
     psic=0.d0
     do ii=1,fcw_number,2
        if(ii==fcw_number) then
           psic(1:fc%nrxxt)=cmplx(fcw_state_r(1:fc%nrxxt,ii),0.d0)
        else
           psic(1:fc%nrxxt)=cmplx(fcw_state_r(1:fc%nrxxt,ii),fcw_state_r(1:fc%nrxxt,ii+1))
        endif
        CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
        if(ii==fcw_number) then
           fcw_state(1:fc%npwt, ii) = psic(fc%nlt(1:fc%npwt))
           if(fc%gstart_t==2) fcw_state(1,ii)=(0.d0,0.d0)
        else
           fcw_state(1:fc%npwt, ii)= 0.5d0*(psic(fc%nlt(igkt(1:fc%npwt)))+conjg( psic(fc%nltm(igkt(1:fc%npwt)))))
           fcw_state(1:fc%npwt, ii+1)= (0.d0,-0.5d0)*(psic(fc%nlt(igkt(1:fc%npwt))) - conjg(psic(fc%nltm(igkt(1:fc%npwt)))))
           if(fc%gstart_t==2) fcw_state(1,ii)=(0.d0,0.d0)
           if(fc%gstart_t==2) fcw_state(1,ii+1)=(0.d0,0.d0)
        endif

     enddo
     if(fc%dual_t==4.d0) then
        fcw_state_n(1:fc%npwt,1:fcw_number)=fcw_state(1:fc%npwt,1:fcw_number)
     else
        call start_clock('fc_merge')
        call reorderwfp (fcw_number,fc%npwt, npw,fcw_state,fcw_state_n, &
             &fc%npwt,npw, fc%ig_l2gt,ig_l2g, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )
        !   call mergewf(fcw_state(:,1),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
        !  call splitwf(fcw_state_n(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
        call stop_clock('fc_merge')
     endif
    
     CALL memstat( kilobytes )
     if(l_verbose) write(stdout,*) 'memory11', kilobytes
 

  !save on file

  iunfcw = find_free_unit()
  CALL diropn( iunfcw, 'fcw', npw*2, exst )
  do ii=1,fcw_number
     CALL davcio( fcw_state_n(1,ii), 2*npw, iunfcw, ii, 1 )
  enddo
  close(iunfcw)

  !write number of states

  if(ionode) then
     open(unit=iunfcw,file=trim(tmp_dir)//trim(prefix)//'.nfcws',status='unknown')
     write(iunfcw,*) fcw_number
     close(iunfcw)
  endif

  CALL diropn( iunfcw, 'fmat',fcw_number, exst )
  do ii=1,nsize
     CALL davcio( fcw_mat(1,ii), fcw_number, iunfcw, ii, 1 )
  enddo
  close(iunfcw)

  if(l_verbose) write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  call deallocate_fft_custom(fc)



  iunrestart0 =  find_free_unit()
  open( unit= iunrestart0, file=trim(tmp_dir)//trim(prefix)//'.restart_fk0_status', status='unknown')
  write(iunrestart0,*) -1
  write(iunrestart0,*) fcw_number
  write(iunrestart0,*) fcw_numberx
  close(iunrestart0)


  deallocate(wv_real_all)
  deallocate(fcw_state_r)
  
  
  if( allocated( state_fc_t ) ) deallocate( state_fc_t )
  deallocate(state_real,state_g,state_real_tmp,state_real_tmp2)
  deallocate(tmp_mat)
  if(allocated(e_fake)) deallocate(e_fake)
  deallocate(fcw_state_n)
  deallocate(evc_g,evc_t)


  if( allocated( state_fc ) ) deallocate( state_fc )
  if( allocated( state_g ) ) deallocate( state_g )
  if( allocated( fcw_state_old_r ) ) deallocate( fcw_state_old_r )
  if( allocated( h_state_fc ) ) deallocate( h_state_fc ) 
  if( allocated( evc_g ) ) deallocate( evc_g )
  if( allocated( evc_t ) ) deallocate( evc_t )
  if( allocated( state_fc_t ) ) deallocate( state_fc_t )
  if( allocated( state_g_t ) ) deallocate( state_g_t ) 
  if( allocated( fcw_state_n ) ) deallocate( fcw_state_n )
  if( allocated( wv_real ) ) deallocate( wv_real )
  if( allocated( state_real ) ) deallocate( state_real )
  if( allocated( wv_real_all ) ) deallocate( wv_real_all )
  if( allocated( state_real_tmp ) ) deallocate( state_real_tmp )
  if( allocated( state_real_tmp2 ) ) deallocate( state_real_tmp2 )
  if( allocated( state_real2 ) ) deallocate( state_real2 )
  if( allocated( omat ) ) deallocate( omat )
  if( allocated( eigen ) ) deallocate( eigen )
  if( allocated( work ) ) deallocate( work )
  if( allocated( tmp_mat ) ) deallocate( tmp_mat )
  if( allocated( omat2 ) ) deallocate( omat2 )
  if( allocated( hmat ) ) deallocate( hmat )
  if( allocated( e_fake ) ) deallocate( e_fake )
  if( allocated( vec_fake ) ) deallocate( vec_fake )
  if( allocated( gap ) ) deallocate( gap )
  if( allocated( hmat_i ) ) deallocate( hmat_i )
  if( allocated( hmat_o ) ) deallocate( hmat_o )
  if( allocated( omat_i ) ) deallocate( omat_i )
  if( allocated( ovec ) ) deallocate( ovec )
  if( allocated( g2kint ) ) deallocate( g2kint )
  !
  if( allocated( iwork ) )   deallocate( iwork )
  if( allocated( ifail ) )   deallocate( ifail )
  if( allocated( isuppz ) )  deallocate( isuppz )
  if( allocated( iclustr ) ) deallocate( iclustr )
  if( allocated( igkt ) )    deallocate( igkt )
  CALL memstat( kilobytes )
  if(l_verbose) write(stdout,*) 'memory12', kilobytes
  if(l_verbose) write(stdout,*) 'memory fcw_state = ',  SIZE( fcw_state ) / 64 , ' kb'
  if(l_verbose) write(stdout,*) 'memory fcw_mat   = ',  SIZE( fcw_mat   ) / 64 , ' kb'
  FLUSH(stdout)

  return
!NOT_TO_BE_INCLUDED_END
end subroutine fake_conduction_real







END MODULE fake_cond_mod
