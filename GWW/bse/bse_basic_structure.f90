MODULE bse_basic_structures
!this module describes the basis structures 
!which are obtained from the DFT and GW code
  USE kinds, ONLY : DP

  REAL(kind=DP), ALLOCATABLE :: vg_q(:) ! contains the elements V(G) of the Coulomb potential obtained upon integration over q

  COMPLEX(kind=DP), ALLOCATABLE :: u_trans(:,:,:)!unitarian transformation from bloch wfcs to wannier'

  TYPE wannier_o
! this structures contains the overlap between the wannier square modules
     integer :: numb_v ! number of valence bands for the two spin channels 
     real(kind=dp),dimension(:,:), pointer ::o(:,:) ! overlap matrix (numb_v*numb_v)   
  END TYPE

  TYPE ii_mat
     integer :: numb_v! number of valence bands for the two spin channels 
     integer :: np_max ! maximum number of overlapping wannier orbitals for a given v
     integer, dimension (:,:), pointer :: iimat(:,:) ! (np_max,numb_v) the rows of this matrix contain for each iv, 
                                                     !the set of jv for which   o_mat(iv,jv)>=s_bse
 
  END TYPE

  TYPE vww_prod
!this type contains the v*wv*wv' products needed for the exchange part of the
!direct interaction term of the excitonic Hamiltonian
     integer :: numb_v! number of valence bands for the two spin channels 
     integer :: npw ! number of plane wave per processor
     integer :: np_max ! maximum number of overlapping wannier orbitals for a given v
     complex(kind=dp), dimension (:,:,:), pointer :: vww(:,:,:) ! v*wv*wv' product in G space (npw,np_max,numb_v)

  END TYPE
  
  TYPE bse_z
! this type contains the z terms to build up the Wc term of the excitonic Hamiltonian
!z_beta_v_v'=(v*phi_beta)*wv*Wv'
     integer :: numb_v! number of valence bands for the two spin channels 
     integer :: np_max ! maximum number of overlapping wannier orbitals for a given v
     integer :: numw_prod ! dimension of the polarizability basis
     real(kind=dp), dimension (:,:,:), pointer :: z(:,:,:) ! v*phi_beta*wv*wv' product (numw_prod,np_max,numb_v)
  END TYPE

  TYPE v_state
! this type contains the valence states wavefunctions and single particle energies
!
  integer :: nspin  ! number of spin channels
  integer :: numb_v(2) ! number valence state
  integer :: npw ! number of g-vectors per processor
  real(kind=dp), dimension (:,:),pointer  :: esp(:,:) ! single particle energies (numb_v,nspin) 
  complex(kind=dp), dimension(:,:,:), pointer :: wfn(:,:,:) ! wave function in G space (npw,numb_v,nspin)  
  integer ::gstart

  END TYPE

  TYPE v_state_r
! this type contains the valence states wfns in real space on the dual grid 
  integer :: nspin  ! number of spin channels
  integer :: numb_v(2) ! number of valence states
  integer :: nrxxt ! number of r points per processor
  real(kind=dp), dimension(:,:,:), pointer :: wfnrt(:,:,:) ! wave function in r-spce (dual grid) (nrxxt,numb_v,nspin)

  END TYPE

  TYPE c_state
! this type contains the valence states wavefunctions and single particle energies
!
  integer :: nspin  ! number of spin channels
  integer :: numb_c ! number valence state
  integer :: npw ! number of g-vectors per processor
  real(kind=dp), dimension (:,:),pointer  :: esp(:) ! single particle energies (numb_c) 
  complex(kind=dp), dimension(:,:), pointer :: wfn(:,:) ! wave function in G space (npw,numb_c)  
  integer ::gstart

  END TYPE

  TYPE c_state_r
! this type contains the valence states wfns in real space on the dual grid 
  integer :: nspin  ! number of spin channels
  integer :: numb_c ! number of valence states
  integer :: nrxxt ! number of r points per processor
  real(kind=dp), dimension(:,:), pointer :: wfnrt(:,:) ! wave function in r-spce (dual grid) (nrxxt,numb_c)

  END TYPE

  CONTAINS

      subroutine initialize_v_state_r(v_wfnr)
      implicit none
      type(v_state_r) :: v_wfnr
      nullify(v_wfnr%wfnrt)
      return
      end subroutine

      subroutine initialize_v_state(v_wfn)
      implicit none
      type(v_state) :: v_wfn
      nullify(v_wfn%wfn)
      nullify(v_wfn%esp)
      return
      end subroutine

      subroutine initialize_c_state_r(c_wfnr)
      implicit none
      type(c_state_r) :: c_wfnr
      nullify(c_wfnr%wfnrt)
      return
      end subroutine

      subroutine initialize_c_state(c_wfn)
      implicit none
      type(c_state) :: c_wfn
      nullify(c_wfn%wfn)
      nullify(c_wfn%esp)
      return
      end subroutine
  
      subroutine initialize_wannier_o(o)
      implicit none
      type(wannier_o) :: o
      nullify(o%o)
      return
      end subroutine
  
      subroutine initialize_imat(iimat)
      implicit none
      type(ii_mat) :: iimat
      nullify(iimat%iimat)
      return
      end subroutine

      subroutine initialize_vww_prod(vww)
      implicit none
      type(vww_prod) :: vww
      nullify(vww%vww)
      return
      end subroutine

      subroutine initialize_bse_z(z)
      implicit none
      type(bse_z) :: z
      nullify(z%z)
      return
      end subroutine

      subroutine free_v_state_r(v_wfnr)
      implicit none
      type(v_state_r) :: v_wfnr
      if(associated(v_wfnr%wfnrt)) deallocate (v_wfnr%wfnrt)
      nullify(v_wfnr%wfnrt)
      return
      end subroutine

      subroutine free_v_state(v_wfn)
      implicit none
      type(v_state) :: v_wfn
      if(associated(v_wfn%wfn)) deallocate (v_wfn%wfn)
      nullify(v_wfn%wfn)
      if(associated(v_wfn%esp)) deallocate (v_wfn%esp)
      nullify(v_wfn%esp)
      return
      end subroutine

      subroutine free_c_state_r(c_wfnr)
      implicit none
      type(c_state_r) :: c_wfnr
      if(associated(c_wfnr%wfnrt)) deallocate (c_wfnr%wfnrt)
      nullify(c_wfnr%wfnrt)
      return
      end subroutine

      subroutine free_c_state(c_wfn)
      implicit none
      type(c_state) :: c_wfn
      if(associated(c_wfn%wfn)) deallocate (c_wfn%wfn)
      nullify(c_wfn%wfn)
      if(associated(c_wfn%esp)) deallocate (c_wfn%esp)
      nullify(c_wfn%esp)
      return
      end subroutine
  
      subroutine free_wannier_o(o)
      implicit none
      type(wannier_o) :: o
      if(associated(o%o)) deallocate (o%o)
      nullify(o%o)
      return
      end subroutine
  
      subroutine free_imat(iimat)
      implicit none
      type(ii_mat) :: iimat
      if(associated(iimat%iimat)) deallocate (iimat%iimat)
      nullify(iimat%iimat)
      return
      end subroutine

      subroutine free_vww_prod(vww)
      implicit none
      type(vww_prod) :: vww
      if(associated(vww%vww)) deallocate (vww%vww)
      nullify(vww%vww)
      return
      end subroutine

      subroutine free_bse_z(z)
      implicit none
      type(bse_z) :: z
      if(associated(z%z)) deallocate (z%z)
      nullify(z%z)
      return
      end subroutine

      subroutine make_v_state(numb_v,v)
      use io_global, ONLY : stdout, ionode 
      USE gvect,                 ONLY : gstart
      USE lsda_mod,              ONLY : nspin
      use wavefunctions_module,  ONLY : evc
      use io_files,  ONLY : prefix, iunwfc, tmp_dir
      USE io_files, ONLY: nwordwfc
      USE wvfct,    ONLY : nbnd, npwx,npw,et
      use mp_world, ONLY : mpime
      USE mp,          ONLY :mp_barrier
      USE mp_world,             ONLY : world_comm

      implicit none

      type(v_state) :: v
      integer :: numb_v(2)
  
      integer :: is,ivmax,iv
      logical :: debug

      debug=.false.
     
      call start_clock('make_v_state')
     
      if(debug) then
         write(*,*) 'make_v_state: in, mpime=',mpime
         ! debug MARGHE
         write(*,*) 'nbnd=', nbnd
         write(*,*) 'numb_v(1)=', numb_v(1)
      endif


      v%nspin=nspin
      v%numb_v(:)=numb_v(:)
      v%npw=npw
      v%gstart=gstart


      allocate( evc( npwx, nbnd ) )
  
      if (nspin==1) then
         ivmax= v%numb_v(1)
      else 
         ivmax=max(v%numb_v(1),v%numb_v(2))
      endif
      


      allocate( v%wfn(v%npw,ivmax,v%nspin))
      allocate( v%esp(ivmax,v%nspin))


      do is=1,nspin
         call davcio(evc,2*nwordwfc,iunwfc,is,-1)
         do iv=1,v%numb_v(is)
            v%wfn(1:v%npw,1:v%numb_v(is),is)=evc(1:v%npw,1:v%numb_v(is))
         enddo  
            v%esp(1:v%numb_v(is),is)=et(1:v%numb_v(is),is)
      enddo

      deallocate(evc)

      if(debug) then
         write(*,*) 'make_v_state: out, mpime=',mpime
      endif

      call mp_barrier( world_comm )
      call stop_clock('make_v_state')

      return
      end subroutine

      subroutine make_c_state(numb_v,c)
      use io_global, ONLY : stdout, ionode 
      USE gvect,                 ONLY : gstart
      USE lsda_mod,              ONLY : nspin
      use wavefunctions_module,  ONLY : evc
      use io_files,  ONLY : prefix, iunwfc, tmp_dir
      USE io_files, ONLY: nwordwfc
      USE wvfct,    ONLY : nbnd, npwx,npw,et
      use mp_world, ONLY : mpime
      USE mp,          ONLY :mp_barrier
      USE mp_world,             ONLY : world_comm

      implicit none

      type(c_state) :: c
      integer :: numb_v(2)
  
      integer :: is,ic
      logical :: debug

      debug=.false.
     
      call start_clock('make_c_state')
     
      if(debug) then
         write(*,*) 'make_c_state: in, mpime=',mpime
         ! debug MARGHE
         write(*,*) 'nbnd=', nbnd
         write(*,*) 'numb_v(1)=', numb_v(1)
      endif


      c%nspin=nspin
      c%numb_c=nbnd-numb_v(1)
      c%npw=npw
      c%gstart=gstart


      allocate( evc( npwx, nbnd ) )
  
!      if (nspin==1) then
!         ivmax= v%numb_v(1)
!      else 
!         ivmax=max(v%numb_v(1),v%numb_v(2))
!      endif
      


      allocate( c%wfn(c%npw,c%numb_c))
      allocate( c%esp(c%numb_c))


      do is=1,nspin
         call davcio(evc,2*nwordwfc,iunwfc,is,-1)
         do ic=1,c%numb_c
            c%wfn(1:c%npw,1:c%numb_c)=evc(1:c%npw,numb_v(is)+1:nbnd)
         enddo  
            c%esp(1:c%numb_c)=et(numb_v(is)+1:nbnd,is)
      enddo

      deallocate(evc)

      if(debug) then
         write(*,*) 'make_c_state: out, mpime=',mpime
      endif

      call mp_barrier( world_comm )
      call stop_clock('make_c_state')

      return
      end subroutine

      subroutine c_times_cstate(v,cstate_in,cstate_out)
      ! this subroutine multiplies each line ic of the c_state vector by the ic real component of the v vector
      use kinds, only:DP
      !use bse_wannier, only: qpe_imin,qpe_imax

      implicit none

      type(c_state),intent(in) :: cstate_in
      type(c_state),intent(out) :: cstate_out

      integer :: ib
      real(kind=DP) :: v(cstate_in%numb_c)

      do ib=1,cstate_in%numb_c
         cstate_out%wfn(1:cstate_out%npw,ib)=cmplx(v(ib),0.d0)* cstate_out%wfn(1:cstate_in%npw,ib)
      enddo
     
      return
      end subroutine

 
      subroutine v_wfng_to_wfnr(vwfng,fc,vwfnr)
     !this subroutine FFT the valence wfns to real space in the dual grid

      USE kinds, ONLY : DP
      USE fft_custom_gwl
      USE bse_wannier, ONLY : dual_bse 
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, nbndx
      USE io_global, ONLY : stdout, ionode, ionode_id
      USE mp_world, ONLY : mpime, nproc,world_comm
      USE mp_wave, ONLY : mergewf,splitwf
      USE mp,             ONLY : mp_sum
      USE gvect
      USE wavefunctions_module, ONLY :  psic

      implicit none

      type(v_state) vwfng
      type(v_state_r) vwfnr
      type(fft_cus) :: fc

      COMPLEX(kind=DP), allocatable :: vwfng_t(:,:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)

      integer :: ii,is
      integer ::ivmax  

      call start_clock('v_wfng_to_wfnr')
  
 
      if (vwfng%nspin==1) then
         ivmax= vwfng%numb_v(1)
      else 
         ivmax=max(vwfng%numb_v(1),vwfng%numb_v(2))
      endif


      allocate(vwfng_t(fc%npwt,ivmax,vwfng%nspin))

      vwfnr%nspin=vwfng%nspin
      vwfnr%nrxxt=fc%nrxxt       
      vwfnr%numb_v=vwfng%numb_v

      allocate(vwfnr%wfnrt(vwfnr%nrxxt,ivmax,vwfnr%nspin))
      
      allocate(evc_g(fc%ngmt_g ))

      if(fc%dual_t==4.d0) then
      do is=1,vwfng%nspin
         vwfng_t(1:fc%npwt,1:vwfng%numb_v(is),is)= vwfng%wfn(1:fc%npwt,1:vwfng%numb_v(is),is)
      enddo
      else
         do is=1,vwfng%nspin
           call reorderwfp_col(vwfng%numb_v(is),vwfng%npw,fc%npwt,vwfng%wfn(1,1,is),vwfng_t(1,1,is),vwfng%npw,&
                 & fc%npwt,ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,world_comm )

           !do ii=1,vwfng%numb_v(is)
           !   call mergewf(vwfng%wfn(:,ii,is),evc_g,vwfng%npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
           !   call splitwf(vwfng_t(:,ii,is),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
           !enddo
         enddo
      endif

      do is=1,vwfng%nspin 
         do ii=1,vwfng%numb_v(is),2
            psic(1:fc%nrxxt)=(0.d0,0.d0)
            if (ii==vwfng%numb_v(is)) then
               psic(fc%nlt(1:fc%npwt))  = vwfng_t(1:fc%npwt,ii,is)
               psic(fc%nltm(1:fc%npwt)) = CONJG( vwfng_t(1:fc%npwt,ii,is) )
            else
               psic(fc%nlt(1:fc%npwt))=vwfng_t(1:fc%npwt,ii,is)+(0.d0,1.d0)*vwfng_t(1:fc%npwt,ii+1,is)
               psic(fc%nltm(1:fc%npwt)) =CONJG(vwfng_t(1:fc%npwt,ii,is))+(0.d0,1.d0)*CONJG(vwfng_t(1:fc%npwt,ii+1,is))
            endif
            CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
            vwfnr%wfnrt(1:fc%nrxxt,ii,is)= DBLE(psic(1:fc%nrxxt))
            if(ii/=vwfng%numb_v(is)) vwfnr%wfnrt(1:fc%nrxxt,ii+1,is)= DIMAG(psic(1:fc%nrxxt))
         enddo
      enddo

      deallocate(evc_g)

      call stop_clock('v_wfng_to_wfnr')

      return
      end subroutine

      subroutine c_wfng_to_wfnr(cwfng,fc,cwfnr)
     !this subroutine FFT the conduction wfns to real space in the dual grid

      USE kinds, ONLY : DP
      USE fft_custom_gwl
      USE bse_wannier, ONLY : dual_bse, num_nbndv 
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, nbndx
      USE io_global, ONLY : stdout, ionode, ionode_id
      USE mp_world, ONLY : mpime, nproc,world_comm
      USE mp_wave, ONLY : mergewf,splitwf
      USE mp,             ONLY : mp_sum
      USE gvect
      USE wavefunctions_module, ONLY :  psic

      implicit none

      type(c_state) cwfng
      type(c_state_r) cwfnr
      type(fft_cus) :: fc

      COMPLEX(kind=DP), allocatable :: cwfng_t(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)

      integer :: ii,is
      integer ::icmax  

      call start_clock('c_wfng_to_wfnr')
  
 
!      if (vwfng%nspin==1) then
!         ivmax= vwfng%numb_v(1)
!      else 
!         ivmax=max(vwfng%numb_v(1),vwfng%numb_v(2))
!      endif
!       icmax=nbnd-num_nbndv(1)


      allocate(cwfng_t(fc%npwt,cwfng%numb_c))

      cwfnr%nrxxt=fc%nrxxt       
      cwfnr%numb_c=cwfng%numb_c

      allocate(cwfnr%wfnrt(cwfnr%nrxxt,cwfnr%numb_c))
      
      allocate(evc_g(fc%ngmt_g ))

      if(fc%dual_t==4.d0) then
         cwfng_t(1:fc%npwt,1:cwfng%numb_c)= cwfng%wfn(1:fc%npwt,1:cwfng%numb_c)
      else
         call reorderwfp_col(cwfng%numb_c,cwfng%npw,fc%npwt,cwfng%wfn(1,1),cwfng_t(1,1),cwfng%npw,&
                 & fc%npwt,ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,world_comm )

      endif

      do ii=1,cwfng%numb_c,2
           psic(1:fc%nrxxt)=(0.d0,0.d0)
           if (ii==cwfng%numb_c) then
              psic(fc%nlt(1:fc%npwt))  = cwfng_t(1:fc%npwt,ii)
              psic(fc%nltm(1:fc%npwt)) = CONJG( cwfng_t(1:fc%npwt,ii) )
           else
              psic(fc%nlt(1:fc%npwt))=cwfng_t(1:fc%npwt,ii)+(0.d0,1.d0)*cwfng_t(1:fc%npwt,ii+1)
              psic(fc%nltm(1:fc%npwt)) =CONJG(cwfng_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(cwfng_t(1:fc%npwt,ii+1))
           endif
           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
           cwfnr%wfnrt(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
           if(ii/=cwfng%numb_c) cwfnr%wfnrt(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
      enddo

      deallocate(evc_g)

      call stop_clock('c_wfng_to_wfnr')

      return
      end subroutine

      subroutine write_wfnr(wfnr)
      ! this subroutines writes on disk the type v_state_r for every processor
!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY :  tmp_dir,prefix
      USE mp_world,  ONLY : mpime
      implicit none

      INTEGER, EXTERNAL :: find_free_unit
      type(v_state_r) wfnr
      INTEGER :: iw, iunw,is
      CHARACTER(5) :: nproc

      iunw=find_free_unit()
      
      write(nproc,'(5i1)') &
           & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)


      open( unit=iunw, file=trim(tmp_dir)//trim(prefix)//'.wfnr_t.'// nproc , status='unknown',form='unformatted')

      write(iunw) wfnr%numb_v  
      write(iunw) wfnr%nspin
      write(iunw) wfnr%nrxxt
      
      do is=1,wfnr%nspin
         do iw=1,wfnr%numb_v(is)
            write(iunw)  wfnr%wfnrt(1:wfnr%nrxxt,iw,is)
         enddo
      enddo 
      close(iunw)
      end subroutine

      subroutine read_wfnr(wfnr)
      ! this subroutines reads from disk the type v_state_r for every processor
!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY : tmp_dir, prefix
      USE mp_world,  ONLY : mpime
      implicit none

      INTEGER, EXTERNAL :: find_free_unit
      type(v_state_r) wfnr
      INTEGER :: iw, iunw,is
      CHARACTER(5) :: nproc

      iunw=find_free_unit()
      
      write(nproc,'(5i1)') &
           & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)


      open( unit=iunw, file=trim(tmp_dir)//trim(prefix)//'.wfnr_t.'// nproc , status='old',form='unformatted')

      read(iunw) wfnr%numb_v  
      read(iunw) wfnr%nspin
      read(iunw) wfnr%nrxxt

      do is=1,wfnr%nspin
         do iw=1,wfnr%numb_v(is)
            read(iunw)  wfnr%wfnrt(1:wfnr%nrxxt,iw,is)
         enddo
      enddo
 
      close(iunw)
      end subroutine
	
      subroutine write_cwfnr(wfnr)
      ! this subroutines writes on disk the type v_state_r for every processor
!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY :  tmp_dir,prefix
      USE mp_world,  ONLY : mpime
      implicit none

      INTEGER, EXTERNAL :: find_free_unit
      type(c_state_r) wfnr
      INTEGER :: iw, iunw,is
      CHARACTER(5) :: nproc

      iunw=find_free_unit()
      
      write(nproc,'(5i1)') &
           & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)


      open( unit=iunw, file=trim(tmp_dir)//trim(prefix)//'.cwfnr_t.'// nproc , status='unknown',form='unformatted')

      write(iunw) wfnr%numb_c  
      write(iunw) wfnr%nrxxt
      
      do iw=1,wfnr%numb_c
         write(iunw)  wfnr%wfnrt(1:wfnr%nrxxt,iw)
      enddo

      close(iunw)
      end subroutine

      subroutine read_cwfnr(wfnr)
      ! this subroutines reads from disk the type v_state_r for every processor
!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY : tmp_dir,prefix
      USE mp_world,  ONLY : mpime
      implicit none

      INTEGER, EXTERNAL :: find_free_unit
      type(c_state_r) wfnr
      INTEGER :: iw, iunw,is
      CHARACTER(5) :: nproc

      iunw=find_free_unit()
      
      write(nproc,'(5i1)') &
           & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)


      open( unit=iunw, file=trim(tmp_dir)//trim(prefix)//'.cwfnr_t.'// nproc , status='old',form='unformatted')

      read(iunw) wfnr%numb_c  
      read(iunw) wfnr%nrxxt

      do iw=1,wfnr%numb_c
         read(iunw)  wfnr%wfnrt(1:wfnr%nrxxt,iw)
      enddo

      close(iunw)
      end subroutine
 

      subroutine read_omat(ispin,o)
      ! this subroutines reads the overlap matrix written by pw4gww
!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY : prefix,tmp_dir
      USE io_global,            ONLY : ionode, ionode_id
      USE mp,                   ONLY : mp_bcast
      USE kinds,                ONLY : DP
      USE mp_world,             ONLY : world_comm
      
      implicit none

      INTEGER, EXTERNAL :: find_free_unit

      type(wannier_o) :: o 
      integer ispin

      integer ii,iunu
      real(kind=DP) :: s_bse

      if(ionode) then
         iunu = find_free_unit()
         if (ispin==1) open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.wbse1',status='old',form='unformatted')
         if (ispin==2) open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.wbse2',status='old',form='unformatted')

         read(iunu) o%numb_v
         read(iunu) s_bse

         allocate(o%o(o%numb_v,o%numb_v))

         do ii=1,o%numb_v
            read(iunu) o%o(1:o%numb_v,ii)
         enddo
         close(iunu)
      endif
      
      CALL mp_bcast(o%numb_v, ionode_id , world_comm)
      if(.not.ionode) then     
        allocate(o%o(o%numb_v,o%numb_v))
      endif
      CALL mp_bcast(o%o, ionode_id, world_comm )

      return
      end subroutine

      subroutine read_iimat(iimat,ispin) 
      ! this subroutines reads the ii matrix written by pw4gww
!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY :  prefix, tmp_dir
      USE io_global,            ONLY : ionode, ionode_id
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE kinds,                ONLY : DP

      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      type(ii_mat) :: iimat
      integer ispin

      real(kind=DP) :: s_bse
      integer       :: iv,iuni
      logical       :: debug

      debug=.false.

      if(ionode) then
        iuni = find_free_unit()
        if (ispin==1) open(unit=iuni,file=trim(tmp_dir)//trim(prefix)//'.iwwbse1',status='old',form='unformatted')
        if (ispin==2) open(unit=iuni,file=trim(tmp_dir)//trim(prefix)//'.iwwbse2',status='old',form='unformatted')
        read(iuni) iimat%numb_v
        read(iuni) s_bse
        read(iuni) iimat%np_max
        if(debug) then
          write(*,*) 'From read_iimat numb_v',iimat%numb_v
          write(*,*) 'From read_iimat s_bse', s_bse
          write(*,*) 'From read_iimat np_max',iimat%np_max
        endif
      endif

      CALL mp_bcast(iimat%numb_v, ionode_id, world_comm )
      CALL mp_bcast(iimat%np_max, ionode_id, world_comm )

      allocate(iimat%iimat(iimat%np_max,iimat%numb_v))

      if(ionode) then
        if(debug) then
           write(*,*) 'iimat matrix'  
        endif
        do iv=1, iimat%numb_v
           read(iuni) iimat%iimat(1:iimat%np_max,iv)
           if(debug) then
              write(*,*) 'iv=',iv, iimat%iimat(1:iimat%np_max,iv)
           endif
        enddo
        close(iuni)
      endif

      CALL mp_bcast(iimat%iimat, ionode_id, world_comm )

      return 
      end subroutine

      subroutine read_vww_prod(ispin,numb_v,npw,np_max,iimat,vww)
      !each processor reads the vww(G) written by pw4gww
      !be careful to check that the iimat that is passed to the subroutine is the related to the correct spin channel

!      USE io_files,             ONLY : find_free_unit, prefix,diropn
      USE io_files,             ONLY : prefix,diropn
      USE io_global, ONLY : stdout, ionode 

      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      type(vww_prod) :: vww
      type(ii_mat)   :: iimat
      integer        :: numb_v,npw,np_max,ispin

      integer iv, ip, iungprod, ii,iundebug,i
      logical exst,debug

      debug=.false.
     
      if(debug) then
         iundebug = find_free_unit()
         open(iundebug,file='vww_bse.dat')
      endif      
 
      vww%numb_v=numb_v
      vww%npw=npw
      vww%np_max=np_max

      allocate(vww%vww(npw,np_max,numb_v))

      vww%vww(1:npw,1:np_max,1:numb_v)=dcmplx(0.d0,0.d0)
      
      iungprod = find_free_unit()
      if (ispin==1)  CALL diropn( iungprod, 'vww_bse1.',npw*2, exst)
      if (ispin==2)  CALL diropn( iungprod, 'vww_bse2.',npw*2, exst)

!      if(debug) then 
!         if(ionode) write(stdout,*) 'Read_vww_prod #1'
!      endif

      ii=0
      do iv=1,numb_v
         do ip=1, np_max 
            if(iimat%iimat(ip,iv)>0) then
!               if(debug) then 
!                  if(ionode) write(stdout,*) 'Read_vww_prod #', ii
!               endif
               ii=ii+1
               call davcio(vww%vww(:,ip,iv),npw*2,iungprod,ii,-1)
               if(debug) then
                  if(ionode) then
                     do i=1,npw
                        write(iundebug,*) vww%vww(i,ip,iv)   
                     enddo
                  endif
               endif
            endif
         enddo
      enddo

      close(iungprod)   
      if (debug) close(iundebug)   
      return
      end subroutine

      subroutine read_z(ispin,iimat,z)
      ! the ionode reads the z matrix and broadcast its value to the rest of the
      ! processors.
      !be careful to check that the iimat that is passed to the subroutine is the related to the correct spin channel


!      USE io_files,             ONLY : find_free_unit, prefix
      USE io_files,             ONLY :  prefix, tmp_dir
      USE io_global,            ONLY : ionode, ionode_id
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world, ONLY : world_comm
      USE kinds,                ONLY : DP
      USE io_global, ONLY : stdout,ionode


      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      type(bse_z)    ::z 
      type(ii_mat)   :: iimat
!      integer        :: numw_prod
      integer        ::ispin

      real(kind=DP) :: s_bse
      integer       :: iv,iunz,ii

      logical debug

      debug=.false.

      if(ionode) then 
         iunz = find_free_unit()
         if(debug) then
           if(ionode) write(stdout,*) 'read_z ',trim(tmp_dir)//trim(prefix)//'.zbse1'
         endif


         if (ispin==1) open(unit=iunz,file=trim(tmp_dir)//trim(prefix)//'.zbse1',status='old',form='unformatted')
         if (ispin==2) open(unit=iunz,file=trim(tmp_dir)//trim(prefix)//'.zbse2',status='old',form='unformatted')
         read(iunz) z%numb_v
         read(iunz) s_bse 
         read(iunz) z%np_max
         read(iunz) z%numw_prod

         if(debug) then
           if(ionode) write(stdout,*) 'z%numb_v=', z%numb_v
           if(ionode) write(stdout,*) 's_bse=',s_bse
           if(ionode) write(stdout,*) 'z%np_max=',z%np_max
           if(ionode) write(stdout,*) 'z%numw_prod=', z%numw_prod
         endif

      endif 

      CALL mp_bcast(z%numb_v, ionode_id, world_comm )
      CALL mp_bcast(z%np_max, ionode_id, world_comm )
      CALL mp_bcast(z%numw_prod, ionode_id, world_comm )
      call mp_barrier(world_comm)
      
      allocate(z%z(z%numw_prod,z%np_max,z%numb_v))
     
      if(ionode) then
         do iv=1, z%numb_v
            do ii=1,z%np_max
               if(debug) then
                 if(ionode) write(stdout,*)'read_z, ii=',ii 
               endif
               if (iimat%iimat(ii,iv)>0) read(iunz) z%z(:,ii,iv)
            enddo
         enddo
      endif

      if(debug) then
        if(ionode) write(stdout,*) 'read_z #1'
      endif
!

      CALL mp_bcast(z%z, ionode_id, world_comm )
      call mp_barrier(world_comm)

      if(debug) then
        if(ionode) write(stdout,*) 'read_z #2'
      endif


      if(ionode) close(iunz)
      FLUSH( stdout )      
 
      return

      end subroutine

END MODULE
