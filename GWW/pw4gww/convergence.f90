
MODULE   convergence_gw

  USE kinds, ONLY : DP

  SAVE


  TYPE vpv
     INTEGER :: nr!numer of points 
     INTEGER(kind=DP), DIMENSION(:), POINTER :: r!points in real space where vPv will be evaluated
                                               !the point is /= 0 if it belongs to this MPI task
                                              !in this case it is the point index on charge grid
     REAL(kind=DP), DIMENSION(:,:), POINTER :: vpvr!vpv(:,r) in real space charge grid, real part
     REAL(kind=DP), DIMENSION(:,:), POINTER :: vpvr_im!vpv(:,r) in real space charge grid, imaginary part 
     REAL(kind=DP) :: freq!frequency at which it is calculated, real part
     REAL(kind=DP) :: freq_im!frequency at which it is calculated, imaginary part  
  END type vpv


  TYPE gzero
     INTEGER :: nr!numer of points  
     INTEGER(kind=DP), DIMENSION(:), POINTER :: r!points in real space where vPv will be evaluated 
                                                 !the point is /= 0 if it belongs to this MPI task
                                                 !in this case it is the point index on charge grid 
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: gf!Green's function (:,nspin,r) in real space charge grid, real part
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: gf_im!Green's function (:,nspin,r) in real space charge grid, imaginary part 
     REAL(kind=DP) :: freq!frequency at which it is calculated, real part
     REAL(kind=DP) :: freq_im!frequency at which it is calculated, imaginary part     
  END type gzero

  TYPE exchange
     INTEGER :: nr!numer of points
     INTEGER(kind=DP), DIMENSION(:), POINTER :: r!points in real space where vPv will be evaluated  
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: x!Exchange function (:,r,nspin) in real space charge grid 
  END type exchange

  TYPE hks
     INTEGER :: nr!numer of points 
     INTEGER(kind=DP), DIMENSION(:), POINTER :: r!points in real space where vPv will be evaluated
     REAL(kind=DP), DIMENSION(:,:), POINTER :: h0!Hamiltonian (:,r) in real space charge grid
     REAL(kind=DP), DIMENSION(:,:), POINTER :: vxc!V_xc(:,r) in real space charge grid 
  END type hks



  TYPE convergence_tests
     INTEGER :: nr!numer of points
     INTEGER :: iband!band index, we put iband=0 if no specific band is associated to the data
     INTEGER :: ispin!spin index,  we put iband=0 if no specific band is associated to the data
     INTEGER,DIMENSION(:), POINTER :: r!points in real space where vPv will be evaluated 
                                               !the point is /= 0 if it belongs to this MPI task 
                                              !in this case it is the point index on charge grid
     INTEGER :: r_coord(3)!coordinates in grid space of point, only one point is assumed
     COMPLEX(kind=DP), POINTER :: freq(:)!frequencies at which it is calculated
     INTEGER :: nf!number of frequencies
     TYPE(vpv),POINTER  :: wapprox(:)!vPV or also  an aproximated W
     TYPE(gzero),POINTER :: g0(:)!Green's functions
     TYPE(exchange) :: xf!exact exchange
     TYPE(hks) :: h0!Hamiltonian KS
     REAL(kind=DP), DIMENSION(:), POINTER :: energy!energy (in Ry) calculated at  points r
  END type convergence_tests
  

  TYPE self_energy
     INTEGER :: s_first_state!these two specify the range of KS states
     
     INTEGER :: s_last_state!
     INTEGER :: s_first_spin, s_last_spin!!these two specify the range of KS spin
     INTEGER :: ngrid!grid points in imaginary frequency
     REAL(kind=DP), POINTER :: freq(:)!frequencies (ngrid)
     INTEGER :: nr!total number of r points
     REAL(kind=DP), POINTER ::   psif(:,:,:) !wavefunction at point dimension (nr,s_first_state:s_last_state,s_first_spin,s_last_spin)
     COMPLEX(kind=DP), POINTER :: self(:,:,:,:)!self_energy(ngrid,nr,s_first_state:s_last_state,s_first_spin,s_last_spin)  
     LOGICAL :: l_first!if true, the self_energy has not yet been created
     REAL(kind=DP), POINTER :: ene_ks(:,:)!KS expectation values (nbnd,nspin)
     REAL(kind=DP), POINTER :: ene_xc(:,:)!DFT-XC expectation values (nbnd,nspin)
     REAL(kind=DP), POINTER :: ene_x(:,:)!exact exchnage expectation values (nbnd,nspin)

  END TYPE self_energy

 


  INTERFACE free_memory
     MODULE PROCEDURE free_vpv
     MODULE PROCEDURE free_convergence_tests
     MODULE PROCEDURE free_gzero
     MODULE PROCEDURE free_exchange
     MODULE PROCEDURE free_hks
     MODULE PROCEDURE free_self_energy
   END INTERFACE free_memory


  INTERFACE initialize_memory
     MODULE PROCEDURE initialize_memory_gzero
     MODULE PROCEDURE initialize_memory_vpv
     MODULE PROCEDURE initialize_memory_self_energy
   END  INTERFACE initialize_memory


 
  CONTAINS

 
    
    SUBROUTINE  initialize_memory_self_energy(se)
      IMPLICIT NONE
      TYPE(self_energy) ::se
      nullify(se%freq,se%psif,se%self,se%ene_ks,se%ene_xc,se%ene_x)
      se%l_first=.true.
    END SUBROUTINE initialize_memory_self_energy

    SUBROUTINE free_self_energy(se)
      IMPLICIT NONE
      TYPE(self_energy) ::se
      if(associated(se%freq)) deallocate(se%freq)
      if(associated(se%self)) deallocate(se%self)
      if(associated(se%psif)) deallocate(se%psif)
      if(associated(se%ene_ks)) deallocate(se%ene_ks)
      if(associated(se%ene_xc)) deallocate(se%ene_xc)
      if(associated(se%ene_x)) deallocate(se%ene_x)
      nullify(se%freq,se%psif,se%self,se%ene_ks,se%ene_xc,se%ene_x)
    END SUBROUTINE free_self_energy

    SUBROUTINE  initialize_memory_gzero(c)
      IMPLICIT NONE
      TYPE(gzero) :: c
      nullify(c%r,c%gf,c%gf_im)
      RETURN
    END SUBROUTINE initialize_memory_gzero

    
    SUBROUTINE  initialize_memory_vpv(c)
      IMPLICIT NONE
      TYPE(vpv) :: c
      nullify(c%r,c%vpvr,c%vpvr_im)
      RETURN
    END SUBROUTINE initialize_memory_vpv

    

    SUBROUTINE  free_hks(c)
      IMPLICIT NONE
      TYPE(hks) :: c
      deallocate(c%r,c%h0,c%vxc)
      RETURN
    END SUBROUTINE free_hks


    
  SUBROUTINE  free_exchange(c)
    IMPLICIT NONE
    TYPE(exchange) :: c
    deallocate(c%r,c%x)
    RETURN
  END SUBROUTINE free_exchange


    

    SUBROUTINE  free_gzero(c)
      IMPLICIT NONE
      TYPE(gzero) :: c
      if(associated(c%r)) deallocate(c%r)
      if(associated(c%gf)) deallocate(c%gf)
      if(associated(c%gf_im)) deallocate(c%gf_im)
      nullify(c%r,c%gf,c%gf_im)
    RETURN
  END SUBROUTINE free_gzero

  SUBROUTINE set_se_energies(se, ene_ks,ene_xc,ene_x)
    USE wvfct,                ONLY : nbnd
    USE lsda_mod,             ONLY : nspin
    
    IMPLICIT NONE

    TYPE(self_energy), INTENT(out) :: se
    REAL(kind=DP), INTENT(in) :: ene_ks(nbnd,nspin)
    REAL(kind=DP), INTENT(in) :: ene_xc(nbnd,nspin)
    REAL(kind=DP), INTENT(in) :: ene_x(nbnd,nspin)

    allocate(se%ene_ks(nbnd,nspin))
    allocate(se%ene_xc(nbnd,nspin))
    allocate(se%ene_x(nbnd,nspin))

    se%ene_ks(1:nbnd,1:nspin)=ene_ks(1:nbnd,1:nspin)
    se%ene_xc(1:nbnd,1:nspin)=ene_xc(1:nbnd,1:nspin)
    se%ene_x(1:nbnd,1:nspin)=ene_x(1:nbnd,1:nspin)

    
    
    
    RETURN
  END SUBROUTINE set_se_energies



  SUBROUTINE write_self_energy(se, ix_start,iy_start,iz_start,nr_counter)
    USE io_global, ONLY :stdout,  ionode
    USE io_files,             ONLY : prefix

    IMPLICIT NONE

    TYPE(self_energy) ::se !object to be written on disk
    INTEGER :: ix_start,iy_start,iz_start! last point coordinates
    INTEGER :: nr_counter !last counter position

    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iun, ii,is

    if(ionode) then
       iun=find_free_unit()
       open( unit= iun, file=trim(prefix)//'.easyself', status='unknown', form='unformatted')
       write(iun) ix_start,iy_start,iz_start
       write(iun) nr_counter
       write(iun) se%s_first_state
       write(iun) se%s_last_state
       write(iun) se%s_first_spin
       write(iun) se%s_last_spin
       write(iun) se%ngrid
       write(iun) se%nr
       write(iun) se%l_first
       write(iun)se%freq(1:se%ngrid)
       do is=se%s_first_spin,se%s_last_spin
          do ii=se%s_first_state,se%s_last_state
             write(iun) se%psif(1:se%nr,ii,is)
          enddo
       enddo
       do is=se%s_first_spin,se%s_last_spin
          do ii=se%s_first_state,se%s_last_state
             write(iun) se%self(1:se%ngrid,1:se%nr,ii,is)
          enddo
       enddo
       close(iun)
    endif
    


  END SUBROUTINE write_self_energy

 SUBROUTINE read_self_energy(se, ix_start,iy_start,iz_start,nr_counter)

    USE io_global, ONLY :stdout,  ionode, ionode_id
    USE io_files,             ONLY : prefix
    USE mp, ONLY : mp_bcast
    USE mp_global, ONLY : world_comm
 

    IMPLICIT NONE

    TYPE(self_energy) ::se !object to be read on disk  
    INTEGER :: ix_start,iy_start,iz_start! last point coordinates
    INTEGER :: nr_counter !last counter position  

    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iun, ii, is

    if(ionode) then
       iun=find_free_unit()
       open( unit= iun, file=trim(prefix)//'.easyself', status='old', form='unformatted')
       read(iun) ix_start,iy_start,iz_start
       read(iun) nr_counter
       read(iun) se%s_first_state
       read(iun) se%s_last_state
       read(iun) se%s_first_spin
       read(iun) se%s_last_spin
       read(iun) se%ngrid
       read(iun) se%nr
       read(iun) se%l_first
    endif
    call mp_bcast(ix_start,ionode_id,world_comm)
    call mp_bcast(iy_start,ionode_id,world_comm)
    call mp_bcast(iz_start,ionode_id,world_comm)
    call mp_bcast(nr_counter,ionode_id,world_comm)
    call mp_bcast(se%s_first_state,ionode_id,world_comm)
    call mp_bcast(se%s_last_state,ionode_id,world_comm)
    call mp_bcast(se%s_first_spin,ionode_id,world_comm)
    call mp_bcast(se%s_last_spin,ionode_id,world_comm)
    call mp_bcast(se%ngrid,ionode_id,world_comm)
    call mp_bcast(se%nr,ionode_id,world_comm)
    call mp_bcast(se%l_first,ionode_id,world_comm)

    allocate(se%freq(se%ngrid))
    allocate(se%psif(1:se%nr,se%s_first_state:se%s_last_state,se%s_first_spin:se%s_last_spin))
    allocate(se%self(1:se%ngrid,1:se%nr,se%s_first_state:se%s_last_state,se%s_first_spin:se%s_last_spin))
    if(ionode) then
       read(iun)se%freq(1:se%ngrid)
       do is=se%s_first_spin,se%s_last_spin
          do ii=se%s_first_state,se%s_last_state
             read(iun) se%psif(1:se%nr,ii,is)
          enddo
       enddo
       do is=se%s_first_spin,se%s_last_spin
          do ii=se%s_first_state,se%s_last_state
             read(iun) se%self(1:se%ngrid,1:se%nr,ii,is)
          enddo
       enddo
       close(iun)
    endif
    call mp_bcast(se%freq,ionode_id,world_comm)
    do is=se%s_first_spin,se%s_last_spin
       do ii=se%s_first_state,se%s_last_state
          call mp_bcast(se%psif(1:se%nr,ii,is),ionode_id,world_comm)
       enddo
    enddo
    do is=se%s_first_spin,se%s_last_spin
       do ii=se%s_first_state,se%s_last_state
          call mp_bcast(se%self(1:se%ngrid,1:se%nr,ii,is),ionode_id,world_comm)
       enddo
    enddo


  END SUBROUTINE read_self_energy




  SUBROUTINE check_normalisation(v_states, iband)

      USE fft_base,             ONLY : dfftp, dffts
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : world_comm,mpime
      USE io_global, ONLY : stdout
      USE wannier_gw

      IMPLICIT NONE
      INTEGER :: iband!band considered for convergence test  
      REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,num_nbndv(1))!valence states in real space 

      INTEGER :: ix,iy,iz, id,nd,ii,ib,jb,nd2
      REAL(kind=DP) ::sca

      do ib=1,10!iband
         do jb=4,4!1,iband            
            do id=1,30
               sca=0.d0
               nd=0
               nd2=0
               do iz=1,dffts%nr3,id
                  do iy=1,dffts%nr2,id
                     do ix=1,dffts%nr1,id
                        if( iz>dffts%nr3p_offset(mpime+1) .and. iz <=( dffts%nr3p_offset(mpime+1)+dffts%my_nr3p)) then
                           ii=(iz-dffts%nr3p_offset(mpime+1)-1)*dffts%nr2*dffts%nr1+&
                                (iy-1)*dffts%nr1+ix
                           if(abs(v_states(ii,jb))> 1.d-0) then
                              sca=sca+v_states(ii,ib)*v_states(ii,jb)
                              nd2=nd2+1
                           endif
                           nd=nd+1
                        endif
                     enddo
                  enddo
               enddo
               call mp_sum(nd,world_comm)
               call mp_sum(nd2,world_comm)
               call mp_sum(sca,world_comm)
               sca=sca/dble(nd)
               write(stdout,*) 'NORMALIZATION, STEP :',ib,jb,id,nd,nd2,sca
            enddo
         enddo
      enddo
    return
  END SUBROUTINE check_normalisation



  SUBROUTINE start_convergence(ct,iband,ispin,v_states,lmax,nx,ny,nz,nf, freq,ks_wfcs)

    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2,bg
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE mp_world,  ONLY : mpime, nproc,world_comm
    USE lanczos
    USE io_files,              ONLY : prefix, tmp_dir
    USE mp_pools,  ONLY : intra_pool_comm
    USE mp_wave, ONLY : splitwf
    USE lsda_mod, ONLY : nspin
    USE scf,       ONLY : rho,rho_core,rhog_core,scf_type,create_scf_type,destroy_scf_type

    IMPLICIT NONE
    TYPE(convergence_tests) :: ct!object for convergence tests
    INTEGER :: iband!band considered for convergence test
    INTEGER :: ispin!spin considered for convergence test
    REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,num_nbnds,nspin)!valence states in real space
    LOGICAL, INTENT(in) :: lmax!if true selec the maximu for v_states otherwise x,y,z
    INTEGER, INTENT(in) :: nx,ny,nz
    INTEGER, INTENT(in) :: nf!number of frequencies 
    COMPLEX(kind=DP), INTENT(in) :: freq(nf)
    COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd,nspin)!KS wavefunctions in G space

    INTEGER :: ip,ir,imax,iw
    REAL(kind=DP) :: charge,charge_max
    LOGICAL :: lfound
    LOGICAL :: l_old
    COMPLEX(kind=DP), ALLOCATABLE :: psi_old(:),psi_old_vpv(:,:,:)
    TYPE(lanczos_chain) :: lc(2),lc_vpv(2)
    REAL(kind=DP), ALLOCATABLE :: head(:),head1(:),head2(:),head3(:)
    INTEGER :: iun,nh
    REAL(kind=DP), ALLOCATABLE :: freqh(:)
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP) :: omegah
!!!!variables for wings                                                                                                                                                                                            
    INTEGER :: n_g, ngm_k
    REAL(kind=DP) :: omega_g
    INTEGER :: ii, ipol
    COMPLEX(kind=DP), ALLOCATABLE :: e_head(:,:,:), e_head_g0(:)
    INTEGER :: npwx_g
    INTEGER, ALLOCATABLE :: k2g_ig_l2g(:)
    REAL(kind=DP), ALLOCATABLE  :: freqs_head(:)
    COMPLEX(kind=DP), ALLOCATABLE :: wing(:)
    TYPE(scf_type) :: aile
    REAL(kind=DP) :: alpha

   
    REAL(kind=DP), ALLOCATABLE :: vc(:,:)
    REAL(kind=DP) :: broadening
    COMPLEX(kind=DP), ALLOCATABLE :: phi_save(:)

    call create_scf_type ( aile , .true. )        
    
    l_wing=.true.
    write(stdout,*) 'Routine calculate_convergence'
    if(.not.l_truncated_coulomb .and. l_wing) then
!here reads wings stuff                                                                                                                                                                                            
       allocate(k2g_ig_l2g(ngm))
       call ktogamma_ig_l2g ( k2g_ig_l2g, at, bg )
       if(ionode) then
          iun =  find_free_unit()
          open( unit= iun, file=trim(tmp_dir)//'/_ph0/'//trim(prefix)//'.e_head', status='old',form='unformatted')
          read(iun) n_g
          read(iun) omega_g
       endif
       call mp_bcast(n_g, ionode_id,world_comm)
       call mp_bcast(omega_g, ionode_id,world_comm)
       allocate(freqs_head(n_g+1))
       if(ionode) then
          read(iun) freqs_head(1:n_g+1)
          read(iun) ngm_k
       endif
       call mp_bcast(freqs_head, ionode_id,world_comm)
       call mp_bcast(ngm_k, ionode_id,world_comm)

       allocate(e_head_g0(ngm_k))
       allocate(e_head(npw, n_g+1,3))
       e_head(:,:,:)= (0.d0,0.d0)
       do ii=1,n_g+1
          do ipol=1,3
             e_head_g0(:)=(0.d0,0.d0)
             if(ionode) read(iun) e_head_g0(1:ngm_k)
             call splitwf(e_head(:, ii,ipol),e_head_g0,npw,k2g_ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
          enddo
       enddo
       if(ionode) close(iun)
       deallocate(e_head_g0)
       allocate(wing(npw))
    else
       allocate(wing(npw))
    endif


    call initialize_memory(lc(1))
    call initialize_memory(lc(2))
    call initialize_memory(lc_vpv(1))
    call initialize_memory(lc_vpv(2))

    allocate(head(nf))
    if(.not.l_truncated_coulomb) then
       if(ionode) then
          allocate(head1(nf),head2(nf),head3(nf))
          allocate(freqh(nf))
          iun =  find_free_unit()
          open( unit= iun, file=trim(tmp_dir)//'/_ph0/'//trim(prefix)//'.head', status='old',form='unformatted')
          read(iun) nh
          read(iun) omegah
          read(iun) freqh(1:nh+1)
          read(iun) head1(1:nh+1)
          read(iun) head2(1:nh+1)
          read(iun) head3(1:nh+1)
          close(iun)
          do iw=1,nf
             head(nf-iw+1)=(1.d0/3.d0)*(head1(iw)+head2(iw)+head3(iw))
             head(nf-iw+1)=head1(iw)!DEBUG
             !if(l_verbose) write(stdout,*) 'HEAD: ', iw, head(nf-iw+1),nf,nh
          enddo
          deallocate(head1,head2,head3,freqh)
       endif
       call mp_bcast(head, ionode_id, world_comm)
    endif
    ct%nr=1
    ct%iband=iband
    ct%ispin=ispin
    allocate(ct%r(1))
    ct%r=0
    ct%nf=nf
    allocate(ct%freq(nf))
    ct%freq(1:nf)=freq(1:nf)
    
    allocate(ct%wapprox(nf))
    allocate(ct%g0(nf))
   
    if(lmax) then
       charge_max=-1.d0 
       do ip=0,nproc-1
          if(mpime==ip) then
             imax=0          
             do ir=1,dffts%nnr
                charge=(v_states(ir,ct%iband,ispin))**2.d0
                if(charge>charge_max) then
                   charge_max=charge
                   imax=ir
                   ct%r(1)=imax
                endif
             enddo
             if(imax>0) then
                lfound=.true.
             else
                lfound=.false.
             endif
             call mp_bcast(lfound,ip,world_comm)
             call mp_bcast(charge_max,ip,world_comm)   
          else
             call mp_bcast(lfound,ip,world_comm)  
             call mp_bcast(charge_max,ip,world_comm)   
             if(lfound) then
                ct%r(1)=0
             endif
             
          endif
       enddo
    else
       if( nz>dffts%nr3p_offset(mpime+1) .and. nz <=( dffts%nr3p_offset(mpime+1)+dffts%my_nr3p)) then
          ct%r(1)=(nz-dffts%nr3p_offset(mpime+1)-1)*dffts%nr2*dffts%nr1+&
               (ny-1)*dffts%nr1+nx


       endif
       ct%iband=0
       ct%ispin=0
       ct%r_coord(1)=nx
       ct%r_coord(2)=ny
       ct%r_coord(3)=nz
    endif
   if(lmax) write(stdout,*) 'MAXIMUM PSI :',ct%iband, sqrt(charge_max)
  
    allocate(psi_old(npw),psi_old_vpv(npw,num_nbndv(1),2))
    !call calculate_vpv_complex(ct%wapprox(1),0.d0,0.d0,ct%r,ct%nr,v_states,.true.,1.d-12,.false.,psi_old_vpv)
  !call calculate_gzero(ct%g0(1),0.d0,ct%r,ct%nr)
    alpha=0.5d0!1.d0 !DEBUG
    allocate(phi_save(npw))
    do iw=1,ct%nf
       if(l_verbose) write(stdout,*) 'IMAGINARY FREQ:',aimag(ct%freq(iw))
       if(iw==1) then
          l_old=.false.
       else
          l_old=.true.
       endif
     
       if(.not.l_truncated_coulomb .and. l_wing) then
          call read_wing ( aile, 1, .true.,1,ct%nf-iw+1 ) 
          wing(1:npw)=(1.d0/3.d0)*(e_head(1:npw,ct%nf-iw+1,1)+e_head(1:npw,ct%nf-iw+1,2)+e_head(1:npw,ct%nf-iw+1,3))
          wing(1:npw)=aile%of_g(1:npw,1)!DEBUG
       endif

       if(aimag(ct%freq(iw))==0.d0) then
          broadening = (aimag(ct%freq(iw))+aimag(ct%freq(iw-1)))/4.d0
          write(stdout,*) 'DEBUG BROADENIG', broadening
       else
          broadening =0.d0
       endif
       call calculate_vpv_complex(ct%wapprox(iw),dble(ct%freq(iw)),aimag(ct%freq(iw))+broadening,ct%r,ct%nr,v_states,.true.,&
        &   easy_w_thrs,lc_vpv,head(iw),wing, alpha,ks_wfcs, l_old,phi_save)
       call calculate_gzero_complex(ct%g0(iw),dble(ct%freq(iw)),aimag(ct%freq(iw))+broadening,ct%r,ct%nr,l_old,psi_old, &
            l_easy_lanczos_g, lc,ks_wfcs)
       if(broadening /= 0.d0) then
          ct%g0(iw)%gf_im=0.d0
       endif

    enddo
    deallocate(phi_save)
    deallocate(psi_old)
    call calculate_x(ct%xf,ct%r,ct%nr,v_states)
    call calculate_hks(ct%h0,ct%r,ct%nr,v_states)
    call free_memory(lc(1))
    call free_memory(lc(2))
    call free_memory(lc_vpv(1))
    call free_memory(lc_vpv(2))
    deallocate(head)
!deallocate wing part                                                                                                                                                                                              
    if(.not.l_truncated_coulomb.and.l_wing) then
       deallocate(k2g_ig_l2g,freqs_head)
       deallocate(e_head)
      
    endif
    deallocate(wing)
    !                
    call  destroy_scf_type (aile )
    return
  END SUBROUTINE start_convergence


  SUBROUTINE calculate_convergence(ct,v_states,se,nr_counter)

    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2,bg
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : prefix, tmp_dir
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE mp_world,  ONLY : mpime, nproc,world_comm
    USE constants, ONLY : rytoev
    USE lsda_mod,      ONLY : nspin
    USE scf,       ONLY : rho,rho_core,rhog_core



    IMPLICIT NONE
    TYPE(convergence_tests) :: ct!object for convergence tests 
    REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,num_nbnds,nspin)!valence states in real space
    TYPE(self_energy) ::se !object for collecting calculated self-energies
    INTEGER, INTENT(inout):: nr_counter!counter for self_energy
   

    INTEGER :: ir,ii,iw,jw,dw,iww,jww,dww,jj,ix,iy,iz,ipol,is
    REAL(kind=DP) :: sca,psif,sca1,sca2,sca3
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iun,iun1,iun2
    REAL(kind=DP), ALLOCATABLE :: vprods(:,:), vprods_im(:,:)
    REAL(kind=DP), ALLOCATABLE :: wterms(:,:), wterms_im(:,:)
    REAL(kind=DP), ALLOCATABLE :: gwmat_rr(:,:),gwmat_ri(:,:),gwmat_ir(:,:),gwmat_ii(:,:)
    CHARACTER(4) :: nfile,nfilex,nfiley,nfilez
    CHARACTER(5) :: nfilel,nfilel2
    COMPLEX(kind=DP), ALLOCATABLE :: cmat(:,:),sigmac(:)
    REAL(kind=DP) emax,freq_im
    INTEGER :: ifirst,ilast,kk,ifirst_spin,ilast_spin
    INTEGER :: passed,nr_calculate
    REAL(kind=DP) :: vpsif, psif_xc,psif_x
    INTEGER :: irr
    CHARACTER(4) :: nfile_orb



    write(stdout,*) 'Routine calculate_convergence'

    if(l_whole_s) then
       ifirst=1
       ilast=num_nbnds
       ifirst_spin=1
       ilast_spin=nspin
    else
       if(ct%iband>0) then
          ifirst=ct%iband
          ilast=ct%iband
          ifirst_spin=ct%ispin
          ilast_spin=ct%ispin
       else
          ifirst=s_first_state
          ilast=s_last_state
          ifirst_spin=s_first_spin
          ilast_spin=s_last_spin
       endif
    endif

    if(se%l_first) then
       se%s_first_state=ifirst
       se%s_last_state=ilast
       se%s_first_spin=ifirst_spin
       se%s_last_spin=ilast_spin

       se%ngrid=2*(ct%nf-1)+1
       
       allocate(se%freq(se%ngrid))
       do jj=-ct%nf+1,ct%nf-1
          se%freq(jj+ct%nf)= abs(aimag(ct%freq(1)))/dble(ct%nf)*jj
       enddo
       if(ct%iband>0) then
          se%nr=s_last_state-s_first_state+1
       else
          if(easy_grid_type==1) then
             if(easy_psi_thrs==0.d0 .and. easy_grid_param(5)>=10000) then
                se%nr=(int((dffts%nr1-easy_grid_param(1)-1)/easy_grid_param(4))+1)* &
      &(int((dffts%nr2-easy_grid_param(2)-1)/easy_grid_param(4))+1)*(int((dffts%nr3-easy_grid_param(3)-1)/easy_grid_param(4))+1)
             else
                nr_calculate=0.d0
                do iz=1+easy_grid_param(3),dffts%nr3,easy_grid_param(4)
                   do iy=1+easy_grid_param(2),dffts%nr2,easy_grid_param(4)
                      do ix=1+easy_grid_param(1),dffts%nr1,easy_grid_param(4)
                         passed=0
                         if(ix<easy_grid_param(5) .and. iy<easy_grid_param(5) .and. iz<easy_grid_param(5)) then
                            if( iz>dffts%nr3p_offset(mpime+1) .and. iz <=( dffts%nr3p_offset(mpime+1)+dffts%my_nr3p)) then
                               ii=(iz-dffts%nr3p_offset(mpime+1)-1)*dffts%nr2*dffts%nr1+&
                                    (iy-1)*dffts%nr1+ix
                               do is=se%s_first_spin,se%s_last_spin
                                  do jj=s_first_state,s_last_state
                                     if(abs(v_states(ii,jj,is))> easy_psi_thrs) passed=1
                                  enddo
                               enddo
                            endif
                         endif
                         call mp_sum(passed,world_comm)
                         if(passed>0) nr_calculate=nr_calculate+1
                      enddo
                   enddo
                enddo
                se%nr=nr_calculate
             endif
          else
             se%nr=easy_grid_param(4)
          endif
          write(stdout,*) 'TOTAL NUMBER OF R POINTS TO BE CALCULATED:', se%nr
       endif
       allocate(se%psif(se%nr,se%s_first_state:se%s_last_state,se%s_first_spin:se%s_last_spin))
       if(  easy_grid_type /= 1 ) then
          allocate(se%self(se%ngrid,se%nr,se%s_first_state:se%s_last_state,se%s_first_spin:se%s_last_spin))
       else
           allocate(se%self(1,1,se%s_first_state:se%s_last_state,se%s_first_spin:se%s_last_spin))
       endif
       se%self=0.d0
       se%l_first=.false.
    endif

    allocate(ct%energy(ct%nr))


    do ir=1,ct%nr!NOTE THAT ONLY 1 POINTS HAS BEEN IMPLEMENTED SO FAR
!!!!!THE FOLLOWING PART WAS TO USE AS A CONVERGENCE ESTIMATOR IN GWL        
!       if(ct%iband>0) then
!          sca=0.d0
!          sca1=0.d0
!          sca2=0.d0
!          sca3=0.d0
!          do ii=1,dffts%nnr
!             sca=sca+ct%g0(1)%gf(ii,ir)*ct%wapprox(1)%vpvr(ii,ir)*v_states(ii,ct%iband)
!             sca1=sca1+ct%xf%x(ii,ir)*v_states(ii,ct%iband)
!             sca2=sca2+ct%h0%h0(ii,ir)*v_states(ii,ct%iband)
!             sca3=sca3+v_states(ii,ct%iband)**2.d0
!          enddo
!          call mp_sum(sca,world_comm)
!          sca=sca/dble((dffts%nr1*dffts%nr2*dffts%nr3))
!          call mp_sum(sca1,world_comm)
!          sca1=sca1/dble((dffts%nr1*dffts%nr2*dffts%nr3))
!          call mp_sum(sca2,world_comm)
!          sca2=sca2/dble((dffts%nr1*dffts%nr2*dffts%nr3))
!          call mp_sum(sca3,world_comm)
!          sca3=sca3/(dffts%nr1*dffts%nr2*dffts%nr3)
!          
!
!          ct%energy(ir)=sca+sca1
!          psif=0.d0
!          if(ct%r(ir)/=0) psif=v_states(ct%r(ir),ct%iband)
!          call mp_sum(psif,world_comm)
!          ct%energy(ir)=ct%energy(ir)/psif
!          
!          write(stdout,*) 'CONVERGED ENERGY FOR TEST:', sca/psif*rytoev,sca1/psif*rytoev,sca2/psif*rytoev,sca3
!          write(stdout,*)  ct%energy(ir)*rytoev,psif
!       endif
!!!!!!!!!!TILL HERE

       if(easy_grid_type==1.or.easy_grid_type==2) then
          write(stdout,*) 'POINT NUMBER', nr_counter
       endif
       do is=ifirst_spin,ilast_spin
          do kk=ifirst,ilast
             do iw=1,ct%nf
                sca=0.d0
                sca1=0.d0
                do ii=1,dffts%nnr
                   sca=sca+ct%g0(iw)%gf(ii,is,ir)*v_states(ii,kk,is)
                   sca1=sca1+ct%g0(iw)%gf_im(ii,is,ir)*v_states(ii,kk,is)
                enddo
                call mp_sum(sca,world_comm)
                sca=sca/dble((dffts%nr1*dffts%nr2*dffts%nr3))
                call mp_sum(sca1,world_comm)
                sca1=sca1/dble((dffts%nr1*dffts%nr2*dffts%nr3))
                psif=0.d0
                if(ct%r(ir)/=0) psif=v_states(ct%r(ir),kk,is)
                call mp_sum(psif,world_comm)
               ! if(l_verbose) write(stdout,*) iw, sca/psif*rytoev,sca1/psif*rytoev
             
             enddo
             do iw=1,ct%nf
                sca=0.d0
                sca1=0.d0
                do ii=1,dffts%nnr
                   sca=sca+ct%wapprox(iw)%vpvr(ii,ir)*v_states(ii,kk,is)
                   sca1=sca1+ct%wapprox(iw)%vpvr_im(ii,ir)*v_states(ii,kk,is)
                enddo
                call mp_sum(sca,world_comm)
                sca=sca/dble((dffts%nr1*dffts%nr2*dffts%nr3))
                call mp_sum(sca1,world_comm)
                sca1=sca1/dble((dffts%nr1*dffts%nr2*dffts%nr3))
                psif=0.d0
                if(ct%r(ir)/=0) psif=v_states(ct%r(ir),kk,is)
                call mp_sum(psif,world_comm)
                !if(l_verbose) write(stdout,*) iw, sca/psif*rytoev,sca1/psif*rytoev
             
             enddo
             if(ionode .and. kk==ct%iband .and. l_verbose) then
                iun=find_free_unit()
                write(nfile,'(4i1)') ct%iband/1000,mod(ct%iband,1000)/100,mod(ct%iband,100)/10,mod(ct%iband,10)
                open( unit= iun, file=trim(prefix)//'.convdata.'//nfile, status='unknown')
                write(iun,*) ct%nf
                do iw=1,ct%nf
                   write(iun,*) ct%wapprox(iw)%freq
                enddo
             endif
             psif=0.d0
             if(ct%r(ir)/=0) psif=v_states(ct%r(ir),kk,is)
             call mp_sum(psif,world_comm)
             
             


             psif_x=0.d0
             do irr=1,dffts%nnr
                psif_x=psif_x+ct%xf%x(irr,1,is)*v_states(irr,kk,is)
             enddo
             call mp_sum(psif_x,world_comm)

             psif_xc=0.d0
             if(ct%r(1)/=0) psif_xc=psif*ct%h0%vxc(ct%r(1),is)
             call mp_sum(psif_xc,world_comm)
             
            
             if(l_whole_s .or. easy_grid_type==1.or.easy_grid_type==2) se%psif(nr_counter,kk,is)=psif
             if( easy_grid_type==1.or.easy_grid_type==2) then
                write(stdout,*) 'WAVE_FUNCTIONS BAND AND VALUE',kk,psif,0.d0,0.d0
             else
                write(stdout,*) 'WAVE_FUNCTIONS OF MAX ', ct%iband, ct%ispin, '  in ', kk,psif,0.d0,0.d0
             endif
             allocate(vprods(dffts%nnr,ct%nf),wterms(dffts%nnr,ct%nf))
             allocate(vprods_im(dffts%nnr,ct%nf),wterms_im(dffts%nnr,ct%nf))
             allocate(gwmat_rr(ct%nf,ct%nf),gwmat_ri(ct%nf,ct%nf),gwmat_ir(ct%nf,ct%nf),gwmat_ii(ct%nf,ct%nf))
          
          
             do iw=1,ct%nf
                vprods(1:dffts%nnr,iw)=ct%g0(iw)%gf(1:dffts%nnr,is,ir)*v_states(1:dffts%nnr,kk,is)
                wterms(1:dffts%nnr,iw)=ct%wapprox(iw)%vpvr(1:dffts%nnr,ir)
                vprods_im(1:dffts%nnr,iw)=ct%g0(iw)%gf_im(1:dffts%nnr,is,ir)*v_states(1:dffts%nnr,kk,is)
                wterms_im(1:dffts%nnr,iw)=ct%wapprox(iw)%vpvr_im(1:dffts%nnr,ir)

             enddo
             
             call DGEMM('T','N',ct%nf,ct%nf,dffts%nnr,1.d0,wterms,dffts%nnr,vprods,dffts%nnr,0.d0,gwmat_rr,ct%nf)
             call mp_sum(gwmat_rr, world_comm)
             call DGEMM('T','N',ct%nf,ct%nf,dffts%nnr,1.d0,wterms_im,dffts%nnr,vprods,dffts%nnr,0.d0,gwmat_ir,ct%nf)
             call mp_sum(gwmat_ir, world_comm)
             call DGEMM('T','N',ct%nf,ct%nf,dffts%nnr,1.d0,wterms,dffts%nnr,vprods_im,dffts%nnr,0.d0,gwmat_ri,ct%nf)
             call mp_sum(gwmat_ri, world_comm)
             call DGEMM('T','N',ct%nf,ct%nf,dffts%nnr,1.d0,wterms_im,dffts%nnr,vprods_im,dffts%nnr,0.d0,gwmat_ii,ct%nf)
             call mp_sum(gwmat_ii, world_comm)

       
             allocate(cmat(ct%nf,ct%nf),sigmac(2*ct%nf))
             do iw=1,ct%nf
                do jw=1,ct%nf
                   if(easy_grid_type/=1.and.easy_grid_type/=2) then
                      sca=(gwmat_rr(iw,jw)-gwmat_ii(iw,jw))/dble((dffts%nr1*dffts%nr2*dffts%nr3))/psif
                      sca1=(gwmat_ir(iw,jw)+gwmat_ri(iw,jw))/dble((dffts%nr1*dffts%nr2*dffts%nr3))/psif
                   else!if equally spaced grid we do not divide for psif as it can be zero
                      sca=(gwmat_rr(iw,jw)-gwmat_ii(iw,jw))/dble((dffts%nr1*dffts%nr2*dffts%nr3))
                      sca1=(gwmat_ir(iw,jw)+gwmat_ri(iw,jw))/dble((dffts%nr1*dffts%nr2*dffts%nr3))
                   endif
                   cmat(iw,jw)=cmplx(sca,sca1)
                   !if(ionode.and. kk==ct%iband .and. l_verbose) write(iun,*) sca,sca1
                enddo
                write(stdout,*) iw,dble(cmat(iw,ct%nf))
             enddo
             if(ionode.and. kk==ct%iband .and. l_verbose) close(iun)
          
             sigmac=0.d0
             do dw=-ct%nf+1,ct%nf
                dww=dw+ct%nf
                do jw=-ct%nf+1,ct%nf
                   iw=dw-jw
                   if(jw>-ct%nf .and. jw <ct%nf) then
                      if(jw>-ct%nf .and.jw <= 0) then 
                         jww=jw+ct%nf
                      elseif(jw>0 .and. jw< ct%nf) then
                         jww=-jw+ct%nf
                      endif
                      if(iw>-ct%nf .and. iw <= 0) then
                         iww=iw+ct%nf
                         sigmac(dww)=sigmac(dww)+cmat(jww,iww)
                      elseif(iw>0 .and. iw< ct%nf) then
                         iww=-iw+ct%nf
                         sigmac(dww)=sigmac(dww)+conjg(cmat(jww,iww))
                      endif
                      
                   endif
                enddo
             enddo
          
      
             emax=abs(aimag(ct%freq(1)))
             
             sigmac=sigmac*(emax*2.d0)/dble(2*ct%nf-1)*(1.d0)/(2.d0*3.1415926d0)

             if(easy_grid_type==0) then
  !MAXVAL VC case
                !sigmac=sigmac*psif*(ene_c)/psifvc!*(-1.d0)

  ! LONG CASE
                !sigmac=(sigmac-psif_xc/psif-psif_x/psif)+se%ene_xc(kk,is)-se%ene_x(kk,is)

             endif

             if(ionode .and. l_verbose ) then
                iun=find_free_unit()
                write(nfile,'(4i1)') ct%iband/1000,mod(ct%iband,1000)/100,mod(ct%iband,100)/10,mod(ct%iband,10)
                open(unit=iun,file=trim(prefix)//'sigmac.'//nfile//'.dat',status='unknown')
                do iw=1,2*ct%nf
                   write(iun,*) 2.d0*emax/dble(2*ct%nf)*iw-emax,dble(sigmac(iw)),aimag(sigmac(iw))
                enddo
                close(iun)
             endif

                
             !if(l_whole_s .or. easy_grid_type==1.or.easy_grid_type==2) then
             if(l_whole_s .or.easy_grid_type==2) then 
                do jj=-ct%nf+1,ct%nf-1
                   se%self(jj+ct%nf,nr_counter,kk,is)=sigmac(jj+ct%nf)
                enddo
             endif

             if(ionode) then
                if(kk==ct%iband) then
             
                   write(nfilel,'(5i1)') &
                  & ct%iband/10000,mod(ct%iband,10000)/1000,mod(ct%iband,1000)/100,mod(ct%iband,100)/10,mod(ct%iband,10)
                   if(is==1) then
                      iun1 = find_free_unit()
                      open( unit=iun1, file=trim(prefix)//'-'//'re_on_im'// nfilel, status='unknown',form='formatted')
                      iun2 = find_free_unit()
                      open( unit=iun2, file=trim(prefix)//'-'//'im_on_im'// nfilel, status='unknown',form='formatted')
                   else
                      iun1 = find_free_unit()
                      open( unit=iun1, file=trim(prefix)//'-'//'re_on_im2'// nfilel, status='unknown',form='formatted')
                      iun2 = find_free_unit()
                      open( unit=iun2, file=trim(prefix)//'-'//'im_on_im2'// nfilel, status='unknown',form='formatted')
                   endif
                   do jj=-ct%nf+1,ct%nf-1
                      freq_im= emax/dble(ct%nf)*jj
                      write(iun1,*) freq_im,0.d0,dble(sigmac(jj+ct%nf)),0.d0
                      write(iun2,*) freq_im,0.d0,aimag(sigmac(jj+ct%nf)),0.d0
                   enddo
                   close(iun1)
                   close(iun2)
                else
                   if(ct%iband>0 ) then
                      write(nfilel,'(5i1)') &
                        & ct%iband/10000,mod(ct%iband,10000)/1000,mod(ct%iband,1000)/100,mod(ct%iband,100)/10,mod(ct%iband,10)
                      write(nfilel2,'(5i1)') &
                           & kk/10000,mod(kk,10000)/1000,mod(kk,1000)/100,mod(kk,100)/10,mod(kk,10)
                      iun1 = find_free_unit()
                      open( unit=iun1, file=trim(prefix)//'-'//'re_on_im_off'// nfilel//'_'//nfilel2, &
                           &       status='unknown',form='formatted')
                      iun2 = find_free_unit()
                      open( unit=iun2, file=trim(prefix)//'-'//'im_on_im_off'// nfilel//'_'//nfilel2, &
                           &       status='unknown',form='formatted')
                   else
                      if(l_verbose) then
                         write(nfilex,'(4i1)') &
                              & ct%r_coord(1)/1000,mod(ct%r_coord(1),1000)/100,mod(ct%r_coord(1),100)/10, mod(ct%r_coord(1),10)
                         write(nfiley,'(4i1)') &
                         & ct%r_coord(2)/1000,mod(ct%r_coord(2),1000)/100,mod(ct%r_coord(2),100)/10, mod(ct%r_coord(2),10)
                         write(nfilez,'(4i1)') &
                         & ct%r_coord(3)/1000,mod(ct%r_coord(3),1000)/100,mod(ct%r_coord(3),100)/10, mod(ct%r_coord(3),10)
 
                         write(nfilel2,'(5i1)') &
                        & kk/10000,mod(kk,10000)/1000,mod(kk,1000)/100,mod(kk,100)/10,mod(kk,10)
                         iun1 = find_free_unit()
                         write(nfile_orb,'(4i1)') &
                              & kk/1000,mod(kk,1000)/100,mod(kk,100)/10, mod(kk,10)
                       
                         if(is==1) then
                            open( unit=iun1, file=trim(prefix)//'-gwl_orbital_1_'//nfile_orb//'/re_on_im_r_'&
                                 &//nfilex//'_'//nfiley//'_'//nfilez//'__'//nfilel2,status='unknown',form='formatted')
                         else
                            open( unit=iun1, file=trim(prefix)//'-gwl_orbital_2_'//nfile_orb//'/re_on_im_r_2'&
                                 &//nfilex//'_'//nfiley//'_'//nfilez//'__'//nfilel2, status='unknown',form='formatted')
                         endif
                     
                         iun2 = find_free_unit()
                     

                         if(is==1) then
                            open( unit=iun2, file=trim(prefix)//'-gwl_orbital_1_'//nfile_orb//'/im_on_im_r_'&
                                  &//nfilex//'_'//nfiley//'_'//nfilez//'__'//nfilel2,     status='unknown',form='formatted')
                         else
                            open( unit=iun2, file=trim(prefix)//'-gwl_orbital_2_'//nfile_orb//'/im_on_im_r_2'&
                                 &//nfilex//'_'//nfiley//'_'//nfilez//'__'//nfilel2,  status='unknown',form='formatted')
                         endif

                      endif
                   endif
                   if(.not.(easy_grid_type==1.or.easy_grid_type==2).or.l_verbose) then
                      do jj=-ct%nf+1,ct%nf-1
                         freq_im= emax/dble(ct%nf)*jj
                         write(iun1,*) freq_im,0.d0,dble(sigmac(jj+ct%nf)),0.d0
                         write(iun2,*) freq_im,0.d0,aimag(sigmac(jj+ct%nf)),0.d0
                      enddo
                      close(iun1)
                      close(iun2)
                   endif
                endif
                
             endif
      
      


             deallocate(cmat,sigmac)
          
          
             deallocate(vprods,vprods_im,wterms,wterms_im)
             deallocate(gwmat_rr,gwmat_ri,gwmat_ir,gwmat_ii)
             
         
          enddo
       enddo
       nr_counter=nr_counter+1
    enddo

   
    return

  END SUBROUTINE calculate_convergence

  SUBROUTINE average_self_energy(se)
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : prefix, tmp_dir
    USE wannier_gw

    IMPLICIT NONE

    TYPE(self_energy) :: se
    INTEGER :: kk,ir,is
    COMPLEX(kind=DP), ALLOCATABLE :: self_ave(:)
    REAL(kind=DP) :: den
    CHARACTER(5) :: nfile
    INTEGER :: iun1,iun2
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: nr_calculate

    allocate(self_ave(se%ngrid))
    
    do is=se%s_first_spin,se%s_last_spin
       do kk=se%s_first_state,se%s_last_state
          self_ave=0.d0
          den=0.d0
          nr_calculate=0
          do ir=1,se%nr
             if(abs(se%psif(ir,kk,is))>easy_psi_thrs) then
                den=den+se%psif(ir,kk,is)**2.d0
                self_ave(1:se%ngrid)=self_ave(1:se%ngrid)+se%self(1:se%ngrid,ir,kk,is)*se%psif(ir,kk,is)
                nr_calculate=nr_calculate+1
             endif
          enddo
          write(stdout,*) 'Averaging over r points',kk,den,dble(nr_calculate)
          if(easy_average_type==0) then 
             self_ave(1:se%ngrid)=self_ave(1:se%ngrid)/den
          else
             self_ave(1:se%ngrid)=self_ave(1:se%ngrid)/dble(nr_calculate)
          endif

          if(ionode) then
             write(nfile,'(5i1)') &
                  & kk/10000,mod(kk,10000)/1000,mod(kk,1000)/100,mod(kk,100)/10,mod(kk,10)
             if(is==1) then
                iun1 = find_free_unit()
                open( unit=iun1, file=trim(prefix)//'-'//'re_on_im'// nfile, status='unknown',form='formatted')
                iun2 = find_free_unit()
                open( unit=iun2, file=trim(prefix)//'-'//'im_on_im'// nfile, status='unknown',form='formatted')
             else
                iun1 = find_free_unit()
                open( unit=iun1, file=trim(prefix)//'-'//'re_on_im2'// nfile, status='unknown',form='formatted')
                iun2 = find_free_unit()
                open( unit=iun2, file=trim(prefix)//'-'//'im_on_im2'// nfile, status='unknown',form='formatted')
             endif
             do ir=1,se%ngrid
             
                write(iun1,*) se%freq(ir),0.d0,dble(self_ave(ir)),0.d0
                write(iun2,*) se%freq(ir),0.d0,aimag(self_ave(ir)),0.d0
             enddo
             close(iun1)
             close(iun2)
          endif


       enddo
    end do
    deallocate(self_ave)
  END SUBROUTINE average_self_energy



    SUBROUTINE solve_off_diagonal(se)
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : prefix, tmp_dir
    USE wannier_gw

    IMPLICIT NONE

    TYPE(self_energy),INTENT(inout) :: se!self_energies object
   

    INTEGER :: kk,ir,it
    COMPLEX(kind=DP), ALLOCATABLE :: self_diag(:,:)
    REAL(kind=DP) :: den
    CHARACTER(5) :: nfile
    INTEGER :: iun1,iun2
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP), ALLOCATABLE :: psi_mat(:,:),sigma_vec(:,:)
    REAL(kind=DP), ALLOCATABLE :: sigma_vec_im(:,:)
    INTEGER :: nn!matrix dimesion
    INTEGER, ALLOCATABLE :: ipiv(:)
    INTEGER :: info

    write(stdout,*) 'Routine solve_off_diagonal'

    nn=s_last_state-s_first_state+1
    if(nn /= se%nr) then
       call errore('solve_off_diagonal','mismatch bands/points ',1)  
    endif
    allocate(self_diag(se%ngrid,s_first_state:s_last_state))
    allocate(psi_mat(nn,nn))
    allocate(sigma_vec(nn,nn))
    allocate(sigma_vec_im(nn,nn))
    allocate(ipiv(nn))
    psi_mat=0.d0

    do it=1,se%ngrid
       do ir=1,nn
          do kk=s_first_state,s_last_state
            sigma_vec(ir,kk) =dble(se%self(it,ir,kk,1))*se%psif(ir,kk,1)
            sigma_vec_im(ir,kk) =aimag(se%self(it,ir,kk,1))*se%psif(ir,kk,1)
         enddo
      enddo
   
      psi_mat(1:nn,1:nn)=se%psif(1:nn,s_first_state:s_last_state,1)
      call DGESV(nn,nn,psi_mat,nn,ipiv,sigma_vec,nn,info)!psi_mat changed in exit
      if(info/=0) then
         call errore('solve_off_diagonal','DGESV error:',info)
      endif
      psi_mat(1:nn,1:nn)=se%psif(1:nn,s_first_state:s_last_state,1)
      call DGESV(nn,nn,psi_mat,nn,ipiv,sigma_vec_im,nn,info)!psi_mat changed in exit         
      if(info/=0) then
         call errore('solve_off_diagonal','DGESV error:',info)
      endif
      do kk=se%s_first_state,se%s_last_state
         self_diag(it,kk)=cmplx(sigma_vec(kk-s_first_state+1,kk-s_first_state+1),&
              sigma_vec_im(kk-s_first_state+1,kk-s_first_state+1))
         write(stdout,*) 'VECTOR',it,kk,sigma_vec(:,kk)
         write(stdout,*) 'VECTOR IM',it,kk,sigma_vec_im(:,kk)
      enddo
     
   enddo

   if(ionode) then
      do  kk=se%s_first_state,se%s_last_state
         write(nfile,'(5i1)') &
              & kk/10000,mod(kk,10000)/1000,mod(kk,1000)/100,mod(kk,100)/10,mod(kk,10)
         iun1 = find_free_unit()
         open( unit=iun1, file=trim(prefix)//'-'//'re_on_im'// nfile, status='unknown',form='formatted')
         iun2 = find_free_unit()
         open( unit=iun2, file=trim(prefix)//'-'//'im_on_im'// nfile, status='unknown',form='formatted')
         do it=1,se%ngrid
            
            write(iun1,*) se%freq(it),0.d0,dble(self_diag(it,kk)),0.d0
            write(iun2,*) se%freq(it),0.d0,aimag(self_diag(it,kk)),0.d0
         enddo
         close(iun1)
         close(iun2)
      enddo
   endif
   


   deallocate(psi_mat,sigma_vec)
   deallocate(sigma_vec_im)
   deallocate(ipiv,self_diag)

  END SUBROUTINE solve_off_diagonal




  SUBROUTINE free_convergence_tests(ct)
    IMPLICIT NONE
    TYPE(convergence_tests) :: ct
    INTEGER :: iw
    do iw=1,ct%nf
       call free_memory(ct%wapprox(iw))
       call free_memory(ct%g0(iw))
    enddo
       
    deallocate(ct%r)
    deallocate(ct%energy)
    CALL free_memory(ct%xf)
    CALL free_memory(ct%h0)
    deallocate(ct%wapprox)
    deallocate(ct%g0)
    RETURN
  END SUBROUTINE free_convergence_tests

  SUBROUTINE  free_vpv(c)
    IMPLICIT NONE
    TYPE(vpv) :: c
    if(associated(c%r)) deallocate(c%r)
    if(associated(c%vpvr)) deallocate(c%vpvr)
    if(associated(c%vpvr_im)) deallocate(c%vpvr_im)
    nullify(c%r,c%vpvr,c%vpvr_im)

    RETURN
  END SUBROUTINE free_vpv
   

  SUBROUTINE calculate_gzero(g0,freq,r,nr)

    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE g_psi_mod,            ONLY : h_diag, s_diag
    USE uspp_init,            ONLY : init_us_2

    IMPLICIT NONE
    TYPE(gzero) :: g0!green's function to be created and initialised
    REAL(kind=DP) :: freq!frequency to be calculated now just 0 implemented
    INTEGER :: r(nr)!list of r points 
    INTEGER :: nr!number of r points 

    INTEGER :: ir,ii,jj,ig,iv
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:,:),psi_g2(:,:)
    INTEGER :: kter
    LOGICAL :: lconv_root,lfirst
    REAL(kind=DP) :: anorm
    EXTERNAL :: hpsi_pw4gww2,cg_psi_pw4gww
    REAL(kind=DP), PARAMETER :: ethr=1d-10
    REAL(kind=DP), ALLOCATABLE :: et0(:,:)

   

    g0%freq=freq
    g0%nr=nr
    allocate(g0%r(nr))
    g0%r(1:nr)=r(1:nr)
    allocate(g0%gf(dffts%nnr,1,nr))

   
    allocate(et0(1,1))
    et0(1,1)=+(et(num_nbndv(1)+1,1)+et(num_nbndv(1),1))/2.d0

    call allocate_bec_type ( nkb, 1, becp)
    IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
    allocate (h_diag(npw, 1),s_diag(npw,1))
    allocate(psi_g2(npw,1),psi_g(npw,1))
    do ig = 1, npw
       g2kin (ig) = ( g (1,ig)**2 + g (2,ig)**2 + g (3,ig)**2 ) * tpiba2
    enddo
    h_diag=0.d0
    do ig = 1, npw
       h_diag(ig,1)=g2kin(ig)
    enddo
    !loop on points                                                                                                                                                                                                    
    do ir=1,nr

       psic(1:dfftp%nnr)=0.d0
       if(g0%r(ir) /= 0) psic(g0%r(ir))=1.d0*dble((dffts%nr1*dffts%nr2*dffts%nr3))
       CALL fwfft ('Wave', psic, dffts)
       psi_g(1:npw,1) = psic(dffts%nl(igk_k(1:npw,1)))
       psi_g2(1:npw,1)=psi_g(1:npw,1)
       call cgsolve_all_gamma (hpsi_pw4gww2,cg_psi_pw4gww,et0,psi_g,psi_g2, &
            h_diag,npw,npw,ethr,1,kter,lconv_root,anorm,1,1)

       psic(:)=(0.d0,0.d0)
       psic(dffts%nl(1:npw))  = psi_g2(1:npw,1)
       psic(dffts%nlm(1:npw)) = CONJG( psi_g2(1:npw,1) )
       CALL invfft ('Wave', psic, dffts)
       g0%gf(1:dffts%nnr,1,ir)= DBLE(psic(1:dffts%nnr))


     enddo

     deallocate(h_diag,s_diag,psi_g2,et0,psi_g)
     call deallocate_bec_type(becp)
      
    return
  END SUBROUTINE calculate_gzero


 SUBROUTINE calculate_gzero_complex(g0,freq,freq_im,r,nr,l_old,psi_old,l_lanczos,lc,ks_wfcs)
!lanczos AND multiple R points NOT IMPLEMENTED YET

    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE mp_global, ONLY : world_comm
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE g_psi_mod,            ONLY : h_diag, s_diag
    USE lanczos 
    USE lsda_mod,             ONLY : nspin, current_spin
    USE uspp_init,            ONLY : init_us_2
    
    IMPLICIT NONE

    TYPE(gzero) :: g0!green's function to be created and initialised 
    REAL(kind=DP) :: freq!frequency to be calculated real part
    REAL(kind=DP) :: freq_im!frequency to be calculated imaginary part   
    INTEGER :: r(nr)!list of r points  
    INTEGER :: nr!number of r points  
    LOGICAL :: l_old!if true in psi_old a previous trial solution is given
    COMPLEX(kind=DP) :: psi_old(npw) !in input: previous solution, in output: this solution in G space
    LOGICAL :: l_lanczos!if true useus lanczos algorithm for matrix inversion
    TYPE(lanczos_chain) :: lc(2)
    COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd,nspin)!KS wavefunctions in G space          


    INTEGER :: ir,ii,jj,ig,iv
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:,:),psi_g2(:,:),psi_g1(:,:),psi_g3(:,:)
    INTEGER :: kter
    LOGICAL :: lconv_root,lfirst
    REAL(kind=DP) :: anorm
    EXTERNAL :: hpsi_pw4gww2,cg_psi_pw4gww_square,hpsi_square
    REAL(kind=DP), PARAMETER :: ethr=1d-12!DEBUG ERA 10
    REAL(kind=DP), ALLOCATABLE :: et_fermi(:,:),et_zero(:,:),et_freq(:),et_freq_set(:)
    COMPLEX(kind=DP) :: freqc
    INTEGER :: is
    REAL(kind=DP) :: sca1,sca2
    REAL(kind=DP) :: smearing=0.0d0
    call start_clock('gzero_complex')

    g0%nr=nr
    allocate(g0%r(nr))
    g0%r(1:nr)=r(1:nr)
    allocate(g0%gf(dffts%nnr,nspin,nr))
    allocate(g0%gf_im(dffts%nnr,nspin,nr))
    g0%freq=freq
    g0%freq_im=freq_im
    freqc=cmplx(freq,freq_im)

    allocate(et_fermi(1,1),et_zero(1,1),et_freq(2),et_freq_set(2))

    et_zero=0.d0
    call allocate_bec_type ( nkb, 1, becp)
    IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
    allocate (h_diag(npw, 1),s_diag(npw,1))
    allocate(psi_g2(npw,1),psi_g(npw,1),psi_g1(npw,1),psi_g3(npw,1))
    do ig = 1, npw
       g2kin (ig) = ( g (1,ig)**2 + g (2,ig)**2 + g (3,ig)**2 ) * tpiba2
    enddo

    do is=1,nspin
       current_spin=is
       et_fermi(1,1)=(et(num_nbndv(is)+1,is)+et(num_nbndv(is),is))/2.d0
!DEBUG INIZIO
       !et_fermi(1,1)=(et(num_nbndv(1)+1,1)+et(num_nbndv(1),1))/2.d0
      ! write(stdout,*) 'DEBUG FERMI ENERGY',  et_fermi(1,1)*13.606, num_nbndv(is),is
       
!DEBUG FINE

       do ig = 1, npw
          h_diag(ig,1)=1.d0*(g2kin(ig)-(freq+et_fermi(1,1)))**2.d0+freq_im**2.d0
       enddo

    !loop on points 
!set frequencies

       et_freq(1)=g0%freq+et_fermi(1,1)+smearing
       et_freq(2)=g0%freq_im
       et_freq_set(1)=et_fermi(1,1)
       et_freq_set(2)=0.d0
       call hpsi_square( npw,psi_g,psi_g1,et_freq,-1,2)
       call hpsi_square( npw,psi_g,psi_g1,et_freq,-1,2)
       !call hpsi_square( npw,psi_g,psi_g1,et_freq,-2,2)!DEBUG
     
    
       do ir=1,nr

          psic(1:dfftp%nnr)=0.d0
          if(g0%r(ir) /= 0) psic(g0%r(ir))=1.d0!*dble((dffts%nr1*dffts%nr2*dffts%nr3))
          CALL fwfft ('Wave', psic, dffts)
          psi_g(1:npw,1) = psic(dffts%nl(igk_k(1:npw,1)))
          !call pv_operator(psi_g(1,1),is,ks_wfcs,.true.)
          if(lc(is)%l_first .and. l_lanczos) then
             !call create_lanczos_chain(lc,hpsi_pw4gww2,1,n_self_lanczos,psi_g(1:npw,1),.true.,et_freq_set,0)
             call create_krylov(lc(is),hpsi_pw4gww2,1,n_self_lanczos,psi_g(1,1),et_freq_set,is)!,ks_wfcs)
             lc(is)%l_first=.false.
          endif
          if(l_lanczos) then
             !call solve_lanczos(lc,psi_g,freqc,psi_g2,psi_g3,.false.,psi_g,psi_g,hpsi_pw4gww2,et_freq_set)
             call solve_krylov(lc(is),psi_g,freqc,psi_g2,psi_g3)
         else
!first calculate the real part
             call h_psi( npw, npw, 1, psi_g, psi_g1 )
             !call pv_operator(psi_g1(1,1),is,ks_wfcs,.true.)
             psi_g1(1:npw,1)=psi_g1(1:npw,1)-(freq+et_fermi(1,1))*psi_g(1:npw,1)
          
             if(.not. l_old) then
                psi_g2(1:npw,1)=psi_g1(1:npw,1)
             else
                psi_g2(1:npw,1)=psi_old(1:npw)
             endif
          
       
             call cgsolve_all_gamma(hpsi_square,cg_psi_pw4gww_square,et_zero,psi_g1,psi_g2, &
                  h_diag,npw,npw,ethr,1,kter,lconv_root,anorm,1,1)
             
       
          endif

          !check for quality of solution
          !call h_psi( npw, npw, 1, psi_g2, psi_g3 )
          !psi_g3(1:npw,1)=psi_g3(1:npw,1)-(freq+et_fermi(1,1))*psi_g2(1:npw,1)
          !psi_g3(1:npw,1)=psi_g3(1:npw,1)-psi_g(1:npw,1)
          !sca1=0.d0
          !do ig=1,npw
          !   sca1=sca1+2.0*dble(psi_g3(ig,1)*conjg(psi_g3(ig,1)))
          !enddo
          !if(gstart==2) sca1=sca1 -dble(psi_g3(1,1)*conjg(psi_g3(1,1)))
          !call mp_sum(sca1,world_comm)
          !sca2=0.d0
          !do ig=1,npw
          !   sca2=sca2+2.0*dble(psi_g(ig,1)*conjg(psi_g(ig,1)))
          !enddo
          !if(gstart==2) sca2=sca2 -dble(psi_g(1,1)*conjg(psi_g(1,1)))
          !call mp_sum(sca2,world_comm)
          !write(stdout,*) 'SOLUTION QUALITY', sca1,sca2,sca1/sca2
          !!
          !call pc_operator2(psi_g2,is,ks_wfcs,.true.)
          psi_old(1:npw)=psi_g2(1:npw,1)
          psic(:)=(0.d0,0.d0)
          psic(dffts%nl(1:npw))  = psi_g2(1:npw,1)
          psic(dffts%nlm(1:npw)) = CONJG( psi_g2(1:npw,1) )
          CALL invfft ('Wave', psic, dffts)
          g0%gf(1:dffts%nnr,is,ir)= DBLE(psic(1:dffts%nnr))*dble((dffts%nr1*dffts%nr2*dffts%nr3))

!then the imaginary part
          if(.not.l_lanczos) then
             psi_g1(1:npw,1)=freq_im*psi_g(1:npw,1)
             psi_g3(1:npw,1)=psi_g1(1:npw,1)
             call cgsolve_all_gamma(hpsi_square,cg_psi_pw4gww_square,et_zero,psi_g1,psi_g3, &
                  h_diag,npw,npw,ethr,1,kter,lconv_root,anorm,1,1)
          endif

          !call pc_operator2(psi_g3,is,ks_wfcs,.true.)
          psic(:)=(0.d0,0.d0)
          psic(dffts%nl(1:npw))  = psi_g3(1:npw,1)
          psic(dffts%nlm(1:npw)) = CONJG( psi_g3(1:npw,1) )
          CALL invfft ('Wave', psic, dffts)
          g0%gf_im(1:dffts%nnr,is,ir)= DBLE(psic(1:dffts%nnr))*dble((dffts%nr1*dffts%nr2*dffts%nr3))


       enddo
    enddo
!DEBUG INIZIO
!     g0%gf_im(1:dffts%nnr,2,1)=  g0%gf_im(1:dffts%nnr,1,1)
!     g0%gf(1:dffts%nnr,2,1)=  g0%gf(1:dffts%nnr,1,1)
!DEBUG FINE
    deallocate(h_diag,s_diag,psi_g2,et_fermi,psi_g,psi_g1,psi_g3)
    deallocate(et_zero,et_freq,et_freq_set)
    call deallocate_bec_type(becp)
    call stop_clock('gzero_complex')
    return
  END SUBROUTINE calculate_gzero_complex



  SUBROUTINE calculate_vpv_complex(c,freq,freq_im,r,nr,v_states,l_w,thrs,lc,&
                                  &  head, wing, alpha,ks_wfcs, l_save,phi_save)


    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE g_psi_mod,            ONLY : h_diag, s_diag
    USE mp_world,  ONLY : world_comm
    USE lanczos
    USE lsda_mod,    ONLY : nspin,  current_spin
    USE io_files,             ONLY :  prefix, diropn
    USE gvecw,     ONLY : ecutwfc
    USE uspp_init, ONLY : init_us_2


    IMPLICIT NONE
    TYPE(vpv) :: c!vpv structure to be created and initialised
    REAL(kind=DP) :: freq!frequency to be calculated now just 0 implemenmted , real part
    REAL(kind=DP) :: freq_im!frequency to be calculated now just 0 implemenmted, imaginary part
    INTEGER :: r(nr)!list of r points 
    INTEGER :: nr!number of r points 
    REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,num_nbnds,nspin)!valence states in real space
    LOGICAL :: l_w!if true calculate W
    REAL(kind=DP) :: thrs!threshold for W calculations, now only for imaginary frequencies
    TYPE(lanczos_chain) :: lc(2) !structure for lanczos chains
    REAL(kind=DP) :: head! in case of extended systems, the head of the RPA symm. diel. matrix
    COMPLEX(kind=DP) :: wing(npw)!wing of the symmetric irriducible dielectric matrix
    REAL(kind=DP) :: alpha !initial parameter for Richardson algorithm
    COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd,nspin)!KS wave functions in real space
    COMPLEX(kind=DP), INTENT(inout) :: phi_save(npw)!previous solutionn to the problem
    LOGICAL, INTENT(in) :: l_save!if true in phi_save previous solution to the problem
    
    INTEGER :: ir,ii,jj,ig,iv
    REAL(kind=DP), ALLOCATABLE :: fac(:)
    REAL(kind=DP):: qq
    COMPLEX(kind=DP), ALLOCATABLE :: verre(:),verre0(:)
    INTEGER :: numv(2)
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g2(:,:,:),psi_g1(:,:),psi_g3(:,:),psi_g4(:,:,:),psi_g5(:,:)
    REAL(kind=DP), ALLOCATABLE :: psi_r(:,:),psi_v(:),psi_v0(:)
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:,:,:)
    INTEGER :: kter
    LOGICAL :: lconv_root,lfirst
    REAL(kind=DP) :: anorm
    EXTERNAL :: hpsi_square,cg_psi_pw4gww_square,hpsi_pw4gww
    COMPLEX(kind=DP), ALLOCATABLE :: opsi(:)
    REAL(kind=DP), PARAMETER :: ethr=1d-10
    REAL(kind=DP) et_freq(2)
    INTEGER :: it
    INTEGER :: max_iter=50!maximum number of iterations (15)
    REAL(kind=DP) :: sca
    COMPLEX(kind=DP) :: freqc
    COMPLEX(kind=DP), ALLOCATABLE :: verreg(:,:,:)
    REAL(kind=DP) :: erre,erre2
    REAL(kind=DP), ALLOCATABLE :: et_scissor(:,:)
    REAL(kind=DP) :: scissor_w=0.d0!DEBUG ATTENZIONE ERA 0.d0
    REAL(kind=DP) :: wing_fact
    COMPLEX(kind=DP), ALLOCATABLE :: phi0(:),phi1(:)
    REAL(kind=DP) :: previous_norm=10.d0
    REAL(kind=DP) :: testa(3)
    REAL(kind=DP), ALLOCATABLE :: norms(:,:),norms_lan(:,:)
    TYPE(lanczos_chain) :: lc_save(2)
    REAL(kind=DP) :: energy0, energy1, energy,energy_old,mod_force, mod_force_old
    COMPLEX(kind=DP), ALLOCATABLE :: force(:),phi_old(:)
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g_phi(:,:)
    REAL(kind=DP) :: alpha_old,par_a,par_b=1.d0,par_c=1.d0, alpha_new
    LOGICAL :: l_updated
    COMPLEX(kind=DP), ALLOCATABLE :: h_diag2(:,:), s_diag2(:,:)
    INTEGER :: is

    LOGICAL :: l_update_basis_w, l_restart
    COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)
    INTEGER :: iungprod
    LOGICAL       :: exst
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP), ALLOCATABLE :: pbr(:),wpbr(:)
    INTEGER :: ifreq, ngm_max
    REAL(kind=DP), ALLOCATABLE :: omat(:,:)

    if(freq_im<0.2d0) then
       l_update_basis_w=.false.
       max_iter=100
    else
       l_update_basis_w=.false.
    endif

    l_restart=.false.
    l_wing=.true.
    previous_norm=10.d0
    energy_old=10.d10
    mod_force_old=10.d0
    alpha_old=alpha
    l_updated=.false.
  
    call start_clock('vpv_complex')

    call initialize_memory(lc_save(1))
    call initialize_memory(lc_save(2))
   
    if(.not.lc(1)%l_first) call copy(lc(1),lc_save(1))
    if(.not.lc(2)%l_first .and. nspin==2) call copy(lc(2),lc_save(2))
  
    numv(1)=num_nbndv(1)
    if(nspin==2) numv(2)=num_nbndv(2)
    allocate(norms(numv(1),nspin),norms_lan(numv(1),nspin))
    allocate(force(npw),phi_old(npw))
    freqc=cmplx(freq,freq_im)

    c%freq=freq
    c%freq_im=freq_im
    allocate(c%r(nr))
    c%r(1:nr)=r(1:nr)
    allocate(c%vpvr(dffts%nnr,nr))
    allocate(c%vpvr_im(dffts%nnr,nr))
    
    if(.not.l_w) max_iter=1

 
    allocate(phi0(npw),phi1(npw))
    allocate(psi_g_phi(npw,numv(1)))


    c%vpvr_im=0.d0

    et_freq(1)=freq
    et_freq(2)=freq_im

!allocate

    allocate(fac(ngm))
    if(l_truncated_coulomb) then
       do ig=1,ngm
          qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
          if (qq > 1.d-8) then
             fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
          else
             fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
          endif
       enddo

       if(gstart==2) fac(1)=0.d0  
       fac(:)=fac(:)/omega
    else
       fac(1:npw)=vg_q(1:npw)
    endif
   
    if(.not.l_truncated_coulomb .and. l_wing) then
       wing(gstart:npw)=wing(gstart:npw)*sqrt( fac(gstart:npw)*omega/e2/fpi)
    endif
    allocate(verre(npw),verre0(npw))
    allocate(verreg(npw,numv(1),nspin))
  
    call allocate_bec_type ( nkb, numv(1), becp)
    IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
    allocate (h_diag(npw, numv(1)),s_diag(npw,numv(1)))
    allocate(psi_g2(npw,numv(1),nspin))
    do ig = 1, npw
       g2kin (ig) = ( g (1,ig)**2 + g (2,ig)**2 + g (3,ig)**2 ) * tpiba2
    enddo
    h_diag=0.d0
    do iv = 1, numv(1)
       do ig = 1, npw
          h_diag(ig,iv)=(g2kin(ig)-(freq+et(iv,1)))**2.d0+freq_im**2.d0
       enddo
    enddo
    allocate(psi_r(dffts%nnr,2),psi_v(dffts%nnr),psi_v0(dffts%nnr))
    allocate(psi_g(npw,numv(1),nspin),psi_g1(npw,numv(1)),psi_g3(npw,numv(1)))
    allocate(opsi(npw))
    allocate(psi_g4(npw,numv(1),nspin),psi_g5(npw,numv(1)))
    if(nspin==2) then
       allocate(h_diag2(npw, numv(2)),s_diag2(npw,numv(2)))
       h_diag2=0.d0
       do iv = 1, numv(2)
          do ig = 1, npw
             h_diag2(ig,iv)=(g2kin(ig)-(freq+et(iv,2)))**2.d0+freq_im**2.d0
          enddo
       enddo
    endif



    allocate(et_scissor(numv(1),1))
    et_scissor(1:numv(1),1)=et(1:numv(1),1)-scissor_w

!loop on points 
    do ir=1,nr
     
       it=0
       do while(it < max_iter) 
          it=it+1
!first real part
       !product with v
          if(it==1) alpha=1.0
          if(it==1.or.l_easy_update_basis_w) then
             psic(1:dfftp%nnr)=0.d0
             if(c%r(ir) /= 0) psic(c%r(ir))=1.!*dble((dffts%nr1*dffts%nr2*dffts%nr3))
             CALL fwfft ('Wave', psic, dffts)
             verre(1:npw) = psic(dffts%nl(igk_k(1:npw,1)))
             if(l_easy_dielectric_constant) verre(gstart:npw)=0.d0
             if(gstart==2) then
                erre=verre(1)
                wing_fact=sqrt(fac(1))
             else
                erre=0.d0
                wing_fact=0.d0
             endif
             call mp_sum(erre, world_comm)
             call mp_sum(wing_fact, world_comm)
             verre0(1:npw)=verre(1:npw)
             verre(1:npw)=sqrt(fac(1:npw))*verre(1:npw)

             
             do is=1,nspin
                do iv=1,numv(is)
                   verreg(1:npw,iv,is)=verre(1:npw)
                enddo
             enddo

             if(it==1) then
                if(.not.l_save) then
                   psic(:)=(0.d0,0.d0)
                   psic(dffts%nl(1:npw))  = verre(1:npw)
                   psic(dffts%nlm(1:npw)) = CONJG( verre(1:npw))
                   phi0(1:npw)=verre0(1:npw)
                else
                   phi0(1:npw)=phi_save(1:npw)
                endif
             else
                psic(:)=(0.d0,0.d0)
                psic(dffts%nl(1:npw))  = phi0(1:npw)*sqrt(fac(1:npw))
                psic(dffts%nlm(1:npw)) = CONJG( phi0(1:npw)*sqrt(fac(1:npw))  )
              
             endif
             CALL invfft ('Wave', psic, dffts)
             psi_v0(1:dffts%nnr)= DBLE(psic(1:dffts%nnr))
          endif

          if(lc(1)%l_first.or. (l_easy_update_basis_w.and.it>1) .or. (l_update_basis_w .and. it==1) ) then!DEBUG ERA 1
             psi_v(1:dffts%nnr)=psi_v0(1:dffts%nnr)
             do is=1,nspin
                do iv=1,numv(is),2
!!product with psi_v 
                   if(iv/=numv(1)) then
                      psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv,is)
                      psi_r(1:dffts%nnr,2)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv+1,is)
                   else
                      psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv,is)
                   endif
!!fourier transfrm to G   
                   if(iv/=numv(is)) then
                      psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),psi_r(1:dffts%nnr,2))
                   else
                      psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),0.d0)
                   endif
                   CALL fwfft ('Wave', psic, dffts)
                   if(iv/=numv(is)) then
                      psi_g(1:npw,iv,is)=0.5d0*(psic(dffts%nl(1:npw))+conjg(psic(dffts%nlm(1:npw))))
                      psi_g(1:npw,iv+1,is)=(0.d0,-0.5d0)*(psic(dffts%nl(1:npw))-conjg(psic(dffts%nlm(1:npw))))
                      if(gstart==2) psi_g(1,iv,is)=dble(psi_g(1,iv,is))
                      if(gstart==2) psi_g(1,iv+1,is)=dble(psi_g(1,iv+1,is))
                   else
                      psi_g(1:npw,iv,is)=psic(dffts%nl(1:npw))
                      if(gstart==2) psi_g(1,iv,is)=dble(psi_g(1,iv,is))
                   endif
              
!!project on conduction manifold   
                   call start_clock('opsi_pc')
                   if(iv/=numv(is)) then
                      call pc_operator2(psi_g(:,iv,is),is,ks_wfcs,.true.)!DEBUG ERA FALSE
                      call pc_operator2(psi_g(:,iv+1,is),is,ks_wfcs,.true.)
                   else
                      call pc_operator2(psi_g(:,iv,is),is,ks_wfcs,.true.)
                   endif
                   call stop_clock('opsi_pc')
                enddo

             enddo

             norms=0.d0
             do is=1,nspin
                do iv=1,numv(is)
                   do ig=1,npw
                      norms(iv,is)=norms(iv,is)+2.d0*dble(conjg(psi_g(ig,iv,is))*psi_g(ig,iv,is))
                   enddo
                   if(gstart==2) norms(iv,is)=norms(iv,is)-dble(conjg(psi_g(1,iv,is))*psi_g(1,iv,is))
                enddo
             enddo
             call mp_sum(norms,world_comm)
           
             if(it>1) then
                do is=1,nspin
                   do iv=1,numv(is)
                      verreg(1:npw,iv,is)=phi0(1:npw)*sqrt(fac(1:npw))
                   enddo
                   if(gstart==2) verreg(1,:,is)=0.d0 
                   call norms_lanczos(lc(is),verreg(1,1,is),norms_lan(1,is))
                enddo
             endif

           
!if required start lanczos chains
            ! if(it>1 .and. l_verbose) write(stdout,*) 'SCARTO MAX', &
            !      &maxval(abs(norms(1:numv(1),1)-norms_lan(1:numv(1),1))/norms(1:numv(1),1)),alpha
            ! if(it>1 .and. l_verbose .and. nspin==2) write(stdout,*) 'SCARTO MAX2', &
            !      &maxval(abs(norms(1:numv(2),2)-norms_lan(1:numv(2),2))/norms(1:numv(2),2)),alpha
             if(it==1 ) then
                
                call start_clock('vpv_lanczos')

                do is=1,nspin
                   et_scissor(1:numv(is),1)=et(1:numv(is),is)-scissor_w
                   current_spin=is
                   evc(1:npw, 1:numv(is))=ks_wfcs(1:npw,1:numv(is),is)
                   call free_memory(lc(is))
!
                   
                   call create_krylov(lc(is),hpsi_pw4gww,numv(is),n_pola_lanczos,psi_g(1,1,is),et_scissor,is)!,ks_wfcs)
                   call product_chain(lc(is),v_states(1,1,is),0)
                   
                   lc(is)%l_first=.false.
                             
                   if(it==1) call copy(lc(is),lc_save(is))
                
                enddo
                
                call stop_clock('vpv_lanczos')
             else 
                do is=1,nspin
                   evc(1:npw, 1:numv(is))=ks_wfcs(1:npw,1:numv(is),is)
                   current_spin=is
                   et_scissor(1:numv(is),1)=et(1:numv(is),is)-scissor_w
                   if((maxval(abs(norms(1:numv(is),is)-norms_lan(1:numv(is),is))/norms(1:numv(is),is)) > easy_w_update_lanczos &
                        &.and. maxval (lc(is)%ns_chain) > 5*n_pola_lanczos) .or. l_update_basis_w) then
!ATTENZIONE ERA 5
                      call free_memory(lc(is))
                      
                      call create_krylov(lc(is),hpsi_pw4gww,numv(is),n_pola_lanczos,psi_g(1,1,is),et_scissor,is)!,ks_wfcs)
                      call product_chain(lc(is),v_states(1,1,is),0)
                      lc(is)%l_first=.false.
                   elseif(.not.l_updated) then
                   
                      do iv=1,numv(is)
                         if(abs(norms(iv,is)-norms_lan(iv,is))/norms(iv,is) > easy_w_update_lanczos) then
                         !   if( l_verbose) write(stdout,*) 'UPDATING LANCZOS', iv,is!DEBUG_NO
                            call update_krylov(lc(is),hpsi_pw4gww,psi_g(1,1,is),et_scissor, 2,iv,is,ks_wfcs)
                            call product_chain_krylov(lc(is),v_states(1,1,is),2,iv)
                         endif
                      enddo
                   endif
                enddo
             endif
          elseif (l_easy_update_basis_w.and.it==1) then
             !call copy(lc_save,lc)
          endif
        
          
          call start_clock('vpv_lanczos')
          
          do is=1,nspin
             do iv=1,numv(is)
                verreg(1:npw,iv,is)=phi0(1:npw)*sqrt(fac(1:npw))
             enddo
             if(gstart==2 ) verreg(1,:,is)=0.d0!ATTENZIONE
          enddo
          !call product_wfc(v_states,numv,psi_g,psi_g_phi)
          
          do is=1,nspin
             if(lc(1)%l_krylov) then
                call solve_krylov(lc(is),verreg(1,1,is),freqc,psi_g2(1,1,is),psi_g4(1,1,is))
             else
               ! call solve_lanczos(lc(is),verreg(1,1,is),freqc,psi_g2,psi_g4,.false.,psi_g,psi_g_phi,hpsi_pw4gww,et_scissor)
             endif
             if(gstart==2) psi_g2(1,1:numv(1),is)=0.d0
          enddo
          if(nspin==1) then
             psi_g2=4.d0*psi_g2
          else
             psi_g2=2.d0*psi_g2
          endif
          psi_g4=0.d0
          call stop_clock('vpv_lanczos')
                
          if(gstart==2) phi0(1)=dble(phi0(1))
          
          
          phi1(1:npw)=alpha*(-phi0(1:npw)+verre(1:npw))
          


          energy0=0.d0
          do ig=1,npw
             energy0=energy0+2.d0*dble(conjg(phi0(ig)*verre(ig)))
          enddo
          if(gstart==2) energy0=energy0-dble(conjg(phi0(1)*verre(1)))
          call mp_sum(energy0,world_comm)
          energy0=energy0*2.d0


          do is=1,nspin
             do iv=1,numv(is)
                phi1(1:npw)=phi1(1:npw)+alpha*(-sqrt(fac(1:npw))*psi_g2(1:npw,iv,is))
             enddo
          enddo
          
          if(gstart==2 .and. l_truncated_coulomb) phi1(1)=0.d0 
          
          if(.not.l_truncated_coulomb.and.l_wing) then
!add wing term
             if(gstart==2) then
                wing_fact=phi0(1)
             else
                wing_fact=0.
             endif
             erre2=0.d0
             do ig=gstart,npw
                erre2=erre2+2.d0*dble(conjg(wing(ig))*phi0(ig))
             enddo

             call mp_sum(wing_fact,world_comm)
             call mp_sum(erre2,world_comm)
             phi1(gstart:npw)=phi1(gstart:npw)-wing(gstart:npw)*wing_fact*alpha
             if(gstart==2) phi1(1)=phi1(1)+phi0(1)*(-head)*alpha!OK
             if(gstart==2) phi1(1)=phi1(1)-erre2*alpha
             if(gstart==2) then
                testa(1)=phi1(1)+phi0(1)
                testa(2)=verre(1)  
                testa(3)=phi0(1)
             else
                testa=0.d0
             endif
             call mp_sum(testa,world_comm)
           !  if(l_verbose) write(stdout,*)' TESTA', testa(1:3),head,erre2,wing_fact!DEBUG_NO
          endif
       
          force(1:npw)=phi1(1:npw)/alpha
          phi1(1:npw)=phi1(1:npw)+phi0(1:npw)
          mod_force=0.d0
          do ig=1,npw
             mod_force=mod_force+2.d0*dble(conjg(force(ig))*(force(ig)))
          enddo
          if(gstart==2) mod_force=mod_force-dble(conjg(force(1))*(force(1)))
          call mp_sum(mod_force,world_comm)

          energy1=0.d0
          do ig=1,npw
             energy1=energy1+2.d0*dble(conjg(phi0(ig))*(-(force(ig)+verre(ig))))
          enddo
          if(gstart==2) energy1=energy1-dble(conjg(phi0(1))*(-(force(1)+verre(1))))
          call mp_sum(energy1,world_comm)

          !energy=-energy0+energy1 
          energy=energy1

        !  if(l_verbose) write(stdout,*) 'ENERGIA :',energy0,energy1,energy!DEBUG_NO

          psic(:)=(0.d0,0.d0)
          psic(dffts%nl(1:npw))  = phi1(1:npw)*sqrt(fac(1:npw))-verre(1:npw)*sqrt(fac(1:npw))
          psic(dffts%nlm(1:npw)) = CONJG( phi1(1:npw)*sqrt(fac(1:npw))-verre(1:npw)*sqrt(fac(1:npw)))
          CALL invfft ('Wave', psic, dffts)
          c%vpvr(1:dffts%nnr,ir)= DBLE(psic(1:dffts%nnr))



         
       !now imaginary part
!TO BE DONE
    
          
          c%vpvr_im(1:dffts%nnr,ir)=0.d0

          
          
          sca=0.d0
          do ig=1,npw
             sca=sca+2.d0*dble((phi0(ig)-phi1(ig))*conjg(phi0(ig)-phi1(ig)))
          enddo
          if(gstart==2) sca=sca-dble((phi0(1)-phi1(1))*conjg(phi0(1)-phi1(1)))
          call mp_sum(sca,world_comm)
          sca=sca*dble((dffts%nr1*dffts%nr2*dffts%nr3))
          
        !  if(l_verbose) write(stdout,*) 'Convergence criteria:', mod_force,it,alpha!DEBUG_NO
          if(mod_force< thrs) then
             if(l_verbose) write(stdout,*) 'Converged with interations :', it
             write(stdout,*) 'Converged with interations :', it, mod_force
             phi_save(1:npw)=phi1(1:npw)
             exit
          !elseif(mod_force>mod_force_old .or. abs(mod_force_old-mod_force)/mod_force_old < 0.1) then
          elseif(mod_force>mod_force_old) then
           if(l_verbose) write(stdout,*) 'ALPHA UPDATE'
        !elseif(energy>1d10) then   
             alpha=alpha*easy_w_update_alpha
             phi0(1:npw)=phi_old(1:npw)
             l_updated=.true.
             mod_force_old=10.d10
             energy_old=10.d0
!DEBUG INIZIO
!             phi_old(1:npw)=phi0(1:npw)
!             phi0(1:npw)=phi1(1:npw)
!DEBUG FINE

       !      (sca>previous_norm) then
       !      
       !      alpha=alpha/2.d0
       !      it=0!START FROM BEGINNING
       !      previous_norm=10.0
       !      if(l_verbose) write(stdout,*) 'NEW ALPHA FOR RICHARDSON', alpha
       !   elseif((previous_norm-sca)/previous_norm<0.25) then
       !      write(stdout,*) 'START ANEW 1'
       !       alpha=alpha/2.d0
       !       it=0
       !       previous_norm=10.0
       !       if(l_verbose) write(stdout,*) 'NEW ALPHA FOR RICHARDSON', alpha
          else
             l_updated=.false.
             if(energy_old<1d0) then
                par_a=energy_old
                par_b=-mod_force_old*2.d0
                par_c=(energy-par_a-par_b*alpha_old)/alpha_old**2.d0
               ! if(l_verbose) write(stdout,*) 'NEW ALPHA FROM MINIMIZATION', -par_b/(2.d0*par_c)!DEBUG_NO
             endif
             previous_norm=sca
            
             mod_force_old=mod_force
             alpha_old=alpha
             phi_old(1:npw)=phi0(1:npw)
             phi0(1:npw)=phi1(1:npw)
             alpha_new= -par_b/(2.d0*par_c)
            if(energy_old<1d0 .and. alpha_new > 0.1d0  .and. alpha_new < 2.d0 .and. mod_force>1d-14)  alpha=alpha_new !DEBUG
            energy_old=energy
         endif
         !if(l_verbose.and. .not.l_truncated_coulomb ) write(stdout,*) 'costante dielettrica',it, 1./( testa(1)/testa(2)),alpha
         if(it==max_iter .and. .not.l_restart) then
             if(l_verbose) write(stdout,*) 'NOT CONVERGED, RETRYING:', freq_im 
             l_update_basis_w=.true.
             it=0
             previous_norm=10.d0
             energy_old=10.d10
             mod_force_old=10.d0
             alpha_old=alpha
             l_updated=.false.
             l_restart=.true.
          endif
       enddo
       if(it==max_iter) write(stdout,*) 'NOT CONVERGED:' ,it,mod_force
       if(l_easy_dielectric_constant .and. .not.l_truncated_coulomb) then
          write(stdout,*) 'COSTANTE DIELETTRICA',1./( testa(1)/testa(2))
       endif 

       c%vpvr(1:dffts%nnr,ir)= c%vpvr(1:dffts%nnr,ir)*dble((dffts%nr1*dffts%nr2*dffts%nr3))
       c%vpvr_im(1:dffts%nnr,ir)= c%vpvr_im(1:dffts%nnr,ir)*dble((dffts%nr1*dffts%nr2*dffts%nr3))

       !add head term
       if(.not.l_truncated_coulomb.and.(.not.l_wing)) then
          psic(:)=(0.d0,0.d0)
          if(gstart==2) psic(dffts%nl(1))  = fac(1)*erre*(1.d0/(head+1.d0)-1.d0)
          if(gstart==2) write(stdout,*) 'DEBUG PSIC', fac(1)*erre*(1.d0/(head+1.d0)-1.d0)
          if(l_verbose) write(stdout,* ) 'HEAD ', head
          CALL invfft ('Wave', psic, dffts)
          c%vpvr(1:dffts%nnr,ir)= c%vpvr(1:dffts%nnr,ir)+ &
               & DBLE(psic(1:dffts%nnr))*dble((dffts%nr1*dffts%nr2*dffts%nr3))
                    
       endif
    
    enddo
    
    !call copy(lc_save,lc)

    deallocate(h_diag,s_diag,opsi)
    deallocate(fac,verre,verre0,psi_g2,psi_r,psi_g,psi_v,psi_v0,psi_g1,psi_g3)
    deallocate(psi_g4,psi_g5)
    call deallocate_bec_type(becp)
    deallocate(verreg)
    deallocate(phi0,phi1)
    deallocate(norms,norms_lan)
    call free_memory(lc_save(1))
    call free_memory(lc_save(2))
    deallocate(force,phi_old)
    deallocate(psi_g_phi)
    deallocate(et_scissor)
    if(nspin==2) then
       deallocate(h_diag2,s_diag2)
    endif

    call stop_clock('vpv_complex')
    return

  
  END SUBROUTINE calculate_vpv_complex



  SUBROUTINE calculate_x(xf,r,nr,v_states)


    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE g_psi_mod,            ONLY : h_diag, s_diag
    USE lsda_mod,             ONLY : lsda, nspin,current_spin,isk

    IMPLICIT NONE
    TYPE(exchange) :: xf!exchange structure to be created and initialised   
    INTEGER :: r(nr)!list of r points  
    INTEGER :: nr!number of r points 
    REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,num_nbnds,nspin)!valence states in real space 

    REAL(kind=DP), ALLOCATABLE :: fac(:)
    REAL(kind=DP):: qq
    INTEGER :: ir,ii,jj,ig,iv,is
    REAL(kind=DP), ALLOCATABLE :: psi_r(:,:),psi_v(:)
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:,:)
    INTEGER :: numv

   
    xf%nr=nr
    allocate(xf%r(nr))
    xf%r(1:nr)=r(1:nr)
    allocate(xf%x(dffts%nnr,nr,nspin))

    xf%x=0.d0


    allocate(fac(ngm))
    allocate(psi_r(dffts%nnr,2),psi_v(dffts%nnr))
    allocate(psi_g(npw,num_nbndv(1)))

    if(l_truncated_coulomb) then
       do ig=1,ngm
          qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
          if (qq > 1.d-8) then
             fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
          else
             fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
          endif
       enddo
       fac(:)=fac(:)/omega
    else
       fac(:)=0.d0
       fac(1:npw)=vg_q(1:npw)
    endif
    do is=1,nspin
       numv=num_nbndv(is)
       do ir=1,nr
         !product with v
         
          psi_v(1:dfftp%nnr)=0.d0
          if(xf%r(ir) /= 0) psi_v(xf%r(ir))=1.d0!*dble((dffts%nr1*dffts%nr2*dffts%nr3))
      
        
          do iv=1,numv,2
!!product with psi_v   
             if(iv/=numv) then
                psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv,is)
                psi_r(1:dffts%nnr,2)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv+1,is)
             else
                psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv,is)
             endif
!!fourier transfrm to G   
             if(iv/=numv) then
                psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),psi_r(1:dffts%nnr,2))
             else
                psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),0.d0)
             endif
             CALL fwfft ('Wave', psic, dffts)
             if(iv/=numv) then
                psi_g(1:npw,iv)=0.5d0*(psic(dffts%nl(1:npw))+conjg(psic(dffts%nlm(1:npw))))
                psi_g(1:npw,iv+1)=(0.d0,-0.5d0)*(psic(dffts%nl(1:npw))-conjg(psic(dffts%nlm(1:npw))))
                if(gstart==2) psi_g(1,iv)=dble(psi_g(1,iv))
                if(gstart==2) psi_g(1,iv+1)=dble(psi_g(1,iv+1))
             else
                psi_g(1:npw,iv)=psic(dffts%nl(1:npw))
                if(gstart==2) psi_g(1,iv)=dble(psi_g(1,iv))
             endif
          enddo
          do iv=1,numv
             psi_g(1:npw,iv)=fac(1:npw)*psi_g(1:npw,iv)
          enddo


          do iv=1,numv,2
             
        !!fourier transform to R space
             psic(:)=(0.d0,0.d0)
             if(iv/=numv) then
                psic(dffts%nl(1:npw))  = psi_g(1:npw,iv)+(0.d0,1.d0)*psi_g(1:npw,iv+1)
                psic(dffts%nlm(1:npw)) = CONJG( psi_g(1:npw,iv) )+(0.d0,1.d0)*conjg(psi_g(1:npw,iv+1))
             else
                psic(dffts%nl(1:npw))  = psi_g(1:npw,iv)
                psic(dffts%nlm(1:npw)) = CONJG( psi_g(1:npw,iv) )
             endif
             CALL invfft ('Wave', psic, dffts)
             if(iv/=numv) then
                psi_r(:,1)= DBLE(psic(:))
                psi_r(:,2)= dimag(psic(:))
             else
                psi_r(:,1)= DBLE(psic(:))
             endif
!!product with psi_v                                                                      
             if(iv/=numv) then
                psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv,is)
                psi_r(1:dffts%nnr,2)=psi_r(1:dffts%nnr,2)*v_states(1:dffts%nnr,iv+1,is)
                xf%x(1:dffts%nnr,ir,is)=xf%x(1:dffts%nnr,ir,is)+psi_r(1:dffts%nnr,1)+psi_r(1:dffts%nnr,2)
             else
                psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv,is)
                xf%x(1:dffts%nnr,ir,is)=xf%x(1:dffts%nnr,ir,is)+psi_r(1:dffts%nnr,1)
             endif

          enddo
       enddo
    enddo
    deallocate(fac,psi_r,psi_v,psi_g)
    RETURN

  END SUBROUTINE calculate_x

SUBROUTINE calculate_hks(h0,r,nr,v_states)

    USE constants, ONLY : e2, pi, tpi, fpi
    USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
    USE fft_base,             ONLY : dfftp, dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE wannier_gw
    USE wavefunctions, ONLY : evc, psic
    USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
    USE gvect
    USE klist,    ONLY : igk_k,xk
    USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE g_psi_mod,            ONLY : h_diag, s_diag
    USE scf,       ONLY : rho, vltot, vrs, rho_core,rhog_core, scf_type
    USE lsda_mod,             ONLY : nspin
    USE uspp_init,            ONLY : init_us_2

    IMPLICIT NONE
    TYPE(hks) :: h0!Hamiltonian function to be created and initialised
    INTEGER :: r(nr)!list of r points   
    INTEGER :: nr!number of r points   
    REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,num_nbndv(1))

    INTEGER :: ir,ii,jj,ig,iv
    COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:,:),psi_g2(:,:)
    INTEGER :: kter
    LOGICAL :: lconv_root,lfirst
    REAL(kind=DP) :: anorm
    EXTERNAL :: hpsi_pw4gww2,cg_psi_pw4gww
    REAL(kind=DP), PARAMETER :: ethr=1d-10
    REAL(kind=DP) :: etxc,vtxc,ehart, charge
    REAL(kind=DP), ALLOCATABLE :: rho_fake_core(:)
    REAL(kind=DP), ALLOCATABLE :: vr(:,:)
 
    allocate(vr(dfftp%nnr,nspin))
 


    h0%nr=nr
    allocate(h0%r(nr))
    h0%r(1:nr)=r(1:nr)
    allocate(h0%h0(dffts%nnr,nr))
    allocate(h0%vxc(dffts%nnr,nspin))

 

    call allocate_bec_type ( nkb, 1, becp)
    IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
    allocate (h_diag(npw, 1),s_diag(npw,1))
    allocate(psi_g2(npw,1),psi_g(npw,1))
    do ig = 1, npw
       g2kin (ig) = ( g (1,ig)**2 + g (2,ig)**2 + g (3,ig)**2 ) * tpiba2
    enddo
    h_diag=0.d0
    do ig = 1, npw
       h_diag(ig,1)=g2kin(ig)
    enddo
    !loop on points                                                                                                                                                                                                                                                                                                          
    do ir=1,nr

       psic(1:dfftp%nnr)=0.d0
       if(h0%r(ir) /= 0) psic(h0%r(ir))=1.d0!*dble((dffts%nr1*dffts%nr2*dffts%nr3))
       !psic(1:dfftp%nnr)=v_states(1:dfftp%nnr,num_nbndv(1))
       CALL fwfft ('Wave', psic, dffts)
       psi_g(1:npw,1) = psic(dffts%nl(igk_k(1:npw,1)))
       call h_psi( npw, npw, 1, psi_g, psi_g2 )

       psic(:)=(0.d0,0.d0)
       psic(dffts%nl(1:npw))  = psi_g2(1:npw,1)
       psic(dffts%nlm(1:npw)) = CONJG( psi_g2(1:npw,1) )
       CALL invfft ('Wave', psic, dffts)
       h0%h0(1:dffts%nnr,ir)= DBLE(psic(1:dffts%nnr))


       
      
       allocate(rho_fake_core(dfftp%nnr))
       rho_fake_core(:)=0.d0
       
       CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vr )
       h0%vxc(1:dffts%nnr,1:nspin)=vr(1:dffts%nnr,1:nspin)
     
       
       deallocate(rho_fake_core)
       

    enddo

     deallocate(h_diag,s_diag,psi_g2,psi_g)
     call deallocate_bec_type(becp)
     deallocate(vr)

    return
  END SUBROUTINE calculate_hks

  subroutine pv_operator(state,ispin,ks_wfcs, l_all)
!this operator project the wavefunction state on the valence                                                               
!subspace, the valence wavefunction are in evc              

    USE io_global,            ONLY : stdout
    USE kinds,    ONLY : DP
    USE gvect
    USE wvfct,    ONLY : npwx, npw, nbnd
    USE wavefunctions, ONLY : evc, psic
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE mp_world, ONLY : world_comm
    USE wannier_gw, ONLY : num_nbndv,num_nbnds
    USE lsda_mod,   ONLy :nspin

    implicit none


    COMPLEX(kind=DP), INTENT(inout) :: state(npw)!state to be projected
    INTEGER, INTENT(in) :: ispin!spin channel 
    COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd,nspin)!KS wavefunctions       
    LOGICAL, INTENT(in) :: l_all! if true project over the entire KS manifold                                                  

    INTEGER :: iv,ig
    REAL(kind=DP), ALLOCATABLE :: prod(:)
    INTEGER :: num_proj

    if(num_nbndv(ispin)==0) return
    if(l_all) then
       num_proj=num_nbndv(ispin)
    else
       num_proj=nbnd
    endif
    allocate(prod(num_proj))
    call dgemm('T','N', num_proj,1,2*npw,2.d0,ks_wfcs(1,1,ispin),2*npwx,state,2*npw,&
         & 0.d0,prod,num_proj)
    do iv=1,num_proj
       if(gstart==2) prod(iv)=prod(iv)-dble(conjg(ks_wfcs(1,iv,ispin))*state(1))
    enddo
    call mp_sum(prod, world_comm)
    call dgemm('N','N',2*npw,1,num_proj,1.d0,ks_wfcs(1,1,ispin),2*npwx,prod,&
         &num_proj,0.d0,state,2*npw)

    deallocate(prod)
    return
  end subroutine pv_operator

  subroutine pc_operator2(state,ispin,ks_wfcs, l_all)
!this operator project the wavefunction state on the conduction                  
!subspace, the valence wavefunction are in evc                                   
!ONLY FOR GAMMA POINT NOW!!!!                                                    
    USE io_global,            ONLY : stdout
    USE kinds,    ONLY : DP
    USE gvect
    USE wvfct,    ONLY : npwx, npw, nbnd
    USE wavefunctions, ONLY : evc, psic
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE mp_world, ONLY : world_comm
    USE wannier_gw, ONLY : num_nbndv,num_nbnds
    USE lsda_mod,   ONLy :nspin

    implicit none


    COMPLEX(kind=DP), INTENT(inout) :: state(npw)!state to be projected            
    INTEGER, INTENT(in) :: ispin!spin channel
    COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd,nspin)!KS wavefunctions
    LOGICAL, INTENT(in) :: l_all! if true project over the entire KS manifold

    INTEGER :: iv,ig
    REAL(kind=DP), ALLOCATABLE :: prod(:)
    INTEGER :: num_proj
    
    if(num_nbndv(ispin)==0) return
    if(l_all) then
       num_proj=num_nbndv(ispin)
    else
       num_proj=nbnd
    endif
    allocate(prod(num_proj))
    call dgemm('T','N', num_proj,1,2*npw,2.d0,ks_wfcs(1,1,ispin),2*npwx,state,2*npw,&
         & 0.d0,prod,num_proj)
    do iv=1,num_proj
       if(gstart==2) prod(iv)=prod(iv)-dble(conjg(ks_wfcs(1,iv,ispin))*state(1))
    enddo
    call mp_sum(prod, world_comm)
    call dgemm('N','N',2*npw,1,num_proj,-1.d0,ks_wfcs(1,1,ispin),2*npwx,prod,&
         &num_proj,1.d0,state,2*npw)
 
    deallocate(prod)
    return
  end subroutine pc_operator2


  subroutine hpsi_pw4gww_krylov( ndim,psi,ppsi,et,is,numv,ks_wfcs)
! ch_psi_all (n, h, ah, e, ik, m)                                                  
    USE kinds,    ONLY : DP
    USE wvfct,    ONLY : npwx, npw, nbnd
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE mp_world, ONLY : mpime, nproc
    USE lsda_mod, ONLY : nspin


    implicit none

    INTEGER, INTENT(in) :: ndim !leading dimension of psi and psip                 
    INTEGER, INTENT(in) :: numv!number of bands                                    
    INTEGER, INTENT(in) ::is!spin index
    COMPLEX(kind=DP), INTENT(inout) :: psi(ndim,numv)
    COMPLEX(kind=DP), INTENT(out) :: ppsi(ndim,numv)
    REAL(kind=DP) ::  et(numv)
    COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd,nspin)!KS wavefunctions

    COMPLEX(kind=DP), ALLOCATABLE :: psi_save(:,:)

    INTEGER :: iv
!apply h_psi                                                                       


    do iv=1,numv
       call pc_operator2(psi(1,iv),is,ks_wfcs,.false.)
    enddo
    call h_psi( ndim, npw, numv, psi, ppsi )
    do iv=1,numv
       ppsi(1:npw,iv)=ppsi(1:npw,iv)-et(iv)*psi(1:npw,iv)
     enddo
    do iv=1,numv
       call pc_operator2(ppsi(1,iv),is,ks_wfcs,.false.)
    enddo


    return

  end subroutine hpsi_pw4gww_krylov




END MODULE convergence_gw



 SUBROUTINE read_wing ( rho, nspin, gamma_only,ipol,iw )
      !
      USE scf,              ONLY : scf_type
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, starting_ns
      USE noncollin_module, ONLY : noncolin, domag
      USE gvect,            ONLY : ig_l2g
      USE io_files,         ONLY : seqopn, prefix, tmp_dir, postfix
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_bands,         ONLY : root_bgrp, intra_bgrp_comm
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      
      USE kinds,       ONLY : DP
      USE io_files,    ONLY : create_directory
      USE io_base,     ONLY : write_rhog, read_rhog
  !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(INOUT)        :: rho
      INTEGER,          INTENT(IN)           :: nspin
      LOGICAL,          INTENT(IN)           :: gamma_only
      INTEGER,          INTENT(IN)           :: ipol!direction
      INTEGER,          INTENT(IN)           :: iw !frequency
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL :: lexist
      INTEGER :: nspin_, iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit
      CHARACTER(5) :: nfile
      CHARACTER :: npol

      write(nfile,'(5i1)') &
               & iw/10000,mod(iw,10000)/1000,mod(iw,1000)/100,mod(iw,100)/10,mod(iw,10)
      write(npol,'(1i1)') ipol

      dirname = TRIM(tmp_dir)//'/_ph0/' // TRIM(prefix) // postfix
      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         nspin_=1
      ELSE
         nspin_=nspin
      ENDIF
      ! read charge density
      CALL read_rhog(TRIM(dirname) // "wing_" // npol // "_" //nfile , &
           root_bgrp, intra_bgrp_comm, &
           ig_l2g, nspin_, rho%of_g, gamma_only )
      IF ( nspin > nspin_) rho%of_r(:,nspin_+1:nspin) = (0.0_dp, 0.0_dp)
     
      !
      RETURN
    END SUBROUTINE read_wing





