!this subroutine performes a GW calculation with EASY strategy
!for states from s_first_state to s_last_state


  SUBROUTINE easy_gw


    USE wannier_gw
    USE wvfct,                ONLY : npw
    USE io_global, ONLY : ionode, stdout
    USE mp, ONLY : mp_barrier, mp_bcast,mp_sum
    USE mp_world, ONLY : world_comm, mpime
    USE io_files, ONLY : prefix, tmp_dir, nwordwfc,iunwfc
    USE fft_base,             ONLY : dffts
    USE convergence_gw
    USE wvfct,                ONLY : nbnd, et, npwx
    USE wavefunctions, ONLY : evc
    USE constants,            ONLY : rytoev
    USE lsda_mod,             ONLY : lsda, nspin,current_spin,isk
    USE io_files,             ONLY : create_directory

    IMPLICIT NONE

    REAL(kind=DP), ALLOCATABLE :: v_states(:,:,:)!valence states in real space    
    INTEGER :: ii,iw,jj,kk
    TYPE(convergence_tests) :: ct!all data for testing convergence of basis sets     
    TYPE(self_energy) :: se!for putting all the calculated data
    COMPLEX(kind=DP), ALLOCATABLE :: freq(:)
    REAL(kind=DP), ALLOCATABLE :: e_xc(:,:),e_h(:,:),e_x(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: ks_wfcs(:,:,:)
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iun
    INTEGER :: ix,iy,iz,ir,im
    INTEGER :: passed
    INTEGER :: nr_counter
    INTEGER :: ix_start,iy_start,iz_start
    REAL(kind=DP) :: disterr,mindisterr,disterrmin
    INTEGER :: is
    INTEGER :: i_restart,i_res
    INTEGER :: n_total, i_first, i_last, n_local, i_stop
    CHARACTER(4) :: nfilex,nfiley,nfilez,nfile_orb
    CHARACTER(5) :: nfilel2
    LOGICAL :: lex, l_old_restart

    INTERFACE
       SUBROUTINE energies_xc( lda, n, m, psi, e_xc, e_h,ispin, v_states )
         USE kinds, ONLY : DP
         USE fft_base,             ONLY : dffts
         USE lsda_mod,             ONLY : nspin
          INTEGER          :: lda, n, m
          COMPLEX(kind=DP) :: psi(lda,m)
          REAL(kind=DP) :: e_xc(m), e_h(m)
          INTEGER, INTENT(in) :: ispin !spin 1,2
          REAL(kind=DP), OPTIONAL :: v_states(dffts%nnr,m, nspin)
        END SUBROUTINE energies_xc
    END INTERFACE


    call start_clock('easy_gw')

    
    allocate(v_states(dffts%nnr,num_nbnds,nspin))
    if(.not.l_truncated_coulomb) call calculate_vg0()
    if(nspin==2) then
       CALL davcio(evc,2*nwordwfc,iunwfc,2,-1)
       call evc_to_real(num_nbnds, v_states(1,1,2))
       CALL davcio(evc,2*nwordwfc,iunwfc,1,-1)
    endif
    call evc_to_real(num_nbnds, v_states(1,1,1))

    allocate(e_xc(nbnd,nspin),e_h(nbnd,nspin),e_x(nbnd,nspin))
    allocate(ks_wfcs(npw,nbnd,nspin))
    do is=1,nspin
       IF (lsda) current_spin  = isk(is)
       if(nspin/=1)  CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)!read wfcs for 
       call  energies_xc( npwx, npw, nbnd, evc, e_xc(:,is),e_h(:,is),is ,v_states)
       ks_wfcs(1:npw,1:nbnd,is)=evc(1:npw,1:nbnd)
    enddo
    CALL dft_exchange(num_nbndv,nbnd,nset,e_x,ks_wfcs)


   
    allocate(freq(n_gauss))
    do iw=1,n_gauss
       freq(iw)=(0.d0,1.d0)*(omega_gauss/dble(n_gauss)*(dble(iw)) -omega_gauss)
       if(abs(aimag(freq(iw))) < 1d-10) freq(iw)=0.d0
    enddo
 
    call initialize_memory(se)
    call set_se_energies(se, et,e_xc,e_x)
    nr_counter=1
    if(easy_grid_type==0) then
       do is=s_first_spin,s_last_spin
          do ii=s_first_state,s_last_state
             
             call start_convergence(ct,ii,is,v_states,.true.,0,0,0,n_gauss,freq,ks_wfcs)
             call calculate_convergence(ct,v_states,se,nr_counter)
             call free_memory(ct)
             
          enddo
       enddo
    elseif(easy_grid_type==1) then
!create directories
       do is=1,nspin
          do ii=s_first_state,s_last_state
             write(nfile_orb,'(4i1)') &
                  & ii/1000,mod(ii,1000)/100,mod(ii,100)/10, mod(ii,10)
             if(is==1) then
                call create_directory(trim(prefix)//'-gwl_orbital_1_'//nfile_orb)
             else
                call create_directory(trim(prefix)//'-gwl_orbital_2_'//nfile_orb)
             endif
          enddo
       enddo
!calculate total number of points
       ix_start=1+easy_grid_param(1)
       iy_start=1+easy_grid_param(2)
       iz_start=1+easy_grid_param(3)

       n_total=0
       do iz=iz_start,dffts%nr3,easy_grid_param(4)
          do iy=iy_start,dffts%nr2,easy_grid_param(4)
             do ix=ix_start,dffts%nr1,easy_grid_param(4)
                passed=0
                if(ix<easy_grid_param(5) .and. iy<easy_grid_param(5) .and. iz<easy_grid_param(5)) then
                   if( iz>dffts%nr3p_offset(mpime+1) .and. iz <=( dffts%nr3p_offset(mpime+1)+dffts%my_nr3p)) then
                      ii=(iz-dffts%nr3p_offset(mpime+1)-1)*dffts%nr2*dffts%nr1+&
                           (iy-1)*dffts%nr1+ix
                      do is=s_first_spin,s_last_spin
                         do jj=s_first_state,s_last_state
                            if(abs(v_states(ii,jj,is))> easy_psi_thrs) passed=1
                         enddo
                      enddo
                   endif
                endif
                call mp_sum(passed,world_comm)
                if(passed>0) then
                   n_total=n_total+1
                endif
             enddo
          enddo
       enddo
       write(stdout,*) 'TOTAL NUMBER OF POINTS:', n_total
       n_local=n_total/easy_split_calc_n
       if(n_local*easy_split_calc_n< n_total) n_local=n_local+1
       i_first=(easy_split_calc_i-1)*n_local+1
       i_last=i_first+n_local-1
       if(i_last>n_total) i_last=n_total

       if(restart_gww>10000) then
          l_old_restart=.false.
          restart_gww=0
       else
          l_old_restart=.true.
       endif

         if(restart_gww>1) then
            i_restart=restart_gww+i_first
            write(stdout,*) 'RESTARTING FROM POINT:', i_restart
         else
            i_restart=i_first
         endif
         i_res=0
         i_stop=i_last
         
         write(stdout,*) 'DOING POINTS RANGE', i_restart,i_stop
          ix_start=1+easy_grid_param(1)
          iy_start=1+easy_grid_param(2)
          iz_start=1+easy_grid_param(3)
          
          do iz=iz_start,dffts%nr3,easy_grid_param(4)
             do iy=iy_start,dffts%nr2,easy_grid_param(4)
                do ix=ix_start,dffts%nr1,easy_grid_param(4)
                   write(stdout,*) 'COORDINATES:',ix,iy,iz
                   passed=0
                !write(stdout,*) 'OFFSET', dffts%nr3p_offset(mpime+1),mpime
                !NOTE OFFSET and MPIME start from 0
                   if(ix<easy_grid_param(5) .and. iy<easy_grid_param(5) .and. iz<easy_grid_param(5)) then
                      if( iz>dffts%nr3p_offset(mpime+1) .and. iz <=( dffts%nr3p_offset(mpime+1)+dffts%my_nr3p)) then
                         ii=(iz-dffts%nr3p_offset(mpime+1)-1)*dffts%nr2*dffts%nr1+&
                              (iy-1)*dffts%nr1+ix
                         do is=s_first_spin,s_last_spin
                            do jj=s_first_state,s_last_state
                               if(abs(v_states(ii,jj,is))> easy_psi_thrs) passed=1
                            enddo
                         enddo
                      end if
                   endif
                   call mp_sum(passed,world_comm)
                   write(stdout,*) 'PASSED', passed
                   if(passed>0) then
                      i_res=i_res+1
                      if(l_old_restart) then
                         if(i_res>=i_restart) then
                            call start_convergence(ct,0,0,v_states,.false.,ix,iy,iz,n_gauss,freq,ks_wfcs)
                            call calculate_convergence(ct,v_states,se,nr_counter)
                            call free_memory(ct)
                         endif
                      else
                         kk=s_first_state
                         write(nfilel2,'(5i1)') &
                              & kk/10000,mod(kk,10000)/1000,mod(kk,1000)/100,mod(kk,100)/10,mod(kk,10)
                           
                         write(nfilex,'(4i1)') &
                              & ix/1000,mod(ix,1000)/100,mod(ix,100)/10, mod(ix,10)
                         write(nfiley,'(4i1)') &
                         & iy/1000,mod(iy,1000)/100,mod(iy,100)/10, mod(iy,10)
                         write(nfilez,'(4i1)') &
                         & iz/1000,mod(iz,1000)/100,mod(iz,100)/10, mod(iz,10)
!call create_directory(trim(prefix)//'-gwl_orbital_1_'//nfile_orb)
                         write(nfile_orb,'(4i1)') &
                              & kk/1000,mod(kk,1000)/100,mod(kk,100)/10, mod(kk,10)

                         inquire( file=trim(prefix)//'-gwl_orbital_1_'//nfile_orb//'/im_on_im_r_'//nfilex//'_'//nfiley&
         &//'_'//nfilez//'__'//nfilel2,exist=lex)
                       !  inquire( file=trim(prefix)//'-'//'im_on_im_r_'//nfilex//'_'//nfiley&
        ! &//'_'//nfilez//'__'//nfilel2,exist=lex)
                         if(.not.lex) then
                            call start_convergence(ct,0,0,v_states,.false.,ix,iy,iz,n_gauss,freq,ks_wfcs)
                            call calculate_convergence(ct,v_states,se,nr_counter)
                            call free_memory(ct)
                         endif
                      endif
                      if(i_res==i_stop) then
                         exit
                      end if
                   endif
            
                enddo
                if(i_res==i_stop) then
                   exit
                end if
             enddo
             if(i_res==i_stop) then
                exit
             end if
          enddo
       end if
  
    
       !if(easy_grid_type==1.or.easy_grid_type==2) call average_self_energy(se)
    !if(easy_grid_type==0 .and. l_whole_s) call solve_off_diagonal(se)
       if(ionode .and. (easy_grid_type==0 .or. i_restart==1)) then
          iun =  find_free_unit()
          open( unit=iun, file=trim(prefix)//'-'//'bands.dat', status='unknown')
          do is=1,nspin
             write(iun,*) num_nbndv(is)
             do ii=1,num_nbnds
                write(iun,*) ii,et(ii,is)*rytoev,0.d0,0.d0,(et(ii,is)-e_xc(ii,is)+e_x(ii,is))*rytoev
             enddo
          enddo
          close(iun)
       endif

    call stop_clock('easy_gw')

    call print_clock('easy_gw')
    call print_clock('gzero_complex')
    call print_clock('vpv_complex')
    call print_clock('vpv_lanczos')
    call print_clock('fft')
    call print_clock('ffts')
    call print_clock('fftw')
    

    call free_memory(se)
    deallocate(e_xc,e_h,e_x)
    deallocate(v_states,ks_wfcs)
    return
  END SUBROUTINE easy_gw
