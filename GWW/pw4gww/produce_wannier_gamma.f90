!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 SUBROUTINE produce_wannier_gamma


      
       USE wannier_gw
       USE wvfct,                ONLY : npw
       USE io_global, ONLY : ionode_id
       USE mp, ONLY : mp_barrier, mp_bcast
       USE mp_world, ONLY : world_comm
       USE io_files, ONLY : prefix, tmp_dir, nwordwfc,iunwfc
       USE wvfct,                ONLY : nbnd, et, npwx
       USE io_global,            ONLY : stdout, ionode
       USE wavefunctions_module, ONLY : evc
       USE exx,      ONLY : ecutfock,vexx,exx_div_check,exx_grid_init,exx_grid_check,exxinit,x_occupation
       USE funct,    ONLY : exx_is_active, dft_is_hybrid,start_exx,stop_exx
       USE wvfct,    ONLY : current_k, et
       USE gvecw,    ONLY : ecutwfc
       USE scf,                  ONLY : scf_type, scf_type_COPY, &
                                   create_scf_type, destroy_scf_type, &
                                   rho, rho_core, rhog_core, &
                                   v, vltot, vrs, kedtau, vnew
       USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon
       USE ldaU,                 ONLY : eth, Hubbard_U, Hubbard_lmax, &
                                   niter_with_fixed_ns, lda_plus_u
 
       USE extfield,             ONLY : tefield, etotefield
       USE lsda_mod,             ONLY : lsda, nspin,current_spin,isk

       USE gvecs,              ONLY : doublegrid
       USE fake_cond_mod
       USE constants,            ONLY : rytoev
       USE fft_base,             ONLY : dfftp
       USE exchange_custom
       USE fft_custom_gwl
       USE klist,                ONLY : nks

       implicit none

       INTEGER, EXTERNAL :: find_free_unit

       REAL(kind=DP), ALLOCATABLE :: e_xc(:,:),e_h(:,:),e_x(:,:)
       INTEGER :: ii,iw,jw,iunuterms,iun,iun2
       REAL(kind=DP), ALLOCATABLE :: tmp_rot(:,:)
       REAL(kind=DP), ALLOCATABLE :: v_states(:,:)!valence states in real space
       COMPLEX(kind=DP), ALLOCATABLE :: o_basis(:,:)!polarization basis 
                                                       !(from diagonalization of O matrix)
       REAL(kind=DP), ALLOCATABLE :: o_mat(:,:)

   !    INTEGER :: fcw_number!number of "producs of fake conduction with valence wannier" states for O matrix method
   !    COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: fcw_state! fcw states for O matrix method
   !    REAL(kind=DP), POINTER, DIMENSION(:,:) :: fcw_mat! "fcw matrix
       REAL(kind=DP), TARGET, ALLOCATABLE :: uterms(:,:), uterms_tmp(:)!matrix for 1/|r-r'| terms between product of wanniers
       REAL(kind=DP) :: charge
       INTEGER :: idumm
       REAL(kind=DP) :: rdumm1,rdumm2,rdumm3
       INTEGER :: is
       COMPLEX(kind=DP), ALLOCATABLE :: ks_wfcs(:,:,:)!Kohn-Sham wfcs (or wannier's) with spin multiplicity
       COMPLEX(kind=DP), ALLOCATABLE :: ks_wfcs_diag(:,:,:)!Kohn-Sham with spin multiplicity
       INTEGER :: num_nbndv_max
       INTEGER :: istate
       INTEGER :: numw_prod_all

       TYPE(exchange_cus) :: exx_cus

!       interface
!          subroutine fake_conduction_wannier(fcw_n,fcw_s,fcw_m,cut,s_cut)

!            USE kinds,                ONLY : DP

!            INTEGER,INTENT(out) :: fcw_n!number of "fake conduction" states for O matrix method
!            COMPLEX(kind=DP), POINTER, DIMENSION(:,:)  :: fcw_s! "fake conduction" states for O matrix method
!            REAL(kind=DP), POINTER, DIMENSION(:,:) :: fcw_m! "fake conduction" matrix
!            REAL(kind=DP), INTENT(in) :: cut!cutoff for planewaves
!            REAL(kind=DP), INTENT(in) :: s_cut!cutoff for orthonormalization
  !          INTEGER, INTENT(in) :: n_fast! number of fast conduction states, 0 = disabled
  !          REAL(kind=DP), INTENT(in) :: o_fast!offset for fast polarizability matrix 0 = disabled
  !          LOGICAL, INTENT(in) :: l_fast!if true fast polarizability matrix calculation
!          end subroutine fake_conduction_wannier
!       end interface
       allocate(e_xc(nbnd,nspin),e_h(nbnd,nspin), e_x(nbnd,nspin))
       allocate(ks_wfcs(npwx,nbnd,nspin))
       allocate(ks_wfcs_diag(npwx,nbnd,nspin))


       call start_clock('produce_wannier')

!setup global cutoff
       ecutoff_global=ecutwfc
       
       if(restart_gww>=2) then
           if(extra_pw_cutoff>0.d0)   call update_numwp(numw_prod, extra_pw_cutoff)
           if(.not.l_truncated_coulomb)  numw_prod=numw_prod+1
       endif

!setup parallel environment
#if !defined(__MPI)
         l_pmatrix=.false.
#endif
#if !defined(__SCALAPACK)
         l_pmatrix=.false.
#endif

         if(l_pmatrix) then
#if defined(__SCALAPACK)
            call blacs_pinfo(p_mpime,p_nproc)
            write(stdout,*) 'PINFO',p_mpime,p_nproc
       !     nprow=int(sqrt(real(p_nproc)))
       !     npcol=p_nproc/nprow
            write(stdout,*) 'NPROW NPCOL', nprow, npcol

            call blacs_get(0,0,icontxt)
            call blacs_gridinit(icontxt,'R',nprow, npcol)
            call blacs_gridinfo(icontxt, nprow, npcol, myrow,mycol)
            write(stdout,*) 'MYROW MYCOL', myrow,mycol
#endif
         endif
         
         if(l_scissor) then
            do is=1,nspin
               et(1:num_nbndv(is),is)=et(1:num_nbndv(is),is)+scissor(1)/rytoev
               et(num_nbndv(is)+1:num_nbnds,is)=et(num_nbndv(is)+1:num_nbnds,is)+scissor(2)/rytoev
            enddo
         endif


         if(l_truncated_coulomb) then
            l_zero=.false.
            l_wing=.false.
         else
            l_zero=.false.
            l_wing=.true.
         endif


  !set ngm max
         call max_ngm_set



         if(l_selfconsistent) then
!NOT_TO_BE_INCLUDED_START
            if(ionode) then
               iun =  find_free_unit()
               open( unit=iun, file='bands.dat', status='old',form='formatted')
               read(iun,*) n_gw_states
            endif
            call mp_bcast(n_gw_states,ionode_id,world_comm)
            allocate(ene_gw(n_gw_states,1))
            if(ionode) then
               do ii=1,n_gw_states
                  read(iun,*) idumm,rdumm1,rdumm2,ene_gw(ii,1),rdumm3
               enddo
               close(iun)
            endif
            call mp_bcast(ene_gw(:,1),ionode_id,world_comm)
            ene_gw(:,1)=ene_gw(:,1)/rytoev
            delta_self=ene_gw(n_gw_states,1)-rdumm1/rytoev!offset for all the conduction states above those calculated
!NOT_TO_BE_INCLUDED_END
         endif


         if( dft_is_hybrid()) then
!NOT_TO_BE_INCLUDED_START                                                                                                                                   
            ecutfock=exchange_fast_dual*ecutwfc
            
            CALL exx_grid_init()
            CALL exx_div_check()
            call stop_exx()
            call  exxinit
            call start_exx()
            current_k= 1
               !the following is very important                                                                                                                            
            if ( exx_is_active())  then
               CALL v_of_rho( rho, rho_core, rhog_core, &
                    ehart, etxc, vtxc, eth, etotefield, charge, v)
               CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
            end if
!NOT_TO_BE_INCLUDED_END                                                                                                                                     
         endif
            

         if(restart_gww <= 0 ) then

!calculate coulomb potential by integration over q   for PBC 
            if(.not.l_truncated_coulomb) call calculate_vg0()


!save KS wfcs if nspin==1
            if(nspin==1) CALL davcio(evc,2*nwordwfc,iunwfc,1,1)



!loop on spin

            do is=1,nspin
               IF (lsda) current_spin  = isk(is)
               if(nspin/=1)  CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)!read wfcs for 
               if(l_verbose) write(stdout,*) 'ATT1'
               FLUSH(stdout)
               CALL wfc_gamma_real(0,is)
               FLUSH(stdout)
               call  energies_xc( npwx, npw, nbnd, evc, e_xc(:,is),e_h(:,is),is )
               if( is == nspin) call write_energies_xc(e_xc)
               FLUSH( stdout )
               call go_wannier(iunwfc,1.d-9,40,num_nbndv(is), 0, is)
               call wfc_gamma_real(0,is)
!if required save MLWF for plotting
               if(l_plot_mlwf) then
                  call write_wfc_plot(0)
               endif
             !write transformation matrix u on file
               if(ionode ) call write_wannier_matrix(e_xc,e_h,is)
               do ii=1,nbnd
                  call mp_barrier( world_comm )
                  call mp_bcast(u_trans(:,ii,is),ionode_id,world_comm)
               enddo
!u_trans TO BE BROADCASTED
            enddo
            deallocate(evc)
               
          else
             if(l_scissor) then
                do is=1,nspin
                   IF (lsda) current_spin  = isk(is)
                   if(nspin/=1)  CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)!read wfcs for                                                                                                                         
                   call  energies_xc( npwx, npw, nbnd, evc, e_xc(:,is),e_h(:,is),is )
                   if( is == nspin) call write_energies_xc(e_xc)
                enddo
             endif
             deallocate(evc)
          endif
          write(stdout,*) 'USE RESTART: 1'
          FLUSH(stdout)

          if(restart_gww <= 1) then
!read coulomb potential   for PBC                   
            if(.not.l_truncated_coulomb) call read_vg0

             allocate( evc( npwx, nbnd ) )
!if required localize the valence wfcs
                if(pmat_type==2 .or. pmat_type==3 .or. pmat_type == 4) then
                   call read_wannier_matrix
                   allocate(tmp_rot(nbnd,nbnd))
                endif
                do is=1,nspin
                   call davcio(evc,2*nwordwfc,iunwfc,is,-1)
                   ks_wfcs_diag(1:npw,1:nbnd,is)=evc(1:npw,1:nbnd)
                   if(pmat_type==2 .or. pmat_type==3 .or. pmat_type == 4) then
                      tmp_rot(:,:)=dble(u_trans(:,:,is))
                      call rotate_wannier_gamma( tmp_rot,1,0)
                   endif
                   ks_wfcs(1:npw,1:nbnd,is)=evc(1:npw,1:nbnd)
                enddo
                if(pmat_type==2 .or. pmat_type==3 .or. pmat_type == 4)  deallocate(tmp_rot)

!if required calculate optimal basis for products of valence wannier fncs times 
!fake conduction states
                if (pmat_type==3 .or. pmat_type==4) then
                   if(l_verbose) write(stdout,*) 'Before fake_conduction_wannier' !ATTENZIONE
                   FLUSH(stdout)
                   !nullify(fcw_state)
                   !call fake_conduction_wannier(fcw_number,fcw_state,fcw_mat,pmat_cutoff,s_pmat)!,n_fast_pmat,off_fast_pmat,l_fast_pola)
                   call start_clock('f_conduction')
                   if(.not.lwannier) then
                      if(.not.l_real) then
                         call fake_conduction_wannier(pmat_cutoff,s_pmat,ks_wfcs,l_frac_occ,ks_wfcs_diag,l_cond_pol_base)
                      else
!NOT_TO_BE_INCLUDED_START
                         call fake_conduction_real(pmat_cutoff,s_pmat,ks_wfcs,l_frac_occ,ks_wfcs_diag,l_cond_pol_base)
!NOT_TO_BE_INCLUDED_END
                      endif
                   else
!NOT_TO_BE_INCLUDED_START
!still to be implemented with spin
                      call fake_conduction_wannier_real(pmat_cutoff,s_pmat)
!NOT_TO_BE_INCLUDED_END
                   endif
                   call stop_clock('f_conduction')
                   if(l_verbose) write(stdout,*) 'After fake_conduction_wannier' !ATTENZIONE
                   FLUSH(stdout)
                   call print_clock('f_conduction')
                   call print_clock('mpsum')
                   call print_clock('fc_optimal')
                   call print_clock('fc_merge')
                   call print_clock('fc_loop')
                   call print_clock('fc_dgemm')
                   deallocate(fcw_state,fcw_mat) 
                endif

!calculates the polarization basis diagonalizing the O matrix
                if(pmat_type==0 .or. pmat_type == 1 .or. pmat_type == 2) then
                   allocate(v_states(dfftp%nnr,num_nbndv(1)))
                else
                   allocate(v_states(1,1))
                endif
                numw_prod_all=numw_prod

                if(extra_pw_cutoff>0.d0)   call update_numwp(numw_prod_all, extra_pw_cutoff)
                allocate(o_basis(npw,numw_prod_all))
                write(stdout,*) 'NUMW_PROD_ALL', numw_prod_all!DEBUG
!calculate array of valence states in real space
                if(l_verbose) write(stdout,*) 'Call evc_to_real'
                FLUSH(stdout)
                if(pmat_type==0 .or. pmat_type == 1 .or. pmat_type == 2) call evc_to_real(num_nbndv(1), v_states)

!allocate and set trial product with random number
                if(l_verbose) write(stdout,*) 'Call o_basis_init'
                FLUSH(stdout)
                !nullify(fcw_state)
                if(pmat_type==0 .or.pmat_type == 1 .or. pmat_type == 2) call o_basis_init(numw_prod,&
             &   o_basis,num_nbndv(1),v_states,pmat_cutoff,pmat_type,fcw_number,fcw_state,fcw_mat,pmat_ethr)
!diagonalize products
                if(l_verbose) write(stdout,*) 'Call o_bands'
                FLUSH(stdout)
                call start_clock('o_bands')
!note: v_states is relevant only for pmat_type ==  0,1,2
                call o_bands(num_nbndv(1), v_states,numw_prod,o_basis,pmat_ethr,pmat_cutoff,pmat_type)
                call stop_clock('o_bands')
!write them to disk
                if(l_v_basis) then
                   call v_basis(numw_prod,o_basis,v_cutoff)
                endif
                if(l_verbose) write(stdout,*) 'Call o_bands write'
                FLUSH(stdout)
!if PBC add also the last one the G=0 element
!NOTE for PBC that numpw BECOMES numpw+1 AT THIS POINT, FOLLOWING ROUTINE
!if required add plane waves to the basis obtained till now
                if(extra_pw_cutoff>0.d0) then
                   call o_extra_pw( o_basis, numw_prod, numw_prod_all,extra_pw_cutoff)
                endif

                call o_basis_write(numw_prod, o_basis,.true.,ecutoff_global,.not.l_truncated_coulomb)
!deallocate arrays
                deallocate(v_states,o_basis)
                deallocate(evc)
                if(.not.l_zero) then
                   call wannier_uterms(nset,.false.,.false.,2, ecutoff_global)
                    FLUSH( stdout )
                    allocate(uterms(numw_prod,numw_prod))
                    if(ionode) then
                       iunuterms =  find_free_unit()
                       open( unit= iunuterms, file=trim(tmp_dir)//trim(prefix)//'.uterms', status='old',form='unformatted')
                    endif
                    allocate(uterms_tmp(numw_prod))
                    do iw=1,numw_prod
                       if(ionode) read(iunuterms) uterms_tmp(1:iw)
                       call mp_bcast(uterms_tmp(1:iw), ionode_id,world_comm)
                       do jw=1,iw
                          uterms(iw,jw)=uterms_tmp(jw)
                          uterms(jw,iw)=uterms(iw,jw)
                       enddo
                    enddo
                    deallocate(uterms_tmp)
                    if(ionode)    close(iunuterms)
                    if(ionode) call write_vpot_matrix(uterms,0)
                    deallocate(uterms)
                    FLUSH( stdout )
                else
                   write(stdout,*) 'NOT LZERO NOT IMPLEMENTED'
                   FLUSH(stdout)
                   stop
                endif
                if(l_verbose) write(stdout,*) 'OUT OF RESTART_GWW1',numw_prod
                FLUSH(stdout)
                call mp_barrier( world_comm )

             endif
             
          if(restart_gww <= 2 ) then
!read coulomb potential   for PBC                                               
              write(stdout,*) 'USE RESTART: 2 LANCZOS RESTART:0'
              FLUSH(stdout)

            if(.not.l_truncated_coulomb) call read_vg0

            if(l_big_system.and.l_list) then
!read list of KS state for which the self-energy will be calculated
               if(ionode) then
                  iun =  find_free_unit()
                  open( unit=iun, file='list_1.dat', status='old')
                  read(iun,*) n_list(1)
                  if(nspin==2) then
                     iun2 =  find_free_unit()
                     open( unit=iun2, file='list_2.dat', status='old')
                     read(iun,*) n_list(2)
                  else
                     n_list(2)=0
                  endif
               endif
               call mp_bcast(n_list,ionode_id,world_comm)
               allocate(i_list(max(n_list(1),n_list(2)),2))
               i_list=0
               if(ionode) then
                  do ii=1,n_list(1)
                     read(iun,*) i_list(ii,1)
                  enddo
                  close(iun)
                  if(nspin==2) then
                     do ii=1,n_list(2)
                        read(iun2,*) i_list(ii,2)
                     enddo
                     close(iun2)
                  endif
               endif
               call mp_bcast(i_list,ionode_id,world_comm)
               s_first_state=1
               s_last_state=num_nbnds
            endif

             do is=1,nspin
                IF (lsda) current_spin  = isk(is)
                allocate( evc( npwx, nbnd ) )
                call davcio(evc,2*nwordwfc,iunwfc,is,-1)
!if required calculate partial occupancies factors
                if(l_frac_occ) then
!NOT_TO_BE_INCLUDED_START
                   call pola_partial(numw_prod,is)
!NOT_TO_BE_INCLUDED_END
                endif

!if EXX is one calculates stuff for Fock operator
                if(dft_is_hybrid()) then
!NOT_TO_BE_INCLUDED_START
                   call  exxinit
                   current_k= 1
!NOT_TO_BE_INCLUDED_END
                endif
             


                call read_wannier_matrix
                allocate(tmp_rot(nbnd,nbnd))
                tmp_rot(1:nbnd,1:nbnd)=dble(u_trans(1:nbnd,1:nbnd,is))
                if(l_t_wannier) call rotate_wannier_gamma( tmp_rot,1,0)
                deallocate(tmp_rot)
             
                if(n_pola_lanczos > numw_prod) n_pola_lanczos=numw_prod
                if(n_self_lanczos > numw_prod) n_self_lanczos=numw_prod
             
                if(n_pola_lanczos_eff == 0) n_pola_lanczos_eff=n_pola_lanczos
                if(n_self_lanczos_eff == 0) n_self_lanczos_eff=n_self_lanczos
             
                if(lanczos_restart <= 0) then
                   call start_clock('pola_basis')
                   if(.not.l_real) then
                      call pola_basis_lanczos(nset,n_pola_lanczos,numw_prod,nsteps_lanczos_pola,is)
                   else
!NOT_TO_BE_INCLUDED_START
                      call pola_basis_lanczos_real(nset,n_pola_lanczos,numw_prod,nsteps_lanczos_pola,is)
!NOT_TO_BE_INCLUDED_END
                   endif
                   call stop_clock('pola_basis')
                endif
                write(stdout,*) 'USE RESTART: 2 LANCZOS_RESTART:1'
                FLUSH(stdout)
              
           
                if(lanczos_restart <= 1) then
                   call start_clock('global_pola')
                   call davcio(evc,2*nwordwfc,iunwfc,is,-1)!re-read for testing and for selfconsistency
                   call global_pola_lanczos(n_pola_lanczos,n_pola_lanczos_eff,s_pola_lanczos,nump_lanczos,&
                        nsteps_lanczos_pola,numw_prod,is,l_ts_eigen)
                   call stop_clock('global_pola')
                endif
                write(stdout,*) 'USE RESTART: 2 LANCZOS_RESTART:2'
                FLUSH(stdout)

                if(lanczos_restart <= 2) then
                   call start_clock('self_basis')
                   call davcio(evc,2*nwordwfc,iunwfc,is,-1)
                   if(.not.l_real) then 
                      call self_basis_lanczos(nset,n_self_lanczos,numw_prod,nsteps_lanczos_self,is,l_full,n_full(is))
                   else
!NOT_TO_BE_INCLUDED_START
                      call self_basis_lanczos_real(nset,n_self_lanczos,numw_prod,nsteps_lanczos_self,is)
!NOT_TO_BE_INCLUDED_END
                   endif
                   call stop_clock('self_basis')
                   call print_clock('self_basis')
                   call print_clock('sl_loop')
                   call print_clock('sl_dgemm')
                   call print_clock('sl_dsyevX')
                   call print_clock('sl_mpbcast')
                   call print_clock('sl_merge')
                endif
                FLUSH( stdout )
           !  ALLOCATE( evc( npwx, nbnd ) )
              !  deallocate(evc)
                write(stdout,*) 'USE RESTART: 2 LANCZOS_RESTART:3'
                FLUSH(stdout)

                if(lanczos_restart <= 3) then
            !    CALL davcio(evc,2*nwordwfc,iunwfc,1,-1)
            !    CALL dft_exchange(num_nbndv,num_nbnds,nset,e_x)
            !    FLUSH( stdout )
                   call mp_barrier( world_comm )
                   call davcio(evc,2*nwordwfc,iunwfc,is,-1)
                   call start_clock('global_self')
                   if(.not.l_big_system) then
                      call global_self_lanczos(n_self_lanczos,n_self_lanczos_eff,s_self_lanczos,nums_lanczos,&
                           nsteps_lanczos_self,numw_prod,s_g_lanczos,is,l_ts_eigen,1,l_full)
                   else
                      if(.not.l_list) then
                         do istate=s_first_state,s_last_state
                            call global_self_lanczos(n_self_lanczos,n_self_lanczos_eff,s_self_lanczos,nums_lanczos,&
                                 nsteps_lanczos_self,numw_prod,s_g_lanczos,is,l_ts_eigen,istate,l_full)
                         enddo
                      else
                         do ii=1,n_list(is)
                            istate=i_list(ii,is)
                            call global_self_lanczos(n_self_lanczos,n_self_lanczos_eff,s_self_lanczos,nums_lanczos,&
                                 nsteps_lanczos_self,numw_prod,s_g_lanczos,is,l_ts_eigen,istate,l_full)
                         enddo
                      endif
                   endif
                   call stop_clock('global_self')
                endif
             
                deallocate(evc)
           !  if(.not.l_truncated_coulomb) call calculate_wing(nset,2)
             enddo
             if(l_big_system.and.l_list) deallocate(i_list)
          endif
     
           
          write(stdout,*) 'USE RESTART: 3 LANCZOS_RESTART /=2,3'
          FLUSH(stdout)

          if(restart_gww <= 3 .and. lanczos_restart /=3 .and. lanczos_restart /=2 ) then !ATTENZIONE RESTART lanczos_restart never been here!
!read coulomb potential   for PBC                                                          
             if(.not.l_truncated_coulomb) call read_vg0

             if(.not.l_truncated_coulomb) call calculate_wing(nset,2)
             FLUSH( stdout )
             ALLOCATE( evc( npwx, nbnd ) )
             do is=1,nspin
                CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)
                ks_wfcs(1:npw,1:nbnd,is)=evc(1:npw,1:nbnd)
             enddo
             CALL dft_exchange(num_nbndv,num_nbnds,nset,e_x,ks_wfcs)
             FLUSH( stdout )
!DEBUG TEST
             if(l_verbose) then
                if(nspin==1) then
                   num_nbndv_max=num_nbndv(1)
                else
                   num_nbndv_max=max(num_nbndv(1),num_nbndv(2))
                endif
                CALL  setup_exx_cus(nspin,num_nbndv_max,num_nbndv,evc, exx_cus, 1.d0, 40.d0, truncation_radius)
                CALL dft_exchange_fast(1,num_nbnds,ks_wfcs(:,:,1),exx_cus)
                write(stdout,*) 'BEFORE periodic_dft_exchange'
                FLUSH(stdout)
                !CALL  periodic_dft_exchange(exx_cus)
                write(stdout,*) 'AFTER periodic_dft_exchange'
                FLUSH(stdout)

                call free_memory_exx_cus(exx_cus)
             endif
           deallocate(evc)
          endif
          write(stdout,*) 'USE RESTART: 4 LANCZOS_RESTART /=2,3'
          FLUSH(stdout)

          if(restart_gww <= 4 .and. l_semicore .and. lanczos_restart /=3 .and. lanczos_restart /=2) then
!NOT_TO_BE_INCLUDED_START
             allocate( evc( npwx, nbnd ) )
             do is=1,nspin
                CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)
                call semicore(n_semicore, num_nbnds,is)
             enddo
             deallocate(evc)
!NOT_TO_BE_INCLUDED_END
          endif
          write(stdout,*) 'USE RESTART: 5 LANCZOS_RESTART /=2,3'
          FLUSH(stdout)

          if(restart_gww <= 5 .and. l_semicore_read .and. lanczos_restart /=3 .and. lanczos_restart /=2) then
!NOT_TO_BE_INCLUDED_START
             if(.not.l_truncated_coulomb) call read_vg0
             allocate( evc( npwx, nbnd ) )
             do is=1,nspin
                CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)
                call semicore_read(num_nbnds,numw_prod,is)
             enddo
             deallocate(evc)
!NOT_TO_BE_INCLUDED_END
          endif

!NOT_TO_BE_INCLUDED_START

          if(restart_gww<=6  .and. l_full) then
             call write_pola_basis(numw_prod,0)
          endif

          if(restart_gww<=6  .and. l_simple) then
             call write_pola_basis(numw_prod,1)
          endif


          if (l_bse) then

          call read_wannier_matrix
          allocate(tmp_rot(nbnd,nbnd))
          allocate( evc( npwx, nbnd ) )

          do is=1,nspin
             allocate(o_mat(num_nbndv(is),num_nbndv(is)))
             call davcio(evc,2*nwordwfc,iunwfc,is,-1)
             
             tmp_rot(:,:)=dble(u_trans(:,:,is))
             call rotate_wannier_gamma( tmp_rot,1,0)
             
             if(.not.l_truncated_coulomb) call read_vg0
             call wannier_bse(is,evc,o_mat)
             
             deallocate(o_mat)
          enddo

          deallocate(tmp_rot)
          deallocate(evc)

          endif
!NOT_TO_BE_INCLUDED_END

          deallocate(e_xc,e_h,e_x)
          deallocate(ks_wfcs,ks_wfcs_diag)
         
          write(stdout,*) 'PW4GWW COMPLETED'
          call stop_clock('produce_wannier')
          call print_clock('produce_wannier')
          call print_clock('f_conduction')
          call print_clock('o_bands')
          call print_clock('pola_basis')
          call print_clock('global_pola')
          call print_clock('self_basis')
          call print_clock('cft3t')
          call print_clock('h_psi')
          call print_clock('fft')
          call print_clock('ffts')
          call print_clock('fftw')
          call print_clock('davcio')
          call print_clock('mpsum')
          call print_clock('global_self')
          call print_clock('lanczos_state')
          FLUSH(stdout)
          return
  END SUBROUTINE produce_wannier_gamma
