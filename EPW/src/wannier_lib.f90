!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Wannier90 v2.0 authors:                                    !
!           Arash A. Mostofi   (Imperial College London)     !
!           Jonathan R. Yates  (University of Oxford)        !
!           Giovanni Pizzi     (EPFL, Switzerland)           !
!           Ivo Souza          (Universidad del Pais Vasco)  !
!                                                            !
! Contributors:                                              !
!          Young-Su Lee        (KIST, S. Korea)              !
!          Matthew Shelley     (Imperial College London)     !
!          Nicolas Poilvert    (Penn State University)       !
!          Raffaello Bianco    (Paris 6 and CNRS)            !
!          Gabriele Sclauzero  (ETH Zurich)                  !
!                                                            !
!  Please cite                                               !
!                                                            !
!  [ref] A. A. Mostofi, J. R. Yates, Y.-S. Lee, I. Souza,    !
!        D. Vanderbilt and N. Marzari, "Wannier90: A Tool    !
!        for Obtaining Maximally Localised Wannier           !
!        Functions", Computer Physics Communications,        !
!        178, 685 (2008)                                     !
!                                                            !
!  in any publications arising from the use of this code.    !
!                                                            !
!  Wannier90 is based on Wannier77, written by N. Marzari,   !
!  I. Souza and D. Vanderbilt. For the method please cite    !
!                                                            !
!  [ref] N. Marzari and D. Vanderbilt,                       !
!        Phys. Rev. B 56 12847 (1997)                        !
!                                                            !
!  [ref] I. Souza, N. Marzari and D. Vanderbilt,             !
!        Phys. Rev. B 65 035109 (2001)                       !
!                                                            !
!                                                            !
! Copyright (C) 2007-13 Jonathan Yates, Arash Mostofi,       !
!                Giovanni Pizzi, Young-Su Lee,               !
!                Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
     num_atoms_loc,atom_symbols_loc,atoms_cart_loc, gamma_only_loc,spinors_loc, &
     nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
     proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)

  use w90_constants
  use w90_parameters
  use w90_io
  use w90_kmesh
 
  implicit none

  character(len=*), intent(in) :: seed__name
  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
  real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
  real(kind=dp), dimension(3,num_kpts_loc), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_tot
  integer, intent(in) :: num_atoms_loc
  character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(kind=dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
  logical, intent(in) :: gamma_only_loc
  logical, intent(in) :: spinors_loc
  integer, intent(out) :: nntot_loc
  integer, dimension(num_kpts_loc,num_nnmax), intent(out) :: nnlist_loc
  integer,dimension(3,num_kpts_loc,num_nnmax), intent(out) :: nncell_loc
  integer, intent(out) :: num_bands_loc
  integer, intent(out) :: num_wann_loc
  real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_site_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_radial_loc
  real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_z_loc
  real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_x_loc
  real(kind=dp), dimension(num_bands_tot), intent(out) :: proj_zona_loc
  integer, dimension(num_bands_tot), intent(out) :: exclude_bands_loc
  integer, dimension(num_bands_tot), optional, intent(out) :: proj_s_loc  
  real(kind=dp), dimension(3,num_bands_tot), optional, intent(out) :: proj_s_qaxis_loc


  real(kind=dp) time0,time1,time2
  character(len=9) :: stat,pos,cdate,ctime
  integer :: ierr
  logical :: wout_found

  time0=io_time()

  library=.true.
!  seedname="wannier"
  seedname=trim(adjustl(seed__name))
  inquire(file=trim(seedname)//'.wout',exist=wout_found)
  if (wout_found) then
     stat='old'
  else
     stat='replace'
  endif
  pos='append'

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.wout',status=trim(stat),position=trim(pos))

  call param_write_header()

  write(stdout,'(/a/)') ' Wannier90 is running in LIBRARY MODE'
  write(stdout,'(a/)') ' Setting up k-point neighbours...'

  ! copy local data into module variables
  mp_grid=mp_grid_loc
  num_kpts=num_kpts_loc
  real_lattice=real_lattice_loc
  recip_lattice=recip_lattice_loc
  allocate ( kpt_latt(3,num_kpts) ,stat=ierr)
  if (ierr/=0) call io_error('Error allocating kpt_latt in wannier_setup')
  kpt_latt=kpt_latt_loc
  num_atoms=num_atoms_loc
  call param_lib_set_atoms(atom_symbols_loc,atoms_cart_loc)
  gamma_only=gamma_only_loc
  spinors=spinors_loc

  ! set num_bands and cell_volume as they are written to output in param_write
  num_bands = num_bands_tot - num_exclude_bands
  call param_read()

  cell_volume = real_lattice(1,1)*(real_lattice(2,2)*real_lattice(3,3)-real_lattice(3,2)*real_lattice(2,3)) +&
                real_lattice(1,2)*(real_lattice(2,3)*real_lattice(3,1)-real_lattice(3,3)*real_lattice(2,1)) +& 
                real_lattice(1,3)*(real_lattice(2,1)*real_lattice(3,2)-real_lattice(3,1)*real_lattice(2,2))
  call param_write()

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to read parameters  ',time1-time0,' (sec)'

  call kmesh_get()


  ! Now we zero all of the local output data, then copy in the data
  ! from the parameters module

  nntot_loc         = 0
  nnlist_loc        = 0   
  nncell_loc        = 0     
  proj_site_loc     = 0.0_dp   
  proj_l_loc        = 0
  proj_m_loc        = 0
  proj_z_loc        = 0.0_dp
  proj_x_loc        = 0.0_dp
  proj_radial_loc   = 0
  proj_zona_loc     = 0.0_dp
  exclude_bands_loc = 0

  nntot_loc       = nntot        
  nnlist_loc(:,1:nntot)   =  nnlist(:,1:nntot)       
  nncell_loc(:,:,1:nntot) =  nncell(:,:,1:nntot)       
  num_bands_loc=num_bands_tot-num_exclude_bands
  num_wann_loc=num_wann
  if(allocated(proj_site)) then
     proj_site_loc(:,1:num_proj)   = proj_site(:,1:num_proj)    
     proj_l_loc(1:num_proj)        = proj_l(1:num_proj)          
     proj_m_loc(1:num_proj)        = proj_m(1:num_proj)           
     proj_z_loc(:,1:num_proj)      = proj_z(:,1:num_proj)     
     proj_x_loc(:,1:num_proj)      = proj_x(:,1:num_proj)       
     proj_radial_loc(1:num_proj)   = proj_radial(1:num_proj)            
     proj_zona_loc(1:num_proj)     = proj_zona(1:num_proj) 
     if(allocated(proj_s) .and. present(proj_s_loc) .and. present(proj_s_qaxis_loc)) then
             proj_s_loc(1:num_proj)     = proj_s(1:num_proj) 
             proj_s_qaxis_loc(:,1:num_proj)   = proj_s_qaxis(:,1:num_proj)    
          end if
       endif
  if(allocated(exclude_bands)) then
     exclude_bands_loc(1:num_exclude_bands) = exclude_bands(1:num_exclude_bands)
  end if

! SP: could not made postproc_setup = .true. from EPW
!  if (postproc_setup) then
   call kmesh_write()
   write(stdout,'(1x,a25,f11.3,a)') 'Time to write kmesh      ',io_time(),' (sec)'
   write(stdout,'(/a)') ' '//trim(seedname)//'.nnkp written.'
!  endif


  call kmesh_dealloc()
  call param_dealloc()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to write kmesh      ',io_time(),' (sec)'

  write(stdout,'(/a/)') ' Finished setting up k-point neighbours.'

  call io_date(cdate,ctime)

  write(stdout,'(2a)') ' Exiting wannier_setup in wannier90 ',ctime

  close(stdout)

  
end subroutine wannier_setup


subroutine wannier_run(seed__name,mp_grid_loc,num_kpts_loc, &
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_loc, &
     num_wann_loc,nntot_loc,num_atoms_loc,atom_symbols_loc, &
     atoms_cart_loc,gamma_only_loc,m_matrix_loc,A_matrix_loc,eigenvalues_loc, &
     u_matrix_loc,u_matrix_opt_loc,lwindow_loc,wann_centres_loc, &
     wann_spreads_loc,spread_loc)


  use w90_constants
  use w90_parameters
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport

  implicit none

  character(len=*), intent(in) :: seed__name
  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
  real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
  real(kind=dp), dimension(3,num_kpts_loc), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_loc
  integer, intent(in) :: num_wann_loc
  integer, intent(in) :: nntot_loc
  integer, intent(in) :: num_atoms_loc
  character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(kind=dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
  logical, intent(in) :: gamma_only_loc
  complex(kind=dp), dimension(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc), intent(in) :: m_matrix_loc
  complex(kind=dp), dimension(num_bands_loc,num_wann_loc,num_kpts_loc), intent(in) :: A_matrix_loc
  real(kind=dp), dimension(num_bands_loc,num_kpts_loc), intent(in) :: eigenvalues_loc
  complex(kind=dp), dimension(num_wann_loc,num_wann_loc,num_kpts_loc), intent(out) :: u_matrix_loc
  complex(kind=dp), dimension(num_bands_loc,num_wann_loc,num_kpts_loc), optional, intent(out) :: u_matrix_opt_loc
  logical, dimension(num_bands_loc,num_kpts_loc), optional, intent(out) :: lwindow_loc
  real(kind=dp), dimension(3,num_wann_loc), optional, intent(out) :: wann_centres_loc
  real(kind=dp), dimension(num_wann_loc), optional, intent(out) :: wann_spreads_loc
  real(kind=dp), dimension(3), optional, intent(out) :: spread_loc

  real(kind=dp) time0,time1,time2
  character(len=9) :: stat,pos,cdate,ctime
  integer :: ierr,loop_k,loop_w
  logical :: wout_found

  integer :: nkp,nn,n,m

  time0=io_time()

  library=.true.
!  seedname="wannier"
  seedname=trim(adjustl(seed__name))
  inquire(file=trim(seedname)//'.wout',exist=wout_found)
  if (wout_found) then
     stat='old'
  else
     stat='replace'
  endif
  pos='append'

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.wout',status=trim(stat),position=trim(pos))

  call io_date(cdate,ctime)

  write(stdout,'(/,2a,/)') ' Resuming Wannier90 at ',ctime

!  call param_write_header

  ! copy local data into module variables
  num_bands=num_bands_loc
  mp_grid=mp_grid_loc
  num_kpts=num_kpts_loc
  real_lattice=real_lattice_loc
  recip_lattice=recip_lattice_loc
  allocate ( kpt_latt(3,num_kpts) ,stat=ierr)
  if (ierr/=0) call io_error('Error allocating kpt_latt in wannier_setup')
  kpt_latt=kpt_latt_loc
  allocate(eigval(num_bands,num_kpts),stat=ierr)
  if (ierr/=0) call io_error('Error allocating eigval in wannier_setup')
  eigval=eigenvalues_loc
  num_atoms=num_atoms_loc
  gamma_only=gamma_only_loc

  call param_lib_set_atoms(atom_symbols_loc,atoms_cart_loc)

  call param_read()

  call param_write()

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to read parameters  ',time1-time0,' (sec)'

  call kmesh_get()

  time2=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to get kmesh        ',time2-time1,' (sec)'

  allocate ( u_matrix( num_wann,num_wann,num_kpts),stat=ierr)
  if (ierr/=0) call io_error('Error in allocating u_matrix in overlap_read')
  
  if (disentanglement) then
     allocate(m_matrix_orig(num_bands,num_bands,nntot,num_kpts),stat=ierr)
     if (ierr/=0) call io_error('Error in allocating m_matrix_orig in overlap_read')
     allocate(a_matrix(num_bands,num_wann,num_kpts),stat=ierr)
     if (ierr/=0) call io_error('Error in allocating a_matrix in overlap_read')
     allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
     if (ierr/=0) call io_error('Error in allocating u_matrix_opt in overlap_read')
  else
     allocate ( m_matrix( num_wann,num_wann,nntot,num_kpts),stat=ierr)
     if (ierr/=0) call io_error('Error in allocating m_matrix in overlap_read')
  endif
    
  if (disentanglement) then
     m_matrix_orig = m_matrix_loc
     a_matrix      = a_matrix_loc
     u_matrix_opt  = cmplx_0
     u_matrix      = cmplx_0
  else
     m_matrix=m_matrix_loc
     u_matrix=a_matrix_loc
  endif

!!$  ! Check Mmn(k,b) is symmetric in m and n for gamma_only case
!!$  if (gamma_only) call overlap_check_m_symmetry()

  if(disentanglement) then
     have_disentangled = .false.
     call dis_main()
     have_disentangled=.true.
     call param_write_chkpt('postdis')
     time1=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time to disentangle      ',time1-time2,' (sec)'     
  else
     if (gamma_only) then
        call overlap_project_gamma()
     else
        call overlap_project()
     endif
     time1=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time to project overlaps ',time1-time2,' (sec)'     
  end if

  if (gamma_only) then
     call wann_main_gamma()
  else
     call wann_main()
  endif

  call param_write_chkpt('postwann')

  time2=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time for wannierise      ',time2-time1,' (sec)'     

  if (wannier_plot .or. bands_plot .or. fermi_surface_plot .or. hr_plot) then
     call plot_main()
     time1=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time for plotting        ',time1-time2,' (sec)'     
  end if

  time2=io_time()
  if (transport) then
     call tran_main()
     time1=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time for transport       ',time1-time2,' (sec)'
  end if

  ! Now we zero all of the local output data, then copy in the data
  ! from the parameters module

  u_matrix_loc=u_matrix
  if(present(u_matrix_opt_loc) .and. present(lwindow_loc)) then
  if(disentanglement) then
     u_matrix_opt_loc=u_matrix_opt
     lwindow_loc=lwindow
  else
     u_matrix_opt_loc=cmplx_0
     do loop_k=1,num_kpts
        do loop_w=1,num_wann
           u_matrix_opt_loc(loop_w,loop_w,loop_k)=cmplx_1
        end do
     end do
     lwindow_loc=.true.
  end if
  end if  


  if(present(wann_centres_loc)) wann_centres_loc=wannier_centres
  if(present(wann_spreads_loc)) wann_spreads_loc=wannier_spreads
  if(present(spread_loc)) then
     spread_loc(1)=omega_total
     spread_loc(2)=omega_invariant
     spread_loc(3)=omega_tilde
  endif
  call hamiltonian_dealloc()
  call overlap_dealloc()
  call kmesh_dealloc()
  call param_dealloc()

  write(stdout,'(1x,a25,f11.3,a)') 'Total Execution Time     ',io_time()-time0,' (sec)'

  if (timing_level>0) call io_print_timings()

  write(stdout,*) 
  write(stdout,'(1x,a)') 'All done: wannier90 exiting'
  close(stdout)



end subroutine wannier_run
