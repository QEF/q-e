!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! pw2wannier was written by Stefano de Gironcoli
! with later additions by
! Jonathan Yates - spinors
! Arash Mostofi - gamma point and transport things
! Timo Thonhauser, Graham Lopez, Ivo Souza
!         uHu, uIu terms for orbital magnetisation
! please send bugs and comments to
! Jonathan Yates and Arash Mostofi
! Takashi Koretsune and Florian Thoele -- noncollinear and USPPs
! Valerio Vitale - Selected columns of density matrix (SCDM)
! Jae-Mo Lihm - SCDM with noncollinear
!
!
! Riccardo De Gennaro: moving the wannier module out of pw2wannier90
!
!
!------------------------------------------------------------------------
module wannier
   USE kinds, only : DP
   !integer, allocatable :: nnb(:)       ! #b  (ik)
   integer              :: nnb          ! #b
   integer, allocatable :: kpb(:,:)     ! k+b (ik,ib)
   integer, allocatable :: g_kpb(:,:,:) ! G_k+b (ipol,ik,ib)
   integer, allocatable :: ig_(:,:)     ! G_k+b (ipol,ik,ib)
   integer, allocatable :: lw(:,:), mw(:,:) ! l and m of wannier (16,n_wannier)
   integer, allocatable :: num_sph(:)   ! num. func. in lin. comb., (n_wannier)
   logical, allocatable :: excluded_band(:)
   ! begin change Lopez, Thonhauser, Souza
   integer  :: iun_nnkp,iun_mmn,iun_amn,iun_band,iun_spn,iun_plot,iun_parity,&
        nnbx,nexband,iun_uhu,&
        iun_uIu,& !ivo
   ! end change Lopez, Thonhauser, Souza
        iun_sHu, iun_sIu ! shc
   integer  :: n_wannier !number of WF
   integer  :: n_proj    !number of projection
   complex(DP), allocatable :: gf(:,:)  ! guding_function(npwx,n_wannier)
   complex(DP), allocatable :: gf_spinor(:,:)
   complex(DP), allocatable :: sgf_spinor(:,:)
   integer               :: ispinw, ikstart, ikstop, iknum
   character(LEN=15)     :: wan_mode    ! running mode
   logical               :: logwann, wvfn_formatted, write_unk, write_eig, &
   ! begin change Lopez, Thonhauser, Souza
                            write_amn,write_mmn,reduce_unk,write_spn,&
                            write_unkg,write_uhu,&
                            write_dmn,read_sym, & !YN
                            write_uIu, spn_formatted, uHu_formatted, uIu_formatted, & !ivo
   ! end change Lopez, Thonhauser, Souza
   ! shc
                            write_sHu, write_sIu, sHu_formatted, sIu_formatted, &
   ! end shc
   ! vv: Begin SCDM keywords
                            scdm_proj
   character(LEN=15)     :: scdm_entanglement
   real(DP)              :: scdm_mu, scdm_sigma
   ! vv: End SCDM keywords
   ! run check for regular mesh
   logical               :: regular_mesh = .true.
   ! input data from nnkp file
   real(DP), allocatable :: center_w(:,:)     ! center_w(3,n_wannier)
   integer,  allocatable :: spin_eig(:)
   real(DP), allocatable :: spin_qaxis(:,:)
   integer, allocatable  :: l_w(:), mr_w(:) ! l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
   integer, allocatable  :: r_w(:)      ! index of radial function (n_wannier) as from table 3.3 of spec.
   real(DP), allocatable :: xaxis(:,:),zaxis(:,:) ! xaxis and zaxis(3,n_wannier)
   real(DP), allocatable :: alpha_w(:)  ! alpha_w(n_wannier) ( called zona in wannier spec)
   !
   real(DP), allocatable :: csph(:,:)    ! expansion coefficients of gf on QE ylm function (16,n_wannier)
   CHARACTER(len=256) :: seedname  = 'wannier'  ! prepended to file names in wannier90
   ! For implementation of wannier_lib
   integer               :: mp_grid(3)            ! dimensions of MP k-point grid
   real(DP)              :: rlatt(3,3),glatt(3,3) ! real and recip lattices (Cartesian co-ords, units of Angstrom)
   real(DP), allocatable :: kpt_latt(:,:)  ! k-points in crystal co-ords. kpt_latt(3,iknum)
   real(DP), allocatable :: atcart(:,:)    ! atom centres in Cartesian co-ords and Angstrom units. atcart(3,nat)
   integer               :: num_bands      ! number of bands left after exclusions
   character(len=3), allocatable :: atsym(:) ! atomic symbols. atsym(nat)
   integer               :: num_nnmax=12
   complex(DP), allocatable :: m_mat(:,:,:,:), a_mat(:,:,:)
   complex(DP), allocatable :: u_mat(:,:,:), u_mat_opt(:,:,:)
   logical, allocatable     :: lwindow(:,:)
   real(DP), allocatable    :: wann_centers(:,:),wann_spreads(:)
   real(DP)                 :: spreads(3)
   real(DP), allocatable    :: eigval(:,:)
   logical                  :: old_spinor_proj  ! for compatability for nnkp files prior to W90v2.0
   integer,allocatable :: rir(:,:)
   logical,allocatable :: zerophase(:,:)
   ! wannier2odd additional variables
   logical :: wannier_plot
   character(len=255) :: wannier_plot_list
   integer, allocatable :: wann_to_plot(:)
   logical :: gamma_trick           ! determines whether or not using SC real wfc (wannier2odd mode)
   logical :: print_rho             ! determines whether or not writing the supercell charge density to file
end module wannier
