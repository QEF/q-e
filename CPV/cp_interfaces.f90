!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! written by Carlo Cavazzoni

!=----------------------------------------------------------------------------=!
   MODULE cp_interfaces
!=----------------------------------------------------------------------------=!

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: bessel2
   PUBLIC :: bessel3
   PUBLIC :: dforce

   PUBLIC :: nlin
   PUBLIC :: nlin_stress
   PUBLIC :: pseudopotential_indexes
   PUBLIC :: compute_dvan
   PUBLIC :: compute_betagx
   PUBLIC :: compute_qradx
   PUBLIC :: interpolate_beta
   PUBLIC :: interpolate_qradb
   PUBLIC :: exact_beta
   PUBLIC :: build_cctab
   PUBLIC :: chkpstab
   PUBLIC :: build_pstab
   PUBLIC :: check_tables
   PUBLIC :: fill_qrl
   PUBLIC :: exact_qradb
   PUBLIC :: compute_xgtab
   PUBLIC :: build_nltab

   PUBLIC :: rhoofr
   PUBLIC :: fillgrad
   PUBLIC :: checkrho
   PUBLIC :: dft_total_charge

   PUBLIC :: writefile
   PUBLIC :: readfile

   PUBLIC :: main_fpmd

   PUBLIC :: runcp_uspp
   PUBLIC :: runcp_uspp_force_pairing

   PUBLIC :: newrho

   PUBLIC :: readempty
   PUBLIC :: writeempty
   PUBLIC :: gram_empty
   PUBLIC :: empty_cp

   PUBLIC :: invfft
   PUBLIC :: fwfft

   PUBLIC :: eigs
   PUBLIC :: fermi_energy
   PUBLIC :: packgam

   PUBLIC :: ortho
   PUBLIC :: ortho_gamma

   PUBLIC :: v2gc
   PUBLIC :: exch_corr_energy
   PUBLIC :: stress_xc
   PUBLIC :: stress_gc

   PUBLIC :: nlrh
   PUBLIC :: nlfh

   PUBLIC :: pstress
   PUBLIC :: pseudo_stress
   PUBLIC :: compute_gagb
   PUBLIC :: stress_har
   PUBLIC :: stress_hartree
   PUBLIC :: add_drhoph
   PUBLIC :: stress_local
   PUBLIC :: stress_kin
   PUBLIC :: stress_nl

   PUBLIC :: interpolate_lambda
   PUBLIC :: update_lambda
   PUBLIC :: elec_fakekine
   PUBLIC :: update_wave_functions
   PUBLIC :: wave_rand_init
   PUBLIC :: kohn_sham
   PUBLIC :: crot
   PUBLIC :: proj

   PUBLIC :: phfacs
   PUBLIC :: strucf

   PUBLIC :: add_core_charge
   PUBLIC :: core_charge_forces

   PUBLIC :: printout_new
   PUBLIC :: printout
   PUBLIC :: print_sfac
   PUBLIC :: open_and_append
   PUBLIC :: cp_print_rho


   PUBLIC :: vofmean
   PUBLIC :: kspotential
   PUBLIC :: vofrhos
   PUBLIC :: vofps
   PUBLIC :: vofloc
   PUBLIC :: force_loc
   PUBLIC :: self_vofhar
   PUBLIC :: localisation
   !
   PUBLIC :: n_atom_wfc
   !
   PUBLIC :: set_eitot
   PUBLIC :: set_evtot
   !
   PUBLIC :: print_projwfc
   PUBLIC :: print_lambda
   !
   PUBLIC :: from_restart_x
   !
   PUBLIC :: move_electrons
   !
   PUBLIC :: compute_stress

   PUBLIC :: protate

   ! ------------------------------------ !


   INTERFACE bessel2
      SUBROUTINE bessel2_x(XG, RW, FINT, LNL, INDL, MMAX)
         USE kinds,     ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN)  :: XG
         REAL(DP), INTENT(IN)  :: RW(:)
         REAL(DP), INTENT(OUT) :: FINT(:,:)
         INTEGER,   INTENT(IN)  :: INDL(:), LNL, MMAX
      END SUBROUTINE bessel2_x
   END INTERFACE

   INTERFACE bessel3
      SUBROUTINE BESSEL3_x(XG, RW, FINT, LNL, INDL, MMAX)
         USE kinds,     ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN)  ::  XG 
         REAL(DP), INTENT(IN)  ::  RW(:)
         REAL(DP), INTENT(OUT)  ::  FINT(:,:)
         INTEGER, INTENT(IN) ::  INDL(:), LNL, MMAX
      END SUBROUTINE BESSEL3_x
   END INTERFACE


   INTERFACE dforce 
      SUBROUTINE dforce_x( i, bec, vkb, c, df, da, v, ldv, ispin, f, n, nspin, v1 )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)    :: i
         REAL(DP)                   :: bec(:,:)
         COMPLEX(DP)                :: vkb(:,:)
         COMPLEX(DP)                :: c(:,:)
         COMPLEX(DP)                :: df(:), da(:)
         INTEGER,     INTENT(IN)    :: ldv
         REAL(DP)                   :: v( ldv, * )
         INTEGER                    :: ispin( : )
         REAL(DP)                   :: f( : )
         INTEGER,     INTENT(IN)    :: n, nspin
         REAL(DP),    OPTIONAL      :: v1( ldv, * )
      END SUBROUTINE dforce_x
   END INTERFACE

   INTERFACE nlin
      SUBROUTINE nlin_x( wsg, wnl )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: wsg( :, : )
         REAL(DP), INTENT(OUT) :: wnl( :, :, :, : )
      END SUBROUTINE nlin_x
   END INTERFACE


   INTERFACE nlin_stress
      SUBROUTINE nlin_stress_x( wnla )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: wnla(:,:,:)
      END SUBROUTINE nlin_stress_x
   END INTERFACE


   INTERFACE pseudopotential_indexes
      SUBROUTINE pseudopotential_indexes_x( nlcc_any )
         IMPLICIT NONE
         LOGICAL, INTENT(OUT) :: nlcc_any
      END SUBROUTINE pseudopotential_indexes_x
   END INTERFACE


   INTERFACE compute_dvan
      SUBROUTINE compute_dvan_x()
         IMPLICIT NONE
      END SUBROUTINE
   END INTERFACE

   INTERFACE compute_betagx
      SUBROUTINE compute_betagx_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE compute_qradx
      SUBROUTINE compute_qradx_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE interpolate_beta
      SUBROUTINE interpolate_beta_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE interpolate_qradb
      SUBROUTINE interpolate_qradb_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE exact_beta
      SUBROUTINE exact_beta_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE build_cctab
      SUBROUTINE build_cctab_x( )
         IMPLICIT NONE
      END SUBROUTINE
   END INTERFACE

   INTERFACE chkpstab
      LOGICAL FUNCTION chkpstab_x(hg, xgtabmax)
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: hg(:)
         REAL(DP), INTENT(IN) :: xgtabmax
      END FUNCTION
   END INTERFACE

   INTERFACE build_pstab
      SUBROUTINE build_pstab_x( )
         IMPLICIT NONE
      END SUBROUTINE
   END INTERFACE

   INTERFACE build_nltab
      SUBROUTINE build_nltab_x( )
         IMPLICIT NONE
      END SUBROUTINE build_nltab_x
   END INTERFACE

   INTERFACE check_tables
      LOGICAL FUNCTION check_tables_x( )
         IMPLICIT NONE
      END FUNCTION check_tables_x
   END INTERFACE

   INTERFACE fill_qrl
      SUBROUTINE fill_qrl_x( is, qrl )
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         INTEGER,  INTENT(IN)  :: is
         REAL(DP), INTENT(OUT) :: qrl( :, :, : )
      END SUBROUTINE
   END INTERFACE

   INTERFACE exact_qradb
      SUBROUTINE exact_qradb_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE compute_xgtab
      SUBROUTINE compute_xgtab_x( xgmin, xgmax, xgtabmax )
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         REAL(DP), INTENT(OUT)  :: xgmax, xgmin, xgtabmax
      END SUBROUTINE
   END INTERFACE


   INTERFACE dft_total_charge
      FUNCTION dft_total_charge_x( c, ngw, fi, n )
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         INTEGER,     INTENT(IN) :: ngw, n
         COMPLEX(DP), INTENT(IN) :: c(:,:)
         REAL (DP),   INTENT(IN) :: fi(:)
         REAL(DP) dft_total_charge_x
      END FUNCTION
   END INTERFACE

  
   INTERFACE rhoofr
      SUBROUTINE rhoofr_fpmd (nfi, tstress, c0, fi, rhor, omega, ekin, dekin)
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         INTEGER,              INTENT(IN) :: nfi
         LOGICAL,              INTENT(IN) :: tstress
         COMPLEX(DP)                     :: c0(:,:)
         REAL(DP),          INTENT(IN) :: fi(:)
         REAL(DP),            INTENT(OUT) :: rhor(:,:)
         REAL(DP),          INTENT(IN) :: omega
         REAL(DP),          INTENT(OUT) :: ekin
         REAL(DP),          INTENT(OUT) :: dekin(6)
      END SUBROUTINE rhoofr_fpmd
      SUBROUTINE rhoofr_cp &
         ( nfi, c, irb, eigrb, bec, rhovan, rhor, rhog, rhos, enl, denl, ekin, dekin )
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         INTEGER nfi
         COMPLEX(DP) c( :, : )
         INTEGER irb( :, : )
         COMPLEX(DP) eigrb( :, : )
         REAL(DP) bec(:,:)
         REAL(DP) rhovan(:, :, : )
         REAL(DP) rhor(:,:)
         COMPLEX(DP) rhog( :, : )
         REAL(DP) rhos(:,:)
         REAL(DP) enl, ekin
         REAL(DP) denl(3,3), dekin(6)
      END SUBROUTINE rhoofr_cp
   END INTERFACE


   INTERFACE fillgrad
      SUBROUTINE fillgrad_x( nspin, rhog, gradr )
         USE kinds,           ONLY: DP         
         USE gvecp,           ONLY: ngm
         USE grid_dimensions, ONLY: nnrx
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nspin
         complex(DP) :: rhog( ngm, nspin )
         real(DP)    :: gradr( nnrx, 3, nspin )
      END SUBROUTINE fillgrad_x
   END INTERFACE


   INTERFACE checkrho
      SUBROUTINE checkrho_x(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
         USE kinds,           ONLY: DP         
         USE grid_dimensions, ONLY: nnrx
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nnr, nspin
         REAL(DP) :: rhor(nnr,nspin), rmin, rmax, rsum, rnegsum
      END SUBROUTINE checkrho_x
   END INTERFACE

   INTERFACE readfile
      SUBROUTINE readfile_cp                                         &
      &     ( flag, ndr,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
      &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
      &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,&
      &       fion, tps, mat_z, occ_f )
         USE kinds,          ONLY : DP
         IMPLICIT NONE
         INTEGER, INTENT(in) :: flag
         integer :: ndr, nfi
         REAL(DP) :: h(3,3), hold(3,3)
         complex(DP) :: c0(:,:), cm(:,:)
         REAL(DP) :: tausm(:,:),taus(:,:), fion(:,:)
         REAL(DP) :: vels(:,:), velsm(:,:) 
         REAL(DP) :: acc(:),lambda(:,:,:), lambdam(:,:,:)
         REAL(DP) :: xnhe0,xnhem,vnhe
         REAL(DP) :: xnhp0(:), xnhpm(:), vnhp(:)
         integer, INTENT(inout) :: nhpcl,nhpdim
         REAL(DP) :: ekincm
         REAL(DP) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
         REAL(DP), INTENT(in) :: ecut, ecutw, delt
         REAL(DP), INTENT(in) :: pmass(:)
         REAL(DP), INTENT(in) :: celldm(6)
         integer, INTENT(in) :: ibrav
         REAL(DP), INTENT(OUT) :: tps
         REAL(DP), INTENT(INOUT) :: mat_z(:,:,:), occ_f(:)
      END SUBROUTINE readfile_cp
      SUBROUTINE readfile_fpmd  &
         ( nfi, trutime, c0, cm, occ, atoms_0, atoms_m, acc, taui, cdmi, ht_m, &
           ht_0, rho, vpot, lambda )
         USE kinds,             ONLY: DP
         USE cell_base,         ONLY: boxdimensions
         USE atoms_type_module, ONLY: atoms_type
         IMPLICIT NONE
         INTEGER,              INTENT(OUT)   :: nfi
         COMPLEX(DP),          INTENT(INOUT) :: c0(:,:), cm(:,:)
         REAL(DP),             INTENT(INOUT) :: occ(:)
         TYPE (boxdimensions), INTENT(INOUT) :: ht_m, ht_0
         TYPE (atoms_type),    INTENT(INOUT) :: atoms_0, atoms_m
         REAL(DP),             INTENT(INOUT) :: rho(:,:)
         REAL(DP),             INTENT(INOUT) :: vpot(:,:)
         REAL(DP),             INTENT(OUT)   :: taui(:,:)
         REAL(DP),             INTENT(OUT)   :: acc(:), cdmi(:)
         REAL(DP),             INTENT(OUT)   :: trutime
         REAL(DP),             INTENT(OUT)   :: lambda(:,:,:)
      END SUBROUTINE readfile_fpmd
   END INTERFACE


   INTERFACE writefile
      SUBROUTINE writefile_cp &
      &     ( ndw,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
      &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
      &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm, &
      &       fion, tps, mat_z, occ_f, rho )
         USE kinds,            ONLY: DP
         implicit none
          integer, INTENT(IN) :: ndw, nfi
         REAL(DP), INTENT(IN) :: h(3,3), hold(3,3)
         complex(DP), INTENT(IN) :: c0(:,:), cm(:,:)
         REAL(DP), INTENT(IN) :: tausm(:,:), taus(:,:), fion(:,:)
         REAL(DP), INTENT(IN) :: vels(:,:), velsm(:,:)
         REAL(DP), INTENT(IN) :: acc(:), lambda(:,:,:), lambdam(:,:,:)
         REAL(DP), INTENT(IN) :: xnhe0, xnhem, vnhe, ekincm
         REAL(DP), INTENT(IN) :: xnhp0(:), xnhpm(:), vnhp(:)
         integer,      INTENT(in) :: nhpcl, nhpdim
         REAL(DP), INTENT(IN) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
         REAL(DP), INTENT(in) :: ecut, ecutw, delt
         REAL(DP), INTENT(in) :: pmass(:)
         REAL(DP), INTENT(in) :: celldm(:)
         REAL(DP), INTENT(in) :: tps
         REAL(DP), INTENT(in) :: rho(:,:)
         integer, INTENT(in) :: ibrav
         REAL(DP), INTENT(in) :: occ_f(:)
         REAL(DP), INTENT(in) :: mat_z(:,:,:)
      END SUBROUTINE writefile_cp
      SUBROUTINE writefile_fpmd  &
         ( nfi, trutime, c0, cm, occ, atoms_0, atoms_m, acc, taui, cdmi, ht_m, &
           ht_0, rho, vpot, lambda, tlast )
         USE kinds,             ONLY: DP
          USE cell_base,         ONLY: boxdimensions
         USE atoms_type_module, ONLY: atoms_type
         IMPLICIT NONE
         INTEGER,              INTENT(IN)    :: nfi
         COMPLEX(DP),          INTENT(IN)    :: c0(:,:), cm(:,:)
         REAL(DP),             INTENT(IN)    :: occ(:)
         TYPE (boxdimensions), INTENT(IN)    :: ht_m, ht_0
         TYPE (atoms_type),    INTENT(IN)    :: atoms_0, atoms_m
         REAL(DP),             INTENT(IN)    :: rho(:,:)
         REAL(DP),             INTENT(INOUT) :: vpot(:,:)
         REAL(DP),             INTENT(IN)    :: taui(:,:)
         REAL(DP),             INTENT(IN)    :: acc(:), cdmi(:)
         REAL(DP),             INTENT(IN)    :: trutime
         REAL(DP),             INTENT(IN)    :: lambda(:,:,:)
         LOGICAL,              INTENT(IN)    :: tlast
      END SUBROUTINE writefile_fpmd
   END INTERFACE
 

   INTERFACE main_fpmd
      SUBROUTINE cpmain_x( tau, fion, etot )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         REAL(DP) :: tau( :, : )
         REAL(DP) :: fion( :, : )
         REAL(DP) :: etot
      END SUBROUTINE cpmain_x
   END INTERFACE


   INTERFACE runcp_uspp
      SUBROUTINE runcp_uspp_x &
         ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, fromscra, restart )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         integer, intent(in) :: nfi
         real(DP) :: fccc, ccc
         real(DP) :: ema0bg(:), dt2bye
         real(DP) :: rhos(:,:)
         real(DP) :: bec(:,:)
         complex(DP) :: c0(:,:), cm(:,:)
         logical, optional, intent(in) :: fromscra
         logical, optional, intent(in) :: restart
      END SUBROUTINE
   END INTERFACE


   INTERFACE runcp_uspp_force_pairing
      SUBROUTINE runcp_uspp_force_pairing_x  &
         ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, intermed, fromscra, &
           restart )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(in) :: nfi
         REAL(DP) :: fccc, ccc
         REAL(DP) :: ema0bg(:), dt2bye
         REAL(DP) :: rhos(:,:)
         REAL(DP) :: bec(:,:)
         COMPLEX(DP) :: c0(:,:), cm(:,:)
         REAL(DP) :: intermed
         LOGICAL, OPTIONAL, INTENT(in) :: fromscra
         LOGICAL, OPTIONAL, INTENT(in) :: restart
      END SUBROUTINE
   END INTERFACE



   INTERFACE newrho
      SUBROUTINE newrho_x( rhor, drho, nfi )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(INOUT) :: rhor(:)
         REAL(DP), INTENT(OUT) ::  drho
         INTEGER, INTENT(IN) :: nfi
      END SUBROUTINE
   END INTERFACE


   INTERFACE readempty
      LOGICAL FUNCTION readempty_x( c_emp, ne )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(OUT) :: c_emp(:,:)
         INTEGER,     INTENT(IN)    :: ne
      END FUNCTION
   END INTERFACE


   INTERFACE writeempty
      SUBROUTINE writeempty_x( c_emp, ne )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: c_emp(:,:)
         INTEGER,     INTENT(IN) :: ne
      END SUBROUTINE
   END INTERFACE


   INTERFACE gram_empty
      SUBROUTINE gram_empty_x  &
         ( tortho, eigr, betae, bec_emp, bec_occ, nkbx, c_emp, c_occ, ngwx, n_emp, n_occ )
         USE kinds,          ONLY: DP
         USE ions_base,      ONLY : nat
         USE uspp,           ONLY : nkb
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nkbx, ngwx, n_emp, n_occ
         COMPLEX(DP)   :: eigr(ngwx,nat)
         REAL(DP)      :: bec_emp( nkbx, n_emp )
         REAL(DP)      :: bec_occ( nkbx, n_occ )
         COMPLEX(DP)   :: betae( ngwx, nkb )
         COMPLEX(DP)   :: c_emp( ngwx, n_emp )
         COMPLEX(DP)   :: c_occ( ngwx, n_occ )
         LOGICAL, INTENT(IN) :: tortho
      END SUBROUTINE
   END INTERFACE


   INTERFACE empty_cp
      SUBROUTINE empty_cp_x ( nfi, c0, v )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         INTEGER,    INTENT(IN) :: nfi
         COMPLEX(DP)            :: c0(:,:)
         REAL(DP)               :: v(:,:)
      END SUBROUTINE
   END INTERFACE


   INTERFACE invfft
      SUBROUTINE invfft_x( grid_type, f, nr1, nr2, nr3, nr1x, nr2x, nr3x, ia )
         USE kinds,      ONLY: DP
         IMPLICIT NONE
         INTEGER,           INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x
         INTEGER, OPTIONAL, INTENT(IN) :: ia
         CHARACTER(LEN=*),  INTENT(IN) :: grid_type
         COMPLEX(DP) :: f(:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE fwfft
      SUBROUTINE fwfft_x( grid_type, f, nr1, nr2, nr3, nr1x, nr2x, nr3x )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x
         CHARACTER(LEN=*), INTENT(IN) :: grid_type
         COMPLEX(DP) :: f(:)
      END SUBROUTINE
   END INTERFACE


   INTERFACE eigs
      SUBROUTINE cp_eigs_x( nfi, lambdap, lambda )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         INTEGER :: nfi
         REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
      END SUBROUTINE
   END INTERFACE

   INTERFACE fermi_energy
      SUBROUTINE fermi_energy_x(eig, occ, wke, ef, qtot, temp, sume)
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         REAL(DP) :: occ(:) 
         REAL(DP) ef, qtot, temp, sume
         REAL(DP) eig(:,:), wke(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE packgam
      SUBROUTINE rpackgam_x( gam, f, aux )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(INOUT)  :: gam(:,:)
         REAL(DP), INTENT(OUT), OPTIONAL :: aux(:)
         REAL(DP), INTENT(IN)  :: f(:)
      END SUBROUTINE
   END INTERFACE


   INTERFACE ortho
      SUBROUTINE ortho_m &
         ( c0, cp, lambda, descla, ccc, nupdwn, iupdwn, nspin )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)    :: descla(:,:)
         INTEGER,     INTENT(IN)    :: nupdwn(:), iupdwn(:), nspin
         COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cp(:,:)
         REAL(DP),    INTENT(INOUT) :: lambda(:,:,:)
         REAL(DP),    INTENT(IN)    :: ccc
      END SUBROUTINE
      SUBROUTINE ortho_cp &
         ( eigr, cp, phi, ngwx, x0, descla, diff, iter, ccc, bephi, becp, nbsp, nspin, nupdwn, iupdwn)
         USE kinds,          ONLY: DP
         USE ions_base,      ONLY: nat
         USE uspp,           ONLY: nkb
         IMPLICIT NONE
         INTEGER,    INTENT(IN)     :: ngwx, nbsp, nspin
         INTEGER,    INTENT(IN)     :: nupdwn( nspin ), iupdwn( nspin )
         INTEGER,     INTENT(IN)    :: descla(:,:)
         COMPLEX(DP) :: cp(ngwx,nbsp), phi(ngwx,nbsp), eigr(ngwx,nat)
         REAL(DP)    :: x0( :, :, : ), diff, ccc
         INTEGER     :: iter
         REAL(DP)    :: bephi(nkb,nbsp), becp(nkb,nbsp)
      END SUBROUTINE
   END INTERFACE

   INTERFACE ortho_gamma
      SUBROUTINE ortho_gamma_x &
         ( iopt, cp, ngwx, phi, becp, qbecp, nkbx, bephi, qbephi, &
           x0, nx0, descla, diff, iter, n, nss, istart )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         INTEGER,  INTENT(IN)  :: iopt
         INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
         INTEGER,  INTENT(IN)  :: n, nss, istart
         COMPLEX(DP) :: phi( ngwx, n ), cp( ngwx, n )
         REAL(DP)    :: bephi( nkbx, n ), becp( nkbx, n )
         REAL(DP)    :: qbephi( nkbx, n ), qbecp( nkbx, n )
         REAL(DP)    :: x0( nx0, nx0 )
         INTEGER,  INTENT(IN)  :: descla( : )
         INTEGER,  INTENT(OUT) :: iter
         REAL(DP), INTENT(OUT) :: diff
      END SUBROUTINE
   END INTERFACE

   INTERFACE v2gc
      SUBROUTINE v2gc_x( v2xc, grho, rhor, vpot )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP) ::  vpot(:,:)
         REAL(DP), intent(in)  ::  v2xc(:,:,:)
         REAL(DP), intent(in)  ::  grho(:,:,:)
         REAL(DP), intent(in)  ::  rhor(:,:)
      END SUBROUTINE
   END INTERFACE
   INTERFACE exch_corr_energy
      SUBROUTINE exch_corr_energy_x(rhoetr, grho, vpot, exc, vxc, v2xc)
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL (DP) :: rhoetr(:,:)
         REAL (DP) :: grho(:,:,:)
         REAL (DP) :: vpot(:,:)
         REAL (DP) :: exc
         REAL (DP) :: vxc
         REAL (DP) :: v2xc(:,:,:)
      END SUBROUTINE
   END INTERFACE
   INTERFACE stress_xc
      SUBROUTINE stress_xc_x &
         ( dexc, strvxc, sfac, vxc, grho, v2xc, gagb, tnlcc, rhocp, box)
         USE kinds,            ONLY: DP
         USE cell_base,        ONLY: boxdimensions
         IMPLICIT NONE
         type (boxdimensions), intent(in) :: box
         LOGICAL :: tnlcc(:)
         COMPLEX(DP) :: vxc(:,:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
         REAL(DP) :: dexc(:), strvxc
         REAL(DP) :: grho(:,:,:)
         REAL(DP) :: v2xc(:,:,:)
         REAL(DP) :: gagb(:,:)
         REAL(DP) :: rhocp(:,:)
      END SUBROUTINE
   END INTERFACE
   INTERFACE stress_gc
      SUBROUTINE stress_gc_x(grho, v2xc, gcpail, omega)
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP) ::  v2xc(:,:,:)
         REAL(DP) ::  grho(:,:,:)
         REAL(DP) ::  gcpail(6)
         REAL(DP) ::  omega
      END SUBROUTINE
   END INTERFACE


   INTERFACE nlrh
      SUBROUTINE nlrh_x( c0, tforce, tstress, fion, bec, becdr, eigr, enl, denl )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP)                 :: eigr(:,:)    
         COMPLEX(DP)                 :: c0(:,:)     
         LOGICAL,     INTENT(IN)     :: tforce     
         LOGICAL,     INTENT(IN)     :: tstress     
         REAL(DP),    INTENT(INOUT)  :: fion(:,:) 
         REAL(DP)                    :: bec(:,:)
         REAL(DP)                    :: becdr(:,:,:)
         REAL(DP),    INTENT(OUT)    :: enl
         REAL(DP),    INTENT(OUT)    :: denl( 6 )
      END SUBROUTINE
   END INTERFACE


   INTERFACE pstress
      SUBROUTINE pstress_x( paiu, desr, dekin, denl, deps, deht, dexc )
         USE kinds,         ONLY: DP
         IMPLICIT NONE
         REAL(DP) :: paiu(3,3)
         REAL(DP) :: desr(6), dekin(6), denl(6), deps(6), deht(6), dexc(6)
      END SUBROUTINE
   END INTERFACE

   INTERFACE pseudo_stress
      SUBROUTINE pseudo_stress_x( deps, epseu, gagb, sfac, dvps, rhoeg, omega )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP),     INTENT(IN)  :: omega
         REAL(DP),     INTENT(OUT) :: deps(:)
         REAL(DP),     INTENT(IN)  :: gagb(:,:)
         COMPLEX(DP),  INTENT(IN)  :: rhoeg(:,:)
         COMPLEX(DP),  INTENT(IN)  :: sfac(:,:)
         REAL(DP),     INTENT(IN)  :: dvps(:,:)
         REAL(DP),     INTENT(IN)  :: epseu
      END SUBROUTINE
   END INTERFACE

   INTERFACE compute_gagb
      SUBROUTINE compute_gagb_x( gagb, gx, ngm, tpiba2 )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,  INTENT(IN)  :: ngm
         REAL(DP), INTENT(IN)  :: gx(:,:)
         REAL(DP), INTENT(OUT) :: gagb(:,:)
         REAL(DP), INTENT(IN)  :: tpiba2
      END SUBROUTINE
   END INTERFACE

   INTERFACE stress_har
      SUBROUTINE stress_har_x(deht, ehr, sfac, rhoeg, gagb, omega )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP),    INTENT(OUT) :: DEHT(:)
         REAL(DP),    INTENT(IN)  :: omega, EHR, gagb(:,:)
         COMPLEX(DP), INTENT(IN)  :: RHOEG(:,:)
         COMPLEX(DP), INTENT(IN)  :: sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE stress_hartree
      SUBROUTINE stress_hartree_x(deht, ehr, sfac, rhot, drhot, gagb, omega )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP),    INTENT(OUT) :: DEHT(:)
         REAL(DP),    INTENT(IN)  :: omega, EHR, gagb(:,:)
         COMPLEX(DP) :: rhot(:)  ! total charge: Sum_spin ( rho_e + rho_I )
         COMPLEX(DP) :: drhot(:,:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE add_drhoph
      SUBROUTINE add_drhoph_x( drhot, sfac, gagb )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(INOUT) :: drhot( :, : )
         COMPLEX(DP), INTENT(IN) :: sfac( :, : )
         REAL(DP),    INTENT(IN) :: gagb( :, : )
      END SUBROUTINE
   END INTERFACE

   INTERFACE stress_local
      SUBROUTINE stress_local_x( deps, epseu, gagb, sfac, rhoe, drhoe, omega )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP),     INTENT(IN)  :: omega
         REAL(DP),     INTENT(OUT) :: deps(:)
         REAL(DP),     INTENT(IN)  :: gagb(:,:)
         COMPLEX(DP),  INTENT(IN)  :: rhoe(:)
         COMPLEX(DP),  INTENT(IN)  :: drhoe(:,:)
         COMPLEX(DP),  INTENT(IN)  :: sfac(:,:)
         REAL(DP),     INTENT(IN)  :: epseu
      END SUBROUTINE
   END INTERFACE

   INTERFACE stress_kin
      SUBROUTINE stress_kin_x(dekin, c0, occ)
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP),    INTENT(OUT) :: dekin(:)
         COMPLEX(DP), INTENT(IN)  :: c0(:,:)
         REAL(DP),    INTENT(IN)  :: occ(:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE stress_nl
      SUBROUTINE stress_nl_x( denl, c0, occ, eigr, bec, enl )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: occ(:)
         COMPLEX(DP), INTENT(IN) :: c0(:,:)
         REAL(DP), INTENT(OUT) :: denl(:)
         COMPLEX(DP), INTENT(IN) :: eigr(:,:)
         REAL(DP) :: bec(:,:)
         REAL(DP), INTENT(IN) :: enl
      END SUBROUTINE
   END INTERFACE


   INTERFACE interpolate_lambda
      SUBROUTINE interpolate_lambda_x( lambdap, lambda, lambdam )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP) :: lambdap(:,:,:), lambda(:,:,:), lambdam(:,:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE update_lambda
      SUBROUTINE update_lambda_x( i, lambda, c0, c2, n, noff, tdist )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: n, noff
         REAL(DP)            :: lambda(:,:)
         COMPLEX(DP)         :: c0(:,:), c2(:)
         INTEGER, INTENT(IN) :: i
         LOGICAL, INTENT(IN) :: tdist   !  if .true. lambda is distributed
      END SUBROUTINE
   END INTERFACE

   INTERFACE elec_fakekine
      SUBROUTINE elec_fakekine_x( ekincm, ema0bg, emass, c0, cm, ngw, n, noff, delt )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         integer, intent(in)      :: ngw    !  number of plane wave coeff.
         integer, intent(in)      :: n      !  number of bands
         integer, intent(in)      :: noff   !  offset for band index
         real(DP), intent(out)    :: ekincm
         real(DP), intent(in)     :: ema0bg( ngw ), delt, emass
         complex(DP), intent(in)  :: c0( ngw, n ), cm( ngw, n )
      END SUBROUTINE
   END INTERFACE

   INTERFACE update_wave_functions
      SUBROUTINE update_wave_functions_x( cm, c0, cp )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN)    :: cp(:,:)
         COMPLEX(DP), INTENT(INOUT) :: c0(:,:)
         COMPLEX(DP), INTENT(OUT)   :: cm(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE wave_rand_init
      SUBROUTINE wave_rand_init_x( cm, n, noff )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)  :: n, noff
         COMPLEX(DP), INTENT(OUT) :: cm(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE kohn_sham
      SUBROUTINE kohn_sham_x( c, ngw, eforces, n, nl, noff )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngw   ! number of plane waves
         INTEGER, INTENT(IN) :: n     ! number of ks states
         INTEGER, INTENT(IN) :: nl    ! local (to the processor) number of states
         INTEGER, INTENT(IN) :: noff  ! offset of the first state in array "c"
         COMPLEX(DP), INTENT(INOUT) ::  c(:,:) 
         COMPLEX(DP) :: eforces(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE crot
      SUBROUTINE crot_gamma ( c0, ngwl, nx, noff, lambda, nrl, eig )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngwl, nx, nrl, noff
         COMPLEX(DP), INTENT(INOUT) :: c0(:,:)
         REAL(DP) :: lambda(:,:)
         REAL(DP) :: eig(:)
      END SUBROUTINE
      SUBROUTINE crot_gamma2 ( c0rot, c0, ngw, n, noffr, noff, lambda, nx, eig )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)    :: ngw, n, nx, noffr, noff
         COMPLEX(DP), INTENT(INOUT) :: c0rot(:,:)
         COMPLEX(DP), INTENT(IN)    :: c0(:,:)
         REAL(DP),    INTENT(IN)    :: lambda(:,:)
         REAL(DP),    INTENT(OUT)   :: eig(:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE proj
      SUBROUTINE proj_gamma( a, b, ngw, n, noff, lambda)
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT( IN )  :: ngw, n, noff
         COMPLEX(DP), INTENT(INOUT) :: a(:,:), b(:,:)
         REAL(DP),    OPTIONAL      :: lambda(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE phfacs
      SUBROUTINE phfacs_x( ei1, ei2, ei3, eigr, mill, taus, nr1, nr2, nr3, nat )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nat
         INTEGER, INTENT(IN) :: nr1, nr2, nr3
         COMPLEX(DP) :: ei1( -nr1 : nr1, nat )
         COMPLEX(DP) :: ei2( -nr2 : nr2, nat )
         COMPLEX(DP) :: ei3( -nr3 : nr3, nat )
         COMPLEX(DP) :: eigr( :, : )
         REAL(DP) :: taus( 3, nat )
         INTEGER      :: mill( :, : )
      END SUBROUTINE
   END INTERFACE

   INTERFACE strucf
      SUBROUTINE strucf_x( sfac, ei1, ei2, ei3, mill, ngm )
         USE kinds,            ONLY: DP
         USE ions_base,        ONLY: nat
         USE grid_dimensions,  ONLY: nr1, nr2, nr3
         IMPLICIT NONE
         COMPLEX(DP) :: ei1( -nr1 : nr1, nat )
         COMPLEX(DP) :: ei2( -nr2 : nr2, nat )
         COMPLEX(DP) :: ei3( -nr3 : nr3, nat )
         INTEGER      :: mill( :, : )
         INTEGER      :: ngm
         COMPLEX(DP), INTENT(OUT) :: sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   
   INTERFACE add_core_charge
      SUBROUTINE add_core_charge_x( rhoetg, rhoetr, sfac, rhoc, nsp)
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         integer :: nsp 
         COMPLEX(DP) :: rhoetg(:)
         REAL(DP)    :: rhoetr(:)
         REAL(DP)    :: rhoc(:,:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE core_charge_forces
      SUBROUTINE core_charge_forces_x &
         ( fion, vxc, rhoc1, tnlcc, atoms, ht, ei1, ei2, ei3 )
         USE kinds,              ONLY: DP
         USE cell_base,          ONLY: boxdimensions
         USE atoms_type_module,  ONLY: atoms_type
         USE grid_dimensions,    ONLY: nr1, nr2, nr3
         USE ions_base,          ONLY: nat
         IMPLICIT NONE
         TYPE (atoms_type), INTENT(IN) :: atoms    !   atomic positions
         TYPE (boxdimensions), INTENT(IN) :: ht    !   cell parameters
         COMPLEX(DP) :: ei1( -nr1:nr1, nat)                  !
         COMPLEX(DP) :: ei2( -nr2:nr2, nat)                  !
         COMPLEX(DP) :: ei3( -nr3:nr3, nat)                  !
         LOGICAL      :: tnlcc(:)                  !   NLCC flags
         REAL(DP)    :: fion(:,:)                 !   ionic forces
         REAL(DP)    :: rhoc1(:,:)                !   derivative of the core charge
         COMPLEX(DP) :: vxc(:,:)                  !   XC potential
      END SUBROUTINE
   END INTERFACE


   INTERFACE printout_new
      SUBROUTINE printout_new_x &
         ( nfi, tfirst, tfilei, tprint, tps, h, stress, tau0, vels, &
           fion, ekinc, temphc, tempp, temps, etot, enthal, econs, econt, &
           vnhh, xnhh0, vnhp, xnhp0, atot, ekin, epot )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nfi
         LOGICAL, INTENT(IN) :: tfirst, tfilei, tprint
         REAL(DP), INTENT(IN) :: tps
         REAL(DP), INTENT(IN) :: h( 3, 3 )
         REAL(DP), INTENT(IN) :: stress( 3, 3 )
         REAL(DP), INTENT(IN) :: tau0( :, : )  ! real positions
         REAL(DP), INTENT(IN) :: vels( :, : )  ! scaled velocities
         REAL(DP), INTENT(IN) :: fion( :, : )  ! real forces
         REAL(DP), INTENT(IN) :: ekinc, temphc, tempp, etot, enthal, econs, econt
         REAL(DP), INTENT(IN) :: temps( : ) ! partial temperature for different ionic species
         REAL(DP), INTENT(IN) :: vnhh( 3, 3 ), xnhh0( 3, 3 ), vnhp( 1 ), xnhp0( 1 )
         REAL(DP), INTENT(IN) :: atot! enthalpy of system for c.g. case
         REAL(DP), INTENT(IN) :: ekin
         REAL(DP), INTENT(IN) :: epot ! ( epseu + eht + exc )
      END SUBROUTINE
   END INTERFACE

   INTERFACE printout
      SUBROUTINE printout_x(nfi, atoms, ekinc, ekcell, tprint, ht, edft)
         USE kinds,             ONLY: DP
         USE atoms_type_module, ONLY: atoms_type
         USE cell_base,         ONLY: boxdimensions
         USE energies,          ONLY: dft_energy_type
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nfi
         TYPE (atoms_type)   :: atoms
         LOGICAL             :: tprint
         type (boxdimensions), intent(in) :: ht
         TYPE (dft_energy_type) :: edft
         REAL(DP) :: ekinc, ekcell
      END SUBROUTINE
   END INTERFACE

   INTERFACE print_sfac
      SUBROUTINE print_sfac_x( rhoe, sfac )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: rhoe(:,:)
         COMPLEX(DP), INTENT(IN) ::  sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE open_and_append
      SUBROUTINE open_and_append_x( iunit, file_name )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: iunit
         CHARACTER(LEN = *), INTENT(IN) :: file_name
      END SUBROUTINE
   END INTERFACE

   INTERFACE cp_print_rho
      SUBROUTINE cp_print_rho_x &
         (nfi, bec, c0, eigr, irb, eigrb, rhor, rhog, rhos, lambdap, lambda, tau0, h )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         INTEGER :: nfi
         INTEGER :: irb(:,:)
         COMPLEX(DP) :: c0( :, : )
         REAL(DP) :: bec( :, : ), rhor( :, : ), rhos( :, : )
         REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
         REAL(DP) :: tau0( :, : ), h( 3, 3 )
         COMPLEX(DP) :: eigrb( :, : ), rhog( :, : )
         COMPLEX(DP) :: eigr( :, : )
      END SUBROUTINE
   END INTERFACE


   INTERFACE vofmean
      SUBROUTINE vofmean_x( sfac, rhops, rhoeg )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         REAL(DP),    INTENT(IN)   :: RHOPS(:,:)
         COMPLEX(DP), INTENT(IN)   :: RHOEG(:)
         COMPLEX(DP), INTENT(IN)   :: sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE kspotential
      SUBROUTINE kspotential_x &
         ( nfi, tprint, tforce, tstress, rhoe, atoms, bec, becdr, eigr, &
          ei1, ei2, ei3, sfac, c0, tcel, ht, fi, vpot, edft )
         USE kinds,             ONLY: DP
         USE energies,          ONLY: dft_energy_type
         USE cell_base,         ONLY: boxdimensions
         USE atoms_type_module, ONLY: atoms_type
         USE wave_types,        ONLY: wave_descriptor 
         IMPLICIT NONE
         INTEGER,              INTENT(IN)    :: nfi
         LOGICAL, INTENT(IN) :: tforce, tstress, tprint
         REAL(DP) :: rhoe(:,:)
         TYPE (atoms_type),    INTENT(INOUT) :: atoms
         REAL(DP) :: bec(:,:)
         REAL(DP) :: becdr(:,:,:)
         COMPLEX(DP) :: eigr(:,:)
         COMPLEX(DP) :: ei1(:,:)
         COMPLEX(DP) :: ei2(:,:)
         COMPLEX(DP) :: ei3(:,:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
         COMPLEX(DP),         INTENT(INOUT) :: c0(:,:)
         LOGICAL   :: tcel
         TYPE (boxdimensions), INTENT(INOUT) ::  ht
         REAL(DP), INTENT(IN) :: fi(:)
         REAL(DP)    :: vpot(:,:)
         TYPE (dft_energy_type) :: edft
      END SUBROUTINE
   END INTERFACE


   INTERFACE vofrhos
      SUBROUTINE vofrhos_x &
         ( tprint, tforce, tstress, rhoe, atoms, vpot, bec, c0, fi, &
           eigr, ei1, ei2, ei3, sfac, box, edft )
         USE kinds,             ONLY: DP
         USE energies,          ONLY: dft_energy_type
         USE cell_base,         ONLY: boxdimensions
         USE atoms_type_module, ONLY: atoms_type
         USE wave_types,        ONLY: wave_descriptor 
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tprint, tforce, tstress
         REAL(DP) :: rhoe(:,:)
         TYPE (atoms_type), INTENT(INOUT) :: atoms
         REAL(DP)            :: vpot(:,:)
         REAL(DP)    :: bec(:,:)
         COMPLEX(DP),   INTENT(IN) :: c0(:,:)
         REAL(DP),    INTENT(IN) :: fi(:)
         COMPLEX(DP) :: eigr(:,:)
         COMPLEX(DP) :: ei1(:,:)
         COMPLEX(DP) :: ei2(:,:)
         COMPLEX(DP) :: ei3(:,:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
         TYPE (boxdimensions),    INTENT(INOUT) :: box 
         TYPE (dft_energy_type) :: edft
      END SUBROUTINE
   END INTERFACE
   

   INTERFACE vofps
      SUBROUTINE vofps_x( eps, vloc, rhoeg, vps, sfac, omega )
         USE kinds,              ONLY: DP 
         IMPLICIT NONE
         REAL(DP),    INTENT(IN)  :: vps(:,:)
         REAL(DP),    INTENT(IN)  :: omega
         COMPLEX(DP), INTENT(OUT) :: vloc(:)
         COMPLEX(DP), INTENT(IN)  :: rhoeg(:)
         COMPLEX(DP), INTENT(IN)  :: sfac(:,:)
         COMPLEX(DP), INTENT(OUT) :: eps
      END SUBROUTINE
   END INTERFACE


   INTERFACE vofloc
      SUBROUTINE vofloc_x( tscreen, ehte, ehti, eh, vloc, rhoeg, &
                     rhops, vps, sfac, omega, screen_coul )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         LOGICAL,     INTENT(IN)    :: tscreen
         REAL(DP),    INTENT(IN)    :: rhops(:,:), vps(:,:)
         COMPLEX(DP), INTENT(INOUT) :: vloc(:)
         COMPLEX(DP), INTENT(IN)    :: rhoeg(:)
         COMPLEX(DP), INTENT(IN)    :: sfac(:,:)
         REAL(DP),    INTENT(OUT)   :: ehte, ehti 
         REAL(DP),    INTENT(IN)    :: omega
         COMPLEX(DP), INTENT(OUT)   :: eh
         COMPLEX(DP), INTENT(IN)    :: screen_coul(:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE force_loc
      SUBROUTINE force_loc_x( tscreen, rhoeg, fion, rhops, vps, eigr, ei1, ei2, ei3, &
                        sfac, omega, screen_coul )
         USE kinds,              ONLY: DP
         USE grid_dimensions,    ONLY: nr1, nr2, nr3
         USE ions_base,          ONLY: nat
         IMPLICIT NONE
         LOGICAL     :: tscreen
         REAL(DP)    :: fion(:,:) 
         REAL(DP)    :: rhops(:,:), vps(:,:)  
         COMPLEX(DP) :: rhoeg(:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
            COMPLEX(DP) :: ei1(-nr1:nr1,nat)
         COMPLEX(DP) :: ei2(-nr2:nr2,nat)
         COMPLEX(DP) :: ei3(-nr3:nr3,nat)
         COMPLEX(DP) :: eigr(:,:)
         REAL(DP)    :: omega
         COMPLEX(DP) :: screen_coul(:)
      END SUBROUTINE
   END INTERFACE
 
   INTERFACE self_vofhar
      SUBROUTINE self_vofhar_x( tscreen, self_ehte, vloc, rhoeg, omega, hmat )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         LOGICAL     :: tscreen
         COMPLEX(DP) :: vloc(:)
         COMPLEX(DP) :: rhoeg(:,:)
         REAL(DP)    :: self_ehte
         REAL(DP), INTENT(IN) :: omega 
         REAL(DP), INTENT(IN) :: hmat( 3, 3 )
      END SUBROUTINE
   END INTERFACE
   
   INTERFACE localisation
      SUBROUTINE localisation_x( wfc, atoms_m, ht)
         USE kinds,              ONLY: DP
         USE cell_base,   ONLY: boxdimensions
         USE atoms_type_module, ONLY: atoms_type
         IMPLICIT NONE 
         COMPLEX(DP), INTENT(IN) :: wfc(:)
         TYPE (atoms_type), INTENT(in) :: atoms_m
         TYPE (boxdimensions), INTENT(in) :: ht
      END SUBROUTINE
   END INTERFACE


   INTERFACE n_atom_wfc
      FUNCTION n_atom_wfc_x()
         INTEGER n_atom_wfc_x
      END FUNCTION
   END INTERFACE

   INTERFACE set_eitot
      SUBROUTINE set_eitot_x( eitot )
         USE kinds, ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: eitot(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE set_evtot
      SUBROUTINE set_evtot_x( c0, ctot, lambda, iupdwn_tot, nupdwn_tot )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN)  :: c0(:,:)
         COMPLEX(DP), INTENT(OUT) :: ctot(:,:)
         REAL(DP),    INTENT(IN)  :: lambda(:,:,:)
         INTEGER,     INTENT(IN)  :: iupdwn_tot(2), nupdwn_tot(2)
      END SUBROUTINE
   END INTERFACE

   INTERFACE print_projwfc
      SUBROUTINE print_projwfc_x ( c0, lambda, eigr, vkb )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN)  :: c0(:,:), eigr(:,:), vkb(:,:)
         REAL(DP),    INTENT(IN)  :: lambda(:,:,:)
      END SUBROUTINE
   END INTERFACE



   INTERFACE from_restart_x
      SUBROUTINE from_restart_x( &
         ht0, htm, phi, c0, cm, lambdap, lambda, lambdam, ei1, ei2, ei3, eigr, &
         sfac, vkb, nkb, bec, dbec, taub, irb, eigrb )
         USE kinds,         ONLY: DP
         USE cell_base,     ONLY: boxdimensions
         IMPLICIT NONE
         TYPE (boxdimensions) :: ht0, htm
         COMPLEX(DP)          :: phi(:,:)
         COMPLEX(DP)          :: c0(:,:)
         COMPLEX(DP)          :: cm(:,:)
         REAL(DP)             :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
         COMPLEX(DP)          :: ei1(:,:)
         COMPLEX(DP)          :: ei2(:,:)
         COMPLEX(DP)          :: ei3(:,:)
         COMPLEX(DP)          :: eigr(:,:)
         COMPLEX(DP)          :: sfac(:,:)
         COMPLEX(DP)          :: vkb(:,:)
         INTEGER, INTENT(IN)  :: nkb
         REAL(DP)             :: bec(:,:)
         REAL(DP)             :: dbec(:,:,:,:)
         REAL(DP)             :: taub(:,:)
         INTEGER              :: irb(:,:)
         COMPLEX(DP)          :: eigrb(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE move_electrons
      SUBROUTINE move_electrons_x( &
         nfi, tfirst, tlast, b1, b2, b3, fion, enthal, enb, enbi, fccc, ccc, dt2bye, stress )
         USE kinds,         ONLY: DP
         IMPLICIT NONE
         INTEGER,  INTENT(IN)    :: nfi
         LOGICAL,  INTENT(IN)    :: tfirst, tlast
         REAL(DP), INTENT(IN)    :: b1(3), b2(3), b3(3)
         REAL(DP)                :: fion(:,:)
         REAL(DP), INTENT(IN)    :: dt2bye
         REAL(DP)                :: fccc, ccc
         REAL(DP)                :: enb, enbi
         REAL(DP)                :: enthal
         REAL(DP)                :: ei_unp
         REAL(DP)                :: stress(3,3)
      END SUBROUTINE
   END INTERFACE

   INTERFACE compute_stress
      SUBROUTINE compute_stress_x( stress, detot, h, omega )
         USE kinds, ONLY : DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: stress(3,3)
         REAL(DP), INTENT(IN)  :: detot(3,3), h(3,3), omega
      END SUBROUTINE
   END INTERFACE

   INTERFACE nlfh
      SUBROUTINE nlfh_x( stress, bec, dbec, lambda )
         USE kinds, ONLY : DP
         IMPLICIT NONE
         REAL(DP), INTENT(INOUT) :: stress(3,3) 
         REAL(DP), INTENT(IN)    :: bec( :, : ), dbec( :, :, :, : )
         REAL(DP), INTENT(IN)    :: lambda( :, :, : )
      END SUBROUTINE
   END INTERFACE

   INTERFACE print_lambda
      SUBROUTINE print_lambda_x( lambda, n, nshow, ccc, iunit )
         USE kinds, ONLY : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: lambda(:,:,:), ccc
         INTEGER, INTENT(IN) :: n, nshow
         INTEGER, INTENT(IN), OPTIONAL :: iunit
      END SUBROUTINE
   END INTERFACE

   INTERFACE protate
      SUBROUTINE protate_x ( c0, bec, c0rot, becrot, ngwl, nss, noff, lambda, nrl, &
                           na, nsp, ish, nh, np_rot, me_rot, comm_rot  )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngwl, nss, nrl, noff
         INTEGER, INTENT(IN) :: na(:), nsp, ish(:), nh(:)
         INTEGER, INTENT(IN) :: np_rot, me_rot, comm_rot  
         COMPLEX(DP), INTENT(IN) :: c0(:,:)
         COMPLEX(DP), INTENT(OUT) :: c0rot(:,:)
         REAL(DP), INTENT(IN) :: lambda(:,:)
         REAL(DP), INTENT(IN) :: bec(:,:)
         REAL(DP), INTENT(OUT) :: becrot(:,:)
      END SUBROUTINE

   END INTERFACE

!=----------------------------------------------------------------------------=!
   END MODULE
!=----------------------------------------------------------------------------=!
