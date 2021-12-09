!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! written by Carlo Cavazzoni

#if defined(__CUDA)
#define DEVICEATTR ,DEVICE
#else
#define DEVICEATTR
#endif

#if defined(__CUDA)
#define PINMEM 
#else
#define PINMEM
#endif

!=----------------------------------------------------------------------------=!
   MODULE cp_interfaces
!=----------------------------------------------------------------------------=!

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: dforce

   PUBLIC :: pseudopotential_indexes
   PUBLIC :: compute_dvan
   PUBLIC :: compute_betagx
   PUBLIC :: compute_qradx
   PUBLIC :: interpolate_beta
   PUBLIC :: interpolate_qradb
   PUBLIC :: exact_beta
   PUBLIC :: build_cctab
   PUBLIC :: build_pstab
   PUBLIC :: check_tables
   PUBLIC :: exact_qradb
   PUBLIC :: compute_xgtab

   PUBLIC :: rhoofr
   PUBLIC :: checkrho

   PUBLIC :: writefile
   PUBLIC :: readfile

   PUBLIC :: runcp_uspp
   PUBLIC :: runcp_uspp_force_pairing

   PUBLIC :: eigs
   PUBLIC :: fermi_energy
   PUBLIC :: packgam

   PUBLIC :: ortho

   PUBLIC :: nlfh
   PUBLIC :: nlfl_bgrp

   PUBLIC :: pseudo_stress
   PUBLIC :: compute_gagb
   PUBLIC :: stress_har
   PUBLIC :: stress_hartree
   PUBLIC :: add_drhoph
   PUBLIC :: stress_local
   PUBLIC :: stress_kin

   PUBLIC :: interpolate_lambda
   PUBLIC :: update_lambda
   PUBLIC :: elec_fakekine
   PUBLIC :: wave_rand_init
   PUBLIC :: crot
   PUBLIC :: proj

   PUBLIC :: phfacs
   PUBLIC :: strucf

   PUBLIC :: printout_new
   PUBLIC :: open_and_append
   PUBLIC :: cp_print_rho


   PUBLIC :: vofmean
   PUBLIC :: vofps
   PUBLIC :: vofloc
   PUBLIC :: force_loc
   PUBLIC :: self_vofhar
   !
   PUBLIC :: set_eitot
   PUBLIC :: set_evtot
   !
   PUBLIC :: move_electrons
   !
   PUBLIC :: compute_stress

   PUBLIC :: protate
   PUBLIC :: c_bgrp_expand
   PUBLIC :: c_bgrp_pack
   PUBLIC :: vofrho
   PUBLIC :: enkin
   PUBLIC :: newinit
   PUBLIC :: prefor
   PUBLIC :: denlcc
   PUBLIC :: dotcsc
   PUBLIC :: nlsm1
   PUBLIC :: nlsm2_bgrp
   PUBLIC :: ennl
   PUBLIC :: calrhovan
   PUBLIC :: calbec
   PUBLIC :: caldbec_bgrp
   PUBLIC :: dennl
   PUBLIC :: nlfq_bgrp
   PUBLIC :: collect_bec
   PUBLIC :: nlsm1us
   PUBLIC :: dbeta_eigr

   ! ------------------------------------ !


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
#if defined (__CUDA)
      SUBROUTINE dforce_gpu_x( i, bec, vkb, c, df, da, v, ldv, ispin, f, n, nspin )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)    :: i
         REAL(DP)                   :: bec(:,:)
         COMPLEX(DP), DEVICE        :: vkb(:,:)
         COMPLEX(DP), DEVICE        :: c(:,:)
         COMPLEX(DP)                :: df(:), da(:)
         INTEGER,     INTENT(IN)    :: ldv
         REAL(DP), DEVICE           :: v( :, : )
         INTEGER                    :: ispin( : )
         REAL(DP)                   :: f( : )
         INTEGER,     INTENT(IN)    :: n, nspin
      END SUBROUTINE dforce_gpu_x
#endif
   END INTERFACE



   INTERFACE pseudopotential_indexes
      SUBROUTINE pseudopotential_indexes_x( )
         IMPLICIT NONE
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

   INTERFACE build_pstab
      SUBROUTINE build_pstab_x( )
         IMPLICIT NONE
      END SUBROUTINE
   END INTERFACE

   INTERFACE check_tables
      LOGICAL FUNCTION check_tables_x( gmax )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: gmax
      END FUNCTION check_tables_x
   END INTERFACE

   INTERFACE exact_qradb
      SUBROUTINE exact_qradb_x( tpre )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: tpre
      END SUBROUTINE
   END INTERFACE

   INTERFACE compute_xgtab
      SUBROUTINE compute_xgtab_x( xgmin, xgmax )
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         REAL(DP), INTENT(OUT)  :: xgmax, xgmin
      END SUBROUTINE
   END INTERFACE


   INTERFACE rhoofr
      SUBROUTINE rhoofr_cp &
         ( nfi, c_bgrp, c_d, bec, dbec, rhovan, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin, tstress, ndwwf )
         USE kinds,      ONLY: DP         
         IMPLICIT NONE
         INTEGER nfi
         COMPLEX(DP) :: c_bgrp( :, : )
         COMPLEX(DP) DEVICEATTR :: c_d( :, : )
         REAL(DP) bec(:,:)
         REAL(DP) dbec(:,:,:,:)
         REAL(DP) rhovan(:, :, : )
         REAL(DP) rhor(:,:)
         REAL(DP) drhor(:,:,:,:)
         COMPLEX(DP) rhog( :, : )
         COMPLEX(DP) drhog( :, :, :, : )
         REAL(DP) rhos(:,:)
         REAL(DP) enl, ekin
         REAL(DP) denl(3,3), dekin(6)
         LOGICAL, OPTIONAL, INTENT(IN) :: tstress
         INTEGER, OPTIONAL, INTENT(IN) :: ndwwf
      END SUBROUTINE rhoofr_cp
      SUBROUTINE rhoofr_host &
         ( nfi, c_bgrp, irb, eigrb, bec, dbec, rhovan, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin, tstress, ndwwf )
         USE kinds,      ONLY: DP
         IMPLICIT NONE
         INTEGER nfi
         COMPLEX(DP) c_bgrp( :, : )
         INTEGER irb( :, : )
         COMPLEX(DP) eigrb( :, : )
         REAL(DP) bec(:,:)
         REAL(DP) dbec(:,:,:,:)
         REAL(DP) rhovan(:, :, : )
         REAL(DP) rhor(:,:)
         REAL(DP) drhor(:,:,:,:)
         COMPLEX(DP) rhog( :, : )
         COMPLEX(DP) drhog( :, :, :, : )
         REAL(DP) rhos(:,:)
         REAL(DP) enl, ekin
         REAL(DP) denl(3,3), dekin(6)
         LOGICAL, OPTIONAL, INTENT(IN) :: tstress
         INTEGER, OPTIONAL, INTENT(IN) :: ndwwf
      END SUBROUTINE rhoofr_host
   END INTERFACE

   INTERFACE checkrho
      SUBROUTINE checkrho_x(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
         USE kinds,           ONLY: DP         
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nnr, nspin
         REAL(DP) :: rhor(nnr,nspin), rmin, rmax, rsum, rnegsum
      END SUBROUTINE checkrho_x
   END INTERFACE

   INTERFACE readfile
      SUBROUTINE readfile_x                                         &
      &     ( flag, h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
      &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
      &       xnhh0,xnhhm,vnhh,velh,fion, tps, mat_z, occ_f )
         USE kinds,          ONLY : DP
         IMPLICIT NONE
         INTEGER, INTENT(in) :: flag
         integer ::  nfi
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
         REAL(DP), INTENT(OUT) :: tps
         REAL(DP), INTENT(INOUT) :: mat_z(:,:,:), occ_f(:)
      END SUBROUTINE readfile_x
   END INTERFACE


   INTERFACE writefile
      SUBROUTINE writefile_x &
      &     ( h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
      &       lambda,lambdam,idesc,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
      &       xnhh0,xnhhm,vnhh,velh, fion, tps, mat_z, occ_f, rho )
         USE kinds,            ONLY: DP
         implicit none
         integer, INTENT(IN) :: nfi
         REAL(DP), INTENT(IN) :: h(3,3), hold(3,3)
         complex(DP), INTENT(IN) :: c0(:,:), cm(:,:)
         REAL(DP), INTENT(IN) :: tausm(:,:), taus(:,:), fion(:,:)
         REAL(DP), INTENT(IN) :: vels(:,:), velsm(:,:)
         REAL(DP), INTENT(IN) :: acc(:), lambda(:,:,:), lambdam(:,:,:)
         INTEGER,  INTENT(IN) :: idesc( :, : )
         REAL(DP), INTENT(IN) :: xnhe0, xnhem, vnhe, ekincm
         REAL(DP), INTENT(IN) :: xnhp0(:), xnhpm(:), vnhp(:)
         integer,      INTENT(in) :: nhpcl, nhpdim
         REAL(DP), INTENT(IN) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
         REAL(DP), INTENT(in) :: tps
         REAL(DP), INTENT(in) :: rho(:,:)
         REAL(DP), INTENT(in) :: occ_f(:)
         REAL(DP), INTENT(in) :: mat_z(:,:,:)
      END SUBROUTINE writefile_x
   END INTERFACE
 
!
   INTERFACE runcp_uspp
      SUBROUTINE runcp_uspp_x &
         ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec_bgrp, c0_bgrp, c0_d, cm_bgrp, cm_d, fromscra, restart, compute_only_gradient )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         integer, intent(in) :: nfi
         real(DP) :: fccc, ccc
         real(DP) :: ema0bg(:), dt2bye
         real(DP) :: rhos(:,:)
         real(DP) :: bec_bgrp(:,:)
         complex(DP) :: c0_bgrp(:,:), cm_bgrp(:,:)
         complex(DP) DEVICEATTR :: c0_d(:,:), cm_d(:,:)
         logical, optional, intent(in) :: fromscra
         logical, optional, intent(in) :: restart, compute_only_gradient
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


   INTERFACE eigs
      SUBROUTINE cp_eigs_x( nfi, lambdap, lambda, idesc )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         INTEGER :: nfi
         REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
         INTEGER, INTENT(IN) :: idesc( :, : )
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
      SUBROUTINE ortho_x &
         ( betae, cp_bgrp, phi_bgrp, x0, idesc, diff, iter, ccc, bephi, becp_bgrp )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         INTEGER,  INTENT(IN) :: idesc( :, : )
         COMPLEX(DP) :: betae( :, : )
         COMPLEX(DP) :: cp_bgrp( :, : ), phi_bgrp( :, : )
         REAL(DP)    :: x0( :, :, : ), diff, ccc
         INTEGER     :: iter
         REAL(DP)    :: bephi(:,:)
         REAL(DP)    :: becp_bgrp(:,:)
#if defined (__CUDA)
         ATTRIBUTES( DEVICE ) :: becp_bgrp, bephi, cp_bgrp, phi_bgrp, betae
#endif
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


   INTERFACE wave_rand_init
      SUBROUTINE wave_rand_init_x( cm, global )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(OUT) :: cm(:,:)
         LOGICAL, OPTIONAL, INTENT(IN) :: global
      END SUBROUTINE
   END INTERFACE

   INTERFACE crot
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
         USE fft_base,         ONLY: dfftp
         IMPLICIT NONE
         COMPLEX(DP) :: ei1( -dfftp%nr1 : dfftp%nr1, nat )
         COMPLEX(DP) :: ei2( -dfftp%nr2 : dfftp%nr2, nat )
         COMPLEX(DP) :: ei3( -dfftp%nr3 : dfftp%nr3, nat )
         INTEGER      :: mill( :, : )
         INTEGER      :: ngm
         COMPLEX(DP), INTENT(OUT) :: sfac(:,:)
      END SUBROUTINE
   END INTERFACE

   
   INTERFACE printout_new
      SUBROUTINE printout_new_x &
         ( nfi, tfirst, tfilei, tprint, tps, h, stress, tau0, vels, &
           fion, ekinc, temphc, tempp, temps, etot, enthal, econs, econt, &
           vnhh, xnhh0, vnhp, xnhp0, vnhe, xnhe0, atot, ekin, epot, print_forces, print_stress,tstdout )
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
         REAL(DP), INTENT(IN) :: vnhh( 3, 3 ), xnhh0( 3, 3 ), vnhp( 1 ), xnhp0( 1 ), vnhe, xnhe0
         REAL(DP), INTENT(IN) :: atot! enthalpy of system for c.g. case
         REAL(DP), INTENT(IN) :: ekin
         REAL(DP), INTENT(IN) :: epot ! ( epseu + eht + exc )
         LOGICAL, INTENT(IN) :: print_forces, print_stress, tstdout
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
      SUBROUTINE force_loc_x( tscreen, rhoeg, fion, rhops, vps, ei1, ei2, ei3, &
                        sfac, omega, screen_coul )
         USE kinds,              ONLY: DP
         USE fft_base,           ONLY: dfftp
         USE ions_base,          ONLY: nat
         IMPLICIT NONE
         LOGICAL     :: tscreen
         REAL(DP)    :: fion(:,:) 
         REAL(DP)    :: rhops(:,:), vps(:,:)  
         COMPLEX(DP) :: rhoeg(:)
         COMPLEX(DP), INTENT(IN) :: sfac(:,:)
         COMPLEX(DP) :: ei1(-dfftp%nr1:dfftp%nr1,nat)
         COMPLEX(DP) :: ei2(-dfftp%nr2:dfftp%nr2,nat)
         COMPLEX(DP) :: ei3(-dfftp%nr3:dfftp%nr3,nat)
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
   
   INTERFACE set_eitot
      SUBROUTINE set_eitot_x( eitot )
         USE kinds, ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: eitot(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE set_evtot
      SUBROUTINE set_evtot_x( c0, ctot, lambda, idesc, iupdwn_tot, nupdwn_tot )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN)  :: c0(:,:)
         COMPLEX(DP), INTENT(OUT) :: ctot(:,:)
         REAL(DP),    INTENT(IN)  :: lambda(:,:,:)
         INTEGER, INTENT(IN) :: idesc(:,:)
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



   INTERFACE move_electrons
      SUBROUTINE move_electrons_x( &
         nfi, tprint, tfirst, tlast, b1, b2, b3, fion, enthal, enb, &
            &  enbi, fccc, ccc, dt2bye, stress,l_cprestart )
         USE kinds,         ONLY: DP
         IMPLICIT NONE
         INTEGER,  INTENT(IN)    :: nfi
         LOGICAL,  INTENT(IN)    :: tprint, tfirst, tlast
         REAL(DP), INTENT(IN)    :: b1(3), b2(3), b3(3)
         REAL(DP)                :: fion(:,:)
         REAL(DP), INTENT(IN)    :: dt2bye
         REAL(DP)                :: fccc, ccc
         REAL(DP)                :: enb, enbi
         REAL(DP)                :: enthal
         REAL(DP)                :: ei_unp
         REAL(DP)                :: stress(3,3)
         LOGICAL, INTENT(in)     :: l_cprestart
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
      SUBROUTINE nlfh_x( stress, bec, dbec, lambda, idesc )
         USE kinds, ONLY : DP
         IMPLICIT NONE
         REAL(DP), INTENT(INOUT) :: stress(3,3) 
         REAL(DP), INTENT(IN)    :: bec( :, : ), dbec( :, :, :, : )
         REAL(DP), INTENT(IN)    :: lambda( :, :, : )
         INTEGER, INTENT(IN) :: idesc(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE nlfl_bgrp
      SUBROUTINE nlfl_bgrp_x( bec_bgrp, becdr_bgrp, lambda, idesc, fion )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         REAL(DP) :: bec_bgrp(:,:), becdr_bgrp(:,:,:)
         REAL(DP), INTENT(IN) :: lambda(:,:,:)
         INTEGER, INTENT(IN) :: idesc(:,:)
         REAL(DP), INTENT(INOUT) :: fion(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE protate
      SUBROUTINE protate_x ( c0, bec, c0rot, becrot, ngwl, nss, noff, lambda, nrl, &
                           ityp, nat, ofsbeta, nh, np_rot, me_rot, comm_rot  )
         USE kinds,            ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngwl, nss, nrl, noff
         INTEGER, INTENT(IN) :: ityp(:), nat, ofsbeta(:), nh(:)
         INTEGER, INTENT(IN) :: np_rot, me_rot, comm_rot  
         COMPLEX(DP), INTENT(IN) :: c0(:,:)
         COMPLEX(DP), INTENT(OUT) :: c0rot(:,:)
         REAL(DP), INTENT(IN) :: lambda(:,:)
         REAL(DP), INTENT(IN) :: bec(:,:)
         REAL(DP), INTENT(OUT) :: becrot(:,:)
      END SUBROUTINE
   END INTERFACE

   INTERFACE c_bgrp_expand
    SUBROUTINE c_bgrp_expand_x( c_bgrp )
      USE kinds,              ONLY: DP
      IMPLICIT NONE
      COMPLEX(DP) :: c_bgrp(:,:)
    END SUBROUTINE c_bgrp_expand_x
#if defined (__CUDA)
    SUBROUTINE c_bgrp_expand_gpu_x( c_bgrp )
      USE kinds,              ONLY: DP
      IMPLICIT NONE
      COMPLEX(DP), DEVICE :: c_bgrp(:,:)
    END SUBROUTINE c_bgrp_expand_gpu_x
#endif
   END INTERFACE
   INTERFACE c_bgrp_pack
    SUBROUTINE c_bgrp_pack_x( c_bgrp )
      USE kinds,              ONLY: DP
      IMPLICIT NONE
      COMPLEX(DP) :: c_bgrp(:,:)
    END SUBROUTINE c_bgrp_pack_x
#if defined (__CUDA)
    SUBROUTINE c_bgrp_pack_gpu_x( c_bgrp )
      USE kinds,              ONLY: DP
      IMPLICIT NONE
      COMPLEX(DP), DEVICE :: c_bgrp(:,:)
    END SUBROUTINE c_bgrp_pack_gpu_x
#endif
   END INTERFACE

   INTERFACE vofrho
      SUBROUTINE vofrho_x( nfi, rhor, drhor, rhog, drhog, rhos, rhoc, tfirst, tlast,           &
     &     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         LOGICAL     :: tlast, tfirst
         INTEGER     :: nfi
         REAL(DP)    :: rhor(:,:), drhor(:,:,:,:), rhos(:,:), fion(:,:)
         REAL(DP)    :: rhoc(:), tau0(:,:)
         ! COMPLEX(DP) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat), ei3(-nr3:nr3,nat)
         COMPLEX(DP) :: ei1(:,:), ei2(:,:), ei3(:,:)
         COMPLEX(DP) :: eigrb(:,:)
         COMPLEX(DP) :: rhog(:,:), drhog(:,:,:,:)
         COMPLEX(DP) :: sfac(:,:)
         INTEGER     :: irb(:,:)
      END SUBROUTINE vofrho_x
   END INTERFACE

   INTERFACE enkin
      FUNCTION enkin_x( c, f, n )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         INTEGER,     INTENT(IN) :: n
         COMPLEX(DP), INTENT(IN) :: c( :, : )
         REAL(DP),    INTENT(IN) :: f( : )
         REAL(DP) :: enkin_x
      END FUNCTION enkin_x
#if defined (__CUDA)
      FUNCTION enkin_gpu_x( c, f, n )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         INTEGER,     INTENT(IN) :: n
         COMPLEX(DP), DEVICE, INTENT(IN) :: c( :, : )
         REAL(DP),    DEVICE, INTENT(IN) :: f( : )
         REAL(DP) :: enkin_gpu_x
      END FUNCTION enkin_gpu_x
#endif
   END INTERFACE 

   INTERFACE newinit
      SUBROUTINE newinit_x( h, iverbosity )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         REAL(DP),    INTENT(IN) :: h( 3, 3 )
         INTEGER,     INTENT(IN) :: iverbosity
      END SUBROUTINE newinit_x
   END INTERFACE

   INTERFACE prefor
      SUBROUTINE prefor_x( eigr, betae )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: eigr( :, : ) 
         COMPLEX(DP), INTENT(OUT) :: betae( :, : )
      END SUBROUTINE prefor_x
   END INTERFACE

   INTERFACE denlcc
      SUBROUTINE denlcc_x( nnr, nspin, vxcr, sfac, drhocg, dcc )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         INTEGER,     INTENT(IN) :: nnr, nspin
         REAL(DP),    INTENT(IN) :: vxcr( :, : )
         COMPLEX(DP), INTENT(IN) :: sfac( :, : ) 
         REAL(DP),    INTENT(IN) :: drhocg( :, : ) 
         REAL(DP), INTENT(OUT) ::  dcc( :, : )
      END SUBROUTINE denlcc_x
   END INTERFACE

   INTERFACE dotcsc
      SUBROUTINE dotcsc_x( betae, cp, ngw, n )
         USE kinds, ONLY: dp
         IMPLICIT NONE
         INTEGER,     INTENT(IN) :: ngw, n
         COMPLEX(DP), INTENT(IN) :: cp(:,:)
         COMPLEX(DP), INTENT(INOUT) :: betae(:,:)
      END SUBROUTINE dotcsc_x
   END INTERFACE

   INTERFACE nlsm1
      SUBROUTINE nlsm1_x ( n, betae, c, becp, pptype_ )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)  :: n
         COMPLEX(DP), INTENT(IN)  :: c( :, : )
         COMPLEX(DP), INTENT(INOUT)  :: betae( :, : )
         REAL(DP),    INTENT(OUT) :: becp( :, : )
         INTEGER,     INTENT(IN), OPTIONAL  :: pptype_
      END SUBROUTINE nlsm1_x 
#if defined (__CUDA)
      SUBROUTINE nlsm1_gpu_x ( n, betae, c, becp, pptype_ )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)  :: n
         COMPLEX(DP), INTENT(IN), DEVICE  :: c( :, : )
         COMPLEX(DP), INTENT(INOUT), DEVICE  :: betae( :, : )
         REAL(DP),    INTENT(OUT), DEVICE :: becp( :, : )
         INTEGER,     INTENT(IN), OPTIONAL  :: pptype_
      END SUBROUTINE nlsm1_gpu_x
#endif
   END INTERFACE

   INTERFACE nlsm2_bgrp
      SUBROUTINE  nlsm2_bgrp_x( ngw, nkb, betae, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)  :: ngw, nkb, nbspx_bgrp, nbsp_bgrp
         COMPLEX(DP), INTENT(IN)  :: betae( :, : ), c_bgrp( :, : )
         REAL(DP),    INTENT(OUT) :: becdr_bgrp( :, :, : )
      END SUBROUTINE nlsm2_bgrp_x
#if defined (__CUDA)
      SUBROUTINE  nlsm2_bgrp_gpu_x( ngw, nkb, betae, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)  :: ngw, nkb, nbspx_bgrp, nbsp_bgrp
         COMPLEX(DP), INTENT(IN), DEVICE :: betae( :, : )
         COMPLEX(DP), INTENT(IN), DEVICE :: c_bgrp( :, : )
         REAL(DP),    INTENT(OUT) PINMEM :: becdr_bgrp( :, :, : )
      END SUBROUTINE nlsm2_bgrp_gpu_x
#endif
   END INTERFACE

   INTERFACE ennl
      SUBROUTINE ennl_x( ennl_val, rhovan, bec_bgrp )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: ennl_val
         REAL(DP), INTENT(OUT) :: rhovan( :, :, : )
         REAL(DP), INTENT(IN)  :: bec_bgrp( :, : )
      END SUBROUTINE ennl_x
   END INTERFACE

   INTERFACE calrhovan
      SUBROUTINE calrhovan_x( rhovan, bec, iwf )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP), INTENT(OUT) :: rhovan( :, :, : )
         REAL(DP), INTENT(IN)  :: bec( :, : )
         INTEGER,  INTENT(IN)  :: iwf
      END SUBROUTINE calrhovan_x
   END INTERFACE

   INTERFACE calbec
      SUBROUTINE calbec_x( n, betae, c, bec, pptype_ )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)    :: n
         REAL(DP),    INTENT(OUT)   :: bec( :, : )
         COMPLEX(DP), INTENT(IN)    :: c( :, : )
         COMPLEX(DP), INTENT(INOUT) :: betae( :, : )
         INTEGER,     INTENT(IN), OPTIONAL  :: pptype_
      END SUBROUTINE calbec_x
#if defined (__CUDA)
      SUBROUTINE calbec_gpu_x( n, betae, c, bec, pptype_ )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)    :: n
         REAL(DP),    INTENT(OUT), DEVICE   :: bec( :, : )
         COMPLEX(DP), INTENT(IN), DEVICE    :: c( :, : )
         COMPLEX(DP), INTENT(INOUT), DEVICE :: betae( :, : )
         INTEGER,     INTENT(IN), OPTIONAL  :: pptype_
      END SUBROUTINE calbec_gpu_x
#endif
   END INTERFACE

   INTERFACE caldbec_bgrp
      SUBROUTINE caldbec_bgrp_x( eigr, c_bgrp, dbec, idesc )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN)  ::  c_bgrp( :, : ), eigr( :, : )
         REAL(DP),    INTENT(OUT) ::  dbec( :, :, :, : )
         INTEGER, INTENT(IN) :: idesc( :, : )
      END SUBROUTINE caldbec_bgrp_x
#if defined (__CUDA)
      SUBROUTINE caldbec_bgrp_gpu_x( eigr, c_bgrp, dbec, idesc )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN), DEVICE ::  eigr( :, : )
         COMPLEX(DP), INTENT(IN), DEVICE ::  c_bgrp( :, : )
         REAL(DP),    INTENT(OUT) PINMEM ::  dbec( :, :, :, : )
         INTEGER, INTENT(IN) :: idesc( :, : )
      END SUBROUTINE caldbec_bgrp_gpu_x
#endif
   END INTERFACE

   INTERFACE dennl
      SUBROUTINE dennl_x( bec_bgrp, dbec, drhovan, denl, idesc )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         REAL(DP),    INTENT(IN)  ::  dbec( :, :, :, : )
         REAL(DP),    INTENT(IN)  ::  bec_bgrp( :, : )
         REAL(DP),    INTENT(OUT) ::  drhovan( :, :, :, :, : )
         REAL(DP),    INTENT(OUT) ::  denl( 3, 3 )
         INTEGER, INTENT(IN) :: idesc( :, : )
      END SUBROUTINE dennl_x
   END INTERFACE

   INTERFACE nlfq_bgrp
      SUBROUTINE nlfq_bgrp_x( c_bgrp, betae, bec_bgrp, becdr_bgrp, fion )
         USE kinds,              ONLY: DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) DEVICEATTR :: c_bgrp( :, : )
         COMPLEX(DP), INTENT(IN) DEVICEATTR :: betae( :, : )
         REAL(DP),    INTENT(IN)  ::  bec_bgrp( :, : )
         REAL(DP),    INTENT(OUT) ::  becdr_bgrp( :, :, : )
         REAL(DP),    INTENT(OUT) ::  fion( :, : )
      END SUBROUTINE nlfq_bgrp_x
   END INTERFACE

   INTERFACE collect_bec
      SUBROUTINE collect_bec_x( bec_repl, bec_dist, idesc, nspin )
         USE kinds,       ONLY : DP
         REAL(DP), INTENT(OUT) :: bec_repl(:,:)
         REAL(DP), INTENT(IN)  :: bec_dist(:,:)
         INTEGER, INTENT(IN)   :: idesc(:,:)
         INTEGER,  INTENT(IN)  :: nspin
      END SUBROUTINE collect_bec_x
   END INTERFACE

   INTERFACE nlsm1us
      SUBROUTINE nlsm1us_x ( n, betae, c, becp )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         INTEGER,     INTENT(IN)  :: n
         COMPLEX(DP) DEVICEATTR , INTENT(INOUT)  :: betae( :, : )
         COMPLEX(DP) DEVICEATTR , INTENT(IN)  :: c( :, : )
         REAL(DP) DEVICEATTR ,    INTENT(OUT) :: becp( :, : )
      END SUBROUTINE nlsm1us_x
   END INTERFACE

   INTERFACE dbeta_eigr
      SUBROUTINE dbeta_eigr_x( dbeigr, eigr )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN)  :: eigr( :, : )
         COMPLEX(DP), INTENT(OUT) :: dbeigr( :, :, :, :)
      END SUBROUTINE dbeta_eigr_x
#if defined (__CUDA)
      SUBROUTINE dbeta_eigr_gpu_x( dbeigr, eigr )
         USE kinds,      ONLY : DP
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN),  DEVICE :: eigr( :, : )
         COMPLEX(DP), INTENT(OUT), DEVICE :: dbeigr( :, :, :, :)
      END SUBROUTINE dbeta_eigr_gpu_x
#endif
   END INTERFACE

!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   END MODULE
!=----------------------------------------------------------------------------=!
