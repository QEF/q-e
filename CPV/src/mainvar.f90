!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
#define PINMEM 
#else
#define PINMEM
#endif
!
!----------------------------------------------------------------------------
MODULE cp_main_variables
  !----------------------------------------------------------------------------
  !! Main global variables for CP and allocation routines.
  !
  USE kinds,             ONLY : DP
  USE xc_lib,            ONLY : xclib_dft_is
  USE metagga_cp,        ONLY : kedtaur, kedtaus, kedtaug
  USE cell_base,         ONLY : boxdimensions
  USE wave_types,        ONLY : wave_descriptor, wave_descriptor_init
  USE energies,          ONLY : dft_energy_type
  USE pres_ai_mod,       ONLY : abivol, abisur, jellium, t_gauss, rho_gaus, &
                                v_vol, posv, f_vol
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... structure factors e^{-ig*R}
  !
  ! ...  G = reciprocal lattice vectors
  ! ...  R_I = ionic positions
  !
  COMPLEX(DP), ALLOCATABLE  PINMEM :: eigr(:,:)
  !! it is \(e^{i G \cdot R_I}\), where \(G\) reciprocal 
  !! lattice vectors and \(R_I\) ionic positions.
#if defined (__CUDA)
  COMPLEX(DP), ALLOCATABLE, DEVICE :: eigr_d(:,:)
  !! GPU double of \(\text{eigr}\)
#endif
  !
  ! ... structure factors (summed over atoms of the same kind)
  !
  ! S( s, G ) = sum_(I in s) exp( i G dot R_(s,I) )
  ! s       = index of the atomic specie
  ! R_(s,I) = position of the I-th atom of the "s" specie
  !
  COMPLEX(DP), ALLOCATABLE:: sfac(:,:)
  !! structure factor
  !
  ! ... indexes, positions, and structure factors for the box grid
  !
  REAL(DP), ALLOCATABLE :: taub(:,:)
  COMPLEX(DP), ALLOCATABLE PINMEM :: eigrb(:,:)
  INTEGER,     ALLOCATABLE :: irb(:,:)
  INTEGER,     ALLOCATABLE :: iabox(:)
  INTEGER :: nabox
  ! 
  ! ... nonlocal projectors:
  ! ...    bec   = scalar product of projectors and wave functions
  ! ...    betae = nonlocal projectors in g space = beta x e^(-ig.R) 
  ! ...    becdr = <betae|g|psi> used in force calculation
  ! ...    rhovan= \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  ! ...    deeq  = \int V_eff(r) q_lm(r) dr
  !
  REAL(DP), ALLOCATABLE :: bephi(:,:)
  !! distributed (orhto group)
  REAL(DP), ALLOCATABLE :: becp_bgrp(:,:)
  !! distributed becp (band group)
  REAL(DP), ALLOCATABLE PINMEM :: bec_bgrp(:,:)
  !! scalar product of projectors and wave functions.
  !! Distributed (band group)
  REAL(DP), ALLOCATABLE :: bec_d(:,:)
  !! GPU double of bec (band group)
  REAL(DP), ALLOCATABLE PINMEM :: becdr_bgrp(:,:,:)
  !! \(\langle\text{betae}|g|\text{psi}\rangle\) used in force calculation,
  !! with \(\text{betae}\) nonlocal projectors in g space.  
  !! Distributed (band group)
  REAL(DP), ALLOCATABLE PINMEM :: dbec(:,:,:,:)
  !! derivative of bec distributed (ortho group).
#if defined (__CUDA)
  REAL(DP), ALLOCATABLE, DEVICE :: dbec_d(:,:,:,:)
  !! GPU double of \(\text{dbec}\) 
  ATTRIBUTES( DEVICE ) :: becp_bgrp, bephi, bec_d
#endif
  !
  REAL(DP), ALLOCATABLE :: ema0bg(:)
  !! mass preconditioning
  !
  REAL(DP), ALLOCATABLE :: lambda(:,:,:)
  !! constraints (lambda at t)
  REAL(DP), ALLOCATABLE :: lambdam(:,:,:)
  !! constraints (lambda at t-dt)
  REAL(DP), ALLOCATABLE :: lambdap(:,:,:)
  !! constraints (lambda at t+dt)
  !
  INTEGER, ALLOCATABLE :: idesc(:,:)
  !! laxlib descriptor of the lambda distribution
  !
  INTEGER, PARAMETER :: nacx = 10
  !! max number of averaged quantities saved to the restart
  REAL(DP) :: acc(nacx)
  REAL(DP) :: acc_this_run(nacx)
  !
  ! cell geometry
  !
  TYPE (boxdimensions) :: htm, ht0, htp  ! cell metrics
  !
  ! charge densities and potentials
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  !! charge density in g space
  REAL(DP),    ALLOCATABLE :: rhor(:,:)
  !! charge density in r space (dense grid)
  REAL(DP),    ALLOCATABLE :: rhos(:,:)
  !! charge density in r space (smooth grid)
  REAL(DP),    ALLOCATABLE :: vpot(:,:)
  !! potential in r space (dense grid)
  !
  ! derivative wrt cell
  !
  COMPLEX(DP), ALLOCATABLE :: drhog(:,:,:,:)
  !! derivative of \(\text{rhog}\) wrt cell
  REAL(DP),    ALLOCATABLE :: drhor(:,:,:,:)
  !! derivative of \(\text{rhor}\) wrt cell
  !
  TYPE (wave_descriptor) :: wfill
  !! wave function descriptor for filled
  !
  TYPE(dft_energy_type) :: edft
  !
  INTEGER :: nfi
  !! counter on the electronic iterations
  INTEGER :: nprint_nfi=-1
  !! counter indicating the last time data have been printed on file
  !! ( prefix.pos, ... ), it is used to avoid printing stuff two times.
  INTEGER :: nfi_run=0
  !! counter on the electronic iterations, for the present run
  INTEGER :: iprint_stdout=1
  !! define how often CP writes verbose information to stdout
  !
  ! working buffers
  !
#if defined (__CUDA)
  REAL(DP), ALLOCATABLE, DEVICE :: nlsm1_wrk_d(:,:)
  !! working buffer for nlsm1 function
  COMPLEX(DP), ALLOCATABLE, DEVICE :: caldbec_wrk_d(:,:,:,:)
  !! working buffer for caldbec_bgrp function
  REAL(DP),    ALLOCATABLE, DEVICE :: caldbec_dwrk_d(:,:)
  !! working buffer for caldbec_bgrp function
#endif
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_mainvar( ngw, ngw_g, ngb, ngs, ng, nr1, nr2, nr3, &
                                 nr1x, nr2x, npl, nnr, nrxxs, nat, nax,  &
                                 nsp, nspin, n, nx, nupdwn, nhsa, &
                                 gstart, nudx, tpre, nbspx_bgrp )
      !------------------------------------------------------------------------
      !! Allocate CP main global variables.
      !
      USE mp_bands,    ONLY: intra_bgrp_comm, me_bgrp
      USE mp,          ONLY: mp_max, mp_min
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER,           INTENT(IN) :: ngw, ngw_g, ngb, ngs, ng, nr1,nr2,nr3, &
                                       nnr, nrxxs, nat, nax, nsp, nspin, &
                                       n, nx, nhsa, nr1x, nr2x, npl
      INTEGER,           INTENT(IN) :: nupdwn(:)
      INTEGER,           INTENT(IN) :: gstart, nudx
      LOGICAL,           INTENT(IN) :: tpre
      INTEGER,           INTENT(IN) :: nbspx_bgrp
      !
      INTEGER  :: iss, ierr, nlam, nrcx
      LOGICAL  :: gzero
      INTEGER  :: np_ortho(2), me_ortho(2), ortho_comm, ortho_comm_id, ortho_cntx
      !
      CALL laxlib_getval( np_ortho = np_ortho, me_ortho = me_ortho, ortho_comm = ortho_comm, &
        ortho_comm_id = ortho_comm_id, ortho_cntx = ortho_cntx )
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      ALLOCATE( eigr( ngw, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate eigr ', ierr )
#if defined (__CUDA)
      ALLOCATE( eigr_d( ngw, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate eigr_d ', ierr )
#endif
      ALLOCATE( sfac( ngs, nsp ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate sfac ', ierr )
      ALLOCATE( eigrb( ngb, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate eigrb ', ierr )
      ALLOCATE( irb( 3, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate irb ', ierr )
      ALLOCATE( iabox( nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate iabox ', ierr )
      nabox = 0
      !
      IF ( xclib_dft_is('meta') ) THEN
         !
         ! ... METAGGA
         !
         ALLOCATE( kedtaur( nnr,   nspin ) )
         ALLOCATE( kedtaus( nrxxs, nspin ) )
         ALLOCATE( kedtaug( ng,    nspin ) )
         !
      ELSE
         !
         ! ... dummy allocation required because this array appears in the
         ! ... list of arguments of some routines
         !
         ALLOCATE( kedtaur( 1, nspin ) )
         ALLOCATE( kedtaus( 1, nspin ) )
         ALLOCATE( kedtaug( 1, nspin ) )
         !
      END IF
      !
      ALLOCATE( ema0bg( ngw ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate ema0bg ', ierr )
      !
      ALLOCATE( rhor( nnr, nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate rhor ', ierr )
      ALLOCATE( vpot( nnr, nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate vpot ', ierr )
      ALLOCATE( rhos( nrxxs, nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate rhos ', ierr )
      ALLOCATE( rhog( ng,    nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate rhog ', ierr )
      IF ( tpre ) THEN
            ALLOCATE( drhog( ng,  nspin, 3, 3 ), STAT=ierr )
            IF( ierr /= 0 ) &
               CALL errore( ' allocate_mainvar ', ' unable to allocate drhog ', ierr )
            ALLOCATE( drhor( nnr, nspin, 3, 3 ), STAT=ierr )
            IF( ierr /= 0 ) &
               CALL errore( ' allocate_mainvar ', ' unable to allocate drhor ', ierr )
      ELSE
            ALLOCATE( drhog( 1, 1, 1, 1 ) )
            ALLOCATE( drhor( 1, 1, 1, 1 ) )
      END IF
!==========================================================================
      !
      !  Compute local dimensions for lambda matrixes
      !

      ALLOCATE( idesc( LAX_DESC_SIZE, nspin ) )
      !
      DO iss = 1, nspin
         CALL laxlib_init_desc( idesc( :, iss ), nupdwn( iss ), nudx, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
      END DO
      !
      nrcx = MAXVAL( idesc( LAX_DESC_NRCX, : ) )
      !
      nlam = 1
      IF( SIZE( idesc, 2 ) < 2 ) THEN
         IF( idesc( LAX_DESC_ACTIVE_NODE, 1 ) > 0 ) &
            nlam = idesc(LAX_DESC_NRCX,1)
      ELSE
         IF( ( idesc( LAX_DESC_ACTIVE_NODE, 1) > 0 ) .OR. ( idesc( LAX_DESC_ACTIVE_NODE, 2 ) > 0 ) ) &
            nlam = MAX( idesc(LAX_DESC_NRCX,1), idesc(LAX_DESC_NRCX,2) )
      END IF

      !
      !
      !  ... End with lambda dimensions
      !
      !
      if ( abivol.or.abisur ) then
         !
         allocate(rho_gaus(nnr))
         allocate(v_vol(nnr))
         if (jellium.or.t_gauss) allocate(posv(3,nr1*nr2*nr3))
         if (t_gauss) allocate(f_vol(3,nax,nsp))
         !
      end if
      !
      ALLOCATE( lambda(  nlam, nlam, nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate lambda ', ierr )
      ALLOCATE( lambdam( nlam, nlam, nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate lambdam ', ierr )
      ALLOCATE( lambdap( nlam, nlam, nspin ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate lambdap ', ierr )
      !
      ! becdr, distributed over row processors of the ortho group
      !
      ALLOCATE( becdr_bgrp( nhsa, nbspx_bgrp, 3 ), STAT=ierr )  
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate becdr_bgrp ', ierr )
      ALLOCATE( bec_bgrp( nhsa, nbspx_bgrp ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate bec_bgrp ', ierr )
#if defined (__CUDA)
      ALLOCATE( bec_d( nhsa, nbspx_bgrp ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate bec_d ', ierr )
#endif
      ALLOCATE( bephi( nhsa, nspin*nrcx ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate becphi ', ierr )
      ALLOCATE( becp_bgrp( nhsa, nbspx_bgrp ), STAT=ierr )  
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate becp_bgrp ', ierr )
      !
      IF ( tpre ) THEN
        ALLOCATE( dbec( nhsa, 2*nrcx, 3, 3 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' allocate_mainvar ', ' unable to allocate dbec ', ierr )
      ELSE
        ALLOCATE( dbec( 1, 1, 1, 1 ) )
      END IF
#if defined (__CUDA)
      IF ( tpre ) THEN
        ALLOCATE( dbec_d( nhsa, 2*nrcx, 3, 3 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' allocate_mainvar ', ' unable to allocate dbec_d ', ierr )
      ELSE
        ALLOCATE( dbec_d( 1, 1, 1, 1 ) )
      END IF
#endif

      gzero =  (gstart == 2)
      !
      CALL wave_descriptor_init( wfill, ngw, ngw_g, nupdwn,  nupdwn, &
            1, 1, nspin, 'gamma', gzero )
      !
      ! allocating working buffers
      !
#if defined (__CUDA)
      ALLOCATE( nlsm1_wrk_d( nhsa, nbspx_bgrp ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate nlsm1_wrk_d ', ierr )
      IF ( tpre ) THEN
        ALLOCATE( caldbec_wrk_d( ngw, nhsa, 3, 3 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' allocate_mainvar ', ' unable to allocate caldbec_wrk_d ', ierr )
        ALLOCATE( caldbec_dwrk_d( nhsa, nbspx_bgrp ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' allocate_mainvar ', ' unable to allocate caldbec_wrk_d ', ierr )
      ELSE
        ALLOCATE( caldbec_wrk_d( 1, 1, 1, 1 ) )
        ALLOCATE( caldbec_dwrk_d( 1, 1 ) )
      END IF
#endif
      !
      RETURN
      !
    END SUBROUTINE allocate_mainvar
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_mainvar()
      !------------------------------------------------------------------------
      !! Deallocate main CP global variables.
      !
      IF( ALLOCATED( eigr ) )    DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) )    DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) )   DEALLOCATE( eigrb )
      IF( ALLOCATED( irb ) )     DEALLOCATE( irb )
      IF( ALLOCATED( iabox ) )     DEALLOCATE( iabox )
      IF( ALLOCATED( rhor ) )    DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) )    DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) )    DEALLOCATE( rhog )
      IF( ALLOCATED( drhog ) )   DEALLOCATE( drhog )
      IF( ALLOCATED( drhor ) )   DEALLOCATE( drhor )
      IF( ALLOCATED( bec_bgrp ) )     DEALLOCATE( bec_bgrp )
      IF( ALLOCATED( bec_d ) )     DEALLOCATE( bec_d )
      IF( ALLOCATED( becdr_bgrp ) )   DEALLOCATE( becdr_bgrp )
      IF( ALLOCATED( bephi ) )   DEALLOCATE( bephi )
      IF( ALLOCATED( becp_bgrp ) )    DEALLOCATE( becp_bgrp )
      IF( ALLOCATED( dbec ) )    DEALLOCATE( dbec )
      IF( ALLOCATED( ema0bg ) )  DEALLOCATE( ema0bg )
      IF( ALLOCATED( lambda ) )  DEALLOCATE( lambda )
      IF( ALLOCATED( lambdam ) ) DEALLOCATE( lambdam )
      IF( ALLOCATED( lambdap ) ) DEALLOCATE( lambdap )
      IF( ALLOCATED( kedtaur ) ) DEALLOCATE( kedtaur )
      IF( ALLOCATED( kedtaus ) ) DEALLOCATE( kedtaus )
      IF( ALLOCATED( kedtaug ) ) DEALLOCATE( kedtaug )
      IF( ALLOCATED( vpot ) )    DEALLOCATE( vpot )
      IF( ALLOCATED( taub ) )    DEALLOCATE( taub )
      IF( ALLOCATED( idesc ) )  DEALLOCATE( idesc )
      !
#if defined (__CUDA)
      IF( ALLOCATED( eigr_d ) )    DEALLOCATE( eigr_d )
      IF( ALLOCATED( dbec_d ) )    DEALLOCATE( dbec_d )
      IF( ALLOCATED( nlsm1_wrk_d ) )  DEALLOCATE( nlsm1_wrk_d )
      IF( ALLOCATED( caldbec_wrk_d ) )  DEALLOCATE( caldbec_wrk_d )
      IF( ALLOCATED( caldbec_dwrk_d ) )  DEALLOCATE( caldbec_dwrk_d )
#endif
      !
      RETURN
      !
    END SUBROUTINE deallocate_mainvar
    !
END MODULE cp_main_variables
