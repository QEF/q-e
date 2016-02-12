!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE cp_main_variables
  !----------------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE funct,             ONLY : dft_is_meta
  USE metagga,           ONLY : kedtaur, kedtaus, kedtaug
  USE cell_base,         ONLY : boxdimensions
  USE wave_types,        ONLY : wave_descriptor, wave_descriptor_init
  USE energies,          ONLY : dft_energy_type
  USE pres_ai_mod,       ONLY : abivol, abisur, jellium, t_gauss, rho_gaus, &
                                v_vol, posv, f_vol
  USE descriptors,       ONLY : la_descriptor
  !
  IMPLICIT NONE
  SAVE
  !
  ! ... structure factors e^{-ig*R}
  !
  ! ...  G = reciprocal lattice vectors
  ! ...  R_I = ionic positions
  !
  COMPLEX(DP), ALLOCATABLE :: eigr(:,:)        ! exp (i G   dot R_I)
  !
  ! ... structure factors (summed over atoms of the same kind)
  !
  ! S( s, G ) = sum_(I in s) exp( i G dot R_(s,I) )
  ! s       = index of the atomic specie
  ! R_(s,I) = position of the I-th atom of the "s" specie
  !
  COMPLEX(DP), ALLOCATABLE:: sfac(:,:)
  !
  ! ... indexes, positions, and structure factors for the box grid
  !
  REAL(DP), ALLOCATABLE :: taub(:,:)
  COMPLEX(DP), ALLOCATABLE :: eigrb(:,:)
  INTEGER,     ALLOCATABLE :: irb(:,:)
  ! 
  ! ... nonlocal projectors:
  ! ...    bec   = scalar product of projectors and wave functions
  ! ...    betae = nonlocal projectors in g space = beta x e^(-ig.R) 
  ! ...    becdr = <betae|g|psi> used in force calculation
  ! ...    rhovan= \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  ! ...    deeq  = \int V_eff(r) q_lm(r) dr
  !
  REAL(DP), ALLOCATABLE :: bephi(:,:)      ! distributed (orhto group)
  REAL(DP), ALLOCATABLE :: becp_bgrp(:,:)  ! distributed becp (band group)
  REAL(DP), ALLOCATABLE :: bec_bgrp(:,:)  ! distributed bec (band group)
  REAL(DP), ALLOCATABLE :: becdr_bgrp(:,:,:)  ! distributed becdr (band group)
  REAL(DP), ALLOCATABLE :: dbec(:,:,:,:)    ! derivative of bec distributed(ortho group) 
  !
  ! ... mass preconditioning
  !
  REAL(DP), ALLOCATABLE :: ema0bg(:)
  !
  ! ... constraints (lambda at t, lambdam at t-dt, lambdap at t+dt)
  !
  REAL(DP), ALLOCATABLE :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
  !
  TYPE(la_descriptor), ALLOCATABLE :: descla(:) ! descriptor of the lambda distribution
                                       ! see descriptors_module
  !
  INTEGER, PARAMETER :: nacx = 10      ! max number of averaged
                                       ! quantities saved to the restart
  REAL(DP) :: acc(nacx)
  REAL(DP) :: acc_this_run(nacx)
  !
  ! cell geometry
  !
  TYPE (boxdimensions) :: htm, ht0, htp  ! cell metrics
  !
  ! charge densities and potentials
  !
  ! rhog  = charge density in g space
  ! rhor  = charge density in r space (dense grid)
  ! rhos  = charge density in r space (smooth grid)
  ! vpot  = potential in r space (dense grid)
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  REAL(DP),    ALLOCATABLE :: rhor(:,:), rhos(:,:)
  REAL(DP),    ALLOCATABLE :: vpot(:,:)
  !
  ! derivative wrt cell
  !
  COMPLEX(DP), ALLOCATABLE :: drhog(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: drhor(:,:,:,:)

  TYPE (wave_descriptor) :: wfill     ! wave function descriptor for filled
  !
  TYPE(dft_energy_type) :: edft
  !
  INTEGER :: nfi             ! counter on the electronic iterations
  INTEGER :: nprint_nfi=-1   ! counter indicating the last time data have been
                             ! printed on file ( prefix.pos, ... ), it is used
                             ! to avoid printing stuff two times .
  INTEGER :: nfi_run=0       ! counter on the electronic iterations,
                             ! for the present run
  INTEGER :: iprint_stdout=1 ! define how often CP writes verbose information to stdout
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_mainvar( ngw, ngw_g, ngb, ngs, ng, nr1, nr2, nr3, &
                                 nr1x, nr2x, npl, nnr, nrxxs, nat, nax,  &
                                 nsp, nspin, n, nx, nupdwn, nhsa, &
                                 gstart, nudx, tpre, nbspx_bgrp )
      !------------------------------------------------------------------------
      !
      USE mp_global,   ONLY: np_ortho, me_ortho, intra_bgrp_comm, ortho_comm, &
                             me_bgrp, ortho_comm_id, ortho_cntx
      USE mp,          ONLY: mp_max, mp_min
      USE descriptors, ONLY: la_descriptor, descla_init
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
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      ALLOCATE( eigr( ngw, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate eigr ', ierr )
      ALLOCATE( sfac( ngs, nsp ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate sfac ', ierr )
      ALLOCATE( eigrb( ngb, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate eigrb ', ierr )
      ALLOCATE( irb( 3, nat ), STAT=ierr )
      IF( ierr /= 0 ) &
         CALL errore( ' allocate_mainvar ', ' unable to allocate irb ', ierr )
      !
      IF ( dft_is_meta() ) THEN
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

      ALLOCATE( descla( nspin ) )
      !
      DO iss = 1, nspin
         CALL descla_init( descla( iss ), nupdwn( iss ), nudx, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
      END DO
      !
      nrcx = MAXVAL( descla( : )%nrcx )
      !
      nlam = 1
      IF( SIZE( descla ) < 2 ) THEN
         IF( descla(1)%active_node > 0 ) &
            nlam = descla(1)%nrcx
      ELSE
         IF( ( descla(1)%active_node > 0 ) .OR. ( descla(2)%active_node > 0 ) ) &
            nlam = MAX( descla(1)%nrcx, descla(2)%nrcx )
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

      gzero =  (gstart == 2)
      !
      CALL wave_descriptor_init( wfill, ngw, ngw_g, nupdwn,  nupdwn, &
            1, 1, nspin, 'gamma', gzero )
      !
      RETURN
      !
    END SUBROUTINE allocate_mainvar
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_mainvar()
      !------------------------------------------------------------------------
      !
      IF( ALLOCATED( eigr ) )    DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) )    DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) )   DEALLOCATE( eigrb )
      IF( ALLOCATED( irb ) )     DEALLOCATE( irb )
      IF( ALLOCATED( rhor ) )    DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) )    DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) )    DEALLOCATE( rhog )
      IF( ALLOCATED( drhog ) )   DEALLOCATE( drhog )
      IF( ALLOCATED( drhor ) )   DEALLOCATE( drhor )
      IF( ALLOCATED( bec_bgrp ) )     DEALLOCATE( bec_bgrp )
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
      IF( ALLOCATED( descla ) )  DEALLOCATE( descla )
      !
      RETURN
      !
    END SUBROUTINE deallocate_mainvar
    !
END MODULE cp_main_variables
