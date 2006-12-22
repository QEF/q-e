!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE cp_main_variables
  !----------------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE parameters,        ONLY : natx, nsx, nacx
  USE control_flags,     ONLY : program_name
  USE funct,             ONLY : dft_is_meta
  USE metagga,           ONLY : kedtaur, kedtaus, kedtaug
  USE cell_base,         ONLY : boxdimensions
  USE wave_types,        ONLY : wave_descriptor, wave_descriptor_init
  USE energies,          ONLY : dft_energy_type
  USE pres_ai_mod,       ONLY : abivol, abisur, jellium, t_gauss, rho_gaus, &
                                v_vol, posv, f_vol
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
  COMPLEX(DP), ALLOCATABLE :: ei1(:,:)         ! exp (i G_x dot x_I)
  COMPLEX(DP), ALLOCATABLE :: ei2(:,:)         ! exp (i G_y dot y_I)
  COMPLEX(DP), ALLOCATABLE :: ei3(:,:)         ! exp (i G_z dot z_I)
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
  REAL(DP), ALLOCATABLE :: bec(:,:), becdr(:,:,:)
  REAL(DP), ALLOCATABLE :: bephi(:,:), becp(:,:)
  !
  ! ... mass preconditioning
  !
  REAL(DP), ALLOCATABLE :: ema0bg(:)
  !
  ! ... constraints (lambda at t, lambdam at t-dt, lambdap at t+dt)
  !
  REAL(DP), ALLOCATABLE :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
  !
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
  ! rhopr   since rhor is overwritten in vofrho,
  !         this array is used to save rhor for restart file
  ! vpot  = potential in r space (dense grid)
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  REAL(DP),    ALLOCATABLE :: rhor(:,:), rhos(:,:)
  REAL(DP),    ALLOCATABLE :: rhopr(:,:)  
  REAL(DP),    ALLOCATABLE :: vpot(:,:)
  !
  TYPE (wave_descriptor) :: wfill     ! wave function descriptor for filled
  !
  TYPE(dft_energy_type) :: edft
  !
  INTEGER :: nfi             ! counter on the electronic iterations
  INTEGER :: nprint_nfi=-1   ! counter indicating the last time data have been
                             ! printed on file ( prefix.pos, ... )
  INTEGER :: nfi_run=0       ! counter on the electronic iterations,
                             ! for the present run
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_mainvar( ngw, ngwt, ngb, ngs, ng, nr1, nr2, nr3, &
                                 nr1x, nr2x, npl, nnr, nnrsx, nat, nax,  &
                                 nsp, nspin, n, nx, n_emp, nupdwn, nhsa, &
                                 gzero, smd )
      !------------------------------------------------------------------------
      !
      INTEGER,           INTENT(IN) :: ngw, ngwt, ngb, ngs, ng, nr1, nr2, nr3, &
                                       nnr, nnrsx, nat, nax, nsp, nspin, &
                                       n, nx, n_emp, nhsa, nr1x, nr2x, npl
      INTEGER,           INTENT(IN) :: nupdwn(:)
      LOGICAL,           INTENT(IN) :: gzero
      LOGICAL, OPTIONAL, INTENT(IN) :: smd
      LOGICAL                       :: nosmd
      !
      INTEGER                       :: nudx
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      nosmd = .TRUE.
      nudx = MAXVAL( nupdwn( 1:nspin ) )
      !
      IF ( PRESENT( smd ) ) THEN
         !
         IF( smd ) nosmd = .FALSE.
         !
      END IF
      !
      ALLOCATE( eigr( ngw, nat ) )
      ALLOCATE( sfac( ngs, nsp ) )
      ALLOCATE( ei1( -nr1:nr1, nat ) )
      ALLOCATE( ei2( -nr2:nr2, nat ) )
      ALLOCATE( ei3( -nr3:nr3, nat ) )
      ALLOCATE( eigrb( ngb, nat ) )
      ALLOCATE( irb( 3, nat ) )
      !
      IF ( dft_is_meta() ) THEN
         !
         ! ... METAGGA
         !
         ALLOCATE( kedtaur( nnr,   nspin ) )
         ALLOCATE( kedtaus( nnrsx, nspin ) )
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
      ALLOCATE( ema0bg( ngw ) )
      !
      ALLOCATE( rhor( nnr, nspin ) )
      !
      IF( program_name == 'CP90' ) THEN
         !
         ALLOCATE( rhopr( nnr,   nspin ) )
         ALLOCATE( rhos( nnrsx, nspin ) )
         ALLOCATE( rhog( ng,    nspin ) )
         !
         if ( abivol.or.abisur ) then
            !
            allocate(rho_gaus(nnr))
            allocate(v_vol(nnr))
            if (jellium.or.t_gauss) allocate(posv(3,nr1*nr2*nr3))
            if (t_gauss) allocate(f_vol(3,natx,nsx))
            !
         end if
         !
         IF ( nosmd ) THEN
            !
            ALLOCATE( lambda(  nudx, nudx, nspin ) )
            ALLOCATE( lambdam( nudx, nudx, nspin ) )
            ALLOCATE( lambdap( nudx, nudx, nspin ) )
            !
         END IF
         !
      ELSE IF( program_name == 'FPMD' ) THEN
         !
         ALLOCATE( vpot( nnr, nspin ) )
         ALLOCATE( lambda(  nudx, nudx, nspin ) )
         !
      END IF
      !
      ALLOCATE( becdr( nhsa, n, 3 ) )
      !
      IF ( nosmd ) ALLOCATE( bec( nhsa, n ) )
      !
      ALLOCATE( bephi( nhsa, n ) )
      ALLOCATE( becp(  nhsa, n ) )
      !
      CALL wave_descriptor_init( wfill, ngw, ngwt, nupdwn,  nupdwn, &
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
      IF( ALLOCATED( ei1 ) )     DEALLOCATE( ei1 )
      IF( ALLOCATED( ei2 ) )     DEALLOCATE( ei2 )
      IF( ALLOCATED( ei3 ) )     DEALLOCATE( ei3 )
      IF( ALLOCATED( eigr ) )    DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) )    DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) )   DEALLOCATE( eigrb )
      IF( ALLOCATED( irb ) )     DEALLOCATE( irb )
      IF( ALLOCATED( rhopr ) )   DEALLOCATE( rhopr )
      IF( ALLOCATED( rhor ) )    DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) )    DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) )    DEALLOCATE( rhog )
      IF( ALLOCATED( bec ) )     DEALLOCATE( bec )
      IF( ALLOCATED( becdr ) )   DEALLOCATE( becdr )
      IF( ALLOCATED( bephi ) )   DEALLOCATE( bephi )
      IF( ALLOCATED( becp ) )    DEALLOCATE( becp )
      IF( ALLOCATED( ema0bg ) )  DEALLOCATE( ema0bg )
      IF( ALLOCATED( lambda ) )  DEALLOCATE( lambda )
      IF( ALLOCATED( lambdam ) ) DEALLOCATE( lambdam )
      IF( ALLOCATED( lambdap ) ) DEALLOCATE( lambdap )
      IF( ALLOCATED( kedtaur ) ) DEALLOCATE( kedtaur )
      IF( ALLOCATED( kedtaus ) ) DEALLOCATE( kedtaus )
      IF( ALLOCATED( kedtaug ) ) DEALLOCATE( kedtaug )
      IF( ALLOCATED( vpot ) )    DEALLOCATE( vpot )
      IF( ALLOCATED( taub ) )    DEALLOCATE( taub )
      !
      RETURN
      !
    END SUBROUTINE deallocate_mainvar
    !
END MODULE cp_main_variables
