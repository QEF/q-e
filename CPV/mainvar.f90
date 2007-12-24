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
  USE parameters,        ONLY : nsx, nacx
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
  INTEGER,  ALLOCATABLE :: descla(:,:) ! descriptor of the lambda distribution
                                       ! see descriptors_module
  INTEGER :: nlax = 0                  ! leading dimension of the distribute (by block) lambda matrix 
  INTEGER :: nlam = 1                  ! dimension of lambda matrix, can be 1 or nlax depending on la_proc
  INTEGER :: nrlx = 0                  ! leading dimension of the distribute (by row  ) lambda matrix
  LOGICAL :: la_proc = .FALSE.         ! indicate if a proc own a block of lambda
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
  ! vpot  = potential in r space (dense grid)
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  REAL(DP),    ALLOCATABLE :: rhor(:,:), rhos(:,:)
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
                                 gzero, nudx, smd )
      !------------------------------------------------------------------------
      !
      USE mp_global,   ONLY: np_ortho, me_ortho, intra_image_comm, ortho_comm, &
                             me_image, ortho_comm_id
      USE mp,          ONLY: mp_max, mp_min
      USE descriptors, ONLY: descla_siz_ , descla_init , nlax_ , la_nrlx_ , lambda_node_
      !
      INTEGER,           INTENT(IN) :: ngw, ngwt, ngb, ngs, ng, nr1, nr2, nr3, &
                                       nnr, nnrsx, nat, nax, nsp, nspin, &
                                       n, nx, n_emp, nhsa, nr1x, nr2x, npl
      INTEGER,           INTENT(IN) :: nupdwn(:)
      LOGICAL,           INTENT(IN) :: gzero
      INTEGER,           INTENT(IN) :: nudx
      LOGICAL, OPTIONAL, INTENT(IN) :: smd
      !
      LOGICAL  :: nosmd
      INTEGER  :: iss
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      nosmd = .TRUE.
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
      ALLOCATE( vpot( nnr, nspin ) )
      ALLOCATE( rhos( nnrsx, nspin ) )
      ALLOCATE( rhog( ng,    nspin ) )
      !
      !  Compute local dimensions for lambda matrixes
      !

      ALLOCATE( descla( descla_siz_ , nspin ) )
      !
      nlax = 0
      nrlx = 0
      DO iss = 1, nspin
         CALL descla_init( descla( :, iss ), nupdwn( iss ), nudx, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
         nlax = MAX( nlax, descla( nlax_ , iss ) )
         nrlx = MAX( nrlx, descla( la_nrlx_ , iss ) )
         IF( descla( lambda_node_ , iss ) > 0 ) la_proc = .TRUE.
      END DO
      !
      nlam = 1
      IF( la_proc ) nlam = nlax
      !
      !  ... End with lambda dimensions
      !
      !
      IF( program_name == 'CP90' ) THEN
         !
         if ( abivol.or.abisur ) then
            !
            allocate(rho_gaus(nnr))
            allocate(v_vol(nnr))
            if (jellium.or.t_gauss) allocate(posv(3,nr1*nr2*nr3))
            if (t_gauss) allocate(f_vol(3,nax,nsx))
            !
         end if
         !
         IF ( nosmd ) THEN
            !
            ALLOCATE( lambda(  nlam, nlam, nspin ) )
            ALLOCATE( lambdam( nlam, nlam, nspin ) )
            ALLOCATE( lambdap( nlam, nlam, nspin ) )
            !
         END IF
         !
      ELSE IF( program_name == 'FPMD' ) THEN
         !
         ALLOCATE( lambda(  nlam, nlam, nspin ) )
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
      IF( ALLOCATED( descla ) )  DEALLOCATE( descla )
      !
      RETURN
      !
    END SUBROUTINE deallocate_mainvar
    !
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE distribute_lambda( lambda_repl, lambda_dist, desc )
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       REAL(DP), INTENT(IN)  :: lambda_repl(:,:)
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, j, ic, ir
       IF( desc( lambda_node_ ) > 0 ) THEN
          ir = desc( ilar_ )       
          ic = desc( ilac_ )       
          DO j = 1, desc( nlac_ )
             DO i = 1, desc( nlar_ )
                lambda_dist( i, j ) = lambda_repl( i + ir - 1, j + ic - 1 )
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_lambda
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE distribute_zmat( zmat_repl, zmat_dist, desc )
       USE descriptors, ONLY: lambda_node_ , la_nrl_ , la_me_ , la_npr_ , la_npc_ , la_n_
       REAL(DP), INTENT(IN)  :: zmat_repl(:,:)
       REAL(DP), INTENT(OUT) :: zmat_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, ii, j, me, np
       me = desc( la_me_ )
       np = desc( la_npc_ ) * desc( la_npr_ )
       IF( desc( lambda_node_ ) > 0 ) THEN
          DO j = 1, desc( la_n_ )
             ii = me + 1
             DO i = 1, desc( la_nrl_ )
                zmat_dist( i, j ) = zmat_repl( ii, j )
                ii = ii + np
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_zmat
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_lambda( lambda_repl, lambda_dist, desc )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       REAL(DP), INTENT(OUT) :: lambda_repl(:,:)
       REAL(DP), INTENT(IN)  :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, j, ic, ir
       lambda_repl = 0.0d0
       IF( desc( lambda_node_ ) > 0 ) THEN
          ir = desc( ilar_ )       
          ic = desc( ilac_ )       
          DO j = 1, desc( nlac_ )
             DO i = 1, desc( nlar_ )
                lambda_repl( i + ir - 1, j + ic - 1 ) = lambda_dist( i, j )
             END DO
          END DO
       END IF
       CALL mp_sum( lambda_repl, intra_image_comm )
       RETURN
    END SUBROUTINE collect_lambda
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_zmat( zmat_repl, zmat_dist, desc )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , la_nrl_ , la_me_ , la_npr_ , la_npc_ , la_n_
       REAL(DP), INTENT(OUT) :: zmat_repl(:,:)
       REAL(DP), INTENT(IN)  :: zmat_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, ii, j, me, np, nrl
       zmat_repl = 0.0d0
       me = desc( la_me_ )
       np = desc( la_npc_ ) * desc( la_npr_ )
       nrl = desc( la_nrl_ )
       IF( desc( lambda_node_ ) > 0 ) THEN
          DO j = 1, desc( la_n_ )
             ii = me + 1
             DO i = 1, nrl
                zmat_repl( ii, j ) = zmat_dist( i, j )
                ii = ii + np
             END DO
          END DO
       END IF
       CALL mp_sum( zmat_repl, intra_image_comm )
       RETURN
    END SUBROUTINE collect_zmat
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE setval_lambda( lambda_dist, i, j, val, desc )
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: i, j
       REAL(DP), INTENT(IN)  :: val
       INTEGER,  INTENT(IN)  :: desc(:)
       IF( desc( lambda_node_ ) > 0 ) THEN
          IF( ( i >= desc( ilar_ ) ) .AND. ( i - desc( ilar_ ) + 1 <= desc( nlar_ ) ) ) THEN
             IF( ( j >= desc( ilac_ ) ) .AND. ( j - desc( ilac_ ) + 1 <= desc( nlac_ ) ) ) THEN
                lambda_dist( i - desc( ilar_ ) + 1, j - desc( ilac_ ) + 1 ) = val
             END IF
          END IF
       END IF
       RETURN
    END SUBROUTINE setval_lambda
    !
    !
END MODULE cp_main_variables
