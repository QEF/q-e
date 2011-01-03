!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
  REAL(DP), ALLOCATABLE :: becp_dist(:,:)  ! distributed becp (ortho group)
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
  INTEGER,  ALLOCATABLE :: descla(:,:) ! descriptor of the lambda distribution
                                       ! see descriptors_module
  INTEGER :: nlax = 0                  ! leading dimension of the distribute (by block) lambda matrix 
  INTEGER :: nlam = 1                  ! dimension of lambda matrix, can be 1 or nlax depending on la_proc
  INTEGER :: nrlx = 0                  ! leading dimension of the distribute (by row  ) lambda matrix
  LOGICAL :: la_proc = .FALSE.         ! indicate if a proc own a block of lambda
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
                             me_bgrp, ortho_comm_id
      USE mp,          ONLY: mp_max, mp_min
      USE descriptors, ONLY: descla_siz_ , descla_init , nlax_ , la_nrlx_ , lambda_node_
      !
      INTEGER,           INTENT(IN) :: ngw, ngw_g, ngb, ngs, ng, nr1,nr2,nr3, &
                                       nnr, nrxxs, nat, nax, nsp, nspin, &
                                       n, nx, nhsa, nr1x, nr2x, npl
      INTEGER,           INTENT(IN) :: nupdwn(:)
      INTEGER,           INTENT(IN) :: gstart, nudx
      LOGICAL,           INTENT(IN) :: tpre
      INTEGER,           INTENT(IN) :: nbspx_bgrp
      !
      INTEGER  :: iss
      LOGICAL  :: gzero
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      ALLOCATE( eigr( ngw, nat ) )
      ALLOCATE( sfac( ngs, nsp ) )
      ALLOCATE( eigrb( ngb, nat ) )
      ALLOCATE( irb( 3, nat ) )
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
      ALLOCATE( ema0bg( ngw ) )
      !
      ALLOCATE( rhor( nnr, nspin ) )
      ALLOCATE( vpot( nnr, nspin ) )
      ALLOCATE( rhos( nrxxs, nspin ) )
      ALLOCATE( rhog( ng,    nspin ) )
      IF ( tpre ) THEN
            ALLOCATE( drhog( ng,  nspin, 3, 3 ) )
            ALLOCATE( drhor( nnr, nspin, 3, 3 ) )
      END IF
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
         if ( abivol.or.abisur ) then
            !
            allocate(rho_gaus(nnr))
            allocate(v_vol(nnr))
            if (jellium.or.t_gauss) allocate(posv(3,nr1*nr2*nr3))
            if (t_gauss) allocate(f_vol(3,nax,nsp))
            !
         end if
         !
      ALLOCATE( lambda(  nlam, nlam, nspin ) )
      ALLOCATE( lambdam( nlam, nlam, nspin ) )
      ALLOCATE( lambdap( nlam, nlam, nspin ) )
      !
      ! becdr, distributed over row processors of the ortho group
      !
      ALLOCATE( becdr_bgrp( nhsa, nbspx_bgrp, 3 ) )  
      !
      ALLOCATE( bec_bgrp( nhsa, nbspx_bgrp ) )
      !
      ALLOCATE( bephi( nhsa, nspin*nlax ) )
      ALLOCATE( becp_dist( nhsa, nlax*nspin ) )  
      !
      IF ( tpre ) THEN
        ALLOCATE( dbec( nhsa, 2*nlax, 3, 3 ) )
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
      IF( ALLOCATED( becp_dist ) )    DEALLOCATE( becp_dist )
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
    SUBROUTINE distribute_bec( bec_repl, bec_dist, desc, nspin )
       USE descriptors, ONLY: lambda_node_ , ilar_ , nlar_ , la_n_ , nlax_
       REAL(DP), INTENT(IN)  :: bec_repl(:,:)
       REAL(DP), INTENT(OUT) :: bec_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:,:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nlax
       !
       IF( desc( lambda_node_ , 1 ) > 0 ) THEN
          !
          bec_dist = 0.0d0
          !
          ir = desc( ilar_ , 1 )
          DO i = 1, desc( nlar_ , 1 )
             bec_dist( :, i ) = bec_repl( :, i + ir - 1 )
          END DO
          !
          IF( nspin == 2 ) THEN
             n     = desc( la_n_ , 1 )  !  number of states with spin 1 ( nupdw(1) )
             nlax  = desc( nlax_ , 1 )   !  array elements reserved for each spin ( bec(:,2*nlax) )
             ir = desc( ilar_ , 2 )
             DO i = 1, desc( nlar_ , 2 )
                bec_dist( :, i + nlax ) = bec_repl( :, i + ir - 1 + n )
             END DO
          END IF
          !
       END IF
       RETURN
    END SUBROUTINE distribute_bec
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
       USE mp_global,   ONLY: intra_bgrp_comm
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
       CALL mp_sum( lambda_repl, intra_bgrp_comm )
       RETURN
    END SUBROUTINE collect_lambda
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_bec( bec_repl, bec_dist, desc, nspin )
       USE mp_global,   ONLY: intra_bgrp_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , ilar_ , nlar_ , la_myc_ , nlax_ , la_n_
       USE io_global, ONLY : stdout
       REAL(DP), INTENT(OUT) :: bec_repl(:,:)
       REAL(DP), INTENT(IN)  :: bec_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:,:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nlax, iss
       !
       bec_repl = 0.0d0
       !
       !  bec is distributed across row processor, the first column is enough
       !
       IF( desc( lambda_node_ , 1 ) > 0 .AND. ( desc( la_myc_ , 1 ) == 0 ) ) THEN
          ir = desc( ilar_ , 1 )
          DO i = 1, desc( nlar_ , 1 )
             bec_repl( :, i + ir - 1 ) = bec_dist( :, i )
          END DO
          IF( nspin == 2 ) THEN
             n  = desc( la_n_ , 1 )   ! number of states with spin==1 ( nupdw(1) )
             nlax = desc( nlax_ , 1 ) ! array elements reserved for each spin ( bec(:,2*nlax) )
             ir = desc( ilar_ , 2 )
             DO i = 1, desc( nlar_ , 2 )
                bec_repl( :, i + ir - 1 + n ) = bec_dist( :, i + nlax )
             END DO
          END IF
       END IF
       !
       CALL mp_sum( bec_repl, intra_bgrp_comm )
       !
       RETURN
    END SUBROUTINE collect_bec
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_zmat( zmat_repl, zmat_dist, desc )
       USE mp_global,   ONLY: intra_bgrp_comm
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
       CALL mp_sum( zmat_repl, intra_bgrp_comm )
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
