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
  USE control_flags,     ONLY : lwfnscf, lwfpbe0, lwfpbe0nscf  ! Lingzhu Kong
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
  REAL(DP),    ALLOCATABLE :: rhopr(:,:)   ! Lingzhu Kong
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
  !==========================================================================
  ! Lingzhu Kong
            
     INTEGER  :: my_nbspx
     INTEGER  :: nord2            ! order of expansion ( points on one side)
     INTEGER  :: lap_neig(3,3)    ! new directions
     REAL(DP) :: lap_dir_step(3)  ! step in the new directions
     INTEGER  :: lap_dir_num      ! number of new directions
     REAL(DP) :: b_lap(6)         ! coefficients of the directions
     INTEGER  :: lap_dir(3)       ! activeness of the new directions

     INTEGER  np_in_sp, np_in_sp2  ! number of grid points in the 1st sphere and the shell between 1st and 2nd sphere

! conversion between 3D index (i,j,k) and 1D index np
     INTEGER,     ALLOCATABLE     :: odtothd_in_sp(:,:)
     INTEGER,     ALLOCATABLE     :: thdtood_in_sp(:,:,:)
     INTEGER,     ALLOCATABLE     :: thdtood(:,:,:)
     REAL(DP),    ALLOCATABLE     :: xx_in_sp(:)
     REAL(DP),    ALLOCATABLE     :: yy_in_sp(:)
     REAL(DP),    ALLOCATABLE     :: zz_in_sp(:)

     REAL(DP),    ALLOCATABLE     :: selfv(:,:,:)
     REAL(DP),    ALLOCATABLE     :: pairv(:,:,:,:)
     REAL(DP),    ALLOCATABLE     :: exx_potential(:, :)

     REAL(DP),    ALLOCATABLE     :: clm(:,:)
     REAL(DP),    ALLOCATABLE     :: coeke(:,:)
     REAL(DP),    ALLOCATABLE     :: vwc(:,:)
     INTEGER  ::  lmax
     INTEGER  ::  n_exx =0
!==========================================================================

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
      USE descriptors, ONLY: la_descriptor, descla_init
!==============================================================================
!Lingzhu Kong
      USE mp_global,               ONLY  : nproc_image
      USE fft_base,                ONLY  : dffts
      USE electrons_base,          ONLY  : nbsp
      USE wannier_base,            ONLY  : neigh, exx_ps_rcut, exx_me_rcut, vnbsp
      USE control_flags,           ONLY  : lwfnscf, lwfpbe0, lwfpbe0nscf
!===============================================================================

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
      ELSE
            ALLOCATE( drhog( 1, 1, 1, 1 ) )
            ALLOCATE( drhor( 1, 1, 1, 1 ) )
      END IF
!=========================================================================
!Lingzhu Kong
      IF ( lwfpbe0 .or. lwfpbe0nscf ) THEN

         lmax=6
         nord2 = 3
         ALLOCATE( clm(0:lmax, 0:lmax) )
         ALLOCATE( coeke(-nord2:nord2, 6))
         ALLOCATE( exx_potential(dffts%nnr,nbsp), stat=ierr)

         CALL getnpinsp(exx_ps_rcut, exx_me_rcut, np_in_sp, np_in_sp2 )

         ALLOCATE( odtothd_in_sp(3, np_in_sp + np_in_sp2 ), stat=ierr )
         ALLOCATE( thdtood_in_sp(nr1, nr2, nr3), stat=ierr )
         ALLOCATE( thdtood(nr1, nr2, nr3), stat=ierr )
         ALLOCATE( xx_in_sp(1:np_in_sp+np_in_sp2), stat=ierr )
         ALLOCATE( yy_in_sp(1:np_in_sp+np_in_sp2), stat=ierr )
         ALLOCATE( zz_in_sp(1:np_in_sp+np_in_sp2), stat=ierr )

         my_nbspx   = nbsp / nproc_image

         IF( MOD(nbsp, nproc_image) /= 0)THEN
            my_nbspx = my_nbspx + 1
         ENDIF
         print *, 'my_nbspx =', my_nbspx
      END IF

      IF ( lwfpbe0nscf .or. lwfnscf ) ALLOCATE( rhopr( nnr, nspin ) )

      IF ( lwfpbe0nscf ) THEN
         ALLOCATE( vwc(3, vnbsp) )
         ALLOCATE( pairv( np_in_sp, 2, neigh, my_nbspx), stat=ierr )
         pairv (:,:,:,:) = 0.d0
      ENDIF

      IF ( lwfpbe0 ) THEN
         ALLOCATE( selfv( np_in_sp, 2,          my_nbspx), stat=ierr )
         ALLOCATE( pairv( np_in_sp, 2, neigh/2, my_nbspx), stat=ierr )
         selfv (:,:,:) = 0.d0
         pairv (:,:,:,:) = 0.d0
      ENDIF
!==========================================================================
      !
      !  Compute local dimensions for lambda matrixes
      !

      ALLOCATE( descla( nspin ) )
      !
      DO iss = 1, nspin
         CALL descla_init( descla( iss ), nupdwn( iss ), nudx, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
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
      ALLOCATE( bephi( nhsa, nspin*nrcx ) )
      ALLOCATE( becp_bgrp( nhsa, nbspx_bgrp ) )  
      !
      IF ( tpre ) THEN
        ALLOCATE( dbec( nhsa, 2*nrcx, 3, 3 ) )
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
!====================================================================
!Lingzhu Kong
      IF ( lwfpbe0 )THEN
         IF( ALLOCATED( selfv ) )          DEALLOCATE( selfv )
      ENDIF

      IF ( lwfpbe0nscf .or. lwfnscf)THEN
         IF( ALLOCATED( rhopr ) )          DEALLOCATE( rhopr )
      ENDIF
         
      IF ( lwfpbe0nscf )THEN
         IF( ALLOCATED( vwc)    )          DEALLOCATE( vwc )
      ENDIF
      IF ( lwfpbe0 .or. lwfpbe0nscf ) THEN
         IF( ALLOCATED( pairv ) )          DEALLOCATE( pairv )
         IF( ALLOCATED( exx_potential ) )  DEALLOCATE( exx_potential )
         IF( ALLOCATED( odtothd_in_sp ) )  DEALLOCATE(odtothd_in_sp )
         IF( ALLOCATED( thdtood_in_sp ) )  DEALLOCATE(thdtood_in_sp )
         IF( ALLOCATED( thdtood  ))        DEALLOCATE(thdtood)
         IF( ALLOCATED( xx_in_sp ))        DEALLOCATE(xx_in_sp )
         IF( ALLOCATED( yy_in_sp ))        DEALLOCATE(yy_in_sp )
         IF( ALLOCATED( zz_in_sp ))        DEALLOCATE(zz_in_sp )
         IF( ALLOCATED( clm )     )        DEALLOCATE(clm)
         IF( ALLOCATED( coeke)    )        DEALLOCATE(coeke)
      END IF
!===================================================================
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
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE distribute_lambda( lambda_repl, lambda_dist, desc )
       USE descriptors
       REAL(DP), INTENT(IN)  :: lambda_repl(:,:)
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, j, ic, ir
       IF( desc%active_node > 0 ) THEN
          ir = desc%ir
          ic = desc%ic
          DO j = 1, desc%nc
             DO i = 1, desc%nr
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
       USE descriptors
       REAL(DP), INTENT(IN)  :: bec_repl(:,:)
       REAL(DP), INTENT(OUT) :: bec_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc(:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nrcx
       !
       IF( desc( 1 )%active_node > 0 ) THEN
          !
          bec_dist = 0.0d0
          !
          ir = desc( 1 )%ir
          DO i = 1, desc( 1 )%nr
             bec_dist( :, i ) = bec_repl( :, i + ir - 1 )
          END DO
          !
          IF( nspin == 2 ) THEN
             n     = desc( 1 )%n  !  number of states with spin 1 ( nupdw(1) )
             nrcx  = desc( 1 )%nrcx   !  array elements reserved for each spin ( bec(:,2*nrcx) )
             ir = desc( 2 )%ir
             DO i = 1, desc( 2 )%nr
                bec_dist( :, i + nrcx ) = bec_repl( :, i + ir - 1 + n )
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
       USE descriptors
       REAL(DP), INTENT(IN)  :: zmat_repl(:,:)
       REAL(DP), INTENT(OUT) :: zmat_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, ii, j, me, np
       me = desc%mype
       np = desc%npc * desc%npr
       IF( desc%active_node > 0 ) THEN
          DO j = 1, desc%n
             ii = me + 1
             DO i = 1, desc%nrl
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
       USE descriptors
       REAL(DP), INTENT(OUT) :: lambda_repl(:,:)
       REAL(DP), INTENT(IN)  :: lambda_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, j, ic, ir
       lambda_repl = 0.0d0
       IF( desc%active_node > 0 ) THEN
          ir = desc%ir
          ic = desc%ic
          DO j = 1, desc%nc
             DO i = 1, desc%nr
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
       USE descriptors
       USE io_global, ONLY : stdout
       REAL(DP), INTENT(OUT) :: bec_repl(:,:)
       REAL(DP), INTENT(IN)  :: bec_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc(:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nrcx, iss
       !
       bec_repl = 0.0d0
       !
       !  bec is distributed across row processor, the first column is enough
       !
       IF( desc( 1 )%active_node > 0 .AND. ( desc( 1 )%myc == 0 ) ) THEN
          ir = desc( 1 )%ir
          DO i = 1, desc( 1 )%nr
             bec_repl( :, i + ir - 1 ) = bec_dist( :, i )
          END DO
          IF( nspin == 2 ) THEN
             n  = desc( 1 )%n   ! number of states with spin==1 ( nupdw(1) )
             nrcx = desc( 1 )%nrcx ! array elements reserved for each spin ( bec(:,2*nrcx) )
             ir = desc( 2 )%ir
             DO i = 1, desc( 2 )%nr
                bec_repl( :, i + ir - 1 + n ) = bec_dist( :, i + nrcx )
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
       USE descriptors
       REAL(DP), INTENT(OUT) :: zmat_repl(:,:)
       REAL(DP), INTENT(IN)  :: zmat_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, ii, j, me, np, nrl
       zmat_repl = 0.0d0
       me = desc%mype
       np = desc%npc * desc%npr
       nrl = desc%nrl
       IF( desc%active_node > 0 ) THEN
          DO j = 1, desc%n
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
       USE descriptors
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: i, j
       REAL(DP), INTENT(IN)  :: val
       TYPE(la_descriptor), INTENT(IN)  :: desc
       IF( desc%active_node > 0 ) THEN
          IF( ( i >= desc%ir ) .AND. ( i - desc%ir + 1 <= desc%nr ) ) THEN
             IF( ( j >= desc%ic ) .AND. ( j - desc%ic + 1 <= desc%nc ) ) THEN
                lambda_dist( i - desc%ir + 1, j - desc%ic + 1 ) = val
             END IF
          END IF
       END IF
       RETURN
    END SUBROUTINE setval_lambda
    !
    !
END MODULE cp_main_variables
