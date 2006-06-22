!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!
!     EMPTY STATES
!

! ---------------------------------------------------------------------- !
      MODULE empty_states
! ---------------------------------------------------------------------- !

        USE kinds
        USE io_files, ONLY: empty_file, emptyunit

        IMPLICIT NONE
        SAVE

        PRIVATE
 
        INTEGER  :: max_emp = 0    !  maximum number of iterations
        REAL(DP) :: ethr_emp       !  threshold for convergence
        REAL(DP) :: delt_emp       !  delt for the empty states updating
        REAL(DP) :: emass_emp      !  fictitious mass for the empty states

        LOGICAL  :: prn_emp       = .FALSE.

        CHARACTER(LEN=256) :: fileempty
        LOGICAL  :: first = .TRUE.

        INTERFACE EMPTY
          MODULE PROCEDURE EMPTY_SD
        END INTERFACE

        PUBLIC :: empty , empty_init , empty_print_info, readempty, empty_cp

! ---------------------------------------------------------------------- !
      CONTAINS
! ---------------------------------------------------------------------- !

        SUBROUTINE empty_print_info(iunit)
          !
          USE electrons_module, ONLY: n_emp
          INTEGER, INTENT(IN) :: iunit
          !
          IF ( n_emp > 0 ) WRITE (iunit,620) n_emp, max_emp, delt_emp
620       FORMAT(3X,'Empty states minimization : states = ',I4, &
             ' maxiter = ',I8,' delt = ',F8.4)
          !
          RETURN
        END SUBROUTINE empty_print_info

!----------------------------------------------------------------------

        SUBROUTINE empty_init( max_emp_ , delt_emp_ , ethr_emp_ )

          INTEGER, INTENT(IN) :: max_emp_
          REAL(DP), INTENT(IN) :: delt_emp_ , ethr_emp_

          delt_emp  = delt_emp_
          emass_emp = 200.0d0
          max_emp   = max_emp_
          ethr_emp  = ethr_emp_

          RETURN
        END SUBROUTINE empty_init

!
!=----------------------------------------------------------------------------=!
!

      LOGICAL FUNCTION readempty( c_emp, ne )

        ! ...   This subroutine reads empty states from unit emptyunit

        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE io_global,          ONLY: stdout, ionode, ionode_id
        USE mp,                 ONLY: mp_bcast, mp_sum
        USE mp_wave,            ONLY: splitwf
        USE io_files,           ONLY: scradir
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw

        IMPLICIT none

        COMPLEX(DP), INTENT(OUT) :: c_emp(:,:)
        INTEGER,     INTENT(IN)    :: ne

        LOGICAL :: exst
        INTEGER :: ierr, ig, i, iss
        INTEGER :: ngw_rd, ne_rd, ngw_l
        INTEGER :: ngw_g

        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        !
        ! ... Subroutine Body
        !
        ngw_g    = ngw
        ngw_l    = ngw
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp(ngw_g) )

        IF( LEN( TRIM( scradir ) ) > 1 ) THEN
          fileempty = TRIM( scradir ) // '/' // TRIM( empty_file )
        ELSE
          fileempty = TRIM( empty_file )
        END IF

        IF ( ionode ) THEN

          INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )

          IF ( EXST ) THEN
            !
            OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
            !
            READ(emptyunit) ngw_rd, ne_rd
            !
            IF( ne > ne_rd ) THEN
              EXST = .false.
              WRITE( stdout,10)  TRIM(fileempty) 
              WRITE( stdout,20)  ngw_rd, ne_rd
              WRITE( stdout,20)  ngw_g, ne
            END IF
            !
          END IF

        END IF

 10     FORMAT('*** EMPTY STATES : wavefunctions dimensions changed  ', A )
 20     FORMAT('*** NGW = ', I8, ' NE = ', I4)

        CALL mp_bcast(exst,   ionode_id, intra_image_comm)
        CALL mp_bcast(ne_rd,  ionode_id, intra_image_comm)
        CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)

        IF ( exst ) THEN

           DO i = 1, MIN( ne, ne_rd )
              IF ( ionode ) THEN
                  READ(emptyunit) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
              END IF
              IF( i <= ne ) THEN
                  CALL splitwf(c_emp(:,i), ctmp, ngw_l, ig_l2g, me_image, &
                       nproc_image, ionode_id, intra_image_comm)
              END IF
           END DO

        END IF

        IF ( ionode .AND. EXST ) THEN
          CLOSE(emptyunit)
        END IF

        readempty = exst
        DEALLOCATE(ctmp)

        RETURN
      END FUNCTION readempty

!
!=----------------------------------------------------------------------------=!
!

      SUBROUTINE writeempty( c_emp, ne )

        ! ...   This subroutine writes empty states to unit emptyunit

        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE mp_wave,            ONLY: mergewf
        USE mp,                 ONLY: mp_sum
        USE io_files,           ONLY: scradir
        USE io_global,          ONLY: ionode, ionode_id, stdout
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw

        COMPLEX(DP), INTENT(IN) :: c_emp(:,:)
        INTEGER,     INTENT(IN) :: ne

        INTEGER :: ig, i, ngw_g, iss, ngw_l
        LOGICAL :: exst
        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        !
        ! ... Subroutine Body
        !
        ngw_g    = ngw
        ngw_l    = ngw
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp( ngw_g ) )

        IF( LEN( TRIM( scradir ) ) > 1  ) THEN
          fileempty = TRIM( scradir ) // '/' // TRIM( empty_file )
        ELSE
          fileempty = TRIM( empty_file )
        END IF

        IF( ionode ) THEN
          OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
          REWIND( emptyunit )
          WRITE (emptyunit)  ngw_g, ne
          WRITE( stdout,10)  TRIM(fileempty)
          WRITE( stdout,20)  ngw_g, ne
        END IF

 10     FORMAT('*** EMPTY STATES : writing wavefunctions  ', A )
 20     FORMAT('*** NGW = ', I8, ' NE = ', I4)

        DO i = 1, ne
           ctmp = 0.0d0
           CALL MERGEWF( c_emp(:,i), ctmp(:), ngw_l, ig_l2g, me_image, &
                         nproc_image, ionode_id, intra_image_comm )
           IF( ionode ) THEN
              WRITE (emptyunit) ( ctmp(ig), ig=1, ngw_g )
           END IF
        END DO

        IF( ionode ) THEN
          CLOSE (emptyunit)
        END IF

        DEALLOCATE(ctmp)

        RETURN
      END SUBROUTINE writeempty

!
!-------------------------------------------------------------------------
   SUBROUTINE gram_empty( tortho, eigr, betae, bec_emp, bec_occ, nkbx, c_emp, c_occ, ngwx, n_emp, n_occ )
!-----------------------------------------------------------------------
!
!     gram-schmidt orthogonalization of the empty states ( c_emp ) 
!     c_emp are orthogonalized among themself and to the occupied states c_occ
!
      USE uspp,           ONLY : nkb, nkbus
      USE cvan,           ONLY : nvb
      USE gvecw,          ONLY : ngw
      USE kinds,          ONLY : DP
      USE mp,             ONLY : mp_sum
      USE mp_global,      ONLY : intra_image_comm
      USE ions_base,      ONLY : nat
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx, n_emp, n_occ
      COMPLEX(DP)   :: eigr(ngwx,nat)
      REAL(DP)      :: bec_emp( nkbx, n_emp )
      REAL(DP)      :: bec_occ( nkbx, n_occ )
      COMPLEX(DP)   :: betae( ngwx, nkb )
      COMPLEX(DP)   :: c_emp( ngwx, n_emp )
      COMPLEX(DP)   :: c_occ( ngwx, n_occ )
      LOGICAL, INTENT(IN) :: tortho
!
      REAL(DP) :: anorm, cscnorm
      REAL(DP), ALLOCATABLE :: csc_emp( : )
      REAL(DP), ALLOCATABLE :: csc_occ( : )
      INTEGER :: i, k, inl
      EXTERNAL cscnorm
!
      ALLOCATE( csc_emp( n_emp ) )
      ALLOCATE( csc_occ( n_occ ) )
      !
      ! orthogonalize empty states to the occupied one and among each other
      !
      DO i = 1, n_emp
         !
         csc_emp = 0.0d0
         csc_occ = 0.0d0
         !
         ! compute scalar product csc_occ(k) = <c_emp(i)|c_occ(k)>
         !
         CALL smooth_csv( c_emp(1,i), c_occ(1,1), ngwx, csc_occ, n_occ )
         !
         IF( .NOT. tortho ) THEN
           !
           ! compute scalar product csc_emp(k) = <c_emp(i)|c_emp(k)>
           !
           CALL smooth_csv( c_emp(1,i), c_emp(1,1), ngwx, csc_emp, i-1 )
           !
           CALL mp_sum( csc_emp, intra_image_comm )
           !
         END IF
         !
         CALL mp_sum( csc_occ, intra_image_comm )
         !
         IF( nvb > 1 ) THEN
            !
            CALL grabec( bec_emp(1,i), nkbx, betae, c_emp(1,i), ngwx )
            !
            CALL mp_sum( bec_emp(1:nkbus,i), intra_image_comm )
            !
            CALL bec_csv( bec_emp(1,i), bec_occ, nkbx, csc_occ, n_occ )
            !
            IF( .NOT. tortho ) THEN
               CALL bec_csv( bec_emp(1,i), bec_emp, nkbx, csc_emp, i-1 )
            END IF
            !
            DO k = 1, n_occ
               DO inl = 1, nkbx
                  bec_emp( inl, i ) = bec_emp( inl, i ) - csc_occ(k) * bec_occ( inl, k )
               END DO
            END DO
            !
            IF( .NOT. tortho ) THEN
               DO k = 1, i-1
                  DO inl = 1, nkbx
                     bec_emp( inl, i ) = bec_emp( inl, i ) - csc_emp(k) * bec_emp( inl, k )
                  END DO
               END DO
            END IF
            !
         END IF
         !
         ! calculate orthogonalized c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k    csv(k)|c_occ(k)>
         !                          c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k<i  csv(k)|c_emp(k)>
         !
         DO k = 1, n_occ
            CALL DAXPY( 2*ngw, -csc_occ(k), c_occ(1,k), 1, c_emp(1,i), 1 )
         END DO
         IF( .NOT. tortho ) THEN
            DO k = 1, i - 1
               CALL DAXPY( 2*ngw, -csc_emp(k), c_emp(1,k), 1, c_emp(1,i), 1 )
            END DO
         END IF
         !
         !
         IF( .NOT. tortho ) THEN
            anorm = cscnorm( bec_emp, nkbx, c_emp, ngwx, i, n_emp )
            !
            CALL DSCAL( 2*ngw, 1.0d0/anorm, c_emp(1,i), 1 )
            !
            IF( nvb > 1 ) THEN
               CALL DSCAL( nkbx, 1.0d0/anorm, bec_emp(1,i), 1 )
            END IF
         END IF
         !
      END DO

      DEALLOCATE( csc_emp )
      DEALLOCATE( csc_occ )
      !
      RETURN
   END SUBROUTINE gram_empty


!
!=======================================================================
!

    SUBROUTINE empty_sd( nfi, tortho, c_occ, wfill, vpot, eigr, bec )

      USE pseudopotential, ONLY: nsanl
      USE cell_base, ONLY: tpiba2, omega
      USE electrons_module, ONLY: n_emp, nupdwn_emp, iupdwn_emp, n_emp_l, ei_emp, eigs
      USE electrons_base, ONLY: nspin, nupdwn
      USE cp_electronic_mass, ONLY: emass
      USE time_step, ONLY: delt
      USE forces, ONLY: dforce_all
      USE orthogonalize, ONLY: ortho
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: nproc_image
      USE check_stop, ONLY: check_stop_now
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE control_flags, ONLY: force_pairing, gamma_only, tsde
      USE reciprocal_vectors, ONLY: g, gx, gstart
      USE uspp_param, ONLY: nhm
      USE pseudopotential, ONLY: nspnl
      USE uspp,             ONLY : nkb, vkb
      USE cp_main_variables,  ONLY: ema0bg
      USE optical_properties, ONLY: opticalp
      USE wave_functions,       ONLY : wave_rand_init, elec_fakekine
      USE wave_base,            ONLY : wave_steepest, wave_verlet, frice
      USE wave_constrains,  ONLY : update_lambda
      USE wave_types, ONLY: wave_descriptor, wave_descriptor_init

      IMPLICIT NONE

! ... ARGUMENTS
!
      COMPLEX(DP), INTENT(INOUT) ::  c_occ(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) ::  wfill
      REAL (DP), INTENT(IN) ::   vpot(:,:)
      LOGICAL, INTENT(IN) :: tortho
      INTEGER, INTENT(IN) :: nfi
      COMPLEX(DP) :: eigr(:,:)
      REAL(DP) :: bec(:,:)
!
! ... LOCALS
!

      INTEGER   ::  i, k, j, iter
      INTEGER   ::  iss, ispin_wfc
      INTEGER   ::  n_occ( nspin )
      INTEGER   ::  ig, iprinte, nrl, jl, ngw, in_emp
      REAL(DP) ::  dek, ekinc, ekinc_old
      REAL(DP) :: ampre
      REAL(DP), ALLOCATABLE :: dt2bye( : )
      REAL(DP), PARAMETER :: small    = 1.0d-14
      REAL(DP), ALLOCATABLE :: bece( :, : )
      COMPLEX(DP), ALLOCATABLE :: eforce(:,:), cp_emp(:,:), c_emp(:,:)
      REAL(DP), ALLOCATABLE :: fi(:,:)
      REAL(DP), ALLOCATABLE :: gam(:,:)
      REAL(DP) :: fccc, ccc
      REAL(DP) :: verl1, verl2
      REAL(DP), ALLOCATABLE :: verl3( : )
      LOGICAL       :: gamma, tlast
      LOGICAL, SAVE :: exst

      TYPE (wave_descriptor) :: wempt  ! wave function descriptor for empty states

!
! ... SUBROUTINE BODY
!    

      ngw       = wfill%ngwl
      gamma     = wfill%gamma
      ampre     = 0.001d0

      CALL wave_descriptor_init( wempt, ngw, wfill%ngwt, nupdwn_emp, nupdwn_emp, &
                                 1, 1, nspin, 'gamma' , wfill%gzero )


      n_occ( 1 ) = wfill%nbt( 1 )
      IF( nspin == 2 ) THEN
         n_occ( 2 ) = wfill%nbt( 2 )
      END IF

      ekinc_old = 1.d+10
      ekinc     = 0.0d0

      ALLOCATE( bece( nkb, n_emp * nspin ) )
      ALLOCATE( eforce( ngw, n_emp ) )
      ALLOCATE( fi( n_emp, nspin ) )
      ALLOCATE( c_emp( ngw, n_emp * nspin ) )
      ALLOCATE( cp_emp( ngw, n_emp * nspin ) )
      ALLOCATE( dt2bye( ngw ) )
      ALLOCATE( verl3( ngw ) )

      IF( ionode ) WRITE( stdout,56)

      exst = readempty( c_emp, n_emp * nspin )

      IF( .NOT. exst ) THEN
         !
         CALL wave_rand_init( c_emp )
         !
         IF ( gstart == 2 ) THEN
            c_emp( 1, : ) = DBLE( c_emp( 1, : ) )
         END IF
         !
         DO iss = 1, wempt%nspin
            ispin_wfc = iss
            IF( force_pairing ) ispin_wfc = 1
            CALL gram_empty( .false., eigr, vkb, bece, bec, nkb, c_emp(:,iupdwn_emp(iss):), &
                              c_occ( :, :, ispin_wfc ), ngw, nupdwn_emp( iss ), nupdwn( iss ) )
         END DO
         !
      END IF

      dt2bye = delt * delt * ema0bg / emass 

      cp_emp = c_emp

      fi     = 2.0d0 / nspin

      IF( tsde ) THEN
          fccc = 1.0d0
      ELSE
          fccc = 1.0d0 / ( 1.0d0 + frice )
      END IF

      verl1 = 2.0d0 * fccc
      verl2 = 1.0d0 - verl1
      verl3 = dt2bye * fccc

      tlast = .false.

      ITERATIONS: DO iter = 1, max_emp

        ekinc = 0.0d0

        IF( iter == max_emp ) tlast = .true.

        SPIN_LOOP: DO iss = 1, nspin

          ispin_wfc = iss
          IF( force_pairing ) ispin_wfc = 1

          in_emp = iupdwn_emp(iss)

          IF( n_emp < 1 ) CYCLE SPIN_LOOP 

          bece = 0.0d0

          CALL nlsm1 ( n_emp, 1, nspnl, eigr, c_emp( 1, in_emp ), bece( 1, in_emp ) )

          ! ...   Calculate | dH / dpsi(j) >
          !
          CALL dforce_all( c_emp(:,in_emp:), fi(:,iss), eforce, &
                             vpot(:,iss), eigr, bece, nupdwn_emp(iss), in_emp )

          DO i = in_emp, in_emp + nupdwn_emp( iss ) - 1
             !
             IF ( tsde ) THEN
                CALL wave_steepest( cp_emp(:,i), c_emp(:,i), verl3, eforce(:,i-in_emp+1) )
             ELSE
                CALL wave_verlet( cp_emp(:,i), c_emp(:,i), verl1, verl2, verl3, eforce(:,i-in_emp+1) )
             END IF
             !
             if ( gstart == 2) THEN
                 cp_emp( 1, i ) = CMPLX( DBLE( cp_emp( 1, i ) ), 0.d0 )
             end if
             !
          END DO

          CALL gram_empty( tortho, eigr, vkb, bece, bec, nkb, cp_emp(:,in_emp:), &
                              c_occ( :, :, ispin_wfc ), ngw, nupdwn_emp( iss ), nupdwn( iss ) )
          !
          ! ...     Calculate Eij = < psi(i) | H | psi(j) > = < psi(i) | dH / dpsi(j) >
          !
          IF( tlast ) THEN
             !
             ALLOCATE( gam( n_emp_l(iss), n_emp ) )
             !
             DO i = 1, nupdwn_emp(iss)
                CALL update_lambda( i, gam, c_emp(:,iupdwn_emp(iss):), wempt, eforce(:,i) )
             END DO
             !
             CALL eigs( n_emp, gam, tortho, fi(:,iss), ei_emp(:,iss) )
             !
             DEALLOCATE( gam )
             !
          END IF

          IF (tortho) THEN

              CALL ortho( iss, c_emp(:,in_emp:), cp_emp(:,in_emp:), wempt)

          END IF

        END DO SPIN_LOOP

        CALL elec_fakekine( ekinc, ema0bg, emass, cp_emp, c_emp, ngw, n_emp * nspin, 1, delt )

        IF( tlast ) THEN
           EXIT ITERATIONS
        END IF

        dek = ekinc - ekinc_old
        IF( ionode ) WRITE( stdout,113) ITER, dek, ekinc

        CALL dswap( 2*SIZE( c_emp ), c_emp, 1, cp_emp, 1 )

        ! ...   check for exit
        !
        IF ( check_stop_now() ) THEN
           tlast = .true. 
        END IF

        ! ...   check for convergence
        !
        IF( ( ekinc <  ethr_emp ) .AND. ( iter > 3 ) ) THEN
           IF( ionode ) WRITE( stdout,112) 
           tlast = .true.
        END IF

        ekinc_old = ekinc

      END DO ITERATIONS


      ! CALL opticalp( nfi, omega, c_occ, wfill, c_emp, wempt, vpot, eigr, bec )

      CALL writeempty( c_emp, n_emp * nspin )

      DEALLOCATE( eforce )
      DEALLOCATE( cp_emp )
      DEALLOCATE( c_emp )
      DEALLOCATE( fi )
      DEALLOCATE( dt2bye )
      DEALLOCATE( verl3 )
      DEALLOCATE( bece )
              
 55   FORMAT(1X,I8,4F12.6)
 56   FORMAT(/,3X,'Empty states minimization starting ')
111   FORMAT(I5,2X,F12.8)
113   FORMAT(I5,2X,2D14.6)
112   FORMAT(/,3X,'Empty states: convergence achieved') 

    RETURN
    END SUBROUTINE empty_sd


!
!
!----------------------------------------------------------------------------
!
!

   SUBROUTINE empty_cp ( nfi, c0, v )

      USE kinds,                ONLY : DP
      USE control_flags,        ONLY : force_pairing, iprsta, tsde
      USE io_global,            ONLY : ionode, stdout
      USE cp_main_variables,    ONLY : eigr, rhos, ema0bg
      USE cell_base,            ONLY : omega, ibrav, h
      USE uspp,                 ONLY : becsum, vkb, nkb
      USE grid_dimensions,      ONLY : nnrx
      USE electrons_base,       ONLY : nbsp, nspin, f, nudx, iupdwn, nupdwn
      USE electrons_module,     ONLY : iupdwn_emp, nupdwn_emp, n_emp, ei_emp
      USE ions_base,            ONLY : nat, nsp
      USE gvecw,                ONLY : ngw
      USE orthogonalize_base,   ONLY : calphi, updatc
      USE electrons_base,       ONLY : nupdwn 
      USE reciprocal_vectors,   ONLY : gzero, gstart
      USE wave_functions,       ONLY : wave_rand_init, elec_fakekine
      USE wave_base,            ONLY : wave_steepest, wave_verlet, frice
      USE pseudopotential,      ONLY : nspnl
      USE cvan,                 ONLY : nvb
      USE cp_electronic_mass,   ONLY : emass
      USE time_step,            ONLY : delt
      USE orthogonalize,        ONLY : ortho
      USE check_stop,           ONLY : check_stop_now
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: nfi
      COMPLEX(DP)            :: c0(:,:,:)
      REAL(DP)               :: v(:,:)
      !
      INTEGER  :: i, iss, j, in, in_emp, iter, iter_ortho
      INTEGER  :: n_occs, n_emps, n_empx
      LOGICAL  :: exst
      !
      REAL(DP) :: fccc, ccc, csv, dt2bye, bigr
      REAL(DP) :: verl1, verl2, verl3
      REAL(DP) :: dek, ekinc, ekinc_old
      !
      REAL(DP),    ALLOCATABLE :: emadt2(:)
      REAL(DP),    ALLOCATABLE :: emaver(:)
      COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
      COMPLEX(DP), ALLOCATABLE :: c0_emp(:,:), cm_emp(:,:), phi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: bec_emp(:,:)
      REAL(DP),    ALLOCATABLE :: bephi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: becp_emp(:,:)
      REAL(DP),    ALLOCATABLE :: bec_occ(:,:)
      REAL(DP),    ALLOCATABLE :: lambda_emp(:,:,:), f_emp(:)
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      !
      ! ...  quick exit if empty states have not to be computed
      !
      IF( n_emp < 1 ) RETURN
      !
      IF( force_pairing ) &
         CALL errore(' empty_cp ', ' force pairing not implemented ', 1 )
      !
      ekinc_old = 1.d+10
      ekinc     = 0.0d0
      !
      n_occs = nupdwn( 1 )
      IF( nspin == 2 ) n_occs = n_occs + nupdwn( 2 )
      !
      n_emps = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_emps = n_emps + nupdwn_emp( 2 )
      !
      n_empx = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_empx = MAX( n_empx, nupdwn_emp( 2 ) )
      !
      ALLOCATE( c0_emp( ngw, n_empx * nspin ) )
      ALLOCATE( cm_emp( ngw, n_empx * nspin ) )
      ALLOCATE( phi_emp( ngw, n_empx * nspin ) )
      ALLOCATE( bec_emp( nkb, n_emps ) )
      ALLOCATE( bec_occ( nkb, n_occs ) )
      ALLOCATE( bephi_emp( nkb, n_emps ) )
      ALLOCATE( becp_emp( nkb, n_emps ) )
      ALLOCATE( lambda_emp( n_empx, n_empx, nspin ) )
      ALLOCATE( f_emp( n_empx * nspin ) )
      ALLOCATE( ispin_emp( n_empx * nspin ) )
      !
      phi_emp    = 0.0d0
      bec_emp    = 0.0d0
      bec_occ    = 0.0d0
      bephi_emp  = 0.0d0
      becp_emp   = 0.0d0
      lambda_emp = 0.0d0
      f_emp      = 2.0d0 / nspin
      !
      ispin_emp( 1:nupdwn_emp( 1 ) ) = 1
      IF( nspin == 2 ) ispin_emp( iupdwn_emp(2) : ) = 2
      !
      IF( ionode ) THEN
         WRITE( stdout,56)
      END IF

      exst = readempty( c0_emp, n_empx * nspin )
      !
      CALL prefor( eigr, vkb )
      !
      DO iss = 1, nspin
         in = iupdwn( iss )
         CALL nlsm1 ( nupdwn( iss ), 1, nvb, eigr, c0( 1, in, 1 ), bec_occ( 1, in ) )
      END DO
      !
      IF( .NOT. exst ) THEN
         !
         ! ...  initial random states orthogonal to filled ones
         !
         CALL wave_rand_init( c0_emp )
         !
         IF ( gzero ) THEN
            c0_emp( 1, : ) = (0.0d0, 0.0d0)
         END IF
         !
         CALL nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp )
         !
         DO iss = 1, nspin
            !
            in     = iupdwn(iss)
            in_emp = iupdwn_emp(iss)
            !
            CALL gram_empty( .false. , eigr, vkb, bec_emp( :, in_emp: ), bec_occ( :, in: ), nkb, &
                                c0_emp( :, in_emp: ), c0( :, in:, 1 ), ngw, nupdwn_emp(iss), nupdwn(iss) )
            !
         END DO
         !
      END IF
      !
      ! WRITE(stdout,*) '------- filled filled -------'
      ! DO i = 1, nupdwn( 1 )
      !    DO j = 1, nupdwn( 1 )
      !       CALL dotcsv( csv, eigr, c0(1,i,1), c0(1,j,1), ngw )
      !       WRITE(stdout,*) i, j, csv
      !    END DO
      ! END DO
       ! WRITE(stdout,*) '------- empty filled -------'
       ! DO i = 1, n_emp 
       !   DO j = 1, nupdwn( 1 )
       !       CALL dotcsv( csv, eigr, c0_emp(1,i), c0(1,j,1), ngw )
       !       WRITE(stdout,*) i, j, csv
       !    END DO
       ! END DO
       ! IF( nspin == 2 ) THEN
       ! DO i = iupdwn_emp(2), iupdwn_emp(2) + nupdwn_emp( 2 ) - 1
       !   DO j = iupdwn(2), iupdwn(2) + nupdwn( 2 ) - 1
       !       CALL dotcsv( csv, eigr, c0_emp(1,i), c0(1,j,1), ngw )
       !       WRITE(stdout,*) i, j, csv
       !    END DO
       ! END DO
       ! END IF
       ! WRITE(stdout,*) '------- empty empty -------'
       ! DO i = 1, n_emp
       !    DO j = 1, n_emp
       !       CALL dotcsv( csv, eigr, c0_emp(1,i,1), c0_emp(1,j,1), ngw )
       !       WRITE(stdout,*) i, j, csv
       !    END DO
       ! END DO
      

      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp )
      !
      ! ...  set verlet variables
      !
      IF( tsde ) THEN
          fccc = 1.0d0
      ELSE   
          fccc = 1.0d0 / ( 1.0d0 + frice )
      END IF
      !
      verl1 = 2.0d0 * fccc
      verl2 = 1.0d0 - verl1
      verl3 = 1.0d0 * fccc
      !
      ALLOCATE( c2( ngw ) )
      ALLOCATE( c3( ngw ) )
      ALLOCATE( emadt2( ngw ) )
      ALLOCATE( emaver( ngw ) )

      dt2bye = delt * delt / emass

      ccc    = fccc   * dt2bye
      emadt2 = dt2bye * ema0bg
      emaver = emadt2 * verl3
 
      cm_emp = c0_emp
      !
      ITERATIONS: DO iter = 1, max_emp

         DO i = 1, n_emps, 2
            !
            CALL dforce( bec_emp, vkb, i, c0_emp(1,i), c0_emp(1,i+1), c2, c3, v, &
                         ispin_emp, f_emp, n_emps, nspin )
            !
            IF( tsde ) THEN
               CALL wave_steepest( cm_emp(:, i  ), c0_emp(:, i  ), emaver, c2 )
               CALL wave_steepest( cm_emp(:, i+1), c0_emp(:, i+1), emaver, c3 )
            ELSE
              CALL wave_verlet( cm_emp(:, i  ), c0_emp(:, i  ), verl1, verl2, emaver, c2 )
              CALL wave_verlet( cm_emp(:, i+1), c0_emp(:, i+1), verl1, verl2, emaver, c3 )
            END IF
            !
            if ( gstart == 2) THEN
               cm_emp(1,  i)=CMPLX(DBLE(cm_emp(1,  i)),0.d0)
               cm_emp(1,i+1)=CMPLX(DBLE(cm_emp(1,i+1)),0.d0)
            end if
            !
         END DO

         DO iss = 1, nspin
            !
            in     = iupdwn(iss)
            in_emp = iupdwn_emp(iss)
            !
            CALL gram_empty( .true. , eigr, vkb, bec_emp( :, in_emp: ), bec_occ( :, in: ), nkb, &
                         cm_emp( :, in_emp: ), c0( :, in:, 1 ), ngw, nupdwn_emp(iss), nupdwn(iss) )
            !
         END DO

         ! ... calphi calculates phi
         ! ... the electron mass rises with g**2
         !
         CALL calphi( c0_emp, ngw, bec_emp, nkb, vkb, phi_emp, n_emps, ema0bg )
         !
         CALL ortho( eigr, cm_emp, phi_emp, ngw, lambda_emp, n_empx, &
                        bigr, iter_ortho, ccc, bephi_emp, becp_emp, n_emps, nspin, &
                        nupdwn_emp, iupdwn_emp )
         !
         DO iss = 1, nspin
            !
            CALL updatc( ccc, n_emps, lambda_emp(:,:,iss), n_empx, phi_emp, ngw, &
                         bephi_emp, nkb, becp_emp, bec_emp, &
                         cm_emp, nupdwn_emp(iss), iupdwn_emp(iss) )
            !
         END DO
         !
         CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, bec_emp )
         !
         CALL elec_fakekine( ekinc, ema0bg, emass, c0_emp, cm_emp, ngw, n_emps, 1, delt )
         !
         CALL dswap( 2*SIZE( c0_emp ), c0_emp, 1, cm_emp, 1 )

         dek = ekinc - ekinc_old
         IF( ionode ) WRITE( stdout,113) ITER, dek, ekinc
      
         ! ...   check for exit
         !
         IF ( check_stop_now() ) THEN
            EXIT ITERATIONS
         END IF
      
         ! ...   check for convergence
         !     
         IF( ( ekinc <  ethr_emp ) .AND. ( iter > 3 ) ) THEN
            IF( ionode ) WRITE( stdout,112)
            EXIT ITERATIONS
         END IF
      
         ekinc_old = ekinc


      END DO ITERATIONS

      ! ...  Compute eigenvalues

      CALL eigs0( ei_emp, .false., nspin, nupdwn_emp, iupdwn_emp, .true., f_emp, &
                  n_empx, lambda_emp, n_empx )

      ! ...   Save emptystates to disk

      CALL writeempty( c0_emp, n_empx * nspin )

      ! 
      DEALLOCATE( ispin_emp )
      DEALLOCATE( lambda_emp )
      DEALLOCATE( f_emp )
      DEALLOCATE( emadt2 )
      DEALLOCATE( emaver )
      DEALLOCATE( c2 )
      DEALLOCATE( c3 )
      DEALLOCATE( c0_emp )
      DEALLOCATE( cm_emp )
      DEALLOCATE( phi_emp )
      DEALLOCATE( bec_emp )
      DEALLOCATE( bec_occ )
      DEALLOCATE( bephi_emp )
      DEALLOCATE( becp_emp )
  
 55   FORMAT(1X,I8,4F12.6)
 56   FORMAT(/,3X,'Empty states minimization starting ')
111   FORMAT(I5,2X,F12.8)
113   FORMAT(I5,2X,2D14.6)
112   FORMAT(/,3X,'Empty states: convergence achieved')

      RETURN
   END SUBROUTINE empty_cp

!
!
! ---------------------------------------------------------------------- !
      END MODULE empty_states
! ---------------------------------------------------------------------- !
