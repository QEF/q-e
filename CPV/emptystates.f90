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


   LOGICAL FUNCTION readempty_x( c_emp, ne )

        ! ...   This subroutine reads empty states from unit emptyunit

        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE io_global,          ONLY: stdout, ionode, ionode_id
        USE mp,                 ONLY: mp_bcast, mp_sum
        USE mp_wave,            ONLY: splitwf
        USE io_files,           ONLY: outdir
        USE io_files,           ONLY: empty_file, emptyunit
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw

        IMPLICIT none

        COMPLEX(DP), INTENT(OUT) :: c_emp(:,:)
        INTEGER,     INTENT(IN)    :: ne

        LOGICAL :: exst
        INTEGER :: ierr, ig, i, iss
        INTEGER :: ngw_rd, ne_rd, ngw_l
        INTEGER :: ngw_g

        CHARACTER(LEN=256) :: fileempty

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

        IF( LEN( TRIM( outdir ) ) > 1 ) THEN
          fileempty = TRIM( outdir ) // '/' // TRIM( empty_file )
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

        readempty_x = exst
        DEALLOCATE(ctmp)

        RETURN
   END FUNCTION readempty_x

!
!=----------------------------------------------------------------------------=!
!

   SUBROUTINE writeempty_x( c_emp, ne )

        ! ...   This subroutine writes empty states to unit emptyunit

        USE xml_io_base
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE mp_wave,            ONLY: mergewf
        USE mp,                 ONLY: mp_sum
        USE io_files,           ONLY: empty_file, emptyunit, outdir
        USE io_global,          ONLY: ionode, ionode_id, stdout
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw

        COMPLEX(DP), INTENT(IN) :: c_emp(:,:)
        INTEGER,     INTENT(IN) :: ne

        INTEGER :: ig, i, ngw_g, iss, ngw_l
        LOGICAL :: exst
        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        CHARACTER(LEN=256) :: fileempty
        !
        ! ... Subroutine Body
        !
        ngw_g    = ngw
        ngw_l    = ngw
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp( ngw_g ) )

        IF( LEN( TRIM( outdir ) ) > 1  ) THEN
          fileempty = TRIM( outdir ) // '/' // TRIM( empty_file )
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
   END SUBROUTINE writeempty_x

!
!-------------------------------------------------------------------------
   SUBROUTINE gram_empty_x  &
      ( tortho, eigr, betae, bec_emp, bec_occ, nkbx, c_emp, c_occ, ngwx, n_emp, n_occ )
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
   END SUBROUTINE gram_empty_x



!
!
!----------------------------------------------------------------------------
!
!

   SUBROUTINE empty_cp_x ( nfi, c0, v )

      USE kinds,                ONLY : DP
      USE control_flags,        ONLY : iprsta, tsde, program_name
      USE io_global,            ONLY : ionode, stdout
      USE cp_main_variables,    ONLY : eigr, ema0bg, collect_lambda
      USE descriptors,          ONLY : descla_siz_ , descla_init
      USE cell_base,            ONLY : omega
      USE uspp,                 ONLY : vkb, nkb
      USE grid_dimensions,      ONLY : nnrx
      USE electrons_base,       ONLY : nbsp, nspin, f, nudx, iupdwn, nupdwn
      USE electrons_module,     ONLY : iupdwn_emp, nupdwn_emp, n_emp, ei_emp, &
                                       max_emp, ethr_emp
      USE ions_base,            ONLY : nat, nsp
      USE gvecw,                ONLY : ngw
      USE orthogonalize_base,   ONLY : calphi, updatc
      USE reciprocal_vectors,   ONLY : gzero, gstart
      USE wave_base,            ONLY : wave_steepest, wave_verlet, frice
      USE cvan,                 ONLY : nvb
      USE cp_electronic_mass,   ONLY : emass
      USE time_step,            ONLY : delt
      USE check_stop,           ONLY : check_stop_now
      USE cp_interfaces,        ONLY : writeempty, readempty, gram_empty, ortho, &
                                       wave_rand_init, elec_fakekine, crot, dforce
      USE mp,                   ONLY : mp_comm_split, mp_comm_free
      USE mp_global,            ONLY : intra_image_comm, me_image
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: nfi
      COMPLEX(DP)            :: c0(:,:)
      REAL(DP)               :: v(:,:)
      !
      INTEGER  :: i, iss, j, in, in_emp, iter, iter_ortho
      INTEGER  :: n_occs, n_emps, n_empx, issw, n
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
      REAL(DP),    ALLOCATABLE :: lambda_rep(:,:)
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      !
      INTEGER, SAVE :: np_emp(2), me_emp(2), emp_comm, color
      INTEGER, SAVE :: desc_emp( descla_siz_ , 2 )
      LOGICAL, SAVE :: first = .true.
      !
      ! ...  quick exit if empty states have not to be computed
      !
      IF( n_emp < 1 ) RETURN
      !
      ekinc_old = 1.d+10
      ekinc     = 0.0d0
      !
      !  Here set the group of processors for empty states
      !
      IF( .NOT. first ) THEN
         CALL mp_comm_free( emp_comm )
      END IF 
      !
      np_emp = 1
      IF( me_image < np_emp(1) * np_emp(2) ) THEN
         color = 1
      ELSE
         color = 0
      END IF
      CALL mp_comm_split( intra_image_comm, color, me_image, emp_comm )
      
      if( me_image <  np_emp(1) * np_emp(2) ) then
          me_emp(1) = me_image / np_emp(1)
          me_emp(2) = MOD( me_image, np_emp(1) ) 
      else
          me_emp(1) = me_image
          me_emp(2) = me_image
      endif
      !
      first = .FALSE.
      !
      !  Done with the group
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
      DO iss = 1, nspin
         CALL descla_init( desc_emp( :, iss ), nupdwn_emp( iss ), n_empx, np_emp, me_emp, emp_comm, color )
      END DO
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
         issw = iupdwn( iss )
         CALL nlsm1 ( nupdwn( iss ), 1, nvb, eigr, c0( 1, issw ), bec_occ( 1, iupdwn( iss ) ) )
      END DO
      !
      IF( .NOT. exst ) THEN
         !
         ! ...  initial random states orthogonal to filled ones
         !
         CALL wave_rand_init( c0_emp, n_empx * nspin, 1 )
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
            issw   = iupdwn( iss )
            !
            CALL gram_empty( .false. , eigr, vkb, bec_emp( :, in_emp: ), bec_occ( :, in: ), nkb, &
                             c0_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss) )
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
            CALL dforce( i, bec_emp, vkb, c0_emp, c2, c3, v, SIZE(v,1), ispin_emp, f_emp, n_emps, nspin )
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
            issw   = iupdwn( iss )
            !
            CALL gram_empty( .true. , eigr, vkb, bec_emp( :, in_emp: ), bec_occ( :, in: ), nkb, &
                        cm_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss) )
            !
         END DO

         ! ... calphi calculates phi
         ! ... the electron mass rises with g**2
         !
         CALL calphi( c0_emp, ngw, bec_emp, nkb, vkb, phi_emp, n_emps, ema0bg )
         !
         CALL ortho( eigr, cm_emp, phi_emp, ngw, lambda_emp, desc_emp, &
                        bigr, iter_ortho, ccc, bephi_emp, becp_emp, n_emps, nspin, &
                        nupdwn_emp, iupdwn_emp )
         !
         DO iss = 1, nspin
            !
            CALL updatc( ccc, n_emps, lambda_emp(:,:,iss), n_empx, phi_emp, ngw, &
                         bephi_emp, nkb, becp_emp, bec_emp, &
                         cm_emp, nupdwn_emp(iss), iupdwn_emp(iss), desc_emp(:,iss) )
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

      ! ...  Compute eigenvalues and bring wave functions on Kohn-Sham orbitals

      ALLOCATE( lambda_rep( n_empx, n_empx ) )
      !
      DO iss = 1, nspin
         i = iupdwn_emp(iss)
         n = nupdwn_emp(iss)
         CALL collect_lambda( lambda_rep, lambda_emp(:,:,iss), desc_emp( :, iss ) )
         CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep, n_empx, ei_emp(:,iss) )
         ei_emp( 1:n, iss ) = ei_emp( 1:n, iss ) / f_emp( i : i + n - 1 )
      END DO
      !
      DEALLOCATE( lambda_rep )

      ! ...   Save emptystates to disk

      CALL writeempty( cm_emp, n_empx * nspin )
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
   END SUBROUTINE empty_cp_x

