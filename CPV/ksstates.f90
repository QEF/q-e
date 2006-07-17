!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


MODULE kohn_sham_states


   USE io_files, ONLY: ksunit, ks_file, ks_emp_file

   IMPLICIT NONE
   SAVE

   PRIVATE

   ! ...   print KS states to file KS.indx_ksout if ksout true
   LOGICAL :: tksout                        
   LOGICAL :: tksout_emp 

   INTEGER, ALLOCATABLE :: indx_ksout(:,:)  ! (state inds, spin indxs)
   INTEGER, ALLOCATABLE :: n_ksout(:)       ! (spin indxs)
   INTEGER, ALLOCATABLE :: indx_ksout_emp(:,:)  ! (state inds, spin indxs)
   INTEGER, ALLOCATABLE :: n_ksout_emp(:)       ! (spin indxs)

   PUBLIC :: ks_states_init, ks_states_closeup
   PUBLIC :: n_ksout, indx_ksout, ks_states, tksout

!  ----------------------------------------------
CONTAINS
!  ----------------------------------------------


   SUBROUTINE ks_states_init( nspin, nprnks, iprnks, nprnks_emp, iprnks_emp )

      INTEGER, INTENT(IN) :: nspin, nprnks(:), nprnks_emp(:)
      INTEGER, INTENT(IN) :: iprnks(:,:)
      INTEGER, INTENT(IN) :: iprnks_emp(:,:)

      INTEGER :: i, ip, k, nstates

      ! ...   Tell the code which Kohn-Sham state should be printed to file
      !
      IF( ALLOCATED( n_ksout    ) ) DEALLOCATE( n_ksout )
      IF( ALLOCATED( indx_ksout ) ) DEALLOCATE( indx_ksout )
      !
      tksout = ANY( nprnks > 0 )
      !
      IF( tksout ) THEN
         nstates = MAXVAL( nprnks )
         ALLOCATE( n_ksout( nspin ) )
         ALLOCATE( indx_ksout( nstates, nspin) )
         n_ksout( 1:nspin ) = nprnks( 1:nspin )
         DO i = 1, nspin
           DO k = 1, nprnks( i )
              indx_ksout( k, i ) = iprnks( k, i )
           END DO
         END DO
      END IF

      IF( ALLOCATED( n_ksout_emp    ) ) DEALLOCATE( n_ksout_emp )
      IF( ALLOCATED( indx_ksout_emp ) ) DEALLOCATE( indx_ksout_emp )
      !
      tksout_emp = ANY( nprnks_emp > 0 )
      !
      IF( tksout_emp ) THEN
         nstates = MAXVAL( nprnks_emp )
         ALLOCATE( n_ksout_emp( nspin ) )
         ALLOCATE( indx_ksout_emp( nstates, nspin ) )
         n_ksout_emp( 1:nspin ) = nprnks_emp( 1:nspin )
         DO i = 1, nspin
            DO k = 1, n_ksout_emp( i )
               indx_ksout_emp( k, i ) = iprnks_emp( k, i )
            END DO
         END DO
      END IF

      RETURN
   END SUBROUTINE ks_states_init

!  ----------------------------------------------

   SUBROUTINE ks_states_closeup()
      IF( ALLOCATED( indx_ksout ) ) DEALLOCATE( indx_ksout )
      IF( ALLOCATED( n_ksout ) ) DEALLOCATE( n_ksout )
      tksout = .FALSE.
      IF( ALLOCATED( indx_ksout_emp ) ) DEALLOCATE( indx_ksout_emp )
      IF( ALLOCATED( n_ksout_emp ) ) DEALLOCATE( n_ksout_emp )
      tksout_emp = .FALSE.
      RETURN
   END SUBROUTINE ks_states_closeup

!  ----------------------------------------------


      SUBROUTINE ks_states(cf, wfill, occ, vpot, eigr, vkb, bec )

        ! ...   declare modules
        USE kinds
        USE mp_global,        ONLY : intra_image_comm
        USE io_global,        ONLY : ionode, stdout
        USE ions_base,        ONLY : nsp
        USE wave_types,       ONLY : wave_descriptor, wave_descriptor_init
        USE wave_functions,   ONLY : kohn_sham
        USE forces
        USE control_flags,    ONLY : force_pairing
        USE electrons_base,   ONLY : nupdwn, iupdwn, nspin, nbsp
        USE electrons_module, ONLY : n_emp, nupdwn_emp, iupdwn_emp, nb_l, n_emp_l
        USE empty_states,     ONLY : readempty
        USE uspp,             ONLY : nkb

        IMPLICIT NONE

        ! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) :: cf(:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wfill
        COMPLEX(DP)  ::  eigr(:,:)
        COMPLEX(DP)  ::  vkb(:,:)
        REAL(DP), INTENT(IN)  ::  occ(:,:), bec(:,:)
        REAL(DP) ::  vpot(:,:)

        ! ...   declare other variables
        INTEGER  :: i, ib, ig, ngw, nb_g, iss, iks
        INTEGER  :: n_emps, in, ns
        LOGICAL  :: tortho = .TRUE.
        LOGICAL  :: exist 

        COMPLEX(DP), ALLOCATABLE :: fforce(:,:)
        COMPLEX(DP), ALLOCATABLE :: eforce(:,:)
        COMPLEX(DP), ALLOCATABLE :: ce(:,:)
        REAL(DP), ALLOCATABLE :: bec_emp(:,:)
        REAL(DP), ALLOCATABLE :: fi(:)

        CHARACTER (LEN=6), EXTERNAL :: int_to_char

        TYPE (wave_descriptor) :: wempt  ! wave function descriptor for empty states

        ! ...   end of declarations

        ngw  = wfill%ngwl

        IF( tksout ) THEN

           ALLOCATE( fforce( ngw, nbsp ) )

           DO iss = 1, nspin

              IF( nupdwn( iss ) > 0 ) THEN

                 CALL dforce_all( cf, occ(:,iss), fforce, vpot(:,iss), vkb, bec, nupdwn(iss), iupdwn(iss) )

              END IF

           END DO

           IF( force_pairing ) THEN 
              DO i = 1, nupdwn(2)
                 fforce(:,i) = occ(i,1) * fforce(:,i) + occ(i,2) * fforce(:,i+iupdwn(2)-1)
              END DO
              DO i = nupdwn(2)+1, nupdwn(1)
                 fforce(:,i) = occ(i,1) * fforce(:,i)
              END DO
              DO i = iupdwn(2), iupdwn(2) + nupdwn(2) - 1
                 fforce(:,i) = fforce(:,i-iupdwn(2)+1)
              END DO
           END IF

           DO iss = 1, nspin
              CALL kohn_sham( cf, ngw, fforce, nupdwn( iss ), nb_l(iss), iupdwn( iss ) )
           END DO

           DEALLOCATE( fforce )

        END IF


        IF( tksout_emp ) THEN

          n_emps = nupdwn_emp( 1 )
          IF( nspin == 2 ) n_emps = n_emps + nupdwn_emp( 2 )

          CALL wave_descriptor_init( wempt, ngw, wfill%ngwt, nupdwn_emp, nupdwn_emp, &
                                 1, 1, nspin, 'gamma' , wfill%gzero )
  
          ALLOCATE( ce( ngw,  n_emp * nspin ) )
          !
          ALLOCATE( bec_emp( nkb, n_emps ) )

          exist = readempty( ce, n_emp * nspin )

          IF( .NOT. exist ) &
             CALL errore( ' ks_states ', ' empty states file not found', 1 )

          CALL nlsm1 ( n_emps, 1, nsp, eigr, ce, bec_emp )

          ALLOCATE( eforce( ngw, SIZE( ce, 2 ) ) )

          DO iss = 1, nspin

            in = iupdwn_emp( iss )
            ns = nupdwn_emp( iss ) 

            IF( ns > 0 ) THEN

              ALLOCATE( fi( ns ) )

              fi = 2.0d0 / nspin

              CALL dforce_all( ce, fi, eforce, vpot(:,iss), vkb, bec_emp, ns, in ) 

              CALL kohn_sham( ce, ngw, eforce, ns, n_emp_l(iss), in )

              DEALLOCATE( fi )

            END IF

          END DO

          DEALLOCATE( eforce )

          DEALLOCATE( bec_emp )

        END IF

        CALL print_all_states( cf, ce )

        IF( tksout_emp ) THEN
          DEALLOCATE( ce )
        END IF

        RETURN
        !
      END SUBROUTINE ks_states

!  ----------------------------------------------

      SUBROUTINE print_all_states( cf, ce )

        USE kinds,            ONLY : DP
        USE mp_global,        ONLY : intra_image_comm
        USE io_global,        ONLY : ionode
        USE io_global,        ONLY : stdout
        USE electrons_module, ONLY : iupdwn_emp, nupdwn_emp
        USE electrons_base,   ONLY : nupdwn, iupdwn, nspin

        IMPLICIT NONE

        ! ...   declare subroutine arguments
        COMPLEX(DP),            INTENT(INOUT) :: cf(:,:), ce(:,:)

        ! ...   declare other variables
        INTEGER ::  i, iss, iks

        CHARACTER(LEN=256) :: file_name
        CHARACTER(LEN=10), DIMENSION(2) :: spin_name
        CHARACTER (LEN=6), EXTERNAL :: int_to_char

        IF( tksout .OR. tksout_emp ) THEN

          IF (ionode) THEN
            WRITE( stdout,*) 
            WRITE( stdout,'( "   Khon Sham state")') 
            WRITE( stdout,'( "   ---------------")') 
          END IF

          IF( nspin == 2 ) THEN
            spin_name(1) = '_UP_'
            spin_name(2) = '_DW_'
          ELSE
            spin_name(1) = '_'
            spin_name(2) = '_'
          END IF

          DO iss = 1, nspin
            IF( tksout ) THEN
              DO i = 1, n_ksout(iss)
                iks = indx_ksout(i, iss)
                IF( ( iks > 0 ) .AND. ( iks <= nupdwn( iss ) ) ) THEN
                  file_name = TRIM( ks_file ) // &
                            & trim(spin_name(iss)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( cf(:,iks+iupdwn(iss)-1), file_name )
                END IF
              END DO
            END IF
            IF( tksout_emp ) THEN
              DO i = 1, n_ksout_emp(iss)
                iks = indx_ksout_emp(i, iss)
                IF( ( iks > 0 ) .AND. ( iks <= nupdwn_emp( iss ) ) ) THEN
                  file_name = TRIM( ks_emp_file ) // &
                            & trim(spin_name(iss)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( ce(:,iupdwn_emp(iss)+iks-1), file_name )
                END IF
              END DO
            END IF
          END DO

        END IF

        RETURN
        ! ...
      END SUBROUTINE print_all_states


!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE print_ks_states( c, file_name )

        USE kinds
        USE mp, ONLY: mp_sum
        USE io_global, ONLY: ionode, ionode_id
        USE io_global, ONLY: stdout
        USE gvecw, ONLY: ngw
        USE fft_base, ONLY: dfftp, dffts
        USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x, nnrx
        USE fft_module, ONLY: invfft
        USE xml_io_base, ONLY: write_rho_xml
        USE mp_global,       ONLY: nproc_image, me_image, intra_image_comm

        IMPLICIT NONE

        COMPLEX(DP),      INTENT(IN) :: c(:)
        CHARACTER(LEN=*), INTENT(IN) :: file_name
        REAL(DP),    ALLOCATABLE :: rpsi2(:)
        COMPLEX(DP), ALLOCATABLE :: psi(:)
        INTEGER   ::  i
        REAL(DP) :: charge

        ALLOCATE( psi( nnrx ) )
        ALLOCATE( rpsi2( nnrx ) )

        CALL c2psi( psi, dffts%nnr, c, c, ngw, 1 )
        CALL invfft( 'Wave', psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

        DO i = 1, nnrx
           rpsi2( i ) = DBLE( psi( i ) )**2
        END DO
        charge = SUM( rpsi2 )

        CALL write_rho_xml( file_name, rpsi2, &
                            nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
        
        CALL mp_sum( charge, intra_image_comm )

        IF ( ionode ) THEN
          WRITE( stdout,'(3X,A15," integrated charge : ",F14.5)')  &
     &      TRIM(file_name), charge / DBLE(nr1*nr2*nr3)
        END IF

        DEALLOCATE( rpsi2, psi )
        ! ...
        RETURN
        ! ...
      END SUBROUTINE print_ks_states

!  ----------------------------------------------
!
END MODULE kohn_sham_states
