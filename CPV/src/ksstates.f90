!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


MODULE kohn_sham_states


   IMPLICIT NONE
   SAVE

   PRIVATE

   ! ...   print KS states to file KS.indx_ksout if ksout true
   LOGICAL :: tksout                        
   CHARACTER(LEN=2 ), PARAMETER :: ks_file       = 'KS'

   INTEGER, ALLOCATABLE :: indx_ksout(:,:)  ! (state inds, spin indxs)
   INTEGER, ALLOCATABLE :: n_ksout(:)       ! (spin indxs)

   PUBLIC :: ks_states_init, ks_states_closeup
   PUBLIC :: n_ksout, indx_ksout, tksout, print_all_states

!  ----------------------------------------------
CONTAINS
!  ----------------------------------------------


   SUBROUTINE ks_states_init( nspin, nprnks, iprnks )

      INTEGER, INTENT(IN) :: nspin, nprnks(:)
      INTEGER, INTENT(IN) :: iprnks(:,:)

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

      RETURN
   END SUBROUTINE ks_states_init

!  ----------------------------------------------

   SUBROUTINE ks_states_closeup()
      IF( ALLOCATED( indx_ksout ) ) DEALLOCATE( indx_ksout )
      IF( ALLOCATED( n_ksout ) ) DEALLOCATE( n_ksout )
      tksout = .FALSE.
      RETURN
   END SUBROUTINE ks_states_closeup

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE print_all_states( ctot, iupdwn_tot, nupdwn_tot )

        USE kinds,            ONLY : DP
        USE mp_bands,         ONLY : intra_bgrp_comm
        USE io_global,        ONLY : ionode, stdout
        USE electrons_base,   ONLY : nupdwn, iupdwn, nspin

        IMPLICIT NONE

        ! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(IN) :: ctot(:,:)
        INTEGER,     INTENT(IN) :: iupdwn_tot(2)
        INTEGER,     INTENT(IN) :: nupdwn_tot(2)

        ! ...   declare other variables
        INTEGER ::  i, iss, iks, itot

        CHARACTER(LEN=256) :: file_name
        CHARACTER(LEN=10), DIMENSION(2) :: spin_name
        CHARACTER (LEN=6), EXTERNAL :: int_to_char

        CALL errore('print_all_states','unused since a long time, please check',1)
        IF( tksout ) THEN

          IF (ionode) THEN
            WRITE( stdout,*) 
            WRITE( stdout,'( "   Kohn Sham state")') 
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
                  itot = iks + iupdwn_tot(iss) - 1 
                  file_name = TRIM( ks_file ) // &
                            & trim(spin_name(iss)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( ctot( :, itot ), file_name )
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
        USE io_global, ONLY: ionode, ionode_id, stdout
        USE fft_base, ONLY: dfftp, dffts
        USE fft_interfaces, ONLY: invfft
        USE fft_helper_subroutines, ONLY: c2psi_gamma
        USE fft_rho, ONLY: rho_r2g
        USE mp_bands, ONLY: intra_bgrp_comm, inter_bgrp_comm, my_bgrp_id,&
             root_bgrp_id, root_bgrp
        USE gvect, ONLY: ig_l2g, mill, ecutrho
        USE io_base, ONLY: write_rhog
 
        IMPLICIT NONE

        COMPLEX(DP),      INTENT(IN) :: c(:)
        CHARACTER(LEN=*), INTENT(IN) :: file_name
        COMPLEX(DP), ALLOCATABLE :: psi(:), rhog(:,:)
        REAL(DP), ALLOCATABLE :: rhor(:,:)
        INTEGER   ::  i
        ! FIXME: reciprocal lattice vectors
        REAL(DP) :: bogus1(3), bogus2(3), bogus3(3)
        REAL(DP) :: charge

        ALLOCATE( rhor(dfftp%nnr,1), psi(dfftp%nnr) )

        CALL c2psi_gamma( dffts, psi, c )
        CALL invfft( 'Wave', psi, dffts )

        ! FIXME: not sure things will work in presence of a double grid
        DO i = 1, dfftp%nnr
           rhor( i,1 ) = DBLE( psi( i ) )**2
        END DO
        charge = SUM( rhor )
        CALL mp_sum( charge, intra_bgrp_comm )

        DEALLOCATE( psi )
        ALLOCATE( rhog(dfftp%ngm,1) )
        CALL rho_r2g( dfftp, rhor, rhog )
        DEALLOCATE ( rhor )

        ! Only the first band group collects and writes
        IF ( my_bgrp_id == root_bgrp_id ) CALL write_rhog &
                ( file_name, root_bgrp, intra_bgrp_comm, &
                bogus1, bogus2, bogus3, .true., &
                mill, ig_l2g, rhog, ecutrho )

        IF ( ionode ) THEN
          WRITE( stdout,'(3X,A15," integrated charge : ",F14.5)')  &
     &      TRIM(file_name), charge / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
        END IF

        DEALLOCATE( rhog )
        ! ...
        RETURN
        ! ...
      END SUBROUTINE print_ks_states

!  ----------------------------------------------
!
END MODULE kohn_sham_states
