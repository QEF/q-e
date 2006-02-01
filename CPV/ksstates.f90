!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!  BEGIN manual

      MODULE kohn_sham_states

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  ----------------------------------------------
!  END manual

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

        PUBLIC :: ks_states_init, kohn_sham, ks_states_closeup
        PUBLIC :: n_ksout, indx_ksout, ks_states, tksout
        PUBLIC :: ks_states_force_pairing

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE ks_states_init( nspin, tprnks, tprnks_emp )

!  (describe briefly what this routine does...)
!  ----------------------------------------------
!  END manual

        INTEGER, INTENT(IN) :: nspin
        LOGICAL, INTENT(IN) :: tprnks(:,:)
        LOGICAL, INTENT(IN) :: tprnks_emp(:,:)

        INTEGER :: i, ip, k, nstates


! ...   Tell the code which Kohn-Sham state should be printed to file
        !
        tksout = .FALSE.
        IF( ALLOCATED( n_ksout    ) ) DEALLOCATE( n_ksout )
        IF( ALLOCATED( indx_ksout ) ) DEALLOCATE( indx_ksout )
        !
        IF( ANY( tprnks ) ) THEN
          tksout = .TRUE.
          ALLOCATE( n_ksout( nspin ) )
          DO i = 1, nspin
            n_ksout(i) = COUNT( tprnks(:,i) )
          END DO
          nstates = MAXVAL(n_ksout)
          ALLOCATE( indx_ksout( nstates, nspin) )
          indx_ksout = 0
          DO i = 1, nspin
            ip = 1
            DO k = 1, SIZE( tprnks, 1)
              IF( tprnks(k,i) ) THEN
                indx_ksout( ip, i ) = k
                ip = ip + 1
              END IF
            END DO
          END DO
        END IF

        tksout_emp = .FALSE.
        IF( ALLOCATED( n_ksout_emp    ) ) DEALLOCATE( n_ksout_emp )
        IF( ALLOCATED( indx_ksout_emp ) ) DEALLOCATE( indx_ksout_emp )
        !
        IF( ANY( tprnks_emp ) ) THEN
          tksout_emp = .TRUE.
          ALLOCATE( n_ksout_emp( nspin ) )
          DO i = 1, nspin
            n_ksout_emp(i) = COUNT( tprnks_emp(:,i) )
          END DO
          nstates = MAXVAL(n_ksout_emp)
          ALLOCATE( indx_ksout_emp( nstates, nspin) )
          indx_ksout_emp = 0
          DO i = 1, nspin
            ip = 1
            DO k = 1, SIZE( tprnks_emp, 1)
              IF( tprnks_emp(k,i) ) THEN
                indx_ksout_emp( ip, i ) = k
                ip = ip + 1
              END IF
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
!  BEGIN manual

      SUBROUTINE kohn_sham(ispin, c, cdesc, eforces )

!  (describe briefly what this routine does...)
!  ----------------------------------------------
!  END manual

! ...   declare modules
        USE kinds
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp, ONLY: mp_sum
        USE io_global, ONLY: ionode
        USE wave_base, ONLY: hpsi
        USE wave_functions, ONLY: crot
        USE wave_constrains, ONLY: update_lambda
        USE wave_types, ONLY: wave_descriptor
        USE electrons_module, ONLY: nb_l
        USE forces
        USE brillouin, ONLY: kpoints, kp

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) ::  c(:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: cdesc
        INTEGER, INTENT(IN) :: ispin

! ...   descriptor of the electronic state distribution
                                       
        COMPLEX(DP) :: eforces(:,:,:)

! ...   declare other variables
        INTEGER ::  ik, ib, nk, ngw, nb_g, nrl, ibl
        COMPLEX(DP), ALLOCATABLE :: cgam(:,:)
        REAL(DP),    ALLOCATABLE :: gam(:,:)
        REAL(DP),    ALLOCATABLE :: eig(:)
        LOGICAL :: tortho = .TRUE.

! ...   end of declarations
!  ----------------------------------------------

          nb_g = cdesc%nbt( ispin )

          IF( nb_g < 1 ) THEN
            
            eforces = 0.0d0

          ELSE

            ! ...     Distribute electronic states cyclically across processors
       
            nrl = nb_l( ispin )

            ALLOCATE( eig(nb_g) )
            IF( kp%gamma_only ) THEN
              ALLOCATE( cgam(1,1), gam(nrl,nb_g) )
            ELSE
              ALLOCATE(cgam(nrl,nb_g), gam(1,1))
            END IF
            DO ik = 1, kp%nkpt
              IF( cdesc%gamma ) THEN
                DO ib = 1, nb_g
                  CALL update_lambda( ib, gam, c(:,:,ik), cdesc, eforces(:,ib,ik) )
                END DO
                CALL crot( ispin, c(:,:,ik), cdesc, gam, eig)
              ELSE
                DO ib = 1, nb_g
                  CALL update_lambda( ib, cgam, c(:,:,ik), cdesc, eforces(:,ib,ik) )
                END DO
                CALL crot( ispin, ik, c(:,:,:), cdesc, cgam, eig)
              END IF
            END DO
            DEALLOCATE(cgam, gam, eig)

          END IF

        RETURN
! ...
      END SUBROUTINE kohn_sham

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE ks_states(cf, wfill, ce, wempt, occ, vpot, eigr, bec )

!  (describe briefly what this routine does...)
!  ----------------------------------------------
!  END manual

! ...   declare modules
        USE kinds
        USE mp_global, ONLY: mpime, nproc, group, root
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE wave_types, ONLY: wave_descriptor
        USE forces
        USE brillouin, ONLY: kpoints, kp
        USE pseudo_projector, ONLY: projector
        USE control_flags, ONLY: timing, force_pairing

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) :: cf(:,:,:,:), ce(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wfill, wempt
        COMPLEX(DP)  ::  eigr(:,:)
        REAL(DP), INTENT(IN)  ::  occ(:,:,:), bec(:,:)
        REAL (DP) ::  vpot(:,:)

! ...   declare other variables
        INTEGER ::  i, ik, ib, nk, ig, ngw, nb_g, nb_l, ispin, nspin, iks
        INTEGER ::  ispin_wfc
        LOGICAL  :: tortho = .TRUE.
        CHARACTER(LEN=4) :: nom
        CHARACTER(LEN=256) :: file_name
        CHARACTER(LEN=10), DIMENSION(2) :: spin_name
        REAL(DP) :: s0, s1, s2, s3

        COMPLEX(DP), ALLOCATABLE :: eforce(:,:,:)
        REAL(DP), ALLOCATABLE :: fi(:,:)

        REAL(DP), EXTERNAL :: cclock
        CHARACTER (LEN=6), EXTERNAL :: int_to_char


! ...   end of declarations
!  ----------------------------------------------

        nk    = wfill%nkl
        nspin = wfill%nspin

        IF( .NOT. wfill%gamma ) &
          CALL errore( ' ks_states ', ' only gamma is implemented ', 1 )

        IF( timing ) s0 = cclock()

        DO ispin = 1, nspin

          ngw  = wfill%ngwl
          nb_l = wfill%nbl( ispin )

          ispin_wfc = ispin
          IF( force_pairing ) ispin_wfc = 1

          IF( nb_l > 0 ) THEN

            ALLOCATE(  eforce( ngw,  nb_l, nk ))

            CALL dforce_all( ispin, cf(:,:,1,ispin_wfc), wfill, occ(:,1,ispin), eforce(:,:,1), &
              vpot(:,ispin), eigr, bec )

            CALL kohn_sham( ispin, cf(:,:,:,ispin_wfc), wfill, eforce )

            DEALLOCATE( eforce )

          END IF

          IF( tksout_emp ) THEN

            ngw  = wempt%ngwl
            nb_l = wempt%nbl( ispin )

            IF( nb_l > 0 ) THEN

              ALLOCATE( fi( nb_l, nk ) )
              DO ik = 1, nk
                fi( 1:nb_l, ik ) = 2.0d0 / nspin
              END DO

              ALLOCATE(  eforce( ngw,  nb_l, nk ))

              CALL dforce_all( ispin, ce(:,:,1,ispin), wempt, fi(:,1), eforce(:,:,1), &
                               vpot(:,ispin), eigr, bec )

              CALL kohn_sham( ispin, ce(:,:,:,ispin), wempt, eforce )

              DEALLOCATE( eforce )
              DEALLOCATE( fi )

            END IF

          END IF

        END DO

        IF( timing ) s1 = cclock()

        IF( tksout .OR. tksout_emp ) THEN
          IF (ionode) THEN
            WRITE( stdout,*) 
            WRITE( stdout,'( "   Khon Sham state")') 
            WRITE( stdout,'( "   ---------------")') 
          END IF

          IF(nspin .EQ. 2) THEN
            spin_name(1) = '_UP.'
            spin_name(2) = '_DW.'
          ELSE
            spin_name(1) = '.'
            spin_name(2) = '.'
          END IF

          DO ispin = 1, nspin
            ispin_wfc = ispin
            IF( force_pairing ) ispin_wfc = 1
            IF( tksout ) THEN
              DO i = 1, n_ksout(ispin)
                iks = indx_ksout(i, ispin)
                IF( ( iks > 0 ) .AND. ( iks <= wfill%nbt( ispin ) ) ) THEN
                  file_name = TRIM( ks_file ) // &
                            & trim(spin_name(ispin)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( cf(:,iks,1,ispin_wfc), file_name )
                END IF
              END DO
            END IF
            IF( tksout_emp ) THEN
              DO i = 1, n_ksout_emp(ispin)
                iks = indx_ksout_emp(i, ispin)
                IF( ( iks > 0 ) .AND. ( iks <= wempt%nbt( ispin ) ) ) THEN
                  file_name = TRIM( ks_emp_file ) // &
                            & trim(spin_name(ispin)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( ce(:,iks,1,ispin), file_name )
                END IF
              END DO
            END IF
          END DO

        END IF

        IF( timing ) THEN
          s2 = cclock()
          IF( ionode ) THEN
            WRITE( stdout,fmt="(3X,'time for KS ortho     = ',F8.2)") (s1-s0)
            WRITE( stdout,fmt="(3X,'time for KS print out = ',F8.2)") (s2-s1)
          END IF
        END IF
! ...
        RETURN
! ...
      END SUBROUTINE ks_states

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE ks_states_force_pairing(cf, wfill, ce, wempt, occ, vpot, eigr, bec )

!  (describe briefly what this routine does...)
!  ----------------------------------------------
!  END manual

! ...   declare modules
        USE kinds
        USE mp_global, ONLY: mpime, nproc, group, root
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE wave_types, ONLY: wave_descriptor
        USE forces
        USE brillouin, ONLY: kpoints, kp
        USE pseudo_projector, ONLY: projector
        USE control_flags, ONLY: timing
        USE electrons_module, ONLY: nupdwn, nspin

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) :: cf(:,:,:,:), ce(:,:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: wfill, wempt
        COMPLEX(DP)  ::  eigr(:,:)
        REAL(DP), INTENT(IN)  ::  occ(:,:,:), bec(:,:)
        REAL (DP) ::  vpot(:,:)

! ...   declare other variables
        INTEGER ::  i, ik, ib, nk, ig, ngw, nb_g, nb_l, iks, nb, ispin
        LOGICAL  :: tortho = .TRUE.
        CHARACTER(LEN=4) :: nom
        CHARACTER(LEN=256) :: file_name
        CHARACTER(LEN=10), DIMENSION(2) :: spin_name
        REAL(DP) :: s0, s1, s2, s3

        COMPLEX(DP), ALLOCATABLE :: eforce(:,:,:,:)
        REAL(DP), ALLOCATABLE :: fi(:,:)

        CHARACTER (LEN=6), EXTERNAL :: int_to_char
        REAL(DP), EXTERNAL :: cclock


! ...   end of declarations
!  ----------------------------------------------

        nk    = wfill%nkl

        IF( .NOT. wfill%gamma ) &
          CALL errore( ' ks_states_force_pairing ', ' only gamma is implemented ', 1 )

        IF( nspin == 1 ) &
          CALL errore(' ks_states_forced_pairing ',' inconsistent nspin ', 1)

        IF( nupdwn(1) < nupdwn(2) ) &
          CALL errore(' ks_states_forced_pairing ',' inconsistent nupdwn ', 1)

        IF( timing ) s0 = cclock()

        ngw  = wfill%ngwl
        nb   = nupdwn(1)

        IF( nb > 0 ) THEN

          ALLOCATE(  eforce( ngw,  nb, 1, 2 ) )

          CALL dforce_all( 1, cf(:,:,1,1), wfill, occ(:,1,1), eforce(:,:,1,1), &
              vpot(:,1), eigr, bec )
          CALL dforce_all( 2, cf(:,:,1,1), wfill, occ(:,1,2), eforce(:,:,1,2), &
              vpot(:,2), eigr, bec )

          DO i = 1, nupdwn(2)
            eforce(:,i,1,1) = occ(i,1,1) * eforce(:,i,1,1) + occ(i,1,2) * eforce(:,i,1,2)
          END DO
          DO i = nupdwn(2)+1, nupdwn(1)
            eforce(:,i,1,1) = occ(i,1,1) * eforce(:,i,1,1)
          END DO

          CALL kohn_sham( 1, cf(:,:,:,1), wfill, eforce(:,:,:,1) )

          DEALLOCATE( eforce )

        END IF

        IF( tksout_emp ) THEN

          ngw  = wempt%ngwl
          nb_l = wempt%nbl( 1 )

          IF( nb_l > 0 ) THEN

            ALLOCATE( fi( nb_l, nk ) )
            DO ik = 1, nk
              fi( 1:nb_l, ik ) = 2.0d0
            END DO

            ALLOCATE(  eforce( ngw,  nb_l, 1, 1 ))

            CALL dforce_all( 1, ce(:,:,1,1), wempt, fi(:,1), eforce(:,:,1,1), vpot(:,1), &
                             eigr, bec )

            CALL kohn_sham( 1, ce(:,:,:,1), wempt, eforce(:,:,:,1) )

            CALL dforce_all( 2, ce(:,:,1,2), wempt, fi(:,1), eforce(:,:,1,1), vpot(:,2), &
                             eigr, bec )

            CALL kohn_sham( 2, ce(:,:,:,2), wempt, eforce(:,:,:,1) )

            DEALLOCATE( eforce )
            DEALLOCATE( fi )

          END IF

        END IF


        IF( timing ) s1 = cclock()

        IF( tksout .OR. tksout_emp ) THEN
          IF (ionode) THEN
            WRITE( stdout,*) 
            WRITE( stdout,'( "   Khon Sham state")') 
            WRITE( stdout,'( "   ---------------")') 
          END IF

          IF(nspin .EQ. 2) THEN
            spin_name(1) = '_UP.'
            spin_name(2) = '_DW.'
          ELSE
            spin_name(1) = '.'
            spin_name(2) = '.'
          END IF

          DO ispin = 1, nspin
            IF( tksout ) THEN
              DO i = 1, n_ksout(ispin)
                iks = indx_ksout(i, ispin)
                IF( ( iks > 0 ) .AND. ( iks <= wfill%nbt( ispin ) ) ) THEN
                  file_name = TRIM( ks_file ) // &
                            & trim(spin_name(ispin)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( cf(:,iks,1,1), file_name )
                END IF
              END DO
            END IF
            IF( tksout_emp ) THEN
              DO i = 1, n_ksout_emp(ispin)
                iks = indx_ksout_emp(i, ispin)
                IF( ( iks > 0 ) .AND. ( iks <= wempt%nbt( ispin ) ) ) THEN
                  file_name = TRIM( ks_emp_file ) // &
                            & trim(spin_name(ispin)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( ce(:,iks,1,ispin), file_name )
                END IF
              END DO
            END IF
          END DO

        END IF

        IF( timing ) THEN
          s2 = cclock()
          IF( ionode ) THEN
            WRITE( stdout,fmt="(3X,'time for KS ortho     = ',F8.2)") (s1-s0)
            WRITE( stdout,fmt="(3X,'time for KS print out = ',F8.2)") (s2-s1)
          END IF
        END IF
! ...
        RETURN
! ...
      END SUBROUTINE ks_states_force_pairing


!  ----------------------------------------------

      SUBROUTINE print_ks_states( psi, file_name )

        USE kinds
        USE fft, ONLY: pw_invfft
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp, ONLY: mp_barrier, mp_sum
        USE io_global, ONLY: ionode, ionode_id
        USE io_global, ONLY: stdout
        USE fft_base, ONLY: dfftp
        USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x, nnrx

        IMPLICIT NONE

        COMPLEX(DP), INTENT(IN) :: psi(:)
        CHARACTER(LEN=*) :: file_name
        COMPLEX(DP), ALLOCATABLE :: zcomp(:)
        REAL(DP), ALLOCATABLE :: rcomp2(:)
        COMPLEX(DP), ALLOCATABLE :: psi2(:)
        INTEGER   ::  i, j, k, istr, izl
        REAL(DP) :: charge
        LOGICAL   :: top

        ALLOCATE( zcomp( nr3 ), rcomp2( nr3 ) )
        ALLOCATE( psi2( nnrx ) )

        CALL pw_invfft( psi2, psi, psi )
        
        INQUIRE( UNIT=ksunit, OPENED=top )
        IF( top ) THEN
          WRITE( stdout,fmt="('** WARNING: Print_ks_states, ksunit already OPENED')")
        END IF
       
        istr = INDEX( file_name, ' ' )
        IF(ionode) THEN
          OPEN(UNIT=ksunit, FILE=file_name(1:istr), STATUS='unknown')
        END IF

        charge = 0.0d0

        izl = dfftp%ipp( mpime + 1 )

        DO i = 1, nr1

          DO j = 1, nr2

            zcomp = 0.0d0

            istr = i + nr1 * ( j - 1 )

            DO k = 1, dfftp%npl 
               zcomp( izl + k ) = psi2( istr + nr1 * nr2 * ( k - 1 ) )
            END DO

            CALL mp_sum( zcomp( 1 : nr3 ) )

            IF ( ionode ) THEN
              rcomp2 = DBLE(zcomp)**2
              WRITE(ksunit, fmt='(F10.5)') ( rcomp2(k), k=1, nr3 )
              charge = charge + SUM(rcomp2)
            END IF

            CALL mp_barrier()

          END DO

        END DO

        IF ( ionode ) THEN
          CLOSE(ksunit)
          WRITE( stdout,'(3X,A15," integrated charge : ",F14.5)')  &
     &      file_name(1:istr), charge / DBLE(nr1*nr2*nr3)
        END IF
        DEALLOCATE(zcomp, rcomp2, psi2)
! ...
        RETURN
! ...
      END SUBROUTINE print_ks_states

!  ----------------------------------------------

      END MODULE kohn_sham_states

