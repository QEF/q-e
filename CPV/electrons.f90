!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE electrons_module
!=----------------------------------------------------------------------------=!
#include "f_defs.h"
        USE kinds
        USE parameters,         ONLY: nspinx
        USE dspev_module,       ONLY: pdspev_drv, dspev_drv
        USE electrons_base,     ONLY: nbnd, nbndx, nbsp, nbspx, nspin, nel, nelt, &
                                      nupdwn, iupdwn, telectrons_base_initval, f, &
                                      nudx
        USE cp_electronic_mass, ONLY: ecutmass => emass_cutoff
        USE cp_electronic_mass, ONLY: emass
        USE cp_electronic_mass, ONLY: emass_precond


        IMPLICIT NONE
        SAVE

        PRIVATE 

!  ...  declare module-scope variables

        LOGICAL :: band_first = .TRUE.

        INTEGER :: n_emp               =  0  ! number of empty states
        INTEGER :: nupdwn_emp(nspinx)  =  0  ! number of empty states
        INTEGER :: iupdwn_emp(nspinx)  =  0  ! number of empty states

        INTEGER :: nb_l(nspinx)    =  0  ! local number of states ( for each spin components )
        INTEGER :: n_emp_l(nspinx) =  0
        !
        INTEGER  :: max_emp = 0    !  maximum number of iterations for empty states
        REAL(DP) :: ethr_emp       !  threshold for convergence
        !
        INTEGER, ALLOCATABLE :: ib_owner(:)
        INTEGER, ALLOCATABLE :: ib_local(:)

        REAL(DP), ALLOCATABLE :: ei(:,:)
        REAL(DP), ALLOCATABLE :: ei_emp(:,:)

!  ...  Fourier acceleration

        LOGICAL :: toccrd = .FALSE.  ! read occupation number from standard input

        PUBLIC :: electrons_setup
        PUBLIC :: bmeshset, occn_info
        PUBLIC :: deallocate_electrons
        PUBLIC :: n_emp, ei_emp, n_emp_l, ib_owner, ib_local, nb_l
        PUBLIC :: ei, nupdwn_emp, iupdwn_emp
        PUBLIC :: print_eigenvalues
        PUBLIC :: max_emp, ethr_emp
        PUBLIC :: empty_print_info, empty_init
 

!
!  end of module-scope declarations
!
!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!



   SUBROUTINE occn_info( occ )
     !
     !   This subroutine prints occupation numbers to stdout
     !
     USE io_global, ONLY: stdout, ionode
     !
     REAL(DP) :: occ(:)
     INTEGER  :: i, iss
     !
     IF( ionode ) THEN
       WRITE( stdout, fmt="(3X,'Occupation number from init')" )
       IF( nspin == 1 ) THEN
         WRITE( stdout, fmt = " (3X, 'nbnd = ', I5 ) " ) nbnd
         WRITE( stdout, fmt = " (3X,10F5.2)" ) ( occ( i ), i = 1, nbnd )
       ELSE
         DO iss = 1, nspin
           WRITE( stdout, fmt = " (3X,'spin = ', I3, ' nbnd = ', I5 ) " ) iss, nupdwn( iss )
           WRITE( stdout, fmt = " (3X,10F5.2)" ) ( occ( i+iupdwn(iss)-1 ), i = 1, nupdwn( iss ) )
         END DO
       END IF
     END IF
     !
     RETURN
   END SUBROUTINE occn_info

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE bmeshset

     !   This subroutine initialize the variables for the 
     !   distribution across processors of the overlap matrixes 
     !   of sizes ( nx, nx )

     USE mp_global, ONLY: me_image, nproc_image

     IMPLICIT NONE

     INTEGER :: i, ierr

     IF( band_first ) THEN
       CALL errore(' bmeshset ',' module not initialized ',0)
     END IF

     DO i = 1, nspin 
       !
       IF( i > nspinx ) CALL errore( ' bmeshset ',' spin too large ', i)
       !
       nb_l( i ) = nupdwn( i ) / nproc_image
       IF( me_image < MOD( nupdwn( i ), nproc_image ) ) nb_l( i ) = nb_l( i ) + 1
       !
       n_emp_l( i ) = nupdwn_emp( i ) / nproc_image
       IF( me_image < MOD( nupdwn_emp( i ), nproc_image ) ) n_emp_l( i ) = n_emp_l( i ) + 1
       !
     END DO

     IF( ALLOCATED( ib_owner ) ) DEALLOCATE( ib_owner )
     ALLOCATE( ib_owner( MAX( n_emp, nbndx ) ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_owner ', ierr)
     IF( ALLOCATED( ib_local ) ) DEALLOCATE( ib_local )
     ALLOCATE( ib_local( MAX( n_emp, nbndx ) ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_local ', ierr)

     !  here define the association between processors and electronic states
     !  round robin distribution is used

     ib_local =  0
     ib_owner = -1
     DO i = 1, MAX( n_emp, nbndx )
       ib_local( i ) = ( i - 1 ) / nproc_image        !  local index of the i-th band 
       ib_owner( i ) = MOD( ( i - 1 ), nproc_image )  !  owner of th i-th band
       IF( me_image <= ib_owner( i ) ) THEN
         ib_local( i ) = ib_local( i ) + 1
       END IF
     END DO

     RETURN
   END SUBROUTINE bmeshset

!  ----------------------------------------------
!
!
!
!  ----------------------------------------------


   SUBROUTINE electrons_setup( n_emp_ , emass_inp, ecutmass_inp )

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: n_emp_
     REAL(DP),  INTENT(IN) :: emass_inp, ecutmass_inp
     INTEGER :: ierr, i
 

     IF( .NOT. telectrons_base_initval ) &
       CALL errore( ' electrons_setup ', ' electrons_base not initialized ', 1 )

     n_emp = n_emp_
     !
     ! assure that the number of empty states is an even number
     !
     n_emp = n_emp + MOD( n_emp, 2 )
     !
     nupdwn_emp(1) = n_emp
     iupdwn_emp(1) = 1

     IF( nspin == 2 ) THEN
        nupdwn_emp(2) = n_emp
        iupdwn_emp(2) = 1 + n_emp
     END IF

     IF( ALLOCATED( ei ) ) DEALLOCATE( ei )
     ALLOCATE( ei( nudx, nspin ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei ',ierr)
     ei = 0.0_DP

     IF( ALLOCATED( ei_emp ) ) DEALLOCATE( ei_emp )
     IF( n_emp > 0 ) THEN
       ALLOCATE( ei_emp( n_emp, nspin ), STAT=ierr)
       IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei_emp ',ierr)
       ei_emp = 0.0_DP
     END IF

     ecutmass = ecutmass_inp
     emass    = emass_inp
     IF ( ecutmass < 0.0_DP ) &
       CALL errore(' electrons ',' ecutmass out of range ' , 0)

     band_first = .FALSE.

     RETURN
   END SUBROUTINE electrons_setup

!----------------------------------------------------------------------

        SUBROUTINE empty_print_info(iunit)
          !
          USE kinds,            ONLY: DP
          INTEGER, INTENT(IN) :: iunit
          !
          IF ( n_emp > 0 ) WRITE (iunit,620) n_emp, max_emp, ethr_emp
620       FORMAT(3X,'Empty states minimization : states = ',I4, &
             ' maxiter = ',I8,' ethr = ',D10.4)
          !
          RETURN
        END SUBROUTINE empty_print_info

!----------------------------------------------------------------------

        SUBROUTINE empty_init( max_emp_ , ethr_emp_ )

          USE kinds,            ONLY: DP

          INTEGER, INTENT(IN) :: max_emp_
          REAL(DP), INTENT(IN) :: ethr_emp_

          max_emp   = max_emp_
          ethr_emp  = ethr_emp_

          RETURN
        END SUBROUTINE empty_init


!  ----------------------------------------------


   SUBROUTINE print_eigenvalues( ei_unit, tfile, nfi, tps )
      !
      use constants,  only : autoev 
      USE io_global,  ONLY : stdout, ionode
      !
      INTEGER,  INTENT(IN) :: ei_unit
      LOGICAL,  INTENT(IN) :: tfile
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps
      !
      INTEGER :: i, j, ik
      !
      IF ( tfile ) THEN
          WRITE(ei_unit,30) nfi, tps
      END IF
      !
      ik = 1
      !
      DO j = 1, nspin
         !
         WRITE( stdout,1002) ik, j
         WRITE( stdout,1004) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         !
         IF( n_emp .GT. 0 ) THEN
            WRITE( stdout,1005) ik, j
            WRITE( stdout,1004) ( ei_emp( i, j ) * autoev , i = 1, n_emp )
            WRITE( stdout,1006) ( ei_emp( 1, j ) - ei( nupdwn(j), j ) ) * autoev
         END IF
         !
         IF( tfile ) THEN
            WRITE(ei_unit,1010) ik, j
            WRITE(ei_unit,1020) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
            IF( n_emp .GT. 0 ) THEN
               WRITE(ei_unit,1011) ik, j
               WRITE(ei_unit,1020) ( ei_emp( i, j ) * autoev , i = 1, n_emp )
               WRITE(ei_unit,1021) ( ei_emp( 1, j ) - ei( nupdwn(j), j ) ) * autoev
            END IF
         END IF
         !
      END DO
      !
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
 1002 FORMAT(/,3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1005 FORMAT(/,3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1004 FORMAT(10F8.2)
 1006 FORMAT(/,3X,'Electronic Gap (eV) = ',F8.2,/)
 1010 FORMAT(3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1011 FORMAT(3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1020 FORMAT(10F8.2)
 1021 FORMAT(3X,'Electronic Gap (eV) = ',F8.2)
 1030 FORMAT(3X,'nfill = ', I4, ', nempt = ', I4, ', kp = ', I3, ', spin = ',I2)
      !
      RETURN
   END SUBROUTINE print_eigenvalues



!  ----------------------------------------------

   SUBROUTINE deallocate_electrons
      INTEGER :: ierr
      IF(ALLOCATED(ei))       THEN
            DEALLOCATE(ei, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ei ',ierr )
      END IF
      IF(ALLOCATED(ei_emp))   THEN
            DEALLOCATE(ei_emp, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ei_emp ',ierr )
      END IF
      IF(ALLOCATED(ib_owner)) THEN
            DEALLOCATE(ib_owner, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ib_owner ',ierr )
      END IF
      IF(ALLOCATED(ib_local)) THEN
            DEALLOCATE(ib_local, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ib_local ',ierr )
      END IF
      RETURN
   END SUBROUTINE deallocate_electrons
        



!=----------------------------------------------------------------------------=!
  END MODULE electrons_module
!=----------------------------------------------------------------------------=!
