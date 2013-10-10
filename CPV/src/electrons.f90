!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE electrons_module
!=----------------------------------------------------------------------------=!
        USE kinds
        USE dspev_module,       ONLY: pdspev_drv, dspev_drv
        USE electrons_base,     ONLY: nbnd, nbndx, nbsp, nbspx, nspin, nel, nelt, &
                                      nupdwn, iupdwn, telectrons_base_initval, f, &
                                      nudx, nupdwn_bgrp, iupdwn_bgrp, nudx_bgrp, &
                                      nbsp_bgrp, nbspx_bgrp, i2gupdwn_bgrp

        USE cp_electronic_mass, ONLY: ecutmass => emass_cutoff, emass, emass_precond


        IMPLICIT NONE
        SAVE

        PRIVATE 

!  ...  declare module-scope variables

        INTEGER, PARAMETER :: nspinx  = 2
        LOGICAL :: band_first = .TRUE.

        INTEGER :: nb_l(nspinx)    =  0  ! local number of states ( for each spin components )
        !
        INTEGER, ALLOCATABLE :: ib_owner(:)
        INTEGER, ALLOCATABLE :: ib_local(:)

        REAL(DP), ALLOCATABLE :: ei(:,:)

!  ...  Fourier acceleration

        LOGICAL :: toccrd = .FALSE.  ! read occupation number from standard input

        PUBLIC :: electrons_setup
        PUBLIC :: bmeshset, occn_info
        PUBLIC :: deallocate_electrons
        PUBLIC :: ib_owner, ib_local, nb_l
        PUBLIC :: ei
        PUBLIC :: print_eigenvalues
        PUBLIC :: distribute_c, collect_c
        PUBLIC :: distribute_b, collect_b
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
    SUBROUTINE distribute_b( b, b_bgrp )
      REAL(DP), INTENT(IN) :: b(:,:)
      REAL(DP), INTENT(OUT) :: b_bgrp(:,:)
      INTEGER :: iss, n1, n2, m1, m2
      DO iss = 1, nspin
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         b_bgrp(:,n1:n2) = b(:,m1:m2)
      END DO
      RETURN
    END SUBROUTINE distribute_b
!
    SUBROUTINE collect_b( b, b_bgrp )
      USE mp_global, ONLY : inter_bgrp_comm
      USE mp,        ONLY : mp_sum
      REAL(DP), INTENT(OUT) :: b(:,:)
      REAL(DP), INTENT(IN)  :: b_bgrp(:,:)
      INTEGER :: iss, n1, n2, m1, m2
      b = 0.0d0
      DO iss = 1, nspin
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         b(:,m1:m2) = b_bgrp(:,n1:n2)
         !write(1000+mpime,*) 'n1, n2 = ', n1, n2 ! debug
         !write(1000+mpime,*) 'm1, m2 = ', m1, m2 ! debug
      END DO
      CALL mp_sum( b, inter_bgrp_comm )
      RETURN
    END SUBROUTINE collect_b


    SUBROUTINE distribute_c( c, c_bgrp )
      COMPLEX(DP), INTENT(IN) :: c(:,:)
      COMPLEX(DP), INTENT(OUT) :: c_bgrp(:,:)
      INTEGER :: iss, n1, n2, m1, m2
      DO iss = 1, nspin
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         c_bgrp(:,n1:n2) = c(:,m1:m2)
      END DO
      RETURN
    END SUBROUTINE distribute_c
!
    SUBROUTINE collect_c( c, c_bgrp )
      USE mp_global, ONLY : inter_bgrp_comm
      USE mp,        ONLY : mp_sum
      COMPLEX(DP), INTENT(OUT) :: c(:,:)
      COMPLEX(DP), INTENT(IN)  :: c_bgrp(:,:)
      INTEGER :: iss, n1, n2, m1, m2
      c = 0.0d0
      DO iss = 1, nspin
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         c(:,m1:m2) = c_bgrp(:,n1:n2)
         !write(1000+mpime,*) 'n1, n2 = ', n1, n2 ! debug
         !write(1000+mpime,*) 'm1, m2 = ', m1, m2 ! debug
      END DO
      CALL mp_sum( c, inter_bgrp_comm )
      RETURN
    END SUBROUTINE collect_c

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE bmeshset

     !   This subroutine initialize the variables for the 
     !   distribution across processors of the overlap matrixes 
     !   of sizes ( nx, nx )

     USE mp_global, ONLY: me_bgrp, nproc_bgrp

     IMPLICIT NONE

     INTEGER :: i, ierr

     IF( band_first ) THEN
       CALL errore(' bmeshset ',' module not initialized ',0)
     END IF

     DO i = 1, nspin 
       !
       IF( i > nspinx ) CALL errore( ' bmeshset ',' spin too large ', i)
       !
       nb_l( i ) = nupdwn( i ) / nproc_bgrp
       IF( me_bgrp < MOD( nupdwn( i ), nproc_bgrp ) ) nb_l( i ) = nb_l( i ) + 1
       !
     END DO

     IF( ALLOCATED( ib_owner ) ) DEALLOCATE( ib_owner )
     ALLOCATE( ib_owner( nbndx ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_owner ', ierr)
     IF( ALLOCATED( ib_local ) ) DEALLOCATE( ib_local )
     ALLOCATE( ib_local( nbndx ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_local ', ierr)

     !  here define the association between processors and electronic states
     !  round robin distribution is used

     ib_local =  0
     ib_owner = -1
     DO i = 1, nbndx
       ib_local( i ) = ( i - 1 ) / nproc_bgrp        !  local index of the i-th band 
       ib_owner( i ) = MOD( ( i - 1 ), nproc_bgrp )  !  owner of th i-th band
       IF( me_bgrp <= ib_owner( i ) ) THEN
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


   SUBROUTINE electrons_setup( emass_inp, ecutmass_inp )

     IMPLICIT NONE
     REAL(DP),  INTENT(IN) :: emass_inp, ecutmass_inp
     INTEGER :: ierr, i
 

     IF( .NOT. telectrons_base_initval ) &
       CALL errore( ' electrons_setup ', ' electrons_base not initialized ', 1 )

     !
     IF( ALLOCATED( ei ) ) DEALLOCATE( ei )
     ALLOCATE( ei( nudx, nspin ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei ',ierr)
     ei = 0.0_DP

     ecutmass = ecutmass_inp
     emass    = emass_inp
     IF ( ecutmass < 0.0_DP ) &
       CALL errore(' electrons ',' ecutmass out of range ' , 0)

     band_first = .FALSE.

     RETURN
   END SUBROUTINE electrons_setup

!----------------------------------------------------------------------

   SUBROUTINE print_eigenvalues( ei_unit, tfile, tstdout, nfi, tps )
      !
      use constants,  only : autoev 
      USE io_global,  ONLY : stdout, ionode
      !
      INTEGER,  INTENT(IN) :: ei_unit
      LOGICAL,  INTENT(IN) :: tfile, tstdout
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
         IF( tstdout ) THEN
            WRITE( stdout,1002) ik, j
            WRITE( stdout,1004) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         END IF
         !
         IF( tfile ) THEN
            WRITE(ei_unit,1010) ik, j
            WRITE(ei_unit,1020) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         END IF
         !
      END DO
      !
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
 1002 FORMAT(/,3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1004 FORMAT(10F8.2)
 1010 FORMAT(3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1020 FORMAT(10F8.2)
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
