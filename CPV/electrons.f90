!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sat Feb 12 11:56:12 MET 2000
!  ----------------------------------------------
!  BEGIN manual

      MODULE electrons_module

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE electron_mass_init(alat, delt, hg, ngw)
!  SUBROUTINE band_init(occ, emp, nk)
!  SUBROUTINE electrons_info(nel_out, nx_out, n_out)
!  SUBROUTINE electrons_print_info(unit)
!  SUBROUTINE electrons_setup(nel_inp, emass_inp, ecutmass_inp,f_inp)
!  SUBROUTINE deallocate_electrons
!  ----------------------------------------------
!  END manual

        USE kinds
        USE parallel_toolkit, ONLY: pdspev_drv, dspev_drv, pzhpev_drv, zhpev_drv
        USE parallel_types, ONLY: processors_grid, descriptor, CYCLIC_SHAPE
        USE descriptors_module, ONLY: desc_init, get_local_dims, &
          get_global_dims, owner_of, local_index

        USE electrons_base, ONLY: nbnd, nbndx, nspin, nel, nelt, nupdwn, iupdwn
        USE cp_electronic_mass, ONLY: ecutmass => emass_cutoff
        USE cp_electronic_mass, ONLY: emass
        USE cp_electronic_mass, ONLY: emass_precond


        IMPLICIT NONE
        SAVE

        PRIVATE 

!  ...  declare module-scope variables

        LOGICAL :: band_first = .TRUE.

        INTEGER :: n_emp  ! number of empty states

        INTEGER :: nb_l, nb_g
        INTEGER :: nb_emp_l, nb_emp_g
        INTEGER, ALLOCATABLE :: ib_owner(:)
        INTEGER, ALLOCATABLE :: ib_local(:)

        TYPE (descriptor) :: occ_desc
        TYPE (descriptor) :: emp_desc

        REAL(dbl), ALLOCATABLE :: f(:,:)
        REAL(dbl), ALLOCATABLE :: ei(:,:,:)
        REAL(dbl), ALLOCATABLE :: ei_emp(:,:,:)
        REAL(dbl), ALLOCATABLE :: pmss(:)

!  ...  Fourier acceleration

        LOGICAL :: toccrd = .FALSE.  ! read occupation number from standard input

        INTERFACE eigs
          MODULE PROCEDURE rceigs
        END INTERFACE

        INTERFACE sumgam
          MODULE PROCEDURE rsumgam, csumgam
        END INTERFACE

        PUBLIC :: electrons_setup, eigs
        PUBLIC :: electron_mass_init, band_init, bmeshset
        PUBLIC :: electrons_info, electrons_print_info
        PUBLIC :: deallocate_electrons, fermi_energy
        PUBLIC :: pmss, n_emp, emass, emp_desc, ei_emp
        PUBLIC :: occ_desc, ei, nspin, nelt, nel, nupdwn
        PUBLIC :: nbnd

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE electron_mass_init( alat, hg, ngw )

!  (describe briefly what this routine does...)
!  Calculate: PMSS = EMASS * (2PI/Alat)^2 * |G|^2 / ECUTMASS 
!  ----------------------------------------------

! ...     declare subroutine arguments
          REAL(dbl), INTENT(IN) :: alat
          REAL(dbl), INTENT(IN) :: hg(:)
          INTEGER,   INTENT(IN) :: ngw
          REAL(dbl) :: tpiba2
          INTEGER :: ierr

          ALLOCATE( pmss( ngw ), STAT=ierr)
          IF( ierr/=0 ) CALL errore( ' electron_mass_init ',' allocating pmss ', ierr)

          tpiba2 = ( 8.d0 * datan( 1.d0 ) / alat ) ** 2

          CALL emass_precond( pmss, hg, ngw, tpiba2, ecutmass )
          pmss = emass / pmss
          
          RETURN 
        END SUBROUTINE electron_mass_init

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE band_init( occ )

     !   This subroutine fill in the input array with the 
     !   occupations values read from input

     USE io_global, ONLY: stdout, ionode

     REAL(dbl) :: occ(:,:,:)
     INTEGER   :: ik, i, nk, ispin

     IF( SIZE( occ, 1 ) < nbnd ) &
       CALL errore(' band_init ',' wrong dimension ', 1)
     IF( SIZE( occ, 3 ) < nspin ) &
       CALL errore(' band_init ',' wrong dimension ', 2)

     nk  = SIZE( occ, 2 )
     occ = 0.0d0
          
     IF( nspin == 1 ) THEN
       DO ik = 1, nk
         occ( 1:nbnd, ik, 1 ) = f( 1:nbnd, 1 )
       END DO
     ELSE
       DO ik = 1, nk
         occ( 1:nupdwn(1), ik, 1 ) = f( 1:nupdwn(1), 1 )
       END DO
       DO ik = 1, nk
         occ( 1:nupdwn(2), ik, 2 ) = f( 1:nupdwn(2), 2 )
       END DO
     END IF

     IF( ionode ) THEN
       WRITE( stdout, fmt="(3X,'Occupation number from init')" )
       IF( nspin == 1 ) THEN
         WRITE( stdout, fmt = " (3X, 'nbnd = ', I5 ) " ) nbnd
         WRITE( stdout, fmt = " (3X,10F5.2)" ) ( occ( i, 1, 1 ), i = 1, nbnd )
       ELSE
         DO ispin = 1, nspin
           WRITE( stdout, fmt = " (3X,'spin = ', I3, ' nbnd = ', I5 ) " ) ispin, nupdwn( ispin )
           WRITE( stdout, fmt = " (3X,10F5.2)" ) ( occ( i, 1, ispin ), i = 1, nupdwn( ispin ) )
         END DO
       END IF
     END IF

     RETURN
   END SUBROUTINE band_init

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE bmeshset

     !   This subroutine initialize the descriptors for the 
     !   distribution across processors of the overlap matrixes 
     !   of sizes ( nx, nx )

     USE mp_global, ONLY: mpime, nproc
     USE mp, ONLY: mp_sum 
     USE bands_mesh, ONLY: bands_grid

     IMPLICIT NONE

     INTEGER :: i, n1, n2, n3, ierr

     IF( band_first ) THEN
       CALL errore(' bmeshset ',' module not initialized ',0)
     END IF

     CALL desc_init( occ_desc, 1, nbnd, 1, 1, 0, 0, 0, 0, 0, 0, bands_grid, &
       CYCLIC_SHAPE, CYCLIC_SHAPE, CYCLIC_SHAPE)
     CALL desc_init( emp_desc, 1, n_emp, 1, 1, 0, 0, 0, 0, 0, 0, bands_grid, &
       CYCLIC_SHAPE, CYCLIC_SHAPE, CYCLIC_SHAPE)

     CALL get_local_dims(occ_desc, nb_l, n2, n3)
     CALL get_global_dims(occ_desc, nb_g, n2, n3)
     CALL get_local_dims(emp_desc, nb_emp_l, n2, n3)
     CALL get_global_dims(emp_desc, nb_emp_g, n2, n3)

     ALLOCATE( ib_owner( nbndx ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_owner ', ierr)
     ALLOCATE( ib_local( nbndx ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_local ', ierr)

     !  here define the association between processors and electronic states
     !  round robin distribution is used

     ib_local =  0
     ib_owner = -1
     DO i = 1, nbndx
       ib_local( i ) = ( i - 1 ) / nproc
       ib_owner( i ) = MOD( ( i - 1 ), nproc )
       IF( mpime <= ib_owner( i ) ) THEN
         ib_local( i ) = ib_local( i ) + 1
       END IF
     END DO

     !  check consistency for fill states

     ierr = 0
     DO i = 1, nbnd
       IF( ib_owner(i) /= owner_of( i, occ_desc, 'R' ) ) THEN
         ierr = 1
       END IF
       IF( mpime == owner_of( i, occ_desc, 'R' ) ) THEN
         IF( ib_local(i) /= local_index(i, occ_desc, 'R' ) ) THEN
           ierr = 2
         END IF
       END IF
     END DO
     CALL mp_sum( ierr )
     IF( ierr /= 0 ) CALL errore( ' bmeshset ',' occ_desc ',ierr)

     !  check consistency for empty states

     ierr = 0
     DO i = 1, n_emp
       IF( ib_owner(i) /= owner_of( i, emp_desc, 'R' ) ) THEN
         ierr = 1
       END IF
       IF( mpime == owner_of( i, emp_desc, 'R' ) ) THEN
         IF( ib_local(i) /= local_index(i, emp_desc, 'R' ) ) THEN
           ierr = 2
         END IF
       END IF
     END DO
     CALL mp_sum( ierr )
     IF( ierr /= 0 )  CALL errore( ' bmeshset ',' emp_desc ',ierr)

     RETURN
   END SUBROUTINE bmeshset

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE electrons_info(nel_out, nx_out, nemp_out, nspin_out)
          INTEGER, INTENT(OUT) :: nel_out, nx_out, nspin_out, nemp_out
            nel_out = nelt
            nx_out  = nbnd
            nemp_out  = n_emp
            nspin_out = nspin
          RETURN
        END SUBROUTINE 

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE electrons_print_info( unit )

          INTEGER, INTENT(IN) :: unit
          INTEGER :: i

          IF( nspin == 1) THEN
            WRITE(unit,6) nelt, nbnd
            WRITE(unit,7) ( f( i, 1 ), i = 1, nbnd )
          ELSE
            WRITE(unit,8) nelt
            WRITE(unit,9) nel(1)
            WRITE(unit,7) ( f( i, 1 ), i = 1, nupdwn(1))
            WRITE(unit,10) nel(2)
            WRITE(unit,7) ( f( i, 2 ), i = 1, nupdwn(2))
          END IF
6         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Number of Electron = ',I5,', of States = ',I5,/ &
                  ,3X,'Occupation numbers :')
7         FORMAT(2X,10F5.2)
8         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Local Spin Density calculation',/ &
                  ,3X,'Number of Electron = ',I5)
9         FORMAT(  3X,'Spins up   = ', I5, ', occupations: ')  
10        FORMAT(  3X,'Spins down = ', I5, ', occupations: ')  

          RETURN
        END SUBROUTINE electrons_print_info

!  ----------------------------------------------
!  ----------------------------------------------
   SUBROUTINE electrons_setup(toccrd_inp, nbnd_ , nelec_ , nelu_ , &
     neld_ , nspin_ , n_emp_ , emass_inp, ecutmass_inp, f_inp, nkp)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

     IMPLICIT NONE
     LOGICAL, INTENT(IN) :: toccrd_inp
     INTEGER, INTENT(IN) :: nbnd_ , nelec_
     INTEGER, INTENT(IN) :: nelu_
     INTEGER, INTENT(IN) :: neld_, nspin_
     INTEGER, INTENT(IN) :: n_emp_

     REAL(dbl),  INTENT(IN) :: emass_inp, ecutmass_inp
     REAL(dbl),  INTENT(IN) :: f_inp(:,:)
     INTEGER, INTENT(IN) :: nkp
     INTEGER :: ierr, i
     CHARACTER(LEN=80) :: msg
 
! ...     end of declarations
!  ----------------------------------------------

     toccrd = toccrd_inp

     nbnd  = nbnd_
     nspin = nspin_
     n_emp = n_emp_

     IF ( nspin /= 1 .AND. nspin /= 2 ) THEN
       WRITE( msg, fmt = '( " invalid nspin: ", I5 )' ) nspin
       CALL errore( ' electrons_setup ', msg, 1 )
     END IF

     IF( toccrd ) THEN

       ! ...   User values

       IF( nspin == 1 ) THEN

         !  here nbnd can be larger than the default

         nel(1)    = nelec_ / 2 + MOD( nelec_ , 2 )
         nel(2)    = nelec_ / 2
         nelt      = nel(1) + nel(2)
         nbnd      = nbnd_
         nupdwn(1) = nbnd_
         nupdwn(2) = nbnd_

       ELSE IF( nspin == 2 ) THEN
       
         IF( nelu_ > 0 .OR. neld_ > 0 ) THEN
           nel(1)  = nelu_
           nel(2)  = neld_
         ELSE
           nel(1)  = nelec_ / 2 + MOD( nelec_ , 2 )
           nel(2)  = nelec_ / 2
         END IF
         nelt  = nel(1) + nel(2)
         nbnd      = nbnd_
         nupdwn(1) = nbnd_
         nupdwn(2) = nbnd_

       END IF

     ELSE

       ! ...   default values

       IF( nspin == 1 ) THEN

         nel(1)  = nelec_ / 2 + MOD( nelec_ , 2 )
         nel(2)  = nelec_ / 2
         nelt    = nel(1) + nel(2)

         IF( nbnd /= ( ( nelt + 1 ) / 2 ) ) THEN
           nbnd = ( ( nelt + 1 ) / 2 )
           CALL infomsg( ' electrons_setup ', ' nbnd modified according to nelt ', nbnd )
         END IF

         nupdwn(1) = nel(1)
         nupdwn(2) = nel(2)

         IF( nbnd /= MAX( nupdwn(1), nupdwn(2) ) ) &
           CALL errore( ' electrons_setup ', ' inconsistent number of states ', 1 )

       ELSE IF( nspin == 2 ) THEN

         IF( nelu_ /= 0 .OR. neld_ /= 0 ) THEN
           nel(1) = nelu_
           nel(2) = neld_
         ELSE
           nel(1)  = nelec_ / 2 + MOD( nelec_ , 2 )
           nel(2)  = nelec_ / 2
         END IF

         nelt    = nel(1) + nel(2)

         nupdwn(1) = nel(1)
         nupdwn(2) = nel(2)

         IF( nbnd /= MAX( nupdwn(1), nupdwn(2) ) ) THEN
           nbnd = MAX( nupdwn(1), nupdwn(2) )
           CALL infomsg( ' electrons_setup ', ' nbnd modified according to nupdwn ', nbnd )
         END IF

       END IF


     END IF


! .. Do some sanity check

     IF( nelt /= nelec_ ) THEN
       CALL errore(' electrons ',' inconsistent nel up, nel down ', 1 )
     END IF

     IF( nelt < 1 ) THEN
       CALL errore(' electrons ',' invalid nelu + neld ', 2 )
     END IF

     IF( nbnd < 1 ) &
       CALL errore(' electrons ',' nbnd out of range ',nbnd)

     IF( nbnd < nupdwn(1) .OR. nbnd < nupdwn(2)) &
       CALL errore(' electrons ',' inconsistent nbnd and nupdwn(1) or nupdwn(2) ',nbnd)

     IF( ( 2 * nbnd ) < nelt ) &
       CALL errore(' electrons ',' TOO FEW STATES ', nbnd )

     nbndx = MAX( nbnd, nupdwn(1), nupdwn(2), n_emp )

     ! .. Set occupation numbers and eingenvalues matrixes

     ALLOCATE( f( nbnd, nspin ), STAT=ierr )
     IF( ierr/=0 ) CALL errore( ' electrons ',' allocating f ',ierr)
     f  = 0.0_dbl

     ALLOCATE( ei( nbnd, nkp, nspin ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei ',ierr)
     ei = 0.0_dbl

     IF( n_emp > 0 ) THEN
       ALLOCATE( ei_emp( n_emp, nkp, nspin ), STAT=ierr)
       IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei_emp ',ierr)
       ei_emp = 0.0_dbl
     END IF

     IF( toccrd ) THEN

       IF( SIZE(f,1) > SIZE(f_inp,1) .OR. SIZE(f,2) > SIZE(f_inp,2) ) & 
         CALL errore(' electrons ',' wrong sizes for f_inp ',1)

       f = f_inp( 1:nbnd , 1:nspin )

     ELSE

       IF( nspin == 1 ) THEN

         DO i = 1, nelt
           f( ( i - 1 ) / 2 + 1, 1 ) = f( ( i - 1 ) / 2 + 1, 1 ) + 1.0d0
         END DO

       ELSE

         f( 1 : nel(1) , 1 ) = 1.0d0
         f( 1 : nel(2) , 2 ) = 1.0d0

       END IF

     END IF

     ecutmass = ecutmass_inp
     emass    = emass_inp
     IF (ecutmass < 0.0_dbl) &
       CALL errore(' electrons ',' ecutmass out of range ' , 0)

     band_first = .FALSE.

     RETURN
   END SUBROUTINE electrons_setup


!=======================================================================

        SUBROUTINE rceigs( nei, gam, cgam, tortho, f, ei, gamma_symmetry )

          USE mp, ONLY: mp_sum
          USE mp_global, ONLY: mpime, nproc, group
          USE energies, only: eig_total_energy
          USE constants, only: au

!OMPUTES:IF (THORTO) 
!           COMPUTES THE EIGENVALUES OF THE COMPLEX HERMITIAN MATRIX GAM
!           THE EIGENVALUES OF GAMMA ARE PRINTED OUT IN ELECTRON VOLTS.
!        ELSE
!           THE EIGENVALUES ARE CALCULATED IN MAIN AS <PSI|H|PSI>, PASSED
!           IN ei() AND PRINTED OUT IN EELECTRON VOLTS.
!        END IF
!
      

! ... ARGUMENTS
          REAL(dbl), INTENT(IN) :: f(:)
          LOGICAL, INTENT(IN) :: tortho, gamma_symmetry
          REAL(dbl), INTENT(INOUT)     ::  gam(:,:)
          COMPLEX(dbl),  INTENT(INOUT) :: cgam(:,:)
          REAL(dbl)  ::  ei(:)
          INTEGER, INTENT(IN) :: nei

      
! ... LOCALS
          INTEGER :: i, nrl, n, ierr
          INTEGER,     ALLOCATABLE :: index(:)
          REAL(dbl),   ALLOCATABLE :: ftmp(:)

          REAL(dbl),   ALLOCATABLE :: vv(:,:)
          REAL(dbl),   ALLOCATABLE :: aux(:)
          REAL(dbl),    ALLOCATABLE :: g(:,:)
          COMPLEX(dbl), ALLOCATABLE :: cg(:,:)
          COMPLEX(dbl), ALLOCATABLE :: caux(:)
!
! ... SUBROUTINE BODY
!    

          IF( nei < 1 ) THEN
            IF( SIZE( ei ) > 1 ) ei = 0.0d0
            RETURN
          END IF

          n   = nei
          IF ( gamma_symmetry ) THEN
            nrl = SIZE( gam, 1)
            IF( n > SIZE( gam, 2 ) ) CALL errore( ' eigs ',' n and gam inconsistent dimensions ',n )
          ELSE
            nrl = SIZE( cgam, 1)
            IF( n > SIZE( cgam, 2 ) ) CALL errore( ' eigs ',' n and cgam inconsistent dimensions ',n )
          END IF

          IF( n < 1 ) CALL errore( ' eigs ',' n wrong value ',n )
          IF( nrl < 1 ) CALL errore( ' eigs ',' nrl wrong value ',nrl )
          IF( n > SIZE( f ) ) CALL errore( ' eigs ',' n and f inconsistent dimensions ',n )

          ALLOCATE( ftmp( n ), STAT=ierr )
          IF( ierr/=0 ) CALL errore( ' eigs ',' allocating ftmp ',ierr )
          ftmp = f( 1:n )
          WHERE ( ftmp < 1.d-6 ) ftmp = 1.d-6

          IF ( gamma_symmetry ) THEN
            ALLOCATE( g(nrl,n), STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' eigs ',' allocating g ',ierr )
            g = gam(1:nrl,1:n)
          ELSE
            ALLOCATE(cg(nrl,n), STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' eigs ',' allocating cg ',ierr )
            cg = cgam(1:nrl,1:n)
          END IF

          IF (tortho) THEN

            IF( gamma_symmetry ) THEN
              IF ( ( nproc < 2 ) .OR. ( n < nproc ) ) THEN
                ALLOCATE( aux( n*(n+1)/2 ), STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' allocating aux ',ierr )
                ! debug
                !WRITE( 6, * ) 
                !WRITE( 6, * ) '  <psi|H|psi> '
                !DO i = 1, SIZE( g, 2 )
                !  WRITE( 6, fmt='(10D12.4)' ) g( :, i ) 
                !END DO
                !WRITE( 6, * ) 
                CALL rpackgam( g, ftmp(:), aux )
                CALL dspev_drv( 'N', 'L', n, aux, ei, g, n )
                DEALLOCATE(aux, STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating aux ',ierr )
              ELSE
                CALL rpackgam( g, ftmp(:) )
                ALLOCATE( vv(nrl,n), STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' allocating vv ',ierr )
                CALL pdspev_drv('N', g, nrl, ei, vv, nrl, nrl, n, nproc, mpime)
                DEALLOCATE( vv, STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating vv ',ierr )
              END IF
            ELSE 
              IF ( (nproc < 2) .OR. (n < nproc) ) THEN
                ALLOCATE(caux(n*(n+1)/2), STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' allocating caux ',ierr )
                CALL cpackgam( cg, ftmp(:), caux )
                CALL zhpev_drv( 'N', 'L', n, caux, ei, cg, n )
                DEALLOCATE(caux, STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating caux ',ierr )
              ELSE
                ALLOCATE(caux(1), STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' allocating caux ',ierr )
                CALL cpackgam( cg, ftmp(:) )
                CALL pzhpev_drv('N', cg, nrl, ei, caux, nrl, nrl, n, nproc, mpime)
                DEALLOCATE(caux, STAT=ierr)
                IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating caux ',ierr )
              END IF
            END IF

          ELSE

            ALLOCATE(index(n), STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' eigs ',' allocating index ',ierr )
            ei = 0.0_dbl
            DO i = 1, n
              IF ( ib_owner(i) == mpime ) THEN
                IF ( gamma_symmetry ) THEN
                  ei(i) = gam(ib_local(i),i) / ftmp(i)
                ELSE
                  ei(i) = REAL(cgam(ib_local(i),i)) / ftmp(i)
                END IF
              END IF
            END DO
            CALL mp_sum(ei,group)
            index(1) = 0
            CALL hpsort(n, ei, index)
            DEALLOCATE(index, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating index ',ierr )

          END IF

          DEALLOCATE(ftmp, STAT=ierr)
          IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating ftmp ',ierr )
          IF (gamma_symmetry) THEN
            DEALLOCATE(g, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating g ',ierr )
            ! gam(1:nrl,1:n) = g
          ELSE
            DEALLOCATE(cg, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' eigs ',' deallocating cg ',ierr )
            ! cgam(1:nrl,1:n) cg
          END IF

          RETURN
        END SUBROUTINE


!  ----------------------------------------------

        SUBROUTINE deallocate_electrons
          INTEGER :: ierr
          IF(ALLOCATED(pmss))     THEN
            DEALLOCATE(pmss, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating pmss ',ierr )
          END IF
          IF(ALLOCATED(f))        THEN
            DEALLOCATE(f, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating f ',ierr )
          END IF
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
        

!  ----------------------------------------------
        SUBROUTINE rsumgam(ib, prod, gam)
          USE mp_global, ONLY: mpime, group, nproc
          USE mp, ONLY: mp_sum
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: ib
          REAL(dbl) :: prod(:), gam(:,:)
          CALL mp_sum(prod, group)
          IF ( ib_owner(ib) == mpime ) gam(ib_local(ib),:) = prod(:)
          RETURN
        END SUBROUTINE
!  ----------------------------------------------

!  ----------------------------------------------
        SUBROUTINE csumgam(ib, cprod, cgam)
          USE mp_global, ONLY: mpime, group, nproc
          USE mp, ONLY: mp_sum
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: ib
          COMPLEX(dbl) :: cprod(:), cgam(:,:)
          CALL mp_sum(cprod, group)
          IF(ib_owner(ib) == mpime ) cgam(ib_local(ib),:) = cprod(:)
          RETURN
        END SUBROUTINE
!  ----------------------------------------------

        SUBROUTINE cpackgam(cgam, f, caux)
          USE mp_global, ONLY: mpime, nproc, group
          USE mp, ONLY: mp_sum
          IMPLICIT NONE
          COMPLEX(dbl), INTENT(INOUT)  :: cgam(:,:)
          COMPLEX(dbl), INTENT(OUT), OPTIONAL :: caux(:)
          REAL(dbl), INTENT(IN)  :: f(:)
          INTEGER n, nrl, i, j, k, jl
          nrl = SIZE(cgam, 1)
          n   = SIZE(cgam, 2)
          IF( PRESENT( caux ) ) THEN
            caux = CMPLX(0.0d0)
            IF( mpime < n ) THEN
              DO i = 1, n
                j = mpime + 1
                DO jl = 1, nrl
                  IF( j >= i ) THEN
                    !   maps (j,i) index to low-tri packed (k) index
                    k = (i-1)*n + j - i*(i-1)/2  
                    caux(k) = cgam(jl,i) / f(j)
                  END IF
                  j = j + nproc
                END DO
              END DO
            END IF
            CALL mp_sum(caux, group)
          ELSE
            IF( mpime < n ) THEN
              DO i = 1, n
                j = mpime + 1
                DO jl = 1, nrl
                  cgam( jl, i ) = cgam( jl, i ) / f(j)
                  j = j + nproc
                END DO
              END DO
            END IF
          END IF
          RETURN
        END SUBROUTINE

        SUBROUTINE rpackgam(gam, f, aux)
          USE mp_global, ONLY: mpime, nproc, group
          USE mp, ONLY: mp_sum
          IMPLICIT NONE
          REAL(dbl), INTENT(INOUT)  :: gam(:,:)
          REAL(dbl), INTENT(OUT), OPTIONAL :: aux(:)
          REAL(dbl), INTENT(IN)  :: f(:)
          INTEGER n, nrl, i, j, k, jl
          nrl = SIZE(gam, 1)
          n   = SIZE(gam, 2)
          IF( PRESENT( aux ) ) THEN
            aux = 0.0d0
            IF( mpime < n ) THEN
              DO i = 1, n
                j = mpime + 1
                DO jl = 1, nrl
                  IF( j >= i ) THEN
                    !   maps (j,i) index to low-tri packed (k) index
                    k = (i-1)*n + j - i*(i-1)/2  
                    aux(k) = gam(jl,i) / f(j)
                  END IF
                  j = j + nproc
                END DO
              END DO
            END IF
            CALL mp_sum(aux, group)
          ELSE
            IF( mpime < n ) THEN
              DO i = 1, n
                j = mpime + 1
                DO jl = 1, nrl
                  gam(jl,i) = gam(jl,i) / f(j)
                  j = j + nproc
                END DO
              END DO
            END IF
          END IF
          RETURN
        END SUBROUTINE

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE fermi_energy(kp, eig, occ, wke, ef, qtot, temp, sume)

!  this routine computes Fermi energy and weights of occupied states
!  using an improved Gaussian-smearing method
!  refs: C.L.Fu and K.M.Ho, Phys.Rev. B28, 5480 (1983)
!        M.Methfessel and A.T.Paxton Phys.Rev. B40 (15 aug. 89).
!
!  taken from APW code by J. Soler and A. Williams (jk+ss)
!  added computation of occupation numbers without k-point weight
!  ----------------------------------------------
!  END manual

      USE brillouin, ONLY: kpoints
      USE io_global, ONLY: stdout

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(dbl) :: occ(:,:,:)
      TYPE (kpoints) :: kp
      REAL(dbl) ef, qtot, temp, sume
      REAL(dbl) eig(:,:,:), wke(:,:,:)
      REAL(dbl), PARAMETER  :: tol = 1.d-10
      INTEGER,   PARAMETER  :: nitmax = 100
      INTEGER ne, nk, nspin

! ... declare functions
      REAL(dbl) stepf

! ... declare other variables
      REAL(dbl) sumq,emin,emax,fac,t,drange
      INTEGER ik,ispin,ie,iter

!  end of declarations
!  ----------------------------------------------

      nspin = SIZE( occ, 3)
      nk    = SIZE( occ, 2)
      ne    = SIZE( occ, 1)
      sumq=0.d0
      sume=0.d0
      emin=eig(1,1,1)
      emax=eig(1,1,1)
      fac=2.d0
      IF (nspin.EQ.2) fac=1.d0

      DO ik=1,nk
        DO ispin=1,nspin
          DO ie=1,ne
            wke(ie,ik,ispin) = kp%weight(ik) * fac
            occ(ie,ik,ispin) = fac
            sumq=sumq+wke(ie,ik,ispin)
            sume=sume+wke(ie,ik,ispin)*eig(ie,ik,ispin)
            emin=MIN(emin,eig(ie,ik,ispin))
            emax=MAX(emax,eig(ie,ik,ispin))
          END DO
        END DO
      END DO
      ef=emax
      IF (abs(sumq-qtot).LT.tol) RETURN
      IF (sumq.LT.qtot) THEN
        WRITE( stdout,*) 'FERMIE: NOT ENOUGH STATES'
        WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
        STOP
      END IF
      t = MAX(temp,1.d-6)
      drange = t * SQRT( - LOG( tol*.01d0) )
      emin = emin - drange
      emax = emax + drange
      DO iter = 1, nitmax
        ef   = 0.5d0 * (emin+emax)
        sumq = 0.d0
        sume = 0.d0
        DO ik = 1, nk
          DO ispin = 1, nspin
            DO ie = 1, ne
              wke(ie,ik,ispin) = fac / 2.d0 * kp%weight(ik) * stepf((eig(ie,ik,ispin)-ef)/t)
              occ(ie,ik,ispin) = fac / 2.d0 * stepf((eig(ie,ik,ispin)-ef)/t)
              sumq = sumq + wke(ie,ik,ispin)
              sume = sume + wke(ie,ik,ispin) * eig(ie,ik,ispin)
            END DO
          END DO
        END DO
        IF (ABS(sumq-qtot).LT.tol) RETURN
        IF (sumq.LE.qtot) emin=ef
        IF (sumq.GE.qtot) emax=ef
      END DO

      WRITE( stdout,*) 'FERMIE: ITERATION HAS NOT CONVERGED.'
      WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
      STOP

      END SUBROUTINE



!  ----------------------------------------------
      END MODULE electrons_module
!  ----------------------------------------------

