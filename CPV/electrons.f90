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
        USE parallel_toolkit,   ONLY: pdspev_drv, dspev_drv, pzhpev_drv, zhpev_drv
        USE electrons_base,     ONLY: nbnd, nbndx, nbsp, nbspx, nspin, nel, nelt, &
                                      nupdwn, iupdwn, telectrons_base_initval, f
        USE cp_electronic_mass, ONLY: ecutmass => emass_cutoff
        USE cp_electronic_mass, ONLY: emass
        USE cp_electronic_mass, ONLY: emass_precond


        IMPLICIT NONE
        SAVE

        PRIVATE 

!  ...  declare module-scope variables

        LOGICAL :: band_first = .TRUE.

        INTEGER :: n_emp       =  0  ! number of empty states

        INTEGER :: nb_l(nspinx)    =  0  ! local number of states ( for each spin components )
        INTEGER :: n_emp_l(nspinx) =  0
        INTEGER, ALLOCATABLE :: ib_owner(:)
        INTEGER, ALLOCATABLE :: ib_local(:)

        REAL(DP), ALLOCATABLE :: ei(:,:,:)
        REAL(DP), ALLOCATABLE :: ei_emp(:,:,:)
        REAL(DP), ALLOCATABLE :: pmss(:)

!  ...  Fourier acceleration

        LOGICAL :: toccrd = .FALSE.  ! read occupation number from standard input

        INTERFACE eigs
          MODULE PROCEDURE rceigs
        END INTERFACE

        PUBLIC :: electrons_setup, eigs, cp_eigs
        PUBLIC :: pmss_init, occn_init, bmeshset, occn_info
        PUBLIC :: deallocate_electrons, fermi_energy
        PUBLIC :: pmss, n_emp, emass, ei_emp, n_emp_l, ib_owner, ib_local, nb_l
        PUBLIC :: ei, nspin, nelt, nupdwn
        PUBLIC :: nbnd
        PUBLIC :: print_eigenvalues

!
!  end of module-scope declarations
!
!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!


   SUBROUTINE pmss_init( ema0bg, ngw )
     !
     !  Calculate: PMSS = EMASS * (2PI/Alat)^2 * |G|^2 / ECUTMASS 
     !
     REAL(DP), INTENT(IN) :: ema0bg(:)
     INTEGER,   INTENT(IN) :: ngw
     INTEGER :: ierr
     !
     ALLOCATE( pmss( ngw ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' pmss_init ',' allocating pmss ', ierr)
     !
     pmss = emass / ema0bg
     !     
     RETURN 
   END SUBROUTINE pmss_init

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE occn_init( occ )

     !   This subroutine fill in the input array with the 
     !   occupations values read from input

     USE io_global, ONLY: stdout, ionode

     REAL(DP) :: occ(:,:,:)
     INTEGER   :: ik, i, nk, ispin

     IF( SIZE( occ, 1 ) < nbnd ) &
       CALL errore(' occn_init ',' wrong dimension ', 1)
     IF( SIZE( occ, 3 ) < nspin ) &
       CALL errore(' occn_init ',' wrong dimension ', 2)

     nk  = SIZE( occ, 2 )
     occ = 0.0d0
     !     
     IF( nspin == 1 ) THEN
       DO ik = 1, nk
         occ( 1:nbnd, ik, 1 ) = f( 1:nbnd )
       END DO
     ELSE IF( nspin == 2 ) THEN
       DO ik = 1, nk
         occ( 1:nupdwn(1), ik, 1 ) = f( 1:nupdwn(1) )
       END DO
       DO ik = 1, nk
         occ( 1:nupdwn(2), ik, 2 ) = f( iupdwn(2) : ( iupdwn(2) + nupdwn(2) - 1 ) )
       END DO
     ELSE
       WRITE( stdout, * ) ' nspin = ', nspin
       CALL errore(' occn_init ',' nspin not implemented ', 3)
     END IF
     !
     RETURN
   END SUBROUTINE occn_init


   SUBROUTINE occn_info( occ )
     !
     !   This subroutine prints occupation numbers to stdout
     !
     USE io_global, ONLY: stdout, ionode
     !
     REAL(DP) :: occ(:,:,:)
     INTEGER   :: ik, i, nk, ispin
     !
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
     !
     RETURN
   END SUBROUTINE occn_info

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE bmeshset

     !   This subroutine initialize the variables for the 
     !   distribution across processors of the overlap matrixes 
     !   of sizes ( nx, nx )

     USE mp_global, ONLY: mpime, nproc

     IMPLICIT NONE

     INTEGER :: i, ierr

     IF( band_first ) THEN
       CALL errore(' bmeshset ',' module not initialized ',0)
     END IF

     DO i = 1, nspin 
       IF( i > nspinx ) CALL errore( ' bmeshset ',' spin too large ', i)
       nb_l( i ) = nupdwn( i ) / nproc
       IF( mpime < MOD( nupdwn( i ), nproc ) ) nb_l( i ) = nb_l( i ) + 1
       n_emp_l( i ) = n_emp / nproc
       IF( mpime < MOD( n_emp, nproc ) ) n_emp_l( i ) = n_emp_l( i ) + 1
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
       ib_local( i ) = ( i - 1 ) / nproc        !  local index of the i-th band 
       ib_owner( i ) = MOD( ( i - 1 ), nproc )  !  owner of th i-th band
       IF( mpime <= ib_owner( i ) ) THEN
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


   SUBROUTINE electrons_setup( n_emp_ , emass_inp, ecutmass_inp, nkp )

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: n_emp_
     REAL(DP),  INTENT(IN) :: emass_inp, ecutmass_inp
     INTEGER, INTENT(IN) :: nkp
     INTEGER :: ierr, i
 

     IF( .NOT. telectrons_base_initval ) &
       CALL errore( ' electrons_setup ', ' electrons_base not initialized ', 1 )

     n_emp = n_emp_

     IF( n_emp > nbndx ) &
       CALL errore( ' electrons_setup ', ' too many empty states ', 1 )

     IF( ALLOCATED( ei ) ) DEALLOCATE( ei )
     ALLOCATE( ei( nbnd, nkp, nspin ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei ',ierr)
     ei = 0.0_DP

     IF( ALLOCATED( ei_emp ) ) DEALLOCATE( ei_emp )
     IF( n_emp > 0 ) THEN
       ALLOCATE( ei_emp( n_emp, nkp, nspin ), STAT=ierr)
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


!  ----------------------------------------------


   SUBROUTINE print_eigenvalues( ei_unit, tfile, nfi, tps )
      !
      use constants,  only : au 
      USE io_global,  ONLY : stdout, ionode
      !
      INTEGER,  INTENT(IN) :: ei_unit
      LOGICAL,  INTENT(IN) :: tfile
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps
      !
      INTEGER :: ik, i, j, nkpt
      !
      nkpt  = 1
      !
      IF ( tfile ) THEN
          WRITE(ei_unit,30) nfi, tps
      END IF
      !
      DO ik = 1, nkpt
         !
         DO j = 1, nspin
            !
            WRITE( stdout,1002) ik, j
            WRITE( stdout,1004) ( ei( i, ik, j ) * au, i = 1, nupdwn(j) )
            !
            IF( n_emp .GT. 0 ) THEN
               WRITE( stdout,1005) ik, j
               WRITE( stdout,1004) ( ei_emp( i, ik, j ) * au , i = 1, n_emp )
               WRITE( stdout,1006) ( ei_emp( 1, ik, j ) - ei( nupdwn(j), ik, j ) ) * au
            END IF
            !
            IF( tfile ) THEN
               WRITE(ei_unit,1010) ik, j
               WRITE(ei_unit,1020) ( ei( i, ik, j ) * au, i = 1, nupdwn(j) )
               IF( n_emp .GT. 0 ) THEN
                  WRITE(ei_unit,1011) ik, j
                  WRITE(ei_unit,1020) ( ei_emp( i, ik, j ) * au , i = 1, n_emp )
                  WRITE(ei_unit,1021) ( ei_emp( 1, ik, j ) - ei( nupdwn(j), ik, j ) ) * au
               END IF
            END IF
            !
         END DO
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
          REAL(DP), INTENT(IN) :: f(:)
          LOGICAL, INTENT(IN) :: tortho, gamma_symmetry
          REAL(DP), INTENT(INOUT)     ::  gam(:,:)
          COMPLEX(DP),  INTENT(INOUT) :: cgam(:,:)
          REAL(DP)  ::  ei(:)
          INTEGER, INTENT(IN) :: nei

      
! ... LOCALS
          INTEGER :: i, nrl, n, ierr
          INTEGER,     ALLOCATABLE :: index(:)
          REAL(DP),   ALLOCATABLE :: ftmp(:)

          REAL(DP),   ALLOCATABLE :: vv(:,:)
          REAL(DP),   ALLOCATABLE :: aux(:)
          REAL(DP),    ALLOCATABLE :: g(:,:)
          COMPLEX(DP), ALLOCATABLE :: cg(:,:)
          COMPLEX(DP), ALLOCATABLE :: caux(:)
!
! ... SUBROUTINE BODY
!    

          IF( nei < 1 ) THEN
            IF( SIZE( ei ) > 1 ) ei = 0.0d0
            RETURN
          END IF

          n   = nei
          nrl = n / nproc
          IF( mpime < MOD( n, nproc ) ) nrl = nrl + 1

          IF ( gamma_symmetry ) THEN
            IF( n > SIZE( gam, 2 ) ) CALL errore( ' eigs ',' n and gam inconsistent dimensions ',n )
          ELSE
            IF( n > SIZE( cgam, 2 ) ) CALL errore( ' eigs ',' n and cgam inconsistent dimensions ',n )
          END IF

          IF( n < 1 ) CALL errore( ' eigs ',' n wrong value ',n )
          IF( n > SIZE( f ) ) CALL errore( ' eigs ',' n and f inconsistent dimensions ',n )
          IF( nrl < 1 ) CALL errore( ' eigs ',' nrl wrong value ',nrl )
          IF( nrl > SIZE( f ) ) CALL errore( ' eigs ',' nrl and f inconsistent dimensions ',n )

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
            ei = 0.0_DP
            DO i = 1, n
              IF ( ib_owner(i) == mpime ) THEN
                IF ( gamma_symmetry ) THEN
                  ei(i) = gam(ib_local(i),i) / ftmp(i)
                ELSE
                  ei(i) = DBLE(cgam(ib_local(i),i)) / ftmp(i)
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
        END SUBROUTINE rceigs


!  ----------------------------------------------

        SUBROUTINE deallocate_electrons
          INTEGER :: ierr
          IF(ALLOCATED(pmss))     THEN
            DEALLOCATE(pmss, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating pmss ',ierr )
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
!
!  tools subroutines
!
!  ----------------------------------------------


        SUBROUTINE cpackgam(cgam, f, caux)
          USE mp_global, ONLY: mpime, nproc, group
          USE mp, ONLY: mp_sum
          IMPLICIT NONE
          COMPLEX(DP), INTENT(INOUT)  :: cgam(:,:)
          COMPLEX(DP), INTENT(OUT), OPTIONAL :: caux(:)
          REAL(DP), INTENT(IN)  :: f(:)
          INTEGER n, nrl, i, j, k, jl
          nrl = SIZE(cgam, 1)
          n   = SIZE(cgam, 2)
          IF( PRESENT( caux ) ) THEN
            caux = CMPLX(0.0d0, 0.d0)
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
        END SUBROUTINE cpackgam

!  ----------------------------------------------

        SUBROUTINE rpackgam(gam, f, aux)
          USE mp_global, ONLY: mpime, nproc, group
          USE mp, ONLY: mp_sum
          IMPLICIT NONE
          REAL(DP), INTENT(INOUT)  :: gam(:,:)
          REAL(DP), INTENT(OUT), OPTIONAL :: aux(:)
          REAL(DP), INTENT(IN)  :: f(:)
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
        END SUBROUTINE rpackgam

!  ----------------------------------------------

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE fermi_energy(eig, occ, wke, ef, qtot, temp, sume)

!  this routine computes Fermi energy and weights of occupied states
!  using an improved Gaussian-smearing method
!  refs: C.L.Fu and K.M.Ho, Phys.Rev. B28, 5480 (1983)
!        M.Methfessel and A.T.Paxton Phys.Rev. B40 (15 aug. 89).
!
!  taken from APW code by J. Soler and A. Williams (jk+ss)
!  added computation of occupation numbers without k-point weight
!  ----------------------------------------------
!  END manual

      USE brillouin, ONLY: kpoints, kp
      USE io_global, ONLY: stdout

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP) :: occ(:,:,:)
      REAL(DP) ef, qtot, temp, sume
      REAL(DP) eig(:,:,:), wke(:,:,:)
      REAL(DP), PARAMETER  :: tol = 1.d-10
      INTEGER,   PARAMETER  :: nitmax = 100
      INTEGER ne, nk, nspin

! ... declare functions
      REAL(DP) stepf

! ... declare other variables
      REAL(DP) sumq,emin,emax,fac,t,drange
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

      END SUBROUTINE fermi_energy

!  ----------------------------------------------
!
!
!
!  ----------------------------------------------

   SUBROUTINE cp_eigs( nfi, bec, c0, irb, eigrb, rhor, rhog, rhos, lambdap, lambda, tau0, h )

     use ensemble_dft, only: tens, ismear, z0, c0diag, becdiag
     use electrons_base, only: nx => nbspx, n => nbsp, ispin, f, nspin
     use electrons_base, only: nel, iupdwn, nupdwn, nudx, nelt
     use energies, only: enl, ekin
     use uspp, only: rhovan => becsum
     use grid_dimensions, only: nnr => nnrx
     use io_global, only: stdout

     IMPLICIT NONE

     INTEGER :: nfi
     INTEGER :: irb(:,:)
     COMPLEX(DP) :: c0( :, :, :, : )
     REAL(DP) :: bec( :, : ), rhor( :, : ), rhos( :, : )
     REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
     REAL(DP) :: tau0( :, : ), h( 3, 3 )
     COMPLEX(DP) :: eigrb( :, : ), rhog( :, : )

     real(DP), allocatable:: rhodip(:,:)
     real(DP) :: dipol( 3 )
     LOGICAL, SAVE :: lprimo
     INTEGER :: i


         if(.not.tens) then
            call eigs0(.false.,nspin,nupdwn,iupdwn,.true.,f,nx,lambda,nudx)
         else
            call eigs0(.false.,nspin,nupdwn,iupdwn,.false.,f,nx,lambdap,nudx)
         endif

         WRITE( stdout,*)

        RETURN
      END SUBROUTINE cp_eigs


!=----------------------------------------------------------------------------=!
  END MODULE electrons_module
!=----------------------------------------------------------------------------=!
