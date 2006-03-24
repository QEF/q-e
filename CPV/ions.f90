!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!  --------------------------------------------------------------------------  !
!  --------------------------------------------------------------------------  !
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Wed Nov 17 07:24:21 MET 1999
!  ----------------------------------------------
!  BEGIN manual

   MODULE ions_module

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE print_scaled_positions(unit)
!  SUBROUTINE displacement(ht)
!  SUBROUTINE set_reference_positions(cdm_ref, tau_ref, atoms, ht)
!  SUBROUTINE ions_print_info(tfor,tsdp,tzerop,tv0rd,nv0rd,nbeg, &
!                             taurdr,iunit)
!  SUBROUTINE deallocate_ions
!  REAL(DP) FUNCTION moveions(tsdp,thdyn,nfi,htm,ht0)
!  SUBROUTINE update_ions
!  SUBROUTINE velocity_scaling(nfi,delt,ht)
!  ----------------------------------------------
!  END manual


! ...   declare modules
        USE kinds
        USE atoms_type_module
        USE ions_base, ONLY: nsp, nax, nat, na, pmass, &
             tions_base_init, ions_cofmass, ions_vel, ions_kinene, &
             ions_thermal_stress, ions_vrescal
        USE io_global, ONLY: stdout

        !  nsp   number of atomic species
        !  nax   maximum number of atoms per specie
        !  nat   total number of atoms
        !  na(:) number of atoms per specie
        !  pmass(:)   mass (converted to a.u.) of ions

        IMPLICIT NONE
        SAVE

        PRIVATE


        PUBLIC :: neighbo, print_scaled_positions, displacement
        PUBLIC :: set_velocities
        PUBLIC :: set_reference_positions, atoms_init
        PUBLIC :: update_ions
        PUBLIC :: velocity_scaling
        PUBLIC :: max_ion_forces, moveions
        PUBLIC :: resort_position
       

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS


!  BEGIN manual -------------------------------------------------------------   
!  ------------------- TEMPLATE OF SUBROUTINE COMMENTS ----------------------   
!     SUBROUTINE pippo(arg1, arg2, ...)                                         
!  Descrive briefly what pippo does and the meaning of arguments                
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE neighbo( taus, na, nsp, rcutg, ht ) 

!  Calculate, for each atom, the neighbouring atom in a radius rcutg            
!    and print the distance                                                     
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ... declare modules
        USE cell_module, ONLY: s_to_r, pbcs

        IMPLICIT NONE

! ... declare subroutine arguments
        REAL(DP), INTENT(IN) :: taus(:,:)   !  scaled positions
        INTEGER,   INTENT(IN) :: na(:)       !  number of atoms x specie
        INTEGER,   INTENT(IN) :: nsp         !  number of specie
        REAL(DP), INTENT(IN) :: rcutg       !  radius to compute neighbour
        REAL(DP), INTENT(IN) :: ht(3,3)     !  system cell
    
! ... declare other variables
        REAL(DP) :: rcutg2, sdist(3), rdist(3), xlm, ylm, zlm, erre2, erre
        INTEGER :: is1, ia, ja, is2
        INTEGER, ALLOCATABLE :: icoor( :, : ), isa( : )
        INTEGER :: is1ia, is2ja

! ... end of declarations
!  ----------------------------------------------

        ALLOCATE( icoor( MAXVAL( na ), nsp ), isa( nsp ) )

        isa( 1 ) = 1
        DO is1 = 2, nsp
          isa( is1 ) = isa( is1 - 1 ) + na( is1 - 1 )
        END DO

        rcutg2 = rcutg**2 

        WRITE( stdout,125)
        WRITE( stdout,126)
        DO is1 = 1, nsp
          DO ia = 1, na(is1)
            icoor(ia,is1) = 0
            is1ia = isa(is1) + ia - 1
            DO is2 = 1, nsp
              DO ja = 1, na(is2)
                is2ja = isa(is2) + ja - 1
                IF(.NOT.(is1 == is2 .AND. ja == ia)) THEN
                  xlm = taus(1,is1ia) - taus(1,is2ja)
                  ylm = taus(2,is1ia) - taus(2,is2ja)
                  zlm = taus(3,is1ia) - taus(3,is2ja)
                  CALL pbcs(xlm,ylm,zlm,sdist(1),sdist(2),sdist(3),1)
                  CALL s_to_r(sdist,rdist,ht)
                  erre2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  IF((erre2 <= rcutg2) .AND. (icoor(ia,is1) == 0)) THEN
                    icoor(ia,is1) = icoor(ia,is1) + 1
                    erre = SQRT(erre2)
                    WRITE( stdout,254) IA,IS1
                    WRITE( stdout,255) ERRE,JA,IS2
                  ELSE IF(erre2 <= rcutg2) THEN
                    icoor(ia,is1) = icoor(ia,is1) + 1
                    erre = SQRT(erre2)
                    WRITE( stdout,255) ERRE,JA,IS2
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO

        DEALLOCATE( icoor, isa )

125     FORMAT(3X,'Neighbours')
126     FORMAT(3X,'----------')
254     FORMAT(2(I3))
255     FORMAT(2X,F12.6,2(I3))

        RETURN
      END SUBROUTINE neighbo


!  BEGIN manual -------------------------------------------------------------   

        SUBROUTINE print_scaled_positions(atoms, unit, string)

!  Print onto "unit" the scaled positions of every atoms                        
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

          INTEGER :: unit
          TYPE (atoms_type) :: atoms
          CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: string
          INTEGER :: is,ia,k
          IF(PRESENT(string)) WRITE(unit,100) string
          DO is = 1, atoms%nsp
            WRITE(unit,1300) is, atoms%na(is)
            DO ia = atoms%isa(is), atoms%isa(is) + atoms%na(is) - 1
              WRITE(unit,555) atoms%label(is), (atoms%taus(k,ia), k = 1,3)
            END DO
          END DO
          RETURN
 100      FORMAT(/,3X,'Scaled atomic positions ',A50)
 1300     FORMAT(3X,'Species ',I3,' atoms = ',I4)
  555     FORMAT(3X, A4, 3(1X,F12.6), 3L6)
        END SUBROUTINE print_scaled_positions


!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE anneions( anner, atoms_m, atoms_0, atoms_p, ht)

!  Descrive briefly what it does
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ... declare modules
      USE cell_module, ONLY: r_to_s, boxdimensions
      USE time_step, ONLY: delt

      IMPLICIT NONE
 
! ... declare subroutine arguments
      REAL(DP) :: anner
      TYPE (boxdimensions), INTENT(IN) :: ht
      TYPE (atoms_type), INTENT(IN) :: atoms_0, atoms_m
      TYPE (atoms_type) :: atoms_p

! ... declare other variables
      REAL(DP) :: fions(3)
      REAL(DP) :: alfap, dt2hbm, dt2
      INTEGER   :: k, j, i, isa
 
! ... end of declarations
!  ----------------------------------------------

      dt2 = delt ** 2
      alfap = .5d0 * SQRT(anner)
      isa   = 0
      DO k = 1, atoms_0%nsp
        dt2hbm = 0.5_DP * dt2 / atoms_0%m(k)
        DO j = 1, atoms_0%na(k)
          isa = isa + 1
          CALL r_to_s( atoms_0%for(:,isa), fions, ht)
          DO i = 1, 3
            IF( atoms_0%mobile(i,isa) > 0 ) THEN
              atoms_p%taus(i,isa) = atoms_0%taus(i,isa) + alfap * atoms_m%taus(i,isa) + dt2hbm * fions(i)
            END IF
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE anneions

!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE displacement(dis, atoms, tau_ref, ht)

!  Calculate the sum of the quadratic displacements of the atoms in the ref.    
!    of cdm respect to the initial positions.                                   
!    tau_ref: initial positions in real units in the ref. of cdm                
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


!  ----------------------------------------------
!  att!     tau_ref: starting position in center-of-mass ref. in real units
!  ----------------------------------------------

! ...   declare modules
        USE cell_module, ONLY: s_to_r, boxdimensions

        IMPLICIT NONE

! ...   declare subroutine arguments
        TYPE (boxdimensions), INTENT(IN) :: ht
        TYPE (atoms_type), INTENT(IN) :: atoms
        REAL (DP), INTENT(OUT) :: dis(:)
        REAL (DP), INTENT(IN)  :: tau_ref(:,:)

! ...   declare other variables
        REAL(DP) :: sdist(3), rdist(3), r2, cdm(3)
        INTEGER   :: i, j, k, isa

! ...   end of declarations
!  ----------------------------------------------

! ...   Compute the current value of cdm "Centro Di Massa"
        CALL ions_cofmass(atoms%taus, atoms%m, atoms%na, atoms%nsp, cdm)
        IF( SIZE( dis ) < atoms%nsp ) &
          CALL errore(' displacement ',' size of dis too small ', 1)
        isa = 0
        DO k = 1, atoms%nsp
          dis(k) = 0.0_DP
          r2     = 0.0_DP
          DO j = 1, atoms%na(k)
            isa = isa + 1
            sdist = atoms%taus(:,isa) - cdm
            CALL s_to_r(sdist, rdist, ht)
            r2 = r2 + SUM( ( rdist(:) - tau_ref(:,isa) )**2 )
          END DO
          dis(k) = dis(k) + r2 / DBLE(atoms%na(k))
        END DO

        RETURN
      END SUBROUTINE displacement



!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE set_velocities(atoms_m, atoms_0, vel_srt, ht_0, delt)

!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

        USE cell_module, ONLY: boxdimensions
        USE cell_base, ONLY: r_to_s, s_to_r
        USE constants, ONLY: angstrom_au

        TYPE (atoms_type) :: atoms_m, atoms_0
        TYPE (boxdimensions) :: ht_0
        REAL(DP), INTENT(IN) :: delt
        REAL(DP), INTENT(IN) :: vel_srt( :, : )
        
        INTEGER :: i, ia, nat
        REAL(DP) :: sv0( 3 )

        nat = atoms_0%nat

        DO ia = 1, nat
          atoms_m%taus( :, ia ) = atoms_0%taus( :, ia )
          atoms_m%taur( :, ia ) = atoms_0%taur( :, ia )
          CALL r_to_s( vel_srt(:,ia), sv0(:), ht_0)
          DO i = 1, 3
            IF( atoms_0%mobile( i, ia ) > 0 ) THEN
              atoms_0%taus( i, ia ) = atoms_0%taus( i, ia ) + sv0( i ) * delt
            END IF
          ENDDO
          CALL s_to_r( atoms_0%taus( :, ia ), atoms_0%taur( :, ia ), ht_0 )
        END DO

        RETURN
      END SUBROUTINE set_velocities



!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE set_reference_positions(cdm_ref, tau_ref, atoms, ht)

!  Calculate the real position of atoms relative to the center of mass (cdm)
!  and store them in tau_ref
!    cdm_ref: initial position of the center of mass (cdm) in scaled units  
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ...   declare modules
        USE cell_module, ONLY: S_TO_R
        USE cell_module, ONLY: boxdimensions

        IMPLICIT NONE

! ...   declare subroutine arguments
        TYPE(boxdimensions), INTENT(IN) :: ht
        TYPE (atoms_type) :: atoms
        REAL(DP) :: cdm_ref(:), tau_ref(:,:)

! ...   declare other variables
        REAL(DP) :: sdist(3)
        INTEGER :: isa

        CALL ions_cofmass(atoms%taus, atoms%m, atoms%na, atoms%nsp, cdm_ref)
        DO isa = 1, atoms%nat
          sdist( 1:3 ) = atoms%taus( 1:3 , isa ) - cdm_ref( 1:3 )
          CALL s_to_r( sdist, tau_ref(:,isa), ht )
        END DO

        RETURN
      END SUBROUTINE set_reference_positions


!   -------------------------------------------------------------   

       SUBROUTINE atoms_init(atoms_m, atoms_0, atoms_p, stau, ind_srt, if_pos, atml, h)

         !   Allocate and fill the three atoms structure using scaled position an cell

         IMPLICIT NONE
!
         TYPE (atoms_type)    :: atoms_0, atoms_p, atoms_m
         REAL(DP), INTENT(IN) :: h( 3, 3 )
         REAL(DP), INTENT(IN) :: stau(:,:)
         CHARACTER(LEN=3), INTENT(IN) :: atml(:)
         INTEGER, INTENT(IN) :: ind_srt( : )
         INTEGER, INTENT(IN) :: if_pos( :, : )
          
         LOGICAL, ALLOCATABLE :: ismb(:,:)
         INTEGER              :: ia, isa
         LOGICAL              :: nofx

         ALLOCATE( ismb( 3, nat ) )

         ismb = .TRUE.
         nofx = .TRUE.
         DO isa = 1, nat
           ia = ind_srt( isa )
           ismb( 1, isa ) = ( if_pos( 1, ia ) /= 0 )
           ismb( 2, isa ) = ( if_pos( 2, ia ) /= 0 )
           ismb( 3, isa ) = ( if_pos( 3, ia ) /= 0 )
           nofx = nofx .AND. ismb( 1, isa )
           nofx = nofx .AND. ismb( 2, isa )
           nofx = nofx .AND. ismb( 3, isa )
         END DO

         CALL atoms_type_init(atoms_m, stau, ismb, atml, pmass, na, nsp, h)
         CALL atoms_type_init(atoms_0, stau, ismb, atml, pmass, na, nsp, h)
         CALL atoms_type_init(atoms_p, stau, ismb, atml, pmass, na, nsp, h)

         CALL print_scaled_positions( atoms_0, stdout, 'from standard input')

         IF( .NOT. nofx ) THEN
           WRITE( stdout, 10 )
 10        FORMAT( /, &
                   3X, 'Position components with 0 are kept fixed', /, &
                   3X, '  ia  x  y  z ' )
           DO isa = 1, nat
             ia = ind_srt( isa )
             WRITE( stdout, 20 ) isa, if_pos( 1, ia ), if_pos( 2, ia ), if_pos( 3, ia )
           END DO
 20        FORMAT( 3X, I4, I3, I3, I3 )
         END IF

         DEALLOCATE( ismb )

         RETURN
       END SUBROUTINE atoms_init

!  --------------------------------------------------------------------------   

   REAL(DP) FUNCTION moveions(TSDP, thdyn, NFI, atoms_m, atoms_0, atoms_p, htm, ht0, vnosep)

      !  Moves the ions

! ... declare modules
      USE constants,          ONLY : amu_au
      USE cell_module,        ONLY : dgcell, r_to_s, s_to_r, boxdimensions
      use control_flags,      ONLY : tnosep, tcap, tcp, tdampions, lconstrain
      use time_step,          ONLY : delt
      use ions_base,          ONLY : fricp, if_pos
      USE constraints_module, ONLY : check_constraint

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL, INTENT(IN) :: tsdp, thdyn
      INTEGER, INTENT(IN) :: nfi
      TYPE (boxdimensions), INTENT(IN)    :: htm
      TYPE (boxdimensions), INTENT(INOUT) :: ht0
      TYPE (atoms_type) :: atoms_m, atoms_0, atoms_p
      REAL(DP), INTENT(IN) :: vnosep

! ... declare other variables
      REAL(DP), DIMENSION(3,3) :: annep, svarps, tmat, svarpd, gcm1, gcdot
      REAL(DP), DIMENSION(3)   :: fions, svel, rvel, stp, tau_diff
      REAL(DP)                 :: const, dumm, vrnos, gfact, ffact, dt2bym
      REAL(DP)                 :: fordt2, dt2, delthal
      INTEGER :: i, j, k, l, is, ia, isa


! ...     end of declarations

      dt2     = delt ** 2
      fordt2  = 4.0d0 * dt2
      delthal = 0.5d0 * delt
 
      ! ...   Determines DGCELL/DT dynamically and GCM1

      IF( thdyn ) THEN
        CALL invmat( 3, ht0%g, gcm1, dumm )
        CALL dgcell( gcdot, htm, ht0, delt )
      END IF

      IF( tnosep ) THEN
        vrnos = vnosep
      ELSE
        vrnos = 0.0d0
      END IF

      !
! ...   Steepest descent of ionic degrees of freedom 

      IF( tdampions ) THEN

          gfact = 1.0_DP / (1.0_DP + fricp )
          isa = 0
          DO is = 1, atoms_0%nsp
            dt2bym = dt2 / atoms_0%m(is)
            DO ia = 1, atoms_0%na(is)
              isa = isa + 1
              CALL r_to_s( atoms_0%for(:,isa), fions, ht0)
              tau_diff = (atoms_0%taus(:,isa) - atoms_m%taus(:,isa))
              DO k = 1, 3
                IF( atoms_0%mobile(k,isa) > 0 ) THEN
                  atoms_p%taus(k,isa) = atoms_m%taus(k,isa) + &
                    2.0d0 * (tau_diff(k) + 0.5d0 * dt2bym * fions(k)) * gfact
                ELSE
                  atoms_p%taus(k,isa) = atoms_m%taus(k,isa) + 2.0d0 * tau_diff(k)
                END IF
              ENDDO
            ENDDO
          ENDDO

      ELSE IF( tsdp ) THEN

          IF(thdyn) THEN
            annep = MATMUL(gcm1,gcdot)
          ELSE
            annep = 0.D0
          END IF

          isa = 0
          DO is = 1, atoms_0%nsp

!....       SVARPS = DT2 * ( M * ( 1 + GCM1*GCDOT ))^-1
            DO j = 1, 3
              DO i = 1, 3
                svarps(i,j) = annep(i,j)
              END DO
              svarps(j,j) = svarps(j,j) + 1.D0
            END DO
            CALL invmat (3, svarps, tmat, dumm)
            svarps = dt2 * tmat / atoms_0%m(is) 

            DO ia = 1, atoms_0%na(is)
              isa = isa + 1
              CALL r_to_s( atoms_0%for(:,isa), fions, ht0)
              DO k = 1, 3 
                IF(  atoms_0%mobile(k,isa) > 0 )THEN
                  dumm = 0.0d0
                  DO l = 1, 3
                    dumm = dumm + svarps(k,l) * fions(l)
                  END DO
                  atoms_p%taus(k,isa) = atoms_0%taus(k,isa) + dumm
                ELSE
                  atoms_p%taus(k,isa) = atoms_0%taus(k,isa)
                END IF
              END DO
            END DO
          END DO

      ELSE  !....   NEWTON DYNAMICS FOR IONIC DEGREES OF FREEDOM

          IF(TNOSEP.OR.thdyn) THEN
    
!....       Determines friction matrix annep according to nose dynamics
! ...       ANNEP = 1 + (GCM1*GCDOT + VRNOS)*DELT/2 
! ...       SVARPD = ANNEP^-1

            DO j = 1, 3
              DO i = 1, 3
                annep(i,j) = 0.0d0
              END DO
              annep(j,j) = 1.0d0 + vrnos * delthal
            END DO
            IF(thdyn) THEN
              DO j = 1, 3
                DO i = 1, 3
                  DO k = 1, 3
                    annep(i,j) = annep(i,j) + gcm1(i,k) * gcdot(k,j) * delthal
                  END DO
                END DO
              END DO
            END IF
            CALL invmat (3, annep,svarpd,dumm)

          ELSE

            svarpd = RESHAPE( (/ 1.0d0, 0.0d0, 0.0d0, &
                                 0.0d0, 1.0d0, 0.0d0, &
                                 0.0d0, 0.0d0, 1.0d0  &
                               /), (/ 3, 3 /) )
          END IF

! ...     Verlet
          isa = 0
          DO is = 1, atoms_0%nsp
            dt2bym = dt2 / atoms_0%m(is)
            DO ia = 1, atoms_0%na(is) 
              isa = isa + 1
              CALL r_to_s( atoms_0%for(:,isa), fions, ht0 )
              tau_diff = (atoms_0%taus(:,isa) - atoms_m%taus(:,isa))
              DO k = 1, 3
                stp(k) = 2.0_DP * tau_diff(k) + dt2bym * fions(k)
              END DO
              DO k = 1, 3
                IF ( atoms_0%mobile(k,isa) > 0 ) THEN
                  ffact = 0.0_DP
                  DO l = 1, 3
                    ffact = ffact + svarpd(k,l) * stp(l)
                  END DO
                  atoms_p%taus(k,isa) = atoms_m%taus(k,isa) + ffact
                ELSE
                  atoms_p%taus(k,isa) = atoms_0%taus(k,isa)
                END IF
              END DO
            END DO
          END DO

        END IF

        IF (tcp)   THEN
          CALL ions_vel( atoms_0%vels, atoms_p%taus, atoms_m%taus, atoms_0%na, atoms_0%nsp, delt )
          CALL ions_kinene( atoms_0%ekint, atoms_0%vels, atoms_0%na, atoms_0%nsp, ht0%hmat, atoms_0%m )
          CALL velocity_scaling(nfi, tcap, atoms_m, atoms_0, atoms_p, delt, ht0) 
        END IF

        IF ( lconstrain ) THEN
           !
           ! ... constraints are imposed here
           !
           DO ia = 1, atoms_p%nat
              !
              CALL s_to_r( atoms_p%taus(:,ia), atoms_p%taur(:,ia), ht0 )
              !
           END DO
           !
           CALL check_constraint( atoms_p%nat, atoms_p%taur, atoms_0%taur, &
                                  atoms_0%for, if_pos, atoms_p%ityp, 1.D0, &
                                  delt, amu_au )
           !
           DO ia = 1, atoms_p%nat
              !
              CALL r_to_s( atoms_p%taur(:,ia), atoms_p%taus(:,ia), ht0 )
              !
           END DO
           !
        END IF
        !
        CALL ions_vel( atoms_0%vels, atoms_p%taus, atoms_m%taus, atoms_0%na, atoms_0%nsp, delt )
        CALL ions_kinene( atoms_0%ekint, atoms_0%vels, atoms_0%na, atoms_0%nsp, ht0%hmat, atoms_0%m )
        !
        moveions = atoms_0%ekint
        !
        RETURN
        !
      END FUNCTION moveions

!  --------------------------------------------------------------------------   

      SUBROUTINE update_ions(atoms_m, atoms_0, atoms_p)

        !  Update ionic positions and velocities in atoms structures

        IMPLICIT NONE
        TYPE(atoms_type) :: atoms_m, atoms_0, atoms_p
        INTEGER :: is, ia, ub

          ub = atoms_m%nat
          atoms_m%taus(1:3,1:ub) = atoms_0%taus(1:3,1:ub)
          atoms_m%vels(1:3,1:ub) = atoms_0%vels(1:3,1:ub)
          atoms_m%for(1:3,1:ub) = atoms_0%for(1:3,1:ub)
          atoms_0%taus(1:3,1:ub) = atoms_p%taus(1:3,1:ub)
          atoms_0%vels(1:3,1:ub) = atoms_p%vels(1:3,1:ub)
          atoms_0%for(1:3,1:ub) = atoms_p%for(1:3,1:ub)

        RETURN
      END SUBROUTINE update_ions



!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE velocity_scaling(nfi, tcap, atoms_m, atoms_0, atoms_p, delt, ht)

!  Descrive briefly what it does
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ... declare modules
      USE cell_module, ONLY: R_TO_S, boxdimensions
      use constants, ONLY: factem
      use io_global, ONLY: ionode
      use control_flags, ONLY: tolp
      use ions_nose, ONLY: tempw

      implicit none
!
! ... ARGUMENTS
!
      integer, INTENT(IN) :: nfi
      type (boxdimensions), intent(in) :: ht
      type (atoms_type) :: atoms_m, atoms_0, atoms_p
      LOGICAL, INTENT(IN) :: tcap
      REAL(DP), INTENT(IN) :: delt
!
! ... Locals
!
      REAL(DP) :: tempp, temp1, temp2
      REAL(DP), ALLOCATABLE :: fions(:,:)
      INTEGER :: isa
!
! ... Subroutine body
!

        tempp   = atoms_0%ekint * factem
        tempp   = tempp * 2.0_DP / atoms_0%doft
        temp1   = tempw + tolp
        temp2   = tempw - tolp
        isa     = 0
 
        IF( (tempp > temp1) .OR. (tempp < temp2) ) THEN
          IF( ionode ) THEN
            WRITE( stdout,400) 
 400        FORMAT(3X,'Rescaling Ionic Velocities')
          END IF
          ALLOCATE(fions(3,atoms_0%nat) )
          DO isa = 1, atoms_0%nat
            CALL r_to_s( atoms_0%for(:,isa), fions(:,isa), ht)
          END DO
          CALL ions_vrescal( tcap, tempw, tempp, atoms_p%taus, atoms_0%taus, atoms_m%taus, &
            atoms_0%na, atoms_0%nsp, fions, atoms_0%mobile, atoms_0%m, delt )
          DEALLOCATE( fions )
        END IF

      RETURN
      END SUBROUTINE velocity_scaling


!  BEGIN manual -------------------------------------------------------------   

      REAL(DP) FUNCTION max_ion_forces( atoms )

!  Descrive briefly what it does
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

        IMPLICIT NONE 
        TYPE (atoms_type) :: atoms
        INTEGER :: ia
        REAL(DP) :: fmax
        fmax = 0.0d0
        DO ia = 1, atoms%nat
          IF( atoms%mobile(1, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(1, ia) ) )
          IF( atoms%mobile(2, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(2, ia) ) )
          IF( atoms%mobile(3, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(3, ia) ) )
        END DO
        max_ion_forces = fmax
        RETURN
      END FUNCTION max_ion_forces

!
!

!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE resort_position( pos, fion, atoms, isrt, ht )

!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   
         
        ! This subroutine copys positions and forces into
        ! array "pos" and "for" using the same atoms sequence
        ! as in the input file

        USE cell_module, ONLY: s_to_r
        USE cell_module, ONLY: boxdimensions, pbcs

        IMPLICIT NONE

        REAL(DP), INTENT(OUT) :: pos(:,:), fion(:,:)
        TYPE (atoms_type), INTENT(IN) :: atoms
        TYPE (boxdimensions), INTENT(IN) :: ht
        INTEGER, INTENT(IN) :: isrt( : )
        INTEGER :: ia, is, isa, ipos

        isa = 0
        DO is = 1, atoms%nsp
          DO ia = 1, atoms%na(is)
            isa  = isa + 1
            ipos = isrt( isa )
            CALL s_to_r( atoms%taus( : , isa ), pos( :, ipos ), ht )
            fion( :, ipos ) = atoms%for( : , isa )
          END DO
        END DO

        RETURN
      END SUBROUTINE resort_position
     


   END MODULE ions_module
