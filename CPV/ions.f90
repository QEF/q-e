!
! Copyright (C) 2002 FPMD group
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
!  SUBROUTINE neighbo(atoms,rcutg,ht)
!  SUBROUTINE print_scaled_positions(unit)
!  SUBROUTINE anneions(ht)
!  SUBROUTINE displacement(ht)
!  SUBROUTINE cdm_displacement(cdm_ref, atoms, ht)
!  SUBROUTINE set_reference_positions(cdm_ref, tau_ref, atoms, ht)
!  SUBROUTINE ions_setup(nax_inp, pos_inp,ipos_inp, &
!                        anne_inp, anner_inp )
!  SUBROUTINE ions_print_info(tfor,tsdp,tzerop,tv0rd,nv0rd,nbeg, &
!                             taurdr,iunit)
!  SUBROUTINE deallocate_ions
!  REAL(dbl) FUNCTION moveions(tsdp,thdyn,nfi,htm2,htm,ht0)
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

! ...   declare module-scope variables

        TYPE (constrains_class) :: constrains

        REAL(dbl), ALLOCATABLE ::  taui(:,:)   
! ...     taui = real ionic positions in the center of mass reference
! ...     system at istep = 0 
! ...     this array is used to compute mean square displacements, 
! ...     it is initialized when NBEG = -1, NBEG = 0 and TAURDR = .TRUE.
! ...     first index: x,y,z, second index: atom sortred by specie with respect input
! ...     this array is saved in the restart file

        REAL(dbl) ::  cdmi(3)
! ...     center of mass reference system (related to the taui)
! ...     this vector is computed when NBEG = -1, NBEG = 0 and TAURDR = .TRUE.
! ...     this array is saved in the restart file

! ...   velocity rescaling
        INTEGER   :: icapstep

! ...   annealing
        LOGICAL   :: tanne
        REAL(dbl) :: anner

! ...   Damped dynamics
        REAL(dbl)  :: gdelt

        LOGICAL   :: tneighbo   ! print neigbours list for each atom
        REAL(dbl) :: neighbo_radius


! ...   Module private control flag
        LOGICAL :: tsetup = .FALSE.


        PUBLIC :: neighbo, print_scaled_positions, anneions, displacement
        PUBLIC :: cdm_displacement, set_velocities
        PUBLIC :: set_reference_positions, ions_setup, atoms_init
        PUBLIC :: ions_print_info, deallocate_ions, constraints_setup
        PUBLIC :: apply_constraints, update_ions
        PUBLIC :: velocity_scaling
        PUBLIC :: cdmi, taui
        PUBLIC :: tneighbo, neighbo_radius
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

      SUBROUTINE neighbo(atoms, rcutg, ht) 

!  Calculate, for each atom, the neighbouring atom in a radius rcutg            
!    and print the distance                                                     
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ... declare modules
        USE cell_module, ONLY: s_to_r
        USE cell_module, ONLY: boxdimensions, pbcs

        IMPLICIT NONE

! ... declare subroutine arguments
        REAL(dbl), INTENT(IN) :: rcutg
        TYPE (boxdimensions), INTENT(IN) :: ht
        TYPE (atoms_type), INTENT(IN) :: atoms
    
! ... declare other variables
        REAL(dbl) :: rcutg2,sdist(3),rdist(3),xlm,ylm,zlm,erre2,erre
        INTEGER :: is1,ia,ja,is2,icoor(atoms%nax,atoms%nsp)
        INTEGER :: is1ia, is2ja

! ... end of declarations
!  ----------------------------------------------

        rcutg2 = rcutg**2 

        WRITE( stdout,125)
        WRITE( stdout,126)
        DO is1 = 1, atoms%nsp
          DO ia = 1, atoms%na(is1)
            icoor(ia,is1) = 0
            is1ia = atoms%isa(is1) + ia - 1
            DO is2 = 1, atoms%nsp
              DO ja = 1, atoms%na(is2)
                is2ja = atoms%isa(is2) + ja - 1
                IF(.NOT.(is1 == is2 .AND. ja == ia)) THEN
                  xlm = atoms%taus(1,is1ia) - atoms%taus(1,is2ja)
                  ylm = atoms%taus(2,is1ia) - atoms%taus(2,is2ja)
                  zlm = atoms%taus(3,is1ia) - atoms%taus(3,is2ja)
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

      RETURN
  125     FORMAT(3X,'Neighbours')
  126     FORMAT(3X,'----------')
  254     FORMAT(2(I3))
  255     FORMAT(2X,F12.6,2(I3))
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
          return
 100      FORMAT(/,3X,'Scaled atomic positions ',A50)
 1300     FORMAT(3X,'Species ',I3,' atoms = ',I4)
  555     FORMAT(3X, A4, 3(1X,F12.6), 3L6)
        end SUBROUTINE print_scaled_positions


!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE anneions(atoms_m, atoms_0, atoms_p, ht)

!  Descrive briefly what it does
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ... declare modules
      USE cell_module, ONLY: r_to_s, boxdimensions
      USE time_step, ONLY: delt

      IMPLICIT NONE
 
! ... declare subroutine arguments
      TYPE (boxdimensions), INTENT(IN) :: ht
      TYPE (atoms_type), INTENT(IN) :: atoms_0, atoms_m
      TYPE (atoms_type) :: atoms_p

! ... declare other variables
      REAL(dbl) :: fions(3)
      REAL(dbl) :: alfap, dt2hbm, dt2
      INTEGER   :: k, j, i, isa
 
! ... end of declarations
!  ----------------------------------------------

      dt2 = delt ** 2
      alfap = .5d0 * SQRT(anner)
      isa   = 0
      DO k = 1, atoms_0%nsp
        dt2hbm = 0.5_dbl * dt2 / atoms_0%m(k)
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
        REAL (dbl), INTENT(OUT) :: dis(:)
        REAL (dbl), INTENT(IN)  :: tau_ref(:,:)

! ...   declare other variables
        REAL(dbl) :: sdist(3), rdist(3), r2, cdm(3)
        INTEGER   :: i, j, k, isa

! ...   end of declarations
!  ----------------------------------------------

! ...   Compute the current value of cdm "Centro Di Massa"
        CALL ions_cofmass(atoms%taus, atoms%m, atoms%na, atoms%nsp, cdm)
        IF( SIZE( dis ) < atoms%nsp ) &
          CALL errore(' displacement ',' size of dis too small ', 1)
        isa = 0
        DO k = 1, atoms%nsp
          dis(k) = 0.0_dbl
          r2     = 0.0_dbl
          DO j = 1, atoms%na(k)
            isa = isa + 1
            sdist = atoms%taus(:,isa) - cdm
            CALL s_to_r(sdist, rdist, ht)
            r2 = r2 + SUM( ( rdist(:) - tau_ref(:,isa) )**2 )
          END DO
          dis(k) = dis(k) + r2 / REAL(atoms%na(k))
        END DO

        RETURN
      ENd  SUBROUTINE displacement

!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE cdm_displacement(cdm_ref, atoms, ht)

!  Calculate the quadratic displacements of the cdm at the current time step
!    with respect to the initial position                                            
!    cdm_ref: initial position of cdm in scaled units                           
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

! ...   declare modules
        USE cell_module, ONLY: S_TO_R
        USE cell_module, ONLY: boxdimensions

        IMPLICIT NONE

! ...   declare subroutine arguments
        TYPE (boxdimensions), INTENT(IN) :: ht
        TYPE (atoms_type), INTENT(IN) :: atoms
        REAL(dbl) :: cdm_ref(3)

! ...   declare other variables
        REAL(dbl)  :: r2, cdm0(3),rcdm0(3), rcdmi(3)

! ...   end of declarations
!  ----------------------------------------------
        CALL ions_cofmass(atoms%taus, atoms%m, atoms%na, atoms%nsp, cdm0)
        CALL s_to_r(cdm0, rcdm0, ht)
        CALL s_to_r(cdm_ref, rcdmi, ht)
        r2 = SUM( (rcdm0(:)-rcdmi(:))**2 )

        WRITE( stdout,1000) R2
1000    FORMAT(/,3X,'Center of mass displacement (a.u.): ',F10.6)

        RETURN
      ENd  SUBROUTINE cdm_displacement


!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE set_velocities(atoms_m, atoms_0, vel_srt, ht_0, delt)

!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

        USE cell_module, ONLY: boxdimensions
        USE cell_base, ONLY: r_to_s, s_to_r
        USE constants, ONLY: angstrom_au

        TYPE (atoms_type) :: atoms_m, atoms_0
        TYPE (boxdimensions) :: ht_0
        REAL(dbl), INTENT(IN) :: delt
        REAL(dbl), INTENT(IN) :: vel_srt( :, : )
        
        INTEGER :: i, ia, nat
        REAL(dbl) :: sv0( 3 )

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
      END SUBROUTINE



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
        REAL(dbl) :: cdm_ref(:), tau_ref(:,:)

! ...   declare other variables
        REAL(dbl) :: sdist(3)
        INTEGER :: isa

        CALL ions_cofmass(atoms%taus, atoms%m, atoms%na, atoms%nsp, cdm_ref)
        DO isa = 1, atoms%nat
          sdist( 1:3 ) = atoms%taus( 1:3 , isa ) - cdm_ref( 1:3 )
          CALL s_to_r( sdist, tau_ref(:,isa), ht )
        END DO

        RETURN
      END subroutine set_reference_positions

!  BEGIN manual -------------------------------------------------------------   

        SUBROUTINE ions_setup(  anne_inp, anner_inp, nconstr_inp, constr_tol_inp, &
           constr_type_inp, constr_dist_inp, constr_inp, gdeltions_inp, tneig, neig_radius)


!  Check the number of atoms and of species. Does something about constrains    
!    (I'll see later what).  
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


! ...     declare modules
          USE constants, ONLY: scmass
          USE io_global, ONLY: ionode
          USE control_flags, ONLY: tdampions
! ...
          IMPLICIT NONE

          LOGICAL, INTENT(IN) :: anne_inp
          REAL(dbl), INTENT(IN) :: anner_inp
          INTEGER, INTENT(IN) :: nconstr_inp
          REAL(dbl), INTENT(IN) :: constr_tol_inp
          INTEGER, INTENT(IN) :: constr_type_inp(:)
          REAL(dbl), INTENT(IN) :: constr_dist_inp(:)
          INTEGER, INTENT(IN) :: constr_inp(:,:)
          REAL(dbl), INTENT(IN) :: gdeltions_inp
          LOGICAL, INTENT(IN) :: tneig
          REAL(dbl), INTENT(IN) :: neig_radius
! ...
! ...     declare other variables
          INTEGER :: is, ia, idd1, idd2, k, istep, isa, ic
          INTEGER :: ierr

! ...     end of declarations
!  ----------------------------------------------
! ...

          IF( tsetup ) THEN
            CALL errore(' ions_setup ', ' module ions already started ', 0)
          END IF

          IF( .NOT. tions_base_init ) THEN
            CALL errore(' ions_setup ', ' ions base module not initializeda ', 0)
          END IF

          tanne = anne_inp
          anner = anner_inp

          tneighbo = tneig
          neighbo_radius = neig_radius

          gdelt = gdeltions_inp
          IF(tdampions .AND. ABS(gdelt+1.0d0) < 1.0d-10 ) THEN
            CALL errore(' ions_setup ', ' gdelt too small ', 1)
          END IF

          if(nax < 1) &
            call errore(' IONS ', ' NAX OUT OF RANGE ',1)
          if(nsp < 1) &
            call errore(' IONS ',' NSP OUT OF RANGE ',1)
          if(nsp > SIZE( na ) ) &
            call errore(' IONS ',' NSP too large, increase NSX parameter ',nsp)

          ALLOCATE( taui(3,nat), STAT=ierr )
          IF( ierr /= 0 ) &
            CALL errore(' ions_setup ',' allocating memory ',ierr)

          taui     = 0.0d0

! ...     Constraints Allocation
          CALL allocate_constrains(constrains, na, constr_inp, &
            constr_dist_inp, constr_type_inp, constr_tol_inp, nconstr_inp)

          tsetup = .TRUE.

          return
        END SUBROUTINE ions_setup

!  BEGIN manual -------------------------------------------------------------   

       SUBROUTINE atoms_init(atoms_m, atoms_0, atoms_p, tau_srt, ind_srt, &
         if_pos, atml, aunits, alat, ht_0)

!  Allocate and fill the three atoms structure using position an cell
!  read from stdin. These initial values could be overwritten in module
!  restart.f90 with the values read from restart file.
!  This subroutine is called from modules init.f90
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

         USE cell_module, ONLY: boxdimensions, s_to_r, r_to_s
         USE constants, ONLY: angstrom_au

         IMPLICIT NONE
!
         TYPE (atoms_type)    :: atoms_0, atoms_p, atoms_m
         TYPE (boxdimensions) :: ht_0
         REAL(dbl), INTENT(IN) :: tau_srt(:,:)
         CHARACTER(LEN=*), INTENT(IN) :: aunits
         CHARACTER(LEN=3), INTENT(IN) :: atml(:)
         REAL(dbl), INTENT(IN) :: alat
         INTEGER, INTENT(IN) :: ind_srt( : )
         INTEGER, INTENT(IN) :: if_pos( :, : )
          
         REAL(dbl), ALLOCATABLE :: stau(:,:)
         LOGICAL  , ALLOCATABLE :: ismb(:,:)
         REAL(dbl) :: p( 3 )
         INTEGER :: ia, isa
         LOGICAL, SAVE :: atoms_init_first = .true.

         IF( .NOT. atoms_init_first ) THEN
           CALL errore(' atoms_init ',' atoms objects already initialized ', 1)
         END IF

         ALLOCATE( stau( 3, nat ), ismb( 3, nat ) )

         DO ia = 1, nat
           CALL r_to_s( tau_srt(:,ia), stau(:,ia), ht_0)
         END DO

         ismb = .TRUE.
         DO isa = 1, nat
           ia = ind_srt( isa )
           IF( if_pos( 1, ia ) == 0 ) ismb( 1, isa ) = .FALSE.
           IF( if_pos( 2, ia ) == 0 ) ismb( 2, isa ) = .FALSE.
           IF( if_pos( 3, ia ) == 0 ) ismb( 3, isa ) = .FALSE.
           ! WRITE(6,*) 'DEBUG atoms_init ', ismb(1:3,isa)
         END DO

         CALL atoms_type_init(atoms_m, stau, ismb, atml, pmass, na, nsp, ht_0%hmat)
         CALL atoms_type_init(atoms_0, stau, ismb, atml, pmass, na, nsp, ht_0%hmat)
         CALL atoms_type_init(atoms_p, stau, ismb, atml, pmass, na, nsp, ht_0%hmat)

         CALL print_scaled_positions( atoms_0, stdout, 'from standard input')

         DEALLOCATE( stau, ismb )

         atoms_init_first = .FALSE.

         RETURN
       END SUBROUTINE atoms_init


!  BEGIN manual -------------------------------------------------------------   

       SUBROUTINE ions_print_info( iunit )

!  Print info about input parameter for ion dynamic
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

         USE io_global, ONLY: ionode
         USE control_flags, ONLY: tranp, amprp, tnosep, tolp, &
             tfor, tsdp, tzerop, tv0rd, taurdr, nv0rd, nbeg, tcp, tcap
         USE ions_base, ONLY: tau_srt, tau_units, if_pos, ind_srt
         USE ions_nose, ONLY: tempw

         integer, intent(in) :: iunit
         integer is, ia, k, ic, isa
         LOGICAL :: ismb( 3 )

          WRITE(iunit,50)
          IF(.NOT.tfor) THEN
            WRITE(iunit,518)
          ELSE
            WRITE(iunit,520)
            IF(tsdp) THEN
              WRITE(iunit,521)
            ELSE
              WRITE(iunit,522)
            END IF
            IF(tzerop) then 
              IF(tv0rd) THEN
                WRITE(iunit,850) nv0rd
              ELSE
                WRITE(iunit,635)
              ENDIF
            ENDIF
          END IF

          DO is = 1, nsp
            IF(tranp(is)) THEN
              WRITE( stdout,510)
              WRITE( stdout,512) is, amprp(is)
            END IF
          END DO

          IF (nbeg < 0 .OR. taurdr) THEN
            WRITE(iunit,660) tau_units
            isa = 0
            DO IS = 1, nsp
              WRITE(iunit,1000) is, na(is), pmass(is)
              DO IA = 1, na(is)
                isa = isa + 1
                WRITE(iunit,1010) ( tau_srt(k,isa), K = 1,3 )
              END DO
            END DO
          ELSE
            WRITE(iunit,661)
          ENDIF

          IF( tfor ) THEN

            IF( ANY( ( if_pos( 1:3, 1:nat ) == 0 )  ) ) THEN

              WRITE(iunit,1020)
              WRITE(iunit,1022)

              DO isa = 1, nat
                ia = ind_srt( isa )
                ismb( 1 ) = ( if_pos(1,ia) /= 0 )
                ismb( 2 ) = ( if_pos(2,ia) /= 0 )
                ismb( 3 ) = ( if_pos(3,ia) /= 0 )
                IF( .NOT. ALL( ismb ) ) THEN
                  WRITE( iunit, 1023 ) isa, ( ismb(k), K = 1, 3 )
                END IF
              END DO

            ELSE

              WRITE(iunit,1021)

            END IF

          END IF


          IF(tfor) THEN
            IF(tcp) THEN
              WRITE( stdout,555) tempw,tolp
            ELSE IF(tcap) THEN
              WRITE( stdout,560) tempw,tolp
            ELSE IF(tnosep) THEN
              WRITE( stdout,595)
            ELSE
              WRITE( stdout,550)
            END IF
          END IF

          IF( constrains%nc  >  0 ) THEN
            IF( ionode ) THEN
              WRITE( stdout,10) constrains%nc, constrains%tolerance
            END IF
            DO ic = 1, constrains%nc
              IF( ionode ) THEN
                WRITE( stdout,30) ic, constrains%tp(ic)%what
                IF( constrains%tp(ic)%what  ==  1 ) THEN
                  WRITE( stdout,40)  constrains%tp(ic)%distance%ia1, &
                      constrains%tp(ic)%distance%is1, &
                      constrains%tp(ic)%distance%ia2, &
                      constrains%tp(ic)%distance%is2
                END IF
              END IF
            END DO
          END IF
    
 1000     FORMAT(3X,'Species ',I3,' atoms = ',I4,' mass (a.u.) = ',F12.2)
 1010     FORMAT(3X,3(1X,F12.6))
 1020     FORMAT(/,3X,'NOT all atoms are allowed to move ')
 1021     FORMAT(/,3X,'All atoms are allowed to move')
 1022     FORMAT(  3X,' indx  ..x.. ..y.. ..z..')
 1023     FORMAT(  3X,I4,3(1X,L5))

 10   FORMAT(/,3X,'Constraints ionic dynamics', /, &
             3X,'--------------------------', /, &
             3X,'number of constraints = ',I3,' tolerance = ',D14.6)
 30   FORMAT(4X,I5,', type ',I2)
 40   FORMAT(4X,'atoms (', I3, ',', I2, ') and (', I3, ',', I2, ')' )
   50 FORMAT(//,3X,'Ions Simulation Parameters',/ &
               ,3X,'--------------------------') 
  115 FORMAT(   3X,'Ewald sum over ',I1,'*',I1,'*',I1,' cells')
  510 FORMAT(   3X,'Initial random displacement of ionic coordinates',/ &
               ,3X,' specie  amplitude')
  512 FORMAT(   3X,I7,2X,F9.6)
  518 FORMAT(   3X,'Ions are not allowed to move')
  520 FORMAT(   3X,'Ions are allowed to move')
  521 FORMAT(   3X,'Ions dynamics with steepest descent')
  522 FORMAT(   3X,'Ions dynamics with newton equations')
  550 FORMAT(   3X,'Ionic temperature is not controlled')
  555 FORMAT(   3X,'Ionic temperature control via ', &
                   'rescaling of velocities :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  560 FORMAT(   3X,'Ionic temperature control via ', &
                   'canonical velocities rescaling :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  595 FORMAT(   3X,'Ionic temperature control via nose thermostat')
  635 FORMAT(   3X,'Zero initial momentum for ions')
  660 FORMAT(   3X,'Ionic position from standard input',/ &
                3X,'( ',A10,' )')
  661 FORMAT(   3X,'Ionic position read from restart file')
  850 FORMAT(   3X,'Initial ion velocities read from unit : ',I4)


         RETURN
       END SUBROUTINE ions_print_info

!  BEGIN manual -------------------------------------------------------------   
!     SUBROUTINE deallocate_ions
!  Deallocate ions input variables
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

        SUBROUTINE deallocate_ions

          INTEGER :: ierr

          IF( .NOT. tsetup ) THEN
            CALL errore(' deallocate_ions ', ' module ions not yet started ', 0)
          END IF
          ierr = 0
          IF( ALLOCATED( taui ) ) DEALLOCATE( taui, STAT = ierr )
          IF( ierr /= 0 ) CALL errore(' deallocate_ions ', ' deallocating memory 1', ierr)
    
          CALL deallocate_constrains(constrains)

          tsetup = .FALSE.

          RETURN
        END SUBROUTINE deallocate_ions


!  BEGIN manual -------------------------------------------------------------   

      REAL(dbl) FUNCTION moveions(TSDP, thdyn, NFI, atoms_m, atoms_0, &
        atoms_p, htm2, htm, ht0, vnosep)

!  Moves the ions
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   



! ... declare modules
      USE cell_module, ONLY: dgcell, r_to_s, s_to_r, boxdimensions
      use control_flags, ONLY: tnosep, tcap, tcp, tdampions
      use time_step, ONLY: delt
      use mp_global, ONLY: mpime

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL, INTENT(IN) :: tsdp, thdyn
      INTEGER, INTENT(IN) :: nfi
      TYPE (boxdimensions), INTENT(IN)    :: htm2, htm
      TYPE (boxdimensions), INTENT(INOUT) :: ht0
      TYPE (atoms_type) :: atoms_m, atoms_0, atoms_p
      REAL(dbl), INTENT(IN) :: vnosep

! ... declare other variables
      REAL(dbl), DIMENSION(3,3) :: annep, svarps, tmat, svarpd, gcm1, gcdot
      REAL(dbl), DIMENSION(3)   :: fions, svel, rvel, stp, tau_diff
      REAL(dbl)                 :: const, dumm, vrnos, gfact, ffact, dt2bym
      REAL(dbl)                 :: fordt2, dt2, delthal
      INTEGER :: i, j, k, l, is, ia, isa


! ...     end of declarations

      dt2     = delt ** 2
      fordt2  = 4.0d0 * dt2
      delthal = 0.5d0 * delt
 
! ...   Determines DGCELL/DT dynamically and GCM1

        IF( thdyn ) THEN
          CALL dgcell( gcm1, gcdot, htm2, htm, ht0, delt )
        END IF

        IF( tnosep ) THEN
          vrnos = vnosep
        ELSE
          vrnos = 0.0d0
        END IF

!....   Steepest descent of ionic degrees of freedom 

        IF( tdampions ) THEN

          gfact = 1.0_dbl / (1.0_dbl + gdelt)
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
                stp(k) = 2.0_dbl * tau_diff(k) + dt2bym * fions(k)
              END DO
              DO k = 1, 3
                IF ( atoms_0%mobile(k,isa) > 0 ) THEN
                  ffact = 0.0_dbl
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

        IF (tanne) THEN
          CALL anneions(atoms_m, atoms_0, atoms_p, ht0)
        END IF

        CALL apply_constraints(ht0, atoms_0, atoms_p, dumm)

        CALL ions_vel( atoms_0%vels, atoms_p%taus, atoms_m%taus, atoms_0%na, atoms_0%nsp, delt )
        CALL ions_kinene( atoms_0%ekint, atoms_0%vels, atoms_0%na, atoms_0%nsp, ht0%hmat, atoms_0%m )
       
        moveions = atoms_0%ekint

      RETURN
      END FUNCTION moveions


!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE constraints_setup(ht,atoms)

!  Fill the constraints structures with the initial positions of atoms 
!  read from stdin
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

        USE cell_module, ONLY: boxdimensions, s_to_r
        USE io_global, ONLY: ionode
        TYPE (boxdimensions), INTENT(IN) :: ht
        TYPE (atoms_type) :: atoms
        REAL(dbl)  :: sr12(3), rr12(3), disto, distn
        INTEGER :: ic, i, is1, ia1, is2, ia2, isa1, isa2
        IF( ionode .AND. (constrains%nc  >  0) ) THEN
          WRITE( stdout,10)
        END IF
 10     FORMAT(/,3X,'Constraints settings')
        DO ic = 1, constrains%nc
          ia1  = constrains%tp(ic)%distance%ia1
          is1  = constrains%tp(ic)%distance%is1
          ia2  = constrains%tp(ic)%distance%ia2
          is2  = constrains%tp(ic)%distance%is2
          isa1 = atoms%isa(is1) + ia1 - 1
          isa2 = atoms%isa(is2) + ia2 - 1
          sr12 = atoms%taus(:,isa2) - atoms%taus(:,isa1)
          CALL s_to_r(sr12, rr12, ht)
          disto = SQRT(rr12(1)**2 + rr12(2)**2 + rr12(3)**2)
          IF( disto  ==  0.0d0 ) THEN
            CALL errore (' constraints setup ', ' zero atomic distance ', ic )
          END IF 
          SELECT CASE (constrains%tp(ic)%what)
            CASE (2)
              distn = constrains%tp(ic)%distance%val 
              atoms%taus(:,isa2) = atoms%taus(:,isa1) + distn * sr12 / disto
            CASE DEFAULT
              constrains%tp(ic)%distance%val = disto
          END SELECT
          IF( ionode ) THEN
            WRITE( stdout,20) ic, ia1, is1, ia2, is2, constrains%tp(ic)%distance%val
          END IF
 20       FORMAT(4X,'c ',I5,' - ',4I4,F12.6)
        END DO
        IF( ionode .AND. (constrains%nc  >  0) ) THEN
          WRITE( stdout,30)
 30     FORMAT(//)
        END IF
        RETURN
      END SUBROUTINE constraints_setup

!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE apply_constraints(ht, atoms_0, atoms_p, wc)

!  Descrive briefly what it does
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   


        USE cell_module, ONLY: boxdimensions, s_to_r, r_to_s

        TYPE (boxdimensions), INTENT(IN) :: ht
        TYPE (atoms_type) :: atoms_0, atoms_p
        REAL(dbl)  :: wc
        LOGICAL :: moved(atoms_0%nax, atoms_0%nsp)
        LOGICAL :: moving(atoms_0%nax,atoms_0%nsp), done
        REAL(dbl)  :: snewa(3,MAX(constrains%nc,1))
        REAL(dbl)  :: snewb(3,MAX(constrains%nc,1))
        REAL(dbl)  :: solda(3,MAX(constrains%nc,1))
        REAL(dbl)  :: soldb(3,MAX(constrains%nc,1))
        REAL(dbl)  :: sab(3), prab(3), rab(3), rd(3), sd(3), rnew(3)
        REAL(dbl)  :: pabsq, rabsq, diffsq, rpab, gab
        REAL(dbl)  :: rma(constrains%nc), rmb(constrains%nc)
        INTEGER, PARAMETER :: itermax = 20
        INTEGER :: ia, ib, is, ic, isa, iter, is1, ia1, is2, ia2
        INTEGER :: isa1, isa2

        wc = 0.0d0
        DO ic = 1, constrains%nc
          ia1  = constrains%tp(ic)%distance%ia1
          is1  = constrains%tp(ic)%distance%is1
          ia2  = constrains%tp(ic)%distance%ia2
          is2  = constrains%tp(ic)%distance%is2
          isa1 = atoms_0%isa(is1) + ia1 - 1
          isa2 = atoms_0%isa(is2) + ia2 - 1
          snewa(:,ic) = atoms_p%taus(:,isa1)
          snewb(:,ic) = atoms_p%taus(:,isa2) 
          solda(:,ic) = atoms_0%taus(:,isa1)
          soldb(:,ic) = atoms_0%taus(:,isa2)
          rma(ic)     = 1.0d0 / atoms_0%m(is1)
          rmb(ic)     = 1.0d0 / atoms_0%m(is2)
        END DO

        moving = .FALSE.
        moved  = .TRUE.
        done   = .FALSE.
        iter   = 1
        DO WHILE(( .NOT. done ) .AND. ( iter  <=  itermax ))
          done = .TRUE.
          DO ic = 1, constrains%nc
            ia1  = constrains%tp(ic)%distance%ia1
            is1  = constrains%tp(ic)%distance%is1
            ia2  = constrains%tp(ic)%distance%ia2
            is2  = constrains%tp(ic)%distance%is2
            IF( moved(ia1,is1) .OR. moved(ia2,is2) ) THEN
              sab(1) = snewa(1,ic) - snewb(1,ic)
              sab(2) = snewa(2,ic) - snewb(2,ic)
              sab(3) = snewa(3,ic) - snewb(3,ic)
              ! CALL pbcs(sab(1),sab(2),sab(3),sab(1),sab(2),sab(3),1)
              CALL s_to_r(sab,prab,ht)
              pabsq  = prab(1)**2 + prab(2)**2 + prab(3)**2
              rabsq  = constrains%tp(ic)%distance%val
              diffsq = rabsq - pabsq
              IF( ABS(diffsq)  >  (rabsq * constrains%tolerance) ) THEN
                sab(1) = solda(1,ic) - soldb(1,ic)
                sab(2) = solda(2,ic) - soldb(2,ic)
                sab(3) = solda(3,ic) - soldb(3,ic)
                ! CALL pbcs(sab(1),sab(2),sab(3),sab(1),sab(2),sab(3),1)
                CALL s_to_r(sab,rab,ht)
                rpab = rab(1) * prab(1) + rab(2) * prab(2) + rab(3) * prab(3)
                ! WRITE( stdout,*) diffsq,  rpab 
                gab = diffsq / ( 2.0d0 * ( rma(ic) + rmb(ic) ) * rpab )
                wc = wc + gab * rabsq
                rd(:) = rab(:) * gab
                CALL r_to_s(rd,sd,ht)
                snewa(:,ic) = snewa(:,ic) + rma(ic) * sd(:)
                snewb(:,ic) = snewb(:,ic) - rmb(ic) * sd(:)
                
                moving(ia1,is1) = .TRUE.
                moving(ia2,is2) = .TRUE.
                done = .FALSE.
              END IF
            END IF
          END DO
          moved = moving
          moving = .FALSE.
          iter = iter + 1
        END DO
        DO ic = 1, constrains%nc
          ia1  = constrains%tp(ic)%distance%ia1
          is1  = constrains%tp(ic)%distance%is1
          ia2  = constrains%tp(ic)%distance%ia2
          is2  = constrains%tp(ic)%distance%is2
          isa1 = atoms_p%isa(is1) + ia1 - 1
          isa2 = atoms_p%isa(is2) + ia2 - 1
          atoms_p%taus(:,isa1) = snewa(:,ic)
          atoms_p%taus(:,isa2) = snewb(:,ic)
        END DO
        RETURN
      END SUBROUTINE apply_constraints


!  BEGIN manual -------------------------------------------------------------   

      SUBROUTINE update_ions(atoms_m, atoms_0, atoms_p)

!  Update ionic positions and velocities in atoms structures
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

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
      REAL(dbl), INTENT(IN) :: delt
!
! ... Locals
!
      REAL(dbl) :: tempp, temp1, temp2
      REAL(dbl), ALLOCATABLE :: fions(:,:)
      INTEGER :: isa
!
! ... Subroutine body
!
        IF(MOD(nfi,icapstep) /= 0) THEN
          RETURN
        END IF

        tempp   = atoms_0%ekint * factem
        tempp   = tempp * 2.0_dbl / atoms_0%doft
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

      REAL(dbl) FUNCTION max_ion_forces( atoms )

!  Descrive briefly what it does
!  --------------------------------------------------------------------------   
!  END manual ---------------------------------------------------------------   

        IMPLICIT NONE 
        TYPE (atoms_type) :: atoms
        INTEGER :: ia
        REAL(dbl) :: fmax
        fmax = 0.0d0
        DO ia = 1, atoms%nat
          IF( atoms%mobile(1, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(1, ia) ) )
          IF( atoms%mobile(2, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(2, ia) ) )
          IF( atoms%mobile(3, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(3, ia) ) )
        END DO
        max_ion_forces = fmax
        RETURN
      END FUNCTION

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

        REAL(dbl), INTENT(OUT) :: pos(:,:), fion(:,:)
        TYPE (atoms_type), INTENT(IN) :: atoms
        TYPE (boxdimensions), INTENT(IN) :: ht
        INTEGER, INTENT(IN) :: isrt( : )
        INTEGER :: ia, is, isa, ipos

        IF( .NOT. tsetup ) THEN 
          CALL errore( ' resort_position ', ' ions module not initialized ', 1 )
        END IF

        isa = 0
        DO is = 1, atoms%nsp
          DO ia = 1, atoms%na(is)
            isa  = isa + 1
            ipos = isrt( isa )
            !WRITE( 6, * ) ' DEBUG resort_posi ipos ', ipos
            !WRITE( 6, * ) atoms%for(1,isa), atoms%for(2,isa), atoms%for(3,isa)
            !WRITE( 6, * ) ' DEBUG resort_posi isa ', isa
            CALL s_to_r( atoms%taus( : , isa ), pos( :, ipos ), ht )
            fion( :, ipos ) = atoms%for( : , isa )
          END DO
        END DO

        RETURN
      END SUBROUTINE
     


   END MODULE ions_module
