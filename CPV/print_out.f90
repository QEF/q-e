!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
   MODULE print_out_module
!=----------------------------------------------------------------------------=!

        USE kinds
        USE io_global, ONLY: ionode, ionode_id, stdout
        USE io_files, ONLY: sfacunit, sfac_file, opt_unit

        IMPLICIT NONE
        SAVE

        PRIVATE
        
        REAL(dbl) :: old_clock_value = -1.0d0

        REAL(dbl) :: timeform = 0.0d0, &
                     timernl = 0.0d0,  &
                     timerho = 0.0d0,  &
                     timevof = 0.0d0,  &
                     timerd = 0.0d0
        REAL(dbl) :: timeorto = 0.0d0, timeloop = 0.0d0
        INTEGER   :: timecnt = 0

        REAL(dbl), EXTERNAL :: cclock


        PUBLIC :: printout, printmain
        PUBLIC :: print_time, print_sfac, printacc, print_legend
        PUBLIC :: cp_print_rho

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!


   SUBROUTINE printout(nfi, atoms, ekinc, ekcell, tprint, toptical, ht, kp, &
                avgs, avgs_run, edft) 

      USE cell_module, only: s_to_r, boxdimensions, press
      use constants, only: factem, au_gpa, au, uma_si, bohr_radius_cm, scmass
      use energies, only: print_energies, dft_energy_type
      use mp_global, only: mpime
      use electrons_module, only: ei, ei_emp, n_emp
      use brillouin, only: kpoints
      use time_step, ONLY: tps
      USE electrons_nose, ONLY: electrons_nose_nrg, xnhe0, vnhe, qne, ekincw
      USE ions_module, ONLY: displacement, cdm_displacement
      USE polarization, ONLY: pdipole, pdipolt, p
      USE optical_properties, ONLY:  WRITE_DIELEC
      USE control_flags, ONLY: tdipole, tnosee, tnosep, tnoseh, iprsta
      USE atoms_type_module, ONLY: atoms_type
      USE sic_module, ONLY: ind_localisation, pos_localisation, nat_localisation, self_interaction
      USE sic_module, ONLY: rad_localisation
      USE ions_base, ONLY: ions_temp, cdmi, taui
      USE ions_nose, ONLY: ndega, ions_nose_nrg, xnhp0, vnhp, qnp, gkbt, kbt, nhpcl, nhpdim, atm2nhp, ekin2nhp, gkbt2nhp
      USE cell_nose, ONLY: cell_nose_nrg, qnh, temph, xnhh0, vnhh
      USE cell_base, ONLY: iforceh
      USE printout_base, ONLY: printout_base_open, printout_base_close, &
            printout_pos, printout_cell, printout_stress
      USE environment, ONLY: start_cclock_val

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nfi
      TYPE (atoms_type) :: atoms
      REAL(dbl) ::   ekinc, ekcell
      LOGICAL   ::   tprint, toptical, ttsic
      type (boxdimensions), intent(in) :: ht
      type (kpoints), intent(in) :: kp
      REAL(dbl) :: avgs(:), avgs_run(:)
      TYPE (dft_energy_type) :: edft
!
! ...
      INTEGER is, ia, k, i, j, ik, isa, iunit, nfill, nempt
      REAL(dbl) :: tau(3), vel(3), stress_tensor(3,3), temps( atoms%nsp )
      REAL(dbl) :: tempp, econs, ettt, out_press, ekinpr, enosee
      REAL(dbl) :: enthal, totalmass, enoseh, temphc, enosep
      REAL(dbl) :: dis(atoms%nsp), h(3,3)
      LOGICAL :: tfile, topen
      INTEGER, SAVE :: nfi_old = -999
      LOGICAL, SAVE :: first = .true.
      CHARACTER(LEN=4), ALLOCATABLE :: labelw( : )
      REAL(dbl), ALLOCATABLE :: tauw( :, : )

! ...   Subroutine Body

      tfile = ( tprint .AND. ( nfi > 0 ) )
      ttsic = ( self_interaction /= 0 )
      
! ...   Calculate Ions temperature tempp (in Kelvin )

      CALL ions_temp( tempp, temps, ekinpr, atoms%vels, atoms%na, atoms%nsp, ht%hmat, atoms%m, ndega, nhpdim, atm2nhp, ekin2nhp )

! ...   Calculate MSD for each specie, starting from taui positions

      CALL displacement(dis, atoms, taui, ht)

! ...   Stress tensor (in GPa) and pressure (in GPa)

      stress_tensor = MATMUL( ht%pail(:,:), ht%a(:,:) ) * au_gpa / ht%deth
      out_press = ( stress_tensor(1,1) + stress_tensor(2,2) + stress_tensor(3,3) ) / 3.0d0

! ...   Enthalpy (in Hartree)

      enthal = edft%etot + press * ht%deth

      IF( tnoseh ) THEN
        enoseh = cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
      ELSE
        enoseh = 0.0d0
      END IF

      if( COUNT( iforceh == 1 ) > 0 ) then
         temphc = 2.0d0 * factem * ekcell / DBLE( COUNT( iforceh == 1 ) )
      else
         temphc = 0.0d0
      endif

      IF( tnosep ) THEN
        enosep = ions_nose_nrg( xnhp0, vnhp, qnp, gkbt2nhp, kbt, nhpcl, nhpdim )
      ELSE
        enosep = 0
      END IF

      IF( tnosee ) THEN
        enosee = electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
      ELSE
        enosee = 0
      END IF


! ...   Energy expectation value for physical ions hamiltonian
! ...   in Born-Oppenheimer approximation

      econs  = atoms%ekint + ekcell + enthal

! ...   Car-Parrinello constant of motion

      ettt   = econs + ekinc + enosee + enosep + enoseh

      IF( nfi > 0 ) THEN

! ...   sum up values to be averaged

        avgs(1) = avgs(1) + ekinc
        avgs(2) = avgs(2)
        avgs(3) = avgs(3)
        avgs(4) = avgs(4) + edft%etot
        avgs(5) = avgs(5) + tempp
        avgs(6) = avgs(6) + enthal
        avgs(7) = avgs(7) + econs
        avgs(8) = avgs(8) + out_press  ! pressure in GPa
        avgs(9) = avgs(9) + ht%deth

! ...   sum up values to be averaged

        avgs_run(1) = avgs_run(1) + ekinc
        avgs_run(2) = avgs_run(2)
        avgs_run(3) = avgs_run(3)
        avgs_run(4) = avgs_run(4) + edft%etot
        avgs_run(5) = avgs_run(5) + tempp
        avgs_run(6) = avgs_run(6) + enthal
        avgs_run(7) = avgs_run(7) + econs
        avgs_run(8) = avgs_run(8) + out_press  ! pressure in GPa
        avgs_run(9) = avgs_run(9) + ht%deth

      END IF

! ...   Check Memory

      IF( iprsta > 1 ) CALL memstat(mpime)

! ...   Print physical variables to fortran units

      IF ( ionode ) THEN

        IF ( tprint .OR. first .OR. ( nfi < 0 ) ) THEN

          ! ...  Write total energy components to standard output

          WRITE( stdout, 5 ) tps
          CALL print_energies( ttsic, edft )

          IF( ttsic ) THEN
            WRITE ( stdout, *) '  Sic_correction: type ', self_interaction
            WRITE ( stdout, *) '  Localisation: ', rad_localisation
            DO i = 1, nat_localisation
              write( stdout, *) '    Atom ', ind_localisation(i), ' : ', pos_localisation(4,i)
            END DO
          END IF

        END IF

        IF ( tprint .OR. ( nfi < 0 ) ) THEN

          IF( tfile ) THEN
            ! ...  Open units 30, 31, ... 40 for simulation output 
            CALL printout_base_open()
            !
          END IF

          ! ...  Write Dielectric tensor and electronic conductivity on unit fort.30

          IF ( toptical .AND. tfile ) THEN
            CALL write_dielec( nfi, tps )
          END IF

          ! ...  Write Polarizability tensor to stdout and fortran unit 32

          IF ( tdipole ) THEN
            WRITE( stdout,19) 
            WRITE( stdout,20) (pdipole(i),i=1,3)
            WRITE( stdout,20) (p(i),i=1,3)
            WRITE( stdout,20) (pdipolt(i),i=1,3)
            IF ( tfile ) THEN 
               WRITE(32,30) nfi, tps
               WRITE(32,20) (pdipole(i),i=1,3)
               WRITE(32,20) (p(i),i=1,3)
               WRITE(32,20) (pdipolt(i),i=1,3)
            END IF
          END IF
          
          ! ...  Write energies, pressure, volume and MSD to unit fort.33

          IF( tfile ) THEN
            WRITE(33,2000) nfi, ekinc, temphc, tempp, edft%etot, enthal, econs, ettt, &
              tps, ht%deth, out_press, (dis(i),i=1,atoms%nsp)
          END IF
 2000 FORMAT(I7, 1X,  F8.5, 1X, F6.1,  1X, F6.1, 1X, F12.5, 1X, F12.5, &
                 1X, F12.5, 1X, F12.5, 1X, F9.2, 1X, F9.2,  1X, F9.2,  &
                 1X, 100F9.4)

          ! ...  Write Positions and velocities to stdout and unit fort.35

          ALLOCATE( labelw( atoms%nat ) )
          ALLOCATE( tauw( 3, atoms%nat ) )
          DO is = 1, atoms%nsp
            DO ia = atoms%isa(is), atoms%isa(is) + atoms%na(is) - 1
              labelw( ia ) = atoms%label(is)
              CALL s_to_r( atoms%taus(:,ia), tauw(:,ia), ht )
            END DO
          END DO
          
          WRITE( stdout,11) 
          CALL printout_pos( stdout, nfi, tauw     , atoms%nat, tps, labelw )
          IF( tfile ) CALL printout_pos( 35, nfi, tauw     , atoms%nat, tps )

          DO ia = 1, atoms%nat
            CALL s_to_r( atoms%vels(:,ia), tauw(:,ia), ht )
          END DO
          WRITE( stdout,12) 
          CALL printout_pos( stdout, nfi, tauw     , atoms%nat, tps, labelw )
          IF( tfile ) CALL printout_pos( 34, nfi, tauw     , atoms%nat, tps )

          WRITE( stdout,13) 
          CALL printout_pos( stdout, nfi, atoms%for, atoms%nat, tps, labelw )
          IF( tfile ) CALL printout_pos( 37, nfi, atoms%for, atoms%nat, tps )

          DEALLOCATE( labelw )
          DEALLOCATE( tauw )

          ! ...  Write to the standard output the center of mass displacement

          CALL cdm_displacement(cdmi, atoms, ht)

          ! ...  Write Cell parameter to unit fort.36 and stdout
          ! ...  Write Stress tensor to unit fort.38 and stdout

          WRITE( stdout, 10 )
          CALL printout_cell( stdout, nfi, ht%a, tps )
          WRITE( stdout, 17 )
          CALL printout_stress( stdout, nfi, stress_tensor, tps )
          IF( tfile ) THEN
            CALL printout_cell( 36, nfi, ht%a, tps )
            CALL printout_stress( 38, nfi, stress_tensor, tps )
          END IF


          ! ...  System density:

          totalmass = 0.0
          DO is = 1, atoms%nsp
            totalmass = totalmass + atoms%m(is) * atoms%na(is) / scmass
          END DO
          WRITE( stdout, fmt='(/,3X,"System Density [g/cm^3] : ",F10.4)' ) &
            totalmass / ht%deth * 11.2061 ! UMA_SI * 1000.0 / BOHR_RADIUS_CM**3 

          ! ...       Write eigenvalues to stdout and unit fort.31

          IF( tprint ) THEN
            !
            nfill = SIZE(ei,1)
            nempt = n_emp
            IF ( tfile ) THEN
              WRITE(31,30) nfi, tps
              WRITE(31,1030) nfill, nempt,  SIZE(ei,2),  SIZE(ei,3)
            END IF
            DO ik = 1, kp%nkpt
              DO j = 1, SIZE( ei, 3 )
                WRITE( stdout,1002) ik, j
                WRITE( stdout,1004) ( ei( i, ik, j ) * au, i = 1, SIZE( ei, 1 ) )
                IF( nempt .GT. 0 ) THEN
                  WRITE( stdout,1005) ik, j
                  WRITE( stdout,1004) ( ei_emp( i, ik, j ) * au , i = 1, SIZE( ei_emp, 1 ) )
                  WRITE( stdout,1006) ( ei_emp( 1, ik, j ) - ei( SIZE(ei,1), ik, j ) ) * au
                END IF
                IF( tfile ) THEN
                  WRITE(31,1010) IK, j
                  WRITE(31,1020) ( ei( i, ik, j ) * au, i = 1, SIZE( ei, 1 ) )
                  IF( nempt .GT. 0 ) THEN
                    WRITE(31,1011) ik, j
                    WRITE(31,1020) ( ei_emp( i, ik, j ) * au , i = 1, SIZE( ei_emp, 1 ) )
                    WRITE(31,1021) ( ei_emp( 1, ik, j ) - ei( SIZE(ei,1), ik, j ) ) * au
                  END IF
                END IF
              END DO
            END DO
            !
          END IF

! ...       Write partial temperature and MSD for each atomic specie tu stdout

          WRITE( stdout, 1944 )
          DO is = 1, atoms%nsp
             WRITE( stdout, 1945 ) is, temps(is), dis(is)
          END DO

          IF( tfile .AND. ( tnosee .OR. tnosep ) ) THEN
            IF(tfile) WRITE(39,*) nfi, enosep, enosee
          END IF

          IF( tfile ) THEN
            ! ...   Close and flush unit 30, ... 40
            CALL printout_base_close()
            !
          END IF

        END IF

        IF( nfi >= 0 ) THEN

          ! ...  Print energies on standard output EVERY MD STEP!  

          IF( tprint .or. first ) WRITE( stdout, 1947 )
          WRITE( stdout, 1948) nfi, ekinc, temphc, tempp, edft%etot, enthal, econs, ettt

        END IF

 1947 FORMAT(//3X,'nfi', 5X,'ekinc', 2X,'temph', 2X,'tempp', 9X,'etot', 7X,'enthal', &
               8X,'econs', 9X,'ettt')
 1946 FORMAT(//6X,       5X,'(AU) ', 3X,'(K )', 9X,'(AU)', &
               2X,' (K) ', 8X,' (AU)', 9X,'(AU)')
 1948 FORMAT(        I6, 1X, F9.5,  1X, F6.1,  1X, F6.1, 1X, F12.5, 1X, F12.5, &
                         1X, F12.5, 1X, F12.5 )

   5  FORMAT(/,3X,'Simulated Time (ps): ',F12.6)
  10  FORMAT(/,3X,'Cell Variables (AU)',/)
  11  FORMAT(/,3X,'Atomic Positions (AU)',/)
  12  FORMAT(/,3X,'Atomic Velocities (AU)',/)
  13  FORMAT(/,3X,'Atomic Forces (AU)',/)
  17  FORMAT(/,3X,'Total Stress (GPa)',/)
  19  FORMAT(/,3X,'Dipole moment (AU)',/)
  20  FORMAT(6X,3(F18.8,2X))
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
  50  FORMAT(6X,3(F18.8,2X))
  293 FORMAT(/,3X,'Atomic Coordinates (AU):')
  253 FORMAT(3F12.5)
  252 FORMAT(3E14.6)
  254 FORMAT(3F14.8)
  255 FORMAT(3X,A3,3F10.4,3E12.4)
 100  FORMAT(3X,A3,3(1X,E14.6))
 101  FORMAT(3X,3(1X,E14.6))
 1002 FORMAT(/,3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1005 FORMAT(/,3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1004 FORMAT(10F8.2)
 1006 FORMAT(/,3X,'Electronic Gap (eV) = ',F8.2,/)
 1010 FORMAT(3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1011 FORMAT(3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1020 FORMAT(10F8.2)
 1021 FORMAT(3X,'Electronic Gap (eV) = ',F8.2)
 1030 FORMAT(3X,'nfill = ', I4, ', nempt = ', I4, ', kp = ', I3, ', spin = ',I2)
 1944 FORMAT(//'   Partial temperatures (for each ionic specie) ', &
             /,'   Species  Temp (K)   MSD (AU)')
 1945 FORMAT(3X,I6,1X,F10.2,1X,F10.4)

      END IF

      IF( nfi >= 0 ) nfi_old = nfi
      first   = .FALSE.

    RETURN
  END SUBROUTINE printout

!=----------------------------------------------------------------------------=!


  SUBROUTINE printmain( tbeg, taurdr, atoms )

      use ions_module, only: print_scaled_positions
      USE atoms_type_module, ONLY: atoms_type

      implicit none

      logical TBEG,TAURDR
      TYPE (atoms_type) :: atoms

      if(ionode) then
        WRITE( stdout,*)
        WRITE( stdout,*) '  ===> MAIN (FROM FILE NDR) <==='
        WRITE( stdout,*)

        IF(.NOT.TBEG) THEN
          WRITE( stdout,*) '  Initial cell (HT0) from restart file NDR'
        ELSE
          WRITE( stdout,*) '  Initial cell (HT0) from input file'
        END IF

        WRITE( stdout,*)
        IF(TAURDR) THEN
          WRITE( stdout,10) 
        ELSE
          WRITE( stdout,9) 
        END IF
        call print_scaled_positions(atoms,  stdout )

      end if

 35   FORMAT(3(1X,F6.2),6X,F6.4,14X,I5)
  9   FORMAT('   SCALED ATOMIC COORDINATES (FROM NDR)  : ')
 10   FORMAT('   SCALED ATOMIC COORDINATES (FROM STDIN): ')
555   FORMAT(3(4X,3(1X,F8.4)))
600   FORMAT(4X,I3,3(2X,F12.8))

    RETURN
  END SUBROUTINE printmain


!=----------------------------------------------------------------------------=!


  SUBROUTINE print_legend( )
    IMPLICIT NONE
      IF( ionode ) THEN
        WRITE( stdout, *) 
        WRITE( stdout, *) '  Short Legend and Physical Units in the Output'
        WRITE( stdout, *) '  ---------------------------------------------'
        WRITE( stdout, *) '  NFI    [int]  - step index'
        WRITE( stdout, *) '  EKINC  [A.U.] - kinetic energy of the fictitious electronic dynamics'
        WRITE( stdout, *) '  TEMPH  [K]    - Temperature of the fictitious cell dynamics'
        WRITE( stdout, *) '  TEMP   [K]    - Ionic temperature'
        WRITE( stdout, *) '  ETOT   [A.U.] - Scf total energy (Kohn-Sham hamiltonian)'
        WRITE( stdout, *) '  ENTHAL [A.U.] - Enthalpy ( ETOT + P * V )'
        WRITE( stdout, *) '  ECONS  [A.U.] - Enthalpy + kinetic energy of ions and cell'
        WRITE( stdout, *) '  ETTT   [A.U.] - Constant of motion for the CP lagrangian'
        WRITE( stdout, *) 
      END IF
    RETURN
  END SUBROUTINE print_legend


!=----------------------------------------------------------------------------=!


    SUBROUTINE print_time( tprint, texit, timeform_ , timernl_ , timerho_ , &
      timevof_ , timerd_ , timeorto_ , timeloop_ , timing )

      USE fft, ONLY: fft_time_stat 
      USE orthogonalize, ONLY: print_ortho_time
      USE potentials, ONLY: print_vofrho_time
      USE stress, ONLY: print_stress_time
      USE printout_base, ONLY: pprefix

      IMPLICIT NONE

      REAL(dbl) :: timeform_ , timernl_ , timerho_ , timevof_ , timerd_
      REAL(dbl) :: timeorto_ , timeloop_
      LOGICAL, INTENT(IN) :: timing, tprint, texit

      REAL(dbl)  :: timeav
      REAL(dbl), SAVE :: timesum = 0.0d0
      INTEGER, SAVE   :: index   = 0
      CHARACTER(LEN=256) :: file_name

      timeform = timeform + timeform_
      timernl  = timernl + timernl_
      timerho  = timerho + timerho_
      timevof  = timevof + timevof_
      timerd   = timerd + timerd_
      timeorto = timeorto + timeorto_
      timeloop = timeloop + timeloop_
      timecnt = timecnt + 1

      ! IF( timing .AND. ( tprint .OR. texit ) ) THEN
      IF( timing ) THEN

        IF( ionode ) THEN

          file_name = trim(pprefix)//'.tmo'
          CALL open_and_append( opt_unit, file_name )
          IF( index == 0 ) WRITE( opt_unit, 930 ) 
          CALL print_ortho_time( opt_unit )
          CLOSE( opt_unit )
 930      FORMAT('     RHOSET    SIGSET      DIAG     TRASF      ITER  BACKTRAS       TOT')

          file_name = trim(pprefix)//'.tmv'
          CALL open_and_append( opt_unit, file_name )
          IF( index == 0 ) WRITE( opt_unit, 940 ) 
          CALL print_vofrho_time( opt_unit )
          CLOSE( opt_unit )
 940      FORMAT('        ESR     FWFFT        XC       HAR    INVFFT    STRESS       TOT')

          file_name = trim(pprefix)//'.tms'
          CALL open_and_append( opt_unit, file_name )
          IF( index == 0 ) WRITE( opt_unit, 920 ) 
          CALL print_stress_time( opt_unit )
          CLOSE( opt_unit )
 920      FORMAT('        EK        EXC       ESR        EH        EL       ENL       TOT')

          file_name = trim(pprefix)//'.tml'
          CALL open_and_append( opt_unit, file_name )
          IF( index == 0 ) WRITE( opt_unit, 910 ) 
          IF( timecnt > 0 ) THEN
            timeform = timeform / timecnt
            timernl  = timernl  / timecnt
            timerho  = timerho  / timecnt
            timevof  = timevof  / timecnt
            timerd   = timerd   / timecnt
            timeorto = timeorto / timecnt
            timeloop = timeloop / timecnt
            WRITE( opt_unit, 999 ) TIMEFORM, TIMERNL, TIMERHO, TIMEVOF, TIMERD, TIMEORTO, TIMELOOP
 999        FORMAT(7(F9.3))
          END IF
          CLOSE( opt_unit )
 910      FORMAT('     FORM     NLRH   RHOOFR   VOFRHO    FORCE    ORTHO     LOOP')

          index   = index + 1
          timesum = timesum + timeloop

        END IF

        timeform = 0.0d0
        timernl  = 0.0d0
        timerho  = 0.0d0
        timevof  = 0.0d0
        timerd   = 0.0d0
        timeorto = 0.0d0
        timeloop = 0.0d0
        timecnt = 0

      END IF

      IF ( texit ) THEN
        IF( index > 0 ) timeav  = timesum / REAL( MAX( index, 1 ) )
        IF (ionode) THEN
          WRITE( stdout,*)
          WRITE( stdout, fmt='(3X,"Execution time statistics (SEC)")')
          WRITE( stdout, fmt='(3X,"Mean time for MD step ..",F12.3)') timeav
          CALL fft_time_stat( index )
        END IF
      ENDIF

      RETURN
    END SUBROUTINE print_time


!=----------------------------------------------------------------------------=!


    SUBROUTINE print_sfac( rhoe, desc, sfac )

      USE mp_global, ONLY: mpime, nproc, group
      USE mp, ONLY: mp_max, mp_get, mp_put
      USE fft, ONLY : pfwfft, pinvfft
      USE charge_types, ONLY: charge_descriptor
      USE reciprocal_vectors, ONLY: ig_l2g, gx, g
      USE gvecp, ONLY: ngm

      TYPE (charge_descriptor), INTENT(IN) :: desc
      REAL(dbl), INTENT(IN) :: rhoe(:,:,:,:)
      COMPLEX(dbl), INTENT(IN) ::  sfac(:,:)

      INTEGER :: nspin, ispin, ip, nsp, ngx_l, ng, is, ig
      COMPLEX(dbl), ALLOCATABLE :: rhoeg(:,:)
      COMPLEX(dbl), ALLOCATABLE :: rhoeg_rcv(:,:)
      REAL   (dbl), ALLOCATABLE :: hg_rcv(:)
      REAL   (dbl), ALLOCATABLE :: gx_rcv(:,:)
      INTEGER     , ALLOCATABLE :: ig_rcv(:)
      COMPLEX(dbl), ALLOCATABLE :: sfac_rcv(:,:)

        nspin = SIZE(rhoe,4)
        nsp   = SIZE(sfac,2)
        ngx_l = ngm
        CALL mp_max(ngx_l, group)
        ALLOCATE(rhoeg(ngm,nspin))
        ALLOCATE(hg_rcv(ngx_l))
        ALLOCATE(gx_rcv(3,ngx_l))
        ALLOCATE(ig_rcv(ngx_l))
        ALLOCATE(rhoeg_rcv(ngx_l,nspin))
        ALLOCATE(sfac_rcv(ngx_l,nsp))
! ...   FFT: rho(r) --> rho(g)
        DO ispin = 1, nspin
          CALL pfwfft(rhoeg(:,ispin),rhoe(:,:,:,ispin))
        END DO
        IF( ionode ) THEN
          OPEN(sfacunit, FILE=TRIM(sfac_file), STATUS='UNKNOWN')
        END IF

        DO ip = 1, nproc
          CALL mp_get(ng, ngm, mpime, ionode_id, ip-1, ip)
          CALL mp_get(hg_rcv(:), g(:), mpime, ionode_id, ip-1, ip)
          CALL mp_get(gx_rcv(:,:), gx(:,:), mpime, ionode_id, ip-1, ip)
          CALL mp_get(ig_rcv(:), ig_l2g(:), mpime, ionode_id, ip-1, ip)
          DO ispin = 1, nspin
            CALL mp_get( rhoeg_rcv(:,ispin), rhoeg(:,ispin), mpime, ionode_id, ip-1, ip)
          END DO
          DO is = 1, nsp
            CALL mp_get( sfac_rcv(:,is), sfac(:,is), mpime, ionode_id, ip-1, ip)
          END DO
          IF( ionode ) THEN
            DO ig = 1, ng
              WRITE(sfacunit,100) ig_rcv(ig), &
                hg_rcv(ig), gx_rcv(1,ig), gx_rcv(2,ig), gx_rcv(3,ig), &
                (sfac_rcv(ig,is),is=1,nsp), &
                (rhoeg_rcv(ig,ispin),ispin=1,nspin)
            END DO
          END IF
 100      FORMAT(1X,I8,1X,4(1X,D13.6),1X,50(2X,2D14.6))
        END DO

        IF( ionode ) THEN
          CLOSE(sfacunit)
        END IF

        DEALLOCATE(rhoeg)
        DEALLOCATE(hg_rcv)
        DEALLOCATE(gx_rcv)
        DEALLOCATE(ig_rcv)
        DEALLOCATE(rhoeg_rcv)
        DEALLOCATE(sfac_rcv)

      RETURN
    END SUBROUTINE print_sfac


!=----------------------------------------------------------------------------=!


    SUBROUTINE printacc( nfi, rhoe, desc, rhoout, atoms, ht, nstep_run, avgs, avgs_run )

      USE charge_density, ONLY: printrho
      USE cell_module, ONLY: boxdimensions
      USE atoms_type_module, ONLY: atoms_type
      USE charge_types, ONLY: charge_descriptor

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nfi, nstep_run
      LOGICAL, INTENT(IN) :: rhoout
      REAL(dbl), intent(in) :: rhoe(:,:,:,:)
      TYPE (charge_descriptor), intent(in) :: desc
      REAL (dbl) :: avgs(:), avgs_run(:)
      TYPE (atoms_type) :: atoms
      TYPE (boxdimensions), intent(in) :: ht
 
      IF ( nfi < 1 ) THEN
        RETURN
      END IF

      IF( rhoout ) THEN
         CALL printrho(nfi, rhoe, desc, atoms, ht)
      END IF

      avgs     = avgs     / REAL( nfi )
      avgs_run = avgs_run / REAL( nstep_run )

      IF( ionode ) THEN
        WRITE( stdout,1949)
        WRITE( stdout,1951) avgs(1), avgs_run(1)
        WRITE( stdout,1954) avgs(4), avgs_run(4)
        WRITE( stdout,1955) avgs(5), avgs_run(5)
        WRITE( stdout,1956) avgs(6), avgs_run(6)
        WRITE( stdout,1957) avgs(7), avgs_run(7)
        WRITE( stdout,1958) avgs(8), avgs_run(8)
        WRITE( stdout,1959) avgs(9), avgs_run(9)
        WRITE( stdout,1990)
 1949   FORMAT(//,3X,'Averaged Physical Quantities',/ &
                ,3X,'             ',' accomulated','    this run')
 1951   FORMAT(3X,'EKINC        ',F12.5,F12.5,' (AU)')
 1954   FORMAT(3X,'TOTEL ENERGY ',F12.5,F12.5,' (AU)')
 1955   FORMAT(3X,'TEMPERATURE  ',F12.5,F12.5,' (K )')
 1956   FORMAT(3X,'ENTHALPY     ',F12.5,F12.5,' (AU)')
 1957   FORMAT(3X,'ECONS        ',F12.5,F12.5,' (AU)')
 1958   FORMAT(3X,'PRESSURE     ',F12.5,F12.5,' (Gpa)')
 1959   FORMAT(3X,'VOLUME       ',F12.5,F12.5,' (AU)')
 1990   FORMAT(/)
      END IF

      RETURN
    END SUBROUTINE printacc

!=----------------------------------------------------------------------------=!


    SUBROUTINE open_and_append( iunit, file_name )
      USE io_global, ONLY: ionode
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iunit
      CHARACTER(LEN = *), INTENT(IN) :: file_name
      INTEGER :: ierr
      IF( ionode ) THEN
        OPEN( UNIT = iunit, FILE = trim( file_name ), &
          STATUS = 'unknown', POSITION = 'append', IOSTAT = ierr)
        IF( ierr /= 0 ) &
          CALL errore( ' open_and_append ', ' opening file '//trim(file_name), 1 )
      END IF
      RETURN
    END SUBROUTINE open_and_append

!=----------------------------------------------------------------------------=!

   SUBROUTINE cp_print_rho(nfi, bec, c0, eigr, irb, eigrb, rhor, rhog, rhos, lambdap, lambda, tau0, h )
   
     use kinds, only: dbl
     use ensemble_dft, only: tens, ismear, z0, c0diag, becdiag, dval, zaux, e0, zx
     use electrons_base, only: nx => nbspx, n => nbsp, ispin => fspin, f, nspin
     use electrons_base, only: nel, iupdwn, nupdwn, nudx, nelt
     use energies, only: enl, ekin
     use ions_base, only: nsp
     use uspp, only: rhovan => becsum
     use grid_dimensions, only: nnr => nnrx
     use io_global, only: stdout
     USE control_flags, ONLY: printwfc, trhor

     IMPLICIT NONE

     INTEGER :: nfi
     INTEGER :: irb(:,:)
     COMPLEX(dbl) :: c0( :, :, :, : )
     REAL(dbl) :: bec( :, : ), rhor( :, : ), rhos( :, : ), lambda( :, : ), lambdap( :, : )
     REAL(dbl) :: tau0( :, : ), h( 3, 3 )
     COMPLEX(dbl) :: eigrb( :, : ), rhog( :, : )
     COMPLEX(dbl) :: eigr( :, : )

     INTEGER :: is, istart, nss, i, j
     LOGICAL, SAVE :: trhor_save

     !se 0 stampa densita' di carica
     if( printwfc == 0 ) then
       call calbec(1,nsp,eigr,c0,bec)
       if(.not.tens) then
         call rhoofr(nfi,c0,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
       else
         write(6,*) 'Print wfc: ', printwfc
         !     calculation of the rotated quantities
         call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)
         !     calculation of rho corresponding to the rotated wavefunctions
         call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
       endif
       call write_rho_xsf(tau0,h,rhor)
#ifdef PARA
       call write_rho(47,nspin,rhor)
#endif            
     else if(printwfc <= n ) then
       !stampa funzione d'onda_n **2
       !usa lo stesso macchinario della occupazione variabile
       !PER ADESSO SOLO NSPIN=1 
       trhor_save=trhor
       write(6,*) 'Plotting band :', printwfc
       do  is=1,nspin
         istart=iupdwn(is)
         nss=nupdwn(is)
         call ddiag(nss,nss,lambda,dval(1),zaux(1,1,is),1)
         do i=1,nss
           e0(i+istart-1)=dval(i)
         enddo
       enddo
       do  is=1,nspin
         nss=nupdwn(is)
         istart=iupdwn(is)
         do i=1,nss
           do j=1,nss
             zx(j,i,is)=zaux(i,j,is)!ATTENZIONE ALLO SPIN
           end do
         enddo
       enddo

       call rotate(zx,c0(:,:,1,1),bec,c0diag,becdiag)
       do i=1,n
         !                  if(i.ne.printwfc) call zero(2*ngw, c0diag(1,i))!ATTENZIONE modifiche temporanee
         if(i.ne.92 .and. i.ne.93 .and. i.ne.94) c0diag(:,i) = 0.0d0
       enddo
       call calbec(1,nsp,eigr,c0diag,becdiag)
       call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
       call write_rho_xsf(tau0,h,rhor)
       trhor=trhor_save
     endif

     RETURN
   END SUBROUTINE cp_print_rho


!=----------------------------------------------------------------------------=!
   END MODULE print_out_module
!=----------------------------------------------------------------------------=!
