!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   SUBROUTINE printout_new_x   &
     ( nfi, tfirst, tfilei, tprint, tps, h, stress, tau0, vels, &
       fion, ekinc, temphc, tempp, temps, etot, enthal, econs, econt, &
       vnhh, xnhh0, vnhp, xnhp0, atot, ekin, epot )
!=----------------------------------------------------------------------------=!

      !
      USE kinds,             ONLY : DP
      USE control_flags,     ONLY : iprint
      USE energies,          ONLY : print_energies, dft_energy_type
      USE printout_base,     ONLY : printout_base_open, printout_base_close, &
                                    printout_pos, printout_cell, printout_stress
      USE constants,         ONLY : au_gpa, amu_si, bohr_radius_cm, &
                                    amu_au, BOHR_RADIUS_ANGS, pi
      USE ions_base,         ONLY : na, nsp, nat, ind_bck, atm, ityp, pmass, &
                                    cdm_displacement, ions_displacement
      USE cell_base,         ONLY : s_to_r, get_volume
      USE efield_module,     ONLY : tefield, pberryel, pberryion, &
                                    tefield2, pberryel2, pberryion2
      USE cg_module,         ONLY : tcg, itercg
      USE sic_module,        ONLY : self_interaction, sic_alpha, sic_epsilon
      USE electrons_module,  ONLY : print_eigenvalues
      USE pres_ai_mod,      ONLY : P_ext, Surf_t, volclu, surfclu, abivol, &
                                   abisur, pvar, n_ele

      USE xml_io_base,       ONLY : save_print_counter
      USE cp_main_variables, ONLY : nprint_nfi
      USE io_files,          ONLY : outdir
      USE control_flags,     ONLY : ndw, tdipole
      USE polarization,      ONLY : print_dipole
      USE io_global,         ONLY : ionode, ionode_id, stdout
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nfi
      LOGICAL, INTENT(IN) :: tfirst, tfilei, tprint
      REAL(DP), INTENT(IN) :: tps
      REAL(DP), INTENT(IN) :: h( 3, 3 )
      REAL(DP), INTENT(IN) :: stress( 3, 3 )
      REAL(DP), INTENT(IN) :: tau0( :, : )  ! real positions
      REAL(DP), INTENT(IN) :: vels( :, : )  ! scaled velocities
      REAL(DP), INTENT(IN) :: fion( :, : )  ! real forces
      REAL(DP), INTENT(IN) :: ekinc, temphc, tempp, etot, enthal, econs, econt
      REAL(DP), INTENT(IN) :: temps( : ) ! partial temperature for different ionic species
      REAL(DP), INTENT(IN) :: vnhh( 3, 3 ), xnhh0( 3, 3 ), vnhp( 1 ), xnhp0( 1 )
      REAL(DP), INTENT(IN) :: atot! enthalpy of system for c.g. case
      REAL(DP), INTENT(IN) :: ekin
      REAL(DP), INTENT(IN) :: epot ! ( epseu + eht + exc )
      !
      REAL(DP) :: stress_gpa( 3, 3 )
      REAL(DP) :: dis( nsp )
      REAL(DP) :: out_press, volume
      REAL(DP) :: totalmass
      INTEGER  :: isa, is, ia, kilobytes
      REAL(DP),         ALLOCATABLE :: tauw( :, : )
      CHARACTER(LEN=3), ALLOCATABLE :: labelw( : )
      LOGICAL  :: tsic, tfile
      LOGICAL, PARAMETER :: nice_output_files=.false.
      !
      ALLOCATE( labelw( nat ) )
      !
      ! avoid double printing to files by refering to nprint_nfi
      !
      tfile = tfilei .and. ( nfi .gt. nprint_nfi )
      !
      CALL memstat( kilobytes )
      !
      IF( ionode .AND. tfile .AND. tprint ) THEN
         CALL printout_base_open()
      END IF
      !
      IF( tprint ) THEN
         IF ( tfile ) THEN
            ! we're writing files, let's save nfi
            CALL save_print_counter( nfi, outdir, ndw )
         ELSE IF ( tfilei ) then
            ! not there yet, save the old nprint_nfi
            CALL save_print_counter( nprint_nfi, outdir, ndw )
         END IF
      END IF
      !
      volume = get_volume( h )
      !
      stress_gpa = stress * au_gpa
      !
      out_press = ( stress_gpa(1,1) + stress_gpa(2,2) + stress_gpa(3,3) ) / 3.0d0
      !
      IF( nfi > 0 ) THEN
         CALL update_accomulators &
              ( ekinc, ekin, epot, etot, tempp, enthal, econs, out_press, volume )
      END IF

      IF( ionode ) THEN
         !
         IF( tprint ) THEN
            !
            tsic = ( self_interaction /= 0 )
            !
            CALL print_energies( tsic, sic_alpha = sic_alpha, sic_epsilon = sic_epsilon )
            !
            CALL print_eigenvalues( 31, tfile, nfi, tps )
            !
            WRITE( stdout, * )
            !
            IF( kilobytes > 0 ) &
               WRITE( stdout, fmt="(3X,'Allocated memory (kb) = ', I9 )" ) kilobytes
            !
            WRITE( stdout, * )
            !
            IF( tdipole ) CALL print_dipole( 32, tfile, nfi, tps )
            !
            CALL printout_cell( stdout, h )
            !
            IF( tfile ) CALL printout_cell( 36, h, nfi, tps )
            !
            !  System density:
            !
            totalmass = 0.0d0
            DO is = 1, nsp
              totalmass = totalmass + pmass(is) * na(is) / amu_au
            END DO
            totalmass = totalmass / volume * 11.2061d0 ! AMU_SI * 1000.0 / BOHR_RADIUS_CM**3 
            WRITE( stdout, fmt='(/,3X,"System Density [g/cm^3] : ",F10.4,/)' ) totalmass

            CALL cdm_displacement( dis(1), tau0 )
            !
            WRITE( stdout,1000) dis(1)
            !
            CALL ions_displacement( dis, tau0 )
            !
            CALL printout_stress( stdout, stress_gpa )
            !
            IF( tfile ) CALL printout_stress( 38, stress_gpa, nfi, tps )
            !
            ! ... write out a standard XYZ file in angstroms
            !
            labelw( ind_bck(1:nat) ) = atm( ityp(1:nat) )
            !
            CALL printout_pos( stdout, tau0, nat, what = 'pos', &
                               label = labelw, sort = ind_bck )
            !
            IF( tfile ) then
               if (.not.nice_output_files) then
                  CALL printout_pos( 35, tau0, nat, nfi = nfi, tps = tps )
               else
                  CALL printout_pos( 35, tau0, nat, what = 'xyz', &
                               nfi = nfi, tps = tps, label = labelw, &
                               fact= BOHR_RADIUS_ANGS ,sort = ind_bck )
               endif
            END IF
            !
            ALLOCATE( tauw( 3, nat ) )
            !
            isa = 0
            !
            DO is = 1, nsp
               !
               DO ia = 1, na(is)
                  !
                  isa = isa + 1
                  !
                  CALL s_to_r( vels(:,isa), tauw(:,isa), h )
                  !
               END DO
               !
            END DO
            !
            WRITE( stdout, * )
            !
            CALL printout_pos( stdout, tauw, nat, &
                               what = 'vel', label = labelw, sort = ind_bck )
            !
            IF( tfile ) then
               if (.not.nice_output_files) then
                  CALL printout_pos( 34, tauw, nat, nfi = nfi, tps = tps )
               else
                  CALL printout_pos( 34, tauw, nat, nfi = nfi, tps = tps, &
                               what = 'vel', label = labelw, sort = ind_bck )
               endif
            END IF
            !
            WRITE( stdout, * )
            !
            CALL printout_pos( stdout, fion, nat, &
                               what = 'for', label = labelw, sort = ind_bck )
            !
            IF( tfile ) then
               if (.not.nice_output_files) then
                  CALL printout_pos( 37, fion, nat, nfi = nfi, tps = tps )
               else
                  CALL printout_pos( 37, fion, nat, nfi = nfi, tps = tps, &
                       what = 'for', label = labelw, sort = ind_bck )
               endif
            END IF
            !
            DEALLOCATE( tauw )
            !
            ! ...       Write partial temperature and MSD for each atomic specie tu stdout
            !
            WRITE( stdout, * ) 
            WRITE( stdout, 1944 )
            !
            DO is = 1, nsp
               WRITE( stdout, 1945 ) is, temps(is), dis(is)
            END DO
            !
            IF( tfile ) WRITE( 33, 2948 ) nfi, ekinc, temphc, tempp, etot, enthal, &
                                          econs, econt, volume, out_press, tps
            IF( tfile ) WRITE( 39, 2949 ) nfi, vnhh(3,3), xnhh0(3,3), vnhp(1), &
                                          xnhp0(1), tps
            !
         END IF
         !
       END IF
      !

       IF( ionode .AND. tfile .AND. tprint ) THEN
         !
         ! ... Close and flush unit 30, ... 40
         !
         CALL printout_base_close()
         !
      END IF
      !
        IF( ( MOD( nfi, iprint ) == 0 ) .OR. tfirst )  THEN
           !
           WRITE( stdout, * )
           WRITE( stdout, 1947 )
           if (abivol.and.pvar) write(stdout,*) 'P = ', P_ext*au_gpa
           !
        END IF
      ! 
      if (abivol) then
         write(stdout,*) nfi, 'ab-initio volume = ', volclu, ' a.u.^3'
         write(stdout,*) nfi, 'PV = ', P_ext*volclu, ' ha'
      end if
      if (abisur) then
         write(stdout,*) nfi, 'ab-initio surface = ', surfclu, ' a.u.^2'
         if (abivol) write(stdout,*) nfi, 'spherical surface = ', &
                 4.d0*pi*(0.75d0*volclu/pi)**(2.d0/3.d0), ' a.u.^2'
         write(stdout,*) nfi, 't*S = ', Surf_t*surfclu, ' ha'
      end if
      if (abivol.or.abisur) write(stdout,*) nfi, &
         ' # of electrons within the isosurface = ', n_ele
      IF( .not. tcg ) THEN
         !
         WRITE( stdout, 1948 ) nfi, ekinc, temphc, tempp, etot, enthal, econs, &
                            econt, vnhh(3,3), xnhh0(3,3), vnhp(1),  xnhp0(1)
      ELSE
         IF ( MOD( nfi, iprint ) == 0 .OR. tfirst ) THEN
            !
            WRITE( stdout, * )
            WRITE( stdout, 255 ) 'nfi','tempp','E','-T.S-mu.nbsp','+K_p','#Iter'
            !
         END IF
         !
         WRITE( stdout, 256 ) nfi, INT( tempp ), etot, atot, econs, itercg
         !
      END IF

      IF( tefield) THEN
         IF(ionode) write(stdout,'( A14,F12.6,2X,A14,F12.6)') 'Elct. dipole 1',-pberryel,'Ionic dipole 1',-pberryion
      ENDIF
      IF( tefield2) THEN
         IF(ionode) write(stdout,'( A14,F12.6,2X,A14,F12.6)') 'Elct. dipole 2',-pberryel2,'Ionic dipole 2',-pberryion2
      ENDIF

      !
      !
      DEALLOCATE( labelw )
      !
255   FORMAT( '     ',A5,A8,3(1X,A12),A6 )
256   FORMAT( 'Step ',I5,1X,I7,1X,F12.5,1X,F12.5,1X,F12.5,1X,I5 )
1000  FORMAT(/,3X,'Center of mass square displacement (a.u.): ',F10.6,/)
1944  FORMAT(//'   Partial temperatures (for each ionic specie) ', &
             /,'   Species  Temp (K)   Mean Square Displacement (a.u.)')
1945  FORMAT(3X,I6,1X,F10.2,1X,F10.4)
1947  FORMAT( 2X,'nfi',4X,'ekinc',2X,'temph',2X,'tempp',8X,'etot',6X,'enthal', &
           & 7X,'econs',7X,'econt',4X,'vnhh',3X,'xnhh0',4X,'vnhp',3X,'xnhp0' )
1948  FORMAT( I5,1X,F8.5,1X,F6.1,1X,F6.1,4(1X,F11.5),4(1X,F7.4) )
2948  FORMAT( I6,1X,F8.5,1X,F6.1,1X,F6.1,4(1X,F11.5),F10.2, F8.2, F8.5 )
2949  FORMAT( I6,1X,4(1X,F7.4), F8.5 )
      !
      RETURN
   END SUBROUTINE printout_new_x
   !  
   !

!=----------------------------------------------------------------------------=!
   SUBROUTINE printout_x(nfi, atoms, ekinc, ekcell, tprint, ht, edft) 
!=----------------------------------------------------------------------------=!

      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: tdipole, tnosee, tnosep, tnoseh, iprsta, iprint, &
                                    toptical, tconjgrad
      use constants,          only: k_boltzmann_au, au_gpa, amu_si, bohr_radius_cm
      use energies,           only: print_energies, dft_energy_type
      use mp_global,          only: me_image, intra_image_comm
      use electrons_module,   only: print_eigenvalues
      use time_step,          ONLY: tps
      USE electrons_nose,     ONLY: electrons_nose_nrg, xnhe0, vnhe, qne, ekincw
      USE sic_module,         ONLY: ind_localisation, pos_localisation, nat_localisation, &
                                    self_interaction, sic_rloc
      !!USE ions_base,        ONLY: ions_displacement, cdm_displacement
      USE ions_base,          ONLY: ions_temp, cdmi, taui, nsp
      USE ions_nose,          ONLY: ndega, ions_nose_nrg, xnhp0, vnhp, qnp, gkbt, &
                                    kbt, nhpcl, nhpdim, atm2nhp, ekin2nhp, gkbt2nhp
      USE cell_nose,          ONLY: cell_nose_nrg, qnh, temph, xnhh0, vnhh
      USE cell_base,          ONLY: iforceh, boxdimensions, s_to_r, press
      USE printout_base,      ONLY: printout_base_open, printout_base_close, &
                                    printout_pos, printout_cell, printout_stress
      USE environment,        ONLY: start_cclock_val
      USE atoms_type_module,  ONLY: atoms_type
      USE cp_interfaces,      ONLY: printout_new

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nfi
      TYPE (atoms_type)   :: atoms
      LOGICAL             :: tprint
      type (boxdimensions), intent(in) :: ht
      TYPE (dft_energy_type) :: edft
      REAL(DP) :: ekinc, ekcell
!
! ...
      INTEGER   :: is, ia, k, i, j, ik, isa, iunit, nfill, nempt
      REAL(DP) :: tau(3), vel(3), stress_tensor(3,3), temps( atoms%nsp )
      REAL(DP) :: tempp, econs, ettt, out_press, ekinpr, enosee
      REAL(DP) :: enthal, totalmass, enoseh, temphc, enosep
      REAL(DP) :: h(3,3)
      REAL(DP) :: epot
      !!REAL(DP) :: dis( nsp )
      LOGICAL   :: tfile, topen, tsic, tfirst
      CHARACTER(LEN=3), ALLOCATABLE :: labelw( : )
      REAL(DP), ALLOCATABLE :: tauw( :, : )
      INTEGER   :: old_nfi = -1

      ! ...   Subroutine Body

      tfile = ( MOD( nfi, iprint ) == 0 )   !  print quantity to trajectory files
      tsic  = ( self_interaction /= 0 )
      
      ! ...   Calculate Ions temperature tempp (in Kelvin )

      CALL ions_temp &
           ( tempp, temps, ekinpr, atoms%vels, atoms%na, atoms%nsp, ht%hmat, &
             atoms%m, ndega, nhpdim, atm2nhp, ekin2nhp )

      ! ...   Stress tensor (in GPa) and pressure (in GPa)

      stress_tensor = MATMUL( ht%pail(:,:), ht%a(:,:) ) * au_gpa / ht%deth
      !
      out_press = ( stress_tensor(1,1) + stress_tensor(2,2) + stress_tensor(3,3) ) / 3.0d0

      ! ...   Enthalpy (in Hartree)

      enthal = edft%etot + press * ht%deth

      IF( tnoseh ) THEN
        enoseh = cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
      ELSE
        enoseh = 0.0d0
      END IF

      if( COUNT( iforceh == 1 ) > 0 ) then
         temphc = 2.0d0 / k_boltzmann_au * ekcell / DBLE( COUNT( iforceh == 1 ) )
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

      epot   = edft%eht + edft%exc + edft%epseu

      ! ...   Print physical variables to fortran units

      tfirst = tprint .OR. tconjgrad
      stress_tensor = stress_tensor / au_gpa

      ALLOCATE( tauw( 3, atoms%nat ) )
      DO is = 1, atoms%nsp
        DO ia = atoms%isa(is), atoms%isa(is) + atoms%na(is) - 1
          CALL s_to_r( atoms%taus(:,ia), tauw(:,ia), ht )
        END DO
      END DO

      CALL  printout_new &
     ( nfi, tfirst, tfile, tprint, tps, ht%hmat, stress_tensor, tauw, atoms%vels, &
       atoms%for, ekinc, temphc, tempp, temps, edft%etot, enthal, econs, ettt, &
       vnhh, xnhh0, vnhp, xnhp0, 0.0d0, edft%ekin, epot )

      DEALLOCATE( tauw )

      old_nfi = nfi

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
 1944 FORMAT(//'   Partial temperatures (for each ionic specie) ', &
             /,'   Species  Temp (K)   MSD (AU)')
 1945 FORMAT(3X,I6,1X,F10.2,1X,F10.4)


    RETURN
  END SUBROUTINE printout_x



!=----------------------------------------------------------------------------=!
  SUBROUTINE print_legend()
!=----------------------------------------------------------------------------=!
    !
    USE io_global, ONLY : ionode, stdout
    !
    IMPLICIT NONE
    !
    IF ( .NOT. ionode ) RETURN
    !
    WRITE( stdout, *) 
    WRITE( stdout, *) '  Short Legend and Physical Units in the Output'
    WRITE( stdout, *) '  ---------------------------------------------'
    WRITE( stdout, *) '  NFI    [int]          - step index'
    WRITE( stdout, *) '  EKINC  [HARTREE A.U.] - kinetic energy of the fictitious electronic dynamics'
    WRITE( stdout, *) '  TEMPH  [K]            - Temperature of the fictitious cell dynamics'
    WRITE( stdout, *) '  TEMP   [K]            - Ionic temperature'
    WRITE( stdout, *) '  ETOT   [HARTREE A.U.] - Scf total energy (Kohn-Sham hamiltonian)'
    WRITE( stdout, *) '  ENTHAL [HARTREE A.U.] - Enthalpy ( ETOT + P * V )'
    WRITE( stdout, *) '  ECONS  [HARTREE A.U.] - Enthalpy + kinetic energy of ions and cell'
    WRITE( stdout, *) '  ECONT  [HARTREE A.U.] - Constant of motion for the CP lagrangian'
    WRITE( stdout, *) 
    !
    RETURN
    !
  END SUBROUTINE print_legend




!=----------------------------------------------------------------------------=!
    SUBROUTINE print_sfac_x( rhoe, sfac )
!=----------------------------------------------------------------------------=!

      USE kinds,              ONLY : DP
      USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
      USE mp,                 ONLY: mp_max, mp_get, mp_put
      USE reciprocal_vectors, ONLY: ig_l2g, gx, g
      USE gvecp,              ONLY: ngm
      USE fft_base,           ONLY : dfftp
      USE cp_interfaces,      ONLY: fwfft
      USE io_global,          ONLY : ionode, ionode_id, stdout
      USE io_files,           ONLY: sfacunit, sfac_file, opt_unit

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: rhoe(:,:)
      COMPLEX(DP), INTENT(IN) ::  sfac(:,:)

      INTEGER :: nspin, ispin, ip, nsp, ngx_l, ng, is, ig
      COMPLEX(DP), ALLOCATABLE :: rhoeg(:,:)
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: rhoeg_rcv(:,:)
      REAL   (DP), ALLOCATABLE :: hg_rcv(:)
      REAL   (DP), ALLOCATABLE :: gx_rcv(:,:)
      INTEGER     , ALLOCATABLE :: ig_rcv(:)
      COMPLEX(DP), ALLOCATABLE :: sfac_rcv(:,:)

        nspin = SIZE(rhoe,2)
        nsp   = SIZE(sfac,2)
        ngx_l = ngm
        CALL mp_max(ngx_l, intra_image_comm)
        ALLOCATE(rhoeg(ngm,nspin))
        ALLOCATE(hg_rcv(ngx_l))
        ALLOCATE(gx_rcv(3,ngx_l))
        ALLOCATE(ig_rcv(ngx_l))
        ALLOCATE(rhoeg_rcv(ngx_l,nspin))
        ALLOCATE(sfac_rcv(ngx_l,nsp))

        ! ...   FFT: rho(r) --> rho(g)
        !
        ALLOCATE( psi( SIZE( rhoe, 1 ) ) )
        !
        DO ispin = 1, nspin
          psi = rhoe(:,ispin)
          CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
          CALL psi2rho( 'Dense', psi, dfftp%nnr, rhoeg(:,ispin), ngm )
        END DO

        DEALLOCATE( psi )

        IF( ionode ) THEN
          OPEN(sfacunit, FILE=TRIM(sfac_file), STATUS='UNKNOWN')
        END IF

        DO ip = 1, nproc_image
          CALL mp_get(ng, ngm, me_image, ionode_id, ip-1, ip, intra_image_comm )
          CALL mp_get(hg_rcv(:), g(:), me_image, ionode_id, ip-1, ip, intra_image_comm )
          CALL mp_get(gx_rcv(:,:), gx(:,:), me_image, ionode_id, ip-1, ip, intra_image_comm )
          CALL mp_get(ig_rcv(:), ig_l2g(:), me_image, ionode_id, ip-1, ip, intra_image_comm )
          DO ispin = 1, nspin
            CALL mp_get( rhoeg_rcv(:,ispin), rhoeg(:,ispin), me_image, ionode_id, ip-1, ip, intra_image_comm )
          END DO
          DO is = 1, nsp
            CALL mp_get( sfac_rcv(:,is), sfac(:,is), me_image, ionode_id, ip-1, ip, intra_image_comm )
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
    END SUBROUTINE print_sfac_x




!=----------------------------------------------------------------------------=!
   SUBROUTINE printacc( )
!=----------------------------------------------------------------------------=!

      USE kinds,               ONLY : DP
      USE cp_main_variables,   ONLY : acc, acc_this_run, nfi, nfi_run
      USE io_global,           ONLY : ionode, stdout

      IMPLICIT NONE
      !
      REAL(DP) :: avgs(9)
      REAL(DP) :: avgs_run(9)
 
      avgs     = 0.0d0
      avgs_run = 0.0d0
      !
      IF ( nfi > 0 ) THEN
         avgs  = acc( 1:9 ) / DBLE( nfi )
      END IF
      !
      IF ( nfi_run > 0 ) THEN
         avgs_run = acc_this_run(1:9) / DBLE( nfi_run )
      END IF

      IF( ionode ) THEN
        WRITE( stdout,1949)
        WRITE( stdout,1951) avgs(1), avgs_run(1)
        WRITE( stdout,1952) avgs(2), avgs_run(2)
        WRITE( stdout,1953) avgs(3), avgs_run(3)
        WRITE( stdout,1954) avgs(4), avgs_run(4)
        WRITE( stdout,1955) avgs(5), avgs_run(5)
        WRITE( stdout,1956) avgs(6), avgs_run(6)
        WRITE( stdout,1957) avgs(7), avgs_run(7)
        WRITE( stdout,1958) avgs(8), avgs_run(8)
        WRITE( stdout,1959) avgs(9), avgs_run(9)
        WRITE( stdout,1990)
 1949   FORMAT(//,3X,'Averaged Physical Quantities',/ &
              ,3X,'                  ',' accomulated','      this run')
 1951   FORMAT(3X,'ekinc         : ',F14.5,F14.5,' (AU)')
 1952   FORMAT(3X,'ekin          : ',F14.5,F14.5,' (AU)')
 1953   FORMAT(3X,'epot          : ',F14.5,F14.5,' (AU)')
 1954   FORMAT(3X,'total energy  : ',F14.5,F14.5,' (AU)')
 1955   FORMAT(3X,'temperature   : ',F14.5,F14.5,' (K )')
 1956   FORMAT(3X,'enthalpy      : ',F14.5,F14.5,' (AU)')
 1957   FORMAT(3X,'econs         : ',F14.5,F14.5,' (AU)')
 1958   FORMAT(3X,'pressure      : ',F14.5,F14.5,' (Gpa)')
 1959   FORMAT(3X,'volume        : ',F14.5,F14.5,' (AU)')
 1990   FORMAT(/)
      END IF

      RETURN
    END SUBROUTINE printacc



!=----------------------------------------------------------------------------=!
    SUBROUTINE open_and_append_x( iunit, file_name )
!=----------------------------------------------------------------------------=!
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
    END SUBROUTINE open_and_append_x




!=----------------------------------------------------------------------------=!
   SUBROUTINE print_projwfc_x ( c0, lambda, eigr, vkb )
!=----------------------------------------------------------------------------=!
      USE kinds,            ONLY: DP
      USE electrons_base,   ONLY: nspin, nbnd, nbsp, iupdwn, nupdwn
      USE electrons_module, ONLY: ei, ei_emp, n_emp, iupdwn_emp, nupdwn_emp
      USE cp_interfaces,    ONLY: set_evtot, set_eitot
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(IN)  :: c0(:,:), eigr(:,:), vkb(:,:)
      REAL(DP),    INTENT(IN)  :: lambda(:,:,:)
      !
      INTEGER  :: nupdwn_tot( 2 ), iupdwn_tot( 2 )
      COMPLEX(DP), ALLOCATABLE :: ctmp(:,:)
      REAL(DP),    ALLOCATABLE :: eitot(:,:)
      !
      nupdwn_tot = nupdwn + nupdwn_emp
      iupdwn_tot(1) = iupdwn(1)
      iupdwn_tot(2) = nupdwn(1) + 1
      !
      ALLOCATE( eitot( nupdwn_tot(1), nspin ) )
      !
      CALL set_eitot( eitot )
      !
      ALLOCATE( ctmp( SIZE( c0, 1 ), nupdwn_tot(1) * nspin ) )
      ! 
      CALL set_evtot( c0, ctmp, lambda, iupdwn_tot, nupdwn_tot )
      !
      CALL projwfc( ctmp, SIZE(ctmp,2), eigr, vkb, nupdwn_tot(1), eitot(1,1)  )
      !
      DEALLOCATE( eitot )
      DEALLOCATE( ctmp )

      RETURN
   END SUBROUTINE


!=----------------------------------------------------------------------------=!
   SUBROUTINE update_accomulators &
      ( ekinc, ekin, epot, etot, tempp, enthal, econs, press, volume )
!=----------------------------------------------------------------------------=!

      USE kinds,               ONLY : DP
      USE cp_main_variables,   ONLY : acc, acc_this_run, nfi_run

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ekinc, ekin, epot, etot, tempp
      REAL(DP), INTENT(IN) :: enthal, econs, press, volume

      nfi_run = nfi_run + 1

      ! ...   sum up values to be averaged

      acc(1) = acc(1) + ekinc
      acc(2) = acc(2) + ekin
      acc(3) = acc(3) + epot
      acc(4) = acc(4) + etot
      acc(5) = acc(5) + tempp
      acc(6) = acc(6) + enthal
      acc(7) = acc(7) + econs
      acc(8) = acc(8) + press  ! pressure in GPa
      acc(9) = acc(9) + volume

      ! ...   sum up values to be averaged

      acc_this_run(1) = acc_this_run(1) + ekinc
      acc_this_run(2) = acc_this_run(2) + ekin
      acc_this_run(3) = acc_this_run(3) + epot
      acc_this_run(4) = acc_this_run(4) + etot
      acc_this_run(5) = acc_this_run(5) + tempp
      acc_this_run(6) = acc_this_run(6) + enthal
      acc_this_run(7) = acc_this_run(7) + econs
      acc_this_run(8) = acc_this_run(8) + press  ! pressure in GPa
      acc_this_run(9) = acc_this_run(9) + volume

      RETURN
   END SUBROUTINE
