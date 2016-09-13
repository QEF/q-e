!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
! modified by Gerardo Ballabio and Carlo Cavazzoni
! on Thu Jul 22 18:57:26 CEST 2004
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! --------------------------------------------------------------------
! this routine writes the crystal structure in XSF, GRD and PDB format
! from a FPMD output files
! --------------------------------------------------------------------
PROGRAM fpmd_postproc

  USE kinds,      ONLY : DP
  USE constants,  ONLY : bohr => BOHR_RADIUS_ANGS
  USE io_files,   ONLY : prefix, iunpun, xmlpun, tmp_dir
  USE mp_global,  ONLY : mp_startup, mp_global_end

  USE iotk_module
  USE xml_io_base

  IMPLICIT NONE

  INTEGER, PARAMETER :: maxsp = 20

  INTEGER                    :: natoms, nsp, na(maxsp), atomic_number(maxsp)
  INTEGER                    :: ounit, cunit, punit, funit, dunit, bunit, ksunit
  INTEGER                    :: nr1, nr2, nr3, ns1, ns2, ns3
  INTEGER                    :: np1, np2, np3, np, ispin
  INTEGER, ALLOCATABLE       :: ityp(:)
  REAL(DP)              :: at(3, 3), atinv(3, 3), ht0(3, 3), h0(3, 3)
  REAL(DP)              :: rhof, rhomax, rhomin, rhoc(6)
  REAL(DP), ALLOCATABLE :: rho_in(:,:,:), rho_out(:,:,:)
  REAL(DP), ALLOCATABLE :: tau_in(:,:), tau_out(:,:)
  REAL(DP), ALLOCATABLE :: sigma(:,:), force(:,:)
  REAL(DP), ALLOCATABLE :: stau0(:,:), svel0(:,:), force0(:,:)

  CHARACTER(len=256) :: filepp, fileout, output, outdir
  CHARACTER(len=256) :: filecel, filepos, filefor, filepdb
  CHARACTER(len=256) :: print_state
  CHARACTER(len=3)   :: atm( maxsp ), lab
  CHARACTER(len=4)   :: charge_density
  LOGICAL            :: lcharge, lforces, ldynamics, lpdb, lrotation
  LOGICAL            :: lbinary, found
  INTEGER            :: nframes
  INTEGER            :: ios, nat, ndr
  INTEGER            :: nproc, mpime, world, root

  REAL(DP) :: x, y, z, fx, fy, fz
  INTEGER  :: i, j, k, n, ix, iy, iz, ierr

  REAL(DP) :: euler(6)

  NAMELIST /inputpp/ prefix, fileout, output, outdir, &
                     lcharge, lforces, ldynamics, lpdb, lrotation, &
                     ns1, ns2, ns3, np1, np2, np3, print_state, &
                     atomic_number, nframes, ndr, charge_density, &
                     lbinary

  ! default values

  dunit = 14
#if defined(__INTEL_COMPILER)
  CALL remove_stack_limit ( )
#endif
  !  initialize mpi
  CALL mp_startup  ( )
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix    = 'cp'
  fileout   = 'out'
  output    = 'xsf'  ! 'grd'
  lcharge   = .false.
  lforces   = .false.
  ldynamics = .false.
  lpdb      = .false.
  lrotation = .false.
  ns1   = 0
  ns2   = 0
  ns3   = 0
  np1   = 1
  np2   = 1
  np3   = 1                  ! 
  nframes = 1                ! number of MD step to be read to buind the trajectory
  ndr   = 51                 ! restart file number
  atomic_number = 1          ! atomic number of the species in the restart file
  charge_density = 'full'    ! specify the component to plot: 'full', 'spin'
  print_state = ' '          ! specify the Kohn-Sham state to plot: 'KS_1'
  lbinary = .TRUE.

  

  call input_from_file()

  ! read namelist
  READ( 5, inputpp, iostat=ios)

  ! set file names
  !
  tmp_dir = TRIM(outdir)
  filecel = TRIM(tmp_dir) // TRIM(prefix) // '.cel'
  filepos = TRIM(tmp_dir) // TRIM(prefix) // '.pos'
  filefor = TRIM(tmp_dir) // TRIM(prefix) // '.for'
  !
  filepdb = TRIM(fileout) // '.pdb'
  !
  ! append extension
  !
  IF (output == 'xsf') THEN
     IF (ldynamics) THEN
        fileout = TRIM(fileout) // '.axsf'
     ELSE
        fileout = TRIM(fileout) // '.xsf'
     END IF
  ELSE IF (output == 'grd') THEN
     fileout = TRIM(fileout) // '.grd'
  ELSE IF (output == 'xyz') THEN
     fileout = TRIM(fileout) // '.xyz'
  END IF

  np = np1 * np2 * np3
  IF (np1 < 1 .OR. np2 < 1 .OR. np3 < 1) THEN
     WRITE(*,*) 'Error: zero or negative replicas not allowed'
     STOP
  END IF


  ! check for wrong input
  IF (ldynamics .AND. nframes < 2) THEN
     WRITE(*,*) 'Error: dynamics requested, but only one frame'
     STOP
  END IF
  IF (.NOT. ldynamics) nframes = 1

  IF (ldynamics .AND. lcharge) THEN
     WRITE(*,*) 'Error: dynamics with charge density not supported'
     STOP
  END IF

  IF (ldynamics .AND. ( print_state /= ' ' ) ) THEN
     WRITE(*,*) 'Error: dynamics with print_state not supported'
     STOP
  END IF

  !
  !  Now read the XML data file
  !

  filepp = restart_dir( tmp_dir, ndr )
  !
  filepp = TRIM( filepp ) // '/' // TRIM(xmlpun)
  !
  CALL iotk_open_read( dunit, file = TRIM( filepp ), BINARY = .FALSE., &
     ROOT = attr, IERR = ierr )
  IF( ierr /= 0 ) CALL errore( ' cppp ', ' Cannot open file '//TRIM(filepp), 1 )

     CALL iotk_scan_begin( dunit, "IONS", FOUND = found )

     IF( .NOT. found ) THEN
        CALL errore( ' cppp ', ' IONS not found in data-file.xml ', 1 )
     END IF

     CALL iotk_scan_dat( dunit, "NUMBER_OF_ATOMS", nat )
     CALL iotk_scan_dat( dunit, "NUMBER_OF_SPECIES", nsp )

     ALLOCATE( ityp( nat * np ) ) ! atomic species

     DO i = 1, nsp
        !
        CALL iotk_scan_begin( dunit, "SPECIE" // TRIM( iotk_index( i ) ), FOUND = found )
        !
        IF( .NOT. found ) THEN
           CALL errore( ' cppp ', "SPECIE" // TRIM( iotk_index( i ) ) // ' not found in data-file.xml ', 1 )
        END IF

        CALL iotk_scan_dat( dunit, "ATOM_TYPE", atm(i) )
        !
        ! CALL iotk_scan_dat( dunit, &
        !                     TRIM( atm(i) )//"_MASS", amass(i), ATTR = attr )
        !
        ! CALL iotk_scan_dat( dunit, &
        !                     "PSEUDO_FOR_" // TRIM( atm(i) ), psfile(i) )
        !
        CALL iotk_scan_end( dunit, "SPECIE" // TRIM( iotk_index( i ) ) )
        !
     END DO

      !
      ! CALL iotk_scan_empty( dunit, "UNITS_FOR_ATOMIC_POSITIONS", attr )
      ! CALL iotk_scan_attr( attr, "UNIT", pos_unit  )
      !
      DO i = 1, nat
         !
         CALL iotk_scan_empty( dunit, "ATOM" // TRIM( iotk_index( i ) ), attr )
         CALL iotk_scan_attr( attr, "SPECIES", lab )
         CALL iotk_scan_attr( attr, "INDEX",   ityp(i) )
         ! CALL iotk_scan_attr( attr, "tau",     tau(:,i) )
         ! CALL iotk_scan_attr( attr, "if_pos",  if_pos(:,i) )
         !
      END DO

     CALL iotk_scan_end( dunit, "IONS" )

     CALL iotk_scan_begin( dunit, "PLANE_WAVES" )
        CALL iotk_scan_empty( dunit, "FFT_GRID", attr )
        CALL iotk_scan_attr( attr, "nr1", nr1 )
        CALL iotk_scan_attr( attr, "nr2", nr2 )
        CALL iotk_scan_attr( attr, "nr3", nr3 )
     CALL iotk_scan_end( dunit, "PLANE_WAVES" )

     ALLOCATE( stau0( 3, nat ) )
     ALLOCATE( svel0( 3, nat ) )
     ALLOCATE( force0( 3, nat ) ) ! forces, atomic units

     CALL iotk_scan_begin( dunit, "TIMESTEPS", attr )
        CALL iotk_scan_begin( dunit, "STEP0" )
           CALL iotk_scan_begin( dunit, "IONS_POSITIONS" )
              CALL iotk_scan_dat( dunit, "stau", stau0 )
              CALL iotk_scan_dat( dunit, "svel", svel0 )
              CALL iotk_scan_dat( dunit, "force", force0 )
           CALL iotk_scan_end( dunit, "IONS_POSITIONS" )
           CALL iotk_scan_begin( dunit, "CELL_PARAMETERS" )
              CALL iotk_scan_dat( dunit, "ht", ht0 )
           CALL iotk_scan_end( dunit, "CELL_PARAMETERS" )
        CALL iotk_scan_end( dunit, "STEP0" )
     CALL iotk_scan_end( dunit, "TIMESTEPS" )
     !
     ispin = 1
     !
     !
  CALL iotk_close_read( dunit )

  !
  !   End of reading from data file
  !

  IF ( nsp > maxsp ) THEN
     WRITE(*,*) 'Error: too many atomic species'
     STOP
  END IF

  natoms = nat
  !
  ! Count atoms in each species
  !
  na = 0
  DO i = 1, nat
     na( ityp( i ) ) = na( ityp( i ) ) + 1                 ! total number of atoms
  END DO

  ! assign species (from input) to each atom
  !
  k = 0
  DO i = 1, nsp
     DO j = 1, na(i)
        k = k + 1
        ityp(k) = atomic_number(i)
     END DO
  END DO


  ! allocate arrays
  ALLOCATE(tau_in(3, nat))                  ! atomic positions, angstroms
  ALLOCATE(tau_out(3, nat * np))            ! replicated positions
  ALLOCATE(sigma(3, nat ) )                 ! scaled coordinates
  !
  IF (lforces) ALLOCATE( force( 3, nat * np ) )

  ! charge density
  IF ( lcharge .OR. print_state /= ' ' ) THEN
     IF (ns1 == 0) ns1 = nr1
     IF (ns2 == 0) ns2 = nr2
     IF (ns3 == 0) ns3 = nr3
     ALLOCATE( rho_in ( nr1, nr2, nr3 ) )     ! original charge density
     ALLOCATE( rho_out( ns1, ns2, ns3 ) )     ! rescaled charge density
  END IF

  ! open output file for trajectories or charge density 
  !
  ounit = 10
  OPEN(ounit, file=fileout, status='unknown')

  ! open Cell trajectory file
  !
  cunit = 11
  OPEN(cunit, file=filecel, status='old')

  ! open Positions trajectory file
  !
  punit = 12
  OPEN(punit, file=filepos, status='old')

  ! open Force trajectory file
  !
  funit = 13
  if (lforces) OPEN(funit, file=filefor, status='old')

  ! open PDB file
  !
  bunit = 15
  OPEN(bunit, file=filepdb, status='unknown')

  ! Unit for KS states
  !
  ksunit = 16

  ! XSF file header
  !
  IF ( output == 'xsf' ) THEN
     IF ( ldynamics ) WRITE(ounit,*) 'ANIMSTEPS', nframes
     WRITE( ounit, * ) 'CRYSTAL'
  END IF

  DO n = 1, nframes
     !
     IF ( ldynamics ) WRITE(*,'("frame",1X,I4)') n

     ! read data from files produced by fpmd
     !
     CALL read_fpmd( lforces, lcharge, lbinary, cunit, punit, funit, dunit, &
                     natoms, nr1, nr2, nr3, ispin, at, tau_in, force, &
                     rho_in, prefix, tmp_dir, ndr, charge_density )

     IF( nframes == 1 ) THEN
        !
        !  use values from the XML file
        !
        IF( lforces ) force( 1:3, 1:nat ) = force0( 1:3, 1:nat ) 
        !
        h0 = TRANSPOSE( ht0 )
        !
        ! from scaled to real coordinates
        !
        tau_in( :, : ) = MATMUL( h0( :, : ), stau0( :, : ) )
        !
        ! convert atomic units to Angstroms
        !
        at     = h0     * bohr
        tau_in = tau_in * bohr
        !
     END IF

     WRITE(*,'(2x,"Cell parameters (Angstroms):")')
     WRITE(*,'(3(2x,f10.6))') ((at(i, j), i=1,3), j=1,3)
     !
     WRITE(*,'(2x,"Atomic coordinates (Angstroms):")')
     WRITE(*,'(3(2x,f10.6))') ((tau_in(i, j), i=1,3), j=1,natoms)

     ! compute scaled coordinates
     !
     CALL inverse( at, atinv )
     sigma(:,:) = MATMUL(atinv(:,:), tau_in(:,:))

     ! compute cell dimensions and Euler angles
     CALL at_to_euler( at, euler )

     IF (lpdb) THEN
        ! apply periodic boundary conditions
        DO i = 1, natoms
           DO j = 1, 3
              sigma(j, i) = sigma(j, i) - FLOOR(sigma(j, i))
           END DO
        END DO
        ! recompute Cartesian coordinates
        tau_in(:,:) = MATMUL(at(:,:), sigma(:,:))
     END IF

     IF (lrotation) THEN
        ! compute rotated cell
        CALL euler_to_at( euler, at )
        ! rotate atomic positions as well
        tau_in(:,:) = MATMUL(at(:,:), sigma(:,:))
     END IF

     ! replicate atoms
     k = 0
     DO ix = 1, np1
        DO iy = 1, np2
           DO iz = 1, np3
              DO j = 1, natoms
                 k = k + 1
                 tau_out(:, k) = tau_in(:, j) + (ix-1) * at(:, 1) + &
                                 (iy-1) * at(:, 2) + (iz-1) * at(:, 3)
                 ityp(k) = ityp(j)
                 IF (lforces) force(:, k) = force(:, j)
              END DO
           END DO
        END DO
     END DO
     natoms = natoms * np

     ! compute supercell
     at(:, 1) = at(:, 1) * np1
     at(:, 2) = at(:, 2) * np2
     at(:, 3) = at(:, 3) * np3
     euler(1) = euler(1) * np1
     euler(2) = euler(2) * np2
     euler(3) = euler(3) * np3

     IF ( lcharge ) &
        CALL scale_charge( rho_in, rho_out, nr1, nr2, nr3, ns1, ns2, ns3, &
                           np1, np2, np3 )

     IF ( output == 'xsf' ) THEN
        ! write data as XSF format
        CALL write_xsf( ldynamics, lforces, lcharge, ounit, n, at, &
                        natoms, ityp, tau_out, force, rho_out, &
                        ns1, ns2, ns3 )
     ELSE IF( output == 'grd' ) THEN
        ! write data as GRD format
        CALL write_grd( ounit, at, rho_out, ns1, ns2, ns3 )
     ELSE IF( output == 'xyz' ) THEN
        ! write data as XYZ format
        CALL write_xyz( ldynamics, lforces, lcharge, ounit, n, at, &
                        natoms, ityp, tau_out, force, rho_out, &
                        ns1, ns2, ns3 )
     END IF

  END DO

  CLOSE(ounit)

  IF ( print_state /= ' ' ) THEN
     !
     CALL read_density( TRIM( print_state ) // '.xml', dunit, nr1, nr2, nr3, rho_in, lbinary )
     CALL scale_charge( rho_in, rho_out, nr1, nr2, nr3, ns1, ns2, ns3, np1, np2, np3 )
     !
     IF (output == 'xsf') THEN
        ! write data as XSF format
        OPEN( unit = ksunit, file = TRIM( print_state ) // '.xsf' )
        WRITE( ksunit, * ) 'CRYSTAL'  !  XSF files need this one line header
        CALL write_xsf( ldynamics, lforces, .true., ksunit, n, at, &
                        natoms, ityp, tau_out, force, rho_out, ns1, ns2, ns3 )
     ELSE IF( output == 'grd' ) THEN
        OPEN( unit = ksunit, file = TRIM( print_state ) // '.grd' )
        CALL write_grd( ksunit, at, rho_out, ns1, ns2, ns3 )
     END IF
     !
     CLOSE( ksunit )
     !
  END IF

  ! write atomic positions as PDB format
  CALL write_pdb( bunit, at, tau_out, natoms, ityp, euler, lrotation )

  ! free allocated resources
  CLOSE(punit)
  CLOSE(cunit)
  IF (lforces) CLOSE(funit)

  DEALLOCATE(tau_in)
  DEALLOCATE(tau_out)
  DEALLOCATE(ityp)
  IF( ALLOCATED( force  ) ) DEALLOCATE(force)
  IF( ALLOCATED( rho_in ) ) DEALLOCATE(rho_in)
  IF( ALLOCATED( rho_out ) ) DEALLOCATE(rho_out)
  DEALLOCATE( stau0 )
  DEALLOCATE( svel0 )
  DEALLOCATE( force0 )

  CALL mp_global_end ()
  STOP
END PROGRAM fpmd_postproc

!
!
!


SUBROUTINE read_fpmd( lforces, lcharge, lbinary, cunit, punit, funit, dunit, &
                      natoms, nr1, nr2, nr3, ispin, at, tau, force, &
                      rho, prefix, tmp_dir, ndr, charge_density )

  USE kinds,      ONLY: DP
  USE constants,  ONLY: bohr => BOHR_RADIUS_ANGS
  USE xml_io_base
  USE iotk_module

  IMPLICIT NONE

  LOGICAL, INTENT(in)   :: lforces, lcharge, lbinary
  INTEGER, INTENT(in)   :: cunit, punit, funit, dunit
  INTEGER, INTENT(in)   :: natoms, nr1, nr2, nr3, ispin, ndr
  REAL(DP), INTENT(out) :: at(3, 3), tau(3, natoms), force(3, natoms)
  REAL(DP), INTENT(out) :: rho(nr1, nr2, nr3)
  CHARACTER(LEN=*), INTENT(IN) :: prefix
  CHARACTER(LEN=*), INTENT(IN) :: tmp_dir
  CHARACTER(LEN=*), INTENT(IN) :: charge_density

  INTEGER  :: i, j, ix, iy, iz
  REAL(DP) :: rhomin, rhomax, rhof
  REAL(DP) :: x, y, z, fx, fy, fz
  CHARACTER(LEN=256) :: filename
  INTEGER       :: n1, n2, n3
  REAL(DP), ALLOCATABLE :: rho_plane(:)

  ! read cell vectors
  ! NOTE: colums are lattice vectors
  !
  READ(cunit,*)
  DO i = 1, 3
     READ(cunit,*) ( at(i, j), j=1,3 )
  END DO
  at(:, :) = at(:, :) * bohr

  ! read atomic coordinates
  READ(punit,*)
  IF (lforces) READ(funit,*)
  DO i = 1, natoms
     ! convert atomic units to Angstroms
     READ(punit,*) x, y, z
     tau(1, i) = x * bohr
     tau(2, i) = y * bohr
     tau(3, i) = z * bohr

     IF (lforces) THEN
        ! read forces
        READ (funit,*) fx, fy, fz
        force(1, i) = fx
        force(2, i) = fy
        force(3, i) = fz
     END IF
  END DO

  IF (lcharge) THEN

     filename = restart_dir( tmp_dir, ndr )
     !
     IF( charge_density == 'spin' ) THEN
        filename = TRIM( filename ) // '/' // 'spin-polarization'
     ELSE
        filename = TRIM( filename ) // '/' // 'charge-density'
     END IF
     !
     !
     IF ( check_file_exst ( TRIM(filename)//'.dat' ) ) THEN
        !
        CALL read_density( TRIM(filename)//'.dat', dunit, nr1, nr2, nr3, rho, lbinary )
        !
     ELSEIF ( check_file_exst ( TRIM(filename)//'.xml' ) ) THEN
        !
        CALL read_density( TRIM(filename)//'.xml', dunit, nr1, nr2, nr3, rho, lbinary )
        !
     ELSE         
        CALL infomsg ('read_fpmd', 'file '//TRIM(filename)//' not found' )
     ENDIF
     !
  END IF

  RETURN
END SUBROUTINE read_fpmd



SUBROUTINE read_density( filename, dunit, nr1, nr2, nr3, rho, lbinary )

   USE kinds,      ONLY: DP
   USE xml_io_base
   USE iotk_module

   IMPLICIT NONE

   LOGICAL, INTENT(in)   :: lbinary
   INTEGER, INTENT(in)   :: dunit
   INTEGER, INTENT(in)   :: nr1, nr2, nr3
   REAL(DP), INTENT(out) :: rho(nr1, nr2, nr3)
   CHARACTER(LEN=*), INTENT(IN) :: filename

   INTEGER  :: ix, iy, iz, ierr
   REAL(DP) :: rhomin, rhomax, rhof
   INTEGER       :: n1, n2, n3
   REAL(DP), ALLOCATABLE :: rho_plane(:)

     !
     WRITE(*,'("Reading density from: ", A80)' ) TRIM( filename ) 
     !
     CALL iotk_open_read( dunit, file = TRIM( filename ) , BINARY = lbinary, ROOT = attr, IERR = ierr )
     !

     CALL iotk_scan_begin( dunit, "CHARGE-DENSITY" )
     CALL iotk_scan_empty( dunit, "INFO", attr )
     CALL iotk_scan_attr( attr, "nr1",   n1 )
     CALL iotk_scan_attr( attr, "nr2",   n2 )
     CALL iotk_scan_attr( attr, "nr3",   n3 )
         !
     ALLOCATE( rho_plane( n1 * n2 ) )

     ! read charge density from file
     ! note: must transpose
     DO iz = 1, n3
       CALL iotk_scan_dat( dunit, "z" // iotk_index( iz ), rho_plane )
       IF( iz <= nr3 ) THEN
          DO iy = 1, MIN( n2, nr2 )
             DO ix = 1, MIN( n1, nr1 )
                rho(ix, iy, iz) = rho_plane( ix + ( iy - 1 ) * n1 )
             END DO
          END DO
       END IF
     END DO

     CALL iotk_scan_end( dunit, "CHARGE-DENSITY" )
     CALL iotk_close_read( dunit )

     rhomin = MINVAL(rho(:,:,:))
     rhomax = MAXVAL(rho(:,:,:))

     ! print some info
     WRITE(*,'(2x,"Density grid:")')
     WRITE(*,'(3(2x,i6))') nr1, nr2, nr3
     WRITE(*,'(2x,"spin = ",A4)') filename
     WRITE(*,'(2x,"Minimum and maximum values:")')
     WRITE(*,'(3(2x,1pe12.4))') rhomin, rhomax

   RETURN
END SUBROUTINE read_density

!
!
!
! compute inverse of 3*3 matrix
!
SUBROUTINE inverse( at, atinv )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  REAL(DP), INTENT(in)  :: at(3, 3)
  REAL(DP), INTENT(out) :: atinv(3, 3)

  REAL(DP) :: det

  atinv(1, 1) = at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2)
  atinv(2, 1) = at(2, 3) * at(3, 1) - at(2, 1) * at(3, 3)
  atinv(3, 1) = at(2, 1) * at(3, 2) - at(2, 2) * at(3, 1)
  atinv(1, 2) = at(1, 3) * at(3, 2) - at(1, 2) * at(3, 3)
  atinv(2, 2) = at(1, 1) * at(3, 3) - at(1, 3) * at(3, 1)
  atinv(3, 2) = at(1, 2) * at(3, 1) - at(1, 1) * at(3, 2)
  atinv(1, 3) = at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2)
  atinv(2, 3) = at(1, 3) * at(2, 1) - at(1, 1) * at(2, 3)
  atinv(3, 3) = at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)

  det = at(1, 1) * atinv(1, 1) + at(1, 2) * atinv(2, 1) + &
        at(1, 3) * atinv(3, 1)
  atinv(:,:) = atinv(:,:) / det;

  RETURN
END SUBROUTINE inverse

! generate cell dimensions and Euler angles from cell vectors
! euler(1:6) = a, b, c, alpha, beta, gamma
! I didn't call the array "celldm" because that could be confusing,
! since in PWscf the convention is different:
! celldm(1:6) = a, b/a, c/a, cos(alpha), cos(beta), cos(gamma)
SUBROUTINE at_to_euler( at, euler )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  REAL(DP), INTENT(in)  :: at(3, 3)
  REAL(DP), INTENT(out) :: euler(6)

  REAL(DP), PARAMETER :: rad2deg = 180.0d0 / 3.14159265358979323846d0
  REAL(DP) :: dot(3, 3)
  INTEGER :: i, j

  DO i = 1, 3
     DO j = i, 3
        dot(i, j) = dot_product(at(:,i), at(:,j))
     END DO
  END DO
  DO i = 1, 3
     euler(i) = sqrt(dot(i, i))
  END DO
  euler(4) = acos(dot(2, 3) / (euler(2) * euler(3))) * rad2deg
  euler(5) = acos(dot(1, 3) / (euler(1) * euler(3))) * rad2deg
  euler(6) = acos(dot(1, 2) / (euler(1) * euler(2))) * rad2deg

  RETURN
END SUBROUTINE at_to_euler

! generate cell vectors back from cell dimensions and Euler angles
! euler(1:6) = a, b, c, alpha, beta, gamma
! here I follow the PDB convention, namely, c is oriented along the z
! axis and b lies in the yz plane, or to put it another way, at is
! lower triangular
SUBROUTINE euler_to_at( euler, at )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  REAL(DP), PARAMETER :: deg2rad = 3.14159265358979323846d0 / 180.0d0

  REAL(DP), INTENT(in)  :: euler(6)
  REAL(DP), INTENT(out) :: at(3, 3)

  REAL(DP) :: cos_ab, cos_ac, cos_bc, temp1, temp2

  cos_bc = COS(euler(4) * deg2rad)
  cos_ac = COS(euler(5) * deg2rad)
  cos_ab = COS(euler(6) * deg2rad)

  temp1 = SQRT(1.0d0 - cos_bc*cos_bc) ! sin_bc
  temp2 = (cos_ab - cos_bc*cos_ac) / temp1

  at(1, 1) = SQRT(1.0d0 - cos_ac*cos_ac - temp2*temp2) * euler(1)
  at(2, 1) = temp2 * euler(1)
  at(3, 1) = cos_ac * euler(1)
  at(1, 3) = 0.0d0
  at(2, 3) = 0.0d0
  at(3, 3) = euler(3)
  at(1, 2) = 0.0d0
  at(2, 2) = temp1 * euler(2)
  at(3, 2) = cos_bc * euler(2)

  RETURN
END SUBROUTINE euler_to_at

! map charge density from a grid to another by linear interpolation
! along the three axes
SUBROUTINE scale_charge( rho_in, rho_out, nr1, nr2, nr3, ns1, ns2, ns3, &
                         np1, np2, np3 )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  INTEGER, INTENT(in)        :: nr1, nr2, nr3, ns1, ns2, ns3, np1, np2, np3
  REAL(DP), INTENT(in)  :: rho_in( nr1, nr2, nr3 )
  REAL(DP), INTENT(out) :: rho_out( ns1, ns2, ns3 )

  INTEGER       :: i, j, k
  INTEGER       :: i0(ns1), j0(ns2), k0(ns3), i1(ns1), j1(ns2), k1(ns3)
  REAL(DP) :: x0(ns1), y0(ns2), z0(ns3), x1(ns1), y1(ns2), z1(ns3)

  ! precompute interpolation data
  DO i = 1, ns1
     CALL scale_linear( i, nr1, ns1, np1, i0(i), i1(i), x0(i), x1(i) )
  END DO
  DO j = 1, ns2
     CALL scale_linear( j, nr2, ns2, np2, j0(j), j1(j), y0(j), y1(j) )
  END DO
  DO k = 1, ns3
     CALL scale_linear( k, nr3, ns3, np3, k0(k), k1(k), z0(k), z1(k) )
  END DO

  ! interpolate linearly along three axes
  DO i = 1, ns1
     DO j = 1, ns2
        DO k = 1, ns3
           rho_out(i, j, k) = &
              rho_in(i1(i), j1(j), k1(k)) * x0(i) * y0(j) * z0(k) + &
              rho_in(i0(i), j1(j), k1(k)) * x1(i) * y0(j) * z0(k) + &
              rho_in(i1(i), j0(j), k1(k)) * x0(i) * y1(j) * z0(k) + &
              rho_in(i1(i), j1(j), k0(k)) * x0(i) * y0(j) * z1(k) + &
              rho_in(i0(i), j0(j), k1(k)) * x1(i) * y1(j) * z0(k) + &
              rho_in(i0(i), j1(j), k0(k)) * x1(i) * y0(j) * z1(k) + &
              rho_in(i1(i), j0(j), k0(k)) * x0(i) * y1(j) * z1(k) + &
              rho_in(i0(i), j0(j), k0(k)) * x1(i) * y1(j) * z1(k)
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE scale_charge

! compute grid parameters for linear interpolation
SUBROUTINE scale_linear( n, nr, ns, np, n0, n1, r0, r1 )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  INTEGER, INTENT(in)        :: n, nr, ns, np
  INTEGER, INTENT(out)       :: n0, n1
  REAL(DP), INTENT(out) :: r0, r1

  ! map new grid point onto old grid
  ! mapping is: 1 --> 1, ns+1 --> (nr*np)+1
  r0 = REAL((n-1) * nr*np, DP) / ns + 1.0d0
  ! indices of neighbors
  n0 = int(r0)
  n1 = n0 + 1
  ! distances from neighbors
  r0 = r0 - n0
  r1 = 1.0d0 - r0
  ! apply periodic boundary conditions
  n0 = MOD(n0 - 1, nr) + 1
  n1 = MOD(n1 - 1, nr) + 1

  RETURN
END SUBROUTINE scale_linear

SUBROUTINE write_xsf( ldynamics, lforces, lcharge, ounit, n, at, &
                      natoms, ityp, tau, force, rho, nr1, nr2, nr3 )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  LOGICAL, INTENT(in)       :: ldynamics, lforces, lcharge
  INTEGER, INTENT(in)       :: ounit, n, natoms, ityp(natoms)
  INTEGER, INTENT(in)       :: nr1, nr2, nr3
  REAL(DP), INTENT(in) :: at(3, 3), tau(3, natoms), force(3, natoms)
  REAL(DP), INTENT(in) :: rho(nr1, nr2, nr3)

  INTEGER :: i, j, ix, iy, iz

  ! write cell
  IF (ldynamics) THEN
     WRITE(ounit,*) 'PRIMVEC', n
  ELSE
     WRITE(ounit,*) 'PRIMVEC'
  END IF
  WRITE(ounit,'(2(3f15.9/),3f15.9)') at
  IF (ldynamics) THEN
     WRITE(ounit,*) 'CONVVEC', n
     WRITE(ounit,'(2(3f15.9/),3f15.9)') at
  END IF

  ! write atomic coordinates (and forces)
  IF (ldynamics) THEN
     WRITE(ounit,*) 'PRIMCOORD', n
  ELSE
     WRITE(ounit,*) 'PRIMCOORD'
  END IF
  WRITE(ounit,*) natoms, 1
  DO i = 1, natoms
     IF (lforces) THEN
        WRITE (ounit,'(i3,3x,3f15.9,1x,3f12.5)') ityp(i), &
              (tau(j, i), j=1,3), (force(j, i), j=1,3)
     ELSE
        WRITE (ounit,'(i3,3x,3f15.9,1x,3f12.5)') ityp(i), &
              (tau(j, i), j=1,3)
     END IF
  END DO

  ! write charge density
  IF (lcharge) THEN
     ! XSF scalar-field header
     WRITE(ounit,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
     WRITE(ounit,'(a)') '3D_PWSCF'
     WRITE(ounit,'(a)') 'DATAGRID_3D_UNKNOWN'

     ! mesh dimensions
     WRITE(ounit,*) nr1, nr2, nr3
     ! origin
     WRITE(ounit,'(3f10.6)') 0.0d0, 0.0d0, 0.0d0
     ! lattice vectors
     WRITE(ounit,'(3f10.6)') ((at(i, j), i=1,3), j=1,3)
     ! charge density
     WRITE(ounit,'(6e13.5)') &
          (((rho(ix, iy, iz), ix=1,nr1), iy=1,nr2), iz=1,nr3)

     WRITE(ounit,'(a)') 'END_DATAGRID_3D'
     WRITE(ounit,'(a)') 'END_BLOCK_DATAGRID_3D'
  END IF

  RETURN
END SUBROUTINE write_xsf


SUBROUTINE write_xyz( ldynamics, lforces, lcharge, ounit, n, at, &
                      natoms, ityp, tau, force, rho, nr1, nr2, nr3 )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  LOGICAL, INTENT(in)       :: ldynamics, lforces, lcharge
  INTEGER, INTENT(in)       :: ounit, n, natoms, ityp(natoms)
  INTEGER, INTENT(in)       :: nr1, nr2, nr3
  REAL(DP), INTENT(in) :: at(3, 3), tau(3, natoms), force(3, natoms)
  REAL(DP), INTENT(in) :: rho(nr1, nr2, nr3)
  INTEGER :: i, j, ix, iy, iz
  CHARACTER*2 :: label(103)
  DATA label /" H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", &
              "Na", "Mg", "Al", "Si", " P", " S", "Cl", "Ar", " K", "Ca", &
              "Sc", "Ti", " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
              "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", &
              "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", &
              "Sb", "Te", " I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
              "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
              "Lu", "Hf", "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
              "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
              "Pa", " U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
              "Md", "No", "Lr"/

  ! write natoms
  write(ounit,*) natoms

  ! write cell
  write(ounit,'(9(F10.4,2X))') at

  ! write atomic coordinates (and forces)
  DO i = 1, natoms
     IF (lforces) THEN
        WRITE (ounit,'(a2,3x,3f15.9,1x,3f12.5)') label(ityp(i)), &
              (tau(j, i), j=1,3), (force(j, i), j=1,3)
     ELSE
        WRITE (ounit,'(a2,3x,3f15.9,1x,3f12.5)') label(ityp(i)), &
              (tau(j, i), j=1,3)
     END IF
  END DO

  RETURN
END SUBROUTINE write_xyz


SUBROUTINE write_grd( ounit, at, rho, nr1, nr2, nr3 )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  INTEGER, INTENT(in)       :: ounit
  INTEGER, INTENT(in)       :: nr1, nr2, nr3
  REAL(DP), INTENT(in) :: at(3, 3), rho(nr1, nr2, nr3)

  INTEGER       :: i, j, k
  REAL(DP) :: euler(6)

  CALL at_to_euler( at, euler )

  WRITE(ounit,*) 'charge density'
  WRITE(ounit,*) '(1p,e12.5)'
  WRITE(ounit,fmt='(6f9.3)') (euler(i), i=1,6)
  WRITE(ounit,fmt='(3i5)') nr1 - 1, nr2 - 1, nr3 - 1
  WRITE(ounit,fmt='(7i5)') 1, 0, 0, 0, nr1 - 1, nr2 - 1, nr3 - 1
  WRITE(ounit,fmt='(1p,e12.5)') (((rho(i, j, k), i=1,nr1), j=1,nr2), k=1,nr3)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE write_pdb( bunit, at, tau, natoms, ityp, euler, lrotation )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  INTEGER, INTENT(in)       :: bunit, natoms
  INTEGER, INTENT(in)       :: ityp(natoms)
  REAL(DP), INTENT(in) :: at(3, 3), tau(3, natoms), euler(6)
  LOGICAL, INTENT(in)       :: lrotation

  INTEGER     :: i, j
  CHARACTER*2 :: label(103)
  DATA label /" H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", &
              "Na", "Mg", "Al", "Si", " P", " S", "Cl", "Ar", " K", "Ca", &
              "Sc", "Ti", " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
              "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", &
              "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", &
              "Sb", "Te", " I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
              "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
              "Lu", "Hf", "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
              "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
              "Pa", " U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
              "Md", "No", "Lr"/

  WRITE(bunit,'("HEADER    PROTEIN")')
  WRITE(bunit,'("COMPND    UNNAMED")')
  WRITE(bunit,'("AUTHOR    GENERATED BY ...")')

  IF (lrotation) &
     WRITE(bunit,'("CRYST1",3F9.3,3F7.2,1X,A10,I3)') euler, "P 1", 1

  DO i = 1, natoms
     WRITE(bunit,'("ATOM  ",I5,1X,A2,3X,2A3,I3,3X,F9.3,2F8.3,2F6.2," ")') &
           i, label(ityp(i)), "UKN", "", 1, (tau(j, i), j=1,3), 1.0d0, 0.0d0
  END DO

  WRITE(bunit,'("MASTER        0    0    0    0    0    0    0    0 ", I4,"    0 ",I4,"    0")') natoms, natoms
  WRITE(bunit,'("END")')

  RETURN
END SUBROUTINE write_pdb

! PDB File Format
!---------------------------------------------------------------------------
!Field |    Column    | FORTRAN |                                         
!  No. |     range    | format  | Description                                   
!---------------------------------------------------------------------------
!   1. |    1 -  6    |   A6    | Record ID (eg ATOM, HETATM)       
!   2. |    7 - 11    |   I5    | Atom serial number                            
!   -  |   12 - 12    |   1X    | Blank                                         
!   3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")   
!   4. |   17 - 17    |   A1    | Alternative location code (if any)            
!   5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue 
!   -  |   21 - 21    |   1X    | Blank                                         
!   6. |   22 - 22    |   A1    | Chain identifier code                         
!   7. |   23 - 26    |   I4    | Residue sequence number                       
!   8. |   27 - 27    |   A1    | Insertion code (if any)                       
!   -  |   28 - 30    |   3X    | Blank                                         
!   9. |   31 - 38    |  F8.3   | Atom's x-coordinate                         
!  10. |   39 - 46    |  F8.3   | Atom's y-coordinate                         
!  11. |   47 - 54    |  F8.3   | Atom's z-coordinate                         
!  12. |   55 - 60    |  F6.2   | Occupancy value for atom                      
!  13. |   61 - 66    |  F6.2   | B-value (thermal factor)                    
!   -  |   67 - 67    |   1X    | Blank                                         
!  14. |   68 - 68    |   I3    | Footnote number                               
!---------------------------------------------------------------------------



