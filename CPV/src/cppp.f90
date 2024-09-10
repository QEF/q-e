!
! Copyright (C) 2002-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! --------------------------------------------------------------------
! this routine writes the crystal structure in XSF, GRD and PDB format
! from CP output files
! --------------------------------------------------------------------
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
! modified by Gerardo Ballabio and Carlo Cavazzoni
! on Thu Jul 22 18:57:26 CEST 2004
! Adapted to new XML format by Paolo Giannozzi, August 2019
!
PROGRAM cp_postproc

  USE kinds,      ONLY : DP
  USE io_files,   ONLY : prefix, tmp_dir, xmlfile
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE matrix_inversion, ONLY : invmat
  USE qes_types_module, ONLY : output_type
  USE qexsd_module,     ONLY : qexsd_readschema
  USE qexsd_copy,       ONLY : qexsd_copy_atomic_species, &
       qexsd_copy_atomic_structure, qexsd_copy_basis_set

  IMPLICIT NONE

  INTEGER, PARAMETER :: maxsp = 20
  !
  ! QE variables, read from files
  !
  TYPE (output_type)    :: output_obj 
  INTEGER, ALLOCATABLE  :: ityp(:)
  INTEGER               :: nat, nsp, ibrav
  CHARACTER(len=6)      :: atm( maxsp )
  REAL(DP)              :: alat, amass(maxsp)
  REAL(DP)              :: at(3, 3)
  REAL(DP), ALLOCATABLE :: tau(:,:), tau_in(:,:), tau_out(:,:)
  REAL(DP), ALLOCATABLE :: sigma(:,:), force(:,:)
  !
  ! local variables
  !
  CHARACTER(len=256) :: filepp, fileout, output, outdir
  CHARACTER(len=256) :: filecel, filepos, filefor, filepdb
  LOGICAL            :: lforces, ldynamics, lpdb, lrotation
  INTEGER            :: ounit, cunit, punit, funit, dunit, bunit, ksunit
  INTEGER            :: natoms, nframes
  INTEGER            :: atomic_number(maxsp)
  !! FIXME: should be deleted from input and replaced by QE function
  !! FIXME: with the same name: atomic_number(atm(ityp(n)))
  INTEGER            :: np1, np2, np3, np
  INTEGER            :: ios, ndr, ierr
  INTEGER, ALLOCATABLE  :: n_atomic(:)
  REAL(DP)           :: atinv(3,3)
  INTEGER            :: i, j, k, n, ix, iy, iz, int_dum

  REAL(DP) :: euler(6)

  NAMELIST /inputpp/ prefix, fileout, output, outdir, &
                     lforces, ldynamics, lpdb, lrotation, &
                     atomic_number, nframes, ndr

  ! default values

  dunit = 14
  !  initialize mpi
  CALL mp_startup  ( )
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix    = 'cp'
  fileout   = 'out'
  output    = 'xsf'
  lforces   = .false.
  ldynamics = .true.
  lpdb      = .false.
  lrotation = .false.
  np1   = 1
  np2   = 1
  np3   = 1                  ! 
  nframes = 1                ! number of MD step to be read to buind the trajectory
  ndr   = 51                 ! restart file number
  atomic_number = 1          ! atomic number of the species in the restart file

  call input_from_file()

  ! read namelist
  READ( 5, inputpp, iostat=ios)
  IF (ios /= 0) THEN
     WRITE(*,*) 'Error reading namelist &inputpp'
     STOP
  END IF

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
  IF ( nframes == 1 ) THEN
     WRITE(*,*) 'Error: single frame not implemented'
     STOP
  END IF
  !
  !  Now read the XML data file
  !

  filepp = xmlfile( ndr )
  CALL qexsd_readschema ( filepp, ierr, output_obj )
  IF( ierr > 0 ) CALL errore(' cppp ', ' Cannot open file '//TRIM(filepp), ierr)
  !
  !   End of reading from data file - now copy to variables
  !
  CALL qexsd_copy_atomic_species ( output_obj%atomic_species, &
       nsp, atm, amass )
  !
  IF ( nsp > maxsp ) THEN
     WRITE(*,*) 'Error: too many atomic species'
     STOP
  END IF
  !
  CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
       atm, nat, tau, ityp, alat, at(:,1), at(:,2), at(:,3), ibrav, int_dum)
  !
  !   End of reading from data file
  !
  ! allocate arrays
  ALLOCATE(tau_in(3, nat))                  ! atomic positions, angstroms
  ALLOCATE(tau_out(3, nat * np))            ! replicated positions
  ALLOCATE(sigma(3, nat ) )                 ! scaled coordinates
  ALLOCATE (n_atomic(nat*np))
  IF (lforces) ALLOCATE( force( 3, nat * np ) )
  !
  ! assign species (from input) to each atom
  !
  DO i = 1, nat
     n_atomic(i) = atomic_number(ityp(i))
  END DO

  ! open output file for trajectories 
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

     ! read data from files produced by cp
     !
     CALL read_cp( lforces, cunit, punit, funit, dunit, &
                     nat, at, tau_in, force, ierr ) 
     IF ( ierr /= 0 ) THEN
         WRITE(*,'("End of file after ",i5," frames")') n
         EXIT
     END IF

     WRITE(*,'(2x,"Cell parameters (Angstroms):")')
     WRITE(*,'(3(2x,f10.6))') ((at(i, j), i=1,3), j=1,3)
     !
     WRITE(*,'(2x,"Atomic coordinates (Angstroms):")')
     WRITE(*,'(3(2x,f10.6))') ((tau_in(i, j), i=1,3), j=1,nat)

     ! compute scaled coordinates
     !
     CALL invmat( 3, at, atinv )
     sigma(:,:) = MATMUL(atinv(:,:), tau_in(:,:))

     ! compute cell dimensions and Euler angles
     CALL at_to_euler( at, euler )

     IF (lpdb) THEN
        ! apply periodic boundary conditions
        DO i = 1, nat
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
              DO j = 1, nat
                 k = k + 1
                 tau_out(:, k) = tau_in(:, j) + (ix-1) * at(:, 1) + &
                                 (iy-1) * at(:, 2) + (iz-1) * at(:, 3)
                 n_atomic(k) = n_atomic(j)
                 IF (lforces) force(:, k) = force(:, j)
              END DO
           END DO
        END DO
     END DO
     natoms = nat * np

     ! compute supercell
     at(:, 1) = at(:, 1) * np1
     at(:, 2) = at(:, 2) * np2
     at(:, 3) = at(:, 3) * np3
     euler(1) = euler(1) * np1
     euler(2) = euler(2) * np2
     euler(3) = euler(3) * np3

     IF ( output == 'xsf' ) THEN
        ! write data as XSF format
        CALL write_xsf( ldynamics, lforces, ounit, n, at, & 
                        natoms, n_atomic, tau_out, force )
     ELSE IF( output == 'xyz' ) THEN
        ! write data as XYZ format
        CALL write_xyz( ldynamics, lforces, ounit, n, at, &
                        natoms, n_atomic, tau_out, force )
     END IF

  END DO

  CLOSE(ounit)

  ! write atomic positions as PDB format
  CALL write_pdb( bunit, at, tau_out, natoms, n_atomic, euler, lrotation )

  ! free allocated resources
  CLOSE(punit)
  CLOSE(cunit)
  IF (lforces) CLOSE(funit)

  DEALLOCATE(tau_in)
  DEALLOCATE(tau_out)
  DEALLOCATE(ityp)
  DEALLOCATE(n_atomic)
  IF( ALLOCATED( force  ) ) DEALLOCATE(force)

  CALL mp_global_end ()
  STOP
END PROGRAM cp_postproc

!
!
!


SUBROUTINE read_cp( lforces, cunit, punit, funit, dunit, &
                      nat, at, tau, force, ierr )

  USE kinds,      ONLY: DP
  USE constants,  ONLY: bohr => BOHR_RADIUS_ANGS
  USE io_files,   ONLY: check_file_exist, restart_dir

  IMPLICIT NONE

  LOGICAL, INTENT(in)   :: lforces
  INTEGER, INTENT(in)   :: cunit, punit, funit, dunit
  INTEGER, INTENT(in)   :: nat
  INTEGER, INTENT(out)  :: ierr
  REAL(DP), INTENT(out) :: at(3, 3), tau(3, nat), force(3, nat)

  INTEGER  :: i, j, ix, iy, iz
  REAL(DP) :: x, y, z, fx, fy, fz
  CHARACTER(LEN=256) :: filename
  INTEGER       :: n1, n2, n3

  ierr = 0
  ! read cell vectors
  ! NOTE: colums are lattice vectors
  !
  READ(cunit,*, end=10, err=10)
  DO i = 1, 3
  READ(cunit,*, end=10, err=10) ( at(i, j), j=1,3 )
  END DO
  at(:, :) = at(:, :) * bohr

  ! read atomic coordinates
  READ(punit,*, end=10, err=10)
  IF (lforces) READ(funit,*, end=10, err=10)
  DO i = 1, nat
     ! convert atomic units to Angstroms
     READ(punit,*, end=10, err=10) x, y, z
     tau(1, i) = x * bohr
     tau(2, i) = y * bohr
     tau(3, i) = z * bohr

     IF (lforces) THEN
        ! read forces
        READ (funit,*, end=10, err=10) fx, fy, fz
        force(1, i) = fx
        force(2, i) = fy
        force(3, i) = fz
     END IF
  END DO
  RETURN
  10 ierr=1
END SUBROUTINE read_cp

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

SUBROUTINE write_xsf( ldynamics, lforces, ounit, n, at, &
                      natoms, ityp, tau, force )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  LOGICAL, INTENT(in)       :: ldynamics, lforces
  INTEGER, INTENT(in)       :: ounit, n, natoms, ityp(natoms)
  REAL(DP), INTENT(in) :: at(3, 3), tau(3, natoms), force(3, natoms)

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

  RETURN
END SUBROUTINE write_xsf


SUBROUTINE write_xyz( ldynamics, lforces, ounit, n, at, &
                      natoms, ityp, tau, force )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  LOGICAL, INTENT(in)       :: ldynamics, lforces
  INTEGER, INTENT(in)       :: ounit, n, natoms, ityp(natoms)
  REAL(DP), INTENT(in) :: at(3, 3), tau(3, natoms), force(3, natoms)
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
  !! FIXME: should be replaced by QE function "atom_name"
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
  !! FIXME: should be replaced by QE function "atom_name"

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

