!
! Copyright (C) 2001-2005 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
!
! Auxiliary functions for the vibrational analysis
!
!-------------------------------------------------------------------------

SUBROUTINE print_matrix (matrix, dim1, dim2, title, filep)
  !
  ! ... modules
  !
  USE kinds,                   ONLY: DP
  !
  IMPLICIT NONE
  !
  ! ... input variables
  !
  INTEGER,             INTENT(IN) :: dim1, dim2, filep
  CHARACTER (len=*),   INTENT(IN) :: title
  REAL      (KIND=DP), INTENT(IN) :: matrix(dim1,dim2)
  !
  ! ... local variables
  !
  INTEGER                         :: i,j
  !
  ! ...
  !
  WRITE (filep,*)
  WRITE (filep,*) TRIM(title)
  WRITE (filep,*)
  DO i=1,dim1
     DO j=1,dim2
        WRITE (filep,FMT='(F12.6,3X)',ADVANCE='NO') matrix(i,j)
     END DO
     WRITE (filep,*)
  END DO
  WRITE (filep,*)
  !
  RETURN
END SUBROUTINE print_matrix
!
!
!  ----------------------------------------------
!
!
SUBROUTINE write_xyz(tau,free_text,position_in_file,file_name)
  !
  ! Printout of the coordinates in an xyz format
  ! tau              - coordinates
  ! free_text        - text for the second line of output
  ! position_in_file - either 'APPEND' (for an xyz animation),
  !                    'REWIND' to overwrite, or 'ASIS'
  !
  ! ... modules
  !
  USE constants,              ONLY : bohr_radius_angs
  USE input_parameters,       ONLY : atom_label
  USE io_files,               ONLY : iunxyz
  USE ions_base,              ONLY : nat, ityp
  USE kinds,                  ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... input variables
  !
  CHARACTER (LEN=*),   INTENT(IN) :: position_in_file, free_text
  CHARACTER (LEN=*),   INTENT(IN) :: file_name
  REAL      (KIND=DP), INTENT(IN) :: tau(3,nat)
  !
  ! ... local variables
  !
  INTEGER                         :: atom
  CHARACTER (LEN=*), PARAMETER    :: xyz_fmt  = "(A2,3(2X,F14.10))"

  OPEN( UNIT = iunxyz, FILE = file_name, &
       STATUS = "UNKNOWN", ACTION = "WRITE", &
       POSITION = TRIM(position_in_file))
  !
  WRITE( UNIT = iunxyz, FMT = '(I5)' ) nat
  WRITE( UNIT = iunxyz, FMT = '(A)'  ) free_text
  !
  DO atom = 1, nat
     !
     WRITE( UNIT = iunxyz, FMT = xyz_fmt ) &
          TRIM( atom_label( ityp( atom ) ) ), &
          tau(1,atom)*bohr_radius_angs, &
          tau(2,atom)*bohr_radius_angs, &
          tau(3,atom)*bohr_radius_angs
     !
  END DO
  !
  CLOSE(iunxyz)
  !
  RETURN
END SUBROUTINE write_xyz
!
!
!  ----------------------------------------------
!
!
SUBROUTINE calculate_dipole (dipole, dipole_moment,tau)
  !
  ! Driver routine for the calculation of the dipole moment.
  ! Currently it uses only subroutine 'poles' from cplib.f90
  ! It could easily be extended to use Wannier functions as well.
  !
  ! ... modules
  !
  USE cp_main_variables,    ONLY : irb, eigrb, bec, rhor, rhog, rhos
  USE electrons_base,       ONLY : nspin
  USE energies,             ONLY : ekin, enl
  USE ions_base,            ONLY : nat
  USE kinds,                ONLY : DP
  USE wavefunctions_module, ONLY : c0
  USE uspp,                 ONLY : becsum
  !
  IMPLICIT NONE
  !
  ! ... input variables
  !
  REAL (KIND=DP),   INTENT(IN)  :: tau(3,nat)
  !
  ! ... output variables
  !
  REAL (KIND=DP),   INTENT(OUT) :: dipole(3), dipole_moment
  !
  ! ... local variables
  !
  REAL (KIND=DP)                :: quadrupole
  INTEGER                       :: nfi=10 ! ... dummy value
  LOGICAL                       :: ion_flag = .TRUE., coc_flag=.TRUE.
  !
  !
  !
  dipole        = 0.0
  dipole_moment = 0.0
  rhor          = 0.0
  rhog          = 0.0
  rhos          = 0.0
  !
  CALL rhoofr(nfi,c0,irb,eigrb,bec,becsum,rhor,rhog,rhos,enl,ekin)
  IF (nspin.EQ.2) THEN
     rhor(:,1)=rhor(:,1)+rhor(:,2)
  END IF
  CALL poles (dipole_moment,dipole,quadrupole,rhor(:,1),&
       ion_flag,tau,coc_flag)
  !
  ! STILL HAVE TO ADD IN THE WANNIER BASED DIPOLE
  ! AS AN ALTERNATIVE WAY OF DOING THE CALCULATION
  !
  RETURN
END SUBROUTINE calculate_dipole
!
!
!  ----------------------------------------------
!
!
SUBROUTINE orthonormalize(M,dim1,dim2)
  !
  ! modified Gram-Schimd - orthogonalizing for 
  ! dim1 vectors of length dim2
  !
  !
  ! ... modules
  !
  USE kinds,             ONLY : DP
  !
  IMPLICIT NONE
  !
  !... input variables
  !
  INTEGER,      INTENT(IN)   :: dim1,dim2         !dimensions of input matrix
  REAL(KIND=DP),INTENT(INOUT):: M(dim1,dim2)      !input matrix
  !
  !... local variables
  !
  REAL (KIND=DP)             :: t(dim2,dim2), z
  INTEGER                    :: i,j,k
  !
  t=0.0
  !
  DO k=1,dim2
     !
     ! ... normalizing vector k
     z=0.0
     DO  i=1,dim1
        z=z+M(i,k)**2
     END DO
     t(k,k)=SQRT(z)
     M(:,k)=M(:,k)/t(k,k)
     !
     ! ... removing the projection of vector k
     ! ... from the remaining vectors k+1...dim2
     !
     DO j=k+1,dim2
        z=0.0
        DO i=1,dim1
           z=z+M(i,j)*M(i,k)
        END DO
        !
        t(k,j)=z
        !
        DO i=1,dim1
           M(i,j)=M(i,j)-t(k,j)*M(i,k)
        END DO
     END DO
     !
  END DO
  !
  RETURN
END SUBROUTINE orthonormalize
!
!
!  ----------------------------------------------
!
!
SUBROUTINE symmetrize_matrix(M, dim)
  !
  ! IMPOSING SYMMETRY ON MATRIX M
  ! M(i,j)==M(j,i)
  !
  ! ... modules
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... input/output variables
  !
  INTEGER,        INTENT(IN)    :: dim
  REAL (KIND=DP), INTENT(INOUT) :: M(dim,dim)
  !
  ! ... local variables
  !
  INTEGER                       :: i,j
  REAL (KIND=DP)                :: sum
  !
  !
  DO i=1,dim
     DO j = i+1,dim
        sum      = (M(i,j)+M(j,i))/2
        M(i,j) = sum
        M(j,i) = sum
     END DO
  END DO
  !
  RETURN
END SUBROUTINE symmetrize_matrix
  !
  !
  !  ----------------------------------------------
  !
  !
SUBROUTINE relax_wavefunction (fion)
  !
  ! A DRIVER ROUTINE FOR WAVE FUNCTION RELAXATION
  !
  ! ... modules
  !
  USE cell_base,            ONLY : b1, b2, b3
  USE cg_module,            ONLY : tcg
  USE control_flags,        ONLY : tprnfor, thdyn
  USE cp_electronic_mass,   ONLY : emass
  USE cp_main_variables,    ONLY : eigr, ei1, ei2, ei3, &
       sfac, irb, eigrb, taub, nfi
  USE efield_module,        ONLY : tefield, efield_update
  USE wave_base,            ONLY : frice
  USE energies,             ONLY : eself, etot
  USE from_scratch_module,  ONLY : from_scratch
  USE gvecs,                ONLY : ngs
  USE ions_positions,       ONLY : tau0
  USE kinds,                ONLY : DP
  USE parameters,           ONLY : natx
  USE phase_factors_module, ONLY : strucf     
  USE reciprocal_vectors,   ONLY : mill_l
  USE time_step,            ONLY : dt2
  !
  IMPLICIT NONE
  !
  ! ... input variables
  !
  REAL (KIND=DP)                ::  fion(3,natx)
  !
  ! ... local variables
  !
  LOGICAL                       :: tfirst, tlast
  REAL (KIND=DP)                :: enthal, fccc
  REAL (KIND=DP)                :: dt2bye, enb, enbi, ccc
  !
  ! ... for smooth restart in the new coordinates
  ! 
  dt2bye = dt2 / emass
  fccc = 1.D0 / ( 1.D0 + frice )
  CALL initbox( tau0, taub, irb )
  CALL phbox( taub, eigrb )
  CALL phfac( tau0, ei1, ei2, ei3, eigr )
  CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
  IF ( thdyn ) CALL formf( tfirst, eself )
  IF (tefield ) CALL efield_update( eigr )
  !
  ! ... relax wavefunction in new position
  !
  IF (tcg) THEN
     tprnfor = .TRUE. ! ... atomic forces are calculated only at the
     !     end of the wavefunction relaxation process
     CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, fion, &
          enthal, enb, enbi, fccc, ccc, dt2bye )
     tprnfor = .FALSE.
  ELSE
     tprnfor = .FALSE. ! ... no need to calculate atomic forces at 
     !     each step of the MD damped dynamics
     CALL cprmain ( tau0, fion, etot )
     !
     ! ... one more scf step to get the forces on the nuclei
     !
     tprnfor = .TRUE.
     CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, fion, &
          enthal, enb, enbi, fccc, ccc, dt2bye )
     tprnfor = .FALSE.
  END IF
  !
  RETURN
END SUBROUTINE relax_wavefunction
