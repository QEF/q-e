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
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
#ifdef DFT_CP
  USE cp_main_variables,    ONLY : irb, eigrb, bec, rhor, &
                                   rhog, rhos
  USE electrons_base,       ONLY : nspin
  USE energies,             ONLY : ekin, enl
  USE wavefunctions_module, ONLY : c0
  USE uspp,                 ONLY : becsum
  USE grid_dimensions,      ONLY : nnrx
  USE charge_density,       ONLY : rhoofr
#endif
#ifdef DFT_PW
  USE scf,                  ONLY : rhor=>rho
  USE lsda_mod,             ONLY : nspin
  USE gvect,                ONLY : nrxx
  USE cell_base,            ONLY : alat
#endif
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
#ifdef DFT_PW
   REAL(DP)                     :: tmp_rho(nrxx,nspin)
#endif
#ifdef DFT_CP
   REAL(DP)                     :: tmp_rho(nnrx,nspin)
#endif
  !
  !
  !
  ! ... initiating variables
  dipole         = 0.0
  dipole_moment  = 0.0
  !
#ifdef DFT_CP
  !
  ! ... CP charge density
  !
  rhog           = 0.0
  rhos           = 0.0
  !
  CALL rhoofr(nfi,c0(:,:,1,1),irb,eigrb,bec,becsum,rhor,rhog,rhos,enl,ekin)
#endif

#ifdef DFT_PW
  !
  ! ... PW charge density
  !
#endif
  !
  tmp_rho = rhor
  IF (nspin.EQ.2) THEN
     tmp_rho(:,1)=tmp_rho(:,1)+tmp_rho(:,2)
  END IF
  !
  ! ... computing dipole
  !
#ifdef DFT_CP
  CALL poles (dipole_moment,dipole,quadrupole,tmp_rho(:,1),&
       ion_flag,tau,coc_flag)
#endif
#ifdef DFT_PW
  CALL poles (dipole_moment,dipole,quadrupole,tmp_rho(:,1),&
       ion_flag,tau*alat,coc_flag)
#endif
  !
  ! STILL HAVE TO ADD IN THE WANNIER BASED DIPOLE
  ! AS AN ALTERNATIVE WAY FOR CALCULATING THE DIPOLE
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
  !
  ! ... modules
  !
  USE parameters,           ONLY : natx
  USE kinds,                ONLY : DP
  !
#ifdef DFT_CP
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
  USE ions_base,            ONLY : ityp, nat
  USE phase_factors_module, ONLY : strucf     
  USE reciprocal_vectors,   ONLY : mill_l
  USE time_step,            ONLY : dt2
  USE io_global,            ONLY : stdout
#endif
  !
#ifdef DFT_PW
  USE force_mod,            ONLY : force
  USE ener,                 ONLY : etot
  USE cell_base,            ONLY : alat
  USE constants,            ONLY : e2
#endif
  IMPLICIT NONE
  !
  ! ... input variables
  !
  REAL (KIND=DP)                ::  fion(3,natx)
  !
#ifdef DFT_CP
  !
  ! ... local variables
  !
  LOGICAL                       :: tfirst, tlast
  REAL (KIND=DP)                :: enthal, fccc
  REAL (KIND=DP)                :: dt2bye, enb, enbi, ccc
  integer                       :: ipol, na
  !
  ! ... for smooth restart in the new coordinates
  ! 
  fion = 0.0
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
  ! ... write on output the forces
  !
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), ( fion(ipol,na), ipol = 1, 3 )
  END DO
9035 FORMAT(5X,'atom ',I3,' type ',I2,'   force = ',3F14.8)

  !
#endif
  !
#ifdef DFT_PW
  fion = 0.0
  force = 0.0
  !call hinit0 ()
  CALL hinit1 ()
  call electrons()
  call forces()
  ! convert to hartree/bohr
  fion = force /e2
  etot = etot / e2
#endif
  !
  RETURN
END SUBROUTINE relax_wavefunction
!
!  ----------------------------------------------
!
subroutine set_guess_wfc ( disp_sign )
  !
  USE kinds,                ONLY : DP
#ifdef DFT_CP
  USE ions_positions,       ONLY : tau0
  USE cp_main_variables,    ONLY : lambda, lambdam, nfi
  USE wavefunctions_module, ONLY : c0, cm
  USE vibrations,           ONLY : ref_c0
#endif
#ifdef DFT_PW
  USE scf,                  ONLY : rho
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, iunoldwfc2, prefix
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE klist,                ONLY : nks
  USE scf,                  ONLY : rho, rho_core, vr
  USE gvect,                ONLY : nrxx, ngm, g, gg, gstart,  & 
                                   nr1, nr2, nr3, nl, &
                                   eigts1, eigts2, eigts3, & 
                                   nrx1, nrx2, nrx3
  USE cell_base,            ONLY : omega, bg, alat
  USE ener,                 ONLY : ehart, etxc, vtxc
  USE extfield,             ONLY : etotefield
  USE vlocal,               ONLY : strf

#endif
  !
  ! ... input variables
  !
  INTEGER, INTENT(IN)           :: disp_sign
  !
  ! ... local variables
  !
  logical                       :: exst
  REAL(kind=DP)                 :: charge
  !
  !
#ifdef DFT_CP
  nfi = 0
#endif
  !
  !

  IF(disp_sign.EQ.-1) THEN
     !
     ! ... use reference wavefunction as initial guess
     !
#ifdef DFT_CP
     c0(:,:,1,1) = ref_c0(:,:,1,1)
#endif
     !
#ifdef DFT_PW
     INQUIRE (file=TRIM(prefix)//'.old2rho', EXIST=exst)
     if ( exst ) then
        !
        ! ... charge-density from file to memory
        !  
        CALL io_pot( - 1, 'old2rho', rho, 1 )
        !
        CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3,   &
             nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
             ehart, etxc, vtxc, etotefield, charge, vr )
        !
        ! .. swapping the wavefunctions
        !
        CALL diropn( iunoldwfc2, 'old2wfc', nwordwfc, exst )
        !
        DO ik = 1, nks
           !
           ! ... "2old"  -> "now"
           !
           CALL davcio( evc, nwordwfc, iunoldwfc2, ik, - 1 )
           CALL davcio( evc, nwordwfc, iunwfc,     ik, + 1 )
           !
        END DO
        !
        CLOSE( UNIT = iunoldwfc2 )
        !
     end if
     history = 1
     call update_pot()
     !
#endif
     !
  ELSE
     !
     ! ... extrapolate wave function coefficients from their value
     ! ... at dispalcement at the negative direction
     !
#ifdef DFT_CP
     c0(:,:,1,1) = 2*ref_c0(:,:,1,1)-c0(:,:,1,1) 
#endif
#ifdef DFT_PW
     history = 2
     call update_pot()
#endif
     !
  END IF
  !
  ! ... set wavefunction velocity to zero
  !
#ifdef DFT_CP
  cm      = c0
  lambdam = lambda
#endif
  !
  return
end subroutine set_guess_wfc
!
!
!  ----------------------------------------------
!
!
SUBROUTINE cofmass( tau, mass, nat, com )
  !
  ! ... calculation of the center of mass
  !
  USE kinds,        only : DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN)  :: tau(3,nat), mass(nat)
  REAL(DP), INTENT(OUT) :: com(3)
  INTEGER,  INTENT(IN)  :: nat

  REAL(DP) :: tmas
  INTEGER :: is, i, ia
  !
  tmas=SUM(mass(1:nat))
  !
  do i=1,3
     com(i)=DOT_PRODUCT(tau(i,1:nat),mass(1:nat))
     com(i)=com(i)/tmas
  end do
  !
  RETURN
END SUBROUTINE cofmass

!-----------------------------------------------------------------------
SUBROUTINE vib_pbc(rin,a1,a2,a3,ainv,rout)
  !-----------------------------------------------------------------------
  !
  !     brings atoms inside the unit cell
  !
  IMPLICIT NONE
  ! input
  REAL(8) rin(3), a1(3),a2(3),a3(3), ainv(3,3)
  ! output
  REAL(8) rout(3)
  ! local
  REAL(8) x,y,z
  !
  ! bring atomic positions to crystal axis
  !
  x = ainv(1,1)*rin(1)+ainv(1,2)*rin(2)+ainv(1,3)*rin(3)
  y = ainv(2,1)*rin(1)+ainv(2,2)*rin(2)+ainv(2,3)*rin(3)
  z = ainv(3,1)*rin(1)+ainv(3,2)*rin(2)+ainv(3,3)*rin(3)
  !
  ! bring x,y,z in the range between -0.5 and 0.5
  !
  x = x - NINT(x)
  y = y - NINT(y)
  z = z - NINT(z)
  !
  ! bring atomic positions back in cartesian axis
  !
  rout(1) = x*a1(1)+y*a2(1)+z*a3(1)
  rout(2) = x*a1(2)+y*a2(2)+z*a3(2)
  rout(3) = x*a1(3)+y*a2(3)+z*a3(3)
  !
  RETURN
END SUBROUTINE vib_pbc
