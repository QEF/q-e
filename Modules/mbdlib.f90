!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#############################################################
! This module computes the Many-Body Dispersion (MBD)
! van der Waals correction to the system.
! The implementation is based on a portable Fortran library
! by Jan Hermann.
! Written by S. Goger (Luxembourg) and H-Y. Ko (Cornell)
!#############################################################

MODULE libmbd_interface

  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout 
  USE tsvdw_module,     ONLY : veff_pub, vfree_pub, vdw_isolated
  USE ions_base,        ONLY : nat, atm, zv, tau, taui, ntyp => nsp, ityp
  USE cell_base,        ONLY : alat, at, bg, omega  ! at: lattice vectors in real space
  USE funct,            ONLY : get_dft_short
  USE control_flags,    ONLY : conv_elec
  USE constants,        ONLY : ry_kbar
  USE mbd,              ONLY : mbd_input_t, mbd_calc_t
  USE input_parameters, ONLY : tforces, tstress

  IMPLICIT NONE

  PUBLIC:: mbd_interface, init_mbd, clean_mbd
  REAL(dp), PUBLIC:: EmbdvdW  ! MBD correction to the energy
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC:: FmbdvdW  ! Ionic force contribs. (-dE/dr)
  REAL(dp), DIMENSION(3, 3),             PUBLIC:: HmbdvdW  ! Cell derivative contribs. (-dE/da)
  REAL(dp), DIMENSION(3, 3),             PUBLIC:: mbdstress  ! Cell derivative contribs. (-dE/da)

  INTEGER:: na
  TYPE(mbd_input_t):: inp
  TYPE(mbd_calc_t):: calc
  REAL(dp), DIMENSION(:), ALLOCATABLE:: ratios
  REAL(dp), DIMENSION(:,:), ALLOCATABLE:: mbd_gradient
  INTEGER:: code, ierr, my_rank, total
  CHARACTER(200):: origin, msg

  CONTAINS

!#############################################################
! This subroutine sets up the library before the first call
!#############################################################
  SUBROUTINE init_mbd()
  IMPLICIT NONE

  !
  ! Allocation of variables that depend on the number of atoms
  !
  ALLOCATE(inp%atom_types(nat))
  IF(.NOT.ALLOCATED(mbd_gradient)) ALLOCATE(mbd_gradient(3, nat))
  IF(tforces .OR. tstress .AND. .NOT.ALLOCATED(FmbdvdW)) ALLOCATE(FmbdvdW(3, nat))
  ALLOCATE(ratios(nat))


  !
  ! Passing atom types and coordinates for LibMBD
  !
  DO na = 1, nat
    inp%atom_types(na) = trim(atm(ityp(na)))
  ENDDO
  inp%coords = tau*alat  ! HK-TODO: this one works for PW (check if it is for CP)

  !
  ! If we pass lattice vectors to the library, it uses the algorithm for
  ! periodic system automatically
  !
  IF( .NOT.vdw_isolated ) THEN
    inp%lattice_vectors = at*alat  ! Lattice vector in real space
  !                            !HK-TODO: cell vectors at (in PW) and h (in CP, note the transpose relationship)
    inp%k_grid = [1, 1, 1] !the k points grid is not needed to be the same as for the PW calculation, but it would help the convergence (TODO)
  ENDIF

  select case (TRIM(get_dft_short()))  ! An empirical factor needs to be set based on the functiona
  CASE ('PBE')
    inp%xc = 'pbe'
  CASE ('PBE0')
    inp%xc = 'pbe0'
  CASE ('HSE')
    inp%xc = 'hse'
  CASE DEFAULT
    ! Block it off since parametrization is not possible
    CALL errore( 'libmbd_interface', 'current xc functional not yet supported for MBD@rsSCS, use PBE, PBE0 or HSE', 1 )
  END SELECT

  CALL calc%init(inp)
  CALL calc%get_exception(code, origin, msg)
  IF (code > 0) THEN
    WRITE( stdout, * ) msg
    CALL errore( 'libmbd_interface', 'Many-Body Dispersion call crashed. This is most likely due to a numerical &
&  error, please check your system carefully.', 1 )
    STOP
  ENDIF


  END SUBROUTINE init_mbd

!#############################################################
! This subroutine calculates the energy and (if needed) forces and stress
!#############################################################
    
  SUBROUTINE mbd_interface()

  IF (.NOT.conv_elec) RETURN ! Wavefunction derivatives are still in progress,
!for now we only can add correction for converged wavefunction

  !
  ! Passing the current parameters to the library
  CALL calc%update_coords(tau*alat)
  DO na = 1, nat
    ratios(na)=veff_pub(na)/vfree_pub(ityp(na))
  ENDDO


  CALL calc%update_vdw_params_from_ratios(ratios)
  IF( .NOT.vdw_isolated ) THEN
    CALL calc%update_lattice_vectors(at*alat)
  ENDIF

  CALL calc%evaluate_vdw_method(EmbdvdW)  !MBD energy
  IF (tforces .OR. tstress) THEN
    CALL calc%get_gradients(mbd_gradient)  
    FmbdvdW = -mbd_gradient  ! Ionic forces with correct sign
  ENDIF
  !

  IF(tforces .OR. tstress .AND. .NOT.vdw_isolated ) THEN
    CALL calc%get_lattice_stress(mbdstress)
    HmbdvdW=-2.0d0*mbdstress/omega ! conversion to rydberg
  ENDIF

  RETURN
  END SUBROUTINE mbd_interface

!#############################################################
! Subroutine to de-allocate internal variables
!#############################################################

  SUBROUTINE clean_mbd()
  IMPLICIT NONE

  CALL calc%destroy()
  IF(ALLOCATED(inp%atom_types)) DEALLOCATE(inp%atom_types)
  IF(ALLOCATED(ratios)) DEALLOCATE(ratios)
  IF(ALLOCATED(mbd_gradient)) DEALLOCATE(mbd_gradient)
  IF(ALLOCATED(veff_pub)) DEALLOCATE(veff_pub)
  IF(ALLOCATED(vfree_pub)) DEALLOCATE(vfree_pub)

  END SUBROUTINE clean_mbd
END MODULE libmbd_interface
