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
  USE ions_base,        ONLY : nat, atm, tau, ityp
  USE cell_base,        ONLY : alat, at, ainv
  USE funct,            ONLY : get_dft_short
  USE control_flags,    ONLY : conv_elec
  USE constants,        ONLY : ry_kbar
#if !defined(__NOMBD)
  USE mbd,              ONLY : mbd_input_t, mbd_calc_t
#endif
  IMPLICIT NONE

  PUBLIC:: mbd_interface, init_mbd, clean_mbd
  REAL(dp), PUBLIC:: EmbdvdW  ! MBD correction to the energy
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC:: FmbdvdW  ! Ionic force contribs. (-dE/dr)
  REAL(dp), DIMENSION(3, 3),             PUBLIC:: HmbdvdW  ! Cell derivative contribs. (-dE/da)
  REAL(dp), DIMENSION(3, 3)                    :: cell_derivs


  INTEGER:: na
  LOGICAL:: do_gradients
#if !defined(__NOMBD)
  TYPE(mbd_input_t):: inp
  TYPE(mbd_calc_t):: calc
#endif
  REAL(dp), DIMENSION(:), ALLOCATABLE:: ratios
  REAL(dp), DIMENSION(:,:), ALLOCATABLE:: mbd_gradient
  INTEGER:: code, ierr, my_rank, total
  CHARACTER(200):: origin, msg

  CONTAINS

!#############################################################
! This subroutine sets up the library before the first call
!#############################################################
SUBROUTINE init_mbd ( nks_start, nk1, nk2, nk3, k1, k2, k3, tprnfor, tstress )
  !
  INTEGER, INTENT(IN) :: nks_start, nk1, nk2, nk3, k1, k2, k3
  LOGICAL, INTENT(IN) :: tprnfor, tstress
  !
  ! Allocation of variables that depend on the number of atoms
  !
#if defined(__NOMBD)
    CALL errore( 'libmbd_interface', 'Many-Body Dispersion not compiled',1)
#else
  ALLOCATE(inp%atom_types(nat))
  !
  EmbdvdW  = 0.0_dp
  do_gradients = tprnfor .OR. tstress
  IF ( do_gradients ) THEN
     !
     IF(.NOT.ALLOCATED(mbd_gradient)) ALLOCATE(mbd_gradient(3, nat))
     !
     IF(.NOT.ALLOCATED(FmbdvdW)) ALLOCATE(FmbdvdW(3, nat))
     !
  END IF
  !
  ALLOCATE(ratios(nat))

  inp%log_level=1
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
    !
    IF ( nks_start == 0 ) THEN
      ! K-point mesh
      inp%k_grid = [nk1, nk2, nk3]
      inp%k_grid_shift = 0.5_DP
      !
      IF (k1 .EQ. 0 .AND. k2 .EQ. 0 .AND. k3 .EQ. 0) &
        CALL infomsg('mbdlib','k-point shift ignored')
      !
    ELSE
      inp%k_grid = [1, 1, 1] !set default k points grid
      inp%k_grid_shift = 0.5_DP ! set default shift
    ENDIF
    !
  ENDIF
  !
  WRITE(stdout, '(5x,"mbdlib: K-point grid set to ",3I3,", shift: ",F4.2)') &
          inp%k_grid, inp%k_grid_shift
  !
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
#endif

  END SUBROUTINE init_mbd

!#############################################################
! This subroutine calculates the energy and (if needed) forces and stress
!#############################################################
    
  SUBROUTINE mbd_interface()

#if !defined(__NOMBD)
  IF (.NOT.conv_elec) RETURN ! Wavefunction derivatives are still in progress,
!for now we only can add correction for converged wavefunction
 CALL infomsg('mbdlib','MBD wavefunction derivatives not yet supported. '//&
  & 'Performing non-self-consistent MBD calculation upon SCF convergence.') 

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
  IF ( do_gradients ) THEN
    CALL calc%get_gradients(mbd_gradient)
    FmbdvdW = -mbd_gradient  ! Ionic forces with correct sign
  ENDIF
  !

  IF( do_gradients .AND. .NOT.vdw_isolated ) THEN
    CALL calc%get_lattice_stress(cell_derivs)
    HmbdvdW=MATMUL(cell_derivs, TRANSPOSE(ainv))
  ENDIF

  RETURN
#endif
  END SUBROUTINE mbd_interface

!#############################################################
! Subroutine to de-allocate internal variables
!#############################################################

  SUBROUTINE clean_mbd()
  IMPLICIT NONE

#if !defined(__NOMBD)
  CALL calc%destroy()
  IF(ALLOCATED(inp%atom_types)) DEALLOCATE(inp%atom_types)
  IF(ALLOCATED(ratios)) DEALLOCATE(ratios)
  IF(ALLOCATED(mbd_gradient)) DEALLOCATE(mbd_gradient)
  IF(ALLOCATED(veff_pub)) DEALLOCATE(veff_pub)
  IF(ALLOCATED(vfree_pub)) DEALLOCATE(vfree_pub)
#endif

  END SUBROUTINE clean_mbd
END MODULE libmbd_interface
