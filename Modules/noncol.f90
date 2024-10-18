!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE noncollin_module
  !
  !! Variables for noncollinear magnetism and spin-orbit interactions
  !
  USE kinds, ONLY : DP
  USE parameters, ONLY : ntypx
  !
  INTEGER :: npol
  !! number of coordinates of wfc
  INTEGER :: report
  !! print the local quantities (magnet. and rho) every #report iterations
  INTEGER :: nspin_lsda = 1
  !! equal to 1 when \(\text{nspin}=1,4\) and to 2 when \(\text{nspin}=2\)
  INTEGER :: nspin_mag = 1
  !! equal to 1 when \(\text{nspin}=1,4\) (\(\text{domag}=\text{FALSE}\)), to 2 when 
  !! \(\text{nspin}=2\) and to 4 when \(\text{nspin}=4\) (\(\text{domag}=\text{TRUE}\))
  INTEGER :: nspin_gga = 1
  !! equal to 1 when \(\text{nspin}=1,4\) (\(\text{domag}=\text{FALSE}\)), equal
  !! to 2 when \(\text{nspin}=2,4\) (\(\text{domag}=\text{TRUE}\)) (needed with GGA)
  INTEGER :: i_cons = 0
  !! indicator for constrained local quantities
  !
  INTEGER, ALLOCATABLE :: pointlist(:)
  !! when spherical (non-overlapping) integration regions are defined around
  !! atoms this index says for each point in the fft grid to which atom
  !! it is assigned (0 if no atom is selected).
  !
  LOGICAL :: noncolin
  !! TRUE if noncollinear magnetism is allowed
  LOGICAL :: domag
  !! TRUE if total magnetization is present, FALSE for nonmagnetic calculation
  LOGICAL :: lsign=.FALSE.
  !! if TRUE use the sign feature to calculate rhoup and rhodw
  INTEGER :: colin_mag = -1
  !! equal to 0 if the system does not have a collinear magnetism
  !! equal to -1 if the collinearity is not checked.
  !! larger than 0 if the system has a collinear magnetism (nspin_mag = 2)
  !! equal to 1 if the symmetries with time-reversal is not detected
  !! equal to 2 if the symmetries with time-reversal is detected

  !
  REAL(DP) :: angle1(ntypx)
  !! define the polar coordinates of the starting magnetization
  !! direction for each atom - first angle
  REAL(DP) :: angle2(ntypx)
  !! starting mag. direction - second angle.
  REAL(DP) :: mcons(3,ntypx)=0.d0
  !! constrained values for local variables
  REAL(DP) :: magtot_nc(3)
  !! total magnetization
  REAL(DP) :: bfield(3)=0.d0
  !! magnetic field used in some cases
  REAL(DP) :: vtcon
  !! contribution of the constraining fields to the total energy
  REAL(DP) :: r_m(ntypx)=0.0d0
  !! radius for local integrations for each type
  REAL(DP) :: lambda
  !! prefactor in the penalty functional for constraints
  !
  REAL(DP), ALLOCATABLE :: factlist(:)
  !! weight factors for local integrations
  REAL(DP), ALLOCATABLE :: m_loc(:,:)
  !! local integrated magnetization
  REAL(DP) :: ux(3)
  !! versor for deciding signs in GGA
  !
  !! Variables needed for calculations with spin-orbit
  !
  LOGICAL :: lspinorb
  !! if .TRUE. this calculation uses spin-orbit interactions
  LOGICAL :: lforcet
  !! if .TRUE. apply Force Theorem to calculate MAE
  LOGICAL :: starting_spin_angle
  !! if .TRUE. the initial wavefunctions are spin-angle functions.
  !
  SAVE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_noncol()
      !------------------------------------------------------------------------
      !! Deallocates arrays related to noncollinear magnetism.
      !
      IF ( ALLOCATED( pointlist) )       DEALLOCATE( pointlist )
      IF ( ALLOCATED( factlist ) )       DEALLOCATE( factlist )
      IF ( ALLOCATED( m_loc    ) )       DEALLOCATE( m_loc )
      !
    END SUBROUTINE deallocate_noncol
    !
END MODULE noncollin_module

