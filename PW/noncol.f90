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
  USE kinds, ONLY : DP
  USE parameters, ONLY : ntypx
  !
  SAVE
  !
  INTEGER :: &
      npol,               & !  number of coordinates of wfc
      report,             & !  print the local quantities (magnet. and rho)
                            !  every #report iterations
      i_cons                !  indicator for contrained local quantities
  !
  INTEGER, ALLOCATABLE :: &
      pointlist(:,:),     & !  points in the integration volume
                            !  around atom na
      pointnum(:)           !  number of such points
  !
  LOGICAL :: &
      noncolin              !  true if noncollinear magnetism is allowed
  !
  REAL (DP) :: &
      angle1(ntypx),       &!  Define the polar coordinates of the starting
      angle2(ntypx),       &!  magnetization's direction for each atom
      mcons(3,ntypx),      &!  constrained values for local variables
      magtot_nc(3),        &!  total magnetization
      bfield(3),           &!  magnetic field used in some cases
      vtcon,               &!  contribution of the constraining fields to
                            !  the total energy
      r_m = 0.0d0,         &!  Radius for local integrations
      lambda                !  prefactor in the penalty functional 
                            !  for constraints
  !
  REAL (DP), ALLOCATABLE :: &
      factlist(:,:),       &! weightenig factors for local integrations
      r_loc(:),            &! local integrated charge 
      m_loc(:,:)            ! local integrated magnetization
     ! mcons(:,:),           ! constrained values for local variables
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_noncol()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( pointlist) )       DEALLOCATE( pointlist )
      IF ( ALLOCATED( factlist ) )       DEALLOCATE( factlist )
      IF ( ALLOCATED( pointnum ) )       DEALLOCATE( pointnum )
      IF ( ALLOCATED( r_loc    ) )       DEALLOCATE( r_loc )
      IF ( ALLOCATED( m_loc    ) )       DEALLOCATE( m_loc )
      !
    END SUBROUTINE deallocate_noncol
    !
END MODULE noncollin_module
