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
      baco_ibm_xlf,       & !  variable used to avoid a bug in mpxlf_r compiler
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
  REAL (KIND=DP) :: &
      angle1(ntypx),       &!  Define the polar coordinates of the starting
      angle2(ntypx),       &!  magnetization's direction for each atom
      mcons(3,ntypx),      &!  constrained values for local variables
      vtcon,               &!  contribution of the constraining fields to
                            !  the total energy
      r_m,                 &!  Radius for local integrations
      lambda                !  prefactor in the penalty functional 
                            !  for constraints
  !
  REAL (KIND=DP), ALLOCATABLE :: &
      factlist(:,:),       &! weightenig factors for local integrations
      r_loc(:),            &! local integrated charge 
      m_loc(:,:)            ! local integrated magnetization
     ! mcons(:,:),           ! constrained values for local variables
  !
END MODULE noncollin_module
