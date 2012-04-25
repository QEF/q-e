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
      nspin_lsda = 1,     & !  =1 when nspin=1,4 =2 when nspin=2 
      nspin_mag = 1,      & !  =1 when nspin=1,4 (domag=.false.), =2 when
                            !   nspin=2, =4 nspin=4 (domag=.true.)
      nspin_gga = 1,      & !  =1 when nspin=1,4 (domag=.false.)   
                            !  =2 when nspin=2,4 (domag=.true.) (needed with gga)
      i_cons = 0            !  indicator for constrained local quantities
  !
  INTEGER, ALLOCATABLE :: &
  !                         !  when spherical (non-overlapping) integration
      pointlist(:)          !  regions are defined around atoms this index
                            !  say for each point in the fft grid to which 
                            !  atom it is assigned (0 if no atom is selected)
  !
  LOGICAL :: &
      noncolin, &           !  true if noncollinear magnetism is allowed
      lsign=.FALSE.         !  if true use the sign feature to calculate
                            !  rhoup and rhodw
  !
  REAL (DP) :: &
      angle1(ntypx),       &!  Define the polar coordinates of the starting
      angle2(ntypx),       &!  magnetization's direction for each atom
      mcons(3,ntypx)=0.d0, &!  constrained values for local variables
      magtot_nc(3),        &!  total magnetization
      bfield(3)=0.d0,      &!  magnetic field used in some cases
      vtcon,               &!  contribution of the constraining fields to
                            !  the total energy
      r_m(ntypx) = 0.0d0,  &!  Radius for local integrations for each type
      lambda                !  prefactor in the penalty functional 
                            !  for constraints
  !
  REAL (DP), ALLOCATABLE :: &
      factlist(:),         &! weight factors for local integrations
      r_loc(:),            &! local integrated charge 
      m_loc(:,:)            ! local integrated magnetization

  REAL(DP) ::     &
     ux(3)                 ! versor for deciding signs in gga
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_noncol()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( pointlist) )       DEALLOCATE( pointlist )
      IF ( ALLOCATED( factlist ) )       DEALLOCATE( factlist )
      IF ( ALLOCATED( r_loc    ) )       DEALLOCATE( r_loc )
      IF ( ALLOCATED( m_loc    ) )       DEALLOCATE( m_loc )
      !
    END SUBROUTINE deallocate_noncol
    !
END MODULE noncollin_module
