!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine simpson (mesh, func, rab, asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integrator for function stored on the
  !     radial logarithmic mesh
  !
  use kinds, ONLY: DP
  implicit none
  integer :: i, mesh

  real(kind=DP) :: rab (mesh), func (mesh), f1, f2, f3, r12, asum
  !     routine assumes that mesh is an odd number so run check
  !     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
  !       write(*,*) '***error in subroutine radlg'
  !       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
  !       stop
  !     endif
  asum = 0.0d0
  r12 = 1.0d0 / 12.0d0

  f3 = func (1) * rab (1) * r12
  do i = 2, mesh - 1, 2
     f1 = f3
     f2 = func (i) * rab (i) * r12
     f3 = func (i + 1) * rab (i + 1) * r12
     asum = asum + 4.0d0 * f1 + 16.0d0 * f2 + 4.0d0 * f3
  enddo

  return
end subroutine simpson

!=----------------------------------------------------------------------------=!

subroutine simpson_cp90( mesh, func, rab, intg )

  implicit none
  integer mesh
  real(kind=8)  func(mesh), rab(mesh), intg

  real(kind=8) c(4)
  integer I

  if ( mesh .lt. 8 ) call errore ('simpson','few mesh points',8)

  c(1) = 109.0 / 48.d0
  c(2) = -5.d0 / 48.d0
  c(3) = 63.d0 / 48.d0
  c(4) = 49.d0 / 48.d0

  intg = ( func(1)*rab(1) + func(mesh  )*rab(mesh  ) )*c(1) &
       + ( func(2)*rab(2) + func(mesh-1)*rab(mesh-1) )*c(2) &
       + ( func(3)*rab(3) + func(mesh-2)*rab(mesh-2) )*c(3) &
       + ( func(4)*rab(4) + func(mesh-3)*rab(mesh-3) )*c(4)
  do i=5,mesh-4
     intg = intg + func(i)*rab(i)
  end do

  return
end subroutine simpson_cp90

!
!===============================================================
!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
      SUBROUTINE simpson_fpmd(n, func, dx, s)

      USE kinds
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  N
      REAL(dbl), INTENT(IN)  ::  func(N), dx
      REAL(dbl), INTENT(OUT) ::  S

      REAL(dbl) :: C1,C2,C3,C4
      PARAMETER(C1=109.d0/48.d0,C2=-5.d0/48.d0, C3=63.d0/48.d0,C4=49.d0/48.d0)
      INTEGER I

! ... Subroutine body

      S =     func(1)*C1
      S = S + func(2)*C2
      S = S + func(3)*C3
      S = S + func(4)*C4
      DO I = 5, (N-5)
         S = S + func(I)
      END DO
      S = S + func(N-4)*C4
      S = S + func(N-3)*C3
      S = S + func(N-2)*C2
      S = S + func(N-1)*C1

      S = S * dx

      RETURN
      END SUBROUTINE simpson_fpmd
