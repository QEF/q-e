!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symv (vect, nsym, s, sname, t_rev, at, bg)
   !--------------------------------------------------------------------
   !
   ! This routine symmetrizes a vector keeping only the component that
   ! remains invariant under the symmetry operations of the crystal.
   ! On input and on output vect is in cartesian coordinates.
   ! The vector is supposed to be axial: inversion does not change it. 
   ! Time reversal changes its sign. Note that only groups compatible with 
   ! a finite magnetization give an nonzero output vector. 
   ! 
   !
   USE kinds
   implicit none
   !
   !    I/O variables first
   !
   integer, intent(in) ::      & ! 
              nsym,            & ! input: the number of symmetries
              t_rev(48),       & ! input: the time reversal informations
              s (3, 3, 48)       ! input: the rotation matrices
   real(DP), intent(in) :: at(3,3), bg(3,3) 
   real(DP), intent(inout) :: vect(3)  ! inp/out: the vector to rotate
   CHARACTER(LEN=45), INTENT(IN) ::  sname(48)   ! name of the symmetries
   !
   !   the local variables
   !
   integer :: isym               ! counter on symmetries

   real(DP) :: work (3), segno

   if (nsym.eq.1) return
!
!  The vector is transformed in crystal axis
!
   work(:) = vect(1)*at(1,:) + vect(2)*at(2,:) + vect(3)*at(3,:)
   vect = work
!
!  It is symmetrized
!
   work = 0.d0
   do isym = 1, nsym
      segno=1.0_DP
      IF (sname(isym)(1:3)=='inv') segno=-1.0_DP
      IF (t_rev(isym)==1) segno=-1.0_DP*segno
      work (:) = work (:) + segno * &
                       s (:, 1, isym) * vect (1) + &
                       s (:, 2, isym) * vect (2) + &
                       s (:, 3, isym) * vect (3)
   enddo
   work=work/nsym
!
!  And back in cartesian coordinates.
!
   vect(:) = work(1) * bg(:,1) + work(2) * bg(:,2) + work(3) * bg(:,3)

   return
end subroutine symv
