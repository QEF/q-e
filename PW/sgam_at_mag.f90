!
! Copyright (C) 2001-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sgam_at_mag ( nrot, s, nat, bg, irt, m_loc, sname, sym, t_rev)
  !-----------------------------------------------------------------------
  !
  !   Find magnetic symmetries, i.e. point-group symmetries that are
  !   also symmetries of the local magnetization - including
  !   rotation + time reversal operations
  !
  USE kinds
  implicit none
  !
  integer, intent(in) :: nrot, nat, s (3, 3, 48), irt (48, nat)
  real(DP), intent(in) :: m_loc(3, nat), bg (3, 3)
  character(len=45), intent(in) :: sname (48)
  !
  ! nrot : order of the parent group (all sym.ops of the lattice)
  ! s    : sym.ops (rotation part only) of the parent group
  ! m_loc: local magnetization, must be invariant under the sym.op.
  ! bg   : basis of the reciprocal-space lattice
  ! sname: point-group symmetry operation name - if it contains an
  !        inversion, the magnetization must be reversed
  !
  logical, intent(inout) :: sym (48)
  !
  ! sym(isym) = .true. if rotation isym is a sym.op. of the crystal
  !                    (i.e. not of the bravais lattice only)
  integer, intent(out) :: t_rev(48)
  !
  ! t_rev(isym) = .true. if sym.op. isym, plus time reversal, is a
  !                      symmetry of the local magnetization
  !
  integer :: na, nb, irot
  logical :: t1, t2
  real(DP) , allocatable ::  mxau(:,:), mrau(:,:)
  ! magnetization and rotated magnetization in crystal axis
  !
  allocate ( mxau(3,nat), mrau(3,nat) )
  !
  !     Compute the local magnetization of each atom in the basis of
  !     the direct lattice vectors
  !
  do na = 1, nat
     mxau (:, na)= bg (1, :) * m_loc (1, na) + &
                   bg (2, :) * m_loc (2, na) + &
                   bg (3, :) * m_loc (3, na)
  enddo
  !
  do irot = 1, nrot
     !
     t_rev(irot) = 0
     !
     if ( sym (irot) ) then
        !
        ! mrau = rotated local magnetization
        !
        do na = 1, nat
            mrau(:,na) = s(1,:,irot) * mxau(1,na) + &
                         s(2,:,irot) * mxau(2,na) + &
                         s(3,:,irot) * mxau(3,na)  
        enddo
        if (sname(irot)(1:3)=='inv') mrau = -mrau
        !
        ! check if this a magnetic symmetry
        !
        t1 = .true.
        t2 = .true.
        do na = 1, nat
           !
           nb = irt (irot,na)
           if ( nb < 1 .or. nb > nat ) call errore ('check_mag_sym', &
               'internal error: out-of-bound atomic index', na) 
           !
           t1 = ( abs(mrau(1,na) - mxau(1,nb)) +       &
                  abs(mrau(2,na) - mxau(2,nb)) +       &
                  abs(mrau(3,na) - mxau(3,nb)) < 1.0D-5 ) .and. t1
           t2 = ( abs(mrau(1,na) + mxau(1,nb))+       &
                  abs(mrau(2,na) + mxau(2,nb))+       &
                  abs(mrau(3,na) + mxau(3,nb)) < 1.0D-5 ) .and. t2
           !
        enddo
        !
        if ( .not.t1 .and. .not.t2 ) then
           ! not a magnetic symmetry
           sym(irot) = .false.
        else if( t2 .and. .not. t1 ) then
           ! magnetic symmetry with time reversal
           t_rev(irot) = 1
        end if
        !
     end if
     !
  enddo
  !
  !   deallocate work space
  !
  deallocate ( mrau, mxau )
  !
  return
END SUBROUTINE sgam_at_mag
