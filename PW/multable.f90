!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine multable (nsym, s, table)
  !-----------------------------------------------------------------------
  !
  !  sets up the multiplication table for a group represented by 3x3
  !  integer matrices and checks that {s} is a group indeed:
  !
  !  table(n,m) = index( s(n)*s(m) )
  !
  use parameters
  implicit none
  !
  !    here the dummy variables
  !
  integer :: nsym, s (3, 3, 48), table (48, 48)
  ! input: the number of symmetry of the
  ! input: the symmetry matrices
  ! output: the multiplication table
  !
  !  and here the local variables
  !
  integer :: irot, jrot, krot, ipol, jpol, kpol, ss (3, 3)
  ! \
  !   counter on rotations
  ! /
  ! \
  !   counters on polarizations
  ! /
  ! buffer multiplication matrix

  logical :: found, smn
  ! if true the table has been set
  ! used to check symmetries
  do irot = 1, nsym
     do jrot = 1, nsym
        !
        do ipol = 1, 3
           ! sets up th
           do jpol = 1, 3
              ! product
              ss (ipol, jpol) = 0
              ! matrix
              do kpol = 1, 3
                 !
                 ss (ipol, jpol) = ss (ipol, jpol) + s (ipol, kpol, jrot) * s ( &
                      kpol, jpol, irot)
                 ! ss=s(j)*s(
                 !
              enddo
              !
           enddo
           !
        enddo
        !
        !     here checks that the input matrices really form a group
        !     and sets the multiplication table
        !
        found = .false.
        do krot = 1, nsym
           smn = .true.
           do ipol = 1, 3
              do jpol = 1, 3
                 smn = smn.and. (s (ipol, jpol, krot) .eq.ss (ipol, jpol) )
              enddo
           enddo
           if (smn) then
              if (found) call error ('Multable', 'Not a group', 1)
              found = .true.
              table (jrot, irot) = krot
           endif
        enddo

     enddo
     if (.not.found) call error ('Multable', ' Not a group', 2)

  enddo
  return
end subroutine multable
