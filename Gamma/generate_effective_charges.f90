!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine generate_effective_charges &
     (nat,nsym,s,irt,at,bg,n_diff_sites,equiv_atoms,has_equivalent,zstar)
  !-----------------------------------------------------------------------
  !
  ! generate all effective charges
  !
#include "machine.h"
  use parameters, only : DP
  implicit none
  integer :: nat, nsym, n_diff_sites, irt(48,nat), equiv_atoms(nat,nat),&
       s(3,3,48), has_equivalent(nat)
  integer :: isym, na, ni, nj, sni, i, j, k, l
  integer :: table(48,48), invs(3,3,48)
  real(kind=DP) :: zstar(3,3,nat), at(3,3), bg(3,3)
  logical :: done(nat), no_equivalent_atoms
  !
  no_equivalent_atoms=.true.
  do na = 1,nat
     no_equivalent_atoms = no_equivalent_atoms .and. has_equivalent(na).eq.0
  end do
  if (no_equivalent_atoms) return
  ! transform to cartesian axis
  do na = 1,nat
     if (has_equivalent(na).eq.0 ) then
        call trntns(zstar(1,1,na),at,bg,-1)
        done(na)=.true.
     else
        zstar(:,:,na) = 0.d0
        done(na)=.false.
     end if
  end do
  !
  ! recalculate S^-1 (once again)
  !
  call multable (nsym,s,table)
  call inverse_s(nsym,s,table,invs)
  !
  do isym = 1,nsym
     do na = 1,n_diff_sites
        ni = equiv_atoms(na,1)
        sni = irt(isym,ni)
        if ( .not.done(sni) ) then
           do i = 1,3
              do j = 1,3
                 do k = 1,3
                    do l = 1,3
                       zstar(i,j,sni) =  zstar(i,j,sni) +  &
                            invs(i,k,isym)*invs(j,l,isym)*zstar(k,l,ni)
                    end do
                 end do
              end do
           end do
           done(sni)=.true.
        end if
     end do
  end do
  ! ritorna ad assi cartesiani
  do na = 1,nat
     call trntns(zstar(1,1,na),at,bg, 1)
  end do
  !
  return
end subroutine generate_effective_charges
