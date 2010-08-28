!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine generate_dynamical_matrix   &
     (nat, nsym, s, invs, irt, at, bg, n_diff_sites, equiv_atoms, &
     has_equivalent, dyn)
  !-----------------------------------------------------------------------
  !
  !  generate the complete dynamical matrix from independent modes only
  !  Input: dyn = irreducible dyn.mat.  Output: dyn = complete dyn.mat.
  !
  USE kinds, only : DP
  implicit none
  integer :: nat, nsym, n_diff_sites, irt(48,nat), invs(48), &
       equiv_atoms(nat,nat), s(3,3,48),  has_equivalent(nat)
  real(DP) :: at(3,3), bg(3,3)
  complex(DP) :: dyn(3*nat,3*nat)
  !
  integer :: isym, na, nb, ni, nj, sni, snj, smu_i, smu_j,  &
       i, j, k, l, mu_k, mu_l
  complex(DP), allocatable :: irreducible_dyn(:,:)
  complex(DP) :: work(3,3)
  logical :: no_equivalent_atoms
  integer, allocatable ::done(:,:)
  !
  no_equivalent_atoms=.true.
  do na = 1,nat
     no_equivalent_atoms = no_equivalent_atoms .and. has_equivalent(na).eq.0
  end do
  if (no_equivalent_atoms) return
  !
  allocate  ( irreducible_dyn( 3*nat, 3*nat))
  call zcopy(3*nat*3*nat,dyn,1,irreducible_dyn,1)
  !
  do na = 1,nat
     if (has_equivalent(na).eq.0 ) then
        do nb = 1,nat
           do i = 1,3
              do j = 1,3
                 work(i,j) = irreducible_dyn(3*(na-1)+i,3*(nb-1)+j)
              end do
           end do
           !
           !  transform to crystal axis
           !
           call trntnsc(work,at,bg,-1)
           do i = 1,3
              do j = 1,3
                 irreducible_dyn(3*(na-1)+i,3*(nb-1)+j) = work(i,j)
              end do
           end do
        end do
     end if
  end do
  !
  allocate  (done( 3*nat, 3*nat))
  do smu_i = 1,3*nat
     do smu_j = 1,3*nat
        dyn(smu_i,smu_j) = (0.d0,0.d0)
        done(smu_i,smu_j)= 0
     end do
  end do
  !
  do isym = 1,nsym
     do na = 1,n_diff_sites
        ni = equiv_atoms(na,1)
        sni = irt(isym,ni)
        do i = 1,3
           smu_i = 3*(sni-1)+i
           do nj = 1,nat
              snj = irt(isym,nj)
              do j = 1,3
                 smu_j = 3*(snj-1)+j
                 if (done(smu_i,smu_j).eq.0) then
                    do k = 1,3
                       mu_k = 3*(ni-1)+k
                       do l = 1,3
                          mu_l = 3*(nj-1)+l
                          dyn(smu_i,smu_j) = dyn(smu_i,smu_j) + &
                               s(i,k,invs(isym)) * s(j,l,invs(isym)) * &
                               irreducible_dyn(mu_k,mu_l)
                          !  rotation matrices are S^-1
                       end do
                    end do
                    done(smu_i,smu_j)=1
                 end if
              end do
           end do
        end do
     end do
  end do
  !
  deallocate(done)
  deallocate(irreducible_dyn)
  !
  do na = 1,nat
     do nb = 1,nat
        do i = 1,3
           do j = 1,3
              work(i,j) = dyn(3*(na-1)+i,3*(nb-1)+j)
           end do
        end do
        !  back to cartesian axes
        call trntnsc(work,at,bg, 1)
        do i = 1,3
           do j = 1,3
              dyn(3*(na-1)+i,3*(nb-1)+j) = work(i,j)
           end do
        end do
     end do
  end do
  !
  return
end subroutine generate_dynamical_matrix
