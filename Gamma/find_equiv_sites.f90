!
subroutine find_equiv_sites (nat,nax,nsym,irt,has_equivalent,     &
     n_diff_sites,n_equiv_atoms,equiv_atoms)
  !
  implicit none
  integer :: nat, nax, nsym, na, nb, ns, n_diff_sites, irt(48,nat),    &
       equiv_atoms(nax,nat), n_equiv_atoms(nax), has_equivalent(nax)
  !
  n_diff_sites = 0
  do na = 1,nat
     has_equivalent(na) = 0
  end do
  !
  do na = 1,nat
     if (has_equivalent(na).eq.0) then
        n_diff_sites = n_diff_sites + 1
        n_equiv_atoms (n_diff_sites) =  1
        equiv_atoms(n_diff_sites,1) = na
        !
        do nb = na+1,nat
           do ns = 1, nsym
              if ( irt(ns,nb) .eq. na) then
                 has_equivalent(nb) = 1
                 n_equiv_atoms (n_diff_sites) =  &
                      n_equiv_atoms (n_diff_sites) + 1
                 equiv_atoms(n_diff_sites, &
                      n_equiv_atoms(n_diff_sites)) = nb
                 go to 10
              end if
           end do
10         continue
        end do
     end if
  end do
  !
  return
end subroutine find_equiv_sites
