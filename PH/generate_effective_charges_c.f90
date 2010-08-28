!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine generate_effective_charges_c &
     (nat, nsym, s, invs, irt, at, bg, n_diff_sites, equiv_atoms, &
      has_equivalent, asr, nasr, zv, ityp, ntyp, atm, zstar)
  !-----------------------------------------------------------------------
  !
  ! generate all effective charges
  !
  USE io_global, ONLY : stdout
  USE kinds, only : DP
  USE symme, only : crys_to_cart
  implicit none
  integer :: nat, nsym, n_diff_sites, irt(48,nat), equiv_atoms(nat,nat),&
       s(3,3,48), invs(48), has_equivalent(nat), nasr
  logical :: asr
  integer :: isym, na, ni, sni, i, j, k, l
  integer :: ityp(nat), ntyp
  real(DP) :: zstar(3,3,nat), at(3,3), bg(3,3), sumz, zv(ntyp)
  logical :: done(nat), no_equivalent_atoms
  character(3) :: atm(ntyp)
  !
  no_equivalent_atoms=.true.
  do na = 1,nat
     no_equivalent_atoms = no_equivalent_atoms .and. has_equivalent(na).eq.0
  end do
  if (no_equivalent_atoms) goto 100
  !
  !  zstar in input is in crystal axis
  !
  do na = 1,nat
     if (has_equivalent(na).eq.0 ) then
        done(na)=.true.
     else
        zstar(:,:,na) = 0.d0
        done(na)=.false.
     end if
  end do
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
                            s(i,k,invs(isym))*s(j,l,invs(isym))*zstar(k,l,ni)
                    end do
                 end do
              end do
           end do
           done(sni)=.true.
        end if
     end do
  end do

100 continue
  !
  ! return to Cartesian axis
  !
  do na = 1,nat
     call crys_to_cart ( zstar(:,:,na) )
  end do
  !
  ! add the diagonal part
  !
  do i = 1, 3
     do na = 1, nat
        zstar(i, i, na) = zstar (i, i, na) + zv (ityp (na) )
     enddo
  enddo
  IF (asr.AND.nasr>0) THEN
     DO i=1,3
        DO j=1,3
           sumz=0.0_DP
           DO na=1,nat
              IF (na.ne.nasr) sumz=sumz+zstar(i,j,na)
           ENDDO
           zstar(i,j,nasr)=-sumz
        ENDDO
     ENDDO
  ENDIF
  !
  return
end subroutine generate_effective_charges_c
