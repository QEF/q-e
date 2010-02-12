!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symmetrize_at(nsym, s, nat, tau, ityp, at, bg, &
           nr1, nr2, nr3, irt, ftau, alat, omega)
  !-----------------------------------------------------------------------
  !
  !     given a point group, this routine finds the subgroup which is
  !     the point group of the crystal under consideration
  !     non symmorphic groups non allowed, provided that fractional
  !     translations are commensurate with the FFT grid
  !
  !     It sets the array sym, which for each operation of the original
  !     point group is true if this operation is also an operation of the
  !     total point group
  !
  USE io_global,  ONLY : stdout
  USE cellmd, ONLY: at_old, lmovecell
  USE symm_base, ONLY: invs
  USE kinds
  implicit none
  !
  !     input variables
  !
  integer :: nsym, s (3, 3, 48), nat, ityp (nat), nr1, nr2, nr3
  real(DP) :: tau (3, nat), at (3, 3), bg (3, 3), alat, omega
  ! nsym : number of symmetry operation of the crystal
  ! s    : symmetry operations of parent group
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  ! nr*  : dimensions of the FFT mesh
  ! tau  : cartesian coordinates of the atoms
  ! at   : basis of the real-space lattice
  ! bg   :  "   "   "  reciprocal-space lattice
  !
  !     output variables
  !
  integer :: irt (48, nat), ftau (3, 48)
  ! irt(isym,na) : sym.op. isym sends atom na into atom irt(isym,na)
  ! ftau(:,isym) : fractional translation associated to sym.op. isym
  !                (in FFT coordinates: crystal axis, multiplied by nr*)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  !    local variables
  !
  integer :: na, icar, ipol, jpol, kpol, lpol, irot
  ! counters
  real(DP) , allocatable :: xau (:,:)
  ! atomic coordinates in crystal axis
  real(DP) :: work, obnr(3), bg_old(3,3), sat(3,3), wrk(3,3), ba(3,3)
  !
  allocate(xau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !

  xau = tau
  tau = 0.d0

  call cryst_to_cart( nat, xau, bg, -1 )

  !
  obnr(1) = 1.d0/DBLE(nr1)
  obnr(2) = 1.d0/DBLE(nr2)
  obnr(3) = 1.d0/DBLE(nr3)
  do irot = 1, nsym
     do na = 1, nat
        do kpol = 1, 3
           work =  s (1, kpol, irot) * xau (1, na) + &
                   s (2, kpol, irot) * xau (2, na) + &
                   s (3, kpol, irot) * xau (3, na) - &
                   ftau(kpol,irot)* obnr(kpol) 
           tau (kpol, irt(irot,na)) = tau (kpol, irt(irot,na)) + work &
                                    - nint(work-xau(kpol,irt(irot,na))) 
        enddo
     enddo
  enddo
  tau (:,:) = tau(:,:)/nsym
  !
  !  If the cell is moving then the lattice vectors has to be
  !  symmetrized as well
  !
  if (lmovecell) then
     CALL recips( at_old(1,1), at_old(1,2), at_old(1,3), &
                  bg_old(1,1), bg_old(1,2), bg_old(1,3) )
     do ipol=1,3
        do jpol=1,3
           ba(ipol,jpol) = bg_old(1,ipol) * at(1,jpol) + &
                           bg_old(2,ipol) * at(2,jpol) + &
                           bg_old(3,ipol) * at(3,jpol) 
        end do
     end do
     at = 0.d0
     !
     !  at(i) = 1/nsym sum_S  at_old(m) S(l,m) <bg_old(l)|at(k)> invS(i,k) 
     !
     do irot=1,nsym
        do icar = 1, 3
           do lpol = 1, 3
              sat(icar,lpol) = at_old(icar,1) * s(lpol,1,irot) &
                             + at_old(icar,2) * s(lpol,2,irot) &
                             + at_old(icar,3) * s(lpol,3,irot) 
           end do
        end do
        do icar = 1, 3
           do kpol =1, 3
              wrk(icar,kpol) = sat(icar,1) * ba(1,kpol) &
                             + sat(icar,2) * ba(2,kpol) &
                             + sat(icar,3) * ba(3,kpol)
           end do
        end do
        do icar = 1, 3
           do ipol =1, 3
              at(icar,ipol) = at(icar,ipol)  &
                            + wrk(icar,1) * s(ipol,1,invs(irot)) &
                            + wrk(icar,2) * s(ipol,2,invs(irot)) &
                            + wrk(icar,3) * s(ipol,3,invs(irot))
           end do
        end do
     end do
     at(:,:) = at(:,:) / nsym
     CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  end if
  !
  !   deallocate work space
  !
  deallocate (xau)
  !
  call cryst_to_cart(nat, tau, at, 1)

  write (stdout,*) " SYMMETRIZED ATOMIC COORDINATES "

  call output_tau(lmovecell, .FALSE.)
  !
  return
end subroutine symmetrize_at

