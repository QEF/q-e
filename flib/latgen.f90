!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
subroutine latgen(ibrav,celldm,a1,a2,a3,omega)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic p (sc)                8  orthorhombic p
  !       2  cubic f (fcc)               9  one face centered orthorhombic
  !       3  cubic i (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal p   11  body centered orthorhombic
  !       5  trigonal r                 12  monoclinic p
  !       6  tetragonal p (st)          13  one face centered monoclinic
  !       7  tetragonal i (bct)         14  triclinic p
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  use kinds, only: DP
  implicit none
  integer ibrav
  real(DP) celldm(6), a1(3), a2(3), a3(3), omega
  !
  real(DP), parameter:: sr2 = 1.414213562373d0, &
                             sr3 = 1.732050807569d0
  integer i,j,k,l,iperm,ir
  real(DP) term, cbya, s, term1, term2, singam, sen
  !
  if(ibrav == 0) go to 100
  !
  !     user-supplied lattice
  !
  do ir=1,3
     a1(ir)=0.d0
     a2(ir)=0.d0
     a3(ir)=0.d0
  end do
  !
  if (ibrav == 1) then
     !
     !     simple cubic lattice
     !
     if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
     !
  else if (ibrav == 2) then
     !
     !     fcc lattice
     !
     if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
     !
  else if (ibrav == 3) then
     !
     !     bcc lattice
     !
     if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     term=celldm(1)/2.d0
     do ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     end do
     a2(1)=-term
     a3(1)=-term
     a3(2)=-term
     !
  else if (ibrav.eq.4) then
     !
     !     hexagonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (3) <= 0.d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
     !
  else if (ibrav == 5) then
     !
     !     trigonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (4) <= - 0.5d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     term1=sqrt(1.d0+2.d0*celldm(4))
     term2=sqrt(1.d0-celldm(4))
     a2(2)=sr2*celldm(1)*term2/sr3
     a2(3)=celldm(1)*term1/sr3
     a1(1)=celldm(1)*term2/sr2
     a1(2)=-a1(1)/sr3
     a1(3)= a2(3)
     a3(1)=-a1(1)
     a3(2)= a1(2)
     a3(3)= a2(3)
     !
  else if (ibrav == 6) then
     !
     !     tetragonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (3) <= 0.d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
     !
  else if (ibrav == 7) then
     !
     !     body centered tetragonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (3) <= 0.d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
     !
  else if (ibrav == 8) then
     !
     !     Simple orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or. &
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
     !
  else if (ibrav == 9) then
     !
     !
     !     One face centered orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or. &
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1) = 0.5 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a3(3) = celldm(1) * celldm(3)
     !
  else if (ibrav == 10) then
     !
     !     All face centered orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.&
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a2(1) = 0.5 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
     !
  else if (ibrav == 11) then
     !
     !     Body centered orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.&
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1) = 0.5 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
     !
  else if (ibrav == 12) then
     !
     !     Simple monoclinic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.   &
         celldm (3) <= 0.d0 .or. abs (celldm (4) ) > 1.d0) &
         call errore ('latgen', 'wrong celldm', ibrav)
     !
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
     !
  else if (ibrav == 13) then
     !
     !     One face centered monoclinic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.   &
         celldm (3) <= 0.d0 .or. abs (celldm (4) ) > 1.d0) &
         call errore ('latgen', 'wrong celldm', ibrav)
     !
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = celldm(1) * celldm(4)
     a1(3) = celldm(1) * sen
     a2(1) = a1(1)
     a2(3) = - a1(3)
     a3(1) = celldm(1) * celldm(2)
     a3(2) = celldm(1) * celldm(3)
     !
  else if (ibrav == 14) then
     !
     !     Triclinic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.  &
         celldm (3) <= 0.d0 .or. abs (celldm (4) ) > 1.d0 .or. &
         abs (celldm (5) ) > 1.d0 .or. abs (celldm (6) ) > 1.d0) &
         call errore ('latgen', 'wrong celldm', ibrav)
     !
     singam=sqrt(1.d0-celldm(6)**2)
     term= sqrt((1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
          -celldm(4)**2-celldm(5)**2-celldm(6)**2)/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
     !
  else
     !
     call errore('latgen',' nonexistent bravais lattice',ibrav)
     !
  end if
  !
  !  calculate unit-cell volume omega
  !
100 omega=0.d0
  s=1.d0
  i=1
  j=2
  k=3
  !
101 do iperm=1,3
     omega=omega+s*a1(i)*a2(j)*a3(k)
     l=i
     i=j
     j=k
     k=l
  end do
!
  i=2
  j=1
  k=3
  s=-s
  if(s < 0.d0) go to 101
  omega=abs(omega)
  return
!
end subroutine latgen
