!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
      subroutine latgen(ibrav,celldm,a1,a2,a3,omega)
!-----------------------------------------------------------------------
!     sets up the crystallographic vectors a1,a2, and a3.
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
!
!     NOTA BENE: the cases ibrav=5, 7, 10 have been modified so as to yield
!     right-handed axis triplets. Boxes for US PPs do not seem to work with
!     the original (left-handed) axis triplets. PG 27-Jan-2001
!
      implicit none
      integer ibrav
      real(kind=8) celldm(6),a1(3),a2(3),a3(3)
!
      integer i,j,k,l,iperm,ir
      real(kind=8) sr2,sr3,term,cbya,s,term1,omega,term2,singam,sen
!
      if(ibrav.eq.0) go to 100
!
      sr2=sqrt(2.d0)
      sr3=sqrt(3.d0)
!
      do ir=1,3
         a1(ir)=0.d0
         a2(ir)=0.d0
         a3(ir)=0.d0
      end do
!
      if (ibrav.eq.1) then
!
         a1(1)=celldm(1)
         a2(2)=celldm(1)
         a3(3)=celldm(1)
!
      else if (ibrav.eq.2) then
!
         term=celldm(1)/2.d0
         a1(1)=-term
         a1(3)=term
         a2(2)=term
         a2(3)=term
         a3(1)=-term
         a3(2)=term
!
      else if (ibrav.eq.3) then
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
         cbya=celldm(3)
         a1(1)=celldm(1)
         a2(1)=-celldm(1)/2.d0
         a2(2)=celldm(1)*sr3/2.d0
         a3(3)=celldm(1)*cbya
!
      else if (ibrav.eq.5) then
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
      else if (ibrav.eq.6) then
!
         cbya=celldm(3)
         a1(1)=celldm(1)
         a2(2)=celldm(1)
         a3(3)=celldm(1)*cbya
!
      else if (ibrav.eq.7) then
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
      else if (ibrav.eq.8) then
!
         a1(1)=celldm(1)
         a2(2)=celldm(1)*celldm(2)
         a3(3)=celldm(1)*celldm(3)
!
      else if (ibrav.eq.9) then
!
         a1(1) = 0.5 * celldm(1)
         a1(2) = a1(1) * celldm(2)
         a2(1) = - a1(1)
         a2(2) = a1(2)
         a3(3) = celldm(1) * celldm(3)
!
      else if (ibrav.eq.10) then
!
         a2(1) = 0.5 * celldm(1)
         a2(2) = a2(1) * celldm(2)
         a1(1) = a2(1)
         a1(3) = a2(1) * celldm(3)
         a3(2) = a2(1) * celldm(2)
         a3(3) = a1(3)
!
      else if (ibrav.eq.11) then
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
      else if (ibrav.eq.12) then
!
         sen=sqrt(1.d0-celldm(4)**2)
         a1(1)=celldm(1)
         a2(1)=celldm(1)*celldm(2)*celldm(4)
         a2(2)=celldm(1)*celldm(2)*sen
         a3(3)=celldm(1)*celldm(3)
!
      else if (ibrav.eq.13) then
!
         sen = sqrt( 1.d0 - celldm(4) ** 2 )
         a1(1) = celldm(1) * celldm(4)
         a1(3) = celldm(1) * sen
         a2(1) = a1(1)
         a2(3) = - a1(3)
         a3(1) = celldm(1) * celldm(2)
         a3(2) = celldm(1) * celldm(3)
!
      else if (ibrav.eq.14) then
!
         singam=sqrt(1.d0-celldm(6)**2)
         term=sqrt((1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
     &     -celldm(4)**2-celldm(5)**2-celldm(6)**2)/(1.d0-celldm(6)**2))
         a1(1)=celldm(1)
         a2(1)=celldm(1)*celldm(2)*celldm(6)
         a2(2)=celldm(1)*celldm(2)*singam
         a3(1)=celldm(1)*celldm(3)*celldm(5)
         a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
         a3(3)=celldm(1)*celldm(3)*term
!
      else
!
         call errore('latgen',' bravais lattice not programmed',ibrav)
!
      end if
!
  100 omega=0.d0
      s=1.d0
      i=1
      j=2
      k=3
!
 101  do iperm=1,3
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
      if(s.lt.0.d0) go to 101
      omega=abs(omega)
      return
!
      end subroutine

!---------------------------------------------------------------------
subroutine latgen_pw (ibrav, celldm, a1, a2, a3)
  !---------------------------------------------------------------------
  !
  !   This routine sets up the crystallographic vectors a1, a2, and a3.
  !   The a's are expressed in units of celldm(1) ( a0 for
  !   cubic lattices). This version contains all 14 Bravais lattices.
  !
  use parameters
  implicit none
  !
  !     First the input variables
  !
  real(kind=DP) :: celldm (6), a1 (3), a2 (3), a3 (3)
  ! input: the dimensions of the lattice
  ! output: the first lattice vector
  ! output: the second lattice vector
  ! output: the third lattice vector
  integer :: ibrav
  ! input: the index of the Bravais lattice
  !
  !    Here the local variables required by the routine
  !
  real(kind=DP) :: sr2, sr3
  ! the square root of 2.0
  ! the square root of 3.0
  parameter (sr2 = 1.414213562373d0, sr3 = 1.732050807569d0)
  real(kind=DP) :: term, term1, term2, cbya, singam, sin
  !\
  ! \
  !   Auxiliary variables
  ! /
  !/
  ! the sine function

  integer :: ipol
  ! counter on the coordinates
  if (ibrav.eq.0) then
     !
     !     user supplied lattice
     !
     if (celldm (1) .le.0.d0) call errore ('latgen', 'wrong celldm', &
          1)
  else
     do ipol = 1, 3
        a1 (ipol) = 0.d0
        a2 (ipol) = 0.d0
        a3 (ipol) = 0.d0
     enddo

  endif
  if (ibrav.eq.1) then
     !
     !     simple cubic lattice
     !
     if (celldm (1) .le.0.d0) call errore ('latgen', 'wrong celldm', &
          1)
     a1 (1) = 1.0d0
     a2 (2) = 1.0d0
     a3 (3) = 1.0d0
  elseif (ibrav.eq.2) then
     !
     !     fcc lattice
     !
     if (celldm (1) .le.0.d0) call errore ('latgen', 'wrong celldm', &
          2)
     term = 1.0d0 / 2.d0
     a1 (1) = - term
     a1 (3) = term
     a2 (2) = term
     a2 (3) = term
     a3 (1) = - term
     a3 (2) = term
  elseif (ibrav.eq.3) then
     !
     !     bcc lattice
     !
     if (celldm (1) .le.0.d0) call errore ('latgen', 'wrong celldm', &
          3)
     term = 1.0d0 / 2.d0
     do 20 ipol = 1, 3
        a1 (ipol) = term
        a2 (ipol) = term
        a3 (ipol) = term
20   end do
     a2 (1) = - term
     a3 (1) = - term
     a3 (2) = - term
  elseif (ibrav.eq.4) then
     !
     !     hexagonal lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (3) .le.0.d0) call errore ( &
          'latgen', 'wrong celldm', 4)
     cbya = celldm (3)
     a1 (1) = 1.0d0
     a2 (1) = - 1.0d0 / 2.d0
     a2 (2) = sr3 / 2.d0
     a3 (3) = cbya
  elseif (ibrav.eq.5) then
     !
     !     trigonal lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (4) .le. - 0.5d0) call errore &
          ('latgen', 'wrong celldm', 5)
     term1 = sqrt (1.d0 + 2.d0 * celldm (4) )
     term2 = sqrt (1.d0 - celldm (4) )
     a1 (2) = sr2 * term2 / sr3
     a1 (3) = term1 / sr3
     a2 (1) = term2 / sr2
     a2 (2) = - a2 (1) / sr3
     a2 (3) = a1 (3)
     a3 (1) = - a2 (1)
     a3 (2) = a2 (2)
     a3 (3) = a1 (3)
  elseif (ibrav.eq.6) then
     !
     !     tetragonal lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (3) .le.0.d0) call errore ( &
          'latgen', 'wrong celldm', 6)
     cbya = celldm (3)
     a1 (1) = 1.0d0
     a2 (2) = 1.0d0
     a3 (3) = cbya
  elseif (ibrav.eq.7) then
     !
     !     body centered tetragonal lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (3) .le.0.d0) call errore ( &
          'latgen', 'wrong celldm', 7)
     cbya = celldm (3)
     a1 (1) = 1.0d0 / 2.d0
     a1 (2) = a1 (1)
     a1 (3) = cbya / 2.d0
     a2 (1) = a1 (1)
     a2 (2) = - a1 (1)
     a2 (3) = a1 (3)
     a3 (1) = - a1 (1)
     a3 (2) = - a1 (1)
     a3 (3) = a1 (3)
  elseif (ibrav.eq.8) then
     !
     !     Simple orthorhombic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0) call errore ('latgen', 'wrong celldm', 8)
     a1 (1) = 1.0d0
     a2 (2) = celldm (2)
     a3 (3) = celldm (3)
  elseif (ibrav.eq.9) then
     !
     !     One face centered orthorhombic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0) call errore ('latgen', 'wrong celldm', 9)
     a1 (1) = 0.5d0
     a1 (2) = 0.5d0 * celldm (2)
     a2 (1) = - a1 (1)
     a2 (2) = a1 (2)

     a3 (3) = celldm (3)
  elseif (ibrav.eq.10) then
     !
     !     All face centered orthorhombic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0) call errore ('latgen', 'wrong celldm', 10)
     a1 (1) = 0.5d0
     a1 (2) = 0.5d0 * celldm (2)
     a2 (1) = a1 (1)
     a2 (3) = 0.5d0 * celldm (3)
     a3 (2) = 0.5d0 * celldm (2)
     a3 (3) = a2 (3)
  elseif (ibrav.eq.11) then
     !
     !     Body centered orthorhombic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0) call errore ('latgen', 'wrong celldm', 11)
     a1 (1) = 0.5d0
     a1 (2) = 0.5d0 * celldm (2)
     a1 (3) = 0.5d0 * celldm (3)
     a2 (1) = - a1 (1)
     a2 (2) = a1 (2)
     a2 (3) = a1 (3)
     a3 (1) = - a1 (1)
     a3 (2) = - a1 (2)
     a3 (3) = a1 (3)
  elseif (ibrav.eq.12) then
     !
     !     Simple monoclinic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0.or.abs (celldm (4) ) .gt.1.d0) call errore ('latgen', &
          'wrong celldm', 12)
     sin = sqrt (1.d0 - celldm (4) **2)
     a1 (1) = 1.0d0
     a2 (1) = celldm (2) * celldm (4)
     a2 (2) = celldm (2) * sin
     a3 (3) = celldm (3)
  elseif (ibrav.eq.13) then
     !
     !     One face centered monoclinic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0.or.abs (celldm (4) ) .gt.1.d0) call errore ('latgen', &
          'wrong celldm', 13)
     sin = sqrt (1.d0 - celldm (4) **2)
     a1 (1) = celldm (4)
     a1 (3) = sin
     a2 (1) = a1 (1)
     a2 (3) = - a1 (3)
     a3 (1) = celldm (2)
     a3 (2) = celldm (3)
  elseif (ibrav.eq.14) then
     !
     !     Triclinic lattice
     !
     if (celldm (1) .le.0.d0.or.celldm (2) .le.0.d0.or.celldm (3) &
          .le.0.d0.or.abs (celldm (4) ) .gt.1.d0.or.abs (celldm (5) ) &
          .gt.1.d0.or.abs (celldm (6) ) .gt.1.d0) call errore ('latgen', &
          'wrong celldm', 14)
     sin = sqrt (1.d0 - celldm (4) **2)
     singam = sqrt (1.d0 - celldm (6) **2)
     term = sqrt ( (1.d0 + 2.d0 * celldm (4) * celldm (5) * celldm ( &
          6) - celldm (4) **2 - celldm (5) **2 - celldm (6) **2) / &
          (1.d0 - celldm (6) **2) )
     a1 (1) = 1.0d0
     a2 (1) = celldm (2) * celldm (6)
     a2 (2) = celldm (2) * singam
     a3 (1) = celldm (3) * celldm (5)
     a3 (2) = celldm (3) * (celldm (4) - celldm (5) * celldm (6) ) &
          / singam
     a3 (3) = celldm (3) * term
  else
     call errore ('latgen', ' wrong ibrav ', ibrav)

  endif
  return
end subroutine latgen_pw
