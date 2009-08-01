!
! Copyright (C) 2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      SUBROUTINE compute_ux(m_loc,ux,nat)
!
!   This subroutine diagonalizes the spin density matrix and gives as output
!   the spin up and spin down components of the charge
!
   USE kinds, ONLY : dp
   USE constants, ONLY: pi, eps12
   USE io_global,  ONLY :  stdout

   implicit none

   INTEGER, INTENT(IN) :: nat    ! number of atoms
   REAL(DP), INTENT(OUT) :: ux(3)     ! fixed direction to calculate signs
   REAL(DP), INTENT(IN) :: m_loc(3,nat)     ! local moments
         
   REAL(DP) :: amag,uxmod !  
   REAL(DP) :: m(3)       ! magnetization

   INTEGER :: i, na        ! counter on polarizations, atoms
   INTEGER :: it,counter(3)

   ux(:)=m_loc(:,1)
   uxmod=SQRT(ux(1)**2+ux(2)**2+ux(3)**2)
   if (uxmod>eps12) ux=ux/uxmod
   do na=2,nat
      amag=SQRT(m_loc(1,na)**2+m_loc(2,na)**2+m_loc(3,na)**2)
      if (amag > eps12) then
         m(:)=m_loc(:,na)/amag
         if (ABS(ux(1)*m(1)+ux(2)*m(2)+ux(3)*m(3)).lt.1.d-3) then
            ux=ux+0.5d0*m
            uxmod=SQRT(ux(1)**2+ux(2)**2+ux(3)**2)
            ux=ux/uxmod
         endif
      endif
   enddo
   WRITE( stdout,'(/,5x,"Fixed quantization axis for GGA: ", 3f12.6)') &
                          ux(1), ux(2), ux(3)
   return
   end subroutine compute_ux
