!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_hub(forceh)
   !----------------------------------------------------------------------
   !
   ! This routine computes the Hubbard contribution to the force. It gives
   ! in output the product (dE_{hub}/dn_{ij}^{alpha})(dn_{ij}^{alpha}
   ! /du(alpha,ipol)) which is the force acting on the atom at tau_{alpha}
   ! (in the unit ceel) along the direction ipol.
   !
#include "machine.h"
   USE kinds, ONLY: DP
   USE basis, ONLY: nat, ityp
   USE brilz, ONLY: at, bg
   USE ldaU,  ONLY: hubbard_lmax, hubbard_l, hubbard_u, hubbard_alpha, ns
   USE lsda_mod, ONLY: nspin
   USE symme,    ONLY: s, nsym, irt
   use io_files, only : prefix, iunocc
#ifdef __PARA
   use para
#endif
   implicit none
   real (kind=DP) :: forceh(3,nat)  ! output: the Hubbard forces

   integer :: alpha, na, nt, is, m1, m2, ipol, ldim

   logical ::  exst

   real (kind=DP), allocatable :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat) ! the derivative of the atomic occupations

   ldim= 2 * Hubbard_lmax + 1
   allocate(dns(ldim,ldim,nspin,nat))
   forceh(:,:) = 0.d0
   dns(:,:,:,:) = 0.d0

#ifdef __PARA
   if (me.eq.1.and.mypool.eq.1) then
#endif
      call seqopn (iunocc, trim(prefix)//'.occup', 'formatted', exst)
      read(iunocc,*) ns
      close(unit=iunocc,status='keep')
#ifdef __PARA
   end if
#endif

   do ipol = 1,3
      do alpha = 1,nat                 ! the displaced atom
         call dndtau(dns,ldim,alpha,ipol)
         do na = 1,nat                 ! the Hubbard atom
            nt = ityp(na)
            if (Hubbard_U(nt).ne.0.d0.or. Hubbard_alpha(nt).ne.0.d0) then
               do is = 1,nspin
                  do m2 = 1,ldim
                     forceh(ipol,alpha) = forceh(ipol,alpha) -  &
                           Hubbard_U(nt) * 0.5d0           * dns(m2,m2,is,na)
                     do m1 = 1,ldim
                        forceh(ipol,alpha) = forceh(ipol,alpha) +    &
                           Hubbard_U(nt) * ns(m2,m1,is,na) * dns(m1,m2,is,na)
                     end do
                  end do
               end do
            end if
         end do
      end do
   end do

   deallocate(dns)
   
   if (nspin.eq.1) forceh(:,:) = 2.d0 * forceh(:,:)
   !
   ! The symmetry matrices are in the crystal basis so...
   ! Transform to crystal axis...
   !
   do na=1, nat
      call trnvect(forceh(1,na),at,bg,-1)
   end do
   !
   ! ...symmetrize...
   !
   call symvect(nat,forceh,nsym,s,irt)
   !
   ! ... and transform back to cartesian axis
   !
   do na=1, nat
      call trnvect(forceh(1,na),at,bg, 1)
   end do

   return
end subroutine force_hub
