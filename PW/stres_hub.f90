!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stres_hub ( sigmah )
   !----------------------------------------------------------------------
   !
   ! This routines computes the Hubbard contribution to the internal stress
   ! tensor. It gives in output the array sigmah(i,j) which corresponds to 
   ! the quantity -(1/\Omega)dE_{h}/d\epsilon_{i,j}
   !
#include "machine.h"
   use pwcom
   use io, only : prefix
#ifdef PARA
   use para
#endif
   implicit none
   real (kind=DP) :: sigmah(3,3)        ! output: the Hubbard stresses

   integer :: ipol, jpol, nworddw, nworddb, na, nt, is,isi, m1,m2,m3,m4
   real (kind=DP) :: omin1, current_sum, inverse_sum, sum, temp, flag
   logical :: exst
   real (kind=DP), allocatable :: dns(:,:,:,:)
   !       dns(nat,nspin,5,5), ! the derivative of the atomic occupations
 
   sigmah(:,:) = 0.d0

   allocate (dns(nat,nspin,5,5))
   dns(:,:,:,:) = 0.d0

#ifdef PARA
   if (me.eq.1.and.mypool.eq.1) then
#endif
   call seqopn(iunocc,trim(prefix)//'.occup','formatted',exst)
   read(iunocc,*) ns
   close(unit=iunocc,status='keep')
#ifdef PARA
   end if
#endif

   ! now we open the files containing dwfc and dbeta

   nworddw = 2*npwx*natomwfc
   nworddb = 2*npwx*nkb

   call diropn(23,'dwat',nworddw,exst)
   call diropn(25,'dbeta',nworddb,exst)

#ifdef DEBUG
   do na=1,nat
      do is=1,nspin
         nt = ityp(na)
         if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
            write (*,'(a,2i3)') 'NS(NA,IS) ', na,is
            write (*,'(5f10.4)') ((ns(na,is,m1,m2),m2=1,5),m1=1,5)
         end if
      end do
   end do
#endif
   omin1 = 1.d0/omega
   do ipol = 1,3
      do jpol = 1,ipol
         call dndepsilon(dns,ipol,jpol)
         do na = 1,nat                 
            nt = ityp(na)
            if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
               do is = 1,nspin
#ifdef DEBUG
                  write (*,'(a,4i3)') 'DNS(IPOL,JPOL,NA,IS) ', ipol,jpol,na,is
                  write (*,'(5f10.4)') ((dns(na,is,m1,m2),m2=1,5),m1=1,5)
#endif
                  do m2 = 1,5
                     sigmah(ipol,jpol) = sigmah(ipol,jpol) - omin1 * &
                           Hubbard_U(nt) * 0.5d0 * dns(na,is,m2,m2) 
                     do m1 = 1,5
                        sigmah(ipol,jpol) = sigmah(ipol,jpol) + omin1 * &
                           Hubbard_U(nt) * ns(na,is,m2,m1) * dns(na,is,m1,m2)
                     end do
                  end do
               end do
            end if
         end do
      end do
   end do
   !
   ! Symmetryze the stress tensor
   !
   !      write(6,*) 'ns =', ns(1,2,1,1)
   do ipol = 1,3
      do jpol = ipol,3
         sigmah(ipol,jpol) = sigmah(jpol,ipol)
      end do
   end do

   call trntns(sigmah,at,bg,-1)
   call symtns(sigmah,nsym,s)
   call trntns(sigmah,at,bg,1)

   close(23,status='delete')
   close(25,status='delete')

   deallocate (dns) 

   return
end  subroutine stres_hub
