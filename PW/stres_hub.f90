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
  USE kinds, ONLY: DP
  USE basis, ONLY: nat, ityp
  USE cell_base, ONLY: omega, at, bg
  USE ldaU,  ONLY: hubbard_lmax, hubbard_l, hubbard_u, hubbard_alpha, ns, &
                   U_projection
  USE lsda_mod, ONLY: nspin
  USE symme,    ONLY: s, nsym
  USE io_files, ONLY : prefix, iunocc
  USE wvfct,    ONLY : gamma_only   
#ifdef DEBUG
  USE io_global,      ONLY : stdout
#endif
#ifdef __PARA
   use para
#endif
   implicit none
   real (kind=DP) :: sigmah(3,3)        ! output: the Hubbard stresses

   integer :: ipol, jpol, na, nt, is,isi, m1,m2,m3,m4
   integer :: ldim
   real (kind=DP) :: omin1, current_sum, inverse_sum, sum, temp, flag
   logical :: exst
   real (kind=DP), allocatable :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat), ! the derivative of the atomic occupations
 
   if (U_projection .ne. "atomic") call errore("stres_hub", &
                   " stress for this U_projection_type not implemented",1)

   if (gamma_only) call errore('stres_hub',&
                   ' LDA+U, stress AND gamma-only not implemented yet',1)

   sigmah(:,:) = 0.d0

   ldim = 2 * Hubbard_lmax + 1
   allocate (dns(ldim,ldim,nspin,nat))
   dns(:,:,:,:) = 0.d0

#ifdef __PARA
   if (me.eq.1.and.mypool.eq.1) then
#endif
   call seqopn(iunocc,trim(prefix)//'.occup','formatted',exst)
   read(iunocc,*) ns
   close(unit=iunocc,status='keep')
#ifdef __PARA
   end if
#endif

#ifdef DEBUG
   do na=1,nat
      do is=1,nspin
         nt = ityp(na)
         if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
            WRITE( stdout,'(a,2i3)') 'NS(NA,IS) ', na,is
            do m1=1,ldim
               WRITE( stdout,'(7f10.4)') (ns(m1,m2,is,na),m2=1,ldim)
            end do
         end if
      end do
   end do
#endif
   omin1 = 1.d0/omega
   do ipol = 1,3
      do jpol = 1,ipol
         call dndepsilon(dns,ldim,ipol,jpol)
         do na = 1,nat                 
            nt = ityp(na)
            if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
               do is = 1,nspin
#ifdef DEBUG
                  WRITE( stdout,'(a,4i3)') 'DNS(IPOL,JPOL,NA,IS) ', ipol,jpol,na,is
                  WRITE( stdout,'(5f10.4)') ((dns(m1,m2,is,na),m2=1,5),m1=1,5)
#endif
                  do m2 = 1, 2 * Hubbard_l(nt) + 1
                     sigmah(ipol,jpol) = sigmah(ipol,jpol) - omin1 * &
                           Hubbard_U(nt) * 0.5d0 * dns(m2,m2,is,na) 
                     do m1 = 1, 2 * Hubbard_l(nt) + 1
                        sigmah(ipol,jpol) = sigmah(ipol,jpol) + omin1 * &
                           Hubbard_U(nt) * ns(m2,m1,is,na) * dns(m1,m2,is,na)
                     end do
                  end do
               end do
            end if
         end do
      end do
   end do
   if (nspin.eq.1) sigmah(:,:) = 2.d0 * sigmah(:,:)

   !
   ! Symmetryze the stress tensor
   !
   do ipol = 1,3
      do jpol = ipol,3
         sigmah(ipol,jpol) = sigmah(jpol,ipol)
      end do
   end do
   
   call trntns(sigmah,at,bg,-1)
   call symtns(sigmah,nsym,s)
   call trntns(sigmah,at,bg,1)

   deallocate (dns)

   return
end  subroutine stres_hub
