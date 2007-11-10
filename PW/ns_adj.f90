!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ns_adj
!-----------------------------------------------------------------------
! This routine tries to suggest to the code the right atomic orbital to 
! localize the charge on.
!
   USE kinds,     ONLY : DP
   USE ions_base, ONLY : nat, ntyp => nsp, ityp
   USE ldaU,      ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, starting_ns
   USE scf,       ONLY : rho
   USE lsda_mod,  ONLY : nspin
   USE io_global, ONLY : stdout
 
   implicit none
   !
   integer, parameter:: ldim=7
   integer :: na,nt,is,m1,m2,majs,mins,adjs,mol(ldim),nel,i,j,l,index(ldim) 
   real(DP) :: totoc, delta,lambda(ldim) 
   complex(DP) :: vet(ldim,ldim), f(ldim,ldim), temp
   logical :: adjust
 
   if (ALL(starting_ns == -1.d0)) return
   write (stdout,*) "Modify starting ns matrices according to input values "
 
   if (2*Hubbard_lmax+1>ldim) call errore('ns_adj',' ldim too small',ldim) 

   do na = 1,nat
      nt = ityp(na)
      if (Hubbard_U(nt).ne.0.d0) then
         do is=1,nspin
            do m1 = 1, 2 * Hubbard_l(nt) + 1
               do m2 = 1, 2 * Hubbard_l(nt) + 1
                  f(m1,m2) = rho%ns(m1,m2,is,na)
               end do
            end do
            call cdiagh(2*Hubbard_l(nt)+1, f, ldim, lambda, vet)
            do i = 1, 2 * Hubbard_l(nt) + 1
               if (starting_ns(i,is,nt) >= 0.d0) lambda(i) = starting_ns(i,is,nt)
            end do
            do m1 = 1,2 * Hubbard_l(nt) + 1
               do m2 = m1, 2 * Hubbard_l(nt) + 1
                  temp = 0.d0
                  do i = 1,2 * Hubbard_l(nt) + 1
                     temp = temp + CONJG(vet(m1,i))*lambda(i)*vet(m2,i)     
                  end do
                  rho%ns(m1,m2,is,na) =  DBLE(temp)
                  rho%ns(m2,m1,is,na) = rho%ns(m1,m2,is,na)
               end do
            end do
         end do
      end if
   end do ! on na
   CALL write_ns
   return
end subroutine ns_adj
