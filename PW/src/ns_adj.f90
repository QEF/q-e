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
   USE noncollin_module, ONLY : noncolin, npol
   USE io_global, ONLY : stdout
 
   implicit none
   !
   integer, parameter:: ldmx=7
   integer :: na, nt, is, m1, m2, ldim, i, j, l 
   real(DP) :: lambda(npol*ldmx) 
   complex(DP) :: vet(npol*ldmx,npol*ldmx), f(npol*ldmx,npol*ldmx), temp
 
   if (ALL(starting_ns == -1.0_dp)) return
   write (stdout,*) "Modify starting ns matrices according to input values "
 
   if (2*Hubbard_lmax+1>ldmx) call errore('ns_adj',' ldmx too small',ldmx) 

   do na = 1, nat
      nt = ityp(na)
      if (Hubbard_U(nt).ne.0.d0) then
       ldim = 2 * Hubbard_l(nt) + 1 

       if (noncolin) then

         do m1 = 1, ldim
            do m2 = 1, ldim
              f(m1, m2)           = rho%ns_nc(m1, m2, 1, na)
              f(m1, ldim+m2)      = rho%ns_nc(m1, m2, 2, na)
              f(ldim+m1, m2)      = rho%ns_nc(m1, m2, 3, na)
              f(ldim+m1, ldim+m2) = rho%ns_nc(m1, m2, 4, na)
            end do
         end do
         call cdiagh( npol*ldim, f, npol*ldmx, lambda, vet)
         j = 0 
         do is = 1, npol
           do i = 1, ldim
             j = j + 1
             if (starting_ns(i,is,nt) >= 0.d0) lambda(j) = starting_ns(i,is,nt)
           enddo
         enddo
         do m1 = 1, npol*ldim
            do m2 = m1, npol*ldim
               temp = 0.d0
               do i = 1, npol*ldim
                  temp = temp + vet(m1,i)*lambda(i)*CONJG(vet(m2,i))     
               end do
               f(m1,m2) =  temp
               f(m2,m1) =  CONJG(temp) 
            end do
         end do
         do m1 = 1, ldim
            do m2 = 1, ldim
               rho%ns_nc(m1, m2, 1, na) = f(m1, m2) 
               rho%ns_nc(m1, m2, 2, na) = f(m1, ldim+m2) 
               rho%ns_nc(m1, m2, 3, na) = f(ldim+m1, m2) 
               rho%ns_nc(m1, m2, 4, na) = f(ldim+m1, ldim+m2) 
            end do
         end do

       else

         do is = 1, nspin
            do m1 = 1, ldim
               do m2 = 1, ldim 
                  f(m1,m2) = rho%ns(m1,m2,is,na)
               enddo
            enddo
            call cdiagh(ldim, f, ldmx, lambda, vet)
            do i = 1, ldim
               if (starting_ns(i,is,nt) >= 0.d0) lambda(i) = starting_ns(i,is,nt)
            enddo
            do m1 = 1, ldim
               do m2 = m1, ldim
                  temp = 0.d0
                  do i = 1, ldim
                     temp = temp + CONJG(vet(m1,i))*lambda(i)*vet(m2,i)
                  enddo
                  rho%ns(m1,m2,is,na) =  DBLE(temp)
                  rho%ns(m2,m1,is,na) = rho%ns(m1,m2,is,na)
               enddo
            enddo
         enddo

       endif 

      endif
   enddo ! on na

   if (noncolin) then
     CALL write_ns_nc
   else
     CALL write_ns
   endif

   ! reset starting_ns so that this step is not repeated
   starting_ns = -1.0_dp

   return
end subroutine ns_adj
