!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dndtau(dns,ldim,alpha,ipol)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
#include "machine.h"
   USE atom, ONLY: nchi, lchi, oc
   USE basis, ONLY: nat, natomwfc, ityp
   USE klist, ONLY: nks, xk
   USE lsda_mod, ONLY: lsda, nspin, current_spin, isk
   USE ldaU, ONLY: swfcatom, Hubbard_l, Hubbard_U, Hubbard_alpha
   USE wavefunctions_module,    ONLY : evc
   USE uspp, ONLY: nkb, vkb
   USE wvfct, ONLY: nbnd, npwx, npw, igk, wg
   use becmod
   use io_files
#ifdef __PARA
   use para
#endif
   implicit none

   integer ::               &
              ik,           & ! counter on k points
              ibnd,         & !    "    "  bands
              is,           & !    "    "  spins
              i, na, nt, n, alpha, ipol, counter, m1,m2, l
   integer :: ldim
   real (kind=DP) :: &
             dns(ldim,ldim,nspin,nat)
   complex (kind=DP) :: ZDOTC, c_one, c_zero
   integer, allocatable :: offset(:)
   ! offset(nat): offset of d electrons of atom d in the natomwfc ordering
   complex (kind=DP), allocatable :: &
                      proj(:,:), wfcatom(:,:), spsi(:,:),dproj(:,:)
   !                  proj(natomwfc,nbnd), wfcatom(npwx,natomwfc),
   !                  spsi(npwx,nbnd), dproj(natomwfc,nbnd)

   CALL start_clock('dndtau')

   allocate ( offset(nat) )
   allocate ( proj(natomwfc,nbnd), dproj(natomwfc,nbnd), &
              spsi(npwx,nbnd), wfcatom(npwx,natomwfc), becp(nkb,nbnd) )

   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   counter = 0
   do na=1,nat
      offset(na) = 0
      nt=ityp(na)
      do n=1,nchi(nt)
         if (oc(n,nt) >= 0.d0) then
            l=lchi(n,nt)
            if (l == Hubbard_l(nt)) offset(na) = counter
            counter = counter + 2 * l + 1
         end if
      end do
   end do

   if(counter.ne.natomwfc)call errore('new_ns','nstart<>counter',1)

   dns(:,:,:,:) = 0.d0
   !
   !    we start a loop on k points
   !
   if (nks.gt.1) rewind (iunigk)

   do ik = 1, nks

      if (lsda) current_spin = isk(ik)
      !
      ! now we need the first derivative of proj with respect to tau(alpha,ipol)
      !

      if (nks.gt.1) read (iunigk) npw, igk

      call davcio(evc,nwordwfc,iunwfc,ik,-1)
      call davcio(swfcatom,nwordatwfc,iunat,ik,-1)
      c_one= (1.d0, 0.d0)
      c_zero = (0.d0, 0.d0)
      call ZGEMM ('C', 'N', natomwfc, nbnd, npw, c_one, &
                            swfcatom, npwx, evc, npwx, c_zero, proj, natomwfc)
#ifdef __PARA
      call reduce(2*natomwfc*nbnd,proj)
#endif

      call init_us_2 (npw,igk,xk(1,ik),vkb)

      call ccalbec(nkb, npwx, npw, nbnd, becp, vkb, evc)

      call s_psi  (npwx, npw, nbnd, evc, spsi )

      call atomic_wfc( ik, wfcatom )

      dproj(:,:) = (0.d0,0.d0)
      call dprojdtau(dproj,wfcatom,spsi,alpha,ipol,offset(alpha))
      !
      ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
      ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
      !
      do na = 1,nat
         nt = ityp(na)
         if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
            do m1 = 1,ldim
               do m2 = m1,ldim
                  do ibnd = 1,nbnd
                     dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                             wg(ibnd,ik) *            &
                                DREAL(  proj(offset(na)+m1,ibnd)  *   &
                                conjg(dproj(offset(na)+m2,ibnd))  +   &
                                       dproj(offset(na)+m1,ibnd)  *   &
                                conjg(proj(offset(na)+m2,ibnd)) )
                  end do
               end do
            end do
         end if
      end do
   end do                 ! on k-points

#ifdef __PARA
   call poolreduce(ldim*ldim*nspin*nat,dns)
#endif
   !
   ! impose hermiticity of dn_{m1,m2}
   !
   do na = 1,nat
      do is = 1,nspin
         do m1 = 1,ldim
            do m2 = m1+1,ldim
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            end do
         end do
      end do
   end do

   deallocate ( proj, dproj, spsi, wfcatom, becp )
   deallocate ( offset )

   CALL stop_clock('dndtau')
   return
end subroutine dndtau

