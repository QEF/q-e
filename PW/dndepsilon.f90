!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dndepsilon ( dns,ipol,jpol )
   !-----------------------------------------------------------------------
   ! This routine computes the derivative of the ns atomic occupations with
   ! respect to the strain epsilon(ipol,jpol) used to obtain the hubbard
   ! contribution to the internal stres tensor.
   !
#include "machine.h"
   use pwcom
   use becmod
   use io
#ifdef PARA
   use para
#endif
   implicit none
   !
   ! I/O variables first
   !
   real (kind=DP) :: dns(nat,nspin,5,5)
   integer :: ipol, jpol
   !
   ! local variable
   !
   integer :: ik,    & ! counter on k points
              ibnd,  & !    "    "  bands
              is,    & !    "    "  spins
              i, na, nt, n, counter, m1, m2, l
   complex (kind=DP) :: ZDOTC

   integer, allocatable :: offset(:)
   ! offset(nat)  ! offset of d electrons of atom d in the natomwfc ordering
   complex (kind=DP), allocatable :: &
                      proj(:,:), wfcatom(:,:), spsi(:,:), dproj(:,:)
   ! proj(natomwfc,nbnd), wfcatom(npwx,natomwfc),
   ! spsi(npwx,nbnd), dproj(natomwfc,nbnd)

   allocate (offset(nat), proj(natomwfc,nbnd), wfcatom(npwx,natomwfc),  &
             spsi(npwx,nbnd), dproj(natomwfc,nbnd) )


   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   counter = 0
   do na=1,nat
      offset(na) = 0
      nt=ityp(na)
      do n=1,nchi(nt)
         if (oc(n,nt).gt.0.d0.or..not.newpseudo(nt)) then
            l=lchi(n,nt)
            if (l.eq.2) offset(na) = counter
            counter = counter + 2 * l + 1
         end if
      end do
   end do

   if(counter.ne.natomwfc) call error('new_ns','nstart<>counter',1)

   dns(:,:,:,:) = 0.d0
   !
   !    we start a loop on k points
   !
   if (nks.gt.1) rewind (iunigk)

   do ik = 1, nks

      if (nks.gt.1) read (iunigk) npw, igk

      !
      ! now we need the first derivative of proj with respect to
      ! epsilon(ipol,jpol)
      !
      call davcio(evc,nwordwfc,iunwfc,ik,-1)
      call init_us_2 (npw,igk,xk(1,ik),vkb)
      call ccalbec(nkb, npwx, npw, nbnd, becp, vkb, evc)

      call s_psi  (npwx, npw, nbnd, evc, spsi )
      call atomic_wfc( ik, wfcatom )

      dproj(:,:) = (0.d0,0.d0)

      call dprojdepsilon(ik,dproj,wfcatom,spsi,ipol,jpol)

      call davcio(swfcatom,nwordatwfc,iunat,ik,-1)

      do ibnd = 1, nbnd
         do i=1,natomwfc
            proj(i,ibnd) = ZDOTC(npw,swfcatom(1,i),1,evc(1,ibnd),1)
         enddo
      enddo

#ifdef PARA
       call reduce(2*natomwfc*nbnd,proj)
#endif
      !
      ! compute the derivative of the occupation numbers (quantities dn(m1,m2))
      ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
      !
      do na = 1,nat
         nt = ityp(na)
         if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
            do m1 = 1,5
               do m2 = m1,5
                  do ibnd = 1,nbnd
                     dns(na,isk(ik),m1,m2) = dns(na,isk(ik),m1,m2) + &
                                             wg(ibnd,ik) *           &
                              DREAL( proj(offset(na)+m1,ibnd) *      &
                              conjg(dproj(offset(na)+m2,ibnd) ) +    &
                                    dproj(offset(na)+m1,ibnd)*       &
                              conjg( proj(offset(na)+m2,ibnd) ) )
                  end do
               end do
            end do
         end if
      end do

   end do                 ! on k-points

#ifdef PARA
   call poolreduce(nat*nspin*25,dns)
#endif
   !
   ! impose hermeticity of dn_{m1,m2}
   !
   do na = 1,nat
      do is = 1,nspin
         do m1 = 1,5
            do m2 = m1+1,5
               dns(na,is,m2,m1) = dns(na,is,m1,m2)
            end do
         end do
      end do
   end do

   deallocate (offset, proj, wfcatom, spsi, dproj )

   return
end subroutine dndepsilon
