!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dprojdtau(dproj,wfcatom,spsi,alpha,ipol,offset)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
   ! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
#include "machine.h"
   use pwcom
   USE wavefunctions,    ONLY : evc
   use becmod
#ifdef __PARA
   use para
#endif
   implicit none
   integer :: &
              alpha,   &! input: the displaced atom
              ipol,    &! input: the component of displacement
              offset    ! input: the offset of the wfcs of the atom "alpha"
   complex (kind=DP) :: &
           wfcatom(npwx,natomwfc), &! input: the atomic wfc
           spsi(npwx,nbnd),        &! input: S|evc>
           dproj(natomwfc,nbnd)     ! output: the derivative of the projection

   integer :: ig, jkb2, na, m1, ibnd, iwf, nt, ib, ih,jh, ldim

   real (kind=DP) :: gvec, a1, a2

   complex (kind=DP):: ZDOTC

   complex (kind=DP), allocatable :: dwfc(:,:), work(:), dbeta(:), &
                                     betapsi(:,:), dbetapsi(:,:), &
                                     wfatbeta(:,:), wfatdbeta(:,:)
   !      dwfc(npwx,ldim),          ! the derivative of the atomic d wfc
   !      work(npwx),            ! the beta function
   !      dbeta(npwx),           ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(natomwfc,nhm),! <wfc|beta>
   !      wfatdbeta(natomwfc,nhm)! <wfc|dbeta>

   nt = ityp(alpha)

   ldim = 2 * Hubbard_l(nt) + 1

   allocate ( dwfc(npwx,ldim), work(npwx), dbeta(npwx), betapsi(nhm,nbnd), &
         dbetapsi(nhm,nbnd), wfatbeta(natomwfc,nhm), wfatdbeta(natomwfc,nhm) )

   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   !

   if (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) then
      do ig = 1,npw
         gvec = g(ipol,igk(ig)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
         do m1 = 1, ldim
            dwfc(ig,m1) = dcmplx(0.d0,-1.d0) * gvec * wfcatom(ig,offset+m1)
         end do
      end do

      call ZGEMM('C','N',ldim, nbnd, npw, (1.d0,0.d0), &
                  dwfc, npwx, spsi, npwx, (0.d0,0.d0), &
                  dproj(offset+1,1), natomwfc)
   end if

#ifdef __PARA
   call reduce(2*natomwfc*nbnd,dproj)
#endif

   jkb2 = 0
   do nt=1,ntyp
      do na=1,nat
         if ( ityp(na) .eq. nt ) then
            do ih=1,nh(nt)
               jkb2 = jkb2 + 1
               if (na.eq.alpha) then
                  do ig = 1, npw
                     gvec = g(ipol,igk(ig)) * tpiba
                     dbeta(ig) = cmplx(0.d0,-1.d0) * vkb(ig,jkb2) * gvec
                     work(ig) = vkb(ig,jkb2)
                  end do
                  do ibnd=1,nbnd
                     dbetapsi(ih,ibnd)= ZDOTC(npw,dbeta,1,evc(1,ibnd),1)
                     betapsi(ih,ibnd) = becp(jkb2,ibnd)
                  end do
                  do iwf=1,natomwfc
                     wfatbeta(iwf,ih) = ZDOTC(npw,wfcatom(1,iwf),1,work,1)
                     wfatdbeta(iwf,ih)= ZDOTC(npw,wfcatom(1,iwf),1,dbeta,1)
                  end do
               end if
            end do
#ifdef __PARA
            call reduce(2*nhm*nbnd,dbetapsi)
            call reduce(2*natomwfc*nhm,wfatbeta)
            call reduce(2*natomwfc*nhm,wfatdbeta)
#endif
            if (na.eq.alpha) then
               do ibnd=1,nbnd
                  do ih=1,nh(nt)
                     do jh=1,nh(nt)
                        do iwf=1,natomwfc
                           dproj(iwf,ibnd) = &
                               dproj(iwf,ibnd) + qq(ih,jh,nt) *         &
                               ( wfatdbeta(iwf,ih)*betapsi(jh,ibnd) +   &
                                  wfatbeta(iwf,ih)*dbetapsi(jh,ibnd) )
                        end do
                     end do
                  end do
               end do
            end if
         end if
      end do
   end do

   deallocate ( dwfc, work, dbeta, betapsi, dbetapsi, wfatbeta, wfatdbeta )

   return
end subroutine dprojdtau
