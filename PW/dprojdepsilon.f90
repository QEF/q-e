!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dprojdepsilon ( ik,dproj,wfcatom,spsi,ipol,jpol )       
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the strain epsilon(i,j)
   ! (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
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
   integer :: ik, ipol, jpol
   complex (kind=DP) :: &
           dproj(natomwfc,nbnd),   &! output: the derivative of the projection
           wfcatom(npwx,natomwfc), &! input: the atomic wfc
           spsi(npwx,nbnd)          ! input: S|evc>
   integer :: i, ig, jkb2, lmax_wfc, na, m1, ibnd, iwf, nt, ib, ih,jh, &
              nworddw, nworddb
   real (kind=DP) :: xyz(3,3), q, eps, a1, a2
   parameter (eps=1.0e-8)
 
   complex (kind=DP) :: ZDOTC
 
   complex (kind=DP), allocatable :: &
           dwfc(:,:), aux(:,:), work(:), dbeta(:,:), aux1(:,:), &
           betapsi(:,:), dbetapsi(:,:), wfatbeta(:,:), wfatdbeta(:,:)

   !       dwfc(npwx,natomwfc),   ! the derivative of the atomic d wfc
   !       aux(npwx,natomwfc),    ! auxiliary array
   !       work(npwx),            ! the beta function
   !       dbeta(npwx,nkb),       ! the derivative of the beta function
   !       aux1(npwx,nkb),        ! auxiliary array
   !       betapsi(nhm,nbnd),     ! <beta|evc>
   !       dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !       wfatbeta(natomwfc,nhm),! <wfc|beta>
   !       wfatdbeta(natomwfc,nhm)! <wfc|dbeta>

   real (kind=DP), allocatable :: gk(:,:), qm1(:) 
   !       gk(3,npwx),
   !       qm1(npwx) 
 
   ! xyz are the three unit vectors in the x,y,z directions
   xyz(:,:) = 0.d0
   do i=1,3
      xyz(i,i) = 1.d0
   end do

   dproj(:,:) = (0.d0,0.d0)
   !      write(6,*) 'dprojde: ik =',ik,' ipol =',ipol,' jpol =',jpol
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}> 
   !
   allocate ( qm1(npwx), gk(3,npwx) )
   allocate ( dwfc(npwx,natomwfc), aux(npwx,natomwfc) )

   nworddw = 2*npwx*natomwfc
   nworddb = 2*npwx*nkb

   lmax_wfc = 0
   do nt=1, ntyp
      do ib=1,nchi(nt)
         lmax_wfc=max(lmax_wfc,lchi(ib,nt))
      end do
   end do

   ! here the derivative of the Bessel function
   if (ipol*jpol.eq.1) then
      call gen_at_dj (ik,natomwfc,lmax_wfc,dwfc)
      call davcio(dwfc,nworddw,23,ik*2-1,1)
   else
      call davcio(dwfc,nworddw,23,ik*2-1,-1)
   end if

   ! and here the derivative of the spherical harmonic
   if (jpol.eq.1) then
      call gen_at_dy (ik,natomwfc,lmax_wfc,xyz(1,ipol),aux)
      call davcio(aux,nworddw,23,ik*2,1)
   else
      call davcio(aux,nworddw,23,ik*2,-1)
   end if

   do ig = 1,npw
      gk(1,ig) = (xk(1,ik)+g(1,igk(ig)))*tpiba
      gk(2,ig) = (xk(2,ik)+g(2,igk(ig)))*tpiba
      gk(3,ig) = (xk(3,ik)+g(3,igk(ig)))*tpiba
      q = sqrt(gk(1,ig)**2+gk(2,ig)**2+gk(3,ig)**2)
      if (q.gt.eps) then
         qm1(ig)=1.d0/q
      else
         qm1(ig)=0.d0
      end if
      a1 = -1.d0*gk(ipol,ig)
      a2 = -1.d0*gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
      do iwf = 1,natomwfc
         dwfc(ig,iwf) = aux(ig,iwf)*a1 + dwfc(ig,iwf)*a2
         if (ipol.eq.jpol) dwfc(ig,iwf) = dwfc(ig,iwf) - wfcatom(ig,iwf)*0.5d0
         do ibnd = 1,nbnd
            dproj(iwf,ibnd) = dproj(iwf,ibnd)+conjg(dwfc(ig,iwf))*spsi(ig,ibnd)
         end do
      end do
   end do

   a1 = 0.d0
   a2 = 0.d0

#ifdef PARA
   call reduce(2*natomwfc*nbnd,dproj)
#endif

   deallocate ( dwfc, aux )
   !
   ! Now the derivatives of the beta functions: we compute the term 
   ! <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   allocate (dbeta(npwx,nkb), aux1(npwx,nkb), work(npwx), &
             dbetapsi(nhm,nbnd), betapsi(nhm,nbnd), wfatbeta(natomwfc,nhm), &
             wfatdbeta(natomwfc,nhm) )

   ! here the derivative of the Bessel function
   if (ipol*jpol.eq.1) then
      call gen_us_dj (ik,dbeta)
      call davcio(dbeta,nworddb,25,ik*2-1,1)
   else
      call davcio(dbeta,nworddb,25,ik*2-1,-1)
   end if

   ! and here the derivative of the spherical harmonic
   if (jpol.eq.1) then
      call gen_us_dy (ik,xyz(1,ipol),aux1)
      call davcio(aux1,nworddb,25,ik*2,1)
   else
      call davcio(aux1,nworddb,25,ik*2,-1)
   end if

   jkb2 = 0
   do nt=1,ntyp
      do na=1,nat
         if ( ityp(na) .eq. nt ) then
            do ih=1,nh(nt)
               jkb2 = jkb2 + 1
               do ig = 1,npw
                  work(ig) = vkb(ig,jkb2)
                  ! now we compute the true dbeta function
                  dbeta(ig,jkb2) = - aux1(ig,jkb2)*gk(jpol,ig) - &
                        dbeta(ig,jkb2) * gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
                  if (ipol.eq.jpol) &
                     dbeta(ig,jkb2) = dbeta(ig,jkb2) - work(ig)*0.5d0         
               end do
               do ibnd = 1,nbnd
                  betapsi(ih,ibnd)= becp(jkb2,ibnd)
                  dbetapsi(ih,ibnd)= ZDOTC(npw,dbeta(1,jkb2),1,evc(1,ibnd),1)
               end do
               do iwf=1,natomwfc
                  wfatbeta(iwf,ih) = ZDOTC(npw,wfcatom(1,iwf),1,work,1)
                  wfatdbeta(iwf,ih)= ZDOTC(npw,wfcatom(1,iwf),1,dbeta(1,jkb2),1)
               end do
            end do
#ifdef PARA
            call reduce(2*nhm*nbnd,dbetapsi)
            call reduce(2*natomwfc*nhm,wfatbeta)
            call reduce(2*natomwfc*nhm,wfatdbeta)
#endif
            do ibnd = 1,nbnd
               do ih=1,nh(nt)
                  do jh = 1,nh(nt)
                     do iwf=1,natomwfc
                        dproj(iwf,ibnd) = dproj(iwf,ibnd) +               &
                                          qq(ih,jh,nt) *                  &
                               ( wfatdbeta(iwf,ih)*betapsi(jh,ibnd) +     &
                                 wfatbeta(iwf,ih)*dbetapsi(jh,ibnd) )     
                     end do
                  end do
               end do
            end do
         end if
      end do
   end do


   deallocate (dbeta, aux1, work, dbetapsi, betapsi, wfatbeta, wfatdbeta )
   deallocate ( qm1, gk )

   return
end subroutine dprojdepsilon
