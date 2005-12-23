!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
module wannier
   USE kinds, only : DP
   integer, allocatable :: nnb(:)       ! #b  (ik)
   integer, allocatable :: kpb(:,:)     ! k+b (ik,ib)
   integer, allocatable :: g_kpb(:,:,:) ! G_k+b (ipol,ik,ib)
   integer, allocatable :: ig_(:,:)     ! G_k+b (ipol,ik,ib)
   integer  :: iun_nnkp, iun_mmn, iun_amn, iun_band, nnbx, n_wannier
   complex(DP), allocatable :: gf(:,:)  ! guding_function(npwx,n_wannier)
   real(DP), allocatable :: rw(:,:)     ! rw(3,n_wannier)
   real(DP), allocatable :: alpha_w(:)  ! alpha_w(n_wannier)
end module wannier
!
!-----------------------------------------------------------------------
program pw2wannier90
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout, ionode
   USE mp_global,  ONLY : mpime, kunit
   USE mp,         ONLY : mp_bcast
   USE cell_base,  ONLY : at, bg
   use lsda_mod,   ONLY : nspin, isk
   use klist,      ONLY : nkstot
   use ktetra,     ONLY : k1, k2, k3, nk1, nk2, nk3
   use io_files,   ONLY : nd_nmbr, prefix, tmp_dir
   !
   implicit none
   integer :: ik, i, kunittmp
   CHARACTER(LEN=4) :: spin_component
   CHARACTER(len=256) :: outdir
   integer :: ispinw
 
   namelist / inputpp / outdir, prefix, spin_component
   !
   call start_postproc (nd_nmbr)
   !
   !   set default values for variables in namelist
   !
   outdir = './'
   prefix = ' '
   spin_component = 'none'
   !
   !     reading the namelist inputpp
   !
   read (5, inputpp)
   !
   !     Check of namelist variables
   !
   tmp_dir = TRIM(outdir) 
   SELECT CASE ( TRIM( spin_component ) )
     CASE ( 'up' )
       ispinw = 1
     CASE ( 'down' )
       ispinw = 2
     CASE DEFAULT
       ispinw = 0
   END SELECT
   !
   !   Now allocate space for pwscf variables, read and check them.
   !
   call read_file  
   call openfil_pp
   !
   call read_nnkp
   !
   call compute_mmn
   !
   call compute_amn
   !
   call write_band
   !
   call stop_pp
   stop
end program pw2wannier90
!
!-----------------------------------------------------------------------
subroutine read_nnkp
   !-----------------------------------------------------------------------
   !
   use kinds, only: DP
   use constants, only: eps8
   use klist, only: nkstot
   use cell_base, only : at, bg
   use gvect, only : g,gg
   use parser, only : find_free_unit
 
   use wannier
   implicit none
   real(DP) :: g_(3), gg_
   integer :: ik, ib, ig, ipol, idum
 
   iun_nnkp = find_free_unit()
 
   open (unit=iun_nnkp, file='wannier.nnkp',form='formatted')
 
   allocate ( nnb(nkstot) )
! get nnbx
   nnbx=0
   do ik=1,nkstot
      read (iun_nnkp,*) nnb(ik)
      nnbx = max (nnbx, nnb(ik) )
      do ib = 1, nnb(ik)
         read(iun_nnkp,*)
      end do
   end do
! allocate
   allocate ( kpb(nkstot,nnbx), g_kpb(3,nkstot,nnbx),ig_(nkstot,nnbx) )
 
! read
   rewind(iun_nnkp)
   do ik=1,nkstot
      read (iun_nnkp,*) idum
      if (idum /= nnb(ik)) call errore('pw2w','ekkekka ',1)
      do ib = 1, nnb(ik)
         read(iun_nnkp,*) idum, kpb(ik,ib), (g_kpb(ipol,ik,ib), ipol =1,3)
         g_(:) = REAL( g_kpb(:,ik,ib) )
         call trnvect (g_, at, bg, 1)
         gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
         ig_(ik,ib) = 0
         ig = 1
         do while  (gg(ig) <= gg_ + eps8) 
            if ( (abs(g(1,ig)-g_(1)) < eps8) .and.  &
                 (abs(g(2,ig)-g_(2)) < eps8) .and.  &
                 (abs(g(3,ig)-g_(3)) < eps8)  ) ig_(ik,ib) = ig
            ig= ig +1
         end do
!         write (*,'(4i6,3f10.6)') ik, ib, kpb(ik,ib), ig_(ik,ib), g_(:)
      end do
   end do
 
   close (iun_nnkp)

   return

end subroutine read_nnkp
!
subroutine compute_mmn
   !-----------------------------------------------------------------------
   !
   use kinds, only: DP
   use wvfct, only : nbnd, npw, npwx, igk, g2kin
   use wavefunctions_module, only : evc, psic
   use gsmooth, only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
   use klist, only : nkstot, xk
   use io_files, only : nwordwfc, iunwfc
   use gvect, only : g, ngm, ecutwfc
   use cell_base, only : tpiba2
   use parser, only : find_free_unit
   use wannier
   implicit none
   integer :: mmn_tot, ik, ikp, ipol, ib, ibnd, jbnd, npwq
   complex(DP), allocatable :: phase(:), aux(:), evcq(:,:)
   integer, allocatable :: igkq(:)
   complex(DP) :: mmn, ZDOTC
   real(DP)::aa

   allocate( phase(nrxxs), aux(npwx), evcq(npwx,nbnd), igkq(npwx) )
   
   iun_mmn = find_free_unit()
   open (unit=iun_mmn, file='wannier.mmn',form='formatted')

   mmn_tot = 0
   do ik=1,nkstot
      mmn_tot = mmn_tot + nnb(ik) * nbnd * nbnd
   end do

   write (*,*) "MMN"
   write (iun_mmn,*) mmn_tot
   do ik=1,nkstot
   write (*,*) ik
! read wfc at k
      call davcio (evc, nwordwfc, iunwfc, ik, -1 )
      call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)

      do ib=1,nnb(ik)
         ikp = kpb(ik,ib)
! read wfc at k+b
         call davcio (evcq, nwordwfc, iunwfc, ikp, -1 )
         call gk_sort (xk(1,ikp), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
! compute the phase
         phase(:) = (0.d0,0.d0)
         if ( ig_(ik,ib)>0) phase( nls(ig_(ik,ib)) ) = (1.d0,0.d0)
         call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
! loops on bands
         do ibnd=1,nbnd
            psic(:) = (0.d0, 0.d0)
            psic(nls (igk (1:npw) ) ) = evc (1:npw, ibnd)
            call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
            psic(1:nrxxs) = psic(1:nrxxs) * phase(1:nrxxs)
            call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
            aux(1:npwq) = psic(nls (igkq(1:npwq) ) )
            aa = 0.d0
            do jbnd=1,nbnd   
               mmn = ZDOTC ( npwq, aux,1,evcq(1,jbnd),1)
               call reduce(mmn,2)
               write (iun_mmn,'(7i5)') ibnd, jbnd, ik, ikp, &
                                       (g_kpb(ipol,ik,ib), ipol=1,3)
               write (iun_mmn,'(2f18.12)') mmn
               aa = aa + abs(mmn)**2
            end do
!            write (*,*) ik,ib,ibnd, aa
         end do
      end do
   end do
   close (iun_mmn)
   return
end subroutine compute_mmn
!
!-----------------------------------------------------------------------
subroutine compute_amn
   !-----------------------------------------------------------------------
   !
   use kinds, only: DP
   use klist, only: nkstot, xk
   use wvfct, only : nbnd, npw, npwx, igk, g2kin
   use wavefunctions_module, only : evc
   use io_files, only : nwordwfc, iunwfc
   use gvect, only : g, ngm, ecutwfc
   use cell_base, only : tpiba2
   use parser, only : find_free_unit
   use wannier
   implicit none
   complex(DP) :: amn, ZDOTC
   integer :: amn_tot, ik, ibnd, iw

   call read_gf_definition

   iun_amn = find_free_unit()
   open (unit=iun_amn, file='wannier.amn',form='formatted')

   amn_tot = nkstot * nbnd * n_wannier

   write (*,*) "AMN"
   write (iun_amn,*) amn_tot
   do ik=1,nkstot
      write (*,*) ik
      call davcio (evc, nwordwfc, iunwfc, ik, -1 )
      call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      call generate_guiding_functions(ik)
      do ibnd=1,nbnd
         do iw=1,n_wannier   
            amn = ZDOTC(npw,evc(1,ibnd),1,gf(1,iw),1)
            call reduce(amn,2)
            write(iun_amn,'(3i5,2f18.12)') ibnd,iw,ik, amn
         end do
      end do
   end do
   close (iun_amn)
   return
end subroutine compute_amn
!
subroutine read_gf_definition
!
! read guiding functions definition
!
   use wvfct, only : npwx
   use wannier
   implicit none 
   integer :: i, iw

   read (*,*) n_wannier
   allocate( rw(3,n_wannier), alpha_w(n_wannier), gf(npwx,n_wannier) )
   do iw=1,n_wannier
      read (*,*) (rw(i,iw),i=1,3), alpha_w(iw)
   end do
! convert them in carthesian coordinated in unit of alat
!   rw(:,:) = rw(:,:) / celldm(1)/0.529177d0
   return
end subroutine read_gf_definition
!
subroutine generate_guiding_functions(ik)
   use constants, only : pi, tpi, fpi
   use wvfct, only : npw, g2kin
   use gvect, only : g
   use cell_base,  ONLY :tpiba2, omega
   use wannier
   implicit none
   integer iw, ig, ik
   real(DP) :: arg, anorm, fac, alpha_w2
   complex(DP) :: ZDOTC

   do iw =1, n_wannier
      fac = ( sqrt(fpi)*alpha_w(iw) )**(1.5) / sqrt(omega)
      alpha_w2 = alpha_w(iw)**2
      do ig=1,npw
         arg = tpi * (g(1,ig)*rw(1,iw)+g(2,ig)*rw(2,iw)+g(3,ig)*rw(3,iw))
         gf(ig,iw) = fac * exp(-0.5d0*alpha_w2*g2kin(ig)*tpiba2) * &
                           CMPLX(cos(arg),-sin(arg))
      end do
      anorm = ZDOTC(npw,gf(1,iw),1,gf(1,iw),1)
      write (*,*) iw, ik, anorm
   end do
   return
end subroutine generate_guiding_functions

subroutine write_band
   use wvfct, only : nbnd, et
   use klist, only : nkstot
   use constants, only: rytoev
   use parser, only : find_free_unit
   use wannier
   implicit none
   integer ik, ibnd
   iun_band = find_free_unit()
   open (unit=iun_band, file='BAND.dat',form='formatted')
   do ik=1,nkstot
      do ibnd=1,nbnd
         write (iun_band,'(2i5,f18.12)') ibnd, ik, et(ibnd,ik)*rytoev
      end do
   end do
   return
end subroutine write_band
