!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
   use wannier
   !
   implicit none
   integer :: ik, i, kunittmp
   CHARACTER(LEN=4) :: spin_component
   CHARACTER(len=256) :: outdir
   ! these are in wannier module.....-> integer :: ispinw, ikstart, ikstop, iknum
 
   namelist / inputpp / outdir, prefix, spin_component, wan_mode, seedname
   !
   call start_postproc (nd_nmbr)
   !
   !   set default values for variables in namelist
   !
   outdir = './'
   prefix = ' '
   seedname = 'wannier'
   spin_component = 'none'
   wan_mode = 'standalone'
   !
   !     reading the namelist inputpp
   !
   read (5, inputpp)
   !
   !     Check of namelist variables
   !
   tmp_dir = TRIM(outdir) 
   !
   !   Now allocate space for pwscf variables, read and check them.
   !
   logwann = .true.
   write(stdout,*)
   write(stdout,*) ' Reading nscf_save data'
   call read_file  
   write(stdout,*)
   !
   SELECT CASE ( TRIM( spin_component ) )
     CASE ( 'up' )
       write(stdout,*) ' Spin CASE ( up )'
       ispinw  = 1
       ikstart = 1
       ikstop  = nkstot/2
       iknum   = nkstot/2
     CASE ( 'down' )
       write(stdout,*) ' Spin CASE ( down )'
       ispinw = 2
       ikstart = nkstot/2 + 1
       ikstop  = nkstot
       iknum   = nkstot/2
     CASE DEFAULT
       write(stdout,*) ' Spin CASE ( default = unpolarized )'
       ispinw = 0
       ikstart = 1
       ikstop  = nkstot
       iknum   = nkstot
   END SELECT
   !
   write(stdout,*)
   write(stdout,*) ' Wannier mode is: ',wan_mode
   write(stdout,*)
   !
   if(wan_mode.eq.'standalone') then
   !
      write(stdout,*) ' -----------------'
      write(stdout,*) ' *** Reading nnkp '
      write(stdout,*) ' -----------------'
      write(stdout,*)
      call read_nnkp
      write(stdout,*) ' Opening pp-files '
      call openfil_pp
      call ylm_expansion
      write(stdout,*)
      write(stdout,*)
      write(stdout,*) ' ---------------'
      write(stdout,*) ' *** Compute  A '
      write(stdout,*) ' ---------------'
      write(stdout,*)
      call compute_amn
      write(stdout,*)
      write(stdout,*) ' ---------------'
      write(stdout,*) ' *** Compute  M '
      write(stdout,*) ' ---------------'
      write(stdout,*) 
      call compute_mmn
      write(stdout,*)
      write(stdout,*) ' ----------------'
      write(stdout,*) ' *** Write bands '
      write(stdout,*) ' ----------------'
      write(stdout,*)
      call write_band
      write(stdout,*)
      write(stdout,*) ' --------------------'
      write(stdout,*) ' *** Write plot info '
      write(stdout,*) ' --------------------'
      write(stdout,*)
      call write_plot
      write(stdout,*)
      write(stdout,*) ' ------------'
      write(stdout,*) ' *** Stop pp '
      write(stdout,*) ' ------------' 
      write(stdout,*)
      call stop_pp
   !
   endif
   !
   if(wan_mode.eq.'library') then
   !
      !call setup_nnkp
      call compute_mmn
      call compute_amn
      call write_band
      !call run_wannier
      call stop_pp
   !
   endif
   !
   if(wan_mode.eq.'wannier2sic') then
   !
      call read_nnkp
      call wan2sic
   !
   endif
   !
   stop
end program pw2wannier90
!
!-----------------------------------------------------------------------
subroutine read_nnkp
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout
   use kinds,     only: DP
   use constants, only : eps6
   use cell_base, only : at, bg, alat
   use gvect,     only : g, gg
   use io_files,  only : find_free_unit
   use klist,     only : nkstot, xk
   use wvfct,     only : npwx, nbnd
   use wannier
   use constants,       only : tpi

   implicit none
   real(DP) :: g_(3), gg_
   integer :: ik, ib, ig, ipol, iw, idum, indexb
   integer numwan, numk, i, j
   real(DP) :: rlatt(3,3), glatt(3,3), xx(3), bohr, box, xnorm, znorm, coseno
   CHARACTER(LEN=80) :: line1, line2

   iun_nnkp = find_free_unit()
   open (unit=iun_nnkp, file=TRIM(seedname)//".nnkp",form='formatted')
   nnbx=0

!   check the information from *.nnkp with the nscf_save data
   write(stdout,*) ' Checking info from wannier.nnkp file' 
   write(stdout,*)
   bohr = 0.5291772108d0

   call scan_file_to('real_lattice')
   do j=1,3
   read(iun_nnkp,*) (rlatt(i,j),i=1,3)
      do i = 1,3
         rlatt(i,j) = rlatt(i,j)/(alat*bohr)
      enddo
   enddo
   do j=1,3
      do i=1,3
         if(abs(rlatt(i,j)-at(i,j)).gt.eps6) then
            write(stdout,*)  ' Something wrong! '
            write(stdout,*)  ' rlatt(i,j) =',rlatt(i,j),  ' at(i,j)=',at(i,j)
            stop
         endif  
      enddo
   enddo
   write(stdout,*) ' - Real lattice is ok'

   call scan_file_to('recip_lattice')
   do j=1,3
      read(iun_nnkp,*) (glatt(i,j),i=1,3)
      do i = 1,3
         glatt(i,j) = (alat*bohr)*glatt(i,j)/tpi
      enddo
   enddo
   do j=1,3
      do i=1,3
         if(abs(glatt(i,j)-bg(i,j)).gt.eps6) then
            write(stdout,*)  ' Something wrong! '
            write(stdout,*)  ' glatt(i,j)=',glatt(i,j), ' bg(i,j)=',bg(i,j)
            stop
         endif
      enddo
   enddo
   write(stdout,*) ' - Reciprocal lattice is ok'

   call scan_file_to('kpoints')
   read(iun_nnkp,*) numk
   if(numk.ne.iknum) then
      write(stdout,*)  ' Something wrong! '
      write(stdout,*)  ' numk=',numk, ' iknum=',iknum
      stop
   endif
   do i=1,numk
      read(iun_nnkp,*) xx(1), xx(2), xx(3)
      CALL cryst_to_cart( 1, xx, bg, 1 )
      if(abs(xx(1)-xk(1,i)).gt.eps6.or. &
         abs(xx(2)-xk(2,i)).gt.eps6.or. &
         abs(xx(3)-xk(3,i)).gt.eps6) then
         write(stdout,*)  ' Something wrong! '
         write(stdout,*) ' k-point ',i,' is wrong'
         write(stdout,*) xx(1), xx(2), xx(3) 
         write(stdout,*) xk(1,i), xk(2,i), xk(3,i)
         stop
      endif
   enddo
   write(stdout,*) ' - K-points are ok'

   call scan_file_to('projections')
   read(iun_nnkp,*) numwan
   n_wannier = numwan
   allocate( center_w(3,n_wannier), alpha_w(n_wannier), gf(npwx,n_wannier), &
             l_w(n_wannier), mr_w(n_wannier), r_w(n_wannier), &
             zaxis(3,n_wannier), xaxis(3,n_wannier), csph(16,n_wannier) )
   write(stdout,'("  - Number of wannier functions is ok (",i3,")")') n_wannier 
   do iw=1,numwan
      read(iun_nnkp,*) (center_w(i,iw), i=1,3), l_w(iw), mr_w(iw), r_w(iw)
      read(iun_nnkp,*) (zaxis(i,iw),i=1,3),(xaxis(i,iw),i=1,3),alpha_w(iw),box
      xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
                                             xaxis(3,iw)*xaxis(3,iw))
      if (xnorm < eps6) call errore ('read_nnkp',' |xaxis| < eps ',1)
      znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
                                             zaxis(3,iw)*zaxis(3,iw))
      if (znorm < eps6) call errore ('read_nnkp',' |zaxis| < eps ',1)
      coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
                                          xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
      if (abs(coseno) > eps6) &
          call errore('read_nnkp',' xaxis and zaxis are not orthogonal !',1)
      if (alpha_w(iw) < eps6) &
          call errore('read_nnkp',' zona value must be positive', 1)
      ! convert wannier center in cartesian coordinates (in unit of alat)
      CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
   enddo
   write(stdout,*) ' - All guiding functions are given '
   !
   call scan_file_to('nnkpts')
   read (iun_nnkp,*) nnb 
   nnbx = max (nnbx, nnb )
!  end of check
   allocate ( kpb(iknum,nnbx), g_kpb(3,iknum,nnbx),ig_(iknum,nnbx) )
!  read data about neighbours
   write(stdout,*)
   write(stdout,*) ' Reading data about k-point neighbours '
   write(stdout,*)
   do ik=1, iknum
      do ib = 1, nnb
         read(iun_nnkp,*) idum, kpb(ik,ib), (g_kpb(ipol,ik,ib), ipol =1,3)
         g_(:) = REAL( g_kpb(:,ik,ib) )
         call trnvect (g_, at, bg, 1)
         gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
         ig_(ik,ib) = 0
         ig = 1
         do while  (gg(ig) <= gg_ + eps6) 
            if ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                 (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                 (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
            ig= ig +1
         end do
      end do
   end do
   write(stdout,*) ' All neighbours are found '
   write(stdout,*)

   call scan_file_to('exclude_bands')
   read (iun_nnkp,*) nexband
   allocate( excluded_band(nbnd) )
   excluded_band(1:nbnd)=.false.
   do i=1,nexband
     read(iun_nnkp,*) indexb
     if (indexb<1 .or. indexb>nbnd) &
          call errore('read_nnkp',' wrong excluded band index ', 1)
     excluded_band(indexb)=.true.
   enddo
    
   close (iun_nnkp)
   return
end subroutine read_nnkp

subroutine scan_file_to (keyword)
   use wannier, only :iun_nnkp
   USE io_global,  ONLY : stdout
   implicit none
   character(len=*) :: keyword
   character(len=80) :: line1, line2
!
! by uncommenting the following line the file scan restarts every time 
! from the beginning thus making the reading independent on the order 
! of data-blocks
!   rewind (iun_nnkp)
!
10 continue
   read(iun_nnkp,*,end=20) line1, line2
   if(line1.ne.'begin')  goto 10
   if(line2.ne.keyword) goto 10
   return
20 write (stdout,*) keyword," data-block missing "
   stop
end subroutine scan_file_to
!
subroutine compute_mmn
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout
   use kinds,           only: DP
   use wvfct,           only : nbnd, npw, npwx, igk, g2kin
   use wavefunctions_module, only : evc, psic
   use gsmooth,         only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
   use klist,           only : nkstot, xk
   use io_files,        only : nwordwfc, iunwfc
   use io_files,        only : find_free_unit
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2, omega, alat, tpiba, at, bg
   USE ions_base,       only : nat, ntyp => nsp, ityp, tau
   use constants,       only : tpi
   use uspp,            only : nkb, vkb
   USE uspp_param,      ONLY : nh, tvanp, lmaxq
   use becmod,          only : becp
   use wannier
   implicit none
   integer :: mmn_tot, ik, ikp, ipol, ib, npwq, i, m, n
   integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nbt
   integer :: ikevc, ikpevcq
   complex(DP), allocatable :: phase(:), aux(:), evcq(:,:), becp2(:,:), Mkb(:,:)
   complex(DP), allocatable :: qb(:,:,:,:), qgm(:)
   real(DP), allocatable    :: qg(:), ylm(:,:), dxk(:,:)
   integer, allocatable     :: igkq(:)
   complex(DP)              :: mmn, ZDOTC, phase1
   real(DP)                 :: aa, arg, g_(3)

   allocate( phase(nrxxs), aux(npwx), evcq(npwx,nbnd), igkq(npwx) )
   
   iun_mmn = find_free_unit()
   open (unit=iun_mmn, file=TRIM(seedname)//".mmn",form='formatted')

   mmn_tot = 0
   do ik=1,iknum
      mmn_tot = mmn_tot + nnb * nbnd * nbnd
   end do
   !
   !   USPP
   !
   CALL init_us_1
   allocate ( becp(nkb,nbnd),becp2(nkb,nbnd))
   !
   !     qb is  FT of Q(r) 
   !
   nbt = 0
   do ik=1, iknum
      nbt = nbt + nnb
   enddo
   !
   allocate( qg(nbt) )
   allocate (dxk(3,nbt))
   !
   ind = 0
   do ik=1,iknum
      do ib=1,nnb
         ind = ind + 1
         ikp = kpb(ik,ib) 
         !
         g_(:) = REAL( g_kpb(:,ik,ib) )
         call trnvect (g_, at, bg, 1)
         dxk(:,ind) = xk(:,ikp) +g_(:) - xk(:,ik) 
         qg(ind) = dxk(1,ind)*dxk(1,ind)+dxk(2,ind)*dxk(2,ind)+dxk(3,ind)*dxk(3,ind)
      enddo
      write (*,'(i3,12f8.4)')  ik, qg((ik-1)*nnb+1:ik*nnb)
   enddo 
   !
   allocate( ylm(nbt,lmaxq*lmaxq), qgm(nbt) )
   allocate( qb (nkb, nkb, ntyp, nbt) )
   !
   call ylmr2 (lmaxq*lmaxq, nbt, dxk, qg, ylm)
   qg(:) = sqrt(qg(:)) * tpiba
   !
   do nt = 1, ntyp
      if (tvanp (nt) ) then 
         do ih = 1, nh (nt)
            do jh = 1, nh (nt)
               CALL qvan2 (nbt, ih, jh, nt, qg, qgm, ylm)
               qb (ih, jh, nt, 1:nbt) = omega * qgm(1:nbt)
            enddo 
         enddo 
      endif 
   enddo 
   !
   deallocate (qg, qgm, ylm )
   !
   write (stdout,*) "MMN"
   write (iun_mmn,*) mmn_tot
   write (iun_mmn,*) nbnd-nexband, iknum, nnb 
   !
   allocate( Mkb(nbnd,nbnd) )
   !
  write(stdout,*) " iknum = ",iknum

   ind = 0
   do ik=1,iknum
   write (stdout,*) ik
      ikevc = ik + ikstart - 1 
      call davcio (evc, nwordwfc, iunwfc, ikevc, -1 )
      call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      !
      !  USPP
      !
      call init_us_2 (npw, igk, xk(1,ik), vkb)
      ! below we compute the product of beta functions with |psi> 
      call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
      !
      !
      !do ib=1,nnb(ik)
      do ib=1,nnb
         ind = ind + 1
         ikp = kpb(ik,ib)
! read wfc at k+b
         ikpevcq = ikp + ikstart - 1
         call davcio (evcq, nwordwfc, iunwfc, ikpevcq, -1 )
         call gk_sort (xk(1,ikp), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
! compute the phase
         phase(:) = (0.d0,0.d0)
         if ( ig_(ik,ib)>0) phase( nls(ig_(ik,ib)) ) = (1.d0,0.d0)
         call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
         !
         !  USPP
         !
         call init_us_2 (npwq, igkq, xk(1,ikp), vkb)
         ! below we compute the product of beta functions with |psi> 
         call ccalbec (nkb, npwx, npwq, nbnd, becp2, vkb, evcq)
         !
         !
         Mkb(:,:) = (0.0d0,0.0d0) 
         ijkb0 = 0
         do nt = 1, ntyp
            if ( tvanp(nt) ) then
               do na = 1, nat
                  !
                  arg = DOT_PRODUCT( dxk(:,ind), tau(:,na) ) * tpi 
                  phase1 = CMPLX ( COS(arg), -SIN(arg) )
                  !
                  if ( ityp(na) == nt ) then
                     do jh = 1, nh(nt)
                        jkb = ijkb0 + jh
                        do ih = 1, nh(nt)
                           ikb = ijkb0 + ih
                           !
                           do m = 1,nbnd
                           do n = 1,nbnd
                           Mkb(m,n) = Mkb(m,n) + &
                                  phase1 * qb(ih,jh,nt,ind) * &
                                  CONJG( becp(ikb,m) ) * becp2(jkb,n) 
                           enddo
                           enddo
                        enddo !ih
                     enddo !jh
                     ijkb0 = ijkb0 + nh(nt)
                  endif  !ityp
               enddo  !nat 
            else  !tvanp
               do na = 1, nat
                  if ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
               enddo
            endif !tvanp
         enddo !ntyp
         !
         !
! loops on bands
         !
         write (iun_mmn,'(7i5)') ik, ikp, (g_kpb(ipol,ik,ib), ipol=1,3)
         
         !
         do m=1,nbnd
            psic(:) = (0.d0, 0.d0)
            psic(nls (igk (1:npw) ) ) = evc (1:npw, m)
            call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
            psic(1:nrxxs) = psic(1:nrxxs) * phase(1:nrxxs)
            call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
            aux(1:npwq) = psic(nls (igkq(1:npwq) ) )
            aa = 0.d0
            !
            !
            do n=1,nbnd   
              !
              mmn = ZDOTC (npwq, aux,1,evcq(1,n),1)
              !
              !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
              !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 > 
              !
              call reduce(2,mmn)
              Mkb(m,n) = mmn + Mkb(m,n)
              !
              aa = aa + abs(mmn)**2
              !
            end do ! n
         end do   ! m

         do n=1,nbnd
            if (excluded_band(n)) cycle
            do m=1,nbnd
               if (excluded_band(m)) cycle
               write (iun_mmn,'(2f18.12)') Mkb(m,n)
            enddo
         enddo

      end do !ikp
   end do  !ik
!
   deallocate (becp, becp2, qb, Mkb, dxk)
!
   close (iun_mmn)
   deallocate( phase, aux, evcq, igkq )
   return
end subroutine compute_mmn
!
!-----------------------------------------------------------------------
subroutine compute_amn
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout
   use kinds,           only : DP
   use klist,           only : nkstot, xk
   use wvfct,           only : nbnd, npw, npwx, igk, g2kin
   use wavefunctions_module, only : evc
   use io_files,        only : nwordwfc, iunwfc
   use io_files,        only : find_free_unit
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2
   use uspp,            only : nkb, vkb
   use becmod,          only : becp
   use wannier
   implicit none
   complex(DP) :: amn, ZDOTC
   complex(DP), allocatable :: sgf(:,:)
   integer :: amn_tot, ik, ibnd, ibnd1, iw,i
   integer :: ikevc

   !call read_gf_definition.....>   this is done at the begin

   iun_amn = find_free_unit()
   open (unit=iun_amn, file=TRIM(seedname)//".amn",form='formatted')
   amn_tot = iknum * nbnd * n_wannier
   write (stdout,*) "AMN"
   write (iun_amn,*) amn_tot
   write (iun_amn,*) nbnd-nexband,  iknum, n_wannier 
   !
   allocate( sgf(npwx,n_wannier))
   allocate ( becp(nkb,n_wannier))
   CALL init_us_1
   !
   do ik=1,iknum
      write (stdout,*) ik
      ikevc = ik + ikstart - 1
      call davcio (evc, nwordwfc, iunwfc, ikevc, -1 )
      call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      call generate_guiding_functions(ik)   ! they are called gf(npw,n_wannier)
      !
      !  USPP
      !
      call init_us_2 (npw, igk, xk (1, ik), vkb)
!      ! below we compute the product of beta functions with trial func.
      call ccalbec (nkb, npwx, npw, n_wannier, becp, vkb, gf)
!      ! and we use it for the product S|trial_func>
      call s_psi (npwx, npw, n_wannier, gf, sgf)  
       ! sgf(:,:) = gf(:,:)
       !
      do iw = 1,n_wannier
         ibnd1 = 0 
         do ibnd = 1,nbnd
            amn = ZDOTC(npw,evc(1,ibnd),1,sgf(1,iw),1) 
            call reduce(2,amn)
            if (excluded_band(ibnd)) cycle
            ibnd1=ibnd1+1
            write(iun_amn,'(3i5,2f18.12)') ibnd1, iw, ik, amn
         end do
      end do
   end do  ! k-points
   deallocate (sgf)
   deallocate (becp)
   deallocate (csph)
   !
   close (iun_amn)
   return
end subroutine compute_amn
!
!
subroutine generate_guiding_functions(ik)
   !
   USE io_global,  ONLY : stdout
   use constants, only : pi, tpi, fpi, eps8
   use wvfct, only : npw, g2kin, igk
   use gvect, only : ig1, ig2, ig3, g
   use cell_base,  ONLY : tpiba2, omega, tpiba
   use wannier
   use klist,      only : xk 
   USE cell_base, ONLY : bg
   implicit none
   integer, parameter :: lmax=3, lmax2=(lmax+1)**2
   integer :: iw, ig, ik, bgtau(3), isph, l, mesh_r
   integer :: lmax_iw, lm, ipol, n1, n2, n3, nr1, nr2, nr3, iig
   real(DP) :: arg, anorm, fac, alpha_w2, yy, alfa
   complex(DP) :: ZDOTC, kphase, lphase, gff, lph
   real(DP), allocatable :: gk(:,:), qg(:), ylm(:,:), radial(:,:)
   complex(DP), allocatable :: sk(:) 
   !
   allocate( gk(3,npw), qg(npw), ylm(npw,lmax2), sk(npw), radial(npw,0:lmax) )
   !
   do ig = 1, npw
      gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
      gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
      gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
      qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
   enddo

   call ylmr2 (lmax2, npw, gk, qg, ylm)
   ! define qg as the norm of (k+g) in a.u.
   qg(:) = sqrt(qg(:)) * tpiba

   do iw = 1, n_wannier
      !
      gf(:,iw) = (0.d0,0.d0)

      call radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial) 

      do lm = 1, lmax2
         if ( abs(csph(lm,iw)) < eps8 ) cycle
         l = int (sqrt( lm-1.d0))
         lphase = (0.d0,-1.d0)**l
         !
         do ig=1,npw
            gf(ig,iw) = gf(ig,iw) + csph(lm,iw) * ylm(ig,lm) * radial(ig,l) * lphase
         end do !ig
      end do ! lm
      do ig=1,npw
         iig = igk(ig)
         arg = ( gk(1,ig)*center_w(1,iw) + gk(2,ig)*center_w(2,iw) + &
                                           gk(3,ig)*center_w(3,iw) ) * tpi
         ! center_w are cartesian coordinates in units of alat 
         sk(ig) = CMPLX(cos(arg), -sin(arg) )
         gf(ig,iw) = gf(ig,iw) * sk(ig) 
      end do
      anorm = REAL(ZDOTC(npw,gf(1,iw),1,gf(1,iw),1))
      call reduce(1,anorm)
      write (stdout,*) ik, iw, anorm
      gf(:,iw) = gf(:,iw) / dsqrt(anorm)
   end do
   !
   deallocate ( gk, qg, ylm, sk, radial)
   return
end subroutine generate_guiding_functions

subroutine write_band
   USE io_global,  ONLY : stdout
   use wvfct, only : nbnd, et
   use klist, only : nkstot
   use constants, only: rytoev
   use io_files, only : find_free_unit
   use wannier
   implicit none
   integer ik, ibnd, ikevc
   iun_band = find_free_unit()
   open (unit=iun_band, file=TRIM(seedname)//".eig",form='formatted')
   do ik=ikstart,ikstop
      ikevc = ik - ikstart + 1
      do ibnd=1,nbnd
         write (iun_band,'(2i5,f18.12)') ibnd, ikevc, et(ibnd,ik)*rytoev
      end do
   end do
   return
end subroutine write_band

subroutine write_plot
   USE io_global,  ONLY : stdout
   use wvfct, only : nbnd, npw, igk, g2kin
   use wavefunctions_module, only : evc, psic
   use io_files, only : find_free_unit, nwordwfc, iunwfc
   use wannier
   use gsmooth,         only : nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
   use klist,           only : nkstot, xk
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2

   implicit none
   integer ik, ibnd, ibnd1, ikevc, i1
#ifdef __PARA
   integer nxxs
   COMPLEX(DP),allocatable :: psic_all(:)

   nxxs = nrx1s * nrx2s * nrx3s
   allocate(psic_all(nxxs) )

#endif

   
   iun_plot = find_free_unit()

   open (unit=iun_plot, file=TRIM(seedname)//".bloch",form='formatted')
   do ik=ikstart,ikstop
      ikevc = ik - ikstart + 1
      call davcio (evc, nwordwfc, iunwfc, ikevc, -1 )
      call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)

      ibnd1 = 0
      do ibnd=1,nbnd
         if (excluded_band(ibnd)) cycle
         ibnd1=ibnd1 + 1
         psic(:) = (0.d0, 0.d0)
         psic(nls (igk (1:npw) ) ) = evc (1:npw, ibnd)
         call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
#ifdef __PARA
         call cgather_smooth(psic,psic_all)
         write (iun_plot,'(2i5,f18.12)') ibnd1, ikevc, psic_all(1:nxxs)
#else
         write (iun_plot,'(2i5,f18.12)') ibnd1, ikevc, psic(1:nrxxs)
#endif
      end do

   end do  !ik
#ifdef __PARA
   deallocate( psic_all )
#endif
   return
end subroutine write_plot

subroutine wan2sic 

  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  use io_files, only : iunwfc, iunatsicwfc, nwordwfc, nwordwann
  USE cell_base, only : omega, tpiba2
  use gvect, only : g, ngm, ecutwfc
  use gsmooth, only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  use wavefunctions_module, only : evc, psic
  use wvfct, only : nbnd, npwx, npw, igk, g2kin
  use klist, only : nkstot, xk, wk
  use wannier

  integer :: i, j, nn, ik, ibnd, iw, ikevc 
  complex(DP), allocatable :: orbital(:,:), orb(:,:), u_matrix(:,:,:) 

  open (20, file = TRIM(seedname)//".dat" , form = 'formatted', status = 'unknown')
  write(stdout,*) ' wannier plot '

  allocate ( u_matrix( n_wannier, n_wannier, nkstot) )
  allocate ( orbital( npwx, n_wannier), orb( nrxxs, n_wannier))

  !
  do i = 1, n_wannier
     do j = 1, n_wannier
        do ik = 1, nkstot
           read (20, * ) u_matrix(i,j,ik)
           !do nn = 1, nnb(ik)
           do nn = 1, nnb
              read (20, * ) ! m_matrix (i,j,nkp,nn)
           enddo
        enddo  !nkp
     enddo !j
  enddo !i
  !
  orb(:,:) = (0.0d0,0.0d0)
  do ik=1,iknum
     ikevc = ik + ikstart - 1
     call davcio (evc, nwordwfc, iunwfc, ikevc, -1)
     call gk_sort (xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
     write(stdout,*) 'npw ',npw
     do iw=1,n_wannier
        do j=1,npw
           orbital(j,iw) = (0.0d0,0.0d0)
           do ibnd=1,n_wannier
              orbital(j,iw) = orbital(j,iw) + u_matrix(iw,ibnd,ik)*evc(j,ibnd)
              write(stdout,*) j, iw, ibnd, ik, orbital(j,iw), &
                              u_matrix(iw,ibnd,ik), evc(j,ibnd)
           enddo !ibnd
        end do  !j
     end do !wannier
     call davcio (orbital, nwordwann, iunatsicwfc, ikevc, +1)
  end do ! k-points

  deallocate ( u_matrix) 
  write(stdout,*) ' dealloc u '
  deallocate (  orbital)
  write(stdout,*) ' dealloc orbital '
  deallocate ( orb )
  write(stdout,*) ' dealloc orb '
  !
end subroutine wan2sic 

subroutine ylm_expansion 
   USE io_global,  ONLY : stdout
   use kinds, ONLY :  DP
   USE random_numbers,       ONLY : rndm
   use wannier
   implicit none
   ! local variables
   integer, parameter :: lmax2=16
   integer ::  lm, i, ir, iw, m
   real(DP) :: capel
   real(DP), allocatable :: r(:,:), rr(:), rp(:,:), ylm_w(:), ylm(:,:), mly(:,:)
   real(DP) :: u(3,3)

   allocate (r(3,lmax2), rp(3,lmax2), rr(lmax2), ylm_w(lmax2))
   allocate (ylm(lmax2,lmax2), mly(lmax2,lmax2) )

   ! generate a set of nr=lmax2 random vectors
   do ir=1,lmax2
      do i=1,3
         r(i,ir) = rndm() -0.5d0
      end do
   end do
   rr(:) = r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:)
   !- compute ylm(ir,lm)
   call ylmr2(lmax2, lmax2, r, rr, ylm)
   !- store the inverse of ylm(ir,lm) in mly(lm,ir)
   call invmat(lmax2, ylm, mly, capel)
   !- check that r points are independent
   call check_inverse(lmax2, ylm, mly)

   do iw=1, n_wannier

      !- define the u matrix that rotate the reference frame
      call set_u_matrix (xaxis(:,iw),zaxis(:,iw),u)
      !- find rotated r-vectors 
      rp(:,:) = matmul ( u(:,:) , r(:,:) )
      !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
      call ylm_wannier(ylm_w,l_w(iw),mr_w(iw),rp,lmax2) 

      csph(:,iw) = matmul (mly(:,:), ylm_w(:))

      write (stdout,*) 
      write (stdout,'(2i4,2(2x,3f6.3))') l_w(iw), mr_w(iw), xaxis(:,iw), zaxis(:,iw)
      write (stdout,'(16i6)')   (lm, lm=1,lmax2)
      write (stdout,'(16f6.3)') (csph(lm,iw), lm=1,lmax2)

   end do
   deallocate (r, rp, rr, ylm_w, ylm, mly )

   return
end subroutine ylm_expansion

subroutine check_inverse(lmax2, ylm, mly)
   use kinds, ONLY :  DP
   use constants, ONLY :  eps8
   implicit none
   ! I/O variables
   integer :: lmax2
   real(DP) :: ylm(lmax2,lmax2), mly(lmax2,lmax2)
   ! local variables
   real(DP), allocatable :: uno(:,:)
   real(DP) :: capel
   integer :: lm
   !
   allocate (uno(lmax2,lmax2) )
   uno = matmul(mly, ylm)
   capel = 0.d0
   do lm = 1, lmax2
      uno(lm,lm) = uno(lm,lm) - 1.d0
   end do
   capel = capel + SUM ( abs(uno(1:lmax2,1:lmax2) ) )
!   write (stdout,*) "capel = ", capel
   if (capel > eps8) call errore('ylm_expansion', &
                    ' inversion failed: r(*,1:nr) are not all independent !!',1)
   deallocate (uno)
   return
end subroutine check_inverse
   
subroutine set_u_matrix(x,z,u)
   use kinds, ONLY :  DP
   use constants, ONLY : eps8
   implicit none
   ! I/O variables
   real(DP) :: x(3),z(3),u(3,3)
   ! local variables
   real(DP) :: xx, zz, y(3), coseno

   xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
   if (xx < eps8) call errore ('set_u_matrix',' |xaxis| < eps ',1)
!   x(:) = x(:)/xx
   zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
   if (zz < eps8) call errore ('set_u_matrix',' |zaxis| < eps ',1)
!   z(:) = z(:)/zz

   coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
   if (abs(coseno) > eps8) call errore('set_u_matrix',' xaxis and zaxis are not orthogonal !',1)

   y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
   y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
   y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz

   u(1,:) = x(:)/xx
   u(2,:) = y(:)
   u(3,:) = z(:)/zz

!   write (stdout,'(3f10.7)') u(:,:)

   return

end subroutine set_u_matrix

subroutine ylm_wannier(ylm,l,mr,r,nr) 
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr) 
! of the spherical harmonic identified  by indices (l,mr) 
! in table 3.1 of the wannierf90 specification.
! 
! No reference to the particular ylm ordering internal to quantum-espresso
! is assumed. 
!
! If ordering in wannier90 code is changed or extended this should be the 
! only place to be modified accordingly
!
   use kinds, ONLY :  DP
   use constants, ONLY : pi, fpi, eps8
   implicit none
! I/O variables
!
   integer :: l, mr, nr
   real(DP) :: ylm(nr), r(3,nr)
!
! local variables
!
   real(DP), external :: s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy, &
                        fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(DP) :: rr, cost, phi
   integer :: ir
   real(DP) :: bs2, bs3, bs6, bs12
   bs2 = 1.d0/sqrt(2.d0)
   bs3=1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
!
   if (l > 3 .OR. l < -5 ) call errore('ylm_wannier',' l out of range ', 1)
   if (l>=0) then
      if (mr < 1 .OR. mr > 2*l+1) call errore('ylm_wannier','mr out of range' ,1)
   else
      if (mr < 1 .OR. mr > abs(l)+1 ) call errore('ylm_wannier','mr out of range',1)
   end if

   do ir=1, nr
      rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      if (rr < eps8) call errore('ylm_wannier',' rr too small ',1)

      cost =  r(3,ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      if (r(1,ir) > eps8) then
         phi = atan( r(2,ir)/r(1,ir) )
      else if (r(1,ir) < -eps8 ) then
         phi = atan( r(2,ir)/r(1,ir) ) + pi
      else
         phi = sign( pi/2.d0,r(2,ir) )
      end if

    
      if (l==0) then   ! s orbital
                    ylm(ir) = s(cost,phi)  
      end if
      if (l==1) then   ! p orbitals
         if (mr==1) ylm(ir) = p_z(cost,phi) 
         if (mr==2) ylm(ir) = px(cost,phi)
         if (mr==3) ylm(ir) = py(cost,phi)
      end if
      if (l==2) then   ! d orbitals
         if (mr==1) ylm(ir) = dz2(cost,phi)
         if (mr==2) ylm(ir) = dxz(cost,phi)
         if (mr==3) ylm(ir) = dyz(cost,phi)
         if (mr==4) ylm(ir) = dx2my2(cost,phi)
         if (mr==5) ylm(ir) = dxy(cost,phi)
      endif
      if (l==3) then   ! f orbitals
         if (mr==1) ylm(ir) = fz3(cost,phi)
         if (mr==2) ylm(ir) = fxz2(cost,phi)
         if (mr==3) ylm(ir) = fyz2(cost,phi)
         if (mr==4) ylm(ir) = fzx2my2(cost,phi)
         if (mr==5) ylm(ir) = fxyz(cost,phi)
         if (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
         if (mr==7) ylm(ir) = fy3x2my2(cost,phi)
      endif
      if (l==-1) then  !  sp hybrids
         if (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) ) 
         if (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) ) 
      end if
      if (l==-2) then  !  sp2 hybrids 
         if (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         if (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         if (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi) 
      end if
      if (l==-3) then  !  sp3 hybrids
         if (mr==1) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+p_z(cost,phi))
         if (mr==2) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-p_z(cost,phi))
         if (mr==3) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-p_z(cost,phi))
         if (mr==4) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+p_z(cost,phi))
      end if
      if (l==-4) then  !  sp3d hybrids
         if (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         if (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         if (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi) 
         if (mr==4) ylm(ir) = bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
         if (mr==5) ylm(ir) =-bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
      end if
      if (l==-5) then  ! sp3d2 hybrids
         if (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5*dx2my2(cost,phi)
         if (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5*dx2my2(cost,phi)
         if (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5*dx2my2(cost,phi)
         if (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5*dx2my2(cost,phi)
         if (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
         if (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
      end if

   end do

   return

end subroutine ylm_wannier

!======== l = 0 =====================================================================
function s(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) :: s, cost,phi
   s = 1.d0/ sqrt(fpi)
   return
end function s
!======== l = 1 =====================================================================
function p_z(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::p_z, cost,phi
   p_z =  sqrt(3.d0/fpi) * cost
   return
end function p_z
function px(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/fpi) * sint * cos(phi)
   return
end function px
function py(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/fpi) * sint * sin(phi)
   return
end function py
!======== l = 2 =====================================================================
function dz2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dz2, cost, phi
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
   return
end function dz2
function dxz(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
   return
end function dxz
function dyz(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
   return
end function dyz
function dx2my2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
   return
end function dx2my2
function dxy(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
   return
end function dxy
!======== l = 3 =====================================================================
function fz3(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fz3, cost, phi
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   return
end function fz3
function fxz2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   return
end function fxz2
function fyz2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   return
end function fyz2
function fzx2my2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   return
end function fzx2my2
function fxyz(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   return
end function fxyz
function fxx2m3y2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   return
end function fxx2m3y2
function fy3x2my2(cost,phi)
   use kinds, ONLY :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   return
end function fy3x2my2
!
!
!-----------------------------------------------------------------------
subroutine radialpart(ng, q, alfa, rvalue, lmax, radial)
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform 
  ! of the radial functions.
  !
  USE kinds,      ONLY : dp
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  !
  implicit none
  ! I/O
  integer :: ng, rvalue, lmax
  real(DP) :: q(ng), alfa, radial(ng,0:lmax)
  ! local variables
  real(DP), parameter :: xmin=-6.d0, dx=0.025d0, rmax=10.d0

  real(DP) :: rad_int, pref, x
  integer :: l, lp1, ir, ig, mesh_r
  real(DP), allocatable :: bes(:), func_r(:), r(:), rij(:), aux(:)

  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
  allocate ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r) )
  allocate ( aux(mesh_r))
  !
  !    compute the radial mesh
  !
  do ir = 1, mesh_r
     x = xmin  + DBLE (ir - 1) * dx 
     r (ir) = exp (x) / alfa
     rij (ir) = dx  * r (ir)
  enddo
  !
  if (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
  if (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * & 
                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
  if (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(3.0d0/2.0d0) * &
                     (1.d0 - 1.5d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
                                           exp(-alfa*r(ir)/3.0d0)
  pref = fpi/sqrt(omega)
  !
  do l = 0, lmax
     do ig=1,ng
       call sph_bes (mesh_r, r(1), q(ig), l, bes)
       aux(:) = bes(:) * func_r(:) * r(:)
       call simpson (mesh_r, aux, rij, rad_int)
       radial(ig,l) = rad_int * pref
     enddo
  enddo

  deallocate (bes, func_r, r, rij, aux )
  return
end subroutine radialpart


