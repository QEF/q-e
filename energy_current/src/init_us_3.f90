!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_3(npw_, xvkb_)
   !----------------------------------------------------------------------
   !
   !   Calculates xbeta functions  with
   !   structure factor, for all atoms, in reciprocal space. On input:
   !      npw_       : number of PWs
   !  On output:
   !      xvkb_(npwx,nkb,3) : beta functions (npw_ <= npwx).
   !
   USE atom, ONLY: rgrid
   USE kinds, ONLY: DP
   USE ions_base, ONLY: nat, ntyp => nsp, ityp, tau
   USE cell_base, ONLY: tpiba
   USE constants, ONLY: tpi, pi, fpi
   USE gvect, ONLY: eigts1, eigts2, eigts3, mill, g, gg, ngl, igtongl, gl, gstart
   USE wvfct, ONLY: npwx, npw, nbnd
   USE us, ONLY: nqx, dq, spline_ps, tab, tab_d2y
   ! USE splinelib
   USE uspp, ONLY: nkb, nhtol, nhtolm, indv, ap, aainit, vkb
   USE uspp_param, ONLY: upf, lmaxkb, nhm, nh
   use zero_mod, ONLY: tabr, tabr_d2y
   use hartree_mod, ONLY: ec_test  
   USE cell_base, ONLY: omega
   use mp, ONLY: mp_sum, mp_min
!
!modules for UT
!
   use splines
   USE wavefunctions, ONLY: psic
   use fft_interfaces, only: invfft, fwfft
   use fft_base, only: dffts
   use atom, ONLY: rgrid
   use cell_base, ONLY: at, alat
   use mp_world, ONLY: mpime
   use io_global, ONLY: ionode
   use mp_pools, only: intra_pool_comm
!
   implicit none
!
   !
   !local variables for SBT
   integer ::nexp, nr, npi, err
   real(DP), allocatable ::rr(:), kk(:), ya(:), yb(:), yc(:)
   real(DP) ::rhomin, rhomax, dr, fixedk, kpmin, skk, cf
!
   TYPE(spline_data) :: spl_beta
   real(DP) :: vint, vint_save, origin
   integer :: i, n1, n2, n3
   integer :: ndm
   INTEGER, INTENT(IN) :: npw_
   COMPLEX(DP), INTENT(OUT) :: xvkb_(npwx, nkb, 3)
   !
   !local variables for UT
   !
   complex(DP) :: vkbr(1:dffts%nnr), vkbr_2(1:dffts%nnr), vkbr_3(1:dffts%nnr, 3)
   real(DP)    :: stampa(nkb, 1:dffts%nnr, 3, 3), stampa_b(nkb, 1:dffts%nnr, 3, 3), ris(4, 4)
   complex(DP) :: add
   integer     :: cont_1, cont_2, cont_3, icont, a, b, ikb, ipol
   integer     :: ir, ir_found
   real(DP)    :: u(3), u_x(3), u_y(3), u_z(3), modulus, value
   integer     :: iqq, ix, iy, iz
!  integer     :: nr3s_end,nr3s_start,ii
   logical     :: l_spline
   !
   !
   !     Local variables
   !
   integer :: i0, i1, i2, i3, ig, l, lm, na, nt, nb, ih, jkb, igl, in, it, ll, lp

   real(DP) :: px, ux, vx, wx, arg, cost, xg
   real(DP), allocatable :: gk(:, :), sg(:), vq(:), ylm(:, :), betagl(:)
   real(DP), allocatable ::beta_save(:), beta_save2(:)
   complex(DP), allocatable ::xvkb1(:, :, :)
   complex(DP), allocatable :: sk(:)
   real(DP), allocatable :: aux(:)
   real(DP), allocatable :: xdata(:)
   integer :: iq

   
! testing variables
   integer, external :: find_free_unit
   integer :: iun
   integer :: ii,nr3s_end, nr3s_start, vkb_pol 
   character(len=20) :: string
   
!
   if (lmaxkb .lt. 0) return
   call start_clock('init_us_3')
   allocate (xvkb1(npw_, nhm, 3))
   allocate (sk(npw_))
   allocate (ylm(npw_, (lmaxkb + 2)**2))
   allocate (betagl(ngl))
   allocate (beta_save(ngl))
   allocate (beta_save2(ngl))
   !
   call ylmr2((lmaxkb + 2)**2, npw_, g(1:3, 1:npw_), gg(1:npw_), ylm)
   call aainit(lmaxkb + 2)
   ndm = MAXVAL(upf(:)%kkbeta)
   allocate (aux(ndm))
!!----------------------------------------------------------------------------------------I
!  !test componenti a g=0 armoniche sferiche
!  l_test=.false.
!  if (l_test) then
!     open(unit=15,file='gzero',status='unknown')
!     do lm=1,(lmaxkb+2)**2
!          write(15,"(I5,F12.7)") lm,ylm(1,lm)
!     end do
!     write(15,*) 'gl:',gl(1)
!     close(15)
!  end if
!!-----------------------------------------------------------------------------------------F
   jkb = 0
   xvkb_ = 0.d0
   cost = SQRT(4.*pi/3.)
   do it = 1, ntyp
      xvkb1(1:npw_, 1:nhm, 1:3) = (0.d0, 0.d0)
      do ih = 1, nh(it)
!gli indici (ih,it) identificano il proiettore che andiamo ad inizializzare.
!
!l è il "momento angolare+1" dello pseudo considerato
         l = nhtol(ih, it) + 1
!
!       if (ionode) print*,"ih,LLLLLLL+1=1",ih,l
!
!lp è l'indice combinato del'armonica sferica riferito allo pseudo che stiamo caricando
         lp = nhtolm(ih, it)
!
!
!       if (ionode) print*,"ih,LPPPPPPP=1",ih,lp
!
!nb è la beta function che dà la dipendenza radiale a questo pseudo
         nb = indv(ih, it)
!
!       if (ionode) print*,"ih,NBBBBBBB=1",ih,nb
!
!in-1,...,nhtol+1=l sono i momenti angolari che contribuiscono a xvkb(:,ikb), ovvero
!quelli che differiscono di 1 dal momento angolare del proiettore.
!in,...,nhtol+2=l+1 sono i "momenti angolari+1" che contribuiscono a xvkb(:,ikb).
         if (l .eq. 1) then
            in = 1
         else
            in = l - 1
         end if
!
!
!        if (ionode) print*,"ih,in=1,2",ih,in,l+1
!ll indicizza i momenti angolari +1 che contribuiscono a xvkb
         do ll = in, l + 1
!igl indicizza le shell. Qui calcoliamo il contributo:
!f_(nb,ll)(q)=\int _0 ^\infty dr r^3 beta_nb(r) j_ll(q.r), con q che varia sulle shell q=1,...,nlg
!Questo dipende solo da nb, ovvero della beta function, e da ll. Inseriamo questo in
!betagl(1:ngl), che si in questa versione si ricalcola e sovrascrive per ogni coppia (ih,ll)
!no e' una scelta ottimanel in quanto diversi ih possiedono la medesima beta function ed
!il medesimo integrale si calcola al momento piu' volte
            do igl = 1, ngl
               xg = sqrt(gl(igl))*tpiba
!xg è il modullo della igl-esima shell
               if (spline_ps) then
                  CALL errore('init_us_3', 'splines not implemented', 1)
!                   if (ll==l) then
!                       betagl(igl) = splint(xdata, tabr(:,nb,it,0), tabr_d2y(:,nb,it,0), xg)
!                   end if
!                   if (ll==l+1) then
!                       betagl(igl) = splint(xdata, tabr(:,nb,it,1), tabr_d2y(:,nb,it,1), xg)
!                   end if
! !nb: se l=1 questa condizione non è mai soddisfatta e siamo sempre in uno dei due casi precedenti
!                   if ((ll==l-1)) then
!                       betagl(igl) = splint(xdata, tabr(:,nb,it,-1), tabr_d2y(:,nb,it,-1), xg)
!                   end if
               else
                  px = xg/dq - int(xg/dq)
                  ux = 1.d0 - px
                  vx = 2.d0 - px
                  wx = 3.d0 - px
                  i0 = INT(xg/dq) + 1
                  i1 = i0 + 1
                  i2 = i0 + 2
                  i3 = i0 + 3
                  if (ll == l) then
                     betagl(igl) = tabr(i0, nb, it, 0)*ux*vx*wx/6.d0 + &
                    &tabr(i1, nb, it, 0)*px*vx*wx/2.d0 - &
                    &tabr(i2, nb, it, 0)*px*ux*wx/2.d0 + &
                    &tabr(i3, nb, it, 0)*px*ux*vx/6.d0
                  end if
                  if (ll == l + 1) then
                     betagl(igl) = tabr(i0, nb, it, 1)*ux*vx*wx/6.d0 + &
                    &tabr(i1, nb, it, 1)*px*vx*wx/2.d0 - &
                    &tabr(i2, nb, it, 1)*px*ux*wx/2.d0 + &
                    &tabr(i3, nb, it, 1)*px*ux*vx/6.d0
                  end if
!nb: se l=1 questa condizione non è mai soddisfatta e siamo sempre in uno dei due casi precedenti
                  if ((ll == l - 1)) then
                     betagl(igl) = tabr(i0, nb, it, -1)*ux*vx*wx/6.d0 + &
                    &tabr(i1, nb, it, -1)*px*vx*wx/2.d0 - &
                    &tabr(i2, nb, it, -1)*px*ux*wx/2.d0 + &
                    &tabr(i3, nb, it, -1)*px*ux*vx/6.d0
                  end if
               end if
            end do
!!------------------------------------------------------------------------------I
!     !parte test
!           l_test=.false.
!           if (l_test) then
!              if ((ih==2).and.(ll==1)) then
!                  beta_save(1:ngl)=betagl(1:ngl)
!              end if
!           end if
!!------------------------------------------------------------------------------F
!Ora dobbiamo eseguire il ciclo su tutte le armoniche sferiche con momento angolare+1= ll
            do lm = (ll - 1)**2 + 1, ll**2
               do ig = gstart, npw_
!nelle seguente espressioni compare (-i)^(ll-1) invece del solito (-i)^(ll) perchè
!il momento angolare che stiamo considerando non è ll ma ll-1. ll è solo un indice.
                  add = ylm(ig, lm)
                  xvkb1(ig, ih, 1) = xvkb1(ig, ih, 1) - &!ap(lm,lp,3)*ylm(ig,lm)
                                     (cmplx(cost*ap(lm, 3, lp)*add*betagl(igtongl(ig)))*((0.d0, -1.d0)**(ll - 1)))
!
                  xvkb1(ig, ih, 2) = xvkb1(ig, ih, 2) - &!ap(lm,lp,4)*ylm(ig,lm)
                                     (cmplx(cost*ap(lm, 4, lp)*add*betagl(igtongl(ig)))*((0.d0, -1.d0)**(ll - 1)))
!
                  xvkb1(ig, ih, 3) = xvkb1(ig, ih, 3) + &!ap(lm,lp,2)*ylm(ig,lm)
                                     (cmplx(cost*ap(lm, 2, lp)*add*betagl(igtongl(ig)))*((0.d0, -1.d0)**(ll - 1)))
!
               end do
            end do
         end do
         if (gstart == 2) then
            xvkb1(1, ih, 1:3) = 0.d0
            if ((lp == 2) .or. (lp == 3) .or. (lp == 4)) then
               do ir = 1, upf(it)%kkbeta
                  aux(ir) = upf(it)%beta(ir, nb)*rgrid(it)%r(ir)*rgrid(it)%r(ir)
               end do
               call simpson(upf(it)%kkbeta, aux, rgrid(it)%rab, vint)
               vint = vint*cost/sqrt(omega)
               if (lp == 2) xvkb1(1, ih, 3) = +vint
               if (lp == 3) xvkb1(1, ih, 1) = -vint
               if (lp == 4) xvkb1(1, ih, 2) = -vint
            end if
!!------------------------------------------------------------------------------I
!     !parte test
!           l_test=.false.
!           if (l_test) then
!              if ((ih==2).and.(ll==1)) then
!                  vint_save=vint
!              end if
!           end if
!!------------------------------------------------------------------------------F=0.d0
         end if
!!------------------------------------------------------------------------------I
!     l_test=.false.
!     if (l_test) then
!         if (ionode) then
!            open(unit=15,file='coerenza',position='append')
!            write(15,*) "IH,LP,IT: ",ih,lp,it,nhm
!            do ig=gstart,npw
!                write(15,*) ig
!                write(15,"(2F12.7)") -ylm(ig,3)*ylm(ig,lp),dble(xvkb1(ig,ih,1))
!                write(15,"(2F12.7)") -ylm(ig,4)*ylm(ig,lp),dble(xvkb1(ig,ih,2))
!                write(15,"(2F12.7)") +ylm(ig,2)*ylm(ig,lp),dble(xvkb1(ig,ih,3))
!            end do
!            close(15)
!         end if
!     end if
!!------------------------------------------------------------------------------F
      end do
!fin qui abbiamo inizializzato xvkb1(npwx,nhm,3) per il tipo di atomo che stiamo considerando.
!Andiamo ora ad inizializzare xvkb(3,npwx,nkb) per tutti gli atomi di questo
!tipo aggiungendo il fattore di struttura.
      do na = 1, nat
!finalmente inizializza xvkb. Prima carica tutti i proiettori
!di tipo 1, poi quelli di tipo 2,... quindi alla fine vkb sarà caricato così:
! atomi con tipo 1 - atomi con tipo 2 - atomi con tipo 3 - .....
         if (ityp(na) .eq. it) then
            do ig = 1, npw_
               sk(ig) = eigts1(mill(1, ig), na)* &
                        &eigts2(mill(2, ig), na)* &
                        &eigts3(mill(3, ig), na)
            end do
            do ih = 1, nh(it)
               jkb = jkb + 1
               do ig = 1, npw_
                  xvkb_(ig, jkb, 1:3) = xvkb1(ig, ih, 1:3)*sk(ig)
               end do
            end do
         endif
      end do
   end do

!!----------------------------------------------------------------------------------------I
!!test per valore nell'origine con calcolo esplicito
!    origin=0.d0
!!!!!!!!!!
!    l_test=.false.
!    if (l_test) then
!       do ig=gstart,npw
!          origin=origin+beta_save(igtongl(ig))
!       end do
!       call mp_sum(origin, intra_pool_comm)
!       origin=origin*2.d0 !perche' e' un calcolo a gamma ed abbiamo solo la meta' dei vettori G
!       origin=origin+vint
!       origin=origin*cost/(4.d0*pi)
!       origin=origin/sqrt(omega)
!       print*, "ORIGIN VALUE:", origin
!       origin=0.d0
!       do n1=-5,5
!         do n2=-5,5
!            do n3=-5,5
!               modulus=alat*sqrt(dble(n1**2+n2**2+n3**2))
!               do ir=1,rgrid(1)%mesh
!                         if (rgrid(1)%r(ir)>modulus) then
!                     ir_found=ir
!                     exit
!                  end if
!                  ir_found=ir
!                      end do
!              value=upf(1)%beta(ir_found,2)/rgrid(1)%r(ir_found)
!              origin=origin+modulus*value/sqrt(12.d0*pi)
!            end do
!         end do
!       end do
!       print*, "ORIGIN VALUE B:", origin
!    end if
!!----------------------------------------------------------------------------------------F
!!----------------------------------------------------------------------------------------I

!!test per xvkb

   if (ec_test) then
     
        if (nat>1) then
           CALL errore('init_us_3', 'Test not working with nat > 1', 1)
        end if
        if  ( sqrt(tau(1,1)**2+tau(2,1)**2+tau(3,1)**2)  > 1.E-3 ) then
           print*, sqrt(tau(1,1)**2+tau(2,1)**2+tau(3,1)**2) 
           CALL errore('init_us_3', 'Test not working if atom is not in the origin', 1)   
        end if
                   
!
!!!!!!!!!!!! First test: we print the x-vkb in real space along the three principal axes of the cell computed in two ways:
!
! 1- fourier transforming x-vkb in real space
! 2- fourier transforming vkb in real space a manually multiplying by x
!
! The test produces different files of name "total_axes_x/y/z_ikb" where axes_x/y/z runs thorugh the three principal directions (axes) of the cell, which
! is supposed to be cubic.
!
! Each file contains 6 records. Three for each x-vbk polarizazion (xvkb, yvkab and zvkab) calculated in way (1) and other three for method 2.
!
        do vkb_pol=1,3
         do ikb=1,nkb   
         
            cont_1=0
            cont_2=0
            cont_3=0

            print*,"Comincia il test per test", vkb_pol
            print*,"Testiamo: ", ikb
            print*,"CONTROLLO",nkb,nh(1),nbnd,nkb
            print*,"CONTROLLO_POS", tau(:,:)*alat

            ! vkbr contains xvkb_ in real space
            psic=0.d0
            psic(dffts%nl(1:npw))=xvkb_(1:npw,ikb,vkb_pol)
            psic(dffts%nlm(1:npw))=CONJG(xvkb_(1:npw,ikb,vkb_pol))
            call invfft ('Wave', psic, dffts)
            vkbr(1:dffts%nnr)=1/sqrt(omega)*psic(1:dffts%nnr)

            psic=0.d0
            psic(dffts%nl(1:npw))=vkb(1:npw,ikb)
            psic(dffts%nlm(1:npw))=CONJG(vkb(1:npw,ikb))
            call invfft ('Wave', psic, dffts)
            vkbr_2(1:dffts%nnr)=1/sqrt(omega)*psic(1:dffts%nnr)
!
            nr3s_start=0
            nr3s_end =0
            do ii=1,mpime + 1
               nr3s_start=nr3s_end+1
               nr3s_end=nr3s_end+dffts%nr3p(ii)
            end do
            do iz=1,dffts%nr3p(mpime+1)     !primo ciclo   sui punti x
               do iy=1,dffts%nr2           !secondo ciclo sui punti x
                  do ix=1,dffts%nr1        !terzo ciclo   sui punti x
                     iqq=(iz-1)*(dffts%nr1x*dffts%nr2x)+(iy-1)*dffts%nr1+ix
                     u_x(1:3)=real(ix-1)/real(dffts%nr1)*at(1:3,1)*alat
                     u_y(1:3)=real(iy-1)/real(dffts%nr2)*at(1:3,2)*alat
                     u_z(1:3)=real(iz+nr3s_start-1-1)/real(dffts%nr3)*at(1:3,3)*alat
                     u(1:3) = u_x(1:3)+u_y(1:3)+u_z(1:3)
                      modulus=sqrt(u(1)**2+u(2)**2+u(3)**2)
!inizializza stampa
                     if ((iz==1).and.(iy==1)) then
                          cont_1=cont_1+1
                          stampa(ikb, cont_1,1, vkb_pol)=dble(vkbr(iqq))
                          stampa_b(ikb, cont_1,1, vkb_pol)=u(vkb_pol)*dble(vkbr_2(iqq))
                     end if
                     if ((ix==1).and.(iz==1)) then!
                          cont_2=cont_2+1
                          stampa(ikb,cont_2,2, vkb_pol)=dble(vkbr(iqq))
                          stampa_b(ikb,cont_2,2, vkb_pol)=u(vkb_pol)*dble(vkbr_2(iqq))
                     end if!
                     if ((ix==1).and.(iy==1)) then
                          cont_3=cont_3+1
                          stampa(ikb, cont_3,3, vkb_pol)=dble(vkbr(iqq))
                          stampa_b(ikb, cont_3,3, vkb_pol)=u(vkb_pol)*dble(vkbr_2(iqq))
                     end if
!!fine
!!                      do ir=1,rgrid(1)%mesh
!!                                 if (rgrid(1)%r(ir)>modulus) then
!!                            ir_found=ir
!!                            exit
!!                         end if
!!                         ir_found=ir
!!                                    end do
!!                     print*,"IRFOUND",ir_found
!!                      value=u(1)/sqrt(4.d0*pi)*upf(1)%beta(ir_found,1)/rgrid(1)%r(ir_found)
!!                      if ((iz==1).and.(ix==1)) then
!!                         if (ionode) then
!!
!!                            open(unit=15,file='plotta_b',access='append')
!!
!! !                           write(15,"(6F12.7,I5)") u(1),u(2),u(3),value,&
!! !&upf(1)%beta(ir_found,1),rgrid(1)%r(ir_found),ir_found
!!
!!                             write(15,"(2F12.7)") u(2),u(2)*dble(vkbr_2(iqq))
!!
!!                            close(15)
!!
!!                            open(unit=15,file='plotta_c',access='append')
!!
!! !                           write(15,"(6F12.7,I5)") u(1),u(2),u(3),value,&
!! !&upf(1)%beta(ir_found,1),rgrid(1)%r(ir_found),ir_found
!!
!!                             write(15,"(2F12.7)") u(2),dble(vkbr_2(iqq))
!!
!!                            close(15)
!!
!!                         end if
!!                      end if
                  end do
               end do
            end do
            
         end do   
        end do     
     
        if (ionode) then
            do ikb=1,nkb
               iun = find_free_unit()
               
               write(string,'(I0)') ikb
               open(unit=iun,file='total_axes_x_'//trim(string),status='unknown')
               do icont=1,cont_1
                   write(iun,"(6F12.7)") stampa(ikb,icont,1, 1),stampa(ikb,icont,1, 2), stampa(ikb,icont,1,3),&
                   stampa_b(ikb,icont,1, 1), stampa_b(ikb,icont,1, 2),stampa_b(ikb,icont,1,3) 
               end do
               close(iun)
               
               write(string,'(I0)') ikb
               open(unit=iun,file='total_axes_y_'//trim(string),status='unknown')
               do icont=1,cont_1
                   write(iun,"(6F12.7)") stampa(ikb,icont,2, 1),stampa(ikb,icont,2, 2), stampa(ikb,icont,2,3),&
                   stampa_b(ikb,icont,2,1), stampa_b(ikb,icont,2,2),stampa_b(ikb,icont,2,3) 
               end do
               close(iun)
               
               write(string,'(I0)') ikb
               open(unit=iun,file='total_axes_z_'//trim(string),status='unknown')
               do icont=1,cont_1
                   write(iun,"(6F12.7)") stampa(ikb,icont,3, 1),stampa(ikb,icont,3,2), stampa(ikb,icont,3,3),&
                   stampa_b(ikb,icont,3,1), stampa_b(ikb,icont,3,2),stampa_b(ikb,icont,3,3) 
               end do
               close(iun)
               
            end do   
        end if
     
     
!!test per xvkb: We evaluate xvkb by transforming vkb in real space and multiplying by x. The zero current should be only slightly changed.

        xvkb_=0.d0
        do ikb=1,nkb
          psic=0.d0
          psic(dffts%nl(1:npw))=vkb(1:npw,ikb)
          psic(dffts%nlm(1:npw))=CONJG(vkb(1:npw,ikb))
          call invfft ('Wave', psic, dffts)
          vkbr_2(1:dffts%nnr)=psic(1:dffts%nnr)

          !vkbr_2 contains the 3D beta function in real space beta(r)

          nr3s_start=0
          nr3s_end =0
          do ii=1,mpime + 1
             nr3s_start=nr3s_end+1
             nr3s_end=nr3s_end+dffts%nr3p(ii)
          end do
          do iz=1,dffts%nr3p(mpime+1)     !primo ciclo   sui punti x
             do iy=1,dffts%nr2           !secondo ciclo sui punti x
                do ix=1,dffts%nr1        !terzo ciclo   sui punti x
                   iqq=(iz-1)*(dffts%nr1x*dffts%nr2x)+(iy-1)*dffts%nr1+ix
                   u_x(1:3)=real(ix-1)/real(dffts%nr1)*at(1:3,1)*alat
                   u_y(1:3)=real(iy-1)/real(dffts%nr2)*at(1:3,2)*alat
                   u_z(1:3)=real(iz+nr3s_start-1-1)/real(dffts%nr3)*at(1:3,3)*alat
                   u(1:3) = u_x(1:3)+u_y(1:3)+u_z(1:3)

                   !using vkbr_2 we muliply bx x in real space: x*beta(r) We multiply by x. Note that we have to be carefully of pbc.  
                   do ipol=1,3
                      if (u(ipol).le.(alat/2.d0)) then
                          vkbr_3(iqq,ipol)=vkbr_2(iqq)*u(ipol)
                      else
                          vkbr_3(iqq,ipol)=vkbr_2(iqq)*(-alat+u(ipol))
                      end if
                   end do
                end do
             end do
          end do

          !We transform x*beta(r) back to reciprocal space.
          do ipol=1,3
             psic=0.d0
             psic(1:dffts%nnr)=vkbr_3(1:dffts%nnr,ipol)
             call fwfft ('Wave', psic, dffts)
             xvkb_(1:npw,ikb,ipol)=psic(dffts%nl(1:npw))
          end do
        end do
   end if
!!----------------------------------------------------------------------------------------F
!!----------------------------------------------------------------------------------------I
!!test per le armoniche sferiche
!     l_test=.false.
!     if (l_test) then
!        open(unit=15,file='spherical_harmonics',status='unknown')
!        do ig=gstart,npw
!                write(15,"(A,I5,2F12.7)") "1",ig,ylm(ig,1),sqrt(1.d0/(4.d0*pi))
!                write(15,"(A,I5,2F12.7)") "2",ig,ylm(ig,2),sqrt(3.d0/(4.d0*pi))*g(3,ig)/sqrt(gg(ig))
!                write(15,"(A,I5,2F12.7)") "3",ig,ylm(ig,3),sqrt(3.d0/(4.d0*pi))*(-1.d0)*g(1,ig)/sqrt(gg(ig))
!                write(15,"(A,I5,2F12.7)") "4",ig,ylm(ig,4),sqrt(3.d0/(4.d0*pi))*(-1.d0)*g(2,ig)/sqrt(gg(ig))
!                write(15,"(A,I5,2F12.7)") "5",ig,ylm(ig,5),sqrt(5.d0/(16.d0*pi))&
!*(2.d0*(g(3,ig)**2)-g(1,ig)**2-g(2,ig)**2)/gg(ig)
!                write(15,"(A,I5,2F12.7)") "6",ig,ylm(ig,6),sqrt(15.d0/(4.d0*pi))*(-1.d0)*g(3,ig)*g(1,ig)/gg(ig)
!                write(15,"(A,I5,2F12.7)") "7",ig,ylm(ig,7),sqrt(15.d0/(4.d0*pi))*(-1.d0)*g(3,ig)*g(2,ig)/gg(ig)
!                write(15,"(A,I5,2F12.7)") "8",ig,ylm(ig,8),sqrt(15.d0/(16.d0*pi))*(+1.d0)*(g(1,ig)**2-g(2,ig)**2)/gg(ig)
!                write(15,"(A,I5,2F12.7)") "9",ig,ylm(ig,9),sqrt(15.d0/(4.d0*pi))*(+1.d0)*g(1,ig)*g(2,ig)/gg(ig)
!        end do
!         close(15)
!     end if
!!----------------------------------------------------------------------------------------F
!!----------------------------------------------------------------------------------------I
!!test per ap
!     l_test=.false.
!     if (l_test) then
!        open(unit=15,file='clebsch_c',status='unknown')
!         do i=1,25
!            write(15,"(' I2-3-4  1 ' ,I5,3F7.4)") i ,ap(i,2,1), ap(i,3,1), ap(i,4,1)
!         enddo
!         do i=1,25
!            write(15,"(' I2-3-4  2 ' ,I5,3F7.4)") i ,ap(i,2,2), ap(i,3,2), ap(i,4,2)
!         enddo
!         do i=1,25
!            write(15,"(' I2-3-4  3 ' ,I5,3F7.4)")  i ,ap(i,2,3), ap(i,3,3), ap(i,4,3)
!         enddo
!         do i=1,25
!            write(15,"(' I2-3-4  4 ' ,I5,3F7.4)")  i ,ap(i,2,4), ap(i,3,4), ap(i,4,4)
!         enddo
!        do a=4,4
!            do b=4,4
!               do ig=1,npw
!                  ris(a,b)=0
!                  do lm=1,9
!                     ris(a,b)=ris(a,b)+ap(lm,a,b)*ylm(ig,lm)
!!                      if (ig==5) &
!!                      write(15,"(3I5,F12.7)") lm,a,b,ap(lm,a,b)
!                  end do
!!                  if (ig==5) &
!!                     write(15,"(3I5,2F12.7)") ig,a,b,ris(a,b),ylm(ig,a)*ylm(ig,b)
!               end do
!            end do
!        end do
!        close(15)
!     end if
!!----------------------------------------------------------------------------------------F

   deallocate (ylm)
   deallocate (sk)
   deallocate (xvkb1)
   deallocate (betagl)
   deallocate (aux)

   call stop_clock('init_us_3')
   return
end subroutine init_us_3

