!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
module init_us_3_mod
contains
   subroutine init_us_3(npw_, xvkb_, tabr, ec_test)
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
      USE ions_base, ONLY: nat, ntyp => nsp, ityp
      USE cell_base, ONLY: tpiba, omega
      USE constants, ONLY: tpi, pi, fpi
      USE gvect, ONLY: eigts1, eigts2, eigts3, mill, g, gg, ngl, igtongl, gl, gstart
      USE wvfct, ONLY: npwx
      use uspp_data, ONLY: dq
      ! USE splinelib
      USE uspp, ONLY: nkb, nhtol, nhtolm, indv, ap, aainit
      USE uspp_param, ONLY: upf, lmaxkb, nhm, nh
      use mp, ONLY: mp_sum, mp_min
!
!modules for UT
!
      !use splines
      use atom, ONLY: rgrid
!
      implicit none
      !TYPE(spline_data) :: spl_beta
      real(DP) :: vint
      integer :: ndm
      INTEGER, INTENT(IN) :: npw_
      COMPLEX(DP), INTENT(OUT) :: xvkb_(npwx, nkb, 3)
      logical, intent(in) :: ec_test
      real(dp), intent(in) :: tabr(:, :, :, :)
      !
      !local variables for UT
      !
      complex(DP) :: add
      integer     :: ir
      !
      !     Local variables
      !
      integer :: i0, i1, i2, i3, ig, l, lm, na, nb, ih, jkb, igl, in, it, ll, lp
      real(DP) :: px, ux, vx, wx, cost, xg
      real(DP), allocatable :: ylm(:, :), betagl(:)
      complex(DP), allocatable ::xvkb1(:, :, :)
      complex(DP), allocatable :: sk(:)
      real(DP), allocatable :: aux(:)
      integer :: ii
      !
      if (lmaxkb .lt. 0) return
      call start_clock('init_us_3')
      allocate (xvkb1(npw_, nhm, 3))
      allocate (sk(npw_))
      allocate (ylm(npw_, (lmaxkb + 2)**2))
      allocate (betagl(ngl))
      ! this can be moved and calculated once per trajectory
      call ylmr2((lmaxkb + 2)**2, npw_, g(1:3, 1:npw_), gg(1:npw_), ylm)
      call aainit(lmaxkb + 2)
      ndm = MAXVAL(upf(:)%kkbeta)
      allocate (aux(ndm))
      jkb = 0
      xvkb_ = 0.d0
      cost = SQRT(4.*pi/3.)
      do it = 1, ntyp
         xvkb1(1:npw_, 1:nhm, 1:3) = (0.d0, 0.d0)
         do ih = 1, nh(it)
!
!indexes (ih,it) identify the projector taken into consideration.
!
!l il the  "angular momentum + 1" of the projector considered.
            l = nhtol(ih, it) + 1
!lp is the combined index of the spherical harmonics of the considered projector,
            lp = nhtolm(ih, it)
!nb is the beta function input of the UPF providing the radial dependence to the projector
            nb = indv(ih, it)
! in-1,...,nhtol+1=l are the angular momenta that contribute to xvkb(:,ikb).
! These corresponds to the range of L in Eq. 41 where the Clebsch-Gordan coefficients are different from zero.
! and differ by at max 1 from the angular momentum of the projector.
! in,...,nhtol+2=l+1 are the "angular momenta + 1" contributing to xvkb(:,ikb).
            if (l .eq. 1) then
               in = 1
            else
               in = l - 1
            end if
!ll indexes the angular moment contributing to xvkb (see previous comments)
            do ll = in, l + 1
!igl indexes the shells of G vectors. Here we evaluate :
!
! f_(nb,ll)(q)=\int _0 ^\infty dr r^3 beta_nb(r) j_ll(q.r), with q varying over the shells q=1,...,nlg [1]
!
! This integral depends only on nb, i.e. on the beta function, and on ll (L in Eq. 41)
! We save its value in betagl(1:ngl), that is here recomputed and overwrite for each pair (ih,ll)
! This is not an optimal choice because various ih correspond to the same beta function and
! the same integral is at the moment recomputed many times. Optimization TBD.
!
! Note that the integrals [1] have already been saved in tabr on a grid of q values. Here we just interpolate.
!
               do igl = 1, ngl
                  xg = sqrt(gl(igl))*tpiba
!xg is the modulus of the igl-shell
                     px = xg/dq - int(xg/dq)
                     ux = 1.d0 - px
                     vx = 2.d0 - px
                     wx = 3.d0 - px
                     i0 = INT(xg/dq) + 1
                     i1 = i0 + 1
                     i2 = i0 + 2
                     i3 = i0 + 3
                     ii = ll - l + 2
                     if (ii < 1 .or. ii > 3) then
                        call ERRORE('init_us_3', 'Internal index error', 1)
                     end if
                     betagl(igl) = tabr(i0, nb, it, ii)*ux*vx*wx/6.d0 + &
                       &tabr(i1, nb, it, ii)*px*vx*wx/2.d0 - &
                       &tabr(i2, nb, it, ii)*px*ux*wx/2.d0 + &
                       &tabr(i3, nb, it, ii)*px*ux*vx/6.d0
               end do
!
! Now we have to cycle over all spherical harmonics with angular momentum + 1 = ll
               do lm = (ll - 1)**2 + 1, ll**2
                  do ig = gstart, npw_
! It the following expressions we have (-i)^(ll-1) instead of (-i)^(ll) because
! The angular momentum considered is not ll but ll-1. ll is just an index corresponding to (angular momentum + 1).
! See also previous comments.
                     add = ylm(ig, lm)
                     xvkb1(ig, ih, 1) = xvkb1(ig, ih, 1) - &!ap(lm,lp,3)*ylm(ig,lm)
                                        (cmplx(cost*ap(lm, 3, lp)*add*betagl(igtongl(ig)))*((0.d0, -1.d0)**(ll - 1)))
                     xvkb1(ig, ih, 2) = xvkb1(ig, ih, 2) - &!ap(lm,lp,4)*ylm(ig,lm)
                                        (cmplx(cost*ap(lm, 4, lp)*add*betagl(igtongl(ig)))*((0.d0, -1.d0)**(ll - 1)))
                     xvkb1(ig, ih, 3) = xvkb1(ig, ih, 3) + &!ap(lm,lp,2)*ylm(ig,lm)
                                        (cmplx(cost*ap(lm, 2, lp)*add*betagl(igtongl(ig)))*((0.d0, -1.d0)**(ll - 1)))
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
            end if
         end do
         ! Up to here we initialized xvkb1(npwx,nhm,3) for the type of atom considered.
         ! We now initialize xvkb(3,npwx,nkb) for all atoms of this type adding the structur factors.
         do na = 1, nat
            !We first load all projectors for atoms of type 1, than type 2,...
            !At the end xvkb will be loaded like this:
            !atoms of type 1 - atoms of type 2 - atoms of type 3 - .....
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
      !!test for xvkb (if needed, only for development)
      if (ec_test) then
         call init_us_3_test(npw_, xvkb_)
      end if
      deallocate (ylm)
      deallocate (sk)
      deallocate (xvkb1)
      deallocate (betagl)
      deallocate (aux)
      call stop_clock('init_us_3')
      return
   end subroutine init_us_3

   subroutine init_us_3_test(npw_, xvkb_)!, rgrid, ntyp, ityp, tau, tpiba, )

      ! Subroutine for tests. Only for development.
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
      USE wvfct, ONLY: npwx, npw, nbnd
      ! USE splinelib
      USE uspp, ONLY: nkb, vkb
      USE uspp_param, ONLY: upf, lmaxkb, nhm, nh
      USE cell_base, ONLY: omega
      use mp, ONLY: mp_sum, mp_min
!
!modules for UT
!
      !use splines
      USE wavefunctions, ONLY: psic
      use fft_interfaces, only: invfft, fwfft
      use fft_base, only: dffts
      use cell_base, ONLY: at, alat
      use mp_world, ONLY: mpime
      use io_global, ONLY: ionode
      use mp_pools, only: intra_pool_comm
!
      implicit none
!
      !TYPE(spline_data) :: spl_beta
      INTEGER, INTENT(IN) :: npw_
      COMPLEX(DP), INTENT(inOUT) :: xvkb_(npwx, nkb, 3)
      !
      !local variables for UT
      !
      integer     :: cont_1, cont_2, cont_3, icont, ikb, ipol
      real(DP)    :: u(3), u_x(3), u_y(3), u_z(3), modulus!, value
      integer     :: iqq, ix, iy, iz

! testing variables used only for ec_test

      complex(DP) :: vkbr(1:dffts%nnr), vkbr_2(1:dffts%nnr), vkbr_3(1:dffts%nnr, 3)
      real(DP)    :: stampa(nkb, 1:dffts%nnr, 3, 3), stampa_b(nkb, 1:dffts%nnr, 3, 3) !, ris(4, 4)
      integer, external :: find_free_unit
      integer :: iun, iun2, iun3
      integer :: ii, nr3s_end, nr3s_start, vkb_pol
      character(len=20) :: string, dimension_string

      if (nat > 1) then
         CALL errore('init_us_3', 'Test not working with nat > 1', 1)
      end if
      if (sqrt(tau(1, 1)**2 + tau(2, 1)**2 + tau(3, 1)**2) > 1.E-3) then
         print *, sqrt(tau(1, 1)**2 + tau(2, 1)**2 + tau(3, 1)**2)
         CALL errore('init_us_3', 'Test not working if atom is not in the origin', 1)
      end if

!
!!!!!!!!!!!! First test: we print the x-vkb in real space along the three principal axes of the cell computed in two ways:
!
! (1) fourier transforming xvkb_ in real space
! (2) fourier transforming vkb in real space a manually multiplying by x or y or z
!
! The test produces different files of name "total_axes_x/y/z_ikb" where axes_x/y/z runs thorugh the three principal directions (axes) of the cell, which
! is supposed to be cubic.
!
! Each file contains 6 records. Three for each xvbk polarizazion (xvkb, yvkb and zvkb) calculated in way (1) and other three for method (2).
!

      write (dimension_string, '(I0," ",I0, " ", I0)') dffts%nr3p(mpime + 1), dffts%nr2, dffts%nr1
      iun3 = find_free_unit()
      open (unit=iun3, file='full_unformatted.xmf')
      write (iun3, '(a)') &
  '<?xml version="1.0" ?><!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>&
  &<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2"><Domain>&
  &    <Grid GridType="Uniform">&
  &      <Topology TopologyType="3DCORECTMesh" Dimensions="' &
  //trim(dimension_string)//'"/>&
  &      <Geometry GeometryType="ORIGIN_DXDYDZ">&
  &        <DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="8" Format="XML">&
  &          0 0 0&
  &        </DataItem>&
  &        <DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="8" Format="XML">&
  &          1 1 1&
  &        </DataItem>&
  &      </Geometry>'

      do vkb_pol = 1, 3
         do ikb = 1, nkb

            cont_1 = 0
            cont_2 = 0
            cont_3 = 0

            print *, "Comincia il test per test", vkb_pol
            print *, "Testiamo: ", ikb
            print *, "CONTROLLO", nkb, nh(1), nbnd, nkb
            print *, "CONTROLLO_POS", tau(:, :)*alat

            ! vkbr contains xvkb_ in real space
            psic = 0.d0
            psic(dffts%nl(1:npw)) = xvkb_(1:npw, ikb, vkb_pol)
            psic(dffts%nlm(1:npw)) = CONJG(xvkb_(1:npw, ikb, vkb_pol))
            call invfft('Wave', psic, dffts)
            vkbr(1:dffts%nnr) = 1/sqrt(omega)*psic(1:dffts%nnr)

            ! vkbr_2 constains vkb in real space
            psic = 0.d0
            psic(dffts%nl(1:npw)) = vkb(1:npw, ikb)
            psic(dffts%nlm(1:npw)) = CONJG(vkb(1:npw, ikb))
            call invfft('Wave', psic, dffts)
            vkbr_2(1:dffts%nnr) = 1/sqrt(omega)*psic(1:dffts%nnr)
!
            nr3s_start = 0
            nr3s_end = 0
            do ii = 1, mpime + 1
               nr3s_start = nr3s_end + 1
               nr3s_end = nr3s_end + dffts%nr3p(ii)
            end do
            iun = find_free_unit()
            write (string, '(I0,"_",I0)') ikb, vkb_pol
            write (iun3, *) '<Attribute Name="A_'//trim(string)//'" Active="1" AttributeType="Scalar" Center="Cell">&
&              <DataItem Dimensions="', dimension_string, '" NumberType="Float" Precision="8"&
&              Format="Binary">full_unformatted_A_', string, '</DataItem>&
&      </Attribute>'
            write (iun3, *) '<Attribute Name="B_'//trim(string)//'" Active="1" AttributeType="Scalar" Center="Cell">&
&              <DataItem Dimensions="', dimension_string, '" NumberType="Float" Precision="8"&
&              Format="Binary">full_unformatted_B_'//trim(string)//'</DataItem>&
&      </Attribute>'

            open (unit=iun, file='full_unformatted_A_'//trim(string), FORM="unformatted", access='stream')
            iun2 = find_free_unit()
            open (unit=iun2, file='full_unformatted_B_'//trim(string), FORM="unformatted", access='stream')
            do iz = 1, dffts%nr3p(mpime + 1)     !primo ciclo   sui punti x
               do iy = 1, dffts%nr2           !secondo ciclo sui punti x
                  do ix = 1, dffts%nr1        !terzo ciclo   sui punti x
                     iqq = (iz - 1)*(dffts%nr1x*dffts%nr2x) + (iy - 1)*dffts%nr1 + ix
                     u_x(1:3) = real(ix - 1)/real(dffts%nr1)*at(1:3, 1)*alat
                     u_y(1:3) = real(iy - 1)/real(dffts%nr2)*at(1:3, 2)*alat
                     u_z(1:3) = real(iz + nr3s_start - 1 - 1)/real(dffts%nr3)*at(1:3, 3)*alat
                     u(1:3) = u_x(1:3) + u_y(1:3) + u_z(1:3)
                     modulus = sqrt(u(1)**2 + u(2)**2 + u(3)**2)
                     write (iun) dble(vkbr(iqq))
                     write (iun2) u(vkb_pol)*dble(vkbr_2(iqq))

                     !init "stampa" variable
                     if ((iz == 1) .and. (iy == 1)) then
                        cont_1 = cont_1 + 1
                        stampa(ikb, cont_1, 1, vkb_pol) = dble(vkbr(iqq))
                        stampa_b(ikb, cont_1, 1, vkb_pol) = u(vkb_pol)*dble(vkbr_2(iqq))
                     end if
                     if ((ix == 1) .and. (iz == 1)) then!
                        cont_2 = cont_2 + 1
                        stampa(ikb, cont_2, 2, vkb_pol) = dble(vkbr(iqq))
                        stampa_b(ikb, cont_2, 2, vkb_pol) = u(vkb_pol)*dble(vkbr_2(iqq))
                     end if!
                     if ((ix == 1) .and. (iy == 1)) then
                        cont_3 = cont_3 + 1
                        stampa(ikb, cont_3, 3, vkb_pol) = dble(vkbr(iqq))
                        stampa_b(ikb, cont_3, 3, vkb_pol) = u(vkb_pol)*dble(vkbr_2(iqq))
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
            close (iun)
            close (iun2)
         end do
      end do
      write (iun3, *) '</Grid></Domain></Xdmf>'
      close (iun3)

      if (ionode) then
         do ikb = 1, nkb
            iun = find_free_unit()

            write (string, '(I0)') ikb
            open (unit=iun, file='total_axes_x_'//trim(string), status='unknown')
            do icont = 1, cont_1
               write (iun, "(6F12.7)") stampa(ikb, icont, 1, 1), stampa(ikb, icont, 1, 2), stampa(ikb, icont, 1, 3), &
                  stampa_b(ikb, icont, 1, 1), stampa_b(ikb, icont, 1, 2), stampa_b(ikb, icont, 1, 3)
            end do
            close (iun)

            write (string, '(I0)') ikb
            open (unit=iun, file='total_axes_y_'//trim(string), status='unknown')
            do icont = 1, cont_1
               write (iun, "(6F12.7)") stampa(ikb, icont, 2, 1), stampa(ikb, icont, 2, 2), stampa(ikb, icont, 2, 3), &
                  stampa_b(ikb, icont, 2, 1), stampa_b(ikb, icont, 2, 2), stampa_b(ikb, icont, 2, 3)
            end do
            close (iun)

            write (string, '(I0)') ikb
            open (unit=iun, file='total_axes_z_'//trim(string), status='unknown')
            do icont = 1, cont_1
               write (iun, "(6F12.7)") stampa(ikb, icont, 3, 1), stampa(ikb, icont, 3, 2), stampa(ikb, icont, 3, 3), &
                  stampa_b(ikb, icont, 3, 1), stampa_b(ikb, icont, 3, 2), stampa_b(ikb, icont, 3, 3)
            end do
            close (iun)

         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Second est per xvkb: We evaluate xvkb by transforming vkb in real space and multiplying by x. In this test
! the xvkb_ is updated. The zero current should be only slightly changed w.r.t. a calculation done with ec_test = .false.
!

      xvkb_ = 0.d0
      do ikb = 1, nkb
         psic = 0.d0
         psic(dffts%nl(1:npw)) = vkb(1:npw, ikb)
         psic(dffts%nlm(1:npw)) = CONJG(vkb(1:npw, ikb))
         call invfft('Wave', psic, dffts)
         vkbr_2(1:dffts%nnr) = psic(1:dffts%nnr)

         !vkbr_2 contains the 3D beta function in real space beta(r)

         nr3s_start = 0
         nr3s_end = 0
         do ii = 1, mpime + 1
            nr3s_start = nr3s_end + 1
            nr3s_end = nr3s_end + dffts%nr3p(ii)
         end do
         do iz = 1, dffts%nr3p(mpime + 1)     !primo ciclo   sui punti x
            do iy = 1, dffts%nr2           !secondo ciclo sui punti x
               do ix = 1, dffts%nr1        !terzo ciclo   sui punti x
                  iqq = (iz - 1)*(dffts%nr1x*dffts%nr2x) + (iy - 1)*dffts%nr1 + ix
                  u_x(1:3) = real(ix - 1)/real(dffts%nr1)*at(1:3, 1)*alat
                  u_y(1:3) = real(iy - 1)/real(dffts%nr2)*at(1:3, 2)*alat
                  u_z(1:3) = real(iz + nr3s_start - 1 - 1)/real(dffts%nr3)*at(1:3, 3)*alat
                  u(1:3) = u_x(1:3) + u_y(1:3) + u_z(1:3)

                  !using vkbr_2 we muliply bx x in real space: x*beta(r) We multiply by x. Note that we have to be carefully of pbc.
                  do ipol = 1, 3
                     if (u(ipol) .le. (alat/2.d0)) then
                        vkbr_3(iqq, ipol) = vkbr_2(iqq)*u(ipol)
                     else
                        vkbr_3(iqq, ipol) = vkbr_2(iqq)*(-alat + u(ipol))
                     end if
                  end do
               end do
            end do
         end do

         !We transform x*beta(r) back to reciprocal space.
         do ipol = 1, 3
            psic = 0.d0
            psic(1:dffts%nnr) = vkbr_3(1:dffts%nnr, ipol)
            call fwfft('Wave', psic, dffts)
            xvkb_(1:npw, ikb, ipol) = psic(dffts%nl(1:npw))
         end do
      end do
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

   end subroutine
end module
