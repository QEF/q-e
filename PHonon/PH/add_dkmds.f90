!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine add_dkmds(ik, uact, jpol, dvkb)
  !--------=========-------------------------------------------------------
  !
  ! This subroutine adds to dvpsi the terms which depend on the augmentation
  ! charge. It assumes that the variable dpqq, has been set. In the noncollinear
  ! and spin_orbit case the variable dpqq_so must be set.
  ! NB: I think this routine is called only for q=0; case q/=0 not implemented
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, tpiba
  USE gvect, ONLY : g
  USE lsda_mod, ONLY: lsda, current_spin, isk, nspin
  USE klist, ONLY : xk, ngk, igk_k
  USE spin_orb, ONLY : lspinorb
  USE uspp, ONLY : nkb, qq, qq_so, vkb
  USE wvfct, ONLY : npwx, nbnd
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module, ONLY : noncolin, npol
  USE wavefunctions_module,    ONLY : evc
  USE uspp_param, only: nh
  USE becmod, ONLY: calbec
  USE phus,   ONLY : alphap

  USE qpoint,     ONLY : ikks, ikqs
  USE lrus,       ONLY : becp1, dpqq, dpqq_so
  USE eqv,        ONLY : dvpsi
  USE control_lr, ONLY : nbnd_occ

  implicit none

  integer, intent(in) :: ik, jpol
  complex(DP), intent(in) :: uact (3 * nat)
  complex(DP), intent(in) :: dvkb (npwx,nkb,3)


  real(DP), parameter :: eps = 1.d-12

  integer :: npw, npwq, ipol, ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, ig, igg, mu, ikq

  logical :: ok

  complex(DP), allocatable :: ps1(:,:), ps2(:,:,:)
  complex(DP), allocatable :: ps1_nc(:,:,:), ps2_nc(:,:,:,:)
  complex(DP), allocatable :: alphadk(:,:,:), becp2(:,:)
  complex(DP), allocatable :: alphadk_nc(:,:,:,:), becp2_nc(:,:,:)
  complex(DP), allocatable :: aux(:), aux1(:,:)


  integer :: i,j,is

#ifdef TIMING_ADD_DKMDS
  call start_clock('add_dkmds')
  call start_clock('add_dkmds2')
#endif
  allocate(aux(npwx))
  allocate(aux1(npwx*npol,nbnd))
  if (nkb.gt.0) then
     if (noncolin) then
        allocate (ps1_nc(nkb,npol,nbnd))
        allocate (ps2_nc(nkb,npol,3,nbnd))
        allocate (alphadk_nc(nkb,npol,nbnd,3))
        allocate (becp2_nc(nkb,npol,nbnd))
     else
        allocate (ps1(nkb,nbnd))
        allocate (ps2(nkb,3,nbnd))
        allocate (alphadk(nkb,nbnd,3))
        allocate (becp2(nkb,nbnd))
     end if
  end if

  if (noncolin) then
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
  else
     ps1 = (0.d0, 0.d0)
     ps2 = (0.d0, 0.d0)
  endif
  !
  !   First we calculate the alphadk = <d/dk d/du beta|psi>
  !   and becp2 = < d/dk beta | psi>
  !
  if (lsda) current_spin = isk (ik)
  npw = ngk(ik)
  ikq = ik
  npwq= ngk(ik)
  if (noncolin) then
     call calbec (npw, dvkb(:,:,jpol), evc, becp2_nc)
  else
     call calbec (npw, dvkb(:,:,jpol), evc, becp2)
  endif

#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds2')
  call start_clock('add_dkmds3')
#endif

  do ipol = 1, 3
     do ibnd = 1, nbnd
        do ig = 1, npw
           aux1 (ig, ibnd) = evc(ig,ibnd) * tpiba * (0.d0,1.d0) * &
                ( xk(ipol,ik) + g(ipol,igk_k(ig,ik)) )
        enddo
        if (noncolin) then
           do ig = 1, npw
              aux1 (ig+npwx, ibnd) = evc(ig+npwx,ibnd)*tpiba*(0.d0,1.d0) * &
                ( xk(ipol,ik) + g(ipol,igk_k(ig,ik)) )
           enddo
        endif
     enddo
     if (noncolin) then
        call calbec(npw, dvkb(:,:,jpol), aux1, alphadk_nc(:,:,:,ipol))
     else
        call calbec(npw, dvkb(:,:,jpol), aux1, alphadk(:,:,ipol))
     endif
  enddo
#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds3')
  call start_clock('add_dkmds4')
#endif

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           mu = 3 * (na - 1)
           if ( abs (uact (mu + 1) ) + &
                abs (uact (mu + 2) ) + &
                abs (uact (mu + 3) ) > eps) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ipol = 1, 3
                       do ibnd=1, nbnd_occ(ik)
                          !
                          ! first we calculate the part coming from the
                          ! overlapp matrix S
                          !
                          if (noncolin) then
                             if (lspinorb) then
                                ps1_nc (ikb,1,ibnd)=ps1_nc(ikb,1,ibnd) +      &
                                     (qq_so(ih,jh,1,nt)*                      &
                                     alphadk_nc(jkb, 1, ibnd, ipol) +         &
                                      qq_so(ih,jh,2,nt)*                      &
                                     alphadk_nc(jkb, 2, ibnd, ipol) )*        &
                                      (0.d0,1.d0)*uact (mu + ipol)
                                ps1_nc (ikb,2,ibnd)=ps1_nc(ikb,2,ibnd) +      &
                                     (qq_so(ih,jh,3,nt)*                      &
                                     alphadk_nc(jkb, 1, ibnd, ipol) +         &
                                      qq_so(ih,jh,4,nt)*                      &
                                     alphadk_nc(jkb, 2, ibnd, ipol) )*        &
                                      (0.d0,1.d0)*uact (mu + ipol)
                                ps2_nc(ikb,1,ipol,ibnd)= &
                                            ps2_nc(ikb,1,ipol,ibnd)+ &
                                   (qq_so(ih,jh,1,nt)*becp2_nc(jkb,1,ibnd)+   &
                                    qq_so(ih,jh,2,nt)*becp2_nc(jkb,2,ibnd))* &
                                     uact (mu + ipol) * tpiba
                                ps2_nc(ikb,2,ipol,ibnd)= &
                                            ps2_nc(ikb,2,ipol,ibnd)+ &
                                   (qq_so(ih,jh,3,nt)*becp2_nc(jkb,1,ibnd)+   &
                                    qq_so(ih,jh,4,nt)*becp2_nc(jkb,2,ibnd))* &
                                     uact (mu + ipol) * tpiba
                                !
                                ! second part
                                !
                                ps1_nc(ikb,1,ibnd)=ps1_nc(ikb,1,ibnd) +     &
                                         (dpqq_so(ih,jh,1,jpol,nt)*      &
                                   alphap(ipol, ik)%nc(jkb,1,ibnd)+  &
                                          dpqq_so(ih,jh,2,jpol,nt)*          &
                                   alphap(ipol, ik)%nc(jkb,2,ibnd) )*&
                                      uact (mu + ipol)
                                ps1_nc(ikb,2,ibnd)=ps1_nc(ikb,2,ibnd) +     &
                                         (dpqq_so(ih,jh,3,jpol,nt)*   &
                                   alphap(ipol, ik)%nc(jkb,1,ibnd)+  &
                                          dpqq_so(ih,jh,4,jpol,nt)*        &
                                   alphap(ipol, ik)%nc(jkb,2,ibnd) )*&
                                      uact (mu + ipol)
                                ps2_nc(ikb,1,ipol,ibnd)= &
                                       ps2_nc(ikb,1,ipol,ibnd) +      &
                                      (dpqq_so(ih,jh,1,jpol,nt)*         &
                                       becp1(ik)%nc(jkb,1,ibnd)+   &
                                       dpqq_so(ih,jh,2,jpol,nt)*            &
                                       becp1(ik)%nc(jkb,2,ibnd))*  &
                                      (0.d0,-1.d0)*uact(mu+ipol)*tpiba
                                ps2_nc(ikb,2,ipol,ibnd)= &
                                       ps2_nc(ikb,2,ipol,ibnd) +      &
                                      (dpqq_so(ih,jh,3,jpol,nt)*          &
                                       becp1(ik)%nc(jkb,1,ibnd)+   &
                                       dpqq_so(ih,jh,4,jpol,nt)*       &
                                       becp1(ik)%nc(jkb,2,ibnd))*  &
                                      (0.d0,-1.d0)*uact(mu+ipol)*tpiba
                             else
                                do is=1,npol
                                   ps1_nc (ikb,is,ibnd)=ps1_nc(ikb,is,ibnd)+  &
                                     (0.d0,1.d0) * qq (ih, jh, nt) *          &
                                     alphadk_nc(jkb, is, ibnd, ipol) *        &
                                      uact (mu + ipol)
                                   ps2_nc(ikb,is,ipol,ibnd)= &
                                          ps2_nc(ikb,is,ipol,ibnd)+ &
                                     qq(ih,jh,nt)*becp2_nc(jkb, is, ibnd)*   &
                                     uact (mu + ipol) * tpiba

                                   ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd) + &
                                            dpqq(ih,jh,jpol,nt) *            &
                                      alphap(ipol, ik)%nc(jkb, is, ibnd)* &
                                      uact (mu + ipol)
                                   ps2_nc(ikb,is,ipol,ibnd)= &
                                          ps2_nc(ikb,is,ipol,ibnd) + &
                                          dpqq(ih,jh,jpol,nt)*(0.d0,-1.d0)* &
                                          becp1(ik)%nc(jkb, is, ibnd)*  &
                                          uact (mu + ipol) * tpiba
                                enddo
                             endif
                          else
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +           &
                                  (0.d0,1.d0) * qq (ih, jh, nt) *          &
                                  alphadk(jkb, ibnd, ipol) *               &
                                  uact (mu + ipol)
                             ps2 (ikb, ipol, ibnd) = ps2 (ikb, ipol, ibnd) +  &
                                  qq (ih, jh, nt) *                           &
                                  becp2(jkb, ibnd) *                          &
                                  uact (mu + ipol) * tpiba
                          !
                          !  and here the part of the matrix K(r)
                          !
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +      &
                                  dpqq(ih,jh,jpol,nt) *               &
                                 alphap(ipol, ik)%k(jkb, ibnd) *  &
                                  uact (mu + ipol)
                             ps2 (ikb, ipol, ibnd) = ps2 (ikb, ipol, ibnd) + &
                                  dpqq(ih,jh,jpol,nt)*(0.d0,-1.d0)*           &
                                  becp1(ik)%k(jkb, ibnd) *                 &
                                  uact (mu + ipol) * tpiba
                          endif
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds4')
  call start_clock('add_dkmds5')
#endif
  !
  !      This term is proportional to beta(k+q+G)
  !
  if (nkb.gt.0) then
     if (noncolin) then
        call zgemm ('N', 'N', npwq, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     else
        call zgemm ('N', 'N', npwq, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
!        dvpsi = matmul(vkb, ps1) + dvpsi
     endif
  endif
#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds5')
  call start_clock('add_dkmds6')
#endif
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        do ibnd = 1, nbnd
           if (noncolin) then
              ok = ok .or. (abs(ps2_nc(ikb,1,ipol,ibnd)).gt.eps )  &
                      .or. (abs(ps2_nc(ikb,2,ipol,ibnd)).gt.eps )
           else
              ok = ok.or. (abs (ps2 (ikb, ipol, ibnd)).gt.eps )
           endif
        enddo
        if (ok) then
           do ig = 1, npw
              igg = igk_k (ig,ikq)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, ik) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              if (noncolin) then
                 dvpsi(1:npw,ibnd) = ps2_nc(ikb,1,ipol,ibnd) * aux(1:npw) +    &
                                   dvpsi(1:npw,ibnd)
                 dvpsi(npwx+1:npwx+npw,ibnd)=ps2_nc(ikb,2,ipol,ibnd)  &
                                * aux(1:npw)+dvpsi(npwx+1:npwx+npw,ibnd)
              else
                 dvpsi(1:npw,ibnd) = ps2(ikb,ipol,ibnd) * aux(1:npw) +    &
                   dvpsi(1:npw,ibnd)
              endif
           enddo
        endif
     enddo
  enddo

  deallocate (aux)
  deallocate(aux1)
  if (noncolin) then
     if (allocated(ps1_nc)) deallocate(ps1_nc)
     if (allocated(ps2_nc)) deallocate(ps2_nc)
     if (allocated(alphadk_nc)) deallocate (alphadk_nc)
     if (allocated(becp2_nc)) deallocate (becp2_nc)
  else
     if (allocated(ps1))     deallocate(ps1)
     if (allocated(ps2))     deallocate(ps2)
     if (allocated(alphadk)) deallocate (alphadk)
     if (allocated(becp2)) deallocate (becp2)
  end if

#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds6')
  call stop_clock('add_dkmds')
#endif
  return

end subroutine add_dkmds
