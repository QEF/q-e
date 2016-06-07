!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

subroutine dft_exchange(nbnd_v,nbnd_s,n_set, e_x,ks_wfcs)
!this subroutine calculates the exchange
!energy term for every state and writes on disk


  
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : prefix, tmp_dir, iunwfc, nwordwfc
  USE mp_global,            ONLY : nproc_pool, me_pool
  USE kinds,    ONLY : DP
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi, RYTOEV
  USE wvfct,     ONLY : npwx, npw, nbnd, wg
  USE gvecw,     ONLY : gcutw
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2,bg
  USE wannier_gw
  USE gvect
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE uspp
  USE uspp_param,           ONLY : lmaxq,upf,nh, nhm
  USE wavefunctions_module, ONLY : psic
 ! USE realus,  ONLY : adduspos_gamma_r
  USE cell_base,            ONLY : at, bg, omega
  USE mp, ONLY : mp_sum, mp_bcast
  USE mp_world, ONLY : world_comm
  USE control_flags,        ONLY : gamma_only
  !USE exx, ONLY : exx_divergence_new, exx_grid_init, yukawa,exx_divergence_old
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE fft_base,             ONLY : dfftp
  USE io_global, ONLY : ionode
  USE lsda_mod,  ONLY : nspin

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

   INTEGER, INTENT(in)  :: nbnd_v(nspin) !number of valence states for both spin channels
   INTEGER, INTENT(in)  :: nbnd_s !number of states considered
   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
   REAL(kind=DP), INTENT(out) :: e_x(nbnd,nspin)!in output exchange energies
   COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npwx,nbnd,nspin)!all kohn sham wavefunctions

   REAL(kind=DP), ALLOCATABLE :: fac(:) 
   REAL(kind=DP) :: qq_fact,exxdiv
   INTEGER :: ig,iiv,iv,jjs,js,hw,ks
   REAL(kind=DP), ALLOCATABLE :: becpr(:,:)
   REAL(kind=DP), ALLOCATABLE :: tmpreal1(:), tmpreal_v(:,:),tmpreal_s(:,:)
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   INTEGER :: jmin,jmax
   COMPLEX(kind=DP), ALLOCATABLE :: prod_g(:),prod_c(:),prod_g2(:,:)
   REAL(kind=DP), ALLOCATABLE :: prod_r(:)
   REAL(kind=DP) :: exc
   INTEGER :: iun
   INTEGER, PARAMETER :: n_int=20
   REAL(kind=DP) :: qx,qy,qz
   INTEGER :: ix,iy,iz,n_int_loc,iunu
   REAL(kind=DP), ALLOCATABLE :: e_x_off(:,:,:)
   COMPLEX(kind=DP) :: c_exc
   INTEGER :: isv

  
   allocate(fac(ngm))
   if(l_whole_s) then
      allocate(e_x_off(nbnd_s,nbnd_s,nspin))
      e_x_off(:,:,:)=0.d0
   endif
!sets factors terms
!sets factors terms
!this has already  been called   call exx_grid_init()
   if(l_truncated_coulomb) then

      do ig=1,ngm
     
         qq_fact = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
         
         if (qq_fact > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq_fact))*(1.d0-dcos(dsqrt(qq_fact)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
            
         endif
         
      end do
      fac(:)=fac(:)/omega
   else 
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   endif

   e_x(:,:)=0.d0
   CALL gk_sort(xk(1,1),ngm,g,gcutw,npw0,igk0,g2kin_bp)

 
  


   allocate(tmpreal1(dfftp%nnr))
   allocate(tmpreal_v(dfftp%nnr,n_set))
   allocate(tmpreal_s(dfftp%nnr,n_set))
   allocate(prod_g(ngm),prod_g2(ngm,nbnd_s))
   allocate(prod_c(dfftp%nnr))
   allocate(prod_r(dfftp%nnr))

!external loop on valence state
   do isv=1,nspin
      do iiv=1,ceiling(real(nbnd_v(isv))/real(n_set))
   !read states and do fourier transform
         do hw=(iiv-1)*n_set+1,min(iiv*n_set,nbnd_v(isv)),2
            psic(:)=(0.d0,0.d0)
            psic(:)=(0.d0,0.d0)
            IF ( hw < min(iiv*n_set,nbnd_v(isv))) then
               psic(nls(1:npw0))  = ks_wfcs(1:npw0,hw,isv) + &
                    ( 0.D0, 1.D0 ) * ks_wfcs(1:npw0,hw+1,isv)
               psic(nlsm(1:npw0)) = CONJG( ks_wfcs(1:npw,hw,isv) - &
                    ( 0.D0, 1.D0 ) * ks_wfcs(1:npw0,hw+1,isv) )
            ELSE
               psic(nls(1:npw0))  = ks_wfcs(1:npw0,hw,isv)
               psic(nlsm(1:npw0)) = CONJG( ks_wfcs(1:npw0,hw,isv) )
            END IF
         
            CALL invfft ('Wave', psic, dffts)
            tmpreal1(1:dfftp%nnr)=dble(psic(1:dfftp%nnr))
            if(doublegrid) then
               call interpolate(tmpreal_v(:,hw-(iiv-1)*n_set),tmpreal1,1)
            else
               tmpreal_v(:,hw-(iiv-1)*n_set)=tmpreal1(:)
            endif
            if ( hw < min(iiv*n_set,nbnd_v(isv))) then
               tmpreal1(1:dfftp%nnr)=aimag(psic(1:dfftp%nnr))
               if(doublegrid) then
                  call interpolate(tmpreal_v(:,hw-(iiv-1)*n_set+1),tmpreal1,1)
               else
                  tmpreal_v(:,hw-(iiv-1)*n_set+1)=tmpreal1(:)
               endif
            endif
            
         enddo

      
         do jjs=1,ceiling(real(nbnd_s)/real(n_set))
   !external loop on states 
      !read states and do fourier transform
            do hw=(jjs-1)*n_set+1,min(jjs*n_set,nbnd_s),2
               psic(:)=(0.d0,0.d0)
               psic(:)=(0.d0,0.d0)
               IF ( hw < min(jjs*n_set,nbnd_s)) then
                  psic(nls(1:npw0))  = ks_wfcs(1:npw0,hw,isv) + &
                       ( 0.D0, 1.D0 ) * ks_wfcs(1:npw0,hw+1,isv)
                  psic(nlsm(1:npw0)) = CONJG( ks_wfcs(1:npw,hw,isv) - &
                       ( 0.D0, 1.D0 ) * ks_wfcs(1:npw0,hw+1,isv) )
               ELSE
                  psic(nls(1:npw0))  = ks_wfcs(1:npw0,hw,isv)
                  psic(nlsm(1:npw0)) = CONJG( ks_wfcs(1:npw0,hw,isv) )
               END IF
                  
               CALL invfft ('Wave', psic, dffts)
               tmpreal1(1:dfftp%nnr)=dble(psic(1:dfftp%nnr))
               if(doublegrid) then
                  call interpolate(tmpreal_s(:,hw-(jjs-1)*n_set),tmpreal1,1)
               else
                  tmpreal_s(:,hw-(jjs-1)*n_set)=tmpreal1(:)
               endif
               if ( hw < min(jjs*n_set,nbnd_s)) then
                  tmpreal1(1:dfftp%nnr)=aimag(psic(1:dfftp%nnr))
                  if(doublegrid) then
                     call interpolate(tmpreal_s(:,hw-(jjs-1)*n_set+1),tmpreal1,1)
                  else
                     tmpreal_s(:,hw-(jjs-1)*n_set+1)=tmpreal1(:)
                  endif
               endif
               
            enddo

      !internal loop on valence states
            do iv=(iiv-1)*n_set+1,min(iiv*n_set,nbnd_v(isv))
               
               jmin=(jjs-1)*n_set+1
               jmax=min(jjs*n_set,nbnd_s)

!for whole X operator for given iv calculate products in real space with all the 
!KS states and store in G space 
               if(l_whole_s) then
!NOT_TO_BE_INCLUDED_START
                  do ks=1,nbnd_s,1
                     psic(:)=(0.d0,0.d0)
                     psic(nls(1:npw0))  = ks_wfcs(1:npw0,ks,isv)
                     psic(nlsm(1:npw0)) = CONJG( ks_wfcs(1:npw0,ks,isv) )
                     CALL invfft ('Wave', psic, dffts)
                     prod_c(1:dfftp%nnr)=dcmplx(dble(psic(1:dfftp%nnr))*tmpreal_v(1:dfftp%nnr,iv-(iiv-1)*n_set)&
                          &  ,0.d0)
                     CALL fwfft ('Dense', prod_c, dfftp)
                     prod_g2(1:ngm,ks)=prod_c(nl(1:ngm))
                  enddo
!NOT_TO_BE_INCLUDED_END
               endif
               
               do js=jmin,jmax
                  !do product in real speace
                  prod_r(:)=tmpreal_v(:,iv-(iiv-1)*n_set)*tmpreal_s(:,js-(jjs-1)*n_set)
                     ! if(okvan) call adduspos_gamma_r & ATTENZIONE
                     ! (iv,js, prod_r(:),1,becpr(:,iv),becpr(:,js))

                  prod_c(:)=dcmplx(prod_r(:),0.d0)
                  CALL fwfft ('Dense', prod_c, dfftp)
                     !go to g_space
                  prod_g(1:ngm)=prod_c(nl(1:ngm))
                  !calculated exchange
                  exc=0.d0
                  do ig=1,ngm
                     exc=exc+2.d0*dble(conjg(prod_g(ig))*prod_g(ig))*fac(ig)*wg(iv,isv)*dble(nspin)/2.d0
                  enddo
                  if(gstart==2) exc=exc-dble(prod_g(1))*dble(prod_g(1))*fac(1)*wg(iv,isv)*dble(nspin)/2.d0
                  call mp_sum(exc,world_comm)
                  exc=-exc
                  e_x(js,isv)=e_x(js,isv)+exc
        
!poor programmer solution for off diagonal terms....
!ONLY FOR NORMCONSERVING PSEUDOS
                  if(l_whole_s) then
!NOT_TO_BE_INCLUDED_START
                     write(stdout,*) 'Call complete X operator part',iv
                     FLUSH(stdout)
                     do ks=1,nbnd_s,1
                        c_exc=(0.d0,0.d0)
                        do ig=1,ngm
                           c_exc=c_exc+conjg(prod_g2(ig,ks))*prod_g(ig)*fac(ig)+&
                                &prod_g2(ig,ks)*conjg(prod_g(ig))*fac(ig)
                        enddo
                        if(gstart==2) c_exc=c_exc-conjg(prod_g2(1,ks))*prod_g(1)*fac(1)
                        call mp_sum(c_exc,world_comm)
                        c_exc=-c_exc
                        e_x_off(ks,js,isv)=e_x_off(ks,js,isv)+dble(c_exc)
                     enddo
!NOT_TO_BE_INCLUDED_END
                  endif
               enddo
            enddo
         enddo
         
      enddo
   enddo!ivv
   do isv=1,nspin
      do iv=1,nbnd_s
         write(stdout,*) 'Exchange energy', iv,isv, e_x(iv,isv)
      enddo
   enddo
!write on file

   if(ionode) then
      iun = find_free_unit()
      open(unit=iun,file=trim(tmp_dir)//trim(prefix)//'.exchange',status='unknown',form='unformatted')
      write(iun) nbnd_s
      do isv=1,nspin
!NOT_TO_BE_INCLUDED_START
         if(l_selfconsistent)  e_x(1:nbnd_s,isv)=0.d0
!NOT_TO_BE_INCLUDED_END
         write(iun) e_x(1:nbnd_s,isv)
      enddo
      close(iun)
   endif

!if required write on disk off-diagonal terms


   if(l_whole_s) then
!NOT_TO_BE_INCLUDED_START
      if(ionode) then
         do iv=1,nbnd_s
            write(stdout,*) 'Exchange energy off', iv, e_x_off(iv,iv,1)
         enddo
      !write on file

         iun = find_free_unit()
         open(unit=iun,file=trim(tmp_dir)//trim(prefix)//'.exchange_off',status='unknown',form='unformatted')
         write(iun) nbnd_s
         do isv=1,nspin
            do js=1,nbnd_s
               write(iun) e_x_off(1:nbnd_s,js,isv)
            enddo
         enddo
         close(iun)
      endif
!NOT_TO_BE_INCLUDED_END
   endif


   deallocate(tmpreal1,tmpreal_s,tmpreal_v)

   deallocate(fac)
   deallocate(prod_c,prod_g,prod_g2)
   deallocate(prod_r)
   if(okvan) deallocate(becpr)
   if(l_whole_s) then
!NOT_TO_BE_INCLUDED_START
      deallocate(e_x_off)
!NOT_TO_BE_INCLUDED_END
   endif

 end subroutine dft_exchange



!----------------------------------------------------------------------
subroutine addus_charge(r_ij,becp_iw,becp_jw)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE gvect,                ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, &
                                   eigts3, mill
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, nkb
  USE uspp_param,           ONLY : lmaxq, upf, nh
  USE wavefunctions_module, ONLY : psic
  USE control_flags ,       ONLY : gamma_only
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft

  !
  implicit none

  COMPLEX(kind=DP), INTENT(inout) :: r_ij(dfftp%nnr)!where to add the us term
  COMPLEX(kind=DP), INTENT(in) ::  becp_iw( nkb)!overlap of wfcs with us  projectors
  COMPLEX(kind=DP), INTENT(in) ::  becp_jw( nkb)!overlap of wfcs with us  projectors



  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, is
  ! counters

  real(DP), allocatable :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics

  complex(DP) :: skk
  complex(DP), allocatable ::  aux (:,:), qgm(:)
  ! work space for rho(G,nspin)
  ! Fourier transform of q
  INTEGER, ALLOCATABLE :: ind_cor(:,:,:)
  INTEGER :: ijkb0, ikb,np

  if (.not.okvan) return


  allocate (aux ( ngm, nspin))
  allocate (qmod( ngm))
  allocate (qgm( ngm))
  allocate (ylmk0( ngm, lmaxq * lmaxq))

  aux (:,:) = (0.d0, 0.d0)
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
  
!found index correspondence


  allocate(ind_cor(ntyp,nat,maxval(nh(1:ntyp))))
  
  ijkb0 = 0
  do  np = 1, ntyp
     if ( upf(np)%tvanp ) then
        do  na = 1, nat
           if ( ityp(na) == np ) then
              do ih = 1, nh(np)
                ikb = ijkb0 + ih
                ind_cor(np,na,ih)=ikb
             enddo
             ijkb0=ijkb0+nh(np)
          endif
       enddo
    else
       do na=1,nat
          if(ityp(na) == np) ijkb0=ijkb0+nh(np)
       enddo
    endif
 enddo
 
 
  





  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        do ih = 1, nh (nt)
           do jh = 1, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    !
                    !  Multiply becsum and qg with the correct structure factor
                    !
                    do is = 1, nspin
                       do ig = 1, ngm
                          skk = eigts1 (mill(1,ig), na) * &
                                eigts2 (mill(2,ig), na) * &
                                eigts3 (mill(3,ig), na)
                          aux(ig,is)=aux(ig,is) + qgm(ig)*skk*&
                                  &conjg(becp_iw(ind_cor(nt,na,ih)))*becp_jw(ind_cor(nt,na,jh))
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo

  deallocate(ind_cor)
  !
  deallocate (ylmk0)
  deallocate (qgm)
  deallocate (qmod)
  !
  !     convert aux to real space and add to the charge density
  !
  do is = 1, nspin!SPIN TO BE IMPLEMENTED YET
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     if (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is))
     CALL invfft ('Dense', psic, dfftp)
     r_ij(:)=r_ij(:)+psic(:)
  enddo
  deallocate (aux)

  return
end subroutine addus_charge

