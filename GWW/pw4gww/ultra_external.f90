! FOR GWW
!
! Author: P. Umari
!
subroutine ultra_external( nbnd_start, nbnd_end, radius, itask)
!this subroutine calculates ultralocalized wanniers
!projecting out the part outside a sphere of radius R
!centered ad the center of each maximally localized generalized wannier
!
!it updates the localziation transform matrix
!and writes wavefunctions on file

!#ifdef __GWW

  USE io_files,             ONLY : find_free_unit,nwordwfc, iunwfc, prefix, diropn
  USE io_global,            ONLY : stdout, ionode_id
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE gvect,                ONLY : gstart
  use mp_global,            ONLY : nproc_pool, me_pool
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx
  USE basis
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE cell_base,            ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE uspp,                 ONLY : okvan,nkb, vkb
  USE uspp_param,           ONLY : lmaxq, nh, nhm
  USE realus,               ONLY : adduspos_real, augmentation_qq
  USE mp_global,            ONLY : intra_image_comm, me_pool
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE mp,                   ONLY : mp_sum
  USE wavefunctions_module, ONLY : psic, evc
  USE ions_base,            ONLY : nat
  USE klist,                ONLY : xk
  USE becmod,               ONLY : calbec

  implicit none

  INTEGER, INTENT(in) :: nbnd_start!first band defining the manifold to be localized
  INTEGER, INTENT(in) :: nbnd_end!last band defining the manifold to be localized
  REAL(kind=DP), INTENT(in) :: radius!localization radius
  INTEGER, INTENT(in) :: itask !if == 1 ultralocalizes {C'}
!--------------------------------
! definition of lecture unit for wannierized evc --> vc1
  integer :: iun_wannier
!---------------------------------
  COMPLEX(kind=DP), ALLOCATABLE :: wfc_mlwf(:,:)
  REAL(kind=DP), ALLOCATABLE :: ultra_trans(:,:)
  INTEGER :: n_man
  INTEGER :: iun_mlwf
  LOGICAL :: exst
  INTEGER :: ii,i,ix,iy,iz, nn, jj, ig,j,k
  COMPLEX(kind=DP), ALLOCATABLE :: op(:)
  INTEGER :: nr3s_start,nr3s_end,nr3_start,nr3_end
  REAL(kind=DP) :: r(3), rd(3), dist, sca
  COMPLEX(kind=DP), ALLOCATABLE :: op_wfc(:)
  REAL(kind=DP), ALLOCATABLE :: u_trans_old(:,:)
  INTEGER :: iunrealwan2
  REAL(kind=DP), ALLOCATABLE :: tmpreal(:),tmpreals(:)
  REAL(kind=DP), ALLOCATABLE :: becp_ext(:,:), qq_op(:,:,:), op_real(:)
  REAL(kind=DP), ALLOCATABLE :: becp_gw2(:,:)

!-------------------------------------------
! delcarion of temporary wannierized evc --> evc1
  COMPLEX(kind=DP), ALLOCATABLE :: evc1(:,:)
!-------------------------------------------

! open


  n_man=nbnd_end-nbnd_start+1
  allocate( wfc_mlwf( npwx, nbnd ) )
  allocate(ultra_trans(n_man,n_man))

!---------------------------------------
! allocation of temporary wannierized evc --> evc1
!!   allocate(evc1(npwx, nbnd))
  allocate(evc1(npw,nbnd_normal))
!----------------------------------------

  if(okvan) then
     allocate(becp_ext(nkb,nbnd))
     allocate(qq_op(nhm,nhm,nat))
  endif

#ifndef __MPI
  nr3s_start=1
  nr3s_end=dffts%nr3
  nr3_start=1
  nr3_end=dfftp%nr3
#else
  nr3s_start=0
  nr3s_end =0
  nr3_start=0
  nr3_end =0
  do i=1,me_pool + 1
     nr3s_start=nr3s_end+1
     nr3s_end=nr3s_end+dffts%npp(i)
     nr3_start=nr3_end+1
     nr3_end=nr3_end+dfftp%npp(i)
  end do
#endif



!read orthonormal wannier
  iun_mlwf=find_free_unit()
  CALL diropn( iun_mlwf, 'wfc_mlwf', nwordwfc, exst )
  CALL davcio(wfc_mlwf,nwordwfc,iun_mlwf,1,-1)
  close(iun_mlwf)

!calculate becs
  if  ( nkb > 0 .and. okvan) then
     CALL init_us_2( npw, igk, xk(1,1), vkb )
     call calbec(npw, vkb, wfc_mlwf, becp_ext, nbnd)
  endif


  allocate(op(dffts%nnr), op_wfc(npwx))
!loop on wavefunctions
  do ii=nbnd_start,nbnd_end

     write(stdout,*) 'II', ii!ATTENZIONE
     call flush_unit(stdout)
 !calculate operator O
     do ix=1,dffts%nr1
        do iy=1,dffts%nr2
           do iz=1,dffts%npp(me_pool+1)
              nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
              r(1)=(dble(ix-1)/dble(dffts%nr1))*at(1,1)*alat
              r(2)=(dble(iy-1)/dble(dffts%nr2))*at(2,2)*alat
              r(3)=(dble(iz+nr3s_start-1-1)/dble(dffts%nr3))*at(3,3)*alat
              rd(:)=r(:)-wannier_centers(:,ii)*alat
              do  i=1,3
                 if(rd(i) > at(i,i)*alat/2.d0) then
                    rd(i)=rd(i)-at(i,i)*alat
                 else if(rd(i) < -at(i,i)*alat/2.d0) then
                    rd(i)=rd(i)+at(i,i)*alat
                 endif
              enddo
              dist=dsqrt(rd(1)**2.d0 + rd(2)**2.d0 + rd(3)**2.d0)
              if(dist >= radius) then
                 op(nn)=(1.d0,0.d0)
              else
                 op(nn)=(0.d0, 0.d0)
              endif
           enddo
        enddo
     enddo
     if(okvan) then
        allocate(op_real(dfftp%nnr))
!calculate operator O
        do ix=1,dfftp%nr1
           do iy=1,dfftp%nr2
              do iz=1,dfftp%npp(me_pool+1)
                 nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                 r(1)=(dble(ix-1)/dble(dfftp%nr1))*at(1,1)*alat
                 r(2)=(dble(iy-1)/dble(dfftp%nr2))*at(2,2)*alat
                 r(3)=(dble(iz+nr3_start-1-1)/dble(dfftp%nr3))*at(3,3)*alat
                 rd(:)=r(:)-wannier_centers(:,ii)*alat
                 do  i=1,3
                    if(rd(i) > at(i,i)*alat/2.d0) then
                       rd(i)=rd(i)-at(i,i)*alat
                    else if(rd(i) < -at(i,i)*alat/2.d0) then
                       rd(i)=rd(i)+at(i,i)*alat
                    endif
                 enddo
                 dist=dsqrt(rd(1)**2.d0 + rd(2)**2.d0 + rd(3)**2.d0)
                 if(dist >= radius) then
                    op_real(nn)=1.d0
                 else
                    op_real(nn)=0.d0
                 endif
              enddo
           enddo
        enddo

        call augmentation_qq(op_real,qq_op)
        deallocate(op_real)
     endif
   !  write(*,*) 'ATTENZIONE1'
     call flush_unit(stdout)
 !calculate product wfc*O in real space and back on G space
        psic(:) = ( 0.D0, 0.D0 )
        psic(nls(igk(1:npw)))  = wfc_mlwf(1:npw,ii)
        psic(nlsm(igk(1:npw))) = CONJG( wfc_mlwf(1:npw,ii))
        CALL invfft ('Wave', psic, dffts)
        op(:)=op(:)*dble(psic(:))
        CALL fwfft ('Wave', op, dffts)
        op_wfc(1:npw)=op(nls(igk(1:npw)))
   ! calculates factors and normalizatio factor
        do jj=nbnd_start,nbnd_end
           sca=0.d0
           do ig=1,npw
              sca=sca+2.d0*dble(conjg(op_wfc(ig))*wfc_mlwf(ig,jj))
           enddo
           if(gstart==2) sca=sca-dble(conjg(op_wfc(1))*wfc_mlwf(1,jj))
           call mp_sum(sca)
!add US part
           if(okvan) call adduspos_real(sca,qq_op,becp_ext(:,ii),becp_ext(:,jj))
         !  write(*,*) 'SCA',ii,jj,sca
           call flush_unit(stdout)
           ultra_trans(ii-nbnd_start+1,jj-nbnd_start+1)=-sca
        enddo
       ! ultra_trans(ii-nbnd_start+1,ii-nbnd_start+1)=ultra_trans(ii-nbnd_start+1,ii-nbnd_start+1)+1.d0
        ultra_trans(ii-nbnd_start+1,ii-nbnd_start+1)=1.d0!ATTENZIONE

!now normalize
        sca=0.d0
        do jj=nbnd_start,nbnd_end
           sca=sca+ultra_trans(ii-nbnd_start+1,jj-nbnd_start+1)**2.d0
        enddo
        sca=sqrt(sca)
        do jj=nbnd_start,nbnd_end
           ultra_trans(ii-nbnd_start+1,jj-nbnd_start+1)=ultra_trans(ii-nbnd_start+1,jj-nbnd_start+1)/sca
        enddo
     enddo
  deallocate(op_wfc)
!update u_trans
 !I have : \Psi_j U_{j,i} = \tilde{w_i}
  allocate(u_trans_old(nbnd_normal,nbnd_normal))
  u_trans_old(:,:)=dble(u_trans(:,:))
  u_trans(nbnd_start:nbnd_end,nbnd_start:nbnd_end)=(0.d0,0.d0)
  do i=nbnd_start,nbnd_end
     do j=nbnd_start,nbnd_end
        do k=nbnd_start,nbnd_end
           u_trans(i,k)=u_trans(i,k)+u_trans_old(j,k)*ultra_trans(i-nbnd_start+1,j-nbnd_start+1)
        enddo
     enddo
  enddo
  deallocate(u_trans_old)
  !write wavefunctions on opportune file on R grid

  iunrealwan2 =  find_free_unit()
  if(itask/=1) then
     CALL diropn( iunrealwan2, 'realwan', dfftp%nnr, exst )
  else
     CALL diropn( iunrealwan2, 'realwan_prim', dfftp%nnr, exst )
  endif
  deallocate(wfc_mlwf,ultra_trans)
!rotate wavefunctions
!    ALLOCATE( evc( npwx, nbnd ) )
    allocate(tmpreal(dfftp%nnr),tmpreals(dffts%nnr))
!    CALL davcio(evc,nwordwfc,iunwfc,1,-1)

    allocate(u_trans_old(nbnd_normal,nbnd_normal))
    u_trans_old(:,:)=dble(u_trans(:,:))
    call rotate_wannier_gamma( u_trans_old,1,1)
    deallocate(u_trans_old)

!--------------------------------------------------
! as old evc has been transformed, now is evc1, by rotate_wannier_gamma
! davcio on evc1 iun_wannier
   iun_wannier = find_free_unit()
   call diropn(iun_wannier,"wfc_w",2*nwordwfc,exst)
   call davcio(evc1,nwordwfc,iun_wannier,1,-1)
   close(iun_wannier)
!---------------------------------------------------

!update becs
    if(okvan) then
       allocate(becp_gw2(nkb,nbnd))
!       call calbec(npw, vkb, evc, becp_gw2, nbnd)
! ----------------------------------------------------
! now is evc1 to be used
       call calbec(npw, vkb, evc1, becp_gw2, nbnd)
!-----------------------------------------------------
       if(itask/=1) then
          becp_gw(:,nbnd_start:nbnd_end)=becp_gw2(:,nbnd_start:nbnd_end)
       else
          becp_gw_c(:,nbnd_start:nbnd_end)=becp_gw2(:,nbnd_start:nbnd_end)
       endif
       deallocate(becp_gw2)
    endif


    do ii = nbnd_start, nbnd_end, 2
       psic(:) = ( 0.D0, 0.D0 )
       IF ( ii < nbnd_end ) THEN
          psic(nls(igk(1:npw)))  = evc1(1:npw,ii) + &
               ( 0.D0, 1.D0 ) * evc1(1:npw,ii+1)
          psic(nlsm(igk(1:npw))) = CONJG( evc1(1:npw,ii) - &
               ( 0.D0, 1.D0 ) * evc1(1:npw,ii+1) )
       ELSE
          psic(nls(igk(1:npw)))  = evc1(1:npw,ii)
          psic(nlsm(igk(1:npw))) = CONJG( evc1(1:npw,ii) )
       END IF
       CALL invfft ('Wave', psic, dffts)

!put on array
       tmpreals(:)=dble(psic(:))
       if(doublegrid) then
          call interpolate(tmpreal,tmpreals,1)
       else
          tmpreal(:)=tmpreals(:)
       endif
!writes on file

!if itask==1 writes {c'} from position 1

       if(itask /= 1) then
          call  davcio( tmpreal,dfftp%nnr,iunrealwan2,ii,1)
       else
          call  davcio( tmpreal,dfftp%nnr,iunrealwan2,ii-nbnd_start+1,1)
       endif
       if ( ii < nbnd_end ) then
!put on array

          tmpreals(:)=aimag(psic(:))
          if(doublegrid) then
             call interpolate(tmpreal,tmpreals,1)
          else
             tmpreal(:)=tmpreals(:)
          endif
!writes on file
!if itask==1 writes {c'} from position 1
          if(itask /= 1) then
             call  davcio( tmpreal,dfftp%nnr,iunrealwan2,ii+1,1)
          else
             call  davcio( tmpreal,dfftp%nnr,iunrealwan2,ii-nbnd_start+1+1,1)
          endif
       endif
    enddo

     deallocate(tmpreal,tmpreals)
!    deallocate(evc)
!---------------------------------------------
! as before evc1 is the one which presently in use
    deallocate(evc1)
!-------------------------------------------
    if(okvan) deallocate(becp_ext)
    if(okvan) deallocate(qq_op)
    close(iunrealwan2)


!#endif
  return

end subroutine ultra_external
