! FOR GWW
!
! Author: P. Umari
!
SUBROUTINE ultralocalization_para(nbndv,nbnd_max,ultra_thr,isubspace,max_array2, itask)

!this subroutine
!1-read R space wanniers
!2-ultralocalized
!3-determine cubic length
!4-calculates A matrix : \tilde{w}_i=A_{i,j}\Psi_j
!4-write on file

!this version writes the ultralocalized wannier functions
!on the big grid R, in order to use davcio for faster i/o
!on parallel architectures

!if alfa is equal to 1000 does not perform ultralocalization


!#ifdef __GWW

  USE io_files,             ONLY : find_free_unit, diropn
  USE io_global,            ONLY : stdout, ionode_id
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  use mp_global,            ONLY : nproc_pool, me_pool
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx
  USE basis
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE cell_base,            ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE uspp,                 ONLY : okvan,nkb
  USE realus,               ONLY : adduspos_gamma_r
  USE mp_global,            ONLY : intra_image_comm, me_pool
  USE fft_base,             ONLY : dfftp, dffts
  USE mp,                   ONLY : mp_bcast, mp_sum

  implicit none

  REAL(kind=DP), INTENT(in) :: ultra_thr!threshold for convergence
  INTEGER, INTENT(in) :: nbndv !number of first bands wich are localized separately
  INTEGER, INTENT(in) :: nbnd_max !number of bands, usually nbnd
  INTEGER, INTENT(in) :: isubspace! == 0 valence subspace; == 1 conduction subspace; == 2 second conduction subspace
  INTEGER, INTENT(in) :: max_array2! max number of states to be considered in the same array
  INTEGER, INTENT(in) :: itask !if == 1 ultralocalizes {C'}

  REAL(kind=DP), ALLOCATABLE ::  tmpreal(:),tmpreali(:,:),tmprealj(:)
  INTEGER :: iw,jw,kw,ir,it, ii,jj
  REAL(kind=DP) :: sca,sca1
  COMPLEX(kind=DP) :: scac1,scac2,scac3
  INTEGER :: iunrealwan
  LOGICAL :: exst
  INTEGER :: ix,iy,iz,nn
  REAL(kind=DP) :: center(3,max_array2),center_old(3,max_array2),rad, radmax
  REAL(kind=DP) :: rdistance
  REAL(kind=DP) :: rx,ry,rz
  REAL(kind=DP), ALLOCATABLE :: loc_mat(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: eigenvector(:,:),eigenvector_old(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: eigenvector2(:)
  INTEGER, PARAMETER :: max_loc_ite=1!5 ATTENZIONE
  INTEGER :: lwork, m,info
  REAL(kind=DP), ALLOCATABLE :: work(:)
  REAL(kind=DP) :: eigen(max_array2),eigenold(max_array2)
  INTEGER, ALLOCATABLE :: iwork(:),ifail(:)
  COMPLEX(kind=DP), ALLOCATABLE :: eigx(:,:),eigy(:,:),eigz(:,:)
  COMPLEX(kind=DP):: center_berry
  INTEGER :: nc1,nc2,nc3,n1,n2,n3
  INTEGER :: nll,no1,no2,no3,nrsmin
  INTEGER :: iunrealwan2
  CHARACTER(5) :: nfile
  REAL(kind=DP) :: norm
  COMPLEX(kind=DP),allocatable :: exp_x(:),exp_y(:),exp_z(:)
  INTEGER :: n_first,n_last, n_bands,jww,kww,iww
  INTEGER, ALLOCATABLE :: min_1(:,:,:),max_1(:,:,:),min_2(:,:,:),max_2(:,:,:)
  LOGICAL :: converged(nbnd_normal)
  REAL(kind=DP), ALLOCATABLE :: sums(:,:,:)
  INTEGER :: iq,ifirst,ilast,iqq
  COMPLEX(kind=DP), ALLOCATABLE :: c_mat(:,:)
  REAL(kind=DP) :: cutoff
  REAL(kind=DP) :: alfa
  REAL(kind=DP), ALLOCATABLE :: becp_gw2(:,:)
  REAL(kind=DP), ALLOCATABLE :: tmp_s(:),tmp_r(:)
  LOGICAL :: is_even
  INTEGER, PARAMETER :: nmaxeig = 4!4 ATTENZIONE
  REAL(kind=DP), PARAMETER :: minover = 0.5d0
  REAL(kind=DP) :: loc_tmp

  REAL(kind=DP), ALLOCATABLE :: eig_set(:), eigvector_tmp(:)
  REAL(kind=DP), ALLOCATABLE :: eigvector_set(:,:)
  INTEGER :: itt, ipp,ipp_exit
  LOGICAL :: good

  INTEGER :: nr3s_start, nr3s_end
  INTEGER :: nr3_start, nr3_end

  INTEGER :: ivv
  INTEGER :: numnbndset, nbnd_start, nbnd_end

  REAL(kind=DP) :: loc_mat_tmp
  REAL(kind=DP) :: scamax
  INTEGER :: i_max, nmaxeig_tmp
  INTEGER :: max_array
  INTEGER :: istop

  max_array=max_array2
#ifndef __PARA
  !dfftp%npp(1)= nr3 ! not needed - PG
  !dffts%npp(1)= dffts%nr3 ! not needed - PG
  nr3s_start=1
  nr3s_end=dffts%nr3
  nr3_start=1
  nr3_end=dfftp%nr3
#else
  nr3s_start=0
  nr3s_end =0
  nr3_start=0
  nr3_end =0
  do ii=1,me_pool + 1
     nr3s_start=nr3s_end+1
     nr3s_end=nr3s_end+dffts%npp(ii)
     nr3_start=nr3_end+1
     nr3_end=nr3_end+dfftp%npp(ii)
  end do
#endif



  allocate(eig_set(nmaxeig))
  allocate(eigvector_set(nbnd_normal,nmaxeig),eigvector_tmp(nbnd_normal))


  allocate(tmpreal(dfftp%nnr))
  allocate(eigenvector2(nbnd_normal),eigenvector_old(nbnd_normal,max_array))
  allocate(iwork(5*nbnd_normal),ifail(nbnd_normal))
  allocate(eigx(nbnd_normal,nbnd_normal),eigy(nbnd_normal,nbnd_normal),eigz(nbnd_normal,nbnd_normal))
  allocate(exp_x(dfftp%nnr),exp_y(dfftp%nnr),exp_z(dfftp%nnr))
  allocate(sums(dfftp%nr1,dfftp%nr2,dfftp%npp(me_pool+1)))
  allocate(tmp_s(dffts%nnr),tmp_r(dfftp%nnr))
  if(okvan) allocate(becp_gw2(nkb,nbnd))

  if(isubspace==0) then
     n_first=1
     n_last=nbndv
     cutoff=cutoff_wsq
     numnbndset=nbndv
  else if(isubspace==1.or.isubspace==2) then
     if(itask/=1) then
        n_first=nbndv+1
        n_last=nbnd_max
        cutoff=cutoff_wsq_c
        numnbndset=n_last-nbndv
     else
        n_first=nbndv+1
        n_last=nbndv+num_nbndc_set
        cutoff=cutoff_wsq_c
        numnbndset=n_last-nbndv
     endif
  else
     write(stdout,*) 'ultralocalization isubspace ILLEGAL'
     stop
  endif
  if(isubspace==1 .and. ultra_alpha_c==1000.d0) numnbndset=nset
  if(isubspace==2) numnbndset=nset
  if((n_last-n_first+1) < nmaxeig) then
     nmaxeig_tmp=n_last-n_first+1
  else
     nmaxeig_tmp=nmaxeig
  endif

  if(isubspace==0) then
     alfa=ultra_alpha_v
  else if(isubspace==1) then
     alfa=ultra_alpha_c
     if(itask==1) alfa=ultra_alpha_c_prim
  else if(isubspace==2) then
     alfa=ultra_alpha_c2
  endif



  allocate(tmpreali(dffts%nnr,numnbndset))
  allocate(tmprealj(dffts%nnr))
  allocate(min_1(dfftp%nr2,dfftp%npp(me_pool+1),numnbndset),min_2(dfftp%nr2,dfftp%npp(me_pool+1),numnbndset))
  allocate(max_1(dfftp%nr2,dfftp%npp(me_pool+1),numnbndset),max_2(dfftp%nr2,dfftp%npp(me_pool+1),numnbndset))
  allocate(eigenvector(nbnd_normal,numnbndset))

  exp_x(:)=(0.d0,0.d0)
  exp_y(:)=(0.d0,0.d0)
  exp_z(:)=(0.d0,0.d0)
  do ix=1,dffts%nr1
     do iy=1,dffts%nr2
        do iz=1,dffts%npp(me_pool+1)
           nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
           exp_x(nn)=exp((0.d0,1.d0)*tpi*dble(ix-1)/dble(dffts%nr1))
           exp_y(nn)=exp((0.d0,1.d0)*tpi*dble(iy-1)/dble(dffts%nr2))
           exp_z(nn)=exp((0.d0,1.d0)*tpi*dble(iz+nr3s_start-1-1)/dble(dffts%nr3))
        enddo
     enddo
  enddo




  do ivv=1,ceiling(real(n_last-n_first+1)/real(numnbndset))
     nbnd_start=(ivv-1)*numnbndset+n_first
     nbnd_end=nbnd_start+numnbndset-1
     if(nbnd_end > n_last) nbnd_end=n_last
     n_bands=nbnd_end-nbnd_start+1

     radmax=no_radius/alat

     nrsmin=min(dfftp%nr1,dfftp%nr2)
     nrsmin=min(dfftp%nr3,nrsmin)
 !if nrsmin is even set to nrsmin -1,
     if(is_even(nrsmin)) then
        nrsmin=nrsmin-1
        write(stdout,*) 'nrsmin set to:', nrsmin
     endif

     write(stdout,*)' ULTRALOCALIZATION PARA', nbnd_start, nbnd_end !ATTENZIONE
     CALL flush_unit( stdout )
!open outpu file

     iunrealwan = find_free_unit()
     CALL diropn( iunrealwan, 'real_whole', dffts%nnr, exst )
     iunrealwan2 =  find_free_unit()

     if(itask/=1) then
        CALL diropn( iunrealwan2, 'realwan', dfftp%nnr, exst )
     else
        CALL diropn( iunrealwan2, 'realwan_prim', dfftp%nnr, exst )
     endif

!read in wave-functions
     do iw=nbnd_start,nbnd_end
        CALL davcio( tmpreali(:,iw-nbnd_start+1),dffts%nnr,iunrealwan,iw,-1)
     enddo
     CLOSE(iunrealwan)
!first valence subspace

     if(alfa /= 1000d0) allocate(loc_mat(n_bands,n_bands,max_array))
     lwork=8*nbnd_normal
     allocate(work(lwork))

     exp_x(:)=(0.d0,0.d0)
     exp_y(:)=(0.d0,0.d0)
     exp_z(:)=(0.d0,0.d0)
     do ix=1,dffts%nr1
        do iy=1,dffts%nr2
           do iz=1,dffts%npp(me_pool+1)
              nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
              exp_x(nn)=exp((0.d0,1.d0)*tpi*dble(ix-1)/dble(dffts%nr1))
              exp_y(nn)=exp((0.d0,1.d0)*tpi*dble(iy-1)/dble(dffts%nr2))
              exp_z(nn)=exp((0.d0,1.d0)*tpi*dble(iz+nr3s_start-1-1)/dble(dffts%nr3))
           enddo
        enddo
     enddo

     do iw=nbnd_start,nbnd_end
        if(alfa/=1000d0) then
           istop=nbnd_end
        else
           istop=iw
        endif
        do jw=iw,istop
           iww=iw-nbnd_start+1
           jww=iw-nbnd_start+1
           eigx(iw,jw)=(0.d0,0.d0)
           eigy(iw,jw)=(0.d0,0.d0)
           eigz(iw,jw)=(0.d0,0.d0)
           tmp_s(:)=tmpreali(:,iww)*tmpreali(:,jww)

           do ir=1,dffts%nnr
              eigx(iw,jw)=eigx(iw,jw)+exp_x(ir)*tmp_s(ir)
              eigy(iw,jw)=eigy(iw,jw)+exp_y(ir)*tmp_s(ir)
              eigz(iw,jw)=eigz(iw,jw)+exp_z(ir)*tmp_s(ir)
           enddo
           eigx(iw,jw)=eigx(iw,jw)/real(dffts%nr1*dffts%nr2*dffts%nr3)
           eigy(iw,jw)=eigy(iw,jw)/real(dffts%nr1*dffts%nr2*dffts%nr3)
           eigz(iw,jw)=eigz(iw,jw)/real(dffts%nr1*dffts%nr2*dffts%nr3)

           call mp_sum(eigx(iw,jw))
           call mp_sum(eigy(iw,jw))
           call mp_sum(eigz(iw,jw))

           eigx(jw,iw)=eigx(iw,jw)
           eigy(jw,iw)=eigy(iw,jw)
           eigz(jw,iw)=eigz(iw,jw)
        enddo
     enddo



     if(okvan) then
        exp_x(:)=(0.d0,0.d0)
        exp_y(:)=(0.d0,0.d0)
        exp_z(:)=(0.d0,0.d0)
        do ix=1,dfftp%nr1
           do iy=1,dfftp%nr2
              do iz=1,dfftp%npp(me_pool+1)
                 nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                 exp_x(nn)=exp((0.d0,1.d0)*tpi*dble(ix-1)/dble(dfftp%nr1))
                 exp_y(nn)=exp((0.d0,1.d0)*tpi*dble(iy-1)/dble(dfftp%nr2))
                 exp_z(nn)=exp((0.d0,1.d0)*tpi*dble(iz+nr3_start-1-1)/dble(dfftp%nr3))
              enddo
           enddo
        enddo

        do iw=nbnd_start,nbnd_end
           if(alfa/=1000d0) then
              istop=nbnd_end
           else
              istop=iw
           endif
           do jw=iw,istop
              iww=iw-nbnd_start+1
              jww=iw-nbnd_start+1
              tmp_r(:)=0.d0

              if(itask == 0) then
                 call adduspos_gamma_r(iw,jw,tmp_r,1,becp_gw(:,iw),becp_gw(:,jw))
              else
                 call adduspos_gamma_r(iw,jw,tmp_r,1,becp_gw_c(:,iw),becp_gw_c(:,jw))
              endif
              scac1=(0.d0,0.d0)
              scac2=(0.d0,0.d0)
              scac3=(0.d0,0.d0)
              do ir=1,dfftp%nnr
                 scac1=scac1+exp_x(ir)*tmp_r(ir)
                 scac2=scac2+exp_y(ir)*tmp_r(ir)
                 scac3=scac3+exp_z(ir)*tmp_r(ir)
              enddo

              call mp_sum(scac1)
              call mp_sum(scac2)
              call mp_sum(scac3)

              eigx(iw,jw)=eigx(iw,jw)+scac1/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
              eigy(iw,jw)=eigy(iw,jw)+scac2/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
              eigz(iw,jw)=eigz(iw,jw)+scac3/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
              eigx(jw,iw)=eigx(iw,jw)
              eigy(jw,iw)=eigy(iw,jw)
              eigz(jw,iw)=eigz(iw,jw)
           enddo
        enddo
     endif


     max_array=min(n_bands,max_array2)


     do iq=1,n_bands/max_array
        ifirst=(iq-1)*max_array+nbnd_start
        ilast=min(iq*max_array+nbnd_start,nbnd_end)
        converged(:)=.false.
        iqq=0
        do iw=ifirst,ilast
           iqq=iqq+1
!read wavefunction

           center(1,iqq)=aimag(log(eigx(iw,iw)))*at(1,1)/tpi
           center(2,iqq)=aimag(log(eigy(iw,iw)))*at(2,2)/tpi
           center(3,iqq)=aimag(log(eigz(iw,iw)))*at(3,3)/tpi
!now loop till convergence
           center_old(:,iqq)=center(:,iqq)
        enddo
        if( alfa /= 1000d0) then
           do it=1,max_loc_ite

 !calculate screen function with center and radius
           iqq=0
           do iw=ifirst,ilast
              iqq=iqq+1
              if(.not.converged(iqq)) then


                 write(stdout,*) 'Ultralocalization preparing distance'
                 CALL flush_unit( stdout )
                 do iy=1,dffts%nr2
                    do iz=1,dffts%npp(me_pool+1)
                       min_1(iy,iz,iqq)=0
                       max_1(iy,iz,iqq)=0
                       min_2(iy,iz,iqq)=0
                       max_2(iy,iz,iqq)=0
                       do ix=1,dffts%nr1
                          rx=rdistance(real(ix-1)*at(1,1)/real(dffts%nr1),center(1,iqq),at(1,1))
                          ry=rdistance(real(iy-1)*at(2,2)/real(dffts%nr2),center(2,iqq),at(2,2))
                          rz=rdistance(real(iz+nr3s_start-1-1)*at(3,3)/real(dffts%nr3),center(3,iqq),at(3,3))
                          if(sqrt(rx**2.d0+ry**2.d0+rz**2.d0) <= radmax) then
                             if(min_1(iy,iz,iqq)==0) min_1(iy,iz,iqq)=ix
                          else
                             if(max_1(iy,iz,iqq)==0 .and. min_1(iy,iz,iqq)/=0)  then
                                max_1(iy,iz,iqq)=ix
                                exit
                             endif
                          endif
                       enddo
                       if(min_1(iy,iz,iqq)/=0 .and. max_1(iy,iz,iqq)==0)  max_1(iy,iz,iqq)=dffts%nr1+1
                       if(min_1(iy,iz,iqq)==1 .and. max_1(iy,iz,iqq)/=(dffts%nr1+1)) then
                          do ix=dffts%nr1,1,-1
                             rx=rdistance(real(ix-1)*at(1,1)/real(dffts%nr1),center(1,iqq),at(1,1))
                             ry=rdistance(real(iy-1)*at(2,2)/real(dffts%nr2),center(2,iqq),at(2,2))
                             rz=rdistance(real(iz+nr3s_start-1-1)*at(3,3)/real(dffts%nr3),center(3,iqq),at(3,3))
                             if(sqrt(rx**2.d0+ry**2.d0+rz**2.d0) <= radmax) then
                                if(max_2(iy,iz,iqq)==0) max_2(iy,iz,iqq)=ix
                             else
                                if(min_2(iy,iz,iqq)==0 .and. max_2(iy,iz,iqq)/=0) then
                                   min_2(iy,iz,iqq)=ix
                                   exit
                                endif
                             endif
                          enddo
                       else
                          max_2(iy,iz,iqq)=0
                          min_2(iy,iz,iqq)=0
                       endif
                       if(max_2(iy,iz,iqq)==(max_1(iy,iz,iqq)-1)) then
                          max_2(iy,iz,iqq)=0
                          min_2(iy,iz,iqq)=0
                       endif
                    enddo
                 enddo
              endif
           enddo
!now calculate localization matrix
           loc_mat(:,:,:)=0.d0
           write(stdout,*) 'Ultralocalization calculate matrix'
           CALL flush_unit( stdout )
           do jw=nbnd_start,nbnd_end
              jww=jw-nbnd_start+1
              do kw=jw,nbnd_end
                 kww=kw-nbnd_start+1
                 sums(:,:,:)=0.d0
                 tmp_s(:)=tmpreali(:,jww)*tmpreali(:,kww)
                 do iy=1,dffts%nr2
                    do iz=1,dffts%npp(me_pool+1)
                       sca=0.d0
                       nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x
                       do ix=1,dffts%nr1
                          !nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
                          nn=nn+1
                          sca=sca+tmp_s(nn)
                          sums(ix,iy,iz)=sca
                       enddo
                    enddo
                 enddo
                 iqq=0
                 do iw=ifirst,ilast
                    iqq=iqq+1
                    if(.not.converged(iqq))then
                       loc_mat_tmp=0.d0
                          do iz=1,dffts%npp(me_pool+1)
                             do iy=1,dffts%nr2
                             if(max_1(iy,iz,iqq)/=0) then
                                if(min_1(iy,iz,iqq)/=1) then
                                   loc_mat_tmp=loc_mat_tmp+&
                                        &sums(max_1(iy,iz,iqq)-1,iy,iz)-sums(min_1(iy,iz,iqq)-1,iy,iz)
                                else
                                   loc_mat_tmp=loc_mat_tmp + sums(max_1(iy,iz,iqq)-1,iy,iz)
                                endif
                             endif
                             if(max_2(iy,iz,iqq)/=0) then
                                loc_mat_tmp=loc_mat_tmp+sums(max_2(iy,iz,iqq),iy,iz)-sums(min_2(iy,iz,iqq),iy,iz)
                             endif
                          enddo
                       enddo
                       loc_mat(jww,kww,iqq)=loc_mat_tmp/real(dffts%nr1*dffts%nr2*dffts%nr3)
                       call mp_sum(loc_mat(jww,kww,iqq))
                       loc_mat(kww,jww,iqq)=loc_mat(jww,kww,iqq)
                    endif
                 enddo
              enddo
           enddo



           if(okvan) then
              iqq=0

              write(stdout,*) 'Ultralocalization preparing distance US'
              CALL flush_unit( stdout )
              do iw=ifirst,ilast
                 iqq=iqq+1
                 if(.not.converged(iqq)) then

                    do iy=1,dfftp%nr2
                       do iz=1,dfftp%npp(me_pool+1)
                          min_1(iy,iz,iqq)=0
                          max_1(iy,iz,iqq)=0
                          min_2(iy,iz,iqq)=0
                          max_2(iy,iz,iqq)=0
                          do ix=1,dfftp%nr1
                             rx=rdistance(real(ix-1)*at(1,1)/real(dfftp%nr1),center(1,iqq),at(1,1))
                             ry=rdistance(real(iy-1)*at(2,2)/real(dfftp%nr2),center(2,iqq),at(2,2))
                             rz=rdistance(real(iz+nr3_start-1-1)*at(3,3)/real(dfftp%nr3),center(3,iqq),at(3,3))
                             if(sqrt(rx**2.d0+ry**2.d0+rz**2.d0) <= radmax) then
                                if(min_1(iy,iz,iqq)==0) min_1(iy,iz,iqq)=ix
                             else
                                if(max_1(iy,iz,iqq)==0 .and. min_1(iy,iz,iqq)/=0)  then
                                   max_1(iy,iz,iqq)=ix
                                   exit
                                endif
                             endif
                          enddo
                          if(min_1(iy,iz,iqq)/=0 .and. max_1(iy,iz,iqq)==0)  max_1(iy,iz,iqq)=dfftp%nr1+1
                          if(min_1(iy,iz,iqq)==1 .and. max_1(iy,iz,iqq)/=(dfftp%nr1+1)) then
                             do ix=dfftp%nr1,1,-1
                                rx=rdistance(real(ix-1)*at(1,1)/real(dfftp%nr1),center(1,iqq),at(1,1))
                                ry=rdistance(real(iy-1)*at(2,2)/real(dfftp%nr2),center(2,iqq),at(2,2))
                                rz=rdistance(real(iz+nr3_start-1-1)*at(3,3)/real(dfftp%nr3),center(3,iqq),at(3,3))
                                if(sqrt(rx**2.d0+ry**2.d0+rz**2.d0) <= radmax) then
                                   if(max_2(iy,iz,iqq)==0) max_2(iy,iz,iqq)=ix
                                else
                                   if(min_2(iy,iz,iqq)==0 .and. max_2(iy,iz,iqq)/=0) then
                                      min_2(iy,iz,iqq)=ix
                                      exit
                                   endif
                                endif
                             enddo
                          else
                             max_2(iy,iz,iqq)=0
                             min_2(iy,iz,iqq)=0
                          endif
                          if(max_2(iy,iz,iqq)==(max_1(iy,iz,iqq)-1)) then
                             max_2(iy,iz,iqq)=0
                             min_2(iy,iz,iqq)=0
                          endif
                       enddo
                    enddo
                 endif
              enddo
!now calculate localization matrix
              write(stdout,*) 'Ultralocalization calculate matrix US'
              CALL flush_unit( stdout )
        do jw=nbnd_start,nbnd_end
           jww=jw-nbnd_start+1
           do kw=jw,nbnd_end
              kww=kw-nbnd_start+1
              sums(:,:,:)=0.d0
              tmp_r(:)=0.d0
              if(itask==0) then
                 call adduspos_gamma_r(jw,kw,tmp_r,1,becp_gw(:,jw),becp_gw(:,kw))
              else
                 call adduspos_gamma_r(jw,kw,tmp_r,1,becp_gw_c(:,jw),becp_gw_c(:,kw))
              endif
              do iy=1,dfftp%nr2
                 do iz=1,dfftp%npp(me_pool+1)
                    sca=0.d0
                    nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x
                    do ix=1,dfftp%nr1
                       nn=nn+1
                       !nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                       sca=sca+tmp_r(nn)
                       sums(ix,iy,iz)=sca
                    enddo
                 enddo
              enddo
              sums(:,:,:)=sums(:,:,:)/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
              iqq=0
              do iw=ifirst,ilast
                 iqq=iqq+1
                 if(.not.converged(iqq))then
                    loc_tmp=0.d0
                       do iz=1,dfftp%npp(me_pool+1)
                          do iy=1,dfftp%nr2
                          if(max_1(iy,iz,iqq)/=0) then
                             if(min_1(iy,iz,iqq)/=1) then
                                loc_tmp=loc_tmp+sums(max_1(iy,iz,iqq)-1,iy,iz)-sums(min_1(iy,iz,iqq)-1,iy,iz)
                             else
                                loc_tmp=loc_tmp+ sums(max_1(iy,iz,iqq)-1,iy,iz)
                             endif
                          endif
                          if(max_2(iy,iz,iqq)/=0) then
                             loc_tmp=loc_tmp+sums(max_2(iy,iz,iqq),iy,iz)-sums(min_2(iy,iz,iqq),iy,iz)
                          endif
                       enddo
                    enddo
                    call mp_sum(loc_tmp)
                    loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq)+loc_tmp
                    loc_mat(kww,jww,iqq)=loc_mat(jww,kww,iqq)
                 endif
              enddo
           enddo
        enddo
     endif








     write(stdout,*) 'Ultralocalization finding eigenvectors'
     CALL flush_unit( stdout )

!loop on ultra-localized wanniers
     iqq=0
     do iw=ifirst,ilast
        iqq=iqq+1
        if(.not.converged(iqq)) then
           loc_mat(iw-nbnd_start+1,iw-nbnd_start+1,iqq) = loc_mat(iw-nbnd_start+1,iw-nbnd_start+1,iqq) + alfa
           if(it==1)  eigenold(iqq)=loc_mat(iw-nbnd_start+1,iw-nbnd_start+1,iqq)

!find upperstate
!for avoiding numerical instabilities only the first processor perform diagonalization
           eig_set(:)=0.d0
           eigvector_set(:,:)=0.d0
           if(me_pool == 0) then
              call dsyevx('V','I','U',n_bands,loc_mat(:,:,iqq),n_bands,0.d0,0.d0,n_bands-nmaxeig_tmp+1,n_bands,0.d0,m,eig_set,&
                   & eigvector_set,nbnd_normal,work,lwork,iwork,ifail,info)
              if(info/=0) then
                 write(stdout,*) 'Error from dsyevx: ', info
                 stop
              endif

           endif

           CALL mp_bcast( eig_set(:), ionode_id )
           CALL mp_bcast( eigvector_set(:,:), ionode_id )
           scamax=0.d0
           i_max=0
           do ipp=nmaxeig_tmp,1,-1
              good=.true.
              eigvector_tmp(:)=eigvector_set(:,ipp)
              sca=0.d0
              do  ii=1,n_bands
                 sca=sca+eigvector_tmp(ii)**2.d0
              enddo
              write(stdout,*) 'Modulus vector', sca
              CALL flush_unit( stdout )
              do itt=1,iqq-1
                 sca=0.d0
                 do  ii=1,n_bands
                    sca=sca+eigenvector(ii,itt)*eigvector_tmp(ii)
                 enddo
                 eigvector_tmp(:)=eigvector_tmp(:)-sca*eigenvector(:,itt)
              enddo
              sca=0.d0
              do  ii=1,n_bands
                 sca=sca+eigvector_tmp(ii)**2.d0
              enddo
              if(sca>=scamax) then
                 scamax=sca
                 i_max=ipp
              endif
              if(sca >=   minover) then
                 ipp_exit=ipp
                 exit
              endif
           enddo
           if(sca < minover) ipp_exit=i_max

          !if(ipp<nmaxeig) write(stdout,*) 'Ultraloc too close', iqq,ipp
           write(stdout,*) 'Ultraloc', iqq,ipp_exit, sca
           if(ipp_exit==0) ipp_exit=1
           eigen(iqq)=eig_set(ipp_exit)
           eigenvector(:,iw-nbnd_start+1)=eigvector_set(:,ipp_exit)
           write(stdout,*) 'Ultra:',iw,it,eigen(iqq),eigenold(iqq)!ATTENZIONE
           if(eigen(iqq)<eigenold(iqq).and.it>1) then
              write(stdout, *) 'Increasing..stopping'
              eigenvector(:,iw-nbnd_start+1)=eigenvector_old(:,iqq)
              center(:,iqq)=center_old(:,iqq)
              converged(iqq)=.true.
           endif
           if(abs(eigen(iqq)-eigenold(iqq))<= ultra_thr) converged(iqq)=.true.
           if( it == max_loc_ite) converged(iqq)=.true.
           if(.not.converged(iqq)) then
              eigenold(iqq)=eigen(iqq)
              eigenvector_old(:,iqq)=eigenvector(:,iw-nbnd_start+1)
              center_old(:,iqq)=center(:,iqq)


!update center
              do ii=nbnd_start,nbnd_end
                 eigenvector2(ii)=(0.d0,0.d0)
                 do jj=nbnd_start,nbnd_end
                    eigenvector2(ii)=eigenvector2(ii)+eigx(ii,jj)*eigenvector(jj-nbnd_start+1,iw-nbnd_start+1)
                 enddo
              enddo
              center_berry=(0.d0,0.d0)
              do ii=nbnd_start,nbnd_end
                 center_berry=center_berry+eigenvector(ii-nbnd_start+1,iw-nbnd_start+1)*eigenvector2(ii)
              enddo
              center(1,iqq)=aimag(log(center_berry))*at(1,1)/tpi

              do ii=nbnd_start,nbnd_end
                 eigenvector2(ii)=(0.d0,0.d0)
                 do jj=nbnd_start,nbnd_end
                    eigenvector2(ii)=eigenvector2(ii)+eigy(ii,jj)*eigenvector(jj-nbnd_start+1,iw-nbnd_start+1)
                 enddo
              enddo
              center_berry=(0.d0,0.d0)
              do ii=nbnd_start,nbnd_end
                 center_berry=center_berry+eigenvector(ii-nbnd_start+1,iw-nbnd_start+1)*eigenvector2(ii)
              enddo
              center(2,iqq)=aimag(log(center_berry))*at(2,2)/tpi

              do ii=nbnd_start,nbnd_end
                 eigenvector2(ii)=(0.d0,0.d0)
                 do jj=nbnd_start,nbnd_end
                    eigenvector2(ii)=eigenvector2(ii)+eigz(ii,jj)*eigenvector(jj-nbnd_start+1,iw-nbnd_start+1)
                 enddo
              enddo
              center_berry=(0.d0,0.d0)
              do ii=nbnd_start,nbnd_end
                 center_berry=center_berry+eigenvector(ii-nbnd_start+1,iw-nbnd_start+1)*eigenvector2(ii)
              enddo
              center(3,iqq)=aimag(log(center_berry))*at(3,3)/tpi
           endif
        endif
     enddo!on iw
  enddo!on iterations it
  endif!alfa/=1000
!write on file
  iqq=0
  do iw=ifirst,ilast
     iqq=iqq+1
     write(stdout,*) 'Ultralocalization overlap with starting wannier:', iw
     CALL flush_unit( stdout )
!construct ultralocalized wannier in real space
     if( alfa /= 1000.d0) then

        tmprealj(:)=0.d0
        do jw=nbnd_start,nbnd_end
           tmprealj(:)=tmprealj(:)+eigenvector(jw-nbnd_start+1,iw-nbnd_start+1)*tmpreali(:,jw-nbnd_start+1)
        enddo

        if(okvan) then
           becp_gw2(:,iw)=0.d0
           do jw=nbnd_start,nbnd_end
              if(itask==0) then
                 becp_gw2(:,iw)=becp_gw2(:,iw)+eigenvector(jw-nbnd_start+1,iw-nbnd_start+1)*becp_gw(:,jw)
              else
                 becp_gw2(:,iw)=becp_gw2(:,iw)+eigenvector(jw-nbnd_start+1,iw-nbnd_start+1)*becp_gw_c(:,jw)
              endif
           enddo
        endif
     else
        tmprealj(:)=tmpreali(:,iw-nbnd_start+1)
        if(okvan) then
           if(itask==0) then
              becp_gw2(:,iw)=becp_gw(:,iw)
           else
              becp_gw2(:,iw)=becp_gw_c(:,iw)
           endif
        endif
     endif

      tmp_s(:)=tmprealj(:)*tmprealj(:)
      if(doublegrid) then
         call interpolate(tmp_r,tmp_s,1)
      else
         tmp_r(:)=tmp_s(:)
      endif
      if(okvan) call adduspos_gamma_r(iw,iw,tmp_r,1,becp_gw2(:,iw),becp_gw2(:,iw))


!determines integer coordinates of center of wannier wfcs

      nc1=nint(center(1,iqq)/(at(1,1))*dble(dfftp%nr1))
      nc2=nint(center(2,iqq)/(at(2,2))*dble(dfftp%nr2))
      nc3=nint(center(3,iqq)/(at(3,3))*dble(dfftp%nr3))

      nc1=nc1+1
      nc2=nc2+1
      nc3=nc3+1

      if(nc1<1) nc1=dfftp%nr1+nc1
      if(nc1>dfftp%nr1) nc1=nc1-dfftp%nr1

      if(nc2<1) nc2=dfftp%nr2+nc2
      if(nc2>dfftp%nr2) nc2=nc2-dfftp%nr2

      if(nc3<1) nc3=dfftp%nr3+nc3
      if(nc3>dfftp%nr3) nc3=nc3-dfftp%nr3


      write(stdout,*)'Wannier :', iw, 'Center :', nc1,nc2,nc3
      CALL flush_unit( stdout )
!determines integer cubic radius

      do nll=0,nrsmin/2
         norm=0.d0
         do ix=-nll,nll
            do iy=-nll,nll
               do iz=-nll,nll

                  n1=nc1+ix
                  if(n1<1) n1=dfftp%nr1+n1
                  if(n1>dfftp%nr1) n1=n1-dfftp%nr1

                  n2=nc2+iy
                  if(n2<1) n2=dfftp%nr2+n2
                  if(n2>dfftp%nr2) n2=n2-dfftp%nr2

                  n3=nc3+iz
                  if(n3<1) n3=dfftp%nr3+n3
                  if(n3>dfftp%nr3) n3=n3-dfftp%nr3

                  if(n3 >= nr3_start .and. n3 <= nr3_end) then
                     nn=(n3-nr3_start)*dfftp%nr1x*dfftp%nr2x+(n2-1)*dfftp%nr1x+n1
                     norm=norm+tmp_r(nn)
                  endif
               enddo
            enddo
         enddo

         norm=norm/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         call mp_sum(norm)
         if(norm >= cutoff) then
            exit
         endif
      enddo
      if(nll > nrsmin/2) nll=nrsmin/2
      write(stdout,*) 'Wannier :',iw,'Norm :',norm,'Nl :', nll
      CALL flush_unit( stdout )
!determines integer origin

      no1=nc1-nll
      if(no1<1) no1=dfftp%nr1+no1
      if(no1>dfftp%nr1) no1=no1-dfftp%nr1

      no2=nc2-nll
      if(no2<1) no2=dfftp%nr2+no2
      if(no2>dfftp%nr2) no2=no2-dfftp%nr2

      no3=nc3-nll
      if(no3<1) no3=dfftp%nr3+no3
      if(no3>dfftp%nr3) no3=no3-dfftp%nr3

!put on array
      if(doublegrid) then
         call interpolate(tmp_r,tmprealj,1)
      else
         tmp_r(:)=tmprealj(:)
      endif
!writes on file
!if itask==1 writes {c'} from position 1
      if(itask /= 1) then
         call  davcio( tmp_r,dfftp%nnr,iunrealwan2,iw,1)
      else
         call  davcio( tmp_r,dfftp%nnr,iunrealwan2,iw-n_first+1,1)
      endif

!put on arrays center and radius

      if(itask/=1) then
         w_centers(1,iw)=nc1
         w_centers(2,iw)=nc2
         w_centers(3,iw)=nc3
         w_radii(iw)=nll
      else
         w_centers_c(1,iw-n_first+1)=nc1
         w_centers_c(2,iw-n_first+1)=nc2
         w_centers_c(3,iw-n_first+1)=nc3
         w_radii_c(iw-n_first+1)=nll
      endif

!writes on file




   enddo!on iw
enddo!on iq



!update becp's
if(okvan) then
   if(itask==0) then
      becp_gw(:,nbnd_start:nbnd_end)=becp_gw2(:,nbnd_start:nbnd_end)
   else
      becp_gw_c(:,nbnd_start:nbnd_end)=becp_gw2(:,nbnd_start:nbnd_end)
   endif
endif

 close(iunrealwan2)


!update u_trans

 if(alfa/= 1000d0) then
    allocate(c_mat(nbnd_normal,nbnd_normal))
    c_mat(:,:)=(0.d0,0.d0)

    do iw=nbnd_start,nbnd_end
       do jw=nbnd_start,nbnd_end
          do kw=nbnd_start,nbnd_end
             c_mat(iw,jw)=c_mat(iw,jw)+eigenvector(kw-nbnd_start+1,iw-nbnd_start+1)*u_trans(kw,jw)
          enddo
       enddo
    enddo

    do iw=nbnd_start,nbnd_end
       do  jw=nbnd_start,nbnd_end
          u_trans(jw,iw)=c_mat(jw,iw)
       enddo
    enddo
    deallocate(c_mat)
 endif
 deallocate(work)
 if(alfa/=1000d0) deallocate(loc_mat)

enddo !on ivv

if(okvan) deallocate(becp_gw2)


 deallocate(eig_set,eigvector_set,eigvector_tmp)
 deallocate(tmpreal)
 deallocate(tmpreali,tmprealj)
 deallocate(eigenvector2,eigenvector_old)
 deallocate(iwork,ifail)
 deallocate(eigx,eigy,eigz)
 deallocate(min_1,min_2,max_1,max_2)
 deallocate(sums)
 deallocate(tmp_s,tmp_r)
 deallocate(exp_x,exp_y,exp_z)
!#endif

END SUBROUTINE ultralocalization_para

