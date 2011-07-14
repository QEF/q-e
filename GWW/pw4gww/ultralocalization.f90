! FOR GWW
!
! Author: P. Umari
!
SUBROUTINE ultralocalization(nbndv,ultra_thr,isubspace,max_array)
!
!this subroutine
!1-read R space wanniers
!2-ultralocalized
!3-determine cubic length
!4-calculates A matrix : \tilde{w}_i=A_{i,j}\Psi_j
!4-write on file

! #ifdef __GWW

  USE io_files,             ONLY : find_free_unit, diropn
  USE io_global,            ONLY : stdout
  USE fft_base,             ONLY : dffts
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE fft_base,             ONLY : dfftp
  use mp_global,            ONLY : nproc_pool, me_pool
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx
  USE basis
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE uspp,                 ONLY : okvan,nkb
  USE realus,               ONLY : adduspos_gamma_r

  implicit none

  REAL(kind=DP), INTENT(in) :: ultra_thr!threshold for convergence
  INTEGER, INTENT(in) :: nbndv !number of first bands wich are localized separately
  INTEGER, INTENT(in) :: isubspace! == 0 valence subspace; == 1 conduction subspace
  INTEGER, INTENT(in) :: max_array! max number of states to be considered in the same array

  REAL(kind=DP), ALLOCATABLE ::  tmpreal(:),tmpreali(:,:),tmprealj(:)
  INTEGER :: iw,jw,kw,ir,it, ii,jj
  REAL(kind=DP) :: sca,sca1
  COMPLEX(kind=DP) :: scac1,scac2,scac3
  INTEGER :: iunrealwan
  LOGICAL :: exst
  INTEGER :: ix,iy,iz,nn
  REAL(kind=DP) :: center(3,max_array),center_old(3,max_array),rad, radmax
  REAL(kind=DP) :: rdistance
  REAL(kind=DP) :: rx,ry,rz
  REAL(kind=DP), ALLOCATABLE :: loc_mat(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: eigenvector(:,:),eigenvector_old(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: eigenvector2(:)
  INTEGER, PARAMETER :: max_loc_ite=5!1 ATTENZIONE
  INTEGER :: lwork, m,info
  REAL(kind=DP), ALLOCATABLE :: work(:)
  REAL(kind=DP) :: eigen(max_array),eigenold(max_array)
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
  LOGICAL :: converged(nbnd)
  REAL(kind=DP), ALLOCATABLE :: sums(:,:,:)
  INTEGER :: iq,ifirst,ilast,iqq
  COMPLEX(kind=DP), ALLOCATABLE :: c_mat(:,:)
  REAL(kind=DP) :: cutoff
  REAL(kind=DP) :: alfa
  REAL(kind=DP), ALLOCATABLE :: becp_gw2(:,:)
  REAL(kind=DP), ALLOCATABLE :: tmp_s(:),tmp_r(:)
  LOGICAL :: is_even
  INTEGER, PARAMETER :: nmaxeig = 1!4 ATTENZIONE
  REAL(kind=DP), PARAMETER :: minover = 0.3d0

  REAL(kind=DP), ALLOCATABLE :: eig_set(:), eigvector_tmp(:)
  REAL(kind=DP), ALLOCATABLE :: eigvector_set(:,:)
  INTEGER :: itt, ipp
  LOGICAL :: good

  allocate(eig_set(nmaxeig))
  allocate(eigvector_set(nbnd,nmaxeig),eigvector_tmp(nbnd))


  allocate(tmpreal(dfftp%nnr))
  allocate(eigenvector2(nbnd),eigenvector_old(nbnd,max_array))
  allocate(iwork(5*nbnd),ifail(nbnd))
  allocate(eigx(nbnd,nbnd),eigy(nbnd,nbnd),eigz(nbnd,nbnd))
  allocate(exp_x(dfftp%nnr),exp_y(dfftp%nnr),exp_z(dfftp%nnr))
  allocate(sums(dfftp%nr1,dfftp%nr2,dfftp%nr3))
  allocate(tmp_s(dffts%nnr),tmp_r(dfftp%nnr))
  if(okvan) allocate(becp_gw2(nkb,nbnd))

  if(isubspace==0) then
    n_first=1
    n_last=nbndv
    n_bands=nbndv
    cutoff=cutoff_wsq
  else if(isubspace==1) then
    n_first=nbndv+1
    n_last=nbnd
    n_bands=nbnd-nbndv
    cutoff=cutoff_wsq_c
  else
    write(stdout,*) 'ultralocalization isubspace ILLEGAL'
    stop
 endif

 alfa=0.01d0

 allocate(tmpreali(dffts%nnr,n_bands))
 allocate(tmprealj(dffts%nnr))
 allocate(min_1(dfftp%nr2,dfftp%nr3,n_bands),min_2(dfftp%nr2,dfftp%nr3,n_bands))
 allocate(max_1(dfftp%nr2,dfftp%nr3,n_bands),max_2(dfftp%nr2,dfftp%nr3,n_bands))
 allocate(eigenvector(nbnd,n_bands))

  radmax=no_radius/alat

  nrsmin=min(dfftp%nr1,dfftp%nr2)
  nrsmin=min(dfftp%nr3,nrsmin)
!if nrsmin is even set to nrsmin -1,
  if(is_even(nrsmin)) then
    nrsmin=nrsmin-1
    write(stdout,*) 'nrsmin set to:', nrsmin
  endif

  write(*,*)' ULTRALOCALIZATION' !ATTENZIONE
!open outpu file

   iunrealwan = find_free_unit()
   CALL diropn( iunrealwan, 'real_whole', dffts%nnr, exst )
!read in wave-functions
  do iw=n_first,n_last
    CALL davcio( tmpreali(:,iw-n_first+1),dffts%nnr,iunrealwan,iw,-1)
  enddo
   CLOSE(iunrealwan)
!first valence subspace

  allocate(loc_mat(n_bands,n_bands,max_array))
  lwork=8*nbnd
  allocate(work(lwork))

  exp_x(:)=(0.d0,0.d0)
  exp_y(:)=(0.d0,0.d0)
  exp_z(:)=(0.d0,0.d0)
  do ix=1,dffts%nr1
     do iy=1,dffts%nr2
        do iz=1,dffts%nr3
           nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
           exp_x(nn)=exp((0.d0,1.d0)*tpi*real(ix)/real(dffts%nr1))
           exp_y(nn)=exp((0.d0,1.d0)*tpi*real(iy)/real(dffts%nr2))
           exp_z(nn)=exp((0.d0,1.d0)*tpi*real(iz)/real(dffts%nr3))
        enddo
     enddo
  enddo

  do iw=n_first,n_last
    do jw=iw,n_last
       iww=iw-n_first+1
       jww=iw-n_first+1
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
          do iz=1,dfftp%nr3
             nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
             exp_x(nn)=exp((0.d0,1.d0)*tpi*real(ix)/real(dfftp%nr1))
             exp_y(nn)=exp((0.d0,1.d0)*tpi*real(iy)/real(dfftp%nr2))
             exp_z(nn)=exp((0.d0,1.d0)*tpi*real(iz)/real(dfftp%nr3))
          enddo
       enddo
    enddo

    do iw=n_first,n_last
      do jw=iw,n_last
         iww=iw-n_first+1
         jww=iw-n_first+1
         tmp_r(:)=0.d0
         call adduspos_gamma_r(iw,jw,tmp_r,1,becp_gw(:,iw),becp_gw(:,jw))
         scac1=(0.d0,0.d0)
         scac2=(0.d0,0.d0)
         scac3=(0.d0,0.d0)
         do ir=1,dfftp%nnr
            scac1=scac1+exp_x(ir)*tmp_r(ir)
            scac2=scac2+exp_y(ir)*tmp_r(ir)
            scac3=scac3+exp_z(ir)*tmp_r(ir)
          enddo
          eigx(iw,jw)=eigx(iw,jw)+scac1/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
          eigy(iw,jw)=eigy(iw,jw)+scac2/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
          eigz(iw,jw)=eigz(iw,jw)+scac3/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
          eigx(jw,iw)=eigx(iw,jw)
          eigy(jw,iw)=eigy(iw,jw)
          eigz(jw,iw)=eigz(iw,jw)
        enddo
      enddo
    endif
  deallocate(exp_x,exp_y,exp_z)


  do iq=1,n_bands/max_array
     ifirst=(iq-1)*max_array+n_first
     ilast=min(iq*max_array+n_first,n_last)
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
     do it=1,max_loc_ite

 !calculate screen function with center and radius
     iqq=0
     do iw=ifirst,ilast
        iqq=iqq+1
        if(.not.converged(iqq)) then

          do iy=1,dffts%nr2
              do iz=1,dffts%nr3
                min_1(iy,iz,iqq)=0
                max_1(iy,iz,iqq)=0
                min_2(iy,iz,iqq)=0
                max_2(iy,iz,iqq)=0
                do ix=1,dffts%nr1
                  nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
                  rx=rdistance(real(ix)*at(1,1)/real(dffts%nr1),center(1,iqq),at(1,1))
                  ry=rdistance(real(iy)*at(2,2)/real(dffts%nr2),center(2,iqq),at(2,2))
                  rz=rdistance(real(iz)*at(3,3)/real(dffts%nr3),center(3,iqq),at(3,3))
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
                    nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
                    rx=rdistance(real(ix)*at(1,1)/real(dffts%nr1),center(1,iqq),at(1,1))
                    ry=rdistance(real(iy)*at(2,2)/real(dffts%nr2),center(2,iqq),at(2,2))
                    rz=rdistance(real(iz)*at(3,3)/real(dffts%nr3),center(3,iqq),at(3,3))
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
      do jw=n_first,n_last
        jww=jw-n_first+1
        do kw=jw,n_last
          kww=kw-n_first+1
          sums(:,:,:)=0.d0
          tmp_s(:)=tmpreali(:,jww)*tmpreali(:,kww)
          do iy=1,dffts%nr2
            do iz=1,dffts%nr3
              sca=0.d0
              do ix=1,dffts%nr1
                nn=(iz-1)*dffts%nr1x*dffts%nr2x+(iy-1)*dffts%nr1x+ix
                sca=sca+tmp_s(nn)
                sums(ix,iy,iz)=sca
              enddo
            enddo
          enddo
          iqq=0
          do iw=ifirst,ilast
            iqq=iqq+1
            if(.not.converged(iqq))then
              do iy=1,dffts%nr2
                do iz=1,dffts%nr3
                   if(max_1(iy,iz,iqq)/=0) then
                      if(min_1(iy,iz,iqq)/=1) then
                        loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq)+&
                                 &sums(max_1(iy,iz,iqq)-1,iy,iz)-sums(min_1(iy,iz,iqq)-1,iy,iz)
                      else
                        loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq) + sums(max_1(iy,iz,iqq)-1,iy,iz)
                      endif
                   endif
                  if(max_2(iy,iz,iqq)/=0) then
                      loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq)+sums(max_2(iy,iz,iqq),iy,iz)-sums(min_2(iy,iz,iqq),iy,iz)
                   endif
                enddo
              enddo
              loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq)/real(dffts%nr1*dffts%nr2*dffts%nr3)
              loc_mat(kww,jww,iqq)=loc_mat(jww,kww,iqq)
            endif
         enddo
       enddo
    enddo



    if(okvan) then
      iqq=0
       do iw=ifirst,ilast
          iqq=iqq+1
          if(.not.converged(iqq)) then

            do iy=1,dfftp%nr2
                do iz=1,dfftp%nr3
                  min_1(iy,iz,iqq)=0
                  max_1(iy,iz,iqq)=0
                  min_2(iy,iz,iqq)=0
                  max_2(iy,iz,iqq)=0
                  do ix=1,dfftp%nr1
                    nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                    rx=rdistance(real(ix)*at(1,1)/real(dfftp%nr1),center(1,iqq),at(1,1))
                    ry=rdistance(real(iy)*at(2,2)/real(dfftp%nr2),center(2,iqq),at(2,2))
                    rz=rdistance(real(iz)*at(3,3)/real(dfftp%nr3),center(3,iqq),at(3,3))
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
                      nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                      rx=rdistance(real(ix)*at(1,1)/real(dfftp%nr1),center(1,iqq),at(1,1))
                      ry=rdistance(real(iy)*at(2,2)/real(dfftp%nr2),center(2,iqq),at(2,2))
                      rz=rdistance(real(iz)*at(3,3)/real(dfftp%nr3),center(3,iqq),at(3,3))
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
        do jw=n_first,n_last
          jww=jw-n_first+1
          do kw=jw,n_last
            kww=kw-n_first+1
            sums(:,:,:)=0.d0
            tmp_r(:)=0.d0
            call adduspos_gamma_r(jw,kw,tmp_r,1,becp_gw(:,jw),becp_gw(:,kw))
            do iy=1,dfftp%nr2
              do iz=1,dfftp%nr3
                sca=0.d0
                do ix=1,dfftp%nr1
                  nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
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
                do iy=1,dfftp%nr2
                  do iz=1,dfftp%nr3
                     if(max_1(iy,iz,iqq)/=0) then
                        if(min_1(iy,iz,iqq)/=1) then
                          loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq)+&
                                   &sums(max_1(iy,iz,iqq)-1,iy,iz)-sums(min_1(iy,iz,iqq)-1,iy,iz)
                        else
                          loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq) + sums(max_1(iy,iz,iqq)-1,iy,iz)
                        endif
                     endif
                      if(max_2(iy,iz,iqq)/=0) then
                        loc_mat(jww,kww,iqq)=loc_mat(jww,kww,iqq)+sums(max_2(iy,iz,iqq),iy,iz)-sums(min_2(iy,iz,iqq),iy,iz)
                     endif
                  enddo
                enddo
                loc_mat(kww,jww,iqq)=loc_mat(jww,kww,iqq)
              endif
           enddo
         enddo
      enddo
    endif











!loop on ultra-localized wanniers
      iqq=0
      do iw=ifirst,ilast
        iqq=iqq+1
        if(.not.converged(iqq)) then
          loc_mat(iw-n_first+1,iw-n_first+1,iqq) = loc_mat(iw-n_first+1,iw-n_first+1,iqq) + alfa
          if(it==1)  eigenold(iqq)=loc_mat(iw-n_first+1,iw-n_first+1,iqq)

!find upperstate

          call dsyevx('V','I','U',n_bands,loc_mat(:,:,iqq),n_bands,0.d0,0.d0,n_bands-nmaxeig+1,n_bands,0.d0,m,eig_set,&
   & eigvector_set,nbnd,work,lwork,iwork,ifail,info)
          do ipp=nmaxeig,1,-1
             good=.true.
             eigvector_tmp(:)=eigvector_set(:,ipp)
             sca=0.d0
             do  ii=1,n_bands
               sca=sca+eigvector_tmp(ii)**2.d0
             enddo
             write(*,*) 'Modulus vector', sca
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

             if(sca >=   minover) exit
          enddo
          !if(ipp<nmaxeig) write(*,*) 'Ultraloc too close', iqq,ipp
           write(*,*) 'Ultraloc', iqq,ipp, sca
          if(ipp==0) ipp=1
          eigen(iqq)=eig_set(ipp)
          eigenvector(:,iw-n_first+1)=eigvector_set(:,ipp)
          if(info/=0) then
             write(stdout,*) 'Error from dsyevx: ', info
             stop
          endif
          write(stdout,*) 'Ultra:',iw,it,eigen(iqq),eigenold(iqq)!ATTENZIONE
          if(eigen(iqq)<eigenold(iqq).and.it>1) then
            write(stdout, *) 'Increasing..stopping'
            eigenvector(:,iw-n_first+1)=eigenvector_old(:,iqq)
            center(:,iqq)=center_old(:,iqq)
            converged(iqq)=.true.
          endif
          if(abs(eigen(iqq)-eigenold(iqq))<= ultra_thr) converged(iqq)=.true.
          if( it == max_loc_ite) converged(iqq)=.true.
          if(.not.converged(iqq)) then
            eigenold(iqq)=eigen(iqq)
            eigenvector_old(:,iqq)=eigenvector(:,iw-n_first+1)
            center_old(:,iqq)=center(:,iqq)


!update center
            do ii=n_first,n_last
              eigenvector2(ii)=(0.d0,0.d0)
              do jj=n_first,n_last
                eigenvector2(ii)=eigenvector2(ii)+eigx(ii,jj)*eigenvector(jj-n_first+1,iw-n_first+1)
              enddo
            enddo
            center_berry=(0.d0,0.d0)
            do ii=n_first,n_last
              center_berry=center_berry+eigenvector(ii-n_first+1,iw-n_first+1)*eigenvector2(ii)
            enddo
            center(1,iqq)=aimag(log(center_berry))*at(1,1)/tpi

            do ii=n_first,n_last
              eigenvector2(ii)=(0.d0,0.d0)
              do jj=n_first,n_last
                eigenvector2(ii)=eigenvector2(ii)+eigy(ii,jj)*eigenvector(jj-n_first+1,iw-n_first+1)
              enddo
            enddo
            center_berry=(0.d0,0.d0)
            do ii=n_first,n_last
              center_berry=center_berry+eigenvector(ii-n_first+1,iw-n_first+1)*eigenvector2(ii)
            enddo
            center(2,iqq)=aimag(log(center_berry))*at(2,2)/tpi

            do ii=n_first,n_last
               eigenvector2(ii)=(0.d0,0.d0)
              do jj=n_first,n_last
                eigenvector2(ii)=eigenvector2(ii)+eigz(ii,jj)*eigenvector(jj-n_first+1,iw-n_first+1)
              enddo
            enddo
            center_berry=(0.d0,0.d0)
            do ii=n_first,n_last
              center_berry=center_berry+eigenvector(ii-n_first+1,iw-n_first+1)*eigenvector2(ii)
            enddo
           center(3,iqq)=aimag(log(center_berry))*at(3,3)/tpi
         endif
       endif
        enddo!on iw
      enddo!on iterations it
!write on file
   iqq=0
    do iw=ifirst,ilast
      iqq=iqq+1
      write(stdout,*) 'Ultralocalization overlap with starting wannier:', iw
!construct ultralocalized wannier in real space

      tmprealj(:)=0.d0
      do jw=n_first,n_last
        tmprealj(:)=tmprealj(:)+eigenvector(jw-n_first+1,iw-n_first+1)*tmpreali(:,jw-n_first+1)
      enddo

      if(okvan) then
        becp_gw2(:,iw)=0.d0
        do jw=n_first,n_last
          becp_gw2(:,iw)=becp_gw2(:,iw)+eigenvector(jw-n_first+1,iw-n_first+1)*becp_gw(:,jw)
        enddo
      endif

      tmp_s(:)=tmprealj(:)*tmprealj(:)
      if(doublegrid) then
        call interpolate(tmp_r,tmp_s,1)
      else
        tmp_r(:)=tmp_s(:)
      endif
      if(okvan) call adduspos_gamma_r(iw,iw,tmp_r,1,becp_gw2(:,iw),becp_gw2(:,iw))


!determines integer coordinates of center of wannier wfcs

      nc1=aint(center(1,iqq)/(at(1,1))*real(dfftp%nr1))
      nc2=aint(center(2,iqq)/(at(2,2))*real(dfftp%nr2))
      nc3=aint(center(3,iqq)/(at(3,3))*real(dfftp%nr3))


      if(nc1<1) nc1=dfftp%nr1+nc1
      if(nc1>dfftp%nr1) nc1=nc1-dfftp%nr1

      if(nc2<1) nc2=dfftp%nr2+nc2
      if(nc2>dfftp%nr2) nc2=nc2-dfftp%nr2

      if(nc3<1) nc3=dfftp%nr3+nc3
      if(nc3>dfftp%nr3) nc3=nc3-dfftp%nr3


      write(stdout,*)'Wannier :', iw, 'Center :', nc1,nc2,nc3

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

                 nn=(n3-1)*dfftp%nr1x*dfftp%nr2x+(n2-1)*dfftp%nr1x+n1
                 norm=norm+tmp_r(nn)
               enddo
             enddo
          enddo
        norm=norm/real(dfftp%nr1*dfftp%nr2*dfftp%nr3)
        if(norm >= cutoff) then
          exit
        endif
      enddo
      if(nll > nrsmin/2) nll=nrsmin/2
          write(stdout,*) 'Wannier :',iw,'Norm :',norm,'Nl :', nll
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

                do ix=0,2*nll
                  do iy=0,2*nll
                    do iz=0,2*nll

                      n1=no1+ix
                      if(n1<1) n1=dfftp%nr1+n1
                      if(n1>dfftp%nr1) n1=n1-dfftp%nr1

                      n2=no2+iy
                      if(n2<1) n2=dfftp%nr2+n2
                      if(n2>dfftp%nr2) n2=n2-dfftp%nr2

                      n3=no3+iz
                      if(n3<1) n3=dfftp%nr3+n3
                      if(n3>dfftp%nr3) n3=n3-dfftp%nr3
                      nn=(n3-1)*dfftp%nr1x*dfftp%nr2x+(n2-1)*dfftp%nr1x+n1
                     tmpreal(iz*(2*nll+1)*(2*nll+1)+iy*(2*nll+1)+ix+1)=&
                &  tmp_r(nn)
                  enddo
                enddo
              enddo

!put on arrays center and radius

              w_centers(1,iw)=nc1
              w_centers(2,iw)=nc2
              w_centers(3,iw)=nc3
              w_radii(iw)=nll


!writes on file


              write(nfile,'(5i1)') iw/10000,mod(iw,10000)/1000,mod(iw,1000)/100,mod(iw,100)/10,mod(iw,10)
              open(unit=iunrealwan2,file='realwan'// nfile,status='unknown', form='unformatted')
              write(iunrealwan2) no1,no2,no3,nll
              write(iunrealwan2) tmpreal(1:(2*nll+1)**3)
              close(iunrealwan2)

  enddo!on iw
  enddo!on iq

!update becp's
  becp_gw(:,n_first:n_last)=becp_gw2(:,n_first:n_last)

  deallocate(eig_set,eigvector_set,eigvector_tmp)
  deallocate(loc_mat)
  deallocate(tmpreal)
  deallocate(tmpreali,tmprealj)
  deallocate(eigenvector2,eigenvector_old)
  deallocate(work,iwork,ifail)
  deallocate(eigx,eigy,eigz)
  deallocate(min_1,min_2,max_1,max_2)
  deallocate(sums)
  deallocate(tmp_s,tmp_r)
  if(okvan) deallocate(becp_gw2)

!update u_trans

  allocate(c_mat(nbnd,nbnd))
  c_mat(:,:)=(0.d0,0.d0)

  do iw=n_first,n_last
    do jw=n_first,n_last
       do kw=n_first,n_last
          c_mat(iw,jw)=c_mat(iw,jw)+eigenvector(kw-n_first+1,iw-n_first+1)*u_trans(kw,jw)
       enddo
    enddo
  enddo

  do iw=n_first,n_last
    do  jw=n_first,n_last
       u_trans(jw,iw)=c_mat(jw,iw)
    enddo
  enddo

  deallocate(eigenvector,c_mat)

! #endif
END SUBROUTINE ultralocalization


FUNCTION rdistance(r1,r2,rl)
!find minimum distance between n1 and n2 with periodicity l

  USE kinds, ONLY : dp


  implicit none

  REAL(kind=DP)  :: rdistance
  REAL(kind=DP)  :: r1,r2,rl
  INTEGER :: i,  imin
  REAL(kind=DP) :: rmin

  rmin=2.d0*rl
  do i=-1,1
     if(abs(r1-(r2+real(i)*rl)) < rmin) then
       rmin = abs(r1-(r2+real(i)*rl))
       imin= i
     endif
  enddo
  rdistance = r1-(r2+real(imin)*rl)

END FUNCTION rdistance


FUNCTION is_even(n)
!true, if n is even, false otherwise

   implicit none

   LOGICAL :: is_even
   INTEGER :: n

   if(mod(n,2)==0) then
     is_even = .true.
   else
     is_even = .false.
   endif

END FUNCTION is_even

