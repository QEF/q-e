! FOR GWW
!
! Author: P. Umari
!
 subroutine product_wannier_para_c(nbndv, lcomplete, lrestart)

!this subroutine
!1-for every ordered couple of wannier wfcs iw,jw
! -formed by conduction prim  and conductions wanniers
!2-determine if they overlap
!3-if yes, calculate product in r space
!4-determine center and radius of overlap 3 lengths
!5-write product on R space on disk
!6-go to G space and save on disk
!7-put on array, which is enlarged

!this version writes the ultralocalized wannier functions
!on the big grid R, in order to use davcio for faster i/o
!on parallel architectures



! #ifdef __GWW

  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : find_free_unit
  USE io_files,             ONLY : nwordwfc
  USE io_files,             ONLY : prefix
  USE io_files,             ONLY : tmp_dir, iunwfc, iunigk, diropn
  USE io_global,            ONLY : stdout, ionode
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  use mp_global,            ONLY : nproc_pool, me_pool
  USE kinds,                ONLY : DP
  USE us
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbndx, ecutwfc
  USE control_flags,        ONLY: gamma_only
  USE gvect
  USE basis
  USE klist
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE ions_base,            ONLY : ityp, tau, nat,ntyp => nsp
  USE uspp,                 ONLY : okvan
  USE wannier_gw
  USE realus,               ONLY : adduspos_gamma_r
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft
  USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum


  implicit none

  INTEGER, INTENT(in) :: nbndv!number of valence subspace bands (and wannier wfcs)
  LOGICAL, INTENT(in) :: lcomplete!if .true. consider all space
  LOGICAL, INTENT(in) :: lrestart!if .true. is a restart calculation and the data on the centers of wannier are
                                 ! not available


  INTEGER :: iunwannier,iungprod,iunprod,iunrprod, iuncouples

  !  --- Internal definitions ---

   REAL(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   REAL(kind=DP), ALLOCATABLE :: tmpspacej(:),tmpspacejs(:)
   COMPLEX(kind=DP), ALLOCATABLE :: warray(:) !temporary array for i/o
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacec(:)
   REAL(kind=DP), ALLOCATABLE :: tmpreal(:),tmpreal2(:)!temporary reading array
   INTEGER :: i,j,k,ig,jg
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   INTEGER :: nt,na, ih, jh, ikb, jkb, ijkb0
   INTEGER :: nmax
   INTEGER :: nll,no1,no2,no3,ndist,nop(3)
   INTEGER :: n1,n2,n3,nn
   REAL(kind=DP) :: norm
   INTEGER :: ix, iy, iz
   INTEGER :: iunrealwan, iunrealwanprim
   CHARACTER(5) :: nfilei,nfilej
   COMPLEX(kind=DP) :: sca,sca1
   INTEGER :: max_nx,max_ny,max_nz
   REAL(kind=DP) :: maxsca,rsca
   INTEGER :: iw,jw
   INTEGER :: rspacel(3)
   LOGICAL :: match
   TYPE(wannier_product) :: wiwj
   INTEGER :: ndistance
   LOGICAL :: exst
   LOGICAL :: is_even

   INTEGER :: nr3_start, nr3_end

   REAL(kind=DP) :: sum_vcvc!for consistency control

   REAL(kind=DP) :: cutoff_product


#ifndef __PARA
  !dfftp%npp(1) = nr3 ! no longer needed - PG
  nr3_start=1
  nr3_end=dfftp%nr3
#else

  nr3_start=0
  nr3_end =0
  do i=1,me_pool + 1
     nr3_start=nr3_end+1
     nr3_end=nr3_end+dfftp%npp(i)
  end do
#endif



  write(stdout,*) 'Routine product_wannier_para: start', lsmallgrid

  if(okvan .and. lsmallgrid) write(stdout,*) 'ATTENTION: USPP AND SMALLGRID'

  allocate(tmpspacei(dfftp%nnr,num_nbndc_set),tmpspacej(dfftp%nnr),tmpreal(dfftp%nnr),tmpreal2(dfftp%nnr))
  allocate(tmpspacejs(dffts%nnr),tmpspacec(dfftp%nnr))

   numw_prod_c=0



   rspacel(1)=dfftp%nr1
   rspacel(2)=dfftp%nr2
   rspacel(3)=dfftp%nr3

   !reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)

   if(.not.lsmallgrid) then
     allocate (warray(ngm))
   else
     allocate (warray(npw0))
   endif


! open files for ouput product in g space and their centers and radii


   iungprod = find_free_unit()
   if(.not.lsmallgrid) then
      CALL diropn( iungprod, 'wiwjwfc_prim', max_ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc_prim', npw0*2, exst )
   endif
   iunrealwan =  find_free_unit()
   CALL diropn( iunrealwan, 'realwan', dfftp%nnr, exst )
   iunrealwanprim =  find_free_unit()
   CALL diropn( iunrealwanprim, 'realwan_prim', dfftp%nnr, exst )
   if(ionode) then
      iunprod = find_free_unit()
      open( unit= iunprod, file=trim(prefix)//'.wiwjprod_prim', status='unknown',form='unformatted')
      iuncouples =  find_free_unit()
      open( unit= iuncouples, file=trim(prefix)//'.ij_prim', status='unknown',form='unformatted')
   endif


!open output file


   sum_vcvc = 0.d0


   do iw=1,num_nbndc_set
!read real wfcs iw
      call davcio( tmpspacei(:,iw),dfftp%nnr,iunrealwanprim,iw,-1)
   enddo

   do jw=nbnd_normal,nbndv+1,-1
!determine cutoff
      if(num_nbnd_first==0) then
         cutoff_product=cutoff_wpr_prim
      else
         if(jw > nbndv+num_nbnd_first) then
            cutoff_product=cutoff_wpr_prim2
         else
            cutoff_product=cutoff_wpr_prim
         endif
      endif

      do iw=num_nbndc_set,1,-1
!determines bottoms and tops of boxes in grid units

         match = .true.

         if(.not.lrestart .and. .not.lcomplete) then

           ! write(stdout,*) 'Center iw', iw, w_centers_c(1:3,iw)
           ! write(stdout,*) 'Center jw', jw, w_centers(1:3,jw)


            do i=1,3
               ndist=ndistance(w_centers_c(i,iw),w_centers(i,jw),rspacel(i))
               if(ndist-w_radii_c(iw)-w_radii(jw) > 0 ) then
                  match = .false.
                  exit
               endif
            enddo
         endif

         if( match) then

           if(ionode) write(stdout,*) 'Matching :', iw,jw
!determine center of product and radii
            if(.not.lrestart) then
               do i=1,3
                  call noverlap(w_centers_c(i,iw),w_centers(i,jw),w_radii_c(iw),w_radii(jw),&
                       &wiwj%center(i),wiwj%radius(i),rspacel(i))
                  if(is_even(rspacel(i)).and.wiwj%radius(i)==(rspacel(i)/2)) &
                       & wiwj%radius(i)=wiwj%radius(i)-1
               enddo

             !  write(stdout,*) 'Center:', wiwj%center(1:3),wiwj%radius(1:3)
            endif
!read wannier function jw

            call davcio( tmpspacej,dfftp%nnr,iunrealwan,jw,-1)



         !   write(stdout,*) 'read wannier'

! calculates product in R space

!determines origin
            if(.not.lrestart) then
               nop(:)= wiwj%center(:)-wiwj%radius(:)
               do i=1,3
                  if(nop(i)<1) nop(i)=rspacel(i)+nop(i)
                  if(nop(i)>rspacel(i)) nop(i)=nop(i)-rspacel(i)
               enddo



            endif

        tmpreal(:)=0.d0
!adds US term if required
        tmpreal2(:)=0.d0
        if(okvan) call adduspos_gamma_r(iw+nbndv,jw,tmpreal2,1,becp_gw_c(:,iw+nbndv),becp_gw(:,jw))




        rsca=0.d0

        if(.not.lcomplete.and..not.lrestart) then
           do ix=0,2*wiwj%radius(1)
              do iy=0,2*wiwj%radius(2)
                 do iz=0,2*wiwj%radius(3)
                    n1=nop(1)+ix
                    if(n1<1) n1=dfftp%nr1+n1
                    if(n1>dfftp%nr1) n1=n1-dfftp%nr1

                    n2=nop(2)+iy
                    if(n2<1) n2=dfftp%nr2+n2
                    if(n2>dfftp%nr2) n2=n2-dfftp%nr2

                    n3=nop(3)+iz
                    if(n3<1) n3=dfftp%nr3+n3
                    if(n3>dfftp%nr3) n3=n3-dfftp%nr3

                    if(n3 >= nr3_start .and. n3 <= nr3_end) then
                       nn=(n3-nr3_start)*dfftp%nr1x*dfftp%nr2x+(n2-1)*dfftp%nr1x+n1
                       if(nn<1 .or. nn > dfftp%nnr)  then
                          CALL errore( 'rsca', 'rsca', nn )
                       endif
                       rsca=rsca+(tmpspacei(nn,iw)*tmpspacej(nn)+tmpreal2(nn))**2.d0
                       tmpreal(nn)=&
                            &tmpspacei(nn,iw)*tmpspacej(nn)+tmpreal2(nn)
                    endif
                 enddo
              enddo
           enddo
        else
           do nn=1,dfftp%nnr
              rsca=rsca+(tmpspacei(nn,iw)*tmpspacej(nn)+tmpreal2(nn))**2.d0
              tmpreal(nn)=&
                   &tmpspacei(nn,iw)*tmpspacej(nn)+tmpreal2(nn)
           enddo
        endif


         rsca=rsca/real(rspacel(1)*rspacel(2)*rspacel(3))


         call mp_sum(rsca)

!check if there is overlap
         if(rsca >= cutoff_product) then
         if(ionode)   write(stdout,*) 'OVERLAP :', iw,jw,rsca


            numw_prod_c=numw_prod_c+1

            wiwj%i=iw
            wiwj%j=jw

            if(ionode) write(iunprod) wiwj
            if(ionode) write(iuncouples) iw, jw


!put on grid for Fourier transform


!bring back to G space on charge grid (4times denser in NC case)
!maybe wfcs grid could be enough?

          if(.not.lsmallgrid) then
            tmpspacec(:)=dcmplx(tmpreal(:),0.d0)
            !ATTENZIONE
            rsca=0.d0
            do nn=1,dfftp%nnr
              rsca=rsca+conjg( tmpspacec(nn))*tmpspacec(nn)
            enddo
            call mp_sum(rsca)
            CALL fwfft ('Dense', tmpspacec, dfftp)
             rsca=0.d0
            do nn=1,ngm
              rsca=rsca+2.d0*conjg( tmpspacec(nl(nn)))*tmpspacec(nl(nn))
            enddo
            if(gamma_only .and. gstart == 2) then
               rsca=rsca-conjg( tmpspacec(nl(1)))*tmpspacec(nl(1))
            endif
            call mp_sum(rsca)
            if(jw>nbndv) sum_vcvc =sum_vcvc+rsca
          else
!if smallgrid and doublegrid, a interpolation to smallgrid is required
            if(lsmallgrid .and. doublegrid) then
              CALL interpolate (tmpspacej, tmpspacejs, -1)
            else
              tmpspacejs(:)=tmpspacej(:)
            endif
            tmpspacec(:)=dcmplx(tmpspacejs(:),0.d0)
            CALL fwfft ('Wave', tmpspacec, dffts)
          endif


!writes on file
          if(.not.lsmallgrid) then
            warray(1:max_ngm)=tmpspacec(nl(1:max_ngm))
            CALL davcio(warray,max_ngm*2,iungprod,numw_prod_c,1)
          else
            warray(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
            CALL davcio(warray,npw0*2,iungprod,numw_prod_c,1)
          endif
          call mp_barrier
        if(ionode)  write(stdout,*) 'written', numw_prod_c

          endif!overlap




        endif!match



     enddo!jw
  enddo!iw


  write(stdout,*) 'SUM VCVC', sum_vcvc

  close(iungprod)
  if(ionode) close(iunprod)
  if(ionode) close(iuncouples)
  close(iunrealwan)
  close(iunrealwanprim)
  deallocate(tmpspacei)
  deallocate(tmpspacej)
  deallocate(tmpreal)
  deallocate(tmpreal2)
  deallocate(tmpspacejs)
  deallocate(tmpspacec)
  deallocate(warray)

! #endif
  return
end subroutine product_wannier_para_c



