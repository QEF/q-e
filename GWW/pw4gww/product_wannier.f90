! FOR GWW
!
! Author: P. Umari
!
SUBROUTINE product_wannier(nbndv)
!
!this subroutine
!1-for every ordered couple of wannier wfcs iw,jw
! -formed by valence times valence and conductions wanniers
!2-determine if they overlap
!3-if yes, calculate product in r space
!4-determine center and radius of overlap 3 lengths
!5-write product on R space on disk
!6-go to G space and save on disk
!7-put on array, which is enlarged

! #ifdef __GWW

  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : find_free_unit
  USE io_files,             ONLY : nwordwfc
  USE io_files,             ONLY : prefix
  USE io_files,             ONLY : tmp_dir, iunwfc, iunigk, diropn
  USE io_global,            ONLY : stdout
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  use mp_global,            ONLY : nproc_pool, me_pool
  USE kinds,                ONLY : DP
  USE us
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE gvect
  USE basis
  USE klist
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE ions_base,            ONLY : ityp, tau, nat,ntyp => nsp
  USE uspp,                 ONLY : okvan
  USE wannier_gw
  USE realus,               ONLY : adduspos_gamma_r

  implicit none

  INTEGER, INTENT(in) :: nbndv!number of valence subspace bands (and wannier wfcs)


  INTEGER :: iunwannier,iungprod,iunprod,iunrprod

  !  --- Internal definitions ---

   REAL(kind=DP), ALLOCATABLE :: tmpspacei(:)
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
   INTEGER :: iunrealwan
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

   REAL(kind=DP) :: sum_vcvc!for consistency control


   write(stdout,*) 'Routine product_wannier: start'

   if(okvan .and. lsmallgrid) write(stdout,*) 'ATTENTION: USPP AND SMALLGRID'

   allocate(tmpspacei(dfftp%nnr),tmpspacej(dfftp%nnr),tmpreal(dfftp%nnr),tmpreal2(dfftp%nnr))
   allocate(tmpspacejs(dffts%nnr),tmpspacec(dfftp%nnr))

   numw_prod=0



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
     CALL diropn( iungprod, 'wiwjwfc', ngm*2, exst )
   else
    CALL diropn( iungprod, 'wiwjwfc', npw0*2, exst )
   endif
   iunprod = find_free_unit()
   open( unit= iunprod, file=trim(prefix)//'.wiwjprod', status='unknown',form='unformatted')



!open output file


   sum_vcvc = 0.d0


   do iw=1,nbndv

!read real wfcs iw


       tmpreal(:) =0.d0

       write(nfilei,'(5i1)') iw/10000,mod(iw,10000)/1000,mod(iw,1000)/100,mod(iw,100)/10,mod(iw,10)
       iunrealwan = find_free_unit()
       open(unit=iunrealwan,file='realwan'// nfilei,status='old', form='unformatted')
       read(iunrealwan) no1,no2,no3,nll
       read(iunrealwan) tmpreal(1:(2*nll+1)**3)
       close(iunrealwan)

       tmpspacei(:)=0.d0

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
               tmpspacei(nn)=tmpreal(iz*(2*nll+1)*(2*nll+1)+iy*(2*nll+1)+ix+1)

             enddo
          enddo
       enddo






      do jw=min(iw,nbndv),nbnd

!determines bottoms and tops of boxes in grid units

        match = .true.

        do i=1,3
           ndist=ndistance(w_centers(i,iw),w_centers(i,jw),rspacel(i))
           if(ndist-w_radii(i)-w_radii(j) > 0 ) then
               match = .false.
               exit
           endif
        enddo

        if( match) then

         write(stdout,*) 'Matching :', iw,jw
!determine center of product and radii
         do i=1,3
           call noverlap(w_centers(i,iw),w_centers(i,jw),w_radii(iw),w_radii(jw),&
            &wiwj%center(i),wiwj%radius(i),rspacel(i))
            if(is_even(rspacel(i)).and.wiwj%radius(i)==(rspacel(i)/2)) &
                & wiwj%radius(i)=wiwj%radius(i)-1
         enddo

         write(stdout,*) 'Center:', wiwj%center(1:3),wiwj%radius(1:3)

!read wannier function jw


         tmpreal(:) =0.d0

         write(nfilej,'(5i1)') jw/10000,mod(jw,10000)/1000,mod(jw,1000)/100,mod(jw,100)/10,mod(jw,10)
         iunrealwan = find_free_unit()
         open(unit=iunrealwan,file='realwan'// nfilej,status='old', form='unformatted')
         read(iunrealwan) no1,no2,no3,nll
         read(iunrealwan) tmpreal(1:(2*nll+1)**3)
         close(iunrealwan)

         tmpspacej(:)=0.d0

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
                 tmpspacej(nn)=tmpreal(iz*(2*nll+1)*(2*nll+1)+iy*(2*nll+1)+ix+1)

               enddo
            enddo
         enddo

! calculates product in R space

!determines origin
        nop(:)= wiwj%center(:)-wiwj%radius(:)
        do i=1,3
           if(nop(i)<1) nop(i)=rspacel(i)+nop(i)
           if(nop(i)>rspacel(i)) nop(i)=nop(i)-rspacel(i)
        enddo

        tmpreal(:)=0.d0
!adds US term if required
        tmpreal2(:)=0.d0
        if(okvan) call adduspos_gamma_r(iw,jw,tmpreal2,1,becp_gw(:,iw),becp_gw(:,jw))


        rsca=0.d0

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

                 nn=(n3-1)*dfftp%nr1x*dfftp%nr2x+(n2-1)*dfftp%nr1x+n1
                 rsca=rsca+(tmpspacei(nn)*tmpspacej(nn)+tmpreal2(nn))**2.d0
                 tmpreal(iz*(2*wiwj%radius(2)+1)*(2*wiwj%radius(1)+1)+iy*(2*wiwj%radius(1)+1)+ix+1)=&
                 &tmpspacei(nn)*tmpspacej(nn)+tmpreal2(nn)

               enddo
            enddo
         enddo

        rsca=rsca/real(rspacel(1)*rspacel(2)*rspacel(3))
!check if there is overlap
        if(rsca >= cutoff_wpr) then
           write(stdout,*) 'OVERLAP :', iw,jw,rsca
           if(jw > nbndv) sum_vcvc = sum_vcvc + rsca

           numw_prod=numw_prod+1

           wiwj%i=iw
           wiwj%j=jw

           write(iunprod) wiwj
!write on R space
           iunrprod = find_free_unit()
           if(.not.lsmallgrid) then
             open(unit=iunrprod,file= 'w_prod_r'//nfilei//nfilej,status='unknown',form='unformatted')
             write(iunrprod) tmpreal(1:(2*wiwj%radius(1)+1)*(2*wiwj%radius(2)+1)*(2*wiwj%radius(3)+1))
             close(iunrprod)
           endif

!put on grid for Fourier transform
          tmpspacej(:)=(0.d0,0.d0)

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

                   nn=(n3-1)*dfftp%nr1x*dfftp%nr2x+(n2-1)*dfftp%nr1x+n1
                   tmpspacej(nn)=tmpreal(iz*(2*wiwj%radius(2)+1)*(2*wiwj%radius(1)+1)+iy*(2*wiwj%radius(1)+1)+ix+1)
                 enddo
              enddo
           enddo

!bring back to G space on charge grid (4times denser in NC case)
!maybe wfcs grid could be enough?

          if(.not.lsmallgrid) then
            tmpspacec(:)=dcmplx(tmpspacej(:),0.d0)
            !ATTENZIONE
            rsca=0.d0
            do nn=1,dfftp%nnr
              rsca=rsca+conjg( tmpspacec(nn))*tmpspacec(nn)
            enddo
            write(*,*) 'RSCA', rsca/(dfftp%nr1*dfftp%nr2*dfftp%nr3)!ATTENZIONE
            !CALL cft3s( tmpspacec,nr1, nr2, nr3, nr1x, nr2x, nr3x, -1 )
             CALL fwfft ('Dense', tmpspacec, dfftp)
             rsca=0.d0
            do nn=1,ngm
              rsca=rsca+2.d0*conjg( tmpspacec(nl(nn)))*tmpspacec(nl(nn))
            enddo
             rsca=rsca-conjg( tmpspacec(nl(1)))*tmpspacec(nl(1))
            write(*,*) 'RSCA', rsca!ATTENZIONE
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
            warray(1:ngm)=tmpspacec(nl(1:ngm))
            CALL davcio(warray,ngm*2,iungprod,numw_prod,1)
          else
            warray(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
            CALL davcio(warray,npw0*2,iungprod,numw_prod,1)
          endif


          endif!overlap




        endif!match



     enddo!jw
  enddo!iw

  write(stdout,*) 'SUM VCVC', sum_vcvc

  close(iungprod)
  close(iunprod)
  deallocate(tmpspacei)
  deallocate(tmpspacej)
  deallocate(tmpreal)
  deallocate(tmpreal2)
  deallocate(tmpspacejs)
  deallocate(tmpspacec)
  deallocate(warray)

! #endif __GWW
END SUBROUTINE product_wannier


FUNCTION ndistance(n1,n2,l)
!find minimum distance between n1 and n2 with periodicity l


  implicit none

  INTEGER :: ndistance
  INTEGER :: n1,n2,l
  INTEGER :: i, nmin, imin

  nmin=2*l
  do i=-1,1
     if(abs(n1-(n2+i*l)) < nmin) then
       nmin = abs(n1-(n2+i*l))
       imin= i
     endif
  enddo
  ndistance = abs(n1-(n2+imin*l))

END FUNCTION ndistance


SUBROUTINE noverlap(n1,n2,r1,r2,c1,r,l)

!find center and radius of overlap with periodicity l

  implicit none

  INTEGER n1,n2,r1,r2,c1,r,l
  INTEGER :: i, nmin, imin
  INTEGER n22

  nmin=2*l
  do i=-1,1
     if(abs(n1-(n2+i*l)) < nmin) then
       nmin = abs(n1-(n2+i*l))
       imin= i
     endif
  enddo
!n22 is not periodic any more
  n22=n2+imin*l
  if(n22 >= n1) then
    c1=nint(real(n22-n1)/real(2)) + n1

    if(c1 < 1) c1=l+c1
    if(c1 > l) c1=c1-l

    r=((n22 + r2)-(n1-r1))/2
    if(r>=l/2) r=l/2
   else
    c1=nint(real((n1-n22))/real(2)) + n22

    if(c1 < 1) c1 = l + c1
    if(c1 > l) c1 = c1 - l

    r=nint(real(((n1+r1)-(n22-r2)))/real(2))
    if(r>=l/2) r=l/2

   endif

END SUBROUTINE noverlap

