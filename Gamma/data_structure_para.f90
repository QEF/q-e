!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine data_structure_para
  !-----------------------------------------------------------------------
  !
  ! distribute columns to processes for parallel fft
  ! columns are sets of g-vectors along z: g(k) = i1*b1+i2*b2+i3*b3 ,
  ! with g^2<gcut and (i1,i2) running over the (xy) plane.
  ! Columns are "active" for a given (i1,i2) if they contain a nonzero
  ! number of wavevectors
  !
#ifdef __PARA
  use para
  use pwcom
  use mp, only: mp_sum
  use mp_global, only: intra_pool_comm
  !
  implicit none
  !
  integer, allocatable :: &
       ngc(:),        &! number of g-vectors per column (dense grid)
       ngcs(:),       &! number of g-vectors per column (smooth grid
       ngcw(:),       &! number of wavefct plane waves per colum
       in1(:),        &! index i for column (i1,i2)
       in2(:),        &! index j for column (i1,i2)
       index(:)        ! used to order column
  integer ic,         &! fft index for this column (dense grid)
       ics,           &! as above for the smooth grid
       icm,           &! fft index for column (-i1,-i2) (dense grid)
       icms,          &! as above for the smooth grid
       ncp_(maxproc), &! number of column per processor (work space)
       ngp(maxproc),  &! number of g-vectors per proc (dense grid)
       ngps(maxproc), &! number of g-vectors per proc (smooth grid)
       ngpw(maxproc)   ! number of wavefct plane waves per proc
  !
  integer np, nps1,           &! counters on planes
       nq, nqs,               &! counters on planes
       max1,min1,max2,min2,   &! aux. variables
       m1, m2, n1, n2, i, mc, &! generic counter
       idum, nct_,            &! check variables
       j,jj,                  &! counters on processors
       n1m1,n2m1,n3m1,        &! nr1-1 and so on
       i1, i2, i3,            &! counters on G space
       good_fft_dimension      ! a function with obvious meaning
  logical has_gzero
  real(kind=8), allocatable :: aux(:)  ! used to order columns
  real(kind=8)  amod, gkcut        ! square modulus of G vectors
  !
  ! set the dimensions of fft arrays
  !
  nrx1 =good_fft_dimension(nr1 )
  nrx2 = nr2
  nrx3 =good_fft_dimension(nr3 )
  !
  nrx1s=good_fft_dimension(nr1s)
  nrx2s=nr2s
  nrx3s=good_fft_dimension(nr3s)
  !
  !     compute number of columns for each processor
  !
  ncplane = nrx1*nrx2
  ncplanes= nrx1s*nrx2s
  !
  !    global variables allocated here
  !
  allocate  (ipc ( ncplane))    
  allocate  (ipcs( ncplanes))    
  !
  !    local variables to be deallocated at the end
  !
  allocate  (ngc ( ncplane))    
  allocate  (ngcs( ncplane))    
  allocate  (ngcw( ncplane))    
  allocate  (in1 ( ncplane))    
  allocate  (in2 ( ncplane))    
  allocate  (index(ncplane))    
  allocate  (aux(  ncplane))    
  !
  ! set the number of plane per process
  !
  if (nr3.lt.nproc) call errore('set_fft_para',                      &
       &                'some processors have no planes ',-1)
  if (nr3s.lt.nproc) call errore('set_fft_para',                     &
       &                'some processors have no smooth planes ',-1)
  !
  if (nproc.eq.1) then
     npp(1) = nr3
     npps(1)= nr3s
  else
     np = nr3/nproc
     nq = nr3 - np*nproc
     nps1 = nr3s/nproc
     nqs = nr3s - nps1*nproc
     do i = 1, nproc
        npp(i) = np
        if (i.le.nq) npp(i) = np + 1
        npps(i) = nps1
        if (i.le.nqs) npps(i) = nps1 + 1
     end do
  end if
  !
  !     Now compute for each point of the big plane how many column have
  !     non zero vectors on the smooth and dense grid
  !
  do mc = 1,ncplane
     ipc(mc) = 0
  end do
  do mc = 1,ncplanes
     ipcs(mc)= 0
  end do
  !
  nct = 0
  ncts= 0
  !
  ! NOTA BENE: the exact limits for a correctly sized FFT grid are:
  ! -nr/2,..,+nr/2  for nr even; -(nr-1)/2,..,+(nr-1)/2  for nr odd.
  ! If the following limits are increased, a slightly undersized fft
  ! grid, with some degree of G-vector refolding, can be used
  ! (at your own risk - a check is done in ggen).
  !
  n1m1=nr1/2
  n2m1=nr2/2
  n3m1=nr3/2
  !
  gkcut = ecutwfc / tpiba2
  do i1=-n1m1,n1m1
     do i2=-n2m1,n2m1
        !
        ! nct counts columns containing G-vectors for the dense grid
        !
        nct=nct+1
        if (nct.gt.ncplane) call errore('set_fft_para','too many columns',1)
        ngc (nct) = 0
        ngcs(nct) = 0
        ngcw(nct) = 0
        !
        do i3 = -n3m1,n3m1
           amod = (bg(1,1)*i1 + bg(1,2)*i2 + bg(1,3)*i3)**2 +  &
                  (bg(2,1)*i1 + bg(2,2)*i2 + bg(2,3)*i3)**2 +  &
                  (bg(3,1)*i1 + bg(3,2)*i2 + bg(3,3)*i3)**2
           if (amod.le.gcutm)  ngc (nct)= ngc (nct)+ 1
           if (amod.le.gcutms) ngcs(nct)= ngcs(nct)+ 1
           if (amod.le.gkcut)  ngcw(nct)= ngcw(nct)+ 1
        enddo
        if (ngc(nct).gt.0) then
           !
           ! this column contains G-vectors
           !
           in1(nct) = i1
           in2(nct) = i2
           if (ngcs(nct).gt.0) then
              !
              ! ncts counts columns contaning G-vectors for the smooth grid
              !
              ncts=ncts+1
              if (ncts.gt.ncplanes) &
                   call errore('set_fft_para','too many columns',2)
           end if
        else
           !
           ! this column has no G-vectors: reset the counter
           !
           nct=nct-1
        end if
     enddo
  end do
  !
  if(nct .eq.0) call errore('set_fft_para','number of column 0', 1)
  if(ncts.eq.0) call errore('set_fft_para','number smooth column 0', 1)
  !
  !   Sort the columns. First the column with the largest number of G
  !   vectors on the wavefunction sphere (active columns),
  !   then on the smooth sphere, then on the big sphere. Dirty trick:
  !
  do mc = 1,nct
     aux(mc)=-(ngcw(mc)*nrx3**2 + ngcs(mc)*nrx3 + ngc(mc))
  end do
  call hpsort(nct,aux,index)
  !
  ! assign columns to processes
  !
  has_gzero=.false.
  do j=1,nproc
     ncp (j) = 0
     ncps(j) = 0
     nkcp(j) = 0
     ngp (j) = 0
     ngps(j) = 0
     ngpw(j) = 0
  end do
  !
  do mc=1, nct
     i = index(mc)
     !
     ! index contains the desired ordering of columns (see above)
     !
     i1=in1(i)
     i2=in2(i)
     !
     if ( i1.lt.0.or.(i1.eq.0.and.i2.lt.0) ) go to 30
     !
     ! only half of the columns, plus column (0,0), are scanned:
     ! column (-i1,-i2) must be assigned to the same proc as column (i1,i2)
     !
     ! ic  :  position, in fft notation, in dense grid, of column ( i1, i2)
     ! icm :      "         "      "          "    "         "    (-i1,-i2)
     ! ics :      "         "      "        smooth "         "    ( i1, i2)
     ! icms:      "         "      "        smooth "         "    (-i1,-i2)
     !
     m1 = i1 + 1
     if (m1.lt.1) m1 = m1 + nr1
     m2 = i2 + 1
     if (m2.lt.1) m2 = m2 + nr2
     ic = m1 + (m2-1)*nrx1
     !
     n1 = -i1 + 1
     if (n1.lt.1) n1 = n1 + nr1
     n2 = -i2 + 1
     if (n2.lt.1) n2 = n2 + nr2
     icm = n1 + (n2-1)*nrx1
     !
     m1 = i1 + 1
     if (m1.lt.1) m1 = m1 + nr1s
     m2 = i2 + 1
     if (m2.lt.1) m2 = m2 + nr2s
     ics = m1 + (m2-1)*nrx1s
     !
     n1 =-i1 + 1
     if (n1.lt.1) n1 = n1 + nr1s
     n2 =-i2 + 1
     if (n2.lt.1) n2 = n2 + nr2s
     icms = n1 + (n2-1)*nrx1s
     !
     jj=1
     if (ngcw(i).gt.0) then
        !
        ! this is an active column: find which processor has currently
        ! the smallest number of plane waves
        !
        do j=1,nproc
           if (ngpw(j).lt.ngpw(jj)) jj = j
        end do
     else
        !
        ! this is an inactive column: find which processor has currently
        ! the smallest number of G-vectors
        !
        do j=1,nproc
           if (ngp(j).lt.ngp(jj)) jj = j
        end do
     end if
     !
     ! jj is the processor to which this column is assigned
     ! use -jj for inactive columns, jj for active columns
     !
     ipc(ic) = -jj
     ncp(jj) = ncp(jj) + 1
     ngp(jj) = ngp(jj) + ngc(i)
     if (ngcs(i).gt.0) then
        ncps(jj)=ncps(jj)+1
        ngps(jj)=ngps(jj)+ngcs(i)
        ipcs(ics)=-jj
     endif
     if (ngcw(i).gt.0) then
        ipcs(ics)=jj
        ipc(ic) = jj
        ngpw(jj)= ngpw(jj) + ngcw(i)
        nkcp(jj)= nkcp(jj) + 1
     endif
     !
     if (i1.eq.0.and.i2.eq.0.and.jj.eq.me) has_gzero = .true.
     !
     ! now assign the (-i1,-i2) column to the same processor
     !
     if (i1.eq.0.and.i2.eq.0) go to 30
     !
     ! do not count twice column (0,0) !
     !
     ipc(icm) = -jj
     ncp(jj) = ncp(jj) + 1
     ngp(jj) = ngp(jj) + ngc(i)
     if (ngcs(i).gt.0) then
        ncps(jj)=ncps(jj)+1
        ngps(jj)=ngps(jj)+ngcs(i)
        ipcs(icms)=-jj
     endif
     if (ngcw(i).gt.0) then
        ipcs(icms)=jj
        ipc(icm) = jj
        ngpw(jj)= ngpw(jj) + ngcw(i)
        nkcp(jj)= nkcp(jj) + 1
     endif
30   continue
  enddo
  !
  ! ipc  is the processor for this column in the dense grid
  ! ipcs is the same, for the smooth grid
  !
  ! write(6,'(" Proc  planes cols    G   planes cols    G    columns  G"/ &
  !         & "       (dense grid)     (smooth grid)   (wavefct grid)'")')
  !do i=1,nproc
  !   write(6,'(i3,2x,3(i5,2i7))') i, npp(i),ncp(i),ngp(i), &
  !        npps(i),ncps(i),ngps(i), nkcp(i), ngpw(i)
  !end do
  !
  ! ngm contains the number of G vectors on this processor (me)
  ! remember that ngp counts all G-vectors, while we want only half
  ! same for ngms (smooth grid)
  !
  if (has_gzero) then
     ngm = (ngp (me)-1)/2 + 1
     ngms= (ngps(me)-1)/2 + 1
  else
     ngm = ngp (me)/2
     ngms= ngps(me)/2
  end if
  !
  ! nxx and nxxs are copies of nrxx and nrxxs, the local fft data size,
  ! to be stored in "parallel" commons. Not a very elegant solution.
  !
  if (nproc.eq.1) then
     nrxx =nrx1*nrx2*nrx3
     nrxxs=nrx1s*nrx2s*nrx3s
  else
     nrxx = max(nrx3*ncp(me), nrx1*nrx2*npp(me))
     nrxxs= max(nrx3s*ncps(me), nrx1s*nrx2s*npps(me))
  end if
  nxx = nrxx
  nxxs= nrxxs
  !
  !   computing the starting column for each processor
  !
  do i=1,nproc
     if(ngpw(i).eq.0) call errore('set_fft_para', &
          &        'some processors have no pencils, not yet implemented',1)
     if (i.eq.1) then
        ncp0(i) = 0
        ncp0s(i)= 0
     else
        ncp0(i) = ncp0 (i-1) + ncp (i-1)
        ncp0s(i)= ncp0s(i-1) + ncps(i-1)
     endif
  enddo
  !
  !  Now compute the arrays ipc and icpl (dense grid):
  !     ipc contain the number of the column for that processor.
  !         zero if the column do not belong to the processor.
  !         Note that input ipc is used and overwritten.
  !     icpl contains the point in the plane for each column
  !
  !- active columns first........
  !
  allocate (icpl( nct))    
  allocate (icpls( ncts))    
  !
  do j=1,nproc
     ncp_(j) = 0
  end do
  do mc =1,ncplane
     if (ipc(mc).gt.0) then
        j = ipc(mc)
        ncp_(j) = ncp_(j) + 1
        icpl(ncp_(j) + ncp0(j)) = mc
        if (j.eq.me) then
           ipc(mc) = ncp_(j)
        else
           ipc(mc) = 0
        end if
     end if
  end do
  !
  !-..... ( intermediate check ) ....
  !
  do j=1,nproc
     if (ncp_(j).ne.nkcp(j))                                        &
          &        call errore('set_fft_para','ncp_(j).ne.nkcp(j)',j)
  end do
  !
  !- ........then the remaining columns
  !
  do mc =1,ncplane
     if (ipc(mc).lt.0) then
        j = -ipc(mc)
        ncp_(j) = ncp_(j) + 1
        icpl(ncp_(j) + ncp0(j)) = mc
        if (j.eq.me) then
           ipc(mc) = ncp_(j)
        else
           ipc(mc) = 0
        end if
     end if
  end do
  !
  !-... ( final check )
  !
  nct_ = 0
  do j=1,nproc
     if (ncp_(j).ne.ncp(j))                                         &
          &        call errore('set_fft_para','ncp_(j).ne.ncp(j)',j)
     nct_ = nct_ + ncp_(j)
  end do
  if (nct_.ne.nct)                                                  &
       &     call errore('set_fft_para','nct_.ne.nct',1)
  !
  !   now compute the arrays ipcs and icpls
  !   (as ipc and icpls, for the smooth grid)
  !
  !   active columns first...
  !
  do j=1,nproc
     ncp_(j) = 0
  end do
  do mc =1,ncplanes
     if (ipcs(mc).gt.0) then
        j = ipcs(mc)
        ncp_(j)=ncp_(j) + 1
        icpls(ncp_(j) + ncp0s(j)) = mc
        if (j.eq.me) then
           ipcs(mc) = ncp_(j)
        else
           ipcs(mc) = 0
        endif
     endif
  enddo
  !
  !-..... ( intermediate check ) ....
  !
  do j=1,nproc
     if (ncp_(j).ne.nkcp(j))                                        &
          &        call errore('set_fft_para','ncp_(j).ne.nkcp(j)',j)
  end do
  !
  !    and then all the others
  !
  do mc =1,ncplanes
     if (ipcs(mc).lt.0) then
        j = -ipcs(mc)
        ncp_(j) = ncp_(j) + 1
        icpls(ncp_(j) + ncp0s(j)) = mc
        if (j.eq.me) then
           ipcs(mc) = ncp_(j)
        else
           ipcs(mc) = 0
        end if
     end if
  end do
  !
  !-... ( final check )
  !
  nct_ = 0
  do j=1,nproc
     if (ncp_(j).ne.ncps(j))                                        &
          &        call errore('set_fft_para','ncp_(j).ne.ncps(j)',j)
     nct_ = nct_ + ncp_(j)
  end do
  if (nct_.ne.ncts)                                                 &
       &     call errore('set_fft_para','nct_.ne.ncts',1)
  !
  deallocate (aux)
  deallocate (index)
  deallocate (in2)
  deallocate (in1)
  deallocate (ngcw)
  deallocate (ngcs)
  deallocate (ngc )
  !
  ngm_l  = ngm
  ngms_l = ngms
  ngm_g  = ngm
  ngms_g = ngms
  call mp_sum( ngm_g , intra_pool_comm )
  call mp_sum( ngms_g, intra_pool_comm )
  !
#endif
  return
end subroutine data_structure_para
