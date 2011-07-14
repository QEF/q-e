!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine init_gper(ik)
!
! - To compute the number of 2D G vectors with |G+k|^2<ecut2d (ngper),
! - the number of shells (ngpsh) for them,
!
  USE io_global,        ONLY :  stdout
  USE noncollin_module, ONLY : npol
  USE cell_base,        ONLY : bg, tpiba, tpiba2
  USE fft_base,         ONLY : dfftp
  USE cond
  implicit none
  integer :: ipol, igper, icount, ik, k, ig, i, il, j, jl, iw
  integer, allocatable :: nshell(:,:)
  real(DP) :: norm2
  real(DP), parameter :: eps=1.d-8
  real(DP), allocatable :: gnorm2(:)

  allocate( gnorm2( nrx * nry ) )
  allocate( nshell( nrx, nry ) )
!
! Compute ngper, ngpsh, and (i,j) --> shell
!
  do i=1, nrx
    do j=1, nry
      nshell(i,j)=0
    enddo
  enddo

  ngper=1
  ngpsh=1
  gnorm2(1)=( (xyk(1,ik)*bg(1,1)+xyk(2,ik)*bg(1,2))**2+ &
              (xyk(1,ik)*bg(2,1)+xyk(2,ik)*bg(2,2))**2 )*tpiba2
  nshell(1,1)=1
  do i=1, nrx
    il=i-1
    if (il.gt.nrx/2) il=il-nrx
    do j=1, nry
      jl=j-1
      if (jl.gt.nry/2) jl=jl-nry
      norm2=(((il+xyk(1,ik))*bg(1,1)+(jl+xyk(2,ik))*bg(1,2))**2+ &
            ((il+xyk(1,ik))*bg(2,1)+(jl+xyk(2,ik))*bg(2,2))**2)* &
                 tpiba2
      if (norm2.le.ecut2d.and.(i*j).ne.1) then
         ngper=ngper+1
         icount=0
         do iw=1, ngpsh
           if (abs(norm2-gnorm2(iw)).gt.eps) then
             icount=icount+1
           else
             nshell(i,j)=iw
           endif
         enddo
         if (icount.eq.ngpsh) then
           ngpsh=ngpsh+1
           gnorm2(ngpsh)=norm2
           nshell(i,j)=ngpsh
         endif
      endif
    enddo
  enddo

!
! Global variables
!
  allocate( gper( 2, ngper ) )
  if (lorb) allocate( nl_2ds( npol*ngper ) )
  if (lorb) allocate( nl_2d( npol*ngper ) )
  allocate( ninsh( ngpsh ) )
  allocate( gnsh( ngpsh ) )
!
! construct the g perpendicular vectors
!

! ninsh() is used as a pointer on the (first-1) element
! in each shell
!
  do i=1, ngpsh
    ninsh(i)=0
  enddo
  do i=1, nrx
    do j=1, nry
      if (nshell(i,j).ne.0) then
        do k=nshell(i,j)+1, ngpsh
           ninsh(k)=ninsh(k)+1
        enddo
      endif
    enddo
  enddo
!
! To form g
!
  do i=1, nrx
    il=i-1
    if (il.gt.nrx/2) il=il-nrx
    do j=1, nry
      jl=j-1
      if (jl.gt.nry/2) jl=jl-nry
      if (nshell(i,j).ne.0) then
         ninsh(nshell(i,j))=ninsh(nshell(i,j))+1
         igper=ninsh(nshell(i,j))
         do ipol=1,2
           gper(ipol,igper)=(il+xyk(1,ik))*bg(ipol,1)+   &
                             (jl+xyk(2,ik))*bg(ipol,2)
         enddo
         if (lorb) THEN
           nl_2ds(igper) = i+(j-1)*nrx
           ig = i
           iw = j
           if (il.lt.0) ig = il+1+dfftp%nr1
           if (jl.lt.0) iw = jl+1+dfftp%nr2
           nl_2d(igper) = ig+(iw-1)*dfftp%nr1
         END IF
      endif
    enddo
  enddo
!
! To set up |g| and number of g in each shell
!
  do i=ngpsh, 2, -1
    igper=ninsh(i)
    gnsh(i)=sqrt(gper(1,igper)**2+gper(2,igper)**2)*tpiba
    ninsh(i)=ninsh(i)-ninsh(i-1)
  enddo
  igper=ninsh(1)
  gnsh(1)=sqrt(gper(1,igper)**2+gper(2,igper)**2)*tpiba
  ninsh(1)=igper

  WRITE( stdout,*) 'ngper, shell number = ', ngper, ngpsh

  deallocate(gnorm2)
  deallocate(nshell)

  return
end subroutine init_gper
