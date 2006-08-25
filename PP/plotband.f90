
!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program read_bands
#include "f_defs.h"

  USE io_global,  ONLY : stdout

  implicit none
  real, allocatable :: e(:,:), k(:,:), e_in(:), kx(:)
  real :: k1(3), k2(3), xk1, xk2, ps
  real, allocatable :: e_rap(:,:), k_rap(:,:)
  integer, allocatable :: nbnd_rapk(:), rap(:,:)
  integer, allocatable :: npoints(:)
  integer :: nks = 0, nbnd = 0, ios, nlines, n,i,j,ni,nf,nl, ierr, ilen
  integer :: nks_rap = 0, nbnd_rap = 0
  logical, allocatable :: high_symmetry(:), is_in_range(:)
  character(len=256) :: filename, filename1
  namelist /plot/ nks, nbnd
  namelist /plot_rap/ nks_rap, nbnd_rap
  integer :: n_interp, init
  real, allocatable :: k_interp(:), e_interp(:), coef_interp(:,:)

  real :: emin = 1.e10, emax =-1.e10, etic, eref, deltaE, Ef

  integer, parameter :: max_lines=99 
  real :: point(max_lines+1), mine
  integer :: nrap(max_lines)
  integer :: ilines, irap, ibnd, ipoint

  real, parameter :: cm=28.453, xdim=15.0*cm, ydim=10.0*cm, &
                     x0=2.0*cm, y0=2.0*cm, eps=1.e-4

  logical :: exist_rap
  logical, allocatable :: todo(:)


  call get_file ( filename )
  open(unit=1,file=filename,form='formatted')
  read (1, plot, iostat=ios)
  !
  if (nks <= 0 .or. nbnd <= 0 .or. ios /= 0) then
     stop 'Error reading file header'
  else
     print '("Reading ",i4," bands at ",i4," k-points")', nbnd, nks
  end if

  filename1=TRIM(filename)//".rap"
  exist_rap=.true.
  open(unit=21,file=filename1,form='formatted',err=100,iostat=ios)

100 if (ios .ne. 0) then
     exist_rap=.false.
  endif
  if (exist_rap) read (21, plot_rap, iostat=ios)
  if (nks_rap.ne.nks.or.nbnd_rap.ne.nbnd.or.ios.ne.0) then
     write(6,'("file with representation but not compatible with bands")')
     exist_rap=.false.
  endif
  !
  allocate (e(nbnd,nks))
  allocate (k(3,nks), e_in(nks), kx(nks), npoints(nks), high_symmetry(nks))

  if (exist_rap) then
     allocate(nbnd_rapk(nks))
     allocate(e_rap(nbnd,nks))
     allocate(rap(nbnd,nks))
     allocate(k_rap(3,nks))
     allocate(todo(nbnd))
  end if

  do n=1,nks
     read(1,*,end=20,err=20) ( k(i,n), i=1,3 )
     read(1,*,end=20,err=20) (e(i,n),i=1,nbnd)
     if (n==1) then
        kx(n) = sqrt (k(1,1)**2 + k(2,1)**2 + k(3,1)**2)
     else
        kx(n) = kx(n-1) + sqrt ( (k(1,n)-k(1,n-1))**2 + &
                                 (k(2,n)-k(2,n-1))**2 + &
                                 (k(3,n)-k(3,n-1))**2 )
     end if

     if (exist_rap) then
        read(21,*,end=20,err=20) (k_rap(i,n),i=1,3)
        read(21,*,end=20,err=20) (rap(i,n),i=1,nbnd)
        if (abs(k(1,n)-k_rap(1,n))+abs(k(2,n)-k_rap(2,n))+  &
            abs(k(3,n)-k_rap(3,n))  > eps ) then
            write(stdout,'("Incompatible k points in rap file")')
            deallocate(nbnd_rapk)
            deallocate(e_rap)
            deallocate(rap)
            deallocate(k_rap)
            deallocate(todo)
            close(unit=21)
            exist_rap=.false.
        end if
     end if
  end do
  close(unit=1)
  if (exist_rap) close(unit=21)

  do n=1,nks
     do i=1,nbnd
        emin = min(emin, e(i,n))
        emax = max(emax, e(i,n))
     end do
  end do
  print '("Range:",2f8.4,"eV  Emin, Emax > ",$)', emin, emax
  read(5,*) emin, emax

  allocate (is_in_range(nbnd))
  is_in_range(:) = .false.
  do i=1,nbnd
     is_in_range(i) = any (e(i,1:nks) >= emin .and. e(i,1:nks) <= emax)
  end do

  do n=1,nks
     if (n==1 .or. n==nks) then
        high_symmetry(n) = .true.
     else
        do i=1,3
           k1(i) = k(i,n)-k(i,n-1)
           k2(i) = k(i,n+1)-k(i,n)
        end do
        ps = ( k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3) ) / &
         sqrt( k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3) ) / &
         sqrt( k2(1)*k2(1) + k2(2)*k2(2) + k2(3)*k2(3) )
        high_symmetry(n) = abs(ps-1.0) .gt.1.0e-4
     end if

     if (high_symmetry(n)) then
        if (n==1) then
           nlines=0
           npoints(1)=1
        else if (n==nks) then
           npoints(nlines+1) = npoints(nlines+1)+1
           nlines=nlines+1
        else
           npoints(nlines+1) = npoints(nlines+1)+1
           nlines=nlines+1
           npoints(nlines+1) = 1
        end if
        point(nlines+1)=n
        WRITE( stdout,'("high-symmetry point: ",3f7.4)') (k(i,n),i=1,3)
     else
        npoints(nlines+1) = npoints(nlines+1)+1
     end if
  end do
  !
  print '("output file (xmgr) > ",$)'
  read(5,'(a)', end=25, err=25)  filename
  if (filename == ' ' ) then
     print '("skipping ...")'
     go to 25
  end if
  if (.not.exist_rap) then
     open (unit=2,file=filename,form='formatted',status='unknown',&
           iostat=ios)  
     ! draw bands
     do i=1,nbnd
        if (is_in_range(i)) then
           if ( mod(i,2) /= 0) then
              write (2,'(2f10.4)') (kx(n), e(i,n),n=1,nks)
           else
              write (2,'(2f10.4)') (kx(n), e(i,n),n=nks,1,-1)
           end if
        end if
     end do
     close (unit = 2)
  else
!
!   In this case we write a diffent file for each line and for each
!   representation. Each file contains the bands of that representation.
!   The file is called filename.#line.#rap
!
!
!   First determine for each line how many representations are there
!
     do ilines=1,nlines
        nrap(ilines)=0
        do ipoint=1,npoints(ilines)-2
           n=point(ilines) + ipoint
           do ibnd=1,nbnd
              nrap(ilines)=max(nrap(ilines),rap(ibnd,n))
           end do
        end do
        write(6,*) 'lines nrap',ilines, nrap(ilines)
     end do
!
!   Then, for each line and for each representation along that line
!
     do ilines=1,nlines
        if (nrap(ilines)==0) then
!
!   Along this line the symmetry decomposition has not been done. 
!   Plot all the bands as in the standard case
!
           if (ilines<10) then
              write(filename1,'(a,".",i1)') TRIM(filename), ilines
           else
              write(filename1,'(a,".",i1)') TRIM(filename), ilines
           endif
           open (unit=2,file=filename1,form='formatted',status='unknown',&
                iostat=ios)  
           ! draw bands
           do i=1,nbnd
              if (is_in_range(i)) then
                 if ( mod(i,2) /= 0) then
                    write (2,'(2f10.4)') (kx(n), e(i,n),n=point(ilines),&
                                                          point(ilines+1))
                 else
                    write (2,'(2f10.4)') (kx(n), e(i,n),n=point(ilines+1), &
                                                          point(ilines),-1 )
                 end if
              end if
           end do
           close (unit = 2)
        endif
        do irap=1, nrap(ilines)
!
!     open a file
!
           if (ilines>99.or.irap>12) then
              write(6,'("too many lines or rap")')
              stop
           endif
           if (ilines < 10) then
              if (irap < 10 ) then
                 write(filename1,'(a,".",i1,".",i1)') TRIM(filename),ilines,irap
              else
                 write(filename1,'(a,".",i1,".",i2)') TRIM(filename),ilines,irap
              endif
           else
              if (irap < 10 ) then
                 write(filename1,'(a,".",i2,".",i1)') TRIM(filename),ilines,irap
              else
                 write(filename1,'(a,".",i2,".",i2)') TRIM(filename),ilines,irap
              endif
           endif
           open (unit=2,file=filename1,form='formatted',status='unknown',&
                 iostat=ios)
!  For each k point along this line selects only the bands which belong
!  to the irap representation
           nbnd_rapk=100000
           do n=point(ilines)+1, point(ilines+1)-1
              nbnd_rapk(n)=0
              do i=1,nbnd
                 if (rap(i,n)==irap) then
                    nbnd_rapk(n) = nbnd_rapk(n) + 1
                    e_rap(nbnd_rapk(n),n)=e(i,n)
                 endif
              enddo
           enddo
!
!   on the two high symmetry points the representation is different. So for each
!   band choose the closest eigenvalue available.
!
           todo=.true.
           do i=1,nbnd_rapk(point(ilines)+1)
              mine=1.d8
              do j=1,nbnd
                 if (abs(e_rap(i,point(ilines)+1)-e(j,point(ilines)))<mine &
                                                           .and. todo(j)) then
                    e_rap(i,point(ilines))=e(j,point(ilines))
                    mine=abs( e_rap(i,point(ilines)+1)-e(j,point(ilines)))
                    todo(j)=.false.
                 end if
              end do
           end do
           todo=.true.
           do i=1,nbnd_rapk(point(ilines+1)-1)
              mine=1.d8
              do j=1,nbnd
                 if (abs(e_rap(i,point(ilines+1)-1)- &
                          e(j,point(ilines+1)))<mine .and. todo(j)) then
                    e_rap(i,point(ilines+1))=e(j,point(ilines+1))
                    mine=abs(e_rap(i,point(ilines+1)-1)-e(j,point(ilines+1)) )
                    todo(j)=.false.
                 end if
              end do
           end do
           do i=1,MINVAL(nbnd_rapk)
              if (is_in_range(i)) then
                 if ( mod(i,2) /= 0) then
                    write (2,'(2f10.4)') ((kx(n), e_rap(i,n)), &
                                        n=point(ilines),point(ilines+1))
                 else
                    write (2,'(2f10.4)') ((kx(n), e_rap(i,n)), &
                                       n=point(ilines+1),point(ilines),-1)
                 end if
              end if
           end do
           if (MINVAL(nbnd_rapk)==0) THEN
              close (unit = 2,status='delete')
           else
              close (unit = 2)
           endif
        end do
     end do
  endif
  print '("bands in xmgr format written to file ",a)', filename
  !
25 continue
  if (exist_rap) then
     deallocate(nbnd_rapk)
     deallocate(e_rap)
     deallocate(rap)
     deallocate(k_rap)
     deallocate(todo)
  endif
  print '("output file (ps) > ",$)'
  read(5,'(a)',end=30,err=30)  filename
  if (filename == ' ' ) then
     print '("stopping ...")'
     go to 30
  end if
  open (unit=1,file=filename,form='formatted',status='unknown',&
       iostat=ios)  
  print '("Efermi > ",$)'
  read(5,*) Ef
  print '("deltaE, reference E (for tics) ",$)'
  read(5,*) deltaE, eref
  !
  write (1,'(a)') '%! PS-Adobe-1.0'
  write (1,*) '/localdict 100 dict def'
  write (1,*) 'localdict begin'
  write (1,*) '% delete next line for insertion in a LaTeX file'
  write (1,*) ' 0 0 moveto'
  write (1,*) 'gsave'
  write (1,*) '/nm  {newpath moveto} def'
  write (1,*) '/riga {newpath moveto lineto stroke} def'
  write (1,*) '/banda {3 1 roll moveto {lineto} repeat stroke} def'
  write (1,*) '/dot {newpath  1 0 360 arc fill} def'
  write (1,*) '/Times-Roman findfont 12 scalefont setfont'
  write (1,*) 'currentpoint translate'
  write (1,*) '% Landscape: uncomment next line'
  write (1,*) ' 90 rotate 0 21 neg 28.451 mul translate 1.5 1.5 scale'
  write (1,*) '% Landscape:   comment next line'
  write (1,*) '% 1.2 1.2 scale'
  write (1,'(2(f8.3,1x)," translate")') x0, y0
  write (1,*) '0 setgray 0.5 setlinewidth'
  ! draw tics on axis
  ni=nint((eref-emin)/deltaE)+1
  nf=nint((emax-eref)/deltaE)+1
  do i=-ni,nf
     etic=eref+i*deltaE
     if (etic >= emin .and. etic <= emax) then
        write (1,'(2(f8.3,1x)," moveto -5 0 rlineto stroke")') &
             0.0,(etic-emin)*ydim/(emax-emin)
        write (1,'(2(f8.3,1x)," moveto (",f4.1,") show")')   &
             -30.,(etic-emin)*ydim/(emax-emin), etic-eref
     end if
  end do
  ! draw the Fermi Energy
  if (Ef > emin .and. Ef < emax) then
     write (1,'("[2 4] 0 setdash newpath ",2(f8.3,1x), " moveto ")') &
          0.0, (Ef-emin)/(emax-emin)*ydim
     write (1,'(2(f8.3,1x)," lineto stroke [] 0 setdash")') &
          xdim, (Ef-emin)/(emax-emin)*ydim
  end if
  ! draw axis and set clipping region
  write (1,*) '1 setlinewidth'
  write (1,'(8(f8.3,1x))') 0.0,0.0,0.0,ydim,xdim,ydim,xdim,0.0 
  write (1,*) 'newpath moveto lineto lineto lineto closepath clip stroke'
  write (1,*) '0.5 setlinewidth'
  ! draw high-symmetry lines
  do n=1,nks
     if (high_symmetry(n)) then
        write (1,'(4(f8.3,1x)," riga")') &
             kx(n)*xdim/kx(nks), 0.0, kx(n)*xdim/kx(nks), ydim
     end if
     do i=1,nbnd
        if (is_in_range(i)) write (1,'(2(f8.3,1x)," dot")' ) &
             kx(n)*xdim/kx(nks), (e(i,n)-emin)*ydim/(emax-emin)
     end do
  end do
  ! draw bands
  allocate (k_interp(4*nks), e_interp(4*nks), coef_interp(nks,4))
  do i=1,nbnd
     if (is_in_range(i)) then
        ! No interpolation:
        !         write (1,'(9(f8.3,1x))') ( kx(n)*xdim/kx(nks), &
        !             (e(i,n)-emin)*ydim/(emax-emin),n=nks,1,-1)
        !         write (1,'(i4," banda")' ) nks-1
        ! Spline interpolation with twice as many points:
        !
        ni=1
        nf=1
        do nl=1,nlines
           ni=nf
           nf=nf + npoints(nl)-1
           n_interp= 2*(nf-ni)+1
           do n=1,n_interp
              k_interp(n)=kx(ni)+(n-1)*(kx(nf)-kx(ni))/(n_interp-1)
           end do
           do n=ni,nf
              e_in(n-ni+1)=e(i,n)
           end do
           call spline_interpol ( kx(ni), e_in, nf-ni+1, &
                k_interp, e_interp, n_interp )
           write (1,'(9(f8.3,1x))') ( k_interp(n)*xdim/kx(nks), &
                (e_interp(n)-emin)*ydim/(emax-emin),n=n_interp,1,-1)
           write (1,'(i4," banda")' ) n_interp-1
        end do
     end if
  end do

  write (1,*) 'grestore'
  write (1,*) '% delete next lines for insertion in a tex file'
  write (1,'(a)') '%%Page'
  write (1,*) 'showpage'
  close (unit=1)
  print '("bands in PostScript format written to file ",a)', filename
30 continue

  stop
20 print '("Error reading k-point # ",i4)', n
  stop

end program read_bands

subroutine spline_interpol (xin, yin, nin, xout, yout, nout)

  ! xin and xout should be in increasing order, with
  ! xout(1) <= xin(1), xout(nout) <= xin(nin)

  implicit none
  integer, intent(in) :: nin, nout
  real, intent(in)  :: xin(nin), yin(nin), xout(nout)
  real, intent(out) :: yout(nout)
  ! work space (automatically allocated)
  real :: d2y(nin)
  real :: dy1, dyn

  dy1 = (yin(2)-yin(1))/(xin(2)-xin(1))
  dyn = 0.0

  call spline( xin, yin, nin, dy1, dyn, d2y)
  call splint( nin, xin, yin, d2y, nout, xout, yout)

  return                                              
end subroutine spline_interpol

subroutine spline(x, y, n, yp1, ypn, d2y)

  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x(n), y(n), yp1, ypn
  real, intent(out):: d2y(n)
  ! work space (automatically allocated)
  real :: work(n)
  integer :: i, k
  real :: sig, p, qn, un

  d2y(1)=-0.5
  work(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

  do i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*d2y(i-1)+2.0
     d2y(i)=(sig-1.0)/p
     work(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*work(i-1))/p
  end do
  qn=0.5
  un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  
  d2y(n)=(un-qn*work(n-1))/(qn*d2y(n-1)+1.0)
  do k=n-1,1,-1
     d2y(k)=d2y(k)*d2y(k+1)+work(k)
  end do

  return
end subroutine spline


subroutine splint (nspline, xspline, yspline, d2y, nfit, xfit, yfit)

  implicit none
  ! input
  integer, intent(in) :: nspline, nfit
  real, intent(in) :: xspline(nspline), yspline(nspline), xfit(nfit), &
       d2y(nspline)
  real, intent(out) :: yfit(nfit)
  integer :: klo, khi, i
  real :: a, b, h 

  klo=1
  do i=1,nfit
     do khi=klo+1, nspline
        if(xspline(khi) >= xfit(i)) then
           if(xspline(khi-1) <= xfit(i)) then
              klo = khi-1
           else
              if (klo == 1 .and. khi-1 == 1) then
                 ! the case xfit(i) < xspline(1) should not happen
                 ! but since it may be due to a numerical artifact
                 ! we just continue
                 print *, '  SPLINT WARNING: xfit(i) < xspline(1)', &
                      xfit(i), xspline(1)  
              else
                 stop '  SPLINT ERROR: xfit not properly ordered'
              end if
           end if
           h= xspline(khi) - xspline(klo) 
           a= (xspline(khi)-xfit(i))/h
           b= (xfit(i)-xspline(klo))/h
           
           yfit(i) = a*yspline(klo) + b*yspline(khi) &
                + ( (a**3-a)*d2y(klo) + (b**3-b)*d2y(khi)  )*h*h/6.0
           go to 10
        end if
     end do

     ! the case xfit(i) > xspline(nspline) should also not happen
     ! but again it may be due to a numerical artifact
     ! A properly chosen extrapolation formula should be used here
     ! (and in the case  xfit(i) < xspline(1) above as well) but
     ! I am too lazy to write one - PG

     print *, '  SPLINT WARNING: xfit(i) > xspline(nspline)', &
                      xfit(i), xspline(nspline)  
     khi = klo+1
     h= xspline(khi) - xspline(klo) 
     a= (xspline(khi)-xfit(i))/h
     b= (xfit(i)-xspline(klo))/h
           
     yfit(i) = a*yspline(klo) + b*yspline(khi) &
          + ( (a**3-a)*d2y(klo) + (b**3-b)*d2y(khi)  )*h*h/6.0
     !
10   continue
  end do
  
  return
end subroutine splint

