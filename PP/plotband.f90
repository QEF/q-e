
program read_bands

  implicit none
  integer, parameter:: maxk=100
  real :: k(3,maxk), e_in(maxk), kx(maxk)
  real, allocatable :: e(:,:)
  real :: k1(3), k2(3), xk1, xk2, ps
  integer :: npoints(maxk)
  integer :: npk, nbands, nlines, n,i,ni,nf,nl, iargc
  logical :: high_symmetry(maxk)
  logical, allocatable :: is_in_range(:)
  character(len=80) :: filename, prgname

  integer, parameter ::  max_interp=4*maxk
  integer :: n_interp, init
  real :: k_interp(max_interp), e_interp(max_interp), coef_interp(maxk,4)

  real :: emin = 1.e10, emax =-1.e10,  etic, eref, deltaE, Ef
  real, parameter :: cm=28.453, xdim=15.0*cm, ydim=10.0*cm, &
                     x0=2.0*cm, y0=2.0*cm


  i=iargc()
  call getarg(0,prgname)
  if (i==0) then
10   print '("input file > ",$)'
     read(5,'(a)',end=10,err=10)  filename
  else if(i==1) then
     call getarg(1,filename)
  else
     print '("usage:",a20," [input file] ")',prgname
  end if
  open(unit=1,file=filename,form='formatted')

  print '("number of bands > ",$)'
  read(5,*) nbands
  allocate (e(nbands,maxk))
  do n=1,maxk
     read(1,*,end=20,err=30) ( k(i,n), i=1,3 )
     read(1,*,end=30,err=30) (e(i,n),i=1,nbands)
     if (n==1) then
        kx(n) = 0.0
     else
        kx(n) = kx(n-1) + sqrt ( (k(1,n)-k(1,n-1))**2 + &
                                 (k(2,n)-k(2,n-1))**2 + &
                                 (k(3,n)-k(3,n-1))**2 )
     end if
  end do
  close(unit=1)
  print '("Warning: max # of k-point (",i3,") read")', maxk
20 npk=n-1
  print '(i3," k-point read")', npk

  do n=1,npk
     do i=1,nbands
        emin = min(emin, e(i,n))
        emax = max(emax, e(i,n))
     end do
  end do
  print '("Range:",2f8.4,"eV  Emin, Emax > ",$)', emin, emax
  read(5,*) emin, emax
  print '("Efermi > ",$)'
  read(5,*) Ef
  print '("deltaE, reference E (for tics) ",$)'
  read(5,*) deltaE, eref

  allocate (is_in_range(nbands))
  is_in_range(:) = .false.
  do i=1,nbands
     is_in_range(i) = any (e(i,1:npk) >= emin .and. e(i,1:npk) <= emax)
  end do

  do n=1,npk
     if (n==1 .or. n==npk) then
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
        else if (n==npk) then
           npoints(nlines+1) = npoints(nlines+1)+1
           nlines=nlines+1
        else
           npoints(nlines+1) = npoints(nlines+1)+1
           nlines=nlines+1
           npoints(nlines+1) = 1
        end if
        write(6,'("high-symmetry point: ",3f7.4)') (k(i,n),i=1,3)
     else
        npoints(nlines+1) = npoints(nlines+1)+1
     end if
  end do

25 continue
  print '("output file > ",$)'
  read(5,'(a)',end=25,err=25)  filename

  open (unit=1,file=filename,form='formatted',status='unknown')

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

  write (1,'(2(f8.3,x)," translate")') x0, y0
  write (1,*) '0 setgray 0.5 setlinewidth'
  ! draw tics on axis
  ni=nint((eref-emin)/deltaE)+1
  nf=nint((emax-eref)/deltaE)+1
  do i=-ni,nf
     etic=eref+i*deltaE
     if (etic >= emin .and. etic <= emax) then
        write (1,'(2(f8.3,x)," moveto -5 0 rlineto stroke")') &
             0.0,(etic-emin)*ydim/(emax-emin)
        write (1,'(2(f8.3,x)," moveto (",f4.1,") show")')   &
             -30.,(etic-emin)*ydim/(emax-emin), etic-eref
     end if
  end do
  ! draw the Fermi Energy
  if (Ef > emin .and. Ef < emax) then
     write (1,'("[2 4] 0 setdash newpath ",2(f8.3,x), " moveto ")') &
          0.0, (Ef-emin)/(emax-emin)*ydim
     write (1,'(2(f8.3,x)," lineto stroke [] 0 setdash")') &
          xdim, (Ef-emin)/(emax-emin)*ydim
  end if
  ! draw axis and set clippping region
  write (1,*) '1 setlinewidth'
  write (1,'(8(f8.3,x))') 0.0,0.0,0.0,ydim,xdim,ydim,xdim,0.0 
  write (1,*) 'newpath moveto lineto lineto lineto closepath clip stroke'
  write (1,*) '0.5 setlinewidth'
  ! draw high-symmetry lines
  do n=1,npk
     if (high_symmetry(n)) then
        write (1,'(4(f8.3,x)," riga")') &
             kx(n)*xdim/kx(npk), 0.0, kx(n)*xdim/kx(npk), ydim
     end if
     do i=1,nbands
        if (is_in_range(i)) write (1,'(2(f8.3,x)," dot")' ) &
             kx(n)*xdim/kx(npk), (e(i,n)-emin)*ydim/(emax-emin)
     end do
  end do
  ! draw bands
  do i=1,nbands
     if (is_in_range(i)) then
        ! No interpolation:
        !         write (1,'(9(f8.3,x))') ( kx(n)*xdim/kx(npk), &
        !             (e(i,n)-emin)*ydim/(emax-emin),n=npk,1,-1)
        !         write (1,'(i4," banda")' ) npk-1
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
           write (1,'(9(f8.3,x))') ( k_interp(n)*xdim/kx(npk), &
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

  stop
30 print '("Error reading k-point # ",i4)', n
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

