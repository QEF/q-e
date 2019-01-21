program grafica


  implicit none

  integer n,idumm,i,j
  real, allocatable :: d(:,:),t(:,:),proj(:)
  real :: x1,x2,s,off,f
  integer :: m, ipro
  real, external :: espandi
  real :: x,rdumm

  open(unit=20, file='qpe_range.dat', status='old')
  read(20,*) n
  allocate (d(n,4))
  allocate(proj(n))
  do i=1,n
     read(20,*) idumm,d(i,1),d(1,2),d(i,3),d(i,4)
  enddo
  close(20)
  write(*,*) 'x1'
  read(*,*) x1
  write(*,*) 'x2'
  read(*,*) x2
  write(*,*) 'sigma'
  read(*,*) s
  write(*,*) 'offset'
  read(*,*) off
  write(*,*) 'number of points'
  read(*,*) m
  write(*,*) 'Use projections file 0=No 1=Yes'
  read(*,*) ipro
  if(ipro==0) then
     proj=1.
  else
     open(unit=19, file='proj.dat', status='unknown')
     do i=1,n
        read(19,*) rdumm,proj(i)
     enddo
     close(19)
  endif
  write(*,*) 'Multiplicative factor'
  read(*,*) f

  open(unit=21, file='grafico_lda.dat', status='unknown')
  open(unit=22, file='grafico_gwp.dat', status='unknown')
  open(unit=23, file='grafico_gw.dat', status='unknown')
  open(unit=24, file='grafico_hf.dat', status='unknown')

  do i=1,m
     x=x1+(x2-x1)/real(m)*real(i)
     do j=1,4
        write(20+j,*) x+off, espandi(n,d(1,j),s,x,proj)*f
     enddo

  enddo
  close(21)
  close(22)
  close(23)
  close(24)
  deallocate(d,proj)
  stop
end program grafica


function espandi(n,d,s,x,proj)

  implicit none

  real :: espandi
  integer::n
  real :: d(n),proj(n)
  real::s
  real ::x


  integer :: i
  real :: xd

  espandi=0.
 
  do i=1,n
     xd=x-d(i)
     espandi=espandi+(1./(s*sqrt(2.*3.1415926)))*exp(-(xd**2.)/(2.*s**2.))*proj(i)
  enddo

  return
end function espandi
