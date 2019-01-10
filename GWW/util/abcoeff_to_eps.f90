program series

  implicit none
  real(kind=8), ALLOCATABLE :: a(:,:),b(:,:)
  integer :: nstep,nn,nint
  real(kind=8) :: emin,emax,delta
  real(kind=8), parameter :: ry2ev=13.6057 
  


  integer :: i,j,k,idumm1, idumm2,nk
  real(kind=8) :: freq, fact
  complex(kind=8) :: z, f0,f1,eps, term


  real(kind=8):: omega=270.0114 
  real(kind=8) :: norm(3), nspin, offset
  


  write(*,*) 'Number of Lanczos steps tot, considered :'
  read(*,*) nstep
  nstep=nstep-1
  read(*,*) nint
  nint=nint-1
  write(*,*) 'Emin Emax Delta(eV) e_steps omega Nk nspin offset'

  read(*,*) emin
  read(*,*) emax
  read(*,*) delta
  read(*,*) nn

  read(*,*) omega
  read(*,*) nk

  read(*,*) nspin
  read(*,*) offset

  emin=emin/ry2ev
  emax=emax/ry2ev
  delta=delta/ry2eV

  allocate(b(nstep,3),a(nstep,3))
  open(unit=20,file='ab_coeff.dat',status='old')
  read(20,*) norm(1:3)
  do i=1,3
     do j=1,nstep
        read(20,*) idumm1,idumm2, a(j,i)
     enddo
  enddo
  do i=1,3
     do j=1,nstep
        read(20,*) idumm1,idumm2, b(j,i)
     enddo
  enddo
  close(20)

  do i=1,3!loop on cart
     if(i==1) then
        open(unit=20,file='lanczos.eps1x.dat',status='unknown')
        open(unit=21,file='lanczos.eps2x.dat',status='unknown')
     else if(i==2) then
        open(unit=20,file='lanczos.eps1y.dat',status='unknown')
        open(unit=21,file='lanczos.eps2y.dat',status='unknown')
     else if (i==3) then
        open(unit=20,file='lanczos.eps1z.dat',status='unknown')
        open(unit=21,file='lanczos.eps2z.dat',status='unknown')
     endif
      write(20,*) '# Energy(eV)   Eps1 Eps1'
      write(21,*) '# Energy(eV)   Eps2 Eps2'
     do j=1,nn!loop on freq
        freq=(emax-emin)/dble(nn)*dble(j)+emin
        z=cmplx(freq,delta)
   term=(0.5d0/b(nint,i))*(freq-a(nint,i)-sqrt((a(nint,i)-freq)**2.d0+4.d0*b(nint,i)**2.0))
        f0=term
        do k=nint-1,1,-1
           f1=z-a(k,i)-b(k,i)**2.d0/f0
           f0=f1
        enddo
        eps=1.d0-4.d0*3.14159264d0*2.d0*norm(i)*nspin/f0/(omega*dble(nk))
    
        write(20,*) freq*ry2ev+offset,real(eps),real(eps)
        write(21,*) freq*ry2ev+offset,aimag(eps),aimag(eps)
     enddo
     
     close(20)
     close(21)
  enddo

  deallocate(a,b)

  stop

end program series
