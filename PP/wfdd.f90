!
! Copyright (C) 2005  MANU/YUDONG WU/NICOLA MARZARI/ROBERTO CAR
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module wanpar
! nw:       the number of the G vectors
! nit:      the number of total iteration during searching
! nsd:      the number of steepest descent iterations
! ibrav:    the structure index, the same as ibrav in CP code.
  integer :: nw,  nit, nsd, ibrav
  logical adapt, restart
! wfdt:     time step during searching
! maxwfdt:  the maximum time step during searching
! b1,b2,b3: the reciprocal lattice
! alat:     the lattice parameter
! a1,a2,a3: the real-space lattice
  real(kind=8) :: wfdt, maxwfdt, b1(3), b2(3), b3(3), alat
  real(kind=8) :: a1(3), a2(3), a3(3), tolw

! wfg:      the G vectors involoved in general symmetry calculation
!           the units are b1, b2, b3.
!           For example: 
!            the ith G vector: wfg(i,1)*b1+wfg(i,2)*b2+wfg(i,3)*b3
  integer, allocatable :: wfg(:,:)

! weight:   the weight of each G vectors
  real(kind=8), allocatable :: weight(:)
!
!       These are the Input variables for Damped Dynamics
!
! q:            imaginary mass of the Unitary Matrix
! dt:           Time Step for damped dynamics
! cgordd:       1=conjugate gradient/SD 
!               any other number = damped dynamics
! fric:         damping coefficient, b/w 0 and 1
! nsteps:       Max No. of MD Steps
  real(kind=8) :: q, dt, fric
  integer :: cgordd, nsteps


end module wanpar

!----------------------------------------------------------------------
program wf
!----------------------------------------------------------------------
!
!    This program works on the overlap matrix calculated 
!        from parallel machine and search the unitary transformation
!        Uall corresponding to the Maximally localized Wannier functions.
!
!    The overlap matrix and lattice information are read from fort.38.
!
!
!    Searching parameters are in the input file:
!
!       cgordd  wfdt   maxwfdt   nit   nsd  q dt fric nsteps
!
!
!    The final unitary matrix Uall is output to fort.39.
!    Some running information is output to fort.24.
!
!                                            Yudong Wu 
!                                            June 28,2001
!
!       This code has been modified to include Damped dynamics to
!       find the maximally localized wannier functions.
!
!                                                Manu
!                                                September 16,2001
!
!
!     copyright MANU/YUDONG WU/NICOLA MARZARI/ROBERTO CAR
!
!----------------------------------------------------------------------
  use wanpar
!
  implicit none
  
  integer :: i, j, inw, n, nspin, nupdwn(2),ini
  complex(kind=8), allocatable :: O(:, :, :), Ospin(:, :, :)
  real(kind=8), allocatable :: Uall(:,:), Uspin(:,:), u1(:,:)

        read (5,*) cgordd, wfdt, maxwfdt, nit, nsd
        read (5,*)  q, dt, fric, adapt, nsteps, tolw
        read (5,*) restart


!
!    input the overlap matrix from fort.38
!
  rewind 38
  read(38, '(i5, 2i2, i3, f9.5)') n, nw, nspin, ibrav, alat
  allocate(wfg(nw, 3), weight(nw), O(nw,n,n), Uall(n,n), u1(n,n))
  if (nspin.eq.2) then
     read(38, '(i5)') nupdwn(1)
  end if
  nupdwn(2)=n-nupdwn(1)
  read(38, *) a1
  read(38, *) a2
  read(38, *) a3
  read(38, *) b1
  read(38, *) b2
  read(38, *) b3
  do inw=1, nw
     read(38, *) wfg(inw, :), weight(inw)
  end do
  do inw=1, nw
     do i=1, n
        do j=1, n
           read(38, *) O(inw, i, j)
        end do
     end do
  end do
  if(restart) then
  do i=1, n
     do j=1, n
        read(39, *) Uall(i, j)
     end do
  end do
  else
  Uall=0.0
     do i=1,n
        Uall(i,i)=1.d0
     end do 
  end if

  rewind 24

  if(cgordd.eq.1) then
    if (nspin.eq.1) then
      call searchwf(n, O, Uall)
    else
!
!   For those spin-polarized calculation,
!    spin up and spin down parts are dealt with seperately
!    and the total unitary matrices are put together.
!
     write(24, *) "spin up:"
     allocate(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
     do i=1, nupdwn(1)
        do j=1, nupdwn(1)
           Uspin(i, j)=Uall(i, j)
           Ospin(:, i, j)=O(:, i, j)
        end do
     end do
     call searchwf(nupdwn(1), Ospin, Uspin)
     do i=1, nupdwn(1)
        do j=1, nupdwn(1)
           Uall(i, j)=Uspin(i, j)
        end do
     end do
     deallocate(Uspin, Ospin)
     write(24, *) "spin down:"
     allocate(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
     do i=1, nupdwn(2)
        do j=1, nupdwn(2)
           Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
           Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
        end do
     end do
     call searchwf(nupdwn(2), Ospin, Uspin)
     do i=1, nupdwn(2)
        do j=1, nupdwn(2)
           Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
        end do
     end do
     deallocate(Uspin, Ospin)
    end if


  else

    if (nspin.eq.1) then
      call ddyn(n,O,Uall)
    else
!
!   For those spin-polarized calculation,
!    spin up and spin down parts are dealt with seperately
!    and the total unitary matrices are put together.
!
     write(24, *) "spin up:"
     allocate(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
     do i=1, nupdwn(1)
        do j=1, nupdwn(1)
           Uspin(i, j)=Uall(i, j)
           Ospin(:, i, j)=O(:, i, j)
        end do
     end do
     call ddyn(nupdwn(1), Ospin, Uspin)
     do i=1, nupdwn(1)
        do j=1, nupdwn(1)
           Uall(i, j)=Uspin(i, j)
        end do
     end do
     deallocate(Uspin, Ospin)
     write(24, *) "spin down:"
     allocate(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
     do i=1, nupdwn(2)
        do j=1, nupdwn(2)
           Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
           Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
        end do
     end do
     call ddyn(nupdwn(2), Ospin, Uspin)
     do i=1, nupdwn(2)
        do j=1, nupdwn(2)
           Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
        end do
     end do
     deallocate(Uspin, Ospin)
    end if
  endif



  rewind 39
  do i=1, n
     do j=1, n
        write(39, *) Uall(i, j)
     end do
  end do

!u1=matmul(Uall,transpose(Uall))

! do i=1, n
!     do j=1, n
!        write(6, *) u1(i, j)
!     end do
!  end do

  deallocate(wfg, weight, O, Uall,u1)

contains
!-------------------------------------------------------------------------
subroutine ddyn(m,Omat,Umat)
!    (m,m) is the size of the matrix Ospin.
!    Ospin is input overlap matrix.
!    Uspin is the output unitary transformation.
!             Rough guess for Uspin can be carried in.
!
!
!                                        MANU
!                                        SEPTEMBER 17, 2001
!-------------------------------------------------------------------------

  use wanpar

!  implicit none
  integer :: f3(nw), f4(nw), i,j,inw
  integer ,intent(in) :: m
  real(kind=8), intent(inout) :: Umat(m,m)
  complex(kind=8), intent(inout) :: Omat(nw,m,m)
  complex(kind=8) :: U2(m,m),U3(m,m)
  integer :: adjust,ini, ierr1,nnn
  real(kind=8), allocatable, dimension(:) :: wr
  real(kind=8), allocatable, dimension(:,:) :: W
  real(kind=8) :: t0, U(m,m), t2
  real(kind=8) :: A(m,m),oldt0,Wm(m,m),U1(m,m)
  real(kind=8) :: Aminus(m,m), Aplus(m,m),f2(3*m-2)
!  real(kind=8) :: Aminus(m,m), Aplus(m,m),f2(4*m)
  real(kind=8) :: temp(m,m)
  complex(kind=8) :: d(m,m), alpha, beta1, ci
  complex(kind=8) :: f1(2*m-1), wp(m*(m+1)/2),z(m,m)
  complex(kind=8), allocatable, dimension(:, :) :: X1
  complex(kind=8), allocatable, dimension(:, :, :) :: Oc
  real(kind=8) , allocatable , dimension(:) :: mt
  real(kind=8), parameter :: autoaf=0.529177d0
  real(kind=8) :: spread, sp, pi2, oldspread
  real(kind=8) :: wfc(3,n), gr(nw,3)

  alpha=(1.d0,0.d0)
  beta1=(0.d0,0.d0)
  ci   =(0.d0,1.d0)
  pi2  =2.d0*3.14159265358979d0

  allocate(mt(nw))
  allocate(X1(m,m))
  allocate(Oc(nw,m,m))

!   fric=friction
  allocate (W(m,m),wr(m))

!   Umat=0.d0
!   do i=1,m
!       Umat(i,i)=1.d0
!   end do

        U2=Umat*alpha

!
! update Oc using the initial guess of Uspin
!
  do inw=1, nw
    X1(:, :)=Omat(inw, :, :)
     U3=beta1
!    call ZGEMUL(U2, m, 'T', X1, m, 'N', U3, m, m,m,m) 
     call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
    X1=beta1
!    call ZGEMUL(U3, m, 'N', U2, m, 'N', X1, m, m,m,m) 
     call ZGEMM ('N','N', m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
    Oc(inw, :, :)=X1(:, :)
  end do

        U2=beta1
        U3=beta1

 oldspread=0.0
  write(24, *) "spread: (unit \AA^2)"
  do i=1, m
     mt=1.d0-DBLE(Oc(:,i,i)*conjg(Oc(:,i,i)))
     sp= (alat*autoaf/pi2)**2*SUM(mt*weight)
     write(24, '(f10.7)') (alat*autoaf/pi2)**2*SUM(mt*weight)
     oldspread=oldspread+sp
  end do

  oldspread=oldspread/m
     write(51, '(f10.7)') oldspread

    oldt0=0.d0
    A=0.d0
    Aminus=A
    temp=Aminus


!        START ITERATIONS HERE

  do ini=1, nsteps

    t0=0.d0     !use t0 to store the value of omega
    do inw=1, nw
       do i=1, m
          t0=t0+DBLE(conjg(Oc(inw, i, i))*Oc(inw, i, i))
       end do
    end do

    write(6, *) t0


        if(ABS(t0-oldt0).lt.tolw) then
           write(6,*) "MLWF Generated at Step",ini
           go to 241
        end if

        if(adapt) then 
          if(oldt0.lt.t0) then
            fric=fric/2.
            A=Aminus
            Aminus=temp
          end if
        end if

!   calculate d(omega)/dA and store result in W
!   this is the force for the damped dynamics
!

    W=0.d0
    do inw=1, nw
       t2=weight(inw)
       do i=1,m
          do j=1,m
             W(i,j)=W(i,j)+t2*DBLE(Oc(inw,i,j)*conjg(Oc(inw,i,i)        &
                  -Oc(inw,j,j))+conjg(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
          end do
       end do
    end do


!   the verlet scheme to calculate A(t+dt)

        Aplus=0.d0

   do i=1,m
     do j=i+1,m
         Aplus(i,j)=Aplus(i,j)+(2*dt/(2*dt+fric))*(2*A(i,j)               &
         -Aminus(i,j)+(dt*dt/q)*W(i,j)) + (fric/(2*dt+fric))*Aminus(i,j)
     enddo
   enddo

        Aplus=Aplus-transpose(Aplus)
        Aplus=(Aplus-A)

    do i=1, m
       do j=i,m 
        wp(i + (j-1)*j/2) = cmplx(0.0d0, Aplus(i,j))
       end do
    end do

    call zhpev('V','U',m,wp,wr,z,m,f1,f2,ierr1)
!    call zhpev(21, wp, wr, z, m, m, f2, 4*m)

    if (ierr1.ne.0) then 
   write(6,*) "failed to diagonalize W!"
    stop
    end if

    d=0.d0
    do i=1, m
       d(i, i)=exp(ci*wr(i)*dt)
    end do      !d=exp(d)

!   U=z*exp(d)*z+
!
     U3=beta1
     call ZGEMM ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)  
     U2=beta1
     call ZGEMM ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
    U=DBLE(U2)
    U2=beta1
    U3=beta1

   temp=Aminus
   Aminus=A
   A=Aplus


!   update Umat
!
    U1=beta1
    call DGEMM ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1 

!   update Oc
!
    U2=Umat*alpha
    U3=beta1
  do inw=1, nw
    X1(:, :)=Omat(inw, :, :)
    call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
    X1=beta1
    call ZGEMM ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
    Oc(inw, :, :)=X1(:, :)
  end do
    U2=beta1
    U3=beta1

        if(ABS(t0-oldt0).ge.tolw.and.ini.ge.nsteps) then
        go to 241
        end if

    oldt0=t0

   end do
241  spread=0.0
  write(24, *) "spread: (unit \AA^2)"
  do i=1, m
     mt=1.d0-DBLE(Oc(:,i,i)*conjg(Oc(:,i,i)))
     sp= (alat*autoaf/pi2)**2*SUM(mt*weight)
     write(24, '(f10.7)') (alat*autoaf/pi2)**2*SUM(mt*weight)
     spread=spread+sp
  end do

  spread=spread/m
     write(51, '(f10.7)') spread


  deallocate(wr, W)


!   output wfc's and spreads of the max. loc. wf's
!
  allocate(wr(nw), W(nw, nw))
  do inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  end do
!
! set up a matrix with the element (i,j) is G_i·G_j·weight(j)
! to check the correctness of choices on G vectors
!
  do i=1, nw
     do j=1, nw
        W(i,j)=SUM(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
     end do
  end do
!  write(24, *) "wannier function centers: (unit:\AA)"
  do i=1, m
     mt=-aimag(log(Oc(:,i,i)))/pi2
     wfc(1, i)=SUM(mt*weight*gr(:,1))
     wfc(2, i)=SUM(mt*weight*gr(:,2))
     wfc(3, i)=SUM(mt*weight*gr(:,3))
     do inw=1, nw
        wr(inw)=SUM(wfc(:,i)*gr(inw,:))-mt(inw)
     end do
     mt=wr
     f3=0
     adjust=0
!
!   balance the phase factor if necessary
!
!     do while(SUM((mt-f3)**2).gt.0.01d0)
!        f4=f3
!        f3=nint(mt-mt(1))
!        if (adjust.gt.200) f3=f3-1
!        if (adjust.gt.100.and.adjust.le.200) f3=f3+1
!        mt=wr+matmul(W, f3)
!        write(6,*) "mt:", mt
!        write(6,*) "f3:", f3
!        adjust=adjust+1
!        if (adjust.gt.300) stop "unable to balance the phase!"
!     end do
     wfc(1,i)=(wfc(1,i)+SUM(mt*weight*gr(:,1)))*alat
     wfc(2,i)=(wfc(2,i)+SUM(mt*weight*gr(:,2)))*alat
     wfc(3,i)=(wfc(3,i)+SUM(mt*weight*gr(:,3)))*alat
  end do
 
!  if (ibrav.eq.1.or.ibrav.eq.6.or.ibrav.eq.8) then
!     do i=1, m
!        if (wfc(1, i).lt.0) wfc(1, i)=wfc(1, i)+a1(1)
!        if (wfc(2, i).lt.0) wfc(2, i)=wfc(2, i)+a2(2)
!        if (wfc(3, i).lt.0) wfc(3, i)=wfc(3, i)+a3(3)
!     end do
!  end if
  do i=1, m
     write(26, '(3f11.6)') wfc(:,i)*autoaf
  end do
 
     write(6,*) "Friction =", fric
     write(6,*) "Mass =", q


  deallocate(wr, W)
 
  return

  end subroutine ddyn

!-----------------------------------------------------------------------
      subroutine searchwf(m, Omat, Umat)
!-----------------------------------------------------------------------
!    (m,m) is the size of the matrix Ospin.
!    Ospin is input overlap matrix.
!    Uspin is the output unitary transformation.
!             Rough guess for Uspin can be carried in.
!
  use wanpar

!
!
!     conjugated gradient to search maximization
!
  implicit none
!
  real(kind=8), parameter :: autoaf=0.529177d0
  integer, intent(in) :: m
  complex(kind=8), intent(in) :: Omat(nw, m, m)
  real(kind=8), intent(inout) :: Umat(m,m)
!
  integer :: i, j, k, l, ig, ierr, ti, tj, tk, inw, ir, adjust
  integer :: f3(nw), f4(nw), istep
  real(kind=8) :: slope, slope2, t1, t2, t3, pi2, mt(nw),t21,temp1,maxdt
  real(kind=8) :: U(m,m), wfc(3, m), Wm(m,m), schd(m,m), f2(3*m-2), gr(nw, 3)
  real(kind=8) :: Uspin2(m,m),temp2,wfdtold,oldt1,t01, d3(m,m), d4(m,m), U1(m,m)
  real(kind=8), allocatable, dimension(:) :: wr
  real(kind=8), allocatable, dimension(:,:) :: W
  complex(kind=8) :: ci, ct1, ct2, ct3, z(m, m), X(m, m), d(m,m), d2(m,m)
  complex(kind=8) :: f1(2*m-1), wp(m*(m+1)/2), Oc(nw, m, m), alpha, beta1
  complex(kind=8) ::  Oc2(nw, m, m),wp1(m*(m+1)/2), X1(m,m), U2(m,m), U3(m,m)

!
  ci=(0.d0,1.d0)
  alpha=(1.0d0, 0.0d0)
  beta1=(0.0d0, 0.0d0)
  pi2=2.d0*3.14159265358979d0
!
  allocate(W(m,m), wr(m))


!  Umat=0.d0
!  do i=1,m
!    Umat(i,i)=1.d0
!  end do
  Oc=beta1
  Oc2=beta1
  X1=beta1
  U2=Umat*alpha

!
! update Oc using the initial guess of Uspin
!
  do inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=beta1
     call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
     X1=beta1
     call ZGEMM ('N','N', m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
     Oc(inw, :, :)=X1(:, :)
  end do

     U2=beta1
     U3=beta1

  W=0.d0
  schd=0.d0
  oldt1=0.d0
  wfdtold=0.d0

  do k=1, nit
    t01=0.d0     !use t1 to store the value of omiga
    do inw=1, nw
       do i=1, m
          t01=t01+DBLE(conjg(Oc(inw, i, i))*Oc(inw, i, i))
       end do
    end do

    write(6,*) t01

    if(ABS(oldt1-t01).lt.tolw) go to 40
    
    oldt1=t01
  
!   calculate d(omiga)/dW and store result in W
!   W should be a real symmetric matrix for gamma-point calculation
!
    Wm=W
    W=0.d0
    do inw=1, nw
       t2=weight(inw)
       do i=1,m
          do j=i+1,m
             W(i,j)=W(i,j)+t2*DBLE(Oc(inw,i,j)*conjg(Oc(inw,i,i)        &
              -Oc(inw,j,j))+conjg(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
          end do
       end do
    end do
    W=W-transpose(W)
  
!   calculate slope=d(omiga)/d(lamda)
    slope=SUM(W**2)
  
!   calculate slope2=d2(omiga)/d(lamda)2
    slope2=0.d0
    do ti=1, m
       do tj=1, m
          do tk=1, m
             t2=0.d0
             do inw=1, nw
                t2=t2+DBLE(Oc(inw,tj,tk)*conjg(Oc(inw,tj,tj)+Oc(inw,tk,tk) &
                          -2.d0*Oc(inw,ti,ti))-4.d0*Oc(inw,ti,tk)          &
                          *conjg(Oc(inw,ti,tj)))*weight(inw)
             end do
             slope2=slope2+W(tk,ti)*W(ti,tj)*2.d0*t2
          end do
       end do
     end do
    slope2=2.d0*slope2
  
!   use parabola approximation. Defined by 1 point and 2 slopes
    if (slope2.lt.0) wfdt=-slope/2.d0/slope2
    if (maxwfdt.gt.0.and.wfdt.gt.maxwfdt) wfdt=maxwfdt
  
    if (k.lt.nsd) then
       schd=W    !use steepest-descent technique

!   calculate slope=d(omiga)/d(lamda)
    slope=SUM(schd**2)

!       schd=schd*maxwfdt
    do i=1, m
       do j=i, m
        wp1(i + (j-1)*j/2) = cmplx(0.0d0, schd(i,j))
       end do
    end do

    call zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
    if (ierr.ne.0) stop 'failed to diagonalize W!'

    else
!
!     construct conjugated gradient
!        d(i)=-g(i)+belta(i)*d(i-1)
!        belta^FR(i)=g(i)t*g(i)/(g(i-1)t*g(i-1))
!        belta^PR(i)=g(i)t*(g(i)-g(i-1))/(g(i-1)t*g(i-1))
!
        call DGEMM ('T','N', m,m,m,alpha,Wm,m,Wm,m,beta1,d3,m)

       
       t1=0.d0
       do i=1, m
          t1=t1+d3(i, i)
       end do
       if (t1.ne.0) then
          d4=(W-Wm)
          call DGEMM ('T','N', m,m,m,alpha,W,m,d4,m,beta1,d3,m)
          t2=0.d0
          do i=1, m
             t2=t2+d3(i, i)
          end do
          t3=t2/t1
          schd=W+schd*t3
       else
          schd=W
       end if
!
!        calculate the new d(Lambda) for the new Search Direction
!        added by Manu. September 19, 2001
!
!   calculate slope=d(omiga)/d(lamda)
    slope=SUM(schd**2)
!------------------------------------------------------------------------
!   schd=schd*maxwfdt
    do i=1, m
       do j=i, m
        wp1(i + (j-1)*j/2) = cmplx(0.0d0, schd(i,j))
       end do
    end do

    call zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
    if (ierr.ne.0) stop 'failed to diagonalize W!'

      maxdt=maxwfdt

11    d=0.d0
    do i=1, m
       d(i, i)=exp(ci*(maxwfdt)*wr(i))
    end do      !d=exp(d)
 
!   U=z*exp(d)*z+
     U3=beta1
     call ZGEMM ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     call ZGEMM ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
     U=DBLE(U2)
     U2=beta1
     U3=beta1
!
!   update Uspin
    U1=beta1
    call DGEMM ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1
 
!    Uspin2=matmul(Uspin, U2)
!
!   update Oc
!
     U2=Umat*alpha
     U3=beta1
     do inw=1, nw
      X1(:,:)=Omat(inw,:,:)
      call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
      X1=beta1
      call ZGEMM ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
      Oc2(inw, :,:)=X(:,:)
     end do
     U2=beta1 
     U3=beta1
!
    t21=0.d0     !use t21 to store the value of omiga
    do inw=1, nw
       do i=1, m
          t21=t21+DBLE(conjg(Oc2(inw, i, i))*Oc2(inw, i, i))
       end do
    end do
 
      temp1=-((t01-t21)+slope*maxwfdt)/(maxwfdt**2)
      temp2=slope
      wfdt=-temp2/(2*temp1)

        if (wfdt.gt.maxwfdt.or.wfdt.lt.0.d0) then
        maxwfdt=2*maxwfdt
        go to 11
        end if

        maxwfdt=maxdt
!
!
!   use parabola approximation. Defined by 2 point and 1 slopes
!    if (slope2.lt.0) wfdt=-slope/2.d0/slope2
!    if (maxwfdt.gt.0.and.wfdt.gt.maxwfdt) wfdt=maxwfdt
!
!    write(6, '(e12.5E2,1x,e11.5E2,1x,f6.2)') slope2, slope, wfdt
!-------------------------------------------------------------------------
!
!      schd is the new searching direction
!
    end if
  
    d=0.d0
    do i=1, m
       d(i, i)=exp(ci*wfdt*wr(i))
    end do          !d=exp(d)

 
!   U=z*exp(d)*z+
!
     U3=beta1
     call ZGEMM ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     call ZGEMM ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
     U=DBLE(U2)
     U2=beta1
     U3=beta1

!   update Uspin
!
    U1=beta1
    call DGEMM ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1

!   update Oc
!
       U2=Umat*alpha
       U3=beta1
     do inw=1, nw
       X1(:, :)=Omat(inw, :, :)
       call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
       X1=beta1
       call ZGEMM ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
       Oc(inw, :, :)=X1(:, :)
     end do
    U2=beta1
    U3=beta1
  end do

40  deallocate(W, wr)

!
! calculate the spread
!
  write(24, *) "spread: (unit \AA^2)"
  do i=1, m
     mt=1.d0-DBLE(Oc(:,i,i)*conjg(Oc(:,i,i)))
     write(24, '(f10.7)') (alat*autoaf/pi2)**2*SUM(mt*weight)
  end do

!
! calculate wannier-function centers
!
  allocate(wr(nw), W(nw, nw))
  do inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  end do
!
! set up a matrix with the element (i,j) is G_i·G_j·weight(j)
! to check the correctness of choices on G vectors
!
  do i=1, nw
     do j=1, nw
        W(i,j)=SUM(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
     end do
  end do
! write(24, *) "wannier function centers: (unit:\AA)"
  do i=1, m
     mt=-aimag(log(Oc(:,i,i)))/pi2
     wfc(1, i)=SUM(mt*weight*gr(:,1))
     wfc(2, i)=SUM(mt*weight*gr(:,2))
     wfc(3, i)=SUM(mt*weight*gr(:,3))
     do inw=1, nw
        wr(inw)=SUM(wfc(:,i)*gr(inw,:))-mt(inw)
     end do
     mt=wr
     wfc(1, i)=(wfc(1,i)+SUM(mt*weight*gr(:,1)))*alat
     wfc(2, i)=(wfc(2,i)+SUM(mt*weight*gr(:,2)))*alat
     wfc(3, i)=(wfc(3,i)+SUM(mt*weight*gr(:,3)))*alat
  end do


  do i=1, m
     write(26, '(3f11.6)') wfc(:,i)*autoaf
  end do
  deallocate(wr, W)
  return
  end subroutine searchwf


 end program wf
