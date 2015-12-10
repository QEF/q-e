!
! Copyright (C) 2005  MANU/YUDONG WU/NICOLA MARZARI/ROBERTO CAR
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE wanpar
! nw:       the number of the G vectors
! nit:      the number of total iteration during searching
! nsd:      the number of steepest descent iterations
! ibrav:    the structure index, the same as ibrav in CP code.
  INTEGER :: nw,  nit, nsd, ibrav
  LOGICAL adapt, restart
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
  INTEGER, ALLOCATABLE :: wfg(:,:)

! weight:   the weight of each G vectors
  real(kind=8), ALLOCATABLE :: weight(:)
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
  INTEGER :: cgordd, nsteps


END MODULE wanpar

!----------------------------------------------------------------------
PROGRAM wfdd
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
  USE wanpar
!
  IMPLICIT NONE

  INTEGER :: i, j, inw, n, nspin, nupdwn(2)
  COMPLEX(kind=8), ALLOCATABLE :: O(:, :, :), Ospin(:, :, :)
  real(kind=8), ALLOCATABLE :: Uall(:,:), Uspin(:,:), u1(:,:)

        READ (5,*) cgordd, wfdt, maxwfdt, nit, nsd
        READ (5,*)  q, dt, fric, adapt, nsteps, tolw
        READ (5,*) restart


!
!    input the overlap matrix from fort.38
!
  REWIND 38
  READ(38, '(i5, 2i2, i3, f9.5)') n, nw, nspin, ibrav, alat
  ALLOCATE(wfg(nw, 3), weight(nw), O(nw,n,n), Uall(n,n), u1(n,n))
  IF (nspin==2) THEN
     READ(38, '(i5)') nupdwn(1)
  ENDIF
  nupdwn(2)=n-nupdwn(1)
  READ(38, *) a1
  READ(38, *) a2
  READ(38, *) a3
  READ(38, *) b1
  READ(38, *) b2
  READ(38, *) b3
  DO inw=1, nw
     READ(38, *) wfg(inw, :), weight(inw)
  ENDDO
  DO inw=1, nw
     DO i=1, n
        DO j=1, n
           READ(38, *) O(inw, i, j)
        ENDDO
     ENDDO
  ENDDO
  IF(restart) THEN
  DO i=1, n
     DO j=1, n
        READ(39, *) Uall(i, j)
     ENDDO
  ENDDO
  ELSE
  Uall=0.0d0
     DO i=1,n
        Uall(i,i)=1.d0
     ENDDO
  ENDIF

  REWIND 24

  IF(cgordd==1) THEN
    IF (nspin==1) THEN
      CALL searchwf(n, O, Uall)
    ELSE
!
!   For those spin-polarized calculation,
!    spin up and spin down parts are dealt with seperately
!    and the total unitary matrices are put together.
!
     WRITE(24, *) "spin up:"
     ALLOCATE(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
     DO i=1, nupdwn(1)
        DO j=1, nupdwn(1)
           Uspin(i, j)=Uall(i, j)
           Ospin(:, i, j)=O(:, i, j)
        ENDDO
     ENDDO
     CALL searchwf(nupdwn(1), Ospin, Uspin)
     DO i=1, nupdwn(1)
        DO j=1, nupdwn(1)
           Uall(i, j)=Uspin(i, j)
        ENDDO
     ENDDO
     DEALLOCATE(Uspin, Ospin)
     WRITE(24, *) "spin down:"
     ALLOCATE(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
     DO i=1, nupdwn(2)
        DO j=1, nupdwn(2)
           Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
           Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
        ENDDO
     ENDDO
     CALL searchwf(nupdwn(2), Ospin, Uspin)
     DO i=1, nupdwn(2)
        DO j=1, nupdwn(2)
           Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
        ENDDO
     ENDDO
     DEALLOCATE(Uspin, Ospin)
    ENDIF


  ELSE

    IF (nspin==1) THEN
      CALL ddyn(n,O,Uall)
    ELSE
!
!   For those spin-polarized calculation,
!    spin up and spin down parts are dealt with seperately
!    and the total unitary matrices are put together.
!
     WRITE(24, *) "spin up:"
     ALLOCATE(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
     DO i=1, nupdwn(1)
        DO j=1, nupdwn(1)
           Uspin(i, j)=Uall(i, j)
           Ospin(:, i, j)=O(:, i, j)
        ENDDO
     ENDDO
     CALL ddyn(nupdwn(1), Ospin, Uspin)
     DO i=1, nupdwn(1)
        DO j=1, nupdwn(1)
           Uall(i, j)=Uspin(i, j)
        ENDDO
     ENDDO
     DEALLOCATE(Uspin, Ospin)
     WRITE(24, *) "spin down:"
     ALLOCATE(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
     DO i=1, nupdwn(2)
        DO j=1, nupdwn(2)
           Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
           Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
        ENDDO
     ENDDO
     CALL ddyn(nupdwn(2), Ospin, Uspin)
     DO i=1, nupdwn(2)
        DO j=1, nupdwn(2)
           Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
        ENDDO
     ENDDO
     DEALLOCATE(Uspin, Ospin)
    ENDIF
  ENDIF



  REWIND 39
  DO i=1, n
     DO j=1, n
        WRITE(39, *) Uall(i, j)
     ENDDO
  ENDDO

!u1=matmul(Uall,transpose(Uall))

! do i=1, n
!     do j=1, n
!        write(6, *) u1(i, j)
!     end do
!  end do

  DEALLOCATE(wfg, weight, O, Uall,u1)

CONTAINS
!-------------------------------------------------------------------------
SUBROUTINE ddyn(m,Omat,Umat)
!    (m,m) is the size of the matrix Ospin.
!    Ospin is input overlap matrix.
!    Uspin is the output unitary transformation.
!             Rough guess for Uspin can be carried in.
!
!
!                                        MANU
!                                        SEPTEMBER 17, 2001
!-------------------------------------------------------------------------

  USE wanpar

  USE constants, ONLY : tpi, autoaf => BOHR_RADIUS_ANGS
!  implicit none
  INTEGER :: f3(nw), f4(nw), i,j,inw
  INTEGER ,INTENT(in) :: m
  real(kind=8), INTENT(inout) :: Umat(m,m)
  COMPLEX(kind=8), INTENT(inout) :: Omat(nw,m,m)
  COMPLEX(kind=8) :: U2(m,m),U3(m,m)
  INTEGER :: adjust,ini, ierr1
  real(kind=8), ALLOCATABLE, DIMENSION(:) :: wr
  real(kind=8), ALLOCATABLE, DIMENSION(:,:) :: W
  real(kind=8) :: t0, U(m,m), t2
  real(kind=8) :: A(m,m),oldt0,Wm(m,m),U1(m,m)
  real(kind=8) :: Aminus(m,m), Aplus(m,m),f2(3*m-2)
!  real(kind=8) :: Aminus(m,m), Aplus(m,m),f2(4*m)
  real(kind=8) :: temp(m,m)
  COMPLEX(kind=8) :: d(m,m), alpha, beta1, ci
  COMPLEX(kind=8) :: f1(2*m-1), wp(m*(m+1)/2),z(m,m)
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:, :) :: X1
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:, :, :) :: Oc
  real(kind=8) , ALLOCATABLE , DIMENSION(:) :: mt
  real(kind=8) :: spread, sp, oldspread
  real(kind=8) :: wfc(3,n), gr(nw,3)

  alpha=(1.d0,0.d0)
  beta1=(0.d0,0.d0)
  ci   =(0.d0,1.d0)

  ALLOCATE(mt(nw))
  ALLOCATE(X1(m,m))
  ALLOCATE(Oc(nw,m,m))

!   fric=friction
  ALLOCATE (W(m,m),wr(m))

!   Umat=0.d0
!   do i=1,m
!       Umat(i,i)=1.d0
!   end do

        U2=Umat*alpha

!
! update Oc using the initial guess of Uspin
!
  DO inw=1, nw
    X1(:, :)=Omat(inw, :, :)
     U3=beta1
!    call ZGEMUL(U2, m, 'T', X1, m, 'N', U3, m, m,m,m)
     CALL zgemm ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
    X1=beta1
!    call ZGEMUL(U3, m, 'N', U2, m, 'N', X1, m, m,m,m)
     CALL zgemm ('N','N', m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
    Oc(inw, :, :)=X1(:, :)
  ENDDO

        U2=beta1
        U3=beta1

 oldspread=0.0d0
  WRITE(24, *) "spread: (unit \AA^2)"
  DO i=1, m
     mt=1.d0-dble(Oc(:,i,i)*conjg(Oc(:,i,i)))
     sp= (alat*autoaf/tpi)**2*sum(mt*weight)
     WRITE(24, '(f10.7)') (alat*autoaf/tpi)**2*sum(mt*weight)
     oldspread=oldspread+sp
  ENDDO

  oldspread=oldspread/m
     WRITE(51, '(f10.7)') oldspread

    oldt0=0.d0
    A=0.d0
    Aminus=A
    temp=Aminus


!        START ITERATIONS HERE

  DO ini=1, nsteps

    t0=0.d0     !use t0 to store the value of omega
    DO inw=1, nw
       DO i=1, m
          t0=t0+dble(conjg(Oc(inw, i, i))*Oc(inw, i, i))
       ENDDO
    ENDDO

    WRITE(6, *) t0


        IF(abs(t0-oldt0)<tolw) THEN
           WRITE(6,*) "MLWF Generated at Step",ini
           GOTO 241
        ENDIF

        IF(adapt) THEN
          IF(oldt0<t0) THEN
            fric=fric/2.d0
            A=Aminus
            Aminus=temp
          ENDIF
        ENDIF

!   calculate d(omega)/dA and store result in W
!   this is the force for the damped dynamics
!

    W=0.d0
    DO inw=1, nw
       t2=weight(inw)
       DO i=1,m
          DO j=1,m
             W(i,j)=W(i,j)+t2*dble(Oc(inw,i,j)*conjg(Oc(inw,i,i)        &
                  -Oc(inw,j,j))+conjg(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
          ENDDO
       ENDDO
    ENDDO


!   the verlet scheme to calculate A(t+dt)

        Aplus=0.d0

   DO i=1,m
     DO j=i+1,m
         Aplus(i,j)=Aplus(i,j)+(2*dt/(2*dt+fric))*(2*A(i,j)               &
         -Aminus(i,j)+(dt*dt/q)*W(i,j)) + (fric/(2*dt+fric))*Aminus(i,j)
     ENDDO
   ENDDO

        Aplus=Aplus-transpose(Aplus)
        Aplus=(Aplus-A)

    DO i=1, m
       DO j=i,m
        wp(i + (j-1)*j/2) = cmplx(0.0d0, Aplus(i,j), kind=8)
       ENDDO
    ENDDO

    CALL zhpev('V','U',m,wp,wr,z,m,f1,f2,ierr1)
!    call zhpev(21, wp, wr, z, m, m, f2, 4*m)

    IF (ierr1/=0) THEN
   WRITE(6,*) "failed to diagonalize W!"
    STOP
    ENDIF

    d=0.d0
    DO i=1, m
       d(i, i)=exp(ci*wr(i)*dt)
    ENDDO      !d=exp(d)

!   U=z*exp(d)*z+
!
     U3=beta1
     CALL zgemm ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     CALL zgemm ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
    U=dble(U2)
    U2=beta1
    U3=beta1

   temp=Aminus
   Aminus=A
   A=Aplus


!   update Umat
!
    U1=beta1
    CALL dgemm ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1

!   update Oc
!
    U2=Umat*alpha
    U3=beta1
  DO inw=1, nw
    X1(:, :)=Omat(inw, :, :)
    CALL zgemm ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
    X1=beta1
    CALL zgemm ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
    Oc(inw, :, :)=X1(:, :)
  ENDDO
    U2=beta1
    U3=beta1

        IF(abs(t0-oldt0)>=tolw.and.ini>=nsteps) THEN
        GOTO 241
        ENDIF

    oldt0=t0

   ENDDO
241  spread=0.0d0
  WRITE(24, *) "spread: (unit \AA^2)"
  DO i=1, m
     mt=1.d0-dble(Oc(:,i,i)*conjg(Oc(:,i,i)))
     sp= (alat*autoaf/tpi)**2*sum(mt*weight)
     WRITE(24, '(f10.7)') (alat*autoaf/tpi)**2*sum(mt*weight)
     spread=spread+sp
  ENDDO

  spread=spread/m
     WRITE(51, '(f10.7)') spread


  DEALLOCATE(wr, W)


!   output wfc's and spreads of the max. loc. wf's
!
  ALLOCATE(wr(nw), W(nw, nw))
  DO inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  ENDDO
!
! set up a matrix with the element (i,j) is G_i * G_j * weight(j)
! to check the correctness of choices on G vectors
!
  DO i=1, nw
     DO j=1, nw
        W(i,j)=sum(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
     ENDDO
  ENDDO
!  write(24, *) "wannier function centers: (unit:\AA)"
  DO i=1, m
     mt=-aimag(log(Oc(:,i,i)))/tpi
     wfc(1, i)=sum(mt*weight*gr(:,1))
     wfc(2, i)=sum(mt*weight*gr(:,2))
     wfc(3, i)=sum(mt*weight*gr(:,3))
     DO inw=1, nw
        wr(inw)=sum(wfc(:,i)*gr(inw,:))-mt(inw)
     ENDDO
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
     wfc(1,i)=(wfc(1,i)+sum(mt*weight*gr(:,1)))*alat
     wfc(2,i)=(wfc(2,i)+sum(mt*weight*gr(:,2)))*alat
     wfc(3,i)=(wfc(3,i)+sum(mt*weight*gr(:,3)))*alat
  ENDDO

!  if (ibrav.eq.1.or.ibrav.eq.6.or.ibrav.eq.8) then
!     do i=1, m
!        if (wfc(1, i).lt.0) wfc(1, i)=wfc(1, i)+a1(1)
!        if (wfc(2, i).lt.0) wfc(2, i)=wfc(2, i)+a2(2)
!        if (wfc(3, i).lt.0) wfc(3, i)=wfc(3, i)+a3(3)
!     end do
!  end if
  DO i=1, m
     WRITE(26, '(3f11.6)') wfc(:,i)*autoaf
  ENDDO

     WRITE(6,*) "Friction =", fric
     WRITE(6,*) "Mass =", q


  DEALLOCATE(wr, W)

  RETURN

  END SUBROUTINE ddyn

!-----------------------------------------------------------------------
  SUBROUTINE searchwf(m, Omat, Umat)
!-----------------------------------------------------------------------
!    (m,m) is the size of the matrix Ospin.
!    Ospin is input overlap matrix.
!    Uspin is the output unitary transformation.
!             Rough guess for Uspin can be carried in.
!
  USE wanpar
  USE constants, ONLY : tpi, autoaf => BOHR_RADIUS_ANGS
!
!
!     conjugated gradient to search maximization
!
  IMPLICIT NONE
!
  INTEGER, INTENT(in) :: m
  COMPLEX(kind=8), INTENT(in) :: Omat(nw, m, m)
  real(kind=8), INTENT(inout) :: Umat(m,m)
!
  INTEGER :: i, j, k, l, ig, ierr, ti, tj, tk, inw, ir
  real(kind=8) :: slope, slope2, t1, t2, t3, mt(nw),t21,temp1,maxdt
  real(kind=8) :: U(m,m), wfc(3, m), Wm(m,m), schd(m,m), f2(3*m-2), gr(nw, 3)
  real(kind=8) :: Uspin2(m,m),temp2,wfdtold,oldt1,t01, d3(m,m), d4(m,m), U1(m,m)
  real(kind=8), ALLOCATABLE, DIMENSION(:) :: wr
  real(kind=8), ALLOCATABLE, DIMENSION(:,:) :: W
  COMPLEX(kind=8) :: ci, ct1, ct2, ct3, z(m, m), X(m, m), d(m,m), d2(m,m)
  COMPLEX(kind=8) :: f1(2*m-1), wp(m*(m+1)/2), Oc(nw, m, m), alpha, beta1
  COMPLEX(kind=8) ::  Oc2(nw, m, m),wp1(m*(m+1)/2), X1(m,m), U2(m,m), U3(m,m)

!
  ci=(0.d0,1.d0)
  alpha=(1.0d0, 0.0d0)
  beta1=(0.0d0, 0.0d0)
!
  ALLOCATE(W(m,m), wr(m))


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
  DO inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=beta1
     CALL zgemm ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
     X1=beta1
     CALL zgemm ('N','N', m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
     Oc(inw, :, :)=X1(:, :)
  ENDDO

     U2=beta1
     U3=beta1

  W=0.d0
  schd=0.d0
  oldt1=0.d0
  wfdtold=0.d0

  DO k=1, nit
    t01=0.d0     !use t1 to store the value of omiga
    DO inw=1, nw
       DO i=1, m
          t01=t01+dble(conjg(Oc(inw, i, i))*Oc(inw, i, i))
       ENDDO
    ENDDO

    WRITE(6,*) t01

    IF(abs(oldt1-t01)<tolw) GOTO 40

    oldt1=t01

!   calculate d(omiga)/dW and store result in W
!   W should be a real symmetric matrix for gamma-point calculation
!
    Wm=W
    W=0.d0
    DO inw=1, nw
       t2=weight(inw)
       DO i=1,m
          DO j=i+1,m
             W(i,j)=W(i,j)+t2*dble(Oc(inw,i,j)*conjg(Oc(inw,i,i)        &
              -Oc(inw,j,j))+conjg(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
          ENDDO
       ENDDO
    ENDDO
    W=W-transpose(W)

!   calculate slope=d(omiga)/d(lamda)
    slope=sum(W**2)

!   calculate slope2=d2(omiga)/d(lamda)2
    slope2=0.d0
    DO ti=1, m
       DO tj=1, m
          DO tk=1, m
             t2=0.d0
             DO inw=1, nw
                t2=t2+dble(Oc(inw,tj,tk)*conjg(Oc(inw,tj,tj)+Oc(inw,tk,tk) &
                          -2.d0*Oc(inw,ti,ti))-4.d0*Oc(inw,ti,tk)          &
                          *conjg(Oc(inw,ti,tj)))*weight(inw)
             ENDDO
             slope2=slope2+W(tk,ti)*W(ti,tj)*2.d0*t2
          ENDDO
       ENDDO
     ENDDO
    slope2=2.d0*slope2

!   use parabola approximation. Defined by 1 point and 2 slopes
    IF (slope2<0) wfdt=-slope/2.d0/slope2
    IF (maxwfdt>0.and.wfdt>maxwfdt) wfdt=maxwfdt

    IF (k<nsd) THEN
       schd=W    !use steepest-descent technique

!   calculate slope=d(omiga)/d(lamda)
    slope=sum(schd**2)

!       schd=schd*maxwfdt
    DO i=1, m
       DO j=i, m
        wp1(i + (j-1)*j/2) = cmplx(0.0d0, schd(i,j), kind=8)
       ENDDO
    ENDDO

    CALL zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
    IF (ierr/=0) STOP 'failed to diagonalize W!'

    ELSE
!
!     construct conjugated gradient
!        d(i)=-g(i)+belta(i)*d(i-1)
!        belta^FR(i)=g(i)t*g(i)/(g(i-1)t*g(i-1))
!        belta^PR(i)=g(i)t*(g(i)-g(i-1))/(g(i-1)t*g(i-1))
!
        CALL dgemm ('T','N', m,m,m,alpha,Wm,m,Wm,m,beta1,d3,m)


       t1=0.d0
       DO i=1, m
          t1=t1+d3(i, i)
       ENDDO
       IF (t1/=0) THEN
          d4=(W-Wm)
          CALL dgemm ('T','N', m,m,m,alpha,W,m,d4,m,beta1,d3,m)
          t2=0.d0
          DO i=1, m
             t2=t2+d3(i, i)
          ENDDO
          t3=t2/t1
          schd=W+schd*t3
       ELSE
          schd=W
       ENDIF
!
!        calculate the new d(Lambda) for the new Search Direction
!        added by Manu. September 19, 2001
!
!   calculate slope=d(omiga)/d(lamda)
    slope=sum(schd**2)
!------------------------------------------------------------------------
!   schd=schd*maxwfdt
    DO i=1, m
       DO j=i, m
        wp1(i + (j-1)*j/2) = cmplx(0.0d0, schd(i,j), kind=8)
       ENDDO
    ENDDO

    CALL zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
    IF (ierr/=0) STOP 'failed to diagonalize W!'

      maxdt=maxwfdt

11    d=0.d0
    DO i=1, m
       d(i, i)=exp(ci*(maxwfdt)*wr(i))
    ENDDO      !d=exp(d)

!   U=z*exp(d)*z+
     U3=beta1
     CALL zgemm ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     CALL zgemm ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
     U=dble(U2)
     U2=beta1
     U3=beta1
!
!   update Uspin
    U1=beta1
    CALL dgemm ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1

!    Uspin2=matmul(Uspin, U2)
!
!   update Oc
!
     U2=Umat*alpha
     U3=beta1
     DO inw=1, nw
      X1(:,:)=Omat(inw,:,:)
      CALL zgemm ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
      X1=beta1
      CALL zgemm ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
      Oc2(inw, :,:)=X(:,:)
     ENDDO
     U2=beta1
     U3=beta1
!
    t21=0.d0     !use t21 to store the value of omiga
    DO inw=1, nw
       DO i=1, m
          t21=t21+dble(conjg(Oc2(inw, i, i))*Oc2(inw, i, i))
       ENDDO
    ENDDO

      temp1=-((t01-t21)+slope*maxwfdt)/(maxwfdt**2)
      temp2=slope
      wfdt=-temp2/(2*temp1)

        IF (wfdt>maxwfdt.or.wfdt<0.d0) THEN
        maxwfdt=2*maxwfdt
        GOTO 11
        ENDIF

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
    ENDIF

    d=0.d0
    DO i=1, m
       d(i, i)=exp(ci*wfdt*wr(i))
    ENDDO          !d=exp(d)


!   U=z*exp(d)*z+
!
     U3=beta1
     CALL zgemm ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     CALL zgemm ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
     U=dble(U2)
     U2=beta1
     U3=beta1

!   update Uspin
!
    U1=beta1
    CALL dgemm ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1

!   update Oc
!
       U2=Umat*alpha
       U3=beta1
     DO inw=1, nw
       X1(:, :)=Omat(inw, :, :)
       CALL zgemm ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
       X1=beta1
       CALL zgemm ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
       Oc(inw, :, :)=X1(:, :)
     ENDDO
    U2=beta1
    U3=beta1
  ENDDO

40  DEALLOCATE(W, wr)

!
! calculate the spread
!
  WRITE(24, *) "spread: (unit \AA^2)"
  DO i=1, m
     mt=1.d0-dble(Oc(:,i,i)*conjg(Oc(:,i,i)))
     WRITE(24, '(f10.7)') (alat*autoaf/tpi)**2*sum(mt*weight)
  ENDDO

!
! calculate wannier-function centers
!
  ALLOCATE(wr(nw), W(nw, nw))
  DO inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  ENDDO
!
! set up a matrix with the element (i,j) is G_i * G_j * weight(j)
! to check the correctness of choices on G vectors
!
  DO i=1, nw
     DO j=1, nw
        W(i,j)=sum(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
     ENDDO
  ENDDO
! write(24, *) "wannier function centers: (unit:\AA)"
  DO i=1, m
     mt=-aimag(log(Oc(:,i,i)))/tpi
     wfc(1, i)=sum(mt*weight*gr(:,1))
     wfc(2, i)=sum(mt*weight*gr(:,2))
     wfc(3, i)=sum(mt*weight*gr(:,3))
     DO inw=1, nw
        wr(inw)=sum(wfc(:,i)*gr(inw,:))-mt(inw)
     ENDDO
     mt=wr
     wfc(1, i)=(wfc(1,i)+sum(mt*weight*gr(:,1)))*alat
     wfc(2, i)=(wfc(2,i)+sum(mt*weight*gr(:,2)))*alat
     wfc(3, i)=(wfc(3,i)+sum(mt*weight*gr(:,3)))*alat
  ENDDO


  DO i=1, m
     WRITE(26, '(3f11.6)') wfc(:,i)*autoaf
  ENDDO
  DEALLOCATE(wr, W)
  RETURN
  END SUBROUTINE searchwf


 END PROGRAM wfdd
