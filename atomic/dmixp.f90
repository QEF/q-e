!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
subroutine dmixp (nmsh,a,b,beta,tr2,in,id,r2,conv,itmax)
  !-------------------------------------------------------------------
  !
  !     dmixp - written by k.c. pandey, apr. 1981. 
  !     (heavily modified version - PG - last modified march 1996)
  !
  !     using anderson's extrapolation scheme
  !     eqs. 4.1 - 4.9, eqs. 4.15 - 4.18 of:
  !     donald g. anderson, j. assoc. computing machinery, 12, 547 (1965)
  !     computes input vector for the next iteration
  !
  !---- input -------
  !     nmsh = size of a and b
  !     a = y(l) = new vector (potential)
  !     b = x(l) = old vector
  !     beta = mixing factor
  !        in simple schemes one uses: new input = (1-beta)*b + beta*a
  !     tr2 = convergence factor
  !        if r2 < tr2  we assume convergence has been achieved
  !        no new input vector is calculated
  !     in = iteration number
  !        even if id = 2 or 3, for in = 1
  !        calculations are done as if id = 1
  !        for in = 2 and id = 3, calculations are done as if id = 2
  !     id = 1,2, or 3
  !        depending on whether we want to use: only a and b; a, b and
  !        the last iteration; a, b, and the last two iterations
  !---- output --------
  !     a = y(l)-x(l)
  !     b = x(l+1) = new input vector (potential)
  !     r2 = mean-square of (a(i)-b(i))**2

  use io_global, only : stdout
  use kinds, only : DP
  implicit none
  integer :: nmsh, id, in, ierr, itmax
  logical :: conv
  real(DP):: a (nmsh), b (nmsh), r2, tr2, beta

  ! detol: if det < d11*d22*detol we assume that the det = 0
  !        in this case only last iteration is used for
  !        computation of the new input vector

  real(DP), parameter:: detol =1.0e-9_DP
  real(DP):: t1,t2, d11,d22,d12, aa2,det,dett,rd1m,rd2m, ddot
  real(DP),allocatable:: c(:), d(:), a1(:), a2(:), b1(:), b2(:)
  integer:: i
  external ddot
  save c, d, a1, a2, b1, b2


  ! a = y(l) = V_out(l), b = x(l) = V_in(l)
  conv=.false.
  call trns(a,b,-1.0_DP,nmsh)
  ! Now: a = r(l) == (V_out-V_in)(l)
  r2 = ddot(nmsh,a,1,a,1)/DBLE(nmsh)

  if(r2.lt.tr2.or.in==itmax) then
     if (r2.lt.tr2) conv=.true.
     ! Convergence achieved or last iteration: Deallocate working space
     if (id.eq.3.and.in.ne.1) then
        deallocate(b2)
        deallocate(a2)
        deallocate(d )
     end if
     if (id.eq.2 .or. id.eq.3.and.in.ne.1) then 
        deallocate(b1)
        deallocate(a1)
        deallocate(c )
     end if
     return
  end if
  !
  if(id.eq.1) then

     call trns(b,a,beta,nmsh)
     ! id=1, simple mixing: b = beta*a+ (1-beta)*b
     return

  else if(id.lt.1.or.id.gt.3.or.in.lt.1) then
     write(stdout,'('' *** stop in dmixp *** id,in = '',2i6 )') id,in
  end if

  if (in.eq.1) then
     ! First iteration: allocate working space
     if (id.eq.2 .or. id.eq.3) then 
        allocate(c(nmsh),stat=ierr)
        allocate(a1(nmsh),stat=ierr)
        allocate(b1(nmsh),stat=ierr)
     end if
     if (id.eq.3) then
        allocate(d(nmsh),stat=ierr)
        allocate(a2(nmsh),stat=ierr)
        allocate(b2(nmsh),stat=ierr)
     end if

     call dcopy(nmsh,a,1,a1,1)
     call dcopy(nmsh,b,1,b1,1)
     ! a1 = r(l-1),  b1 = V_in(l-1)
     call trns(b,a,beta,nmsh)
     ! Simple mixing: b = beta*a+ (1-beta)*b
     return
  end if

  do i=1,nmsh
     c(i)=a1(i)-a(i)
  end do
  ! Now: c = r(l-1)-r(l)

  d11 =ddot(nmsh,c,1,c,1)
  rd1m=ddot(nmsh,a,1,c,1)

  if(id.eq.3.and.in.gt.2) then

     do i=1,nmsh
        d(i)=a2(i)-a(i)
     end do
     ! Now: d = r(l-2)-r(l)

     d22 =ddot(nmsh,d,1,d,1)
     d12 =ddot(nmsh,c,1,d,1)
     rd2m=ddot(nmsh,a,1,d,1)
     aa2=d11*d22
     det=aa2-d12*d12
     dett=det/aa2
     if(abs(dett).lt.detol) then
        t1=-rd1m/d11
        t2=0.0_dp
     else
        t1=(-rd1m*d22+rd2m*d12)/det
        t2=( rd1m*d12-rd2m*d11)/det
     end if
     ! Write the new input vector on c
     do i=1,nmsh
        c(i) = (1.0_dp-t1-t2)*b(i) + ( t1*b1(i) + t2*b2(i) ) &
             + beta * ( a(i) + t1*c(i) + t2*d(i) )
     end do

  else

     ! t1 = theta1   (eq 4.8)
     t1 = -rd1m/d11
     ! Write the new input vector on c
     do i=1,nmsh
        c(i) = (1.0_dp-t1)*b(i) + t1*b1(i) + beta*(a(i)+t1*c(i))
     end do

  end if

  ! Save a1=r(l-1), a2=r(l-2), b1=V_in(l-1), b2=V_in(l-2), for next iteration

  if(id.eq.3) then
     call dcopy(nmsh,b1,1,b2,1)
     call dcopy(nmsh,a1,1,a2,1)
  end if
  call dcopy(nmsh,a,1,a1,1)
  call dcopy(nmsh,b,1,b1,1)

  ! Copy the new input vector on b

  call dcopy(nmsh,c,1,b,1)

  return
end subroutine dmixp
!
!     trns calculates a(i)=a(i)+c*b(i) for i=1,n
!     b(i) and n remain unchanged
!
subroutine trns(a,b,c,n)

  use kinds, only : DP
  implicit none
  integer:: n,i
  real(DP):: a(n),b(n),c
  do i=1,n
     a(i)=a(i)+c*b(i)
  end do
  return
end subroutine trns
