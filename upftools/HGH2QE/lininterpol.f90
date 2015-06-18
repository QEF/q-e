subroutine interpolate(znuc,mesh2,amesh2,rr2,rhoatom2)
implicit real*8 (a-h,o-z)
dimension :: rr2(mesh2),rhoatom2(mesh2)
allocatable rr1(:),rhoatom1(:)
 z=znuc
 xmin = -7.0d0
 amesh1=0.0125d0  !! Modify this value to increase the no of grid points
 rmax =100.0d0
 mesh1 = 1 + (log(z*rmax)-xmin)/amesh1
 mesh1 = (mesh1/2)*2+1 ! mesh is odd (for historical reasons?)
 allocate(rr1(mesh1),rhoatom1(mesh1))
 do i=1,mesh1
   rr1(i) = exp (xmin+(i-1)*amesh1)/z
 end do
open(12,file='rhoatom.dat')
 ir=1
 do i=1,(mesh1/4)+1
  if ((ir+3).le.mesh1) then
    !write(*,*)"ir begin",ir
    read(12,*)(rhoatom1(j),j=ir,ir+3)
    !write(*,*)i,(rhoatom1(j),j=ir,ir+3)
    ir=j
    !write(*,*)"ir",ir
  else if ((ir+2).le.mesh1) then
    !write(*,*)"ir begin",ir
    read(12,*)(rhoatom1(j),j=ir,ir+2)
    !write(*,*)i,(rhoatom1(j),j=ir,ir+2)
    ir=j
    !write(*,*)"ir",ir
  else if ((ir+1).le.mesh1) then
    !write(*,*)"ir begin",ir
    read(12,*)(rhoatom1(j),j=ir,ir+1)
    !write(*,*)i,(rhoatom1(j),j=ir,ir+1)
    ir=j
    !write(*,*)"ir",ir
  else if (ir.le.mesh1) then
    !write(*,*)"ir begin",ir
    read(12,*)(rhoatom1(j),j=ir,ir+0)
    !write(*,*)i,(rhoatom1(j),j=ir,ir+0)
    ir=j
    !write(*,*)"ir",ir
  end if
 end do
close (12)


do i=1,mesh1
write(25,*)rr1(i),rhoatom1(i)
end do

rhoatom2=0.d0
!!do i2=1,mesh2
!! LINEAR INTERPOLATION
! do i1=1,mesh1
!  if((rr2(i2).lt.rr1(i1)))rlow=rr1(i1-1)
!  if((rr2(i2).lt.rr1(i1)))ilow=i1-1
!  if((rr2(i2).lt.rr1(i1)))exit
! end do
! do i1=mesh1,1
!  if((rr2(i2).gt.rr1(i1)))rhigh=rr1(i1+1)
!  if((rr2(i2).gt.rr1(i1)))ihigh=i1+1
!  if((rr2(i2).gt.rr1(i1)))exit
! end do
! dx1=rhigh-rlow
! dx2=rr2(i2)-rlow
! drhoatom=rhoatom1(ihigh)-rhoatom1(ilow)
! rhoatom2(i2)=rhoatom1(ilow)+drhoatom*dx2/dx1
! if (rhoatom2(i2).lt.0.d0)rhoatom2(i2)=0.d0

!! LAGRANGE INTERPOLATION
! do i1=1,mesh1
! do j1=1,mesh1
!   if (i1.ne.j1) then
!    dr2=rr2(i2)-rr1(j1)
!    dr1=rr1(i1)-rr1(j1)
!    dfact=dr2/dr1
    !!write(*,*)"dr1,dr2",dr1,dr2
!    rhoatom2(i2)=rhoatom2(i2)+rhoatom1(i1)*dfact
!   end if
! end do
! end do

!! SPLINE INTERPOLATION
do i1=2,mesh1
  do i2=1,mesh2
    if((rr2(i2).le.rr1(i1)).and.(rr2(i2).ge.rr1(i1-1))) then
      r1=rr1(i1-1)
      r2=rr1(i1)
      if((i1-1).eq.1)  then
         s1=rhoatom1(i1-1)/rr1(i1-1)
      else
         dy=rhoatom1(i1-1)-rhoatom1(i1-2)
         dx=rr1(i1-1)-rr1(i1-2)
         s1=dy/dx
      end if
      dy=rhoatom1(i1)-rhoatom1(i1-1)
      dx=rr1(i1)-rr1(i1-1)
      s2=dy/dx
      dx1=rr1(i1)-rr1(i1-1)
      dy1=rhoatom1(i1)-rhoatom1(i1-1)
      t=(rr2(i2)-rr1(i1))/dx1
      a=s1*dx1-dy1
      b=-s2*dx1+dy1
      rhoatom2(i2)=(1.d0-t)*rhoatom1(i1-1)+t*rhoatom1(i1)+t*(1.d0-t)*(a*(1.d0-t)+b*t)
    end if
 end do 
end do
end subroutine interpolate
