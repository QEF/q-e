!
!---------------------------------------------------------------
      subroutine lderivps
!---------------------------------------------------------------
!
!  numerical integration of the radial schroedinger equation 
!  computing logaritmic derivatives
!  thresh dermines the absolute accuracy for the eigenvalue
! 
!  nonlocal projectors are allowed
!
!
!
use ld1inc
implicit none

integer  ::       &
         lam           ! the angular momentum

real(kind=dp) ::  &
         ze2,     &    ! the nuclear charge in Ry units
         jam,     &    ! the total angular momentum
         e             ! the eigenvalue

integer ::      &
         ikrld, &       ! index of matching radius
         nc,    &       ! counter on logaritmic derivatives
         nbf,   &       ! number of b functions
         n,ie          ! generic counter 

real(kind=dp) :: &
         lamsq,            & ! combined angular momentum
         b(0:3),c(4),      & ! used for starting guess of the solution 
         b0e, rr1,rr2,     & ! auxiliary
         xl1, x4l6, ddx12, &
         x6l12, x8l20,     & !
         cb(ndm),j1(ndm),d(4),delta, compute_log, &
         eta(ndm,nwfx),bm(nwfx),cm(nwfx,nwfx), &
         vaux(ndm),        &
         int_0_inf_dr,     &  ! integral function
         coef(nwfx),xc(4), &  !
         dlchis(ndm,nwfsx)    ! the logarithmic derivatives

real(kind=dp),allocatable :: &
         aux(:),     &   ! the square of the wavefunction
         al(:)           ! the known part of the differential equation

integer :: &
         ib,jb,iib,jjb, &  ! counters on beta functions
         nst,nstop,     &  ! auxiliary for integrals
         ierr,          &  ! used to control allocation
         ios,           &  ! used for I/O control
         is, ind           ! counters on index

character(len=80) :: flld

if (nld.eq.0) return
if (nld.gt.nwfsx) call errore('lderivps','nld is too large',1)

allocate(al(mesh),stat=ierr)
allocate(aux(mesh),stat=ierr)

ze2=0.d0

ikrld= 1
do n=1,mesh
   if (r(n).lt.rlderiv) ikrld=n
enddo
write(6,'(5x,''Computing logaritmic derivative in'',f10.5)') &
             (r(ikrld)+r(ikrld+1))*0.5d0
if (ikrld.gt.mesh-1) &
    call errore('lderiv','ikrld is wrong',1)
do is=1,nspin
   do nc=1,nld
      if (rel.lt.2) then
         lam=nc-1
         jam=0.d0
      else
         lam=nc/2
         if (mod(nc,2)==0) jam=lam-0.5d0
         if (mod(nc,2)==1) jam=lam+0.5d0
      endif
      xl1=lam+1
      x4l6=4*lam+6
      x6l12=6*lam+12
      x8l20=8*lam+20
      ddx12=dx*dx/12.d0
      nst=(lam+1)**2  
      nbf=nbeta
      if (pseudotype.eq.1) then
          if (rel.eq.2) then
             if (abs(jam-lam+0.5d0).lt.1.d-2.or. lam==0 ) then
                ind=1
             else
                ind=2
             endif
             do n=1,mesh
                vpstot(n,is)=vpstot(n,is)+vnlo(n,lam,ind)
                vaux(n)=vnlo(n,lam,ind)
             enddo
          else
             do n=1,mesh
                vpstot(n,is)=vpstot(n,is)+vnl(n,lam)
                vaux(n)=vnl(n,lam)
             enddo
          endif
          nbf=0.d0
       endif
 
       do n=1,4
          al(n)=vpstot(n,is)-ze2/r(n)
       enddo
       call series(al,r,r2,b)
    
       npte= (emaxld-eminld)/deld 
       npte=npte+1
       if (npte.gt.ndm) &
          call errore('lderivps ', 'npte is too large ',1)
       do ie=1,npte
          e=eminld+deld*(ie-1.d0)
          lamsq=(lam+0.5d0)**2
!
!     b) find the value of solution s in the first two points
!
          b0e=b(0)-e
          c(1)=0.5*ze2/xl1
          c(2)=(c(1)*ze2+b0e)/x4l6
          c(3)=(c(2)*ze2+c(1)*b0e+b(1))/x6l12
          c(4)=(c(3)*ze2+c(2)*b0e+c(1)*b(1)+b(2))/x8l20
          rr1=(1.d0+r(1)*(c(1)+r(1)* &
                (c(2)+r(1)*(c(3)+r(1)*c(4)))))*r(1)**(lam+1)
          rr2=(1.d0+r(2)*(c(1)+r(2)* &
                 (c(2)+r(2)*(c(3)+r(2)*c(4)))))*r(2)**(lam+1)
          aux(1)=rr1/sqr(1)
          aux(2)=rr2/sqr(2)

          do n=1,mesh
             al(n)=( (vpstot(n,is)-e)*r2(n) + lamsq )*ddx12
             al(n)=1.d0-al(n)
          enddo

          call integrate_outward (lam,jam,e,mesh,ndm,dx,r,r2,sqr,al, &
                     b,aux,betas,ddd,qq,nbf,nwfsx,lls,jjs,ikrld+5)

!
!    compute the logaritmic derivative and save in dlchi
!            
         do n=-3,3
            aux(ikrld+n)= aux(ikrld+n)*sqr(ikrld+n)
         enddo

         dlchis(ie,nc)=compute_log(aux(ikrld-3),r(ikrld),dx)
      enddo
      if (pseudotype.eq.1) then
         do n=1,mesh
            vpstot(n,is)=vpstot(n,is)-vaux(n)
         enddo
      endif
   enddo

   if (file_logderps .ne. ' ') then
       flld=file_logderps
       if (is.eq.2) flld=trim(file_logderps)//'.01'
       open(unit=25,file=flld, status='unknown', iostat=ios, err=300 )
300    call errore('lderivps','opening file '//flld, abs(ios))

      do ie=1,npte
         e= eminld+deld*(ie-1)
         write(25,'(10f14.6)') e, (dlchis(ie,nc),nc=1,nld)
      enddo
      close(unit=25)
   endif
enddo

deallocate(aux)
deallocate(al)

return
end
