!
!---------------------------------------------------------------
      subroutine lderiv
!---------------------------------------------------------------
!
!  numerical integration of the radial schroedinger equation 
!  computing logarithmic derivatives
!  thresh dermines the absolute accuracy for the eigenvalue
!
!
!
use ld1inc

implicit none

integer ::        &
         lam         ! the angular momentum

real(kind=dp) ::  &
         ze2,     &  ! the nuclear charge in Ry units
         e           ! the eigenvalue

integer ::       &
         ikrld,  &   ! index of matching radius
         nc,     &   ! counter on logarithmic derivatives
         idum,   &   ! integer variable for lschps
         is,     &   ! counter on spin
         ierr,   &   ! used for allocation control
         ios,    &   ! used for I/O control
         n,ie        ! generic counter 

real(kind=dp) ::           &
         rab(ndm),         & ! derivative of the radial mesh
         dlchi(ndm,nwfsx), & ! the logarithmic derivative
         j,                & ! total angular momentum for log_der
         compute_log 

real(kind=dp) ::        &
         aux(ndm),       &  ! the square of the wavefunction
         aux_dir(ndm,2)     ! the square of the wavefunction

character(len=80) :: flld   ! auxiliary variable

if (nld.eq.0) return
if (nld.gt.nwfsx) &
   call errore('lderiv','nld is too large',1)

ze2=-zed*2.d0

ikrld= 1
do n=1,mesh
   if (r(n).lt.rlderiv) ikrld=n
enddo
write(6,'(5x,''Computing logarithmic derivative in'',f10.5)') &
         (r(ikrld)+r(ikrld+1))*0.5d0

if (ikrld.gt.mesh-1) &
      call errore('lderiv','ikrld is wrong',1)

do is=1,nspin
   do nc=1,nld
      if (rel.lt.2) then
         lam=nc-1
         j=0.d0
      else
         lam=nc/2
         if (mod(nc,2)==0) j=lam-0.5d0
         if (mod(nc,2)==1) j=lam+0.5d0
      endif
      npte= (emaxld-eminld)/deld 
      npte=npte+1
      if (npte.gt.ndm) &
          call errore('lderiv','npte is too large',1)
      do ie=1,npte
         e=eminld+deld*(ie-1.d0)
!
!    integrate outward up to ikrld+1
!
         if (rel.eq.1) then
            call lschps(3,zed,exp(dx),dx,mesh,idum,ikrld+5, &
                       1,lam,e,aux,r,vpot(1,is)) 
         else if (rel.eq.2) then
            rab=r*dx
            call dir_outward(ndm,ikrld+5,lam,j,e,dx,&
                             aux_dir,r,rab,vpot(1,is))
            aux(:)=aux_dir(:,2)
         else
            call intref(lam,e,ikrld+5,dx,r,r2,sqr, &
                    vpot(1,is),ze2,aux)
         endif
!
!    compute the logarithmic derivative and save in dlchi
!            
         dlchi(ie,nc)=compute_log(aux(ikrld-3),r(ikrld),dx)
      enddo
   enddo
   if (file_logderae .ne. ' ') then
       flld=file_logderae
       if (is.eq.2) flld=trim(file_logderae)//'.01'
       open(unit=25,file=flld, status='unknown', iostat=ios, &
                        err=300 )
300    call errore('lderiv','opening file '//file_core,abs(ios))
       do ie=1,npte
          e= eminld+deld*(ie-1)
          write(25,'(10f14.6)') e, (dlchi(ie,nc),nc=1,nld)
       enddo
       close(unit=25)
   endif
enddo
 
return
end
