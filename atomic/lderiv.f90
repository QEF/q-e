!
!---------------------------------------------------------------
subroutine lderiv
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation 
  !  computing logarithmic derivatives for Coulomb potential
  !
  !
  use ld1inc

  implicit none

  integer ::        &
       lam,    &   ! the angular momentum
       ikrld,  &   ! index of matching radius
       nc,     &   ! counter on logarithmic derivatives
       idum,   &   ! integer variable for lschps
       is,     &   ! counter on spin
       ierr,   &   ! used for allocation control
       ios,    &   ! used for I/O control
       n,ie        ! generic counter 

  real(kind=dp) ::           &
       rab(ndm),         & ! derivative of the radial mesh
       aux(ndm),         & ! the square of the wavefunction
       aux_dir(ndm,2),   & ! the square of the wavefunction
       ze2,              & ! the nuclear charge in Ry units
       e,                & ! the eigenvalue
       j                   ! total angular momentum for log_der

  real(kind=dp), external ::           &
      compute_log 

  real(kind=dp), allocatable ::        &
       dlchi(:, :)         ! the logarithmic derivative

  character(len=256) :: flld   ! auxiliary variable

  if (nld == 0 .or. file_logderae == ' ') return
  if (nld > nwfsx) call errore('lderiv','nld is too large',1)

  ze2=-zed*2.0_dp

  do n=1,mesh
     if (r(n) > rlderiv) go to 10
  enddo
  call errore('lderiv','wrong rlderiv?',1)
10 ikrld = n-1
  write(6,'(5x,''Computing logarithmic derivative in'',f10.5)') &
       (r(ikrld)+r(ikrld+1))*0.5_dp

  npte= (emaxld-eminld)/deld + 1
  allocate ( dlchi(npte, nld) )

  do is=1,nspin
     do nc=1,nld
        if (rel < 2) then
           lam=nc-1
           j=0.0_dp
        else
           lam=nc/2
           if (mod(nc,2)==0) j=lam-0.5_dp
           if (mod(nc,2)==1) j=lam+0.5_dp
        endif
        do ie=1,npte
           e=eminld+deld*(ie-1.0_dp)
           !
           !    integrate outward up to ikrld+1
           !
           if (rel == 1) then
              call lschps(3,zed,exp(dx),dx,mesh,idum,ikrld+5, &
                   1,lam,e,aux,r,vpot(1,is)) 
           else if (rel == 2) then
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
           dlchi(ie, nc) = compute_log(aux(ikrld-3),r(ikrld),dx)
        enddo
     enddo
     flld=file_logderae
     if (is == 2) flld=trim(file_logderae)//'.01'
     open(unit=25,file=flld, status='unknown', iostat=ios, &
          err=300 )
300  call errore('lderiv','opening file '//flld,abs(ios))
     do ie=1,npte
        e= eminld+deld*(ie-1)
        write(25,'(10f14.6)') e, (dlchi(ie,nc),nc=1,nld)
     enddo
     close(unit=25)

  enddo
  deallocate (dlchi)
  return
end subroutine lderiv
