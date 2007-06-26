!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine lderiv
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation 
  !  computing logarithmic derivatives for Coulomb potential
  !
  !
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
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

  real(DP) ::           &
       aux(ndm),         & ! the square of the wavefunction
       aux_dir(ndm,2),   & ! the square of the wavefunction
       ze2,              & ! the nuclear charge in Ry units
       e,                & ! the eigenvalue
       j                   ! total angular momentum for log_der

  real(DP), external ::           &
      compute_log 

  real(DP), allocatable ::        &
       dlchi(:, :)         ! the logarithmic derivative

  character(len=256) :: flld


  if (nld == 0 .or. file_logder == ' ') return
  if (nld > nwfsx) call errore('lderiv','nld is too large',1)

  ze2=-zed*2.0_dp

  do n=1,mesh
     if (r(n) > rlderiv) go to 10
  enddo
  call errore('lderiv','wrong rlderiv?',1)
10 ikrld = n-1
  write(stdout,'(5x,''Computing logarithmic derivative in'',f10.5)') &
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
              call dir_outward(ndm,ikrld+5,lam,j,e,dx,&
                   aux_dir,r,rab,vpot(1,is))
              aux(:)=aux_dir(:,1)
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
      
     if (nspin == 2 .and. is == 1) then
        flld = trim(file_logder)//'up'
     else if (nspin == 2 .and. is == 2) then
        flld = trim(file_logder)//'dw'
     else
        flld = trim(file_logder)
     end if
     if (ionode) &
        open(unit=25, file=flld, status='unknown', iostat=ios, err=300 )
300  call mp_bcast(ios,ionode_id)
     call errore('lderivps','opening file '//flld, abs(ios))
     if (ionode) then
        do ie=1,npte
           e= eminld+deld*(ie-1)
           write(25,'(10f14.6)') e, (dlchi(ie,nc),nc=1,nld)
        enddo
        close(unit=25)
     endif
  enddo
  deallocate (dlchi)
  return
end subroutine lderiv
