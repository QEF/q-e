!-----------------------------------------------------------------------
subroutine lr_calculate_spectrum()
  !---------------------------------------------------------------------
  ! ... calculates the spectrum
  ! ... by solving tridiagonal problem for each value of omega
  !---------------------------------------------------------------------
  !
  use io_global,            only : stdout
  use constants,            only : pi, rytoev
  use kinds,                only : dp
  use cell_base,            only : omega
  use lr_variables,         only : n_ipol, ipol, alpha_store, beta_store, gamma_store, zeta_store,&
                                   norm0, itermax!, broadening
  USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout

  !
  implicit none
  !
  integer :: i, info, ip, ip2
  real(kind=dp) :: omeg, omegmax, delta_omeg,broadening
  complex(kind=dp) :: omeg_c
  complex(kind=dp) :: green, eps
  complex(kind=dp), allocatable :: a(:), b(:), c(:), r(:,:)
  !
  complex(kind=dp), external :: zdotc
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_calculate_spectrum>")')
  endif
  omegmax=2.50d0
  delta_omeg=0.001d0
  omeg=0.0d0
  broadening=0.002d0
  !
  allocate(a(itermax))
  allocate(b(itermax-1))
  allocate(c(itermax-1))
  allocate(r(n_ipol,itermax))
  !
  a(:) = (0.0d0,0.0d0)
  b(:) = (0.0d0,0.0d0)
  c(:) = (0.0d0,0.0d0)
  r(:,:) = (0.0d0,0.0d0)
  !
  write(stdout,'(/,5X,"Spectrum calculation")')
  !
  !   Start the omega loop
  !
  do while(omeg<omegmax)
     !
     omeg_c = cmplx(omeg,broadening,dp)
     !
     do ip=1, n_ipol
        !
        a(:) = omeg_c
        !
        do i=1,itermax-1
           !
           b(i)=cmplx(-beta_store(ip,i),0.0d0,dp)
           c(i)=cmplx(-gamma_store(ip,i),0.0d0,dp)
           !
        end do
        !
        r(ip,:) =(0.0d0,0.0d0)
        r(ip,1)=(1.0d0,0.0d0)
        !
        call zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info) !r=>(w-T^itermax)^-1|1,0,...,0| =(w-T^itermax)^-1 e_1
        if(info /= 0) write(stdout,*) " unable to solve tridiagonal system 1"
        !
     end do
     !
     if ( n_ipol==1 ) then
        ! 
        green=ZDOTC(itermax,zeta_store(1,1,:),1,r(1,:),1) !green=>zeta^T * (w-T^itermax)^-1 * e_1  = u^T * V^itermax * (w-T^itermax)^-1 * e_1
        green=green*cmplx(norm0(1),0.0d0,dp)
        !
        eps=(1.d0,0.d0)-(32.d0*pi/omega)*green
        !
        write(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
            ipol, ipol, rytoev*omeg, dble(green), aimag(green)
        write(stdout,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
            ipol, ipol, rytoev*omeg, dble(eps), aimag(eps)            
        !
     else if ( n_ipol==3 ) then
        !
        do ip=1,n_ipol
           !
           do ip2=1,n_ipol
              !
              green=ZDOTC(itermax,zeta_store(ip,ip2,:),1,r(ip,:),1)
              green=green*cmplx(norm0(ip),0.0d0,dp)
              !
              eps=(1.d0,0.d0)-(32.d0*pi/omega)*green
              !
              write(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, rytoev*omeg, dble(green), aimag(green)
              write(stdout,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, rytoev*omeg, dble(eps), aimag(eps)
              !
           end do
           !
        end do
        !
     end if
     !
     omeg=omeg+delta_omeg
     !
  enddo
  !
  deallocate(a)
  deallocate(b)
  deallocate(c)
  deallocate(r)
  !
  return
  !
end subroutine lr_calculate_spectrum
!-----------------------------------------------------------------------
