!-----------------------------------------------------------------------
program lr_calculate_spectrum
  !---------------------------------------------------------------------
  ! ... calculates the spectrum
  ! ... by solving tridiagonal problem for each value of omega
  !---------------------------------------------------------------------
  !
  use io_files,            only : tmp_dir, prefix,trimcheck,nd_nmbr
  USE global_version,      ONLY : version_number
  USE io_global,           ONLY : stdout,ionode
  USE environment,           ONLY: environment_start
  
  implicit none
  !
  integer, parameter :: dp=kind(0.d0)
  real(dp), parameter :: pi=3.14d0
  real(dp), parameter :: ry=13.6056981d0 
  real(dp) :: omega, q
  integer :: n_ipol, ipol
  real(dp), allocatable, dimension(:,:) :: beta_store, gamma_store
  complex(dp), allocatable, dimension(:,:,:) :: zeta_store
  real(dp) :: norm0(3)
  integer :: itermax, itermax0, itermax_actual
  !
  integer :: i, info, ip, ip2, counter
  integer :: ios
  integer :: sym_op
  real(kind=dp) :: omeg, omegmax, delta_omeg, z1,z2
  real(kind=dp) :: average(3), av_amplitude(3), epsil
  real (kind=dp) :: alpha_temp
  real (kind=dp) :: f_sum 
  complex(kind=dp) :: omeg_c
  complex(kind=dp) :: green(3,3), eps(3,3)
  complex(kind=dp), allocatable :: a(:), b(:), c(:), r(:,:)
  logical :: parallel, skip, exst
  character(len=60) :: terminator
  
  character(len=256) :: outdir, filename
  
  !
  complex(kind=dp), external :: zdotc
  character(len=6), external :: int_to_char
  !
  !DEBUGGING
  real(kind=dp)  :: test
  !
  !
  !User controlled variable initialisation
  !
  namelist / lr_input / itermax, itermax0, itermax_actual, terminator,&
                      & omegmax, delta_omeg, omeg, parallel, ipol, outdir, prefix,&
                      & epsil, sym_op
  !
  prefix = 'pwscf'
  outdir = './'
  itermax=1000
  itermax0=1000
  !itermax_actual=1000
  terminator="no"
  omegmax=2.50d0
  delta_omeg=0.001d0
  omeg=0.0d0
  epsil=0.02
  parallel=.true.
  ipol=1
  sym_op=0
  !
  !
  f_sum=0.0d0
  !
  !DEBUG
  test=0.0d0
  !
  
  CALL environment_start ( 'TDDFPT_PP' )

! The code starts here

if (ionode) then !No need for parallelization in this code
  
  call input_from_file()
  read (5, lr_input, iostat = ios)
  
  if (itermax0 < 151 .and. trim(terminator).ne."no") then
   write(*,*) "Itermax0 is less than 150, no terminator scheme can be used!"
   terminator="no"
  endif
  
  outdir = trimcheck(outdir)
  tmp_dir = outdir
  if (ipol < 4) then
    n_ipol=1
  else
    n_ipol=3
    ipol = 1
  endif
  
  ! Polarization symmetry
  if ( .not. sym_op == 0 ) then
   if (sym_op == 1) then
    write(stdout,'(5x,"All polarization axes will be considered to be equal.")')
    n_ipol=3
    ipol=1
   else
    write(stdout,'(5x,"Not supported yet")')
   endif
  endif
  ! Terminator Scheme
  if (trim(terminator)=="no") then
     !
     itermax0=itermax
     !
  end if

  !
  !
  !print *,"n_ipol=",n_ipol
  !
  !Initialisation of coefficients
  !
  allocate(beta_store(n_ipol,itermax))
  allocate(gamma_store(n_ipol,itermax))
  allocate(zeta_store(n_ipol,n_ipol,itermax))
  !
  !beta_store=0.d0
  !gamma_store=0.d0
  !zeta_store=(0.d0,0.d0)
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
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (sym_op == 0) then
  do ip=1,n_ipol
   ! Read the coefficents 
    if (n_ipol==3) filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ip))
    if (n_ipol==1) filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ipol))
    filename = trim(tmp_dir) // trim(filename)
   !
   inquire (file = filename, exist = exst)
   !
   if (.not.exst) then
      !
      call errore("tddfpt_calculate_spectrum", "Error reading file",1)
      !WRITE( *,*) "WARNING: " // trim(filename) // " does not exist"
      !stop
      !
   end if
   
   !
   open (158, file = filename, form = 'formatted', status = 'old')
   !
   read(158,*) itermax_actual
   write(stdout,'(/5X,"Reading ",I6," Lanczos steps for direction ",I1)') itermax_actual, ip
   write(stdout,'(5X,I6," steps will be considered")') itermax0
   if (itermax0 > itermax_actual .or. itermax0 > itermax) then
    call errore("tddfpt_calculate_spectrum", "Error in Itermax0",1)
   endif
   !
   read(158,*) norm0(ip)
   !print *, "norm0(", ip,")=",norm0(ip) 
   !
   do i=1,itermax0
      !
      read(158,*) beta_store(ip,i)
      read(158,*) gamma_store(ip,i)
      read(158,*) zeta_store (ip,:,i) 
      !
    !  print *, "ip=",ip,"i=",i,"beta_store=",beta_store(ip,i),"gamma_store=",gamma_store(ip,i),"zeta_store=",zeta_store (ip,:,i)
   end do
   !
   close(158)
   beta_store(ip,itermax0+1:)=0.d0
   gamma_store(ip,itermax0+1:)=0.d0
   zeta_store(ip,:,itermax0+1:)=(0.d0,0.d0)
  
  enddo
 else if (sym_op==1) then
  filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ipol))
  filename = trim(tmp_dir) // trim(filename)
  !
  inquire (file = filename, exist = exst)
  !
  if (.not.exst) then
     !
     call errore("tddfpt_calculate_spectrum", "Error reading file",1)
     !WRITE( *,*) "ERROR: " // trim(filename) // " does not exist"
     !stop
     !
  end if
   !
   open (158, file = filename, form = 'formatted', status = 'old')
   !
   read(158,*) itermax_actual
   write(stdout,'(/5X,"Reading ",I6," Lanczos steps for direction ",I1)') itermax_actual, ip
   write(stdout,'(5X,I6," steps will be considered")') itermax0
  if (itermax0 > itermax_actual .or. itermax0 > itermax) then
   call errore("tddfpt_calculate_spectrum", "Error in Itermax0",1)
  endif

   !
   read(158,*) norm0(1)
   !
   norm0(2)=norm0(1)
   norm0(3)=norm0(1)
   do i=1,itermax0
      !
      read(158,*) beta_store(1,i)
      beta_store(2,i)=beta_store(1,i)
      beta_store(3,i)=beta_store(1,i)
      read(158,*) gamma_store(1,i)
      gamma_store(2,i)=gamma_store(1,i)
      gamma_store(3,i)=gamma_store(1,i)
      read(158,*) zeta_store (1,1,i)
      zeta_store (2,2,i)=zeta_store (1,1,i)
      zeta_store (3,3,i)=zeta_store (1,1,i)
      zeta_store (1,2,i)=(0.0d0,0.0d0)
      zeta_store (1,3,i)=(0.0d0,0.0d0)
      zeta_store (2,1,i)=(0.0d0,0.0d0)
      zeta_store (2,3,i)=(0.0d0,0.0d0)
      zeta_store (3,1,i)=(0.0d0,0.0d0)
      zeta_store (3,2,i)=(0.0d0,0.0d0)
      !
   end do
   !
   close(158)
   beta_store(ip,itermax0+1:)=0.d0
   gamma_store(ip,itermax0+1:)=0.d0
   zeta_store(ip,:,itermax0+1:)=(0.d0,0.d0)

 endif



 
  !
  ! 
  !  Terminatore
  !
  !
skip=.false.
if (trim(terminator).ne."no") then
  !
  average=0.d0
  av_amplitude=0.d0
  !
  do ip=1,n_ipol
     !
     write(stdout,'(/5x,"Polarization direction:",I1)') ip
     counter=0
     !
     do i=151,itermax0
        !
        if (skip .eqv. .true.) then 
           skip=.false.
           cycle
        end if
        !
        if (mod(i,2)==1) then
           !
           if ( i.ne.151 .and. abs( beta_store(ip,i)-average(ip)/counter ) > 2.d0 ) then
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.true.
              !
           else
              !
              average(ip)=average(ip)+beta_store(ip,i)
              av_amplitude(ip)=av_amplitude(ip)+beta_store(ip,i)
              counter=counter+1
              !print *, "t1 ipol",ip,"av_amp",av_amplitude(ip)
              !
           end if
           !
        else
           !
           if ( i.ne.151 .and. abs( beta_store(ip,i)-average(ip)/counter ) > 2.d0 ) then
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.true.
              !
           else
              !
              average(ip)=average(ip)+beta_store(ip,i)
              av_amplitude(ip)=av_amplitude(ip)-beta_store(ip,i)
              counter=counter+1
              !print *, "t2 ipol",ip,"av_amp",av_amplitude(ip)
              !
           end if
           !
        end if
        !
     end do
     !
     average(ip)=average(ip)/counter
     av_amplitude(ip)=av_amplitude(ip)/counter     
     !print *, "t3 ipol",ip,"av_amp",av_amplitude(ip)
     !
     write(stdout,'(5x,"Average =",3F15.8)') average(ip)
     write(stdout,'(5x,"Average oscillation amplitude =",F15.8)') av_amplitude(ip)
  end do
  !
  if (trim(terminator)=="constant") av_amplitude=0 
  !
  !
  do ip=1,n_ipol
     !
     do i=itermax0,itermax
        !
        if (mod(i,2)==1) then
           !
           beta_store(ip,i)=average(ip)+av_amplitude(ip)
           gamma_store(ip,i)=average(ip)+av_amplitude(ip)
           !
        else
           !
           beta_store(ip,i)=average(ip)-av_amplitude(ip)
           gamma_store(ip,i)=average(ip)-av_amplitude(ip)
           !
        end if
        !
     end do
     !
  end do
  !
end if
  !
  
 do ip=1,n_ipol
  ! Read the coefficents
   if (n_ipol==3) filename = trim(prefix) // ".beta_term." // trim(int_to_char(ip))
   if (n_ipol==1) filename = trim(prefix) // ".beta_term." // trim(int_to_char(ipol))
   filename = trim(tmp_dir) // trim(filename)
  !
   open (17, file = filename, status = 'unknown')
  !
  do i=1,itermax
     !
     write(17,'(i5,2x,e21.15)') i,beta_store(ip,i)
     !
  end do 
  !
  close(17)
 enddo
  !
  !
  !  Spectrum calculation 
  !
  !
  write (stdout,'(/5x,"Data ready, starting to calculate observables")')
  filename = trim(prefix) // ".plot"
  write (stdout,'(/5x,"Output file name: ",A20)') filename
  write (stdout,'(5x,"alpha:absorption coefficient")')
  write (stdout,'(5x,"CHI:susceptibility tensor")')
  write (stdout,'(5x,"Energy unit in output file is eV")')
  open(17,file=filename,status="unknown")
  !   Start the omega loop
  !
  do while(omeg<omegmax)
     !
     omeg_c = cmplx(omeg,epsil,dp)
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
        call zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info) !|w_t|=(w-L) |1,0,0,...,0|
        ! This solves 
        if(info /= 0) write(*,*) " unable to solve tridiagonal system 1"
        !
     end do
     !
     if ( n_ipol==1 ) then
        ! 
        green(1,1)=ZDOTC(itermax,zeta_store(1,1,:),1,r(1,:),1) !green=<zeta_store|r>=<d0psi|evc1>|r>
        !green is <r_i|(w-L)^-1|[r_j,evc0]> 
!        green=sum(conjg(zeta_store(1,1,:))*r(1,:))
        green(1,1)=green(1,1)*cmplx(norm0(1),0.0d0,dp) !green/=scaling factor
        !
        !eps(1,1)=(1.d0,0.d0)-(32.d0*pi/omega)*green(1,1) !This
        !
        write(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') & !chi is the susceptibility
            ipol, ipol, ry*omeg, dble(green(1,1)), aimag(green(1,1))
!        write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') & !epsilon is the polarizability
!            ipol, ipol, ry*omeg, dble(eps), aimag(eps)            
        !
     else if ( n_ipol==3 ) then
        !
        do ip=1,n_ipol
           !
           do ip2=1,n_ipol
              !
              green(ip,ip2)=ZDOTC(itermax,zeta_store(ip,ip2,:),1,r(ip,:),1)
              green(ip,ip2)=green(ip,ip2)*cmplx(norm0(ip),0.0d0,dp)
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              write(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, ry*omeg, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           end do
           !
        end do
        !
        !These are the absorbtion coefficient
        !
        ! Dario's interpretation of the absorbtion coefficient breaks the sum rule
        !alpha_temp= -omeg*ry*aimag(green(1,1)+green(2,2)+green(3,3))/3.d0 
        ! 
        alpha_temp= -omeg*aimag(green(1,1)+green(2,2)+green(3,3))/3.d0
        write(17,'(5x,"alpha",2x,3(e21.15,2x))') &
            omeg*ry, alpha_temp
        !
        !if (is_peak(omeg,alpha_temp)) write(stdout,'(5x,"Possible resonance in the vicinity of ",F15.8," Ry")') omeg-4.0d0*delta_omeg
        if (is_peak(omeg,alpha_temp)) write(stdout,'(5x,"Possible resonance in the vicinity of ",F15.8," Ry")') omeg-2.0d0*delta_omeg
        f_sum=f_sum+integrator(omeg,alpha_temp)
     end if
     !
     omeg=omeg+delta_omeg
     !
  enddo
  !
close(17)
  !
  if ( n_ipol==3 )  write(stdout,'(5x,"Integral of absorbtion coefficient =",F15.8)') f_sum
  !omeg=0.0d0
  !do while (omeg<1)
    !omeg=omeg+0.000001d0
    !test=test+integrator(omeg,cos(omeg))
  !enddo
  !print *, "test=",test,"real=",sin(omeg),"difference=",abs(test-sin(omeg))
  if (allocated(beta_store)) deallocate(beta_store)
  if (allocated(gamma_store)) deallocate(gamma_store)
  if (allocated(zeta_store)) deallocate(zeta_store)
  !
  deallocate(a)
  deallocate(b)
  deallocate(c)
  deallocate(r)
  !
endif
contains
 LOGICAL FUNCTION is_peak(omeg,alpha)
 !
 ! A simple algorithm for detecting peaks 
 ! Increments of omega between alpha steps should be constant
 ! omega must increase monothonically
 ! no checks performed!
 ! OBM 2010
 !
     IMPLICIT NONE
     !Input and output
     real(kind=dp),intent(in) :: omeg, alpha !x and y
     !Internal
     real(kind=dp),save :: omeg_save = 0.0d0, thm1,h2m1,first_der_save=9.0d99
     real(kind=dp),save :: alpha_save(3) = 0.0d0
     integer, save :: current_iter = 0
     logical, save :: trigger=.true.
     real(kind=dp) :: first_der, second_der
     
     is_peak=.false.
     !counter 
     !Rotate the variables
     if (current_iter < 3) then
      current_iter = current_iter + 1
      omeg_save = omeg
      alpha_save(current_iter) = alpha
      return
     else
      if (current_iter == 3) then
       current_iter = current_iter + 1
       thm1=(omeg-omeg_save)
       h2m1=1.0d0/(thm1*thm1) !for second derivative
       thm1=0.5d0/thm1        !for first derivative
       !thm1=0.083333333333333d0/thm1        !for first derivative
      endif
      !alpha_save(1)=alpha_save(2) !t-2h
      !alpha_save(2)=alpha_save(3) !t-h
      !alpha_save(3)=alpha_save(4) !t
      !alpha_save(4)=alpha_save(5) !t+h
      !alpha_save(5)=alpha         !t+2h
      alpha_save(1)=alpha_save(2)  !t-h
      alpha_save(2)=alpha_save(3)  !t
      alpha_save(3)=alpha          !t+h
     endif
     !The derivatives
     first_der = (alpha_save(3)-alpha_save(1))*thm1
     second_der = (alpha_save(3)-2.0d0*alpha_save(2)+alpha_save(1))*h2m1 ! second derivative corresponds to t, 3 steps before
     !first_der = (-alpha_save(5)+8.0d0*(alpha_save(4)-alpha_save(2))+alpha_save(1))*thm1 !first derivative corresponds to t, 3 steps before
     !second_der = (alpha_save(4)-2.0d0*alpha_save(3)+alpha_save(2))*h2m1 ! second derivative corresponds to t, 3 steps before
     !Decide
     !print *,"w",omeg-0.25d0/thm1,"f=",abs(first_der),"s=",second_der
     !print *,"w",omeg-0.5d0/thm1,"f=",abs(first_der),"s=",second_der
     !if (abs(first_der) < 1.0d-8 .and. second_der < 0 ) is_peak=.true.
     if (second_der < 0) then
      if (trigger) then
        if (abs(first_der) <abs(first_der_save)) then
         first_der_save = first_der
         return
        else
         is_peak=.true.
         trigger=.false.
         return
        endif 
      endif
     else
       first_der_save=9.0d99
       trigger=.true.
     endif
     !
     return
 END FUNCTION is_peak
!------------------------------------------------
 REAL(kind=dp) FUNCTION integrator(omeg,alpha)
!this function calculates the integral every three points using Simpson's rule
  IMPLICIT NONE
  !Input and output
  real(kind=dp),intent(in) :: omeg, alpha !x and y
  !internal
  integer, save :: current_iter = 1
  real(kind=dp),save :: omeg_save = 0.0d0,dh=0.0, alpha_save(2)
! 
  integrator=0.0d0
     if (current_iter < 3) then
      if (current_iter == 2) dh=0.16666666666666666667D0*(omeg-omeg_save)
      omeg_save = omeg
      alpha_save(current_iter) = alpha
      current_iter = current_iter + 1
      return
     else
      !simpsons rule \int (x-h) (x+h) f(x) dx ~ 1h/6 (f(x-h) + 4f(x) + f(x+h)) 
      integrator = dh*(alpha_save(1) + 4.0d0*alpha_save(2) + alpha)
      alpha_save(1)=alpha_save(2)
      alpha_save(2)=alpha
      current_iter = current_iter + 1
     endif

 END FUNCTION integrator

end program lr_calculate_spectrum
!-----------------------------------------------------------------------

