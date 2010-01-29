!-----------------------------------------------------------------------
program lr_calculate_spectrum
  !---------------------------------------------------------------------
  ! ... calculates the spectrum
  ! ... by solving tridiagonal problem for each value of omega
  !---------------------------------------------------------------------
  !
  use io_files,            only : tmp_dir, prefix,trimcheck,nd_nmbr
  USE global_version,      ONLY : version_number
  USE io_global,           ONLY : stdout
  USE environment,           ONLY: environment_start
  
  implicit none
  !
  integer, parameter :: dp=kind(0.d0)
  real(dp), parameter :: pi=3.14d0
  real(dp), parameter :: ry=13.6d0 
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
  !CALL startup (nd_nmbr, "TDDFPT_PP", "1.0")
  CALL environment_start ( 'TDDFPT_PP' )


  
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
    write(*,*) "All polarization axes will be considered to be equal."
    n_ipol=3
    ipol=1
   else
    write(*,*) "Not supported yet"
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
      WRITE( *,*) "WARNING: " // trim(filename) // " does not exist"
      stop
      !
   end if
   
   !
   open (158, file = filename, form = 'formatted', status = 'old')
   !
   read(158,*) itermax_actual
   print *, "Reading", itermax_actual, " Lanczos steps for direction", ip
   print *, itermax0, " steps will be considered"
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
      read(158,*) zeta_store (ip,:,i) !warning, in the old part, the imaginary part was set to zero forcibly. inquire (see Dario's thesis Improving the numerical efficency part)
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
     WRITE( *,*) "ERROR: " // trim(filename) // " does not exist"
     stop
     !
  end if
   !
   open (158, file = filename, form = 'formatted', status = 'old')
   !
   read(158,*) itermax_actual
  print *, "Reading", itermax_actual, " Lanczos steps in direction", ipol
  print *, itermax0, " steps will be considered"
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
  end do
  !
  write(*,*) "average =",average
  if (trim(terminator)=="constant") av_amplitude=0 
  write(*,*) "average oscillation amplitude =",av_amplitude
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
  write (*,*) "Output energy unit is in eV" 
  filename = trim(prefix) // ".plot"
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
        !green is apparently <r_i|(w-L)^-1|[r_j,evc0]> 
        !(does this mean alpha_ij in the thesis is actually chi?)
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
        write(17,'(5x,"alpha",2x,3(e21.15,2x))') &
            ry*omeg, -ry*omeg*aimag(green(1,1)+green(2,2)+green(3,3))/3.d0
        !
     end if
     !
     omeg=omeg+delta_omeg
     !
  enddo
  !
close(17)
  !
  if (allocated(beta_store)) deallocate(beta_store)
  if (allocated(gamma_store)) deallocate(gamma_store)
  if (allocated(zeta_store)) deallocate(zeta_store)
  !
  deallocate(a)
  deallocate(b)
  deallocate(c)
  deallocate(r)
  !
end program lr_calculate_spectrum
!-----------------------------------------------------------------------

