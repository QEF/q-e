!-----------------------------------------------------------------------
program lr_calculate_spectrum
  !---------------------------------------------------------------------
  ! ... calculates the spectrum
  ! ... by solving tridiagonal problem for each value of omega
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : dp
  USE constants,            ONLY : pi,rytoev,evtonm,rytonm
  use io_files,            only : tmp_dir, prefix,trimcheck,nd_nmbr
  USE global_version,      ONLY : version_number
  USE io_global,           ONLY : stdout,ionode
  USE environment,           ONLY: environment_start,environment_end
  USE mp_global,        ONLY : mp_startup,mp_global_end,mp_barrier
  
  implicit none
  !
  !Constants
  !
!  integer, parameter :: dp=kind(0.d0)
!  real(dp), parameter :: pi=3.14159265d0
!  real(dp), parameter :: ry=13.6056981d0
!  real(dp), parameter :: ev2nm=1239.84172
!  real(dp), parameter :: ry2nm=91.1266519
  !
  !User controlled variables
  !
  real(dp) :: omega(3)
  integer :: ipol
  integer :: itermax, itermax0, itermax_actual
  integer :: sym_op
  integer :: verbosity
  real(kind=dp) :: start,end,increment
  real(kind=dp) :: omegmax,delta_omeg
  character(len=60) :: extrapolation
  character(len=256) :: outdir, filename
  integer :: units
  !
  !General use variables & counters
  !
  integer :: n_ipol
  real(dp), allocatable, dimension(:,:) :: beta_store, gamma_store
  complex(dp), allocatable, dimension(:,:,:) :: zeta_store
  real(dp) :: norm0(3)
  integer :: i,j, info, ip, ip2, counter
  real(kind=dp) :: average(3), av_amplitude(3), epsil
  integer :: ios
  real (kind=dp) :: alpha_temp(3),scale,wl
  real (kind=dp) :: f_sum 
  complex(kind=dp) :: omeg_c
  complex(kind=dp) :: green(3,3), eps(3,3)
  complex(kind=dp), allocatable :: a(:), b(:), c(:), r(:,:)
  logical :: skip, exst
  real(kind=dp) :: omeg, z1,z2
  real(DP) :: degspin 
  ! 
  !
  !For perceived color analysis
  !
  real(dp),parameter :: vis_start=0.116829041,vis_start_wl=780
  real(dp),parameter :: vis_end=0.239806979,vis_end_wl=380
  real(dp) :: perceived_red=0.0d0,perceived_green=0.0d0,perceived_blue=0.0d0
  real(dp) :: perceived_renorm
  integer  :: perceived_itermax,perceived_iter
  real(dp), allocatable :: perceived_intensity(:)
  real(dp), allocatable :: perceived_evaluated(:)
  logical  :: do_perceived 
  !
  !Subroutines etc.
  !
  complex(kind=dp), external :: zdotc
  character(len=6), external :: int_to_char
  !
  !DEBUGGING
  real(kind=dp)  :: test

  !
  !User controlled variable initialisation
  !
  namelist / lr_input / itermax, itermax0, itermax_actual, extrapolation,&
                      & end, increment, start, ipol, outdir, prefix,&
                      & epsil, sym_op, verbosity, units,omeg,omegmax,delta_omeg


  !
  !Initialization of system variables
  !
  !for the time being, degspin is set by hand
  degspin=2.d0
  !
  prefix = 'pwscf'
  outdir = './'
  itermax=1000
  itermax0=1000
  !itermax_actual=1000
  extrapolation="no"
  end=2.50d0
  increment=0.001d0
  start=0.0d0
  epsil=0.02
  ipol=1
  sym_op=0
  verbosity=0
  units=0
  omeg=-1
  delta_omeg=-1
  omegmax=-1
  !
  !Other initialisation
  !
  f_sum=0.0d0
  !
  !DEBUG
  test=0.0d0
  !

#ifdef __PARA
  CALL mp_startup ( )
if (ionode) then 
   write(*,*) "Warning: Only a single cpu will be used!"
endif
#endif

  CALL environment_start ( 'TDDFPT_PP' )

! The code starts here

if (ionode) then !No need for parallelization in this code
  
  call input_from_file()
  read (5, lr_input, iostat = ios)
  
  if (itermax0 < 151 .and. trim(extrapolation).ne."no") then
   write(*,*) "Itermax0 is less than 150, no extrapolation scheme can be used!"
   extrapolation="no"
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
    call errore("tddfpt_pp","Unsupported symmetry operation",1)
   if (sym_op == 1) then
    write(stdout,'(5x,"All polarization axes will be considered to be equal.")')
    n_ipol=3
    ipol=1
   else
    call errore("tddfpt_pp","Unsupported symmetry operation",1)
   endif
  endif
  ! Terminator Scheme
  if (trim(extrapolation)=="no") then
     !
     itermax0=itermax
     !
  end if
  !
  !Check the unit system used
  !
  if (units < 0 .or. units >2) then 
   call errore("tddfpt_pp","Unsupported unit system",1)
  endif
  if ( units /= 0 .and. verbosity > 4) then
    write(stdout,'(5x,"Verbosity this high is not supported when non-default units are used")')
    verbosity=4
  endif
  !
  if (omeg>0 .or. omegmax>0 .or. delta_omeg>0) then
    write(stdout,'(5x,"Warning, omeg, omegmax and delta_omeg depreciated, use start,end,increment instead")')
    start=omeg
    end=omegmax
    increment=delta_omeg
    units = 0
  endif
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
  call read_b_g_z_file()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call extrapolate()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
  !
  !
  !  Spectrum calculation 
  !
  !
  write (stdout,'(/5x,"Data ready, starting to calculate observables")')
  write (stdout,'(/5x,"Broadening = ",f15.8," Ry")') epsil
  filename = trim(prefix) // ".plot"
  write (stdout,'(/5x,"Output file name: ",A20)') filename
  write(stdout,'(/,5x,"chi_i_j: dipole polarizability tensor in units of e^2*a_0^2/energy")')
  if (n_ipol == 3) then
     write(stdout,'(/,5x,"S: oscillator strength in units of 1/energy")')
     write(stdout,'(/,5x,"S(\hbar \omega) = 2m/( 3 \pi e^2 \hbar) \omega sum_j chi_j_j")')
     write(stdout,'(/,5x,"S(\hbar \omega) satisfies the f-sum rule: \int_0^\infty dE S(E) = N_el ")')
  else
   write (stdout,'(/,5x,"Insufficent info for S")')
  endif
  if (units == 0) then
   write (stdout,'(/,5x,"Functions are reported in \hbar.\omega Energy unit is (Ry)")')
  else if (units == 1) then
   write (stdout,'(/,5x,"Functions are reported in \hbar.\omega Energy unit is (eV)")')
  else if (units == 2) then
   write (stdout,'(/,5x,"Functions are reported in (nm), Energy unit is (eV) ")')
  endif
  !
  ! The static polarizability
  !
  if (verbosity>0) then
  write (stdout,'(/,5x,"Static dipole polarizability Tensor:")')
  call calc_chi(0.0d0,epsil,green(:,:)) 
    do ip=1,n_ipol
        !
        do ip2=1,n_ipol
              !
              write(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e21.15," + i",e21.15)') &
                  ip2, ip, dble(green(ip,ip2)), aimag(green(ip,ip2))
              !
           end do
          !
      end do
  endif

!!!! The output file:
  open(17,file=filename,status="unknown")



!  The perceived color analysis uses the perception fit from the following program:
!      RGB VALUES FOR VISIBLE WAVELENGTHS   by Dan Bruton (astro@tamu.edu)
!
!
! Lets see if the environment is suitable for perceived color analysis
!

  if (verbosity > 2) then
   if (units == 0 .and. start<vis_start .and. end>vis_end .and. n_ipol == 3) then
    write (stdout,'(/,5x,"Will attempt to calculate perceived color")')
    do_perceived=.true.
    perceived_iter=1
    perceived_itermax=int((vis_end-vis_start)/increment)
    allocate(perceived_intensity(perceived_itermax))
    allocate(perceived_evaluated(perceived_itermax))
    perceived_intensity(:)=0.0d0
    perceived_evaluated(:)=-1.0d0
    perceived_renorm=-9999999.0d0
   elseif (units == 2 .and. start<vis_end_wl .and. end>vis_start_wl .and. n_ipol == 3) then
    write (stdout,'(/,5x,"Will attempt to calculate perceived color")')
    do_perceived=.true.
    perceived_iter=1
    perceived_itermax=int((vis_start_wl-vis_end_wl)/increment)
    allocate(perceived_intensity(perceived_itermax))
    allocate(perceived_evaluated(perceived_itermax))
    perceived_intensity(:)=0.0d0
    perceived_evaluated(:)=-1.0d0
    perceived_renorm=-9999999.0d0
   elseif (verbosity>2) then
    write (stdout,'(/,5x,"Will not calculate perceived color")')
    do_perceived=.false.
   endif
  endif
 !
 ! Header of the output plot file
 !
  if (units == 0) then
   write (17,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (Ry) Re(chi) (e^2*a_0^2/Ry) Im(chi) (e^2*a_0^2/Ry) ")')
  else if (units == 1) then
   write (17,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (eV) Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV) ")')
  else if (units == 2) then
   write (17,'("#Chi is reported as CHI_(i)_(j) wavelength (nm) Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV) ")')
  endif
  if (n_ipol == 3) then
    write(17,'("# S(E) satisfies the sum rule ")' )
  endif

  !
  !   Start the omega loop
  !
  !Units conversion and omega history
  omega(1)=omega(2)
  omega(2)=omega(3)
  if (units == 0) then
   omega(3)=start
  else if (units == 1) then
   omega(3)=start/rytoev
  else if (units == 2) then
   omega(3)=rytonm/start
  endif
   
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIRST STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (verbosity > 0 .and. n_ipol == 3) then ! In order to gain speed, I perform first term seperately
    ! 
    call calc_chi(omega(3),epsil,green(:,:))
    if (units == 1 .or. units == 2) then
     green(:,:)=green(:,:)/rytoev
    endif

    do ip=1,n_ipol
        !
        do ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              write(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           end do
          !
      end do

    alpha_temp(3)= omega(3)*aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0) 
    if (units == 1 .or. units == 2) then
     alpha_temp(3)=alpha_temp(3)*rytoev
    endif
    !alpha is ready
    write(17,'(5x,"S(E)=",2x,2(e21.15,2x))') &
            start, alpha_temp(3)
    f_sum=0.3333333333333333d0*increment*alpha_temp(3)
    start=start+increment
  endif
!!!!!!!!!!!!!!!!!!OMEGA LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do while(start<end) 
   !Units conversion and omega history
    omega(1)=omega(2)
    omega(2)=omega(3)
    if (units == 0) then
      omega(3)=start
    else if (units == 1) then
     omega(3)=start/rytoev
    else if (units == 2) then
     omega(3)=rytonm/start
    endif
     !
     call calc_chi(omega(3),epsil,green(:,:))
    if (units == 1 .or. units == 2) then
     green(:,:)=green(:,:)/rytoev
    endif
     !
     do ip=1,n_ipol
        !
        do ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              write(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           end do
          !
      end do
      if (n_ipol==3) then
        !
        !These are the absorbtion coefficient
        !
        alpha_temp(1)=alpha_temp(2)
        alpha_temp(2)=alpha_temp(3)
        alpha_temp(3)= omega(3)*aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0) 
        if (units == 1 .or. units == 2) then
         alpha_temp(3)=alpha_temp(3)*rytoev
        endif
        !alpha is ready
        write(17,'(5x,"S(E)=",2x,2(e21.15,2x))') &
            start, alpha_temp(3)
        !
        if (verbosity > 0 ) then
         if ( is_peak(omega(3),alpha_temp(3))) &
            write(stdout,'(5x,"Possible peak at ",F15.8," Ry; Intensity=",E11.2)') omega(1),alpha_temp(1)
         f_sum=f_sum+integrator(increment,alpha_temp(3))
         if ( omega(3)<vis_end .and. omega(3)>vis_start .and. do_perceived ) then
             perceived_intensity(perceived_iter)=alpha_temp(3)
             perceived_evaluated(perceived_iter)=omega(3)
             perceived_iter=perceived_iter+1
             if (alpha_temp(3) > perceived_renorm) perceived_renorm=alpha_temp(3) !Renormalization to 1
         endif
        endif
     end if
     !
     start=start+increment
     !
  enddo
  ! 
! In order to gain speed, I perform last term seperately
 if (verbosity > 0 .and. n_ipol == 3) then     
  !Units conversion
    if (units == 0) then
      omega(3)=start
    else if (units == 1) then
     omega(3)=start/rytoev
    else if (units == 2) then
     omega(3)=rytonm/start
    endif

    call calc_chi(omega(3),epsil,green(:,:)) 
    if (units == 1 .or. units == 2) then
     green(:,:)=green(:,:)/rytoev
    endif
    do ip=1,n_ipol
        !
        do ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              write(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           end do
          !
      end do
    alpha_temp(3)= omega(3)*aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0) 
    if (units == 1 .or. units == 2) then
     alpha_temp(3)=alpha_temp(3)*rytoev
    endif
    !alpha is ready
    write(17,'(5x,"S(E)=",2x,2(e21.15,2x))') &
            start, alpha_temp(3)

    f_sum=f_sum+0.3333333333333333d0*increment*alpha_temp(3)
  endif

close(17)
  !
  if ( n_ipol==3 .and. verbosity >4 )  then 
   !S(w)=2m_e/(pi e^2 hbar)
   write(stdout,'(5x,"Integral of absorbtion coefficient ",F15.8)') f_sum
  endif
! The perceived color analysis!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (allocated(perceived_intensity)) then
   write(stdout,'(5x,"Perceived color analysis is experimental")')
    perceived_intensity(:)=perceived_intensity(:)/perceived_renorm
    perceived_intensity(:)=1.0d0-perceived_intensity(:) !inverse spectrum
    filename = trim(prefix) // "-spectra.ppm"
    open(UNIT=20,FILE=filename,STATUS='UNKNOWN')
    write(20, '(A2)') 'P6'
    write(20, '(I0,'' '',I0)') perceived_itermax, int(perceived_itermax/8)
    write(20, '(A)') '255'
    do j=1, int(perceived_itermax/8)
       do i=1,perceived_itermax
        if (perceived_evaluated(i)<0.0d0) then
          write(20, '(3A1)', advance='no') achar(0), achar(0), achar(0)
         else
          wl=91.1266519/perceived_evaluated(i) !hc/hbar.omega=lambda (hbar.omega in rydberg units)
          !
          !WARNING alpha_temp duty change: now contains R G and B
          !
          call wl_to_color(wl,alpha_temp(1),alpha_temp(2),alpha_temp(3))
          !Now the intensities
          !First the degradation toward the end
          if (wl >700) then
            scale=.3+.7* (780.-wl)/(780.-700.)
          else if (wl<420.) then
            scale=.3+.7*(wl-380.)/(420.-380.)
         else
            scale=1.
         endif
         alpha_temp(:)=scale*alpha_temp(:)
         !Then the data from absorbtion spectrum
         alpha_temp(:)=perceived_intensity(i)*alpha_temp(:)
         !The perceived color can also be calculated here
         if (j==1) then
          perceived_red=perceived_red+alpha_temp(1)
          perceived_green=perceived_green+alpha_temp(2)
          perceived_blue=perceived_blue+alpha_temp(3)
         endif
         if (alpha_temp(1)>1.0d0) print *,alpha_temp(1)
          write(20, '(3A1)', advance='no') achar(int(255*alpha_temp(1))), &
                                           achar(int(255*alpha_temp(2))), &
                                           achar(int(255*alpha_temp(3)))
        endif
     end do
    end do
    close(20)
    !Now lets write a file with perceived color
    perceived_red=perceived_red/(1.0d0*perceived_itermax)
    perceived_green=perceived_green/(1.0d0*perceived_itermax)
    perceived_blue=perceived_blue/(1.0d0*perceived_itermax)
    write(stdout,'(5x,"Perceived R G B ",3(F15.8,1X))') perceived_red,perceived_green,perceived_blue
    filename = trim(prefix) // "-perceived.ppm"
    open(UNIT=20,FILE=filename,STATUS='UNKNOWN')
    write(20, '(A2)') 'P6'
    write(20, '(I0,'' '',I0)') 180,180
    write(20, '(A)') '255'
    do j=1, 180
     do i=1, 180
          write(20, '(3A1)', advance='no') achar(int(255*perceived_red)), &
                                            achar(int(255*perceived_green)), &
                                            achar(int(255*perceived_blue))

     enddo
    enddo
    close(20)
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (verbosity > 9) then
   start=0.0
   f_sum=0.0
   do while (start<end)
     f_sum=f_sum+integrator(increment,start)
     start=start+increment
   enddo
   f_sum=f_sum+0.3333333333333333d0*increment*start
   write(stdout,'(5x,"Integral test:",F15.8,"actual: ",F15.8:)') f_sum,0.5*start*start
  endif


  !
  !Deallocations
  !
  if (allocated(perceived_intensity)) deallocate(perceived_intensity)
  if (allocated(perceived_evaluated)) deallocate(perceived_evaluated)
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
  CALL environment_end( 'TDDFPT_PP' )
  !
endif
#ifdef __PARA
  CALL mp_barrier ()
  CALL mp_global_end ()
#endif

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
 REAL(kind=dp) FUNCTION integrator(dh,alpha)
!this function calculates the integral every three points using Simpson's rule
  IMPLICIT NONE
  !Input and output
  real(kind=dp),intent(in) :: dh, alpha !x and y
  !internal
  !integer, save :: current_iter = 1
  logical,save :: flag=.true.
  !real(kind=dp),save :: omeg_save = 0.0d0,dh=0.0d0, alpha_save(2)
! 
  integrator=0.0d0
  !COMPOSITE SIMPSON INTEGRATOR, (precision level ~ float)
  ! \int a b f(x) dx = ~ h/3 (f(a) + \sum_odd-n 2*f(a+n*h) + \sum_even-n 4*f(a+n*h) +f(b))
  if (flag) then !odd steps
   integrator=dh*1.33333333333333333D0*alpha
   flag = .false.
   return
  endif
  if (.not. flag) then !even steps
   integrator=dh*0.66666666666666666D0*alpha
   flag = .true.
   return
  endif
 END FUNCTION integrator
!------------------------------------------------
subroutine read_b_g_z_file()
!Reads the coefficients from the designated file
IMPLICIT NONE
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
end subroutine read_b_g_z_file
!------------------------------------------------
subroutine extrapolate()
!
!This subroutine applies the "extrapolation" scheme for extrapolating the reduced matrix
!
IMPLICIT NONE
  !
  ! 
  !  Terminatore
  !
  !
skip=.false.
if (trim(extrapolation).ne."no") then
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
  if (trim(extrapolation)=="constant") av_amplitude=0 
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
 if (verbosity > -1) then
    ! Write all the coefficients in a file for detailed post-processing
   do ip=1,n_ipol
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
 endif
end subroutine extrapolate
!------------------------------------------------
subroutine calc_chi(freq,broad,chi)
! Calculates the susceptibility, 
IMPLICIT NONE
real(kind=dp), intent(in) :: freq
real(kind=dp), intent(in) :: broad
complex(kind=dp), intent(out) :: chi(:,:)

     omeg_c = cmplx(freq,broad,dp)
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
        r(ip,:) =cmplx(0.0d0,0.0d0,dp)
        r(ip,1)=cmplx(1.0d0,0.0d0,dp)
        !
        call zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info) !|w_t|=(w-L) |1,0,0,...,0|
        if(info /= 0) &
         call errore("tddfpt_pp", "Unable to solve tridiagonal system",1)

        !p=-div.rho'
        !p= chi . E
        ! Thus
        !chi = - <zeta|w_t>
        ! Notice that brodening has a positive sign, thus the abs. coefficient is Im(tr(chi)) not -Im(Tr(chi)) as usual
        do ip2=1,n_ipol
              chi(ip,ip2)=ZDOTC(itermax,zeta_store(ip,ip2,:),1,r(ip,:),1)
              chi(ip,ip2)=chi(ip,ip2)*cmplx(norm0(ip),0.0d0,dp) 
              ! The response charge density is defined as 2.*evc0*q, see Eq. (43) in JCP 128, 154105 (2008). The dipole is therefore 
              ! given by 2.*degspin* zeta^T * (w-T^itermax)^-1 * e_1. See also Eq. (15) in that paper.
              ! the minus sign accounts for the negative electron charge (perturbation is -e E x, rather than E x)
              chi(ip,ip2)=chi(ip,ip2)*cmplx(-2.d0*degspin, 0.d0, dp)
 
        enddo
    enddo
end subroutine calc_chi
!------------------------------------------------
subroutine wl_to_color(wavelength,red,green,blue)
! Gives the colour intensity of a given wavelength in terms of red green and blue
IMPLICIT NONE
real(kind=dp), intent(in) :: wavelength
real(kind=dp), intent(out) :: red,green,blue

 if ((wavelength>=380.).and.(wavelength<=440.)) then
              red = -1.*(wavelength-440.)/(440.-380.)
              green = 0.
              blue = 1.
 endif
 if ((wavelength>=440.).and.(wavelength<=490.)) then
   red = 0.
   green = (wavelength-440.)/(490.-440.)
   blue = 1.
 endif
 if ((wavelength>=490.).and.(wavelength<=510.)) then 
   red = 0.
   green = 1.
   blue = -1.*(wavelength-510.)/(510.-490.)
 endif
 if ((wavelength>=510.).and.(wavelength<=580.)) then 
   red = (wavelength-510.)/(580.-510.)
   green = 1.
   blue = 0.
 endif
 if ((wavelength>=580.).and.(wavelength<=645.)) then
   red = 1.
   green = -1.*(wavelength-645.)/(645.-580.)
   blue = 0.
 endif
 if ((wavelength>=645.).and.(wavelength<=780.)) then
   red = 1.
   green = 0.
   blue = 0.
 endif
end subroutine wl_to_color

!------------------------------------------------

end program lr_calculate_spectrum
!-----------------------------------------------------------------------

