!-----------------------------------------------------------------------
PROGRAM lr_calculate_spectrum
  !---------------------------------------------------------------------
  ! ... calculates the spectrum
  ! ... by solving tridiagonal problem for each value of omega
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : dp
  USE constants,           ONLY : pi,rytoev,evtonm,rytonm
  USE io_files,            ONLY : tmp_dir, prefix,nd_nmbr
  USE global_version,      ONLY : version_number
  USE io_global,           ONLY : stdout,ionode, ionode_id
  USE environment,         ONLY: environment_start,environment_end
  USE mp_global,           ONLY : mp_startup,mp_global_end
  USE mp_world,            ONLY : world_comm
  USE mp,                  ONLY : mp_bcast, mp_barrier

  IMPLICIT NONE
  !
  CHARACTER(len=256), EXTERNAL :: trimcheck
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
  INTEGER :: ipol
  INTEGER :: itermax, itermax0, itermax_actual
  INTEGER :: sym_op
  INTEGER :: verbosity
  real(kind=dp) :: start,end,increment
  real(kind=dp) :: omegmax,delta_omeg
  CHARACTER(len=60) :: extrapolation,td
  CHARACTER(len=256) :: outdir, filename,eign_file
  INTEGER :: units
  !
  !General use variables & counters
  !
  INTEGER :: n_ipol
  real(dp), ALLOCATABLE, DIMENSION(:,:) :: beta_store, gamma_store
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: zeta_store
  real(dp) :: norm0(3) 
  INTEGER :: i,j, info, ip, ip2, counter
  real(kind=dp) :: average(3), av_amplitude(3), epsil
  INTEGER :: ios
  real (kind=dp) :: alpha_temp(3),scale,wl
  real (kind=dp) :: f_sum
  COMPLEX(kind=dp) :: omeg_c
  COMPLEX(kind=dp) :: green(3,3), eps(3,3)
  COMPLEX(kind=dp), ALLOCATABLE :: a(:), b(:), c(:), r(:,:)
  LOGICAL :: skip, exst
  real(kind=dp) :: omeg, z1,z2
  real(DP) :: degspin
  !
  !
  !For perceived color analysis
  !
  real(dp),PARAMETER :: vis_start=0.116829041,vis_start_wl=780
  real(dp),PARAMETER :: vis_end=0.239806979,vis_end_wl=380
  real(dp) :: perceived_red=0.0d0,perceived_green=0.0d0,perceived_blue=0.0d0
  real(dp) :: perceived_renorm
  INTEGER  :: perceived_itermax,perceived_iter
  real(dp), ALLOCATABLE :: perceived_intensity(:)
  real(dp), ALLOCATABLE :: perceived_evaluated(:)
  LOGICAL  :: do_perceived
  !
  !Subroutines etc.
  !
  COMPLEX(kind=dp), EXTERNAL :: zdotc
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  !DEBUGGING
  real(kind=dp)  :: test

  !
  !User controlled variable initialisation
  !
  NAMELIST / lr_input / itermax, itermax0, itermax_actual, extrapolation,&
                      & end, increment, start, ipol, outdir, prefix,&
                      & epsil, sym_op, verbosity, units,omeg,omegmax,&
                      & delta_omeg, td,eign_file


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
  td='lanczos'
  eign_file = 'pwscf.eigen'
  !
  !Other initialisation
  !
  f_sum=0.0d0
  !
  !DEBUG
  test=0.0d0
  !

#ifdef __MPI
  CALL mp_startup ( )
IF (ionode) THEN
   WRITE(*,*) "Warning: Only a single cpu will be used!"
ENDIF
#endif

  CALL environment_start ( 'TDDFPT_PP' )

! The code starts here

  ios=0
  IF (ionode) THEN ! No need for parallelization in this code
     CALL input_from_file()
     READ (5, lr_input, iostat = ios)
  ENDIF

  CALL mp_bcast ( ios, ionode_id , world_comm )
  CALL errore ('lr_readin', 'reading lr_input namelist', abs (ios) )

  if(trim(td)=="davidson" .or. trim(td)=='david') then
    if(ionode) call spectrum_david()
    goto  555
  endif

  IF (ionode) THEN
  
    IF (itermax0 < 151 .and. trim(extrapolation)/="no") THEN
       WRITE(*,*) "Itermax0 is less than 150, no extrapolation scheme can be used!"
       extrapolation="no"
    ENDIF

    outdir = trimcheck(outdir)
    tmp_dir = outdir
    IF (ipol < 4) THEN
      n_ipol=1
    ELSE
      n_ipol=3
      ipol = 1
    ENDIF

    ! Polarization symmetry
    IF ( .not. sym_op == 0 ) THEN
    ! CALL errore("tddfpt_pp","Unsupported symmetry operation",1)
    IF (sym_op == 1) THEN
      WRITE(stdout,'(5x,"All polarization axes will be considered to be equal.")')
      n_ipol=3
      ipol=1
    ELSE
     CALL errore("tddfpt_pp","Unsupported symmetry operation",1)
    ENDIF

  ENDIF
  ! Terminator Scheme
  IF (trim(extrapolation)=="no") THEN
     !
     itermax=itermax0
     !
  ENDIF
  !
  !Check the unit system used
  !
  IF (units < 0 .or. units >2) THEN
   CALL errore("tddfpt_pp","Unsupported unit system",1)
  ENDIF
  IF ( units /= 0 .and. verbosity > 4) THEN
    WRITE(stdout,'(5x,"Verbosity this high is not supported when non-default units are used")')
    verbosity=4
  ENDIF
  !
  IF (omeg>0 .or. omegmax>0 .or. delta_omeg>0) THEN
    WRITE(stdout,'(5x,"Warning, omeg, omegmax and delta_omeg depreciated, use start,end,increment instead")')
    start=omeg
    end=omegmax
    increment=delta_omeg
    units = 0
  ENDIF
  !
  !Initialisation of coefficients
  !
  ALLOCATE(beta_store(n_ipol,itermax))
  ALLOCATE(gamma_store(n_ipol,itermax))
  ALLOCATE(zeta_store(n_ipol,n_ipol,itermax))
  !
  !beta_store=0.d0
  !gamma_store=0.d0
  !zeta_store=(0.d0,0.d0)
  !
  ALLOCATE(a(itermax))
  ALLOCATE(b(itermax-1))
  ALLOCATE(c(itermax-1))
  ALLOCATE(r(n_ipol,itermax))
  !
  a(:) = (0.0d0,0.0d0)
  b(:) = (0.0d0,0.0d0)
  c(:) = (0.0d0,0.0d0)
  r(:,:) = (0.0d0,0.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL read_b_g_z_file()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL extrapolate()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !
  !
  !  Spectrum calculation
  !
  !
  WRITE (stdout,'(/5x,"Data ready, starting to calculate observables")')
  WRITE (stdout,'(/5x,"Broadening = ",f15.8," Ry")') epsil
  filename = trim(prefix) // ".plot"
  WRITE (stdout,'(/5x,"Output file name: ",A20)') filename
  WRITE(stdout,'(/,5x,"chi_i_j: dipole polarizability tensor in units of e^2*a_0^2/energy")')
  IF (n_ipol == 3) THEN
     WRITE(stdout,'(/,5x,"S: oscillator strength in units of 1/energy")')
     WRITE(stdout,'(/,5x,"S(\hbar \omega) = 2m/( 3 \pi e^2 \hbar) \omega sum_j chi_j_j")')
     WRITE(stdout,'(/,5x,"S(\hbar \omega) satisfies the f-sum rule: \int_0^\infty dE S(E) = N_el ")')
  ELSE
   WRITE (stdout,'(/,5x,"Insufficent info for S")')
  ENDIF
  IF (units == 0) THEN
   WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega Energy unit is (Ry)")')
  ELSEIF (units == 1) THEN
   WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega Energy unit is (eV)")')
  ELSEIF (units == 2) THEN
   WRITE (stdout,'(/,5x,"Functions are reported in (nm), Energy unit is (eV) ")')
  ENDIF
  !
  ! The static polarizability
  !
  IF (verbosity>0) THEN
  WRITE (stdout,'(/,5x,"Static dipole polarizability Tensor:")')
  CALL calc_chi(0.0d0,epsil,green(:,:))
    DO ip=1,n_ipol
        !
        DO ip2=1,n_ipol
           !
           IF(n_ipol == 3) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e21.15," + i",e21.15)') &
                  ip2, ip, dble(green(ip,ip2)), aimag(green(ip,ip2))
           !
           IF(n_ipol == 1) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e21.15," + i",e21.15)') &
                   ipol, ipol, dble(green(ip,ip2)), aimag(green(ip,ip2))
           !
           ENDDO
          !
      ENDDO
  ENDIF

!!!! The output file:
  OPEN(17,file=filename,status="unknown")



!  The perceived color analysis uses the perception fit from the following program:
!      RGB VALUES FOR VISIBLE WAVELENGTHS   by Dan Bruton (astro@tamu.edu)
!
!
! Lets see if the environment is suitable for perceived color analysis
!

  do_perceived=.false.

  IF (verbosity > 2) THEN
   IF (units == 0 .and. start<vis_start .and. end>vis_end .and. n_ipol == 3) THEN
    WRITE (stdout,'(/,5x,"Will attempt to calculate perceived color")')
    do_perceived=.true.
    perceived_iter=1
    perceived_itermax=int((vis_end-vis_start)/increment)
    ALLOCATE(perceived_intensity(perceived_itermax))
    ALLOCATE(perceived_evaluated(perceived_itermax))
    perceived_intensity(:)=0.0d0
    perceived_evaluated(:)=-1.0d0
    perceived_renorm=-9999999.0d0
   ELSEIF (units == 2 .and. start<vis_end_wl .and. end>vis_start_wl .and. n_ipol == 3) THEN
    WRITE (stdout,'(/,5x,"Will attempt to calculate perceived color")')
    do_perceived=.true.
    perceived_iter=1
    perceived_itermax=int((vis_start_wl-vis_end_wl)/increment)
    ALLOCATE(perceived_intensity(perceived_itermax))
    ALLOCATE(perceived_evaluated(perceived_itermax))
    perceived_intensity(:)=0.0d0
    perceived_evaluated(:)=-1.0d0
    perceived_renorm=-9999999.0d0
   ELSEIF (verbosity>2) THEN
    WRITE (stdout,'(/,5x,"Will not calculate perceived color")')
    do_perceived=.false.
   ENDIF
  ENDIF
 !
 ! Header of the output plot file
 !
  IF (units == 0) THEN
   WRITE (17,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (Ry) Re(chi) (e^2*a_0^2/Ry) Im(chi) (e^2*a_0^2/Ry) ")')
  ELSEIF (units == 1) THEN
   WRITE (17,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (eV) Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV) ")')
  ELSEIF (units == 2) THEN
   WRITE (17,'("#Chi is reported as CHI_(i)_(j) wavelength (nm) Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV) ")')
  ENDIF
  IF (n_ipol == 3) THEN
    WRITE(17,'("# S(E) satisfies the sum rule ")' )
  ENDIF

  !
  !   Start the omega loop
  !
  !Units conversion and omega history
  omega(1)=omega(2)
  omega(2)=omega(3)
  IF (units == 0) THEN
   omega(3)=start
  ELSEIF (units == 1) THEN
   omega(3)=start/rytoev
  ELSEIF (units == 2) THEN
   omega(3)=rytonm/start
  ENDIF

  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIRST STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (verbosity > 0 .and. n_ipol == 3) THEN ! In order to gain speed, I perform first term separately
    !
    CALL calc_chi(omega(3),epsil,green(:,:))
    IF (units == 1 .or. units == 2) THEN
     green(:,:)=green(:,:)/rytoev
    ENDIF

    DO ip=1,n_ipol
        !
        DO ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           ENDDO
          !
      ENDDO

    alpha_temp(3)= omega(3)*aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
    IF (units == 1 .or. units == 2) THEN
     alpha_temp(3)=alpha_temp(3)*rytoev
    ENDIF
    !alpha is ready
    WRITE(17,'(5x,"S(E)=",2x,2(e21.15,2x))') &
            start, alpha_temp(3)
    f_sum=0.3333333333333333d0*increment*alpha_temp(3)
    start=start+increment
 ENDIF
!!!!!!!!!!!!!!!!!!OMEGA LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO WHILE(start<end)
   !Units conversion and omega history
    omega(1)=omega(2)
    omega(2)=omega(3)
    IF (units == 0) THEN
      omega(3)=start
    ELSEIF (units == 1) THEN
     omega(3)=start/rytoev
    ELSEIF (units == 2) THEN
     omega(3)=rytonm/start
    ENDIF
     !
     CALL calc_chi(omega(3),epsil,green(:,:))
    IF (units == 1 .or. units == 2) THEN
     green(:,:)=green(:,:)/rytoev
    ENDIF
     !
     DO ip=1,n_ipol
        !
        DO ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
           IF(n_ipol == 3) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
           IF(n_ipol == 1) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ipol, ipol, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           ENDDO
          !
      ENDDO
      IF (n_ipol==3) THEN
        !
        !These are the absorbtion coefficient
        !
        alpha_temp(1)=alpha_temp(2)
        alpha_temp(2)=alpha_temp(3)
        alpha_temp(3)= omega(3)*aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
        IF (units == 1 .or. units == 2) THEN
         alpha_temp(3)=alpha_temp(3)*rytoev
        ENDIF
        !alpha is ready
        WRITE(17,'(5x,"S(E)=",2x,2(e21.15,2x))') &
            start, alpha_temp(3)
        !
        IF (verbosity > 0 ) THEN
         IF ( is_peak(omega(3),alpha_temp(3))) &
            WRITE(stdout,'(5x,"Possible peak at ",F15.8," Ry; Intensity=",E11.2)') omega(1),alpha_temp(1)
         f_sum=f_sum+integrator(increment,alpha_temp(3))
         IF ( omega(3)<vis_end .and. omega(3)>vis_start .and. do_perceived ) THEN
             perceived_intensity(perceived_iter)=alpha_temp(3)
             perceived_evaluated(perceived_iter)=omega(3)
             perceived_iter=perceived_iter+1
             IF (alpha_temp(3) > perceived_renorm) perceived_renorm=alpha_temp(3) !Renormalization to 1
         ENDIF
        ENDIF
     ENDIF
     !
     start=start+increment
     !
  ENDDO
  !
! In order to gain speed, I perform last term seperately
 IF (verbosity > 0 .and. n_ipol == 3) THEN
  !Units conversion
    IF (units == 0) THEN
      omega(3)=start
    ELSEIF (units == 1) THEN
     omega(3)=start/rytoev
    ELSEIF (units == 2) THEN
     omega(3)=rytonm/start
    ENDIF

    CALL calc_chi(omega(3),epsil,green(:,:))
    IF (units == 1 .or. units == 2) THEN
     green(:,:)=green(:,:)/rytoev
    ENDIF
    DO ip=1,n_ipol
        !
        DO ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
!              write(*,'(5x,"eps_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
!                  ip2, ip, ry*omeg, dble(eps), aimag(eps)
              !
           ENDDO
          !
      ENDDO
    alpha_temp(3)= omega(3)*aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
    IF (units == 1 .or. units == 2) THEN
     alpha_temp(3)=alpha_temp(3)*rytoev
    ENDIF
    !alpha is ready
    WRITE(17,'(5x,"S(E)=",2x,2(e21.15,2x))') &
            start, alpha_temp(3)

    f_sum=f_sum+0.3333333333333333d0*increment*alpha_temp(3)
  ENDIF

CLOSE(17)
  !
  IF ( n_ipol==3 .and. verbosity >4 )  THEN
   !S(w)=2m_e/(pi e^2 hbar)
   WRITE(stdout,'(5x,"Integral of absorbtion coefficient ",F15.8)') f_sum
  ENDIF
! The perceived color analysis!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (allocated(perceived_intensity)) THEN
   WRITE(stdout,'(5x,"Perceived color analysis is experimental")')
    perceived_intensity(:)=perceived_intensity(:)/perceived_renorm
    perceived_intensity(:)=1.0d0-perceived_intensity(:) !inverse spectrum
    filename = trim(prefix) // "-spectra.ppm"
    OPEN(UNIT=20,FILE=filename,STATUS='UNKNOWN')
    WRITE(20, '(A2)') 'P6'
    WRITE(20, '(I0,'' '',I0)') perceived_itermax, int(perceived_itermax/8)
    WRITE(20, '(A)') '255'
    DO j=1, int(perceived_itermax/8)
       DO i=1,perceived_itermax
        IF (perceived_evaluated(i)<0.0d0) THEN
          WRITE(20, '(3A1)', advance='no') achar(0), achar(0), achar(0)
         ELSE
          wl=91.1266519/perceived_evaluated(i) !hc/hbar.omega=lambda (hbar.omega in rydberg units)
          !
          !WARNING alpha_temp duty change: now contains R G and B
          !
          CALL wl_to_color(wl,alpha_temp(1),alpha_temp(2),alpha_temp(3))
          !Now the intensities
          !First the degradation toward the end
          IF (wl >700) THEN
            scale=.3+.7* (780.-wl)/(780.-700.)
          ELSEIF (wl<420.) THEN
            scale=.3+.7*(wl-380.)/(420.-380.)
         ELSE
            scale=1.
         ENDIF
         alpha_temp(:)=scale*alpha_temp(:)
         !Then the data from absorbtion spectrum
         alpha_temp(:)=perceived_intensity(i)*alpha_temp(:)
         !The perceived color can also be calculated here
         IF (j==1) THEN
          perceived_red=perceived_red+alpha_temp(1)
          perceived_green=perceived_green+alpha_temp(2)
          perceived_blue=perceived_blue+alpha_temp(3)
         ENDIF
         IF (alpha_temp(1)>1.0d0) PRINT *,alpha_temp(1)
          WRITE(20, '(3A1)', advance='no') achar(int(255*alpha_temp(1))), &
                                           achar(int(255*alpha_temp(2))), &
                                           achar(int(255*alpha_temp(3)))
        ENDIF
     ENDDO
    ENDDO
    CLOSE(20)
    !Now lets write a file with perceived color
    perceived_red=perceived_red/(1.0d0*perceived_itermax)
    perceived_green=perceived_green/(1.0d0*perceived_itermax)
    perceived_blue=perceived_blue/(1.0d0*perceived_itermax)
    WRITE(stdout,'(5x,"Perceived R G B ",3(F15.8,1X))') perceived_red,perceived_green,perceived_blue
    filename = trim(prefix) // "-perceived.ppm"
    OPEN(UNIT=20,FILE=filename,STATUS='UNKNOWN')
    WRITE(20, '(A2)') 'P6'
    WRITE(20, '(I0,'' '',I0)') 180,180
    WRITE(20, '(A)') '255'
    DO j=1, 180
     DO i=1, 180
          WRITE(20, '(3A1)', advance='no') achar(int(255*perceived_red)), &
                                            achar(int(255*perceived_green)), &
                                            achar(int(255*perceived_blue))

     ENDDO
    ENDDO
    CLOSE(20)
   ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (verbosity > 9) THEN
   start=0.0
   f_sum=0.0
   DO WHILE (start<end)
     f_sum=f_sum+integrator(increment,start)
     start=start+increment
   ENDDO
   f_sum=f_sum+0.3333333333333333d0*increment*start
   WRITE(stdout,'(5x,"Integral test:",F15.8,"actual: ",F15.8:)') f_sum,0.5*start*start
  ENDIF


  !
  !Deallocations
  !
  IF (allocated(perceived_intensity)) DEALLOCATE(perceived_intensity)
  IF (allocated(perceived_evaluated)) DEALLOCATE(perceived_evaluated)
  !
  IF (allocated(beta_store)) DEALLOCATE(beta_store)
  IF (allocated(gamma_store)) DEALLOCATE(gamma_store)
  IF (allocated(zeta_store)) DEALLOCATE(zeta_store)
  !
  DEALLOCATE(a)
  DEALLOCATE(b)
  DEALLOCATE(c)
  DEALLOCATE(r)
  !
  CALL environment_end( 'TDDFPT_PP' )
  !
ENDIF
555 print *, "Calculation is finished."
#ifdef __MPI
  CALL mp_barrier (world_comm)
  CALL mp_global_end ()
#endif

CONTAINS
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
     real(kind=dp),INTENT(in) :: omeg, alpha !x and y
     !Internal
     real(kind=dp),SAVE :: omeg_save = 0.0d0, thm1,h2m1,first_der_save=9.0d99
     real(kind=dp),SAVE :: alpha_save(3) = 0.0d0
     INTEGER, SAVE :: current_iter = 0
     LOGICAL, SAVE :: trigger=.true.
     real(kind=dp) :: first_der, second_der

     is_peak=.false.
     !counter
     !Rotate the variables
     IF (current_iter < 3) THEN
      current_iter = current_iter + 1
      omeg_save = omeg
      alpha_save(current_iter) = alpha
      RETURN
     ELSE
      IF (current_iter == 3) THEN
       current_iter = current_iter + 1
       thm1=(omeg-omeg_save)
       h2m1=1.0d0/(thm1*thm1) !for second derivative
       thm1=0.5d0/thm1        !for first derivative
       !thm1=0.083333333333333d0/thm1        !for first derivative
      ENDIF
      !alpha_save(1)=alpha_save(2) !t-2h
      !alpha_save(2)=alpha_save(3) !t-h
      !alpha_save(3)=alpha_save(4) !t
      !alpha_save(4)=alpha_save(5) !t+h
      !alpha_save(5)=alpha         !t+2h
      alpha_save(1)=alpha_save(2)  !t-h
      alpha_save(2)=alpha_save(3)  !t
      alpha_save(3)=alpha          !t+h
     ENDIF
     !The derivatives
     first_der = (alpha_save(3)-alpha_save(1))*thm1
     second_der = (alpha_save(3)-2.0d0*alpha_save(2)+alpha_save(1))*h2m1 ! second derivative corresponds to t, 3 steps before
     !first_der = (-alpha_save(5)+8.0d0*(alpha_save(4)-alpha_save(2))+alpha_save(1))*thm1 !first derivative corresponds to t, 3 steps before
     !second_der = (alpha_save(4)-2.0d0*alpha_save(3)+alpha_save(2))*h2m1 ! second derivative corresponds to t, 3 steps before
     !Decide
     !print *,"w",omeg-0.25d0/thm1,"f=",abs(first_der),"s=",second_der
     !print *,"w",omeg-0.5d0/thm1,"f=",abs(first_der),"s=",second_der
     !if (abs(first_der) < 1.0d-8 .and. second_der < 0 ) is_peak=.true.
     IF (second_der < 0) THEN
      IF (trigger) THEN
        IF (abs(first_der) <abs(first_der_save)) THEN
         first_der_save = first_der
         RETURN
        ELSE
         is_peak=.true.
         trigger=.false.
         RETURN
        ENDIF
      ENDIF
     ELSE
       first_der_save=9.0d99
       trigger=.true.
     ENDIF
     !
     RETURN
 END FUNCTION is_peak
!------------------------------------------------
 REAL(kind=dp) FUNCTION integrator(dh,alpha)
!this function calculates the integral every three points using Simpson's rule
  IMPLICIT NONE
  !Input and output
  real(kind=dp),INTENT(in) :: dh, alpha !x and y
  !internal
  !integer, save :: current_iter = 1
  LOGICAL,SAVE :: flag=.true.
  !real(kind=dp),save :: omeg_save = 0.0d0,dh=0.0d0, alpha_save(2)
!
  integrator=0.0d0
  !COMPOSITE SIMPSON INTEGRATOR, (precision level ~ float)
  ! \int a b f(x) dx = ~ h/3 (f(a) + \sum_odd-n 2*f(a+n*h) + \sum_even-n 4*f(a+n*h) +f(b))
  IF (flag) THEN !odd steps
   integrator=dh*1.33333333333333333D0*alpha
   flag = .false.
   RETURN
  ENDIF
  IF (.not. flag) THEN !even steps
   integrator=dh*0.66666666666666666D0*alpha
   flag = .true.
   RETURN
  ENDIF
 END FUNCTION integrator
!------------------------------------------------
SUBROUTINE read_b_g_z_file()
!Reads the coefficients from the designated file
IMPLICIT NONE
IF (sym_op == 0) THEN
  DO ip=1,n_ipol
   ! Read the coefficents
    IF (n_ipol==3) filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ip))
    IF (n_ipol==1) filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ipol))
    filename = trim(tmp_dir) // trim(filename)
   !
   INQUIRE (file = filename, exist = exst)
   !
   IF (.not.exst) THEN
      !
      CALL errore("tddfpt_calculate_spectrum", "Error reading file",1)
      !WRITE( *,*) "WARNING: " // trim(filename) // " does not exist"
      !stop
      !
   ENDIF

   !
   OPEN (158, file = filename, form = 'formatted', status = 'old')
   !
   READ(158,*) itermax_actual
   WRITE(stdout,'(/5X,"Reading ",I6," Lanczos steps for direction ",I1)') itermax_actual, ip
   WRITE(stdout,'(5X,I6," steps will be considered")') itermax0
   IF (itermax0 > itermax_actual .or. itermax0 > itermax) THEN
    CALL errore("tddfpt_calculate_spectrum", "Error in Itermax0",1)
   ENDIF
   !
   READ(158,*) norm0(ip)
   !print *, "norm0(", ip,")=",norm0(ip)
   !
   DO i=1,itermax0
      !
      READ(158,*) beta_store(ip,i)
      READ(158,*) gamma_store(ip,i)
      READ(158,*) zeta_store (ip,:,i)
      !
    !  print *, "ip=",ip,"i=",i,"beta_store=",beta_store(ip,i),"gamma_store=",gamma_store(ip,i),"zeta_store=",zeta_store (ip,:,i)
   ENDDO
   !
   CLOSE(158)
   beta_store(ip,itermax0+1:)=0.d0
   gamma_store(ip,itermax0+1:)=0.d0
   zeta_store(ip,:,itermax0+1:)=(0.d0,0.d0)

  ENDDO
 ELSEIF (sym_op==1) THEN
  filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ipol))
  filename = trim(tmp_dir) // trim(filename)
  !
  INQUIRE (file = filename, exist = exst)
  !
  IF (.not.exst) THEN
     !
     CALL errore("tddfpt_calculate_spectrum", "Error reading file",1)
     !WRITE( *,*) "ERROR: " // trim(filename) // " does not exist"
     !stop
     !
  ENDIF
   !
   OPEN (158, file = filename, form = 'formatted', status = 'old')
   !
   READ(158,*) itermax_actual
   WRITE(stdout,'(/5X,"Reading ",I6," Lanczos steps for direction ",I1)') itermax_actual, ipol
   WRITE(stdout,'(5X,I6," steps will be considered")') itermax0
  IF (itermax0 > itermax_actual .or. itermax0 > itermax) THEN
   CALL errore("tddfpt_calculate_spectrum", "Error in Itermax0",1)
  ENDIF

   !
   READ(158,*) norm0(1)
   !print *, "norm0(", ip,")=",norm0(ip)
   !
   norm0(2)=norm0(1)
   norm0(3)=norm0(1)
   DO i=1,itermax0
      !
      READ(158,*) beta_store(1,i)
      beta_store(2,i)=beta_store(1,i)
      beta_store(3,i)=beta_store(1,i)
      READ(158,*) gamma_store(1,i)
      gamma_store(2,i)=gamma_store(1,i)
      gamma_store(3,i)=gamma_store(1,i)
      READ(158,*) zeta_store (1,1,i)
      zeta_store (2,2,i)=zeta_store (1,1,i)
      zeta_store (3,3,i)=zeta_store (1,1,i)
      zeta_store (1,2,i)=(0.0d0,0.0d0)
      zeta_store (1,3,i)=(0.0d0,0.0d0)
      zeta_store (2,1,i)=(0.0d0,0.0d0)
      zeta_store (2,3,i)=(0.0d0,0.0d0)
      zeta_store (3,1,i)=(0.0d0,0.0d0)
      zeta_store (3,2,i)=(0.0d0,0.0d0)
      !
   ENDDO
   !
   CLOSE(158)
   beta_store(:,itermax0+1:)=0.d0
   gamma_store(:,itermax0+1:)=0.d0
   zeta_store(:,:,itermax0+1:)=(0.d0,0.d0)

 ENDIF
END SUBROUTINE read_b_g_z_file
!------------------------------------------------
SUBROUTINE extrapolate()
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
IF (trim(extrapolation)/="no") THEN
  !
  average=0.d0
  av_amplitude=0.d0
  !
  DO ip=1,n_ipol
     !
     WRITE(stdout,'(/5x,"Polarization direction:",I1)') ip
     counter=0
     !
     DO i=151,itermax0
        !
        IF (skip .eqv. .true.) THEN
           skip=.false.
           CYCLE
        ENDIF
        !
        IF (mod(i,2)==1) THEN
           !
           IF ( i/=151 .and. abs( beta_store(ip,i)-average(ip)/counter ) > 2.d0 ) THEN
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.true.
              !
           ELSE
              !
              average(ip)=average(ip)+beta_store(ip,i)
              av_amplitude(ip)=av_amplitude(ip)+beta_store(ip,i)
              counter=counter+1
              !print *, "t1 ipol",ip,"av_amp",av_amplitude(ip)
              !
           ENDIF
           !
        ELSE
           !
           IF ( i/=151 .and. abs( beta_store(ip,i)-average(ip)/counter ) > 2.d0 ) THEN
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.true.
              !
           ELSE
              !
              average(ip)=average(ip)+beta_store(ip,i)
              av_amplitude(ip)=av_amplitude(ip)-beta_store(ip,i)
              counter=counter+1
              !print *, "t2 ipol",ip,"av_amp",av_amplitude(ip)
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
     average(ip)=average(ip)/counter
     av_amplitude(ip)=av_amplitude(ip)/counter
     !print *, "t3 ipol",ip,"av_amp",av_amplitude(ip)
     !
     WRITE(stdout,'(5x,"Average =",3F15.8)') average(ip)
     WRITE(stdout,'(5x,"Average oscillation amplitude =",F15.8)') av_amplitude(ip)
  ENDDO
  !
  IF (trim(extrapolation)=="constant") av_amplitude=0
  !
  !
  DO ip=1,n_ipol
     !
     DO i=itermax0,itermax
        !
        IF (mod(i,2)==1) THEN
           !
           beta_store(ip,i)=average(ip)+av_amplitude(ip)
           gamma_store(ip,i)=average(ip)+av_amplitude(ip)
           !
        ELSE
           !
           beta_store(ip,i)=average(ip)-av_amplitude(ip)
           gamma_store(ip,i)=average(ip)-av_amplitude(ip)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
ENDIF
  !
 IF (verbosity > -1) THEN
    ! Write all the coefficients in a file for detailed post-processing
   DO ip=1,n_ipol
     IF (n_ipol==3) filename = trim(prefix) // ".beta_term." // trim(int_to_char(ip))
     IF (n_ipol==1) filename = trim(prefix) // ".beta_term." // trim(int_to_char(ipol))
     filename = trim(tmp_dir) // trim(filename)
    !
     OPEN (17, file = filename, status = 'unknown')
    !
    DO i=1,itermax
       !
       WRITE(17,'(i5,2x,e21.15)') i,beta_store(ip,i)
       !
    ENDDO
    !
    CLOSE(17)
   ENDDO
 ENDIF
END SUBROUTINE extrapolate
!------------------------------------------------
SUBROUTINE calc_chi(freq,broad,chi)
! Calculates the susceptibility,
IMPLICIT NONE
real(kind=dp), INTENT(in) :: freq
real(kind=dp), INTENT(in) :: broad
COMPLEX(kind=dp), INTENT(out) :: chi(:,:)

     omeg_c = cmplx(freq,broad,dp)
     !
     DO ip=1, n_ipol
        !
        a(:) = omeg_c
        !
        DO i=1,itermax-1
           !
           b(i)=cmplx(-beta_store(ip,i),0.0d0,dp)
           c(i)=cmplx(-gamma_store(ip,i),0.0d0,dp)
           !
        ENDDO
        !
        r(ip,:) =cmplx(0.0d0,0.0d0,dp)
        r(ip,1)=cmplx(1.0d0,0.0d0,dp)
        !
        CALL zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info) !|w_t|=(w-L) |1,0,0,...,0|
        IF(info /= 0) &
         CALL errore("tddfpt_pp", "Unable to solve tridiagonal system",1)

        !p=-div.rho'
        !p= chi . E
        ! Thus
        !chi = - <zeta|w_t>
        ! Notice that brodening has a positive sign, thus the abs. coefficient is Im(tr(chi)) not -Im(Tr(chi)) as usual
        DO ip2=1,n_ipol
              chi(ip,ip2)=zdotc(itermax,zeta_store(ip,ip2,:),1,r(ip,:),1)
              chi(ip,ip2)=chi(ip,ip2)*cmplx(norm0(ip),0.0d0,dp)
              ! The response charge density is defined as 2.*evc0*q, see Eq. (43) in JCP 128, 154105 (2008). The dipole is therefore
              ! given by 2.*degspin* zeta^T * (w-T^itermax)^-1 * e_1. See also Eq. (15) in that paper.
              ! the minus sign accounts for the negative electron charge (perturbation is -e E x, rather than E x)
              chi(ip,ip2)=chi(ip,ip2)*cmplx(-2.d0*degspin, 0.d0, dp)

        ENDDO
    ENDDO
END SUBROUTINE calc_chi
!------------------------------------------------
SUBROUTINE wl_to_color(wavelength,red,green,blue)
! Gives the colour intensity of a given wavelength in terms of red green and blue
IMPLICIT NONE
real(kind=dp), INTENT(in) :: wavelength
real(kind=dp), INTENT(out) :: red,green,blue

 IF ((wavelength>=380.).and.(wavelength<=440.)) THEN
              red = -1.*(wavelength-440.)/(440.-380.)
              green = 0.
              blue = 1.
 ENDIF
 IF ((wavelength>=440.).and.(wavelength<=490.)) THEN
   red = 0.
   green = (wavelength-440.)/(490.-440.)
   blue = 1.
 ENDIF
 IF ((wavelength>=490.).and.(wavelength<=510.)) THEN
   red = 0.
   green = 1.
   blue = -1.*(wavelength-510.)/(510.-490.)
 ENDIF
 IF ((wavelength>=510.).and.(wavelength<=580.)) THEN
   red = (wavelength-510.)/(580.-510.)
   green = 1.
   blue = 0.
 ENDIF
 IF ((wavelength>=580.).and.(wavelength<=645.)) THEN
   red = 1.
   green = -1.*(wavelength-645.)/(645.-580.)
   blue = 0.
 ENDIF
 IF ((wavelength>=645.).and.(wavelength<=780.)) THEN
   red = 1.
   green = 0.
   blue = 0.
 ENDIF
END SUBROUTINE wl_to_color

  subroutine spectrum_david()
  
    implicit none
    real(dp) :: energy,chi(4)
    integer :: ieign, nstep, istep
    real(dp) :: frequency, temp
    real(dp), allocatable :: absorption(:,:)

    open(18,file=trim(eign_file),action='read')

    nstep=(end-start)/increment+1
    allocate(absorption(nstep,5)) ! Column 1: Energy; 2: Toal; 3,4,5: X,Y,Z
    absorption(:,:)=0.0d0

    read(18,*)  ! Jump to the second line
522 read(18,*,END=521)  energy,chi
      frequency=start
      istep=1
      do while( .not. istep .gt. nstep )

        absorption(istep,1)=frequency
        temp=frequency-energy
        temp=epsil/(temp**2+epsil**2)
        absorption(istep,2)=absorption(istep,2)+chi(1)*temp
        absorption(istep,3)=absorption(istep,3)+chi(2)*temp
        absorption(istep,4)=absorption(istep,4)+chi(3)*temp
        absorption(istep,5)=absorption(istep,5)+chi(4)*temp
        istep=istep+1
        frequency=frequency+increment
      enddo
    goto 522

521 close(18)

    filename=trim(prefix)//".plot"
    OPEN(17,file=filename,status="unknown")
    write(17,'("#",2x,"Energy(Ry)",10x,"total",13x,"X",13x,"Y",13x,"Z")')
    write(17,'("#  Broadening is: ",5x,F10.7,5x,"Ry")') epsil
    istep=1
    do while( .not. istep .gt. nstep )
      write(17,'(5E20.8)') absorption(istep,1),absorption(istep,1)*absorption(istep,2),&
                           absorption(istep,1)*absorption(istep,3),absorption(istep,1)*&
                           absorption(istep,4),absorption(istep,1)*absorption(istep,5)
      istep=istep+1
    enddo

     print *, "   The spectrum is in file: ", filename

    close(17)
    return
  end subroutine spectrum_david

!------------------------------------------------

END PROGRAM lr_calculate_spectrum
!-----------------------------------------------------------------------

