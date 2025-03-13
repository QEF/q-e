!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM lr_calculate_spectrum
  !---------------------------------------------------------------------
  !
  ! Calculates the spectrum by solving tridiagonal problem for each value 
  ! of the frequency omega
  !
  ! Modified by Osman Baris Malcioglu (2008)
  ! Modified by Xiaochuan Ge (2013)
  ! Modified by Iurii Timrov (2015)  
  ! Modified by Tommaso Gorni (2022)  
  !
  USE kinds,               ONLY : dp
  USE constants,           ONLY : pi,rytoev,evtonm,rytonm
  USE io_files,            ONLY : tmp_dir, prefix
  USE io_global,           ONLY : stdout,ionode, ionode_id
  USE environment,         ONLY : environment_start,environment_end
  USE mp_global,           ONLY : mp_startup,mp_global_end, my_image_id
  USE mp_world,            ONLY : world_comm
  USE mp,                  ONLY : mp_bcast, mp_barrier

  IMPLICIT NONE
  !
  CHARACTER(len=256), EXTERNAL :: trimcheck
  !
  ! User controlled variables
  !
  REAL(dp) :: omega(3)
  INTEGER :: ipol
  INTEGER :: itermax, itermax0, itermax_actual
  INTEGER :: sym_op
  INTEGER :: verbosity
  REAL(kind=dp) :: start,end,increment
  REAL(kind=dp) :: omegmax,delta_omeg
  CHARACTER(len=60) :: extrapolation, td
  CHARACTER(len=256) :: outdir, filename, filename1, &
                      & eign_file, tmp_dir_lr
  INTEGER :: units
  LOGICAL :: eels
  LOGICAL :: magnons
  !
  ! General use variables & counters
  !
  INTEGER :: n_ipol, i,j, info, ip, ip2, counter, ios, n_op
  REAL(dp) :: norm0(3), factor_eels, volume, &
            & alat, q1, q2, q3, modulus_q, nelec, &
            & average(3), av_amplitude(3), epsil, &
            & alpha_temp(3), scale, wl, f_sum, &
            & omeg, z1,z2, degspin, integration_function, start_save
  COMPLEX(kind=dp) :: omeg_c
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: beta_store, gamma_store
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: alpha_store
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: gamma_magnons_store
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: zeta_store
  COMPLEX(kind=dp) :: green(3,3), &  ! susceptibility chi
                      eps(3,3),   &  ! dielectric function
                      epsm1(3,3)     ! inverse dielectric function
  COMPLEX(kind=dp), ALLOCATABLE :: a(:), b(:), c(:), r(:,:)
  LOGICAL :: skip, exst
  !
  COMPLEX(dp) :: a_av(3), c_av(3), a_ampli(3), c_ampli(3)
  !
  ! For perceived color analysis
  !
  REAL(dp),PARAMETER :: vis_start    = 0.116829041, &
                      & vis_start_wl = 780
  REAL(dp),PARAMETER :: vis_end      = 0.239806979, &
                      & vis_end_wl   = 380
  REAL(dp) :: perceived_red   = 0.0d0, &
            & perceived_green = 0.0d0, &
            & perceived_blue  = 0.0d0
  REAL(dp) :: perceived_renorm
  INTEGER  :: perceived_itermax,perceived_iter
  REAL(dp), ALLOCATABLE :: perceived_intensity(:), &
                         & perceived_evaluated(:)
  LOGICAL  :: do_perceived
  !
  ! Subroutines etc.
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! User controlled variable initialisation
  !
  NAMELIST / lr_input / itermax, itermax0, itermax_actual, extrapolation,&
                      & end, increment, start, ipol, outdir, prefix,&
                      & epsil, sym_op, verbosity, units, omeg, omegmax,&
                      & delta_omeg, td, eign_file, eels, magnons
  !
  ! Checking for the path to the output directory.
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  !
  ! Initialization of system variables
  !
  ! degspin is read from file in read_b_g_z_file
  ! degspin = 2.d0
  !
  prefix = 'pwscf'
  itermax = 1000
  itermax0 = 1000
  extrapolation = "no"
  end = 2.50d0
  increment = 0.001d0
  start = 0.0d0
  epsil = 0.02
  ipol = 1
  sym_op = 0
  verbosity = 0
  units = 0
  omeg = -1
  delta_omeg = -1
  omegmax = -1
  td = 'lanczos'
  eels = .false.
  f_sum = 0.0d0
  eign_file = 'pwscf.eigen' 
  magnons = .false.
  !
#if defined(__MPI)
  CALL mp_startup ( )
  IF (ionode) THEN
     WRITE(*,*) "Warning: Only a single CPU will be used!"
  ENDIF
#endif
  !
  CALL environment_start ( 'TDDFPT_PP' )
  !
  ! No need for parallelization in this code
  !
  ios = 0
  !
  IF (ionode) THEN
     ! 
     CALL input_from_file()
     READ (5, lr_input, iostat = ios)
     !
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id , world_comm )
  CALL errore ('lr_calculate_spectrum', 'reading lr_input namelist', abs (ios) )
  !
  IF (ionode) THEN
     !
     IF (trim(td)=="davidson" .or. trim(td)=='david' .and. eels) &
           & CALL errore ('lr_calculate_spectrum', 'EELS + Davidson is not supported!', abs (ios) )
     !
     ! Davidson case 
     !
     IF (trim(td)=="davidson" .or. trim(td)=='david') THEN
        !
        CALL spectrum_david()
        GOTO  555
        !
     ENDIF
     !
     ! Lanczos case
     !
     IF (itermax0 < 151 .and. trim(extrapolation)/="no") THEN
        WRITE(*,*) "Itermax0 is less than 150, no extrapolation scheme can be used!"
        extrapolation="no"
     ENDIF
     !
     IF (omeg>0 .or. omegmax>0 .or. delta_omeg>0) THEN
        !
        WRITE(stdout,'(5x,"Warning, omeg, omegmax and delta_omeg depreciated, &
                        &  use start,end,increment instead")')
        ! 
        start = omeg
        end = omegmax
        increment = delta_omeg
        units = 0
        !
     ENDIF
     !
     ! Call the proper driver
     !
     IF (eels) THEN
       !
       CALL compute_eels_spectrum()
       !
     ELSEIF (magnons) THEN
       !
       CALL compute_magnon_spectrum()
       !
     ELSE
       !
       CALL compute_optical_spectrum()
       !
     ENDIF
     !
     CALL environment_end( 'TDDFPT_PP' )
     !
  ENDIF
  !
555 IF (trim(td)=="davidson" .or. trim(td)=='david') &
              & print *, "Calculation is finished."
  !
#if defined(__MPI)
  CALL mp_barrier (world_comm)
  CALL mp_global_end ()
#endif
  !
  STOP
  !
CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE compute_optical_spectrum()
  !-----------------------------------------------------------------------------
  !
  ! Compute the susceptibility in the optical case (dipole-dipole response function)
  ! 
  IMPLICIT NONE
     
     ! Definition of the temporary directory tmp_dir,
     ! where to read the data produced by TDDFPT
     !
     outdir = trimcheck(outdir)
     tmp_dir = outdir
     !
     IF (ipol < 4) THEN
       n_ipol=1
     ELSE
       n_ipol=3
       ipol = 1
     ENDIF
     !
     ! Polarization symmetry
     !
     IF ( .not. sym_op == 0 ) THEN
        !
        IF (sym_op == 1) THEN
           WRITE(stdout,'(5x,"All polarization axes will be considered to be equal.")')
           n_ipol = 3
           ipol = 1
        ELSE
           CALL errore("lr_calculate_spectrum","Unsupported symmetry operation",1)
        ENDIF
        !
     ENDIF
     !
     ! Terminator scheme
     !
     IF (trim(extrapolation)=="no") THEN
        !
        itermax = itermax0
        !
     ENDIF
     !
     ! Check the units (Ry, eV, nm)
     !
     IF (units < 0 .or. units >3) CALL errore("lr_calculate_spectrum","Unsupported unit system",1)
     !
     IF (units==3) CALL errore("lr_calculate_spectrum","meV unit is supported only for magnon= .true.",1)
     !
     IF ( units /= 0 .and. verbosity > 4) THEN
        WRITE(stdout,'(5x,"Such a high verbosity is not supported when &
                        & non-default units are used")')
        verbosity = 4
     ENDIF
     !
     ! Initialisation of coefficients
     !
     ALLOCATE(beta_store(n_ipol,itermax))
     ALLOCATE(gamma_store(n_ipol,itermax))
     ALLOCATE(zeta_store(n_ipol,n_ipol,itermax))
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
     !
     ! Read beta, gamma, and zeta coefficients
     !
     CALL read_b_g_z_file()
     !
     ! Optional: use an extrapolation scheme
     !
     CALL extrapolate()
     !
     !  Spectrum calculation
     !
     WRITE (stdout,'(/5x,"Data ready, starting to calculate observables...")')
     WRITE (stdout,'(/5x,"Broadening = ",f15.8," Ry")') epsil
     !
     filename = trim(prefix) // ".plot_chi.dat"
     !
     WRITE (stdout,'(/5x,"Output file name: ",A)') trim(filename)
     !
     filename1 = trim(prefix) // ".plot_S.dat"
     WRITE (stdout,'(/5x,"Output file name for the oscillator strength S &
                     & ",A)') trim(filename1)
     !
     WRITE(stdout,'(/,5x,"chi_i_j: dipole polarizability tensor &
                                     & in units of e^2*a_0^2/energy")')
     !
     IF (n_ipol == 3) THEN
       WRITE(stdout,'(/,5x,"S: oscillator strength in units of 1/energy")')
       WRITE(stdout,'(/,5x,"S(\hbar \omega) = 2m/( 3 \pi e^2 \hbar) &
                              & \omega sum_j chi_j_j")')
       WRITE(stdout,'(/,5x,"S(\hbar \omega) satisfies the f-sum rule: &
                              & \int_0^\infty dE S(E) = N_el ")')
     ELSE
       WRITE (stdout,'(/,5x,"Insufficent info for S")')
     ENDIF
     !
     ! Units
     !
     IF (units == 0) THEN      ! Ry
        WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega &
                              & Energy unit is (Ry)")')
     ELSEIF (units == 1) THEN  ! eV
        WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega &
                              & Energy unit is (eV)")')
     ELSEIF (units == 2) THEN  ! nm
        WRITE (stdout,'(/,5x,"Functions are reported in (nm), &
                              & Energy unit is (eV) ")')
     ENDIF
     !
     ! The static dipole polarizability / static charge-density susceptibility
     !
     IF (verbosity>0) THEN
        !
        WRITE (stdout,'(/,5x,"Static dipole polarizability tensor:")')
        !
        CALL calc_chi(0.0d0,epsil,green(:,:))
        !
        DO ip=1,n_ipol
          DO ip2=1,n_ipol
             !
             IF (n_ipol == 3) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e21.15," + i",e21.15)') &
                                         & ip2, ip, dble(green(ip,ip2)), aimag(green(ip,ip2))
             !
             IF (n_ipol == 1) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e21.15," + i",e21.15)') &
                                         & ipol, ipol, dble(green(ip,ip2)), aimag(green(ip,ip2))
            !
          ENDDO
        ENDDO
        !
     ENDIF
     !
     ! Open the output file
     !
     OPEN(17,file=filename,status="unknown")
     !
     IF (n_ipol==3) THEN
        !
        OPEN(18,file=filename1,status="unknown")
        !
        IF (units == 0) THEN
           WRITE(18,'("#",10x,"\hbar \omega(Ry)",4x,"Oscillator strength")')
        ELSEIF (units == 1) THEN
           WRITE(18,'("#",10x,"\hbar \omega(eV)",4x,"Oscillator strength")')
        ELSEIF (units == 2) THEN
           WRITE(18,'("#",10x,"wavelength(nm)",4x,"Oscillator strength")')
        ENDIF
        !
     ENDIF
     !
     !-----------------------------------------------------------------------------!
     !                          PERCEIVED COLOR ANALYSIS                           !
     !-----------------------------------------------------------------------------!
     !
     ! The perceived color analysis uses the perception fit from the following program:
     ! RGB VALUES FOR VISIBLE WAVELENGTHS   by Dan Bruton (astro@tamu.edu)
     ! Let's see if the environment is suitable for perceived color analysis
     ! This is needed for optics, not for EELS.
     !
     do_perceived = .false.
     !
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
        WRITE (17,'("#",16x,"\hbar \omega(Ry)",5x,"Re(chi) (e^2*a_0^2/Ry)",x,"Im(chi) (e^2*a_0^2/Ry)")')
     ELSEIF (units == 1) THEN
        WRITE (17,'("#",16x,"\hbar \omega(eV)",5x,"Re(chi) (e^2*a_0^2/eV)",x,"Im(chi) (e^2*a_0^2/eV)")')
     ELSEIF (units == 2) THEN
        WRITE (17,'("#",16x,"wavelength(nm)",5x,"Re(chi) (e^2*a_0^2/eV)",x,"Im(chi) (e^2*a_0^2/eV)")')
     ENDIF
     !
     ! Start a loop on frequency
     !
     ! Units conversion and omega history
     !
     omega(1) = omega(2)
     omega(2) = omega(3)
     !
     IF (units == 0) THEN
        omega(3) = start
     ELSEIF (units == 1) THEN
        omega(3) = start/rytoev
     ELSEIF (units == 2) THEN
        omega(3) = rytonm/start
     ENDIF
     !
     !-------------------------------------------------------------!
     !                          FIRST STEP                         !
     !-------------------------------------------------------------!
     !
     ! In order to gain speed, we perform the first step seperately
     !
     IF (n_ipol == 3) THEN 
        !
        ! Calculation of the susceptibility
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           !
           green(:,:) = green(:,:)/rytoev
           !
        ENDIF
        !
        DO ip=1,n_ipol
           DO ip2=1,n_ipol
              !
              WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
              !
           ENDDO
        ENDDO
        !
        alpha_temp(3) = omega(3) * aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
        !
        IF (units == 1 .or. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready
        !
        WRITE(18,'(5x,2x,2(e21.15,2x))') start, alpha_temp(3)
        !
        ! This is for the f-sum rule
        !
        f_sum = increment*alpha_temp(3)/3.0d0
        !
        start = start + increment
        !
     ENDIF
     !
     !--------------------------------------------------------------!
     !                       OMEGA LOOP                             !
     !--------------------------------------------------------------!
     !
     DO WHILE (start < end)
        !
        ! Units conversion and omega history
        !
        omega(1) = omega(2)
        omega(2) = omega(3)
        !
        IF (units == 0) THEN
           omega(3) = start
        ELSEIF (units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF (units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility for a given frequency omega.
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           !
           green(:,:) = green(:,:)/rytoev
           !
        ENDIF
        !
        ! Writing of chi
        !
        DO ip=1,n_ipol
          DO ip2=1,n_ipol
             !
             IF (n_ipol == 3) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                           & ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
             IF (n_ipol == 1) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                           & ipol, ipol, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
             !
          ENDDO
        ENDDO
      !
      IF (n_ipol==3) THEN
        !
        ! Calculation of the absorption coefficient
        !
        alpha_temp(1) = alpha_temp(2)
        alpha_temp(2) = alpha_temp(3)
        alpha_temp(3) = omega(3) * aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
        !
        IF (units == 1 .or. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready
        !
        WRITE(18,'(5x,2x,2(e21.15,2x))') start, alpha_temp(3)
        !
        IF (verbosity > 0 ) THEN
           IF ( is_peak(omega(3),alpha_temp(3))) &
               WRITE(stdout,'(5x,"Possible peak at ",F15.8," Ry; &
                                  & Intensity=",E11.2)') omega(1),alpha_temp(1)
               !
               ! f-sum rule
               !
               f_sum = f_sum + integrator(increment,alpha_temp(3))
               !
               ! Perceived color analysis
               !
               IF ( omega(3)<vis_end .and. omega(3)>vis_start .and. do_perceived ) THEN
                  !
                  perceived_intensity(perceived_iter) = alpha_temp(3)
                  perceived_evaluated(perceived_iter) = omega(3)
                  perceived_iter = perceived_iter + 1
                  !
                  ! Renormalization to 1
                  !  
                  IF (alpha_temp(3) > perceived_renorm) perceived_renorm = alpha_temp(3) 
                  !
               ENDIF
           ENDIF
        ENDIF
        !
        start = start + increment
        !
     ENDDO
     !
     !------------------------------------------------------------------!
     !                      END OF OMEGA LOOP                           !
     !------------------------------------------------------------------!
     !
     !------------------------------------------------------------------!
     !                          LAST STEP                               !
     !------------------------------------------------------------------!
     !
     ! In order to gain speed, we perform the last step seperately
     !
     IF (n_ipol == 3) THEN
        !
        ! Units conversion
        !
        IF (units == 0) THEN
           omega(3) = start
        ELSEIF (units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF (units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           green(:,:) = green(:,:)/rytoev
        ENDIF
        !
        DO ip=1,n_ipol
          DO ip2=1,n_ipol
             !
             WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
             !
          ENDDO
        ENDDO
        !
        ! Absorption coefficient
        !
        alpha_temp(3) = omega(3) * aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
        !
        IF (units == 1 .or. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready 
        !
        WRITE(18,'(5x,2x,2(e21.15,2x))') start, alpha_temp(3)
        !
        ! alpha is ready
        !
        f_sum = f_sum + increment*alpha_temp(3)/3.0d0
        !
     ENDIF
     !
     CLOSE(17)
     IF (n_ipol==3) CLOSE(18)
     !
     IF ( n_ipol==3 .and. verbosity >4 ) THEN
        !
        ! S(w)=2m_e/(pi e^2 hbar)
        WRITE(stdout,'(5x,"Integral of absorbtion coefficient ",F15.8)') f_sum
        !     
     ENDIF
     !
     !------------------------------------------------------------------------!
     !                      Perceived color analysis                          !
     !------------------------------------------------------------------------!
     !
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
     !
     !----------------------------------------------------------------------!
     !                   End of perceived color analysis                    !
     !----------------------------------------------------------------------!
     !
     ! f-sum rule
     !
     IF (verbosity > 9) THEN
        !
        start = 0.0
        f_sum = 0.0
        !
        DO WHILE (start<end)
           !
           f_sum = f_sum + integrator(increment,start)
           start = start + increment
           !
        ENDDO
        ! 
        f_sum = f_sum + increment*start/3.0d0
        !
        WRITE(stdout,'(5x,"Integral test:",F15.8,"actual: ",F15.8:)') &
                                             f_sum,0.5*start*start
        !
     ENDIF
     !
     ! Deallocations
     !
     IF (allocated(perceived_intensity)) DEALLOCATE(perceived_intensity)
     IF (allocated(perceived_evaluated)) DEALLOCATE(perceived_evaluated)
     !
     IF (allocated(beta_store)) DEALLOCATE(beta_store)
     IF (allocated(gamma_store)) DEALLOCATE(gamma_store)
     IF (allocated(zeta_store)) DEALLOCATE(zeta_store)
     !
     IF (allocated(alpha_store)) DEALLOCATE(alpha_store)
     IF (allocated(gamma_magnons_store)) DEALLOCATE(gamma_magnons_store)
     !
     DEALLOCATE(a)
     DEALLOCATE(b)
     DEALLOCATE(c)
     DEALLOCATE(r)

     RETURN
!-----------------------------------------------------------------------
END SUBROUTINE compute_optical_spectrum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE compute_eels_spectrum()
  !-------------------------------------------------------------------------
  !
  ! Compute the susceptibility in the EELS case (charge-charge response function)
  !
  IMPLICIT NONE

     outdir = trimcheck(outdir)
     tmp_dir = outdir
     !
     tmp_dir_lr = TRIM (tmp_dir) // 'tmp_eels/'
     tmp_dir = tmp_dir_lr
     !
     n_ipol = 1
     ipol = 1
     !
     ! Terminator scheme
     !
     IF (trim(extrapolation)=="no") THEN
        !
        itermax = itermax0
        !
     ENDIF
     !
     ! Check the units (Ry, eV, nm)
     !
     IF (units < 0 .or. units >3) CALL errore("lr_calculate_spectrum","Unsupported unit system",1)
     !
     IF (units == 3) CALL errore("lr_calculate_spectrum","meV unit is supported only for magnon= .true.",1) 
     !
     IF ( units /= 0 .and. verbosity > 4) THEN
        WRITE(stdout,'(5x,"Such a high verbosity is not supported when &
                        & non-default units are used")')
        verbosity = 4
     ENDIF
     !
     ! Initialisation of coefficients
     !
     ALLOCATE(beta_store(n_ipol,itermax))
     ALLOCATE(gamma_store(n_ipol,itermax))
     ALLOCATE(zeta_store(n_ipol,n_ipol,itermax))
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
     !
     ! Read beta, gamma, and zeta coefficients
     !
     CALL read_b_g_z_file()
     !
     ! Optional: use an extrapolation scheme
     !
     CALL extrapolate()
     !
     !  Spectrum calculation
     !
     WRITE (stdout,'(/5x,"Data ready, starting to calculate observables...")')
     WRITE (stdout,'(/5x,"Broadening = ",f15.8," Ry")') epsil
     !
     filename = trim(prefix) // ".plot_chi.dat"
     !
     WRITE (stdout,'(/5x,"Output file name for the susceptibility: ",A)') &
            trim(filename)
     !
     filename1 = trim(prefix) // ".plot_eps.dat"
     !
     ! Units
     !
     IF (units == 0) THEN      ! Ry
        WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega &
                              & Energy unit is (Ry)")')
     ELSEIF (units == 1) THEN  ! eV
        WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega &
                              & Energy unit is (eV)")')
     ELSEIF (units == 2) THEN  ! nm
        CALL errore("lr_calculate_spectrum", "Unsupported units (nm).",1) 
     ENDIF
     !
     ! The static dipole polarizability / static charge-density susceptibility
     !
     IF (verbosity>0) THEN
        !
        WRITE (stdout,'(/,5x,"Static charge-density susceptibility:")')
        !
        CALL calc_chi(0.0d0,epsil,green(:,:))
        !
        DO ip=1,n_ipol
          DO ip2=1,n_ipol
             !
             WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e21.15," + i",e21.15)') &
                       & ipol, ipol, dble(green(ip,ip2)), aimag(green(ip,ip2))
            !
          ENDDO
        ENDDO
        !
     ENDIF
     !
     ! Open the output file
     !
     OPEN(17,file=filename,status="unknown")
     !
     IF (units==1) THEN
        !
        ! TODO: Make an implementation also for units=0.
        !
        WRITE (stdout,'(/5x,"Output file name for the inverse and direct &
                      &dielectric function: ",A)') trim(filename1)
        !
        OPEN(18,file=filename1,status="unknown")
        !
        WRITE(18,'("#",8x,"\hbar \omega(eV)",11x,"Re(1/eps)",12x, &
                      & "-Im(1/eps)",15x,"Re(eps)",16x,"Im(eps)")')
        !
        ! Calculation of the inverse dielectric function at finite q.
        !
        !     1/eps  =  1 + (4*pi/q^2) * chi [1/eV]
        !
        ! Hartree atomic units: \hbar=1, e^2=1, m=1 
        ! 1/2 Hartree = 1 Ry = 13.6057 eV
        !
        ! q = (2*pi/a) (q1, q2, q3)
        !
        modulus_q = (2.0d0*pi/alat) * sqrt( (q1)**2 + (q2)**2 + (q3)**2 )
        !
        factor_eels = (4.0d0*pi/(modulus_q**2)) * (2.0d0*rytoev/volume)
        !
        start_save = start
        !
     ELSE
        !
        WRITE (stdout,'(/5x,"Inverse and direct dielectric function ", &
                            "requires in the input: units = 1.")')
        !
     ENDIF
     !
     ! Header of the output plot file
     !
     IF (units == 0) THEN
        WRITE (17,'("#",16x,"\hbar \omega(Ry)",5x,"Re(chi) (e^2*a_0^2/Ry)",x,"Im(chi) (e^2*a_0^2/Ry)")')
     ELSEIF (units == 1) THEN
        WRITE (17,'("#",16x,"\hbar \omega(eV)",5x,"Re(chi) (e^2*a_0^2/eV)",x,"Im(chi) (e^2*a_0^2/eV)")')
     ELSEIF (units == 2) THEN
        WRITE (17,'("#",16x,"wavelength(nm)",5x,"Re(chi) (e^2*a_0^2/eV)",x,"Im(chi) (e^2*a_0^2/eV)")')
     ENDIF
     !
     ! Start a loop on frequency
     !
     ! Units conversion and omega history
     !
     omega(1) = omega(2)
     omega(2) = omega(3)
     !
     IF (units == 0) THEN
        omega(3) = start
     ELSEIF (units == 1) THEN
        omega(3) = start/rytoev
     ELSEIF (units == 2) THEN
        omega(3) = rytonm/start
     ENDIF
     !
     !--------------------------------------------------------------!
     !                       OMEGA LOOP                             !
     !--------------------------------------------------------------!
     !
     DO WHILE (start < end)
        !
        ! Units conversion and omega history
        !
        omega(1) = omega(2)
        omega(2) = omega(3)
        !
        IF (units == 0) THEN
           omega(3) = start
        ELSEIF (units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF (units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility for a given frequency omega.
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           !
           green(:,:) = green(:,:)/rytoev
           !
        ENDIF
        !
        ! Writing of chi
        !
        DO ip=1,n_ipol
          DO ip2=1,n_ipol
             !
             WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
             & ipol, ipol, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
             !
          ENDDO
        ENDDO
      !
      ! EELS: writing of 1/eps(q)
      !
      IF (units==1) THEN
         !
         ! Calculation of the inverse dielectric function at finite q.
         !
         !     1/eps  =  1 + (4*pi/q^2) * chi [1/eV]
         !
         epsm1(1,1) = cmplx( 1.0d0 + factor_eels*dble(green(1,1)), factor_eels*aimag(green(1,1)), kind=dp ) 
         !
         ! Calculation of the direct macroscopic dielectric function
         !
         !     eps = 1/(epsm1)
         !
         eps(1,1)   = cmplx( dble(epsm1(1,1))/(dble(epsm1(1,1))**2 + aimag(epsm1(1,1))**2), & 
                            -aimag(epsm1(1,1))/(dble(epsm1(1,1))**2 + aimag(epsm1(1,1))**2), kind=dp )
         !
         !                            frequency       Re(1/eps)       -Im(1/eps)           Re(eps)        Im(eps)
         WRITE(18,'(5x,5(e21.15,2x))')  start,    dble(epsm1(1,1)), -aimag(epsm1(1,1)), dble(eps(1,1)), aimag(eps(1,1))
         !
         ! The f-sum rule (see Eq.(6) in Comput. Phys. Commun. 196, 460 (2015)).
         ! The f_sum will give the number of valence (and semicore) electrons
         ! in the unit cell.  
         !
         ! Convert the frequency from eV to Hartree
         start = start/(rytoev*2.0d0)
         increment = increment/(rytoev*2.0d0) 
         !
         integration_function = -aimag(epsm1(1,1)) * start * volume/(2.0d0*pi**2) 
         !
         f_sum = f_sum + integrator(increment,integration_function)
         !
         ! Convert the frequency back from Hartree to eV
         start = start*(rytoev*2.0d0)
         increment = increment*(rytoev*2.0d0)
         !
      ENDIF
      !
      start = start + increment
      !
     ENDDO
     !
     !------------------------------------------------------------------!
     !                      END OF OMEGA LOOP                           !
     !------------------------------------------------------------------!
     !
     CLOSE(17)
     IF (units==1) CLOSE(18)
     !
     IF (units==1) THEN
        !
        WRITE(stdout,'(/5x,"The f-sum rule is given by Eq.(6) in Comput. Phys. Commun. 196, 460 (2015).")') 
        WRITE(stdout,'(5x,"Integration in the range from",1x,f6.2,1x,"to",1x,f6.2,1x,"eV."/, &
                     & 5x,"The number of valence (and semicore) electrons in the unit cell:",1x,f6.2)') start_save, end, f_sum
        WRITE(stdout,'(5x,"The exact number of electrons:",1x,f6.2)') nelec
        WRITE(stdout,'(5x,"The violation of the f-sum rule:",1x,f6.2,1x,"%")') 100*abs(f_sum-nelec)/nelec
        !     
     ENDIF
     !
     ! Deallocations
     !
     IF (allocated(alpha_store)) DEALLOCATE(alpha_store)
     IF (allocated(beta_store)) DEALLOCATE(beta_store)
     IF (allocated(gamma_store)) DEALLOCATE(gamma_store)
     IF (allocated(zeta_store)) DEALLOCATE(zeta_store)
     !
     DEALLOCATE(a)
     DEALLOCATE(b)
     DEALLOCATE(c)
     DEALLOCATE(r)

       RETURN
!-----------------------------------------------------------------------
END SUBROUTINE compute_eels_spectrum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE compute_magnon_spectrum()
  !-------------------------------------------------------------------------
  !
  ! Compute the susceptibility in the magnon case (magnetization-magnetization response function)
  !
  IMPLICIT NONE

     REAL(DP) :: hbarw

     ! Check the units (Ry, eV, nm)
     !
     IF (units < 0 .or. units >3) CALL errore("lr_calculate_spectrum","Unsupported unit system",1)
     !
     IF (units /= 3) CALL errore("lr_calculate_spectrum","only meV unit=3 is supported for magnon= .true.",1)

     outdir = trimcheck(outdir)
     tmp_dir = outdir
     !
     tmp_dir_lr = TRIM (tmp_dir) // 'tmp_magnons/'
     tmp_dir = tmp_dir_lr
     !
     IF (ipol < 4) THEN
       n_ipol=1
     ELSE
       n_ipol=3
       ipol = 1
     ENDIF
     !
     n_op = 3
     !
     ! Terminator scheme
     !
     IF (trim(extrapolation)=="no") THEN
        !
        itermax = itermax0
        !
     ENDIF

     !
     ! Initialisation of coefficients
     !
     ALLOCATE(beta_store(n_ipol,itermax))
     ALLOCATE(zeta_store(n_ipol,n_op,itermax))
     ALLOCATE(alpha_store(n_ipol,itermax))
     ALLOCATE(gamma_magnons_store(n_ipol, itermax))
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
     !
     ! Read beta, gamma, and zeta coefficients
     !
     CALL read_b_g_z_file()
     !
     ! Optional: use an extrapolation scheme
     !
     CALL extrapolate()
     !
     !  Spectrum calculation
     !
     WRITE (stdout,'(/5x,"Data ready, starting to calculate observables...")')
     WRITE (stdout,'(/5x,"Broadening = ",f15.8," meV")') epsil
     !
     filename = trim(prefix) // ".plot_chi.dat"
     !
     WRITE (stdout,'(/5x,"Output file name: ",A)') trim(filename)
     !
     WRITE(stdout,'(/,5x,"chi_i_j: magnetization-magnetization tensor &
                        & in units of mu_B^2 / meV")')
     !
     ! Open the output file
     !
     OPEN(17,file=filename,status="unknown")
     !
     !
     WRITE(17,'("#",2x,"Chi is reported as CHI_(i)_(j) \hbar \omega (meV) &
                  & Re(chi) (mu_B^2/meV) Im(chi) (mu_B^2/meV) ")')
     !
     !--------------------------------------------------------------!
     !                       OMEGA LOOP                             !
     !--------------------------------------------------------------!
     !

     ! calc_chi needs variables in Rydberg 
     !
     epsil = epsil/rytoev/1000.d0

     DO WHILE (start < end)
        !
        ! Units conversion and omega history
        !
        hbarw = start/rytoev/1000.d0

        !
        ! Calculation of the susceptibility for a given frequency omega.
        !
        CALL calc_chi(hbarw,epsil,green(:,:))
        !
        green(:,:) = green(:,:)/rytoev/1000.d0
        !
        ! Writing of chi
        !
        DO ip=1,n_ipol
          DO ip2=1,n_op
             !
             IF (n_ipol == 3) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                           & ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
             IF (n_ipol == 1) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e21.15,2x))') &
                           & ip2, ipol, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
             !
          ENDDO
        ENDDO
        !
        start = start + increment
        !
     ENDDO
     !
     !------------------------------------------------------------------!
     !                      END OF OMEGA LOOP                           !
     !------------------------------------------------------------------!
     !
     CLOSE(17)
     !
     ! Deallocations
     !
     IF (allocated(alpha_store)) DEALLOCATE(alpha_store)
     IF (allocated(beta_store)) DEALLOCATE(beta_store)
     IF (allocated(gamma_magnons_store)) DEALLOCATE(gamma_magnons_store)
     IF (allocated(zeta_store)) DEALLOCATE(zeta_store)
     !
     DEALLOCATE(a)
     DEALLOCATE(b)
     DEALLOCATE(c)
     DEALLOCATE(r)

       RETURN
!-----------------------------------------------------------------------
END SUBROUTINE compute_magnon_spectrum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE read_b_g_z_file()
  !------------------------------------------------------------------------
  !
  ! This subroutine reads the coefficients from the file.
  !
  IMPLICIT NONE
  !
  IF (sym_op == 0) THEN
     !
     DO ip = 1, n_ipol
        !
        IF (eels) THEN
          filename = trim(prefix) // ".beta_gamma_z." // trim("dat")
        ELSE
          IF (n_ipol==3) filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ip))
          IF (n_ipol==1) filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ipol))
        ENDIF
        !
        filename = trim(tmp_dir) // trim(filename)
        !
        INQUIRE (file = filename, exist = exst)
        !
        IF (.not.exst) CALL errore("read_b_g_z_file", "Error reading file: " &
                                   // trim(filename), 1)
        !
        OPEN (158, file = filename, form = 'formatted', status = 'old')
        ! 
        READ(158,*) itermax_actual
        !
        IF (eels) THEN
           WRITE(stdout,'(/5X,"Reading ",I6," Lanczos steps")') itermax_actual
        ELSE
           WRITE(stdout,'(/5X,"Reading ",I6," Lanczos steps for direction ",I1)') &
                                                       itermax_actual, ip
        ENDIF
        WRITE(stdout,'(5X,I6," steps will be considered",/)') itermax0
        !
        IF (itermax0 > itermax_actual .or. itermax0 > itermax) THEN
           CALL errore("read_b_g_z_file", "Error in itermax0",1)
        ENDIF
        !
        ! Read the norm
        !
        READ(158,*) norm0(ip)
        !
        ! Read the degenaracy wrt spin
        !
        READ(158,*) degspin
        !
        ! Read the lattice parameter
        !
        READ(158,*) alat
        !
        ! Read the unit-cell volume
        !
        READ(158,*) volume
        !
        ! Number of valence (and semicore electrons) in the unit cell
        !
        READ(158,*) nelec
        !
        ! Read the components of the transferred momentum
        !
        READ(158,*) q1
        READ(158,*) q2
        READ(158,*) q3
        !
        ! Read the coefficients
        !
        !
        IF (magnons) THEN
           !
           READ(158,*) alpha_store(ip,1)
           READ(158,*) zeta_store(ip,:,1)
           !
           DO i = 2, itermax0
              !
              READ(158,*) alpha_store(ip,i)
              READ(158,*) beta_store(ip,i)
              READ(158,*) gamma_magnons_store(ip,i)
              READ(158,*) zeta_store (ip,:,i)
              !
           ENDDO
           !
           CLOSE(158)
           !
           IF (itermax > itermax0) THEN
              alpha_store(ip,itermax0+1:)         = (0.0d0, 0.0d0)
              beta_store(ip,itermax0+1:)          =  0.0d0
              gamma_magnons_store(ip,itermax0+1:) = (0.0d0, 0.0d0)
              zeta_store(ip,:,itermax0+1:)        = (0.0d0, 0.0d0)
           ENDIF
           !
        ELSE
           !
           DO i = 1, itermax0
              !
              READ(158,*) beta_store(ip,i)
              READ(158,*) gamma_store(ip,i)
              READ(158,*) zeta_store (ip,:,i)
              !
           ENDDO
           !
           CLOSE(158)
           !
           beta_store(ip,itermax0+1:) = 0.d0
           gamma_store(ip,itermax0+1:) = 0.d0
           zeta_store(ip,:,itermax0+1:) = (0.d0,0.d0)
           !
        ENDIF
        !
     ENDDO
     !
  ELSEIF (sym_op==1) THEN
     !
     filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(ipol))
     filename = trim(tmp_dir) // trim(filename)
     !
     INQUIRE (file = filename, exist = exst)
     !
     IF (.not.exst) CALL errore("read_b_g_z_file", "Error reading file",1)
     !
     OPEN (158, file = filename, form = 'formatted', status = 'old')
     !
     READ(158,*) itermax_actual
     !
     WRITE(stdout,'(/5X,"Reading ",I6," Lanczos steps for direction ",I1)') &
                                                        itermax_actual, ipol
     WRITE(stdout,'(5X,I6," steps will be considered")') itermax0
     !
     IF (itermax0 > itermax_actual .or. itermax0 > itermax) THEN
        CALL errore("read_b_g_z_file", "Error in Itermax0",1)
     ENDIF
     !
     ! Read the norm
     !
     READ(158,*) norm0(1)
     !   
     ! Read the degenaracy wrt spin
     !
     READ(158,*) degspin
     !
     ! Read the lattice parameter
     !
     READ(158,*) alat
     !
     ! Read the unit-cell volume
     !
     READ(158,*) volume
     !
     ! Read the components of the transferred momentum
     !
     READ(158,*) q1
     READ(158,*) q2
     READ(158,*) q3
     !
     norm0(2) = norm0(1)
     norm0(3) = norm0(1)
     !
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
     !
     beta_store(:,itermax0+1:)=0.d0
     gamma_store(:,itermax0+1:)=0.d0
     zeta_store(:,:,itermax0+1:)=(0.d0,0.d0)
     !
  ENDIF
  !
  RETURN
  !
!-----------------------------------------------------------------------
END SUBROUTINE read_b_g_z_file
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE extrapolate()
  !-----------------------------------------------------------------------
  !
  ! This subroutine applies the "extrapolation" scheme 
  ! for extrapolating the reduced matrix.
  !
  IMPLICIT NONE
  !
  !  Terminatore
  !
  skip = .false.
  !
  IF (trim(extrapolation)/="no") THEN
   !
   average = 0.d0
   av_amplitude = 0.d0
   !
   IF (magnons) THEN
      a_av    = (0.0d0, 0.0d0)
      c_av    = (0.0d0, 0.0d0)
      a_ampli = (0.0d0, 0.0d0)
      c_ampli = (0.0d0, 0.0d0)
   ENDIF
   !
   DO ip=1,n_ipol
     !
     IF (.not.eels) WRITE(stdout,'(/5x,"Polarization direction:",I1)') ip
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
              average(ip) = average(ip) + beta_store(ip,i)
              av_amplitude(ip) = av_amplitude(ip) + beta_store(ip,i)
              !
              IF (magnons) THEN
                 a_av(ip) = a_av(ip) + alpha_store(ip,i)
                 c_av(ip) = c_av(ip) + gamma_magnons_store(ip,i)
                 a_ampli(ip) = a_ampli(ip) + alpha_store(ip,i)
                 c_ampli(ip) = c_ampli(ip) + gamma_magnons_store(ip,i)
              ENDIF
              !
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
              !
              IF (magnons) THEN
                 a_av(ip) = a_av(ip) + alpha_store(ip,i)
                 c_av(ip) = c_av(ip) + gamma_magnons_store(ip,i)
                 a_ampli(ip) = a_ampli(ip) - alpha_store(ip,i)
                 c_ampli(ip) = c_ampli(ip) - gamma_magnons_store(ip,i)
              ENDIF
              !
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
     !
     IF (magnons) THEN
        a_av(ip) = a_av(ip)/counter
        c_av(ip) = c_av(ip)/counter
        a_ampli(ip) = a_ampli(ip)/counter
        c_ampli(ip) = c_ampli(ip)/counter
     ENDIF
     !
     !print *, "t3 ipol",ip,"av_amp",av_amplitude(ip)
     !
     WRITE(stdout,'(5x,"Lanczos coefficients:")')
     WRITE(stdout,'(5x,"Average =",3F15.8)') average(ip)
     WRITE(stdout,'(5x,"Average oscillation amplitude =",F15.8)') av_amplitude(ip)
     !
     IF (magnons) THEN
        WRITE(stdout,'(5x,"Average alpha =",2F15.8)') a_av(ip)
        WRITE(stdout,'(5x,"Average alpha oscillation amplitude =",2F15.8)') a_ampli(ip)
        WRITE(stdout,'(5x,"Average gamma =",2F15.8)') c_av(ip)
        WRITE(stdout,'(5x,"Average gamma oscillation amplitude =",2F15.8)') c_ampli(ip)
     ENDIF
     !
   ENDDO
   !

   IF (magnons) THEN
      !
      IF (trim(extrapolation)=="constant") THEN  
         av_amplitude=0.0d0
         a_ampli =(0.0d0, 0.0d0)
         c_ampli =(0.0d0, 0.0d0)
      ELSEIF (trim(extrapolation)=="osc") THEN
         c_av    = cmplx(average,0.0d0,dp)
         c_ampli = cmplx(av_amplitude,0.0d0,dp)
         a_av    = (0.0d0, 0.0d0)
         a_ampli = (0.0d0, 0.0d0)
      ENDIF
      !
      WRITE(stdout,'(5x,"Extrapolation = ",A24)') trim(extrapolation)
      DO ip = 1, n_ipol
         WRITE(stdout,'(5x," ip = ",I3)') ip
         WRITE(stdout,'(5x,"Average beta =",F15.8)') average(ip)
         WRITE(stdout,'(5x,"Average beta oscillation amplitude =",F15.8)') av_amplitude(ip)
         WRITE(stdout,'(5x,"Average alpha =",2F15.8)') a_av(ip)
         WRITE(stdout,'(5x,"Average alpha oscillation amplitude =",2F15.8)') a_ampli(ip)
         WRITE(stdout,'(5x,"Average gamma =",2F15.8)') c_av(ip)
         WRITE(stdout,'(5x,"Average gamma oscillation amplitude =",2F15.8)') c_ampli(ip)
      ENDDO
      !
   ELSEIF (trim(extrapolation)=="constant") THEN
      av_amplitude=0
   ENDIF
   !
   !
   IF (magnons) THEN
      !
      DO ip=1,n_ipol
         !
         DO i=itermax0+1,itermax
            !
            IF (mod(i,2)==1) THEN
               !
               alpha_store(ip,i) = a_av(ip)         + a_ampli(ip)
               beta_store(ip,i)  = average(ip)      + av_amplitude(ip)
               gamma_magnons_store(ip,i) = c_av(ip) + c_ampli(ip)
               !
            ELSE
               !
               alpha_store(ip,i) = a_av(ip)         - a_ampli(ip)
               beta_store(ip,i)  = average(ip)      - av_amplitude(ip)
               gamma_magnons_store(ip,i) = c_av(ip) - c_ampli(ip)
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
   ELSE
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
  ENDIF
  !
  !IF (verbosity > -1) THEN
      !
      ! Write all the coefficients in a file for a detailed post-processing
      !
      DO ip = 1, n_ipol
         !
         IF (eels) THEN
           filename = trim(prefix) // ".beta_term.dat"
         ELSE
           IF (n_ipol==3) filename = trim(prefix) // ".beta_term." // trim(int_to_char(ip))
           IF (n_ipol==1) filename = trim(prefix) // ".beta_term." // trim(int_to_char(ipol))
         ENDIF
         !
         filename = trim(tmp_dir) // trim(filename)
         !
         OPEN (17, file = filename, status = 'unknown')
         !
         IF (magnons) THEN
            !
            DO i = 1,itermax 
               WRITE(17,'(i5,2x,e21.15)') i,beta_store(ip,i)
            ENDDO
            !
         ELSE
            !
            ! Write even iterations
            !
            DO i = 1, itermax0
               IF (mod(i,2)==0) WRITE(17,'(i5,2x,e21.15)') i, beta_store(ip,i)
            ENDDO
            !
            WRITE(17,*) '  '
            !
            ! Write odd iterations
            !
            DO i = 1, itermax0
               IF (mod(i,2)==1) WRITE(17,'(i5,2x,e21.15)') i, beta_store(ip,i)
            ENDDO
            !
         ENDIF
         !
         CLOSE(17)
         !
         IF (magnons) THEN
            !
            ! alpha coefficients
            !
            IF (n_ipol==3) filename = trim(prefix) // ".alpha_term." // trim(int_to_char(ip))
            IF (n_ipol==1) filename = trim(prefix) // ".alpha_term." // trim(int_to_char(ipol))
            !
            filename = trim(tmp_dir) // trim(filename)
            !
            OPEN (17, file = filename, status = 'unknown')
            !
            DO i = 1,itermax
               WRITE(17,'(i5,2x,e21.15,2x,e21.15)') i, dble(alpha_store(ip,i)), aimag(alpha_store(ip,i))
            ENDDO
            !
            CLOSE(17)
            !
            ! gamma coefficients
            !
            IF (n_ipol==3) filename = trim(prefix) // ".gamma_term." // trim(int_to_char(ip))
            IF (n_ipol==1) filename = trim(prefix) // ".gamma_term." // trim(int_to_char(ipol))
            !
            filename = trim(tmp_dir) // trim(filename)
            !
            OPEN (17, file = filename, status = 'unknown')
            !
            DO i = 1,itermax
               WRITE(17,'(i5,2x,e21.15,2x,e21.15)') i, dble(gamma_magnons_store(ip,i)), aimag(gamma_magnons_store(ip,i))
            ENDDO
            !
            CLOSE(17)
            !
         ENDIF
         !
      ENDDO
      !
  !ENDIF
  !
  RETURN
  !
!-----------------------------------------------------------------------
END SUBROUTINE extrapolate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE calc_chi(freq,broad,chi)
  !-----------------------------------------------------------------------------
  !
  ! This subroutine Calculates the susceptibility.
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), INTENT(in) :: freq
  REAL(kind=dp), INTENT(in) :: broad
  COMPLEX(kind=dp), INTENT(out) :: chi(:,:)
  !
  omeg_c = cmplx(freq,broad,dp)
  !
  DO ip =1, n_ipol
     !
     IF (magnons) THEN
        a(1) = omeg_c - alpha_store(ip,1)
        !
        DO i = 1,itermax-1
           !
           a(i+1) = omeg_c - alpha_store(ip,i+1)
           b(i) = cmplx(-beta_store(ip,i+1),0.0d0,dp)
           c(i) = -gamma_magnons_store(ip,i+1)
           !
        ENDDO
     ELSE
        a(:) = omeg_c
        !
        DO i = 1,itermax-1
           !
           b(i) = cmplx(-beta_store(ip,i),0.0d0,dp)
           c(i) = cmplx(-gamma_store(ip,i),0.0d0,dp)
           !
        ENDDO
     ENDIF
     !
     r(ip,:) = (0.0d0,0.0d0)
     r(ip,1) = (1.0d0,0.0d0)
     !
     ! |w_t|=(w-L) |1,0,0,...,0|
     ! 
     CALL zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info)
     !
     IF (info /= 0) CALL errore("calc_chi", "Unable to solve tridiagonal system",1)
     !
     ! p=-div.rho'
     ! p= chi . E
     ! Thus, chi = - <zeta|w_t>
     !
     ! Notice that brodening has a positive sign, 
     ! thus the abs. coefficient is Im(tr(chi)) not -Im(Tr(chi)) as usual
     ! 
     IF (magnons) THEN
        DO ip2 = 1,n_op
           !
           chi(ip,ip2) = dot_product(zeta_store(ip,ip2,:),r(ip,:))
           !
           ! Multiplication with a norm
           !
           chi(ip,ip2) = chi(ip,ip2) * cmplx(norm0(ip),0.0d0,dp)
           !
           ! The response charge density is defined as 2*evc0*q, see Eq. (43) in
           ! JCP 128, 154105 (2008).
           ! Therefore, the dipole is given by 2*degspin* zeta^T *
           ! (w-T^itermax)^-1 * e_1. See also Eq. (15) in that paper.
           ! Optics: The minus sign accounts for the negative electron charge
           ! (perturbation is -e E x, rather than E x)
           !
           !
           chi(ip,ip2) = chi(ip,ip2) * cmplx(degspin, 0.d0, dp)
           !
        ENDDO
     ELSE
        DO ip2 = 1,n_ipol
           !
           chi(ip,ip2) = dot_product(zeta_store(ip,ip2,:),r(ip,:))
           !
           ! Multiplication with a norm
           !
           chi(ip,ip2) = chi(ip,ip2) * cmplx(norm0(ip),0.0d0,dp)
           !
           ! The response charge density is defined as 2*evc0*q, see Eq. (43) in
           ! JCP 128, 154105 (2008). 
           ! Therefore, the dipole is given by 2*degspin* zeta^T *
           ! (w-T^itermax)^-1 * e_1. See also Eq. (15) in that paper.
           ! Optics: The minus sign accounts for the negative electron charge
           ! (perturbation is -e E x, rather than E x)
           !
           ! 
           IF (eels) THEN
              chi(ip,ip2) = chi(ip,ip2) * cmplx( 2.d0*degspin, 0.d0, dp)
           ELSE
              chi(ip,ip2) = chi(ip,ip2) * cmplx(-2.d0*degspin, 0.d0, dp)
           ENDIF
           !
        ENDDO
     ENDIF 
     !
  ENDDO
  !
  RETURN
  !
!-----------------------------------------------------------------------
END SUBROUTINE calc_chi
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE wl_to_color(wavelength,red,green,blue)
  !----------------------------------------------------------------------------
  !
  ! Gives the colour intensity of a given wavelength 
  ! in terms of RGB (red, green and blue).
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), INTENT(in) :: wavelength
  REAL(kind=dp), INTENT(out) :: red,green,blue
  !
  IF ((wavelength>=380.).and.(wavelength<=440.)) THEN
     red = -1.*(wavelength-440.)/(440.-380.)
     green = 0.
     blue = 1.
  ENDIF
  !
  IF ((wavelength>=440.).and.(wavelength<=490.)) THEN
     red = 0.
     green = (wavelength-440.)/(490.-440.)
     blue = 1.
  ENDIF 
  !
  IF ((wavelength>=490.).and.(wavelength<=510.)) THEN
     red = 0.
     green = 1.
     blue = -1.*(wavelength-510.)/(510.-490.)
  ENDIF
  !
  IF ((wavelength>=510.).and.(wavelength<=580.)) THEN
     red = (wavelength-510.)/(580.-510.)
     green = 1.
     blue = 0.
  ENDIF
  !
  IF ((wavelength>=580.).and.(wavelength<=645.)) THEN
     red = 1.
     green = -1.*(wavelength-645.)/(645.-580.)
     blue = 0.
  ENDIF
  !
  IF ((wavelength>=645.).and.(wavelength<=780.)) THEN
     red = 1.
     green = 0.
     blue = 0.
  ENDIF
  !
  RETURN
  !
!-----------------------------------------------------------------------
END SUBROUTINE wl_to_color
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE spectrum_david()
  !-----------------------------------------------------------------
  !
  ! This routine is used to compute an absorption spectrum
  ! when the Davidson algorithm is used.
  !
  ! Written by X. Ge (2013)
  !
  IMPLICIT NONE
  !
  REAL(dp) :: energy, frequency, temp, chi(4)
  INTEGER :: ieign, nstep, istep
  REAL(dp), ALLOCATABLE :: absorption(:,:)
  !
  OPEN(18,file=trim(eign_file),action='read')
  !
  nstep = (end-start)/increment+1
  !
  ! Column 1: Energy; 2: Total; 3,4,5: X,Y,Z
  !   
  ALLOCATE(absorption(nstep,5)) 
  absorption(:,:) = 0.0d0
  !
  READ(18,*)  ! Jump to the second line
522 READ(18,*,END=521)  energy, chi
  !   
  frequency = start
  !
  istep = 1
  !
  DO WHILE( .not. istep .gt. nstep )
     !
     absorption(istep,1) = frequency
     !
     temp = frequency - energy
     temp = epsil/(temp**2 + epsil**2)
     !
     absorption(istep,2) = absorption(istep,2) + chi(1)*temp
     absorption(istep,3) = absorption(istep,3) + chi(2)*temp
     absorption(istep,4) = absorption(istep,4) + chi(3)*temp
     absorption(istep,5) = absorption(istep,5) + chi(4)*temp
     !
     istep = istep + 1
     !
     frequency = frequency + increment
     !
  ENDDO
  !
  GOTO 522
  !
521 CLOSE(18)
  !
  filename = trim(prefix)//".plot.dat"
  !
  OPEN(17,file=filename,status="unknown")
  !
  write(17,'("#",7x,"Energy(Ry)",12x,"Total",17x,"X",18x,"Y",19x,"Z")')  
  !    
  istep = 1
  !
  do while( .not. istep .gt. nstep )
      write(17,'(5E20.8)') absorption(istep,1),absorption(istep,1)*absorption(istep,2),&
                           absorption(istep,1)*absorption(istep,3),absorption(istep,1)*&
                           absorption(istep,4),absorption(istep,1)*absorption(istep,5)
      istep = istep+1
  enddo
  !
  PRINT *, "   The spectrum is in file: ", filename
  !
  CLOSE(17)
  !
  RETURN
  !
!-----------------------------------------------------------------------
END SUBROUTINE spectrum_david
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
LOGICAL FUNCTION is_peak(omeg,alpha)
  !-------------------------------------------------------------------------
  !
  ! A simple algorithm for detecting peaks.
  ! Increments of omega between alpha steps should be constant
  ! omega must increase monothonically
  ! no checks performed!
  ! OBM 2010
  !
  IMPLICIT NONE
  !Input and output
  REAL(kind=dp),INTENT(in) :: omeg, alpha !x and y
  !
  ! Local variables
  !
  REAL(kind=dp),SAVE :: omeg_save = 0.0d0, &
                      & thm1, h2m1,&
                      & first_der_save=9.0d99
  REAL(kind=dp),SAVE :: alpha_save(3) = 0.0d0
  INTEGER, SAVE :: current_iter = 0
  LOGICAL, SAVE :: trigger=.true.
  REAL(kind=dp) :: first_der, second_der
  !
  is_peak = .false.
  ! counter
  ! Rotate the variables
  !   
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
  !
  !The derivatives
  first_der = (alpha_save(3)-alpha_save(1))*thm1
  second_der = (alpha_save(3)-2.0d0*alpha_save(2)+alpha_save(1))*h2m1 
  ! second derivative corresponds to t, 3 steps before
  !first_der = (-alpha_save(5)+8.0d0*(alpha_save(4)-alpha_save(2))+alpha_save(1))*thm1 
  !first derivative corresponds to t, 3 steps before
  !second_der = (alpha_save(4)-2.0d0*alpha_save(3)+alpha_save(2))*h2m1 
  ! second derivative corresponds to t, 3 steps before
  !Decide
  !print *,"w",omeg-0.25d0/thm1,"f=",abs(first_der),"s=",second_der
  !print *,"w",omeg-0.5d0/thm1,"f=",abs(first_der),"s=",second_der
  !if (abs(first_der) < 1.0d-8 .and. second_der < 0 ) is_peak=.true.
  !  
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
  !
!-----------------------------------------------------------------------
END FUNCTION is_peak
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
REAL(kind=dp) FUNCTION integrator(dh,alpha)
  !------------------------------------------------------------------------
  !
  ! This function calculates an integral every 
  ! three points, using the Simpson's rule.
  !
  IMPLICIT NONE
  !Input and output
  REAL(kind=dp),INTENT(in) :: dh, alpha !x and y
  !internal
  LOGICAL,SAVE :: flag=.true.
  !
  ! COMPOSITE SIMPSON INTEGRATOR, (precision level ~ float)
  ! \int a b f(x) dx = ~ h/3 (f(a) + \sum_odd-n 2*f(a+n*h) + \sum_even-n 4*f(a+n*h) +f(b))
  !
  integrator = 0.0d0
  !
  IF (flag) THEN 
     ! odd steps
     integrator = (4.0d0/3.0d0)*dh*alpha
     flag = .false.
  ELSE
     ! even steps
     integrator = (2.0d0/3.0d0)*dh*alpha
     flag = .true.
  ENDIF
  !
  RETURN
  !
!-----------------------------------------------------------------------
END FUNCTION integrator
!-----------------------------------------------------------------------


END PROGRAM lr_calculate_spectrum
!-----------------------------------------------------------------------

