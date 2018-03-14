!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! by P. Bonfa', F. Affinito and C. Cavazzoni, Cineca
!  & S. de Gironcoli, SISSA

module timers
  save
  LOGICAL :: ignore_time = .true.  ! This is used to avoid collection of initialization times
  REAL*8  :: times(20) = 0.d0      ! Array hosting various timers
end module

program test
  !! This mini-app provides a tool for testing and benchmarking the FFT drivers
  !! contained in the FFTXlib.
  !!
  !! The mini-app mimics the workload of vloc_psi (a charge-density transformation
  !! from complex-to-real space and back and contribution to h).
  !!
  !! To compile the test program, once you have properly edit the make.sys file
  !! included in the FFTXlib and type:
  !!
  !!      make TEST
  !!      
  !! N.B.: do not run the make command alone, otherwise the FFT times will
  !!       not be present in the final summary.
  !!
  !! Then you can run your FFT tests using command like:
  !!
  !!      mpirun -np 4 ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4 -gamma .true.
  !!
  !! or, in case of serial build
  !!
  !!      ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4
  !!
  !! Command line arguments:
  !!
  !!-ecutwfc  Plane wave energy cut off
  !!
  !!-ecutrho  Potential energy cut off (not used)
  !!
  !!-alat     Lattice parameter
  !!
  !!-nbnd     Number of bands (fft cycles)
  !!
  !!-ntg      Number of task groups
  !!
  !!-gamma    Enables gamma point trick. Should be about 2 times faster.
  !!
  !!-av1  x y z    First lattice vector, in atomic units. N.B.: when using -av1, -alat is ignored!
  !!
  !!-av2  x y z    Second lattice vector, in atomic units. N.B.: when using -av2, -alat is ignored!
  !!
  !!-av3  x y z    Third lattice vector, in atomic units. N.B.: when using -av3, -alat is ignored!
  !!
  !!-kmax kx ky kz    Reciprocal lattice vector inside the BZ with maximum norm. Used to calculate max(|G+K|). (2pi/a)^2 units.
  !!
  !! Timings of different stages of execution are provided at the end of the
  !! run.
  !! In the present version, a preliminar implementation with non-blocking MPI
  !! calls as been implemented. This version requires the precompilation flags
  !! -D__NON_BLOCKING_SCATTER
  !!
  USE fft_types
  USE stick_base
  USE fft_parallel
  USE fft_support
  USE fft_helper_subroutines
  USE fft_interfaces, ONLY:fwfft, invfft
  USE timers
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor) :: dfftp, dffts, dfft3d
  !
  TYPE(sticks_map) :: smap
  INTEGER :: nx = 128
  !! grid points along x (modified after)
  INTEGER :: ny = 128
  !! grid points along y (modified after)
  INTEGER :: nz = 256
  !! grid points along z (modified after)
  !
  INTEGER :: mype, npes, comm, root
  !! MPI handles
  INTEGER :: ntgs
  !! number of taskgroups
  INTEGER :: nbnd
  !! number of bands
  LOGICAL :: iope
  !! I/O process
  INTEGER :: ierr, i, j, ncount, ib
  INTEGER :: incr=1, right_nr3
  INTEGER :: ngw_, ngm_, ngs_
  REAL*8  :: gcutm, gkcut, gcutms
  REAL*8  :: ecutm, ecutw, ecutms
  REAL*8  :: ecutrho, dual, kmax(3)
  !! cut-off for density
  REAL*8  :: ecutwfc
  !! cut-off for the wave-function
  REAL*8  :: tpiba, alat, alat_in
  !! lattice parameters
  REAL*8  :: time(100)
  REAL*8  :: my_time(100)
  REAL*8  :: time_min(100)
  REAL*8  :: time_max(100)
  REAL*8  :: time_avg(100)
  REAL*8  :: wall
  REAL*8  :: wall_avg
  !
  LOGICAL :: gamma_only = .false.
  LOGICAL :: use_tg
  !! if calculations require only gamma point
  REAL*8  :: at(3, 3), bg(3, 3)
  REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP
  !
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  COMPLEX(DP), ALLOCATABLE :: psic(:)
  COMPLEX(DP), ALLOCATABLE :: psi(:, :)
  !! fake argument returned by the FFT
  REAL(DP), ALLOCATABLE :: v(:)
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:, :)
  !! array representing the potential
  INTEGER :: ngms, ngsx, ngms_g
  INTEGER, ALLOCATABLE :: mill(:, :), ig_l2g(:)
  REAL(DP), ALLOCATABLE :: g(:, :), gg(:)
  INTEGER :: ngm, ngmx, ngm_g, gstart
  !
  !
  integer :: nargs
  CHARACTER(LEN=80) :: arg
  !
#if defined(_OPENMP)
  INTEGER :: PROVIDED
#endif
  !
  !   ........
  !
  !   default parameter (32 water molecules)
  !
  kmax    = HUGE(0.d0)   !not set
  ecutwfc = 80.0d0
  ecutrho = 0.d0
  alat_in = 18.65
  ntgs    = 1
  nbnd    = 1
  !
  at(1, :) = (/0.5d0, 1.0d0, 0.0d0/)
  at(2, :) = (/0.5d0, 0.0d0, 0.5d0/)
  at(3, :) = (/0.0d0, 0.5d0, 1.5d0/)
  !
  DO i = 1, iargc()
    CALL getarg(i, arg)
    IF ((trim(arg) == "-h") .or. (trim(arg) == "--help")) THEN
      print *, "Example usage:  mpirun -np 4 ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4 -gamma .true."
      print *, ""
      print *, "For detailed information see test.f90"
      STOP
    END IF
  END DO
 
  nargs = command_argument_count()
  do i = 1, nargs - 1
    CALL get_command_argument(i, arg)
    IF (TRIM(arg) == '-ecutrho') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) ecutrho
    END IF
    IF (TRIM(arg) == '-ecutwfc') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) ecutwfc
    END IF
    IF (TRIM(arg) == '-alat') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) alat_in
    END IF
    IF (TRIM(arg) == '-av1') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) at(1, 1)
      CALL get_command_argument(i + 2, arg)
      READ (arg, *) at(2, 1)
      CALL get_command_argument(i + 3, arg)
      READ (arg, *) at(3, 1)
      alat_in = 1.0
    END IF
    IF (TRIM(arg) == '-av2') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) at(1, 2)
      CALL get_command_argument(i + 2, arg)
      READ (arg, *) at(2, 2)
      CALL get_command_argument(i + 3, arg)
      READ (arg, *) at(3, 2)
      alat_in = 1.0
    END IF
    IF (TRIM(arg) == '-av3') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) at(1, 3)
      CALL get_command_argument(i + 2, arg)
      READ (arg, *) at(2, 3)
      CALL get_command_argument(i + 3, arg)
      READ (arg, *) at(3, 3)
      alat_in = 1.0
    END IF
    IF (TRIM(arg) == '-kmax') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) kmax(1)
      CALL get_command_argument(i + 2, arg)
      READ (arg, *) kmax(2)
      CALL get_command_argument(i + 3, arg)
      READ (arg, *) kmax(3)
      alat_in = 1.0
    END IF
    IF (TRIM(arg) == '-ntg') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) ntgs
    END IF
    IF (TRIM(arg) == '-nbnd') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) nbnd
    END IF
    IF (TRIM(arg) == '-gamma') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) gamma_only
    END IF
  end do
  if (ecutrho == 0.d0) ecutrho = 4.0d0*ecutwfc

#if defined(__MPI)

#if defined(_OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
  CALL MPI_Init(ierr)
#endif
  CALL mpi_comm_rank(MPI_COMM_WORLD, mype, ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  comm = MPI_COMM_WORLD
  root = 0
  IF (mype == root) THEN
    iope = .true.
  ELSE
    iope = .false.
  ENDIF

#else

  mype = 0
  npes = 1
  comm = 0
  ntgs = 1
  root = 0
  iope = .true.

#endif
  !
  !  Broadcast input parameter first
  !
  CALL MPI_BCAST(ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(alat_in, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(ntgs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(nbnd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !
  !
  ! --------  INITIALIZE DIMENSIONS AND DESCRIPTORS
  !
  !
  ecutw  = ecutwfc
  ! dual
  ecutm  = ecutrho
  ecutms = ecutrho
  !
  at = at*alat_in
  !
  alat = SQRT ( at(1,1)**2+at(2,1)**2+at(3,1)**2 )
  !
  at(:,:) = at(:,:) / alat
  !
  tpiba = 2.0d0*pi/alat
  ! 
  call recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
  !
  !
  IF (gamma_only) kmax = 0.d0
  !
  IF (ALL(kmax < HUGE(0.d0))) THEN
     gkcut = sqrt ( sum(kmax(1:3)**2) ) + sqrt (ecutw/ tpiba**2)
  ELSE
     ! use max(bg)/2 as an estimate of the largest k-point
     gkcut = 0.5d0 * max ( &
        sqrt (sum(bg (1:3, 1)**2) ), &
        sqrt (sum(bg (1:3, 2)**2) ), &
        sqrt (sum(bg (1:3, 3)**2) ) ) + sqrt (ecutw/ tpiba**2)
  END IF
  !
  gkcut = gkcut**2  ! wave function cut-off
  !
  dual = ecutrho / ecutwfc 
  gcutm = dual * ecutwfc / tpiba**2  ! potential cut-off
  !
  IF ( dual > 4.000000001_dp ) THEN
     gcutms = 4.D0 * ecutwfc / tpiba**2 ! smooth mesh cut-off
  ELSE
     gcutms = gcutm                     ! smooth mesh cut-off
  END IF
  !
  if( mype == 0 ) then

    write (*, *) '+-----------------------------------+'
    write (*, *) '|               QE FFT              |'
    write (*, *) '|          testing & timing         |'
    write (*, *) '|         by Carlo Cavazzoni        |'
    write (*, *) '+-----------------------------------+'
    write (*, *)
    write (*, *) 'alat    = ', alat
    write (*, *) 'Ecutwfc = ', ecutw
    write (*, *) 'Ecutrho = ', ecutm
    write (*, *) 'Ecuts   = ', ecutms
    write (*, *) 'Gcutrho = ', SQRT(gcutm)
    write (*, *) 'Gcuts   = ', SQRT(gcutms)
    write (*, *) 'Gcutwfc = ', SQRT(gkcut)
    write (*, *) 'Num bands      = ', nbnd
    write (*, *) 'Num procs      = ', npes
    write (*, *) 'Num Task Group = ', ntgs
    write (*, *) 'Gamma trick    = ', gamma_only
  end if
  !
  nx = 2*int(sqrt(gcutm)*sqrt(at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2)) + 1
  ny = 2*int(sqrt(gcutm)*sqrt(at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2)) + 1
  nz = 2*int(sqrt(gcutm)*sqrt(at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2)) + 1
  !
  !
  if (mype == 0) then
    write (*, *) 'nx = ', nx, ' ny = ', ny, ' nz = ', nz
  end if
  !
  IF (gamma_only) incr = 2
  dffts%has_task_groups = (ntgs > 1)
  use_tg = dffts%has_task_groups
  !
  dffts%rho_clock_label='ffts' ; dffts%wave_clock_label='fftw'
  CALL fft_type_init(dffts, smap, "wave", gamma_only, .true., comm, at, bg, gkcut, gcutms/gkcut, nyfft=ntgs)
  dfftp%rho_clock_label='fft' 
  CALL fft_type_init(dfftp, smap, "rho", gamma_only, .true., comm, at, bg, gcutm, 4.d0, nyfft=ntgs)
  !
  CALL fft_base_info(mype == 0, dffts, dfftp)
  if (mype == 0) then
    write (*, *) 'dffts:  nr1 = ', dffts%nr1, ' nr2 = ', dffts%nr2, ' nr3 = ', dffts%nr3
    write (*, *) '        nr1x= ', dffts%nr1x, ' nr2x= ', dffts%nr2x, ' nr3x= ', dffts%nr3x
  end if
  !
  ngw_ = dffts%nwl(dffts%mype + 1)
  ngs_ = dffts%ngl(dffts%mype + 1)
  ngm_ = dfftp%ngl(dfftp%mype + 1)
  !
  IF (gamma_only) THEN
    ngw_ = (ngw_ + 1)/2
    ngs_ = (ngs_ + 1)/2
    ngm_ = (ngm_ + 1)/2
  END IF
  !
  ngms = ngs_
  CALL MPI_ALLREDUCE(ngms, ngsx, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(ngms, ngms_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  ngm = ngm_
  CALL MPI_ALLREDUCE(ngm, ngmx, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(ngm, ngm_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  !
  ! --------  ALLOCATE
  !
  ALLOCATE (psic(dffts%nnr))
  ALLOCATE (psi(ngms, nbnd))
  ALLOCATE (hpsi(ngms, nbnd))
  ALLOCATE (v(dffts%nnr))
  ALLOCATE (mill(3, ngm))
  ALLOCATE (g(3, ngm))
  ALLOCATE (gg(ngm))
  ALLOCATE (ig_l2g(ngm))
  !
  ! --------  GENERATE G-VECTORS
  !
  ! smallmem = .flase.
  IF( .false. ) THEN
     CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
          g, gg, mill, ig_l2g, gstart, .TRUE. )
  ELSE
     CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart, .FALSE.  )
  END IF
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
  !
  ! --------  RESET TIMERS
  !
  time = 0.0d0
  my_time = 0.0d0
  time_min = 0.0d0
  time_max = 0.0d0
  time_avg = 0.0d0
  !
  ! Test FFT for wave functions - First calls may be biased by MPI and FFT initialization
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  IF (use_tg) THEN
    ALLOCATE (tg_psic(dffts%nnr_tg))
    CALL invfft('tgWave', tg_psic, dffts)
    DEALLOCATE (tg_psic)
  ELSE
    CALL invfft('Wave', psic, dffts)
  END IF
  !
  IF (use_tg) THEN
    ALLOCATE (tg_psic(dffts%nnr_tg))
    CALL fwfft('tgWave', tg_psic, dffts)
    DEALLOCATE (tg_psic)
  ELSE
    CALL fwfft('Wave', psic, dffts)
  END IF
  ! Now for real,
  !
  ! --------  RECORD TIMES
  !
  wall = MPI_WTIME()
  ignore_time = .false.
  !
  ! --------  INITIALIZE WAVE FUNCTIONS psi
  !
  psi = 0.0d0
  !
  ! --------  INITIALIZE POTENTIAL v
  !
  v = 1.0d0
  !
  IF (use_tg) THEN
    !
    ALLOCATE (tg_v(dffts%nnr_tg))
    ALLOCATE (tg_psic(dffts%nnr_tg))
    !
    CALL tg_gather(dffts, v, tg_v)
    !      incr is already set to 2 for gamma_only
    incr = incr*fftx_ntgrp(dffts)
    !
  ENDIF

  !
  ! Execute FFT calls once more and Take time
  !
  ncount = 0
  !
  !
  !
  DO ib = 1, nbnd, incr
    !
    time(1) = MPI_WTIME()
    !
    IF (use_tg) THEN
      !
      call prepare_psi_tg(ib, nbnd, ngms, psi, tg_psic, dffts, gamma_only)
      time(2) = MPI_WTIME()
      !
      CALL invfft('tgWave', tg_psic, dffts); 
      time(3) = MPI_WTIME()
      !
      CALL tg_get_group_nr3(dffts, right_nr3)
      !
      DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
        tg_psic(j) = tg_psic(j)*tg_v(j)
      ENDDO
      !
      time(4) = MPI_WTIME()
      !
      CALL fwfft('tgWave', tg_psic, dffts); 
      time(5) = MPI_WTIME()
      !
      CALL accumulate_hpsi_tg(ib, nbnd, ngms, hpsi, tg_psic, dffts, gamma_only)
      time(6) = MPI_WTIME()
    ELSE
      !
      call prepare_psi(ib, nbnd, ngms, psi, psic, dffts, gamma_only)
      time(2) = MPI_WTIME()
      !
      CALL invfft('Wave', psic, dffts); time(3) = MPI_WTIME()
      !
      DO j = 1, dffts%nnr
        psic(j) = psic(j)*v(j)
      ENDDO
      time(4) = MPI_WTIME()
      !
      CALL fwfft('Wave', psic, dffts); 
      time(5) = MPI_WTIME()
      !
      CALL accumulate_hpsi(ib, nbnd, ngms, hpsi, psic, dffts, gamma_only)
      time(6) = MPI_WTIME()
      !
    ENDIF
    !
    do i = 2, 6
      my_time(i) = my_time(i) + (time(i) - time(i - 1))
    end do
    !
    ncount = ncount + 1
    !
  enddo
  !
  wall = MPI_WTIME() - wall

  ! --------  DEALLOCATE
  !
  DEALLOCATE (psic)
  DEALLOCATE (psi)
  DEALLOCATE (hpsi)
  DEALLOCATE (v)
  DEALLOCATE (mill)
  DEALLOCATE (g)
  DEALLOCATE (gg)
  DEALLOCATE (ig_l2g)
  IF (use_tg) THEN
    !
    DEALLOCATE (tg_psic)
    DEALLOCATE (tg_v)
    !
  ENDIF
  ! --------------------
  !
#if defined(__MPI)
  CALL MPI_ALLREDUCE(my_time, time_min, 10, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(my_time, time_max, 10, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(my_time, time_avg, 10, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(wall, wall_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  time_min = time
  time_max = time
  time_avg = time
#endif

  time_avg = time_avg / npes
  wall_avg = wall_avg / npes

  if( mype == 0 ) then
    
    write(*,*) '**** QE 3DFFT Timing ****'
    write(*,*) 'grid size = ', dffts%nr1, dffts%nr2, dffts%nr3
    write(*,*) 'num proc  = ', npes
    write(*,*) 'num band  = ', nbnd
    write(*,*) 'num task group  = ', ntgs
    write(*,*) 'num fft cycles  = ', ncount

    write(*,100) 
    write(*,1) 
    write(*,100) 
    write(*,2) time_min(2), time_max(2), time_avg(2)
    write(*,3) time_min(3), time_max(3), time_avg(3)
    write(*,4) time_min(4), time_max(4), time_avg(4)
    write(*,5) time_min(5), time_max(5), time_avg(5)
    write(*,6) time_min(6), time_max(6), time_avg(6)
    write(*,7) wall 
    write(*,100) 

100 FORMAT(' +--------------------+----------------+-----------------+----------------+' )
1   FORMAT(' |FFT TEST subroutine |  sec. min      | sec. max        | sec.  avg      |' )
2   FORMAT(' |prepare_psi         | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3,    ' |' )
3   FORMAT(' |invfft              | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3,    ' |' )
4   FORMAT(' |workload            | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,   ' |')
5   FORMAT(' |fwfft               | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,   ' |')
6   FORMAT(' |accumulate_hpsi     | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,   ' |')
7   FORMAT(' |wall time           | ',    D14.5, ' |')

  end if
  ! now print FFT clocks
  call print_clock(mype, npes, ncount)

  CALL fft_type_deallocate(dffts)
  CALL fft_type_deallocate(dfftp)
  CALL fft_type_deallocate(dfft3d)
  CALL sticks_map_deallocate( smap )

#if defined(__MPI)
  CALL mpi_finalize(ierr)
#endif


contains



!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !     first the input variables
  !
  implicit none
  real(DP) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
  ! input: first direct lattice vector
  ! input: second direct lattice vector
  ! input: third direct lattice vector
  ! output: first reciprocal lattice vector
  ! output: second reciprocal lattice vector
  ! output: third reciprocal lattice vector
  !
  !   then the local variables
  !
  real(DP) :: den, s
  ! the denominator
  ! the sign of the permutations
  integer :: iperm, i, j, k, l, ipol
  ! counter on the permutations
  !\
  !  Auxiliary variables
  !/
  !
  ! Counter on the polarizations
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.d0
100 do iperm = 1, 3
     den = den + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s
  if (s.lt.0.d0) goto 100
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
     b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
     b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
     l = i
     i = j
     j = k
     k = l
  enddo
  return
end subroutine recips


  SUBROUTINE ggen ( dfftp, gamma_only, at, bg,  gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart, no_global_sort )
    !----------------------------------------------------------------------
    !
    !     This routine generates all the reciprocal lattice vectors
    !     contained in the sphere of radius gcutm. Furthermore it
    !     computes the indices nl which give the correspondence
    !     between the fft mesh points and the array of g vectors.
    !
    USE fft_types, ONLY: fft_stick_index, fft_type_descriptor
    USE fft_ggen, ONLY : fft_set_nl
    USE fft_param
    !USE mp, ONLY: mp_rank, mp_size, mp_sum
    !USE constants, ONLY : eps8
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor),INTENT(INOUT) :: dfftp
    LOGICAL,  INTENT(IN) :: gamma_only
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3), gcutm
    INTEGER, INTENT(IN) :: ngm_g
    INTEGER, INTENT(INOUT) :: ngm
    REAL(DP), INTENT(OUT) :: g(:,:), gg(:)
    INTEGER, INTENT(OUT) :: mill(:,:), ig_l2g(:), gstart
    !  if no_global_sort is present (and it is true) G vectors are sorted only
    !  locally and not globally. In this case no global array needs to be
    !  allocated and sorted: saves memory and a lot of time for large systems.
    !
    LOGICAL, INTENT(IN) :: no_global_sort
    !
    !     here a few local variables
    !
    REAL(DP) :: tx(3), ty(3), t(3)
    REAL(DP), ALLOCATABLE :: tt(:)
    INTEGER :: ngm_save, n1, n2, n3, ngm_offset, ngm_max, ngm_local
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing only g vectors for the current processor
    INTEGER, ALLOCATABLE :: mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors
    ! (replicated data). When no_global_sort is present and .true.,
    ! only g-vectors for the current processor are stored
    INTEGER, ALLOCATABLE :: igsrt(:), g2l(:)
    !
    INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw, ierr
    INTEGER :: istart, jstart, kstart
    INTEGER :: mype, npe
    LOGICAL :: global_sort, is_local
    INTEGER, ALLOCATABLE :: ngmpe(:)
    !
    ! The 'no_global_sort' is not optional in this case.
    ! This differs from the version present in QE distribution.
    global_sort = .NOT. no_global_sort
    !
    IF( .NOT. global_sort ) THEN
       ngm_max = ngm
    ELSE
       ngm_max = ngm_g
    END IF
    !
    ! save current value of ngm
    !
    ngm_save  = ngm
    !
    ngm = 0
    ngm_local = 0
    !
    !    set the total number of fft mesh points and and initial value of gg
    !    The choice of gcutm is due to the fact that we have to order the
    !    vectors after computing them.
    !
    gg(:) = gcutm + 1.d0
    !
    !    and computes all the g vectors inside a sphere
    !
    ALLOCATE( mill_unsorted( 3, ngm_save ) )
    ALLOCATE( igsrt( ngm_max ) )
    ALLOCATE( g2l( ngm_max ) )
    ALLOCATE( g2sort_g( ngm_max ) )
    !
    g2sort_g(:) = 1.0d20
    !
    ! allocate temporal array
    !
    ALLOCATE( tt( dfftp%nr3 ) )
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfftp%nr1-1)/2
    nj = (dfftp%nr2-1)/2
    nk = (dfftp%nr3-1)/2
    !
    ! gamma-only: exclude space with x < 0
    !
    IF ( gamma_only ) THEN
       istart = 0
    ELSE
       istart = -ni
    ENDIF
    !
    iloop: DO i = istart, ni
       !
       ! gamma-only: exclude plane with x = 0, y < 0
       !
       IF ( gamma_only .and. i == 0 ) THEN
          jstart = 0
       ELSE
          jstart = -nj
       ENDIF
       !
       tx(1:3) = i * bg(1:3,1)
       !
       jloop: DO j = jstart, nj
          !
          IF ( .NOT. global_sort ) THEN
             IF ( fft_stick_index( dfftp, i, j ) == 0 ) CYCLE jloop
             is_local = .TRUE.
          ELSE
             IF ( dfftp%lpara .AND. fft_stick_index( dfftp, i, j ) == 0) THEN
                is_local = .FALSE.
             ELSE
                is_local = .TRUE.
             END IF
          END IF
          !
          ! gamma-only: exclude line with x = 0, y = 0, z < 0
          !
          IF ( gamma_only .and. i == 0 .and. j == 0 ) THEN
             kstart = 0
          ELSE
             kstart = -nk
          ENDIF
          !
          ty(1:3) = tx(1:3) + j * bg(1:3,2)
          !
          !  compute all the norm square
          !
          DO k = kstart, nk
             !
             t(1) = ty(1) + k * bg(1,3)
             t(2) = ty(2) + k * bg(2,3)
             t(3) = ty(3) + k * bg(3,3)
             tt(k-kstart+1) = t(1)**2 + t(2)**2 + t(3)**2
          ENDDO
          !
          !  save all the norm square within cutoff
          !
          DO k = kstart, nk
             IF (tt(k-kstart+1) <= gcutm) THEN
                ngm = ngm + 1
                IF (ngm > ngm_max) CALL fftx_error__ ('ggen 1', 'too many g-vectors', ngm)
                IF ( tt(k-kstart+1) > eps8 ) THEN
                   g2sort_g(ngm) = tt(k-kstart+1)
                ELSE
                   g2sort_g(ngm) = 0.d0
                ENDIF
                IF (is_local) THEN
                  ngm_local = ngm_local + 1
                  mill_unsorted( :, ngm_local ) = (/ i,j,k /)
                  g2l(ngm) = ngm_local
                ELSE
                  g2l(ngm) = 0
                ENDIF
             ENDIF
          ENDDO
       ENDDO jloop
    ENDDO iloop
    IF (ngm  /= ngm_max) &
         CALL fftx_error__ ('ggen', 'g-vectors missing !', abs(ngm - ngm_max))
    !
    igsrt(1) = 0
    IF( .NOT. global_sort ) THEN
       CALL hpsort_eps( ngm, g2sort_g, igsrt, eps8 )
    ELSE
       CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
    END IF
    DEALLOCATE( g2sort_g, tt )
    
    IF( .NOT. global_sort ) THEN
       !
       ! compute adeguate offsets in order to avoid overlap between
       ! g vectors once they are gathered on a single (global) array
       !
       !mype = mp_rank( dfftp%comm )
       !npe  = mp_size( dfftp%comm )
       CALL MPI_COMM_RANK(dfftp%comm, mype, ierr)
       CALL MPI_COMM_SIZE(dfftp%comm, npe, ierr)
       ALLOCATE( ngmpe( npe ) )
       ngmpe = 0
       ngmpe( mype + 1 ) = ngm
       !CALL mp_sum( ngmpe, dfftp%comm )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, ngmpe, 1, MPI_INTEGER, MPI_SUM, dfftp%comm, ierr )
       ngm_offset = 0
       DO ng = 1, mype
          ngm_offset = ngm_offset + ngmpe( ng )
       END DO
       DEALLOCATE( ngmpe )
       !
    END IF
    
    ngm = 0
    !
    ngloop: DO ng = 1, ngm_max
       !
       IF (g2l(igsrt(ng))>0) THEN
          ! fetch the indices
          i = mill_unsorted(1, g2l(igsrt(ng)))
          j = mill_unsorted(2, g2l(igsrt(ng)))
          k = mill_unsorted(3, g2l(igsrt(ng)))
          !
          ngm = ngm + 1
          !
          !  Here map local and global g index !!! N.B: :
          !  the global G vectors arrangement depends on the number of processors
          !
          IF( .NOT. global_sort ) THEN
             ig_l2g( ngm ) = ng + ngm_offset
          ELSE
             ig_l2g( ngm ) = ng
          END IF
       
          g(1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
          gg(ngm) = sum(g(1:3, ngm)**2)
       ENDIF
    ENDDO ngloop

    DEALLOCATE( igsrt, g2l )

    IF (ngm /= ngm_save) &
         CALL fftx_error__ ('ggen', 'g-vectors (ngm) missing !', abs(ngm - ngm_save))
    !
    !     determine first nonzero g vector
    !
    IF (gg(1).le.eps8) THEN
       gstart=2
    ELSE
       gstart=1
    ENDIF
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    CALL fft_set_nl( dfftp, at, g, mill )
    !
  END SUBROUTINE ggen
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms, &
       gs, ggs )
    !-----------------------------------------------------------------------
    !
    ! Initialize number and indices of  g-vectors for a subgrid,
    ! for exactly the same ordering as for the original FFT grid
    !
    !--------------------------------------------------------------------
    !
    USE fft_types, ONLY: fft_stick_index, fft_type_descriptor
    USE fft_ggen, ONLY : fft_set_nl
    USE fft_param
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: gamma_only
    TYPE (fft_type_descriptor), INTENT(INOUT) :: dffts
    ! primitive lattice vectors
    REAL(dp), INTENT(IN) :: at(3,3)
    ! G-vectors in FFT grid
    REAL(dp), INTENT(IN) :: g(:,:), gg(:)    
    ! Miller indices for G-vectors of FFT grid
    INTEGER, INTENT(IN) :: mill(:,:)
    ! cutoff for subgrid
    REAL(DP), INTENT(IN):: gcutms
    ! Local number of G-vectors in subgrid
    INTEGER, INTENT(OUT):: ngms
    ! Optionally: G-vectors and modules
    REAL(DP), INTENT(OUT), POINTER, OPTIONAL:: gs(:,:), ggs(:)
    !
    INTEGER :: i, ng, ngm
    !
    ngm  = SIZE(gg)
    ngms = dffts%ngm
    IF ( ngms > ngm  ) CALL fftx_error__ ('ggens','wrong  number of G-vectors',1)
    !
    IF ( PRESENT(gs) ) ALLOCATE ( gs(3,ngms) )
    IF ( PRESENT(ggs)) ALLOCATE ( ggs(ngms) )
    ng = 0
    DO i = 1, ngm
       IF ( gg(i) > gcutms ) exit
       IF ( PRESENT(gs) ) gs (:,i) = g(:,i)
       IF ( PRESENT(ggs)) ggs(i) = gg(i)
       ng = i
    END DO
    IF ( ng /= ngms ) CALL fftx_error__ ('ggens','mismatch in number of G-vectors',2)
    !
    CALL fft_set_nl ( dffts, at, g )
    !
  END SUBROUTINE ggens
  !
  SUBROUTINE fft_base_info( ionode, dffts, dfftp )
     USE fft_types, ONLY: fft_type_descriptor
     implicit none
     LOGICAL, INTENT(IN) :: ionode
     TYPE(fft_type_descriptor), INTENT(IN) :: dfftp, dffts
     !
     !  Display fft basic information
     !
     IF (ionode) THEN
        WRITE( *,*)
        IF ( dfftp%nproc > 1 ) THEN
           WRITE( stdout, '(5X,"Parallelization info")')
        ELSE
           WRITE( stdout, '(5X,"G-vector sticks info")')
        ENDIF
        WRITE( *, '(5X,"--------------------")')
        WRITE( *, '(5X,"sticks:   dense  smooth     PW", &
                       & 5X,"G-vecs:    dense   smooth      PW")') 
        IF ( dfftp%nproc > 1 ) THEN
           WRITE( *,'(5X,"Min",4X,2I8,I7,12X,2I9,I8)') &
              minval(dfftp%nsp), minval(dffts%nsp), minval(dffts%nsw), &
              minval(dfftp%ngl), minval(dffts%ngl), minval(dffts%nwl)
           WRITE( *,'(5X,"Max",4X,2I8,I7,12X,2I9,I8)') &
              maxval(dfftp%nsp), maxval(dffts%nsp), maxval(dffts%nsw), &
              maxval(dfftp%ngl), maxval(dffts%ngl), maxval(dffts%nwl)
        END IF
        WRITE( *,'(5X,"Sum",4X,2I8,I7,12X,2I9,I8)') &
           sum(dfftp%nsp), sum(dffts%nsp), sum(dffts%nsw), &
           sum(dfftp%ngl), sum(dffts%ngl), sum(dffts%nwl)
     ENDIF
     
     IF(ionode) WRITE( *,*)
     
     RETURN
  END SUBROUTINE fft_base_info
end program test

subroutine start_clock(label)
  use timers
  use mpi, ONLY:MPI_WTIME
  implicit none
  character(len=*) :: label
  if (ignore_time) RETURN
  select case (label)
  case ("cft_1z")
    times(1) = times(1) - MPI_WTIME()
  case ("cft_2xy")
    times(2) = times(2) - MPI_WTIME()
  case ("cgather")
    times(3) = times(3) - MPI_WTIME()
  case ("cgather_grid")
    times(4) = times(4) - MPI_WTIME()
  case ("cscatter_grid")
    times(5) = times(5) - MPI_WTIME()
  case ("cscatter_sym")
    times(6) = times(6) - MPI_WTIME()
  case ("fft")
    times(7) = times(7) - MPI_WTIME()
  case ("fft_scatt_tg")
    times(8) = times(8) - MPI_WTIME()
  case ("fft_scatt_xy")
    times(9) = times(9) - MPI_WTIME()
  case ("fft_scatt_yz")
    times(10) = times(10) - MPI_WTIME()
  case ("fftb")
    times(11) = times(11) - MPI_WTIME()
  case ("fftc")
    times(12) = times(12) - MPI_WTIME()
  case ("fftcw")
    times(13) = times(13) - MPI_WTIME()
  case ("ffts")
    times(14) = times(14) - MPI_WTIME()
  case ("fftw")
    times(15) = times(15) - MPI_WTIME()
  case ("rgather_grid")
    times(16) = times(16) - MPI_WTIME()
  case ("rscatter_grid")
    times(17) = times(17) - MPI_WTIME()
  case ("fft_scatter") !alt version compatibility
    times(18) = times(18) - MPI_WTIME()
  case ("ALLTOALL") !alt version compatibility
    times(19) = times(19) - MPI_WTIME()
  case default
    write (*, *) "Error, label not found", label
  end select
end subroutine

subroutine stop_clock(label)
  use timers
  use mpi, ONLY:MPI_WTIME
  implicit none
  character(len=*) :: label
  if (ignore_time) RETURN
  select case (label)
  case ("cft_1z")
    times(1) = times(1) + MPI_WTIME()
  case ("cft_2xy")
    times(2) = times(2) + MPI_WTIME()
  case ("cgather")
    times(3) = times(3) + MPI_WTIME()
  case ("cgather_grid")
    times(4) = times(4) + MPI_WTIME()
  case ("cscatter_grid")
    times(5) = times(5) + MPI_WTIME()
  case ("cscatter_sym")
    times(6) = times(6) + MPI_WTIME()
  case ("fft")
    times(7) = times(7) + MPI_WTIME()
  case ("fft_scatt_tg")
    times(8) = times(8) + MPI_WTIME()
  case ("fft_scatt_xy")
    times(9) = times(9) + MPI_WTIME()
  case ("fft_scatt_yz")
    times(10) = times(10) + MPI_WTIME()
  case ("fftb")
    times(11) = times(11) + MPI_WTIME()
  case ("fftc")
    times(12) = times(12) + MPI_WTIME()
  case ("fftcw")
    times(13) = times(13) + MPI_WTIME()
  case ("ffts")
    times(14) = times(14) + MPI_WTIME()
  case ("fftw")
    times(15) = times(15) + MPI_WTIME()
  case ("rgather_grid")
    times(16) = times(16) + MPI_WTIME()
  case ("rscatter_grid")
    times(17) = times(17) + MPI_WTIME()
  case ("fft_scatter") !alt version compatibility
    times(18) = times(18) + MPI_WTIME()
  case ("ALLTOALL") !alt version compatibility
    times(19) = times(19) + MPI_WTIME()
  case default
    write (*, *) "Error, label not found", label
  end select
end subroutine
!
subroutine print_clock(mype, npes, ncount)
  use timers
  use mpi
  implicit none
  integer, intent(in) :: mype, npes, ncount
  REAL*8  :: time_min(20)
  REAL*8  :: time_max(20)
  REAL*8  :: time_avg(20)
  integer :: ierr

#if defined(__MPI)
  CALL MPI_ALLREDUCE(times, time_min, 20, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(times, time_max, 20, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(times, time_avg, 20, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  time_avg = time_avg/npes
#else
  time_min(:) = times(:)
  time_max(:) = times(:)
  time_avg(:) = times(:)
#endif

  if (mype == 0) then
    write (*, 10100)
    write (*, 101)
    write (*, 10100)
    if (times(1) > 0.d0) write (*, 102) time_min(1), time_max(1), time_avg(1)
    if (times(2) > 0.d0) write (*, 103) time_min(2), time_max(2), time_avg(2)
    if (times(3) > 0.d0) write (*, 104) time_min(3), time_max(3), time_avg(3)
    if (times(4) > 0.d0) write (*, 105) time_min(4), time_max(4), time_avg(4)
    if (times(5) > 0.d0) write (*, 106) time_min(5), time_max(5), time_avg(5)
    if (times(6) > 0.d0) write (*, 107) time_min(6), time_max(6), time_avg(6)
    if (times(7) > 0.d0) write (*, 108) time_min(7), time_max(7), time_avg(7)
    if (times(8) > 0.d0) write (*, 109) time_min(8), time_max(8), time_avg(8)
    if (times(9) > 0.d0) write (*, 1010) time_min(9), time_max(9), time_avg(9)
    if (times(10) > 0.d0) write (*, 1011) time_min(10), time_max(10), time_avg(10)
    if (times(11) > 0.d0) write (*, 1012) time_min(11), time_max(11), time_avg(11)
    if (times(12) > 0.d0) write (*, 1013) time_min(12), time_max(12), time_avg(12)
    if (times(13) > 0.d0) write (*, 1014) time_min(13), time_max(13), time_avg(13)
    if (times(14) > 0.d0) write (*, 1015) time_min(14), time_max(14), time_avg(14)
    if (times(15) > 0.d0) write (*, 1016) time_min(15), time_max(15), time_avg(15)
    if (times(16) > 0.d0) write (*, 1017) time_min(16), time_max(16), time_avg(16)
    if (times(17) > 0.d0) write (*, 1018) time_min(17), time_max(17), time_avg(17)
    if (times(18) > 0.d0) write (*, 1019) time_min(18), time_max(18), time_avg(18)
    if (times(19) > 0.d0) write (*, 1020) time_min(19), time_max(19), time_avg(19)
    write (*, 10100)
  end if
10100 FORMAT(' +--------------------+----------------+-----------------+----------------+' )
101   FORMAT(' |FFT subroutine      |  sec. min      | sec. max        | sec.  avg      |' )
102   FORMAT(' |cft_1z              | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3,   ' |' )
103   FORMAT(' |cft_2xy             | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3,   ' |' )
104   FORMAT(' |cgather             | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
105   FORMAT(' |cgather_grid        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
106   FORMAT(' |cscatter_grid       | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
107   FORMAT(' |cscatter_sym        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
108   FORMAT(' |fft                 | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
109   FORMAT(' |fft_scatt_tg        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1010  FORMAT(' |fft_scatt_xy        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1011  FORMAT(' |fft_scatt_yz        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1012  FORMAT(' |fftb                | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1013  FORMAT(' |fftc                | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1014  FORMAT(' |fftcw               | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1015  FORMAT(' |ffts                | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1016  FORMAT(' |fftw                | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1017  FORMAT(' |rgather_grid        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1018  FORMAT(' |rscatter_grid       | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1019  FORMAT(' |fft_scatter         | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
1020  FORMAT(' |ALLTOALL            | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
 
end subroutine

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine hpsort_eps(n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  USE fft_param
  implicit none  
  !-input/output variables
  integer, intent(in) :: n  
  integer, intent(inout) :: ind (*)  
  real(DP), intent(inout) :: ra (*)
  real(DP), intent(in) :: eps
  !-local variables
  integer :: i, ir, j, l, iind  
  real(DP) :: rra  
  ! initialize index array
  if (ind (1) .eq.0) then  
     do i = 1, n  
        ind (i) = i  
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    if ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    else  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       if ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       endif
    endif
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    do while ( j .le. ir )  
       if ( j .lt. ir ) then  
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then  
             if (ra(j).lt.ra(j+1)) j = j + 1
          else
             ! this means ra(j) == ra(j+1) within tolerance
             if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
       endif
       ! demote rra
       if ( abs(rra - ra(j)).ge.eps ) then  
          if (rra.lt.ra(j)) then
             ra (i) = ra (j)  
             ind (i) = ind (j)  
             i = j  
             j = j + j  
          else
             ! set j to terminate do-while loop
             j = ir + 1  
          end if
       else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind (j) ) then
             ra (i) = ra (j)
             ind (i) = ind (j)
             i = j
             j = j + j
          else
             ! set j to terminate do-while loop
             j = ir + 1
          endif
       end if
    enddo
    ra(i) = rra
    ind(i) = iind

  end do sorting
  !
end subroutine hpsort_eps

subroutine prepare_psi_tg(ibnd, nbnd, ngms, psi, tg_psi, dffts, gamma_only)
  USE fft_param
  USE fft_types
  USE fft_helper_subroutines
  implicit none
  integer, intent(in) :: ibnd, nbnd, ngms
  TYPE(fft_type_descriptor), intent(in) :: dffts
  complex(DP) :: tg_psi(dffts%nnr_tg)
  complex(DP) :: psi(ngms, nbnd)
  logical, intent(in) :: gamma_only
  integer ioff, idx, j, ntgrp, right_nnr
!
   tg_psi(:) = ( 0.D0, 0.D0 )
   ioff   = 0
   !
   CALL tg_get_nnr( dffts, right_nnr )
   ntgrp = fftx_ntgrp(dffts)
   !
   IF (gamma_only) THEN
      DO idx = 1, 2*ntgrp, 2
         !
         ! ... 2*ntgrp ffts at the same time
         !
         IF( idx + ibnd - 1 < nbnd ) THEN
            DO j = 1, ngms
               tg_psi(dffts%nl (j)+ioff)=     psi(j,idx+ibnd-1)+&
                    (0.0d0,1.d0) * psi(j,idx+ibnd)
               tg_psi(dffts%nlm(j)+ioff)=CONJG(psi(j,idx+ibnd-1) -&
                    (0.0d0,1.d0) * psi(j,idx+ibnd) )
            END DO
         ELSE IF( idx + ibnd - 1 == nbnd ) THEN
            DO j = 1, ngms
               tg_psi(dffts%nl (j)+ioff)=       psi(j,idx+ibnd-1)
               tg_psi(dffts%nlm(j)+ioff)=CONJG( psi(j,idx+ibnd-1) )
            END DO
         END IF
      
         ioff = ioff + right_nnr
      
      END DO
   ELSE
      !
      DO idx = 1, ntgrp
         IF( idx + ibnd - 1 <= nbnd ) THEN
!$omp parallel do
            DO j = 1, ngms
                          ! here we forget about igk
               tg_psi(dffts%nl(j) +ioff) =  psi(j,idx+ibnd-1)
            ENDDO
!$omp end parallel do
         ENDIF
         ioff = ioff + right_nnr
      ENDDO
   END IF
end subroutine prepare_psi_tg

subroutine prepare_psi( ibnd, nbnd, ngms, psi, psic, dffts, gamma_only)
   USE fft_param
   USE fft_types
   USE fft_helper_subroutines
   implicit none
   integer, intent(in) :: ibnd, nbnd, ngms
   TYPE(fft_type_descriptor), intent(in) :: dffts
   complex(DP) :: psic( dffts%nnr )
   complex(DP) :: psi( ngms, nbnd )
   logical, intent(in) :: gamma_only
   integer :: j
   psic(:) = (0.d0, 0.d0)
   IF (gamma_only) THEN
     IF (ibnd < nbnd) THEN
         ! two ffts at the same time
         DO j = 1, ngms
           psic(dffts%nl (j))=      psi(j,ibnd) + (0.0d0,1.d0)*psi(j,ibnd+1)
           psic(dffts%nlm(j))=conjg(psi(j,ibnd) - (0.0d0,1.d0)*psi(j,ibnd+1))
         ENDDO
     ELSE
         DO j = 1, ngms
           psic (dffts%nl (j)) =       psi(j, ibnd)
           psic (dffts%nlm(j)) = conjg(psi(j, ibnd))
         ENDDO
     ENDIF
   ELSE
      DO j = 1, ngms
              ! here we forget about igk
         psic (dffts%nl (j)) = psi(j, ibnd)
      END DO
   END IF
end subroutine prepare_psi


subroutine accumulate_hpsi( ibnd, nbnd, ngms, hpsi, psic, dffts, gamma_only)
   USE fft_types
   USE fft_param
   USE fft_helper_subroutines
   implicit none
   integer, intent(in) :: ibnd, nbnd, ngms
   TYPE(fft_type_descriptor) :: dffts
   complex(DP) :: psic( dffts%nnr )
   complex(DP) :: hpsi( ngms, nbnd )
   logical, intent(in) :: gamma_only
   integer j
   complex(DP) :: fp, fm
   !
   !
   !   addition to the total product
   !
   IF (gamma_only) THEN
      IF (ibnd < nbnd) THEN
          ! two ffts at the same time
          DO j = 1, ngms
            fp = (psic (dffts%nl(j)) + psic (dffts%nlm(j)))*0.5d0
            fm = (psic (dffts%nl(j)) - psic (dffts%nlm(j)))*0.5d0
            hpsi (j, ibnd)   = hpsi (j, ibnd)   + &
                                cmplx( dble(fp), aimag(fm),kind=DP)
            hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + &
                                cmplx(aimag(fp),- dble(fm),kind=DP)
          ENDDO
      ELSE
          DO j = 1, ngms
            hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (dffts%nl(j))
          ENDDO
      ENDIF
    ELSE
       DO j = 1, ngms                                ! here we forget about igk_k
          hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (dffts%nl(j))
       ENDDO
    END IF
    !
end subroutine accumulate_hpsi
subroutine accumulate_hpsi_tg( ibnd, nbnd, ngms, hpsi, tg_psic, dffts, gamma_only)
   USE fft_types
   USE fft_param
   USE fft_helper_subroutines
   implicit none
   integer, intent(in) :: ibnd, nbnd, ngms
   TYPE(fft_type_descriptor) :: dffts
   complex(DP) :: tg_psic( dffts%nnr_tg  )
   complex(DP) :: hpsi( ngms, nbnd )
   logical, intent(in) :: gamma_only
   integer ioff, idx, j, right_inc
   complex(DP) :: fp, fm
   !
   !   addition to the total product
   !
   ioff   = 0
   !
   CALL tg_get_recip_inc( dffts, right_inc )
   !
   IF (gamma_only) THEN
      DO idx = 1, 2*fftx_ntgrp(dffts), 2
         !
         IF( idx + ibnd - 1 < nbnd ) THEN
            DO j = 1, ngms
               fp= ( tg_psic( dffts%nl(j) + ioff ) +  &
                     tg_psic( dffts%nlm(j) + ioff ) ) * 0.5d0
               fm= ( tg_psic( dffts%nl(j) + ioff ) -  &
                     tg_psic( dffts%nlm(j) + ioff ) ) * 0.5d0
               hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                      cmplx( dble(fp), aimag(fm),kind=DP)
               hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + &
                                      cmplx(aimag(fp),- dble(fm),kind=DP)
            ENDDO
         ELSEIF( idx + ibnd - 1 == nbnd ) THEN
            DO j = 1, ngms
               hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                       tg_psic( dffts%nl(j) + ioff )
            ENDDO
         ENDIF
         !
         ioff = ioff + right_inc
         !
      ENDDO
   ELSE
      DO idx = 1, fftx_ntgrp(dffts)
         !
         IF( idx + ibnd - 1 <= nbnd ) THEN
!$omp parallel do
            DO j = 1, ngms
               hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                  tg_psic( dffts%nl(j) + ioff )
                  ! we forgot about igk above
            ENDDO
!$omp end parallel do
         ENDIF
         !
         ioff = ioff + right_inc
         !
      ENDDO
   END IF 
   !
end subroutine accumulate_hpsi_tg
