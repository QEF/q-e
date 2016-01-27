program test
  USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_deallocate
  USE stick_set, ONLY: pstickset
  IMPLICIT NONE
#ifdef __MPI
  include 'mpif.h'
#endif
  TYPE(fft_dlay_descriptor) :: dfftp, dffts, dfft3d 
  INTEGER, PARAMETER :: nx = 128
  INTEGER, PARAMETER :: ny = 128
  INTEGER, PARAMETER :: nz = 256
  !
  INTEGER :: mype, npes, comm, ntgs, root
  LOGICAL :: iope
  INTEGER :: ierr
  INTEGER :: stdout
  INTEGER :: ngw_ , ngm_ , ngs_
  REAL*8  :: bg(3,3)
  REAL*8  :: gcutm, gkcut, gcutms
  LOGICAL :: gamma_only
  !
#if defined(__OPENMP)
  INTEGER :: PROVIDED
#endif
  !
  !   ........
  !
#ifdef __MPI
#if defined(__OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
  CALL MPI_Init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  comm = MPI_COMM_WORLD
  ntgs = 1
  root = 0
  IF(mype==root) THEN
     iope = .true.
  ELSE
     iope = .false.
  ENDIF
#endif
#else
  mype = 0
  npes = 1
  comm = 0
  ntgs = 1
  root = 0
  iope = .true.
#endif
  !
  dffts%nr1 = nx 
  dffts%nr2 = ny
  dffts%nr3 = nz
  !
  dffts%nr1x = nx 
  dffts%nr2x = ny
  dffts%nr3x = nz
  !
  dfftp%nr1 = nx 
  dfftp%nr2 = ny
  dfftp%nr3 = nz
  !
  dfftp%nr1x = nx 
  dfftp%nr2x = ny
  dfftp%nr3x = nz
  !
  dfft3d%nr1 = nx 
  dfft3d%nr2 = ny
  dfft3d%nr3 = nz
  !
  dfft3d%nr1x = nx 
  dfft3d%nr2x = ny
  dfft3d%nr3x = nz

  gamma_only = .true.
  stdout     = 6

  bg      = 0.0d0
  bg(1,1) = 1.0d0
  bg(2,2) = 1.0d0
  bg(3,3) = 0.5d0

  gcutm  = 53.1
  gcutms = 53.1
  gkcut  = 26.5

  CALL pstickset( gamma_only, bg, gcutm, gkcut, gcutms, &
        dfftp, dffts, ngw_ , ngm_ , ngs_ , mype, root, &
        npes, comm, ntgs, iope, stdout, dfft3d )

  !

  CALL fft_dlay_deallocate( dffts )
  CALL fft_dlay_deallocate( dfftp )
  CALL fft_dlay_deallocate( dfft3d )
  
  call fftx_error__( ' test ', ' FFTXlib test program ', 1 )
#ifdef __MPI
  CALL mpi_finalize(ierr)
#endif
end program test
