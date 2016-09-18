!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! by F. Affinito and C. Cavazzoni, Cineca
!  & S. de Gironcoli, SISSA

program test 
  !! This mini-app provides a tool for testing and benchmarking the FFT drivers
  !! contained in the FFTXlib.
  !!
  !! The mini-app mimics a charge-density transformation from complex-to-real
  !! space and back.
  !!
  !! To compile the test program, once you have properly edit the make.sys file 
  !! included in the FFTXlib and type:
  !!
  !!      make TEST
  !!
  !! Then you can run your FFT tests using command like:
  !!
  !!      mpirun -np 4 ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4
  !!
  !! or, in case of serial build
  !!
  !!      ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4
  !!
  !! Command line arguments:
  !! 
  !!-ecutwfc  Plane wave energy cut off
  !!
  !!-alat     Lattice parameter
  !!
  !!-nbnd     Number of bands (fft cycles)
  !!
  !!-ntg      Number of task groups
  !!
  !! Timings of different stages of execution are provided at the end of the
  !! run.
  !! In the present version, a preliminar implementation with non-blocking MPI
  !! calls as been implemented. This version requires the precompilation flags
  !! -D__NON_BLOCKING_SCATTER and -D__DOUBLE_BUFFER
  !!  
  USE fft_types
  USE stick_base
  USE task_groups
  USE fft_parallel
  USE fft_support
  IMPLICIT NONE
#if defined(__MPI)
  include 'mpif.h'
  include 'fft_param.f90'
  INTEGER, ALLOCATABLE :: req_p(:),req_u(:)
#endif
  TYPE(fft_type_descriptor) :: dfftp, dffts, dfft3d
  TYPE(task_groups_descriptor) :: dtgs
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
  INTEGER :: ierr, i, ncount, ib, ireq, nreq, ipsi, iloop
  INTEGER :: ngw_ , ngm_ , ngs_
  REAL*8  :: gcutm, gkcut, gcutms
  REAL*8  :: ecutm, ecutw, ecutms
  REAL*8  :: ecutrho
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
  REAL*8  :: tmp1(10000),tmp2(10000)
  !
  LOGICAL :: gamma_only
   !! if calculations require only gamma point
  REAL*8  :: at(3,3), bg(3,3)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  !
  COMPLEX(DP), ALLOCATABLE :: psis(:,:)
   !! fake wave-function to be (anti-)transformed 
  COMPLEX(DP), ALLOCATABLE :: aux(:)
   !! fake argument returned by the FFT 
  !
  integer :: nargs
  CHARACTER(LEN=80) :: arg
  !
#if defined(__OPENMP)
  INTEGER :: PROVIDED
#endif
  !
  !   ........
  !
  !   default parameter (32 water molecules)
  !
  ecutwfc = 80.0d0
  ecutrho = 0.d0
  alat_in = 18.65
  ntgs    = 1
  nbnd    = 1
  !
  nargs = command_argument_count()
  do i = 1, nargs - 1
     CALL get_command_argument(i, arg)
     IF( TRIM( arg ) == '-ecutrho' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) ecutrho
     END IF 
     IF( TRIM( arg ) == '-ecutwfc' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) ecutwfc
     END IF 
     IF( TRIM( arg ) == '-alat' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) alat_in
     END IF 
     IF( TRIM( arg ) == '-ntg' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) ntgs
     END IF 
     IF( TRIM( arg ) == '-nbnd' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) nbnd
     END IF 
  end do
  if (ecutrho == 0.d0) ecutrho = 4.0d0 * ecutwfc
  
#if defined(__MPI)

#if defined(__OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
  CALL MPI_Init(ierr)
#endif
  CALL mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  comm = MPI_COMM_WORLD
  root = 0
  IF(mype==root) THEN
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
  CALL MPI_BCAST(ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(alat_in, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(ntgs,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(nbnd,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  !
  ecutw  = ecutwfc
  ! dual
  ecutm  = ecutrho
  ecutms = ecutrho
  !
  at(1,:) = (/ 0.5d0 , 1.0d0, 0.0d0 /)
  at(2,:) = (/ 0.5d0 , 0.0d0, 0.5d0 /)
  at(3,:) = (/ 0.0d0 , 0.5d0, 1.5d0 /)
  !
  at = at * alat_in
  !
  alat = SQRT ( at(1,1)**2+at(2,1)**2+at(3,1)**2 )
  !
  tpiba = 2.0d0 * pi / alat 
  !
  gcutm = ecutm    / tpiba**2  ! potential cut-off
  gcutms= ecutms   / tpiba**2  ! smooth mesh cut-off
  gkcut = ecutw    / tpiba**2  ! wave function cut-off
  !
  if( mype == 0 ) then

    write(*,*) '+-----------------------------------+'
    write(*,*) '|               QE FFT              |'
    write(*,*) '|          testing & timing         |'
    write(*,*) '|         by Carlo Cavazzoni        |'
    write(*,*) '+-----------------------------------+'
    write(*,*) 
    write(*,*) 'alat    = ', alat
    write(*,*) 'Ecutwfc = ', ecutw
    write(*,*) 'Ecutrho = ', ecutm
    write(*,*) 'Ecuts   = ', ecutms
    write(*,*) 'Gcutrho = ', SQRT(gcutm)
    write(*,*) 'Gcuts   = ', SQRT(gcutms)
    write(*,*) 'Gcutwfc = ', SQRT(gkcut)
    write(*,*) 'Num bands      = ', nbnd
    write(*,*) 'Num procs      = ', npes
    write(*,*) 'Num Task Group = ', ntgs
  end if
  !
  at = at / alat
  !
  call recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  !
  nx = 2 * int ( sqrt (gcutm) * sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
  ny = 2 * int ( sqrt (gcutm) * sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
  nz = 2 * int ( sqrt (gcutm) * sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
  !
  if( mype == 0 ) then
    write(*,*) 'nx = ', nx,' ny = ', ny, ' nz = ', nz
  end if
  !
  gamma_only = .true.
  CALL fft_type_init( dffts, smap, "wave", gamma_only, .true., comm, at, bg, gkcut )
  CALL fft_type_init( dfftp, smap, "rho", gamma_only, .true., comm, at, bg,  gcutm )
  CALL fft_type_init( dfft3d, smap, "wave", gamma_only, .false., comm, at, bg, gkcut)
  !
  if( mype == 0 ) then
    write(*,*) 'dffts:  nr1 = ', dffts%nr1 ,' nr2 = ', dffts%nr2 , ' nr3 = ', dffts%nr3
    write(*,*) '        nr1x= ', dffts%nr1x,' nr2x= ', dffts%nr2x, ' nr3x= ', dffts%nr3x
  end if

  CALL task_groups_init( dffts, dtgs, ntgs )
  ngw_ = dffts%nwl( dffts%mype + 1 )
  ngs_ = dffts%ngl( dffts%mype + 1 )
  ngm_ = dfftp%ngl( dfftp%mype + 1 )
  IF( gamma_only ) THEN
     ngw_ = (ngw_ + 1)/2
     ngs_ = (ngs_ + 1)/2
     ngm_ = (ngm_ + 1)/2
  END IF

  ALLOCATE( psis( dtgs%tg_nnr * dtgs%nogrp, 2 ) )
  ALLOCATE( req_p(nbnd) )
  ALLOCATE( req_u(nbnd) )
  ALLOCATE( aux( dtgs%tg_nnr * dtgs%nogrp ) )

  time = 0.0d0
  my_time = 0.0d0
  time_min = 0.0d0
  time_max = 0.0d0
  time_avg = 0.0d0

  !
  ! Test FFT for wave functions - First calls may be biased by MPI and FFT initialization
  !
  aux = 0.0d0
  aux(1) = 1.0d0
  aux(2) = 0.7d0
  aux(3) = 0.1d0

  !!write (*,*) (aux(i),i=1,5)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  CALL pack_group_sticks( aux, psis(:,1), dtgs )
  CALL fw_tg_cft3_z( psis(:,1), dffts, aux, dtgs )
  CALL fw_tg_cft3_scatter( psis(:,1), dffts, aux, dtgs )
  CALL fw_tg_cft3_xy( psis(:,1), dffts, dtgs )

  !!write (*,*) (psis(i,1),i=1,5)

  CALL bw_tg_cft3_xy( psis(:,1), dffts, dtgs )
  CALL bw_tg_cft3_scatter( psis(:,1), dffts, aux, dtgs )
  CALL bw_tg_cft3_z( psis(:,1), dffts, aux, dtgs )
  CALL unpack_group_sticks( psis(:,1), aux, dtgs )

  !!write (*,*) (aux(i),i=1,5)

  !
  ! Execute FFT calls once more and Take time
  !
  ncount = 0
  !
  wall = MPI_WTIME() 
  !
#if defined(__DOUBLE_BUFFER)
  ireq = 1
  ipsi = MOD( ireq + 1, 2 ) + 1 
  !
  CALL pack_group_sticks_i( aux, psis(:, ipsi ), dtgs, req_p( ireq ) )
  !
  nreq = 0
  DO ib = 1, nbnd, 2*dtgs%nogrp 
    nreq = nreq + 1
  END DO
  ! 
  DO ib = 1, nbnd, 2*dtgs%nogrp 
 
     ireq = ireq + 1

     aux = 0.0d0
     aux(1) = 1.0d0

     time(1) = MPI_WTIME()

     IF( ireq <= nreq ) THEN
        ipsi = MOD( ireq + 1, 2 ) + 1 
        CALL pack_group_sticks_i( aux, psis(:,ipsi), dtgs, req_p(ireq) )
     END IF

     ipsi = MOD(ipsi-1,2)+1 

     CALL MPI_WAIT( req_p( ireq - 1 ),MPI_STATUS_IGNORE)

     time(2) = MPI_WTIME()

     CALL fw_tg_cft3_z( psis( :, ipsi ), dffts, aux, dtgs )
     time(3) = MPI_WTIME()
     CALL fw_tg_cft3_scatter( psis( :, ipsi ), dffts, aux, dtgs )
     time(4) = MPI_WTIME()
     CALL fw_tg_cft3_xy( psis( :, ipsi ), dffts, dtgs )
     time(5) = MPI_WTIME()
     !
     tmp1=1.d0
     tmp2=0.d0
     !
     do iloop = 1,10 
       CALL DAXPY(10000, pi*iloop, tmp1, 1, tmp2, 1)
     end do 
     !
     time(6) = MPI_WTIME()
     CALL bw_tg_cft3_xy( psis( :, ipsi ), dffts, dtgs )
     time(7) = MPI_WTIME()
     CALL bw_tg_cft3_scatter( psis( :, ipsi ), dffts, aux, dtgs )
     time(8) = MPI_WTIME()
     CALL bw_tg_cft3_z( psis( :, ipsi ), dffts, aux, dtgs )
     time(9) = MPI_WTIME()
     !
     CALL unpack_group_sticks( psis( :, ipsi ), aux, dtgs )
     !
     time(10) = MPI_WTIME()
     !
     do i = 2, 10
        my_time(i) = my_time(i) + (time(i) - time(i-1))
     end do
     !
     ncount = ncount + 1
     !
  enddo
#else
  ipsi = 1 
  ! 
  DO ib = 1, nbnd, 2*dtgs%nogrp 
 
     aux = 0.0d0
     aux(1) = 1.0d0

     time(1) = MPI_WTIME()

     CALL pack_group_sticks( aux, psis(:,ipsi), dtgs )

     time(2) = MPI_WTIME()

     CALL fw_tg_cft3_z( psis( :, ipsi ), dffts, aux, dtgs )
     time(3) = MPI_WTIME()
     CALL fw_tg_cft3_scatter( psis( :, ipsi ), dffts, aux, dtgs )
     time(4) = MPI_WTIME()
     CALL fw_tg_cft3_xy( psis( :, ipsi ), dffts, dtgs )
     time(5) = MPI_WTIME()
     !
     tmp1=1.d0
     tmp2=0.d0
     do iloop = 1,10 
       CALL DAXPY(10000, pi*iloop, tmp1, 1, tmp2, 1)
     end do 
     !
     time(6) = MPI_WTIME()
     CALL bw_tg_cft3_xy( psis( :, ipsi ), dffts, dtgs )
     time(7) = MPI_WTIME()
     CALL bw_tg_cft3_scatter( psis( :, ipsi ), dffts, aux, dtgs )
     time(8) = MPI_WTIME()
     CALL bw_tg_cft3_z( psis( :, ipsi ), dffts, aux, dtgs )
     time(9) = MPI_WTIME()

     CALL unpack_group_sticks( psis( :, ipsi ), aux, dtgs )

     time(10) = MPI_WTIME()

     do i = 2, 10
        my_time(i) = my_time(i) + (time(i) - time(i-1))
     end do

     ncount = ncount + 1

  enddo
#endif

  wall = MPI_WTIME() - wall

  DEALLOCATE( psis, aux )

  CALL fft_type_deallocate( dffts )
  CALL fft_type_deallocate( dfftp )
  CALL fft_type_deallocate( dfft3d )
  CALL task_groups_deallocate( dtgs )

  if( ncount > 0 ) then
     my_time = my_time / DBLE(ncount)
  endif

!write(*,*)my_time(2), my_time(3), my_time(4)

#if defined(__MPI)
  CALL MPI_ALLREDUCE( my_time, time_min, 10, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( my_time, time_max, 10, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( my_time, time_avg, 10, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( wall, wall_avg,       1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
#else
  time_min = time
  time_max = time
#endif

!write(*,*)time_min(2), time_min(3), time_min(4)

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
    write(*,7) time_min(7), time_max(7), time_avg(7)
    write(*,8) time_min(8), time_max(8), time_avg(8)
    write(*,9) time_min(9), time_max(9), time_avg(9)
    write(*,10) time_min(10), time_max(10), time_avg(10)
    write(*,11) wall 
    write(*,100) 

100 FORMAT(' +--------------------+----------------+-----------------+----------------+' )
1   FORMAT(' |FFT subroutine      |  sec. min      | sec. max        | sec.  avg      |' )
2   FORMAT(' |pack_group_sticks/w | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3,   ' |' )
3   FORMAT(' |fw_tg_cft3_z        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3,   ' |' )
4   FORMAT(' |fw_tg_cft3_scatter  | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
5   FORMAT(' |fw_tg_cft3_xy       | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
6   FORMAT(' |workload            | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
7   FORMAT(' |bw_tg_cft3_xy       | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
8   FORMAT(' |bw_tg_cft3_scatter  | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
9   FORMAT(' |bw_tg_cft3_z        | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
10  FORMAT(' |unpack_group_sticks | ',    D14.5, ' | ',   D14.3,  '  | ', D14.3 ,  ' |')
11  FORMAT(' |wall time           | ',    D14.5, ' |')


  end if

  
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


end program test


subroutine start_clock( label )
implicit none
character(len=*) :: label
end subroutine

subroutine stop_clock( label )
implicit none
character(len=*) :: label
end subroutine
