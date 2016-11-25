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
  INTEGER :: ierr, i, j, ncount, ib, ireq, nreq, ipsi, iloop
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
  LOGICAL :: gamma_only
   !! if calculations require only gamma point
  REAL*8  :: at(3,3), bg(3,3)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  !
  COMPLEX(DP), ALLOCATABLE :: psis(:,:)
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
   !! fake wave-function to be (anti-)transformed 
  COMPLEX(DP), ALLOCATABLE :: aux(:)
   !! fake argument returned by the FFT 
  COMPLEX(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)
   !! array representing the potential
  INTEGER, ALLOCATABLE :: nls( : ), nlsm( : )
  INTEGER :: ngms, ngsx, ngms_g
  INTEGER, ALLOCATABLE :: mill(:,:), nl(:), nlm(:), ig_l2g(:)
  COMPLEX(DP), ALLOCATABLE :: g(:,:), gg(:)
  INTEGER :: ngm, ngmx, ngm_g, gstart
  !
  INTEGER :: ibnd, ioff, idx, incr, m
  COMPLEX(DP) :: fp, fm

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
  !
  ! --------  INITIALIZE DIMENSIONS AND DESCRIPTORS
  !
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
  CALL fft_type_init( dffts, smap, "wave", gamma_only, .true., comm, at, bg, gkcut, gcutms/gkcut, ntgs )
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

  ngms = ngs_
  CALL MPI_ALLREDUCE( ngms, ngsx, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( ngms, ngms_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  ngm = ngm_
  CALL MPI_ALLREDUCE( ngm, ngmx, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( ngm, ngm_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )


  ! --------  ALLOCATE


  ALLOCATE( psis( dtgs%tg_nnr * dtgs%nogrp, 2 ) )
  ALLOCATE( psi( ngms, nbnd ) )
  ALLOCATE( req_p(nbnd) )
  ALLOCATE( req_u(nbnd) )
  ALLOCATE( aux( dtgs%tg_nnr * dtgs%nogrp ) )
  ALLOCATE( tg_v( dtgs%tg_nnr * dtgs%nogrp ) )
  ALLOCATE( hpsi( nbnd, dtgs%tg_nnr * dtgs%nogrp ) )
  ALLOCATE( nls( ngms ) )
  ALLOCATE( nlsm( ngms ) )
  ALLOCATE( nl( ngm ) )
  ALLOCATE( nlm( ngm ) )
  ALLOCATE( mill( 3, ngm ) )
  ALLOCATE( g( 3, ngm ) )
  ALLOCATE( gg( ngm ) )
  ALLOCATE( ig_l2g( ngm ) )


  ! --------  GENERATE G-VECTORS


  call ggen ( gamma_only, at, bg, .true., ngm, ngms, ngm_g, ngms_g, mill, &
&                    nl, nls, nlm, nlsm, gg, g, ig_l2g, gstart, gcutm, gcutms, dfftp, dffts )


  ! --------  RESET TIMERS

  time = 0.0d0
  my_time = 0.0d0
  time_min = 0.0d0
  time_max = 0.0d0
  time_avg = 0.0d0

  !
  ! Test FFT for wave functions - First calls may be biased by MPI and FFT initialization
  !
  aux = 0.0d0
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  CALL pack_group_sticks( aux, psis(:,1), dtgs )
  CALL fw_tg_cft3_z( psis(:,1), dffts, aux, dtgs )
  CALL fw_tg_cft3_scatter( psis(:,1), dffts, aux, dtgs )
  CALL fw_tg_cft3_xy( psis(:,1), dffts, dtgs )
  CALL bw_tg_cft3_xy( psis(:,1), dffts, dtgs )
  CALL bw_tg_cft3_scatter( psis(:,1), dffts, aux, dtgs )
  CALL bw_tg_cft3_z( psis(:,1), dffts, aux, dtgs )
  CALL unpack_group_sticks( psis(:,1), aux, dtgs )

  !
  ! --------  INITIALIZE WAVE FUNCTIONS psi

  psi = 0.0d0

  ! --------  INITIALIZE POTENTIAL tg_v

  tg_v = 1.0d0

  ! --------  INITIALIZE FORCE hpsi

  hpsi = 0.0d0

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
  CALL prepare_psi( 1, nbnd, ngms, psi, aux, nls, nlsm, dtgs)
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
        call  prepare_psi( ib+1, nbnd, ngms, psi, aux, nls, nlsm, dtgs)
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
     DO j = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( dffts%mype + 1 )
        psis (j,ipsi) = psis (j,ipsi) * tg_v(j)
     ENDDO
     !
     time(6) = MPI_WTIME()
     CALL bw_tg_cft3_xy( psis( :, ipsi ), dffts, dtgs )
     time(7) = MPI_WTIME()
     CALL bw_tg_cft3_scatter( psis( :, ipsi ), dffts, aux, dtgs )
     time(8) = MPI_WTIME()
     CALL bw_tg_cft3_z( psis( :, ipsi ), dffts, aux, dtgs )
     time(9) = MPI_WTIME()
     !
     IF(ireq == 2)THEN
       CALL unpack_group_sticks_i( psis( :, ipsi ), aux, dtgs , req_u(ireq) )
     ELSE
       CALL MPI_WAIT(req_u(ireq-1), MPI_STATUS_IGNORE, ierr )
       ipsi = MOD( ireq + 1, 2 ) + 1 ! ireq = 2, ipsi = 2; ireq = 3, ipsi = 1
       CALL unpack_group_sticks_i( psis( :, ipsi ), aux, dtgs, req_u(ireq) )
     ENDIF

     call accumulate_hpsi( ib, nbnd, ngms, hpsi, aux, nls, nlsm, dtgs, dffts)

     time(10) = MPI_WTIME()
     !
     do i = 2, 10
        my_time(i) = my_time(i) + (time(i) - time(i-1))
     end do
     !
     ncount = ncount + 1
     !
  enddo
  CALL MPI_WAIT(req_u(ireq),MPI_STATUS_IGNORE,ierr)

#else

  ipsi = 1 
  ! 
  DO ib = 1, nbnd, 2*dtgs%nogrp 
 
     aux = 0.0d0
     aux(1) = 1.0d0

     call  prepare_psi( ib, nbnd, ngms, psi, aux, nls, nlsm, dtgs)

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
     DO j = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( dffts%mype + 1 )
        psis (j,ipsi) = psis (j,ipsi) * tg_v(j)
     ENDDO
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
     !
     do i = 2, 10
        my_time(i) = my_time(i) + (time(i) - time(i-1))
     end do
     !
     ncount = ncount + 1
     !
  enddo
#endif

  wall = MPI_WTIME() - wall

  DEALLOCATE( psis, aux )


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

  CALL fft_type_deallocate( dffts )
  CALL fft_type_deallocate( dfftp )
  CALL fft_type_deallocate( dfft3d )
  CALL task_groups_deallocate( dtgs )
  
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
!
   !-----------------------------------------------------------------------
   SUBROUTINE ggen ( gamma_only, at, bg, no_global_sort, ngm, ngms, ngm_g, ngms_g, mill, &
&                    nl, nls, nlm, nlsm, gg, g, ig_l2g, gstart, gcutm, gcutms, dfftp, dffts )
   !----------------------------------------------------------------------
   !
   !     This routine generates all the reciprocal lattice vectors
   !     contained in the sphere of radius gcutm. Furthermore it
   !     computes the indices nl which give the correspondence
   !     between the fft mesh points and the array of g vectors.
   !
   USE fft_types
   USE fft_param
   !
   IMPLICIT NONE
   !
   REAL(DP), PARAMETER :: eps8  = 1.0E-8_DP
   !
   LOGICAL,  INTENT(IN) :: gamma_only
   REAL(DP), INTENT(IN) :: at(3,3), bg(3,3), gcutm, gcutms
   LOGICAL,  INTENT(IN) :: no_global_sort
   !  if no_global_sort is present (and it is true) G vectors are sorted only
   !  locally and not globally. In this case no global array needs to be
   !  allocated and sorted: saves memory and a lot of time for large systems.
   INTEGER :: ngm, ngms, ngm_g, ngms_g, gstart
   INTEGER :: mill(3,ngm), nlm(ngm), nlsm(ngms), ig_l2g(ngm_g)
   INTEGER :: nl(ngm), nls(ngms)
   REAL(DP) :: gg(ngm), g(3,ngm)
   TYPE(fft_type_descriptor) :: dfftp, dffts
   !
   !     here a few local variables
   !
   REAL(DP) ::  t (3), tt
   INTEGER :: ngm_save, ngms_save, n1, n2, n3, n1s, n2s, n3s, ngm_offset, ngm_max, ngms_max
   INTEGER :: ierr
   !
   REAL(DP), ALLOCATABLE :: g2sort_g(:)
   ! array containing all g vectors, on all processors: replicated data
   ! when no_global_sort is present (and it is true) only g vectors for the current processor are stored
   INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
   ! array containing all g vectors generators, on all processors: replicated data
   ! when no_global_sort is present (and it is true) only g vectors for the current processor are stored
   INTEGER, ALLOCATABLE :: igsrt(:)
   !
   INTEGER :: m1, m2, mc
   INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw
   INTEGER :: mype, npe, comm
   LOGICAL :: global_sort
   INTEGER, ALLOCATABLE :: ngmpe(:)
   !
   global_sort = .NOT. no_global_sort
   !
   comm = dfftp%comm
   mype = dfftp%mype
   npe = dfftp%nproc
   !
   IF( .NOT. global_sort ) THEN
      ALLOCATE( ngmpe( npe ) )
      ngmpe = 0
      ngm_max = ngm
      ngms_max = ngms
   ELSE
      ngm_max = ngm_g
      ngms_max = ngms_g
   END IF
   !
   ! save current value of ngm and ngms
   !
   ngm_save  = ngm
   ngms_save = ngms
   !
   ngm = 0
   ngms = 0
   !
   ! counters
   !
   !    set the total number of fft mesh points and and initial value of gg
   !    The choice of gcutm is due to the fact that we have to order the
   !    vectors after computing them.
   !
   gg(:) = gcutm + 1.d0
   !
   !    and computes all the g vectors inside a sphere
   !
   ALLOCATE( mill_g( 3, ngm_max ),mill_unsorted( 3, ngm_max ) )
   ALLOCATE( igsrt( ngm_max ) )
   ALLOCATE( g2sort_g( ngm_max ) )
   !
   g2sort_g(:) = 1.0d20
   !
   ! max miller indices (same convention as in module stick_set)
   !
   ni = (dfftp%nr1-1)/2
   nj = (dfftp%nr2-1)/2
   nk = (dfftp%nr3-1)/2
   !
   iloop: DO i = -ni, ni
      !
      ! gamma-only: exclude space with x < 0
      !
      IF ( gamma_only .and. i < 0) CYCLE iloop
      jloop: DO j = -nj, nj
         !
         ! gamma-only: exclude plane with x = 0, y < 0
         !
         IF ( gamma_only .and. i == 0 .and. j < 0) CYCLE jloop

         IF( .NOT. global_sort ) THEN
            m1 = mod (i, dfftp%nr1) + 1
            IF (m1 < 1) m1 = m1 + dfftp%nr1
            m2 = mod (j, dfftp%nr2) + 1
            IF (m2 < 1) m2 = m2 + dfftp%nr2
            mc = m1 + (m2 - 1) * dfftp%nr1x
            IF ( dfftp%isind ( mc ) == 0) CYCLE jloop
         END IF

         kloop: DO k = -nk, nk
            !
            ! gamma-only: exclude line with x = 0, y = 0, z < 0
            !
            IF ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) CYCLE kloop
            t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
            !tt = sum(t(:)**2)
            tt = t(1)**2+t(2)**2+t(3)**2
            IF (tt <= gcutm) THEN
               ngm = ngm + 1
               IF (tt <= gcutms) ngms = ngms + 1
               IF (ngm > ngm_max) CALL fftx_error__ ('ggen 1', 'too many g-vectors', ngm)
               mill_unsorted( :, ngm ) = (/ i,j,k /)
               IF ( tt > eps8 ) THEN
                  g2sort_g(ngm) = tt
               ELSE
                  g2sort_g(ngm) = 0.d0
               ENDIF
            ENDIF
         ENDDO kloop
      ENDDO jloop
   ENDDO iloop

   IF( .NOT. global_sort ) THEN
      ngmpe( mype + 1 ) = ngm
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, ngmpe, 1, MPI_INTEGER, MPI_SUM, comm, ierr )
   END IF
   IF (ngm  /= ngm_max) &
         CALL fftx_error__ ('ggen', 'g-vectors missing !', abs(ngm - ngm_max))
   IF (ngms /= ngms_max) &
         CALL fftx_error__ ('ggen', 'smooth g-vectors missing !', abs(ngms - ngms_max))

   igsrt(1) = 0
   IF( .NOT. global_sort ) THEN
      CALL hpsort_eps( ngm, g2sort_g, igsrt, eps8 )
   ELSE
      CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
   END IF
   mill_g(1,:) = mill_unsorted(1,igsrt(:))
   mill_g(2,:) = mill_unsorted(2,igsrt(:))
   mill_g(3,:) = mill_unsorted(3,igsrt(:))
   DEALLOCATE( g2sort_g, igsrt, mill_unsorted )

   IF( .NOT. global_sort ) THEN
      ! compute adeguate offsets in order to avoid overlap between
      ! g vectors once they are gathered on a single (global) array
      !
      ngm_offset = 0
      DO ng = 1, mype
         ngm_offset = ngm_offset + ngmpe( ng )
      END DO
   END IF

   ngm = 0
   ngms = 0
   !
   ngloop: DO ng = 1, ngm_max

      i = mill_g(1, ng)
      j = mill_g(2, ng)
      k = mill_g(3, ng)

#if defined(__MPI)
      IF( global_sort ) THEN
         m1 = mod (i, dfftp%nr1) + 1
         IF (m1 < 1) m1 = m1 + dfftp%nr1
         m2 = mod (j, dfftp%nr2) + 1
         IF (m2 < 1) m2 = m2 + dfftp%nr2
         mc = m1 + (m2 - 1) * dfftp%nr1x
         IF ( dfftp%isind ( mc ) == 0) CYCLE ngloop
      END IF
#endif

      ngm = ngm + 1

      !  Here map local and global g index !!!
      !  N.B. the global G vectors arrangement depends on the number of processors
      !
      IF( .NOT. global_sort ) THEN
         ig_l2g( ngm ) = ng + ngm_offset
      ELSE
         ig_l2g( ngm ) = ng
      END IF

      g (1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
      gg (ngm) = sum(g (1:3, ngm)**2)

      IF (gg (ngm) <= gcutms) ngms = ngms + 1
      IF (ngm > ngm_save) CALL fftx_error__ ('ggen 2', 'too many g-vectors', ngm)
   ENDDO ngloop

   IF (ngm /= ngm_save) &
      CALL fftx_error__ ('ggen', 'g-vectors (ngm) missing !', abs(ngm - ngm_save))
   IF (ngms /= ngms_save) &
      CALL fftx_error__ ('ggen', 'g-vectors (ngms) missing !', abs(ngm - ngms_save))
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
   DO ng = 1, ngm
      n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      mill (1,ng) = n1 - 1
      n1s = n1
      IF (n1<1) n1 = n1 + dfftp%nr1
      IF (n1s<1) n1s = n1s + dffts%nr1

      n2 = nint (sum(g (:, ng) * at (:, 2))) + 1
      mill (2,ng) = n2 - 1
      n2s = n2
      IF (n2<1) n2 = n2 + dfftp%nr2
      IF (n2s<1) n2s = n2s + dffts%nr2

      n3 = nint (sum(g (:, ng) * at (:, 3))) + 1
      mill (3,ng) = n3 - 1
      n3s = n3
      IF (n3<1) n3 = n3 + dfftp%nr3
      IF (n3s<1) n3s = n3s + dffts%nr3

      IF (n1>dfftp%nr1 .or. n2>dfftp%nr2 .or. n3>dfftp%nr3) &
         CALL fftx_error__('ggen','Mesh too small?',ng)

#if defined (__MPI) && !defined (__USE_3D_FFT)
      nl (ng) = n3 + ( dfftp%isind (n1 + (n2 - 1) * dfftp%nr1x) - 1) * dfftp%nr3x
      IF (ng <= ngms) &
         nls (ng) = n3s + ( dffts%isind (n1s+(n2s-1)*dffts%nr1x) - 1 ) * dffts%nr3x
#else
      nl (ng) = n1 + (n2 - 1) * dfftp%nr1x + (n3 - 1) * dfftp%nr1x * dfftp%nr2x
      IF (ng <= ngms) &
         nls (ng) = n1s + (n2s - 1) * dffts%nr1x + (n3s - 1) * dffts%nr1x * dffts%nr2x
#endif
   ENDDO
   !
   DEALLOCATE( mill_g )

   IF ( gamma_only) CALL index_minusg( ngm, ngms, nlm, nlsm, mill, dfftp, dffts )

   IF( ALLOCATED( ngmpe ) ) DEALLOCATE( ngmpe )

   END SUBROUTINE ggen

   !
   !-----------------------------------------------------------------------
   SUBROUTINE index_minusg( ngm, ngms, nlm, nlsm, mill, dfftp, dffts )
   !----------------------------------------------------------------------
   !
   !     compute indices nlm and nlms giving the correspondence
   !     between the fft mesh points and -G (for gamma-only calculations)
   !
   USE fft_types
   !
   IMPLICIT NONE
   !
   TYPE(fft_type_descriptor) :: dfftp, dffts
   INTEGER :: ngm, ngms
   INTEGER :: mill(3,ngm), nlm(ngm), nlsm(ngms)
   !
   INTEGER :: n1, n2, n3, n1s, n2s, n3s, ng
   !
   DO ng = 1, ngm
      n1 = -mill (1,ng) + 1
      n1s = n1
      IF (n1 < 1) THEN
         n1 = n1 + dfftp%nr1
         n1s = n1s + dffts%nr1
      END IF

      n2 = -mill (2,ng) + 1
      n2s = n2
      IF (n2 < 1) THEN
         n2 = n2 + dfftp%nr2
         n2s = n2s + dffts%nr2
      END IF
      n3 = -mill (3,ng) + 1
      n3s = n3
      IF (n3 < 1) THEN
         n3 = n3 + dfftp%nr3
         n3s = n3s + dffts%nr3
      END IF

      IF (n1>dfftp%nr1 .or. n2>dfftp%nr2 .or. n3>dfftp%nr3) THEN
         CALL fftx_error__('index_minusg','Mesh too small?',ng)
      ENDIF

#if defined (__MPI) && !defined (__USE_3D_FFT)
      nlm(ng) = n3 + (dfftp%isind (n1 + (n2 - 1) * dfftp%nr1x) - 1) * dfftp%nr3x
      IF (ng<=ngms) &
         nlsm(ng) = n3s + (dffts%isind (n1s+(n2s-1) * dffts%nr1x) - 1) * dffts%nr3x
#else
      nlm(ng) = n1 + (n2 - 1) * dfftp%nr1x + (n3 - 1) * dfftp%nr1x * dfftp%nr2x
      IF (ng<=ngms) &
         nlsm(ng) = n1s + (n2s - 1) * dffts%nr1x + (n3s-1) * dffts%nr1x * dffts%nr2x
#endif
   ENDDO

   END SUBROUTINE index_minusg
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine hpsort_eps (n, ra, ind, eps)
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
    ra (i) = rra  
    ind (i) = iind  

  end do sorting    
  !
end subroutine hpsort_eps


subroutine prepare_psi( ib, nbnd, ngms, psi, tg_psic, nls, nlsm, dtgs)
   USE task_groups
   USE fft_param
   implicit none
   integer, intent(in) :: ib, nbnd, ngms
   TYPE(task_groups_descriptor), intent(in) :: dtgs
   complex(DP) :: tg_psic( dtgs%tg_nnr * dtgs%nogrp )
   complex(DP) :: psi( ngms, nbnd )
   integer, intent(in) :: nls(ngms), nlsm(ngms)
   integer ioff, n, ibnd, m, idx, j
   !
   ibnd = ib
   m = nbnd
   n = ngms

        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        DO idx = 1, 2*dtgs%nogrp, 2
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 tg_psic(nls (j)+ioff) =        psi(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd)
                 tg_psic(nlsm(j)+ioff) = conjg( psi(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 tg_psic(nls (j)+ioff) =        psi(j,idx+ibnd-1)
                 tg_psic(nlsm(j)+ioff) = conjg( psi(j,idx+ibnd-1) )
              ENDDO
           ENDIF
           ioff = ioff + dtgs%tg_nnr
        END DO

   return
end subroutine prepare_psi

subroutine accumulate_hpsi( ib, nbnd, ngms, hpsi, tg_psic, nls, nlsm, dtgs, dffts)
   USE task_groups
   USE fft_types
   USE fft_param
   implicit none
   integer, intent(in) :: ib, nbnd, ngms
   integer, intent(in) :: nls(ngms), nlsm(ngms)
   TYPE(task_groups_descriptor), intent(in) :: dtgs
   TYPE(fft_type_descriptor) :: dffts
   complex(DP) :: tg_psic( dtgs%tg_nnr * dtgs%nogrp )
   complex(DP) :: hpsi( ngms, nbnd )
   integer ioff, n, ibnd, m, idx, j
   complex(DP) :: fp, fm
   !
   ibnd = ib
   m = nbnd
   n = ngms
        !
        ioff   = 0
        !
        DO idx = 1, 2*dtgs%nogrp, 2
           !
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 fp= ( tg_psic( nls(j) + ioff ) +  &
                       tg_psic( nlsm(j) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(j) + ioff ) -  &
                       tg_psic( nlsm(j) + ioff ) ) * 0.5d0
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                        cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + &
                                        cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                         tg_psic( nls(j) + ioff )
              ENDDO
           ENDIF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( dffts%mype + 1 )
           !
        ENDDO
        !
end subroutine accumulate_hpsi
