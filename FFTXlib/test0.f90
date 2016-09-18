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
  USE fft_types, ONLY: fft_type_descriptor, fft_type_deallocate
  USE fft_interfaces
  USE fft_parallel
  USE fft_scalar
  USE fft_support
  IMPLICIT NONE
#if defined(__MPI)
  include 'mpif.h'
#endif
  include 'fft_param.f90'
  TYPE(fft_type_descriptor) :: dfftp, dffts, dfft3d
  INTEGER :: nx = 128
  INTEGER :: ny = 128
  INTEGER :: nz = 256
  !
  INTEGER :: mype, npes, comm, ntgs, root, nbnd
  LOGICAL :: iope
  INTEGER :: ierr, i, j, k, ncount, ib, ireq, nreq, ipsi, iloop
  INTEGER :: stdout
  INTEGER :: ngw_ , ngm_ , ngs_
  REAL*8  :: gcutm, gkcut, gcutms
  REAL*8  :: ecutm, ecutw, ecutms
  REAL*8  :: ecutwfc, ecutrho
  REAL*8  :: tpiba, alat, alat_in
  REAL*8  :: time(100), my_time(100), time_min(100), time_max(100), time_avg(100)
  REAL*8  :: wall
  REAL*8  :: wall_avg
  !
  REAL*8  :: tmp1(10000),tmp2(10000)
  !
  LOGICAL :: gamma_only
  REAL*8  :: at(3,3), bg(3,3)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  !
  COMPLEX(DP), ALLOCATABLE :: psis(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  COMPLEX(DP) :: f_aux, ff(5)
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
#if defined(__MPI)
  CALL MPI_BCAST(ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(alat_in, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(ntgs,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(nbnd,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
#endif
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
  dffts%nr1   = good_fft_order( nx )
  dffts%nr2   = good_fft_order( ny )
  dffts%nr3   = good_fft_order( nz )
  dffts%nr1x  = good_fft_dimension( dffts%nr1 )
  dffts%nr2x  = dffts%nr2
  dffts%nr3x  = good_fft_dimension( dffts%nr3 )
  !
  if( mype == 0 ) then
    write(*,*) 'dffts:  nr1 = ', dffts%nr1 ,' nr2 = ', dffts%nr2 , ' nr3 = ', dffts%nr3
    write(*,*) '        nr1x= ', dffts%nr1x,' nr2x= ', dffts%nr2x, ' nr3x= ', dffts%nr3x
  end if

  dfftp%nr1   = good_fft_order( nx )
  dfftp%nr2   = good_fft_order( ny )
  dfftp%nr3   = good_fft_order( nz )
  dfftp%nr1x  = good_fft_dimension( dfftp%nr1 )
  dfftp%nr2x  = dfftp%nr2
  dfftp%nr3x  = good_fft_dimension( dfftp%nr3 )
  !
  dfft3d%nr1   = good_fft_order( nx )
  dfft3d%nr2   = good_fft_order( ny )
  dfft3d%nr3   = good_fft_order( nz )
  dfft3d%nr1x  = good_fft_dimension( dfft3d%nr1 )
  dfft3d%nr2x  = dfft3d%nr2
  dfft3d%nr3x  = good_fft_dimension( dfft3d%nr3 )

  gamma_only = .true.
  stdout     = 6
  

  CALL pstickset( gamma_only, bg, gcutm, gkcut, gcutms, &
        dfftp, dffts, ngw_ , ngm_ , ngs_ , mype, root, &
        npes, comm, ntgs, iope, stdout, dfft3d )


!  write (6,'(25i5)') dffts%isind

  ALLOCATE( psis( dffts%tg_nnr * dffts%nogrp, 2 ) )
  ALLOCATE( aux( dffts%tg_nnr * dffts%nogrp ) )

  time = 0.0d0
  my_time = 0.0d0
  time_min = 0.0d0
  time_max = 0.0d0
  time_avg = 0.0d0

  !
  ! Test FFT for wave functions - First calls may be biased by MPI and FFT initialization
  !
  if( mype == 0 ) then
     write (*,*) 'Define a function in Reciprocal space such that it contains'
     write (*,*) '         f(1,1,1) = (1.0,0.0)    | a constant term'
     write (*,*) '         f(2,1,1) = (0.d0,0.5d0) | something varying along x '
     write (*,*) '         f(1,2,1) = (0.d0,0.3d0) | something varying along y '
     write (*,*) '         f(1,1,2) = (0.d0,0.7d0) | something varying along z '
  end if
  aux = 0.0d0
  f_aux = (1.0,0.0)   ;  call put_f_of_G(f_aux,1,1,1,aux,dffts) ! constant
  f_aux = (0.d0,0.5d0);  call put_f_of_G(f_aux,2,1,1,aux,dffts) ! something varying along x
  f_aux = (0.d0,0.3d0);  call put_f_of_G(f_aux,1,2,1,aux,dffts) ! something varying along y
  f_aux = (0.d0,0.7d0);  call put_f_of_G(f_aux,1,1,2,aux,dffts) ! something varying along z

  if( mype == 0 ) write (*,*) 'function in Reciprocal space '
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) = get_f_of_G(i,j,k,aux,dffts)
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do

#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

  call invfft ('Dense',aux,dffts)

  if( mype == 0 ) write (*,*) 'function in Real space (i,j,k)'
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) = get_f_of_R(i,j,k,aux,dffts)
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do

  call fwfft  ('Dense',aux,dffts)

  if( mype == 0 ) write (*,*) 'function in Reciprocal space '
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) = get_f_of_G(i,j,k,aux,dffts)
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do

  !
  ! Execute FFT calls once more, this time as Wave, and Take time
  !
  !
  if( mype == 0 ) write (*,*) ' Execute FFT calls once more, this time as Wave !'

#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  wall = MPI_WTIME() 
#endif

  call invfft ('Wave',aux,dffts)

  if( mype == 0 ) write (*,*) 'function in Real space (i,j,k)'
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) = get_f_of_R(i,j,k,aux,dffts)
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do

  call fwfft  ('Wave',aux,dffts)

  if( mype == 0 ) write (*,*) 'function in Reciprocal space '
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) = get_f_of_G(i,j,k,aux,dffts)
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do

  DEALLOCATE( aux )

  if( mype == 0 ) write (*,*) ' Execute FFT calls once more, this time as with cft3ds !'

  ALLOCATE( aux (dffts%nr1x * dffts%nr2x * dffts%nr3x ) )

#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  wall = MPI_WTIME() 
#endif

  aux = 0.0d0
  f_aux = (1.0,0.0)   ; i=1;j=1;k=1;  aux( i+dffts%nr1x *(j-1) + dffts%nr1x*dffts%nr2x*(k-1) ) = f_aux
  f_aux = (0.d0,0.5d0); i=2;j=1;k=1;  aux( i+dffts%nr1x *(j-1) + dffts%nr1x*dffts%nr2x*(k-1) ) = f_aux
  f_aux = (0.d0,0.3d0); i=1;j=2;k=1;  aux( i+dffts%nr1x *(j-1) + dffts%nr1x*dffts%nr2x*(k-1) ) = f_aux 
  f_aux = (0.d0,0.7d0); i=1;j=1;k=2;  aux( i+dffts%nr1x *(j-1) + dffts%nr1x*dffts%nr2x*(k-1) ) = f_aux 
  
  CALL cfft3ds( aux, dffts%nr1, dffts%nr2, dffts%nr3, &
                   dffts%nr1x,dffts%nr2x,dffts%nr3x, 1, &
                   dffts%isind, dffts%iplw )

  if( mype == 0 ) write (*,*) 'function in Real space (i,j,k)'
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) =aux (i+dffts%nr1x *(j-1) + dffts%nr1x*dffts%nr2x*(k-1) )
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do

  CALL cfft3ds( aux, dffts%nr1, dffts%nr2, dffts%nr3, &
                   dffts%nr1x,dffts%nr2x,dffts%nr3x, -1, &
                   dffts%isind, dffts%iplw )

  if( mype == 0 ) write (*,*) 'function in Reciprocal space '
  do k =1, 5
     if( mype == 0 ) write (*,*) 'k = ',k
     do j =1,5
        do i =1,5
           ff(i) =aux (i+dffts%nr1x *(j-1) + dffts%nr1x*dffts%nr2x*(k-1) )
        end do
        if( mype == 0 ) write (*,'(5("(",2f10.6,")",3x))') (ff(i),i=1,5)
     end do
  end do
#if defined(__MPI)
  wall = MPI_WTIME() - wall
#endif

  DEALLOCATE( psis, aux )

  CALL fft_type_deallocate( dffts )
  CALL fft_type_deallocate( dfftp )
  CALL fft_type_deallocate( dfft3d )

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
