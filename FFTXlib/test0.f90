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
  USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_deallocate
  USE stick_set, ONLY: pstickset
  USE fft_interfaces
  USE fft_parallel
  USE fft_support
  IMPLICIT NONE
#ifdef __MPI
  include 'mpif.h'
  include 'fft_param.f90'
  INTEGER, ALLOCATABLE :: req_p(:),req_u(:)
#endif
  TYPE(fft_dlay_descriptor) :: dfftp, dffts, dfft3d
  INTEGER :: nx = 128
  INTEGER :: ny = 128
  INTEGER :: nz = 256
  !
  INTEGER :: mype, npes, comm, ntgs, root, nbnd
  LOGICAL :: iope
  INTEGER :: ierr, i, ncount, ib, ireq, nreq, ipsi, iloop
  INTEGER :: stdout
  INTEGER :: ngw_ , ngm_ , ngs_
  REAL*8  :: gcutm, gkcut, gcutms
  REAL*8  :: ecutm, ecutw, ecutms
  REAL*8  :: ecutrho
  REAL*8  :: ecutwfc
  REAL*8  :: tpiba, alat, alat_in
  REAL*8  :: tempo(100)
  REAL*8  :: tempo_mio(100)
  REAL*8  :: tempo_min(100)
  REAL*8  :: tempo_max(100)
  REAL*8  :: tempo_avg(100)
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
  COMPLEX(DP) :: f_aux
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
  
#ifdef __MPI

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

  ALLOCATE( psis( dffts%tg_nnr * dffts%nogrp, 2 ) )
  ALLOCATE( req_p(nbnd) )
  ALLOCATE( req_u(nbnd) )
  ALLOCATE( aux( dffts%tg_nnr * dffts%nogrp ) )

  tempo = 0.0d0
  tempo_mio = 0.0d0
  tempo_min = 0.0d0
  tempo_max = 0.0d0
  tempo_avg = 0.0d0

  !
  ! Test FFT for wave functions - First calls may be biased by MPI and FFT initialization
  !
  aux = 0.0d0
  f_aux = (1.0,0.0)
  call put_f_of_G(f_aux,1,1,1,aux,dffts)
  f_aux = (0.d0,0.7d0)
  call put_f_of_G(f_aux,1,1,2,aux,dffts)
  f_aux = (0.d0,0.5d0)
  call put_f_of_G(f_aux,2,1,1,aux,dffts)

  write (*,*) (aux(i),i=1,5)

  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

  call invfft ('Dense',aux,dffts)
  write (*,*) (aux(i),i=1,5)
  write (3,*) aux
  write (*,*) get_f_of_R(1,1,1,aux,dffts),get_f_of_R(2,1,1,aux,dffts),get_f_of_R(3,1,1,aux,dffts),get_f_of_R(4,1,1,aux,dffts)
  call fwfft  ('Dense',aux,dffts)
  write (*,*) (aux(i),i=1,5)
  write (4,*) aux
  write (*,*) get_f_of_G(1,1,1,aux,dffts),get_f_of_G(2,1,1,aux,dffts),get_f_of_G(3,1,1,aux,dffts),get_f_of_G(4,1,1,aux,dffts)
  !
  ! Execute FFT calls once more and Take time
  !
  !

  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  wall = MPI_WTIME() 
  call invfft ('Dense',aux,dffts)
  call fwfft  ('Dense',aux,dffts)
  wall = MPI_WTIME() - wall

  DEALLOCATE( psis, aux )

  CALL fft_dlay_deallocate( dffts )
  CALL fft_dlay_deallocate( dfftp )
  CALL fft_dlay_deallocate( dfft3d )

  if( ncount > 0 ) then
     tempo_mio = tempo_mio / DBLE(ncount)
  endif

!write(*,*)tempo_mio(2), tempo_mio(3), tempo_mio(4)

#ifdef __MPI
  CALL MPI_ALLREDUCE( tempo_mio, tempo_min, 10, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( tempo_mio, tempo_max, 10, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( tempo_mio, tempo_avg, 10, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( wall, wall_avg,       1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
#else
  tempo_min = tempo
  tempo_max = tempo
#endif

!write(*,*)tempo_min(2), tempo_min(3), tempo_min(4)

  tempo_avg = tempo_avg / npes
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
    write(*,2) tempo_min(2), tempo_max(2), tempo_avg(2)
    write(*,3) tempo_min(3), tempo_max(3), tempo_avg(3)
    write(*,4) tempo_min(4), tempo_max(4), tempo_avg(4)
    write(*,5) tempo_min(5), tempo_max(5), tempo_avg(5)
    write(*,6) tempo_min(6), tempo_max(6), tempo_avg(6)
    write(*,7) tempo_min(7), tempo_max(7), tempo_avg(7)
    write(*,8) tempo_min(8), tempo_max(8), tempo_avg(8)
    write(*,9) tempo_min(9), tempo_max(9), tempo_avg(9)
    write(*,10) tempo_min(10), tempo_max(10), tempo_avg(10)
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

  
#ifdef __MPI
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
