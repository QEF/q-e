!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
Module para_const
  Integer, parameter  :: maxproc  = 128  
  ! maximum number of processors
End Module para_const
!
Module pfft
  !
  ! parallel fft information for the dense grid
  !
  Use para_const
  Integer :: npp (maxproc), ncp (maxproc), ncp0 (maxproc), &
       ncplane, nct, nxx
  ! npp    : number of plane per processor
  ! ncp    : number of (density) columns per proc
  ! ncp0   : starting column for each processor
  ! ncplane: number of columns in a plane
  ! nct    : total number of non-zero columns
  ! nxx    : local fft data size
  Integer, Pointer :: ipc(:), icpl(:)
  ! ipc (ncplane): index saying which proc owns columns in a plane
  ! icpl(nct)    : index relating columns and pos. in the plan
End Module pfft
!
Module pffts
  !
  ! parallel fft information for the smooth grid
  !
  Use para_const
  Integer :: nkcp (maxproc)
  ! nkcp   : number of (wfs) columns per processor
  Integer :: npps (maxproc), ncps (maxproc), ncp0s (maxproc), &
       ncplanes, ncts, nxxs 
  ! npps   : number of plane per processor
  ! ncps   : number of (density) columns per proc
  ! ncp0s  : starting column for each processor
  ! ncplanes:number of columns in a plane
  ! ncts   : total number of non-zero columns
  ! nxxs   : local fft data size
  Integer , Pointer :: ipcs(:), icpls(:)
  ! ipcs(ncplanes) : index saying which proc owns columns in a plane
  ! icpls(ncts)    : index relating columns and positions in the plan
End Module pffts
!
Module para
  Use pfft
  Use pffts
  Use mp_global, only: nproc
  ! number of processors= # of tasks
  Integer :: MPI_COMM_POOL, MPI_COMM_ROW  
  ! comunicator handle intra-pool
  !      "        "    inter-pool
  !
  ! general parallel information
  !
  Integer :: npool, nprocp, mypool, me, kunit  
  ! number of pools
  ! number of processors in this task pool
  ! identifier of this task pool
  ! identifier of this task within his pool
  ! granularity of k-point distribution
End Module para
!
