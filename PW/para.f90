!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE para_const
  !
  SAVE
  !
  INTEGER, PARAMETER :: &
      maxproc  = 128     !  maximum number of processors
END MODULE para_const
!
!
MODULE pfft
  !
  ! ... parallel fft information for the dense grid
  !
  USE para_const
  !
  SAVE
  !
  INTEGER :: &
      npp(maxproc),     &!  number of plane per processor
      ncp(maxproc),     &!  number of (density) columns per proc
      ncp0(maxproc),    &!  starting column for each processor
      ncplane,          &!  number of columns in a plane
      nct,              &!  total number of non-zero columns
      nxx                !  local fft data size
  ! 
  ! ... Warning these arrays moved to dfftp and dffts
  !
  ! INTEGER, POINTER :: &
  !   ipc(:),           &!  index saying which proc owns columns in a plane
  !   icpl(:)            !  index relating columns and pos. in the plan
  !
END MODULE pfft
!
!
MODULE pffts
  !
  ! parallel fft information for the smooth grid
  !
  USE para_const
  !
  SAVE
  !
  INTEGER :: &
       nkcp(maxproc)     !  number of (wfs) columns per processor
  INTEGER :: &
       npps(maxproc),   &!  number of plane per processor
       ncps(maxproc),   &!  number of (density) columns per proc
       ncp0s(maxproc),  &!  starting column for each processor
       ncplanes,        &!  number of columns in a plane
       ncts,            &!  total number of non-zero columns
       nxxs              !  local fft data size
  ! 
  ! ... Warning these arrays moved to dfftp and dffts
  !       
  ! INTEGER , POINTER :: &
  !     ipcs(:),        &!  index saying which proc owns columns in a plane
  !     icpls(:)         !  index relating columns and positions in the plan   
  !
END MODULE pffts
!
!
MODULE para
  !
  USE pfft
  USE pffts
  USE mp_global, ONLY :  nproc
  !
  SAVE
  !
  ! ... number of processors= # of tasks
  !
  INTEGER :: &
      MPI_COMM_POOL,    &!  comunicator handle intra-pool
      MPI_COMM_ROW       !       "        "    inter-pool
  !
  ! ... general parallel information
  !
  INTEGER :: &
      npool,            &!  number of pools
      nprocp,           &!  number of processors in this task pool
      mypool,           &!  identifier of this task pool
      me,               &!  identifier of this task within his pool
      kunit              !  granularity of k-point distribution
  !
END MODULE para
