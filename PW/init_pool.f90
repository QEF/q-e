!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_pool()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initialize the pool :  MPI division in pools and images
  !
  USE para,      ONLY : me, mypool, npool, nprocp
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : mpime, me_image, my_image_id, nproc, &
                        nproc_image, nproc_pool, nimage, me_pool, my_pool_id, &
                        intra_image_comm, inter_image_comm, &
                        intra_pool_comm, inter_pool_comm
  USE mp_global, ONLY : mp_global_group_start
  USE parallel_include   
  !
  IMPLICIT NONE
  !
  INTEGER :: ierr = 0, numtask, taskid
  !
  !
  !
  ! ... communicators are allocated
  !
  ALLOCATE( inter_pool_comm( 0 : nimage - 1 ) )
  ALLOCATE( intra_pool_comm( 0 : nimage - 1 ) )  
  !
  ! ... and initialized
  !
  inter_pool_comm(:) = 0 
  intra_pool_comm(:) = 0
  !
#if defined (__PARA)
  !  
  ! ... set "my_image_id", "mypool" and reset "me"
  !
  ! ... my_image_id = 0 : nimage
  ! ... mypool      = 1 : npool       ==>    nimage * npool * nprocp = nproc
  ! ... me          = 1 : nprocp 
  !
  ! ... number of cpus per image
  !
  nproc_image = nproc / nimage
  !
  IF ( MOD( nproc, nimage ) /= 0 ) &
     CALL errore( 'startup', 'nproc /= nproc_image * nimage', 1 )  
  !
  my_image_id = INT( REAL( mpime ) / REAL( nproc_image ) )  
  me          = MOD( mpime, nproc_image )
  me_image    = me
  !
  CALL mp_barrier()
  !
  CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, &
                       my_image_id, mpime, intra_image_comm, ierr )
  !
  CALL errore( 'init_pool', 'intra_image_comm is wrong', ierr )
  !
  CALL mp_barrier()                       
  !     
  CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, me, mpime, inter_image_comm, ierr )  
  !
  CALL errore( 'init_pool', 'inter_image_comm is wrong', ierr )
  !
  ! ... number of cpus per pool
  !
  nproc_pool = nproc_image / npool
  !
  IF ( MOD( nproc, npool ) /= 0 ) &
     CALL errore( 'startup', 'nproc /= nproc_pool * npool', 1 )  
  !
  mypool = INT( REAL( me ) / REAL( nproc_pool ) )    
  me     = MOD( me, nproc_pool )
  !
  ! ... This is added for compatibility with PVM notations
  ! ... parent process (source) will have me=1 - child process me=2,...,NPROC
  !
  me     = me + 1
  mypool = mypool + 1
  !
  CALL mp_barrier( intra_image_comm )
  !
  CALL MPI_COMM_SPLIT( intra_image_comm, mypool, &
                       me_image, intra_pool_comm(my_image_id), ierr )
  !
  CALL errore( 'init_pool', 'intra_pool_comm is wrong', ierr )
  !
  CALL mp_barrier( intra_image_comm )
  !
  CALL MPI_COMM_SPLIT( intra_image_comm, me, &
                       me_image, inter_pool_comm(my_image_id), ierr )
  !
  call errore( 'init_pool', 'inter_pool_comm is wrong', ierr )
  !
  ! ... compatibility with old PWscf routines
  !
  nprocp = nproc_pool
  !  
  ! ... Initialize globally accessible pool variables
  !
  ! ... me      =>  me_pool + 1
  ! ... mypool  =>  my_pool_id + 1
  ! ... nprocp  =>  nproc_pool
  ! 
  CALL mp_global_group_start( ( me - 1 ), ( mypool - 1 ), nprocp, npool )
  !                            
#endif
  !
#if defined (__NEW_PARALLEL_DEBUG)  
  PRINT *, ""  
  PRINT *, "MPIME       = ", MPIME
  PRINT *, "ME          = ", ME
  PRINT *, "MY_IMAGE_ID = ", MY_IMAGE_ID
  PRINT *, "MYPOOL      = ", MYPOOL 
  PRINT *, ""
  CALL MPI_COMM_RANK( intra_image_comm, taskid, ierr )
  CALL MPI_COMM_SIZE( intra_image_comm, numtask, ierr )
  PRINT *, "intra_image_comm : ", taskid, numtask
  CALL MPI_COMM_RANK( inter_image_comm, taskid, ierr )
  CALL MPI_COMM_SIZE( inter_image_comm, numtask, ierr )
  PRINT *, "inter_image_comm : ", taskid, numtask  
  CALL MPI_COMM_RANK( intra_pool_comm(my_image_id), taskid, ierr )
  CALL MPI_COMM_SIZE( intra_pool_comm(my_image_id), numtask, ierr )
  PRINT *, "intra_pool_comm : ", taskid, numtask
  CALL MPI_COMM_RANK( inter_pool_comm(my_image_id), taskid, ierr )
  CALL MPI_COMM_SIZE( inter_pool_comm(my_image_id), numtask, ierr )
  PRINT *, "inter_pool_comm : ", taskid, numtask  
  PRINT *, ""  
  ierr = 666
  IF ( me_image == 0 ) ierr = my_image_id
  PRINT *, "1:MPIME = ", MPIME, " IERR = ", IERR
  CALL mp_bcast( ierr, 0, inter_pool_comm(my_image_id) )
  IF ( ierr /= 666 ) ierr = ierr + my_pool_id
  PRINT *, "2:MPIME = ", MPIME, " IERR = ", IERR
  CALL mp_bcast( ierr, 0, intra_pool_comm(my_image_id) )
  IF ( ierr /= 666 ) ierr = ierr + me_pool
  PRINT *, "3:MPIME = ", MPIME, " IERR = ", IERR
  ierr = 999
  IF ( me_image == 0 ) ierr = 100 + my_image_id
  PRINT *, "4:MPIME = ", MPIME, " IERR = ", IERR
  CALL mp_bcast( ierr, 0, intra_pool_comm(my_image_id) )
  IF ( ierr /= 999 ) ierr = ierr + me_pool
  PRINT *, "5:MPIME = ", MPIME, " IERR = ", IERR
  CALL mp_bcast( ierr, 0, inter_pool_comm(my_image_id) )
  IF ( ierr /= 999 ) ierr = ierr + my_pool_id
  PRINT *, "6:MPIME = ", MPIME, " IERR = ", IERR
  PRINT *, ""    
  !
  CALL stop_pw( .FALSE. )
  !
#endif    
  !
  RETURN
  !
END SUBROUTINE init_pool
