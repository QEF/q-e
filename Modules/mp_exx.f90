!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_exx
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_bcast, mp_size, mp_rank, mp_comm_split
  USE mp_bands, ONLY : nbgrp, nproc_bgrp, me_bgrp, root_bgrp, my_bgrp_id, &
       inter_bgrp_comm, intra_bgrp_comm
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Band groups (processors within a pool of bands)
  ! ... Subdivision of pool group, used for parallelization over bands
  !
  INTEGER :: negrp       = 1  ! number of band groups
  INTEGER :: nproc_egrp  = 1  ! number of processors within a band group
  INTEGER :: me_egrp     = 0  ! index of the processor within a band group
  INTEGER :: root_egrp   = 0  ! index of the root processor within a band group
  INTEGER :: my_egrp_id  = 0  ! index of my band group
  INTEGER :: inter_egrp_comm  = 0  ! inter band group communicator
  INTEGER :: intra_egrp_comm  = 0  ! intra band group communicator  
  !
  ! ... "task" groups (for band parallelization of FFT)
  !
  INTEGER :: ntask_groups = 1  ! number of proc. in an orbital "task group"
  !
  ! should the old band parallelization method be used?
  !
  LOGICAL :: use_old_exx = .FALSE.
  !
  ! variables for the pair parallelization
  !
  INTEGER :: max_pairs ! maximum pairs per band group
  INTEGER, ALLOCATABLE :: egrp_pairs(:,:,:) ! pairs for each band group
  INTEGER, ALLOCATABLE :: band_roots(:) ! root for each band
  LOGICAL, ALLOCATABLE :: contributed_bands(:,:) ! bands for which the bgroup has a pair
  INTEGER, ALLOCATABLE :: nibands(:) ! number of bands for which the bgroup has a pair
  INTEGER, ALLOCATABLE :: ibands(:,:) ! bands for which the bgroup has a pair
  INTEGER :: iexx_start = 0              ! starting band index used in bgrp parallelization
  INTEGER :: iexx_end = 0                ! ending band index used in bgrp parallelization
  INTEGER, ALLOCATABLE :: iexx_istart(:) ! starting band inded for the outer loop
  INTEGER, ALLOCATABLE :: iexx_iend(:)  ! ending band index used in the outer loop
  INTEGER, ALLOCATABLE :: all_start(:)
  INTEGER, ALLOCATABLE :: all_end(:)
  INTEGER :: max_contributors
  !
  ! flag for whether the exx part of the calculation is in progress
  !
  INTEGER :: exx_mode = 0
  !
  ! maximum number of bands for psi
  !
  INTEGER :: max_ibands
  !
  INTEGER :: jblock
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_exx( nband_, ntg_, parent_comm )
    !---------------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into nband_ pools
    ! ... Requires: nband_, read from command line
    ! ...           parent_comm, typically processors of a k-point pool
    ! ...           (intra_pool_comm)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nband_, parent_comm
    INTEGER, INTENT(IN), OPTIONAL :: ntg_
    !
    INTEGER :: parent_nproc = 1, parent_mype = 0
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... nband_ must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    negrp = nband_
    !
    IF ( negrp < 1 .OR. negrp > parent_nproc ) CALL errore( 'mp_start_bands',&
                          'invalid number of band groups, out of range', 1 )
    IF ( MOD( parent_nproc, negrp ) /= 0 ) CALL errore( 'mp_start_bands', &
        'n. of band groups  must be divisor of parent_nproc', 1 )
    ! 
    ! ... Set number of processors per band group
    !
    nproc_egrp = parent_nproc / negrp
    !
    ! ... set index of band group for this processor   ( 0 : negrp - 1 )
    !
    my_egrp_id = parent_mype / nproc_egrp
    !
    ! ... set index of processor within the image ( 0 : nproc_image - 1 )
    !
    me_egrp    = MOD( parent_mype, nproc_egrp )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the intra_egrp_comm communicator is created
    !
    CALL mp_comm_split( parent_comm, my_egrp_id, parent_mype, intra_egrp_comm )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the inter_egrp_comm communicator is created                     
    !     
    CALL mp_comm_split( parent_comm, me_egrp, parent_mype, inter_egrp_comm )  
    !
    IF ( PRESENT(ntg_) ) THEN
       ntask_groups = ntg_
    END IF
    !
    ! turn on the old band group parallelization method
    !
    IF(use_old_exx) THEN
       nbgrp = negrp
       nproc_bgrp = nproc_egrp
       me_bgrp = me_egrp
       root_bgrp = root_egrp
       my_bgrp_id = my_egrp_id
       inter_bgrp_comm = inter_egrp_comm
       intra_bgrp_comm = intra_egrp_comm
       negrp = 1
    END IF
#endif
    RETURN
    !
  END SUBROUTINE mp_start_exx
  !
  SUBROUTINE init_index_over_band(comm,nbnd,m)
    !
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: comm, nbnd
    INTEGER, INTENT(IN) :: m
    
    INTEGER :: npe, myrank, rest, k, rest_i, k_i
    INTEGER :: i, j, ipair, iegrp, root
    INTEGER :: ibnd, npairs, ncontributing
    INTEGER :: n_underloaded ! number of band groups that are under max load
    INTEGER :: pair_bands(nbnd,nbnd)

    jblock = 7
    
    max_ibands = CEILING(float(nbnd)/float(negrp))+2
    IF (ALLOCATED(all_start)) THEN
       DEALLOCATE( all_start, all_end )
       DEALLOCATE( iexx_istart, iexx_iend )
    END IF
    ALLOCATE( all_start(negrp) )
    ALLOCATE( all_end(negrp) )
    ALLOCATE( iexx_istart(negrp) )
    ALLOCATE( iexx_iend(negrp) )
    
    myrank = mp_rank(comm)
    npe = mp_size(comm)

    rest = mod(nbnd, npe)
    k = int(nbnd/npe)
    
    IF ( k >= 1 ) THEN
       IF ( rest > myrank ) THEN
          iexx_start = (myrank)*k + (myrank+1)
          iexx_end = (myrank+1)*k + (myrank+1)
       ELSE
          iexx_start = (myrank)*k + rest + 1
          iexx_end = (myrank+1)*k + rest
       END IF
    ELSE
       IF(my_egrp_id+1.le.m) THEN
          iexx_start = my_egrp_id+1
          iexx_end = my_egrp_id+1
       ELSE
          iexx_start = 0
          iexx_end = 0
       END IF
    END IF

    !determine iexx_start and iexx_end for all of the other bands
    all_start = 0
    all_end = 0
    DO i=1, negrp
       IF ( k >= 1 ) THEN
          IF ( rest > i-1 ) THEN
             all_start(i) = (i-1)*k + i
             all_end(i) = i*k + i
          ELSE
             all_start(i) = (i-1)*k + rest + 1
             all_end(i) = i*k + rest
          END IF
       ELSE
          IF(i.le.m) THEN
             all_start(i) = i
             all_end(i) = i
          ELSE
             all_start(i) = 0
             all_end(i) = 0
          END IF
       END IF
    END DO

    !determine the first and last indices for the outer loop
    rest_i = mod(m, npe)
    k_i = int(m/npe)
    DO i=1, negrp
       IF ( k_i >= 1 ) THEN
          IF ( rest_i > i-1 ) THEN
             iexx_istart(i) = (i-1)*k_i + i
             iexx_iend(i) = i*k_i + i
          ELSE
             iexx_istart(i) = (i-1)*k_i + rest_i + 1
             iexx_iend(i) = i*k_i + rest_i
          END IF
       ELSE
          IF(i.le.m) THEN
             iexx_istart(i) = i
             iexx_iend(i) = i
          ELSE
             iexx_istart(i) = 0
             iexx_iend(i) = 0
          END IF
       END IF
    END DO
    max_pairs = CEILING(REAL(nbnd*m)/REAL(negrp))
    n_underloaded = MODULO(max_pairs*negrp-nbnd*m,negrp)
    !
    ! allocate arrays
    !
    IF (allocated(egrp_pairs)) THEN
       DEALLOCATE(egrp_pairs)
       DEALLOCATE(band_roots)
       DEALLOCATE(contributed_bands)
       DEALLOCATE(nibands)
       DEALLOCATE(ibands)
    END IF
    !
    IF (.not.allocated(egrp_pairs)) THEN
       ALLOCATE(egrp_pairs(2,max_pairs,negrp))
       ALLOCATE(band_roots(m))
       ALLOCATE(contributed_bands(nbnd,negrp))
       ALLOCATE(nibands(negrp))
       ALLOCATE(ibands(nbnd,negrp))
    END IF
    !
    ! assign the pairs for each band group
    !
    pair_bands = 0
    egrp_pairs = 0
    j = 1
    DO iegrp=1, negrp
       npairs = max_pairs
       IF (iegrp.le.n_underloaded) npairs = npairs - 1
       DO ipair=1, npairs
          !
          ! get the first value of i for which the (i,j) pair has not been 
          ! assigned yet
          !
          i = 1
          DO WHILE (pair_bands(i,j).gt.0)
             i = i + 1
             IF(i.gt.m) exit
          END DO
          IF (i.le.m) THEN
             pair_bands(i,j) = iegrp
          END IF
          egrp_pairs(1,ipair,iegrp) = i
          egrp_pairs(2,ipair,iegrp) = j
          
          j = j + 1
          IF (j.gt.nbnd) j = 1
       END DO
       
    END DO
    !
    ! determine the bands for which this band group will calculate a pair
    !
    contributed_bands = .FALSE.
    DO iegrp=1, negrp
       npairs = max_pairs
       IF (iegrp.le.n_underloaded) npairs = npairs - 1
       DO ipair=1, npairs
          contributed_bands(egrp_pairs(1,ipair,iegrp),iegrp) = .TRUE.
       END DO
    END DO
    nibands = 0
    ibands = 0
    DO iegrp=1, negrp
       DO i=1, nbnd
          IF (contributed_bands(i,iegrp)) THEN
             nibands(iegrp) = nibands(iegrp) + 1
             ibands(nibands(iegrp),iegrp) = i
          END IF
       END DO
    END DO
    !
    ! determine the maximum number of contributing egrps for any band
    !
    max_contributors = 0
    DO i=1, nbnd
       !
       ! determine the number of sending egrps
       !
       ncontributing = 0
       DO iegrp=1, negrp
          IF(contributed_bands(i,iegrp)) THEN
             ncontributing = ncontributing + 1
          END IF
       END DO
       IF(ncontributing.gt.max_contributors) THEN
          max_contributors = ncontributing
       END IF
    END DO
    !
    ! create a list of the roots for each band
    !
    DO i=1, m
       DO iegrp=1, negrp
          IF(iexx_iend(iegrp).ge.i) THEN
             band_roots(i) = iegrp - 1
             exit
          END IF
       END DO       
    END DO
    !
  END SUBROUTINE init_index_over_band
  !
  SUBROUTINE set_egrp_indices(nbnd, ib_start, ib_end)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(OUT) :: ib_start, ib_end

    INTEGER :: rest, nbnd_per_bgrp

    rest = mod ( nbnd, negrp )
    nbnd_per_bgrp = int( nbnd / negrp ) 

    IF (rest > my_egrp_id) THEN 
       ib_start =  my_egrp_id    * (nbnd_per_bgrp+1) + 1
       ib_end   = (my_egrp_id+1) * (nbnd_per_bgrp+1) 
    ELSE
       ib_start =  my_egrp_id    * nbnd_per_bgrp + rest + 1
       ib_end   = (my_egrp_id+1) * nbnd_per_bgrp + rest 
    ENDIF

  END SUBROUTINE set_egrp_indices

  INTEGER FUNCTION egrp_start(nbnd)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd

    INTEGER :: rest, nbnd_per_bgrp

    rest = mod ( nbnd, negrp )
    nbnd_per_bgrp = int( nbnd / negrp ) 

    IF (rest > my_egrp_id) THEN 
       egrp_start =  my_egrp_id    * (nbnd_per_bgrp+1) + 1
    ELSE
       egrp_start =  my_egrp_id    * nbnd_per_bgrp + rest + 1
    ENDIF

  END FUNCTION egrp_start

  INTEGER FUNCTION egrp_end(nbnd)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd

    INTEGER :: rest, nbnd_per_bgrp

    rest = mod ( nbnd, negrp )
    nbnd_per_bgrp = int( nbnd / negrp ) 

    IF (rest > my_egrp_id) THEN 
       egrp_end   = (my_egrp_id+1) * (nbnd_per_bgrp+1) 
    ELSE
       egrp_end   = (my_egrp_id+1) * nbnd_per_bgrp + rest 
    ENDIF

  END FUNCTION egrp_end

END MODULE mp_exx
!
!     
