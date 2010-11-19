!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE data_structure( lgamma )
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays
  ! (both the smooth and the hard mesh)
  ! In the parallel case, it distributes columns to processes, too
  !
  USE io_global,  ONLY : stdout
  USE fft_base,   ONLY : dfftp, dffts
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : bg, tpiba, tpiba2
  USE klist,      ONLY : xk, nks
  USE gvect,      ONLY : ngm, ngm_g, gcutm, ecutwfc
  USE grid_dimensions, ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE smooth_grid_dimensions, &
                  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE gsmooth,    ONLY : ngms, ngms_g, gcutms
  USE mp,         ONLY : mp_sum, mp_max
  USE mp_global,  ONLY : intra_pool_comm, nproc_pool, me_pool, my_image_id, &
                         nogrp, nproc, inter_pool_comm,  use_task_groups
  USE stick_base
  USE fft_scalar, ONLY : good_fft_dimension
  USE fft_types,  ONLY : fft_dlay_allocate, fft_dlay_set, fft_dlay_scalar
  !
  USE task_groups,   ONLY : task_groups_init
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in) :: lgamma
  INTEGER :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !

  REAL(DP) :: amod
  ! modulus of G vectors

  INTEGER, ALLOCATABLE :: stw(:,:)
  ! sticks maps

  INTEGER :: ub(3), lb(3)
  ! upper and lower bounds for maps

  REAL(DP) :: gkcut
  ! cut-off for the wavefunctions

  INTEGER  :: ncplane, nxx
  INTEGER  :: ncplanes, nxxs

#ifdef __PARA
  INTEGER, ALLOCATABLE :: st(:,:), sts(:,:)
  ! sticks maps

  INTEGER, ALLOCATABLE :: ngc (:), ngcs (:), ngkc (:)
  INTEGER  ::  ncp (nproc), nct, nkcp (nproc), ncts, ncps(nproc)
  INTEGER  ::  ngp (nproc), ngps(nproc), ngkp (nproc), ncp_(nproc),&
       i, j, jj, idum

  !      nxx              !  local fft data dim
  !      ncplane,        &!  number of columns in a plane
  !      nct,            &!  total number of non-zero columns
  !      ncp(nproc),     &!  number of (density) columns per proc


  LOGICAL :: tk = .true.
  ! map type: true for full space sticks map, false for half space sticks map
  INTEGER, ALLOCATABLE :: in1(:), in2(:), idx(:)
  ! sticks coordinates

  !
  !  Subroutine body
  !

  tk = .not. lgamma

  !
  ! set the values of fft arrays
  !
  nr1x  = good_fft_dimension (nr1)
  nr1sx = good_fft_dimension (nr1s)
  nr2x  = nr2          ! nr2x is there just for compatibility
  nr2sx = nr2s
  nr3x  = good_fft_dimension (nr3)
  nr3sx = good_fft_dimension (nr3s)

  ! compute number of points per plane
  ncplane  = nr1x * nr2x
  ncplanes = nr1sx * nr2sx

  !
  ! check the number of plane per process
  !
  IF ( nr3 < nproc_pool ) &
    CALL infomsg ('data_structure', 'some processors have no planes ')

  IF ( nr3s < nproc_pool ) &
    CALL infomsg ('data_structure', 'some processors have no smooth planes ')

  !
  ! compute gkcut calling an internal procedure
  !
  CALL calculate_gkcut()

#ifdef DEBUG
  WRITE( stdout, '(5x,"ecutrho & ecutwfc",2f12.2)') tpiba2 * gcutm, &
       tpiba2 * gkcut
#endif

  !
  !     Now compute for each point of the big plane how many column have
  !     non zero vectors on the smooth and thick mesh
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1

  ub =  (/  n1,  n2,  n3 /)
  lb =  (/ -n1, -n2, -n3 /)

  ALLOCATE( stw ( lb(1) : ub(1), lb(2) : ub(2) ) )
  ALLOCATE( st  ( lb(1) : ub(1), lb(2) : ub(2) ) )
  ALLOCATE( sts ( lb(1) : ub(1), lb(2) : ub(2) ) )

!
! ...     Fill in the stick maps, for given g-space base (b1,b2,b3)
! ...     and cut-offs
! ...     The value of the element (i,j) of the map ( st ) is equal to the
! ...     number of G-vector belonging to the (i,j) stick.
!
  CALL sticks_maps( tk, ub, lb, bg(:,1), bg(:,2), bg(:,3), gcutm, gkcut, gcutms, st, stw, sts )

  nct  = count( st  > 0 )
  ncts = count( sts > 0 )

  IF ( nct > ncplane )    &
     &    CALL errore('data_structure','too many sticks',1)

  IF ( ncts > ncplanes )  &
     &    CALL errore('data_structure','too many sticks',2)

  IF ( nct  == 0 ) &
     &    CALL errore('data_structure','number of sticks 0', 1)

  IF ( ncts == 0 ) &
     &    CALL errore('data_structure','number smooth sticks 0', 1)

  !
  !   local pointers deallocated at the end
  !
  ALLOCATE( in1( nct ), in2( nct ) )
  ALLOCATE( ngc( nct ), ngcs( nct ), ngkc( nct ) )
  ALLOCATE( idx( nct ) )

!
! ...     initialize the sticks indexes array ist
! ...     nct counts columns containing G-vectors for the dense grid
! ...     ncts counts columns contaning G-vectors for the smooth grid
!

  CALL sticks_countg( tk, ub, lb, st, stw, sts, in1, in2, ngc, ngkc, ngcs )

  CALL sticks_sort( ngc, ngkc, ngcs, nct, idx )

  CALL sticks_dist( tk, ub, lb, idx, in1, in2, ngc, ngkc, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngkp, ngps, st, stw, sts )

  CALL sticks_pairup( tk, ub, lb, idx, in1, in2, ngc, ngkc, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngkp, ngps, st, stw, sts )

  !  set the total number of G vectors

  IF( tk ) THEN
    ngm  = ngp ( me_pool + 1 )
    ngms = ngps( me_pool + 1 )
  ELSE
    IF( st( 0, 0 ) == ( me_pool + 1 ) ) THEN
      ngm  = ngp ( me_pool + 1 ) / 2 + 1
      ngms = ngps( me_pool + 1 ) / 2 + 1
    ELSE
      ngm  = ngp ( me_pool + 1 ) / 2
      ngms = ngps( me_pool + 1 ) / 2
    ENDIF
  ENDIF

  CALL fft_dlay_allocate( dfftp, nproc_pool, nr1x,  nr2x  )
  CALL fft_dlay_allocate( dffts, nproc_pool, nr1sx, nr2sx )

  !  here set the fft data layout structures for dense and smooth mesh,
  !  according to stick distribution

  CALL fft_dlay_set( dfftp, &
       tk, nct, nr1, nr2, nr3, nr1x, nr2x, nr3x, (me_pool+1), &
       nproc_pool, nogrp, ub, lb, idx, in1(:), in2(:), ncp, nkcp, ngp, ngkp, st, stw)
  CALL fft_dlay_set( dffts, &
       tk, ncts, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, (me_pool+1), &
       nproc_pool, nogrp, ub, lb, idx, in1(:), in2(:), ncps, nkcp, ngps, ngkp, sts, stw)

  !  if tk = .FALSE. only half reciprocal space is considered, then we
  !  need to correct the number of sticks

  IF( .not. tk ) THEN
    nct  = nct*2  - 1
    ncts = ncts*2 - 1
  ENDIF

  !
  ! set the number of plane per process
  !

  ! npp ( 1 : nproc_pool ) = dfftp%npp ( 1 : nproc_pool )
  ! npps( 1 : nproc_pool ) = dffts%npp ( 1 : nproc_pool )

  IF ( dfftp%nnp /= ncplane )    &
     &    CALL errore('data_structure','inconsistent plane dimension on dense grid', abs(dfftp%nnp-ncplane) )

  IF ( dffts%nnp /= ncplanes )    &
     &    CALL errore('data_structure','inconsistent plane dimension on smooth grid', abs(dffts%nnp-ncplanes) )

  WRITE( stdout, '(/5x,"Planes per process (thick) : nr3 =", &
       &        i4," npp = ",i4," ncplane =",i6)') nr3, dfftp%npp (me_pool + 1) , ncplane

  IF ( nr3s /= nr3 ) WRITE( stdout, '(5x,"Planes per process (smooth): nr3s=",&
       &i4," npps= ",i4," ncplanes=",i6)') nr3s, dffts%npp (me_pool + 1) , ncplanes

  WRITE( stdout,*)
  WRITE( stdout,'(5X,                                                     &
    & "Proc/  planes cols     G    planes cols    G      columns  G"/5X,  &
    & "Pool       (dense grid)       (smooth grid)      (wavefct grid)")')
  DO i=1,nproc_pool
    WRITE( stdout,'(5x,i4,1x,2(i5,i7,i9),i7,i9)') i, dfftp%npp(i), ncp(i), ngp(i), &
      &        dffts%npp(i), ncps(i), ngps(i), nkcp(i), ngkp(i)
  ENDDO
  IF ( nproc_pool > 1 ) THEN
      WRITE( stdout,'(5x,"tot",2x,2(i5,i7,i9),i7,i9)') &
      sum(dfftp%npp(1:nproc_pool)), sum(ncp(1:nproc_pool)), sum(ngp(1:nproc_pool)), &
      sum(dffts%npp(1:nproc_pool)),sum(ncps(1:nproc_pool)),sum(ngps(1:nproc_pool)),&
      sum(nkcp(1:nproc_pool)), sum(ngkp(1:nproc_pool))
  ENDIF
  WRITE( stdout,*)

  DEALLOCATE( stw, st, sts, in1, in2, idx, ngc, ngcs, ngkc )

  !
  !   ncp0 = starting column for each processor
  !

  ! ncp0( 1:nproc_pool )  = dfftp%iss( 1:nproc_pool )
  ! ncp0s( 1:nproc_pool ) = dffts%iss( 1:nproc_pool )

  !
  !  array ipc and ipcl ( ipc contain the number of the
  !                       column for that processor or zero if the
  !                       column do not belong to the processor,
  !                       ipcl contains the point in the plane for
  !                       each column)
  !
  !  ipc ( 1:ncplane )    = >  dfftp%isind( 1:ncplane )
  !  icpl( 1:nct )        = >  dfftp%ismap( 1:nct )

  !  ipcs ( 1:ncplanes )  = >  dffts%isind( 1:ncplanes )
  !  icpls( 1:ncts )      = >  dffts%ismap( 1:ncts )

  nrxx  = dfftp%nnr
  nrxxs = dffts%nnr

  nrxx = max( nrxx, nrxxs )
  !
  ! nxx is just a copy
  !
  nxx   = nrxx
  nxxs  = nrxxs


#else

  nr1x = good_fft_dimension (nr1)
  nr1sx = good_fft_dimension (nr1s)
  !
  !     nr2x and nr3x are there just for compatibility
  !
  nr2x = nr2
  nr3x = nr3

  nrxx = nr1x * nr2x * nr3x
  nr2sx = nr2s
  nr3sx = nr3s
  nrxxs = nr1sx * nr2sx * nr3sx

  ! nxx is just a copy
  !
  nxx   = nrxx
  nxxs  = nrxxs

  CALL fft_dlay_allocate( dfftp, nproc_pool, max(nr1x, nr3x),  nr2x  )
  CALL fft_dlay_allocate( dffts, nproc_pool, max(nr1sx, nr3sx), nr2sx )

  CALL calculate_gkcut()

  !
  !     compute the number of g necessary to the calculation
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1

  ngm = 0
  ngms = 0

  ub =  (/  n1,  n2,  n3 /)
  lb =  (/ -n1, -n2, -n3 /)
!
  ALLOCATE( stw ( lb(2):ub(2), lb(3):ub(3) ) )
  stw = 0

  DO i1 = - n1, n1
     !
     ! Gamma-only: exclude space with x<0
     !
     IF (lgamma .and. i1 < 0) GOTO 10
     !
     DO i2 = - n2, n2
        !
        ! Gamma-only: exclude plane with x=0, y<0
        !
        IF(lgamma .and. i1 == 0.and. i2 < 0) GOTO 20
        !
        DO i3 = - n3, n3
           !
           ! Gamma-only: exclude line with x=0, y=0, z<0
           !
           IF(lgamma .and. i1 == 0 .and. i2 == 0 .and. i3 < 0) GOTO 30
           !
           amod = (i1 * bg (1, 1) + i2 * bg (1, 2) + i3 * bg (1, 3) ) **2 + &
                  (i1 * bg (2, 1) + i2 * bg (2, 2) + i3 * bg (2, 3) ) **2 + &
                  (i1 * bg (3, 1) + i2 * bg (3, 2) + i3 * bg (3, 3) ) **2
           IF (amod <= gcutm)  ngm  = ngm  + 1
           IF (amod <= gcutms) ngms = ngms + 1
           IF (amod <= gkcut ) THEN
              stw( i2, i3 ) = 1
              IF (lgamma) stw( -i2, -i3 ) = 1
           ENDIF
30         CONTINUE
        ENDDO
20      CONTINUE
     ENDDO
10   CONTINUE
  ENDDO

  CALL fft_dlay_scalar( dfftp, ub, lb, nr1, nr2, nr3, nr1x, nr2x, nr3x, stw )
  CALL fft_dlay_scalar( dffts, ub, lb, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, stw )

  DEALLOCATE( stw )

#endif

  IF( nxx < dfftp%nnr ) &
     CALL errore( ' data_structure ', ' inconsistent value for nxx ', abs( nxx - dfftp%nnr ) )
  IF( nxxs /= dffts%nnr ) &
     CALL errore( ' data_structure ', ' inconsistent value for nxxs ', abs( nxxs - dffts%nnr ) )
  !
  !     compute the global number of g, i.e. the sum over all processors
  !     within a pool
  !
  ngm_g  = ngm  ; CALL mp_sum( ngm_g , intra_pool_comm )
  ngms_g = ngms ; CALL mp_sum( ngms_g, intra_pool_comm )

  IF( use_task_groups ) THEN
    !
    !  Initialize task groups.
    !  Note that this call modify dffts adding task group data.
    !
    CALL task_groups_init( dffts )
    !
  ENDIF


CONTAINS

  SUBROUTINE calculate_gkcut()

    INTEGER :: kpoint

    IF (nks == 0) THEN
       !
       ! if k-points are automatically generated (which happens later)
       ! use max(bg)/2 as an estimate of the largest k-point
       !
       gkcut = 0.5d0 * max ( &
          sqrt (sum(bg (1:3, 1)**2) ), &
          sqrt (sum(bg (1:3, 2)**2) ), &
          sqrt (sum(bg (1:3, 3)**2) ) )
    ELSE
       gkcut = 0.0d0
       DO kpoint = 1, nks
          gkcut = max (gkcut, sqrt ( sum(xk (1:3, kpoint)**2) ) )
       ENDDO
    ENDIF
    gkcut = (sqrt (ecutwfc) / tpiba + gkcut)**2
    !
    ! ... find maximum value among all the processors
    !
    CALL mp_max (gkcut, inter_pool_comm )

  END SUBROUTINE calculate_gkcut


END SUBROUTINE data_structure

