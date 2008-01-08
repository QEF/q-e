!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine data_structure( lgamma )
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
  USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                         ngm, ngm_l, ngm_g, gcutm, ecutwfc
  USE gsmooth,    ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs, &
                         ngms, ngms_l, ngms_g, gcutms
  USE mp,         ONLY : mp_sum, mp_max
  USE mp_global,  ONLY : intra_pool_comm, nproc_pool, me_pool, my_image_id, &
                         nogrp, nproc, inter_pool_comm
  USE stick_base
  USE fft_scalar, ONLY : good_fft_dimension
  USE fft_types,  ONLY : fft_dlay_allocate, fft_dlay_set, fft_dlay_scalar
  !
  USE control_flags, ONLY : use_task_groups
  USE task_groups,   ONLY : task_groups_init
  !
  implicit none
  logical, intent(in) :: lgamma
  integer :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !

  real(DP) :: amod
  ! modulus of G vectors

  integer, allocatable :: stw(:,:)
  ! sticks maps

  integer :: ub(3), lb(3)  
  ! upper and lower bounds for maps

  real(DP) :: gkcut
  ! cut-off for the wavefunctions

  integer  :: ncplane, nxx
  integer  :: ncplanes, nxxs

#ifdef __PARA
  integer, allocatable :: st(:,:), sts(:,:) 
  ! sticks maps

  integer, allocatable :: ngc (:), ngcs (:), ngkc (:)
  integer  ::  ncp (nproc), nct, nkcp (nproc), ncts, ncps(nproc)
  integer  ::  ngp (nproc), ngps(nproc), ngkp (nproc), ncp_(nproc),&
       i, j, jj, idum

  !      nxx                !  local fft data dim
  !      ncplane,          &!  number of columns in a plane
  !      nct,              &!  total number of non-zero columns
  !      ncp(nproc),     &!  number of (density) columns per proc


  logical :: tk = .TRUE.   
  ! map type: true for full space sticks map, false for half space sticks map
  integer, allocatable :: in1(:), in2(:), idx(:)
  ! sticks coordinates

  !
  !  Subroutine body
  !

  tk = .NOT. lgamma

  !
  ! set the values of fft arrays
  !

  nrx1  = good_fft_dimension (nr1)
  nrx1s = good_fft_dimension (nr1s)
  !
  ! nrx2 is there just for compatibility
  !
  nrx2  = nr2
  nrx2s = nr2s
  nrx3  = good_fft_dimension (nr3)
  nrx3s = good_fft_dimension (nr3s)
  !
  !     compute number of columns per plane for each processor
  !
  ncplane  = nrx1 * nrx2
  ncplanes = nrx1s * nrx2s
  !
  !
  ! check the number of plane per process
  !
  if ( nr3 < nproc_pool ) &
    call infomsg ('data_structure', 'some processors have no planes ')

  if ( nr3s < nproc_pool ) &
    call infomsg ('data_structure', 'some processors have no smooth planes ')

  !
  ! compute gkcut calling an internal procedure
  !

  call calculate_gkcut()  

#ifdef DEBUG
  WRITE( stdout, '(5x,"ecutrho & ecutwfc",2f12.2)') tpiba2 * gcutm, &
       tpiba2 * gkcut
#endif
  !
  !
  !     Now compute for each point of the big plane how many column have
  !     non zero vectors on the smooth and thick mesh
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1
  !
  ub =  (/  n1,  n2,  n3 /)
  lb =  (/ -n1, -n2, -n3 /)
!
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

  nct  = COUNT( st  > 0 )
  ncts = COUNT( sts > 0 )

  if ( nct > ncplane )    &
     &    call errore('data_structure','too many sticks',1)

  if ( ncts > ncplanes )  &
     &    call errore('data_structure','too many sticks',2)

  if ( nct  == 0 ) &
     &    call errore('data_structure','number of sticks 0', 1)

  if ( ncts == 0 ) &
     &    call errore('data_structure','number smooth sticks 0', 1)

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
    END IF
  END IF

  CALL fft_dlay_allocate( dfftp, nproc_pool, nrx1,  nrx2  )
  CALL fft_dlay_allocate( dffts, nproc_pool, nrx1s, nrx2s )

  !  here set the fft data layout structures for dense and smooth mesh,
  !  according to stick distribution

  CALL fft_dlay_set( dfftp, &
       tk, nct, nr1, nr2, nr3, nrx1, nrx2, nrx3, (me_pool+1), &
       nproc_pool, nogrp, ub, lb, idx, in1(:), in2(:), ncp, nkcp, ngp, ngkp, st, stw)
  CALL fft_dlay_set( dffts, &
       tk, ncts, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, (me_pool+1), &
       nproc_pool, nogrp, ub, lb, idx, in1(:), in2(:), ncps, nkcp, ngps, ngkp, sts, stw)

  !  if tk = .FALSE. only half reciprocal space is considered, then we
  !  need to correct the number of sticks

  IF( .NOT. tk ) THEN
    nct  = nct*2  - 1
    ncts = ncts*2 - 1
  END IF

  !
  ! set the number of plane per process
  !

  ! npp ( 1 : nproc_pool ) = dfftp%npp ( 1 : nproc_pool )
  ! npps( 1 : nproc_pool ) = dffts%npp ( 1 : nproc_pool )

  if ( dfftp%nnp /= ncplane )    &
     &    call errore('data_structure','inconsistent plane dimension on dense grid', ABS(dfftp%nnp-ncplane) )

  if ( dffts%nnp /= ncplanes )    &
     &    call errore('data_structure','inconsistent plane dimension on smooth grid', ABS(dffts%nnp-ncplanes) )

  WRITE( stdout, '(/5x,"Planes per process (thick) : nr3 =", &
       &        i3," npp = ",i3," ncplane =",i5)') nr3, dfftp%npp (me_pool + 1) , ncplane

  if ( nr3s /= nr3 ) WRITE( stdout, '(5x,"Planes per process (smooth): nr3s=",&
       &i3," npps= ",i3," ncplanes=",i5)') nr3s, dffts%npp (me_pool + 1) , ncplanes

  WRITE( stdout,*)
  WRITE( stdout,'(5X,                                                     &
    & "Proc/  planes cols     G    planes cols    G      columns  G"/5X,  &
    & "Pool       (dense grid)       (smooth grid)      (wavefct grid)")')
  do i=1,nproc_pool
    WRITE( stdout,'(5x,i3,2x,2(i5,i7,i9),i7,i9)') i, dfftp%npp(i), ncp(i), ngp(i), &
      &        dffts%npp(i), ncps(i), ngps(i), nkcp(i), ngkp(i)
  end do
  IF ( nproc_pool > 1 ) THEN
      WRITE( stdout,'(5x,"tot",2x,2(i5,i7,i9),i7,i9)') &
      SUM(dfftp%npp(1:nproc_pool)), SUM(ncp(1:nproc_pool)), SUM(ngp(1:nproc_pool)), &
      SUM(dffts%npp(1:nproc_pool)),SUM(ncps(1:nproc_pool)),SUM(ngps(1:nproc_pool)),&
      SUM(nkcp(1:nproc_pool)), SUM(ngkp(1:nproc_pool))
  END IF
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

  !
  ! nxx is just a copy 
  !
  nxx   = nrxx
  nxxs  = nrxxs


#else

  nrx1 = good_fft_dimension (nr1)
  nrx1s = good_fft_dimension (nr1s)
  !
  !     nrx2 and nrx3 are there just for compatibility
  !
  nrx2 = nr2
  nrx3 = nr3

  nrxx = nrx1 * nrx2 * nrx3
  nrx2s = nr2s
  nrx3s = nr3s
  nrxxs = nrx1s * nrx2s * nrx3s

  ! nxx is just a copy 
  !
  nxx   = nrxx
  nxxs  = nrxxs

  CALL fft_dlay_allocate( dfftp, nproc_pool, MAX(nrx1, nrx3),  nrx2  )
  CALL fft_dlay_allocate( dffts, nproc_pool, MAX(nrx1s, nrx3s), nrx2s )

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

  do i1 = - n1, n1
     !
     ! Gamma-only: exclude space with x<0
     !
     if (lgamma .and. i1 < 0) go to 10 
     !
     do i2 = - n2, n2
        !
        ! Gamma-only: exclude plane with x=0, y<0
        !
        if(lgamma .and. i1 == 0.and. i2 < 0) go to 20
        !
        do i3 = - n3, n3
           !
           ! Gamma-only: exclude line with x=0, y=0, z<0
           !
           if(lgamma .and. i1 == 0 .and. i2 == 0 .and. i3 < 0) go to 30
           !
           amod = (i1 * bg (1, 1) + i2 * bg (1, 2) + i3 * bg (1, 3) ) **2 + &
                  (i1 * bg (2, 1) + i2 * bg (2, 2) + i3 * bg (2, 3) ) **2 + &
                  (i1 * bg (3, 1) + i2 * bg (3, 2) + i3 * bg (3, 3) ) **2
           if (amod <= gcutm)  ngm  = ngm  + 1
           if (amod <= gcutms) ngms = ngms + 1
           if (amod <= gkcut ) then
              stw( i2, i3 ) = 1
              if (lgamma) stw( -i2, -i3 ) = 1
           end if
30         continue
        enddo
20      continue
     enddo
10   continue
  enddo

  call fft_dlay_scalar( dfftp, ub, lb, nr1, nr2, nr3, nrx1, nrx2, nrx3, stw )
  call fft_dlay_scalar( dffts, ub, lb, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, stw )

  deallocate( stw )

#endif

  IF( nxx /= dfftp%nnr ) &
     CALL errore( ' data_structure ', ' inconsistent value for nxx ', ABS( nxx - dfftp%nnr ) )
  IF( nxxs /= dffts%nnr ) &
     CALL errore( ' data_structure ', ' inconsistent value for nxxs ', ABS( nxxs - dffts%nnr ) )
  !
  !     compute the global number of g, i.e. the sum over all processors
  !     within a pool
  !
  ngm_l  = ngm
  ngms_l = ngms
  ngm_g  = ngm
  ngms_g = ngms
  call mp_sum( ngm_g , intra_pool_comm )
  call mp_sum( ngms_g, intra_pool_comm )

  IF( use_task_groups ) THEN
    !
    !  Initialize task groups.
    !  Note that this call modify dffts adding task group data.
    !
    CALL task_groups_init( dffts )
    !
  END IF


  return

contains

  subroutine calculate_gkcut()
    !
    integer :: kpoint
    !
    if (nks == 0) then
       !
       ! if k-points are automatically generated (which happens later)
       ! use max(bg)/2 as an estimate of the largest k-point
       !
       gkcut = sqrt (ecutwfc) / tpiba + 0.5d0 * max ( &
          sqrt (bg (1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2), &
          sqrt (bg (1, 2) **2 + bg (2, 2) **2 + bg (3, 2) **2), &
          sqrt (bg (1, 3) **2 + bg (2, 3) **2 + bg (3, 3) **2) )
    else
       gkcut = 0.0d0
       do kpoint = 1, nks
          gkcut = max (gkcut, sqrt (ecutwfc) / tpiba + sqrt ( &
             xk (1, kpoint) **2 + xk (2, kpoint) **2 + xk (3, kpoint) **2) )
       enddo
    endif
    gkcut = gkcut**2 
    !
    ! ... find maximum value among all the processors
    !
    call mp_max (gkcut, inter_pool_comm )
    !
    return
  end subroutine calculate_gkcut


end subroutine data_structure

