!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine data_structure_para
  !-----------------------------------------------------------------------
  !
  ! distribute columns to processes for parallel fft
  ! columns are sets of g-vectors along z: g(k) = i1*b1+i2*b2+i3*b3 ,
  ! with g^2<gcut and (i1,i2) running over the (xy) plane.
  ! Columns are "active" for a given (i1,i2) if they contain a nonzero
  ! number of wavevectors
  !

#ifdef __PARA

  use para
  use pwcom
  use mp, only: mp_sum
  use mp_global, only: intra_pool_comm, me_pool, nproc_pool
  use fft_scalar, only: good_fft_dimension
  use fft_types
  use stick_base
  !
  implicit none
  !
  integer, allocatable :: &
       ngc(:),        &! number of g-vectors per column (dense grid)
       ngcs(:),       &! number of g-vectors per column (smooth grid
       ngcw(:),       &! number of wavefct plane waves per colum
       in1(:),        &! index i for column (i1,i2)
       in2(:),        &! index j for column (i1,i2)
       index(:)        ! used to order column
  integer ic,         &! fft index for this column (dense grid)
       ics,           &! as above for the smooth grid
       icm,           &! fft index for column (-i1,-i2) (dense grid)
       icms,          &! as above for the smooth grid
       ncp_(maxproc), &! number of column per processor (work space)
       ngp(maxproc),  &! number of g-vectors per proc (dense grid)
       ngps(maxproc), &! number of g-vectors per proc (smooth grid)
       ngpw(maxproc)   ! number of wavefct plane waves per proc
  !
  integer np, nps1,           &! counters on planes
       nq, nqs,               &! counters on planes
       max1,min1,max2,min2,   &! aux. variables
       m1, m2, n1, n2, i, mc, &! generic counter
       idum, nct_,            &! check variables
       j,jj,                  &! counters on processors
       n1m1,n2m1,n3m1,        &! nr1-1 and so on
       i1, i2, i3              ! counters on G space
  logical has_gzero
  real(kind=8), allocatable :: aux(:)  ! used to order columns
  real(kind=8)  amod, gkcut        ! square modulus of G vectors

  integer, allocatable :: st(:,:), stw(:,:), sts(:,:)
  ! sticks maps, define the stick projections in the 2D xy plane
  integer :: ub(3), lb(3)
  ! upper and lower bounds for maps
  logical :: tk = .FALSE.

  !
  ! set the dimensions of fft arrays
  !
  nrx1 =good_fft_dimension( nr1 )
  nrx2 = nr2
  nrx3 =good_fft_dimension( nr3 )
  !
  nrx1s=good_fft_dimension( nr1s )
  nrx2s=nr2s
  nrx3s=good_fft_dimension( nr3s )
  !
  !     compute number of columns for each processor
  !
  ncplane = nrx1*nrx2
  ncplanes= nrx1s*nrx2s
  !
  !    local variables to be deallocated at the end
  !
  allocate  (ngc ( ncplane))    
  allocate  (ngcs( ncplane))    
  allocate  (ngcw( ncplane))    
  allocate  (in1 ( ncplane))    
  allocate  (in2 ( ncplane))    
  allocate  (index(ncplane))    
  allocate  (aux(  ncplane))    
  !
  ! set the number of plane per process
  !
  if ( nr3  < nproc_pool ) call errore('set_fft_para',                      &
       &                'some processors have no planes ',-1)
  if ( nr3s < nproc_pool ) call errore('set_fft_para',                     &
       &                'some processors have no smooth planes ',-1)
  !
  !     Now compute for each point of the big plane how many column have
  !     non zero vectors on the smooth and dense grid
  !
  ! NOTA BENE: the exact limits for a correctly sized FFT grid are:
  ! -nr/2,..,+nr/2  for nr even; -(nr-1)/2,..,+(nr-1)/2  for nr odd.
  ! If the following limits are increased, a slightly undersized fft
  ! grid, with some degree of G-vector refolding, can be used
  ! (at your own risk - a check is done in ggen).
  !
  n1m1=nr1/2+1
  n2m1=nr2/2+1
  n3m1=nr3/2+1
  !
  ub = (/  n1m1,  n2m1,  n3m1 /)
  lb = (/ -n1m1, -n2m1, -n3m1 /)
  !
  ALLOCATE( stw ( lb(1) : ub(1), lb(2) : ub(2) ) )
  ALLOCATE( st  ( lb(1) : ub(1), lb(2) : ub(2) ) )
  ALLOCATE( sts ( lb(1) : ub(1), lb(2) : ub(2) ) )
  !
  gkcut = ecutwfc / tpiba2

  CALL sticks_maps( tk, ub, lb, bg(:,1), bg(:,2), bg(:,3), &
    gcutm, gkcut, gcutms, st, stw, sts )

  nct  = COUNT( st  > 0 )
  ncts = COUNT( sts > 0 )

  if (nct > ncplane)    &
     &    call errore('data_structure','too many sticks',1)

  if (ncts > ncplanes)  &
     &    call errore('data_structure','too many sticks',2)

  if (nct == 0) &
     &    call errore('data_structure','number of sticks 0', 1)

  if (ncts == 0) &
     &    call errore('data_structure','number smooth sticks 0', 1)

  CALL sticks_countg( tk, ub, lb, st, stw, sts, in1, in2, ngc, ngcw, ngcs )

  CALL sticks_sort( ngc, ngcw, ngcs, nct, index )

  CALL sticks_dist( tk, ub, lb, index, in1, in2, ngc, ngcw, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngpw, ngps, st, stw, sts )

  CALL sticks_pairup( tk, ub, lb, index, in1, in2, ngc, ngcw, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngpw, ngps, st, stw, sts )

  !
  ! ngm contains the number of G vectors on this processor (me)
  ! remember that ngp counts all G-vectors, while we want only half
  ! same for ngms (smooth grid)
  !
  IF( st( 0, 0 ) == ( me_pool + 1 ) ) THEN
    ngm  = ngp ( me_pool + 1 ) / 2 + 1
    ngms = ngps( me_pool + 1 ) / 2 + 1
  ELSE
    ngm  = ngp ( me_pool + 1 ) / 2 
    ngms = ngps( me_pool + 1 ) / 2 
  END IF
  !

  CALL fft_dlay_allocate( dfftp, nproc_pool, nrx1,  nrx2  )
  CALL fft_dlay_allocate( dffts, nproc_pool, nrx1s, nrx2s )

  CALL fft_dlay_set( dfftp, &
       tk, nct, nr1, nr2, nr3, nrx1, nrx2, nrx3, ( me_pool + 1 ), &
       nproc_pool, ub, lb, index, in1(:), in2(:), ncp, nkcp, ngp, ngpw, st, stw)
  CALL fft_dlay_set( dffts, &
       tk, ncts, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, ( me_pool + 1 ), &
       nproc_pool, ub, lb, index, in1(:), in2(:), ncps, nkcp, ngps, ngpw, sts, stw)

  IF( .NOT. tk ) THEN
    nct  = nct*2  - 1
    ncts = ncts*2 - 1
  END IF

  npp ( 1 : nproc_pool ) = dfftp%npp ( 1 : nproc_pool )
  npps( 1 : nproc_pool ) = dffts%npp ( 1 : nproc_pool )

  write(6,*)
  write(6,'(                                                        &
    & '' Proc/  planes cols    G   planes cols    G    columns  G''/    &
    & '' Pool       (dense grid)      (smooth grid)   (wavefct grid)'')')
  do i = 1, nproc_pool
    write(6,'(i3,2x,3(i5,2i7))') i, npp(i), ncp(i), ngp(i),          &
      &        npps(i), ncps(i), ngps(i), nkcp(i), ngpw(i)
  end do
  write(6,'(i3,2x,3(i5,2i7))') 0, SUM( npp( 1:nproc_pool ) ), SUM( ncp( 1:nproc_pool ) ), &
    &   SUM( ngp( 1:nproc_pool ) ), SUM( npps( 1:nproc_pool ) ), SUM( ncps( 1:nproc_pool ) ), &
    &   SUM( ngps( 1:nproc_pool ) ), SUM( nkcp( 1:nproc_pool ) ), SUM( ngpw( 1:nproc_pool ) )
  write(6,*)

  !
  ! nxx and nxxs are copies of nrxx and nrxxs, the local fft data size,
  ! to be stored in "parallel" commons. Not a very elegant solution.
  !

  nrxx  = dfftp%nnr 
  nrxxs = dffts%nnr 
  nxx   = nrxx
  nxxs  = nrxxs

  !
  ncp0( 1:nproc_pool ) = dfftp%iss( 1:nproc_pool )
  !   ipc ( 1:ncplane )    = dfftp%isind( 1:ncplane )
  !   icpl( 1:nct )        = dfftp%ismap( 1:nct )

  ncp0s( 1:nproc_pool ) = dffts%iss( 1:nproc_pool )
  !   ipcs ( 1:ncplanes )   = dffts%isind( 1:ncplanes )
  !   icpls( 1:ncts )       = dffts%ismap( 1:ncts )

  !
  deallocate (aux)
  deallocate (index)
  deallocate (in2)
  deallocate (in1)
  deallocate (ngcw)
  deallocate (ngcs)
  deallocate (ngc )
  deallocate ( stw, st, sts )

  !
  ngm_l  = ngm
  ngms_l = ngms
  ngm_g  = ngm
  ngms_g = ngms
  call mp_sum( ngm_g , intra_pool_comm )
  call mp_sum( ngms_g, intra_pool_comm )
  !
#endif

  return
end subroutine data_structure_para
