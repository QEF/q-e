!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine data_structure( lgamma )
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays.
  ! In the parallel case distributes columns to processes, too
  ! This version computes also the smooth and hard mesh
  !
#include "machine.h"
  use pwcom
#ifdef __PARA
  use para
#endif
  use mp, only: mp_sum
  use mp_global, only: intra_pool_comm, nproc_pool, me_pool
  use stick_base
  use fft_scalar, only: good_fft_dimension
  use fft_types, only: fft_dlay_allocate, fft_dlay_set, fft_dlay_scalar
  !
  implicit none
  logical, intent(in) :: lgamma
  integer :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !

  real(kind=DP) :: amod
  ! modulus of G vectors

  integer, allocatable :: st(:,:), stw(:,:), sts(:,:) 
  ! sticks maps

  integer :: ub(3), lb(3)  
  ! upper and lower bounds for maps

  real(kind=DP) :: gkcut
  ! cut-off for the wavefunctions

  integer :: np, nps1, nq, nqs, max1, min1, max2, min2, kpoint, m1, &
       m2, i, mc, nct_, ic, ics

#ifdef __PARA
  ! counters on planes

  integer, allocatable :: ngc (:), ngcs (:), ngkc (:)
  integer  ::  ngp (maxproc), ngps(maxproc), ngkp (maxproc), ncp_(maxproc),&
       j, jj, idum
  ! counters on planes
  ! indices for including meshes
  ! counter on k points
  ! generic counters
  ! check variables
  ! number of columns per plane
  ! number of columns per plane (smooth part)
  ! number of columns per plane (hard part)
  ! from thick plane to column list
  ! from smooth plane to column list
  ! number of g per processor
  ! number of column per processor
  ! counter on processors
  ! counter on processors
  ! used for swap

  logical :: tk = .TRUE.   
  ! map type: true for full space sticks map, false for half space sticks map
  integer, allocatable :: in1(:), in2(:), index(:)
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
  !    global pointers passed with commons are allocate here
  !
  allocate (ipc( ncplane))    
  allocate (ipcs( ncplanes))    
  !
  ! set the number of plane per process
  !
  if (nr3.lt.nprocp) call errore ('data_structure', &
       'some processors have no planes ',  - 1)

  if (nr3s.lt.nprocp) call errore ('data_structure', &
       'some processors have no smooth planes ',  - 1)
  if (nprocp.eq.1) then
     npp (1) = nr3
     npps (1) = nr3s
  else
     np = nr3 / nprocp
     nq = nr3 - np * nprocp
     nps1 = nr3s / nprocp
     nqs = nr3s - nps1 * nprocp
     do i = 1, nprocp
        npp (i) = np
        if (i.le.nq) npp (i) = np + 1
        npps (i) = nps1
        if (i.le.nqs) npps (i) = nps1 + 1
     enddo
  endif
  write (6, '(/5x,"Planes per process (thick) : nr3 =", &
       &        i3," npp = ",i3," ncplane =",i5)') nr3, npp (me) , ncplane

  if (nr3s.ne.nr3) write (6, '(/5x,"Planes per process (smooth): nr3s=",&
       &i3," npps= ",i3," ncplanes=",i5)') &
       nr3s, npps (me) , ncplanes

  !if (nks.eq.0) then
  !   !
  !   ! if k-points are automatically generated (which happens later)
  !   ! use max(bg)/2 as an estimate of the largest k-point
  !   !
  !   gkcut = sqrt (ecutwfc) / tpiba + 0.5d0 * max (sqrt (bg ( &
  !        1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2), sqrt (bg (1, 2) ** &
  !        2 + bg (2, 2) **2 + bg (3, 2) **2), sqrt (bg (1, 3) **2 + bg ( &
  !        2, 3) **2 + bg (3, 3) **2) )
  !else
  !   gkcut = 0.0d0
  !   do kpoint = 1, nks
  !      gkcut = max (gkcut, sqrt (ecutwfc) / tpiba + sqrt (xk ( &
  !           1, kpoint) **2 + xk (2, kpoint) **2 + xk (3, kpoint) **2) )
  !   enddo
  !endif
!
!  gkcut = gkcut * gkcut

  call calculate_gkcut()  ! call the internal procedure

  ! find maximum among the nodes
  call poolextreme (gkcut, + 1)
  !
#ifdef DEBUG
  write (6, '(5x,"ecutrho & ecutwfc",2f12.2)') tpiba2 * gcutm, &
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
  lb(1) = -n1
  lb(2) = -n2
  lb(3) = -n3
  ub(1) =  n1
  ub(2) =  n2
  ub(3) =  n3
!
  ALLOCATE( stw ( lb(1):ub(1), lb(2):ub(2) ) )
  ALLOCATE( st  ( lb(1):ub(1), lb(2):ub(2) ) )
  ALLOCATE( sts ( lb(1):ub(1), lb(2):ub(2) ) )

  
  ipc  = 0
  ic   = 0
  ipcs = 0

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
  ALLOCATE( index( nct ) )

!
! ...     initialize the sticks indexes array ist
! ...     nct counts columns containing G-vectors for the dense grid
! ...     ncts counts columns contaning G-vectors for the smooth grid
!

  CALL sticks_countg( tk, ub, lb, st, stw, sts, in1, in2, ngc, ngkc, ngcs )

  CALL sticks_sort( ngc, ngkc, ngcs, nct, index )

  CALL sticks_dist( tk, ub, lb, index, in1, in2, ngc, ngkc, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngkp, ngps, st, stw, sts )

  CALL sticks_pairup( tk, ub, lb, index, in1, in2, ngc, ngkc, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngkp, ngps, st, stw, sts )

  CALL fft_dlay_allocate( dfftp, nproc_pool, nrx1,  nrx2  )
  CALL fft_dlay_allocate( dffts, nproc_pool, nrx1s, nrx2s )

  CALL fft_dlay_set( dfftp, &
       tk, nct, nr1, nr2, nr3, nrx1, nrx2, nrx3, (me_pool+1), &
       nproc_pool, ub, lb, index, in1(:), in2(:), ncp, nkcp, ngp, ngkp, st, stw)
  CALL fft_dlay_set( dffts, &
       tk, nct, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx2s, (me_pool+1), &
       nproc_pool, ub, lb, index, in1(:), in2(:), ncps, nkcp, ngps, ngkp, sts, stw)

  !
  !  Set the the correspondence between planes in the smooth thick mesh
  !
  if (mod (nr1s, 2) .eq.0) then
     max1 = nr1s / 2
     min1 = nr1 - nr1s / 2 + 1
  else
     max1 = (nr1s - 1) / 2 + 1
     min1 = nr1 - (nr1s - 1) / 2 + 1
  endif

  if (mod (nr2s, 2) .eq.0) then
     max2 = nr2s / 2
     min2 = nr2 - nr2s / 2 + 1
  else
     max2 = (nr2s - 1) / 2 + 1
     min2 = nr2 - (nr2s - 1) / 2 + 1
  endif


  DO mc = 1, nct

    i = index( mc )

    i1 = in1( i )
    i2 = in2( i )

    m1 = mod ( i1 , nr1 ) + 1
    if ( m1 < 1 ) m1 = m1 + nr1
    m2 = mod ( i2 , nr2 ) + 1
    if ( m2 < 1 ) m2 = m2 + nr2
    ic = m1 + (m2 - 1) * nrx1 

    IF( st( i1, i2 ) > 0 .AND. stw( i1, i2 ) > 0 ) THEN
      ipc( ic  ) = st(  i1,  i2 )
    ELSE IF( st( i1, i2 ) > 0 ) THEN
      ipc( ic  ) = -st(  i1,  i2 )
    END IF

    IF( sts( i1, i2 ) > 0 ) THEN

      if ( m1 <= max1 ) then
         n1 = m1
      elseif ( m1 >= min1 ) then
         n1 = nr1s - (nr1 - m1)
      else
         call errore ('data_structure', 'something wrong with n1', 1)
      endif

      if ( m2 <= max2 ) then
         n2 = m2
      elseif ( m2 >= min2 ) then
         n2 = nr2s - (nr2 - m2)
      else
         call errore ('data_structure', 'something wrong with n2', 1)
      endif

      !   check that the indices are within bounds
      !
      if (n1.lt.1.or.n1.gt.nr1s.or.n2.lt.1.or.n2.gt.nr2s) &
         call errore ('data_structure', 'something wrong with n1,n2', 1)

      ics = n1 + (n2 - 1) * nrx1s

      IF( stw( i1, i2 ) > 0 ) THEN
        ipcs( ics ) = sts(  i1,  i2 )
      ELSE
        ipcs( ics ) = -sts(  i1,  i2 )
      END IF

    END IF

  END DO

  write(6,*)
  write(6,'(                                                        &
    & '' Proc/  planes cols    G   planes cols    G    columns  G''/    &
    & '' Pool       (dense grid)      (smooth grid)   (wavefct grid)'')')
  do i=1,nproc_pool
    write(6,'(i3,2x,3(i5,2i7))') i, npp(i), ncp(i), ngp(i),          &
      &        npps(i), ncps(i), ngps(i), nkcp(i), ngkp(i)
  end do
  write(6,'(i3,2x,3(i5,2i7))') 0, SUM(npp(1:nproc)), SUM(ncp(1:nproc)), &
    &   SUM(ngp(1:nproc)), SUM(npps(1:nproc)), SUM(ncps(1:nproc)), &
    &   SUM(ngps(1:nproc)), SUM(nkcp(1:nproc)), SUM(ngkp(1:nproc))
  write(6,*)


  DEALLOCATE( stw, st, sts, in1, in2, index, ngc, ngcs, ngkc )

  !
  !   computing the starting column for each processor
  !   and the total number of G vectors
  !
  ngm = ngp (me)
  ngms = ngps (me)

  do i = 1, nprocp

     if (ngkp (i) .eq.0) call errore ('data_structure', &
          'some processors have no pencils, not yet implemented', 1)
     if (i.eq.1) then
        ncp0 (i) = 0
        ncp0s (i) = 0
     else
        ncp0 (i) = ncp0 (i - 1) + ncp (i - 1)
        ncp0s (i) = ncp0s (i - 1) + ncps (i - 1)
     endif

  enddo

  allocate (icpl( nct))    

  allocate (icpls( ncts))    
  do j = 1, nprocp
     ncp_ (j) = 0
  enddo
  !
  !  Now compute the array ipc and ipcl ( ipc contain the number of the
  !                        column for that processor or zero if the
  !                        column do not belong to the processor,
  !                        ipcl contains the point in the plane for
  !                        each column)
  !
  !- columns with non zero ngkc first........
  !
  do mc = 1, ncplane
     if (ipc (mc) .gt.0) then
        j = ipc (mc)
        ncp_ (j) = ncp_ (j) + 1
        icpl (ncp_ (j) + ncp0 (j) ) = mc
        if (j.eq.me) then
           ipc (mc) = ncp_ (j)
        else
           ipc (mc) = 0
        endif
     endif
  enddo
  !
  !-..... ( intermediate check ) ....
  !
  do j = 1, nprocp
     if (ncp_ (j) .ne.nkcp (j) ) then
        write (6, * ) 'ncp_(j).ne.nkcp(j)', j, ncp_ (j) , nkcp (j)
     endif
  enddo
  !
  !- ........then the remaining columns
  !
  do mc = 1, ncplane
     if (ipc (mc) .lt.0) then
        j = - ipc (mc)
        ncp_ (j) = ncp_ (j) + 1
        icpl (ncp_ (j) + ncp0 (j) ) = mc
        if (j.eq.me) then
           ipc (mc) = ncp_ (j)
        else
           ipc (mc) = 0
        endif
     endif
  enddo
  !
  !-... ( final check )
  !
  nct_ = 0
  do j = 1, nprocp
     if (ncp_ (j) .ne.ncp (j) ) then
        write (6, * ) 'ncp_(j).ne.ncp(j)', j, ncp_ (j) , ncp (j)
     endif
     nct_ = nct_ + ncp_ (j)
  enddo
  if (nct_.ne.nct) then
     write (6, * ) 'nct_.ne.nct', nct_, nct
  endif
  !
  !   here compute the two arrays ipcs and ipcls for the smooth mesh
  !
  do j = 1, nprocp
     ncp_ (j) = 0
  enddo
  !
  !    the columns of the wavefunction sphere first
  !
  do mc = 1, ncplanes
     if (ipcs (mc) .gt.0) then
        j = ipcs (mc)
        ncp_ (j) = ncp_ (j) + 1
        icpls (ncp_ (j) + ncp0s (j) ) = mc
        if (j.eq.me) then
           ipcs (mc) = ncp_ (j)
        else
           ipcs (mc) = 0
        endif
     endif
  enddo
  !
  !    and then all the others
  !
  do mc = 1, ncplanes
     if (ipcs (mc) .lt.0) then
        j = - ipcs (mc)
        ncp_ (j) = ncp_ (j) + 1
        icpls (ncp_ (j) + ncp0s (j) ) = mc
        if (j.eq.me) then
           ipcs (mc) = ncp_ (j)
        else
           ipcs (mc) = 0
        endif
     endif

  enddo
  if (nprocp.eq.1) then
     nrxx = nrx1 * nrx2 * nrx3
     nrxxs = nrx1s * nrx2s * nrx3s
  else
     nrxx = max (nrx3 * ncp (me), nrx1 * nrx2 * npp (me) )
     nrxxs = max (nrx3s * ncps (me), nrx1s * nrx2s * npps (me) )
  endif
  !
  ! nxx is just a copy in the parallel commons of nrxx
  !
  nxx = nrxx

  nxxs = nrxxs

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

  lb(1) = -n1
  lb(2) = -n2
  lb(3) = -n3
  ub(1) =  n1
  ub(2) =  n2
  ub(3) =  n3
!
  ALLOCATE( stw ( lb(2):ub(2), lb(3):ub(3) ) )
  stw = 0

  do i1 = - n1, n1
     do i2 = - n2, n2
        do i3 = - n3, n3
           amod = (i1 * bg (1, 1) + i2 * bg (1, 2) + i3 * bg (1, 3) ) **2 + &
                (i1 * bg (2, 1) + i2 * bg (2, 2) + i3 * bg (2, 3) ) **2 + (i1 * &
                bg (3, 1) + i2 * bg (3, 2) + i3 * bg (3, 3) ) **2
           if (amod.le.gcutm) ngm = ngm + 1
           if (amod.le.gcutms) ngms = ngms + 1
           if (amod <= gkcut ) stw( i2, i3 ) = 1
        enddo
     enddo
  enddo

  call fft_dlay_scalar( dfftp, ub, lb, nr1, nr2, nr3, nrx1, nrx2, nrx3, stw )
  call fft_dlay_scalar( dffts, ub, lb, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, stw )

  deallocate( stw )

#endif

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

  return

contains

  subroutine calculate_gkcut()
    if (nks.eq.0) then
       !
       ! if k-points are automatically generated (which happens later)
       ! use max(bg)/2 as an estimate of the largest k-point
       !
       gkcut = sqrt (ecutwfc) / tpiba + 0.5d0 * max (sqrt (bg ( &
          1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2), sqrt (bg (1, 2) ** &
          2 + bg (2, 2) **2 + bg (3, 2) **2), sqrt (bg (1, 3) **2 + bg ( &
          2, 3) **2 + bg (3, 3) **2) )
    else
       gkcut = 0.0d0
       do kpoint = 1, nks
          gkcut = max (gkcut, sqrt (ecutwfc) / tpiba + sqrt (xk ( &
               1, kpoint) **2 + xk (2, kpoint) **2 + xk (3, kpoint) **2) )
       enddo
    endif
    gkcut = gkcut * gkcut
    return
  end subroutine calculate_gkcut


end subroutine data_structure

