!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine data_structure  
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays.
  ! In the parallel case distributes columns to processes, too
  ! This version computes also the smooth and hard mesh
  !
#include "machine.h"
  use pwcom  
#ifdef PARA
  use para
  use allocate
#endif
  use mp, only: mp_sum
  use mp_global, only: intra_pool_comm
  !
  implicit none
  integer :: n1, n2, n3, i1, i2, i3  
  ! counters on G space
  !
  integer :: good_fft_dimension  
  ! a function with obvious meaning

  real(kind=DP) :: amod  
  ! modulus of G vectors
#ifdef PARA
  ! counters on planes

  integer :: np, nps1, nq, nqs, max1, min1, max2, min2, kpoint, m1, &
       m2, i, mc, nct_ 
  integer, pointer :: ngc (:), ngcs (:), ngkc (:), &
       ic (:), ics (:)
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

  real(kind=DP) :: gkcut  
  ! cut-off for the wavefunctions
  !
  ! set the values of fft arrays
  !
  nrx1 = good_fft_dimension (nr1)  
  nrx1s = good_fft_dimension (nr1s)  
  !
  ! nrx2 is there just for compatibility
  !
  nrx2 = nr2  
  nrx2s = nr2s  
  nrx3 = good_fft_dimension (nr3)  
  nrx3s = good_fft_dimension (nr3s)  
  !
  !     compute number of columns per plane for each processor
  !
  ncplane = nrx1 * nrx2  

  ncplanes = nrx1s * nrx2s  
  !
  !    global pointers passed with commons are allocate here
  !
  call mallocate(ipc, ncplane)  
  call mallocate(ipcs, ncplanes)  
  !
  !   local pointers deallocated at the end
  !
  call mallocate(ngc,   ncplane)  
  call mallocate(ngcs,  ncplane)  
  call mallocate(ngkc,  ncplane)  
  call mallocate(ic,   ncplane)  
  call mallocate(ics,  ncplane)  
  !
  ! set the number of plane per process
  !
  if (nr3.lt.nprocp) call error ('data_structure', &
       'some processors have no planes ',  - 1)

  if (nr3s.lt.nprocp) call error ('data_structure', &
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
  do mc = 1, ncplane  
     ngc (mc) = 0  
     ngcs (mc) = 0  
     ngkc (mc) = 0  
     ipc (mc) = 0  
     ic (mc) = 0  
  enddo
  do mc = 1, ncplanes  
     ipcs (mc) = 0  
  enddo
  do i1 = - n1, n1  
     m1 = mod (i1, nr1) + 1  
     if (m1.lt.1) m1 = m1 + nr1  
     do i2 = - n2, n2  
        m2 = mod (i2, nr2) + 1  
        if (m2.lt.1) m2 = m2 + nr2  
        mc = m1 + (m2 - 1) * nrx1  
        if (mc.lt.1.or.mc.gt.ncplane) call error ('data_structure', &
             'mc is wrong', 1)
        do i3 = - n3, n3  
           amod = (bg (1, 1) * i1 + bg (1, 2) * i2 + bg (1, 3) * i3) **2 + &
                (bg (2, 1) * i1 + bg (2, 2) * i2 + bg (2, 3) * i3) **2 + (bg (3, &
                1) * i1 + bg (3, 2) * i2 + bg (3, 3) * i3) **2
           if (amod.le.gcutm) ngc (mc) = ngc (mc) + 1  
           if (amod.le.gcutms) ngcs (mc) = ngcs (mc) + 1  
           if (amod.le.gkcut) ngkc (mc) = ngkc (mc) + 1  
        enddo
     enddo
  enddo
  !
  !   Now compute the total number of nonzero columns
  !
  nct = 0  
  ncts = 0  
  do mc = 1, ncplane  
     if (ngcs (mc) .gt.0) ncts = ncts + 1  
     if (ngc (mc) .gt.0) then  
        nct = nct + 1  
        ngc (nct) = ngc (mc)  
        ngcs (nct) = ngcs (mc)  
        ngkc (nct) = ngkc (mc)  
        ic (nct) = mc  
        !           write (*,*) nct, ngc(nct), mc
     endif
  enddo
  do mc = nct + 1, ncplane  
     ngc (mc) = 0  
     ngcs (mc) = 0  
     ngkc (mc) = 0  

  enddo
  if (nct.eq.0) call error ('data_structure', 'number of column 0', 1)
  if (ncts.eq.0) call error ('data_structure', 'number smooth column 0', 1)
  !
  !   Now sort the columns. First the column with the largest number of G
  !   vectors on the wavefunction sphere, then on the smooth sphere,
  !   then on the big sphere
  !
  do i = 1, nct  
     do j = i + 1, nct  

        if ( (ngkc (i) .lt.ngkc (j) ) .or. (ngkc (i) .eq.ngkc (j) &
             .and.ngcs (i) .lt.ngcs (j) ) .or. (ngkc (i) .eq.ngkc (j) &
             .and.ngcs (i) .eq.ngcs (j) .and.ngc (i) .lt.ngc (j) ) ) then
           idum = ngkc (i)  
           ngkc (i) = ngkc (j)  

           ngkc (j) = idum  
           idum = ngc (i)  
           ngc (i) = ngc (j)  

           ngc (j) = idum  
           idum = ngcs (i)  
           ngcs (i) = ngcs (j)  

           ngcs (j) = idum  
           idum = ic (i)  
           ic (i) = ic (j)  

           ic (j) = idum  
        endif
     enddo
  enddo
  !
  !  Set the ics index with the correspondence between planes in the smooth
  !  thick mesh
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
  do i = 1, nct  
     if (ngcs (i) .gt.0) then  
        !
        ! this column contain smooth vectors, find the indices on the thick
        ! mesh i1,i2, and the corresponding indices on the smooth mesh n1,n2
        !
        i1 = mod (ic (i) - 1, nrx1) + 1  
        i2 = (ic (i) - i1) / nrx1 + 1  
        if (i1.le.max1) then  
           n1 = i1  
        elseif (i1.ge.min1) then  
           n1 = nr1s - (nr1 - i1)  
        else  
           call error ('data_structure', 'something wrong with n1', 1)  
        endif
        if (i2.le.max2) then  
           n2 = i2  
        elseif (i2.ge.min2) then  
           n2 = nr2s - (nr2 - i2)  
        else  
           call error ('data_structure', 'something wrong with n2', 1)  
        endif
        !
        !   check that the indices are within bounds
        !
        if (n1.lt.1.or.n1.gt.nr1s.or.n2.lt.1.or.n2.gt.nr2s) &
             call error ('data_structure', 'something wrong with n1,n2', 1)
        ics (i) = n1 + (n2 - 1) * nrx1s  
     else  
        ics (i) = 0  
     endif
  enddo
  !
  ! assign columns to processes, set the ipc index
  !
  do j = 1, nprocp  
     ncp (j) = 0  
     ncps (j) = 0  
     nkcp (j) = 0  
     ngp (j) = 0  
     ngps (j) = 0  
     ngkp (j) = 0  
  enddo
  do i = 1, nct  
     jj = 1  
     if (ngkc (i) .gt.0) then  
        do j = 1, nprocp  
           if (ngkp (j) .lt.ngkp (jj) ) jj = j  
        enddo
     else  
        do j = 1, nprocp  
           if (ngp (j) .lt.ngp (jj) ) jj = j  
        enddo
     endif
     ipc (ic (i) ) = - jj  
     ncp (jj) = ncp (jj) + 1  
     ngp (jj) = ngp (jj) + ngc (i)  
     if (ngcs (i) .gt.0) then  
        ncps (jj) = ncps (jj) + 1  
        ngps (jj) = ngps (jj) + ngcs (i)  
        ipcs (ics (i) ) = - jj  
     endif
     if (ngkc (i) .gt.0) then  
        ipcs (ics (i) ) = jj  
        ipc (ic (i) ) = jj  
        ngkp (jj) = ngkp (jj) + ngkc (i)  
        nkcp (jj) = nkcp (jj) + 1  
     endif
  enddo
  !
  call mfree (ics)  
  call mfree (ic)  
  call mfree (ngkc)  
  call mfree (ngcs)  
  call mfree (ngc)  
  !
  !   computing the starting column for each processor
  !   and the total number of G vectors
  !
  ngm = ngp (me)  
  ngms = ngps (me)  
  do i = 1, nprocp  
     !        write(6,*) i, ncp(i), ngp(i), nkcp(i), ngkp(i)
     if (ngkp (i) .eq.0) call error ('data_structure', &
          'some processors have no pencils, not yet implemented', 1)
     if (i.eq.1) then  
        ncp0 (i) = 0  
        ncp0s (i) = 0  
     else  
        ncp0 (i) = ncp0 (i - 1) + ncp (i - 1)  
        ncp0s (i) = ncp0s (i - 1) + ncps (i - 1)  
     endif


  enddo
  call mallocate(icpl, nct)  

  call mallocate(icpls, ncts)  
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
  !
  !     compute the number of g necessary to the calculation
  !
  n1 = nr1 + 1  
  n2 = nr2 + 1  

  n3 = nr3 + 1  
  ngm = 0  
  ngms = 0  
  do i1 = - n1, n1  
     do i2 = - n2, n2  
        do i3 = - n3, n3  
           amod = (i1 * bg (1, 1) + i2 * bg (1, 2) + i3 * bg (1, 3) ) **2 + &
                (i1 * bg (2, 1) + i2 * bg (2, 2) + i3 * bg (2, 3) ) **2 + (i1 * &
                bg (3, 1) + i2 * bg (3, 2) + i3 * bg (3, 3) ) **2
           if (amod.le.gcutm) ngm = ngm + 1  
           if (amod.le.gcutms) ngms = ngms + 1
        enddo
     enddo

  enddo
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
end subroutine data_structure

