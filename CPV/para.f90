!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "../include/machine.h"

module para_mod

  USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_allocate, fft_dlay_deallocate, &
                       fft_dlay_set

  integer maxproc, ncplanex
  parameter (maxproc=64, ncplanex=37000)
  
  character(len=3) :: node
! node:    node number, useful for opening files
  integer nproc, me, mygroup
! nproc:   number of processors
! me:      number of this processor
!
! parallel fft information for the dense grid
!
! npp:     number of plane per processor                      
! n3:      n3(me)+1 = first  plane on proc. me            
! ncp:     number of (density) columns per proc  
! ncp0:    starting column for each processor
! ncplane: number of columns in a plane                   
! nct:     total number of non-zero columns               
! nnr_:    local fft data size                            
! ipc:     index saying which proc owns columns in a plane
! icpl:    index relating columns and pos. in the plane   
!
! n3 -> dfftp%ipp
! ncplane -> dfftp%nnp
! ncp  -> dfftp%nsp
! ncp0 -> dfftp%iss
! npp  -> dfftp%npp
! ipc  -> dfftp%isind
! icpl -> dfftp%ismap
! nnr_ -> dfftp%nnr
!
!  integer  npp(maxproc), n3(maxproc), ncp(maxproc), ncp0(maxproc), &
!           ncplane, nct, nnr_, ipc(ncplanex), icpl(ncplanex)
!
! parallel fft information for the smooth mesh
!
! npps:    number of plane per processor
! ncps:    number of (density) columns per proc 
! ncpw:    number of (wfs) columns per processor
! ncps0:   starting column for each processor
! ncplanes:number of columns in a plane (smooth)
! ncts:    total number of non-zero columns
! nnrs_:   local fft data size
! ipcs:    saying which proc owns columns in a plane
! icpls:   index relating columns and pos. in the plane 
!
! ncpw -> dffts%ncpw
! n3s -> dffts%ipp
! ncplanes -> dffts%nnp
! ncps  -> dffts%nsp
! ncps0 -> dffts%iss
! npps  -> dffts%npp
! ipcs  -> dffts%isind
! icpls -> dffts%ismap
! nnrs_ -> dffts%nnr
!
!  integer npps(maxproc), ncps(maxproc), ncpw(maxproc), ncp0s(maxproc), &
!       ncplanes, ncts,  nnrs_, ipcs(ncplanex), icpls(ncplanex)

  TYPE ( fft_dlay_descriptor ) :: dfftp  ! fft descriptor for potentials
  TYPE ( fft_dlay_descriptor ) :: dffts  ! fft descriptor for smooth mesh

!  PRIVATE :: ipcs, icpls, ipc, icpl, ncp0s, ncp0, npp, npps, ncp, ncps, &
!    ncplane, ncplanes, n3, ncpw, nnr_, nnrs_, nct, ncts


end module para_mod
!
!-----------------------------------------------------------------------
      subroutine startup
!-----------------------------------------------------------------------
!
!  This subroutine initializes the Message Passing environment. 
!  The number of processors "nproc" returned by this routine is 
!  determined by the way the process was started. 
!  This is configuration- and machine-dependent. For inter-
!  active execution it is determined by environment variable MP_PROCS
!
      use para_mod
      use mp, only: mp_start, mp_env
!
      implicit none

      integer ierr
!
      call mp_start()
      call mp_env( nproc, me, mygroup )
      me = me + 1
!
!
! parent process (source) will have me=1 - child process me=2,...,NPROC
! (for historical reasons: MPI uses 0,...,NPROC-1 instead )
!
      if ( nproc > maxproc) &
     &   call errore('startup',' too many processors ',nproc)
!
      if (me.lt.10) then
         write(node,'(i1,2x)') me
      else if (me.lt.100) then
         write(node,'(i2,1x)') me
      else if (me.lt.1000) then
         write(node,'(i3)') me
      else
         call errore('startup','wow, >1000 nodes !!',nproc)
      end if
!
! only the first processor writes
!
      if ( me == 1 ) then
         write(6,'(/5x,''Parallel version (MPI)'')')
         write(6,'(5x,''Number of processors in use:   '',i4)') nproc
      else
         open(6,file='/dev/null',status='unknown')
!
! useful for debugging purposes 
!         open(6,file='out.'//node,status='unknown')
      end if
!
      return
      end

#if defined __PARA

!
!----------------------------------------------------------------------
      subroutine read_rho(unit,nspin,rhor)
!----------------------------------------------------------------------
!
! read from file rhor(nnr,nspin) on first node and distribute to other nodes
!
      use para_mod
      use parm
      implicit none
      include 'mpif.h'
      integer unit, nspin
      real(kind=8) rhor(nnr,nspin)
!
      integer ir, is
      integer root, proc, ierr, n, displs(nproc), sendcount(nproc)
      real(kind=8), allocatable:: rhodist(:)
!
!
      if (me.eq.1) allocate(rhodist(nr1x*nr2x*nr3x))
      root = 0
      do proc=1,nproc
         sendcount(proc) =  dfftp%nnp * ( dfftp%npp(proc) )
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + sendcount(proc-1)
         end if
      end do
      do is=1,nspin
!
! read the charge density from unit "unit" on first node only
!
         if (me.eq.1) read(unit) (rhodist(ir),ir=1,nr1x*nr2x*nr3x)
!
! distribute the charge density to the other nodes
!
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         call mpi_scatterv(rhodist, sendcount, displs, MPI_REAL8,       &
     &                     rhor(1,is),sendcount(me),   MPI_REAL8,       &
     &                     root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('mpi_scatterv','ierr<>0',ierr)
!
! just in case: set to zero unread elements (if any)
!
         do ir=sendcount(me)+1,nnr
            rhor(ir,is)=0.d0
         end do
      end do
      if (me.eq.1) deallocate(rhodist)
!
      return
      end subroutine read_rho
!
!----------------------------------------------------------------------
      subroutine write_rho(unit,nspin,rhor)
!----------------------------------------------------------------------
!
! collect rhor(nnr,nspin) on first node and write to file
!
      use para_mod
      use parm
      implicit none
      include 'mpif.h'
      integer unit, nspin
      real(kind=8) rhor(nnr,nspin)
!
      integer ir, is
      integer root, proc, ierr, displs(nproc), recvcount(nproc)
      real(kind=8), allocatable:: rhodist(:)
!
!
      if (me.eq.1) allocate(rhodist(nr1x*nr2x*nr3x))
!
      root = 0
      do proc=1,nproc
         recvcount(proc) =  dfftp%nnp  * ( dfftp%npp(proc) )
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + recvcount(proc-1)
         end if
      end do
!
      do is=1,nspin
!
! gather the charge density on the first node
!
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         call mpi_gatherv (rhor(1,is), recvcount(me), MPI_REAL8,        &
     &                     rhodist,recvcount, displs, MPI_REAL8,        &
     &                     root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('mpi_gatherv','ierr<>0',ierr)
!
! write the charge density to unit "unit" from first node only
!
         if (me.eq.1) write(unit) (rhodist(ir),ir=1,nr1x*nr2x*nr3x)
      end do
      if (me.eq.1) deallocate(rhodist)
!
      return
      end subroutine write_rho
!
!
!-----------------------------------------------------------------------
      subroutine set_fft_para( b1, b2, b3, gcut, gcuts, gcutw,          &
     &                        nr1, nr2, nr3, nr1s, nr2s, nr3s, nnr,     &
     &                        nr1x,nr2x,nr3x,nr1sx,nr2sx,nr3sx,nnrs )
!-----------------------------------------------------------------------
!
! distribute columns to processes for parallel fft
! based on code written by Stefano de Gironcoli for PWSCF
! columns are sets of g-vectors along z: g(k) = i1*b1+i2*b2+i3*b3 , 
! with g^2<gcut and (i1,i2) running over the (xy) plane.
! Columns are "active" for a given (i1,i2) if they contain a nonzero
! number of wavevectors
!
      use para_mod
      use stick_base
!
      implicit none
      real(kind=8) b1(3), b2(3), b3(3), gcut, gcuts, gcutw
      integer nr1, nr2, nr3, nr1s, nr2s, nr3s, nnr,                     &
     &        nr1x,nr2x,nr3x,nr1sx,nr2sx,nr3sx,nnrs
!
      integer ngc(ncplanex),&! number of g-vectors per column (dense grid)
     &        ngcs(ncplanex),&! number of g-vectors per column (smooth grid
     &        ngcw(ncplanex),&! number of wavefct plane waves per colum
     &        in1(ncplanex),&! index i for column (i1,i2)
     &        in2(ncplanex),&! index j for column (i1,i2)
     &        ic,           &! fft index for this column (dense grid)
     &        ics,          &! as above for the smooth grid
     &        icm,          &! fft index for column (-i1,-i2) (dense grid)
     &        icms,         &! as above for the smooth grid
     &        index(ncplanex),&! used to order column
     &        ncp_(maxproc),&! number of column per processor (work space)
     &        ngp(maxproc), &! number of g-vectors per proc (dense grid)
     &        ngps(maxproc),&! number of g-vectors per proc (smooth grid)
     &        ngpw(maxproc)  ! number of wavefct plane waves per proc
!
      integer np, nps1,            &! counters on planes 
     &        nq, nqs,             &! counters on planes
     &        max1,min1,max2,min2, &! aux. variables
     &        m1, m2, n1, n2, i, mc,&! generic counter
     &        idum, nct_,          &! check variables 
     &        j,jj,                &! counters on processors
     &        n1m1,n2m1,n3m1,      &! nr1-1 and so on
     &        i1, i2, i3,          &! counters on G space
     &        good_fft_dimension    ! a function with obvious meaning
      real(kind=8)                                                      &
     &        aux(ncplanex), &! used to order columns
     &        amod            ! modulus of G vectors
!
      integer  ncp0(maxproc), ipc(ncplanex), icpl(ncplanex)
      integer  ncp0s(maxproc), ipcs(ncplanex), icpls(ncplanex)
      integer  npp(maxproc), npps(maxproc)
      integer  ncp(maxproc), ncps(maxproc)
      integer  n3(maxproc)
      integer  ncplane, ncplanes
      integer  ncpw(maxproc)
      integer  nnr_, nnrs_, nct, ncts
!
      integer, allocatable :: st(:,:), stw(:,:), sts(:,:) ! sticks maps
      integer :: ub(3), lb(3)  ! upper and lower bounds for maps
      integer :: nctw
      logical :: tk = .FALSE.
!
      call tictac(27,0)
!
! set the dimensions of fft arrays
!
      nr1x =good_fft_dimension(nr1 )
      nr2x = nr2
      nr3x =good_fft_dimension(nr3 )
!
      nr1sx=good_fft_dimension(nr1s)
      nr2sx=nr2s
      nr3sx=good_fft_dimension(nr3s)
!
!     compute number of columns for each processor
!
      ncplane = nr1x*nr2x
      ncplanes= nr1sx*nr2sx
      if (ncplane.gt.ncplanex .or. ncplanes.gt.ncplanex)                &
     &     call errore('set_fft_para','ncplanex too small',ncplane)
!
! set the number of plane per process
!
      if (nr3.lt.nproc) call errore('set_fft_para',                      &
     &                'some processors have no planes ',-1)      
      if (nr3s.lt.nproc) call errore('set_fft_para',                     &
     &                'some processors have no smooth planes ',-1)      
!
      if (nproc.eq.1) then
         npp(1) = nr3
         npps(1)= nr3s
      else
         np = nr3/nproc
         nq = nr3 - np*nproc
         nps1 = nr3s/nproc
         nqs = nr3s - nps1*nproc
         do i = 1, nproc
            npp(i) = np
            if (i.le.nq) npp(i) = np + 1
            npps(i) = nps1
            if (i.le.nqs) npps(i) = nps1 + 1
         end do
      end if
      n3(1)= 0
      do i = 2, nproc
         n3(i)=n3(i-1)+npp(i-1)
      end do
!
!     Now compute for each point of the big plane how many column have
!     non zero vectors on the smooth and dense grid
!
      ipc  = 0
      ipcs = 0
!
      nct = 0
      ncts= 0
!
! NOTA BENE: the exact limits for a correctly sized FFT grid are:
! -nr/2,..,+nr/2  for nr even; -(nr-1)/2,..,+(nr-1)/2  for nr odd.
! If the following limits are increased, a slightly undersized fft 
! grid, with some degree of G-vector refolding, can be used
! (at your own risk - a check is done in ggen).
!
      n1m1=nr1/2
      n2m1=nr2/2
      n3m1=nr3/2

      lb(1) = -n1m1
      lb(2) = -n2m1
      lb(3) = -n3m1
      ub(1) =  n1m1
      ub(2) =  n2m1
      ub(3) =  n3m1
!
      ALLOCATE( stw ( lb(1):ub(1), lb(2):ub(2) ) )
      ALLOCATE( st  ( lb(1):ub(1), lb(2):ub(2) ) )
      ALLOCATE( sts ( lb(1):ub(1), lb(2):ub(2) ) )

!
! ...     Fill in the stick maps, for given g-space base (b1,b2,b3)
! ...     and cut-offs
! ...     The value of the element (i,j) of the map ( st ) is equal to the
! ...     number of G-vector belonging to the (i,j) stick.
!

      CALL sticks_maps( tk, ub, lb, b1, b2, b3, gcut, gcutw, gcuts, st, stw, sts )

      nct  = COUNT( st  > 0 )
      nctw = COUNT( stw > 0 )
      ncts = COUNT( sts > 0 )

      if ( nct > ncplane )    &
     &    call errore('set_fft_para','too many sticks',1)

      if ( ncts > ncplanes )  &
     &    call errore('set_fft_para','too many sticks',2)

      if ( nct == 0 ) & 
     &    call errore('set_fft_para','number of sticks 0', 1)

      if ( ncts == 0 ) &
     &    call errore('set_fft_para','number smooth sticks 0', 1)

!
! ...     initialize the sticks indexes array ist
! ...     nct counts columns containing G-vectors for the dense grid
! ...     ncts counts columns contaning G-vectors for the smooth grid
!

      CALL sticks_countg( tk, ub, lb, st, stw, sts, in1, in2, ngc, ngcw, ngcs )

!
!   Sort the columns. First the column with the largest number of G
!   vectors on the wavefunction sphere (active columns), 
!   then on the smooth sphere, then on the big sphere. Dirty trick:
!

      CALL sticks_sort( ngc, ngcw, ngcs, nct, index )

!
!   assign columns to processes
!

      CALL sticks_dist( tk, ub, lb, index, in1, in2, ngc, ngcw, ngcs, nct, &
                ncp, ncpw, ncps, ngp, ngpw, ngps, st, stw, sts )

      CALL sticks_pairup( tk, ub, lb, index, in1, in2, ngc, ngcw, ngcs, nct, &
                ncp, ncpw, ncps, ngp, ngpw, ngps, st, stw, sts )

      ! WRITE(*,*) ' DEBUG : ', COUNT( st > 0 ), COUNT( sts > 0 ), COUNT( stw > 0 )


      CALL fft_dlay_allocate( dfftp, nproc, nr1x, nr2x )
      CALL fft_dlay_allocate( dffts, nproc, nr1sx, nr2sx )


      CALL fft_dlay_set( dfftp, tk, nct, nr1, nr2, nr3, nr1x, nr2x, nr3x, me, &
                nproc, ub, lb, index, in1, in2, ncp, ncpw, ngp, ngpw, st, stw)
      CALL fft_dlay_set( dffts, tk, ncts, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, me, &
                nproc, ub, lb, index, in1, in2, ncps, ncpw, ngps, ngpw, sts, stw)


      DO mc = 1, nct

         i = index(mc)

         i1 = in1(i)
         i2 = in2(i)

         if ( i1.lt.0.or.(i1.eq.0.and.i2.lt.0) ) go to 29
!
! only half of the columns, plus column (0,0), are scanned:
! column (-i1,-i2) must be assigned to the same proc as column (i1,i2)
!
! ic  :  position, in fft notation, in dense grid, of column ( i1, i2)
! icm :      "         "      "          "    "         "    (-i1,-i2)
! ics :      "         "      "        smooth "         "    ( i1, i2)
! icms:      "         "      "        smooth "         "    (-i1,-i2)
!
         m1 = i1 + 1
         if (m1.lt.1) m1 = m1 + nr1
         m2 = i2 + 1
         if (m2.lt.1) m2 = m2 + nr2
         ic = m1 + (m2-1)*nr1x
!
         n1 = -i1 + 1
         if (n1.lt.1) n1 = n1 + nr1
         n2 = -i2 + 1
         if (n2.lt.1) n2 = n2 + nr2
         icm = n1 + (n2-1)*nr1x
!
         m1 = i1 + 1
         if (m1.lt.1) m1 = m1 + nr1s
         m2 = i2 + 1
         if (m2.lt.1) m2 = m2 + nr2s
         ics = m1 + (m2-1)*nr1sx
!
         n1 =-i1 + 1
         if (n1.lt.1) n1 = n1 + nr1s
         n2 =-i2 + 1
         if (n2.lt.1) n2 = n2 + nr2s
         icms = n1 + (n2-1)*nr1sx

         IF( st( i1, i2 ) > 0 .AND. stw( i1, i2 ) > 0 ) THEN
           ipc( ic  ) = st(  i1,  i2 )
           ipc( icm ) = st( -i1, -i2 )
         ELSE IF( st( i1, i2 ) > 0 ) THEN
           ipc( ic  ) = -st(  i1,  i2 )
           ipc( icm ) = -st( -i1, -i2 )
         END IF

         IF( sts( i1, i2 ) > 0 .AND. stw( i1, i2 ) > 0 ) THEN
           ipcs( ics  ) = sts(  i1,  i2 )
           ipcs( icms ) = sts( -i1, -i2 )
         ELSE IF( sts( i1, i2 ) > 0 ) THEN
           ipcs( ics  ) = -sts(  i1,  i2 )
           ipcs( icms ) = -sts( -i1, -i2 )
         END IF

29       CONTINUE

      END DO
      
      DEALLOCATE( st, stw, sts )

      IF( .NOT. tk ) THEN
        nct  = nct*2  - 1
        ncts = ncts*2 - 1
      END IF

!
! ipc  is the processor for this column in the dense grid
! ipcs is the same, for the smooth grid
!
      write(6,'(                                                        &
     & '' Proc  planes cols    G   planes cols    G    columns  G''/    &
     & ''       (dense grid)     (smooth grid)   (wavefct grid)'')')
      do i=1,nproc
         write(6,'(i3,2x,3(i5,2i7))') i, npp(i),ncp(i),ngp(i),          &
     &        npps(i),ncps(i),ngps(i), ncpw(i), ngpw(i)
      end do
      write(6,'(i3,2x,3(i5,2i7))') 0, SUM(npp(1:nproc)), SUM(ncp(1:nproc)), &
        SUM(ngp(1:nproc)), SUM(npps(1:nproc)), SUM(ncps(1:nproc)), &
        SUM(ngps(1:nproc)), SUM(ncpw(1:nproc)), SUM(ngpw(1:nproc))
!
! nnr_ and nnrs_ are copies of nnr and nnrs, the local fft data size,
! to be stored in "parallel" commons. Not a very elegant solution.
!
      if ( nproc == 1 ) then
         nnr =nr1x*nr2x*nr3x
         nnrs=nr1sx*nr2sx*nr3sx
      else
         nnr = max(nr3x*ncp(me), nr1x*nr2x*npp(me))
         nnrs= max(nr3sx*ncps(me), nr1sx*nr2sx*npps(me))
      end if
      nnr_= nnr
      nnrs_= nnrs
!
!   computing the starting column for each processor
!
      do i=1,nproc
         if(ngpw(i).eq.0)                                               &
     &        call errore('set_fft_para',                                &
     &        'some processors have no pencils, not yet implemented',1)
         if (i.eq.1) then 
            ncp0(i) = 0
            ncp0s(i)= 0
         else
            ncp0(i) = ncp0 (i-1) + ncp (i-1)
            ncp0s(i)= ncp0s(i-1) + ncps(i-1)
         endif
      enddo
!
!  Now compute the arrays ipc and icpl (dense grid):
!     ipc contain the number of the column for that processor.
!         zero if the column do not belong to the processor.
!         Note that input ipc is used and overwritten.
!     icpl contains the point in the plane for each column
!
!- active columns first........
!
      do j=1,nproc
         ncp_(j) = 0
      end do
      do mc =1,ncplane
         if (ipc(mc).gt.0) then
            j = ipc(mc)
            ncp_(j) = ncp_(j) + 1
            icpl(ncp_(j) + ncp0(j)) = mc 
            if (j.eq.me) then
               ipc(mc) = ncp_(j)
            else
               ipc(mc) = 0
            end if 
         end if
      end do
!
!-..... ( intermediate check ) ....
!
      do j=1,nproc
         if (ncp_(j).ne.ncpw(j))                                        &
     &        call errore('set_fft_para','ncp_(j).ne.ncpw(j)',j)
      end do
!
!- ........then the remaining columns
!
      do mc =1,ncplane
         if (ipc(mc).lt.0) then
            j = -ipc(mc)
            ncp_(j) = ncp_(j) + 1
            icpl(ncp_(j) + ncp0(j)) = mc 
            if (j.eq.me) then
               ipc(mc) = ncp_(j)
            else
               ipc(mc) = 0
            end if 
         end if
      end do
!
!-... ( final check )
!
      nct_ = 0
      do j=1,nproc
         if (ncp_(j).ne.ncp(j))                                         &
     &        call errore('set_fft_para','ncp_(j).ne.ncp(j)',j)
         nct_ = nct_ + ncp_(j)
      end do
      if (nct_.ne.nct)                                                  &
     &     call errore('set_fft_para','nct_.ne.nct',1)
!
!   now compute the arrays ipcs and icpls 
!   (as ipc and icpls, for the smooth grid)
!
!   active columns first...
!
      do j=1,nproc
         ncp_(j) = 0
      end do
      do mc =1,ncplanes
         if (ipcs(mc).gt.0) then
            j = ipcs(mc)
            ncp_(j)=ncp_(j) + 1
            icpls(ncp_(j) + ncp0s(j)) = mc
            if (j.eq.me) then
               ipcs(mc) = ncp_(j)
            else
               ipcs(mc) = 0
            endif
         endif
      enddo
!
!-..... ( intermediate check ) ....
!
      do j=1,nproc
         if (ncp_(j).ne.ncpw(j))                                        &
     &        call errore('set_fft_para','ncp_(j).ne.ncpw(j)',j)
      end do
!
!    and then all the others
!
      do mc =1,ncplanes
         if (ipcs(mc).lt.0) then
            j = -ipcs(mc)
            ncp_(j) = ncp_(j) + 1
            icpls(ncp_(j) + ncp0s(j)) = mc
            if (j.eq.me) then
               ipcs(mc) = ncp_(j)
            else
               ipcs(mc) = 0
            end if
         end if
      end do
!
!-... ( final check )
!
      nct_ = 0
      do j=1,nproc
         if (ncp_(j).ne.ncps(j))                                        &
     &        call errore('set_fft_para','ncp_(j).ne.ncps(j)',j)
         nct_ = nct_ + ncp_(j)
      end do
      if (nct_.ne.ncts)                                                 &
     &     call errore('set_fft_para','nct_.ne.ncts',1)
      call tictac(27,1)

!      do i = 1, nproc
!        write(6,fmt="('DEBUG fft_setup ',3I5 )" ) i, npp(i), dfftp%npp(i)
!        write(6,fmt="('DEBUG fft_setup ',3I5 )" ) i, npps(i), dffts%npp(i)
!      end do
!      write(6,fmt="('DEBUG fft_setup ',3I9 )" ) nnr_, dfftp%nnr
!      write(6,fmt="('DEBUG fft_setup ',3I9 )" ) nnrs_, dffts%nnr
!      write(6,fmt="('DEBUG fft_setup ',3I9 )" ) nct, dfftp%nst
!      write(6,fmt="('DEBUG fft_setup ',3I9 )" ) ncts, dffts%nst
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine cfftp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,sign)
!----------------------------------------------------------------------
!
!   sign = +-1 : parallel 3d fft for rho and for the potential
!
!   sign = +1 : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!               fft along z using pencils (cft_1)
!               transpose across nodes    (fft_scatter)
!                  and reorder
!               fft along y and x         (cft_2)
!   sign = -1 : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega 
!               fft along x and y         (cft_2)
!               transpose across nodes    (fft_scatter)
!                  and reorder
!               fft along z using pencils (cft_1)
!
!   based on code written by Stefano de Gironcoli for PWSCF
!
      use para_mod
      use work_fft
      use fft_base, only: fft_scatter
!
      implicit none
      integer nr1,nr2,nr3,nr1x,nr2x,nr3x, sign, nppx
      complex(kind=8)  f( dfftp%nnr )
      integer  mc, i, j, ii
!
! the following is needed if the fft is distributed over only one processor
! for the special case nx3.ne.n3. Not an elegant solution, but simple, fast,
! and better than the preceding one that did not work in some cases. Note 
! that fft_scatter does nothing if nproc=1. PG
!

      if ( nr1  /= dfftp%nr1  ) call errore(' cfftp ',' wrong dims ', 1)
      if ( nr2  /= dfftp%nr2  ) call errore(' cfftp ',' wrong dims ', 2)
      if ( nr3  /= dfftp%nr3  ) call errore(' cfftp ',' wrong dims ', 3)
      if ( nr1x /= dfftp%nr1x ) call errore(' cfftp ',' wrong dims ', 4)
      if ( nr2x /= dfftp%nr2x ) call errore(' cfftp ',' wrong dims ', 5)
      if ( nr3x /= dfftp%nr3x ) call errore(' cfftp ',' wrong dims ', 6)

      if (nproc.eq.1) then
         nppx=nr3x
      else
         nppx=dfftp%npp(me)
      end if

      if (sign.eq.1) then

         call cft_1(f,dfftp%nsp(me),nr3,nr3x,sign,aux)
         call fft_scatter(aux,nr3x,dfftp%nnr,f,dfftp%nsp,dfftp%npp,sign)
         call zero(2*dfftp%nnr,f)
         do i=1,dfftp%nst
            mc = dfftp%ismap( i )
            do j=1,dfftp%npp(me)
               f(mc+(j-1)*dfftp%nnp) = aux(j + (i-1)*nppx)
            end do
         end do
!
         call cft_2(f,dfftp%npp(me),nr1,nr2,nr1x,nr2x,sign)
!
      else if (sign.eq.-1) then
!
         call cft_2(f,dfftp%npp(me),nr1,nr2,nr1x,nr2x,sign)
!
         do i=1,dfftp%nst
            mc = dfftp%ismap( i )
            do j=1,dfftp%npp(me)
               aux(j + (i-1)*nppx) = f(mc+(j-1)*dfftp%nnp)
            end do
         end do
         call fft_scatter(aux,nr3x,dfftp%nnr,f,dfftp%nsp,dfftp%npp,sign)
         call cft_1(aux,dfftp%nsp(me),nr3,nr3x,sign,f)
      else
          call errore('cftp','not allowed',abs(sign))
      end if
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine cfftps (f,nr1,nr2,nr3,nr1x,nr2x,nr3x,sign)
!----------------------------------------------------------------------
!
!   sign = +-1 : parallel 3d fft for rho and for the potential
!   sign = +-2 : parallel 3d fft for wavefunctions
!
!   sign = + : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!              fft along z using pencils        (cft_1)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along y (using planes) and x (cft_2)
!   sign = - : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega 
!              fft along x and y(using planes)  (cft_2)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along z using pencils        (cft_1)
!
!   The array "planes" signals whether a fft is needed along y :
!     planes(i)=0 : column f(i,*,*) empty , don't do fft along y
!     planes(i)=1 : column f(i,*,*) filled, fft along y needed
!   "empty" = no active components are present in f(i,*,*) 
!             after (sign>0) or before (sign<0) the fft on z direction
!
!   Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
!   and all planes(i) are set to 1
!
!   based on code written by Stefano de Gironcoli for PWSCF
!
      use para_mod
      use work_fft
      use fft_base, only: fft_scatter
!
      implicit none
      integer nr1,nr2,nr3,nr1x,nr2x,nr3x,sign
      complex(kind=8)  f( dffts%nnr )
      integer  mc, i, j, ii, proc, k, nppx
      integer planes(nr1x)
!
! see comments in cfftp for the logic (or lack of it) of the following
!
      if ( nr1  /= dffts%nr1  ) call errore(' cfftps ',' wrong dims ', 1)
      if ( nr2  /= dffts%nr2  ) call errore(' cfftps ',' wrong dims ', 2)
      if ( nr3  /= dffts%nr3  ) call errore(' cfftps ',' wrong dims ', 3)
      if ( nr1x /= dffts%nr1x ) call errore(' cfftps ',' wrong dims ', 4)
      if ( nr2x /= dffts%nr2x ) call errore(' cfftps ',' wrong dims ', 5)
      if ( nr3x /= dffts%nr3x ) call errore(' cfftps ',' wrong dims ', 6)

      if ( nproc == 1 ) then
         nppx = dffts%nr3x
      else
         nppx = dffts%npp(me)
      end if

      if ( sign > 0 ) then
         if ( sign /= 2 ) then
            call cft_1s(f,dffts%nsp(me),nr3,nr3x,sign,aux)
            call fft_scatter( aux, nr3x, dffts%nnr, f, dffts%nsp, dffts%npp, sign)
            call zero(2*dffts%nnr,f)
            do i = 1, dffts%nst
               mc = dffts%ismap( i )
               do j = 1, dffts%npp(me)
                  f( mc + (j-1) * dffts%nnp ) = aux( j + (i-1) * nppx)
               end do
            end do
            planes = dffts%iplp
         else
            call cft_1s(f,dffts%nsw(me),nr3,nr3x,sign,aux)
            call fft_scatter( aux, nr3x, dffts%nnr, f, dffts%nsw, dffts%npp, sign)
            call zero( 2*dffts%nnr, f )
            ii = 0
            do proc=1,nproc
               do i=1,dffts%nsw(proc)
                  mc = dffts%ismap( i + dffts%iss(proc) )
                  ii = ii + 1 
                  do j=1,dffts%npp(me)
                     f(mc+(j-1)*dffts%nnp) = aux(j + (ii-1)*nppx)
                  end do
               end do
            end do
            planes = dffts%iplw
         end if
!
         call cft_2s(f,dffts%npp(me),nr1,nr2,nr1x,nr2x,sign,planes)
!
      else
!
         if (sign.ne.-2) then
            planes = dffts%iplp
         else
            planes = dffts%iplw
         endif
!
         call cft_2s(f,dffts%npp(me),nr1,nr2,nr1x,nr2x,sign,planes)
!
         if (sign.ne.-2) then
            do i=1,dffts%nst
               mc = dffts%ismap( i )
               do j=1,dffts%npp(me)
                  aux(j + (i-1)*nppx) = f(mc+(j-1)*dffts%nnp)
               end do
            end do
            call fft_scatter(aux,nr3x,dffts%nnr,f,dffts%nsp,dffts%npp,sign)
            call cft_1s(aux,dffts%nsp(me),nr3,nr3x,sign,f)
         else
            ii = 0
            do proc=1,nproc
               do i=1,dffts%nsw(proc)
                  mc = dffts%ismap( i + dffts%iss(proc) )
                  ii = ii + 1 
                  do j=1,dffts%npp(me)
                     aux(j + (ii-1)*nppx) = f(mc+(j-1)*dffts%nnp)
                  end do
               end do
            end do
            call fft_scatter(aux,nr3x,dffts%nnr,f,dffts%nsw,dffts%npp,sign)
            call cft_1s(aux,dffts%nsw(me),nr3,nr3x,sign,f)
         end if
      end if
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine cfftpb(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,sign)
!----------------------------------------------------------------------
!
!   not-so-parallel 3d fft for box grid, implemented only for sign=1
!   G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!   The array f (overwritten on output) is NOT distributed:
!   a copy is present on each processor.
!   The fft along z  is done on the entire grid.
!   The fft along xy is done only on planes that have components on the
!   dense grid for each processor. Note that the final array will no
!   longer be the same on all processors.
!
      use para_mod
      use parm, only: nr3
!
      implicit none
      integer nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,sign
      complex(kind=8) f(nr1bx,nr2bx,nr3bx)
!
      integer ir3, ibig3, imin3, imax3, np3
!
      call parabox(nr3b,irb3,nr3,imin3,imax3)
      np3=imax3-imin3+1
! np3 is the number of planes to be transformed
      if (np3.le.0) return
      call cft_b(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,imin3,imax3,sign)
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine parabox(nr3b,irb3,nr3,imin3,imax3)
!----------------------------------------------------------------------
!
! find if box grid planes in the z direction have component on the dense
! grid on this processor, and if, which range imin3-imax3
!
      use para_mod
! input
      integer nr3b,irb3,nr3
! output
      integer imin3,imax3
! local
      integer ir3, ibig3
!
      imin3=nr3b
      imax3=1
      do ir3=1,nr3b
         ibig3=1+mod(irb3+ir3-2,nr3)
         if(ibig3.lt.1.or.ibig3.gt.nr3)                                 &
     &        call errore('cfftpb','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         if (ibig3.gt.0.and.ibig3.le.dfftp%npp(me)) then
            imin3=min(imin3,ir3)
            imax3=max(imax3,ir3)
         end if
      end do
!
      return
      end
!
!------------------------------------------------------------------------
      subroutine fft_scatter2                                           &
     &     (nproc, me, f_in, nr3x, nnr_, f_aux, ncp_, npp_, sign)
!------------------------------------------------------------------------
!
! transpose the fft grid across nodes
! a) From columns to planes (sign > 0)
!
!    "columns" (or "pencil") representation:
!    processor "me" has ncp_(me) contiguous columns along z
!    Each column is nr3x=nr3+1 elements for a fft of order nr3
!    (the additional element is added to reduce memory conflicts)
!
!    The transpose take places in two steps:
!    1) on each processor the columns are divided into slices along z
!       that are stored contiguously. On processor "me", slices for 
!       processor "proc" are npp_(proc)*ncp_(me) big
!    2) all processors communicate to exchange slices
!       (all columns with z in the slice belonging to "me"
!        must be received, all the others must be sent to "proc")
!    Finally one gets the "planes" representation:
!    processor "me" has npp_(me) complete xy planes
!
!  b) From planes to columns (sign < 0)
!
!  Quite the same in the opposite direction
!
!  The output is overwritten on f_in ; f_aux is used as work space
!
!  based on code written by Stefano de Gironcoli for PWSCF
!
      use para_mod, only: maxproc
      implicit none
!
      include 'mpif.h'
      integer  nproc, me, nr3x, nnr_, sign, ncp_(nproc), npp_(nproc)
      real(kind=8)     f_in(2*nnr_), f_aux(2*nnr_)
!
      integer  dest, from, k, offset1(maxproc),                         &
     &         sendcount(maxproc), sdispls(maxproc),                    &
     &         recvcount(maxproc), rdispls(maxproc),                    &
     &         proc, ierr
!
!
      if (nproc.eq.1) return
      call tictac(28,0)
!
! sendcount(proc): amount of data processor "me" must send to processor proc
! recvcount(proc): amount of data processor "me" must receive from      proc
!
      do proc = 1, nproc
         sendcount(proc) = 2*npp_(proc)*ncp_(me)
         recvcount(proc) = 2*npp_(me)*ncp_(proc)
      end do
!
! offset1(proc) is used to locate the slices to be sent to proc
! sdispls(proc)+1 is the beginning of data that must be sent to proc
! rdispls(proc)+1 is the beginning of data that must be received from proc
!
      offset1(1) = 1
      sdispls(1)=0
      rdispls(1)=0
      do proc = 2, nproc
         offset1(proc) = offset1(proc-1) + 2 * npp_(proc-1)
         sdispls(proc) = sdispls(proc-1) + sendcount(proc-1)
         rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
      end do
!
      if(sign.gt.0) then
!
! "forward" scatter from columns to planes
!
! step one: store contiguously the slices
!
         do proc = 1, nproc
            from = offset1(proc)
            dest = 1 + sdispls(proc)
            do k = 1, ncp_(me)
               call DCOPY ( 2 * npp_(proc),                              &
     &              f_in ( from + 2*(k-1)*nr3x) , 1,                    &
     &              f_aux( dest + 2*(k-1)*npp_(proc) ) , 1 )
            end do
         end do
!
! maybe useless; ensures that no garbage is present in the output
!
         call zero (2*nnr_,f_in)
!
! step two: communication
!
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         call mpi_alltoallv(f_aux,sendcount,sdispls,MPI_REAL8,          &
     &                      f_in ,recvcount,rdispls,MPI_REAL8,          &
     &                      MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('fft_scatter','ierr<>0',ierr)
!
      else
!
!  step two: communication
!
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         call mpi_alltoallv(f_in ,recvcount,rdispls,MPI_REAL8,          &
     &                      f_aux,sendcount,sdispls,MPI_REAL8,          &
     &                      MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('fft_scatter','ierr<>0',ierr)
!
!  step one: store contiguously the columns
!
         call zero (2*nnr_,f_in)
!
         do proc = 1, nproc
            from = 1 + sdispls(proc)
            dest = offset1(proc)
            do k = 1, ncp_(me)
               call DCOPY ( 2 * npp_(proc),                              &
     &              f_aux( from + 2*(k-1)*npp_(proc) ) , 1 ,            &
     &              f_in ( dest + 2*(k-1)*nr3x) , 1 )
            end do
         end do
!
      end if
      call tictac(28,1)
!
      return
      end 
!
!-----------------------------------------------------------------------
      subroutine reduce(size,ps)
!-----------------------------------------------------------------------
!
!     sums a distributed variable s(size) over the processors.
!     This version uses a fixed-length buffer of appropriate (?) size
!
      use para_mod
!
      implicit none
      integer size
      real(kind=8)  ps(size)
!
      include 'mpif.h'
      integer ierr, n, nbuf
      integer, parameter:: MAXB=10000
      real(kind=8) buff(MAXB)
!
      if (nproc.le.1) return
      if (size.le.0) return
      call tictac(29,0)
!
!  syncronize processes
!
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) call errore('reduce','error in barrier',ierr)
!
      nbuf=size/MAXB
!
      do n=1,nbuf
         call mpi_allreduce (ps(1+(n-1)*MAXB), buff, MAXB,              &
     &        MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0)                                                 &
     &        call errore('reduce','error in allreduce1',ierr)
         call DCOPY(MAXB,buff,1,ps(1+(n-1)*MAXB),1)
      end do
!
!    possible remaining elements < maxb
!
      if (size-nbuf*MAXB.gt.0) then
          call mpi_allreduce (ps(1+nbuf*MAXB), buff, size-nbuf*MAXB,    &
     &          MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
          if (ierr.ne.0)                                                &
     &         call errore('reduce','error in allreduce2',ierr)
          call DCOPY(size-nbuf*MAXB,buff,1,ps(1+nbuf*MAXB),1)
      endif
      call tictac(29,1)
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine print_para_times
!----------------------------------------------------------------------
!
      use para_mod
      use timex_mod
!
      implicit none
      include 'mpif.h'
      integer i, ierr
      real(kind=8) mincpu(maxclock), maxcpu(maxclock)
!
!  syncronize processes
!
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if (ierr.ne.0)                                                    &
     &     call errore('print_all_times','error in barrier',ierr)
!
      call mpi_allreduce (cputime, mincpu, maxclock,                    &
     &     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
      if (ierr.ne.0)                                                    &
     &     call errore('print_para_times','error in minimum',ierr)
      call mpi_allreduce (cputime, maxcpu, maxclock,                    &
     &     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (ierr.ne.0)                                                    &
     &     call errore('print_para_times','error in maximum',ierr)
!
      write(6,*)
      write(6,*) ' routine     calls       cpu time        elapsed'
      write(6,*) '             node0  node0,  min,  max     node0'
      write(6,*)
      do i=1, maxclock
         if (ntimes(i).gt.0) write(6,30) routine(i),                    &
     &        ntimes(i), cputime(i), mincpu(i),maxcpu(i), elapsed(i)
      end do
 30   format(a10,i7,4f8.1)
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine write_wfc(unit,c)
!----------------------------------------------------------------------
!
! collect wavefunctions on first node and write to file
!
      use gvec
      use gvecs
      use elct
      use para_mod
      use parms
      use work1
!
      implicit none
      include 'mpif.h'
      integer unit
      complex(kind=8) c(ngw,nx)
      complex(kind=8), pointer:: psis(:)
!
      integer i, ii, ig, proc, ierr, ntot, ncol, mc
      integer nmin(3), nmax(3), n1,n2,nzx,nz,nz_
      integer root, displs(nproc), recvcount(nproc)
      complex(kind=8), allocatable:: psitot(:), psiwr(:,:,:)
!
! nmin, nmax are the bounds on (i,j,k) indexes of wavefunction G-vectors
!
      call nrbounds(ngw,nr1s,nr2s,nr3s,in1p,in2p,in3p,nmin,nmax)
!
! nzx is the maximum length of a column along z
!
      nzx=nmax(3)-nmin(3)+1
!
      root = 0
! root is the first node
      ntot = 0
      do proc=1,nproc
         recvcount(proc) = dffts%nsw(proc)*nzx
!
! recvcount(proc) = size of data received from processor proc
!                   (number of columns times length of each column)
!
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + recvcount(proc-1)
         end if
!
! displs(proc) is the position of data received from processor proc
!
         ntot = ntot + recvcount(proc)
!
! ntot = total size of gathered data
!
      end do
!
! allocate the needed work spaces
!
      psis=> wrk1
      call zero(2*nnrs,psis)
      if (me.eq.1) then
         allocate(psitot(ntot))
         allocate(psiwr(nmin(3):nmax(3),nmin(1):nmax(1),nmin(2):nmax(2)))
         write(unit) n, nmin, nmax
      end if
!
! fill array psis with c_i(G) (as packed columns along z)
!
      do i=1,n
         do ig=1,ngw
!
! ncol+1 is the index of the column            
!
            ncol=(nps(ig)-1)/nr3sx
!
! nz_ is the z component in FFT style (refolded between 1 and nr3s)
!
            nz_ =nps(ig)-ncol*nr3sx
!
! nz is the z component in "natural" style (between nmin(3) and nmax(3))
!
            nz  =nz_-1
            if (nz.ge.nr3s/2) nz=nz-nr3s
!
! dffts%nsw(me) columns along z are stored in contiguous order on each node 
!
            psis(nz-nmin(3)+1+ncol*nzx)=c(ig,i)
         end do
!
! gather all psis arrays on the first node, in psitot
!
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('write_wfc','mpi_barrier 2',ierr)
         call mpi_gatherv (psis, recvcount(me),     MPI_DOUBLE_COMPLEX, &
    &                      psitot,recvcount, displs,MPI_DOUBLE_COMPLEX, &
    &                 root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('write_wfc','mpi_gatherv',ierr)
!
! now order the columns with a node-number-independent ordering
! We use n1,n2,nz for ordering, where G=n1*g(1)+n2*g(2)+nz*g(3)
! and g(1), g(2), g(3) are basis vectors of the reciprocal lattice.
! Note that the entire column is written, even if there are no
! wavefunction components associated to it. This is a waste of 
! disk space but it is simpler to implement than other schemes.
!
         if (me.eq.1) then
            ncol=0
            do proc=1,nproc
               do ii=1,dffts%nsw(proc)
                  ncol=ncol+1
                  ! mc=icpls(ii+ncp0s(proc))
                  mc = dffts%ismap( ii + dffts%iss(proc) )
!
! mc is the position in the xy plane of this column in FFT style
! we need to calculate "natural" indexes n1,n2, centered on the 
! origin, for the smallest possible box.
!
                  n2=(mc-1)/nr1sx
                  n1= mc-1-n2*nr1sx
                  if (n2.ge.nr2s/2) n2=n2-nr2s
                  if (n1.le.nmax(1).and.(n1.ne.0.or.n2.ge.0)) then
!
! NB: n1.gt.nmax(1) correspond to negative indexes that should be
!     refolded into negative n1. However these are absent because
!     only half of the G sphere is considered. The other condition
!     excludes the case n1=0, n2<0 that is also not present 
!
                     do nz=nmin(3),nmax(3)
                        psiwr(nz,n1,n2)=psitot(nz-nmin(3)+1+(ncol-1)*nzx)
                     end do
                  end if
               end do
            end do
!
! write the node-number-independent array 
!
            write(unit) psiwr
!
         end if
      end do
!
      if (me.eq.1) then
         deallocate(psiwr)
         deallocate(psitot)
      end if
!
      return
!
      end subroutine write_wfc
!
!----------------------------------------------------------------------
      subroutine write_wfc2(unit,c)
!----------------------------------------------------------------------
!
! collect wavefunctions on first node and write to file
! This subroutine differs from write_wfc, since it writes
! to the disks only the wavefunctions components less than
! the cut-off, and following a given order, independent
! from the number of processors.
! Each wave function is collected (with the call to mergewf) on the
! master node, according to the array ig_l2g, and then written to 
! the disk
!
      use gvec
      use gvecs
      use elct
      use para_mod
      use parms
      use work1
      use mp_wave, ONLY: mergewf
      use mp, ONLY: mp_sum
      use parallel_include
!
      implicit none
      integer unit
      complex(kind=8) c(ngw,nx)
      complex(kind=8), allocatable :: ctot(:)
!
      integer i, ig, mpime, proc, ierr, ntot, ncol, mc, ngwt, root

      ngwt = ngw
      root = 0
      mpime = me - 1
      CALL mp_sum( ngwt )

      ALLOCATE( ctot( ngwt ) )
      DO i = 1, nx
        CALL mergewf(c(:,i), ctot(:), ngw, ig_l2g, mpime, nproc, root)
        IF( mpime == root ) THEN
          WRITE(unit) ( ctot(ig), ig = 1, ngwt )
        END IF
      END DO
      return
      end subroutine write_wfc2

!
!----------------------------------------------------------------------
      subroutine read_wfc2(unit,c)
!----------------------------------------------------------------------
!
! collect wavefunctions on first node and write to file
!
      use gvec
      use gvecs
      use elct
      use para_mod
      use parms
      use work1
      use mp_wave, ONLY: splitwf
      use mp, ONLY: mp_sum
      use parallel_include
!
      implicit none
      integer unit
      complex(kind=8) c(ngw,nx)
      complex(kind=8), allocatable :: ctot(:)
!
      integer i, ig, mpime, proc, ierr, ntot, ncol, mc, ngwt, root

      ngwt = ngw
      root = 0
      mpime = me - 1
      CALL mp_sum( ngwt )

      ALLOCATE( ctot( ngwt ) )
      DO i = 1, nx
        IF( mpime == root ) THEN
          READ(unit) ( ctot(ig), ig = 1, ngwt )
        END IF
        CALL splitwf(c(:,i), ctot(:), ngw, ig_l2g, mpime, nproc, root)
      END DO
      return
      end subroutine read_wfc2


!
!----------------------------------------------------------------------
      subroutine read_wfc(unit,c)
!----------------------------------------------------------------------
!
! read wavefunctions from file and distribute to all nodes
!
      use gvec
      use gvecs
      use elct
      use para_mod
      use parms
      use work1
!
      implicit none
      include 'mpif.h'
      integer unit
      complex(kind=8) c(ngw,nx)
      complex(kind=8), pointer:: psis(:)
!
      integer i, ii, ig, proc, ierr, ntot, ncol, mc, nr
      integer nmin0(3), nmax0(3), nmin(3), nmax(3), n1,n2,nzx,nz,nz_
      integer root, displs(nproc), sendcount(nproc)
      complex(kind=8), allocatable:: psitot(:), psird(:,:,:)
!
! nmin, nmax are the bounds on (i,j,k) indexes of wavefunction G-vectors
!
      call nrbounds(ngw,nr1s,nr2s,nr3s,in1p,in2p,in3p,nmin,nmax)
!
! nzx is the maximum length of a column along z
!
      nzx=nmax(3)-nmin(3)+1
!
      root = 0
! root is the first node
      ntot = 0
      do proc=1,nproc
         sendcount(proc) =  dffts%nsw(proc)*nzx
!
! sendcount(proc) = size of data send to processor proc
!                   (number of columns times length of each column)
!
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + sendcount(proc-1)
         end if
!
! displs(proc) is the position of data sent to processor proc
!
         ntot = ntot + sendcount(proc)
!
! ntot = total size of gathered data
!
      end do
!
! allocate the needed work spaces
!
      psis=> wrk1
      call zero(2*nnrs,psis)
      if (me.eq.1) then
         allocate(psitot(ntot))
!
! read the smallest FFT box that contains all wavefunction G-vectors
!
         read(unit,end=10,err=10) nr, nmin0, nmax0
! check
         if (nmin(1).ne.nmin0(1) .or. nmin(2).ne.nmin0(2) .or.          &
     &       nmin(3).ne.nmin0(3) .or. nmax(1).ne.nmax0(1) .or.          &
     &       nmax(2).ne.nmax0(2) .or. nmax(3).ne.nmax0(3) ) then
            write(6,*) 'read  nmin, nmax =',nmin, nmax
            write(6,*) 'found nmin, nmax =',nmin0, nmax0
            call errore('read_wfc','wavefunction mismatch',1)
         end if
         if (nr.lt.n) call errore('read_wfc','not enough wavefcts',nr)
         allocate(psird(nmin(3):nmax(3),nmin(1):nmax(1),nmin(2):nmax(2)))
      end if
!
      do i=1,n
!
! read the node-number-independent array 
!
         if (me.eq.1) then
            read(unit,end=20,err=20) psird
!
! reorder as contiguously stored columns along z. See comments in 
! write_wfc for the logic (or lack thereof) of the storage
!
            ncol=0
            do proc=1,nproc
               do ii=1,dffts%nsw(proc)
                  ncol=ncol+1
                  ! mc=icpls(ii+ncp0s(proc))
                  mc = dffts%ismap( ii + dffts%iss(proc) )
                  n2=(mc-1)/nr1sx
                  n1= mc-1-n2*nr1sx
                  if (n2.ge.nr2s/2) n2=n2-nr2s
                  if (n1.le.nmax(1).and.(n1.ne.0.or.n2.ge.0)) then
                     do nz=nmin(3),nmax(3)
                        psitot(nz-nmin(3)+1+(ncol-1)*nzx) = psird(nz,n1,n2)
                     end do
                  end if
               end do
            end do
         end if
!
! distribute the array psitot on all nodes (in psis)
!
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('write_wfc','mpi_barrier 2',ierr)
         call mpi_scatterv(psitot,sendcount,displs,MPI_DOUBLE_COMPLEX,  &
    &                      psis  ,sendcount(me),   MPI_DOUBLE_COMPLEX,  &
    &                 root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('write_wfc','mpi_scatter',ierr)
!
! fill c_i(G) (see write_wfc for the logic-or lack thereof-of ordering)
!
         do ig=1,ngw
            ncol=(nps(ig)-1)/nr3sx
            nz_ =nps(ig)-ncol*nr3sx
            nz  =nz_-1
            if (nz.ge.nr3s/2) nz=nz-nr3s
            c(ig,i)=psis(nz-nmin(3)+1+ncol*nzx)
         end do
      end do
!
      if (me.eq.1) then
         deallocate(psird)
         deallocate(psitot)
      end if
!
      return
 10   call errore('read_wfc','file missing or wrong',ierr)
 20   call errore('read_wfc','wavefunction missing or wrong',ierr)
!
      end subroutine read_wfc
!
!-----------------------------------------------------------------------
      subroutine readpfile                                              &
     &     ( flag,ndr,h,hold,nfi,c0,cm,tau0,taum,vel,velm,acc,          &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh)
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
! iflag=-1 : read only h and hold (for use in variable-cell)
! iflag= 0 : read h, hold, c0
! iflag=+1 : read everything
!
      use para_mod
      use mp
      use elct, only: n, nx, ngw, ng0
      use ions_module, only: nsp, na, natx
      use parameters, only: nacx
!
      implicit none
      integer flag, ndr, nfi
      real(kind=8) h(3,3), hold(3,3)
      complex(kind=8) c0(ngw,n), cm(ngw,n)
      real(kind=8) taum(3,natx,nsp),tau0(3,natx,nsp)
      real(kind=8) vel(3,natx,nsp), velm(3,natx,nsp)
      real(kind=8) acc(nacx),lambda(nx,nx), lambdam(nx,nx)
      real(kind=8) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
      real(kind=8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
!
      integer i, ia, is, j, root, ierr
!
! Only the first node reads
!
      if (flag.eq.-1) then
         write(6,'((a,i3,a))') ' ### reading from file ',ndr,' only h  ##'
      else if (flag.eq.0) then
         write(6,'((a,i3,a))') ' ### reading from file ',ndr,' only c0  ##'
      else
         write(6,'((a,i3))') ' ### reading from file ',ndr
      end if
!
      if (me.eq.1) then
         open (unit=ndr,status='old',form='unformatted')
         read(ndr,end=10,err=10) h, hold
      end if
!
! h and hold are needed here by variable-cell calculation
!
      root = 0
      call mpi_barrier( MPI_COMM_WORLD, ierr)
      call mp_bcast(h,   root)
      call mp_bcast(hold,root)
      if (flag.eq.-1) then
         if (me.eq.1) close (unit=ndr)
         return
      end if
!
! now read and distribute the wave functions
!
      ! call read_wfc(ndr,c0)
      call read_wfc2(ndr,c0)
!
      if (flag.eq.0) then
         if (me.eq.1) close (unit=ndr)
         return
      end if
!
      ! call read_wfc(ndr,cm)
      call read_wfc2(ndr,cm)
!
! arrays whose dimensions exceed what is actually read are explicitely set to zero
! in order to prevent potential problems with uninitialized trailing elements
!
      tau0=0.d0
      taum=0.d0
      vel =0.d0
      velm=0.d0
      acc =0.d0
      lambda =0.d0
      lambdam=0.d0
!
! read all other variables
!
      if (me.eq.1) then
         read(ndr,end=10,err=10) nfi
         read(ndr,end=10,err=10) (((tau0(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         read(ndr,end=10,err=10) (((taum(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         read(ndr,end=10,err=10) (((vel (i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         read(ndr,end=10,err=10) (((velm(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         read(ndr,end=10,err=10) acc
         read(ndr,end=10,err=10) ((lambda(i,j),i=1,n),j=1,n)
         read(ndr,end=10,err=10) ((lambdam(i,j),i=1,n),j=1,n)
         read(ndr,end=10,err=10) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
         read(ndr,end=10,err=10) xnhh0,xnhhm,vnhh,velh
         close (unit=ndr)
      end if
!
! broadcast variables to the other nodes
!
      call mpi_barrier( MPI_COMM_WORLD, ierr)
      call mp_bcast( nfi, root)
      call mp_bcast(tau0, root)
      call mp_bcast(taum, root)
      call mp_bcast(vel , root)
      call mp_bcast(velm, root)
      call mp_bcast(acc , root)
      call mp_bcast(lambda, root)
      call mp_bcast(lambdam,root)
      call mp_bcast(xnhe0 , root)
      call mp_bcast(xnhem , root)
      call mp_bcast(vnhe  , root)
      call mp_bcast(xnhp0 , root)
      call mp_bcast(xnhpm , root)
      call mp_bcast(vnhp  , root)
      call mp_bcast(ekincm, root)
      call mp_bcast(xnhh0, root)
      call mp_bcast(xnhhm, root)
      call mp_bcast(vnhh , root)
      call mp_bcast(velh , root)
!
      return
 10   call errore('readpfile','end of file detected',1)
      end
!-----------------------------------------------------------------------
      subroutine writepfile                                             &
     &     ( ndw,h,hold,nfi,c0,cm,tau0,taum,vel,velm,acc,               &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh)
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      use para_mod
      use elct, only: n, nx, ngw, ng0
      use ions_module, only: nsp, na, natx
      use parameters, only: nacx
!
      implicit none
      integer ndw, nfi
      real(kind=8) h(3,3), hold(3,3)
      complex(kind=8) c0(ngw,n), cm(ngw,n)
      real(kind=8) taum(3,natx,nsp),tau0(3,natx,nsp)
      real(kind=8) vel(3,natx,nsp), velm(3,natx,nsp)
      real(kind=8) acc(nacx),lambda(nx,nx), lambdam(nx,nx)
      real(kind=8) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
      real(kind=8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
!
      include 'mpif.h'
      integer i, ia, is, j
!
! Only the first node writes
!
      if (me.eq.1) then
         open (unit=ndw,status='unknown',form='unformatted')
!
! h and hold are needed here by variable-cell calculation
!
         write(ndw) h, hold
      end if
!
! now write the wave functions
!
      ! call write_wfc(ndw,c0)
      ! call write_wfc(ndw,cm)
      call write_wfc2(ndw,c0)
      call write_wfc2(ndw,cm)
!
! write all other variables
!
      if (me.eq.1) then
         write(ndw) nfi
         write(ndw) (((tau0(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         write(ndw) (((taum(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         write(ndw) (((vel (i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         write(ndw) (((velm(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         write(ndw) acc
         write(ndw) ((lambda(i,j),i=1,n),j=1,n)
         write(ndw) ((lambdam(i,j),i=1,n),j=1,n)
         write(ndw) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
         write(ndw) xnhh0,xnhhm,vnhh,velh
         close (unit=ndw)
      end if
      return
      end

!
!----------------------------------------------------------------------
      subroutine nrbounds(ngw,nr1s,nr2s,nr3s,in1p,in2p,in3p,nmin,nmax)
!----------------------------------------------------------------------
!
! find the bounds for (i,j,k) indexes of all wavefunction G-vectors
! The (i,j,k) indexes are defined as: G=i*g(1)+j*g(2)+k*g(3)
! where g(1), g(2), g(3) are basis vectors of the reciprocal lattice
!
      implicit none
! input
      integer ngw,nr1s,nr2s,nr3s,in1p(ngw),in2p(ngw),in3p(ngw)
! output
      integer nmin(3), nmax(3)
! local
      include 'mpif.h'
      integer nmin0(3), nmax0(3), ig, ierr
!
!
      nmin0(1)=  nr1s
      nmax0(1)= -nr1s
      nmin0(2)=  nr2s
      nmax0(2)= -nr2s
      nmin0(3)=  nr3s
      nmax0(3)= -nr3s
!
      do ig=1,ngw
         nmin0(1) = min(nmin0(1),in1p(ig))
         nmin0(2) = min(nmin0(2),in2p(ig))
         nmin0(3) = min(nmin0(3),in3p(ig))
         nmax0(1) = max(nmax0(1),in1p(ig))
         nmax0(2) = max(nmax0(2),in2p(ig))
         nmax0(3) = max(nmax0(3),in3p(ig))
      end do
!
! find minima and maxima for the FFT box across all nodes
!
      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if (ierr.ne.0) call errore('nrbounds','mpi_barrier 1',ierr)
      call mpi_allreduce (nmin0, nmin, 3, MPI_INTEGER, MPI_MIN,         &
     &           MPI_COMM_WORLD, ierr)
      if (ierr.ne.0) call errore('nrbounds','mpi_allreduce min',ierr)
      call mpi_allreduce (nmax0, nmax, 3, MPI_INTEGER, MPI_MAX,         &
     &           MPI_COMM_WORLD, ierr)
      if (ierr.ne.0) call errore('nrbounds','mpi_allreduce max',ierr)

      return
      end subroutine nrbounds

#endif
