!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!----------------------------------------------------------------------
      subroutine dipol_matrix(tau0,h,rho, dipol)
!----------------------------------------------------------------------
! This subroutine calculates the dipole element between two states, 
! the states are given from rho(r)=Psi_1(r)*Psi_2(r) 
! The density is assumed to be a real function

      use parameters, only: nsx, natx
      use ions_base, only: nsp, na, pmass
      use para_mod
      use grid_dimensions, only: nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr => nnrx
      use io_global, only: ionode
      use mp, only: mp_bcast
   

      implicit none

#ifdef __PARA
      include 'mpif.h'
#endif

      integer specie(100) !atomic species to modify
      integer i, j, k, ir, ia, is, ityp(natx), nat00, tp(3), isa
      integer ir1, ir2, ir3, ip1, ip2, ip3, ipn
      real(8) tau0(3,natx), tau00(3,natx), rho(nnr)
      real(8) h(3,3), l1, l2, l3, shift(3), maxr, minr, cm(3)
      real(8) ll1, ll2, ll3, tot_m
      real(8) rho_aux
      real(8) dipol(3) ! dipole vector: Int_R rho(r)r

#ifdef __PARA
      integer ip, ierr, incr(nproc), displs(nproc)
      real(8), allocatable:: rhow(:)
#endif


      dipol(:)=0.d0


      l1 = h(1,1) + h(2,1) + h(3,1)
      l2 = h(1,2) + h(2,2) + h(3,2)
      l3 = h(1,3) + h(2,3) + h(3,3)
      ll1 = dsqrt(h(1,1)**2+h(1,2)**2+h(1,3)**2)
      ll2 = dsqrt(h(2,1)**2+h(2,2)**2+h(2,3)**2)
      ll3 = dsqrt(h(3,1)**2+h(3,2)**2+h(3,3)**2)
      nat00 = 0
      tot_m = 0.d0
      do i = 1,3
         cm(i) = 0.d0
      end do
      isa = 0
      do is = 1,nsp
         do ia = 1,na(is)
            tot_m = tot_m + pmass(is)
            nat00 = nat00 + 1
            isa = isa + 1
            do i = 1,3
               tau00(i,nat00) = tau0(i,isa)
               cm(i) = cm(i) + tau00(i,nat00)*pmass(is)
            end do
            ityp(nat00) = is
         end do
      end do
      do i = 1,3
         cm(i) = cm(i)/tot_m
      end do

! to center the plot of the charge density at the center of the unit cell where also
! the center of mass of the system is moved

      shift(1) = 0.5d0*l1 - cm(1)
      tp(1) = nint(shift(1)*DBLE(nr1)/ll1)
      shift(1) = 0.d0 !DBLE(tp(1))*ll1/DBLE(nr1)
      shift(2) = 0.5d0*l2 - cm(2)
      tp(2) = nint(shift(2)*DBLE(nr2)/ll2)
      shift(2) = 0.d0 !DBLE(tp(2))*ll2/DBLE(nr2)
      shift(3) = 0.5d0*l3 - cm(3)
      tp(3) = nint(shift(3)*DBLE(nr3)/ll3)
      shift(3) = 0.d0 !DBLE(tp(3))*ll3/DBLE(nr3)

#ifdef __PARA
    

      if (me.eq.1) allocate(rhow(nr1x*nr2x*nr3x))
      do ip=1,nproc
         incr(ip) = dfftp%nnp * ( dfftp%npp(ip) )
         if (ip.eq.1) then
            displs(ip)=0
         else
            displs(ip)=displs(ip-1) + incr(ip)
         end if
      end do
      call mpi_barrier ( MPI_COMM_WORLD, ierr)
      call mpi_gatherv (rho, incr(me), MPI_REAL8,                       &
     &                  rhow,incr, displs, MPI_REAL8,                   &
     &                     0, MPI_COMM_WORLD, ierr)
      if (ierr.ne.0) call errore('mpi_gatherv','ierr<>0',ierr)

! in parallel execution, only the first nodes writes

      if (me.eq.1) then
#endif

      maxr = 0.d0
      minr = 10.d0
      do ir3 = 1,nr3
         if ((ir3-tp(3)).le.0) then
            ip3 = (ir3-tp(3))+nr3
         else
            ip3 = (ir3-tp(3))
         end if
         do ir2 = 1,nr2
            if ((ir2-tp(2)).le.0) then
               ip2 = (ir2-tp(2))+nr2
            else
               ip2 = (ir2-tp(2))
            end if
            do ir1 = 1,nr1
               if ((ir1-tp(1)).le.0) then
                  ip1 = (ir1-tp(1))+nr1
               else
                  ip1 = (ir1-tp(1))
               end if
               ir = ir1 + (ir2-1)*nr1 + (ir3-1)*nr2*nr1
               ipn = ip1 + (ip2-1)*nr1 + (ip3-1)*nr2*nr1
#ifdef __PARA
               rho_aux = rhow(ir) !rhow(ipn)
#else
               rho_aux = rho(ir) !rho(ipn)
#endif
               dipol(1)=dipol(1)+rho_aux*DBLE(ir1)*ll1/DBLE(nr1)
               dipol(2)=dipol(2)+rho_aux*DBLE(ir2)*ll2/DBLE(nr2)
               dipol(3)=dipol(3)+rho_aux*DBLE(ir3)*ll3/DBLE(nr3)
            end do
         end do
      end do
#ifdef __PARA
      deallocate (rhow)
#endif
  

      close(40)
#ifdef __PARA
      end if
#endif

      return
    end subroutine dipol_matrix
