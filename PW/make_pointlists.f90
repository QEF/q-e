!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!--------------------------------------------------------------------------
      subroutine make_pointlists
!--------------------------------------------------------------------------
!
! This initialization is needed in order to integrate charge (or
! magnetic moment) in a sphere around the atomic positions.
! This can be used to simply monitor these quantities during the scf
! cycles or in order to calculate constrains on these quantities.
!
! In the input the integration radius r_m can be given, otherwise it is
! calculated here. The integration is a sum over all points in real
! space with the weight 1, if they are closer than r_m to an atom
! and    1 - (distance-r_m)/(0.2*r_m) if r_m<distance<1.2*r_m            
!

      USE kinds,      ONLY : dp
      USE io_global,  ONLY : stdout
      USE ions_base,  ONLY : nat, tau
      USE cell_base,  ONLY : at
      USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrxx
      USE noncollin_module
      USE para

      implicit none
      integer index0,index,indproc,iat,ir,iat1
      integer i,j,k,i0,j0,k0,ipol,ishift(3)

      real(kind=dp) :: posi(3),distance,shift(3),scalprod, distmin

      if (.not.(noncolin)) return
      WRITE( stdout,*) "  Generating pointlists ..."

! First, the real-space position of every point ir is needed ...

! In the parallel case, find the index-offset to account for the planes
! treated by other procs

      index0 = 0
#ifdef __PARA
      do indproc=1,me-1
         index0 = index0 + nrx1*nrx2*npp(indproc)
      enddo
      
#endif

! Check the minimum distance between two atoms in the system

      distmin = 1.d0

      do iat = 1,nat
         do iat1 = iat,nat
               
! posi is the position of a second atom
            do i = -1,1
               do j = -1,1
                  do k = -1,1

                     distance = 0.d0
                     do ipol = 1,3
                        posi(ipol) = tau(ipol,iat1) + real(i)*at(ipol,1) &
     &                       +real(j)*at(ipol,2) + real(k)*at(ipol,3)
                        distance = distance + (posi(ipol)-tau(ipol,iat)) &
     &                       **2.
                     enddo

                     distance = sqrt(distance)
                     if ((distance.lt.distmin).and.(distance.gt.1.d-8)) &
     &                    distmin = distance

                  enddo ! k
               enddo ! j
            enddo ! i

         enddo                  ! iat1
      enddo                     ! iat
                  
      if ((distmin.lt.(2.d0*r_m*1.2d0)).or.(r_m.lt.1.d-8)) then
! Set the radius r_m to a value a little smaller than the minimum
! distance divided by 2*1.2 (so no point in space can belong to more
! than one atom)
         r_m = 0.5d0*distmin/1.2d0 * 0.99d0
         WRITE( stdout,*) "    new r_m : ",r_m
      endif


! Now, make for every atom a list of points which are in their
! integration sphere, as well as a list of weights.
! This also works in the parallel case.
      
      do iat = 1,nat

         pointnum(iat) = 0
         
         do ir=1,nrxx
            index = index0 + ir - 1

            k0 = index/(nrx1*nrx2)
            index = index - (nrx1*nrx2) * k0
            j0 = index / nrx1
            index = index - nrx1*j0
            i0 = index
            
            do i = i0-nr1,i0+nr1, nr1
               do j = j0-nr2, j0+nr2, nr2
                  do k = k0-nr3, k0+nr3, nr3

                     do ipol=1,3
                        posi(ipol) = real(i)/real(nr1) * at(ipol &
                             ,1) +real(j)/real(nr2) * at(ipol,2) &
                             +real(k)/real(nr3) * at(ipol,3)
                  
                        posi(ipol) = posi(ipol) - tau(ipol,iat)
                     enddo

                     distance = sqrt(posi(1)**2.+posi(2)**2.+posi(3 &
                          )**2.)
                     
               
                     if (distance.le.r_m) then
                        pointnum(iat) = pointnum(iat) + 1
                        factlist(pointnum(iat),iat) = 1.d0
                        pointlist(pointnum(iat),iat) = ir
                  
                     else if (distance.le.1.2*r_m) then
                        pointnum(iat) = pointnum(iat) + 1
                        factlist(pointnum(iat),iat) = 1.d0 - (distance &
                             -r_m)/(0.2*r_m)
                        pointlist(pointnum(iat),iat) = ir

                     endif


                  enddo         ! k
               enddo            ! j
            enddo               ! i

         enddo                  ! ir


      enddo                     ! ipol
           
      end

