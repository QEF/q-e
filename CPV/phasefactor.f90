!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"


!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE phase_factors_module
!=----------------------------------------------------------------------------=!

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE phfacs(eigr,atoms)
!  SUBROUTINE strucf(atoms,eigr,gv)
!  ----------------------------------------------

!  END manual

! ...   declare modules
        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: strucf

! ...   --+ end of module-scope declarations +--

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
!  BEGIN manual

        SUBROUTINE phfacs( ei1, ei2, ei3, eigr, mill, taus, nat )

!  this routine computes the phase factors
!
!    ei1(ix,ia) = exp (i ix G_1 dot R(ia))
!    ei2(iy,ia) = exp (i iy G_2 dot R(ia))
!    ei3(iz,ia) = exp (i iz G_3 dot R(ia))
!
!    eigr(ig,ia) = exp(i G dot R(ia)) =
!                = ei1(ix,ia) * ei2(iy,ia) * ei3(iz,ia)
!
!    G_1,G_2,G_3 = reciprocal lattice generators
!
!    ia = index of ion
!    ig = index of G vector
!    ix,iy,iz = Miller indices
!  ----------------------------------------------
 
!  END manual

          USE constants,       ONLY: tpi
          USE grid_dimensions, ONLY: nr1, nr2, nr3

! ...     declare subroutine arguments
          INTEGER, INTENT(IN) :: nat
          COMPLEX(dbl) :: ei1( -nr1 : nr1, nat )
          COMPLEX(dbl) :: ei2( -nr2 : nr2, nat )
          COMPLEX(dbl) :: ei3( -nr3 : nr3, nat )
          COMPLEX(dbl) :: eigr( :, : )
          REAL(dbl) :: taus( 3, nat )
          INTEGER      :: mill( :, : )

! ...     declare other variables
          COMPLEX(dbl) :: ctep1, ctep2, ctep3, ctem1, ctem2, ctem3
          REAL(dbl) :: ar1, ar2, ar3
          INTEGER :: i, j, k, isa
          INTEGER :: ig, ig1, ig2, ig3, ngw

! ...     --+ end of declarations +--

          DO isa = 1, nat

! ...         Miller index = 0: exp(i 0 dot R(ia)) = 1

              ei1( 0, isa ) = CMPLX( 1.d0, 0.d0 )
              ei2( 0, isa ) = CMPLX( 1.d0, 0.d0 )
              ei3( 0, isa ) = CMPLX( 1.d0, 0.d0 )

! ...         let R_1,R_2,R_3 be the direct lattice generators,
! ...         G_1,G_2,G_3 the reciprocal lattice generators
! ...         by definition G_i dot R_j = 2 pi delta_{ij}
! ...         ionic coordinates are in units of R_1,R_2,R_3
! ...         then G_i dot R(ia) = 2 pi R(ia)_i

              ar1 = tpi * taus(1,isa)  ! G_1 dot R(ia)
              ar2 = tpi * taus(2,isa)  ! G_2 dot R(ia)
              ar3 = tpi * taus(3,isa)  ! G_3 dot R(ia)

! ...         Miller index = 1: exp(i G_i dot R(ia))

              ctep1 = CMPLX( cos( ar1 ), -sin( ar1 ) )
              ctep2 = CMPLX( cos( ar2 ), -sin( ar2 ) )
              ctep3 = CMPLX( cos( ar3 ), -sin( ar3 ) )

! ...         Miller index = -1: exp(-i G_i dot R(ia))

              ctem1 = CONJG(ctep1)
              ctem2 = CONJG(ctep2)
              ctem3 = CONJG(ctep3)

! ...         Miller index > 0: exp(i N G_i dot R(ia)) =
! ...           = exp(i G_i dot R(ia)) exp(i (N-1) G_i dot R(ia))
! ...         Miller index < 0: exp(-i N G_i dot R(ia)) =
! ...           = exp(-i G_i dot R(ia)) exp(-i (N-1) G_i dot R(ia))

              DO i = 1, nr1 - 1
                ei1(  i, isa ) = ei1(  i - 1, isa ) * ctep1
                ei1( -i, isa ) = ei1( -i + 1, isa ) * ctem1
              END DO
              DO j = 1, nr2 - 1
                ei2(  j, isa ) = ei2(  j - 1, isa ) * ctep2
                ei2( -j, isa ) = ei2( -j + 1, isa ) * ctem2
              END DO
              DO k = 1, nr3 - 1
                ei3(  k, isa ) = ei3(  k - 1, isa ) * ctep3
                ei3( -k, isa ) = ei3( -k + 1, isa ) * ctem3
              END DO

          END DO

          ngw = SIZE( eigr, 1 )
          IF( ngw > SIZE( mill, 2 ) ) THEN
            CALL errore(' phfacs ',' eigr inconsisten size ',ngw)
          END IF

          DO ig = 1, ngw
            ig1 = mill( 1, ig )
            ig2 = mill( 2, ig )
            ig3 = mill( 3, ig )
            DO i = 1, nat
              eigr( ig, i ) = ei1( ig1, i ) * ei2( ig2, i ) * ei3( ig3, i )
            END DO
          END DO

          RETURN
        END SUBROUTINE phfacs


!=----------------------------------------------------------------------------=!
!  BEGIN manual

      SUBROUTINE strucf( sfac, atoms, eigr, ei1, ei2, ei3, gv )

!  this routine computes the structure factors
!
!    sfac(ig,is) = (sum over ia) exp(i G dot R(ia)) =
!                  (sum over ia) ei1(ix,ia) * ei2(iy,ia) * ei3(iz,ia)
!
!    ei1(ix,ia) = exp (i ix G_1 dot R(ia))
!    ei2(iy,ia) = exp (i iy G_2 dot R(ia))
!    ei3(iz,ia) = exp (i iz G_3 dot R(ia))
!
!    G_1,G_2,G_3 = reciprocal lattice generators
!
!    ia = index of ion (running over ions of species is)
!    ig = index of G vector
!    is = index of atomic species
!    ix,iy,iz = Miller indices of G vector
!  ----------------------------------------------
!  END manual

      USE atoms_type_module, ONLY: atoms_type
      USE cp_types, ONLY: recvecs
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE ions_base, ONLY: nat

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (atoms_type) :: atoms
      COMPLEX(dbl) :: eigr(:,:)
      COMPLEX(dbl) :: ei1( -nr1 : nr1, nat )
      COMPLEX(dbl) :: ei2( -nr2 : nr2, nat )
      COMPLEX(dbl) :: ei3( -nr3 : nr3, nat )
      TYPE (recvecs) :: gv
      COMPLEX(dbl), INTENT(OUT) :: sfac(:,:)

! ... declare other variables
      INTEGER :: is, ig, ia, ig1, ig2, ig3, isa

! ...     --+ end of declarations +--
      
      CALL phfacs( ei1, ei2, ei3, eigr, gv%mill, atoms%taus, atoms%nat )

      DO ig = 1, gv%ng_l
        ig1 = gv%mill( 1, ig ) 
        ig2 = gv%mill( 2, ig ) 
        ig3 = gv%mill( 3, ig )
        isa = 1
        DO is = 1, atoms%nsp
          sfac( is, ig ) = CMPLX (0.0_dbl, 0.0_dbl)
          DO ia = 1, atoms%na(is)
            sfac( is, ig ) = sfac( is, ig ) + &
              ei1( ig1, isa ) * ei2( ig2, isa ) * ei3( ig3, isa )
            isa = isa + 1
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE strucf

!=----------------------------------------------------------------------------=!
      END MODULE phase_factors_module
!=----------------------------------------------------------------------------=!


!-----------------------------------------------------------------------
      subroutine phfac( tau0, ei1, ei2, ei3, eigr)
!-----------------------------------------------------------------------
!  this subroutine generates the complex matrices ei1, ei2, and ei3
!  used to compute the structure factor and forces on atoms :
!     ei1(n1,ia,is) = exp(-i*n1*b1*tau(ia,is)) -nr1<n1<nr1
!  and similar definitions for ei2 and ei3 ; and :
!     eigr(n,ia,is) = ei1*ei2*ei3 = exp(-i g*tau(ia,is))
!  The value of n1,n2,n3 for a vector g is supplied by arrays mill_l
!  calculated in ggen .
!
      use constants,          only: pi, fpi
      use parameters,         only: natx, nsx
      use control_flags,      only: iprint, iprsta
      use io_global,          only: stdout
      use ions_base,          only: nsp, na, nat
      use cell_base,          only: ainv
      use grid_dimensions,    only: nr1, nr2, nr3
      use gvecw,              only: ngw
      use reciprocal_vectors, only: mill_l
!
      implicit none
      real(kind=8) tau0(3,natx)
!
      complex(kind=8) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),      &
     &                ei3(-nr3:nr3,nat), eigr(ngw,nat)
!
      integer i,j,k, ia, is, ig, isa
      real(kind=8) taup(3), s, ar1,ar2,ar3
      complex(kind=8) ctep1,ctep2,ctep3,ctem1,ctem2,ctem3
!
      if(nr1.lt.3) call errore(' phfac ',' nr1 too small ',nr1)
      if(nr2.lt.3) call errore(' phfac ',' nr1 too small ',nr2)
      if(nr3.lt.3) call errore(' phfac ',' nr1 too small ',nr3)
!
      if(iprsta.gt.3) then
         WRITE( stdout,*) ' phfac: tau0 '
         WRITE( stdout,*) ( ( tau0(i,isa), i=1, 3 ), isa=1, SUM(na(1:nsp)) )
      endif
      isa = 0
      do is=1,nsp
         do ia=1,na(is)
            isa = isa + 1
            do i=1,3
               s=0.d0
               do j=1,3
                  s=s+ainv(i,j)*tau0(j,isa)
               end do
               taup(i)=s
!
! tau0=x1*a1+x2*a2+x3*a3 => taup(1)=x1=tau0*b1 and so on
!
            end do
!
            ar1=2.d0*pi*taup(1)
            ctep1=cmplx(cos(ar1),-sin(ar1))
            ctem1=conjg(ctep1)
            ei1( 0,isa)=cmplx(1.d0,0.d0)
            ei1( 1,isa)=ctep1
            ei1(-1,isa)=ctem1
            do i=2,nr1-1
               ei1( i,isa)=ei1( i-1,isa)*ctep1
               ei1(-i,isa)=ei1(-i+1,isa)*ctem1
            end do
!
            ar2=2.d0*pi*taup(2)
            ctep2=cmplx(cos(ar2),-sin(ar2))
            ctem2=conjg(ctep2)
            ei2( 0,isa)=cmplx(1.d0,0.d0)
            ei2( 1,isa)=ctep2
            ei2(-1,isa)=ctem2
            do j=2,nr2-1
               ei2( j,isa)=ei2( j-1,isa)*ctep2
               ei2(-j,isa)=ei2(-j+1,isa)*ctem2
            end do
!
            ar3=2.d0*pi*taup(3)
            ctep3=cmplx(cos(ar3),-sin(ar3))
            ctem3=conjg(ctep3)
            ei3( 0,isa)=cmplx(1.d0,0.d0)
            ei3( 1,isa)=ctep3
            ei3(-1,isa)=ctem3
            do k=2,nr3-1
               ei3( k,isa)=ei3( k-1,isa)*ctep3
               ei3(-k,isa)=ei3(-k+1,isa)*ctem3
            end do
!
         end do
      end do
!
      if(iprsta.gt.4) then
         WRITE( stdout,*)
         if(nsp.gt.1) then
            isa = 0
            do is=1,nsp
               WRITE( stdout,'(33x,a,i4)') ' ei1, ei2, ei3 (is)',is
               do ig=1,4
                  WRITE( stdout,'(6f9.4)')                                    &
     &                 ei1(ig,1+isa),ei2(ig,1+isa),ei3(ig,1+isa)
               end do
               WRITE( stdout,*)
               isa = isa + na(is)
            end do
         else
            do ia=1,na(1)
               WRITE( stdout,'(33x,a,i4)') ' ei1, ei2, ei3 (ia)',ia
               do ig=1,4
                  WRITE( stdout,'(6f9.4)')                                    &
     &                 ei1(ig,ia),ei2(ig,ia),ei3(ig,ia)
               end do
               WRITE( stdout,*)
            end do
         endif
      endif
!
!     calculation of eigr(g,ia,is)=e^(-ig.r(ia,is))
!
      do isa=1,nat
            do ig=1,ngw
               i = mill_l(1,ig)
               j = mill_l(2,ig)
               k = mill_l(3,ig)
               eigr(ig,isa) = ei1(i,isa) * ei2(j,isa) * ei3(k,isa)
            end do
      end do
!
      return
      end


!-------------------------------------------------------------------------
      subroutine strucf (ei1,ei2,ei3,sfac)
!-----------------------------------------------------------------------
! computes the structure factor sfac(ngs,nsp) and returns it in "pseu"
!
!
      use gvecs,              only: ngs
      use constants,          only: pi, fpi
      use grid_dimensions,    only: nr1, nr2, nr3
      use ions_base,          only: nsp, na, nat
      use reciprocal_vectors, only: mill_l
!
      implicit none
      complex(kind=8) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),          &
     &           ei3(-nr3:nr3,nat), sfac(ngs,nsp)
      integer is, ig, ia, i, j, k, isa
!
      call start_clock( 'strucf' )
      isa = 1
      do is=1,nsp
         do ig=1,ngs
            sfac(ig,is)=(0.,0.)
         end do
         do ia=1,na(is)
            do ig=1,ngs 
               i = mill_l( 1, ig )
               j = mill_l( 2, ig )
               k = mill_l( 3, ig )
               sfac(ig,is)=sfac(ig,is) + ei1(i,isa) *ei2(j,isa) *ei3(k,isa)
            end do
            isa = isa + 1
         end do
      end do
!     
      call stop_clock( 'strucf' )
      return
      end subroutine


!-----------------------------------------------------------------------
      subroutine phbox(taub,eigrb)
!-----------------------------------------------------------------------
!     calculates the phase factors for the g's of the little box
!     eigrt=exp(-i*g*tau) .
!     Uses the same logic for fast calculation as in phfac (see below)
!        
      use io_global,     only: stdout
      use control_flags, only: iprint, iprsta
      use constants,     only: pi, fpi
      use parameters,    only: natx, nsx
      use ions_base,     only: nsp, na, nat
      use gvecb,         only: ngb, mill_b
      use cell_base,     only: ainv
      use small_box,     only: ainvb
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
!                 
      implicit none    
      real(kind=8) :: taub(3,natx)
      complex(kind=8) :: eigrb(ngb,nat)
! local           
      integer :: i,j,k, is, ia, ig, isa
      real(kind=8) taup(3),s,ar1,ar2,ar3
      complex(kind=8), allocatable:: ei1b(:,:), ei2b(:,:), ei3b(:,:)
      complex(kind=8) ctep1,ctep2,ctep3,ctem1,ctem2,ctem3
!
      if(nr1b.lt.3) call errore(' phbox ',' nr1b too small ',nr1b)
      if(nr2b.lt.3) call errore(' phbox ',' nr2b too small ',nr2b)
      if(nr3b.lt.3) call errore(' phbox ',' nr3b too small ',nr3b)
!
      allocate(ei1b(-nr1b:nr1b,nat))
      allocate(ei2b(-nr2b:nr2b,nat))
      allocate(ei3b(-nr3b:nr3b,nat))
!
      if(iprsta.gt.3) then
         WRITE( stdout,*) ' phbox: taub '
         WRITE( stdout,*) ( (taub(i,isa), i=1, 3 ), isa=1, SUM(na(1:nsp)) )
      endif

      do isa=1,nat

            do i=1,3
               s=0.d0
               do j=1,3
                  s = s + ainvb(i,j)*taub(j,isa)
               end do
               taup(i)=s
            end do
!
            ar1=2.d0*pi*taup(1)
            ctep1=cmplx(cos(ar1),-sin(ar1))
            ctem1=conjg(ctep1)
            ei1b( 0,isa)=cmplx(1.d0,0.d0)
            ei1b( 1,isa)=ctep1
            ei1b(-1,isa)=ctem1
            do i=2,nr1b-1
               ei1b( i,isa)=ei1b( i-1,isa)*ctep1
               ei1b(-i,isa)=ei1b(-i+1,isa)*ctem1
            end do
!
            ar2=2.d0*pi*taup(2)
            ctep2=cmplx(cos(ar2),-sin(ar2))
            ctem2=conjg(ctep2)
            ei2b( 0,isa)=cmplx(1.d0,0.d0)
            ei2b( 1,isa)=ctep2
            ei2b(-1,isa)=ctem2
            do j=2,nr2b-1
               ei2b( j,isa)=ei2b( j-1,isa)*ctep2
               ei2b(-j,isa)=ei2b(-j+1,isa)*ctem2
            end do
!
            ar3=2.d0*pi*taup(3)
            ctep3=cmplx(cos(ar3),-sin(ar3))
            ctem3=conjg(ctep3)
            ei3b( 0,isa)=cmplx(1.d0,0.d0)
            ei3b( 1,isa)=ctep3
            ei3b(-1,isa)=ctem3
            do k=2,nr3b-1
               ei3b( k,isa)=ei3b( k-1,isa)*ctep3
               ei3b(-k,isa)=ei3b(-k+1,isa)*ctem3
            end do
!
      end do
!
!     calculation of eigrb(g,ia,is)=e^(-ig.r(ia,is))
!
      do isa=1,nat
         do ig=1,ngb
            i = mill_b(1,ig)
            j = mill_b(2,ig)
            k = mill_b(3,ig)
            eigrb(ig,isa) = ei1b(i,isa) * ei2b(j,isa) * ei3b(k,isa)
         end do
      end do
!
      if(iprsta.gt.4) then
         WRITE( stdout,*)
         if(nsp.gt.1) then
            isa = 0
            do is=1,nsp
               WRITE( stdout,'(33x,a,i4)') ' ei1b, ei2b, ei3b (is)',is
               do ig=1,4
                  WRITE( stdout,'(6f9.4)')                                    &
     &                 ei1b(ig,1+isa),ei2b(ig,1+isa),ei3b(ig,1+isa)
               end do
               WRITE( stdout,*)
               isa = isa + na(is)
            end do
         else
            do ia=1,na(1)
               WRITE( stdout,'(33x,a,i4)') ' ei1b, ei2b, ei3b (ia)',ia
               do ig=1,4
                  WRITE( stdout,'(6f9.4)')                                    &
     &                 ei1b(ig,ia),ei2b(ig,ia),ei3b(ig,ia)
               end do
               WRITE( stdout,*)
            end do
         endif
      endif
!
      deallocate(ei3b)
      deallocate(ei2b)
      deallocate(ei1b)
!
      return
      end subroutine

