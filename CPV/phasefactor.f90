!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


!=----------------------------------------------------------------------------=!
   SUBROUTINE phfacs_x( ei1, ei2, ei3, eigr, mill, taus, nr1, nr2, nr3, nat )
!=----------------------------------------------------------------------------=!

      !  this routine computes the phase factors
      !
      !    ei1(ix,ia) = exp ( -i ix G_1 dot R(ia))
      !    ei2(iy,ia) = exp ( -i iy G_2 dot R(ia))
      !    ei3(iz,ia) = exp ( -i iz G_3 dot R(ia))
      !
      !    eigr(ig,ia) = exp( -i G dot R(ia)) =
      !                = ei1(ix,ia) * ei2(iy,ia) * ei3(iz,ia)
      !
      !    G_1,G_2,G_3 = reciprocal lattice generators
      !
      !    ia = index of ion
      !    ig = index of G vector
      !    ix,iy,iz = Miller indices
      !  ----------------------------------------------
 
          USE kinds,           ONLY: DP
          USE constants,       ONLY: tpi

          IMPLICIT NONE

          ! ...     declare subroutine arguments

          INTEGER, INTENT(IN) :: nat
          INTEGER, INTENT(IN) :: nr1, nr2, nr3
          COMPLEX(DP) :: ei1( -nr1 : nr1, nat )
          COMPLEX(DP) :: ei2( -nr2 : nr2, nat )
          COMPLEX(DP) :: ei3( -nr3 : nr3, nat )
          COMPLEX(DP) :: eigr( :, : )
          REAL(DP) :: taus( 3, nat )
          INTEGER      :: mill( :, : )

          ! ...     declare other variables

          COMPLEX(DP) :: ctep1, ctep2, ctep3, ctem1, ctem2, ctem3
          REAL(DP) :: ar1, ar2, ar3
          INTEGER :: i, j, k, isa
          INTEGER :: ig, ig1, ig2, ig3, ngw

          ! ...     --+ end of declarations +--

          if(nr1 < 3) call errore(' phfacs ',' nr1 too small ',1)
          if(nr2 < 3) call errore(' phfacs ',' nr2 too small ',1)
          if(nr3 < 3) call errore(' phfacs ',' nr3 too small ',1)

          DO isa = 1, nat

              ! ... Miller index = 0: exp(i 0 dot R(ia)) = 1

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

              ! ...         Miller index = 1: exp(-i G_i dot R(ia))

              ctep1 = CMPLX( cos( ar1 ), -sin( ar1 ) )
              ctep2 = CMPLX( cos( ar2 ), -sin( ar2 ) )
              ctep3 = CMPLX( cos( ar3 ), -sin( ar3 ) )

              ! ...         Miller index = -1: exp(-i G_im dot R(ia)) = exp(i G_i dot R(ia))

              ctem1 = CONJG(ctep1)
              ctem2 = CONJG(ctep2)
              ctem3 = CONJG(ctep3)

              ! ...         Miller index > 0: exp(i N G_i dot R(ia)) =
              ! ...           = exp(i G_i dot R(ia)) exp(i (N-1) G_i dot R(ia))
              ! ...         Miller index < 0: exp(-i N G_i dot R(ia)) =
              ! ...           = exp(-i G_i dot R(ia)) exp(-i (N-1) G_i dot R(ia))

              DO i = 1, nr1
                ei1(  i, isa ) = ei1(  i - 1, isa ) * ctep1
                ei1( -i, isa ) = ei1( -i + 1, isa ) * ctem1
              END DO
              DO j = 1, nr2
                ei2(  j, isa ) = ei2(  j - 1, isa ) * ctep2
                ei2( -j, isa ) = ei2( -j + 1, isa ) * ctem2
              END DO
              DO k = 1, nr3
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
      END SUBROUTINE phfacs_x


!=----------------------------------------------------------------------------=!
   SUBROUTINE strucf_x( sfac, ei1, ei2, ei3, mill, ngm )
!=----------------------------------------------------------------------------=!

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


      USE kinds,            ONLY: DP
      USE ions_base,        ONLY: nat, na, nsp
      use grid_dimensions,  only: nr1, nr2, nr3

      IMPLICIT NONE

      ! ... declare subroutine arguments
      !
      COMPLEX(DP) :: ei1( -nr1 : nr1, nat )
      COMPLEX(DP) :: ei2( -nr2 : nr2, nat )
      COMPLEX(DP) :: ei3( -nr3 : nr3, nat )
      INTEGER      :: mill( :, : )
      INTEGER      :: ngm
      COMPLEX(DP), INTENT(OUT) :: sfac(:,:)

      ! ... declare other variables
      !
      INTEGER :: is, ig, ia, ig1, ig2, ig3, isa

      call start_clock( 'strucf' )

      DO ig = 1, ngm
        ig1 = mill( 1, ig ) 
        ig2 = mill( 2, ig ) 
        ig3 = mill( 3, ig )
        isa = 1
        DO is = 1, nsp
          sfac( ig, is ) = CMPLX (0.0d0, 0.0d0)
          DO ia = 1, na(is)
            sfac( ig, is ) = sfac( ig, is ) + &
              ei1( ig1, isa ) * ei2( ig2, isa ) * ei3( ig3, isa )
            isa = isa + 1
          END DO
        END DO
      END DO

      call stop_clock( 'strucf' )

      RETURN
   END SUBROUTINE strucf_x



!-----------------------------------------------------------------------

   subroutine phfac_x( tau0, ei1, ei2, ei3, eigr)

      !-----------------------------------------------------------------------
      !  this subroutine generates the complex matrices ei1, ei2, and ei3
      !  used to compute the structure factor and forces on atoms :
      !     ei1(n1,ia,is) = exp(-i*n1*b1*tau(ia,is)) -nr1<n1<nr1
      !  and similar definitions for ei2 and ei3 ; and :
      !     eigr(n,ia,is) = ei1*ei2*ei3 = exp(-i g*tau(ia,is))
      !  The value of n1,n2,n3 for a vector g is supplied by arrays mill_l
      !  calculated in ggen .
      !
      use kinds,              only: DP
      use control_flags,      only: iprsta
      use io_global,          only: stdout
      use ions_base,          only: nsp, na, nat
      use cell_base,          only: ainv, r_to_s
      use grid_dimensions,    only: nr1, nr2, nr3
      use reciprocal_vectors, only: mill_l
      use gvecw,              only: ngw
      use cp_interfaces,      only: phfacs
!
      implicit none
      real(DP) tau0(3,nat)
!
      complex(DP) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),      &
     &                ei3(-nr3:nr3,nat), eigr(ngw,nat)
!
      integer :: i, isa
      real(DP), allocatable :: taus(:,:)
!
      allocate( taus(3,nat) )
!
      if(iprsta.gt.3) then
         WRITE( stdout,*) ' phfac: tau0 '
         WRITE( stdout,*) ( ( tau0(i,isa), i=1, 3 ), isa=1, nat )
      endif
      CALL r_to_s( tau0, taus, na, nsp, ainv )
      CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )

      deallocate( taus )
!
      return
   end subroutine phfac_x


!-----------------------------------------------------------------------
      SUBROUTINE phbox( taub, eigrb, ainvb )
!-----------------------------------------------------------------------
!     calculates the phase factors for the g's of the little box
!     eigrt=exp(-i*g*tau) .
!     Uses the same logic for fast calculation as in phfac (see below)
!        
      USE kinds,         only: DP
      use io_global,     only: stdout
      use control_flags, only: iprsta
      use ions_base,     only: nsp, na, nat
      use gvecb,         only: ngb, mill_b
      use cell_base,     only: r_to_s
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
      use cp_interfaces, only: phfacs
!                 
      IMPLICIT NONE    
      REAL(DP)    :: taub(3,nat)
      COMPLEX(DP) :: eigrb(ngb,nat)
      REAL(DP)    :: ainvb(3,3)
! local           
      integer :: i,j,k, is, ia, ig, isa
      complex(8), allocatable:: ei1b(:,:), ei2b(:,:), ei3b(:,:)
      real(8), allocatable :: taus(:,:)
!
      allocate(ei1b(-nr1b:nr1b,nat))
      allocate(ei2b(-nr2b:nr2b,nat))
      allocate(ei3b(-nr3b:nr3b,nat))
      allocate( taus( 3, nat ) )
!
      if(iprsta.gt.3) then
         WRITE( stdout,*) ' phbox: taub '
         WRITE( stdout,*) ( (taub(i,isa), i=1, 3 ), isa=1, nat )
      endif
      CALL r_to_s( taub, taus, na, nsp, ainvb )
      CALL phfacs( ei1b, ei2b, ei3b, eigrb, mill_b, taus, nr1b, nr2b, nr3b, nat )
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
      deallocate( taus )
!
      RETURN
      END SUBROUTINE phbox

