!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   subroutine phfac_x( tau0, ei1, ei2, ei3, eigr)
!=----------------------------------------------------------------------------=!
      !
      !  this subroutine generates the complex matrices ei1, ei2, and ei3
      !  used to compute the structure factor and forces on atoms :
      !     ei1(n1,ia,is) = exp(-i*n1*b1*tau(ia,is)) -nr1<n1<nr1
      !  and similar definitions for ei2 and ei3 ; and :
      !     eigr(n,ia,is) = ei1*ei2*ei3 = exp(-i g*tau(ia,is))
      !  The value of n1,n2,n3 for a vector g is supplied by array mill
      !  calculated in ggen .
      !
      use kinds,              only: DP
      use control_flags,      only: iverbosity
      use io_global,          only: stdout
      use ions_base,          only: nsp, na, nat
      use cell_base,          only: ainv, r_to_s
      use fft_base,           only: dfftp
      use gvect, only: mill
      use gvecw,              only: ngw
      use cp_interfaces,      only: phfacs
!
      implicit none
      real(DP) tau0(3,nat)
!
      complex(DP) ei1(-dfftp%nr1:dfftp%nr1,nat), ei2(-dfftp%nr2:dfftp%nr2,nat),      &
     &                ei3(-dfftp%nr3:dfftp%nr3,nat), eigr(ngw,nat)
!
      integer :: i, isa
      real(DP), allocatable :: taus(:,:)
!
      allocate( taus(3,nat) )
!
      if( iverbosity > 2 ) then
         WRITE( stdout,*) ' phfac: tau0 '
         WRITE( stdout,*) ( ( tau0(i,isa), i=1, 3 ), isa=1, nat )
      endif
      CALL r_to_s( tau0, taus, na, nsp, ainv )
      CALL phfacs( ei1, ei2, ei3, eigr, mill, taus, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat )

      deallocate( taus )
!
      return
   end subroutine phfac_x


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

              ei1( 0, isa ) = ( 1.d0, 0.d0 )
              ei2( 0, isa ) = ( 1.d0, 0.d0 )
              ei3( 0, isa ) = ( 1.d0, 0.d0 )

              ! ...         let R_1,R_2,R_3 be the direct lattice generators,
              ! ...         G_1,G_2,G_3 the reciprocal lattice generators
              ! ...         by definition G_i dot R_j = 2 pi delta_{ij}
              ! ...         ionic coordinates are in units of R_1,R_2,R_3
              ! ...         then G_i dot R(ia) = 2 pi R(ia)_i

              ar1 = tpi * taus(1,isa)  ! G_1 dot R(ia)
              ar2 = tpi * taus(2,isa)  ! G_2 dot R(ia)
              ar3 = tpi * taus(3,isa)  ! G_3 dot R(ia)

              ! ...         Miller index = 1: exp(-i G_i dot R(ia))

              ctep1 = CMPLX( cos( ar1 ), -sin( ar1 ) ,kind=DP)
              ctep2 = CMPLX( cos( ar2 ), -sin( ar2 ) ,kind=DP)
              ctep3 = CMPLX( cos( ar3 ), -sin( ar3 ) ,kind=DP)

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
            CALL errore(' phfacs ',' inconsistent size for eigr ',ngw)
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
      use fft_base,         only: dfftp

      IMPLICIT NONE

      ! ... declare subroutine arguments
      !
      COMPLEX(DP) :: ei1( -dfftp%nr1 : dfftp%nr1, nat )
      COMPLEX(DP) :: ei2( -dfftp%nr2 : dfftp%nr2, nat )
      COMPLEX(DP) :: ei3( -dfftp%nr3 : dfftp%nr3, nat )
      INTEGER      :: mill( :, : )
      INTEGER      :: ngm
      COMPLEX(DP), INTENT(OUT) :: sfac(:,:)

      ! ... declare other variables
      !
      INTEGER :: is, ig, ia, ig1, ig2, ig3, isa

      call start_clock( 'strucf' )

!$omp parallel do default(shared), private(ig1,ig2,ig3,isa,is,ia)
      DO ig = 1, ngm
        ig1 = mill( 1, ig ) 
        ig2 = mill( 2, ig ) 
        ig3 = mill( 3, ig )
        isa = 1
        DO is = 1, nsp
          sfac( ig, is ) = (0.0d0, 0.0d0)
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


