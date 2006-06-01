!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!=----------------------------------------------------------------------------=!
   MODULE forces
!=----------------------------------------------------------------------------=!

       USE kinds
       USE cell_base, ONLY: tpiba2

       IMPLICIT NONE
       SAVE
 
       PRIVATE

! ... i^l imaginary unit to the angular momentum
       COMPLEX(DP), PARAMETER :: cimgl(0:3) = (/ (1.0d0,0.0d0), &
         (0.0d0,1.0d0), (-1.0d0,0.0d0), (0.0d0,-1.0d0) /)
       COMPLEX(DP), PARAMETER :: czero = (0.0_DP,0.0_DP)
       REAL(DP),    PARAMETER :: rzero =  0.0_DP

       PUBLIC :: dforce, dforce_all

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!



    SUBROUTINE dforce1( co, ce, dco, dce, fio, fie, hg, v, psi_stored )

      USE fft_base,   ONLY: dffts
      USE gvecw,      ONLY: ngw
      USE fft_module, ONLY: fwfft, invfft

      IMPLICIT NONE

      ! ... declare subroutine arguments
      COMPLEX(DP), INTENT(OUT) :: dco(:), dce(:)
      COMPLEX(DP), INTENT(IN)  :: co(:), ce(:)
      REAL(DP),    INTENT(IN)  :: fio, fie
      REAL(DP),    INTENT(IN)  :: v(:)
      REAL(DP),    INTENT(IN)  :: hg(:)
      COMPLEX(DP), OPTIONAL    :: psi_stored(:)

      ! ... declare other variables
      !
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP) :: fp, fm, aro, are
      REAL(DP)    :: fioby2, fieby2, arg
      INTEGER      :: ig

      !  end of declarations

      ALLOCATE( psi( SIZE(v) ) )

      IF( PRESENT( psi_stored ) ) THEN
        psi = psi_stored * CMPLX(v, 0.0d0)
      ELSE
        CALL c2psi( psi, dffts%nnr, co, ce, ngw, 2 )
        CALL invfft( 'Wave', psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )
        psi = psi * CMPLX(v, 0.0d0)
      END IF

      CALL fwfft( 'Wave', psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )
      CALL psi2c( psi, dffts%nnr, dco, dce, ngw, 2 )

      DEALLOCATE(psi)

      fioby2   = fio * 0.5
      fieby2   = fie * 0.5

      DO ig = 1, SIZE(co)
        fp = dco(ig) + dce(ig)
        fm = dco(ig) - dce(ig)
        aro = CMPLX(  DBLE(fp), AIMAG(fm) )
        are = CMPLX( AIMAG(fp), -DBLE(fm))
        arg = tpiba2 * hg(ig)
        dco(ig) = -fioby2 * (arg * co(ig) + aro)
        dce(ig) = -fieby2 * (arg * ce(ig) + are)
      END DO

    RETURN
    END SUBROUTINE dforce1


!=----------------------------------------------------------------------------=!


    SUBROUTINE dforce2_bec( fio, fie, df, da, eigr, beco, bece )

        !  this routine computes:
        !  the generalized force df=CMPLX(dfr,dfi) acting on the i-th
        !  electron state at the ik-th point of the Brillouin zone
        !  represented by the vector c=CMPLX(cr,ci)
        !  ----------------------------------------------

      USE ions_base,       ONLY: na
      USE pseudopotential, ONLY: nspnl
      USE uspp_param,      only: nh
      USE uspp,            only: nhtol, nhtolm, indv, beta, dvan
      use cvan,            only: ish


      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(DP), INTENT(IN) :: eigr(:,:)
      REAL(DP), INTENT(IN) :: fio, fie
      COMPLEX(DP)  :: df(:), da(:)
      REAL(DP), INTENT(IN) :: beco(:)
      REAL(DP), INTENT(IN) :: bece(:)

! ... declare other variables
      COMPLEX(DP), ALLOCATABLE :: temp(:,:)
      REAL(DP) :: t1
      REAL(DP) :: sgn
      INTEGER   :: l, is, ig, ngw, iv, inl, isa

!  end of declarations
!  ----------------------------------------------

      ngw  = SIZE(df)
      ALLOCATE(temp(ngw,2))

      isa = 1
      
      DO is = 1, nspnl
        !
        DO iv = 1, nh( is )
          !
          l   = nhtol ( iv, is )
          inl = ish(is) + (iv-1) * na(is) + 1

          sgn = 1.0d0
          IF( MOD( l, 2 ) /= 0 ) sgn = -1.0d0   !  ( -1)^l
          
          t1= - fio * dvan( iv, iv, is ) * sgn
          !
          CALL DGEMV('N', 2*ngw, na(is), t1, eigr(1,isa), &
               2*SIZE(eigr,1), beco( inl ), 1, rzero, temp(1,1), 1)
          !
          CALL ZSCAL( ngw, cimgl(l), temp(1,1), 1)
          !
          t1= - fie * dvan( iv, iv, is ) * sgn
          CALL DGEMV('N', 2*ngw, na(is), t1, eigr(1,isa), &
               2*SIZE(eigr,1), bece( inl ), 1, rzero, temp(1,2), 1)
          !
          CALL ZSCAL( ngw, cimgl(l), temp(1,2), 1)
          !
          DO ig=1,ngw
             df(ig) = df(ig) + temp(ig,1) * beta(ig,iv,is)
          END DO
          DO ig=1,ngw
             da(ig) = da(ig) + temp(ig,2) * beta(ig,iv,is)
          END DO
        END DO
        !
        isa = isa + na( is )
        !
      END DO

      DEALLOCATE(temp)

    RETURN
    END SUBROUTINE dforce2_bec



!=----------------------------------------------------------------------------=!

     

    SUBROUTINE dforce( ib, iss, c, f, df, da, v, eigr, bec, nupdwn, iupdwn )
       !
       USE reciprocal_vectors, ONLY: ggp, g, gx
       !
       IMPLICIT NONE
       !
       INTEGER,      INTENT(IN) :: ib, iss     ! band and spin index
       COMPLEX(DP), INTENT(IN)  :: c(:,:)
       COMPLEX(DP), INTENT(OUT) :: df(:), da(:)
       REAL (DP),   INTENT(IN)  :: v(:), bec(:,:), f(:)
       COMPLEX(DP), INTENT(IN)  :: eigr(:,:)
       INTEGER,      INTENT(IN) :: nupdwn(:), iupdwn(:)
       !
       COMPLEX(DP), ALLOCATABLE :: dum( : )   
       !
       INTEGER :: ig, in
       !
       IF( ib > nupdwn( iss ) ) CALL errore( ' dforce ', ' ib out of range ', 1 )
       !
       in = iupdwn( iss ) + ib - 1 
       !
       IF( ib == nupdwn( iss ) ) THEN
          !
          ALLOCATE( dum( SIZE( da ) ) )
          !
          CALL dforce1( c(:,ib), c(:,ib), df, dum, f(ib), f(ib), ggp, v )
          !
          CALL dforce2_bec( f(ib), f(ib), df , dum , eigr, bec( :, in ), bec( :, in ) )
          !
          DEALLOCATE( dum )
          !
       ELSE
          !
          CALL dforce1( c(:,ib), c(:,ib+1), df, da, f(ib), f(ib+1), ggp, v )
          !
          CALL dforce2_bec( f(ib), f(ib+1), df, da, eigr, bec( :, in ), bec( :, in+1 ) )
          !
       END IF
       !
       return
    END SUBROUTINE dforce


!  ----------------------------------------------
  

    SUBROUTINE dforce_all( iss, c, f, cgrad, vpot, eigr, bec, nupdwn, iupdwn )
        !
        IMPLICIT NONE

        INTEGER,               INTENT(IN)    :: iss
        COMPLEX(DP),           INTENT(INOUT) :: c(:,:)
        REAL(DP),              INTENT(IN)    :: vpot(:), f(:)
        COMPLEX(DP),           INTENT(OUT)   :: cgrad(:,:)
        COMPLEX(DP),           INTENT(IN)    :: eigr(:,:)
        REAL(DP),              INTENT(IN)    :: bec(:,:)
        INTEGER,               INTENT(IN)    :: nupdwn(:), iupdwn(:)
       
        INTEGER :: ib
        !
        IF( nupdwn( iss ) > 0 ) THEN
           !
           !   Process two states at the same time
           !
           DO ib = 1, nupdwn( iss )-1, 2
              CALL dforce( ib, iss, c, f, cgrad(:,ib), cgrad(:,ib+1), &
                  vpot, eigr, bec, nupdwn, iupdwn )
           END DO
           !
           IF( MOD( nupdwn( iss ), 2 ) /= 0 ) THEN
              ib = nupdwn( iss )
              CALL dforce( ib, iss, c, f, cgrad(:,ib), cgrad(:,ib), &
                  vpot, eigr, bec, nupdwn, iupdwn )
           END IF
           !
        END IF
        !
        RETURN
    END SUBROUTINE dforce_all



 END MODULE forces
