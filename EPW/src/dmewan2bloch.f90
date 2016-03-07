  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dmewan2bloch ( nbnd, nrr, irvec, ndegen, xk, cuf, dmef, etf, etf_ks)
  !--------------------------------------------------------------------------
  !
  !  From the Dipole in Wannier representation, find the corresponding
  !  Dipole in Bloch representation for a given k point
  !  
  !  input  : number of bands nbnd
  !           number of WS vectors, coordinates and degeneracy 
  !           Dipole in Wannier representation  cdmew(3,nbnd,nbnd,nrr)
  !           kpoint coordinate xk(3)
  !
  !  output : interpolated dipole matrix elements (dmef)
  !
  !  JN, EK 09/2010
  !
  !--------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : cdmew
  USE epwcom,        ONLY : eig_read
  USE constants_epw, ONLY : twopi, ci
  !
  implicit none
  !
  !  input/output variables
  !
  INTEGER, INTENT( IN ) :: nbnd, nrr 
  ! number of bands (possibly of the optimal subspace)
  ! kpoint number for the interpolation
  INTEGER, DIMENSION( 1:3, 1:nrr ), INTENT( IN ) :: irvec
  INTEGER, DIMENSION( 1:nrr      ), INTENT( In ) :: ndegen
  ! record length and unit for direct write of rotation matrix
  ! degeneracy
  REAL   ( kind=DP ), DIMENSION( 1:3                  ), INTENT( IN  ) :: xk
  !  kpoint coordinates for the interpolation
  COMPLEX( kind=DP ), DIMENSION( 1:nbnd, 1:nbnd       ), INTENT( IN  ) :: cuf
  ! Rotation matrix, fine mesh 
  COMPLEX( kind=DP ), DIMENSION( 1:3, 1:nbnd, 1:nbnd  ), INTENT( OUT ) :: dmef
  ! interpolated dipole matrix elements in Bloch basis, fine mesh
  REAL   ( kind=DP ), DIMENSION( 1:nbnd               ), INTENT( IN  ) :: etf
  REAL   ( kind=DP ), DIMENSION( 1:nbnd               ), INTENT( IN  ) :: etf_ks
  ! 
  ! local variables
  !
  INTEGER :: ir, ibnd, jbnd, j 
  REAL   ( kind=DP ) :: rdotk
  COMPLEX( kind=DP ) :: cfac, congj_cuf, cdtmp1, cdtmp2, cdtmp3
  COMPLEX( kind=DP ), DIMENSION( 1:3, 1:nbnd, 1:nbnd ) :: cdmef
  COMPLEX( kind=DP ), DIMENSION( 1:3, 1:nbnd, 1:nbnd ) :: cdmef_tmp
  ! dipole matrix elements in Bloch basis, fine mesh
  !
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform to fine k and k+q meshes
  !----------------------------------------------------------
  !
  !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
  !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
  !
  !  H~_k   is chf ( nbnd, nbnd, 2*ik-1 )
  !  H~_k+q is chf ( nbnd, nbnd, 2*ik   )
  !
  ir = 1
  rdotk = twopi * ( xk( 1 ) * irvec( 1, ir ) + xk( 2 ) * irvec( 2, ir ) + xk(3 ) * irvec( 3, ir ) )
  cfac = exp( ci*rdotk ) / ndegen(ir)
  DO jbnd = 1, nbnd
     DO ibnd = 1, nbnd
        !
        cdmef (1,ibnd,jbnd) = cfac * cdmew (1, ibnd,jbnd, ir )
        cdmef (2,ibnd,jbnd) = cfac * cdmew (2, ibnd,jbnd, ir )
        cdmef (3,ibnd,jbnd) = cfac * cdmew (3, ibnd,jbnd, ir )
        ! 
        cdmef_tmp( 1, ibnd, jbnd ) = 0.0_dp
        cdmef_tmp( 2, ibnd, jbnd ) = 0.0_dp
        cdmef_tmp( 3, ibnd, jbnd ) = 0.0_dp
        !
        dmef( 1, ibnd, jbnd ) = 0.0_dp
        dmef( 2, ibnd, jbnd ) = 0.0_dp
        dmef( 3, ibnd, jbnd ) = 0.0_dp
        !
     ENDDO
  ENDDO
  !
  DO ir = 2, nrr
     !
     rdotk = twopi * ( xk( 1 ) * irvec( 1, ir ) + xk( 2 ) * irvec( 2, ir ) + xk( 3 ) * irvec( 3, ir ) )
     cfac = Exp( ( 0.0_dp, 1.0_dp ) * rdotk ) / ndegen(ir)
     DO jbnd = 1, nbnd
        DO ibnd = 1, nbnd
           cdmef (1,ibnd,jbnd) = cdmef (1,ibnd,jbnd) + cfac * cdmew (1, ibnd,jbnd, ir )
           cdmef (2,ibnd,jbnd) = cdmef (2,ibnd,jbnd) + cfac * cdmew (2, ibnd,jbnd, ir )
           cdmef (3,ibnd,jbnd) = cdmef (3,ibnd,jbnd) + cfac * cdmew (3, ibnd,jbnd, ir )
        ENDDO
     ENDDO
     !
  ENDDO
  !
  ! pmn(k) = U p~ U^dagger
  ! cuf,  passed from hamwan2bloch.
  !
  DO j = 1, nbnd
     DO jbnd = 1, nbnd
        congj_cuf = Conjg( cuf( jbnd, j ) )
        DO ibnd = 1, nbnd
           cdmef_tmp( 1, ibnd, jbnd ) = cdmef_tmp( 1, ibnd, jbnd ) + cdmef( 1, ibnd, j ) * congj_cuf
           cdmef_tmp( 2, ibnd, jbnd ) = cdmef_tmp( 2, ibnd, jbnd ) + cdmef( 2, ibnd, j ) * congj_cuf
           cdmef_tmp( 3, ibnd, jbnd ) = cdmef_tmp( 3, ibnd, jbnd ) + cdmef( 3, ibnd, j ) * congj_cuf
        ENDDO
     ENDDO
  ENDDO
  ! 
  DO j = 1, nbnd
     DO jbnd = 1, nbnd
        cdtmp1 = cdmef_tmp( 1, j, jbnd )
        cdtmp2 = cdmef_tmp( 2, j, jbnd )
        cdtmp3 = cdmef_tmp( 3, j, jbnd )
        DO ibnd = 1, nbnd
           dmef( 1, ibnd, jbnd ) =  dmef( 1, ibnd, jbnd ) + cuf( ibnd, j ) * cdtmp1
           dmef( 2, ibnd, jbnd ) =  dmef( 2, ibnd, jbnd ) + cuf( ibnd, j ) * cdtmp2
           dmef( 3, ibnd, jbnd ) =  dmef( 3, ibnd, jbnd ) + cuf( ibnd, j ) * cdtmp3
        ENDDO
     ENDDO
  ENDDO
  !
  ! Satisfy
  ! Phys. Rev. B 62, 4927–4944 (2000) , Eq. (30)
  !
  IF (eig_read) THEN
     DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
           IF (abs(etf_ks(ibnd) - etf_ks(jbnd)) .gt. 1.d-4) THEN
              dmef(:, ibnd, jbnd) = dmef(:,ibnd, jbnd) * &
                   ( etf(ibnd)    - etf(jbnd) )/ &
                   ( etf_ks(ibnd) - etf_ks(jbnd) )
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  END SUBROUTINE dmewan2bloch
  !
