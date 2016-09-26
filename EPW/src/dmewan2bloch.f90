  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dmewan2bloch ( nbnd, nrr, irvec, ndegen, xk, cuf, dmef, etf, etf_ks, cfac)
  !--------------------------------------------------------------------------
  !!
  !!  From the Dipole in Wannier representation, find the corresponding
  !!  Dipole in Bloch representation for a given k point
  !!  
  !!  input  : number of bands nbnd
  !!           number of WS vectors, coordinates and degeneracy 
  !!           Dipole in Wannier representation  cdmew(3,nbnd,nbnd,nrr)
  !!           kpoint coordinate xk(3)
  !!
  !!  output : interpolated dipole matrix elements (dmef)
  !!
  !!  SP 09/2016: Optimization
  !!  JN, EK 09/2010
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : cdmew
  USE epwcom,        ONLY : eig_read
  USE constants_epw, ONLY : twopi, ci, cone, czero
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nbnd 
  !! number of bands (possibly of the optimal subspace)
  INTEGER, INTENT (in) :: nrr 
  !! kpoint number for the interpolation
  INTEGER, DIMENSION( 1:3, 1:nrr ), INTENT (in) :: irvec
  !! record length and unit for direct write of rotation matrix
  INTEGER, DIMENSION( 1:nrr      ), INTENT (in) :: ndegen
  !! Number of degeneracies
  !
  REAL(kind=DP), INTENT (in) :: xk (3)
  !!  kpoint coordinates for the interpolation
  REAL(kind=DP), INTENT (in) :: etf (nbnd)
  !! Eigenenergies on the fine grid
  REAL(kind=DP), INTENT (in) :: etf_ks (nbnd) 
  !! 
  COMPLEX(kind=DP), INTENT (in) :: cuf (nbnd, nbnd)
  !! Rotation matrix, fine mesh 
  COMPLEX(kind=DP), INTENT (out) :: dmef (3, nbnd, nbnd)
  !! interpolated dipole matrix elements in Bloch basis, fine mesh
  COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
  ! 
  ! local variables
  !
  INTEGER :: ir, ibnd, jbnd, j 
  REAL   ( kind=DP ) :: rdotk
  COMPLEX( kind=DP ) :: congj_cuf (nbnd, nbnd)
  COMPLEX( kind=DP ) :: cdmef (3, nbnd, nbnd)
  COMPLEX( kind=DP ) :: cdmef_tmp (3, nbnd, nbnd)
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
  ! Initialization
  cdmef_tmp (:,:,:) = czero
  dmef (:,:,: ) = czero
  cdmef (:,:,:) = czero
  !   
  ! SUM on ir of: cdmef (1,ibnd,jbnd) = cdmef (1,ibnd,jbnd) + cfac( ir) * cdmew (1, ibnd,jbnd, ir )
  CALL zgemv('n', 3*(nbnd**2), nrr, cone, cdmew(:,:,:,:), 3*(nbnd**2), cfac(:), 1, cone, cdmef(:,:,:), 1  )
  ! 
  !
  ! pmn(k) = U p~ U^dagger
  ! cuf,  passed from hamwan2bloch.
  !
  congj_cuf = Conjg(Transpose(cuf))
  DO j=1, 3
    CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cdmef(j,:,:), &
         nbnd, congj_cuf(:,:), nbnd, czero, cdmef_tmp(j,:,:), nbnd)
  ENDDO
  !
  DO j=1, 3
    CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:,:), &
               nbnd, cdmef_tmp(j,:,:), nbnd, czero, dmef(j,:,:), nbnd)
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
