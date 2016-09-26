  ! 
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !   
  !--------------------------------------------------------------------------
  subroutine hamwan2bloch ( nbnd, nrr, cuf, eig, chw, cfac)
  !--------------------------------------------------------------------------
  !
  !  From the Hamiltonian in Wannier representation, find the corresponding
  !  Hamiltonian in Bloch representation for a given k point
  !  
  !  input  : number of bands nbnd
  !           number of WS vectors, coordinates and degeneracy 
  !           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
  !           kpoint coordinate xk(3)
  !
  !  output : rotation matrix cuf (nbnd, nbnd)
  !           interpolated hamiltonian eigenvalues eig(nbnd)
  !
  !  SP [optimization]
  !
  !  Feliciano Giustino, UCB
  !
  !--------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  !
  implicit none
  !
  !  input variables
  !
  INTEGER, INTENT (in) :: nbnd
  !! number of bands (possibly of the optimal subspace)
  INTEGER, INTENT (in) :: nrr
  !! number of WS points
  ! 
  REAL(kind=DP), INTENT (out) :: eig (nbnd)
  !! interpolated hamiltonian eigenvalues for this kpoint 
  ! 
  COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
  !! Exponential factor
  COMPLEX(kind=DP), INTENT (in) :: chw ( nbnd, nbnd, nrr)
  !! Hamiltonian in Bloch basis, fine mesh
  COMPLEX(kind=DP), INTENT (out) :: cuf(nbnd, nbnd)
  !! Rotation matrix, fine mesh
  !
  ! variables for lapack ZHPEVX 
  !
  INTEGER :: neig, info, ifail( nbnd ), iwork( 5*nbnd )
  REAL(kind=DP) :: w( nbnd )
  REAL(kind=DP) :: rwork( 7*nbnd )
  COMPLEX(kind=DP) :: champ( nbnd*(nbnd+1)/2 )
  COMPLEX(kind=DP) :: cwork( 2*nbnd )
  COMPLEX(kind=DP) :: cz( nbnd, nbnd)
  !
  ! work variables 
  !
  INTEGER :: ibnd, jbnd
  COMPLEX(kind=DP) :: chf(nbnd, nbnd) 
  !
  CALL start_clock('HamW2B')
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
  chf (:,:) = czero
  !
  ! Previous implementation
  !DO ir = 1, nrr
  !   !
  !   rdotk = twopi * ( xk(1)*irvec(1,ir) + xk(2)*irvec(2,ir) + xk(3)*irvec(3,ir))
  !   cfac = exp( ci*rdotk ) / ndegen(ir)
  !   !
  !   ! SP : Significantly faster !
  !   !chf (:,:) = chf (:,:) + cfac * chw (:,:, ir )
  !   CALL zaxpy(nbnd**2,cfac, chw (1,1, ir ),1, chf (1,1),1)
  !   !
  !ENDDO
  ! New one
  CALL zgemv('n', nbnd**2, nrr, cone, chw, nbnd**2, cfac, 1, cone, chf, 1 )
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jbnd = 1, nbnd
   DO ibnd = 1, jbnd
      champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
      ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) * 0.5d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nbnd, cwork, &
               rwork, iwork, ifail, info)
  !
  ! rotation matrix and Ham eigenvalues 
  ! [in Ry, mind when comparing with wannier code]
  ! 
  cuf = conjg( transpose ( cz ) )
  eig = w 
  !
  CALL stop_clock('HamW2B')
  !
  end subroutine hamwan2bloch
  !
