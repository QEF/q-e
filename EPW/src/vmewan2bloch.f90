  ! 
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !
  !--------------------------------------------------------------------------
  subroutine vmewan2bloch ( nbnd, nrr, irvec, ndegen, xk, cuf, vmef, et, et_ks, chw)
  !--------------------------------------------------------------------------
  !
  !  From the Velocity matrix elements in Wannier representation, find the corresponding
  !  MEs in Bloch representation for a given k point
  !  
  !  input  : nbnd, nrr, irvec, ndegen, xk, cuf, et
  !
  !  output : vmef; velocity matrix elements on the fine mesh
  !
  !  Adapted from hamwan2bloch by Jesse Noffsinger and Emmanouil Kioupakis
  !
  !--------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  use elph2,         ONLY : cvmew !, chw
  use cell_base,     ONLY : at, celldm
  USE epwcom,        ONLY : eig_read
  USE constants_epw, ONLY : twopi, ci, czero
  !
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nrr, irvec (3, nrr), ndegen (nrr), ipol
  ! number of bands (possibly of the optimal subspace)
  ! kpoint number for the interpolation
  ! record length and unit for direct write of rotation matrix
  ! number of WS points, crystal coordinates, degeneracy
  !
  ! Hamiltonian in wannier basis
  !
  real(kind=DP) :: xk (3), et(nbnd),et_ks(nbnd),irvec_tmp(3)
  ! kpoint coordinates for the interpolation
  !
  ! output variables
  !  
  complex(kind=DP) :: chf_a(3,nbnd, nbnd), chf_a_tmp(nbnd, nbnd)
  complex(kind=DP) :: vmef (3,nbnd,nbnd)
  ! interpolated hamiltonian eigenvalues for this kpoint 
  complex(kind=DP) :: cuf(nbnd, nbnd), chw(nbnd, nbnd, nrr)
  ! Rotation matrix, fine mesh 
  !
  !
  complex(kind=DP) :: cvmef(3,nbnd, nbnd), cvmef_tmp(nbnd, nbnd)
  ! Hamiltonian in Bloch basis, fine mesh
  integer :: ibnd, jbnd, ir
  real(kind=DP) :: rdotk
  complex(kind=DP) :: cfac
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
  !  
  cvmef (:,:,:) = czero
  !
  DO ir = 1, nrr
     !
     rdotk = twopi * dot_product( xk, dble(irvec( :, ir) ))
     cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
     DO ipol = 1,3
        cvmef (ipol,:,:) = cvmef (ipol,:,:) + cfac * cvmew (ipol, :,:, ir )
     ENDDO
     !
  ENDDO
  !
  !
  ! vmn(k) = U v(amn)~ U^dagger
  !cuf,  passed from hamwan2bloch.
  DO ipol = 1,3
     cvmef_tmp(:,:) = matmul( cvmef(ipol,:,:) ,  conjg(transpose(cuf(:,:)))  )
     vmef(ipol, :,:) = matmul(cuf(:,:) , cvmef_tmp(:,:) )
  ENDDO
  !
  !
  !
  !  get k-derivative of the Hamiltonian in the Wannier gauge
  !
  chf_a (:,:,:) = czero
  !
  DO ir = 1, nrr
     !
     rdotk = twopi * dot_product( xk, dble(irvec( :, ir) ))
     cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
     irvec_tmp(:) = celldm(1) * matmul ( at, dble(irvec(:,ir)) )
     DO ipol = 1, 3
        chf_a (ipol,:,:) = chf_a (ipol,:,:) + &
             ci * irvec_tmp( ipol ) * cfac * chw (:,:, ir )
     ENDDO
     !
  ENDDO
  !
  ! H'mn(k) = U H'~ U^dagger
  ! cuf,  passed from hamwan2bloch.
  DO ipol = 1,3
     chf_a_tmp(:,:) = matmul( chf_a(ipol,:,:) ,  conjg(transpose(cuf(:,:)))  )
     chf_a(ipol, :,:) = matmul(cuf(:,:) , chf_a_tmp(:,:) )
  ENDDO
  !
  !
  DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
        vmef (:,ibnd,jbnd) = chf_a(:, ibnd, jbnd) - &
             ci * (et_ks(jbnd) - et_ks(ibnd) ) *  vmef(:,ibnd,jbnd)
     ENDDO
  ENDDO
  !
  IF (eig_read) THEN
     DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
           IF (abs(et_ks(ibnd) - et_ks(jbnd)) .gt. 1.d-4) THEN
              vmef(:, ibnd, jbnd) = vmef(:,ibnd, jbnd) * &
                   ( et(ibnd)    - et(jbnd) )/ &
                   ( et_ks(ibnd) - et_ks(jbnd) )
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  end subroutine vmewan2bloch
  !
