  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dynwan2bloch ( nmodes, nrr, irvec, ndegen, xxq, cuf, eig)
  !--------------------------------------------------------------------------
  !!
  !!
  !!  WARNING: this SUBROUTINE is identical to hamwan2bloch.f90, except
  !!           that here rdw is a real array, not a complex one. This is
  !!           required to obtain proper phonon dispersion interpolation
  !!           and corresponds to the reality of the interatomic force
  !!           constants
  !!
  !! -------------------------------------------------------------------------
  !!
  !!  From the Hamiltonian in Wannier representation, find the corresponding
  !!  Hamiltonian in Bloch representation for a given k point
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg
  USE phcom,     ONLY : nq1, nq2, nq3
  USE ions_base, ONLY : amass, tau, nat, ityp
  USE elph2,     ONLY : rdw, epsi, zstar
  USE epwcom,    ONLY : lpolar
  USE constants_epw, ONLY : twopi, ci, czero
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nmodes
  !! number of modes (possibly of the optimal subspace)
  INTEGER, INTENT (in) :: nrr
  !! kpoint number for the interpolation
  INTEGER, INTENT (in) :: irvec (3, nrr)
  !! record length and unit for direct write of rotation matrix
  INTEGER, INTENT (in) :: ndegen (nrr)
  !! number of WS points, crystal coordinates, degeneracy
  !
  REAL(kind=DP), INTENT (in) :: xxq (3)
  !! kpoint coordinates for the interpolation
  REAL(kind=DP), INTENT (out) :: eig (nmodes)
  !! interpolated hamiltonian eigenvalues for this kpoint
  COMPLEX(kind=DP), INTENT (out) :: cuf(nmodes, nmodes)
  !! Rotation matrix, fine mesh 
  !
  ! work variables
  ! variables for lapack ZHPEVX
  integer :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
  real(kind=DP) :: w( nmodes ), rwork( 7*nmodes )
  complex(kind=DP) :: champ( nmodes*(nmodes+1)/2 ), &
    cwork( 2*nmodes ), cz( nmodes, nmodes)
  !
  real(kind=DP) :: xq (3)
  complex(kind=DP) :: chf(nmodes, nmodes)
  ! Hamiltonian in Bloch basis, fine mesh
  integer :: imode, jmode, ir, na, nb
  real(kind=DP) :: rdotk, massfac
  complex(kind=DP) :: cfac
  !
  CALL start_clock ( 'DynW2B' )
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform to fine k and k+q meshes
  !----------------------------------------------------------
  !
  !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
  !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
  !
  !  H~_k   is chf ( nmodes, nmodes, 2*ik-1 )
  !  H~_k+q is chf ( nmodes, nmodes, 2*ik   )
  !
  xq = xxq
  chf (:,:) = czero
  !
  DO ir = 1, nrr
    !
    rdotk = twopi * dot_product( xq, dble(irvec( :, ir) ))
    cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
    chf = chf + cfac * rdw (:,:, ir )
    !
  ENDDO
  !
  ! bring xq in cart. coordinates (needed for rgd_blk call)
  CALL cryst_to_cart (1, xq, bg, 1)
  !
  !  add the long-range term to D(q)
  IF (lpolar) THEN
    ! xq has to be in 2pi/a     
    CALL rgd_blk (nq1,nq2,nq3,nat,chf,xq,tau,epsi,zstar,+1.d0)
    !
  ENDIF
  !
  !  divide by the square root of masses 
  !
  DO na = 1, nat
    DO nb = 1, nat
      massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
      !
      chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
         chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) * massfac
      ! 
    ENDDO
  ENDDO
  !
  ! bring xq back to crystal coordinates
  CALL cryst_to_cart (1, xq, at, -1)
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jmode = 1, nmodes
   DO imode = 1, jmode
      champ (imode + (jmode - 1) * jmode/2 ) = &
      ( chf ( imode, jmode) + conjg ( chf ( jmode, imode) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nmodes, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nmodes, cwork, &
               rwork, iwork, ifail, info)
  !
  ! rotation matrix and Ham eigenvalues
  ! [in Ry, mind when comparing with wannier code]
  !
  cuf = cz
  eig = w
  !
  CALL stop_clock ( 'DynW2B' )
  !   
  END SUBROUTINE dynwan2bloch

