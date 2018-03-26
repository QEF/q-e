  !
  ! Copyright (C) 2010-2017 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2001-2008 Quantum-Espresso group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
!-------------------------------------------------------------------------------
SUBROUTINE rgd_blk_epw_fine_mem(imode,nq1,nq2,nq3,q,uq,epmat,nmodes,epsil,zeu,bmat,signe)
!-------------------------------------------------------------------------------
  !!
  !! Compute the long range term for the e-ph vertex
  !! to be added or subtracted from the vertex
  !!
  !! The long-range part can be computed using Eq. (4) of PRL 115, 176401 (2015).
  !! The sum over G is converged using the Ewald summation technique (see for example 
  !! F.2, p.500 in Martin Electronic structure book) where the Ewald factor is ((q+G)**2)/alph/4.0_DP.
  !!
  !! Technical note: From the solution of the Poisson equation, there is an additional factor 
  !! e^{-i(q+G)\tau_\kappa} with respect to Eq. (4) of PRL 115, 176401 (2015).
  !! The full equation can be found in Eq. (S4) of the supplemental materials of PRL 115, 176401 (2015).
  !! 
  !! The final implemented formula is:  
  !!
  !! $$ g_{mn\nu}^{\mathcal L}({\bf k},{\bf q) = i\frac{4\pi e^2}{\Omega} \sum_{\kappa}
  !!   \left(\frac{\hbar}{2 {M_\kappa \omega_{{\bf q}\nu}}}\right)^{\!\!\frac{1}{2}}
  !!   \sum_{{\bf G}\ne -{\bf q}} e^{-({\bf q}+{\bf G})^2/4\alpha}
  !! \frac{ ({\bf q}+{\bf G})\cdot{\bf Z}^*_\kappa \cdot {\bf e}_{\kappa\nu}({\bf q}) } 
  !!  {({\bf q}+{\bf G})\cdot\bm\epsilon^\infty\!\cdot({\bf q}+{\bf G})}\,
  !!   \left[ U_{{\bf k}+{\bf q}}\:U_{{\bf k}}^{\dagger} \right]_{mn} $$
  !!
  !! 10/2016 - SP: Optimization  
  !!
  USE kinds,         ONLY : dp
  USE cell_base,     ONLY : bg, omega, alat
  USE ions_base,     ONLY : tau, nat
  USE constants_epw, ONLY : twopi, fpi, e2, ci, czero
  USE epwcom,        ONLY : shortrange, nbndsub
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nq1
  !! Coarse q-point grid 
  INTEGER, INTENT (in) :: nq2
  !! Coarse q-point grid 
  INTEGER, INTENT (in) :: nq3
  !! Coarse q-point grid 
  INTEGER, INTENT (in) :: nmodes
  !! Max number of modes
  ! 
  REAL (kind=DP), INTENT (in) :: q(3)
  !! q-vector from the full coarse or fine grid.
  REAL (kind=DP), INTENT (in) :: epsil(3,3)
  !! dielectric constant tensor
  REAL (kind=DP), INTENT (in) :: zeu(3,3,nat)
  !! effective charges tensor
  REAL (kind=DP), INTENT (in) :: signe
  !! signe=+/-1.0 ==> add/subtract long range term
  ! 
  COMPLEX (kind=DP), INTENT (in) :: uq(nmodes, nmodes)
  !! phonon eigenvec associated with q
  COMPLEX (kind=DP), INTENT (inout) :: epmat(nbndsub,nbndsub)
  !! e-ph matrix elements 
  COMPLEX (kind=DP), INTENT (in) :: bmat(nbndsub,nbndsub) 
  !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
  !
  ! work variables
  !
  REAL(kind=DP) :: qeq,     &! <q+G| epsil | q+G>
       arg, zaq, g1, g2, g3, gmax, alph, geg
  INTEGER :: na, ipol, m1,m2,m3, imode
  COMPLEX(kind=DP) :: fac, facqd, facq
  COMPLEX(kind=DP) :: epmatl(nbndsub,nbndsub)
  !
  IF (abs(signe) /= 1.0) &
       CALL errore ('rgd_blk',' wrong value for signe ',1)
  !
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  fac = signe*e2*fpi/omega * ci
  !
  epmatl(:,:) = czero   
  !
  DO m1 = -nq1,nq1
    DO m2 = -nq2,nq2
      DO m3 = -nq3,nq3
      !
      g1 = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3) + q(1)
      g2 = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3) + q(2)
      g3 = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3) + q(3)
      !
      qeq = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3 )+      &
             g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3 )+      &
             g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3 )) !*twopi/alat
      !
      IF (qeq > 0.0_DP .and. qeq/alph/4.0_DP < gmax ) THEN
        !
        qeq=qeq*twopi/alat
        facqd = fac*exp(-qeq/alph/4.0d0)/qeq !/(two*wq)
        !
        DO na = 1,nat
          arg = -twopi* ( g1*tau(1,na)+ g2*tau(2,na)+ g3*tau(3,na) )
          facq = facqd * CMPLX(cos(arg),sin(arg),kind=DP)
          DO ipol=1,3
            zaq=g1*zeu(1,ipol,na)+g2*zeu(2,ipol,na)+g3*zeu(3,ipol,na)
            !
            CALL zaxpy(nbndsub**2,facq * zaq * uq(3*(na-1)+ipol,imode), bmat(:,:),1, epmat(:,:),1)
            CALL zaxpy(nbndsub**2,facq * zaq * uq(3*(na-1)+ipol,imode), bmat(:,:),1, epmatl(:,:),1)
            !
          ENDDO !ipol
        ENDDO !nat
      ENDIF
      !
      ENDDO
    ENDDO
  ENDDO
  !
  ! In case we want only the short-range we do
  ! g_s = sqrt(g*g - g_l*g_l)
  ! 
  ! Important notice: It is possible that (g*g - g_l*g_l) < 0, in which 
  ! case the sqrt will give an pure imaginary number. If it is positive we 
  ! will get a pure real number.
  ! In any case, when g_s will be squared both will become real numbers. 
  IF (shortrange) THEN
    !epmat = ZSQRT(epmat*conjg(epmat) - epmatl*conjg(epmatl))
    epmat = SQRT(epmat*conjg(epmat) - epmatl*conjg(epmatl))
  ENDIF        
  !
  !
END SUBROUTINE rgd_blk_epw_fine_mem
