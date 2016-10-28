  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2001-2008 Quantum-Espresso group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !-----------------------------------------------------------------------
  SUBROUTINE rgd_blk (nq1,nq2,nq3,nat,dyn,q,tau,epsil,zeu,signe)
  !-----------------------------------------------------------------------
  !! This is adapted from QE PH/rigid.f90 
  !!
  !! compute the rigid-ion (long-range) term for q 
  !! The long-range term used here, to be added to or subtracted from the
  !! dynamical matrices, is exactly the same of the formula introduced
  !! in:
  !! X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is 
  !! implemented: the Ewald parameter alpha must be large enough to 
  !! have negligible r-space contribution
  !!
  USE kinds,         ONLY : dp
  USE constants_epw, ONLY : pi, fpi, e2
  USE cell_base,     ONLY : bg, omega
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nq1
  !! Coarse q-point grid 
  INTEGER, INTENT (in) :: nq2
  !! Coarse q-point grid 
  INTEGER, INTENT (in) :: nq3
  !! Coarse q-point grid 
  INTEGER, INTENT (in) :: nat
  !! Number of atoms
  ! 
  REAL (kind=DP), INTENT (in) :: q(3)
  !! q-vector from the full coarse or fine grid.
  REAL (kind=DP), INTENT (in) :: epsil(3,3)
  !! dielectric constant tensor
  REAL (kind=DP), INTENT (in) :: zeu(3,3,nat)
  !! effective charges tensor
  REAL (kind=DP), INTENT (in) :: signe
  !! signe=+/-1.0 ==> add/subtract rigid-ion term
  REAL (kind=DP), INTENT (in) :: tau(3,nat)
  !! Atomic positions
  ! 
  COMPLEX (kind=DP), INTENT (inout) :: dyn(3*nat,3*nat)
  !! Dynamical matrix
  !
  ! local variables
  real(DP):: geg                    !  <q+G| epsil | q+G>
  integer :: na,nb, i,j, m1, m2, m3
  integer :: nrx1, nrx2, nrx3
  real(DP) :: alph, fac,g1,g2,g3, facgd, arg, gmax
  real(DP) :: zag(3),zbg(3),zcg(3), fnat(3)
  complex(dp) :: facg
  real(DP) :: eps=1.0d-6
  !
  ! alph is the Ewald parameter, geg is an estimate of G^2
  ! such that the G-space sum is convergent for that alph
  ! very rough estimate: geg/4/alph > gmax = 14 
  ! (exp (-14) = 10^-6)
  !
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  !
  ! Estimate of nrx1,nrx2,nrx3 generating all vectors up to G^2 < geg
  ! Only for dimensions where periodicity is present, e.g. if nr1=1 
  ! and nr2=1, then the G-vectors run along nr3 only.
  ! (useful if system is in vacuum, e.g. 1D or 2D)
  !
  IF (nq1 == 1) THEN 
    nrx1=0
  ELSE
    nrx1 = int ( sqrt (geg) / &
                 sqrt (bg (1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2) ) + 1
  ENDIF
  IF (nq2 == 1) THEN 
    nrx2=0
  ELSE
    nrx2 = int ( sqrt (geg) / &
                 sqrt (bg (1, 2) **2 + bg (2, 2) **2 + bg (3, 2) **2) ) + 1
  ENDIF
  IF (nq3 == 1) THEN 
    nrx3=0
  ELSE
    nrx3 = int ( sqrt (geg) / &
                 sqrt (bg (1, 3) **2 + bg (2, 3) **2 + bg (3, 3) **2) ) + 1
  ENDIF
  !
  IF (abs(abs(signe) - 1.0) > eps) &
       CALL errore ('rgd_blk',' wrong value for signe ',1)
  !
  fac = signe*e2*fpi/omega
  DO m1 = -nrx1,nrx1
    DO m2 = -nrx2,nrx2
      DO m3 = -nrx3,nrx3
        !
        g1 = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3)
        g2 = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3)
        g3 = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3)
        !
        geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
               g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
               g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
        !
        IF (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) THEN
          !
          facgd = fac*exp(-geg/alph/4.0d0)/geg
          !
          DO na = 1,nat
            zag(:)=g1*zeu(1,:,na)+g2*zeu(2,:,na)+g3*zeu(3,:,na)
            fnat(:) = 0.d0
            DO nb = 1,nat
              arg = 2.d0*pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                              g2 * (tau(2,na)-tau(2,nb))+             &
                              g3 * (tau(3,na)-tau(3,nb)))
              zcg(:) = g1*zeu(1,:,nb) + g2*zeu(2,:,nb) + g3*zeu(3,:,nb)
              fnat(:) = fnat(:) + zcg(:)*cos(arg)
            ENDDO
            DO j=1,3 
              DO i=1,3 
                dyn( (na-1)*3+i,(na-1)*3+j ) = dyn((na-1)*3+i,(na-1)*3+j) - facgd * &
                                 zag(i) * fnat(j) 
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        g1 = g1 + q(1)
        g2 = g2 + q(2)
        g3 = g3 + q(3)
        !
        geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
               g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
               g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
        !
        IF (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) THEN
          !
          facgd = fac*exp(-geg/alph/4.0d0)/geg
          !
          DO nb = 1,nat
            zbg(:)=g1*zeu(1,:,nb)+g2*zeu(2,:,nb)+g3*zeu(3,:,nb)
            DO na = 1,nat
              zag(:)=g1*zeu(1,:,na)+g2*zeu(2,:,na)+g3*zeu(3,:,na)
              arg = 2.d0*pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                              g2 * (tau(2,na)-tau(2,nb))+             &
                              g3 * (tau(3,na)-tau(3,nb)))
              !
              facg = facgd * CMPLX(cos(arg),sin(arg),kind=DP)
              DO j=1,3 
                DO i=1,3 
                  dyn( (na-1)*3+i,(nb-1)*3+j ) = dyn((na-1)*3+i,(nb-1)*3+j) + facg * &
                                   zag(i) * zbg(j)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      !
      ENDDO
    ENDDO
  ENDDO
  !
END SUBROUTINE rgd_blk
!
!
!-------------------------------------------------------------------------------
SUBROUTINE rgd_blk_epw(nq1,nq2,nq3,q,uq,epmat,nmodes,epsil,zeu,bmat,signe)
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
  USE kinds,         ONLY : dp
  USE cell_base,     ONLY : bg, omega, alat
  USE ions_base,     ONLY : tau, nat
  USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, cone, two, ryd2mev
  USE epwcom,        ONLY : shortrange
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
  COMPLEX (kind=DP), INTENT (inout) :: epmat(nmodes)
  !! e-ph matrix elements 
  COMPLEX (kind=DP), INTENT (in) :: bmat 
  !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
  !
  ! work variables
  !
  real(DP) :: qeq,     &! <q+G| epsil | q+G>
       arg, zaq,       &
       g1, g2, g3, gmax, alph, geg
  integer :: na, ipol, im, m1,m2,m3!, nrx1,nrx2,nrx3
  complex(dp) :: fac, facqd, facq, epmatl(nmodes), matsq
  !
  IF (abs(signe) /= 1.0) &
       CALL errore ('rgd_blk',' wrong value for signe ',1)
  !
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  fac = signe*e2*fpi/omega * ci
  !
 ! IF (nq1 == 1) THEN 
 !    nrx1=0
 ! ELSE
 !    nrx1 = int ( sqrt (geg) / &
 !                 sqrt (bg (1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2) ) + 1
 ! ENDIF
 ! IF (nq2 == 1) THEN 
 !    nrx2=0
 ! ELSE
 !    nrx2 = int ( sqrt (geg) / &
 !                 sqrt (bg (1, 2) **2 + bg (2, 2) **2 + bg (3, 2) **2) ) + 1
 ! ENDIF
 ! IF (nq3 == 1) THEN 
 !    nrx3=0
 ! ELSE
 !    nrx3 = int ( sqrt (geg) / &
 !                 sqrt (bg (1, 3) **2 + bg (2, 3) **2 + bg (3, 3) **2) ) + 1
 ! ENDIF
  !
  epmatl(:) = czero   
  !
  !DO m1 = -nrx1,nrx1
  ! TO be test
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
            epmat = epmat + facq * zaq * uq(3*(na-1)+ipol,:)*bmat
            epmatl = epmatl + facq * zaq * uq(3*(na-1)+ipol,:)*bmat
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
END SUBROUTINE rgd_blk_epw
!
!
!-------------------------------------------------------------------------------
SUBROUTINE rgd_blk_epw_fine(nq1,nq2,nq3,q,uq,epmat,nmodes,epsil,zeu,bmat,signe)
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
  USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, cone, two, ryd2mev
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
  COMPLEX (kind=DP), INTENT (inout) :: epmat(nbndsub,nbndsub,nmodes)
  !! e-ph matrix elements 
  COMPLEX (kind=DP), INTENT (in) :: bmat(nbndsub,nbndsub) 
  !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
  !
  ! work variables
  !
  REAL(kind=DP) :: qeq,     &! <q+G| epsil | q+G>
       arg, zaq, g1, g2, g3, gmax, alph, geg
  INTEGER :: na, ipol, im, m1,m2,m3, nrx1,nrx2,nrx3, imode
  COMPLEX(kind=DP) :: fac, facqd, facq,  matsq
  COMPLEX(kind=DP) :: epmatl(nbndsub,nbndsub,nmodes)
  !
  IF (abs(signe) /= 1.0) &
       CALL errore ('rgd_blk',' wrong value for signe ',1)
  !
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  fac = signe*e2*fpi/omega * ci
  !
  epmatl(:,:,:) = czero   
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
            DO imode=1, nmodes
              CALL zaxpy(nbndsub**2,facq * zaq * uq(3*(na-1)+ipol,imode), bmat(:,:),1, epmat(:,:,imode),1)
              CALL zaxpy(nbndsub**2,facq * zaq * uq(3*(na-1)+ipol,imode), bmat(:,:),1, epmatl(:,:,imode),1)
            ENDDO
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
END SUBROUTINE rgd_blk_epw_fine
!
!-----------------------------------------------------------------------
SUBROUTINE compute_umn_f (nbnd, cufkk, cufkq, bmatf)
!-----------------------------------------------------------------------
  !!
  !! Calculates $$ U_{k+q} U_k^\dagger = <\Psi_{mk+q}|e^{i{q+G}r}|\Psi_{nk}> $$ in the approximation q+G->0
  !! on the fine grids.
  !!
  USE kinds,     ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  implicit none
  !
  !  input variables
  integer :: nbnd 
  ! number of bands (possibly in the optimal subspace)
  ! number of kpoint blocks, total
  complex(kind=DP) :: cufkk (nbnd, nbnd), cufkq (nbnd, nbnd)
  ! rotation matrix U(k)
  ! rotation matrix U(k+q)
  !
  !  output variables
  complex(kind=DP) :: bmatf (nbnd, nbnd)
  ! overlap wfcs in Bloch representation, fine grid
  !
  !  every pool works with its own subset of k points on the fine grid
  !
  bmatf = czero
  !
  !  U_q^ (k') U_k^dagger (k')
  !
  CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, cufkq, &
       nbnd, cufkk, nbnd, czero, bmatf, nbnd)
  !
  !bmatf = bmatf / dble(nkstot)
  !
  !
END SUBROUTINE compute_umn_f
!
!-----------------------------------------------------------------------
SUBROUTINE compute_umn_c (nbnd, nbndsub, nks, cuk, cukq, bmat)
!-----------------------------------------------------------------------
  !!
  !! Calculates $$ U_{k+q} U_k^\dagger = <\Psi_{mk+q}|e^{i(q+G)r}|\Psi_{nk}> $$ in the approximation q+G->0
  !! on the coarse grids. 
  !!
  USE kinds,     ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE io_global, ONLY : stdout
  implicit none
  !
  INTEGER, INTENT (in) :: nbnd
  !! number of bands
  INTEGER, INTENT (in) :: nks
  !! number of kpoint blocks, per pool
  INTEGER, INTENT (in) :: nbndsub
  !! Number of band on the subspace of Wannier
  ! 
  COMPLEX(kind=DP), INTENT (in) :: cuk (nbnd,nbndsub,nks)
  !! rotation matrix U(k)
  COMPLEX(kind=DP), INTENT (in) :: cukq (nbnd,nbndsub,nks)
  !! rotation matrix U(k+q)
  COMPLEX(kind=DP), INTENT (out) :: bmat (nbnd, nbnd, nks)
  !! overlap wfcs in Bloch representation, fine grid
  !
  !  work variables 
  integer :: ik
  !
  !  every pool works with its own subset of k points on the fine grid
  !
  bmat = czero
  !
  !  U_q^ (k') U_k^dagger (k')
  !
  DO ik=1,nks
     CALL zgemm ('n', 'c', nbnd, nbnd, nbndsub, cone, cukq(:,:,ik), &
          nbnd, cuk(:,:,ik), nbnd, czero, bmat(:,:,ik), nbnd)
  ENDDO
  !
  !WRITE(stdout,'(5x,a)') 'UMN calculated'
  !WRITE(stdout,*)
  !
  !
END SUBROUTINE compute_umn_c
!




