!
! Copyright (C) 2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3ionq (nat, ntyp, ityp, zv, tau, alat, omega, q, at, &
     bg, g, gg, ngm, gcutm, nmodes, u, ug0, npert_1, npert_f, q0mode, &
     d3dyn)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the contribution of the ions to the third order derivative
  ! of the total energy. Both the real and reciprocal space terms are included.
  !
  ! This version of the routine is general, i.e. it can compute D3^ewald(q1,q2,q3) with
  ! the only condition q1+q2+q3 = 0. Notice however, that only the case q1=q, q2=-q, q3=0
  ! has been extensively tested.
  !
  ! Written in February 2010 by L.Paulatto, T.Wassmann and  M.Lazzeri
  !
  ! The exact mechanism of this subroutine is quite complicated, a LaTeX form of all
  ! implemented formulas is reported here for reference and future extensions.
  ! Note that unit-of-measure dependent factors are missing (they can be derived from the code).
  !
! \begin{eqnarray*}
! atom1 & = & \{s_{1}(atom\_index),\tau_{s1}(position),Z_{s1}(charge)\}
! perturbation\_\nu_{1} & = & \{\alpha(cartensian\_direction),s_{1}(atom\_displaced)\}\end{eqnarray*}
! \begin{eqnarray*}
! D_{\nu1,\nu2,\nu3}^{3} & = & \delta_{s3,s1}Z_{s1}Z_{s2}F_{\alpha\beta\gamma}(q_{2},\tau_{s1}-\tau_{s2})
!  & + & \delta_{s1,s2}Z_{s2}Z_{s3}F_{\alpha\beta\gamma}(q_{3},\tau_{s2}-\tau_{s3})
!  & + & \delta_{s2,s3}Z_{s3}Z_{s1}F_{\alpha\beta\gamma}(q_{1},\tau_{s3}-\tau_{s1})
!  & - & \delta_{s1,s2,s3}Z_{s3}\sum_{s'}Z_{s'}F_{\alpha\beta\gamma}(0,\tau_{s3}-\tau_{s'})\end{eqnarray*}
! \begin{eqnarray*}
! F_{\alpha\beta\gamma}(q,\tau) & = & \frac{4\pi e^{2}}{\Omega}e^{i(G+q)\tau}
! \sum_{G}i(G+q)_{\alpha}(G+q)_{\beta}(G+q)_{\gamma}\frac{e^{-(G+q)^{2}/4\eta^{2}}}{(G+q)^{2}}
!  &  & -e^{2}\sum_{R}e^{iqR}\left.\frac{d^{3}f}{dx_{\alpha}dx_{\beta}dx_{\gamma}}\right|_{x=|\tau-R|}\end{eqnarray*}
! \begin{eqnarray*}
! \frac{d^{3}f(x)}{dx_{\alpha}dx_{\beta}dx_{\gamma}} & = &
! (\delta_{\alpha\beta}x_{\gamma}+\delta_{\alpha\gamma}x_{\beta}+\delta_{\beta\gamma}x_{\alpha})f_{1}(x)
!  &  & +x_{\alpha}x_{\beta}x_{\gamma}f_{3}(x)\end{eqnarray*}
! \begin{eqnarray*}
! f_{1}(x) &=& \frac{3erfc(\eta x)+a(\eta x)(3+2x^{2}\eta^{2})}{x^{5}}
! f_{3}(x) &=& -\frac{15erfc(\eta x)+a(\eta x)(15+10\eta^{2}x^{2}+4\eta^{4}x^{4})}{x^{7}}
! a(\xi) &=& \frac{2\xi}{\sqrt{\pi}}e^{-\xi^{2}}
! \end{eqnarray*}
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE constants, ONLY : e2, tpi, fpi, eps16, eps8
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  !  I/O variables
  INTEGER,INTENT(IN) :: nat, &        ! number of atoms
                        ntyp, &       ! number of types of atoms
                        ngm, &        ! number of G vectors
                        ityp (nat), & ! type of each atom
                        nmodes, &     ! number of modes
                        npert_1, &    ! only compute perturbations ...
                        npert_f       ! ... npert_1 < n < npert_f

  REAL (DP),INTENT(IN) :: tau (3, nat), & ! positions of the atoms
                          g (3, ngm), &   ! coordinates of g vectors
                          gg (ngm), &     ! modulus of g vectors
                          zv (ntyp), &    ! charge of each type
                          at (3, 3), &    ! direct lattice vectors
                          bg (3, 3), &    ! reciprocal lattice vectors
                          omega, &        ! volume of the unit cell
                          alat, &         ! length scale
                          gcutm, &        ! cut-off of g vectors
                          q (3)           ! q vector of perturbation -> D3(q,-q,0)

  COMPLEX (DP), INTENT(IN)  :: u (3*nat, nmodes), & ! pattern of the modes
                               ug0 (3*nat, nmodes)  ! pattern of the modes (q=0)

  COMPLEX (DP), INTENT(INOUT) :: d3dyn (3*nat, nmodes, 3*nat) ! derivative of the dyn. matrix

  LOGICAL, INTENT(IN) :: q0mode (300) ! if .true. this mode is to be computed
  ! Actually: all the modes between npert_1 and npert_f are always computed,
  ! but only the ones in q0mode are added to the dynamical matrix
  !
  !   Local variables
  !
  REAL(DP) :: q1(3),q2(3),q3(3) ! three q-vectors of the perturbations
  ! these will become INPUT parameters in future versions,
  ! at the moment it is always q1=q, q2=-q, q3=0

  REAL(DP),PARAMETER :: gamma(3) = (/ 0._dp, 0._dp, 0._dp /)
  INTEGER :: nu_1, nu_2, nu_3, &       ! perturbation indexes
             a_1, a_2, a_3, &          ! xyz indexes
             na_1, na_2, na_3, na_p,&  ! atom indexes
             nc_3cart,na_1cart,nb_2cart! additional indexes for changing to irrep. basis

  REAL(DP)::     alpha, eta, &         ! dumping factor of ewald sum, eta=sqrt(alpha)
                 upperbound, charge,  &! total charge in the cell
                 dtau(3)               ! aux: tau_s1 - tau_s2
  INTEGER :: abc(3)                    ! aux: {\alpha,\beta,\gamma}
  REAL (DP), EXTERNAL ::  qe_erfc
  COMPLEX (DP), ALLOCATABLE :: d3dion (:,:,:), d3dy2 (:,:,:) ! workspace
  COMPLEX (DP) :: work                                       ! more workspace
  !
  ! Undefine the following macros to esclude one of the terms
#define _D3_EWALD_G_SPACE
#define _D3_EWALD_REAL_SPACE
  !
  ! Temporary solution: this choice of q1,q2 and q3 reproduces the
  ! results of the previous code, minus a bug
  q1 = 0._dp
  q2 = q      ! GOOD FOR G-SPACE
  q3 = -q
  ! This alternative choice of q1,q2 and q3 reproduces the "wrong" value of the
  ! real-space term in the old code (only substantial for alpha < 1.0)
 !q1 =  q
 !q2 = -q     ! GOOD FOR R-SPACE
 !q3 =  0._dp
  !
  charge = SUM(zv(ityp(1:nat)))
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is an estimate of the error in the sum over G
  ! (empirical trust!)
  !
  upperbound = 1._dp
  alpha = 2.9_dp
  DO WHILE(upperbound > 1.e-9_dp)
     alpha = alpha - 0.1d0
     IF (alpha <= 0._dp) CALL errore ('d3ion', 'optimal alpha not found', 1)
     upperbound = 2 * charge**2 * SQRT(2 * alpha / tpi) &
                 * qe_erfc( SQRT((tpi/alat)**2 * gcutm / 4 / alpha) )
  ENDDO
  !
  eta = SQRT(alpha)
  WRITE( stdout, '(/5x,"Alpha used in Ewald sum = ",f6.2)') alpha
  !
  ALLOCATE  (d3dion( 3 * nat, nmodes, 3 * nat))
  d3dion (:,:,:) = (0.d0, 0.d0)
  !
  DO na_1 = 1,nat
  loop_a : &
  DO a_1  = 1,3
  nu_1 = a_1 + (na_1-1)*3
  !
  ! Inefficient but simple way to do only a subset of the perturbations
  ! (note: when nu_1 > npert_f BREAK would work as well)
  IF (nu_1 < npert_1 .or. nu_1 > npert_f) THEN
    CYCLE loop_a
  ENDIF
    !
    DO na_2 = 1,nat
    DO a_2  = 1,3
    nu_2 = a_2 + (na_2-1)*3
      !
      DO na_3 = 1,nat
      DO a_3  = 1,3
      nu_3 = a_3 + (na_3-1)*3
        !
        ! abc (read alpha-beta-gamma) is a list of the polarization
        ! for the three modes involved
        abc = (/ a_1,a_2,a_3 /)
        !
        ! delta_s1,s3
        IF (na_1==na_3) THEN
          dtau = tau(:,na_2) - tau(:,na_1)         ! tau_s2 - tau_s1
          work = zv(ityp(na_1)) * zv(ityp(na_2)) & ! z_s1 * z_s2
                * F_abc(q2,dtau,abc,eta)
          !
          d3dion(nu_1, nu_2, nu_3) = d3dion(nu_1, nu_2, nu_3) &
                                    + work
        ENDIF
        !
        ! delta_s1,s2
        IF (na_1==na_2) THEN
          dtau = tau(:,na_3) - tau(:,na_2)         ! tau_s3 - tau_s2
          work = zv(ityp(na_2)) * zv(ityp(na_3)) & ! z_s2 * z_s3
                * F_abc(q3,dtau,abc,eta)
          !
          d3dion(nu_1, nu_2, nu_3) = d3dion(nu_1, nu_2, nu_3) &
                                    + work
        ENDIF
        !
        ! delta_s2,s3
        IF (na_2==na_3) THEN
          dtau = tau(:,na_1) - tau(:,na_3)         ! tau_s1 - tau_s3
          work = zv(ityp(na_3)) * zv(ityp(na_1)) & ! z_s3 * z_s1
                * F_abc(q1,dtau,abc,eta)
          !
          d3dion(nu_1, nu_2, nu_3) = d3dion(nu_1, nu_2, nu_3) &
                                    + work
        ENDIF
        !
        ! delta_s1,s3,s3
        IF (na_1==na_2.and.na_2==na_3) THEN
          DO na_p = 1,nat
            dtau = tau(:,na_3) - tau(:,na_p)         ! tau_s3 - tau_sp
            work = zv(ityp(na_3)) * zv(ityp(na_p)) & ! z_s3 * z_sp
                  * F_abc(gamma,dtau,abc,eta)
            !
            d3dion(nu_1, nu_2, nu_3) = d3dion(nu_1, nu_2, nu_3) &
                                      + work
          ENDDO
        ENDIF
        !
      ENDDO !a_3
      ENDDO !na_3
      !
    ENDDO !a_2
    ENDDO !na_2
    !
  ENDDO loop_a !a_1
  ENDDO        !na_1
  !
#ifdef __MPI
  ! in the parallel case, recollect the modes
  CALL mp_sum( d3dion, intra_pool_comm )
  CALL mp_sum( d3dion, inter_pool_comm )
#endif
  !
  ! The dynamical matrix was computed in cartesian axis, now it is
  ! put on the basis of the modes; d3dy2 used as working array
  !
  ALLOCATE(d3dy2( 3*nat, nmodes, 3*nat))
  d3dy2 (:,:,:) = (0.d0, 0.d0)
  DO nu_3 = npert_1, npert_f
    !
    IF (q0mode (nu_3) ) THEN
      !
      DO nu_1 = 1, 3 * nat
      DO nu_2 = 1, 3 * nat
        !
        work = (0.d0, 0.d0)
        !
        DO nc_3cart = 1, 3 * nat
        DO na_1cart = 1, 3 * nat
        DO nb_2cart = 1, 3 * nat
          work = work + ug0 (nc_3cart, nu_3) &
                * CONJG(u (na_1cart, nu_1) ) &
                * d3dion (nc_3cart, na_1cart, nb_2cart) &
                * u (nb_2cart, nu_2)
        ENDDO
        ENDDO
        ENDDO
        !
        d3dy2 (nu_3, nu_1, nu_2) = work
        !
      ENDDO
      ENDDO
      !
    ENDIF
    !
  ENDDO
  !
#ifdef __MPI
  CALL mp_sum ( d3dy2, inter_pool_comm )
#endif
  !
  ! For debugging purposes (to be removed), the Ewald contribution
  ! can be dumped to file (uncomment the lines that apply).
  ! 1. using internal debugging subroutine
  !   CALL writed3dyn_5(d3dy2,'d3qewald',-1)
  ! 2. using iotk
  !  CALL iotk_write_dat(1077, 'd3ionq', d3dy2)
  ! 3. by hand, the old way
!     open(unit=1077, file='d3ionq-n.xml', action='write', status='unknown')
!     do a_1 = 1,3*nat
!     do a_2 = 1,3*nat
!     do a_3 = 1,3*nat
!       write(1077, '(3i4,2f32.16)') a_1, a_2, a_3, d3dy2(a_1,a_2,a_3)
!     enddo
!     enddo
!     enddo
!     close(1077)
  !
  ! Add the Ewald term to the rest of D3 matrix
  d3dyn = d3dyn+d3dy2
  !
  DEALLOCATE (d3dion, d3dy2)
  !
  RETURN

!-----------------------------------------------------------------------
CONTAINS
    !-------------------------------------------------------------------
    !
    ! dumping factor of Ewald sum
    ! 2/sqrt(pi) eta*x exp(-eta**2 x**2)
    !-----------------------------------------------------------------------
    FUNCTION a_fct(xeta)
        !-------------------------------------------------------------------
        USE constants, ONLY : sqrtpm1 ! 1/sqrt(pi)
        IMPLICIT NONE
        REAL(DP) :: a_fct
        REAL(DP),INTENT(IN) :: xeta
        a_fct = 2*sqrtpm1*xeta*exp(-(xeta)**2)
        ! note: 2*sqrtpm1 == 2/sqrt(pi) == sqrt (8.d0 / tpi) <- from old code
    END FUNCTION
    !
    ! Used by d3f_abc, it's (related to) the second derivative of erfc function
    ! f1
    !-----------------------------------------------------------------------
    FUNCTION d2f_fct(xx, eta)
        !-------------------------------------------------------------------
        IMPLICIT NONE
        REAL(DP) :: d2f_fct
        REAL(DP),INTENT(IN) :: xx, eta
        REAL(DP)            :: xeta
        REAL(DP), EXTERNAL  :: qe_erfc
        xeta = xx*eta
        !
        d2f_fct = 3._dp*qe_erfc(xeta) + a_fct(xeta)*(3._dp + 2*(xeta**2))
        d2f_fct = d2f_fct/xx**5
    END FUNCTION
    !
    ! Used by d3f_abc, it's (related to) the third derivative of erfc function
    ! f3
    !-----------------------------------------------------------------------
    FUNCTION d3f_fct(xx, eta)
        !-------------------------------------------------------------------
        IMPLICIT NONE
        REAL(DP) :: d3f_fct
        REAL(DP),INTENT(IN) :: xx, eta
        REAL(DP)            :: xeta, xeta2
        REAL(DP), EXTERNAL  :: qe_erfc
        xeta  = xx*eta
        xeta2 = xeta**2
        d3f_fct = 15._dp*qe_erfc(xeta) &
                 + a_fct(xeta)*(15._dp + 10._dp*xeta2 + 4*(xeta2**2))
        d3f_fct = -d3f_fct/xx**7
    END FUNCTION
    !
    ! Used for real-space term
    ! d3f(x)/dx_a dx_b dx_c
    !-----------------------------------------------------------------------
    FUNCTION d3f_abc(x, xx, abc, eta)
        !-------------------------------------------------------------------
        IMPLICIT NONE
        REAL(DP) :: d3f_abc
        REAL(DP),INTENT(IN) :: x(3), xx, eta
        INTEGER,INTENT(IN)  :: abc(3)
        !
        REAL(DP) :: delta3 ! delta_{a,b} x_c + delta_{a,c} x_b + delta_{b,c} x_a
        REAL(DP) :: xa_xb_xc ! x_a * x_b * x_c
        !
        d3f_abc=0._dp
        !
        !
        delta3 = 0._dp
        IF(abc(1)==abc(2)) delta3 = delta3 + x(abc(3))
        IF(abc(2)==abc(3)) delta3 = delta3 + x(abc(1))
        IF(abc(3)==abc(1)) delta3 = delta3 + x(abc(2))
        delta3 = delta3*alat
        !
        IF( ABS(delta3) > eps16) THEN
            d3f_abc = d3f_abc + delta3*d2f_fct(xx, eta)
        ENDIF
        !
        !
        xa_xb_xc = x(abc(1))*x(abc(2))*x(abc(3))*alat**3
        !
        IF( ABS(xa_xb_xc) > eps16) THEN
            d3f_abc = d3f_abc + xa_xb_xc*d3f_fct(xx, eta)
        ENDIF
        !
    END FUNCTION
    !
    !
    !-----------------------------------------------------------------------
    FUNCTION F_abc(q,tau,abc,eta)
        !-------------------------------------------------------------------
        USE constants, ONLY : tpi, fpi, e2, eps8
        USE mp_global, ONLY : nproc_image, me_image, intra_image_comm
        IMPLICIT NONE
        COMPLEX(DP)         :: F_abc
        REAL(DP),INTENT(IN) :: q(3), tau(3), eta
        INTEGER, INTENT(IN) :: abc(3)
        COMPLEX(DP),PARAMETER :: ii   = (0._dp, 1._dp), &
                                 zero = (0._dp, 0._dp), &
                                 one  = (1._dp, 0._dp)
        !
        REAL(DP) :: prefG, facq  ! prefactors for G-space term
        REAL(DP) :: Gpq_abc
        REAL(DP) :: Gpq_tau
        INTEGER  :: ng
        !
        INTEGER,PARAMETER :: mxr = 100     ! max number of neighbours
        REAL(DP)    :: r (3,mxr), r2 (mxr) ! shells of neighbours (r and r**2)
        REAL(DP)    :: rr                  ! sqrt(r2)*alat
        REAL(DP)    :: rmax                ! radius containg the shells of ngbrs
        INTEGER     :: nrm, nr             ! number of neighbours in teh shell, and their index
        INTEGER     :: nr_s, nr_e, mykey   ! used to parallelize r-space sum
        COMPLEX(DP) :: facr
        REAL(DP)    :: qdr  ! q*g
        REAL(DP)    :: gtq2 ! (g+q)**2 (atomic units)
        !
        ! First part: the reciprocal space term
        !
        F_abc  = zero
        prefG = fpi * e2 * (tpi/alat)**3 / omega
        !
#ifdef _D3_EWALD_G_SPACE
        !
        sum_on_G : &
        DO ng = 1, ngm
            !
            Gpq_abc =  ( g(abc(1), ng) + q(abc(1)) ) &
                      * ( g(abc(2), ng) + q(abc(2)) ) &
                      * ( g(abc(3), ng) + q(abc(3)) )
            !
            ! Skip null terms
            IF (ABS(Gpq_abc) < eps8) &
                CYCLE sum_on_G
            !
            gtq2 = (  (g(1, ng) + q(1)) **2 &
                    + (g(2, ng) + q(2)) **2 &
                    + (g(3, ng) + q(3)) **2 ) * (tpi/alat) **2
            !
            facq = Gpq_abc * prefG * EXP( - gtq2 / eta**2 / 4._dp) / gtq2
            !
            Gpq_tau = tpi *(  ( g(1, ng) + q(1) ) * tau(1) &
                            + ( g(2, ng) + q(2) ) * tau(2) &
                            + ( g(3, ng) + q(3) ) * tau(3) )
            !
            F_abc = F_abc - ii*facq* EXP(ii*Gpq_tau)
            !
        ENDDO sum_on_G
        !
#endif
!         print*, "  nrm",nrm
#ifdef _D3_EWALD_REAL_SPACE
        !
        ! Second part: the real space term
        !
        rmax   = 5.d0 / eta / alat
        CALL rgen (tau, rmax, mxr, at, bg, r, r2, nrm)
        ! note: r = R - tau : R is a real-space cell vector
        !
        ! In some cases the real-space term does not include any term
        IF( nrm>0 ) THEN
          !
          ! Parallelize the real space sum, it will hardly give any performance
          ! improvement, but cannot hurt (alternatively this term must be computed
          ! by one processor only, i.e. ionode)
          CALL block_distribute( nrm, me_image, nproc_image, nr_s, nr_e, mykey )
          !
          ! If we have more CPUs than nrm some will do nothing
          IF(mykey==0)THEN
            sum_on_R : &
            DO nr = nr_s, nr_e
                rr   = SQRT(r2(nr)) * alat
                qdr = tpi * (  q (1) * (r(1, nr) + tau (1)) &
                            + q (2) * (r(2, nr) + tau (2)) &
                            + q (3) * (r(3, nr) + tau (3)) )
                !
                IF (ABS(qdr) < eps16) THEN
                    facr = - e2*one
                ELSE
                    facr = - e2*CMPLX(cos(qdr), sin(qdr), kind=DP)
                ENDIF
                !
                F_abc = F_abc + facr*d3f_abc(r(1:3,nr),rr,abc,eta)
                !
            ENDDO sum_on_R
          ENDIF
          !
        ENDIF
        !
#endif
        !
        RETURN
        !
    END FUNCTION F_abc


END SUBROUTINE d3ionq



