MODULE coulomb
  !
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER       :: nq1=1, nq2=1, nq3=1
  REAL(DP)      :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP)      :: eps  = 1.d-6
  REAL(DP)      :: exxdiv     = 0._dp
  REAL(DP)      :: exxdiv_eps = 0._dp
  CHARACTER(32) :: exxdiv_treatment  = ' '
  LOGICAL       :: use_regularization = .FALSE.
  ! ... x_gamma_extrapolation
  LOGICAL       :: x_gamma_extrapolation =.FALSE.
  LOGICAl       :: on_double_grid =.FALSE.
  REAL(DP)      :: grid_factor = 1.d0 !8.d0/7.d0
  !
  REAL(DP)      :: yukawa = 0._dp
  REAL(DP)      :: erfc_scrlen = 0._dp
  REAL(DP)      :: erf_scrlen = 0._dp
  REAL(DP)      :: gau_scrlen = 0.d0
  REAL(DP)      :: eps_mat(3,3)
  !
  CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE setup_coulomb ( )
    !--------------------------------------------------------------------
    !
    !	
    !!  This routine initialize the coulomb kernel according to the 
    !!  Gygi-Balderschi scheme. This is compatible with peridic 
    !!  systems and finite sytems calculation. For finite system 
    !!  one can also (should) use assume_isolated (see kcw_readin) 
    !!  Adapted from exx_divergence () inside PW. Generalized to the
    !!  screened Coulomb eps^-1 * vc adpting the madelung constant 
    !!  formula in Rurali and Cartoixa Nano Letters 9. 975 (2009)
    !
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : mp1, mp2, mp3, l_vcut, eps_inf, calculation
    USE martyna_tuckerman,    ONLY : do_comp_mt
    USE cell_base,            ONLY : omega
    !
    IMPLICIT NONE 
    !
    LOGICAL :: exst, skip_eps
    INTEGER :: nqs
    !
    CALL start_clock( 'Coulomb setup' )
    !
    nq1=mp1; nq2=mp2; nq3=mp3
    nqs=nq1*nq2*nq3
    !
    eps_mat = 0.D0
    eps_mat (1,1) = 1.D0; eps_mat (2,2) = 1.D0; eps_mat (3,3) = 1.D0
    !
    yukawa = 0._dp
    erfc_scrlen = 0._dp
    erf_scrlen = 0._dp
    gau_scrlen = 0.d0
    exxdiv_treatment = 'none'
    x_gamma_extrapolation = .false.
    ! NOTABENE: The Screened contribution does not work with x_gamma_extrapolation
    !           Leave x_gamma_extrapolation to FALSE unless you know what you are doing
    !
    skip_eps = .FALSE. 
    IF (l_vcut) THEN
      !
      use_regularization = .true. 
      exxdiv_treatment = 'gb' ! For now only GB scheme
      !
      ! Set the dielectric tensor. Priority is given in this order
      !   1) if a file eps.dat exist, read and set up the dielectric tensor
      INQUIRE( file="eps.dat", exist=exst )
      IF ( exst) THEN  
         !
         WRITE (stdout,'(/,5X, "INFO: Dielectric tensor read from file eps.dat")')  
         CALL read_eps (eps_mat)
         eps_inf = (eps_mat(1,1)+eps_mat(2,2)+eps_mat(3,3))/3.D0
         WRITE (stdout,'(  5X, "INFO: average macroscopic eps", 1F12.6)')  eps_inf
         !
         ! 2) the average epsilon_inf is given in input
      ELSE IF ( ABS(eps_inf -1.D0) .gt. 1E-6) THEN
         WRITE (stdout,'(/,5X, "INFO: average macroscopic eps from input", 1F12.6)')  eps_inf
         eps_mat (1,1) = eps_inf; eps_mat (2,2) = eps_inf; eps_mat (3,3) = eps_inf
         !
         ! 3) No information about the dielectric properties: set eps=I
      ELSE 
         WRITE (stdout,'(/,5X, "INFO: Dielectric tensor = NOT SPECIFIED")')
         WRITE (stdout,'(  5X, "      NO Correction for the Screened Coulomb")')
         skip_eps = .TRUE.
         !
      ENDIF
      !
    ENDIF
    !
    CALL divergence ()
    IF (skip_eps) exxdiv_eps = 0.D0
    ! Compute the divergence and set the q+G=0 term (exxdiv and exxdiv_eps, 
    ! to be used to correct the bare and screened coulomb. USed in bare_pot
    ! and screen_coeff. 
    !
    IF (.NOT. do_comp_mt) THEN
      !
      WRITE (stdout,'(/,5X, "INFO: Coulomb q+G=0 treatment:")') 
      WRITE (stdout,'(  5X, "INFO: Divergence         ", 3x, 1A8)'     )  exxdiv_treatment
      WRITE( stdout,'(  5X, "INFO: q-grid dimension   ", 3x, 3I4)'     )  nq1, nq2, nq3
      WRITE( stdout,'(  5X, "INFO: cell volume        ", 3x, 1F20.12)' )  omega
      WRITE (stdout,'(  5X, "INFO: Gamma Extrapolation", 3x, 1L5 )'    )  x_gamma_extrapolation
      !
      IF ( x_gamma_extrapolation ) THEN
         !
         WRITE( stdout, '(5X, "INFO: extrapolation q->0 dealt with 8/7 -1/7 trick ")')
         grid_factor = 8.d0 / 7.d0
         !
      ELSE
         !
         WRITE( stdout, '(5X, "INFO: extrapolation q->0 term not estimated")' )
         grid_factor = 1.d0
         !
      ENDIF
      !
      WRITE (stdout,'(  5X,    "INFO: Bare Coulomb q+G=0     ", 3x, 1ES15.5 )')  exxdiv
      IF ( (calculation == "screen" .AND. l_vcut) .OR. calculation == 'cc' ) THEN
         !
         IF (.NOT. skip_eps) THEN 
           WRITE (stdout,'(  5X, "INFO: Epsilon infinity       ", 3x, 1F12.6  )')  eps_inf
           WRITE (stdout,'(  5X, "INFO: Dielectric tensor      ", 3x, 3F12.6  )') (eps_mat(:,1))
           WRITE (stdout,'(  5X, "                             ", 3x, 3F12.6  )') (eps_mat(:,2))
           WRITE (stdout,'(  5X, "                             ", 3x, 3F12.6  )') (eps_mat(:,3))
         ENDIF
         WRITE (stdout,'(  5X, "INFO: Screened Coulomb q+G=0 ", 3x, 1F20.12 )')  exxdiv_eps
         IF (.NOT. skip_eps) THEN 
            WRITE (stdout,'(  5X, "      Isotropic estimate     ", 3x, 1F20.12 )')  exxdiv/eps_inf
         ELSE
            WRITE (stdout,'(  5X, "      Isotropic estimate     ", 3x, 1F20.12 )')  0.D0
         ENDIF
         WRITE (stdout,'(  5X, "INFO: uPi(q+G=0) estimation  ", 3X, 1F20.12 )') -exxdiv/omega/nqs
         WRITE (stdout,'(  5X, "INFO: rPi(q+G=0) estimation  ", 3X, 1F20.12 )') -(exxdiv_eps)/omega/nqs
         !
      ENDIF
      !
    ELSE 
      !
      WRITE (stdout,'(/,5X, "INFO: Coulomb kernel treated according to MT scheme")') 
      !
    ENDIF 
    !
    CALL stop_clock( 'Coulomb setup' )
    !  
    RETURN
    !
  END SUBROUTINE setup_coulomb
  !
  !----------------------------------------------------------------------
  SUBROUTINE g2_convolution( ngm, g, xk, xkq, fac )
    !--------------------------------------------------------------------
    !! This routine calculates the 1/|r-r'| part of the exact exchange
    !! expression in reciprocal space (the G^-2 factor).
    !! It then regularizes it according to the specified recipe.
    !
    USE kinds,       ONLY : DP
    USE cell_base,   ONLY : at, tpiba2
    USE constants,   ONLY : fpi, e2, pi
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: ngm
    !! Number of G vectors
    REAL(DP), INTENT(IN) :: g(3,ngm)
    !! Cartesian components of G vectors
    REAL(DP), INTENT(IN) :: xk(3)
    !! current k vector
    REAL(DP), INTENT(IN) :: xkq(3)
    !! current q vector
    REAL(DP), INTENT(INOUT) :: fac(ngm)
    !! Calculated convolution
    !
    ! ... local variables
    !
    INTEGER :: ig !Counters
    REAL(DP) :: q(3), qq, x
    REAL(DP) :: grid_factor_track(ngm), qq_track(ngm)
    REAL(DP) :: nqhalf_dble(3)
    LOGICAL :: odg(3)
    !
    ! Now the Coulomb potential that are computed on the fly
    !
    nqhalf_dble(1:3) = (/ DBLE(nq1)*0.5_DP, DBLE(nq2)*0.5_DP, DBLE(nq3)*0.5_DP /)
    !
    ! Set the grid_factor_track and qq_track
    !
    IF ( x_gamma_extrapolation ) THEN
!$omp parallel do default(shared), private(ig,q,x,odg)
       DO ig = 1, ngm
          q(:)= xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = SUM(q(:)**2) * tpiba2
          x = (q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nqhalf_dble(1)
          odg(1) = ABS(x-NINT(x)) < eps
          x = (q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nqhalf_dble(2)
          odg(2) = ABS(x-NINT(x)) < eps
          x = (q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nqhalf_dble(3)
          odg(3) = ABS(x-NINT(x)) < eps
          IF( ALL( odg(:) ) ) THEN
             grid_factor_track(ig) = 0._DP ! on double grid
          ELSE
             grid_factor_track(ig) = grid_factor ! not on double grid
          ENDIF
       ENDDO
!$omp end parallel do
    ELSE
!$omp parallel do default(shared), private(ig,q)
       DO ig = 1, ngm
          q(:) = xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = SUM(q(:)**2) * tpiba2
       ENDDO
       grid_factor_track = 1._DP
       !
    ENDIF
    ! ... The big loop
    !
!$omp parallel do default(shared), private(ig,qq)
    DO ig = 1, ngm
      !
      qq = qq_track(ig)
      !
      IF(gau_scrlen > 0) THEN
         fac(ig) = e2*((pi/gau_scrlen)**(1.5_DP))*EXP(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
         !
      ELSEIF (qq > eps_qdiv) THEN
         !
         IF ( erfc_scrlen > 0  ) THEN
            fac(ig) = e2*fpi/qq*(1._DP-EXP(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
         ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig) = e2*fpi/qq*(EXP(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
         ELSE
            fac(ig) = e2*fpi/( qq + yukawa ) * grid_factor_track(ig) ! as HARTREE
         ENDIF
         !
      ELSE
         !
         fac(ig) = - exxdiv ! or rather something ELSE (see F.Gygi)
         !
         IF (yukawa>0._DP .AND. .NOT.x_gamma_extrapolation) fac(ig) = fac(ig) + &
                                                            e2*fpi/( qq + yukawa )
         !
         IF (erfc_scrlen>0._DP .AND. .NOT.x_gamma_extrapolation) fac(ig) = fac(ig) + &
                                                                  e2*pi/(erfc_scrlen**2)
         !
      ENDIF
      !
    ENDDO
!$omp end parallel do
    !
  END SUBROUTINE g2_convolution
  !
  !-----------------------------------------------------------------------
  SUBROUTINE divergence()
      !-----------------------------------------------------------------------
    !
    ! ...  This function calculates the G=0 term of the Coulomb potential
    ! ...  in presence  of a unitary charge (q=1). Used for counter-charge
    ! ...  corrections.
    !
    USE kinds,              ONLY : DP
    USE constants,          ONLY : fpi, e2
    USE cell_base,          ONLY : bg, at, alat, omega
    USE gvect,              ONLY : ngm, g
    USE gvecw,              ONLY : gcutw
    USE control_flags,      ONLY : gamma_only
    USE mp_global,          ONLY : intra_pool_comm
    USE mp,                 ONLY : mp_sum
    USE control_kcw,        ONLY : eps_inf
    !
    !
    IMPLICIT NONE
    !
    INTEGER  :: nqs                   ! number of points in the q-grid
    LOGICAL  :: on_double_grid=.false.
    INTEGER  :: iq1,iq2,iq3, ig
    REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)
    INTEGER  :: nqq, iq
    REAL(DP) :: aa, dq
    !
    REAL(DP) :: q_eps_q
    REAL(DP) :: aa_eps, div_eps, det_eps, alpha_eps, dq_eps
    !
    CALL start_clock( 'exx_div' )
    !
    IF ( .NOT. use_regularization ) THEN
       exxdiv     = 0._DP
       exxdiv_eps = 0._DP
       RETURN
    ENDIF
    !
    nqs = nq1 * nq2 * nq3
    tpiba2 = ( fpi / 2.d0 / alat ) ** 2
    alpha  = 10.d0 / gcutw
    alpha_eps = alpha/eps_inf
    !
    dq1= 1.d0/DBLE(nq1)
    dq2= 1.d0/DBLE(nq2) 
    dq3= 1.d0/DBLE(nq3)
    !
    div = 0.d0
    div_eps = 0.d0
    !
    DO iq1 = 1, nq1
      DO iq2 = 1, nq2
        DO iq3 = 1, nq3
          !
          xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                  bg(:,2) * (iq2-1) * dq2 + &
                  bg(:,3) * (iq3-1) * dq3
          !
          DO ig = 1, ngm
            !
            q(1) = xq(1) + g(1,ig)
            q(2) = xq(2) + g(2,ig)
            q(3) = xq(3) + g(3,ig)
            qq = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
            q_eps_q = DOT_PRODUCT( q(:), MATMUL( eps_mat(:,:), q(:) ) )
            !
            IF ( x_gamma_extrapolation ) THEN
              !
              on_double_grid = .true.
              x = 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
              on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
              x = 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
              on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
              x = 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
              on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
              !
            ENDIF
            !
            IF ( .NOT. on_double_grid ) THEN
              IF ( qq > 1.d-8 ) THEN
                div = div + EXP( - alpha * qq ) / ( qq + yukawa / tpiba2 ) &
                                                    * grid_factor
                div_eps = div_eps + EXP( - alpha_eps * q_eps_q ) / ( q_eps_q + yukawa / tpiba2 ) &
                                                    * grid_factor
              ENDIF
            ENDIF
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum( div, intra_pool_comm )
    CALL mp_sum( div_eps, intra_pool_comm )
    !
    IF ( gamma_only ) div = 2.d0 * div
    IF ( gamma_only ) div_eps = 2.d0 * div_eps
    !
    IF ( .NOT. x_gamma_extrapolation ) THEN
      IF ( yukawa < 1.d-8) THEN
        div = div - alpha
        div_eps = div_eps - alpha_eps
      ELSE
        div = div + tpiba2 / yukawa
        div_eps = div_eps + tpiba2 / yukawa
      ENDIF
    ENDIF
    !
    div = div * e2 * fpi / tpiba2 / nqs
    div_eps = div_eps * e2 * fpi / tpiba2 / nqs
    !
    alpha = alpha / tpiba2
    alpha_eps = alpha_eps / tpiba2
    !
    nqq = 100000
    dq = 5.0d0 / SQRT( alpha ) / nqq
    dq_eps = 5.0d0 / SQRT( alpha_eps ) / nqq
    aa = 0.d0
    aa_eps = 0.d0
    !
    DO iq = 0, nqq
      !
      q_ = dq * ( iq + 0.5d0 )
      qq = q_ * q_
      aa = aa - EXP( - alpha * qq ) * yukawa / ( qq + yukawa ) * dq
      q_ = dq_eps * ( iq + 0.5d0 )
      qq = q_ * q_
      aa_eps = aa_eps - EXP( - alpha_eps * q_eps_q ) * yukawa / ( q_eps_q + yukawa ) * dq_eps
      !
    ENDDO
    !
    det_eps = eps_mat(1,1)*(eps_mat(2,2)*eps_mat(3,3)-eps_mat(2,3)*eps_mat(3,2)) & 
             -eps_mat(2,2)*(eps_mat(2,1)*eps_mat(3,3)-eps_mat(2,3)*eps_mat(3,1)) &
             +eps_mat(3,3)*(eps_mat(2,1)*eps_mat(3,2)-eps_mat(2,2)*eps_mat(3,1)) 
    !
    aa = aa * 8.d0 / fpi
    aa_eps = aa_eps * 8.d0 / fpi
    aa = aa + 1.d0 / SQRT( alpha * 0.25d0 * fpi )
    aa_eps = aa_eps + 1.d0 / SQRT( alpha_eps * 0.25d0 * fpi * det_eps)
    !  
    div = div - e2 * omega * aa
    div_eps = div_eps - e2 * omega * aa_eps
    !
    exxdiv = div * nqs
    exxdiv_eps = div_eps * nqs
    !
    CALL stop_clock( 'exx_div' )
    !
    RETURN
    !
  END SUBROUTINE divergence 
  !
  !---------------------------------------------------------------------
  SUBROUTINE read_eps (eps_mat) 
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL (DP) :: eps_mat (3,3)
    INTEGER   :: i, j, ios
    !
    OPEN (765, file = 'eps.dat')
    DO j = 1, 3
      READ (765, *, IOSTAT=ios) (eps_mat(j,i), i =1,3)
      IF (ios /= 0) THEN 
       WRITE(stdout,'(/, 5X,A, I5)') "ERROR: Somethng wrong reading eps.dat", ios
       STOP
      ENDIF
    ENDDO
    CLOSE (765)
    !
    RETURN
    !
  END SUBROUTINE read_eps
  !
END MODULE
