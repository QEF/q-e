  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2016-2018 Samuel Ponce'
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE transport_iter
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutine linked with self-consistent electronic transport  
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE iterativebte( iter, iq, ef0, efcb, error_h, error_el, first_cycle, first_time ) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine computes the scattering rate with the iterative BTE
    !!  (inv_tau).
    !!  The fine k-point and q-point grid have to be commensurate. 
    !!  The k-point grid uses crystal symmetry to decrease computational cost.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : alat, at, omega, bg
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : fsthick, & 
                              eps_acustic, degaussw, nstemp, & 
                              system_2d, int_mob, ncarrier, restart, restart_freq,&
                              mp_mesh_k, nkf1, nkf2, nkf3, vme
    USE pwcom,         ONLY : ef 
    USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, wkf, dmef, vmef, & 
                              wf, wqf, xkf, epf17, nqtotf, nkqtotf, inv_tau_all, xqf, & 
                              F_current, Fi_all, F_SERTA, F_currentcb, Fi_allcb, &
                              F_SERTAcb, inv_tau_allcb
    USE transportcom,  ONLY : transp_temp, mobilityh_save, mobilityel_save, lower_bnd, &
                              ixkqf_tr, s_BZtoIBZ_full
    USE constants_epw, ONLY : zero, one, two, pi, kelvin2eV, ryd2ev, & 
                              electron_SI, bohr2ang, ang2cm, hbarJ, eps6
    USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, world_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE symm_base,     ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE superconductivity, ONLY : kpmq_map
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT (INOUT) :: first_time
    LOGICAL, INTENT (INOUT) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(IN) :: iter
    !! Iteration number
    INTEGER, INTENT(IN) :: iq
    !! Q-point index
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp. Could be unused (0).
    REAL(KIND=DP), INTENT(out) :: error_h
    !! Error on the hole mobility made in the last iterative step.
    REAL(KIND=DP), INTENT(out) :: error_el
    !! Error on the electron mobility made in the last iterative step.
    !
    ! Local variables
    INTEGER :: i, iiq
    !! Cartesian direction index 
    INTEGER :: j
    !! Cartesian direction index 
    INTEGER :: ij
    !! Cartesian direction index 
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: imode
    !! Local mode index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ipool
    !! Index of the pool
    INTEGER :: nkq
    !! Index of the pool the the k+q point is
    INTEGER :: nkq_abs
    !! Index of the k+q point from the full grid. 
    INTEGER :: BZtoIBZ(nkf1*nkf2*nkf3)
    !! Map between the full uniform k-grid and the IBZ
    INTEGER :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
    !! Save the symmetry operation that brings BZ k into IBZ
    INTEGER :: nkqtotf_tmp
    ! 
    REAL(KIND=DP) :: tau
    !! Relaxation time
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND=DP) :: ekq
    !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
    REAL(KIND=DP) :: g2
    !! Electron-phonon matrix elements squared (g2 is Ry^2) 
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: w0g1
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q}] $$ 
    REAL(KIND=DP) :: w0g2 
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q}] $$
    REAL(KIND=DP) :: inv_wq 
    !! Inverse phonon frequency. Defined for efficiency reasons.
    REAL(KIND=DP) :: inv_etemp
    !! Invese temperature inv_etemp = 1/etemp. Defined for efficiency reasons.
    REAL(KIND=DP) :: g2_tmp 
    !! Used to set component to 0 if the phonon freq. is too low. This is defined
    !! for efficiency reasons as if statement should be avoided in inner-most loops.
    REAL(KIND=DP) :: inv_degaussw
    !! 1.0/degaussw. Defined for efficiency reasons. 
    REAL(KIND=DP) :: wq
    !! Phonon frequency $$\omega_{q\nu}$$ on the fine grid.  
    REAL(KIND=DP) :: wgq
    !! Bose-Einstein occupation function $$n_{q\nu}$$
    REAL(KIND=DP) :: fmkq
    !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
    REAL(KIND=DP) :: trans_prob
    !! Transition probability function
    REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
    !! Electronic velocity $$v_{n\mathbf{k}}$$
    REAL(KIND=DP) :: dfnk
    !! Derivative Fermi distribution $$-df_{nk}/dE_{nk}$$
    REAL(KIND=DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND=DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND=DP) :: mobility
    !! Sum of the diagonalized mobilities [cm^2/Vs] 
    REAL(KIND=DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND=DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND=DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs]
    REAL(KIND=DP) :: tdf_sigma(9)
    !! Transport distribution function
    REAL(KIND=DP) :: Sigma(9,nstemp)
    !! Electrical conductivity
    REAL(KIND=DP) :: sigma_up(3,3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND=DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND=DP) :: sigma_vect(3,3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND=DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(kind=DP) :: xkf_all(3,nkqtotf)
    !! Collect k-point coordinate (and k+q) from all pools in parallel case
    REAL(kind=DP) :: xkf_red(3,nkqtotf/2)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(kind=DP) :: xxq(3)
    !! Current q-point 
    REAL(kind=DP) :: xkk(3)
    !! Current k-point on the fine grid
    REAL(kind=DP) :: Fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    !
    !
    REAL(KIND=DP), EXTERNAL :: DDOT
    !! Dot product function
    REAL(KIND=DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function  
    REAL(kind=DP) :: xkf_tmp (3, nkqtotf)
    !! Temporary k-point coordinate (dummy variable)
    REAL(kind=DP) :: wkf_tmp(nkqtotf)
    !! Temporary k-weights (dummy variable)
    ! 
    inv_cell = 1.0d0/omega
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF ( system_2d ) &
       inv_cell = inv_cell * at(3,3) * alat
    !
    ! 
    ! Gather all the k-point coordinate from all the pools
    xkf_all(:,:) = zero 
    xkf_red(:,:) = zero 
    ! 
#ifdef __MPI
    ! 
    CALL poolgather2 ( 3, nkqtotf, nkqf, xkf, xkf_all) 
#else
    !
    xkf_all = xkf
    !
#endif 
    ! 
    IF (first_time) THEN
      first_time = .FALSE.
      DO itemp = 1, nstemp
        !   
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          ! 
          IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
            DO ibnd = 1, ibndmax-ibndmin+1
              !
              ! vkk(3,nbnd) - velocity for k
              IF ( vme ) THEN
                ! vmef is in units of Ryd * bohr
                vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
              ELSE
                ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                ! 1/m  = 2 in Rydberg atomic units
                ! dmef is in units of 1/a.u. (where a.u. is bohr)
                ! v_(k,i) is in units of Ryd * a.u.
                vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
              ENDIF
              ! 
              ! The inverse of SERTA 
              tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
              F_SERTA(:,ibnd,ik+lower_bnd-1,itemp) = vkk(:,ibnd) * tau
              !
              ! In this case we are also computing the scattering rate for another Fermi level position
              ! This is used to compute both the electron and hole mobility at the same time.  
              IF ( ABS(efcb(itemp)) > eps6 ) THEN
                tau = one / inv_tau_allcb(1,ibnd,ik+lower_bnd-1)
                F_SERTAcb(:,ibnd,ik+lower_bnd-1,itemp) = vkk(:,ibnd) * tau
              ENDIF
              !
            ENDDO ! ibnd
          ENDIF
        ENDDO ! ik 
      ENDDO ! itemp
      ! 
      Fi_all = F_SERTA
      CALL mp_sum( Fi_all, world_comm )
      IF ( ABS(efcb(1)) > eps6 ) THEN
        Fi_allcb = F_SERTAcb 
        CALL mp_sum( Fi_allcb, world_comm )
      ENDIF
      ! 
      IF (mp_mesh_k) THEN
        IF ( .not. ALLOCATED(ixkqf_tr) ) ALLOCATE(ixkqf_tr(nkf,nqtotf))
        IF ( .not. ALLOCATED(s_BZtoIBZ_full) ) ALLOCATE(s_BZtoIBZ_full(3,3,nkf,nqtotf))
        ixkqf_tr(:,:) = 0
        s_BZtoIBZ_full(:,:,:,:) = 0
        ! 
        IF ( mpime .eq. ionode_id ) THEN
          ! 
          CALL set_sym_bl( )
          !
          BZtoIBZ(:) = 0
          s_BZtoIBZ(:,:,:) = 0 
          ! What we get from this call is BZtoIBZ
          CALL kpoint_grid_epw ( nrot, time_reversal, .false., s, t_rev, bg, nkf1*nkf2*nkf3, &
                     nkf1,nkf2,nkf3, nkqtotf_tmp, xkf_tmp, wkf_tmp,BZtoIBZ,s_BZtoIBZ)
          ! 
          DO ik = 1, nkqtotf/2
            ikk = 2 * ik - 1
            xkf_red(:,ik) = xkf_all(:,ikk)
          ENDDO 
          ! 
        ENDIF ! mpime
        CALL mp_bcast( xkf_red, ionode_id, inter_pool_comm )
        CALL mp_bcast( s_BZtoIBZ, ionode_id, inter_pool_comm )
        CALL mp_bcast( BZtoIBZ, ionode_id, inter_pool_comm )
        ! 
        DO ik = 1, nkf
          !
          DO iiq=1, nqtotf
            ! 
            CALL kpmq_map( xkf_red(:,ik+lower_bnd-1), xqf (:, iiq), +1, nkq_abs )
            ! 
            ! We want to map k+q onto the full fine k and keep the symm that bring
            ! that point onto the IBZ one.
            s_BZtoIBZ_full(:,:,ik,iiq) = s_BZtoIBZ(:,:,nkq_abs)  
            !
            ixkqf_tr(ik,iiq) = BZtoIBZ(nkq_abs) 
            ! 
          ENDDO ! q-loop
        ENDDO ! k-loop
      ENDIF ! mp_mesh_k
      ! 
    ENDIF ! first_time
    ! 
    !print*,'Start iterative -----------------', iq
    !print*,'Fi_all ',SUM(Fi_all)
    !print*,'Fi_allcb ',SUM(Fi_allcb)
    !print*,'F_current ',sum(F_current)
    !print*,'F_currentcb ',sum(F_currentcb)
    !print*,'F_SERTA ',sum(F_SERTA)
    !print*,'F_SERTAcb ',sum(F_SERTAcb)
    !
    ! In the case of a restart do not add the first step
    IF (first_cycle) THEN
      first_cycle = .FALSE.
    ELSE
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        inv_etemp = 1.0/etemp
        inv_degaussw = 1.0/degaussw
        ! 
        IF(mp_mesh_k) THEN ! Use IBZ k-point grid
          !
          DO ik = 1, nkf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            ! 
            ! We are not consistent with ef from ephwann_shuffle but it should not 
            ! matter if fstick is large enough.
            IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
                 ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
              !
              DO imode = 1, nmodes
                !
                ! the phonon frequency and bose occupation
                wq = wf (imode, iq)
                !
                ! SP : Avoid if statement in inner loops
                ! the coupling from Gamma acoustic phonons is negligible
                IF ( wq .gt. eps_acustic ) THEN
                  inv_wq = 1.0/( two * wq )
                  wgq    = wgauss( -wq*inv_etemp, -99)
                  wgq    = wgq / ( one - two * wgq )
                  g2_tmp = 1.0
                ELSE
                  inv_wq = 0.0
                  wgq    = 0.0
                  g2_tmp = 0.0
                ENDIF
                !
                DO ibnd = 1, ibndmax-ibndmin+1
                  !
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  !
                  DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  energy and fermi occupation at k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0(itemp)
                    fmkq = wgauss( -ekq*inv_etemp, -99)
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    g2 = (abs(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                    !
                    ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                    w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                    w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                    !
                    trans_prob = pi * wqf(iq) * g2 * & 
                                 ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                    !
                    CALL cryst_to_cart(1,Fi_all(:,jbnd,ixkqf_tr(ik,iq),itemp),at,-1)
  
                    CALL dgemv( 'n', 3, 3, 1.d0,&
                        REAL(s_BZtoIBZ_full(:,:,ik,iq), kind=DP), 3, Fi_all(:,jbnd,ixkqf_tr(ik,iq),itemp),1 ,0.d0 , Fi_rot(:), 1 )       
                    CALL cryst_to_cart(1,Fi_all(:,jbnd,ixkqf_tr(ik,iq),itemp),bg,1)
                    CALL cryst_to_cart(1,Fi_rot,bg,1)
                    ! 
                    F_current(:,ibnd,ik+lower_bnd-1,itemp) = F_current(:,ibnd,ik+lower_bnd-1,itemp) +&
                                 two * trans_prob * Fi_rot
                    ! 
                  ENDDO !jbnd
                  !
                  IF ( ABS(efcb(itemp)) > eps6 ) THEN
                    ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                    !
                    DO jbnd = 1, ibndmax-ibndmin+1
                      !
                      !  energy and fermi occupation at k+q
                      ekq = etf (ibndmin-1+jbnd, ikq) - efcb(itemp)
                      fmkq = wgauss( -ekq*inv_etemp, -99)
                      !
                      ! here we take into account the zero-point sqrt(hbar/2M\omega)
                      ! with hbar = 1 and M already contained in the eigenmodes
                      ! g2 is Ry^2, wkf must already account for the spin factor
                      !
                      g2 = (abs(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                      !
                      ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                      w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                      w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                      !
                      trans_prob = pi * wqf(iq) * g2 * &
                                   ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                      !
                      CALL cryst_to_cart(1,Fi_allcb(:,jbnd,ixkqf_tr(ik,iq),itemp),at,-1)

                      CALL dgemv( 'n', 3, 3, 1.d0,&
                          REAL(s_BZtoIBZ_full(:,:,ik,iq), kind=DP), 3, Fi_allcb(:,jbnd,ixkqf_tr(ik,iq),itemp),1 ,0.d0 , Fi_rot(:), 1 )
                      CALL cryst_to_cart(1,Fi_allcb(:,jbnd,ixkqf_tr(ik,iq),itemp),bg,1)
                      CALL cryst_to_cart(1,Fi_rot,bg,1)
                      ! 
                      F_currentcb(:,ibnd,ik+lower_bnd-1,itemp) = F_currentcb(:,ibnd,ik+lower_bnd-1,itemp) +&
                                   two * trans_prob * Fi_rot
                      ! 
                    ENDDO !jbnd
                    ! 
                  ENDIF !  efcb
                  !
                ENDDO !ibnd
                !
              ENDDO !imode
              !
            ENDIF ! endif  fsthick
            !
          ENDDO ! end loop on k
          ! 
        ELSE ! Now the case with FULL k-point grid. 
          ! We need to recast xkf_all with only the full k point (not all k and k+q)
          DO ik = 1, nkqtotf/2
            ikk = 2 * ik - 1
            xkf_red(:,ik) = xkf_all(:,ikk)
          ENDDO
          ! We do some code dupplication wrt to above to avoid branching in a loop.
          ! 
          DO ik = 1, nkf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            ! 
            ! We need to find F_{mk+q}^i (Fi_all). The grids need to be commensurate !
            !CALL ktokpmq ( xk (:, ik), xq, +1, ipool, nkq, nkq_abs )
            xxq = xqf (:, iq)
            xkk = xkf (:, ikk)
            CALL cryst_to_cart (1, xkk, bg, +1)
            CALL cryst_to_cart (1, xxq, bg, +1)
  
            !xkq = xkk + xxq
            !
            ! Note: In this case, Fi_all contains all the k-point across all pools. 
            ! Therefore in the call below, ipool and nkq are dummy variable.
            ! We only want the global index for k+q ==> nkq_abs  
            CALL ktokpmq_fine ( xkf_red ,xkk, xxq, +1, ipool, nkq, nkq_abs )
            ! 
            ! We are not consistent with ef from ephwann_shuffle but it should not 
            ! matter if fstick is large enough.
            IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
                 ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
              !
              DO imode = 1, nmodes
                !
                ! the phonon frequency and bose occupation
                wq = wf (imode, iq)
                !
                ! SP : Avoid if statement in inner loops
                ! the coupling from Gamma acoustic phonons is negligible
                IF ( wq .gt. eps_acustic ) THEN
                  inv_wq = 1.0/( two * wq )
                  wgq    = wgauss( -wq*inv_etemp, -99)
                  wgq    = wgq / ( one - two * wgq )
                  g2_tmp = 1.0
                ELSE
                  inv_wq = 0.0
                  wgq    = 0.0
                  g2_tmp = 0.0
                ENDIF
                !
                DO ibnd = 1, ibndmax-ibndmin+1
                  !
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  !
                  DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  energy and fermi occupation at k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0(itemp)
                    fmkq = wgauss( -ekq*inv_etemp, -99)
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    g2 = (abs(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                    !
                    ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                    w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                    w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                    !
                    trans_prob = pi * wqf(iq) * g2 * &
                                 ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                    !
                    ! IBTE
                    F_current(:,ibnd,ik+lower_bnd-1,itemp) = F_current(:,ibnd,ik+lower_bnd-1,itemp) +&
                                                          two * trans_prob * Fi_all(:,jbnd,nkq_abs,itemp)
                    ! 
                  ENDDO !jbnd
                  !
                  IF ( ABS(efcb(itemp)) > eps6 ) THEN
                    !  energy at k (relative to Ef)
                    ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                    !
                    DO jbnd = 1, ibndmax-ibndmin+1
                      !
                      !  energy and fermi occupation at k+q
                      ekq = etf (ibndmin-1+jbnd, ikq) - efcb(itemp)
                      fmkq = wgauss( -ekq*inv_etemp, -99)
                      !
                      ! here we take into account the zero-point sqrt(hbar/2M\omega)
                      ! with hbar = 1 and M already contained in the eigenmodes
                      ! g2 is Ry^2, wkf must already account for the spin factor
                      !
                      g2 = (abs(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                      !
                      ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                      w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                      w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                      !
                      trans_prob = pi * wqf(iq) * g2 * &
                                   ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                      !
                      ! IBTE
                      F_currentcb(:,ibnd,ik+lower_bnd-1,itemp) = F_currentcb(:,ibnd,ik+lower_bnd-1,itemp) +&
                                                            two * trans_prob * Fi_allcb(:,jbnd,nkq_abs,itemp)
                      ! 
                    ENDDO !jbnd
                    !
                  ENDIF
                  ! 
                ENDDO !ibnd
                !
              ENDDO !imode
              !
            ENDIF ! endif  fsthick
            !
          ENDDO ! end loop on k
          !
        ENDIF ! mp_mesh_k
        !
      ENDDO ! itemp 
      !  
      ! Creation of a restart point
      IF (restart) THEN
        IF ( MOD(iq,restart_freq) == 0 .or. iq == nqtotf ) THEN
          WRITE(stdout, '(a)' ) '     Creation of a restart point'
          ! 
          ! The mp_sum will aggreage the results on each k-points. 
          CALL mp_sum( F_current, inter_pool_comm )
          !
          IF ( ABS(efcb(1)) > eps6 ) THEN
            CALL mp_sum( F_currentcb, inter_pool_comm)
            CALL F_write(iter, iq, nqtotf, nkqtotf/2, error_h, error_el, .TRUE.)
          ELSE
            CALL F_write(iter, iq, nqtotf, nkqtotf/2, error_h, error_el, .FALSE.)
          ENDIF
          ! 
        ENDIF
      ENDIF
      !
    ENDIF ! first_cycle
    !  
    ! 
    ! The k points are distributed among pools: here we collect them
    !
    IF ( iq .eq. nqtotf ) THEN
      !
      DO itemp = 1, nstemp
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
            DO ibnd = 1, ibndmax-ibndmin+1
              tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
              F_current(:,ibnd,ik+lower_bnd-1,itemp) = F_SERTA(:,ibnd,ik+lower_bnd-1,itemp) +&
                                                    tau * F_current(:,ibnd,ik+lower_bnd-1,itemp)
              IF ( ABS(efcb(itemp)) > eps6 ) THEN
                tau = one / inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1)
                F_currentcb(:,ibnd,ik+lower_bnd-1,itemp) = F_SERTAcb(:,ibnd,ik+lower_bnd-1,itemp) +&
                                                    tau * F_currentcb(:,ibnd,ik+lower_bnd-1,itemp) 
              ENDIF 
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      !
      CALL mp_sum( F_current, inter_pool_comm )
      IF ( ABS(efcb(1)) > eps6 ) CALL mp_sum( F_currentcb, inter_pool_comm )
      !
      ! The next Fi is equal to the current Fi+1 F_current. 
      Fi_all = F_current
      F_current = zero
      ! 
      IF ( ABS(efcb(1)) > eps6 ) THEN
        Fi_allcb = F_currentcb
        F_currentcb = zero
      ENDIF
      !
      ! From the F, we compute the HOLE conductivity
      IF (int_mob .OR. (ncarrier < -1E5)) THEN
        WRITE(stdout,'(/5x,a)') repeat('=',67)
        WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
        WRITE(stdout,'(5x,a/)') repeat('=',67)
        ! 
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp) 
          ! 
          IF ( itemp .eq. 1 ) THEN
            Sigma(:,:)   = zero
            tdf_sigma(:) = zero
          ENDIF
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
              !  
              DO ibnd = 1, ibndmax-ibndmin+1
                ! This selects only valence bands for hole conduction
                IF (etf (ibndmin-1+ibnd, ikk) < ef0(itemp) ) THEN
                  IF ( vme ) THEN 
                    vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
                  ELSE
                    vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
                  ENDIF
                  ! 
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      tdf_sigma(ij) = vkk(i,ibnd) * Fi_all(j,ibnd,ik+lower_bnd-1,itemp)
                    ENDDO
                  ENDDO
                  ! 
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  !  
                  ! derivative Fermi distribution
                  ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                  dfnk = w0gauss( ekk / etemp, -99 ) / etemp          
                  !
                  ! electrical conductivity
                  Sigma(:,itemp) = Sigma(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
                ENDIF
              ENDDO ! iband
            ENDIF ! fsthick
          ENDDO ! ik
          !
          ! The k points are distributed among pools: here we collect them
          !
          CALL mp_sum( Sigma(:,:), inter_pool_comm )
          CALL mp_barrier(inter_pool_comm)
          !
          carrier_density = 0.0
          ! 
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only valence bands for hole conduction
              IF (etf (ibndmin-1+ibnd, ikk) < ef0(itemp) ) THEN
                !  energy at k (relative to Ef)
                ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                fnk = wgauss( -ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf(ikk) * (1.0d0 - fnk )
              ENDIF
            ENDDO
          ENDDO
          ! 
          CALL mp_sum( carrier_density, inter_pool_comm )
          CALL mp_barrier(inter_pool_comm)
          ! 
          sigma_up(:,:) = zero
          sigma_up(1,1) = Sigma(1,itemp)
          sigma_up(1,2) = Sigma(2,itemp)
          sigma_up(1,3) = Sigma(3,itemp)
          sigma_up(2,1) = Sigma(4,itemp)
          sigma_up(2,2) = Sigma(5,itemp)
          sigma_up(2,3) = Sigma(6,itemp)
          sigma_up(3,1) = Sigma(7,itemp)
          sigma_up(3,2) = Sigma(8,itemp)
          sigma_up(3,3) = Sigma(9,itemp)
          !
          ! Diagonalize the conductivity matrix
          CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
          !
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          ! carrier_density in cm^-1
          carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)         
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                        ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
          ! 
          error_h = ABS(mobility-mobilityh_save(itemp))
          mobilityh_save(itemp) = mobility
          WRITE(stdout,'(47x, 1E16.4, a)') error_h, '     Err'
          !
        ENDDO ! itemp
      ENDIF ! holes mobility
      ! From the F, we compute the ELECTRON conductivity
      IF (int_mob .OR. (ncarrier > 1E5)) THEN
        WRITE(stdout,'(/5x,a)') repeat('=',67)
        WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
        WRITE(stdout,'(5x,a/)') repeat('=',67)
        ! 
        Sigma(:,:)   = zero
        tdf_sigma(:) = zero 
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp)
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
              DO ibnd = 1, ibndmax-ibndmin+1
                IF ( ABS(efcb(itemp)) < eps6 ) THEN ! Case with 1 Fermi level
                  ! This selects only conduction bands for electron conduction
                  IF (etf (ibndmin-1+ibnd, ikk) > ef0(itemp) ) THEN
                    IF ( vme ) THEN 
                      vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
                    ELSE
                      vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
                    ENDIF
                    ! 
                    ij = 0
                    DO j = 1, 3
                      DO i = 1, 3
                        ij = ij + 1
                        tdf_sigma(ij) = vkk(i,ibnd) * Fi_all(j,ibnd,ik+lower_bnd-1,itemp)
                      ENDDO
                    ENDDO
                    ! 
                    !  energy at k (relative to Ef)
                    ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                    !  
                    ! derivative Fermi distribution
                    ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                    dfnk = w0gauss( ekk / etemp, -99 ) / etemp          
                    !
                    ! electrical conductivity
                    Sigma(:,itemp) = Sigma(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
                    !
                  ENDIF
                ELSE ! In this case we have 2 Fermi level
                  IF (etf (ibndmin-1+ibnd, ikk) > efcb(itemp) ) THEN
                    IF ( vme ) THEN
                      vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
                    ELSE
                      vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
                    ENDIF
                    ! 
                    ij = 0
                    DO j = 1, 3
                      DO i = 1, 3
                        ij = ij + 1
                        tdf_sigma(ij) = vkk(i,ibnd) * Fi_all(j,ibnd,ik+lower_bnd-1,itemp)
                      ENDDO
                    ENDDO
                    ! 
                    !  energy at k (relative to Ef)
                    ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                    !  
                    ! derivative Fermi distribution
                    ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                    dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                    !
                    ! electrical conductivity
                    Sigma(:,itemp) = Sigma(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
                    !
                  ENDIF ! efcb
                ENDIF 
              ENDDO ! iband
            ENDIF ! fsthick
          ENDDO ! ik
          !
          ! The k points are distributed among pools: here we collect them
          !
          CALL mp_sum( Sigma(:,:), inter_pool_comm )
          CALL mp_barrier(inter_pool_comm)
          ! 
          carrier_density = 0.0
          ! 
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only valence bands for hole conduction
              IF ( ABS(efcb(itemp)) < eps6 ) THEN
                IF (etf (ibndmin-1+ibnd, ikk) > ef0(itemp) ) THEN
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  fnk = wgauss( -ekk / etemp, -99)
                  ! The wkf(ikk) already include a factor 2
                  carrier_density = carrier_density + wkf(ikk) * fnk
                ENDIF
              ELSE
                IF (etf (ibndmin-1+ibnd, ikk) > efcb(itemp) ) THEN
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                  fnk = wgauss( -ekk / etemp, -99)
                  ! The wkf(ikk) already include a factor 2
                  carrier_density = carrier_density + wkf(ikk) * fnk
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          ! 
          CALL mp_sum( carrier_density, inter_pool_comm )
          CALL mp_barrier(inter_pool_comm)
          sigma_up(:,:) = zero
          sigma_up(1,1) = Sigma(1,itemp)
          sigma_up(1,2) = Sigma(2,itemp)
          sigma_up(1,3) = Sigma(3,itemp)
          sigma_up(2,1) = Sigma(4,itemp)
          sigma_up(2,2) = Sigma(5,itemp)
          sigma_up(2,3) = Sigma(6,itemp)
          sigma_up(3,1) = Sigma(7,itemp)
          sigma_up(3,2) = Sigma(8,itemp)
          sigma_up(3,3) = Sigma(9,itemp)
          !
          ! Diagonalize the conductivity matrix
          CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
          !
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          ! carrier_density in cm^-1
          carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)         
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                           ef0(itemp)*ryd2ev, carrier_density, mobility_xx, '  x-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
          error_el = ABS(mobility-mobilityel_save(itemp))
          mobilityel_save(itemp) = mobility
          WRITE(stdout,'(47x, 1E16.4, a)') error_el, '     Err'
          ! 
        ENDDO ! itemp
      ENDIF ! Electron mobility
      !
    ENDIF ! iq == nq
    !
    RETURN
    !
    ! ---------------------------------------------------------------------------
    END SUBROUTINE iterativebte
    !----------------------------------------------------------------------------
    ! 
  END MODULE transport_iter
