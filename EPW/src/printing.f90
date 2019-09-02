  !                                                                            
  ! Copyright (C) 2010-2019 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  MODULE printing
  !----------------------------------------------------------------------
  !! 
  !! This module contains various printing routines
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE k_avg(F_out, vkk_all, nb_sp, xkf_sp)
    !-----------------------------------------------------------------------
    !! 
    !! This routines enforces symmetry.
    !! Averages points which leaves the k-point unchanged by symmetry
    !!   e.g. k=[1,1,1] and q=[1,0,0] with the symmetry that change x and y gives 
    !!        k=[1,1,1] and q=[0,1,0].
    !! 
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nstemp
    USE elph2,         ONLY : nkqtotf, ibndmax, ibndmin, nkf, nbndfst
    USE cell_base,     ONLY : bg, at
    USE constants_epw, ONLY : eps6, zero
    USE symm_base,     ONLY : s, nrot
    USE division,      ONLY : fkbounds
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : world_comm
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nb_sp
    !! Lenght of xkf_sp
    INTEGER, INTENT(in) :: xkf_sp(49, nb_sp)
    !! Special points indexes and symmetries
    REAL(KIND = DP), INTENT(inout) :: F_out(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i
    REAL(KIND = DP), INTENT(inout) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity of k
    ! 
    ! Local variables
    LOGICAL :: special_map(nkf)
    !! Special mapping
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: lower_bnd
    !! Lower bounds index after k para
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: sp
    !! Local index
    INTEGER ::nb
    !! Local index
    LOGICAL :: special
    !! Local logical
    INTEGER :: counter_average
    !! Local counter
    INTEGER :: index_sp(nkf)
    !! Index of special points
    REAL(KIND = DP) :: xkk_cart (3)
    !! k-point coordinate in Cartesian unit
    REAL(KIND = DP) :: mean(nbndfst)
    !! Mean of the velocities
    REAL(KIND = DP) :: mean_pop(nbndfst)
    !! Mean of the populations
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: S_vkk(3, nbndfst)
    !! Rotated vector
    REAL(KIND = DP) :: S_F_out(3, nbndfst)
    !! Rotated vector
    REAL(KIND = DP) :: tmp_vkk(3, nbndfst)
    !! Temporary vector
    REAL(KIND = DP) :: tmp_F_out(3, nbndfst)
    !! Temporary vector
    REAL(KIND = DP) :: F_out_loc(3, nbndfst, nktotf, nstemp)
    !! Local F_out where the k-points have been spread
    REAL(KIND = DP) :: vkk_all_loc(3, nbndfst, nktotf)
    !! Local velocity where the k-points have been spread
    ! 
    ! Split the k-point across cores
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    ! 
    F_out_loc(:, :, :, :) = zero
    vkk_all_loc(:, :, :) = zero
    special_map(:) = .FALSE. 
    index_sp(:) = 0
    DO ik = 1, nkf
      DO sp = 1,nb_sp
        IF (ik + lower_bnd - 1 == xkf_sp(1, sp)) THEN
          special_map(ik) = .TRUE. 
          index_sp(ik) = sp
        ENDIF
      ENDDO
    ENDDO ! ik
    !  
    DO itemp = 1, nstemp
      DO ik = 1, nkf
        IF (special_map(ik)) THEN
          counter_average = 0
          tmp_vkk = zero
          tmp_F_out = zero
          DO nb = 1,nrot
            IF (index_sp(ik) > 0) THEN 
              IF (xkf_sp(nb + 1, index_sp(ik)) > 0) THEN
                counter_average = counter_average + 1
                sa(:, :) = DBLE(s(:, :, xkf_sp(nb + 1, index_sp(ik))))
                sb       = MATMUL(bg, sa)
                sr(:, :) = MATMUL(at, TRANSPOSE(sb))
                sr       = TRANSPOSE(sr)
                DO ibnd = 1, nbndfst
                  CALL dgemv('n', 3, 3, 1.d0, sr, 3, vkk_all(:, ibnd, ik + lower_bnd - 1), 1, 0.d0, S_vkk(:, ibnd), 1)
                  CALL dgemv('n', 3, 3, 1.d0, sr, 3, F_out(:, ibnd, ik + lower_bnd - 1, itemp), 1, 0.d0, S_F_out(:, ibnd), 1)
                  tmp_vkk(:, ibnd) = tmp_vkk(:, ibnd) + S_vkk(:, ibnd)
                  tmp_F_out(:, ibnd) = tmp_F_out(:, ibnd) + S_F_out(:, ibnd)
                ENDDO ! ibnd
              ENDIF
            ENDIF
          ENDDO ! sp
          DO ibnd = 1, nbndfst
            vkk_all_loc(:, ibnd, ik + lower_bnd - 1) = tmp_vkk(:, ibnd) / DBLE(counter_average)
            F_out_loc(:, ibnd, ik + lower_bnd - 1, itemp) = tmp_F_out(:, ibnd) / DBLE(counter_average)
          ENDDO
          ! 
        ELSE ! not a special point 
          DO ibnd = 1, nbndfst
            vkk_all_loc(:, ibnd, ik + lower_bnd - 1) = vkk_all(:, ibnd, ik + lower_bnd - 1)
            F_out_loc(:, ibnd, ik + lower_bnd - 1, itemp) = F_out(:, ibnd, ik + lower_bnd - 1, itemp) 
          ENDDO 
        ENDIF ! special
      ENDDO! ik
    ENDDO! itemp
    ! 
    ! Gather from all cores
    CALL mp_sum(vkk_all_loc, world_comm)
    CALL mp_sum(F_out_loc, world_comm)
    ! 
    F_out = F_out_loc
    vkk_all = vkk_all_loc
    !  
    !-----------------------------------------------------------------------
    END SUBROUTINE k_avg
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE print_gkk(iq)
    !-----------------------------------------------------------------------
    !! 
    !! Print the |g| vertex for all n,n' and modes in meV and do average
    !! on degenerate states. 
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, nbndfst, &
                              nkf, epf17, xkf, nkqtotf, wf, adapt_smearing
    USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, eps8
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE division,      ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-point index 
    !
    ! Local variables 
    ! 
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! K-point index
    INTEGER :: ikq
    !! K+q-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: pbnd
    !! Band index
    INTEGER :: nu
    !! Mode index
    INTEGER :: mu
    !! Mode index
    INTEGER :: n
    !! Number of modes
    !
    REAL(KIND = DP) :: xkf_all(3, nkqtotf)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP) :: etf_all(nbndsub, nkqtotf)
    !! Collect eigenenergies from all pools in parallel case
    REAL(KIND = DP) :: wq
    !! Phonon frequency 
    REAL(KIND = DP) :: w_1
    !! Temporary phonon freq. 1
    REAL(KIND = DP) :: w_2
    !! Temporary phonon freq. 2
    REAL(KIND = DP) :: gamma
    !! Temporary electron-phonon matrix element
    REAL(KIND = DP) :: ekk
    !! Eigenenergies at k
    REAL(KIND = DP) :: ekq
    !! Eigenenergies at k+q
    REAL(KIND = DP) :: g2
    !! Temporary electron-phonon matrix element square
    REAL(KIND = DP), ALLOCATABLE :: epc(:, :, :, :)
    !! g vectex accross all pools 
    REAL(KIND = DP), ALLOCATABLE :: epc_sym(:, :, :)
    !! Temporary g-vertex for each pool 
    ! 
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    ! 
    ALLOCATE(epc(nbndfst, nbndfst, nmodes, nktotf))
    ALLOCATE(epc_sym(nbndfst, nbndfst, nmodes))
    ! 
    epc(:, :, :, :)  = zero
    epc_sym(:, :, :) = zero
    !
    ! First do the average over bands and modes for each pool
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      ikq = ikk + 1
      ! 
      DO nu = 1, nmodes
        wq = wf (nu, iq)
        DO ibnd = 1, nbndfst
          DO jbnd = 1, nbndfst
            gamma = (ABS(epf17(jbnd, ibnd, nu, ik)))**two
            IF (wq > 0.d0) THEN
              gamma = gamma / (two * wq)
            ELSE
              gamma = 0.d0
            ENDIF
            gamma = SQRT(gamma)
            ! gamma = |g| [Ry]
            epc(ibnd, jbnd, nu, ik + lower_bnd - 1) = gamma
          ENDDO ! jbnd
        ENDDO   ! ibnd        
      ENDDO ! loop on modes
      !
      !  Here we "SYMMETRIZE": actually we simply take the averages over
      !  degenerate states, it is only a convention because g is gauge-dependent!
      !
      ! first the phonons
      DO ibnd = 1, nbndfst
        DO jbnd = 1, nbndfst
          DO nu = 1, nmodes
            w_1 = wf(nu, iq)
            g2 = 0.d0
            n  = 0
            DO mu = 1, nmodes
              w_2 = wf(mu, iq)
              IF (ABS(w_2 - w_1) < eps8) THEN
                n = n + 1
                g2 = g2 + epc(ibnd, jbnd, mu, ik + lower_bnd - 1) * epc(ibnd, jbnd, mu, ik + lower_bnd - 1)
              ENDIF
            ENDDO
            g2 = g2 / FLOAT(n)
            epc_sym(ibnd, jbnd, nu) = SQRT(g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:, :, :, ik + lower_bnd - 1) = epc_sym
      ! Then the k electrons
      DO nu = 1, nmodes
        DO jbnd = 1, nbndfst
          DO ibnd = 1, nbndfst
            w_1 = etf(ibndmin - 1 + ibnd, ikk)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, nbndfst
              w_2 = etf(ibndmin - 1 + pbnd, ikk)
              IF (ABS(w_2-w_1) < eps8) THEN
                n = n + 1
                g2 = g2 + epc(pbnd, jbnd, nu, ik + lower_bnd - 1) * epc(pbnd, jbnd, nu, ik + lower_bnd - 1)
              ENDIF
            ENDDO
            g2 = g2 / FLOAT(n)
            epc_sym(jbnd, ibnd, nu) = SQRT(g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:, :, :, ik + lower_bnd - 1) = epc_sym
      !
      ! and finally the k+q electrons
      DO nu = 1, nmodes
        DO ibnd = 1, nbndfst
          DO jbnd = 1, nbndfst
            w_1 = etf(ibndmin - 1 + jbnd, ikq)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, nbndfst
              w_2 = etf(ibndmin - 1 + pbnd, ikq)
              IF (ABS(w_2 - w_1) < eps8) then
                n = n + 1
                g2 = g2 + epc(ibnd, pbnd, nu, ik + lower_bnd - 1) * epc(ibnd, pbnd, nu, ik + lower_bnd - 1)
              ENDIF
            ENDDO
            g2 = g2 / FLOAT(n)
            epc_sym(ibnd, jbnd, nu) = SQRT(g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:, :, :, ik + lower_bnd - 1) = epc_sym
      ! 
    ENDDO ! k-points
    ! 
    ! We need quantity from all the pools
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    !
#if defined(__MPI)
    !
    ! Note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    CALL mp_sum(epc, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    !
#endif
    !
    ! Only master writes
    IF (mpime == ionode_id) THEN
      !
      WRITE(stdout, '(5x,a)') ' Electron-phonon vertex |g| (meV)'
      !
      WRITE(stdout, '(/5x,"iq = ",i7," coord.: ", 3f12.7)') iq, xqf(:, iq)
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout,'(5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
        WRITE(stdout, '(5x,a)') ' ibnd     jbnd     imode   enk[eV]    enk+q[eV]  omega(q)[meV]   |g|[meV]'
        WRITE(stdout,'(5x,a)') REPEAT('-',78)
        !
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) 
          DO jbnd = 1, nbndfst
            ekq = etf_all(ibndmin - 1 + jbnd, ikq) 
            DO nu = 1, nmodes
              WRITE(stdout,'(3i9,3f12.4,1e20.10)') ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, nu, ryd2ev * ekk, ryd2ev * ekq, &
                                           ryd2mev * wf(nu, iq), ryd2mev * epc(ibnd, jbnd, nu, ik)
            ENDDO
          ENDDO  
          !
        ENDDO
        WRITE(stdout,'(5x,a/)') REPEAT('-',78)
        !
      ENDDO
    ENDIF ! master node
    ! 
    DEALLOCATE(epc)
    DEALLOCATE(epc_sym)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_gkk
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_serta_sym(F_SERTA, BZtoIBZ, s_BZtoIBZ, BZtoIBZ_mat, vkk_all, etf_all, wkf_all, ef0) 
    !-----------------------------------------------------------------------
    !!
    !!  This routine prints the SERTA mobility using k-point symmetry.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, omega, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf, nbndfst 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              electron_SI, bohr2ang, ang2cm, hbarJ
    USE symm_base,     ONLY : s, nrot
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm, my_pool_id
    USE mp,            ONLY : mp_sum
    USE kinds_epw,     ONLY : SIK2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER(SIK2), INTENT(in) :: s_BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! Corresponding symmetry matrix
    INTEGER, INTENT(in) :: BZtoIBZ_mat(nrot, nktotf)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(KIND = DP), INTENT(in) :: F_SERTA(3, nbndfst, nktotf, nstemp)  
    !! SERTA solution 
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level 
    ! 
    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik 
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: nb 
    !! Number of points equivalent by sym from BZ to IBTZ
    INTEGER :: ikbz
    !! Index of full BZ points
    INTEGER :: ij  
    !! Combined x,y,z index
    INTEGER :: i,j
    !! Cartesian index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma(3, 3, nstemp)
    !! Electrical conductivity
    REAL(KIND = DP) :: Fi_cart(3)
    !! Cartesian Fi_all 
    REAL(KIND = DP) :: Fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    REAL(KIND = DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(KIND = DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate 
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP) :: mobility
    !! Sum of the diagonal mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_od
    !! Sum of the off-diagonal mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs]
    REAL(KIND = DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND = DP) :: sigma_vect(3, 3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    REAL(KIND = DP) :: Fi_check(3, nstemp)
    !! Sum rule on population
    REAL(KIND = DP) :: sfac
    !! Spin factor
    ! 
    IF (noncolin) THEN
      sfac = 1.0
    ELSE
      sfac = 2.0
    ENDIF
    Fi_check(:, :) = zero
    ! 
    ! Hole
    IF (ncarrier < -1E5) THEN
      sigma(:, :, :) = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Hole density  Population SR                Hole mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [h per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) < ef0(itemp)) THEN
              !
              vk_cart(:) = vkk_all(:, ibnd, ik)
              Fi_cart(:) = F_SERTA(:, ibnd, ik, itemp)
              !
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb = 1, nrot
                IF (BZtoIBZ_mat(nb, ik) > 0) THEN 
                  ikbz = BZtoIBZ_mat(nb, ik)
                  ! 
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:, :) = DBLE(s(:, :, s_BZtoIBZ(ikbz)))
                  sb        = MATMUL(bg, sa)
                  sr (:, :) = MATMUL(at, TRANSPOSE(sb))
                  sr        = TRANSPOSE(sr) 
                  !
                  CALL dgemv( 'n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1 ,0.d0 , v_rot(:), 1 )
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0, sr, 3, Fi_cart(:), 1 ,0.d0 , Fi_rot(:), 1 )
                  !
                  DO j = 1, 3
                    DO i = 1, 3
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      sigma(i, j, itemp) = sigma(i, j, itemp) - (v_rot(j) * Fi_rot(i)) * sfac / (nkf1 * nkf2 * nkf3)
                    ENDDO
                  ENDDO
                  !
                  Fi_check(:, itemp) = Fi_check(:, itemp) + Fi_rot(:) * sfac / (nkf1 * nkf2 * nkf3)
                ENDIF ! BZ 
              ENDDO ! ikb
              ! 
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for hole conduction
            IF (etf_all(ibnd, ik) < ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk)
            ENDIF
          ENDDO
        ENDDO
        ! 
        ! Print the resulting mobility
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp))
        ! 
      ENDDO
      ! 
    ENDIF
    ! 
    ! Now electron mobility
    IF (ncarrier > 1E5) THEN
      ! Needed because of residual values from the hole above
      sigma(:, :, :)    = zero
      ! 
      WRITE(stdout, '(/5x,a)') REPEAT('=',71)
      WRITE(stdout, '(5x,"  Temp     Fermi   Elec density  Population SR                Elec mobility ")')
      WRITE(stdout, '(5x,"   [K]      [eV]     [cm^-3]      [e per cell]                  [cm^2/Vs]")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              vk_cart(:) = vkk_all(:, ibnd, ik)
              Fi_cart(:) = F_SERTA(:, ibnd, ik, itemp) 
              ! 
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb = 1, nrot
                IF (BZtoIBZ_mat(nb, ik) > 0) THEN
                  ikbz = BZtoIBZ_mat(nb, ik)
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:, :) = DBLE(s(:, :, s_BZtoIBZ(ikbz)))
                  sb        = MATMUL(bg, sa)
                  sr (:, :) = MATMUL(at, TRANSPOSE(sb))
                  sr        = TRANSPOSE(sr)
                  CALL dgemv( 'n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1 ,0.d0 , v_rot(:), 1)
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0, sr, 3, Fi_cart(:), 1 ,0.d0 , Fi_rot(:), 1)
                  !
                  DO j = 1, 3
                    DO i = 1, 3
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      sigma(i, j, itemp) = sigma(i, j, itemp) - (v_rot(j) * Fi_rot(i)) * sfac / (nkf1 * nkf2 * nkf3)
                    ENDDO
                  ENDDO
                  !
                  Fi_check(:, itemp) = Fi_check(:, itemp) + Fi_rot(:) * sfac / (nkf1 * nkf2 * nkf3)
                ENDIF ! BZ
              ENDDO ! ikb
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for electron conduction
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        ! 
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp)) 
        !
      ENDDO
      ! 
    ENDIF
    !  
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_serta_sym
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE print_serta(F_SERTA, vkk_all, etf_all, wkf_all, ef0) 
    !-----------------------------------------------------------------------
    !!
    !!  This routine prints the SERTA mobility without k-point symmetry
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, omega, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf, nbndfst 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                               electron_SI, bohr2ang, ang2cm, hbarJ
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm, my_pool_id
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: F_SERTA(3, nbndfst, nktotf, nstemp)  
    !! SERTA solution 
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level 
    ! 
    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik 
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: ij  
    !! Combined x,y,z index
    INTEGER :: i, j
    !! Cartesian index
    ! 
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma(3, 3, nstemp)
    !! Electrical conductivity
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP) :: mobility
    !! Sum of the diagonal mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_od
    !! Sum of the off-diagonal mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs]
    REAL(KIND = DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND = DP) :: sigma_vect(3, 3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND = DP) :: Fi_check(3, nstemp)
    !! Sum rule on population
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    REAL(KIND = DP) :: sfac
    !! Spin factor
    !
    Fi_check(:, :) = zero
    IF (noncolin) THEN
      sfac = 1.0
    ELSE
      sfac = 2.0
    ENDIF
    ! 
    ! Hole
    IF (ncarrier < -1E5) THEN
      sigma(:, :, :) = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Hole density  Population SR                Hole mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [h per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) < ef0(itemp)) THEN
              DO j = 1, 3
                DO i = 1, 3
                  sigma(i, j, itemp) = sigma(i, j, itemp) - vkk_all(j, ibnd, ik) * F_SERTA(i, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              Fi_check(:, itemp) = Fi_check(:, itemp) + F_SERTA(:, ibnd, ik, itemp) * sfac / (nkf1 * nkf2 * nkf3)
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for hole conduction
            IF (etf_all(ibnd, ik ) < ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk)
            ENDIF
          ENDDO
        ENDDO
        ! 
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp)) 
        !  
      ENDDO
      ! 
    ENDIF
    ! 
    ! Now electron mobility
    IF (ncarrier > 1E5) THEN
      ! Needed because of residual values from the hole above
      sigma(:, :, :) = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Elec density  Population SR                Elec mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [e per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              DO j = 1, 3
                DO i = 1, 3
                  sigma(i, j, itemp) = sigma(i, j, itemp) - vkk_all(j, ibnd, ik) * F_SERTA(i, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              Fi_check(:, itemp) = Fi_check(:, itemp) + F_SERTA(:, ibnd, ik, itemp) * sfac / (nkf1 * nkf2 * nkf3)
              ! 
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for electron conduction
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        !
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp))
        ! 
      ENDDO
      ! 
    ENDIF
    
    RETURN
    !
    END SUBROUTINE print_serta
    !-----------------------------------------------------------------------
    !  
    !-----------------------------------------------------------------------
    SUBROUTINE print_mob_sym(F_out, BZtoIBZ, s_BZtoIBZ, BZtoIBZ_mat, vkk_all, etf_all, wkf_all, ef0, av_mob) 
    !-----------------------------------------------------------------------
    !!
    !!  This routine prints the iterative mobility using k-point symmetry ( electron or hole )
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf, nbndfst 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              electron_SI, bohr2ang, ang2cm, hbarJ, eps80
    USE symm_base,     ONLY : s, nrot
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm
    USE mp,            ONLY : mp_sum
    USE kinds_epw,     ONLY : SIK2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER(SIK2), INTENT(in) :: s_BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! Corresponding symmetry matrix
    INTEGER, INTENT(in) :: BZtoIBZ_mat(nrot, nktotf)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(KIND = DP), INTENT(in) :: F_out(3, nbndfst, nktotf, nstemp)  
    !! SERTA solution 
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level 
    REAL(KIND = DP), INTENT(inout) :: av_mob(nstemp)
    !! Maximum error for all temperature

    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik 
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: nb 
    !! Number of points equivalent by sym from BZ to IBTZ
    INTEGER :: ikbz
    !! Index of full BZ points
    INTEGER :: ij  
    !! Combined x,y,z index
    INTEGER :: i,j
    !! Cartesian index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma(3, 3, nstemp)
    !! Electrical conductivity
    REAL(KIND = DP) :: Fi_cart(3)
    !! Cartesian Fi_all 
    REAL(KIND = DP) :: Fi_check(3, nstemp)
    !! Check that \Sum_k Fi = (0,0,0)
    REAL(KIND = DP) :: Fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    REAL(KIND = DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(KIND = DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate 
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP) :: mobility
    !! Sum of the diagonal mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_od
    !! Sum of the off-diagonal mobilities [cm^2/Vs]
    REAL(KIND = DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs]
    REAL(KIND = DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND = DP) :: sigma_vect(3, 3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    REAL(KIND = DP) :: sfac
    !! Spin factor
    ! 
    Fi_check(:, :) = zero
    ! 
    IF (noncolin) THEN
      sfac = 1.0
    ELSE
      sfac = 2.0
    ENDIF
    ! 
    IF (ncarrier < -1E5) THEN ! If true print hole
      sigma(:, :, :) = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Hole density  Population SR                Hole mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [h per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all (ibnd, ik) < ef0(itemp)) THEN
              vk_cart(:) = vkk_all(:, ibnd, ik)
              Fi_cart(:) = F_out(:, ibnd, ik, itemp)
              ! 
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb = 1, nrot
                IF (BZtoIBZ_mat(nb, ik) > 0) THEN
                  ikbz = BZtoIBZ_mat(nb, ik)
                  ! 
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:, :) = DBLE(s(:, :, s_BZtoIBZ(ikbz)))
                  sb        = MATMUL(bg, sa)
                  sr (:, :) = MATMUL(at, TRANSPOSE(sb))
                  sr        = TRANSPOSE(sr)
                  ! 
                  CALL dgemv('n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1, 0.d0, v_rot(:), 1)
                  ! 
                  CALL dgemv('n', 3, 3, 1.d0, sr, 3, Fi_cart(:), 1, 0.d0, Fi_rot(:), 1)
                  !
                  DO j = 1, 3
                    DO i = 1, 3
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      sigma(i, j, itemp) = sigma(i, j, itemp) - (v_rot(j) * Fi_rot(i)) * sfac / (nkf1 * nkf2 * nkf3)
                    ENDDO
                  ENDDO
                  !
                  Fi_check(:, itemp) = Fi_check(:, itemp) + Fi_rot(:) * sfac / (nkf1 * nkf2 * nkf3) 
                ENDIF ! BZ
              ENDDO ! ikb
              ! 
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for hole conduction
            IF (etf_all(ibnd, ik) < ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
        ! 
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp))
        ! 
      ENDDO
      ! 
    ENDIF
    ! Now electron mobility
    IF (ncarrier > 1E5) THEN
      ! Needed because of residual values from the hole above
      sigma(:, :, :) = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Elec density  Population SR                Elec mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [e per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              vk_cart(:) = vkk_all(:, ibnd, ik)
              Fi_cart(:) = F_out(:, ibnd, ik, itemp)
              ! 
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb = 1, nrot
                IF (BZtoIBZ_mat(nb, ik) > 0) THEN
                  ikbz = BZtoIBZ_mat(nb, ik)
                  ! 
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:, :) = DBLE(s(:, :, s_BZtoIBZ(ikbz)))
                  sb        = MATMUL(bg, sa)
                  sr (:, :) = MATMUL(at, TRANSPOSE (sb))
                  sr        = TRANSPOSE(sr)
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1, 0.d0 ,v_rot(:), 1 )
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0, sr, 3, Fi_cart(:), 1, 0.d0, Fi_rot(:), 1 )
                  !
                  DO j = 1, 3
                    DO i = 1, 3
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      sigma(i, j, itemp) = sigma(i, j, itemp) - (v_rot(j) * Fi_rot(i)) * sfac / (nkf1 * nkf2 * nkf3)  
                    ENDDO
                  ENDDO 
                  !
                  Fi_check(:, itemp) = Fi_check(:, itemp) + Fi_rot(:) * sfac / (nkf1 * nkf2 * nkf3) 
                ENDIF ! BZ
              ENDDO ! ikb
              ! 
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for electron conduction
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        ! 
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp)) 
        ! 
      ENDDO
      ! 
    ENDIF 
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_mob_sym
    !-----------------------------------------------------------------------
    !  
    !-----------------------------------------------------------------------
    SUBROUTINE print_mob(F_out, vkk_all, etf_all, wkf_all, ef0, av_mob) 
    !-----------------------------------------------------------------------
    !!
    !!  This routine prints the iterative mobility without k-point symmetry ( electron or hole )
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : omega 
    USE epwcom,        ONLY : ncarrier, nstemp, nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf, nbndfst 
    USE transportcom,  ONLY : transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              electron_SI, bohr2ang, ang2cm, hbarJ, eps80
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: F_out(3, nbndfst, nktotf, nstemp)  
    !! SERTA solution 
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level 
    REAL(KIND = DP), INTENT(inout) :: av_mob(nstemp)
    !! Maximum error for all temperature

    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik 
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: nb 
    !! Number of points equivalent by sym from BZ to IBTZ
    INTEGER :: ij  
    !! Combined x,y,z index
    INTEGER :: i,j
    !! Cartesian index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma(3, 3, nstemp)
    !! Electrical conductivity
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP) :: Fi_check(3, nstemp)
    !! Check that \Sum_k Fi = (0,0,0)
    REAL(KIND = DP) :: mobility
    !! Sum of the diagonal mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_od
    !! Sum of the off-diagonal mobilities [cm^2/Vs]
    REAL(KIND = DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs]
    REAL(KIND = DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND = DP) :: sigma_vect(3, 3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    REAL(KIND = DP) :: sfac
    !! Spin factor    
    ! 
    Fi_check(:, :) = zero
    ! 
    IF (noncolin) THEN
      sfac = 1.0
    ELSE
      sfac = 2.0
    ENDIF
    !
    IF (ncarrier < -1E5) THEN ! If true print hole
      sigma(:, :, :)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Hole density  Population SR                Hole mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [h per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) < ef0(itemp)) THEN
              DO j = 1, 3
                DO i = 1, 3
                  sigma(i, j, itemp) = sigma(i, j, itemp) - vkk_all(j, ibnd, ik) * F_out(i, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              Fi_check(:, itemp) = Fi_check(:, itemp) + F_out(:, ibnd, ik, itemp) * sfac / (nkf1 * nkf2 * nkf3)
              ! 
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for hole conduction
            IF (etf_all(ibnd, ik) < ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
        ! 
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp))
        ! 
      ENDDO
      ! 
    ENDIF
    ! Now electron mobility
    IF (ncarrier > 1E5) THEN
      ! Needed because of residual values from the hole above
      sigma(:, :, :)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') REPEAT('=',71)
      WRITE(stdout,'(5x,"  Temp     Fermi   Elec density  Population SR                Elec mobility ")')
      WRITE(stdout,'(5x,"   [K]      [eV]     [cm^-3]      [e per cell]                  [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',71)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        DO ik = 1,  nktotf
          DO ibnd = 1, nbndfst
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              DO j = 1, 3
                DO i = 1, 3
                  sigma(i, j, itemp) = sigma(i, j, itemp) - vkk_all(j, ibnd, ik) * F_out(i, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              Fi_check(:, itemp) = Fi_check(:, itemp) + F_out(:, ibnd, ik, itemp) * sfac / (nkf1 * nkf2 * nkf3) 
              ! 
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for electron conduction
            IF (etf_all(ibnd, ik) > ef0(itemp)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        !
        CALL prtmob(sigma(:, :, itemp), carrier_density, Fi_check(:, itemp), ef0(itemp)) 
        ! 
      ENDDO
      ! 
    ENDIF 
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_mob
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE prtmob(sigma, carrier_density, Fi_check, ef0) 
    !-----------------------------------------------------------------------
    !! 
    !! This routine print the mobility in a nice format and in proper units.
    !! 
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : omega
    USE constants_epw, ONLY : zero, kelvin2eV, ryd2ev, eps80, &
                              electron_SI, bohr2ang, ang2cm, hbarJ
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: sigma(3, 3)
    !! Conductivity tensor
    REAL(KIND = DP), INTENT(in) :: carrier_density
    !! Carrier density in cm^-3
    REAL(KIND = DP), INTENT(in) :: Fi_check(3)
    !! Integrated population vector
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi-level 
    ! 
    ! Local variables
    REAL(KIND = DP) :: mobility(3, 3)
    !! mobility 
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    ! 
    inv_cell = 1.0d0 / omega
    ! carrier_density in cm^-1
    carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
    IF (ABS(carrier_density) < eps80) CALL errore('print_mob_sym', 'The carrier density is 0', 1)
    ! 
    mobility(:, :) = (sigma(:, :) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
    ! 
    WRITE(stdout, '(5x, 1f8.3, 1f9.4, 1E14.5, 1E14.5, 3E16.6)') etemp * ryd2ev / kelvin2eV, ef0 * ryd2ev, &
          carrier_density, SUM(Fi_check(:)), mobility(1, 1), mobility(1, 2), mobility(1, 3)
    WRITE(stdout,'(50x, 3E16.6)') mobility(2, 1), mobility(2, 2), mobility(2, 3) 
    WRITE(stdout,'(50x, 3E16.6)') mobility(3, 1), mobility(3, 2), mobility(3, 3)
    ! 
    !-----------------------------------------------------------------------
    END SUBROUTINE prtmob
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE print_clock_epw
    !-----------------------------------------------------------------------
    !!  
    !! Adapted from the code PH/print_clock_ph - Quantum-ESPRESSO group          
    !!  
    USE io_global,     ONLY : stdout
    USE uspp,          ONLY : nlcc_any
    !  
    IMPLICIT NONE
    !
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x, a)') 'Unfolding on the coarse grid'
    CALL print_clock('dvanqq2')
    CALL print_clock('elphon_wrap')
    WRITE(stdout, '(5x)')
    CALL print_clock('ELPHWAN')
    WRITE(stdout, '(5x, a)') 'INITIALIZATION: '
    CALL print_clock('epq_init')
    WRITE(stdout, '(5x)')
    CALL print_clock('epq_init')
    IF (nlcc_any) CALL print_clock('set_drhoc')
    CALL print_clock('init_vloc')
    CALL print_clock('init_us_1')
    CALL print_clock('newd')
    CALL print_clock('dvanqq2')
    CALL print_clock('drho')
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x)')
    !
    ! Electron-Phonon interpolation 
    !
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x, a)') 'Electron-Phonon interpolation'   
    CALL print_clock('ephwann')
    CALL print_clock('ep-interp')
    CALL print_clock('PH SELF-ENERGY')
    CALL print_clock('ABS SPECTRA')
    CALL print_clock('crys_cart')
    WRITE(stdout, '(5x)')
    CALL print_clock('load data')
    CALL print_clock('Ham: step 1')
    CALL print_clock('Ham: step 2')
    CALL print_clock('Ham: step 3')
    CALL print_clock('Ham: step 4')
    CALL print_clock('ep: step 1')
    CALL print_clock('ep: step 2')
    CALL print_clock('ep: step 3')
    CALL print_clock('ep: step 4')
    CALL print_clock('DynW2B')
    CALL print_clock('HamW2B')
    CALL print_clock('ephW2Bp')
    CALL print_clock('ephW2B')
    CALL print_clock('print_ibte')
    CALL print_clock('vmewan2bloch')
    CALL print_clock('vmewan2blochp')
    CALL print_clock('ephW2Bp1')
    CALL print_clock('ephW2Bp2')
    CALL print_clock('kpoint_paral') 
    !
    ! Eliashberg equations
    WRITE(stdout, '(5x)')
    CALL print_clock('ELIASHBERG')
    !
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x, a)') 'Total program execution'
    CALL print_clock('EPW') 
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_clock_epw
    !-----------------------------------------------------------------------
    ! 
  !-------------------------------------------------------------------------
  END MODULE printing
  !-------------------------------------------------------------------------
