  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !                                                                            
  !----------------------------------------------------------------------
  MODULE io_transport
  !----------------------------------------------------------------------
  !! 
  !! This module contains various writing or reading routines related to transport. 
  !! Most of them are for restart purposes. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_ibte(iqq, iq, totq, ef0, efcb, first_cycle, ind_tot, ind_totcb, &
                          lrepmatw2_restart, lrepmatw5_restart, ctype) 
    !-----------------------------------------------------------------------
    !!
    !! This subroutine computes the transition probability and the scattering rates.
    !! Only the elements larger than threshold are saved on file. 
    !!
    USE kinds,         ONLY : DP, i4b
    USE cell_base,     ONLY : omega
    USE io_global,     ONLY : stdout
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, fsthick, eps_acustic, degaussw, nkf1, nkf2, nkf3,   &  
                              nstemp, scattering_serta, scattering_0rta, shortrange, nqf1, &
                              restart, restart_freq, restart_filq, vme, ncarrier, nqf2, nqf3
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, dmef, vmef, wf, wqf,      & 
                              epf17, nkqtotf, inv_tau_all, inv_tau_allcb, adapt_smearing, &
                              xqf, xkf, wkf, dmef, vmef, nqf, eta, transp_temp, lower_bnd,&
                              nbndfst, nktotf
    USE constants_epw, ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev, eps4, eps8, & 
                              eps6, eps10, bohr2ang, ang2cm
    USE io_files,      ONLY : prefix, diropn, tmp_dir
    USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,     ONLY : world_comm, my_pool_id, npool
    USE io_global,     ONLY : ionode_id
    USE io_var,        ONLY : iunepmat, iunepmatcb, iufilibtev_sup, iunrestart, iuntau,   &
                              iunsparseq, iunsparsek, iunsparsei, iunsparsej, iunsparset, &
                              iunsparseqcb, iunsparsekcb, iunsparseicb, iunsparsejcb,     &
                              iunsparsetcb, iuntaucb
    USE elph2,         ONLY : lrepmatw2_merge, lrepmatw5_merge, threshold
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, MPI_INTEGER8, &
                                 MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_INTEGER4, &
                                 MPI_MODE_CREATE, MPI_INFO_NULL, MPI_MODE_WRONLY, MPI_OFFSET
#endif
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart 
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in selecq
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! Current position inside the file during writing
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! Current position inside the file during writing (electron)
#if defined(__MPI)  
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of element written to file 
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of element written to file 
#else
    INTEGER, INTENT(inout) :: ind_tot
    !! Total number of element written to file 
    INTEGER, INTENT(inout) :: ind_totcb
    !! Total number of element written to file 
#endif   
    INTEGER, INTENT(in) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp. Could be unused (0).
    !
    ! Local variables
    CHARACTER(LEN = 256) :: filint
    !! Name of the file
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates  
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
    !! Index over temperature range
    INTEGER :: ifil
    !! Index over the files
    INTEGER :: nu, mu
    !! Index for modes
    INTEGER :: pbnd
    !! Index for bands
    INTEGER :: ierr
    !! Error
    INTEGER :: ipool
    !! Pool index
    INTEGER :: ind(npool)
    !! Nb of Matrix elements that are non-zero 
    INTEGER :: indcb(npool)
    !! Nb of Matrix elements that are non-zero in the cb
    INTEGER(KIND = i4b) :: sparse_q(nbndfst * nbndfst * nstemp * nkf)
    !! Index of q-points for mapping 
    INTEGER(KIND = i4b) :: sparse_k(nbndfst * nbndfst * nstemp * nkf)
    !! Index of k-points for mapping 
    INTEGER(KIND = i4b) :: sparse_i(nbndfst * nbndfst * nstemp * nkf)
    !! Index of i-bands for mapping
    INTEGER(KIND = i4b) :: sparse_j(nbndfst * nbndfst * nstemp * nkf)
    !! Index of j-bands for mapping
    INTEGER(KIND = i4b) :: sparse_t(nbndfst * nbndfst * nstemp * nkf)
    !! Index of temperature for mapping
    INTEGER(KIND = i4b) :: sparsecb_q(nbndfst * nbndfst * nstemp * nkf)
    !! Index of q-points for cb for mapping 
    INTEGER(KIND = i4b) :: sparsecb_k(nbndfst * nbndfst * nstemp * nkf)
    !! Index of k-points for cb for mapping 
    INTEGER(KIND = i4b) :: sparsecb_i(nbndfst * nbndfst * nstemp * nkf)
    !! Index of i-band for cb for mapping 
    INTEGER(KIND = i4b) :: sparsecb_j(nbndfst * nbndfst * nstemp * nkf)
    !! Index of j-band for cb for mapping 
    INTEGER(KIND = i4b) :: sparsecb_t(nbndfst * nbndfst * nstemp * nkf)
    !! Index of temeprature for cb for mapping 
    REAL(KIND = DP) :: tmp
    !! Temporary variable
    REAL(KIND = DP) :: tmp2 
    !! Temporary variable
    REAL(KIND = DP) :: dfnk
    !! Derivative of f_nk with respect to \varepsilon_nk
    REAL(KIND = DP) :: ekk2
    !! Temporary variable to the eigenenergies for the degenerate average  
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: ekq
    !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared (g2 is Ry^2) 
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: w0g1
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q}] $$ 
    REAL(KIND = DP) :: w0g2 
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q}] $$
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! Inverse phonon frequency. Defined for efficiency reasons.
    REAL(KIND = DP) :: inv_etemp
    !! Invese temperature inv_etemp = 1/etemp. Defined for efficiency reasons.
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! Used to set component to 0 if the phonon freq. is too low. This is defined
    !! for efficiency reasons as if statement should be avoided in inner-most loops.
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency $$\omega_{q\nu}$$ on the fine grid.  
    REAL(KIND = DP) :: wgq(nmodes)
    !! Bose-Einstein occupation function $$n_{q\nu}$$
    REAL(KIND = DP) :: weight
    !! Self-energy factor 
    REAL(KIND = DP) :: fmkq
    !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
    REAL(KIND = DP) :: vkk(3, nbndfst)
    !! Electronic velocity $$v_{n\mathbf{k}}$$
    REAL(KIND = DP) :: trans_prob(nbndfst * nbndfst * nstemp * nkf)
    !! Temporary array to store the scattering rates
    REAL(KIND = DP) :: trans_probcb(nbndfst * nbndfst * nstemp * nkf)
    !! Temporary array to store the scattering rates
    REAL(KIND = DP) :: wkf_all(nktotf)
    !! Weights from all the cores
    REAL(KIND = DP) :: vkk_all(3, nbndfst, nktotf)
    !! Velocities from all the cores
    REAL(KIND = DP) :: inv_eta(nmodes, nbndfst, nktotf)
    !! Inverse of the eta for speed purposes 
    REAL(KIND = DP) :: etf_all(nbndfst, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP) :: epf2_deg(nbndfst, nbndfst, nmodes)
    !! Epc in degeneracies
    REAL(KIND = DP) :: w_1
    !! Temporary electronic energy
    REAL(KIND = DP) :: w_2
    !! Temporary electronic energy
    REAL(KIND = DP) :: carrier_density
    !! Carrier concentration 
    REAL(KIND = DP) :: fnk 
    !! Fermi-Dirac
    REAL(KIND = DP) :: inv_cell
    !! cell volume
    REAL(KIND = DP) :: inv_tau_all_MPI(nbndfst, nktotf, nstemp)
    !! Auxiliary variables
    REAL(KIND = DP) :: inv_tau_allcb_MPI(nbndfst, nktotf, nstemp)
    !! Auxiliary variables
    REAL(KIND = DP), EXTERNAL :: DDOT
    !! Dot product function
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    !  
    inv_cell = 1.0d0 / omega
    ! 
    IF (iqq == 1) THEN
      !
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Scattering rate for IBTE")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Restart and restart_freq inputs deactivated (restart point at every q-points).")')
      WRITE(stdout, '(5x,"No intermediate mobility will be shown.")')
      !
      IF (fsthick < 1.d3) THEN
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a)' ) 'This is computed with respect to the fine Fermi level ',ef * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a,f10.6,a)' ) 'Only states between ', (ef - fsthick) * ryd2ev, ' eV and ', &
                (ef + fsthick) * ryd2ev, ' eV will be included'
        WRITE(stdout, '(5x,a/)')
      ENDIF
      ! 
      ! We save matrix elements larger than threshold defined in ephwann_shuffle
      WRITE(stdout,'(5x,a,1E20.12)') 'Save matrix elements larger than threshold: ', threshold
      WRITE(stdout,'(5x," ")')
      !
    ENDIF
    ! 
    ! In the case of a restart do not add the first step
    IF (first_cycle) THEN
      first_cycle = .FALSE.
    ELSE
      ! 
      ! To avoid if branching in the loop
      inv_eta(:, :, :) = zero 
      IF (adapt_smearing) THEN
        DO ik = 1, nkf
          DO ibnd = 1, nbndfst
            DO imode = 1, nmodes
              inv_eta(imode, ibnd, ik) = 1.0d0 / (SQRT(2d0) * eta(imode, ibnd, ik))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ik = 1, nkf
          DO ibnd = 1, nbndfst
            DO imode = 1, nmodes
              inv_eta(imode, ibnd, ik) = 1.0d0 / degaussw
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      ! 
      ! Average the el-ph matrix elements on degenerate bands and phonon modes. 
      ! This is important to ensure that the mobility tensor perfectly respects crystal symmetry. 
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! 
        wkf_all(ik + lower_bnd - 1) = wkf(ikk)
        !
        ! Average over the k electrons
        DO nu = 1, nmodes
          DO jbnd = 1, nbndfst
            DO ibnd = 1, nbndfst
              w_1 = etf(ibndmin - 1 + ibnd, ikk)
              g2  = zero
              n   = 0
              DO pbnd = 1, nbndfst
                w_2 = etf(ibndmin - 1 + pbnd, ikk)
                IF (ABS(w_2-w_1) < eps6) THEN
                  n = n + 1
                  g2 = g2 + ABS(epf17(jbnd, pbnd, nu, ik))**two
                ENDIF
              ENDDO
              epf2_deg(jbnd, ibnd, nu) = SQRT(g2 / FLOAT(n))
            ENDDO
          ENDDO
        ENDDO
        epf17(:, :, :, ik) = epf2_deg(:, :, :)
        !
        ! Average over the k+q electrons
        DO nu = 1, nmodes
          DO jbnd = 1, nbndfst
            DO ibnd = 1, nbndfst
              w_1 = etf(ibndmin - 1 + jbnd, ikq)
              g2 = 0.d0
              n  = 0
              DO pbnd = 1, nbndfst
                w_2 = etf(ibndmin - 1 + pbnd, ikq)
                IF (ABS(w_2 - w_1) < eps6) THEN
                  n = n + 1
                  g2 = g2 + ABS(epf17(pbnd, ibnd, nu, ik))**two
                ENDIF
              ENDDO
              epf2_deg(jbnd, ibnd, nu) = g2 / FLOAT(n) 
            ENDDO
          ENDDO
        ENDDO
        ! 
        ! Note that we already took the square above
        epf17(:, :, :, ik) = epf2_deg(:, :, :)
      ENDDO ! ik
      ! 
      trans_prob(:)    = zero
      sparse_q(:)      = zero
      sparse_k(:)      = zero
      sparse_i(:)      = zero
      sparse_j(:)      = zero
      sparse_t(:)      = zero
      trans_probcb(:)  = zero
      sparsecb_q(:)    = zero
      sparsecb_k(:)    = zero
      sparsecb_i(:)    = zero
      sparsecb_j(:)    = zero
      sparsecb_t(:)    = zero
      etf_all(:, :)    = zero
      vkk_all(:, :, :) = zero
      ind(:)           = 0
      indcb(:)         = 0
      ! loop over temperatures
      DO itemp = 1, nstemp
        !
        ! Define the inverse so that we can efficiently multiply instead of dividing
        etemp = transp_temp(itemp)
        inv_etemp = 1.0 / etemp
        !
        ! Now pre-treat phonon modes for efficiency for this specific current q-point.
        ! Treat phonon frequency and Bose occupation
        wq(:) = zero
        DO imode = 1, nmodes
          wq(imode) = wf(imode, iq)
          IF (wq(imode) > eps_acustic) THEN
            g2_tmp(imode) = 1.0d0
            wgq(imode)    = wgauss(-wq(imode) * inv_etemp, -99)
            wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
            inv_wq(imode) =  1.0d0 / (two * wq(imode))
          ELSE
            g2_tmp(imode) = 0.0
            wgq(imode)    = 0.0
            inv_wq(imode) = 0.0
          ENDIF
        ENDDO
        !  
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
            IF (ctype == 0 .OR. ctype == -1) THEN ! hole          
              DO ibnd = 1, nbndfst
                ! Energy at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                ! Fermi-Dirac occupation $f_{nk}$
                fnk = wgauss(-ekk * inv_etemp, -99)
                !
                ! This is to know if we need to store the data 
                ! derivative Fermi distribution -df_nk/dE_nk = (f_nk)*(1-f_nk)/ (k_B T) 
                dfnk = w0gauss(ekk * inv_etemp, -99 ) * inv_etemp
                !  
                DO jbnd = 1, nbndfst
                  !
                  ! Energy and fermi occupation at k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0(itemp)
                  ! Fermi-Dirac occupation $f_{mk+q}$ 
                  fmkq = wgauss(-ekq * inv_etemp, -99)
                  !
                  ! We perform a sum over the modes  
                  tmp  = zero
                  tmp2 = zero
                  DO imode = 1, nmodes
                    !
                    ! Here we take into account the zero-point SQRT(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    ! Note that epf17 has already been squared above during averaging. 
                    g2 = epf17(jbnd, ibnd, imode, ik) * inv_wq(imode) * g2_tmp(imode)
                    !
                    ! delta[E_k - E_k+q + w_q] 
                    w0g1 = w0gauss((ekk - ekq + wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                    ! delta[E_k - E_k+q - w_q]
                    w0g2 = w0gauss((ekk - ekq - wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                    !
                    ! Transition probability - See Eq. 41 of arXiv:1908.01733 (2019).
                    ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                    ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                    !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                    !
                    ! This is summed over modes
                    tmp = tmp + two * pi * wqf(iq) * g2 * ((fnk + wgq(imode)) * w0g2 + (one - fnk + wgq(imode)) * w0g1)
                    ! 
                    ! Here we compute the scattering rate - Eq. 64 of arXiv:1908.01733 (2019).
                    ! inv_tau = (2 pi/hbar) * (k+q-point weight) * g2 * 
                    ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] }                  
                    tmp2 = tmp2 + two * pi * wqf(iq) * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                    !
                  ENDDO !imode
                  inv_tau_all(ibnd, ik + lower_bnd - 1, itemp) = inv_tau_all(ibnd, ik + lower_bnd - 1, itemp) + tmp2 
                  !  
                  ! Only save the elements that really contribute
                  ! The check is made on the SERTA mobility - See Eq. 44 of arXiv:1908.01733 (2019). 
                  IF (ABS(tmp2 * dfnk) > threshold) THEN
                    ! 
                    ind(my_pool_id + 1) = ind(my_pool_id + 1) + 1 
                    trans_prob(ind(my_pool_id + 1)) = tmp
                    sparse_q(ind(my_pool_id + 1)) = iq
                    sparse_k(ind(my_pool_id + 1)) = ik + lower_bnd - 1
                    sparse_i(ind(my_pool_id + 1)) = ibnd
                    sparse_j(ind(my_pool_id + 1)) = jbnd 
                    sparse_t(ind(my_pool_id + 1)) = itemp
                    !  
                  ENDIF 
                ENDDO !jbnd
              ENDDO ! ibnd
            ENDIF ! ctype
            !
            ! In this case we are also computing the scattering rate for another Fermi level position
            ! This is used to compute both the electron and hole mobility at the same time.  
            IF (ctype == 0 .OR. ctype == 1) THEN
              ! 
              DO ibnd = 1, nbndfst
                ! Energy at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - efcb(itemp) 
                ! Fermi-Diract distribution $f_{nk}$
                fnk = wgauss(-ekk * inv_etemp, -99)   
                ! 
                ! Derivative Fermi distribution (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                dfnk = w0gauss(ekk * inv_etemp, -99) * inv_etemp
                !
                DO jbnd = 1, nbndfst
                  !
                  !  energy and fermi occupation at k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - efcb(itemp)
                  fmkq = wgauss(-ekq * inv_etemp, -99)
                  ! 
                  tmp  = zero
                  tmp2 = zero
                  DO imode = 1, nmodes
                    ! Same as above but for conduction bands (electrons)
                    g2 = epf17(jbnd, ibnd, imode, ik) * inv_wq(imode) * g2_tmp(imode)
                    w0g1 = w0gauss((ekk - ekq + wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                    w0g2 = w0gauss((ekk - ekq - wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                    tmp = tmp  + two * pi * wqf(iq) * g2 * ((fnk + wgq(imode)) * w0g2 + (one - fnk + wgq(imode)) * w0g1)
                    tmp2 = tmp2 + two * pi * wqf(iq) * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                  ENDDO ! imode
                  inv_tau_allcb(ibnd, ik + lower_bnd - 1, itemp) = inv_tau_allcb(ibnd, ik + lower_bnd - 1, itemp) + tmp2
                  ! 
                  IF (ABS(tmp2 * dfnk) > threshold) THEN
                    indcb (my_pool_id + 1) = indcb(my_pool_id + 1) + 1
                    trans_probcb(indcb(my_pool_id + 1)) = tmp
                    sparsecb_q(indcb(my_pool_id + 1)) = iq
                    sparsecb_k(indcb(my_pool_id + 1)) = ik + lower_bnd - 1
                    sparsecb_i(indcb(my_pool_id + 1)) = ibnd
                    sparsecb_j(indcb(my_pool_id + 1)) = jbnd
                    sparsecb_t(indcb(my_pool_id + 1)) = itemp
                  ENDIF
                ENDDO !jbnd
              ENDDO !ibnd
            ENDIF ! ctype
          ENDIF ! endif fsthick
        ENDDO ! end loop on k
      ENDDO ! itemp
      ! If the q-point is taken, write on file
      CALL mp_sum(ind, world_comm)
      CALL mp_sum(indcb, world_comm)
      ! 
      ! SP - IBTE only with if EPW compiled with MPI
      IF (SUM(ind) > 0) THEN
        ! 
        IF (my_pool_id == 0 ) ind_tot = ind_tot + SUM(ind)
#if defined(__MPI)
        CALL MPI_BCAST(ind_tot, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
#endif
!       WRITE(stdout,'(a,i9,E22.8)') '     Total number of element written ',ind_tot
        IF (ind(my_pool_id + 1) > 0) THEN 
          WRITE(iunepmat) trans_prob(1:ind(my_pool_id + 1))
          ! The flush is crucial otherwise restart wont work correctly. 
          FLUSH(iunepmat)
          DO ifil = 1, ind(my_pool_id + 1)
            WRITE(iunsparseq) sparse_q(ifil)
            WRITE(iunsparseq) sparse_k(ifil)
            WRITE(iunsparseq) sparse_i(ifil)
            WRITE(iunsparseq) sparse_j(ifil)
            WRITE(iunsparseq) sparse_t(ifil)
          ENDDO
          FLUSH(iunsparseq)
        ENDIF
        ! 
        ! Offset for the next q iteration
        lrepmatw2_merge = lrepmatw2_merge + ind(my_pool_id + 1)
      ENDIF 
      IF (SUM(indcb) > 0) THEN
        ! 
        IF (my_pool_id == 0 ) ind_totcb = ind_totcb + SUM(indcb)
#if defined(__MPI)
        CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
#endif
        ! 
        IF (indcb(my_pool_id + 1) > 0) THEN
          WRITE(iunepmatcb) trans_probcb(1:indcb(my_pool_id + 1))
          FLUSH(iunepmatcb)
          DO ifil = 1, indcb(my_pool_id + 1)
            WRITE(iunsparseqcb) sparsecb_q(ifil)
            WRITE(iunsparseqcb) sparsecb_k(ifil)
            WRITE(iunsparseqcb) sparsecb_i(ifil)
            WRITE(iunsparseqcb) sparsecb_j(ifil)
            WRITE(iunsparseqcb) sparsecb_t(ifil)
          ENDDO
          FLUSH(iunsparseqcb)
        ENDIF
        !
        ! Offset for the next q iteration
        lrepmatw5_merge = lrepmatw5_merge + indcb(my_pool_id + 1)
        ! 
      ENDIF ! indcb
      ! 
      ! Save to file restart information in formatted way for possible restart
      lrepmatw2_restart(:) = 0
      lrepmatw5_restart(:) = 0
      lrepmatw2_restart(my_pool_id + 1) = lrepmatw2_merge 
      lrepmatw5_restart(my_pool_id + 1) = lrepmatw5_merge 
      call mp_sum(lrepmatw2_restart, world_comm)
      call mp_sum(lrepmatw5_restart, world_comm)
      !
      inv_tau_all_MPI = inv_tau_all
      inv_tau_allcb_MPI = inv_tau_allcb
      CALL mp_sum(inv_tau_all_MPI, world_comm)
      CALL mp_sum(inv_tau_allcb_MPI, world_comm)
      !
      IF (my_pool_id == 0) THEN
        OPEN(UNIT = iunrestart, FILE = 'restart_ibte.fmt')
        WRITE(iunrestart, *) iqq
        WRITE(iunrestart, *) ind_tot
        WRITE(iunrestart, *) ind_totcb
        WRITE(iunrestart, *) npool
        DO ipool = 1, npool
          WRITE(iunrestart, *) lrepmatw2_restart(ipool)
        ENDDO
        DO ipool = 1, npool 
          WRITE(iunrestart, *) lrepmatw5_restart(ipool)
        ENDDO
        CLOSE(iunrestart)
        ! 
        OPEN(UNIT = iuntau, FORM = 'unformatted', FILE = 'inv_tau_tmp')
        WRITE(iuntau) inv_tau_all_MPI
        CLOSE(iuntau)
        ! 
        OPEN(UNIT = iuntaucb, FORM = 'unformatted', FILE = 'inv_taucb_tmp')
        WRITE(iuntaucb) inv_tau_allcb_MPI
        CLOSE(iuntaucb)
      ENDIF    
      ! 
    ENDIF ! first_cycle
    ! 
    IF (iqq == totq) THEN
      wkf_all(:) = zero
      ! Computes the k-velocity
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ! 
        wkf_all(ik + lower_bnd -1 ) = wkf(ikk) 
        ! 
        DO ibnd = 1, nbndfst
          IF (vme) THEN
            vkk_all(:, ibnd, ik + lower_bnd - 1) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
          ELSE
            vkk_all(:,ibnd, ik + lower_bnd -1 ) = 2.0 * REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
          ENDIF  
          etf_all(ibnd, ik + lower_bnd - 1) = etf(ibndmin - 1 + ibnd, ikk)
        ENDDO
      ENDDO 
      CALL mp_sum(vkk_all, world_comm) 
      CALL mp_sum(etf_all, world_comm) 
      CALL mp_sum(wkf_all, world_comm)
      CALL mp_sum(inv_tau_all, world_comm)
      CALL mp_sum(inv_tau_allcb, world_comm)
      ! 
      IF (my_pool_id == 0) THEN
        ! Now write total number of q-point inside and k-velocity
        OPEN(iufilibtev_sup, FILE = 'IBTEvel_sup.fmt', FORM = 'formatted')
        WRITE(iufilibtev_sup, '(a)') '# Number of elements in hole and electrons  '
        WRITE(iufilibtev_sup, '(2i16)') ind_tot, ind_totcb
        WRITE(iufilibtev_sup, '(a)') '# itemp    ef0    efcb'
        DO itemp = 1, nstemp
          WRITE(iufilibtev_sup, '(i8,2E22.12)') itemp, ef0(itemp), efcb(itemp)
        ENDDO
        WRITE(iufilibtev_sup, '(a)') '# ik  ibnd      velocity (x,y,z)              eig     weight '
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            WRITE(iufilibtev_sup, '(i8,i6,5E22.12)') ik, ibnd, vkk_all(:, ibnd, ik), etf_all(ibnd, ik), wkf_all(ik)
          ENDDO
        ENDDO
        CLOSE(iufilibtev_sup)
        ! 
        ! Save the inv_tau and inv_tau_all on file (formatted)
        OPEN(iufilibtev_sup, FILE = 'inv_tau.fmt', FORM = 'formatted')
        WRITE(iufilibtev_sup, '(a)') '# Hole relaxation time  '
        WRITE(iufilibtev_sup, '(a)') '# itemp    kpt      ibnd    energy [Ry]   relaxation time [?]'
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              WRITE(iufilibtev_sup, '(i5,i8,i6,2E22.12)') itemp, ik, ibnd, etf_all(ibnd, ik), inv_tau_all(ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilibtev_sup)       
        ! 
        ! Save the inv_tau and inv_tau_all on file (formatted)
        OPEN(iufilibtev_sup, FILE = 'inv_taucb.fmt', FORM = 'formatted')
        WRITE(iufilibtev_sup, '(a)') '# Hole relaxation time  '
        WRITE(iufilibtev_sup, '(a)') '# itemp    kpt      ibnd    energy [Ry]   relaxation time [?]'
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              WRITE(iufilibtev_sup, '(i5,i8,i6,2E22.12)') itemp, ik, ibnd, etf_all(ibnd, ik), inv_tau_allcb(ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilibtev_sup)
        ! 
      ENDIF ! master
      ! 
      ! Now print the carrier density for checking
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        IF (ncarrier < 0.0) THEN ! VB
          DO ik = 1, nkf
            DO ibnd = 1, nbndfst
              ! This selects only valence bands for hole conduction
              IF (etf_all(ibnd, ik + lower_bnd - 1 ) < ef0(itemp)) THEN
                ! Energy at k (relative to Ef)
                ekk = etf_all(ibnd, ik + lower_bnd - 1 ) - ef0(itemp)
                fnk = wgauss(-ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf_all(ik + lower_bnd - 1) * (1.0d0 - fnk)
              ENDIF
            ENDDO
          ENDDO
          CALL mp_sum(carrier_density, world_comm)
          carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6)') etemp * ryd2ev / kelvin2eV, ef0(itemp) * ryd2ev,  carrier_density
        ELSE ! CB
          DO ik = 1, nkf
            DO ibnd = 1, nbndfst
              ! This selects only valence bands for hole conduction
              IF (etf_all (ibnd, ik + lower_bnd - 1) > efcb(itemp)) THEN
                !  energy at k (relative to Ef)
                ekk = etf_all(ibnd, ik + lower_bnd - 1) - efcb(itemp)
                fnk = wgauss(-ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf_all(ik + lower_bnd - 1) *  fnk
              ENDIF
            ENDDO
          ENDDO
          CALL mp_sum(carrier_density, world_comm)
          carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6)') etemp * ryd2ev / kelvin2eV, efcb(itemp) * ryd2ev,  carrier_density
        ENDIF ! ncarrier
      ENDDO
    ENDIF ! iqq
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE print_ibte
    !-----------------------------------------------------------------------
    !  
    !----------------------------------------------------------------------------
    SUBROUTINE fin_write(iter, F_in, av_mob_old, elec)
    !----------------------------------------------------------------------------
    !!
    !! Writes the F without magnetic field for restart
    !! 
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilFi_all
    USE io_files,      ONLY : diropn
    USE epwcom,        ONLY : nstemp
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE elph2,         ONLY : ibndmax, ibndmin, nkqtotf, nktotf, nbndfst
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iter
    !! Iteration number
    REAL(KIND = DP), INTENT(in) :: f_in(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i  
    REAL(KIND = DP), INTENT(in) :: av_mob_old(nstemp)
    !! Error in the hole mobility
    LOGICAL, INTENT(in) :: elec
    !! IF true we do electron mobility, if false the hole one. 
    ! 
    ! Local variable
    LOGICAL :: exst
    !! File exist
    INTEGER :: i
    !! Running index for the vector
    INTEGER :: lfi_all
    !! Length of the vector
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: idir
    !! Direction index
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
    exst = .FALSE.
    aux(:) = zero
    IF (mpime == ionode_id) THEN
      !
      lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
      ! First element is the iteration number
      aux(1) = iter
      ! 
      i = 1
      DO itemp = 1, nstemp
        i = i + 1  
        ! Value of the previous h mobility (used for error evaluation)
        aux(i) = av_mob_old(itemp)
      ENDDO
      ! 
      i = 1 + nstemp 
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO idir = 1, 3
              i = i +1
              aux(i) = f_in(idir, ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! 
      ! Electron mobility
      IF (elec) THEN
        CALL diropn(iufilFi_all, 'Fin_restartcb', lfi_all, exst)
        CALL davcio(aux, lfi_all, iufilFi_all, 1, +1)
      ELSE 
        CALL diropn(iufilFi_all, 'Fin_restart', lfi_all, exst)
        CALL davcio( aux, lfi_all, iufilFi_all, 1, +1)
      ENDIF
      CLOSE(iufilFi_all)
      !
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE fin_write
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE fin_read(iter, F_in, av_mob_old, elec)
    !----------------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilFi_all
    USE epwcom,    ONLY : nstemp
    USE constants_epw, ONLY : zero
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, nbndfst, nktotf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iter
    !! Iteration number
    REAL(KIND = DP), INTENT(inout) :: f_in(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i  
    REAL(KIND = DP), INTENT(inout) :: av_mob_old(nstemp)
    !! Error in the hole mobility
    LOGICAL, INTENT(in) :: elec
    !! IF true we do electron mobility, if false the hole one. 
    !
    ! Local variable
    CHARACTER(LEN = 256) :: name1
    !! Variable name
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Running index for the vector
    INTEGER :: lfi_all
    !! Length of the vector
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: idir
    !! Direction index
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      ! 
      ! First inquire if the file exists
      IF (elec) THEN
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restartcb1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restartcb'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        ! 
        IF (exst) THEN ! read the file
          !
          lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
          CALL diropn(iufilFi_all, 'Fin_restartcb', lfi_all, exst)
          CALL davcio(aux, lfi_all, iufilFi_all, 1, -1)
          !
          ! First element is the iteration number
          iter = INT(aux(1)) 
          !
          i = 1
          DO itemp = 1, nstemp
            i = i + 1  
            ! Last value of hole mobility 
            av_mob_old(itemp) = aux(i)
          ENDDO
          ! 
          i = 1 + nstemp 
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                DO idir = 1, 3
                  i = i + 1
                  f_in(idir, ibnd, ik, itemp) = aux(i)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilFi_all)
        ENDIF
      ELSE ! hole
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restart1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restart'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        ! 
        IF (exst) THEN ! read the file
          !
          lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
          CALL diropn(iufilFi_all, 'Fin_restart', lfi_all, exst)
          CALL davcio(aux, lfi_all, iufilFi_all, 1, -1)
          !
          ! First element is the iteration number
          iter = INT(aux(1))
          !
          i = 1
          DO itemp = 1, nstemp
            i = i + 1
            ! Last value of hole mobility 
            av_mob_old(itemp) = aux(i)
          ENDDO
          ! 
          i = 1 + nstemp
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, (nbndfst)
                DO idir = 1, 3
                  i = i + 1
                  f_in(idir, ibnd, ik, itemp) = aux(i)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilFi_all)
        ENDIF
      ENDIF
    ENDIF ! mpime
    ! 
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iter,       ionode_id, world_comm)
      CALL mp_bcast(F_in,       ionode_id, world_comm)
      CALL mp_bcast(av_mob_old, ionode_id, world_comm)
      WRITE(stdout, '(a,i10)' ) '     Restart from iter: ', iter
    ENDIF ! exists
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE fin_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_merge_parallel()
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_var,           ONLY : iunepmat_merge, iunepmat, iunepmatcb_merge,              &
                                 iunepmatcb, iunsparseq_merge, iunsparsek_merge,          &
                                 iunsparsej_merge,iunsparset_merge, iunepmatcb_merge,     &
                                 iunsparseqcb_merge, iunsparsekcb_merge, iunsparsei_merge,&
                                 iunsparseicb_merge, iunsparsejcb_merge, iunsparsetcb_merge
    USE mp_global,        ONLY : my_pool_id, npool, world_comm 
    USE io_files,         ONLY : tmp_dir, prefix 
    USE mp,               ONLY : mp_sum, mp_barrier
    USE io_global,        ONLY : stdout
    USE elph2,            ONLY : lrepmatw2_merge, lrepmatw5_merge
    USE epwcom,           ONLY : int_mob, carrier, ncarrier
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE,MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION,          &
                                 MPI_STATUS_IGNORE, MPI_INTEGER
#endif
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! Pool number, character
    CHARACTER(LEN = 256) :: dirname(2)
    !! Name of the directory to hold files
    CHARACTER(LEN = 256) :: filename(6)
    !! Name of the files to merge files
    CHARACTER(LEN = 256) :: path_to_files(2)
    !! Name of the path to files
    INTEGER :: i2, i4, i5, i6
    !! Indexes to loop over file sizes
    INTEGER :: ipool
    !! Process index
    INTEGER :: lrepmatw2_tot(npool)
    !! Lenght of each file
    INTEGER :: lrepmatw5_tot(npool)
    !! Lenght of each file
    INTEGER :: ich
    !! Loop over directories
    INTEGER :: ifil
    !! Index over the files
    INTEGER :: io_u(6)
    !! Input output units
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: sparse(:, :)
    !! Vaariable for reading and writing the files
    INTEGER, ALLOCATABLE :: sparsecb(:, :)
    !! Vaariable for reading and writing the files
    REAL(KIND = DP), ALLOCATABLE :: trans_prob(:)
    !! Variable for reading and writing trans_prob
    REAL(KIND = DP), ALLOCATABLE :: trans_probcb(:)
    !! Variable for reading and writing trans_prob
#if defined(__MPI)
    INTEGER (KIND = MPI_OFFSET_KIND) :: lsize
    !! Size of what we write
    INTEGER (KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset while writing scattering to files
    !
    IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 0.0))) THEN
      !
      ALLOCATE(trans_prob(lrepmatw2_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating trans_prob', 1)
      ALLOCATE(sparse(5, lrepmatw2_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating sparse', 1)
      !
      io_u(1) = iunepmat_merge
      io_u(2) = iunsparseq_merge
      io_u(3) = iunsparsek_merge
      io_u(4) = iunsparsei_merge
      io_u(5) = iunsparsej_merge
      io_u(6) = iunsparset_merge
      !
      dirname(1) = 'Fepmatkq1'
      dirname(2) = 'Fsparse'
      !
      filename(1) = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkq1' 
      filename(2) = 'sparseq'
      filename(3) = 'sparsek'
      filename(4) = 'sparsei'
      filename(5) = 'sparsej'
      filename(6) = 'sparset'
      !
      path_to_files(1) = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkq1' // '_'
      path_to_files(2) = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_'
      !
      lrepmatw2_tot = 0
      lrepmatw2_tot(my_pool_id + 1) = lrepmatw2_merge
      CALL mp_sum(lrepmatw2_tot, world_comm)
      DO ich = 1, 6
        CALL mp_barrier(world_comm)
        CALL MPI_FILE_OPEN(world_comm, filename(ich), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,io_u(ich), ierr)
      ENDDO
      !
      DO ich = 1, 2 
        ! Read files per processor
        WRITE(my_pool_id_ch, "(I0)") my_pool_id
        filint = TRIM(path_to_files(ich)) // TRIM(my_pool_id_ch)
        OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'old', FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
        IF (ich == 1) THEN
          DO i2 = 1, lrepmatw2_merge
            READ(iunepmat) trans_prob(i2)
          ENDDO
        ELSE
          DO i2 = 1, lrepmatw2_merge
            DO ifil = 1, 5
              READ(iunepmat) sparse(ifil, i2)
            ENDDO
          ENDDO
        ENDIF
        CLOSE(iunepmat, STATUS = 'delete')
        IF (ich == 1) THEN
          lrepmatw = INT(SUM(lrepmatw2_tot(1:my_pool_id + 1)) - lrepmatw2_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
          & 8_MPI_OFFSET_KIND 
          lsize = INT(lrepmatw2_merge, KIND = MPI_OFFSET_KIND) 
          CALL MPI_FILE_WRITE_AT(io_u(1), lrepmatw, trans_prob(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        ELSE
          DO ifil = 1, 5
            lrepmatw = INT(SUM(lrepmatw2_tot(1:my_pool_id + 1)) - lrepmatw2_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
            & 4_MPI_OFFSET_KIND 
            lsize = INT(lrepmatw2_merge, KIND = MPI_OFFSET_KIND) 
            CALL MPI_FILE_WRITE_AT(io_u(ifil+1), lrepmatw, sparse(ifil,:), lsize, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          ENDDO
        ENDIF
      ENDDO
      !
      DO ich = 1, 6
        CALL MPI_FILE_CLOSE(io_u(ich), ierr)
      ENDDO
      !
      DEALLOCATE(trans_prob, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating trans_prob', 1)
      DEALLOCATE(sparse, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating sparse', 1)
      !
    ENDIF
    IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 0.0))) THEN
      !
      ALLOCATE(trans_probcb(lrepmatw5_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating trans_probcb', 1)
      ALLOCATE(sparsecb(5, lrepmatw5_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating sparsecb', 1)
      !
      io_u(1) = iunepmatcb_merge
      io_u(2) = iunsparseqcb_merge
      io_u(3) = iunsparsekcb_merge
      io_u(4) = iunsparseicb_merge
      io_u(5) = iunsparsejcb_merge
      io_u(6) = iunsparsetcb_merge
      !
      dirname(1) = 'Fepmatkqcb1'
      dirname(2) = 'Fsparsecb'
      !
      filename(1) = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkqcb1' 
      filename(2) = 'sparseqcb'
      filename(3) = 'sparsekcb'
      filename(4) = 'sparseicb'
      filename(5) = 'sparsejcb'
      filename(6) = 'sparsetcb'
      !
      path_to_files(1) = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkqcb1' // '_'
      path_to_files(2) = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparsecb' // '_'
      !
      lrepmatw5_tot = 0
      lrepmatw5_tot(my_pool_id + 1) = lrepmatw5_merge
      CALL mp_sum(lrepmatw5_tot, world_comm)
      DO ich = 1, 6
        CALL mp_barrier(world_comm)
        CALL MPI_FILE_OPEN(world_comm, filename(ich), MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, io_u(ich), ierr)
      ENDDO
      !
      DO ich = 1, 2
        ! Read files per processor
        WRITE(my_pool_id_ch, "(I0)") my_pool_id
        filint = TRIM(path_to_files(ich)) // TRIM(my_pool_id_ch)
        OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
        IF (ich == 1) THEN
          DO i2 = 1, lrepmatw5_merge
            READ(iunepmatcb) trans_probcb(i2)
          ENDDO
        ELSE
          DO i2 = 1, lrepmatw5_merge
            DO ifil = 1 ,5
              READ(iunepmatcb) sparsecb(ifil, i2)
            ENDDO
          ENDDO
        ENDIF
        CLOSE(iunepmatcb, STATUS = 'delete')
        IF (ich == 1) THEN
          lrepmatw = INT(SUM(lrepmatw5_tot(1:my_pool_id + 1)) - lrepmatw5_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
          & 8_MPI_OFFSET_KIND 
          lsize = INT(lrepmatw5_merge, KIND = MPI_OFFSET_KIND) 
          CALL MPI_FILE_WRITE_AT(io_u(1), lrepmatw, trans_probcb(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        ELSE
          DO ifil = 1, 5
            lrepmatw = INT(SUM(lrepmatw5_tot(1:my_pool_id + 1)) - lrepmatw5_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
            & 4_MPI_OFFSET_KIND 
            lsize = INT(lrepmatw5_merge, KIND = MPI_OFFSET_KIND) 
            CALL MPI_FILE_WRITE_AT(io_u(ifil + 1), lrepmatw, sparsecb(ifil, :), lsize, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          ENDDO
        ENDIF
      ENDDO
      !
      DO ich = 1, 6
        CALL MPI_FILE_CLOSE(io_u(ich), ierr)
      ENDDO
      ! 
      DEALLOCATE(trans_probcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating trans_probcb', 1)
      DEALLOCATE(sparsecb, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating sparsecb', 1)
      !
    ENDIF ! in all other cases it is still to decide which files to open
    ! 
#endif
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_merge_parallel
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_open(ind_tot, ind_totcb, lrepmatw2_restart, lrepmatw5_restart)
    !----------------------------------------------------------------------------
    ! 
    ! This routine opens all the files needed to save scattering rates for the IBTE.
    ! 
    USE kinds,            ONLY : DP
    USE io_files,         ONLY : tmp_dir, prefix, create_directory, delete_if_present
    USE io_var,           ONLY : iunepmat, iunsparseq, iunsparsek,                 &
                                 iunsparsei, iunsparsej, iunsparset, iunsparseqcb, &
                                 iunsparsekcb, iunsparseicb, iunsparsejcb,         &
                                 iunsparsetcb, iunepmatcb, iunrestart
    USE mp_global,        ONLY : world_comm, my_pool_id, npool
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE elph2,            ONLY : lrepmatw2_merge, lrepmatw5_merge
    USE epwcom,           ONLY : int_mob, carrier, ncarrier
    USE io_global,        ONLY : ionode_id
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND
#endif
    ! 
    IMPLICIT NONE
    !  
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! To restart opening files
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! To restart opening files
#if defined(__MPI)
    INTEGER (KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER (KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of component for the conduction band
#else
    INTEGER, INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER, INTENT(inout) :: ind_totcb
    !! Total number of component for conduction band
#endif  
    ! 
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! my_pool_id in character type
    CHARACTER(LEN = 256) :: dirname(2), dirnamecb(2)
    !! Name of the directory to hold files
    LOGICAL :: exst
    !! Logical for existence of files
    LOGICAL :: exst2
    !! Logical for existence of files
    INTEGER :: ierr
    !! Error index
    INTEGER :: ilrep
    !! index to loop over the reading elements
    INTEGER :: dummy_int
    !! Dummy INTEGER for reading
    INTEGER :: ipool
    !! Pool index
    INTEGER (KIND = 8) :: position_byte
    !! Position in the file in byte
    REAL (KIND = DP) :: dummy_real
    !! Dummy variable for reading
    !
    WRITE(my_pool_id_ch, "(I0)") my_pool_id
    !
    dirname(1)   = 'Fepmatkq1'
    dirname(2)   = 'Fsparse'
    dirnamecb(1) = 'Fepmatkqcb1'
    dirnamecb(2) = 'Fsparsecb'
    ! 
    INQUIRE(FILE = 'restart_ibte.fmt', EXIST = exst)
    !
    IF (my_pool_id == ionode_id) THEN
      IF (exst) THEN
        OPEN(UNIT = iunrestart, FILE = 'restart_ibte.fmt', STATUS = 'old')
        READ (iunrestart,*) 
        READ (iunrestart,*) 
        READ (iunrestart,*) 
        READ (iunrestart,*) 
        DO ipool = 1, npool
          READ (iunrestart,*) lrepmatw2_restart(ipool)
        ENDDO
        DO ipool = 1, npool
          READ (iunrestart,*) lrepmatw5_restart(ipool)
        ENDDO
        CLOSE(iunrestart)
      ENDIF
    ENDIF
    CALL mp_bcast(exst, ionode_id, world_comm )
    CALL mp_bcast(lrepmatw2_restart, ionode_id, world_comm )
    CALL mp_bcast(lrepmatw5_restart, ionode_id, world_comm )
    !
    ! The restart_ibte.fmt exist - we try to restart
    IF (exst) THEN
      ! Hole
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 1E5))) THEN
        !
        filint = './' // ADJUSTL(TRIM(dirname(1))) // '/'//TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw2_restart(my_pool_id + 1) - 1) * 8 + 1  
            READ(iunepmat, POS=position_byte) dummy_real 
          ENDIF
        ELSE 
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fepmatkq1 folder', 1) 
        ENDIF
        ! 
        filint = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN 
            position_byte = (5 * lrepmatw2_restart(my_pool_id + 1) - 1) * 4 + 1
            READ(iunsparseq, POS = position_byte) dummy_int
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fsparse folder', 1) 
        ENDIF
        !
      ENDIF ! Hole
      ! Electron
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 1E5))) THEN
        !
        filint = './' // ADJUSTL(TRIM(dirnamecb(1))) // '/' // TRIM(prefix) // '.epmatkqcb1'//'_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw5_restart(my_pool_id + 1) - 1) * 8 + 1  
            READ(iunepmatcb, POS = position_byte) dummy_real
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fepmatkqcb1 folder', 1)
        ENDIF
        ! 
        filint = './' // ADJUSTL(TRIM(dirnamecb(2))) // '/' // 'sparsecb' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (5 * lrepmatw5_restart(my_pool_id + 1) - 1) * 4 + 1 
            READ(iunsparseqcb, POS = position_byte) dummy_int
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fsparse folder', 1)
        ENDIF
        !
      ENDIF ! electron
      lrepmatw2_merge = lrepmatw2_restart(my_pool_id + 1)
      lrepmatw5_merge = lrepmatw5_restart(my_pool_id + 1)
      !  
    ELSE ! no restart file present
      ! Hole
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 1E5))) THEN
        ! 
        CALL create_directory(ADJUSTL(TRIM(dirname(1))))
        CALL create_directory(ADJUSTL(TRIM(dirname(2))))
        ! 
        filint = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF 
        ! 
        filint = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF 
        ! 
      ENDIF ! Hole
      ! Electron
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 1E5))) THEN
        ! 
        CALL create_directory(ADJUSTL(TRIM(dirnamecb(1))))
        CALL create_directory(ADJUSTL(TRIM(dirnamecb(2))))        
        ! 
        filint = './' // ADJUSTL(TRIM(dirnamecb(1))) // '/' // TRIM(prefix) // '.epmatkqcb1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !  
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        ! 
        filint = './' // ADJUSTL(TRIM(dirnamecb(2))) // '/' // 'sparsecb' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        !  
      ENDIF !electron 
      lrepmatw2_merge = 0
      lrepmatw5_merge = 0
      ! 
    ENDIF ! restart 
    !
    ind_tot   = 0
    ind_totcb = 0
    lrepmatw2_restart(:) = 0
    lrepmatw5_restart(:) = 0
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_open
    !----------------------------------------------------------------------------    
    !
    !----------------------------------------------------------------------------
    SUBROUTINE scattering_write(itemp, etemp, ef0, etf_all)
    !----------------------------------------------------------------------------
    !! 
    !! Write scattering rates
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilscatt_rate
    USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, inv_tau_all, nbndfst, nktotf
    USE epwcom,    ONLY : nbndsub, nstemp
    USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                              meV2invps, eps4
    USE mp,        ONLY : mp_barrier
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndsub, nkqtotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file.
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: temp
    !! Temporary file name used to write scattering rate to file. 
    !
    WRITE(stdout, '(/5x,"Writing scattering rate to file"/)')
    !
    IF (mpime == ionode_id) THEN
      !
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1,'(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1,'(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1,'(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate, FILE = name1, FORM = 'formatted')
      WRITE(iufilscatt_rate, '(a)') '# Inverse scattering time (ps)'
      WRITE(iufilscatt_rate, '(a)') '#      ik       ibnd                 E(ibnd)    scattering rate(1/ps)'
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        DO ibnd = 1, nbndfst
          !
          ! note that ekk does not depend on q
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
          !
          WRITE(iufilscatt_rate, '(i9,2x)', ADVANCE = 'no') ik
          WRITE(iufilscatt_rate, '(i9,2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
          WRITE(iufilscatt_rate, '(E22.14)', ADVANCE = 'no') ryd2ev * ekk
          WRITE(iufilscatt_rate, '(E26.16E3)') ryd2mev * meV2invps * inv_tau_all(itemp, ibnd, ik)
          !
        ENDDO
        !
      ENDDO
      !
      CLOSE(iufilscatt_rate)
    ENDIF
    !CALL mp_barrier(inter_pool_comm)
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE scattering_write
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE scattering_read(etemp, ef0, etf_all, inv_tau_all)
    !----------------------------------------------------------------------------
    !!
    !! Read scattering files
    !!  
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilscatt_rate
    USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, nktotf, nbndfst
    USE epwcom,    ONLY : nbndsub, nstemp
    USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                              meV2invps, eps4
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(out) :: etf_all(nbndsub, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), INTENT(out) :: inv_tau_all(nstemp, nbndfst, nktotf)
    !! Inverse scattering rates
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file. 
    CHARACTER(LEN = 256) :: dummy1
    !! Dummy variable to store the text of the scattering_rate file 
    INTEGER :: ik
    !! K-point index
    INTEGER :: ik_tmp
    !! K-point index read from file
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ibnd_tmp
    !! Local band index read from file
    INTEGER :: ios
    !! Status of reading file
    REAL(KIND = DP) :: temp
    !! Temporary file name used to write scattering rate to file. 
    ! 
    WRITE(stdout,'(/5x,"Reading scattering rate from file"/)')
    !
    IF (mpime == ionode_id) THEN
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1, '(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1, '(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1, '(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate, FILE = name1, STATUS = 'old', IOSTAT = ios)
      WRITE(stdout,'(a16,a22)') '     Open file: ',name1   
      ! There are two comment line at the beginning of the file
      READ(iufilscatt_rate, *) dummy1
      READ(iufilscatt_rate, *) dummy1
      !
      DO ik = 1, nktotf
        !
        DO ibnd = 1, nbndfst
          !
          READ(iufilscatt_rate, *) ik_tmp, ibnd_tmp, etf_all(ibndmin - 1 + ibnd, ik), inv_tau_all(1, ibnd, ik) 
          inv_tau_all(1, ibnd, ik) = inv_tau_all(1, ibnd, ik) / (ryd2mev * meV2invps) 
          !
          ! Check that the file corresponds to the run we are making
          IF (ABS(ibnd_tmp - ibndmin - ibnd + 1) > 0)  CALL errore('scattering_read', &
            'Band read from the scattering_rate file do not match current calculation ', 1)
          ! 
        ENDDO
        ! Check that the file corresponds to the run we are making
        IF (ABS(ik_tmp - ik) > 0)  CALL errore('scattering_read', &
          'k-point read from the scattering_rate file do not match current calculation ', 1)
        !
      ENDDO
      !
      etf_all = etf_all / ryd2ev
      etf_all = etf_all + ef0
      !
      CLOSE(iufilscatt_rate)
    ENDIF
    CALL mp_bcast(etf_all, ionode_id, world_comm)
    CALL mp_bcast(inv_tau_all, ionode_id, world_comm)
    ! 
    WRITE(stdout,'(/5x,"Scattering rate read from file"/)')
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE scattering_read
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE electron_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !! 
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : ibndmax, ibndmin, lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufilsigma_all
    USE io_files,  ONLY : diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: sigmar_all(nbndfst, nktotf)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: sigmai_all(nbndfst, nktotf)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: zi_all(nbndfst, nktotf)
    !! Z parameter of electron-phonon self-energy accross all pools
    ! 
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf + 2)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      !
      lsigma_all = 3 * nbndfst * nktotf + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          i = i + 1
          aux(i) = sigmar_all(ibnd, ik)
        ENDDO
      ENDDO
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          i = i + 1
          aux(i) = sigmai_all(ibnd, ik)
        ENDDO
      ENDDO
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          i = i + 1
          aux(i) = zi_all(ibnd, ik)
        ENDDO
      ENDDO
      CALL diropn(iufilsigma_all, 'sigma_restart', lsigma_all, exst)
      CALL davcio(aux, lsigma_all, iufilsigma_all, 1, +1)
      CLOSE(iufilsigma_all)
    ENDIF
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN 
      sigmar_all(:, 1:lower_bnd - 1) = zero
      sigmai_all(:, 1:lower_bnd - 1) = zero
      zi_all(:, 1:lower_bnd - 1) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      sigmar_all(:, upper_bnd + 1:nktotf) = zero
      sigmai_all(:, upper_bnd + 1:nktotf) = zero
      zi_all(:, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE electron_write
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE electron_read(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !! 
    !! Self-energy reading
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : ibndmax, ibndmin, lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufilsigma_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE constants_epw, ONLY :  zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: sigmar_all(nbndfst, nktotf)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: sigmai_all(nbndfst, nktotf)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: zi_all(nbndfst, nktotf)
    !! Z parameter of electron-phonon self-energy accross all pools
    ! 
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf + 2)
    !! Vector to store the array
    ! 
    CHARACTER(LEN = 256) :: name1
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.sigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.sigma_restart'
#endif    
      INQUIRE(FILE = name1, EXIST = exst)
      ! 
      IF (exst) THEN ! read the file
        !
        lsigma_all = 3 * nbndfst * nktotf + 2
        CALL diropn(iufilsigma_all, 'sigma_restart', lsigma_all, exst)
        CALL davcio(aux, lsigma_all, iufilsigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('electron_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        ! 
        i = 2
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            sigmar_all(ibnd, ik) = aux(i)
          ENDDO
        ENDDO
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            sigmai_all(ibnd, ik) = aux(i)
          ENDDO
        ENDDO
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            zi_all(ibnd, ik) = aux(i)
          ENDDO
        ENDDO
        CLOSE(iufilsigma_all)
      ENDIF
    ENDIF
    ! 
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(sigmar_all, ionode_id, world_comm)
      CALL mp_bcast(sigmai_all, ionode_id, world_comm)
      CALL mp_bcast(zi_all, ionode_id, world_comm)
      ! 
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        sigmar_all(:, 1:lower_bnd - 1) = zero
        sigmai_all(:, 1:lower_bnd - 1) = zero
        zi_all(:, 1:lower_bnd - 1) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        sigmar_all(:, upper_bnd + 1:nktotf) = zero
        sigmai_all(:, upper_bnd + 1:nktotf) = zero
        zi_all(:, upper_bnd + 1:nktotf) = zero
      ENDIF
      ! 
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE electron_read
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE tau_write(iqq, totq, nktotf, second)
    !----------------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nstemp
    USE io_global, ONLY : meta_ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : diropn
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE constants_epw, ONLY : zero
    !! 
    !! Write scattering rates
    !! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! q-point from the selected ones within the fstick window. 
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    LOGICAL, INTENT(in) :: second
    !! IF we have two Fermi level
    ! 
    ! Local variable
    LOGICAL :: exst
    !! Does the file exists
    INTEGER :: i
    !! Running index for the vector
    INTEGER :: itemp
    !! Running index for the temperature
    INTEGER :: ltau_all
    !! Length of the vector
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    REAL(KIND = DP) :: aux(2 * nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array inv_tau_all and zi_all
    !
    IF (mpime == meta_ionode_id) THEN
      !
      ltau_all = 2 * nstemp * (nbndfst) * nktotf + 2
      ! First element is the iteration number
      aux(1) = REAL(iqq - 1, KIND = DP)   ! -1 because we will start at the next one. 
      aux(2) = REAL(totq, KIND = DP)
      i = 2
      ! 
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = inv_tau_all(itemp, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i +1
            aux(i) = zi_allvb(itemp, ibnd, ik) 
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iufiltau_all, 'tau_restart', ltau_all, exst)
      CALL davcio(aux, ltau_all, iufiltau_all, 1, +1 )
      CLOSE(iufiltau_all)
      ! 
      IF (second) THEN
        ! First element is the iteration number
        aux(1) = iqq - 1   ! -1 because we will start at the next one. 
        aux(2) = totq
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              aux(i) = inv_tau_allcb(itemp, ibnd, ik)
            ENDDO
          ENDDO
        ENDDO
        !
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              aux(i) = zi_allcb(itemp, ibnd, ik)     
            ENDDO
          ENDDO
        ENDDO
        ! 
        CALL diropn(iufiltau_all, 'tau_restart_CB', ltau_all, exst)
        CALL davcio(aux, ltau_all, iufiltau_all, 1, +1)
        CLOSE(iufiltau_all)   
      ENDIF
      ! 
    ENDIF
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) inv_tau_all(:, :, 1:lower_bnd - 1) = zero
    IF (upper_bnd < nktotf) inv_tau_all(:, :, upper_bnd + 1:nktotf) = zero
    IF (second) THEN
      IF (lower_bnd > 1) inv_tau_allcb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) inv_tau_allcb(:, :, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! Same for the Znk factor
    IF (lower_bnd > 1) zi_allvb(:, :, 1:lower_bnd - 1) = zero
    IF (upper_bnd < nktotf) zi_allvb(:, :, upper_bnd + 1:nktotf) = zero
    IF (second) THEN
      IF (lower_bnd > 1) zi_allcb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) zi_allcb(:, :, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE tau_write
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE tau_read(iqq, totq, nktotf, second)
    !----------------------------------------------------------------------------
    !!
    !! Scattering read
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, meta_ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE epwcom,    ONLY : nstemp
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : world_comm
    USE mp_world,  ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point from selecq.fmt
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    LOGICAL, INTENT(in) :: second
    !! IF we have two Fermi level
    ! 
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    CHARACTER(LEN = 256) :: name1
    !! Name of the file
    INTEGER :: i
    !! Iterative index
    INTEGER :: itemp
    !! Iterative temperature
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ltau_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    REAL(KIND = DP) :: aux(2 * nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array
    ! 
    IF (mpime == meta_ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart'
#endif 
      INQUIRE(FILE = name1, EXIST = exst)
      ! 
      IF (exst) THEN ! read the file
        !
        ltau_all = 2 * nstemp * nbndfst * nktotf + 2
        CALL diropn(iufiltau_all, 'tau_restart', ltau_all, exst)
        CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('tau_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        ! 
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              inv_tau_all(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        ! 
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              zi_allvb(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO 
        CLOSE(iufiltau_all)
      ENDIF
      ! 
      IF (second) THEN
        ! First inquire if the file exists
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart_CB1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart_CB'
#endif 
        INQUIRE(FILE = name1, EXIST = exst)
        ! 
        IF (exst) THEN ! read the file
          !
          ltau_all = nstemp * nbndfst * nktotf + 2
          CALL diropn(iufiltau_all, 'tau_restart_CB', ltau_all, exst)
          CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
          !
          ! First element is the iteration number
          iqq = INT(aux(1))
          iqq = iqq + 1 ! we need to start at the next q
          nqtotf_read = INT(aux(2))
          IF (nqtotf_read /= totq) CALL errore('tau_read',&
            &'Error: The current total number of q-point is not the same as the read one. ', 1)
          ! 
          i = 2
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                i = i + 1
                inv_tau_allcb(itemp, ibnd, ik) = aux(i)
              ENDDO
            ENDDO
          ENDDO
          ! 
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                i = i + 1
                zi_allcb(itemp, ibnd, ik) = aux(i)
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufiltau_all)
          WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau_CB: ', iqq, '/', totq 
        ENDIF
      ENDIF ! second
    ENDIF
    ! 
    CALL mp_bcast(exst, meta_ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, meta_ionode_id, world_comm)
      CALL mp_bcast(inv_tau_all, meta_ionode_id, world_comm)
      CALL mp_bcast(zi_allvb, meta_ionode_id, world_comm)
      IF (second) CALL mp_bcast(inv_tau_allcb, meta_ionode_id, world_comm)
      IF (second) CALL mp_bcast(zi_allcb, meta_ionode_id, world_comm)
      ! 
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1)      inv_tau_all(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) inv_tau_all(:, :, upper_bnd + 1:nktotf) = zero
      IF (lower_bnd > 1)      zi_allvb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) zi_allvb(:, :, upper_bnd + 1:nktotf) = zero
      !  
      IF (second) THEN
        ! Make everythin 0 except the range of k-points we are working on
        IF (lower_bnd > 1)      inv_tau_allcb(:, :, 1:lower_bnd - 1) = zero
        IF (upper_bnd < nktotf) inv_tau_allcb(:, :, upper_bnd + 1:nktotf) = zero
        IF (lower_bnd > 1)      zi_allcb(:, :, 1:lower_bnd - 1) = zero
        IF (upper_bnd < nktotf) zi_allcb(:, :, upper_bnd + 1:nktotf) = zero
      ENDIF 
      ! 
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau: ', iqq, '/', totq
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE tau_read
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE merge_read(nktotf, nqtotf_new, inv_tau_all_new)
    !----------------------------------------------------------------------------
    !!
    !! File merging
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : ibndmax, ibndmin, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : tmp_dir, diropn
    USE epwcom,    ONLY : nstemp, restart_filq
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    INTEGER, INTENT(out) :: nqtotf_new
    !! Total number of q-points
    REAL(KIND = DP), INTENT(inout) :: inv_tau_all_new(nstemp, nbndfst, nktotf)
    !! Scattering rate read from file restart_filq
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name of the file 
    LOGICAL :: exst
    !! Does the variable exist
    INTEGER :: i, iq, ios
    !! Iterative index
    INTEGER :: itemp
    !! Iterative temperature
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ltau_all
    !! Length of the vector
    INTEGER(KIND = 8) :: unf_recl
    !! Record length unit
    REAL(KIND = DP) :: aux(nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array 
    REAL(KIND = DP) :: dummy
    !! Test what the record length is
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
      name1 = TRIM(tmp_dir) // TRIM(restart_filq)
      INQUIRE(FILE = name1, EXIST = exst)
      ! 
      IF (exst) THEN ! read the file
        !
        ltau_all = nstemp * nbndfst * nktotf + 2
        !CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
        ! 
        INQUIRE(IOLENGTH = unf_recl) dummy  
        unf_recl = unf_recl * INT(ltau_all, KIND = KIND(unf_recl))
        OPEN(UNIT = iufiltau_all, FILE = restart_filq, IOSTAT = ios, FORM ='unformatted', &
             STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
        !  
        CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
        !
        ! First element is the iteration number
        iq = INT(aux(1))
        iq = iq + 1 ! we need to start at the next q
        nqtotf_new = INT(aux(2))
        ! 
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i +1
              inv_tau_all_new(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufiltau_all)
      ENDIF
    ENDIF
    ! 
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(nqtotf_new, ionode_id, world_comm)
      CALL mp_bcast(inv_tau_all_new, ionode_id, world_comm)
      ! 
      WRITE(stdout, '(a,a)' ) '     Correctly read file ', restart_filq
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE merge_read
    !----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io_transport
  !------------------------------------------------------------------------------
