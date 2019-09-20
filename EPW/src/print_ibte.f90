  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
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
!     WRITE(stdout,'(a,i9,E22.8)') '     Total number of element written ',ind_tot
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
