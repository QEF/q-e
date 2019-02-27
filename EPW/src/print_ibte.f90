  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2016-2018 Samuel Ponce'
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE print_ibte( iqq, iq, totq, ef0, efcb, first_cycle, ind_tot, ind_totcb, &
                         lrepmatw2, lrepmatw4, lrepmatw5, lrepmatw6 ) 
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine computes the scattering rate (inv_tau)
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP, i4b
  USE cell_base,     ONLY : omega
  USE io_global,     ONLY : stdout
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, eps_acustic, degaussw, & 
                            nstemp, scattering_serta, scattering_0rta, shortrange, &
                            restart, restart_freq, restart_filq, vme, ncarrier
  USE pwcom,         ONLY : ef
  USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, dmef, vmef, wf, wqf, & 
                            epf17, nkqtotf, inv_tau_all, inv_tau_allcb, &
                            xqf, zi_allvb, zi_allcb, xkf, wkf, dmef, vmef, nqf
  USE transportcom,  ONLY : transp_temp, lower_bnd
  USE constants_epw, ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev, & 
                            eps6, eps10, bohr2ang, ang2cm, eps4, eps8
  USE io_files,      ONLY : prefix, diropn, tmp_dir
  USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_global,     ONLY : world_comm, my_pool_id, npool
  USE io_global,     ONLY : ionode_id
  USE io_epw,        ONLY : iunepmat, iunepmatcb, iufilibtev_sup, iunrestart, &
                            iunsparseq, iunsparsek, iunsparsei, iunsparsej, iunsparset, &
                            iunsparseqcb, iunsparsekcb, iunsparseicb, iunsparsejcb, iunsparsetcb 
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, MPI_INTEGER8, &
                               MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_INTEGER4, &
                               MPI_MODE_CREATE, MPI_INFO_NULL, MPI_MODE_WRONLY, MPI_OFFSET
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (INOUT) :: first_cycle
  !! Use to determine weather this is the first cycle after restart 
  INTEGER, INTENT(IN) :: iqq
  !! Q-point index from selecq
  INTEGER, INTENT(IN) :: iq
  !! Q-point index
  INTEGER, INTENT(IN) :: totq
  !! Total number of q-points in selecq
  REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp. Could be unused (0).
#if defined(__MPI)  
  INTEGER(kind=MPI_OFFSET_KIND), INTENT(INOUT) :: ind_tot
  !! Total number of element written to file 
  INTEGER(kind=MPI_OFFSET_KIND), INTENT(INOUT) :: ind_totcb
  !! Total number of element written to file 
  INTEGER (kind=MPI_OFFSET_KIND), INTENT(inout) :: lrepmatw2
  !! Offset for that core
  INTEGER (kind=MPI_OFFSET_KIND), INTENT(inout) :: lrepmatw4
  !! Offset for that core
  INTEGER (kind=MPI_OFFSET_KIND), INTENT(inout) :: lrepmatw5
  !! Offset for that core
  INTEGER (kind=MPI_OFFSET_KIND), INTENT(inout) :: lrepmatw6
  !! Offset for that core
#else
  INTEGER, INTENT(INOUT) :: ind_tot
  !! Total number of element written to file 
  INTEGER, INTENT(INOUT) :: ind_totcb
  !! Total number of element written to file 
  INTEGER, INTENT (inout) :: lrepmatw2
  !! Offset for that core
  INTEGER, INTENT (inout) :: lrepmatw4
  !! Offset for that core
  INTEGER, INTENT (inout) :: lrepmatw5
  !! Offset for that core
  INTEGER, INTENT (inout) :: lrepmatw6
#endif   
  !
  ! Local variables
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
  INTEGER :: k_all(npool)
#if defined(__MPI)
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw
  !! Offset for that core
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw3
  !! Offset for that core
  INTEGER (kind=MPI_OFFSET_KIND) :: lsize
  !! Offset to tell where to start reading the file
#else
  INTEGER :: lrepmatw
  INTEGER :: lrepmatw3
  INTEGER :: lsize
  !! Offset to tell where to start reading the file
#endif
  CHARACTER (len=256) :: filint
  INTEGER :: ierr
  INTEGER :: ind(npool)
  INTEGER :: indcb(npool)
  INTEGER(kind=i4b) :: sparse_q((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparse_k((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparse_i((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparse_j((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparse_t((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparsecb_q((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparsecb_k((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparsecb_i((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparsecb_j((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  INTEGER(kind=i4b) :: sparsecb_t((ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf )
  !
  REAL(kind=DP) :: tmp
  !! Temporary variable
  REAL(kind=DP) :: dfnk
  !! Derivative of f_nk with respect to \varepsilon_nk
  REAL(kind=DP) :: ekk2
  !! Temporary variable to the eigenenergies for the degenerate average  
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
  REAL(kind=DP) :: weight
  !! Self-energy factor 
  REAL(KIND=DP) :: fmkq
  !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
  REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
  !! Electronic velocity $$v_{n\mathbf{k}}$$
  !REAL(kind=DP) :: trans_prob(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nstemp, nkf)
  REAL(kind=DP) :: trans_prob( (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) * nstemp * nkf)
  !! Temporary array to store the scattering rates
  !REAL(kind=DP) :: trans_probcb(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nstemp, nkf)
  REAL(kind=DP) :: trans_probcb( (ibndmax-ibndmin+1)*(ibndmax-ibndmin+1)* nstemp* nkf )
  !! Temporary array to store the scattering rates
  REAL(kind=DP) :: zi_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the zi
  REAL(KIND=DP), ALLOCATABLE :: inv_tau_all_new (:,:,:)
  !! New scattering rates to be merged
  REAL(KIND=DP) :: xkf_all(3,nkqtotf/2)
  !! k-points coordinate from all cores 
  REAL(KIND=DP) :: wkf_all(nkqtotf/2)
  !! Weights from all the cores
  REAL(KIND=DP) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
  !! Velocities from all the cores
  !
  REAL(KIND=DP) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  REAL(KIND=DP), EXTERNAL :: DDOT
  !! Dot product function
  REAL(KIND=DP), EXTERNAL :: efermig
  !! Function that returns the Fermi energy
  REAL(KIND=DP), EXTERNAL :: wgauss
  !! Compute the approximate theta function. Here computes Fermi-Dirac 
  REAL(KIND=DP), EXTERNAL :: w0gauss
  !! The derivative of wgauss:  an approximation to the delta function  
  REAL(kind=DP) :: carrier_density, fnk, inv_cell
  !  
  inv_cell = 1.0d0/omega
  ! 
  IF ( iqq == 1 ) THEN
    !
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Scattering rate for IBTE")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    WRITE(stdout,'(5x,"restart and restart_freq inputs deactivated (restart point at every q-points).")')
    WRITE(stdout,'(5x,"No intermediate mobility will be shown.")')
    !
    IF ( fsthick < 1.d3 ) THEN
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(5x,a,f10.6,a)' ) 'This is computed with respect to the fine Fermi level ',ef * ryd2ev, ' eV'
      WRITE(stdout, '(5x,a,f10.6,a,f10.6,a)' ) 'Only states between ',(ef-fsthick) * ryd2ev, ' eV and ',&
              (ef+fsthick) * ryd2ev, ' eV will be included'
      WRITE(stdout,'(5x,a/)')
    ENDIF
    lrepmatw = 0
    lrepmatw2 = 0
    lrepmatw4 = 0
    lrepmatw5 = 0
    lrepmatw6 = 0
    !
  ENDIF
  ! 
  ! In the case of a restart do not add the first step
  IF (first_cycle) THEN
    first_cycle = .FALSE.
  ELSE
    ! 
    trans_prob(:) = zero
    sparse_q(:)   = zero
    sparse_k(:)   = zero
    sparse_i(:)   = zero
    sparse_j(:)   = zero
    sparse_t(:)   = zero
    trans_probcb(:) = zero
    sparsecb_q(:) = zero
    sparsecb_k(:) = zero
    sparsecb_i(:) = zero
    sparsecb_j(:) = zero
    sparsecb_t(:) = zero

    !trans_probcb(:) = zero
    etf_all(:,:) = zero
    vkk_all(:,:,:) = zero
    ! local index for each q-point
    ind(:) = 0
    indcb(:) = 0
    ! 
    ! loop over temperatures
    DO itemp = 1, nstemp
      xkf_all(:,:) = 0.0d0
      wkf_all(:) = 0.0d0
      !
      etemp = transp_temp(itemp)
      !
      ! SP: Define the inverse so that we can efficiently multiply instead of
      ! dividing
      !
      inv_etemp = 1.0/etemp
      inv_degaussw = 1.0/degaussw
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
          
          xkf_all(:, ik+lower_bnd - 1 ) = xkf(:,ikk)
          wkf_all(ik+lower_bnd - 1 ) = wkf(ikk)
          ! 
          DO ibnd = 1, ibndmax-ibndmin+1
            !  energy at k (relative to Ef)
            ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
            !
            ! This is to know if we need to store the data 
            ! derivative Fermi distribution
            ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
            dfnk = w0gauss( ekk*inv_etemp, -99 )*inv_etemp
            !  
            DO jbnd = 1, ibndmax-ibndmin+1
              !
              !  energy and fermi occupation at k+q
              ekq = etf (ibndmin-1+jbnd, ikq) - ef0(itemp)
              fmkq = wgauss( -ekq*inv_etemp, -99)
              !
              ! We perform a sum over the modes  
              tmp = zero
              DO imode = 1, nmodes
                !
                ! the phonon frequency and bose occupation
                wq = wf (imode, iq)
                !
                ! SP : Avoid if statement in inner loops
                ! the coupling from Gamma acoustic phonons is negligible
                IF ( wq .gt. eps_acustic ) THEN
                  g2_tmp = 1.0
                  wgq = wgauss( -wq*inv_etemp, -99)
                  wgq = wgq / ( one - two * wgq )
                  ! SP : Define the inverse for efficiency
                  inv_wq =  1.0/(two * wq) 
                ELSE
                  g2_tmp = 0.0
                  wgq = 0.0
                  inv_wq = 0.0
                ENDIF
                !
                ! here we take into account the zero-point sqrt(hbar/2M\omega)
                ! with hbar = 1 and M already contained in the eigenmodes
                ! g2 is Ry^2, wkf must already account for the spin factor
                !
                ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                ! need to treat it like the normal g with abs(g).
                IF ( shortrange .AND. ( abs(xqf (1, iq))> eps8 .OR. abs(xqf (2, iq))> eps8 &
                   .OR. abs(xqf (3, iq))> eps8 )) THEN
                  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                  !     number, in which case its square will be a negative number. 
                  g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp, KIND=DP )
                ELSE
                  g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                ENDIF
                !
                ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                !
                ! transition probability 
                ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                !
                ! This is summed over modes
                tmp = tmp  + pi * wqf(iq) * g2 * ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                !
              ENDDO !imode
              ! Only save the onces that really contribute
              IF (ABS(tmp * dfnk) > 1d-40 ) THEN
                !IF (ind(my_pool_id+1)<10) print*,'ik ibnd jbnd ',ik, ibnd, jbnd
                !IF (ind(my_pool_id+1)<10) print*,'VB ekk ',ekk
                !IF (ind(my_pool_id+1)<10) print*,'VB ekq ',ekq
                !IF (ind(my_pool_id+1)<10) print*,'VB f ',fmkq
                !IF (ind(my_pool_id+1)<10) print*,'VB g1 ',w0g1
                !IF (ind(my_pool_id+1)<10) print*,'VB g2 ',w0g2
                !IF (ind(my_pool_id+1)<10) print*,'VB tmp ',tmp
                !IF (ind (my_pool_id+1)==0) print*,tmp
                !IF (ind (my_pool_id+1)==0) print*,wqf(iq)
                !IF (ind (my_pool_id+1)==0) print*, g2
                !IF (ind (my_pool_id+1)==0) print*, fmkq
                !IF (ind (my_pool_id+1)==0) print*, wgq
                !IF (ind (my_pool_id+1)==0) print*, w0g1
                !IF (ind (my_pool_id+1)==0) print*, w0g2
                !
                ind (my_pool_id+1) = ind (my_pool_id+1) + 1 
                ! 
                trans_prob( ind(my_pool_id+1) ) = tmp
                sparse_q  ( ind(my_pool_id+1) ) = iq
                sparse_k  ( ind(my_pool_id+1) ) = ik+lower_bnd-1
                !print*,'index sp_k ',ind(my_pool_id+1), sparse_k  (ind(my_pool_id+1) )
                sparse_i  ( ind(my_pool_id+1) ) = ibnd
                sparse_j  ( ind(my_pool_id+1) ) = jbnd 
                sparse_t  ( ind(my_pool_id+1) ) = itemp
               ! IF (ik+lower_bnd-1 == 2) print*,ind (my_pool_id+1)+ind_tot, tmp
               ! print*,ind (my_pool_id+1)+ind_tot, iq, ik, ibnd, jbnd, itemp
                !  
              ENDIF 
              !
            ENDDO !jbnd
          ENDDO ! ibnd
          !
          ! In this case we are also computing the scattering rate for another Fermi level position
          ! This is used to compute both the electron and hole mobility at the same time.  
          IF ( ABS(efcb(itemp)) > eps4 ) THEN
            ! 
            DO ibnd = 1, ibndmax-ibndmin+1
              !
              !  energy at k (relative to Ef)
              ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
              ! This is to know if we need to store the data 
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk*inv_etemp, -99 )*inv_etemp
              !
              DO jbnd = 1, ibndmax-ibndmin+1
                !
                !  energy and fermi occupation at k+q
                ekq = etf (ibndmin-1+jbnd, ikq) - efcb(itemp)
                fmkq = wgauss( -ekq*inv_etemp, -99)
                ! 
                tmp = zero
                DO imode = 1, nmodes
                  !
                  ! the phonon frequency and bose occupation
                  wq = wf (imode, iq)
                  !
                  ! SP : Avoid if statement in inner loops
                  ! the coupling from Gamma acoustic phonons is negligible
                  IF ( wq .gt. eps_acustic ) THEN
                    g2_tmp = 1.0
                    wgq = wgauss( -wq*inv_etemp, -99)
                    wgq = wgq / ( one - two * wgq )
                    ! SP : Define the inverse for efficiency
                    inv_wq =  1.0/(two * wq)
                  ELSE
                    g2_tmp = 0.0
                    wgq = 0.0
                    inv_wq = 0.0
                  ENDIF
                  !
                  ! here we take into account the zero-point sqrt(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                  ! need to treat it like the normal g with abs(g).
                  IF ( shortrange .AND. ( abs(xqf (1, iq))> eps8 .OR. abs(xqf (2, iq))> eps8 &
                     .OR. abs(xqf (3, iq))> eps8 )) THEN
                    ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                    !     number, in which case its square will be a negative number. 
                    g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp, KIND=DP)
                  ELSE
                    g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                  ENDIF
                  !
                  ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                  w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                  w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                  !
                  ! transition probability 
                  ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                  ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                  !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                  !
                  ! This is summed over modes
                  tmp = tmp  + pi * wqf(iq) * g2 * ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                  !  
                ENDDO ! imode
                ! 
                ! Only save the onces that really contribute
                IF (ABS(tmp * dfnk) > 1d-40 ) THEN
                  indcb (my_pool_id+1) = indcb (my_pool_id+1) + 1
                  ! 
                  trans_probcb  ( indcb(my_pool_id+1) ) = tmp
                  sparsecb_q    ( indcb(my_pool_id+1) ) = iq
                  sparsecb_k    ( indcb(my_pool_id+1) ) = ik+lower_bnd-1
                  sparsecb_i    ( indcb(my_pool_id+1) ) = ibnd
                  sparsecb_j    ( indcb(my_pool_id+1) ) = jbnd
                  sparsecb_t    ( indcb(my_pool_id+1) ) = itemp
                  !  
                ENDIF

                ! 
              ENDDO !jbnd
              !
            ENDDO !ibnd
            ! 
          ENDIF ! ABS(efcb) < eps
          !
          !print*,'trans_probcb(jbnd,ibnd,itemp,k_all(my_pool_id+1)) ',trans_probcb(1,1,1,k_all(my_pool_id+1))
        ENDIF ! endif  fsthick
        !
      ENDDO ! end loop on k
    ENDDO ! itemp
    ! If the q-point is taken, write on file
    CALL mp_sum( ind,    world_comm )
    CALL mp_sum( indcb,    world_comm )
    ! 
    ! SP - IBTE only with if EPW compiled with MPI
#if defined(__MPI)
    IF ( sum(ind) > 0 ) THEN
      ! 
      IF ( my_pool_id == 0 ) ind_tot = ind_tot + SUM(ind)
      !CALL mp_bcast (ind_tot, ionode_id, world_comm)  
      CALL MPI_BCAST( ind_tot, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
! Can be uncommented
!      WRITE(stdout,'(a,i9,E22.8)') '     Total number of element written ',ind_tot
      ! 
      ! Size of what we write
      lsize =  INT( ind(my_pool_id+1), kind = MPI_OFFSET_KIND )

      ! Offset where we need to start writing (we increment for each q-points)
      lrepmatw = lrepmatw2 + &
            INT( SUM( ind(1:my_pool_id+1)) - ind(my_pool_id+1), kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      lrepmatw3 = lrepmatw4 + &
            INT( SUM( ind(1:my_pool_id+1)) - ind(my_pool_id+1), kind = MPI_OFFSET_KIND ) * 4_MPI_OFFSET_KIND
      ! 
      CALL MPI_FILE_SEEK(iunepmat, lrepmatw,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunepmat, trans_prob, lsize, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      !  
      CALL MPI_FILE_SEEK (iunsparseq, lrepmatw3,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparseq, sparse_q, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      !  
      CALL MPI_FILE_SEEK (iunsparsek, lrepmatw3,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparsek, sparse_k, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      !  
      CALL MPI_FILE_SEEK (iunsparsei, lrepmatw3,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparsei, sparse_i, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      !  
      CALL MPI_FILE_SEEK (iunsparsej, lrepmatw3,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparsej, sparse_j, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      !  
      CALL MPI_FILE_SEEK (iunsparset, lrepmatw3,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparset, sparse_t, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      ! 
      ! Offset for the next q iteration
      lrepmatw2 = lrepmatw2 + INT( SUM( ind(:) ), kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      lrepmatw4 = lrepmatw4 + INT( SUM( ind(:) ), kind = MPI_OFFSET_KIND ) * 4_MPI_OFFSET_KIND
      ! 
      ! now write in the support file
      CALL mp_sum(xkf_all, world_comm)
      CALL mp_sum(wkf_all, world_comm)
      ! 
    ENDIF 
    IF ( sum(indcb) > 0 ) THEN
      ! 
      IF ( my_pool_id == 0 ) ind_totcb = ind_totcb + SUM(indcb)
      !CALL mp_bcast (ind_totcb, ionode_id, world_comm)
      CALL MPI_BCAST( ind_totcb, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
!      WRITE(stdout,'(a,i9,E22.8)') '     Total number of element written in electron ',ind_totcb
      ! 
      ! Size of what we write
      lsize =  INT( indcb(my_pool_id+1), kind = MPI_OFFSET_KIND )

      ! Offset where we need to start writing (we increment for each q-points)
      lrepmatw = lrepmatw5 + &
            INT( SUM( indcb(1:my_pool_id+1)) - indcb(my_pool_id+1), kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      lrepmatw3 = lrepmatw6 + &
            INT( SUM( indcb(1:my_pool_id+1)) - indcb(my_pool_id+1), kind = MPI_OFFSET_KIND ) * 4_MPI_OFFSET_KIND
      ! 
      CALL MPI_FILE_SEEK(iunepmatcb, lrepmatw,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunepmatcb, trans_probcb, lsize, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )
      !  
      CALL MPI_FILE_SEEK (iunsparseqcb, lrepmatw3,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparseqcb, sparsecb_q, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )
      !  
      CALL MPI_FILE_SEEK (iunsparsekcb, lrepmatw3,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparsekcb, sparsecb_k, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )
      !  
      CALL MPI_FILE_SEEK (iunsparseicb, lrepmatw3,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparseicb, sparsecb_i, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )
      !  
      CALL MPI_FILE_SEEK (iunsparsejcb, lrepmatw3,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparsejcb, sparsecb_j, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )
      !  
      CALL MPI_FILE_SEEK (iunsparsetcb, lrepmatw3,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_WRITE(iunsparsetcb, sparsecb_t, lsize, MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )
      ! 
      ! Offset for the next q iteration
      lrepmatw5 = lrepmatw5 + INT( SUM( indcb(:) ), kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      lrepmatw6 = lrepmatw6 + INT( SUM( indcb(:) ), kind = MPI_OFFSET_KIND ) * 4_MPI_OFFSET_KIND
      ! 
    ENDIF ! indcb
#endif
    ! 
    ! Save to file restart information in formatted way for possible restart
    IF (my_pool_id == 0) THEN
      OPEN(unit=iunrestart,file='restart_ibte.fmt')
      WRITE (iunrestart,*) iqq
      WRITE (iunrestart,*) ind_tot
      WRITE (iunrestart,*) ind_totcb
      WRITE (iunrestart,*) lrepmatw2
      WRITE (iunrestart,*) lrepmatw4
      WRITE (iunrestart,*) lrepmatw5
      WRITE (iunrestart,*) lrepmatw6
      CLOSE(iunrestart)
    ENDIF    
    ! 
  ENDIF ! first_cycle
  ! 
  IF ( iqq == totq ) THEN
    wkf_all(:) = zero
    ! Computes the k-velocity
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ! 
      wkf_all( ik+lower_bnd -1 ) = wkf(ikk) 
      ! 
      DO ibnd = 1, ibndmax-ibndmin+1
        IF ( vme ) THEN
          vkk_all(:, ibnd, ik + lower_bnd - 1) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
        ELSE
          vkk_all(:,ibnd, ik + lower_bnd -1 ) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
        ENDIF  
        etf_all(ibnd, ik+lower_bnd-1) = etf(ibndmin-1+ibnd, ikk)
      ENDDO
    ENDDO 
    CALL mp_sum ( vkk_all, world_comm ) 
    CALL mp_sum ( etf_all, world_comm ) 
    CALL mp_sum ( wkf_all, world_comm )
    ! 
    IF ( my_pool_id == 0 ) THEN

      ! Now write total number of q-point inside and k-velocity
      !
      OPEN(iufilibtev_sup,file='IBTEvel_sup.fmt', form='formatted')
      WRITE(iufilibtev_sup,'(a)') '# Number of elements in hole and electrons  '
      WRITE(iufilibtev_sup,'(2i16)') ind_tot, ind_totcb
      WRITE(iufilibtev_sup,'(a)') '# itemp    ef0    efcb'
      DO itemp=1, nstemp
        WRITE(iufilibtev_sup,'(i8,2E22.12)') itemp, ef0(itemp), efcb(itemp)
      ENDDO
      WRITE(iufilibtev_sup,'(a)') '# ik  ibnd      velocity (x,y,z)              eig     weight '
      DO ik = 1, nkqtotf/2
        DO ibnd = 1, ibndmax-ibndmin+1
          WRITE(iufilibtev_sup,'(i8,i6,5E22.12)') ik, ibnd, vkk_all(:,ibnd,ik), etf_all(ibnd, ik), wkf_all(ik)
        ENDDO
      ENDDO
      CLOSE(iufilibtev_sup)
      ! 
    ENDIF ! master
    ! 
    ! Now print the carrier density for checking
    DO itemp=1, nstemp
      etemp = transp_temp(itemp)
      carrier_density = 0.0
      ! 
      IF ( ncarrier < 0.0 ) THEN ! VB
        DO ik = 1, nkf
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibnd, ik+lower_bnd-1 ) < ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik+lower_bnd-1 ) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik+lower_bnd-1 ) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
        CALL mp_sum( carrier_density, world_comm )
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6)') etemp *ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density
      ELSE ! CB
        DO ik = 1, nkf
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibnd, ik+lower_bnd-1 ) > efcb(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik+lower_bnd-1 ) - efcb(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik+lower_bnd-1 ) *  fnk
            ENDIF
          ENDDO
        ENDDO
        CALL mp_sum( carrier_density, world_comm )
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6)') etemp *ryd2ev/kelvin2eV, &
                      efcb(itemp)*ryd2ev,  carrier_density
      ENDIF ! ncarrier
    ENDDO
    ! 
  ENDIF ! iqq
  !
  RETURN
  !
  END SUBROUTINE print_ibte
  !-----------------------------------------------------------------------
