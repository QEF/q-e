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
    SUBROUTINE print_gkk ( iq )
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
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                              nkf, epf17, xkf, nkqtotf, wf
    USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, eps8
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE division,      ONLY : fkbounds
    !
    implicit none
    !
    INTEGER, INTENT (in) :: iq
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
    REAL(kind=DP) :: xkf_all(3, nkqtotf)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(kind=DP) :: etf_all(nbndsub,nkqtotf)
    !! Collect eigenenergies from all pools in parallel case
    REAL(kind=DP) :: wq
    !! Phonon frequency 
    REAL(kind=DP) :: w_1
    !! Temporary phonon freq. 1
    REAL(kind=DP) :: w_2
    !! Temporary phonon freq. 2
    REAL(kind=DP) :: gamma
    !! Temporary electron-phonon matrix element
    REAL(kind=DP) :: ekk
    !! Eigenenergies at k
    REAL(kind=DP) :: ekq
    !! Eigenenergies at k+q
    REAL(kind=DP) :: g2
    !! Temporary electron-phonon matrix element square
    REAL(kind=DP), ALLOCATABLE :: epc(:,:,:,:)
    !! g vectex accross all pools 
    REAL(kind=DP), ALLOCATABLE :: epc_sym(:,:,:)
    !! Temporary g-vertex for each pool 
    ! 
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds( nkqtotf/2, lower_bnd, upper_bnd )
    ! 
    ALLOCATE ( epc (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkqtotf/2) )
    ALLOCATE ( epc_sym (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes) )
    ! 
    epc(:,:,:,:)   = zero
    epc_sym(:,:,:) = zero
    !
    ! First do the average over bands and modes for each pool
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      ikq = ikk + 1
      ! 
      DO nu = 1, nmodes
        wq = wf (nu, iq)
        !DO ibnd = ibndmin, ibndmax
        DO ibnd = 1, ibndmax-ibndmin+1
          DO jbnd = 1, ibndmax-ibndmin+1
            gamma = ( ABS(epf17 (jbnd, ibnd, nu, ik)) )**two
            IF (wq > 0.d0) then
                gamma = gamma / (two * wq)
            ELSE
                gamma = 0.d0
            ENDIF
            gamma = sqrt(gamma)
            ! gamma = |g| [Ry]
            epc(ibnd,jbnd,nu,ik+lower_bnd-1) = gamma
          ENDDO ! jbnd
        ENDDO   ! ibnd        
      ENDDO ! loop on modes
      !
      !  Here we "SYMMETRIZE": actually we simply take the averages over
      !  degenerate states, it is only a convention because g is gauge-dependent!
      !
      ! first the phonons
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
          DO nu = 1, nmodes
            w_1 = wf(nu,iq)
            g2 = 0.d0
            n  = 0
            DO mu = 1, nmodes
              w_2 = wf(mu,iq)
              IF ( abs(w_2-w_1) < eps8 ) THEN
                n = n + 1
                g2 = g2 + epc(ibnd,jbnd,mu,ik+lower_bnd-1)*epc(ibnd,jbnd,mu,ik+lower_bnd-1)
              ENDIF
            ENDDO
            g2 = g2 / float(n)
            epc_sym (ibnd, jbnd, nu) = sqrt (g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:,:,:,ik+lower_bnd-1) = epc_sym
      ! Then the k electrons
      DO nu = 1, nmodes
        DO jbnd = 1, ibndmax-ibndmin+1
          DO ibnd = 1, ibndmax-ibndmin+1
            w_1 = etf (ibndmin-1+ibnd, ikk)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, ibndmax-ibndmin+1
              w_2 = etf (ibndmin-1+pbnd, ikk)
              IF ( abs(w_2-w_1) < eps8 ) THEN
                n = n + 1
                g2 = g2 + epc(pbnd,jbnd,nu,ik+lower_bnd-1)*epc(pbnd,jbnd,nu,ik+lower_bnd-1)
              ENDIF
            ENDDO
            g2 = g2 / float(n)
            epc_sym (ibnd, jbnd, nu) = sqrt (g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:,:,:,ik+lower_bnd-1) = epc_sym
      !
      ! and finally the k+q electrons
      DO nu = 1, nmodes
        DO ibnd = 1, ibndmax-ibndmin+1
          DO jbnd = 1, ibndmax-ibndmin+1
            w_1 = etf (ibndmin-1+jbnd, ikq)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, ibndmax-ibndmin+1
              w_2 = etf(ibndmin-1+pbnd, ikq)
              IF ( abs(w_2-w_1) < eps8 ) then
                n = n + 1
                g2 = g2 + epc(ibnd,pbnd,nu,ik+lower_bnd-1)*epc(ibnd,pbnd,nu,ik+lower_bnd-1)
              ENDIF
            ENDDO
            g2 = g2 / float(n)
            epc_sym (ibnd, jbnd, nu) = sqrt (g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:,:,:,ik+lower_bnd-1) = epc_sym
      ! 
    ENDDO ! k-points
    ! 
    ! We need quantity from all the pools
    xkf_all(:,:) = zero
    etf_all(:,:) = zero
    !
#if defined(__MPI)
    !
    ! Note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
    CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
    CALL mp_sum( epc, inter_pool_comm )
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
    IF (mpime.eq.ionode_id) THEN
      !
      WRITE(stdout, '(5x,a)') ' Electron-phonon vertex |g| (meV)'
      !
      WRITE(stdout, '(/5x,"iq = ",i7," coord.: ", 3f12.7)') iq, xqf (:, iq)
      DO ik = 1, nkqtotf/2
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout,'(5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:,ikk)
        WRITE(stdout, '(5x,a)') ' ibnd     jbnd     imode   enk[eV]    enk+q[eV]  omega(q)[meV]   |g|[meV]'
        WRITE(stdout,'(5x,a)') repeat('-',78)
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          ekk = etf_all (ibndmin-1+ibnd, ikk) 
          DO jbnd = 1, ibndmax-ibndmin+1
            ekq = etf_all (ibndmin-1+jbnd, ikq) 
            DO nu = 1, nmodes
              WRITE(stdout,'(3i9,3f12.4,1e20.10)') ibndmin-1+ibnd, ibndmin-1+jbnd, nu, ryd2ev * ekk, ryd2ev * ekq, &
                                           ryd2mev * wf(nu,iq), ryd2mev * epc(ibnd,jbnd,nu,ik)
            ENDDO
          ENDDO  
          !
          !
        ENDDO
        WRITE(stdout,'(5x,a/)') repeat('-',78)
        !
      ENDDO
    ENDIF ! master node
    ! 
    DEALLOCATE ( epc )
    DEALLOCATE ( epc_sym )
    !
    END SUBROUTINE print_gkk
    ! 
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_serta_sym( F_SERTA, BZtoIBZ, s_BZtoIBZ, BZtoIBZ_mat, vkk_all, etf_all, wkf_all, ef0) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine prints the SERTA mobility using k-point symmetry.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, omega, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              electron_SI, bohr2ang, ang2cm, hbarJ
    USE symm_base,     ONLY : nrot
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm, my_pool_id
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: BZtoIBZ(nkf1*nkf2*nkf3)
    !! BZ to IBZ mapping
    INTEGER, INTENT(in) :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
    !! Corresponding symmetry matrix
    INTEGER, INTENT(in) :: BZtoIBZ_mat(nrot,nkqtotf/2)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(kind=DP), INTENT(in) :: F_SERTA(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)  
    !! SERTA solution 
    REAL(KIND=DP), INTENT(IN) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
    !! Velocity
    REAL(KIND=DP), INTENT(IN) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
    !! Eigenenergies of k
    REAL(KIND=DP), INTENT(IN) :: wkf_all(nkqtotf/2)
    !! Weight of k
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! The Fermi level 

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
    ! 
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: tdf_sigma(9)
    !! Transport distribution function
    REAL(KIND=DP) :: Sigma(9,nstemp)
    !! Electrical conductivity
    REAL(kind=DP) :: Fi_cart(3)
    !! Cartesian Fi_all 
    REAL(kind=DP) :: Fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    REAL(kind=DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(kind=DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(kind=DP) :: sa(3,3)
    !! Symmetry matrix in crystal
    REAL(kind=DP) :: sb(3,3)
    !! Symmetry matrix (intermediate step)
    REAL(kind=DP) :: sr(3,3)
    !! Symmetry matrix in cartesian coordinate 
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
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
    REAL(KIND=DP) :: sigma_up(3,3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND=DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND=DP) :: sigma_vect(3,3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND=DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    !
    inv_cell = 1.0d0/omega
    ! 
    ! Hole
    IF (ncarrier < -1E5) THEN
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) < ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              vk_cart(:) = vkk_all(:,ibnd, ik)
              Fi_cart(:) = F_SERTA(:, ibnd, ik, itemp)
              !
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb=1, nrot
                IF (BZtoIBZ_mat(nb,ik) > 0 ) THEN 
                  ikbz = BZtoIBZ_mat(nb,ik)
                  ! 
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                  sb = matmul ( bg, sa )
                  sr (:,:) = matmul ( at, transpose (sb) )
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, Fi_cart(:),1 ,0.d0 , Fi_rot(:), 1 )
                  !
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      IF (noncolin) THEN
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                        !tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                      ELSE
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                        !tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF ! BZ
              ENDDO ! ikb
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik ) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      !CALL mp_sum( Sigma(:,:), world_comm )
      ! 
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibnd, ik ) < ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik ) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik ) * (1.0d0 - fnk )
            ENDIF
          ENDDO
          !IF (my_pool_id == 0 ) write(990,*)ik, etf_all(1,ik + lower_bnd - 1), carrier_density
          !IF (my_pool_id == 1 ) write(991,*)ik, etf_all(1,ik + lower_bnd - 1), carrier_density
        ENDDO
        !CALL mp_sum( carrier_density, world_comm )
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
      ENDDO
      ! 
    ENDIF
    ! 
    ! Now electron mobility
    IF (ncarrier > 1E5) THEN
      ! Needed because of residual values from the hole above
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              vk_cart(:) = vkk_all(:,ibnd, ik)
              Fi_cart(:) = F_SERTA(:, ibnd, ik, itemp) 
              ! 
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb=1, nrot
                IF (BZtoIBZ_mat(nb,ik) > 0 ) THEN
                  ikbz = BZtoIBZ_mat(nb,ik)
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                  sb = matmul ( bg, sa )
                  sr (:,:) = matmul ( at, transpose (sb) )
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, Fi_cart(:),1 ,0.d0 , Fi_rot(:), 1 )
                  !
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      IF (noncolin) THEN
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                        !tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                      ELSE
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                        !tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF ! BZ
              ENDDO ! ikb
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      !CALL mp_sum( Sigma(:,:), world_comm )
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for electron conduction
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        !CALL mp_sum( carrier_density, world_comm )
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
      ENDDO
      ! 
    ENDIF
    
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_serta_sym
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE print_serta( F_SERTA, vkk_all, etf_all, wkf_all, ef0) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine prints the SERTA mobility
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, omega, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                               electron_SI, bohr2ang, ang2cm, hbarJ
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm, my_pool_id
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(kind=DP), INTENT(in) :: F_SERTA(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)  
    !! SERTA solution 
    REAL(KIND=DP), INTENT(IN) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
    !! Velocity
    REAL(KIND=DP), INTENT(IN) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
    !! Eigenenergies of k
    REAL(KIND=DP), INTENT(IN) :: wkf_all(nkqtotf/2)
    !! Weight of k
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! The Fermi level 

    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik 
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: ij  
    !! Combined x,y,z index
    INTEGER :: i,j
    !! Cartesian index
    ! 
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: tdf_sigma(9)
    !! Transport distribution function
    REAL(KIND=DP) :: Sigma(9,nstemp)
    !! Electrical conductivity
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
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
    REAL(KIND=DP) :: sigma_up(3,3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND=DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND=DP) :: sigma_vect(3,3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND=DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 

    !
    inv_cell = 1.0d0/omega
    ! 
    ! Hole
    IF (ncarrier < -1E5) THEN
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) < ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              ! 
              ij = 0
              DO j = 1, 3
                DO i = 1, 3
                  ij = ij + 1
                  tdf_sigma(ij) = vkk_all(i, ibnd, ik) * F_SERTA(j, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik ) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibnd, ik ) < ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik ) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik ) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
      ENDDO
      ! 
    ENDIF
    ! 
    ! Now electron mobility
    IF (ncarrier > 1E5) THEN
      ! Needed because of residual values from the hole above
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              ij = 0
              DO j = 1, 3
                DO i = 1, 3
                  ij = ij + 1
                  tdf_sigma(ij) = vkk_all(i, ibnd, ik) * F_SERTA(j, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      !CALL mp_sum( Sigma(:,:), world_comm )
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for electron conduction
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        !CALL mp_sum( carrier_density, world_comm )
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
      ENDDO
      ! 
    ENDIF
    
    RETURN
    !
    END SUBROUTINE print_serta
    !-----------------------------------------------------------------------
    !  
    !-----------------------------------------------------------------------
    SUBROUTINE print_mob_sym( F_out, BZtoIBZ, s_BZtoIBZ, BZtoIBZ_mat, vkk_all, etf_all, wkf_all, ef0, av_mob) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine prints the mobility using k-point symmetry ( electron or hole )
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, omega, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              electron_SI, bohr2ang, ang2cm, hbarJ
    USE symm_base,     ONLY : nrot
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: BZtoIBZ(nkf1*nkf2*nkf3)
    !! BZ to IBZ mapping
    INTEGER, INTENT(in) :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
    !! Corresponding symmetry matrix
    INTEGER, INTENT(in) :: BZtoIBZ_mat(nrot,nkqtotf/2)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(kind=DP), INTENT(in) :: F_out(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)  
    !! SERTA solution 
    REAL(KIND=DP), INTENT(IN) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
    !! Velocity
    REAL(KIND=DP), INTENT(IN) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
    !! Eigenenergies of k
    REAL(KIND=DP), INTENT(IN) :: wkf_all(nkqtotf/2)
    !! Weight of k
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! The Fermi level 
    REAL(KIND=DP), INTENT(INOUT) :: av_mob(nstemp)
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
    ! 
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: tdf_sigma(9)
    !! Transport distribution function
    REAL(KIND=DP) :: Sigma(9,nstemp)
    !! Electrical conductivity
    REAL(kind=DP) :: Fi_cart(3)
    !! Cartesian Fi_all 
    REAL(kind=DP) :: Fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    REAL(kind=DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(kind=DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(kind=DP) :: sa(3,3)
    !! Symmetry matrix in crystal
    REAL(kind=DP) :: sb(3,3)
    !! Symmetry matrix (intermediate step)
    REAL(kind=DP) :: sr(3,3)
    !! Symmetry matrix in cartesian coordinate 
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
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
    REAL(KIND=DP) :: sigma_up(3,3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND=DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND=DP) :: sigma_vect(3,3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND=DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 

    !
    inv_cell = 1.0d0/omega
    ! 
    IF (ncarrier < -1E5) THEN ! If true print hole
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) < ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              vk_cart(:) = vkk_all(:,ibnd, ik)
              Fi_cart(:) = F_out(:, ibnd, ik, itemp)
              ! 
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb=1, nrot
                IF (BZtoIBZ_mat(nb,ik) > 0 ) THEN
                  ikbz = BZtoIBZ_mat(nb,ik)
                  ! 
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                  sb = matmul ( bg, sa )
                  sr (:,:) = matmul ( at, transpose (sb) )
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, Fi_cart(:),1 ,0.d0 , Fi_rot(:), 1 )
                  !
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      IF (noncolin) THEN
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                      ELSE
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF ! BZ
              ENDDO ! ikb
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      !CALL mp_sum( Sigma(:,:), world_comm )
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibnd, ik) < ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
        !CALL mp_sum( carrier_density, world_comm )
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        av_mob(itemp) = mobility
      ENDDO
      ! 
    ! Now electron mobility
    ENDIF
    IF (ncarrier > 1E5 ) THEN
      ! Needed because of residual values from the hole above
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              vk_cart(:) = vkk_all(:, ibnd, ik)
              Fi_cart(:) = F_out(:, ibnd, ik, itemp)
              ! 
              ! Loop on the point equivalent by symmetry in the full BZ
              DO nb=1, nrot
                IF (BZtoIBZ_mat(nb,ik) > 0 ) THEN
                  ikbz = BZtoIBZ_mat(nb,ik)
                  ! 
                  ! Transform the symmetry matrix from Crystal to
                  ! cartesian
                  sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                  sb = matmul ( bg, sa )
                  sr (:,:) = matmul ( at, transpose (sb) )
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                  ! 
                  CALL dgemv( 'n', 3, 3, 1.d0,&
                    sr, 3, Fi_cart(:),1 ,0.d0 , Fi_rot(:), 1 )
                  !
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      ! The factor two in the weight at the end is to
                      ! account for spin
                      IF (noncolin) THEN
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                      ELSE
                        tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * Fi_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF ! BZ
              ENDDO ! ikb
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      !CALL mp_sum( Sigma(:,:), world_comm )
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for electron conduction
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
        !CALL mp_sum( carrier_density, world_comm )
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        av_mob(itemp) = mobility
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
    SUBROUTINE print_mob( F_out, vkk_all, etf_all, wkf_all, ef0, av_mob) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine prints the mobility ( electron or hole )
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : at, omega, bg
    USE epwcom,        ONLY : int_mob, ncarrier, nstemp, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, ibndmax, ibndmin, nkqtotf 
    USE transportcom,  ONLY : lower_bnd, transp_temp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                               electron_SI, bohr2ang, ang2cm, hbarJ
    USE noncollin_module, ONLY : noncolin
    USE mp_global,     ONLY : world_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(kind=DP), INTENT(IN) :: F_out(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)  
    !! SERTA solution 
    REAL(KIND=DP), INTENT(IN) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
    !! Velocity
    REAL(KIND=DP), INTENT(IN) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
    !! Eigenenergies of k
    REAL(KIND=DP), INTENT(IN) :: wkf_all(nkqtotf/2)
    !! Weight of k
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! The Fermi level 
    REAL(KIND=DP), INTENT(INOUT) :: av_mob(nstemp)
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
    ! 
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: tdf_sigma(9)
    !! Transport distribution function
    REAL(KIND=DP) :: Sigma(9,nstemp)
    !! Electrical conductivity
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$ 
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
    REAL(KIND=DP) :: sigma_up(3,3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND=DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND=DP) :: sigma_vect(3,3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND=DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 

    !
    inv_cell = 1.0d0/omega
    ! 
    IF (ncarrier < -1E5) THEN ! If true print hole
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) < ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              ij = 0
              DO j = 1, 3
                DO i = 1, 3
                  ij = ij + 1
                  tdf_sigma(ij) = vkk_all(i, ibnd, ik) * F_out(j, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibnd, ik) < ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        av_mob(itemp) = mobility
      ENDDO
      ! 
    ! Now electron mobility
    ENDIF
    IF (ncarrier > 1E5 ) THEN
      ! Needed because of residual values from the hole above
      Sigma(:,:)    = zero
      ! 
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        DO ik=1,  nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              tdf_sigma(:) = zero
              ij = 0
              DO j = 1, 3
                DO i = 1, 3
                  ij = ij + 1
                  tdf_sigma(ij) = vkk_all(i, ibnd, ik) * F_out(j, ibnd, ik, itemp) * wkf_all(ik)
                ENDDO
              ENDDO
              ! 
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              !  
              ! derivative Fermi distribution
              ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
              dfnk = w0gauss( ekk / etemp, -99 ) / etemp
              !
              ! electrical conductivity
              Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:)
              !
            ENDIF ! if below Fermi level
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp      
      ! 
      DO itemp=1, nstemp
        etemp = transp_temp(itemp)
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for electron conduction
            IF (etf_all (ibnd, ik) > ef0(itemp) ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ik) * fnk
            ENDIF
          ENDDO
        ENDDO
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
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm)**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density = carrier_density * inv_cell * ( bohr2ang * ang2cm)**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev/kelvin2eV, &
                      ef0(itemp)*ryd2ev,  carrier_density, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        av_mob(itemp) = mobility
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
    ! Adapted from the code PH/print_clock_ph - Quantum-ESPRESSO group          
    !-----------------------------------------------------------------------
    SUBROUTINE print_clock_epw
    !-----------------------------------------------------------------------
    !  
    !  12/2009 The clock data could be better organized
    !  
    USE io_global,     ONLY : stdout
    USE uspp,          ONLY : nlcc_any
    !  
    implicit none
    !
    WRITE( stdout, * )
    WRITE( stdout,  * ) '    Unfolding on the coarse grid'
    CALL print_clock ('dvanqq2')
    CALL print_clock ('elphon_wrap')
    WRITE( stdout, * )
    CALL print_clock ('ELPHWAN')
    WRITE( stdout,  * ) '    INITIALIZATION: '
    CALL print_clock ('epq_init')
    WRITE( stdout, * )
    CALL print_clock ('epq_init')
    IF (nlcc_any) call print_clock ('set_drhoc')
    CALL print_clock ('init_vloc')
    CALL print_clock ('init_us_1')
    CALL print_clock ('newd')
    CALL print_clock ('dvanqq')
    CALL print_clock ('drho')
    WRITE( stdout, * )
    !
    ! Electron-Phonon interpolation 
    !
    WRITE( stdout, * )
    WRITE( stdout,  * ) '    Electron-Phonon interpolation'   
    CALL print_clock ('ephwann')
    CALL print_clock ('ep-interp')
    CALL print_clock('PH SELF-ENERGY')
    CALL print_clock('ABS SPECTRA')
    CALL print_clock ('crys_cart')
    WRITE( stdout, * )
    CALL print_clock ('load data')
    CALL print_clock ('Ham: step 1')
    CALL print_clock ('Ham: step 2')
    CALL print_clock ('Ham: step 3')
    CALL print_clock ('Ham: step 4')
    CALL print_clock ('ep: step 1')
    CALL print_clock ('ep: step 2')
    CALL print_clock ('ep: step 3')
    CALL print_clock ('ep: step 4')
    CALL print_clock ('DynW2B')
    CALL print_clock ('HamW2B')
    CALL print_clock ('ephW2Bp')
    CALL print_clock ('ephW2Bp1')
    CALL print_clock ('ephW2Bp2') 
    !
    ! Eliashberg equations
    WRITE( stdout, * )
    CALL print_clock( 'ELIASHBERG' )
    !
    WRITE( stdout, * )
    WRITE( stdout,  * ) '    Total program execution'
    CALL print_clock ('EPW') 
    !
    END SUBROUTINE print_clock_epw
    ! 
  END MODULE printing
