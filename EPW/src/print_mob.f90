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
  SUBROUTINE print_serta( F_SERTA, BZtoIBZ, s_BZtoIBZ, vkk_all, etf_all, wkf_all, ef0) 
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine prints the mobility
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
  INTEGER, INTENT(in) :: BZtoIBZ(nkf1*nkf2*nkf3)
  !! BZ to IBZ mapping
  INTEGER, INTENT(in) :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
  !! Corresponding symmetry matrix
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
            ! Loop on full BZ
            nb = 0
            DO ikbz=1, nkf1*nkf2*nkf3
              IF (BZtoIBZ(ikbz) == ik) THEN
                nb = nb + 1
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
            ! Loop on full BZ
            nb = 0
            DO ikbz=1, nkf1*nkf2*nkf3
              IF (BZtoIBZ(ikbz) == ik) THEN
                nb = nb + 1
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
  END SUBROUTINE print_serta
  !-----------------------------------------------------------------------
  !  
  !-----------------------------------------------------------------------
  SUBROUTINE print_mob( F_out, BZtoIBZ, s_BZtoIBZ, vkk_all, etf_all, wkf_all, ef0, av_mob) 
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
  INTEGER, INTENT(in) :: BZtoIBZ(nkf1*nkf2*nkf3)
  !! BZ to IBZ mapping
  INTEGER, INTENT(in) :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
  !! Corresponding symmetry matrix
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
            ! Loop on full BZ
            nb = 0
            DO ikbz=1, nkf1*nkf2*nkf3
              IF (BZtoIBZ(ikbz) == ik) THEN
                nb = nb + 1
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
            ! Loop on full BZ
            nb = 0
            DO ikbz=1, nkf1*nkf2*nkf3
              IF (BZtoIBZ(ikbz) == ik) THEN
                nb = nb + 1
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
  END SUBROUTINE print_mob
  !-----------------------------------------------------------------------
