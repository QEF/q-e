  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !       
  !-----------------------------------------------------------------------
  SUBROUTINE transport_coeffs (ef0,efcb)
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine computes the transport coefficients
  !!
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, meta_ionode_id
  USE cell_base, ONLY : alat, at, omega
  USE io_files,  ONLY : prefix 
  USE io_epw,    ONLY : iufilsigma 
  USE epwcom,    ONLY : nbndsub, fsthick, & 
                        system_2d, nstemp, &
                        int_mob, ncarrier, scatread, &
                        iterative_bte
  USE pwcom,     ONLY : ef 
  USE elph2,     ONLY : ibndmax, ibndmin, etf, nkf, wkf, dmef, & 
                        inv_tau_all, nkqtotf, Fi_all, inv_tau_allcb, &
                        zi_allvb, zi_allcb
  USE transportcom,  ONLY : transp_temp
  USE constants_epw, ONLY : zero, one, bohr2ang, ryd2ev, electron_SI, &
                            kelvin2eV, hbar, Ang2m, hbarJ, ang2cm, czero
  USE mp,        ONLY : mp_sum
  USE mp_global, ONLY : world_comm
  USE mp_world,  ONLY : mpime
  !
  IMPLICIT NONE
  ! 
  REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp (could be 0)
  !
  ! Local variables
  INTEGER :: i
  !! Cartesian direction index 
  INTEGER :: j
  !! Cartesian direction index 
  INTEGER :: ij
  !! Cartesian coupled index for matrix. 
  INTEGER :: ik
  !! K-point index
  INTEGER :: ikk
  !! Odd index to read etf
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: itemp
  !! Temperature index
  INTEGER :: lower_bnd
  !! Lower bounds index after k or q paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k or q paral
  ! 
  REAL(KIND=DP) :: ekk
  !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
 ! REAL(KIND=DP) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND=DP) :: dfnk
  !! Derivative Fermi distribution $$-df_{nk}/dE_{nk}$$
  REAL(KIND=DP) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND=DP) :: tau 
  !! Relaxation time
  REAL(KIND=DP) :: conv_factor1
  !! Conversion factor for the conductivity 
  REAL(KIND=DP) :: inv_cell 
  !! Inverse of the volume in [Bohr^{-3}]
!  REAL(KIND=DP) :: tdf_factor(9)
  !! Transport distribution function factor 
  REAL(KIND=DP) :: carrier_density
  !! Carrier density [nb of carrier per unit cell]
  REAL(KIND=DP) :: carrier_density_prt
  !! Carrier density [nb of carrier per unit cell] in cm^-3 unit
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
  REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
  !! Electron velocity vector for a band. 
  REAL(KIND=DP) :: Sigma(9,nstemp)
  !! Conductivity matrix in vector form
  REAL(KIND=DP) :: SigmaZ(9,nstemp)
  !! Conductivity matrix in vector form with Znk
  REAL(KIND=DP) :: Sigma_m(3,3,nstemp)
  !! Conductivity matrix
  REAL(KIND=DP) :: sigma_up(3,3)
  !! Conductivity matrix in upper-triangle
  REAL(KIND=DP) :: sigma_eig(3)
  !! Eigenvalues from the diagonalized conductivity matrix
  REAL(KIND=DP) :: sigma_vect(3,3)
  !! Eigenvectors from the diagonalized conductivity matrix
  REAL(KIND=DP) :: Znk
  !! Real Znk from \lambda_nk (called zi_allvb or zi_allcb)
  REAL(KIND=DP) :: tdf_sigma(9)
  !! Temporary file
  REAL(kind=DP), PARAMETER :: eps = 1.d-4
  !! Tolerence
  REAL(KIND=DP), EXTERNAL :: wgauss
  !! Compute the approximate theta function. Here computes Fermi-Dirac 
  REAL(KIND=DP), EXTERNAL :: w0gauss
  !! The derivative of wgauss:  an approximation to the delta function 
  REAL(KIND=DP), EXTERNAL :: efermig
  !! Function that returns the Fermi energy
  CHARACTER (len=256) :: filsigma
  !! File for the conductivity  
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  COMPLEX(kind=DP), ALLOCATABLE :: dmef_all(:,:,:,:)
  !! dipole matrix elements on the fine mesh among all pools
  REAL(DP), ALLOCATABLE :: tdf_sigma_m(:,:,:,:)
  !! transport distribution function
  REAL(DP), ALLOCATABLE :: wkf_all(:)
  !! k-point weight on the full grid across all pools
  !
  inv_cell = 1.0d0/omega
  ! for 2d system need to divide by area (vacuum in z-direction)
  IF ( system_2d ) &
     inv_cell = inv_cell * at(3,3) * alat

  ! 
  ! We can read the scattering rate from files. 
  IF ( scatread ) THEN
    conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
    !
    ! Compute the Fermi level 
    DO itemp = 1, nstemp
      ! 
      etemp = transp_temp(itemp)
      ! 
      ! Lets gather the velocities from all pools
#ifdef __MPI
      IF ( .not. ALLOCATED(dmef_all) )  ALLOCATE( dmef_all(3,nbndsub,nbndsub,nkqtotf) )
      IF ( .not. ALLOCATED(wkf_all) )  ALLOCATE( wkf_all(nkqtotf) )
      wkf_all(:) = zero
      dmef_all(:,:,:,:) = czero
      CALL poolgather2 ( 1, nkqtotf, 2*nkf, wkf, wkf_all  )
      CALL poolgatherc4 ( 3, nbndsub, nbndsub, nkqtotf, 2*nkf, dmef, dmef_all )
#else
      dmef_all = dmef
#endif     
      ! 
      ! In this case, the sum over q has already been done. It should therefore be ok 
      ! to do the mobility in sequential. Each cpu does the same thing below
      ALLOCATE ( etf_all ( nbndsub, nkqtotf/2 ) )
      !
      CALL scattering_read(etemp, ef0(itemp), etf_all, inv_tau_all)
      ! 
      ! This is hole mobility. ----------------------------------------------------
      IF (int_mob .OR. (ncarrier < -1E5)) THEN
        IF (itemp == 1) THEN
          WRITE(stdout,'(/5x,a)') repeat('=',67)
          WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
          WRITE(stdout,'(5x,a/)') repeat('=',67)
        ENDIF
        !      
        IF ( itemp .eq. 1 ) THEN        
          IF ( .not. ALLOCATED(tdf_sigma_m) )  ALLOCATE( tdf_sigma_m(3,3,ibndmax-ibndmin+1,nkqtotf) )
          tdf_sigma_m(:,:,:,:) = zero
          Sigma_m(:,:,:)   = zero
        ENDIF
        !
        DO ik = 1, nkqtotf/2 
          !DBSP
          !write(*,*)'ik ',ik
          !write(*,*)'SUM(inv_tau_all) ',SUM(inv_tau_all(:,:,ik))
          !write(*,*)'Sigma_m(:) before ',SUM(Sigma_m)
          !write(*,*)'minval ( abs(etf_all (:, ik) - ef ) )',minval ( abs(etf_all(:, ik) - ef ) )
          !write(*,*)'fsthick ',fsthick
          ikk = 2 * ik - 1
          ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
          IF ( minval ( abs(etf_all (:, ik) - ef ) ) < fsthick ) THEN
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only valence bands for hole conduction
              IF (etf_all (ibndmin-1+ibnd, ik) < ef0(itemp)  ) THEN
                vkk(:,ibnd) = 2.0 * REAL (dmef_all (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                ! We take itemp = 1 only !!!!
                tau = one / inv_tau_all(1,ibnd,ik)
                ekk = etf_all (ibndmin-1+ibnd, ik) -  ef0(itemp)
                ! 
                DO j = 1, 3
                  DO i = 1, 3
                    tdf_sigma_m(i,j,ibnd,ik) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                  ENDDO
                ENDDO
                !
                ! derivative Fermi distribution
                dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                !
                ! electrical conductivity matrix
                Sigma_m(:,:,itemp) = Sigma_m(:,:,itemp) +  wkf_all(ikk) * dfnk * tdf_sigma_m(:,:,ibnd,ik)
                !
              ENDIF ! valence bands
            ENDDO ! ibnd
          ENDIF ! fstick
          !write(*,*)'Sigma_m(:) ',SUM(Sigma_m), 'ef ',ef, 'ef0 ',ef0
        ENDDO ! ik
        ! 
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          ikk = 2 * ik - 1
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only valence bands for hole conduction
            IF (etf_all (ibndmin-1+ibnd, ik) < ef0(itemp)  ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibndmin-1+ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ikk) * (1.0d0 - fnk )
            ENDIF
          ENDDO
        ENDDO
        ! 
        ! Diagonalize the conductivity matrix
        CALL rdiagh(3,Sigma_m(:,:,itemp),3,sigma_eig,sigma_vect)
        ! 
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        !
      ENDIF ! int_mob .OR. (ncarrier < -1E5)
      ! 
      ! This is electron mobility. ----------------------------------------------------
      IF (int_mob .OR. (ncarrier > 1E5)) THEN
        IF (itemp == 1) THEN
          WRITE(stdout,'(/5x,a)') repeat('=',67)
          WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Electron density [cm^-3]  Electron mobility [cm^2/Vs]")')
          WRITE(stdout,'(5x,a/)') repeat('=',67)
        ENDIF
        !      
        IF ( itemp .eq. 1 ) THEN
          IF ( .not. ALLOCATED(tdf_sigma_m) )  ALLOCATE( tdf_sigma_m(3,3,ibndmax-ibndmin+1,nkqtotf) )
          tdf_sigma_m(:,:,:,:) = zero
          Sigma_m(:,:,:)   = zero
        ENDIF
        !
        DO ik = 1, nkqtotf/2
          ikk = 2 * ik - 1
          ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
          IF ( minval ( abs(etf_all (:, ik) - ef ) ) < fsthick ) THEN
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only conduction bands for electron conduction
              IF (etf_all (ibndmin-1+ibnd, ik) > ef0(itemp)  ) THEN
                vkk(:,ibnd) = 2.0 * REAL (dmef_all (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                tau = one / inv_tau_all(1,ibnd,ik)
                ekk = etf_all (ibndmin-1+ibnd, ik) -  ef0(itemp)
                ! 
                DO j = 1, 3
                  DO i = 1, 3
                    tdf_sigma_m(i,j,ibnd,ik) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                  ENDDO
                ENDDO
                !
                ! derivative Fermi distribution
                dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                !
                ! electrical conductivity matrix
                Sigma_m(:,:,itemp) = Sigma_m(:,:,itemp) +  wkf_all(ikk) * dfnk * tdf_sigma_m(:,:,ibnd,ik)
                !
              ENDIF ! valence bands
            ENDDO ! ibnd
          ENDIF ! fstick
        ENDDO ! ik
        ! 
        carrier_density = 0.0
        ! 
        DO ik = 1, nkqtotf/2
          ikk = 2 * ik - 1
          DO ibnd = 1, ibndmax-ibndmin+1
            ! This selects only conduction bands for electron conduction
            IF (etf_all (ibndmin-1+ibnd, ik) > ef0(itemp)  ) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all (ibndmin-1+ibnd, ik) - ef0(itemp)
              fnk = wgauss( -ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density = carrier_density + wkf_all(ikk) * fnk 
            ENDIF
          ENDDO
        ENDDO
        ! 
        ! Diagonalize the conductivity matrix
        CALL rdiagh(3,Sigma_m(:,:,itemp),3,sigma_eig,sigma_vect)
        ! 
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        !
      ENDIF ! int_mob .OR. (ncarrier > 1E5)
      ! 
    ENDDO ! itemp
    !
  ELSE ! Case without reading the scattering rates from files.
    !
    ! This is hole mobility. In the case of intrinsic mobilities we can do both
    ! electron and hole mobility because the Fermi level is the same. This is not
    ! the case for doped mobilities.
    ! 
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds( nkqtotf/2, lower_bnd, upper_bnd )
    !DBSP
    !print*,'inv_tau_all ',SUM(inv_tau_all(:,:,:)) 
    !print*,'zi_allvb ',SUM(zi_allvb(:,:,:)) 
    ! 
    IF (int_mob .OR. (ncarrier < -1E5)) THEN
      ! 
      DO itemp = 1, nstemp
        !
        etemp = transp_temp(itemp)
         
        !DBSP
        !write(stdout,*)'etemp ',etemp 
        !write(stdout,*)'inv_tau_all ', SUM(inv_tau_all(itemp,:,:))
        !write(stdout,*)'inv_tau_all ', SUM(inv_tau_all(:,:,:))
        !
        IF ( itemp .eq. 1 ) THEN 
          !
          ! tdf_sigma_ij(ibnd,ik) = v_i(ik,ibnd) * v_j(ik,ibnd) * tau(ik,ibnd)
          ! i,j - cartesian components and ij combined (i,j) index
          ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
          ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
          ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
          ! this can be reduced to 6 if we take into account symmetry xy=yx, ...
          tdf_sigma(:)  = zero
          Sigma(:,:)    = zero
          SigmaZ(:,:)   = zero
          !
        ENDIF
        !
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          !
          ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
          IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
            !
            ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
            ! 1/m  = 2 in Rydberg atomic units
            ! dmef is in units of 1/a.u. (where a.u. is bohr)
            ! v_(k,i) is in units of Rydberg * a.u.
            !
            DO ibnd = 1, ibndmax-ibndmin+1
              !
              ! This selects only valence bands for hole conduction
              IF (etf (ibndmin-1+ibnd, ikk) < ef0(itemp) ) THEN 
                !
                ! vkk(3,nbnd) - velocity for k
                vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                !
                !  energy at k (relative to Ef)
                ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                !
                tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
                !
                ij = 0
                DO j = 1, 3
                  DO i = 1, 3
                    ij = ij + 1
                    tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                  ENDDO
                ENDDO
                !
                ! derivative Fermi distribution
                ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                !
                ! electrical conductivity
                Sigma(:,itemp) = Sigma(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
                !
                ! Now do the same but with Znk multiplied
                ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                Znk = one / ( one + zi_allvb (itemp,ibnd,ik+lower_bnd-1) )
                tau = one / ( Znk * inv_tau_all(itemp,ibnd,ik+lower_bnd-1) )
                ij = 0
                DO j = 1, 3
                  DO i = 1, 3
                    ij = ij + 1
                    tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                  ENDDO
                ENDDO
                SigmaZ(:,itemp) = SigmaZ(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
 
                !print*,'itemp ik ibnd ',itemp, ik, ibnd
                !print*,'Sigma ',Sigma(:,itemp)
                !print*,'SigmaZ ',SigmaZ(:,itemp)
                !print*,'Znk ',Znk
                !
              ENDIF
              !
            ENDDO ! ibnd
            !
          ENDIF ! endif  fsthick
          !
        ENDDO ! end loop on k
        !
        ! The k points are distributed among pools: here we collect them
        !
        CALL mp_sum( Sigma(:,itemp),  world_comm )
        CALL mp_sum( SigmaZ(:,itemp), world_comm )
        !DBSP
        !write(stdout,*) 'ef0(itemp) ',ef0(itemp)    
        !write(stdout,*) 'Sigma ',SUM(Sigma(:,itemp))    
        !
      ENDDO ! nstemp
      !
      IF (mpime .eq. meta_ionode_id) THEN
        filsigma = TRIM(prefix) // '_elcond_h'
        OPEN(iufilsigma, file = filsigma, form = 'formatted')
        WRITE(iufilsigma,'(a)') "# Electrical conductivity in 1/(Ohm * m)"
        WRITE(iufilsigma,'(a)') "#         Ef(eV)         Temp(K)        Sigma_xx        Sigma_xy        Sigma_xz" // & 
                                                 "       Sigma_yx         Sigma_yy        Sigma_yz " // &
                                                 "        Sigma_xz        Sigma_yz        Sigma_zz"
      ENDIF
      !
      conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
      !
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      ! 
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        ! Sigma in units of 1/(a.u.) is converted to 1/(Ohm * m)
        IF (mpime.eq. meta_ionode_id) THEN
          WRITE(iufilsigma,'(11E16.8)') ef0(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                       conv_factor1 * Sigma(:,itemp) * inv_cell
        ENDIF
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
        CALL mp_sum( carrier_density, world_comm )
        !
        ! Diagonalize the conductivity matrix
        ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
        ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
        ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
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
        CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
        ! 
        !Sigma_diag = (Sigma(1,itemp)+Sigma(5,itemp)+Sigma(9,itemp))/3
        !Sigma_offdiag = (Sigma(2,itemp)+Sigma(3,itemp)+Sigma(4,itemp)+&
        !                 Sigma(6,itemp)+Sigma(7,itemp)+Sigma(8,itemp))/6
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
        carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg' 
        ! 
        ! Now do Znk ----------------------------------------------------------
        sigma_up(:,:) = zero
        sigma_up(1,1) = SigmaZ(1,itemp)
        sigma_up(1,2) = SigmaZ(2,itemp)
        sigma_up(1,3) = SigmaZ(3,itemp)
        sigma_up(2,1) = SigmaZ(4,itemp)
        sigma_up(2,2) = SigmaZ(5,itemp)
        sigma_up(2,3) = SigmaZ(6,itemp)
        sigma_up(3,1) = SigmaZ(7,itemp)
        sigma_up(3,2) = SigmaZ(8,itemp)
        sigma_up(3,3) = SigmaZ(9,itemp)
        CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        ! carrier_density in cm^-1
! DBSP - Z-factor
!        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
!                ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
!        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis [Z]'
!        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis [Z]'
!        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg [Z]'

        ! 
      ENDDO ! nstemp
      !
      IF (mpime .eq. meta_ionode_id) CLOSE(iufilsigma)
      !
    ENDIF ! Hole mob
    ! 
    ! Now the electron conduction and mobilities
    ! 
    IF (int_mob .OR. (ncarrier > 1E5)) THEN
      DO itemp = 1, nstemp
        !
        etemp = transp_temp(itemp)
        IF ( itemp .eq. 1 ) THEN
          tdf_sigma(:)  = zero
          Sigma(:,:)    = zero
          SigmaZ(:,:)   = zero
        ENDIF
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
            IF ( ABS(efcb(itemp)) < eps ) THEN  
              DO ibnd = 1, ibndmax-ibndmin+1
                ! This selects only cond bands for electron conduction
                IF (etf (ibndmin-1+ibnd, ikk) > ef0(itemp) ) THEN
                  vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    ENDDO
                  ENDDO
                  dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                  Sigma(:,itemp) = Sigma(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
                  !
                  ! Now do the same but with Znk multiplied
                  ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                  Znk = one / ( one + zi_allvb (itemp,ibnd,ik+lower_bnd-1) )
                  tau = one / ( Znk * inv_tau_all(itemp,ibnd,ik+lower_bnd-1) )
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    ENDDO
                  ENDDO
                  SigmaZ(:,itemp) = SigmaZ(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)
                ENDIF
              ENDDO 
            ELSE ! In this case we have 2 Fermi levels
              DO ibnd = 1, ibndmax-ibndmin+1
                ! This selects only cond bands for hole conduction
                IF (etf (ibndmin-1+ibnd, ikk) > efcb(itemp) ) THEN
                  vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                  tau = one / inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1)
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    ENDDO
                  ENDDO
                  dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                  Sigma(:,itemp) = Sigma(:,itemp) +  wkf(ikk) * dfnk * tdf_sigma(:)
                  !
                  ! Now do the same but with Znk multiplied
                  ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                  Znk = one / ( one + zi_allcb (itemp,ibnd,ik+lower_bnd-1) )
                  tau = one / ( Znk * inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1) )
                  ij = 0
                  DO j = 1, 3
                    DO i = 1, 3
                      ij = ij + 1
                      tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    ENDDO
                  ENDDO
                  SigmaZ(:,itemp) = SigmaZ(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)                  
                ENDIF 
              ENDDO ! ibnd
            ENDIF ! etcb
          ENDIF ! endif  fsthick
        ENDDO ! end loop on k
        CALL mp_sum( Sigma(:,itemp),  world_comm )
        CALL mp_sum( SigmaZ(:,itemp), world_comm )
        ! 
      ENDDO ! nstemp
      IF (mpime .eq. meta_ionode_id) THEN
        filsigma = TRIM(prefix) // '_elcond_e'
        OPEN(iufilsigma, file = filsigma, form = 'formatted')
        WRITE(iufilsigma,'(a)') "# Electrical conductivity in 1/(Ohm * m)"
        WRITE(iufilsigma,'(a)') "#         Ef(eV)         Temp(K)        Sigma_xx        Sigma_xy        Sigma_xz" // &
                                                 "       Sigma_yx         Sigma_yy        Sigma_yz " // &
                                                "        Sigma_xz        Sigma_yz        Sigma_zz"
      ENDIF
      !
      conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      DO itemp = 1, nstemp
        etemp = transp_temp(itemp)
        IF (mpime .eq. meta_ionode_id) THEN
          ! Sigma in units of 1/(a.u.) is converted to 1/(Ohm * m)
          IF ( ABS(efcb(itemp)) < eps ) THEN 
            WRITE(iufilsigma,'(11E16.8)') ef0(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                         conv_factor1 * Sigma(:,itemp) * inv_cell
          ELSE
            WRITE(iufilsigma,'(11E16.8)') efcb(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                         conv_factor1 * Sigma(:,itemp) * inv_cell
          ENDIF
        ENDIF
        carrier_density = 0.0
        ! 
        DO ik = 1, nkf
          DO ibnd = 1, ibndmax-ibndmin+1
            ikk = 2 * ik - 1
            ! This selects only conduction bands for electron conduction
            IF ( ABS(efcb(itemp)) < eps ) THEN 
              IF (etf (ibndmin-1+ibnd, ikk) > ef0(itemp) ) THEN
                ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                fnk = wgauss( -ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf(ikk) * fnk
              ENDIF
            ELSE
              IF (etf (ibndmin-1+ibnd, ikk) > efcb(itemp) ) THEN
                ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                fnk = wgauss( -ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf(ikk) * fnk
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        CALL mp_sum( carrier_density, world_comm )
        ! Diagonalize the conductivity matrix
        ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
        ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
        ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
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
        CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
        ! 
        !Sigma_diag = (Sigma(1,itemp)+Sigma(5,itemp)+Sigma(9,itemp))/3
        !Sigma_offdiag = (Sigma(2,itemp)+Sigma(3,itemp)+Sigma(4,itemp)+&
        !                 Sigma(6,itemp)+Sigma(7,itemp)+Sigma(8,itemp))/6
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        !
     
        ! carrier_density in cm^-1
        carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
        IF ( ABS(efcb(itemp)) < eps ) THEN
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                   ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
        ELSE
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                   efcb(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
        ENDIF
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
        ! Issue warning if the material is anisotropic
       ! IF (Sigma_offdiag > 0.1*Sigma_diag) THEN
       !   WRITE(stdout,'(5x,a,1f10.5,a)') 'Warning: Sigma_offdiag = ',(Sigma_offdiag*100)/Sigma_diag, '% of Sigma_diag'
       ! ENDIF
        ! Now do the mobility with Znk factor ----------------------------------------------------------
        sigma_up(:,:) = zero
        sigma_up(1,1) = SigmaZ(1,itemp)
        sigma_up(1,2) = SigmaZ(2,itemp)
        sigma_up(1,3) = SigmaZ(3,itemp)
        sigma_up(2,1) = SigmaZ(4,itemp)
        sigma_up(2,2) = SigmaZ(5,itemp)
        sigma_up(2,3) = SigmaZ(6,itemp)
        sigma_up(3,1) = SigmaZ(7,itemp)
        sigma_up(3,2) = SigmaZ(8,itemp)
        sigma_up(3,3) = SigmaZ(9,itemp)
        CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
        !
! DBSP - Z-factor
!        IF ( ABS(efcb(itemp)) < eps ) THEN
!          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
!                                                   ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
!        ELSE
!          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
!                                                   efcb(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
!        ENDIF
!        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis [Z]'
!        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis [Z]'
!        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg [Z]'

        ! 
      ENDDO ! nstemp
      WRITE(stdout,'(5x)')
      WRITE(stdout,'(5x,"Note: Mobility are sorted by ascending values and might not correspond to the expected (x,y,z) axis.")')
      !
      IF (mpime .eq. meta_ionode_id) CLOSE(iufilsigma)
      ! 
    ENDIF ! Electron mobilities
  ENDIF ! scatread
  ! 
  ! IF IBTE we want the SRTA solution to be the first iteration of IBTE
  IF (iterative_bte) THEN
    Fi_all(:,:,:) = zero
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
        DO ibnd = 1, ibndmax-ibndmin+1
          vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd,ikk))
          tau = one / inv_tau_all(1,ibnd,ik+lower_bnd-1)
          Fi_all(:,ibnd,ik+lower_bnd-1) = vkk(:,ibnd) * tau
        ENDDO
      ENDIF
    ENDDO ! kpoints
    CALL mp_sum( Fi_all, world_comm )
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE transport_coeffs
  !--------------------------------------------------------------------------
  ! 
