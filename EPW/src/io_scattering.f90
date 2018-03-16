  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------------
  SUBROUTINE F_write(iter, iq, nqtotf, nktotf, error_h, error_el)
  !----------------------------------------------------------------------------
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE io_epw,       ONLY : iufilFi_all
  USE io_files,     ONLY : diropn
  USE mp,           ONLY : mp_barrier
  USE mp_global,    ONLY : inter_pool_comm
  USE mp_world,     ONLY : mpime
  USE io_global,    ONLY : ionode_id
  USE elph2,        ONLY : F_current, ibndmax, ibndmin
  USE transportcom, ONLY : lower_bnd, upper_bnd, mobilityh_save, mobilityel_save
  USE constants_epw, ONLY : upper_bnd < nktotf ) F_current(:,:,upper_bnd+1:nktotf) = zero
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE F_write
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE F_read(iter, iq, nqtotf, nktotf, error_h, error_el)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilFi_all
  USE constants_epw, ONLY : zero
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm, root_pool
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  USE elph2,        ONLY : F_current, ibndmax, ibndmin, Fi_all
  USE transportcom, ONLY : lower_bnd, upper_bnd, mobilityh_save, mobilityel_save
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter
  !! Iteration number
  INTEGER, INTENT(INOUT) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  REAL(kind=DP), INTENT(INOUT) :: error_h
  !! Error in the hole mobility
  REAL(kind=DP), INTENT(INOUT) :: error_el
  !! Error in the electron mobility

  !
  ! Local variable
  LOGICAL :: exst
  !! 
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
  INTEGER :: nqtotf_read
  !! Total number of q-point read
  ! 
  CHARACTER (len=256) :: name1
 
  REAL(KIND=DP) :: aux ( 3 * (ibndmax-ibndmin+1) * nktotf + 7 )
  !! Vector to store the array
  !
  IF (mpime.eq.ionode_id) THEN
    ! 
    ! First inquire if the file exists
#if defined(__MPI)
    name1 = trim(tmp_dir) // trim(prefix) // '.F_restart1'
#else
    name1 = trim(tmp_dir) // trim(prefix) // '.F_restart'
#endif
    INQUIRE(file = name1, exist=exst)
    ! 
    IF (exst) THEN ! read the file
      !
      lfi_all = 3* (ibndmax-ibndmin+1) * nktotf +3
      CALL diropn (iufilFi_all, 'F_restart', lfi_all, exst)
      CALL davcio ( aux, lfi_all, iufilFi_all, 1, -1 )
      !
      ! First element is the iteration number
      iter = aux(1)
      ! Current iteration number
      iq = aux(2)
      iq = iq + 1 ! we need to start at the next q
      ! Total number of q-points
      nqtotf_read = aux(3)
      ! Last value of hole mobility
      mobilityh_save = aux(4) 
      ! Last value of electron mobility
      mobilityel_save = aux(5) 
      ! Error in hole mobility
      error_h = aux(6) 
      ! Error in electron mobility
      error_el = aux(7) 
      ! This is the error of the previous iteration. Therefore when you restart
      ! from a converged one, you want to be finished. 
      ! This small substraction wont affect anything. 
      error_h = error_h -0.5E-2
      error_el = error_el -0.5E-2
      !
      IF ( nqtotf_read /= nqtotf) CALL errore('io_scattering',&
        &'Error: The current total number of q-point is not the same as the read one. ',1)
      !
      i = 7
      DO ik=1, nktotf
        DO ibnd=1, (ibndmax-ibndmin+1)
          DO idir=1,3
            i = i +1
            F_current(idir,ibnd, ik) = aux(i)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iufilFi_all)
    ENDIF
  ENDIF
  ! 
  CALL mp_bcast (exst, ionode_id, inter_pool_comm)
  CALL mp_bcast (exst, root_pool, intra_pool_comm)
  !
  IF (exst) THEN
    CALL mp_bcast (iter, ionode_id, inter_pool_comm)
    CALL mp_bcast (iter, root_pool, intra_pool_comm)
    CALL mp_bcast (iq, ionode_id, inter_pool_comm)
    CALL mp_bcast (iq, root_pool, intra_pool_comm)
    CALL mp_bcast (F_current, ionode_id, inter_pool_comm)
    CALL mp_bcast (F_current, root_pool, intra_pool_comm)
    CALL mp_bcast (mobilityh_save, ionode_id, inter_pool_comm)
    CALL mp_bcast (mobilityh_save, root_pool, intra_pool_comm)
    CALL mp_bcast (mobilityel_save, ionode_id, inter_pool_comm)
    CALL mp_bcast (mobilityel_save, root_pool, intra_pool_comm)
    CALL mp_bcast (error_h, ionode_id, inter_pool_comm)
    CALL mp_bcast (error_h, root_pool, intra_pool_comm)
    CALL mp_bcast (error_el, ionode_id, inter_pool_comm)
    CALL mp_bcast (error_el, root_pool, intra_pool_comm)
    ! 
    ! The Fi on the commensurate full k-grid is equal to the full
    ! F_current from the previous read run.
    Fi_all = F_current
    CALL mp_bcast (Fi_all, ionode_id, inter_pool_comm)
    CALL mp_bcast (Fi_all, root_pool, intra_pool_comm)
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1 ) F_current(:,:,1:lower_bnd-1) = zero
    IF (upper_bnd < nktotf ) F_current(:,:,upper_bnd+1:nktotf) = zero
    ! 
    WRITE(stdout, '(a,i10,a,i10,a,i10)' ) '     Restart from iter: ',iter,' and iq: ',iq,'/',nqtotf
  ENDIF
  !
  !----------------------------------------------------------------------------
  END SUBROUTINE F_read
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE scattering_write(itemp, etemp, ef0, etf_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, inv_tau_all
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
  INTEGER, INTENT(IN) :: itemp
  !! Temperature index
  REAL(KIND=DP), INTENT(IN) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND=DP), INTENT(IN) :: etf_all(nbndsub, nkqtotf)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
!!  REAL(KIND=DP), INTENT(IN) :: inv_tau_all(nstemp,ibndmax-ibndmin+1, nkqtotf/2)
  !! Total scattering rate on all k-points collected from all pools in parallel case
  ! 
  ! Local variables
  INTEGER :: ik
  !! K-point index
  INTEGER :: ikk
  !! Odd index to read etf
  INTEGER :: ikq
  !! Even k+q index to read etf
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: nrec
  !! Record index
  REAL(KIND=DP) :: ekk
  !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
  REAL(KIND=DP) :: temp
  !! Temporary file name used to write scattering rate to file. 
  !
  CHARACTER (len=256) :: name1
  !! Name used to write scattering rates to file. 
  !
  WRITE(stdout,'(/5x,"Writing scattering rate to file"/)')
  !
  IF (mpime.eq.ionode_id) THEN
    !
    ! Write to file
    temp = etemp * ryd2ev / kelvin2eV
    IF ( temp .lt. 10.d0 - eps4 ) THEN
      WRITE(name1,'(a18,f4.2)') 'scattering_rate_00', temp
    ELSEIF ( temp .ge. 10.d0 - eps4 .AND. temp .lt. 100.d0 -eps4 ) THEN
      WRITE(name1,'(a17,f5.2)') 'scattering_rate_0', temp
    ELSEIF ( temp .ge. 100.d0 -eps4 ) THEN
      WRITE(name1,'(a16,f6.2)') 'scattering_rate_', temp
    ENDIF
    OPEN(iufilscatt_rate,file=name1, form='formatted')
    WRITE(iufilscatt_rate,'(a)') '# Inverse scattering time (ps)'
    WRITE(iufilscatt_rate,'(a)') '#      ik       ibnd                 E(ibnd)    scattering rate(1/ps)'
    !
    DO ik = 1, nkqtotf/2
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      DO ibnd = 1, ibndmax-ibndmin+1
        !
        ! note that ekk does not depend on q
        ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0(itemp)
        !
        WRITE(iufilscatt_rate,'(i9,2x)',advance='no') ik
        WRITE(iufilscatt_rate,'(i9,2x)',advance='no') ibndmin-1+ibnd
        WRITE(iufilscatt_rate,'(E22.14)',advance='no') ryd2ev * ekk
        !WRITE(iufilscatt_rate,'(E24.14)') ryd2mev * meV2invps * inv_tau_all(itemp,ibnd,ik)
        WRITE(iufilscatt_rate,'(E26.16E3)') ryd2mev * meV2invps * inv_tau_all(itemp,ibnd,ik)
        !
      ENDDO
      !
    ENDDO
    !
    CLOSE(iufilscatt_rate)
  ENDIF
  CALL mp_barrier(inter_pool_comm)
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE scattering_write
  !----------------------------------------------------------------------------
  ! 
  !----------------------------------------------------------------------------
  SUBROUTINE scattering_read(etemp, ef0, etf_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, inv_tau_all
  USE epwcom,    ONLY : nbndsub, nstemp
  USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                            meV2invps, eps4, zero
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : inter_pool_comm, root_pool, intra_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(IN) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND=DP), INTENT(IN) :: ef0
  !! Fermi level for the temperature itemp
  REAL(KIND=DP), INTENT(out) :: etf_all(bndsub, nkqtotf/2)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
!  REAL(KIND=DP), INTENT(out) :: inv_tau_all(nstemp,ibndmax-ibndmin+1, nkqtotf/2)
  !! Total scattering rate on all k-points collected from all pools in parallel case
  ! 
  ! Local variables
  INTEGER :: ik
  !! K-point index
  INTEGER :: ik_tmp
  !! K-point index read from file
  INTEGER :: ikk
  !! Odd index to read etf
  INTEGER :: ikq
  !! Even k+q index to read etf
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: ibnd_tmp
  !! Local band index read from file
  INTEGER :: nrec
  !! Record index
  INTEGER :: ios
  !! Status of reading file
  REAL(KIND=DP) :: ekk
  !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
  REAL(KIND=DP) :: temp
  !! Temporary file name used to write scattering rate to file. 
  !
  CHARACTER (len=256) :: name1
  !! Name used to write scattering rates to file. 
  CHARACTER (len=256) :: dummy1
  !! Dummy variable to store the text of the scattering_rate file 
  ! 
  WRITE(stdout,'(/5x,"Reading scattering rate from file"/)')
  !
  IF (mpime.eq.ionode_id) THEN
    ! Write to file
    temp = etemp * ryd2ev / kelvin2eV
    IF ( temp .lt. 10.d0 - eps4 ) THEN
      WRITE(name1,'(a18,f4.2)') 'scattering_rate_00', temp
    ELSEIF ( temp .ge. 10.d0 - eps4 .AND. temp .lt. 100.d0 -eps4 ) THEN
      WRITE(name1,'(a17,f5.2)') 'scattering_rate_0', temp
    ELSEIF ( temp .ge. 100.d0 -eps4 ) THEN
      WRITE(name1,'(a16,f6.2)') 'scattering_rate_', temp
    ENDIF
    OPEN(iufilscatt_rate,file=name1, status='old',iostat=ios)
    WRITE(stdout,'(a16,a22)'),'     Open file: ',name1   
    ! There are two comment line at the beginning of the file
    READ(iufilscatt_rate,*) dummy1
    READ(iufilscatt_rate,*) dummy1
    !
    DO ik = 1, nkqtotf/2
      !
      DO ibnd = 1, ibndmax-ibndmin+1
        !
        READ(iufilscatt_rate,*) ik_tmp, ibnd_tmp, &
             etf_all(ibndmin-1+ibnd, ik), inv_tau_all(1,ibnd,ik) 
        !etf_all (ibndmin-1+ibnd, ik) =  etf_all (ibndmin-1+ibnd, ik)/ryd2ev
        !etf_all (ibndmin-1+ibnd, ik) =  etf_all (ibndmin-1+ibnd, ik) + ef0
        inv_tau_all(1,ibnd,ik) = inv_tau_all(1,ibnd,ik)/(ryd2mev * meV2invps) 
     
        !
        ! Check that the file corresponds to the run we are making
        IF (ABS(ibnd_tmp - ibndmin - ibnd +1 ) > 0)  CALL errore('io_scattering', &
          'Band read from the scattering_rate file do not match current calculation ',1)
        ! 
      ENDDO
      ! Check that the file corresponds to the run we are making
      IF (ABS(ik_tmp-ik) > 0)  CALL errore('io_scattering', &
        'k-point read from the scattering_rate file do not match current calculation ',1)
      !
    ENDDO
    !
    etf_all = etf_all/ryd2ev
    etf_all = etf_all + ef0
    !
    CLOSE(iufilscatt_rate)
    ! 
  ENDIF
  CALL mp_bcast (etf_all, ionode_id, inter_pool_comm)
  CALL mp_bcast (etf_all, root_pool, intra_pool_comm)

  CALL mp_bcast (inv_tau_all, ionode_id, inter_pool_comm)
  CALL mp_bcast (inv_tau_all, root_pool, intra_pool_comm)

  CALL mp_barrier(inter_pool_comm)
  ! 
  WRITE(stdout,'(/5x,"Scattering rate read from file"/)')
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE scattering_read
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE electron_write(iq,nqtotf,nktotf,sigmar_all,sigmai_all,zi_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all
  USE io_epw,    ONLY : iufilsigma_all
  USE io_files,  ONLY : diropn
  USE epwcom,    ONLY : nstemp
  USE constants_epw, ONLY : kelvin2eV, &
                            meV2invps, eps4, zero
  USE transportcom, ONLY : lower_bnd, upper_bnd
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  ! Local variable
  LOGICAL :: exst
  !
  INTEGER, INTENT(IN) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  REAL(KIND=DP), INTENT(INOUT) :: sigmar_all(ibndmax-ibndmin+1, nktotf)
  !! Real part of the electron-phonon self-energy accross all pools
  REAL(KIND=DP), INTENT(INOUT) :: sigmai_all(ibndmax-ibndmin+1, nktotf)
  !! Imaginary part of the electron-phonon self-energy accross all pools
  REAL(KIND=DP), INTENT(INOUT) :: zi_all(ibndmax-ibndmin+1, nktotf)
  !! Z parameter of electron-phonon self-energy accross all pools
  ! 
  ! Local variables
  INTEGER :: i
  !! Iterative index
  INTEGER :: ik
  !! K-point index
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: lsigma_all
  !! Length of the vector
  REAL(KIND=DP) :: aux ( 3 * (ibndmax-ibndmin+1) * nktotf + 2 )
  !! Vector to store the array
  !
  IF (mpime.eq.ionode_id) THEN
    !
    lsigma_all = 3 * (ibndmax-ibndmin+1) * nktotf +2
    ! First element is the current q-point
    aux(1) = iq -1 ! we need to start at the next q
    ! Second element is the total number of q-points
    aux(2) = nqtotf
    !
    i = 2
    ! 
    DO ik=1, nktotf
      DO ibnd=1, (ibndmax-ibndmin+1)
        i = i +1
        aux(i) = sigmar_all(ibnd, ik)
      ENDDO
    ENDDO
    DO ik=1, nktotf
      DO ibnd=1, (ibndmax-ibndmin+1)
        i = i +1
        aux(i) = sigmai_all(ibnd, ik)
      ENDDO
    ENDDO
    DO ik=1, nktotf
      DO ibnd=1, (ibndmax-ibndmin+1)
        i = i +1
        aux(i) = zi_all(ibnd, ik)
      ENDDO
    ENDDO
    CALL diropn (iufilsigma_all, 'sigma_restart', lsigma_all, exst)
    CALL davcio ( aux, lsigma_all, iufilsigma_all, 1, +1 )
    CLOSE(iufilsigma_all)
  ENDIF
  ! 
  ! Make everythin 0 except the range of k-points we are working on
  IF (lower_bnd > 1 ) THEN 
    sigmar_all(:,1:lower_bnd-1) = zero
    sigmai_all(:,1:lower_bnd-1) = zero
    zi_all(:,1:lower_bnd-1) = zero
  ENDIF
  IF (upper_bnd < nktotf ) THEN
    sigmar_all(:,upper_bnd+1:nktotf) = zero
    sigmai_all(:,upper_bnd+1:nktotf) = zero
    zi_all(:,upper_bnd+1:nktotf) = zero
  ENDIF
  !CALL mp_barrier(inter_pool_comm)
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE electron_write
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE electron_read(iq,nqtotf,nktotf,sigmar_all,sigmai_all,zi_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all
  USE io_epw,    ONLY : iufilsigma_all
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE epwcom,    ONLY : nstemp
  USE constants_epw, ONLY : kelvin2eV, &
                            meV2invps, eps4, zero
  USE transportcom, ONLY : lower_bnd, upper_bnd
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm, root_pool
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  ! Local variable
  LOGICAL :: exst
  !
  INTEGER, INTENT(INOUT) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  REAL(KIND=DP), INTENT(OUT) :: sigmar_all(ibndmax-ibndmin+1, nktotf)
  !! Real part of the electron-phonon self-energy accross all pools
  REAL(KIND=DP), INTENT(OUT) :: sigmai_all(ibndmax-ibndmin+1, nktotf)
  !! Imaginary part of the electron-phonon self-energy accross all pools
  REAL(KIND=DP), INTENT(OUT) :: zi_all(ibndmax-ibndmin+1, nktotf)
  !! Z parameter of electron-phonon self-energy accross all pools
  ! 
  ! Local variables
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
  REAL(KIND=DP) :: aux ( 3 * (ibndmax-ibndmin+1) * nktotf + 2 )
  !! Vector to store the array
  ! 
  CHARACTER (len=256) :: name1
  !
  IF (mpime.eq.ionode_id) THEN
    !
    ! First inquire if the file exists
#if defined(__MPI)
    name1 = trim(tmp_dir) // trim(prefix) // '.sigma_restart1'
#else
    name1 = trim(tmp_dir) // trim(prefix) // '.sigma_restart'
#endif    
    INQUIRE(file = name1, exist=exst)
    ! 
    IF (exst) THEN ! read the file
      !
      lsigma_all = 3 * (ibndmax-ibndmin+1) * nktotf +2
      CALL diropn (iufilsigma_all, 'sigma_restart', lsigma_all, exst)
      CALL davcio ( aux, lsigma_all, iufilsigma_all, 1, -1 )
      !
      ! First element is the iteration number
      iq = aux(1)
      iq = iq + 1 ! we need to start at the next q
      nqtotf_read = aux(2)
      !print*, 'iq',iq
      !print*, 'nqtotf_read ',nqtotf_read
      IF ( nqtotf_read /= nqtotf) CALL errore('io_scattering',&
        &'Error: The current total number of q-point is not the same as the read one. ',1)
      ! 
      i = 2
      DO ik=1, nktotf
        DO ibnd=1, (ibndmax-ibndmin+1)
          i = i +1
          sigmar_all(ibnd, ik) = aux(i)
        ENDDO
      ENDDO
      DO ik=1, nktotf
        DO ibnd=1, (ibndmax-ibndmin+1)
          i = i +1
          sigmai_all(ibnd, ik) = aux(i)
        ENDDO
      ENDDO
      DO ik=1, nktotf
        DO ibnd=1, (ibndmax-ibndmin+1)
          i = i +1
          zi_all(ibnd, ik) = aux(i)
        ENDDO
      ENDDO
      CLOSE(iufilsigma_all)
    ENDIF
  ENDIF
  ! 
  CALL mp_bcast (exst, ionode_id, inter_pool_comm)
  CALL mp_bcast (exst, root_pool, intra_pool_comm)  
  !
  IF (exst) THEN
    CALL mp_bcast (iq, ionode_id, inter_pool_comm)
    CALL mp_bcast (iq, root_pool, intra_pool_comm)
    CALL mp_bcast (sigmar_all, ionode_id, inter_pool_comm)
    CALL mp_bcast (sigmar_all, root_pool, intra_pool_comm)
    CALL mp_bcast (sigmai_all, ionode_id, inter_pool_comm)
    CALL mp_bcast (sigmai_all, root_pool, intra_pool_comm)
    CALL mp_bcast (zi_all, ionode_id, inter_pool_comm)
    CALL mp_bcast (zi_all, root_pool, intra_pool_comm)
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1 ) THEN
      sigmar_all(:,1:lower_bnd-1) = zero
      sigmai_all(:,1:lower_bnd-1) = zero
      zi_all(:,1:lower_bnd-1) = zero
    ENDIF
    IF (upper_bnd < nktotf ) THEN
      sigmar_all(:,upper_bnd+1:nktotf) = zero
      sigmai_all(:,upper_bnd+1:nktotf) = zero
      zi_all(:,upper_bnd+1:nktotf) = zero
    ENDIF
    ! 
    WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ',(iq-1),'/',nqtotf
  ENDIF
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE electron_read
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE tau_write(iq,nqtotf,nktotf,second)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nstemp
  USE io_global, ONLY : stdout, meta_ionode_id
  USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb
  USE io_epw,    ONLY : iufiltau_all
  USE io_files,  ONLY : diropn
  USE mp,        ONLY : mp_barrier
  USE mp_world,  ONLY : mpime
  USE constants_epw, ONLY : zero
  USE transportcom, ONLY : lower_bnd, upper_bnd
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(INOUT) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  LOGICAL, INTENT(IN) :: second
  !! IF we have two Fermi level
  ! 
  ! Local variable
  LOGICAL :: exst
  !! 
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
  INTEGER :: idir
  !! Direction index
  ! 
  REAL(KIND=DP) :: aux ( 2 * nstemp * (ibndmax-ibndmin+1) * nktotf +2 )
  !! Vector to store the array inv_tau_all and zi_all
  !
  IF (mpime .eq. meta_ionode_id) THEN
    !
    ltau_all = 2 * nstemp * (ibndmax-ibndmin+1) * nktotf +2
    ! First element is the iteration number
    aux(1) = iq -1   ! -1 because we will start at the next one. 
    aux(2) = nqtotf
    i = 2
    ! 
    DO itemp=1, nstemp
      DO ik=1, nktotf
        DO ibnd=1, (ibndmax-ibndmin+1)
          i = i +1
          aux(i) = inv_tau_all(itemp,ibnd, ik)
        ENDDO
      ENDDO
    ENDDO
    !
    DO itemp=1, nstemp
      DO ik=1, nktotf
        DO ibnd=1, (ibndmax-ibndmin+1)
          i = i +1
          aux(i) = zi_allvb(itemp,ibnd, ik) 
        ENDDO
      ENDDO
    ENDDO
    CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
    CALL davcio ( aux, ltau_all, iufiltau_all, 1, +1 )
    CLOSE(iufiltau_all)
    ! 
    IF (second) THEN
      ! First element is the iteration number
      aux(1) = iq -1   ! -1 because we will start at the next one. 
      aux(2) = nqtotf
      i = 2
      ! 
      DO itemp=1, nstemp
        DO ik=1, nktotf
          DO ibnd=1, (ibndmax-ibndmin+1)
            i = i +1
            aux(i) = inv_tau_allcb(itemp,ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      DO itemp=1, nstemp
        DO ik=1, nktotf
          DO ibnd=1, (ibndmax-ibndmin+1)
            i = i +1
            aux(i) = zi_allcb(itemp,ibnd, ik)     
          ENDDO
        ENDDO
      ENDDO
      ! 
      CALL diropn (iufiltau_all, 'tau_restart_CB', ltau_all, exst)
      CALL davcio ( aux, ltau_all, iufiltau_all, 1, +1 )
      CLOSE(iufiltau_all)   
    ENDIF
    ! 
  ENDIF
  ! 
  ! Make everythin 0 except the range of k-points we are working on
  IF (lower_bnd > 1 ) inv_tau_all(:,:,1:lower_bnd-1) = zero
  IF (upper_bnd < nktotf ) inv_tau_all(:,:,upper_bnd+1:nktotf) = zero
  IF (second) THEN
    IF (lower_bnd > 1 ) inv_tau_allcb(:,:,1:lower_bnd-1) = zero
    IF (upper_bnd < nktotf ) inv_tau_allcb(:,:,upper_bnd+1:nktotf) = zero
  ENDIF
  ! Same for the Znk factor
  IF (lower_bnd > 1 ) zi_allvb(:,:,1:lower_bnd-1) = zero
  IF (upper_bnd < nktotf ) zi_allvb(:,:,upper_bnd+1:nktotf) = zero
  IF (second) THEN
    IF (lower_bnd > 1 ) zi_allcb(:,:,1:lower_bnd-1) = zero
    IF (upper_bnd < nktotf ) zi_allcb(:,:,upper_bnd+1:nktotf) = zero
  ENDIF
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE tau_write
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  SUBROUTINE tau_read(iq,nqtotf,nktotf,second)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, meta_ionode_id
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb
  USE io_epw,    ONLY : iufiltau_all
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE epwcom,    ONLY : nstemp
  USE constants_epw, ONLY : kelvin2eV, &
                            meV2invps, eps4, zero
  USE transportcom, ONLY : lower_bnd, upper_bnd
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : world_comm
  USE mp_world,  ONLY : mpime
  !
  IMPLICIT NONE
  !
  ! Local variable
  LOGICAL :: exst
  !
  INTEGER, INTENT(INOUT) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  LOGICAL, INTENT(IN) :: second
  !! IF we have two Fermi level
  ! 
  ! Local variables
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
  REAL(KIND=DP) :: aux ( 2 * nstemp * (ibndmax-ibndmin+1) * nktotf + 2 )
  !! Vector to store the array
  ! 
  CHARACTER (len=256) :: name1
  !
  IF (mpime .eq. meta_ionode_id) THEN
    !
    ! First inquire if the file exists
#if defined(__MPI)
    name1 = trim(tmp_dir) // trim(prefix) // '.tau_restart1'
#else
    name1 = trim(tmp_dir) // trim(prefix) // '.tau_restart'
#endif 
    INQUIRE(file = name1, exist=exst)
    ! 
    IF (exst) THEN ! read the file
      !
      ltau_all = 2 * nstemp * (ibndmax-ibndmin+1) * nktotf + 2
      CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
      CALL davcio ( aux, ltau_all, iufiltau_all, 1, -1 )
      !
      ! First element is the iteration number
      iq = aux(1)
      iq = iq + 1 ! we need to start at the next q
      nqtotf_read = aux(2)
      !print*, 'iq',iq
      !print*, 'nqtotf_read ',nqtotf_read
      IF ( nqtotf_read /= nqtotf) CALL errore('io_scattering',&
        &'Error: The current total number of q-point is not the same as the read one. ',1)
      ! 
      i = 2
      DO itemp=1, nstemp
        DO ik=1, nktotf
          DO ibnd=1, (ibndmax-ibndmin+1)
            i = i +1
            inv_tau_all(itemp,ibnd, ik) = aux(i)
          ENDDO
        ENDDO
      ENDDO
      ! 
      DO itemp=1, nstemp
        DO ik=1, nktotf
          DO ibnd=1, (ibndmax-ibndmin+1)
            i = i +1
            zi_allvb(itemp,ibnd, ik) = aux(i)
          ENDDO
        ENDDO
      ENDDO 
      CLOSE(iufiltau_all)
    ENDIF
    ! 
    IF (second) THEN
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = trim(tmp_dir) // trim(prefix) // '.tau_restart_CB1'
#else
      name1 = trim(tmp_dir) // trim(prefix) // '.tau_restart_CB'
#endif 
      INQUIRE(file = name1, exist=exst)
      ! 
      IF (exst) THEN ! read the file
        !
        ltau_all = nstemp * (ibndmax-ibndmin+1) * nktotf +2
        CALL diropn (iufiltau_all, 'tau_restart_CB', ltau_all, exst)
        CALL davcio ( aux, ltau_all, iufiltau_all, 1, -1 )
        !
        ! First element is the iteration number
        iq = aux(1)
        iq = iq + 1 ! we need to start at the next q
        nqtotf_read = aux(2)
        IF ( nqtotf_read /= nqtotf) CALL errore('io_scattering',&
          &'Error: The current total number of q-point is not the same as the read one. ',1)
        ! 
        i = 2
        DO itemp=1, nstemp
          DO ik=1, nktotf
            DO ibnd=1, (ibndmax-ibndmin+1)
              i = i +1
              inv_tau_allcb(itemp,ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        ! 
        DO itemp=1, nstemp
          DO ik=1, nktotf
            DO ibnd=1, (ibndmax-ibndmin+1)
              i = i +1
              zi_allcb(itemp,ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufiltau_all)
        WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau_CB: ',iq,'/',nqtotf 
      ENDIF
      ! 
    ENDIF ! second
    !
  ENDIF
  ! 
  CALL mp_bcast (exst, meta_ionode_id, world_comm)
  !
  IF (exst) THEN
    CALL mp_bcast (iq,          meta_ionode_id, world_comm)
    CALL mp_bcast (inv_tau_all, meta_ionode_id, world_comm)
    CALL mp_bcast (zi_allvb,    meta_ionode_id, world_comm)
    IF (second) CALL mp_bcast (inv_tau_allcb, meta_ionode_id, world_comm)
    IF (second) CALL mp_bcast (zi_allcb, meta_ionode_id, world_comm)
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1 )      inv_tau_all(:,:,1:lower_bnd-1) = zero
    IF (upper_bnd < nktotf ) inv_tau_all(:,:,upper_bnd+1:nktotf) = zero
    IF (lower_bnd > 1 )      zi_allvb(:,:,1:lower_bnd-1) = zero
    IF (upper_bnd < nktotf ) zi_allvb(:,:,upper_bnd+1:nktotf) = zero
    !  
    IF (second) THEN
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1 )      inv_tau_allcb(:,:,1:lower_bnd-1) = zero
      IF (upper_bnd < nktotf ) inv_tau_allcb(:,:,upper_bnd+1:nktotf) = zero
      IF (lower_bnd > 1 )      zi_allcb(:,:,1:lower_bnd-1) = zero
      IF (upper_bnd < nktotf ) zi_allcb(:,:,upper_bnd+1:nktotf) = zero
    ENDIF 
    ! 
    WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau: ',iq,'/',nqtotf
  ENDIF
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE tau_read
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  SUBROUTINE merge_read(nktotf, nqtotf_new, inv_tau_all_new)
  !----------------------------------------------------------------------------
  !
#if defined(__SX6)
#  define DIRECT_IO_FACTOR 1
#else
#  define DIRECT_IO_FACTOR 8 
#endif
  ! 
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin
  USE io_epw,    ONLY : iufiltau_all
  USE io_files,  ONLY : tmp_dir, diropn
  USE epwcom,    ONLY : nstemp, restart_filq
  USE constants_epw, ONLY : kelvin2eV, &
                            meV2invps, eps4
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm, root_pool
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  ! Local variable
  LOGICAL :: exst
  !
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  INTEGER, INTENT(OUT) :: nqtotf_new
  !! Total number of q-points
  REAL(KIND=DP), INTENT(INOUT) :: inv_tau_all_new(nstemp, ibndmax-ibndmin+1, nktotf)
  !! Scattering rate read from file restart_filq
  ! 
  ! Local variables
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
  INTEGER(kind=8) :: unf_recl
  !! 
  REAL(KIND=DP) :: aux ( nstemp * (ibndmax-ibndmin+1) * nktotf + 2 )
  !! Vector to store the array 
  CHARACTER (len=256) :: name1 
  ! 
  !
  IF (mpime.eq.ionode_id) THEN
    !
    ! First inquire if the file exists
    name1 = trim(tmp_dir) // trim(restart_filq)
    INQUIRE(file = name1, exist=exst)
    ! 
    IF (exst) THEN ! read the file
      !
      ltau_all = nstemp * (ibndmax-ibndmin+1) * nktotf +2
      !CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
      ! 
      unf_recl = DIRECT_IO_FACTOR * int(ltau_all, kind=kind(unf_recl))
      open (unit = iufiltau_all, file = restart_filq, iostat = ios, form ='unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
      !  
      CALL davcio ( aux, ltau_all, iufiltau_all, 1, -1 )
      !
      ! First element is the iteration number
      iq = aux(1)
      iq = iq + 1 ! we need to start at the next q
      nqtotf_new = aux(2)
      ! 
      i = 2
      DO itemp=1, nstemp
        DO ik=1, nktotf
          DO ibnd=1, (ibndmax-ibndmin+1)
            i = i +1
            inv_tau_all_new(itemp,ibnd, ik) = aux(i)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iufiltau_all)
    ENDIF
  ENDIF
  ! 
  CALL mp_bcast (exst, ionode_id, inter_pool_comm)
  CALL mp_bcast (exst, root_pool, intra_pool_comm)
  !
  IF (exst) THEN
    CALL mp_bcast (nqtotf_new, ionode_id, inter_pool_comm)
    CALL mp_bcast (nqtotf_new, root_pool, intra_pool_comm)
    CALL mp_bcast (inv_tau_all_new, ionode_id, inter_pool_comm)
    CALL mp_bcast (inv_tau_all_new, root_pool, intra_pool_comm)
    ! 
    WRITE(stdout, '(a,a)' ) '     Correctly read file ',restart_filq
  ENDIF
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE merge_read
  !----------------------------------------------------------------------------
