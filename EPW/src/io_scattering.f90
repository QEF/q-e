  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------------
  SUBROUTINE F_write(Fi_all, nbnd, nkq, iter)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilFi_all
  USE io_files,  ONLY : diropn
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nbnd
  !! Number of bands
  INTEGER, INTENT(IN) :: nkq
  !! Number of k-points
  INTEGER, INTENT(IN) :: iter
  !! Iteration number
  REAL(KIND=DP), INTENT(IN) :: Fi_all(3,nbnd, nkq)
  !! Total scattering rate on all k-points collected from all pools in parallel case
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
  ! 
  REAL(KIND=DP) :: aux ( 3*nbnd* nkq +1 )
  !! Vector to store the array
  !
  WRITE(stdout,'(/5x,"Writing Fi_all to file"/)')
  !
  IF (mpime.eq.ionode_id) THEN
    !
    lfi_all = 3*nbnd* nkq +1
    ! First element is the iteration number
    aux(1) = iter
    i = 1
    ! 
    DO ik=1, nkq
      DO ibnd=1, nbnd
        DO idir=1,3
          i = i +1
          aux(i) = Fi_all(idir,ibnd, ik)
        ENDDO
      ENDDO
    ENDDO
    CALL diropn (iufilFi_all, 'Fi_all', lfi_all, exst)
    CALL davcio ( aux, lfi_all, iufilFi_all, 1, +1 )
    CLOSE(iufilFi_all)
  ENDIF
  CALL mp_barrier(inter_pool_comm)
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE F_write
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE F_read(Fi_all, nbnd, nkq, iter)
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
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nbnd
  !! Number of bands
  INTEGER, INTENT(IN) :: nkq
  !! Number of k-points
  INTEGER, INTENT(OUT) :: iter
  !! Iteration number
  REAL(KIND=DP), INTENT(OUT) :: Fi_all(3,nbnd, nkq)
  !! Total scattering rate on all k-points collected from all pools in parallel case
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
  ! 
  CHARACTER (len=256) :: name1
 
  REAL(KIND=DP) :: aux ( 3*nbnd* nkq +1 )
  !! Vector to store the array
  !
  IF (mpime.eq.ionode_id) THEN
    ! 
    ! First inquire if the file exists
    name1 = trim(tmp_dir) // trim(prefix) // '.Fi_all1'
    INQUIRE(file = name1, exist=exst)
    ! 
    IF (exst) THEN ! read the file
      WRITE(stdout,'(/5x,"Restart iterative BTE: reading Fi_all from file"/)')
      !
      lfi_all = 3*nbnd* nkq + 1
      CALL diropn (iufilFi_all, 'Fi_all', lfi_all, exst)
      CALL davcio ( aux, lfi_all, iufilFi_all, 1, -1 )
      !
      ! First element is the iteration number
      iter = aux(1) 
      i = 1
      DO ik=1, nkq
        DO ibnd=1, nbnd
          DO idir=1,3
            i = i +1
            Fi_all(idir,ibnd, ik) = aux(i)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iufilFi_all)
    ENDIF
  ENDIF
  !
  CALL mp_bcast (Fi_all, ionode_id, inter_pool_comm)
  CALL mp_bcast (Fi_all, root_pool, intra_pool_comm)
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE F_read
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE scattering_write(itemp, etemp, ef0, etf_all, inv_tau_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf
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
  REAL(KIND=DP), INTENT(IN) :: inv_tau_all(nstemp,ibndmax-ibndmin+1, nkqtotf/2)
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
        WRITE(iufilscatt_rate,'(E24.14)') ryd2mev * meV2invps * inv_tau_all(itemp,ibnd,ik)
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
  SUBROUTINE scattering_read(etemp, ef0, etf_all, inv_tau_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf
  USE epwcom,    ONLY : nbndsub, nstemp
  USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                            meV2invps, eps4
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
  REAL(KIND=DP), INTENT(out) :: etf_all(nbndsub, nkqtotf/2)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  REAL(KIND=DP), INTENT(out) :: inv_tau_all(nstemp,ibndmax-ibndmin+1, nkqtotf/2)
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
  ! DBSP
  !print*,'ionode_id ',ionode_id
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
  USE elph2,     ONLY : ibndmax, ibndmin
  USE io_epw,    ONLY : iufilsigma_all
  USE io_files,  ONLY : diropn
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
  ! Local variable
  LOGICAL :: exst
  !
  INTEGER, INTENT(IN) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  REAL(KIND=DP), INTENT(IN) :: sigmar_all(ibndmax-ibndmin+1, nktotf)
  !! Real part of the electron-phonon self-energy accross all pools
  REAL(KIND=DP), INTENT(IN) :: sigmai_all(ibndmax-ibndmin+1, nktotf)
  !! Imaginary part of the electron-phonon self-energy accross all pools
  REAL(KIND=DP), INTENT(IN) :: zi_all(ibndmax-ibndmin+1, nktotf)
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
    aux(1) = iq
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

  CALL mp_barrier(inter_pool_comm)
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
  USE elph2,     ONLY : ibndmax, ibndmin
  USE io_epw,    ONLY : iufilsigma_all
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE epwcom,    ONLY : nbndsub, nstemp
  USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
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
    name1 = trim(tmp_dir) // trim(prefix) // '.sigma_restart1'
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
    WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ',iq,'/',nqtotf
  ENDIF
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE electron_read
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE tau_write(iq,nqtotf,nktotf,inv_tau_all)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nstemp
  USE io_global, ONLY : stdout
  USE elph2,     ONLY : ibndmax, ibndmin
  USE io_epw,    ONLY : iufiltau_all
  USE io_files,  ONLY : diropn
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(INOUT) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  ! 
  REAL(KIND=DP), INTENT(IN) :: inv_tau_all(nstemp,ibndmax-ibndmin+1,nktotf)
  !! Total scattering rate on all k-points collected from all pools in parallel case
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
  REAL(KIND=DP) :: aux ( nstemp * (ibndmax-ibndmin+1) * nktotf +1 )
  !! Vector to store the array
  !
  IF (mpime.eq.ionode_id) THEN
    !
    ltau_all = nstemp * (ibndmax-ibndmin+1) * nktotf +2
    ! First element is the iteration number
    aux(1) = iq
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
    CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
    CALL davcio ( aux, ltau_all, iufiltau_all, 1, +1 )
    CLOSE(iufiltau_all)
  ENDIF
  CALL mp_barrier(inter_pool_comm)
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE tau_write
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  SUBROUTINE tau_read(iq,nqtotf,nktotf,inv_tau_all)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilscatt_rate
  USE elph2,     ONLY : ibndmax, ibndmin
  USE io_epw,    ONLY : iufiltau_all
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE epwcom,    ONLY : nbndsub, nstemp
  USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
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
  INTEGER, INTENT(INOUT) :: iq
  !! Current q-point
  INTEGER, INTENT(IN) :: nqtotf
  !! Total number of q-points
  INTEGER, INTENT(IN) :: nktotf
  !! Total number of k-points
  REAL(KIND=DP), INTENT(INOUT) :: inv_tau_all(nstemp,ibndmax-ibndmin+1,nktotf)
  !! Scattering rate accross all pools
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
  REAL(KIND=DP) :: aux ( nstemp * (ibndmax-ibndmin+1) * nktotf + 2 )
  !! Vector to store the array
  ! 
  CHARACTER (len=256) :: name1
  !
  IF (mpime.eq.ionode_id) THEN
    !
    ! First inquire if the file exists
    name1 = trim(tmp_dir) // trim(prefix) // '.tau_restart1'
    INQUIRE(file = name1, exist=exst)
    ! 
    IF (exst) THEN ! read the file
      !
      ltau_all = nstemp * (ibndmax-ibndmin+1) * nktotf +2
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
      CLOSE(iufiltau_all)
    ENDIF
  ENDIF
  ! 
  CALL mp_bcast (exst, ionode_id, inter_pool_comm)
  CALL mp_bcast (exst, root_pool, intra_pool_comm)  
  !
  IF (exst) THEN
    CALL mp_bcast (iq, ionode_id, inter_pool_comm)
    CALL mp_bcast (iq, root_pool, intra_pool_comm)
    CALL mp_bcast (inv_tau_all, ionode_id, inter_pool_comm)
    CALL mp_bcast (inv_tau_all, root_pool, intra_pool_comm)
    ! 
    WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau: ',iq,'/',nqtotf
  ENDIF
  ! 
  !----------------------------------------------------------------------------
  END SUBROUTINE tau_read
  !----------------------------------------------------------------------------


