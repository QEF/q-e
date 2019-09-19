  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE plot_band()
  !-----------------------------------------------------------------------
  !!
  !! This SUBROUTINE writes output files for phonon dispersion and band structure 
  !! RM : this SUBROUTINE should be tested
  !! SP : Modified so that it works with the current plotband.x of QE 5
  !!
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, filqf, filkf
  USE elph2,         ONLY : etf, nkf, nqtotf, wf, xkf, xqf, nkqtotf, nktotf
  USE constants_epw, ONLY : ryd2mev, ryd2ev
  USE io_var,        ONLY : iufilfreq, iufileig
  USE elph2,         ONLY : nkqf
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id
  USE poolgathering, ONLY : poolgather2
  !
  IMPLICIT NONE
  !
  ! Local variables
  INTEGER :: ik
  !! Global k-point index
  INTEGER :: ikk
  !! Index for the k-point
  INTEGER :: ikq
  !! Index for the q-point
  INTEGER :: ibnd
  !! Band index
  INTEGER :: imode
  !! Mode index
  INTEGER :: iq
  !! Global q-point index
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: dist
  !! Distance from Gamma
  REAL(KIND = DP) :: dprev
  !! Previous distance
  REAL(KIND = DP) :: dcurr
  !! Current distance
  REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
  !! K-points on the full k grid (all pools)
  REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
  !! Eigenenergies on the full k grid (all pools)
  !
  IF (filqf /= ' ') THEN
    ! 
    IF (my_pool_id == ionode_id) THEN
      !
      OPEN(iufilfreq, FILE = "phband.freq", FORM = 'formatted')
      WRITE(iufilfreq, '(" &plot nbnd=",i4,", nks=",i6," /")') nmodes, nqtotf
      !
      ! crystal to cartesian coordinates
      CALL cryst_to_cart(nqtotf, xqf, bg, 1)
      !
      dist = 0.d0
      dprev = 0.d0
      dcurr = 0.d0
      DO iq = 1, nqtotf
        !
        IF (iq /= 1) THEN  
          dist = SQRT((xqf(1, iq) - xqf(1, iq - 1)) * (xqf(1, iq) - xqf(1, iq - 1)) & 
                    + (xqf(2, iq) - xqf(2, iq - 1)) * (xqf(2, iq) - xqf(2, iq - 1)) & 
                    + (xqf(3, iq) - xqf(3, iq - 1)) * (xqf(3, iq) - xqf(3, iq - 1)))
        ELSE 
          dist = 0.d0
        ENDIF
        dcurr = dprev + dist
        dprev = dcurr
        WRITE(iufilfreq, '(10x,3f10.6)') xqf(:, iq)
        WRITE(iufilfreq, '(1000f14.4)') (wf(imode, iq) * ryd2mev, imode = 1, nmodes)
        !
      ENDDO
      CLOSE(iufilfreq)
      !
      ! back from cartesian to crystal coordinates
      CALL cryst_to_cart(nqtotf, xqf, at, -1)
      !
    ENDIF
  ENDIF ! filqf
  ! 
  IF (filkf /= ' ') THEN
    !
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
    ENDDO
    !
    ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_band', 'Error allocating xkf_all', 1)
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_band', 'Error allocating etf_all', 1)
    !
#if defined(__MPI)
    CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    CALL mp_barrier(inter_pool_comm)
#else    
    !
    xkf_all = xkf
    etf_all = etf
#endif
    !
    IF (my_pool_id == ionode_id) THEN
      !
      OPEN(iufileig, FILE = "band.eig", FORM = 'formatted')
      WRITE(iufileig, '(" &plot nbnd=",i4,", nks=",i6," /")') nbndsub, nktotf
      !
      ! crystal to cartesian coordinates
      CALL cryst_to_cart(nkqtotf, xkf_all, bg, 1)
      !
      dist = 0.d0
      dprev = 0.d0
      dcurr = 0.d0
      DO ik = 1, nktotf
         !
         ikk = 2 * ik - 1
         ikq = ikk + 1
         !
         IF (ikk /= 1) THEN
            dist = SQRT((xkf_all(1, ikk) - xkf_all(1, ikk - 2)) * (xkf_all(1, ikk) - xkf_all(1, ikk - 2)) &
                      + (xkf_all(2, ikk) - xkf_all(2, ikk - 2)) * (xkf_all(2, ikk) - xkf_all(2, ikk - 2)) &
                      + (xkf_all(3, ikk) - xkf_all(3, ikk - 2)) * (xkf_all(3, ikk) - xkf_all(3, ikk - 2)))
         ELSE
            dist = 0.d0
         ENDIF
         dcurr = dprev + dist
         dprev = dcurr
         WRITE(iufileig, '(10x,3f10.6)') xkf_all(:, ikk)
         WRITE(iufileig, '(1000f20.12)') (etf_all(ibnd, ikk) * ryd2ev, ibnd = 1, nbndsub)
         !
      ENDDO
      CLOSE(iufileig)
      !
      ! back from cartesian to crystal coordinates
      CALL cryst_to_cart(nkqtotf, xkf_all, at, -1)
      !
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_band', 'Error deallocating xkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_band', 'Error deallocating etf_all', 1)
    !
  ENDIF ! filkf
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE plot_band
  !-----------------------------------------------------------------------
