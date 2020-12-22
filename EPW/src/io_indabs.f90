  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !----------------------------------------------------------------------
  MODULE io_indabs
  !----------------------------------------------------------------------
  !!
  !! This module contains the writing and reading routines relates to indabs.
  !! This is bascially for restart
  !! 12/20/2020 Implemented by Xiao Zhang
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !----------------------------------------------------------------------
    SUBROUTINE indabs_write(iq, totq, nktotf, epsilon2_abs, epsilon2_abs_lorenz)
    !----------------------------------------------------------------------
    !!
    !! Write indirect optical spectra
    !!
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : epsilon2_abs, epsilon2_abs_lorenz
    USE io_var,    ONLY : iuindabs
    USE io_files,  ONLY : diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    USE epwcom,    ONLY : nstemp, omegamin, omegamax, omegastep, &
                          neta, nomega
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
!    INTEGER, INTENT(in) :: nomega
    !! Number of energt points
    INTEGER :: leps_all
    !! Length of the vector
!    INTEGER :: neta = 9
    !! Number of broadenings
    REAL(KIND = DP), INTENT(inout) :: epsilon2_abs(3, nomega, neta, nstemp)
    !! Imaginary part of the dielectric function
    REAL(KIND = DP), INTENT(inout) :: epsilon2_abs_lorenz(3, nomega, neta, nstemp)
    !! Imaginary part of the dielectric function, Lorenzian broadening
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: iw
    !! energy index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ieta
    !! eta index
    INTEGER :: ipol
    !! Polarization index
    REAL(KIND = DP) :: aux(2 * 3 * nomega * neta * nstemp + 2)
    !! vector to store the array
    !
    IF (mpime = ionode_id) THEN
      !
      leps_all = 2 * 3 * nomega * neta * nstemp + 2
      aux(1) = REAL(iq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO itemp = 1, nstemp
        DO ipol = 1, 3
          DO ieta = 1, neta
            DO iw = 1, nomega
              i = i + 1
                aux(i) = epsilon2_abs(ipol, iw, ieta, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ipol = 1, 3
          DO ieta = 1, neta
            DO iw = 1, nomega
              i = i + 1
                aux(i) = epsilon2_abs_lorenz(ipol, iw, ieta, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iuindabs, 'indabs_restart', leps_all, exst)
      CALL davcio(aux, leps_all, iuindabs, 1, +1)
      CLOSE(iuindabs)
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE indabs_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE indabs_read(iq, totq, nktotf, epsilon2_abs, epsilon2_abs_lorenz)
    !----------------------------------------------------------------------
    !!
    !! Read indirect optical spectra
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode_id
    USE elph2,     ONLY : epsilon2_abs, epsilon2_abs_lorenz
    USE io_var,    ONLY : iuindabs
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE epwcom,    ONLY : nstemp, omegamin, omegamax, omegastep, &
                          neta, nomega
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
!    INTEGER, INTENT(in) :: nomega
    !! Number of energt points
    INTEGER :: leps_all
    !! Length of the vector
!    INTEGER :: neta = 9
    !! Number of broadenings
    REAL(KIND = DP), INTENT(inout) :: epsilon2_abs(3, nomega, neta, nstemp)
    !! Imaginary part of the dielectric function
    REAL(KIND = DP), INTENT(inout) :: epsilon2_abs_lorenz(3, nomega, neta, nstemp)
    !! Imaginary part of the dielectric function, Lorenzian broadening
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: iw
    !! energy index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ieta
    !! eta index
    INTEGER :: ipol
    !! Polarization index
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    REAL(KIND = DP) :: aux(2 * 3 * nomega * neta * nstemp + 2)
    !! vector to store the array
    !
    !
    CHARACTER(LEN = 256) :: name1
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.indabs_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.indabs_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        leps_all = 2 * 3 * nomega * neta * nstemp + 2
        CALL diropn(iuindabs, 'indabs_restart', leps_all, exst)
        CALL davcio(aux, leps_all, iuindabs_all, 1, -1)
        !
        !
        ! First element is the iteration number
        iq = INT(aux(1))
        iq = iq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('indabs_read', &
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ipol = 1, 3
            DO ieta = 1, neta
              DO iw = 1, nomega
                i = i + 1
                  epsilon2_abs(ipol, iw, ieta, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ipol = 1, 3
            DO ieta = 1, neta
              DO iw = 1, nomega
                i = i + 1
                  epsilon2_abs_lorenz(ipol, iw, ieta, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iuindabs)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(epsilon2_abs, ionode_id, world_comm)
      CALL mp_bcast(epsilon2_abs_lorenz, ionode_id, world_comm)
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !----------------------------------------------------------------------
    END SUBROUTINE indabs_read
    !----------------------------------------------------------------------
    !
  !----------------------------------------------------------------------
  END MODULE io_indabs
  !----------------------------------------------------------------------
