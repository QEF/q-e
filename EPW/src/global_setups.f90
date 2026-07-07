  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2024 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE global_setups
  !----------------------------------------------------------------------
  !!
  !! Setup routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !--------------------------------------------------------------------------
    SUBROUTINE use_wannier_setup()
    !--------------------------------------------------------------------------
    !!
    !! Setup routine for Wannier interpolation
    !!
    !! 10/2025 Zhe Liu initial implementation
    !!
    !--------------------------------------------------------------------------
    USE kinds,              ONLY : DP
    USE input,              ONLY : lsda
    USE noncollin_module,   ONLY : noncolin
    USE global_var,         ONLY : spin_fac
    !
    IMPLICIT NONE
    !
    IF (noncolin .OR. TRIM(lsda) /= 'none') THEN
      spin_fac = 1.0d0
    ELSE
      spin_fac = 2.0d0
    ENDIF
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE use_wannier_setup
    !--------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------
    SUBROUTINE lsda_setup(nbndsub_aux, nbndskip_aux, nelec_aux, etf_aux)
    !---------------------------------------------------------------------------------
    !! This subroutine sets the number of skipped bands, the number of wannierized 
    !! bands and interpolates the eigenvalues for the second spin channel in an LSDA 
    !! calculation.
    !!
    !! The main reason for this is to be able to determine the Fermi level of the system 
    !! when one uses a collinear magnetic calculation. The total numeber of electrons cannot
    !! be calculated without the hamiltonian of both spin channels and then performing the 
    !! integral over the density of states. In that regard this forces that the workflow 
    !! of an LSDA calculation requires a previous unfolding of the Hamiltonian
    !! into the Wannier basis for both spin channels. Once this is done, one is able to 
    !! interpolate any quantity for a given spin channel. 
    !!
    !! There has been several extra changes in the main code to allow for the LSDA 
    !! calculations in EPW. For a more detailed information about the changes please check 
    !! the epw developers manual
    !
    USE kinds,            ONLY : DP
    USE pwcom,            ONLY : nelec
    USE cell_base,        ONLY : at, bg
    USE ions_base,        ONLY : nat
    USE input,            ONLY : mp_mesh_k, lsda
    USE ep_constants,     ONLY : zero, czero, twopi, ci
    USE io_global,        ONLY : ionode_id, stdout
    USE io_var,           ONLY : crystal, epwdata, iuwigner
    USE global_var,       ONLY : xkf, wkf, nkqf, bztoibz, xkf_irr, wkf_irr, s_bztoibz, &
                                 spin_fac
    USE wannier2bloch,    ONLY : hamwan2bloch
    USE wigner,           ONLY : wigner_divide_ndegen
    USE bzgrid,           ONLY : loadkmesh_para
    USE mp,               ONLY : mp_bcast
    USE mp_global,        ONLY : my_pool_id
    USE mp_world,         ONLY : mpime, world_comm
    USE symmetry,         ONLY : kpoints_time_reversal_init
    
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: nbndsub_aux
    !! Bands of the second spin channel
    INTEGER, INTENT(inout) :: nbndskip_aux
    !! Skipped bands on the second spin channel
    REAL(KIND = DP), INTENT(inout) :: nelec_aux
    !! Number of electrons in the second spin channel (LSDA case)
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: etf_aux(:, :)
    !! Interpolated eigenvalues for the second spin channel
    !
    ! Local variables
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: spin_channel
    !! Current spin channel in the run 
    CHARACTER(LEN = 256) :: lsda_opposite
    !! Opposite spin channel in the run 
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! Counter on coarse k-point grid
    INTEGER :: ikk
    !! Counter on k-point when you have paired k and q
    INTEGER :: ikq
    !! Paired counter so that q is adjacent to its k
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: jbnd
    !! Counter on band
    INTEGER :: irk
    !! counter on real space vectors (LSDA) case
    INTEGER :: ierr
    !! Error status
    INTEGER :: nrr_k_aux
    !! Number of skipped bands in the second spin channel (LSDA case)
    INTEGER :: dims_aux
    !! Dimensions of the degeneracy array for the second spin channel (LSDA case)
    INTEGER :: dims_aux2
    !! Dimensions of the degeneracy array for the second spin channel (LSDA case)
    INTEGER, ALLOCATABLE :: irvec_k_aux(:, :)
    !! Wigner-Seitz vectors for the second spin channel (LSDA case)
    INTEGER, ALLOCATABLE :: ndegen_k_aux(:, :, :)
    !! K-point degeneracy for the second spin channel (LSDA case)
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: xxk(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: dummy(3)
    !! Dummy variable
    REAL(KIND = DP), ALLOCATABLE :: dummy2(:, :)
    !! dummy array dim 2
    REAL(KIND = DP), ALLOCATABLE :: dummy3(:, :, :)
    !! dummy array dim 3
    REAL(KIND = DP), ALLOCATABLE :: irvec_r(:, :)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), ALLOCATABLE :: rdotk(:)
    !! $r\cdot k$
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: chw_aux(:, :, :)
    !! Hammiltonian in the wannier basis for the second spin channel 
    !
    spin_channel = lsda
    ! We store the original value of lsda and then change it to the other spin channel
    ! so that we can read and perform the interpolation of the eigen energies on the fine mesh 
    ! this is important to determine the Fermi level in the LSDA case on the fine mesh since 
    ! both spin channels are needed.
    SELECT CASE(spin_channel)
    CASE('up')
      lsda_opposite = 'down'
    CASE('down')
      lsda_opposite = 'up'
    END SELECT
    ! We need to read the fmt files from the other spin channel in order to compute the Fermi level
    IF (mpime == ionode_id) THEN
      ! Obtain the number of electrons and the number of skipped bands in the second spin channel
      filint = 'crystal.fmt'
      IF (TRIM(lsda_opposite) == 'down') filint = 'crystal.down.fmt'
      INQUIRE(FILE = filint, EXIST = exst)
      IF (.NOT. exst) CALL errore('lsda_setup','Error '// TRIM(filint)//' not found', 1)
      OPEN(UNIT = crystal, FILE = TRIM(filint), STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('lsda_setup', 'error opening ' // TRIM(filint), 1)
      READ(crystal,*) ik
      READ(crystal,*) ik
      READ(crystal,*) nelec_aux, nbndskip_aux
      CLOSE(crystal)
      ! Sanity check
      !! Change to crystal and epwdata 
      IF (nelec /= nelec_aux) CALL errore('lsda_setup','Error total number of electrons read in ' & 
                                          //TRIM(filint)//' file is inconsistent with nelec', 1)
      ! Obtain the number of bands and the nrr_k 
      filint = 'epwdata.fmt'
      IF (TRIM(lsda_opposite)=='down') filint = 'epwdata.down.fmt'
      INQUIRE(FILE = filint, EXIST = exst)
      IF (.NOT. exst) CALL errore('lsda_setup','Error '// TRIM(filint)//' not found', 1)
      ! Allocate dummy variables to read epw file
      ALLOCATE(dummy3(3, 3, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup','Error allocating dummy3', 1)
      ALLOCATE(dummy2(3, 3), STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating dummy2', 1)
      ! Open epwdata file
      OPEN(UNIT = epwdata, FILE = TRIM(filint), STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('lsda_setup', 'error opening ' // TRIM(filint), 1)
      READ(epwdata,*) dummy(1)
      READ(epwdata,*) nbndsub_aux, nrr_k_aux, ik, ikk, ikq
      READ(epwdata,*) dummy3, dummy2
      ! Deallocate unnecesary variables
      DEALLOCATE(dummy2, STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating dummy2', 1)
      DEALLOCATE(dummy3, STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating dummy3', 1)
      ! Keep epwdata file open
    ENDIF
    ! Broadcast to all nodes 
    CALL mp_bcast(nbndsub_aux, ionode_id, world_comm)
    CALL mp_bcast(nbndskip_aux, ionode_id, world_comm)
    CALL mp_bcast(nrr_k_aux, ionode_id, world_comm)
    CALL mp_bcast(nelec_aux, ionode_id, world_comm)
    ! Allocate the Hammiltonian for the second spin channel
    ALLOCATE(chw_aux(nbndsub_aux, nbndsub_aux, nrr_k_aux), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating cfac_aux', 1)
    chw_aux(:, :, :) = czero
    ! Read the Hammiltonian
    IF (mpime == ionode_id) THEN
      DO ibnd = 1, nbndsub_aux
        DO jbnd = 1, nbndsub_aux
          DO irk = 1, nrr_k_aux
            READ (epwdata,*) chw_aux(ibnd, jbnd, irk)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(epwdata)
      ! Cose epwdata file
    ENDIF
    ! Broadcast the Hammiltonian
    CALL mp_bcast(chw_aux, ionode_id, world_comm)   
    ! Allocate wigner seitz vectors
    ALLOCATE(irvec_k_aux(3, nrr_k_aux), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating irvec_k_aux', 1)
    ! Read wigner seitz data file
    IF (mpime == ionode_id) THEN
      ! Obtain the wigner seitz vectors 
      filint = 'wigner.fmt'
      IF (TRIM(lsda_opposite) == 'down') filint = 'wigner.down.fmt'
      INQUIRE(FILE = filint, EXIST = exst)
      IF (.NOT. exst) CALL errore('lsda_setup','Error '// TRIM(filint) // ' not found', 1)  
      OPEN(NEWUNIT=iuwigner, FILE=TRIM(filint), ACTION='read', STATUS='old') 
      IF (ios /= 0) CALL errore ('lsda_setup', 'error opening ' // TRIM(filint), 1)
      READ(iuwigner, *) nrr_k_aux, ikk, ikq, dims_aux, dims_aux2
      ! Allocate dummy integer
    ENDIF
    CALL mp_bcast(dims_aux, ionode_id, world_comm)
    CALL mp_bcast(dims_aux2, ionode_id, world_comm)
    !! Allocate ndegen_k_aux
    ALLOCATE(ndegen_k_aux(nrr_k_aux, dims_aux, dims_aux), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup','Error allocating ndegen_k_aux', 1)
    !
    IF (mpime == ionode_id) THEN
      !
      DO irk = 1, nrr_k_aux
        READ (iuwigner, *) irvec_k_aux(:, irk), dummy(1)
        DO ibnd = 1, dims_aux
          READ(iuwigner, *) ndegen_k_aux(irk, ibnd, :)
        ENDDO
      ENDDO
      ! Close and Deallocate unnecesary variables
      CLOSE(iuwigner)
    ENDIF
    ! Broadcast wigner seitz vectors and degeneracies
    CALL mp_bcast(irvec_k_aux, ionode_id, world_comm)
    CALL mp_bcast(ndegen_k_aux, ionode_id, world_comm)
    ! Divide the hamiltonian by the k degeneracy
    CALL wigner_divide_ndegen(chw_aux, 1, nbndsub_aux, nrr_k_aux, 1, ndegen_k_aux, dims_aux)
    !
    DEALLOCATE(ndegen_k_aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating ndegen_k_aux', 1)
    ! Load the fine k mesh for this second spin channel
    CALL kpoints_time_reversal_init()
    !
    CALL loadkmesh_para()
    !
    ! Allocate the phase factors, ukk, wigner seitz vector and eiegen energies for the second spin channel
    ALLOCATE(cfac(nrr_k_aux), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating cfac', 1)
    ALLOCATE(rdotk(nrr_k_aux), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating rdotk', 1)
    ALLOCATE(cufkk(nbndsub_aux, nbndsub_aux))
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating cufkk', 1)
    ALLOCATE(etf_aux(nbndsub_aux, nkqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating etf_aux', 1)
    ALLOCATE(irvec_r(3, nrr_k_aux), STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error allocating irvec_r', 1)
    ! Zeroing everything - initialization is important !
    cfac(:)        = czero
    cufkk(:, :)    = czero
    rdotk(:)       = zero
    etf_aux(:, :)  = zero
    irvec_r = REAL(irvec_k_aux, KIND = DP)
    xxq = 0.d0
    DO ik = 1, nkqf
      !
      xxk = xkf(:, ik)
      !
      IF (2 * (ik / 2) == ik) THEN
        !
        !  this is a k+q point : redefine as xkf (:, ik-1) + xxq
        !
        CALL cryst_to_cart(1, xxq, at, -1)
        xxk = xkf(:, ik - 1) + xxq
        CALL cryst_to_cart(1, xxq, bg, 1)
        !
      ENDIF
      !
      ! Compute the phase factor for the second spin channel
      CALL DGEMV('t', 3, nrr_k_aux, twopi, irvec_r, 3, xxk, 1, 0.0_DP, rdotk, 1)
      cfac(:) = EXP(ci * rdotk(:))
      ! Compute the eigen energies for the second spin channel
      CALL hamwan2bloch(nbndsub_aux, nrr_k_aux, cufkk, etf_aux(:, ik), chw_aux, cfac)
    ENDDO
    ! Dealocate unnecesary variables
    IF (my_pool_id == ionode_id .AND. mp_mesh_k) THEN
      DEALLOCATE(xkf_irr, STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating xkf_irr', 1)
      DEALLOCATE(wkf_irr, STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating wkf_irr', 1)
    ENDIF
    DEALLOCATE(cufkk, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating cufkk', 1)
    DEALLOCATE(chw_aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating chw_aux', 1)
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating cfac', 1)
    DEALLOCATE(rdotk, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating rdotk', 1)
    DEALLOCATE(irvec_k_aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating irvec_k_aux', 1)
    DEALLOCATE(irvec_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating irvec_r', 1)
    DEALLOCATE(wkf, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating wkf', 1)
    DEALLOCATE(xkf, STAT = ierr)
    IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating xkf', 1)
    IF (mp_mesh_k) THEN
      DEALLOCATE(bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating bztoibz',1)
      DEALLOCATE(s_bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('lsda_setup', 'Error deallocating s_bztoibz', 1)
    ENDIF
    !
    END SUBROUTINE lsda_setup
    !-----------------------------------------------------------------------------
  !----------------------------------------------------------------------
  END MODULE global_setups
  !----------------------------------------------------------------------
