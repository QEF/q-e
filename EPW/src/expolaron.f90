!
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!-----------------------------------------------------------------------
MODULE expolaron
!-----------------------------------------------------------------------
  !!
  !! This module contains variables and subroutines of excitonic polaron
  !! Authored by Zhenbang Dai, with a lot of copy/paste from the polaron module
  !!
  USE kinds,  ONLY : dp
  !
  IMPLICIT NONE
  !! pool-dependent Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: Hamil(:, :), Hamil_work(:)
  !! Expolaron Hamiltonian and the work array for diagonalization
  COMPLEX(DP), ALLOCATABLE :: eigVec(:, :)
  !! Polaron eigenvectors
  COMPLEX(DP),  ALLOCATABLE :: Bmat(:,:), BtimesG(:,:), dtau(:,:)
  !! Bqv, Bqv*G, and polaron displacements
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE run_explrn(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
  !-----------------------------------------------------------------------
  !! Main driver to run excitonic polaron calculations
    USE input,         ONLY : explrn, plot_explrn_e, plot_explrn_h, nqc1, nqc2, nqc3, &
                              step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : wf, nktotf, nqtotf
    USE parallelism,   ONLY : fkbounds
    USE modes,         ONLY : nmodes
    USE io_global,     ONLY : stdout, ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT (in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT (in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT (in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT (in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT (in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT (in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT (in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), INTENT (in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    INTEGER :: lower_bnd, upper_bnd
    !! Lower and upper bound of the k parallelization
    INTEGER :: ikminusiexq
    !! Index of k+Q
    INTEGER :: imode
    !! Index of phonon mode
    INTEGER :: nktotf_temp
    !! Temporary variable to store the number of k points in the uniform grid. nktotf may become smaller
    REAL(KIND = DP), ALLOCATABLE :: wf_temp(:, :)
    !! Temporary array to store the phonon frequencies in the uniform grid
    !
    nktotf_temp = nktotf
    !
    ALLOCATE(wf_temp(nmodes, nqtotf))
    !
    wf_temp(:,:) = wf(:,:)
    !
    IF(explrn .AND. (.NOT. plot_explrn_e) .AND. (.NOT. plot_explrn_h)) THEN
      nktotf = nqc1 * nqc2 * nqc3 / (step_k1_explrn * step_k2_explrn * step_k3_explrn)!! ZD: This is because in this run, EPW no longer read pwscf results
      CALL read_exband()
      CALL read_phband()
      CALL read_ph_eig()
      !! Parallelize the construction of ex-plrn Hamiltonian by q=Q-Q'!!
      CALL fkbounds(nktotf, lower_bnd, upper_bnd) !! Maybe redundent, but just to make sure
      WRITE(stdout, '(5x, a, 2I10)') 'Check some params: nktotf, nmodes', nktotf, nmodes
      CALL read_ex_ph(lower_bnd, upper_bnd)
      CALL solve_for_explrn(nrr_q, ndegen_q, irvec_q, rws, nrws)
      ALLOCATE(wf(nmodes, nqtotf))
    !
    ELSEIF (explrn .AND. (plot_explrn_e .OR. plot_explrn_h)) THEN
      !
      nktotf = nqc1 * nqc2 * nqc3 / (step_k1_explrn * step_k2_explrn * step_k3_explrn)
      CALL plot_explrn_densities()
      !
    ENDIF
    !
    nktotf = nktotf_temp
    wf(:,:) = wf_temp
    !
    DEALLOCATE(wf_temp)
    !
  !-----------------------------------------------------------------------
  END SUBROUTINE run_explrn
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_exband()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, nkc1, nkc2, nkc3, step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : eigval_ex, nktotf
    USE ep_constants,  ONLY : ryd2ev
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE      
    !
    CHARACTER(LEN = 256)        :: filename
    !! File name
    INTEGER                     :: iexq
    !! Index of the exciton momentum
    INTEGER                     :: imode
    !! Index of the phonon mode
    INTEGER                     :: ierr
    !! Error status
    INTEGER                     :: ik1, ik2, ik3
    !! Index of Qx, Qy, and Qz
    INTEGER                     :: iq
    !! Index of the exciton momentum if some Q points are skipped
    REAL(KIND=DP), ALLOCATABLE  :: eigval_ex_temp(:,:)
    !! Temporary array for storing the BSE eigenvalues
    !
    filename = './G_full_epmatq/exband.fmt'
    !
    ALLOCATE(eigval_ex(negnv_explrn, nktotf))
    ALLOCATE(eigval_ex_temp(negnv_explrn, nkc1*nkc2*nkc3))
    !
    eigval_ex = 0.d0
    eigval_ex_temp = 0.d0
    !
    ! Formatted reading
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
      ! Convert from eV to Ry
      DO iexq=1, nkc1*nkc2*nkc3
        READ(iuexphg, *) eigval_ex_temp(:, iexq)
        DO imode=1,negnv_explrn
          eigval_ex_temp(imode, iexq) = eigval_ex_temp(imode, iexq) / ryd2ev 
        ENDDO
      ENDDO
      !
      CLOSE(iuexphg)
      !
      iq = 1
      DO ik1 = 1, nkc1, step_k1_explrn
        DO ik2 = 1, nkc2, step_k2_explrn 
          DO ik3 = 1, nkc3, step_k3_explrn
            iexq = (ik1 - 1) * nkc2 * nkc3 + (ik2 - 1) * nkc3 + ik3
            eigval_ex(:, iq) = eigval_ex_temp(:, iexq)
            iq = iq + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    ! End of formatted reading
    CALL mp_bcast(eigval_ex, ionode_id, inter_pool_comm)
    !
    DEALLOCATE(eigval_ex_temp)
    !
    CALL mp_barrier(inter_pool_comm)
  !-----------------------------------------------------------------------
  END SUBROUTINE read_exband
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_phband()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, nkc1, nkc2, nkc3, step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : wf, nktotf
    USE ep_constants,  ONLY : ryd2ev
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE      
    !
    CHARACTER(LEN = 256)        :: filename
    !! File name
    INTEGER                     :: iexq
    !! Index of the phonon momentum
    INTEGER                     :: imode
    !! Index of the phonon mode
    INTEGER                     :: ierr
    !! Error status
    INTEGER                     :: ik1, ik2, ik3
    !! Index of qx, qy, and qz
    INTEGER                     :: iq
    !! Index of the phonon momentum if some q points are skipped
    REAL(KIND=DP), ALLOCATABLE  :: wf_temp(:,:)
    !! Temporary array for storing the phonon eigenvalues
    !
    filename = './G_full_epmatq/phband.fmt'
    !
    DEALLOCATE(wf)
    ALLOCATE(wf(nmodes, nktotf))
    ALLOCATE(wf_temp(nmodes, nkc1*nkc2*nkc3))
    wf = 0.d0
    wf_temp = 0.d0
    !
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
      !
      DO iexq=1,nkc1*nkc2*nkc3
        READ(iuexphg, *) wf_temp(:, iexq)
      ENDDO
      CLOSE(iuexphg)
      !
      iq = 1
      DO ik1 = 1, nkc1, step_k1_explrn
        DO ik2 = 1, nkc2, step_k2_explrn 
          DO ik3 = 1, nkc3, step_k3_explrn
            iexq = (ik1 - 1) * nkc2 * nkc3 + (ik2 - 1) * nkc3 + ik3
            wf(:, iq) = wf_temp(:, iexq)
            iq = iq + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    CALL mp_bcast(wf, ionode_id, inter_pool_comm)
    !
    DEALLOCATE(wf_temp)
    !
    CALL mp_barrier(inter_pool_comm)
  !-----------------------------------------------------------------------
  END SUBROUTINE read_phband
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_ph_eig()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, &
                              step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : wf, nktotf, ph_eig
    USE ep_constants,  ONLY : ryd2ev, ci
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE      
    !
    CHARACTER(LEN = 256)        :: filename
    !! File name
    CHARACTER(LEN = 256)        :: tp
    !! Temporary for convering integer to character
    INTEGER                     :: iexq
    !! Index of phonon q
    INTEGER                     :: imode, imode1
    !! Index of the phonon modes
    INTEGER                     :: ierr
    !! Error status
    REAL(DP)                    :: real_temp, imag_temp 
    !! Temporary variable for reading real and imaginary part of a number
    INTEGER                     :: ik1, ik2, ik3
    !! Index of qx, qy, qz
    INTEGER                     :: iq
    !! Index of the phonon momentum if some q points are skipped
    COMPLEX(KIND=DP), ALLOCATABLE  :: ph_eig_temp(:,:,:)
    !! Temporary array for storing the phonon eigenvectors
    !
    filename = './G_full_epmatq/phband.fmt'
    !
    ALLOCATE(ph_eig(nmodes, nmodes, nktotf))
    ALLOCATE(ph_eig_temp(nmodes, nmodes, nkc1*nkc2*nkc3))
    ph_eig = 0.d0
    ph_eig_temp = 0.d0
    !
    IF(my_pool_id .EQ. ionode_id) THEN
      DO iq=1, nkc1*nkc2*nkc3
        WRITE(tp,"(I10)") iq - 1  
        filename = './G_full_epmatq/ph_eigvec_' // trim(adjustl(tp)) // '.dat'
        OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
        !
        DO imode=1,nmodes
          DO imode1=1,nmodes
            READ(iuexphg, *) real_temp, imag_temp
            ph_eig_temp(imode1, imode, iq) = real_temp + ci * imag_temp
          ENDDO
        ENDDO
        !
        CLOSE(iuexphg)
      ENDDO
      !
      iq = 1
      DO ik1 = 1, nkc1, step_k1_explrn
        DO ik2 = 1, nkc2, step_k2_explrn 
          DO ik3 = 1, nkc3, step_k3_explrn
            iexq = (ik1 - 1) * nkc2 * nkc3 + (ik2 - 1) * nkc3 + ik3
            ph_eig(:, :, iq) = ph_eig_temp(:, :, iexq)
            iq = iq + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    CALL mp_bcast(ph_eig, ionode_id, inter_pool_comm)
    !
    DEALLOCATE(ph_eig_temp)
    !
    CALL mp_barrier(inter_pool_comm)
  !-----------------------------------------------------------------------
  END SUBROUTINE read_ph_eig
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_ex_ph(lower_bnd, upper_bnd)
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, g_scale_plrn, nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, &
                              step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : G_full_epmat, wf, nktotf
    USE ep_constants,  ONLY : ryd2ev, ci
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE      
    !
    INTEGER, INTENT(in)         :: lower_bnd, upper_bnd
    !! Lower and upper bound of the q in G_ss1nu(Q,q)
    CHARACTER(LEN = 256)        :: filename
    !! File name
    CHARACTER(LEN = 256)        :: tp
    !! Temporary for converting integer to character
    INTEGER                     :: nqloc
    !! Number of local q points
    INTEGER                     :: iq, iq_f
    !! Index of q in the coarser q grid because of skipping and full q grid
    INTEGER                     :: iexq, iexq_f
    !! Index of q in the coarser Q grid because of skipping and full q grid  
    INTEGER                     :: sbnd, s1bnd
    !! Index of exciton bands
    INTEGER                     :: imode
    !! Index of phonon modes
    INTEGER                     :: iexq_d, imode_d, sbnd_d, s1bnd_d
    !! Dummy variables to store uesless data while reading
    INTEGER                     :: ierr
    !! Error status
    REAL(kind=DP)               :: wf_d
    !! Dummy variables to store uesless data while reading
    REAL(kind=DP)               :: temp_real_G, temp_imag_G
    !! Temporary variable to store the real and imaginary part of ex-ph matrix
    INTEGER                     :: ik1, ik2, ik3
    !! Index of Qx, Qy, Qz
    COMPLEX(KIND=DP), ALLOCATABLE  :: G_full_epmat_temp(:,:,:,:)
    !! Temporary array for storing the ex-ph matrix elements
    !
    nqloc = upper_bnd - lower_bnd + 1
    filename = './G_full_epmatq/phband.freq'
    !
    ALLOCATE(G_full_epmat(negnv_explrn, negnv_explrn, nktotf, nmodes, nqloc))
    ALLOCATE(G_full_epmat_temp(negnv_explrn, negnv_explrn, nkc1*nkc2*nkc3, nmodes))
    ! Notice the index change
    G_full_epmat = 0.d0
    G_full_epmat_temp = 0.d0
    !
    DO iq=1,nqloc
      CALL iq_to_iqf(iq + lower_bnd - 1, nqc1, nqc2, nqc3, &
      step_k1_explrn, step_k2_explrn, step_k3_explrn, iq_f)
      !Find the index of iq in the finer coarse grid
      WRITE(tp,"(I10)") iq_f - 1
      filename = './G_full_epmatq/G_full_epmatq_' // trim(adjustl(tp)) // '.dat'
      OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
      !
      DO iexq_f=1,nkc1*nkc2*nkc3
        DO sbnd=1,negnv_explrn
          DO s1bnd=1,negnv_explrn
            DO imode=1,nmodes
              READ(iuexphg, *) sbnd_d, s1bnd_d, iexq_d, imode_d, wf_d, &
              temp_real_G, temp_imag_G
              !
              G_full_epmat_temp(sbnd, s1bnd, iexq_f, imode) = (temp_real_G + ci * temp_imag_G) / ryd2ev
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      CLOSE(iuexphg)
      !
      iexq = 1
      DO ik1 = 1, nkc1, step_k1_explrn
        DO ik2 = 1, nkc2, step_k2_explrn 
          DO ik3 = 1, nkc3, step_k3_explrn
            iexq_f = (ik1 - 1) * nkc2 * nkc3 + (ik2 - 1) * nkc3 + ik3
            G_full_epmat(:,:,iexq,:,iq) = G_full_epmat_temp(:,:,iexq_f,:) * g_scale_plrn
            iexq = iexq + 1
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    DEALLOCATE(G_full_epmat_temp)
  !-----------------------------------------------------------------------
  END SUBROUTINE read_ex_ph
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE initialize_Asq()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, init_plrn
    USE global_var,    ONLY : G_full_epmat, wf, eigval_ex, Asq, nktotf
    USE ep_constants,  ONLY : ryd2ev, ci
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE      
    !
    CHARACTER(LEN = 256)        :: filename
    !! File name
    INTEGER                     :: iq, imode
    !! Index of exciton Q and exciton band
    INTEGER                     :: ierr
    !! Error status
    INTEGER                     :: numk, nummode
    !! Number of Q points and exciton bands
    INTEGER                     :: temp_iexq, temp_imode
    !! Dummy variables for reading purpose only
    REAL(DP)                    :: Asqnorm 
    !! Norm of Asq
    REAL(DP)                    :: Asq_real, Asq_imag
    !! Real and imaginary part of Asq
    !
    ALLOCATE(Asq(negnv_explrn, nktotf))
    !
    IF(init_plrn .EQ. 4) THEN
      WRITE(stdout, '(5x, a)') 'Initializing Asq from a previous Asq_disp.plrn file.'
      filename = 'Asq_disp.plrn'
      IF(my_pool_id .EQ. ionode_id) THEN
        OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
        READ(iuexphg, '(5I10)') temp_iexq, temp_iexq, temp_iexq, numk, nummode
        !
        DO iq=1,numk
          DO imode=1, negnv_explrn
            ! READ(iuexphg, '(2I10, 2f15.7)') temp_iexq, temp_imode, &
            ! Asq_real, Asq_imag
            READ(iuexphg, *) temp_iexq, temp_imode, &
            Asq_real, Asq_imag
            Asq(imode, iq) = Asq_real + ci * Asq_imag
          ENDDO
        ENDDO 
        !
        CLOSE(iuexphg)
        CALL normalize_Asq(Asq)
      ENDIF
      !
      CALL mp_bcast(Asq, ionode_id, inter_pool_comm)
      !
    ELSEIF(init_plrn .EQ. 5) THEN
      WRITE(stdout, '(5x, a)') 'Initializing Asq with a constant.'
      Asq = 1.0
      CALL normalize_Asq(Asq)
      ! Check whether it is correctly normalized
      Asqnorm = 0.0
      DO iq=1, nktotf
        DO imode=1,negnv_explrn
          Asqnorm = Asqnorm + ABS(Asq(imode, iq)) ** 2
        ENDDO
      ENDDO
      WRITE(stdout, '(5x, a, f10.4)') 'The norm of the initialized Asq: ', Asqnorm / nktotf
    ELSEIF(init_plrn .EQ. 6) THEN
      !! Note that this initialization will be implemented in prepare_bmat
      WRITE(stdout, '(5x, a)') 'Initializing Asq with given displacement patterns by dtau_disp.plrn.'
      Asq = 1.0
      CALL normalize_Asq(Asq)
      ! Check whether it is correctly normalized
      Asqnorm = 0.0
      DO iq=1, nktotf
        DO imode=1,negnv_explrn
          Asqnorm = Asqnorm + ABS(Asq(imode, iq)) ** 2
        ENDDO
      ENDDO
      WRITE(stdout, '(5x, a, f10.4)') 'The norm of the initialized Asq: ', Asqnorm / nktotf
    ELSE
      WRITE(stdout, '(5x, a)') 'Initialization not supported yet!'
    ENDIF
  !-----------------------------------------------------------------------  
  END SUBROUTINE initialize_Asq
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE prepare_bmat(iscf)
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, init_plrn, nqc1, nqc2, nqc3, &
                              step_k1_explrn, step_k2_explrn, step_k3_explrn, only_pos_modes_explrn
    USE global_var,    ONLY : G_full_epmat, wf, eigval_ex, Asq, nktotf
    USE ep_constants,  ONLY : ryd2ev, eps6, ci
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    USE parallelism,   ONLY : fkbounds
    USE exphonon,      ONLY : indexing_q2
    !
    IMPLICIT NONE      
    !
    CHARACTER(LEN = 256)        :: filename
    !! File name
    INTEGER, INTENT(in)         :: iscf
    !! Whether this is first iteration in the explrn scf calculations
    INTEGER                     :: nqloc
    !! Number of local q points because of skipping
    INTEGER                     :: numk, nummode
    !! Number of total q points and phonon modes
    INTEGER                     :: iq
    !! Index of q point
    INTEGER                     :: iexq, iexqplusq
    !! Index of Q point and Q+q point
    INTEGER                     :: imode, sbnd, s1bnd
    !! Index of phonon bands and exciton bands
    INTEGER                     :: lower_bnd, upper_bnd
    !! Lower and upper bound of the local phono q
    INTEGER                     :: ierr
    !! Error status
    INTEGER                     :: temp_nqc, temp_iq, temp_imode
    !! Dummy variable
    REAL(DP)                    :: disp_real, disp_imag
    !! Real and imaginary part of the atomic displacement
    !
    Bmat = 0.d0
    !
    IF (init_plrn == 6 .AND. iscf == 1) THEN
      filename = 'dtau_disp.plrn'
      ! 
      IF(my_pool_id .EQ. ionode_id) THEN
        OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
        READ(iuexphg, '(5I10)') temp_nqc, temp_nqc, temp_nqc, numk, nummode
        ! 
        DO iq=1,numk
          DO imode=1, nmodes
            ! READ(iuexphg, '(2I10, 2f15.7)') temp_iexq, temp_imode, &
            ! Asq_real, Asq_imag
            READ(iuexphg, *) temp_iq, temp_imode, &
            disp_real, disp_imag
            dtau(imode, iq) = disp_real + ci * disp_imag
          ENDDO
        ENDDO     
        CLOSE(iuexphg)
      ENDIF
      CALL mp_bcast(dtau, ionode_id, inter_pool_comm)
      !
      CALL compute_disp(dtau, Bmat, -1) 
      !
    ELSE
      CALL fkbounds(nktotf, lower_bnd, upper_bnd)
      nqloc = upper_bnd - lower_bnd + 1
      DO iq=1,nqloc 
        !$omp parallel default(shared) private(imode, iexq, iexqplusq, sbnd, s1bnd)
        !$omp do schedule(dynamic)
        DO imode=1,nmodes
          DO iexq=1,nktotf
            CALL indexing_q2(iq+lower_bnd-1, iexq, &
                 nqc1/step_k1_explrn, nqc2/step_k2_explrn, nqc3/step_k3_explrn, 1, iexqplusq)
            DO sbnd=1,negnv_explrn
              DO s1bnd=1,negnv_explrn
                Bmat(imode, iq+lower_bnd-1) = Bmat(imode, iq+lower_bnd-1) + &
                                              CONJG(Asq(s1bnd, iexq)) * Asq(sbnd, iexqplusq) * &
                                              CONJG(G_full_epmat(sbnd, s1bnd, iexq, imode, iq))
              ENDDO
            ENDDO
          ENDDO
          !
          IF(only_pos_modes_explrn) THEN
            ! IF(ABS(wf(imode, iq+lower_bnd-1)) > eps6) THEN
            IF( ((imode <= 3) .AND. (iq + lower_bnd - 1 == 1))  &
                .OR. ABS(wf(imode, iq+lower_bnd-1)) < eps6 * 10.d0) THEN
              Bmat(imode, iq+lower_bnd-1) = 0
            ELSE
              Bmat(imode, iq+lower_bnd-1) = Bmat(imode, iq+lower_bnd-1) / nktotf / wf(imode, iq+lower_bnd-1)
            ENDIF 
          ELSE
            IF(ABS(wf(imode, iq+lower_bnd-1)) < eps6 * 10.d0) THEN 
              Bmat(imode, iq+lower_bnd-1) = 0
            ELSE
              Bmat(imode, iq+lower_bnd-1) = Bmat(imode, iq+lower_bnd-1) / nktotf / ABS(wf(imode, iq+lower_bnd-1))
            ENDIF 
          ! 
          ENDIF
        ENDDO
        !$omp end do
        !$omp end parallel
      ENDDO
      !
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(Bmat, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
   ENDIF
  !-----------------------------------------------------------------------
  END SUBROUTINE prepare_bmat
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE prepare_btimesg()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, init_plrn, nqc1, nqc2, nqc3, &
                              step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : G_full_epmat, wf, eigval_ex, nktotf
    USE ep_constants,  ONLY : ryd2ev, eps6
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    USE parallelism,   ONLY : fkbounds
    USE exphonon,      ONLY : indexing_q2
    !
    IMPLICIT NONE      
    !
    INTEGER                     :: nqloc
    !! Number of local qpts
    INTEGER                     :: iq1, iq2, iq
    !! Index of Q1, Q2, and Q2-Q1        
    INTEGER                     :: sbnd, s1bnd
    !! Index of exction bands     
    INTEGER                     :: imode
    !! Index of phonon bands
    INTEGER                     :: idx1, idx2
    !! Index of the specific element in the flattened matrix
    INTEGER                     :: lower_bnd, upper_bnd
    !! Lower and upper bound of the local q points
    BtimesG = 0.d0
    !
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    nqloc = upper_bnd - lower_bnd + 1
    !
    DO iq2=1,nqloc
      !$omp parallel default(shared) private(imode, iq, iq1, sbnd, s1bnd, idx1, idx2)
      !$omp do schedule(dynamic)
      DO iq1=1,nktotf
        !iq and iq1 correspond to Q and Q' in the long paper Eq.(24), iq2=iq-iq1
        CALL indexing_q2(iq1, iq2+lower_bnd-1, &
             nqc1/step_k1_explrn, nqc2/step_k2_explrn, nqc3/step_k3_explrn, 1, iq)
        !
        DO sbnd=1,negnv_explrn
          idx1 = (iq - 1)* negnv_explrn + sbnd
          DO s1bnd=1,negnv_explrn
            idx2 = (iq1 - 1) * negnv_explrn + s1bnd
            DO imode=1,nmodes
              BtimesG(idx1, idx2) = BtimesG(idx1, idx2) + Bmat(imode, iq2+lower_bnd-1) * &
                                    G_full_epmat(sbnd, s1bnd, iq1, imode, iq2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !$omp end do
      !$omp end parallel
    ENDDO
    !
    CALL mp_barrier(inter_pool_comm)
    CALL mp_sum(BtimesG, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
  !-----------------------------------------------------------------------
  END SUBROUTINE prepare_btimesg
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE prepare_ex_hamil()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, init_plrn
    USE global_var,    ONLY : G_full_epmat, wf, eigval_ex, nktotf
    USE ep_constants,  ONLY : ryd2ev, eps6
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE      

    INTEGER                     :: iq1, iq2, iq
    !! Index of Q1, Q2, and Q2-Q1        
    INTEGER                     :: sbnd, s1bnd
    !! Index of exction bands     
    INTEGER                     :: idx1, idx2
    !! Index of the specific element in the flattened matrix
    INTEGER                     :: ierr
    !! Error status
    !
    ! Diagonalization is only done in the root rank
    ! Now only do direct diagonalization
    !
    IF(my_pool_id .EQ. ionode_id) THEN
      ALLOCATE(Hamil(negnv_explrn*nktotf,negnv_explrn*nktotf), STAT=ierr)
      IF(ierr/=0) CALL errore('prepare_ex_hamil','Error in allocating Hamil', 1)
      Hamil = 0.d0
      !$omp parallel default(shared) private(iq, iq1, sbnd, s1bnd, idx1, idx2)
      !$omp do schedule(dynamic)
      DO iq=1,nktotf
        DO sbnd=1,negnv_explrn
          idx1 = (iq - 1) * negnv_explrn + sbnd
          Hamil(idx1, idx1) = eigval_ex(sbnd, iq)
          DO iq1=1,nktotf
            DO s1bnd=1,negnv_explrn
              idx2 = (iq1 - 1) * negnv_explrn + s1bnd
              Hamil(idx1,idx2) = Hamil(idx1,idx2) - BtimesG(idx1, idx2) / nktotf * 2
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !$omp end do
      !$omp end parallel
    ENDIF
    !
    CALL mp_barrier(inter_pool_comm)
  !-----------------------------------------------------------------------
  END SUBROUTINE prepare_ex_hamil
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE compute_tot_energy()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, niter_plrn, ethrdg_plrn
    USE global_var,    ONLY : G_full_epmat, eigval_ex, wf, Asq , nktotf,&
                              explrn_eigval, explrn_eigval_old, explrn_etot
    USE ep_constants,  ONLY : ryd2ev, eps6
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE
    !
    INTEGER                     :: iq
    !! Index of phonon q point
    INTEGER                     :: imode
    !! Index of the phonon modes
    !
    IF(my_pool_id .EQ. ionode_id) THEN
      explrn_etot = 0.0
      DO iq=1,nktotf
        DO imode=1,nmodes
          explrn_etot = explrn_etot + ABS(Bmat(imode, iq)) ** 2 * wf(imode, iq)
        ENDDO
      ENDDO
      explrn_etot = explrn_etot / nktotf + explrn_eigval
    ENDIF
  !-----------------------------------------------------------------------
  END SUBROUTINE compute_tot_energy
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE solve_for_explrn(nrr_q, ndegen_q, irvec_q, rws, nrws)
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn, niter_plrn, ethrdg_plrn, &
                              nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : G_full_epmat, eigval_ex, wf, Asq , nktotf, &
                              explrn_eigval, explrn_eigval_old, explrn_etot, &
                              ph_eig
    USE ep_constants,  ONLY : ryd2ev, eps6
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER,  INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    INTEGER                     :: dimn
    !! Dimension of the flattened exciton polaron Hamiltonian
    INTEGER                     :: iq1, iq2, iq3
    !! Index of Qx, Qy, Qz
    INTEGER                     :: iexq
    !! Index of Q
    INTEGER                     :: imode
    !! Index of the exciton band    
    INTEGER                     :: iscf
    !! Index of the iteration in the exciton polaron scf calculations
    INTEGER                     :: idx1, idx2
    !! Index of the specific element in the flattened matrix
    INTEGER                     :: neig
    !! Number of eigenvectors after matrix diagonalization
    INTEGER                     :: ierr
    !! Error status
    INTEGER                     :: idx_min(2)
    !! Index of the minimum element of a 2D array
    REAL(DP)                    :: plrn_error
    !! The error between two iterations in polaron scf
    INTEGER :: info
    !! "0" successful exit, "<0" i-th argument had an illegal value, ">0" i eigenvectors failed to converge.
    INTEGER :: ifail(negnv_explrn * nktotf)
    !! Contains the indices of the eigenvectors that failed to converge
    INTEGER :: iwork(5 * negnv_explrn * nktotf)
    !! Work array for matrix diagonalization
    REAL(KIND = DP) :: rwork(7 * negnv_explrn * nktotf)
    !! Real work array
    COMPLEX(KIND = DP) :: cwork(2 * negnv_explrn * nktotf)
    !! Complex work array
    REAL(KIND = DP) :: w1(negnv_explrn * nktotf)
    !! Arrays to store matrix eigenvalues
    !
    dimn = negnv_explrn * nktotf
    !
    IF(my_pool_id == ionode_id) THEN
      ALLOCATE(eigVec(dimn, dimn))
    ENDIF 
    !
    ALLOCATE(dtau(nmodes, nktotf))    
    ALLOCATE(Bmat(nmodes, nktotf)) 
    ALLOCATE(BtimesG(negnv_explrn*nktotf,negnv_explrn*nktotf))
    dtau = 0.d0
    Bmat = 0.d0
    BtimesG = 0.d0
    !
    CALL initialize_Asq()
    !
    explrn_eigval = 1000
    explrn_eigval_old = 1000
    plrn_error=1000
    WRITE(stdout, '(5x, a)') 'Start the excitonic polaron self-consistent calculation.'
    WRITE(stdout, '(5x, a)') 'Iteration    Eigenvalue (eV)    Total Energy (eV)    Error (eV).'
    DO iscf=1,niter_plrn
      CALL prepare_bmat(iscf)
      CALL prepare_btimesg()
      CALL prepare_ex_hamil()
      !
      IF(my_pool_id == ionode_id) THEN
        ! Diagonalize the polaron Hamiltonian
        eigVec = 0.d0
        ALLOCATE(Hamil_work(dimn*(dimn+1)/2))
        Hamil_work = 0.d0
        DO idx2=1,dimn
          DO idx1=1,idx2
            !Hamil_work(idx1 + (idx2 - 1) * idx2 / 2) = Hamil(idx1, idx2)
            Hamil_work(idx1 + (idx2 - 1) * idx2 / 2) = ( Hamil(idx1, idx2) + &
                                              CONJG( Hamil(idx2, idx1) ) ) / 2.0
          ENDDO
        ENDDO
        !
        DEALLOCATE(Hamil)
        !
        CALL ZHPEVX('V', 'A', 'U', dimn, Hamil_work, 0.0, 0.0, 0, 0, -1.0, neig, w1, eigVec, dimn, cwork, &
            rwork, iwork, ifail, info)
        !
        IF(info .NE. 0) THEN
          CALL errore('expoloarn', 'Error in diagoanlizing the polaron Hamiltonian.',1)
        ENDIF
        explrn_eigval_old = explrn_eigval
        explrn_eigval = w1(1)
        !
        DO iexq=1,nktotf
          DO imode=1,negnv_explrn
            Asq(imode, iexq) = eigVec(imode+(iexq-1)*negnv_explrn, 1)
          ENDDO
        ENDDO
        CALL normalize_Asq(Asq)
        DEALLOCATE(Hamil_work) 
        ! Compute the polaron total energy 
        CALL compute_tot_energy()
        !
      ENDIF
      !
      CALL mp_barrier(inter_pool_comm)
      CALL mp_bcast(Asq, ionode_id, inter_pool_comm)
      CALL mp_bcast(explrn_etot, ionode_id, inter_pool_comm)
      CALL mp_bcast(explrn_eigval, ionode_id, inter_pool_comm)
      CALL mp_bcast(explrn_eigval_old, ionode_id, inter_pool_comm)
      !
      plrn_error = ABS(explrn_eigval - explrn_eigval_old)
      !
      IF(iscf .EQ. 1) THEN
        WRITE(stdout, '(5x, I5, 2F20.8)') iscf, explrn_eigval * ryd2ev, &
                                        (explrn_etot - MINVAL(eigval_ex))*ryd2ev
      ELSE
        WRITE(stdout, '(5x, I5, 3F20.8)') iscf, explrn_eigval * ryd2ev, &
                                       (explrn_etot - MINVAL(eigval_ex))*ryd2ev, plrn_error * ryd2ev
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      !
      IF( (plrn_error * ryd2ev * 1000 <= ethrdg_plrn) .OR. &
           iscf == niter_plrn ) THEN
        WRITE(stdout, '(5x, a, 3F20.8)') &
        'Convergence achieved. The converged eigenvalue and total energy are: ', &
        explrn_eigval  * ryd2ev, (explrn_etot - MINVAL(eigval_ex))*ryd2ev
        !
        idx_min = MINLOC(eigval_ex)
        !
        iq3 = MOD(idx_min(2) - 1, nqc3 / step_k3_explrn)
        iq2 = MOD( (idx_min(2) - iq3 - 1) / (nqc3 / step_k3_explrn), nqc2 /step_k2_explrn)
        iq1 = (idx_min(2) - 1) / (nqc2 * nqc3 / step_k2_explrn / step_k3_explrn)
        !
        WRITE(stdout, '(5x, a, 3F20.8)') &
        'Q point holding the lowest-energy exciton: ', &
        (iq1 + 0.0) / (nqc1 / step_k1_explrn), &
        (iq2 + 0.0) / (nqc2 / step_k2_explrn), &
        (iq3 + 0.0) / (nqc3 / step_k3_explrn)  
        !
        WRITE(stdout, '(5x, a, 1F20.8)') &
        'Total energy referenced to the lowest-energy exciton at Q=0: ', &
        (explrn_etot - (eigval_ex(1,1)))*ryd2ev
        !
        CALL compute_disp(Bmat, dtau, 1)
        !
        CALL save_solution(Asq, nktotf, negnv_explrn, 'Asq.plrn')
        ! Save the previous Bmat so that Asq can be reproduced
        ! CALL save_solution(Bmat, nktotf, nmodes, 'bmat.plrn')
        IF(my_pool_id == ionode_id) CALL write_explrn_bmat(Bmat, 'bmat.plrn', wf, nktotf)
        !    
        CALL save_solution(dtau, nktotf, nmodes, 'dtau.plrn')
        !
        IF(my_pool_id == ionode_id) THEN
          CALL write_explrn_dtau_xsf('dtau.plrn.xsf', nqc1 / step_k1_explrn, nqc2 / step_k2_explrn, &
              nqc3 / step_k3_explrn, .TRUE.)
        ENDIF
        !
        CALL interp_explrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
        ! Note that at this point the wf has been changed to correspond to the fine grid
        CALL mp_barrier(inter_pool_comm)
        EXIT
      ENDIF
    ENDDO
    !
    DEALLOCATE(eigval_ex)
    DEALLOCATE(wf)
    DEALLOCATE(ph_eig)
    DEALLOCATE(G_full_epmat)
    DEALLOCATE(Asq)
    DEALLOCATE(Bmat)
    DEALLOCATE(BtimesG)
    DEALLOCATE(dtau)
    IF(my_pool_id == ionode_id) DEALLOCATE(eigVec)
  !-----------------------------------------------------------------------  
  END SUBROUTINE solve_for_explrn
  !-----------------------------------------------------------------------
  !----------------------------Some utilities-----------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE write_explrn_dtau_xsf(filename, nqc1, nqc2, nqc3, write_force)
  !-----------------------------------------------------------------------
    !! Adapted from write_plrn_dtau_xsf from polaron.f90 !!
    USE ep_constants,   ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE io_global,      ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,       ONLY : world_comm
    USE mp,             ONLY : mp_sum, mp_bcast
    USE modes,          ONLY : nmodes
    USE global_var,     ONLY : nktotf
    USE ions_base,      ONLY : nat, amass, ityp, tau, atm, nsp, na, ntypx
    USE cell_base,      ONLY : at, alat
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! File name
    INTEGER, INTENT(in) :: nqc1, nqc2, nqc3    
    !! Dimensions of the phonon q grid
    LOGICAL, INTENT(in) :: write_force
    !! Whether to write force vectors
    INTEGER :: nktotf_temp
    !! Temporary variable to store nktotf since it might change
    INTEGER :: ierr
    !! Error status
    INTEGER :: nat_all
    !! Number of atoms in the supercell
    INTEGER :: iRp
    !! Index of the unit cell in the supercell
    INTEGER :: iatm
    !! Index of the atom in the unit cell
    INTEGER :: iatm_all
    !! Index of the atom in the supercell
    INTEGER :: ika
    !! Index of the atom and its cartesian direction
    INTEGER :: iq, iq1, iq2, iq3
    !! Index of q, qx, qy, qz
    INTEGER :: Rp_vec(1:3)
    !! Coordinate of unit cell in the supercell
    INTEGER, ALLOCATABLE :: elements(:)
    !! Index of species of the atom
    INTEGER :: species(50)
    !! Arrays to store all species index
    INTEGER :: n_grid(3)
    !! The dimension of the real-space grid
    INTEGER :: grid_start(3)
    !! Starting points of the real-space grid
    INTEGER :: grid_end(3)
    !! Ending points of the real-space grid
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell lattice parameters
    REAL(KIND = DP) :: shift(1:3)
    !! The shift vector that brings the atoms from unit cell to supercell
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! The positions of the atoms in the supercell
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! The displacements of each atom after polaron formation
    REAL(DP), ALLOCATABLE        ::  Rp(:,:)
    !! Coordinate of all unit cells in the supercell
    !
    nktotf_temp = nktotf
    ! total number of atoms is
    ! (number of atoms in unit cell) x (number of cells in the supercell)
    nktotf =  nqc1 * nqc2 * nqc3 
    nat_all = nat * nktotf
    !
    ALLOCATE(atoms(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating atoms', 1)
    ALLOCATE(elements(nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating elements', 1)
    ALLOCATE(displacements(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating displacements', 1)
    ALLOCATE(Rp(nktotf, 3))
    !
    iq = 1
    DO iq1=1,nqc1 
      DO iq2=1,nqc2
        DO iq3=1,nqc3 
          Rp(iq, :) = (/iq1-1.0, iq2-1.0, iq3-1.0/)
          ! WRITE(stdout, '(a, 6F20.10)') "Check grids: ", xqc_cryst(iq, :), Rp(iq, :) 
          iq = iq + 1
        ENDDO
      ENDDO
    ENDDO
    atoms = zero
    elements = 0
    displacements = 0.0
    !
    cell(1:3, 1) = at(1:3, 1) * nqc1 
    cell(1:3, 2) = at(1:3, 2) * nqc2 
    cell(1:3, 3) = at(1:3, 3) * nqc3 
    !
    ! Determine species from .cube files 
    CALL read_species_wannier_cube(species, n_grid, grid_start, grid_end)
    !
    iatm_all = 1
    DO iRp = 1, nktotf
      Rp_vec(1:3) = Rp(iRp, :)
      DO iatm = 1, nat        
        ika = (iatm - 1) * 3 + 1
        shift(1:3) = Rp_vec(1) * at(1:3, 1) + Rp_vec(2) * at(1:3, 2) + Rp_vec(3) * at(1:3, 3)        
        elements(iatm_all) = species(ityp(iatm))
        atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
        !
        IF(write_force) displacements(1:3, iatm_all) = REAL(dtau(ika:ika + 2, iRp))
        !
        iatm_all = iatm_all + 1
      ENDDO
    ENDDO
    !
    cell = cell * alat
    atoms = atoms * alat
    !
    IF(write_force) THEN
      CALL write_explrn_xsf_file(filename, cell * bohr2ang, elements, &
          atoms * bohr2ang, nat_all, displacements * bohr2ang)
    ELSE
      CALL write_explrn_xsf_file(filename, cell * bohr2ang, elements, &
          atoms * bohr2ang, nat_all)
    ENDIF
    !
    nktotf = nktotf_temp
    !
    DEALLOCATE(atoms, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating atoms', 1)
    DEALLOCATE(elements, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating elements', 1)
    DEALLOCATE(displacements, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating displacements', 1)
    !
    DEALLOCATE(Rp)
  !-----------------------------------------------------------------------
  END SUBROUTINE write_explrn_dtau_xsf
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE write_explrn_xsf_file(filename, cell, elements, atoms, natm, forces)
  !-----------------------------------------------------------------------
    !! Adapted from write_xsf_file from polaron.f90
    USE io_var,        ONLY : iuindabs, iuexphg
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in):: filename
    !! File name
    INTEGER, INTENT(in) :: elements(:)
    !! Index of species of the atom
    REAL(KIND = DP), INTENT(in) :: cell(3, 3)
    !! Supercell lattice parameters
    REAL(KIND = DP), INTENT(in) :: atoms(:, :)
    !! The positions of the atoms in the supercell
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: forces(:, :)
    !! The atomic displacements
    INTEGER, INTENT(in) :: natm
    !! Number of atoms in the supercell
    INTEGER :: iatm
    !! Index of one atom in the supercell
    !
    OPEN(UNIT = iuexphg, FILE = TRIM(filename), FORM = 'formatted', STATUS = 'unknown')
    !
    WRITE(iuexphg, '(a)') '#'
    WRITE(iuexphg, '(a)') '# Generated by the EPW polaron code'
    WRITE(iuexphg, '(a)') '#'
    WRITE(iuexphg, '(a)') '#'
    WRITE(iuexphg, '(a)') 'CRYSTAL'
    WRITE(iuexphg, '(a)') 'PRIMVEC'
    WRITE(iuexphg, '(3f12.7)') cell(1:3, 1)
    WRITE(iuexphg, '(3f12.7)') cell(1:3, 2)
    WRITE(iuexphg, '(3f12.7)') cell(1:3, 3)
    WRITE(iuexphg, '(a)') 'PRIMCOORD'
    ! The second number is always 1 for PRIMCOORD coordinates,
    ! according to http://www.xcrysden.org/doc/XSF.html
    WRITE (iuexphg, '(2I10)')  natm, 1
    !
    DO iatm = 1, natm
      IF (PRESENT(forces)) THEN
        WRITE(iuexphg,'(I3, 3x, 3f15.9, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm), forces(1:3, iatm)
      ELSE
        WRITE(iuexphg,'(I3, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm)
      ENDIF
    ENDDO
    !
    CLOSE(iuexphg)
  !-----------------------------------------------------------------------
  END SUBROUTINE write_explrn_xsf_file
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_species_wannier_cube(species, n_grid, &
    grid_start_min, grid_end_max)
  !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    !! Read the the species from Wannier function from prefix_0000n.cube file
    !! Adapted from read_wannier_cube from polaron.f90
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone
    USE io_var,        ONLY : iun_plot
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    USE input,         ONLY : nbndsub
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: species(50)
    !! Atomic species in unit cell
    INTEGER, INTENT(out) :: n_grid(3)
    !! The dimension of the real-space grid
    INTEGER, INTENT(out) :: grid_start_min(3)
    !! Min value of the starting points of the real-space grid
    INTEGER, INTENT(out) :: grid_end_max(3)
    !! Max value of the ending points of the real-space grid
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! Wannier functions in real space grid
    CHARACTER(LEN = 60) :: wancube
    !! File name of the cube files that store Wannier functions
    CHARACTER(LEN = 60) :: temp_str
    !! Temporary character variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd
    !! Index of the Wannier functions
    INTEGER :: ie
    !! Index of species
    INTEGER :: idir
    !! Index of cartesian directions
    INTEGER :: i_species
    !! Index of species (ZD: whats the difference between ie?)
    INTEGER :: iline
    !! Auxiliary counter
    INTEGER :: nAtoms
    !! Number of atoms contained in the cube file
    INTEGER :: grid_start(3)
    !! Starting points of the real-space grid
    INTEGER :: grid_end(3)
    !! Ending points of the real-space grid
    REAL(KIND = DP) :: rtempvec(4)
    !! Temporary vector to be read from .cube file
    REAL(KIND = DP) :: norm
    !! Norm square of the Wannier functions
    ! find the max and min of real space grid of Wannier functions of all Wannier orbitals
    !
    grid_start_min(:) = 100000
    grid_end_max(:) = -100000
    DO ibnd = 1, nbndsub
      WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
      OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
      READ(iun_plot, *) temp_str !, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str
      READ(iun_plot, *) n_grid, grid_start, grid_end
      DO idir = 1, 3
        IF (grid_start_min(idir) >= grid_start(idir)) grid_start_min(idir) = grid_start(idir)
        IF (grid_end_max(idir) <= grid_end(idir))   grid_end_max(idir) = grid_end(idir)
      ENDDO
      CLOSE(iun_plot)
    ENDDO
    !
    ! Read the xth Wannier functions from prefix_0000x.cube in ionode
    ! and broadcast to all nodes
    ALLOCATE(wann_func(grid_start_min(1):grid_end_max(1), &
      grid_start_min(2):grid_end_max(2), &
      grid_start_min(3):grid_end_max(3), nbndsub), STAT = ierr)
    !
    IF (ierr /= 0) CALL errore('read_wannier_cube', 'Error allocating wann_func', 1)
    wann_func = zero
    species = 0
    DO ibnd = 1, nbndsub
      WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
      OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
      READ(iun_plot, *) temp_str
      READ(iun_plot, *) n_grid, grid_start, grid_end
      READ(iun_plot, *) nAtoms, rtempvec(1:3)
      !
      DO iline = 1, 3
        READ(iun_plot, '(8A)') temp_str
      ENDDO
      ie = 1
      DO iline = 1, nAtoms
        READ(iun_plot, '(i4, 4f13.5)') i_species, rtempvec
        IF (iline == 1 ) THEN
          species(ie) = i_species
          ie = ie + 1
        ELSE IF (species(ie - 1) /= i_species) THEN
          species(ie) = i_species
          ie = ie + 1
        ENDIF
      ENDDO      
      CLOSE(iun_plot)
    ENDDO
    DEALLOCATE(wann_func)
  !-----------------------------------------------------------------------
  END SUBROUTINE read_species_wannier_cube
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_explrn_wannier_cube(wann_func, species, n_grid, &
    grid_start_min, grid_end_max)
  !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    !! Read the nth Wannier function from prefix_0000n.cube file
    !! Adapted from read_wannier_cube from polaron.f90
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone, ci
    USE io_var,        ONLY : iun_plot
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    USE input,         ONLY : nbndsub
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: species(50)
    !! Atomic species
    INTEGER, INTENT(out) :: n_grid(3)
    !! Number of points in real space grid where Wannier functions are written
    INTEGER, INTENT(out) :: grid_start_min(3)
    !! Initial grid point within this pool
    INTEGER, INTENT(out) :: grid_end_max(3)
    !! Final grid point within this pool
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(out) :: wann_func(:, :, :, :)
    !! Wannier functions
    REAL(KIND = DP), ALLOCATABLE :: wann_func_r(:, :, :, :)
    !! Real part of the Wannier functions
    REAL(KIND = DP), ALLOCATABLE :: wann_func_i(:, :, :, :)
    !! Imaginary-part of the Wannier functions
    CHARACTER(LEN = 60) :: wancube
    !! Name of file containing Wannier function
    CHARACTER(LEN = 60) :: temp_str
    !! Temporary string
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: ie
    !! Atomic species counter
    INTEGER :: idir
    !! Cartesian direction counter
    INTEGER :: i_species
    !! Atomic species index
    INTEGER :: nbnd
    !! Number of Wannier funcions in which polaron wave function is expanded
    INTEGER :: iline
    !! Counter along line
    INTEGER :: nAtoms
    !! Total number of atoms
    INTEGER :: nxx, nyy, nzz
    !! Number of grid points in Cartesian directions
    INTEGER :: n_len_z
    !! Number of grid points within this pool
    INTEGER :: grid_start(3)
    !! Initial grid point within this loop
    INTEGER :: grid_end(3)
    !! Final grid point within this loop
    REAL(KIND = DP) :: rtempvec(4)
    !! Temporary vector to be read from .cube file
    REAL(KIND = DP) :: norm
    !! Norm square of the Wannier function in the real space grid
    !
    ! find the max and min of real space grid of Wannier functions of all Wannier orbitals
    IF (ionode) THEN
      grid_start_min(:) = 100000
      grid_end_max(:) = -100000
      DO ibnd = 1, nbndsub
        WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str !, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        DO idir = 1, 3
         IF (grid_start_min(idir) >= grid_start(idir)) grid_start_min(idir) = grid_start(idir)
          IF (grid_end_max(idir) <= grid_end(idir))   grid_end_max(idir) = grid_end(idir)
        ENDDO
        CLOSE(iun_plot)
      ENDDO
    ENDIF
    !
    CALL mp_bcast(n_grid,         meta_ionode_id, world_comm)
    CALL mp_bcast(grid_start_min, meta_ionode_id, world_comm)
    CALL mp_bcast(grid_end_max,   meta_ionode_id, world_comm)
    !
    ! Read the xth Wannier functions from prefix_0000x.cube in ionode
    ! and broadcast to all nodes
    ALLOCATE(wann_func(grid_start_min(1):grid_end_max(1), &
             grid_start_min(2):grid_end_max(2), &
             grid_start_min(3):grid_end_max(3), nbndsub), STAT = ierr)
    !
    IF (ierr /= 0) CALL errore('read_wannier_cube', 'Error allocating wann_func', 1)
    !
    ALLOCATE(wann_func_r(grid_start_min(1):grid_end_max(1), &
             grid_start_min(2):grid_end_max(2), &
             grid_start_min(3):grid_end_max(3), nbndsub), STAT = ierr)
    !
    IF (ierr /= 0) CALL errore('read_wannier_cube', 'Error allocating wann_func', 1)
    !
    ALLOCATE(wann_func_i(grid_start_min(1):grid_end_max(1), &
             grid_start_min(2):grid_end_max(2), &
             grid_start_min(3):grid_end_max(3), nbndsub), STAT = ierr)
    !
    IF (ierr /= 0) CALL errore('read_wannier_cube', 'Error allocating wann_func', 1)
    !
    wann_func = czero
    wann_func_r = zero
    wann_func_i = zero
    species = 0
    IF (ionode) THEN
      DO ibnd = 1, nbndsub
        ! Read the real part
        WRITE(wancube, "(a, '_', i5.5, 'r.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        READ(iun_plot, *) nAtoms, rtempvec(1:3)
        !
        DO iline = 1, 3
          READ(iun_plot, '(8A)') temp_str
        ENDDO
        ie = 1
        DO iline = 1, nAtoms
          READ(iun_plot, '(i4, 4f13.5)') i_species, rtempvec
          IF (iline == 1 ) THEN
            species(ie) = i_species
            ie = ie + 1
          ELSE IF (species(ie - 1) /= i_species) THEN
            species(ie) = i_species
            ie = ie + 1
          ENDIF
        ENDDO
        n_len_z = grid_end(3) - grid_start(3) + 1
        !
        DO nxx = grid_start(1), grid_end(1)
          DO nyy = grid_start(2), grid_end(2)
            DO nzz = grid_start(3), grid_end(3), 6
              IF (grid_end(3) - nzz < 6) THEN
                READ(iun_plot, *) wann_func_r(nxx, nyy, nzz:grid_end(3) - 1, ibnd)
              ELSE
                READ(iun_plot, '(6E13.5)') wann_func_r(nxx, nyy, nzz:nzz + 5, ibnd)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iun_plot)
        !
        ! Read the imaginary part
        WRITE(wancube, "(a, '_', i5.5, 'i.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        READ(iun_plot, *) nAtoms, rtempvec(1:3)
        !
        DO iline = 1, 3
          READ(iun_plot, '(8A)') temp_str
        ENDDO
        ie = 1
        DO iline = 1, nAtoms
          READ(iun_plot, '(i4, 4f13.5)') i_species, rtempvec
          IF (iline == 1 ) THEN
            species(ie) = i_species
            ie = ie + 1
          ELSE IF (species(ie - 1) /= i_species) THEN
            species(ie) = i_species
            ie = ie + 1
          ENDIF
        ENDDO
        n_len_z = grid_end(3) - grid_start(3) + 1
        !
        DO nxx = grid_start(1), grid_end(1)
          DO nyy = grid_start(2), grid_end(2)
            DO nzz = grid_start(3), grid_end(3), 6
              IF (grid_end(3) - nzz < 6) THEN
                READ(iun_plot, *) wann_func_i(nxx, nyy, nzz:grid_end(3) - 1, ibnd)
              ELSE
                READ(iun_plot, '(6E13.5)') wann_func_i(nxx, nyy, nzz:nzz + 5, ibnd)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iun_plot)
        !
        ! Wannier function is not well normalized
        ! Normalize here will make the calculations with Wannier functions easier
        wann_func(:, :, :, ibnd) = wann_func_r(:, :, :, ibnd) + ci * wann_func_i(:, :, :, ibnd)
        norm = SUM(CONJG(wann_func(:, :, :, ibnd)) * wann_func(:, :, :, ibnd))
        wann_func(:, :, :, ibnd) = wann_func(:, :, :, ibnd) / SQRT(norm)
      ENDDO
    ENDIF
    CALL mp_bcast(wann_func, meta_ionode_id, world_comm)
    CALL mp_bcast(species, meta_ionode_id, world_comm)
    !
    DEALLOCATE(wann_func_r)
    DEALLOCATE(wann_func_i)
  !-----------------------------------------------------------------------
  END SUBROUTINE read_explrn_wannier_cube
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE compute_disp(mat_in, mat_out, itype, nrr_q, ndegen_q, irvec_q, rws, nrws, xxq)
  !-----------------------------------------------------------------------
    !! Adapted from the plrn_bmat_tran in the polaron.f90
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE ions_base,     ONLY : amass, ityp, atm, nsp
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : wf, nktotf, ph_eig, nqtotf
    USE ep_constants,  ONLY : ryd2ev, eps6, eps8, twopi, ci, cone, two, czero
    USE input,         ONLY : nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn, only_pos_modes_explrn
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    USE parallelism,   ONLY : fkbounds
    USE wannier2bloch, ONLY : dynwan2bloch
    !
    IMPLICIT NONE
    !
    INTEGER,            INTENT(in)          :: itype
    !! itype=1 -> Bqv to dtau, itype=-1 -> dtau -> Bqv
    COMPLEX(KIND = DP), INTENT(in)          :: mat_in(:, :)
    !! Input matrix
    INTEGER, INTENT(in), OPTIONAL           :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in), OPTIONAL           :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in), OPTIONAL           :: irvec_q(:, :)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in), OPTIONAL           :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in), OPTIONAL   :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    REAL(KIND = DP), OPTIONAL,  INTENT(in)  :: xxq(:,:)
    !! The list of q points for interpolation
    COMPLEX(KIND = DP), INTENT(out)         :: mat_out(:, :)
    !! Output matrix    
    INTEGER                                 :: iq, iq1, iq2, iq3
    !! Index of q, qx, qy, qz    
    INTEGER                                 :: inu
    !! Index of phonon modes
    INTEGER                                 :: ierr
    !! Error status
    INTEGER                                 :: lower_bnd, upper_bnd
    !! Lower and upper bound for parallelization
    INTEGER                                 :: start_modes
    !! Starting phonon modes when transforming from Bqv to dtau
    INTEGER                                 :: ika
    !! (kappa, alpha)
    INTEGER                                 :: iRp
    !! Index of unit cell in the supercell
    INTEGER                                 :: ina
    !! Index of atoms in the unit cell
    INTEGER                                 :: nqtot
    !! Number of q points to be considered. May be different than uniform grid if a qlist exists
    REAL(DP)                                :: Rp(nktotf, 3)
    !! List of unit cell vectors within supercell
    REAL(DP), ALLOCATABLE                   :: xqc_cryst(:, :) 
    !! Crystal coordinate of all q piints
    COMPLEX(KIND = DP)                      :: expTable(3)
    !! Table for exponential factors
    COMPLEX(KIND = DP)                      :: ctemp, dtemp
    !! Dummy variables
    COMPLEX(KIND = DP), ALLOCATABLE         :: uf(:, :, :)
    !! Phonon eigenvectors of all q points 
    REAL(KIND = DP)                         :: w2(nmodes)
    !! Phonon frequencies for a given q point
    REAL(KIND = DP)                         :: xxq_r(3)
    !! The coordinate of a single q point in the q list
    !
    mat_out = 0.d0
    !
    IF(PRESENT(xxq)) THEN
      nqtot = nqtotf
    ELSE
      nqtot = nktotf
    ENDIF
    !
    ALLOCATE(xqc_cryst(nqtot, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_disp', 'Error allocating xqc_cryst', 1)
    !
    ALLOCATE(uf(nmodes, nmodes, nqtot), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_disp', 'Error allocating uf', 1)
    !
    xqc_cryst = 0.d0
    uf = czero
    !
    IF(PRESENT(xxq)) THEN
      DEALLOCATE(wf)
      ALLOCATE(wf(nmodes, nqtot))
      wf = 0.d0
      !
      DO iq = 1, nqtotf
        xqc_cryst(iq, :) = xxq(:, iq)        
        xxq_r = xxq(:, iq)
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf(:,:,iq), w2, .false.)
        !
        DO inu = 1, nmodes
          IF(w2(inu) >=0.d0) THEN
            wf(inu, iq) =  DSQRT(ABS(w2(inu)))
          ELSE
            wf(inu, iq) = -DSQRT(ABS(w2(inu)))
          ENDIF
        ENDDO
      ENDDO
      !
      iq = 1
      DO iq1=1,nqc1 / step_k1_explrn
        DO iq2=1,nqc2 / step_k2_explrn
          DO iq3=1,nqc3 / step_k3_explrn
            Rp(iq, :) = (/iq1-1.0, iq2-1.0, iq3-1.0/)
            iq = iq + 1
          ENDDO
        ENDDO
      ENDDO
      !
    ELSE
      uf(:,:,:) = ph_eig(:,:,:)
      iq = 1
      DO iq1=1,nqc1 / step_k1_explrn
        DO iq2=1,nqc2 / step_k2_explrn
          DO iq3=1,nqc3 / step_k3_explrn
            xqc_cryst(iq, :) = (/(iq1 - 1.0)/nqc1*step_k1_explrn,&
               (iq2 - 1.0)/nqc2*step_k2_explrn, (iq3 - 1.0)/nqc3*step_k3_explrn/)
            Rp(iq, :) = (/iq1-1.0, iq2-1.0, iq3-1.0/)
            ! WRITE(stdout, '(a, 6F20.10)') "Check grids: ", xqc_cryst(iq, :), Rp(iq, :) 
            iq = iq + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    DO iq = 1, nqtot ! iq -> q
      expTable(1:3) = EXP(twopi * ci * xqc_cryst(iq, 1:3))
      start_modes = 1
      DO inu = start_modes, nmodes ! inu -> nu
        ! WRITE(stdout, *) "Checking nktotf, lower_bnd, upper_bnd, nmodes: ", nktotf, lower_bnd, upper_bnd, nmodes
        ! IF (wf(inu, iq) < eps6 * 1e10 .OR. ((inu <= 3) .AND. (iq == 1))) CYCLE ! cycle zero and imaginary frequency modes
        IF (only_pos_modes_explrn .EQV. .true.) THEN
          IF (wf(inu, iq) < eps6 * 10.d0 ) CYCLE ! cycle zero and imaginary frequency modes
        ELSE 
          IF (ABS(wf(inu, iq)) < eps6 * 10.d0 ) CYCLE ! cycle zero and imaginary frequency modes 
        ENDIF
        !
        DO ika = 1, nmodes ! ika -> kappa alpha
          ina = (ika - 1) / 3 + 1
          ctemp = DSQRT(two / (ABS(wf(inu, iq)) * amass(ityp(ina))))
          ! Parallel run, only calculate the local cell ip
          ! Note that, ip_end obtained from fkbounds should be included
          ! If you have 19 kpts and 2 pool,
          ! lower_bnd= 1 and upper_bnd=10 for the first pool
          ! lower_bnd= 1 and upper_bnd=9 for the second pool
          DO iRp = lower_bnd, upper_bnd !, (nqf1_p + 1)/2
            ! D_{\kappa\alpha\nu,p}(q) = e_{\kappa\alpha,\nu}(q) \exp(iq\cdot R_p)
            !dtemp = uf_q(ika, inu, iq) * PRODUCT(expTable(1:3)**Rp_vec(1:3))
            dtemp = uf(ika, inu, iq) * PRODUCT(expTable(1:3)**Rp(iRp, :))
            IF (itype == 1) THEN ! Bqv -> dtau
              ! \Delta \tau_{\kappa\alpha p} = -\frac{1}{N_p} \sum_{q\nu} C_{\kappa\nu q} D_{\kappa\alpha\nu q}  B_{q\nu}
              ! Dtau(iRp, ika) = Dtau(iRp, ika) + conjg(B(iq, inu)) * ctemp * dtemp
              mat_out(ika, iRp) = mat_out(ika, iRp) -  cone / REAL(nktotf, dp) * dtemp * ctemp &
                                 * mat_in(inu, iq)
            ELSE IF(itype == -1) THEN
              !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
              mat_out(inu, iq) = mat_out(inu, iq) - CONJG(dtemp) / ctemp * mat_in(ika, iRp) !JLB: dtau should be real but just in case
              !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! sum all the cell index ip
    CALL mp_sum(mat_out, inter_pool_comm)
    !
    DEALLOCATE(uf)
    DEALLOCATE(xqc_cryst)
  !-----------------------------------------------------------------------
  END SUBROUTINE compute_disp
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE interp_explrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
  !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    !! Interpolate bmat and write to bmat.band.plrn,
    !! especially used to visualize phonon contribution to ex-polaron in band-mode
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : xqf, wf, nqtotf, nktotf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : czero, ci
    USE io_global,     ONLY : stdout, ionode
    USE io_var,        ONLY : iuexphg
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id, ionode_id
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)         :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in)         :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in)         :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER,  INTENT(in)        :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    INTEGER                     :: ierr
    !! Error status
    INTEGER                     :: itemp
    !! Dummy variable
    INTEGER                     :: iq
    !! Index of unit cell in supercell
    INTEGER                     :: imode
    !! Index of atom in one unit cell
    REAL(KIND = DP)             :: disp_real, disp_imag
    !! Real and imaginary part of the displacement
    !
    ! read dtau.plrn, get the displacement.
    DEALLOCATE(dtau)
    ALLOCATE(dtau(nmodes, nktotf))
    dtau = 0.d0
    !
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE='dtau.plrn', FORM='formatted')
      READ(iuexphg, '(5I10)') itemp, itemp, itemp, itemp, itemp
      !
      DO iq=1, nktotf
        DO imode=1, nmodes
          ! READ(iuexphg, '(2I10, 2f15.7)') temp_iexq, temp_imode, &
          ! Asq_real, Asq_imag
          READ(iuexphg, *) itemp, itemp, disp_real, disp_imag
          dtau(imode, iq) = disp_real + ci * disp_imag
        ENDDO
      ENDDO     
      CLOSE(iuexphg)
    ENDIF
    CALL mp_bcast(dtau, ionode_id, inter_pool_comm)
    !
    DEALLOCATE(Bmat)
    ALLOCATE(Bmat(nmodes, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
    Bmat = czero
    !
    CALL compute_disp(dtau, Bmat, -1, nrr_q, ndegen_q, irvec_q, rws, nrws, xqf)
    !
    IF (ionode) CALL write_explrn_bmat(Bmat, 'bmat.band.plrn', wf)
    !
    ! DEALLOCATE(dtau, STAT = ierr)
    ! IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau', 1)
    ! DEALLOCATE(Bmat, STAT = ierr)
    ! IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating Bmat', 1)
    !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END SUBROUTINE interp_explrn_bq
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE plot_explrn_densities()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg, iunukk
    USE ions_base,     ONLY : amass, ityp, atm, nsp
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nbndsub, nqc1, nqc2 ,nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE global_var,    ONLY : wf, nktotf, ph_eig, Asq, nktotf, A_cv_all, A_cvq
    USE ep_constants,  ONLY : ryd2ev, eps6, twopi, ci, cone, two, czero
    USE input,         ONLY : nqc1, nqc2, nqc3, nbndc_explrn, nbndv_explrn, negnv_explrn, filukk
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    USE parallelism,   ONLY : fkbounds
    USE exphonon,      ONLY : indexing_q2, read_Acv_all
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256)             :: filename
    !! File name
    INTEGER                          :: nqloc, nummode, temp_iexq, temp_imode
    !! Dummy variable for reading files
    INTEGER                          :: ik, iq
    !! Counters for k and q
    INTEGER                          :: iq_f
    !! Index of q in the finer coarse grid
    INTEGER                          :: nktotf_f
    !! The real number of kpts in the fine grid without skipping
    INTEGER                          :: ikplusq, ikminusq
    !! Index of k+q and k-q
    INTEGER                          :: iex
    !! Index of exciton band
    INTEGER                          :: imode
    !! Index of exciton band. The same as iex
    INTEGER                          :: ibnd, jbnd
    !! Counter for electronic bands
    INTEGER                          :: iq1, iq2, iq3
    !! Index of qx, qy, qz
    INTEGER                          :: ic, iv
    !! Index of conduction band and valence band
    INTEGER                          :: im
    !! Index of Wannier functions
    INTEGER                          :: iRp
    !! Index of unit cell in the supercell
    INTEGER                          :: ierr, ios
    !! Error status
    INTEGER                          :: numk
    !! Number of Q points and 
    INTEGER                          :: icomp
    !! Flattened index of im and iRp
    INTEGER                          :: lower_bnd, upper_bnd
    !! Lower and upper bound for parallelization
    REAL(DP)                         :: Asq_real, Asq_imag
    !! Real and imaginary part of Asq
    COMPLEX(KIND = DP), ALLOCATABLE  :: cu_big(:,:,:)
    !! U(k) matrix for all k-points
    COMPLEX(KIND = DP), ALLOCATABLE  :: A_mp_vk(:, :, :, :)
    !! A_mp_vk matrix for plotting explrn wfcs
    COMPLEX(KIND = DP), ALLOCATABLE  :: A_mp_ck(:, :, :, :)
    !! A_mp_ck matrix for plotting explrn wfcs
    REAL(DP), ALLOCATABLE            :: xqc_cryst(:, :)
    !! Crystal coordinate of all q piints
    REAL(DP), ALLOCATABLE            :: Rp(:, :)
    !! List of unit cell vectors within supercell
    COMPLEX(KIND = DP)               :: expTable(3)
    !! Table for exponential factors
    COMPLEX(KIND = DP)               :: ctemp, dtemp
    !! Dummy vriables
    !
    ! Reading Asq !! 
    ALLOCATE(Asq(negnv_explrn, nktotf))
    !
    WRITE(stdout, '(5x, a)') 'Reading Asq from the previously computed Asq.plrn file.'
    filename = 'Asq.plrn'
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
      READ(iuexphg, '(5I10)') nqloc, nqloc, nqloc, numk, nummode ! nqloc are just dummies
      !
      DO iq=1,numk
        DO imode=1, negnv_explrn
          ! READ(iuexphg, '(2I10, 2f15.7)') temp_iexq, temp_imode, &
          ! Asq_real, Asq_imag
          READ(iuexphg, *) temp_iexq, temp_imode, &
          Asq_real, Asq_imag
          Asq(imode, iq) = Asq_real + ci * Asq_imag
        ENDDO
      ENDDO 
      !
      CLOSE(iuexphg)
      CALL normalize_Asq(Asq)
    ENDIF
    !
    CALL mp_bcast(Asq, ionode_id, inter_pool_comm)
    !
    ! Reading BSE eigenvectors 
    nktotf_f = nqc1 * nqc2 * nqc3
    !
    WRITE(stdout, '(5x, a)') 'Reading all BSE eigenvectors.'
    IF(my_pool_id == ionode_id) THEN
      ALLOCATE(A_cv_all(nbndv_explrn,nbndc_explrn,nktotf_f,negnv_explrn,nktotf_f),STAT=ierr)
      A_cv_all = 0.d0
      CALL read_Acv_all()
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    ! Reading the U matrix from the Wannierization 
    ! Adapted from loadumat.f90 
    !
    WRITE(stdout, '(5x, a)') 'Reading the coarse-grid U matrix.'
    ALLOCATE(cu_big(nbndsub, nbndsub, nktotf_f))
    cu_big = czero
    IF (my_pool_id .EQ. ionode_id) THEN
      !
      ! First proc read rotation matrix (coarse mesh) from file
      !
      OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
      IF (ios /=0) CALL errore('loadumat', 'error opening ukk file', iunukk)
      !
      ! dummy operation for skipping unnecessary data (ibndstart and ibndend) here
      !
      READ(iunukk, *) ibnd, jbnd
      DO ibnd = 1, nbndsub
        READ(iunukk, *) jbnd
      ENDDO
      !
      DO ik = 1, nktotf_f
        ! Since we only care about the coarse-grid quantities, we impose in the input file 
        ! that nbndep = nbndsub
        DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
            READ(iunukk, *) cu_big(ibnd, jbnd, ik)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iunukk)
      !
    ENDIF 
    !
    CALL mp_bcast(cu_big, ionode_id, world_comm)
    !
    ! End of reading BSE eigenvectors
    ! Calculate the A_mp^ck and A_mp^vk coefficients
    !
    WRITE(stdout, '(5x, a)') 'Calculating the A_mp_vk and A_mp_ck coefficients.'
    ALLOCATE(A_cvq(nbndv_explrn,nbndc_explrn,nktotf_f,negnv_explrn),STAT=ierr)
    ALLOCATE(A_mp_vk(nbndsub, nktotf_f, nbndv_explrn, nktotf_f))
    ALLOCATE(A_mp_ck(nbndsub, nktotf_f, nbndc_explrn, nktotf_f))
    !
    A_mp_vk = 0.d0
    A_mp_ck = 0.d0
    !
    ! Determine the q-grid and Rp grid
    ALLOCATE(xqc_cryst(nktotf_f, 3))
    !
    ALLOCATE(Rp(nktotf_f, 3))
    !
    iq = 1
    DO iq1=1,nqc1
      DO iq2=1,nqc2
        DO iq3=1,nqc3
          xqc_cryst(iq, :) = (/(iq1 - 1.0)/nqc1, (iq2 - 1.0)/nqc2, (iq3 - 1.0)/nqc3/)
          Rp(iq, :) = (/iq1-1.0, iq2-1.0, iq3-1.0/)
          ! WRITE(stdout, '(a, 6F20.10)') "Check grids: ", xqc_cryst(iq, :), Rp(iq, :) 
          iq = iq + 1
        ENDDO
      ENDDO
    ENDDO
    !
    CALL fkbounds(nktotf_f * nbndsub, lower_bnd, upper_bnd)
    !
    DO iq = 1, nktotf
      CALL iq_to_iqf(iq, nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn, iq_f)
      !
      A_cvq = 0.d0
      IF(my_pool_id == ionode_id) THEN
        ! WRITE(stdout, '(a, 2I10)') 'Checking iq and iq_f: ', iq, iq_f 
        A_cvq(:, :, :, :) = A_cv_all(:,:,:,:, iq_f)
      ENDIF
      CALL mp_bcast(A_cvq, ionode_id, inter_pool_comm)
      !
      ! A_mp_vk
      !
      DO iv = 1, nbndv_explrn 
        DO ik = 1, nktotf_f
          CALL indexing_q2(ik, iq_f, nqc1, nqc2, nqc3, 1, ikplusq)
          expTable(1:3) = EXP(twopi * ci * xqc_cryst(ikplusq, 1:3))
          !
          DO icomp = lower_bnd, upper_bnd
            iRp = (icomp - 1)/ nbndsub + 1
            im = MOD(icomp - 1, nbndsub) + 1
            DO ic = nbndv_explrn + 1, nbndsub
              DO iex = 1, negnv_explrn
                A_mp_vk(im, iRp, iv, ik) = A_mp_vk(im, iRp, iv, ik) + &
                Asq(iex, iq) * A_cvq(iv, ic-nbndv_explrn, ik, iex) * &
                CONJG(cu_big(ic, im, ikplusq)) * & ! U^\dagger_mck+q 
                PRODUCT(expTable ** Rp(iRp,1:3))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      ! A_mp_ck
      !
      DO ic = 1, nbndc_explrn 
        DO ik = 1, nktotf_f
          CALL indexing_q2(iq_f, ik, nqc1, nqc2, nqc3, -1, ikminusq)
          expTable(1:3) = EXP(-twopi * ci * xqc_cryst(ikminusq, 1:3))
          DO icomp = lower_bnd, upper_bnd
            iRp = (icomp - 1)/ nbndsub + 1
            im = MOD(icomp - 1, nbndsub) + 1
            DO iv = 1, nbndv_explrn
              DO iex = 1, negnv_explrn
                A_mp_ck(im, iRp, ic, ik) = A_mp_ck(im, iRp, ic, ik) + &
                Asq(iex, iq) * A_cvq(iv, ic, ikminusq, iex) * &
                cu_big(iv, im, ikminusq) * & ! U_vmk-q 
                PRODUCT(expTable ** Rp(iRp,1:3))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum(A_mp_vk, inter_pool_comm)
    CALL mp_sum(A_mp_ck, inter_pool_comm)
    !
    WRITE(stdout, '(5x, a)') 'Saving the A_mp_vk and A_mp_ck coefficients.'
    ! Save the A_mp_vk and A_mp_ck coeff
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE='A_mp_vk.plrn', FORM='formatted')
      WRITE(iuexphg, '(a)') 'im, iRp, iv, ik, coeffs'
      !
      DO im = 1, nbndsub
        DO iRp = 1, nktotf_f
          DO iv = 1, nbndv_explrn
            DO ik = 1, nktotf_f
              WRITE(iuexphg, '(4I10, 2E20.10)') im, iRp, iv, ik, A_mp_vk(im, iRp, iv, ik)
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
      !
      CLOSE(iuexphg)
      !
      OPEN(UNIT=iuexphg, FILE='A_mp_ck.plrn', FORM='formatted')
      WRITE(iuexphg, '(a)') 'im, iRp, ic, ik, coeffs'
      !
      DO im = 1, nbndsub
        DO iRp = 1, nktotf_f
          DO ic = 1, nbndc_explrn
            DO ik = 1, nktotf_f
              WRITE(iuexphg, '(4I10, 2E20.10)') im, iRp, ic, ik, A_mp_ck(im, iRp, ic, ik)
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
      !
      CLOSE(iuexphg)
    ENDIF
    !
    ! End of calculating the A_mp^ck and A_mp^vk coefficients
    !
    DEALLOCATE(A_cvq)
    DEALLOCATE(xqc_cryst)
    DEALLOCATE(Rp)
    DEALLOCATE(cu_big)
    IF(my_pool_id == ionode_id) DEALLOCATE(A_cv_all)
    !
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout, '(5x, a)') 'Calculating the charge densities.'
    ! Compute the charge densities 
    CALL write_explrn_wavefunction(A_mp_vk, A_mp_ck)
    ! Deallocate all the arrays 
    DEALLOCATE(Asq)
  !-----------------------------------------------------------------------
  END SUBROUTINE plot_explrn_densities
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE write_explrn_wavefunction(A_mp_vk, A_mp_ck)
  !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    !! Write real space wavefunction to file.
    !! Adapted from write_real_space_wavefunction from polaron.f90
    !-----------------------------------------------------------------------------------
    USE ep_constants, ONLY : zero, czero, cone, twopi, ci, bohr2ang
    USE input,        ONLY : nbndsub, step_wf_grid_plrn, nbndc_explrn, nbndv_explrn, &
                              nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn, &
                              plot_explrn_e, plot_explrn_h
    USE global_var,    ONLY : nktotf
    USE modes,         ONLY : nmodes
    USE io_var,        ONLY : iun_plot, iuexphg
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id, ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(INOUT) :: A_mp_vk(:,:,:,:)
    !! A_mp_vk matrix for plotting explrn wfcs
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(INOUT) :: A_mp_ck(:,:,:,:)
    !! A_mp_ck matrix for plotting explrn wfcs
    COMPLEX(KIND = DP), ALLOCATABLE                :: A_mp_vk_sub(:,:,:) ! 
    !! For memory saving (nbndsub, nktotf, upper_bnd-lower_bnd+1)
    COMPLEX(KIND = DP), ALLOCATABLE                :: A_mp_ck_sub(:,:,:) 
    !! For memory saving (nbndsub, nktotf, upper_bnd-lower_bnd+1)
    CHARACTER(LEN = 60)                            :: e_plrn_file
    !! File name of the electron density
    CHARACTER(LEN = 60)                            :: h_plrn_file
    !! File name of the hole density
    CHARACTER(LEN = 60)                            :: filename
    !! File name of dtau.plrn
    INTEGER                                        :: ierr
    !! Error status
    INTEGER                                        :: nqf_p(3)
    !! Dimension of the q grid
    INTEGER                                        :: ibnd
    !! Index for Wannier orbitals
    INTEGER                                        :: itemp
    !! Dummy variable
    INTEGER                                        :: indexkn1
    !! Combined band and k-point index
    INTEGER                                        :: nxx, nyy, nzz
    !! Number of grid points in real space supercell along cartesian directions
    INTEGER                                        :: ig_vec(1:3)
    !! Supercell latice vector coordinates
    INTEGER                                        :: iRp
    !! Lattice vector counter
    INTEGER                                        :: n_grid(3)
    !! Number of grid points in cell where Wannier functions are written
    INTEGER                                        :: grid_start(3)
    !! Initial grid point within this pool
    INTEGER                                        :: grid_end(3)
    !! Final grid point within this pool
    INTEGER                                        :: n_grid_super(3)
    !! Number of grid points in supercellcell where polaron wave function is written
    INTEGER                                        :: rpc(1:3)
    !! Grid point coordinates for lattice vectors
    INTEGER                                        :: species(50)
    !! Atomic species counter
    INTEGER                                        :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER                                        :: shift(1:3)
    !! Shift vector coordinates
    INTEGER                                        :: ishift
    !! Shift counter
    REAL(KIND = DP)                                :: orig(3)
    !! Supercell origin coordinates
    REAL(KIND = DP)                                :: cell(3, 3)
    !! Supercell lattice vectors
    COMPLEX(KIND = DP), ALLOCATABLE                :: wann_func(:, :, :, :)
    !! Wannier function in real space
    COMPLEX(KIND = DP)                             :: ctemp(1:3)
    !! Prefactor for polaron center calculation
    COMPLEX(KIND = DP)                             :: b_vec(1:3)
    !! Prefactor for polaron center calculation
    COMPLEX(KIND = DP), ALLOCATABLE                :: cvec(:)
    !! Polaron wave function at each grid point
    REAL(KIND = DP),    ALLOCATABLE                :: cvec2(:)
    !! Polaron wave function magnitude squared at each grid point
    REAL(DP)                                       :: disp_real, disp_imag
    !! Read and imaginary part of the atomic displacement
    INTEGER                                        :: ik, iq
    !! Counter for k and q
    INTEGER                                        :: imode
    !! Counter for phonon mode
    INTEGER                                        :: iq1, iq2, iq3
    !! Index of qx, qy, qz
    INTEGER                                        :: ic, iv
    !! Index of conduction band and valence band
    INTEGER                                        :: lower_bnd, upper_bnd
    !! Lower and upper bound for parallelization
    INTEGER                                        :: icomp
    !! Combined index for im and iRp
    REAL(DP), ALLOCATABLE                          :: Rp(:,:)
    !! List of all unit cells in the supercell
    REAL(DP)                                       :: progress(11)
    !! Array to store the progress of the calculation
    REAL(DP)                                       :: counting, counting_tot
    !! Counters for computing current progress
    INTEGER                                        :: current_pr
    !! Counters for computing current progress
    INTEGER                                        :: nktotf_f
    !! Number of grid points in the uniform grid
    !
    ! Only use a block of the code to save some memory
    nktotf_f = nqc1 * nqc2 * nqc3
    !
    CALL fkbounds(nktotf_f * nbndv_explrn, lower_bnd, upper_bnd)
    ALLOCATE(A_mp_vk_sub(nbndsub, nktotf_f, upper_bnd-lower_bnd+1))
    DO icomp = lower_bnd, upper_bnd 
      ik = (icomp - 1)/ nbndv_explrn + 1
      iv = MOD(icomp - 1, nbndv_explrn) + 1
      A_mp_vk_sub(:,:,icomp-lower_bnd+1) = A_mp_vk(:,:,iv, ik)
    ENDDO
    DEALLOCATE(A_mp_vk)
    !
    CALL fkbounds(nktotf_f * nbndc_explrn, lower_bnd, upper_bnd)
    ALLOCATE(A_mp_ck_sub(nbndsub, nktotf_f, upper_bnd-lower_bnd+1))
    DO icomp = lower_bnd, upper_bnd 
      ik = (icomp - 1)/ nbndc_explrn + 1
      ic = MOD(icomp - 1, nbndc_explrn) + 1
      A_mp_ck_sub(:,:,icomp-lower_bnd+1) = A_mp_ck(:,:,ic, ik)
    ENDDO
    DEALLOCATE(A_mp_ck)
    !
    nqf_p(1:3) = (/nqc1, nqc2, nqc3/)
    !
    ALLOCATE(Rp(nktotf_f, 3))
    !
    iq = 1
    DO iq1=1,nqc1
      DO iq2=1,nqc2
        DO iq3=1,nqc3
          Rp(iq, :) = (/iq1-1.0, iq2-1.0, iq3-1.0/)
          ! WRITE(stdout, '(a, 6F20.10)') "Check grids: ", xqc_cryst(iq, :), Rp(iq, :) 
          iq = iq + 1
        ENDDO
      ENDDO
    ENDDO
    !
    ! read dtau.plrn, get the displacement.
    !
    ALLOCATE(dtau(nmodes, nktotf))
    dtau = 0.d0
    !
    filename = 'dtau.plrn'
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
      READ(iuexphg, '(5I10)') itemp, itemp, itemp, itemp, itemp
      !
      DO iq=1, nktotf
        DO imode=1, nmodes
          READ(iuexphg, *) itemp, itemp, &
          disp_real, disp_imag
          dtau(imode, iq) = disp_real + ci * disp_imag
        ENDDO
      ENDDO     
      CLOSE(iuexphg)
    ENDIF
    CALL mp_bcast(dtau, ionode_id, inter_pool_comm)
    !
    ! read cube files for the real-space Wannier function Wm(r)
    !
    CALL read_explrn_wannier_cube(wann_func, species, &
       n_grid, grid_start, grid_end)
    !
    cell(1:3, 1) = at(1:3, 1) * nqf_p(1) * alat
    cell(1:3, 2) = at(1:3, 2) * nqf_p(2) * alat
    cell(1:3, 3) = at(1:3, 3) * nqf_p(3) * alat
    !
    orig(1:3) = zero
    n_grid_super(1:3) = nqf_p(1:3) * n_grid(1:3)
    !
    progress = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0/)
    !
    b_vec(1:3) = twopi * ci / REAL(n_grid_super(1:3))
    !
    ALLOCATE(cvec(1:n_grid_super(1)), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating cvec', 1)
    ALLOCATE(cvec2(1:n_grid_super(1)), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating cvec2', 1)
    !
    IF(plot_explrn_e) THEN
      e_plrn_file = 'e_psir_plrn.xsf'
      ! Write the file head including information of structures,
      ! using the same format of
      IF (ionode) THEN
        IF(step_k1_explrn /= 1 .OR. step_k2_explrn /= 1 .OR. step_k3_explrn /= 1) THEN
          CALL write_explrn_dtau_xsf(e_plrn_file, nqc1, nqc2, nqc3, .FALSE.)
        ELSE
          CALL write_explrn_dtau_xsf(e_plrn_file, nqc1, nqc2, nqc3, .TRUE.)
        ENDIF
      ENDIF
      !
      ! Write e_psir_plrn.xsf
      IF (ionode) THEN
        OPEN(UNIT = iuexphg, FILE = TRIM(e_plrn_file), POSITION='APPEND')
        WRITE(iuexphg, '(/)')
        WRITE(iuexphg, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
        WRITE(iuexphg, '(3i6)')  n_grid_super / step_wf_grid_plrn
        WRITE(iuexphg, '(3f12.6)') zero, zero, zero
        WRITE(iuexphg, '(3f12.7)') cell(1:3, 1) * bohr2ang
        WRITE(iuexphg, '(3f12.7)') cell(1:3, 2) * bohr2ang
        WRITE(iuexphg, '(3f12.7)') cell(1:3, 3) * bohr2ang
      ENDIF
      !
      ! Parallelize over kpts and bands
      CALL fkbounds(nktotf_f * nbndv_explrn, lower_bnd, upper_bnd)
      !
      ctemp = czero
      current_pr = 1
      DO nzz = 1, n_grid_super(3), step_wf_grid_plrn
        DO nyy = 1, n_grid_super(2), step_wf_grid_plrn
          ! Check progress 
          counting =  REAL(nzz, KIND=DP) * nyy
          counting_tot = n_grid_super(2) * n_grid_super(3)
          !
          IF ( counting / counting_tot > progress(current_pr)) THEN
            WRITE(stdout, '(5x, a, 1F5.2, a)') 'Current progress: ', &
            counting / counting_tot * 100, '%'
            current_pr = current_pr + 1
          ENDIF
          !
          cvec2 = 0.d0
          DO icomp = lower_bnd, upper_bnd 
            ik = (icomp - 1)/ nbndv_explrn + 1
            iv = MOD(icomp - 1, nbndv_explrn) + 1
            !
            cvec = czero
            DO iRp = 1, nktotf_f !-nqf_p(3)/2, (nqf_p(3)+1)/2  
              Rp_vec(1:3) = Rp(iRp, 1:3)
              rpc(1:3) = Rp_vec(1:3) * n_grid(1:3) !- (n_grid_super(1:3))/2
              !
              DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
                ! To make sure that all the nonzero points are included,
                ! we need to try from -1 to 1 neighbor supercells
                DO ishift = 1, 27
                  shift(1:3) = index_explrn_shift(ishift)
                  ig_vec(1:3) = (/nxx, nyy, nzz/) - rpc(1:3) + shift(1:3) * n_grid_super(1:3)
                  IF (ALL(ig_vec(1:3) <= grid_end(1:3)) .AND. &
                    ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                    DO ibnd = 1, nbndsub !TODO change to nbndsub
                      indexkn1 = (iRp - 1) * nbndsub + ibnd
                      !
                      cvec(nxx) = cvec(nxx)  +  A_mp_vk_sub(ibnd, iRp, icomp-lower_bnd+1) &
                                  * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                    ENDDO !ibnd
                  ENDIF
                ENDDO ! iscx
              ENDDO ! ipx
            ENDDO ! nxx
            !
            DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
              cvec2(nxx) = cvec2(nxx) + ABS(cvec(nxx)) ** 2 
            ENDDO
            !
          ENDDO !ik, iv
          !
          CALL mp_sum(cvec2, inter_pool_comm)
          !
          IF(ionode) THEN
            WRITE (iuexphg, '(5e13.5)', ADVANCE='yes') ABS(cvec2(::step_wf_grid_plrn))
          ENDIF
        ENDDO
      ENDDO
      !
      IF (ionode) THEN
        WRITE (iuexphg, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
        CLOSE(iuexphg)
        WRITE(stdout, "(5x, '|\Psi(r)|^2 of electron written to file.')")
      ENDIF
      !
    ENDIF
    !
    ! Write h_psir_plrn.xsf
    IF(plot_explrn_h) THEN
      h_plrn_file = 'h_psir_plrn.xsf'
      ! Write the file head including information of structures,
      ! using the same format of
      IF (ionode) THEN
        IF(step_k1_explrn /= 1 .OR. step_k2_explrn /= 1 .OR. step_k3_explrn /= 1) THEN
          CALL write_explrn_dtau_xsf(h_plrn_file, nqc1, nqc2, nqc3, .FALSE.)
        ELSE
          CALL write_explrn_dtau_xsf(h_plrn_file, nqc1, nqc2, nqc3, .TRUE.)
        ENDIF
      ENDIF
      !
      IF (ionode) THEN
        OPEN(UNIT = iuexphg, FILE = TRIM(h_plrn_file), POSITION='APPEND')
        WRITE(iuexphg, '(/)')
        WRITE(iuexphg, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
        WRITE(iuexphg, '(3i6)')  n_grid_super / step_wf_grid_plrn
        WRITE(iuexphg, '(3f12.6)') zero, zero, zero
        WRITE(iuexphg, '(3f12.7)') cell(1:3, 1) * bohr2ang
        WRITE(iuexphg, '(3f12.7)') cell(1:3, 2) * bohr2ang
        WRITE(iuexphg, '(3f12.7)') cell(1:3, 3) * bohr2ang
      ENDIF
      !
      CALL fkbounds(nktotf_f * nbndc_explrn, lower_bnd, upper_bnd)
      !
      ctemp = czero
      current_pr = 1
      DO nzz = 1, n_grid_super(3), step_wf_grid_plrn
        DO nyy = 1, n_grid_super(2), step_wf_grid_plrn
          counting =  REAL(nzz, KIND=DP) * nyy
          counting_tot = n_grid_super(2) * n_grid_super(3)
          !
          IF ( counting / counting_tot > progress(current_pr)) THEN
            WRITE(stdout, '(5x, a, 1F5.2, a)') 'Current progress: ', &
            counting / counting_tot * 100, '%'
            current_pr = current_pr + 1
          ENDIF
          !
          cvec2 = 0.d0
          DO icomp = lower_bnd, upper_bnd 
            ik = (icomp - 1)/ nbndc_explrn + 1
            ic = MOD(icomp - 1, nbndc_explrn) + 1
            !  
            cvec = czero
            DO iRp = 1, nktotf_f          
              Rp_vec(1:3) = Rp(iRp, 1:3)
              rpc(1:3) = Rp_vec(1:3) * n_grid(1:3) !- (n_grid_super(1:3))/2
              !
              DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
                !
                ! To make sure that all the nonzero points are included,
                ! we need to try from -1 to 1 neighbor supercells
                DO ishift = 1, 27
                  shift(1:3) = index_explrn_shift(ishift)
                  ig_vec(1:3) = (/nxx, nyy, nzz/) - rpc(1:3) + shift(1:3) * n_grid_super(1:3)
                  IF (ALL(ig_vec(1:3) <= grid_end(1:3)) .AND. &
                    ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                    DO ibnd = 1, nbndsub !TODO change to nbndsub
                      indexkn1 = (iRp - 1) * nbndsub + ibnd
                      !TODO eigvec_wan(indexkn1, 1) should be eigvec_wan(indexkn1, iplrn)
                      ! cvec(nxx) = cvec(nxx)  +  A_mp_ck(ibnd, iRp, ic, ik) * CONJG(wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd))
                      cvec(nxx) = cvec(nxx)  +  A_mp_ck_sub(ibnd, iRp, icomp-lower_bnd+1) &
                                  * CONJG(wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd))
                    ENDDO !ibnd
                  ENDIF
                ENDDO ! iscx
              ENDDO ! ipx
            ENDDO ! nxx
            !
            DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
              cvec2(nxx) = cvec2(nxx) + ABS(cvec(nxx)) ** 2 
            ENDDO
            !
          ENDDO !ik, ic
          !
          CALL mp_sum(cvec2, inter_pool_comm)
          !
          IF(ionode) THEN
            WRITE (iuexphg, '(5e13.5)', ADVANCE='yes') ABS(cvec2(::step_wf_grid_plrn))
          ENDIF
        ENDDO
      ENDDO
      !
      IF (ionode) THEN
        WRITE (iuexphg, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
        CLOSE(iuexphg)
        WRITE(stdout, "(5x, '|\Psi(r)|^2 of hole written to file.')")
      ENDIF
    ENDIF
    !
    DEALLOCATE(wann_func, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(cvec , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(cvec2 , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    !
    DEALLOCATE(dtau)
    DEALLOCATE(A_mp_ck_sub)
    DEALLOCATE(A_mp_vk_sub)
    DEALLOCATE(Rp)
  !-----------------------------------------------------------------------
  END SUBROUTINE write_explrn_wavefunction
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE save_solution(solutions, numk, nummode, filename)
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE ep_constants,  ONLY : ryd2ev, eps6
    USE input,         ONLY : nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool,my_pool_id
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)            :: numk
    !! Number of BZ grid points
    INTEGER, INTENT(in)            :: nummode
    !! Number of modes 
    COMPLEX(DP), INTENT(in)        :: solutions(nummode, numk)
    !! Solutin arrays
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! File name    
    INTEGER                        :: iexq
    !! Index of grid point
    INTEGER                        :: imode
    !! Index of mode
    !
    IF(my_pool_id .EQ. ionode_id) THEN
      OPEN(UNIT=iuexphg, FILE=filename, FORM='formatted')
      WRITE(iuexphg, '(5I10)') nqc1 / step_k1_explrn, nqc2 / step_k2_explrn, &
                               nqc3 / step_k3_explrn, numk, nummode
      !
      DO iexq=1,numk
        DO imode=1,nummode
          WRITE(iuexphg, '(2I10, 2E20.10)') iexq, imode, solutions(imode, iexq)
        ENDDO
      ENDDO 
      !
      CLOSE(iuexphg)
    ENDIF
  !-----------------------------------------------------------------------  
  END SUBROUTINE save_solution
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE write_explrn_bmat(eigvec_wan, filename, etf_all, nktot)
  !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    !! Write Bmat and phonon frequency to filename
    !!
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev
    USE global_var,    ONLY : nkf, nqtotf
    USE input,         ONLY : nstate_plrn, nqf1, nqf2, nqf3, scell_mat_plrn, &
                              nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    USE ions_base,     ONLY : amass, ityp
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in)        :: filename
    !! File name
    COMPLEX(KIND = DP), INTENT(in)        :: eigvec_wan(:, :)
    !! Bmat or dtau to be written
    REAL(KIND = DP), INTENT(in), OPTIONAL :: etf_all(:, :)
    !! Phonon frequencies
    INTEGER, INTENT(in), OPTIONAL         :: nktot
    !! Number of q points
    INTEGER                               :: ik
    !! Index of BZ grid point
    INTEGER                               :: ibnd
    !! Index of band
    INTEGER                               :: nqtot
    !! Actual total number of q points
    !
    OPEN(UNIT = iuexphg, FILE = TRIM(filename))
    !
    IF(PRESENT(nktot)) THEN
      WRITE(iuexphg, '(5I10)') nqc1/step_k1_explrn, nqc2/step_k2_explrn, &
                                     nqc3/step_k3_explrn, nktot, nmodes
      nqtot = nktot
    ELSE
      WRITE(iuexphg, '(5I10)') nqc1, nqc2, &
                                     nqc3, nqtotf, nmodes
      nqtot = nqtotf
    ENDIF
    !
    DO ik = 1, nqtot
      DO ibnd = 1, nmodes ! p
        IF (PRESENT(etf_all)) THEN ! \kappa, \alpha
          !!JLB: Changed format for improved accuracy
          WRITE(iuexphg, '(2I5, 4ES18.10)') ik, ibnd, etf_all(ibnd, ik) * ryd2mev, &
             eigvec_wan(ibnd, ik), ABS(eigvec_wan(ibnd, ik))
        ELSE
          !JLB: Changed format for improved accuracy
          WRITE(iuexphg, '(2ES18.10)') eigvec_wan(ibnd, ik)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(iuexphg)
  !-----------------------------------------------------------------------
  END SUBROUTINE write_explrn_bmat
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE normalize_Asq(Asq)
  !-----------------------------------------------------------------------  
    USE modes,         ONLY : nmodes
    USE input,         ONLY : negnv_explrn
    USE global_var,    ONLY : nktotf
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(inout)  :: Asq(:,:)
    !! Array to be normalized
    INTEGER                     :: iq, imode
    !! Counter for q and mode
    INTEGER                     :: ierr
    !! Error status
    REAL(DP)                    :: Asqnorm 
    !! Norm of Asq
    !
    Asqnorm = 0.0
    DO iq=1,nktotf
      DO imode=1,negnv_explrn
        Asqnorm = Asqnorm + ABS(Asq(imode, iq)) ** 2
      ENDDO
    ENDDO
    !
    Asq(:,:) = Asq(:,:) / (Asqnorm**0.5) * (nktotf**0.5)
    !
  !-----------------------------------------------------------------------
  END SUBROUTINE normalize_Asq
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE iq_to_iqf(iq, nqc1, nqc2, nqc3, step_k1_explrn, step_k2_explrn, step_k3_explrn, iq_f)
  !-----------------------------------------------------------------------
    !  
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: iq
    !! Index of q in the coarser coarse grid
    INTEGER, INTENT(in)  :: nqc1, nqc2, nqc3
    !! Dimensions of the full coarse q grid 
    INTEGER, INTENT(in)  :: step_k1_explrn, step_k2_explrn, step_k3_explrn
    !! Skipping steps
    INTEGER, INTENT(out) :: iq_f
    !! Index of q in the full coarse grid
    INTEGER              :: iq1, iq2, iq3
    !! Index of qx, qy, qz
    INTEGER              :: iq1_f, iq2_f, iq3_f
    !! Index of qx, qy, qz in the full coarse grid
    !
    ! Find the index of iq in the finer coarse grid
    iq3 = mod(iq - 1, nqc3 / step_k3_explrn)
    iq2 = mod( (iq - 1 - iq3) / (nqc3 / step_k3_explrn), nqc2 / step_k2_explrn )
    iq1 = (iq - 1) / (nqc3 * nqc2 / step_k3_explrn / step_k2_explrn)
    !
    iq1_f = iq1 * step_k1_explrn 
    iq2_f = iq2 * step_k2_explrn 
    iq3_f = iq3 * step_k3_explrn 
    !
    iq_f = iq1_f * nqc2 * nqc3 + iq2_f * nqc3 + iq3_f + 1
  !-----------------------------------------------------------------------  
  END SUBROUTINE iq_to_iqf
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  FUNCTION index_explrn_shift(ishift)
  !-----------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ishift
    !! Count for shift
    INTEGER             :: index_explrn_shift(1:3)
    !! Shift 
    !
    index_explrn_shift(1) = (ishift - 1) / 9 - 1
    index_explrn_shift(2) = MOD(ishift - 1, 9) / 3 - 1
    index_explrn_shift(3) = MOD(ishift - 1, 3) - 1
    !
    IF (ANY(index_explrn_shift < -1) .OR. ANY(index_explrn_shift > 1)) THEN
      CALL errore('index_explrn_shift', 'index_explrn_shift not correct!', 1)
    ENDIF
  !-----------------------------------------------------------------------
  END FUNCTION index_explrn_shift
  !-----------------------------------------------------------------------
!-----------------------------------------------------------------------
END MODULE
!-----------------------------------------------------------------------
