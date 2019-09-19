  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !                                                                            
  !----------------------------------------------------------------------
  MODULE io_epw
  !----------------------------------------------------------------------
  !! 
  !! This module contains various writing or reading routines to files. 
  !! Most of them are for restart purposes. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE Fin_write(iter, F_in, av_mob_old, elec)
    !----------------------------------------------------------------------------
    !!
    !! Writes the F without magnetic field for restart
    !! 
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilFi_all
    USE io_files,      ONLY : diropn
    USE epwcom,        ONLY : nstemp
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE elph2,         ONLY : ibndmax, ibndmin, nkqtotf, nktotf, nbndfst
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iter
    !! Iteration number
    REAL(KIND = DP), INTENT(in) :: F_in(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i  
    REAL(KIND = DP), INTENT(in) :: av_mob_old(nstemp)
    !! Error in the hole mobility
    LOGICAL, INTENT(in) :: elec
    !! IF true we do electron mobility, if false the hole one. 
    ! 
    ! Local variable
    LOGICAL :: exst
    !! File exist
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
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
    exst = .FALSE.
    aux(:) = zero
    IF (mpime == ionode_id) THEN
      !
      lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
      ! First element is the iteration number
      aux(1) = iter
      ! 
      i = 1
      DO itemp = 1, nstemp
        i = i + 1  
        ! Value of the previous h mobility (used for error evaluation)
        aux(i) = av_mob_old(itemp)
      ENDDO
      ! 
      i = 1 + nstemp 
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO idir = 1, 3
              i = i +1
              aux(i) = F_in(idir, ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! 
      ! Electron mobility
      IF (elec) THEN
        CALL diropn(iufilFi_all, 'Fin_restartcb', lfi_all, exst)
        CALL davcio(aux, lfi_all, iufilFi_all, 1, +1)
      ELSE 
        CALL diropn(iufilFi_all, 'Fin_restart', lfi_all, exst)
        CALL davcio( aux, lfi_all, iufilFi_all, 1, +1)
      ENDIF
      CLOSE(iufilFi_all)
      !
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE Fin_write
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE Fin_read(iter, F_in, av_mob_old, elec)
    !----------------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilFi_all
    USE epwcom,    ONLY : nstemp
    USE constants_epw, ONLY : zero
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, nbndfst, nktotf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iter
    !! Iteration number
    REAL(KIND = DP), INTENT(inout) :: F_in(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i  
    REAL(KIND = DP), INTENT(inout) :: av_mob_old(nstemp)
    !! Error in the hole mobility
    LOGICAL, INTENT(in) :: elec
    !! IF true we do electron mobility, if false the hole one. 
    !
    ! Local variable
    CHARACTER(LEN = 256) :: name1
    !! Variable name
    LOGICAL :: exst
    !! Does the file exist
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
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      ! 
      ! First inquire if the file exists
      IF (elec) THEN
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restartcb1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restartcb'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        ! 
        IF (exst) THEN ! read the file
          !
          lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
          CALL diropn(iufilFi_all, 'Fin_restartcb', lfi_all, exst)
          CALL davcio(aux, lfi_all, iufilFi_all, 1, -1)
          !
          ! First element is the iteration number
          iter = INT(aux(1)) 
          !
          i = 1
          DO itemp = 1, nstemp
            i = i + 1  
            ! Last value of hole mobility 
            av_mob_old(itemp) = aux(i)
          ENDDO
          ! 
          i = 1 + nstemp 
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                DO idir = 1, 3
                  i = i + 1
                  F_in(idir, ibnd, ik, itemp) = aux(i)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilFi_all)
        ENDIF
      ELSE ! hole
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restart1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restart'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        ! 
        IF (exst) THEN ! read the file
          !
          lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
          CALL diropn(iufilFi_all, 'Fin_restart', lfi_all, exst)
          CALL davcio(aux, lfi_all, iufilFi_all, 1, -1)
          !
          ! First element is the iteration number
          iter = INT(aux(1))
          !
          i = 1
          DO itemp = 1, nstemp
            i = i + 1
            ! Last value of hole mobility 
            av_mob_old(itemp) = aux(i)
          ENDDO
          ! 
          i = 1 + nstemp
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, (nbndfst)
                DO idir = 1, 3
                  i = i + 1
                  F_in(idir, ibnd, ik, itemp) = aux(i)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilFi_all)
        ENDIF
      ENDIF
    ENDIF ! mpime
    ! 
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iter,       ionode_id, world_comm)
      CALL mp_bcast(F_in,       ionode_id, world_comm)
      CALL mp_bcast(av_mob_old, ionode_id, world_comm)
      WRITE(stdout, '(a,i10)' ) '     Restart from iter: ', iter
    ENDIF ! exists
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE Fin_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_merge_parallel()
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_var,           ONLY : iunepmat_merge, iunepmat, iunepmatcb_merge,              &
                                 iunepmatcb, iunsparseq_merge, iunsparsek_merge,          &
                                 iunsparsej_merge,iunsparset_merge, iunepmatcb_merge,     &
                                 iunsparseqcb_merge, iunsparsekcb_merge, iunsparsei_merge,&
                                 iunsparseicb_merge, iunsparsejcb_merge, iunsparsetcb_merge
    USE mp_global,        ONLY : my_pool_id, npool, world_comm 
    USE io_files,         ONLY : tmp_dir, prefix 
    USE mp,               ONLY : mp_sum, mp_barrier
    USE io_global,        ONLY : stdout
    USE elph2,            ONLY : lrepmatw2_merge, lrepmatw5_merge
    USE epwcom,           ONLY : int_mob, carrier, ncarrier
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE,MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION,          &
                                 MPI_STATUS_IGNORE, MPI_INTEGER
#endif
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! Pool number, character
    CHARACTER(LEN = 256) :: dirname(2)
    !! Name of the directory to hold files
    CHARACTER(LEN = 256) :: filename(6)
    !! Name of the files to merge files
    CHARACTER(LEN = 256) :: path_to_files(2)
    !! Name of the path to files
    INTEGER :: i2, i4, i5, i6
    !! Indexes to loop over file sizes
    INTEGER :: ipool
    !! Process index
    INTEGER :: lrepmatw2_tot(npool)
    !! Lenght of each file
    INTEGER :: lrepmatw5_tot(npool)
    !! Lenght of each file
    INTEGER :: ich
    !! Loop over directories
    INTEGER :: ifil
    !! Index over the files
    INTEGER :: io_u(6)
    !! Input output units
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: sparse(:, :)
    !! Vaariable for reading and writing the files
    INTEGER, ALLOCATABLE :: sparsecb(:, :)
    !! Vaariable for reading and writing the files
    REAL(KIND = DP), ALLOCATABLE :: trans_prob(:)
    !! Variable for reading and writing trans_prob
    REAL(KIND = DP), ALLOCATABLE :: trans_probcb(:)
    !! Variable for reading and writing trans_prob
#if defined(__MPI)
    INTEGER (KIND = MPI_OFFSET_KIND) :: lsize
    !! Size of what we write
    INTEGER (KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset while writing scattering to files
    !
    IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 0.0))) THEN
      !
      ALLOCATE(trans_prob(lrepmatw2_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating trans_prob', 1)
      ALLOCATE(sparse(5, lrepmatw2_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating sparse', 1)
      !
      io_u(1) = iunepmat_merge
      io_u(2) = iunsparseq_merge
      io_u(3) = iunsparsek_merge
      io_u(4) = iunsparsei_merge
      io_u(5) = iunsparsej_merge
      io_u(6) = iunsparset_merge
      !
      dirname(1) = 'Fepmatkq1'
      dirname(2) = 'Fsparse'
      !
      filename(1) = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkq1' 
      filename(2) = 'sparseq'
      filename(3) = 'sparsek'
      filename(4) = 'sparsei'
      filename(5) = 'sparsej'
      filename(6) = 'sparset'
      !
      path_to_files(1) = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkq1' // '_'
      path_to_files(2) = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_'
      !
      lrepmatw2_tot = 0
      lrepmatw2_tot(my_pool_id + 1) = lrepmatw2_merge
      CALL mp_sum(lrepmatw2_tot, world_comm)
      DO ich = 1, 6
        CALL mp_barrier(world_comm)
        CALL MPI_FILE_OPEN(world_comm, filename(ich), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,io_u(ich), ierr)
      ENDDO
      !
      DO ich = 1, 2 
        ! Read files per processor
        WRITE(my_pool_id_ch, "(I0)") my_pool_id
        filint = TRIM(path_to_files(ich)) // TRIM(my_pool_id_ch)
        OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'old', FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
        IF (ich == 1) THEN
          DO i2 = 1, lrepmatw2_merge
            READ(iunepmat) trans_prob(i2)
          ENDDO
        ELSE
          DO i2 = 1, lrepmatw2_merge
            DO ifil = 1, 5
              READ(iunepmat) sparse(ifil, i2)
            ENDDO
          ENDDO
        ENDIF
        CLOSE(iunepmat, STATUS = 'delete')
        IF (ich == 1) THEN
          lrepmatw = INT(SUM(lrepmatw2_tot(1:my_pool_id + 1)) - lrepmatw2_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
          & 8_MPI_OFFSET_KIND 
          lsize = INT(lrepmatw2_merge, KIND = MPI_OFFSET_KIND) 
          CALL MPI_FILE_WRITE_AT(io_u(1), lrepmatw, trans_prob(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        ELSE
          DO ifil = 1, 5
            lrepmatw = INT(SUM(lrepmatw2_tot(1:my_pool_id + 1)) - lrepmatw2_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
            & 4_MPI_OFFSET_KIND 
            lsize = INT(lrepmatw2_merge, KIND = MPI_OFFSET_KIND) 
            CALL MPI_FILE_WRITE_AT(io_u(ifil+1), lrepmatw, sparse(ifil,:), lsize, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          ENDDO
        ENDIF
      ENDDO
      !
      DO ich = 1, 6
        CALL MPI_FILE_CLOSE(io_u(ich), ierr)
      ENDDO
      !
      DEALLOCATE(trans_prob, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating trans_prob', 1)
      DEALLOCATE(sparse, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating sparse', 1)
      !
    ENDIF
    IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 0.0))) THEN
      !
      ALLOCATE(trans_probcb(lrepmatw5_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating trans_probcb', 1)
      ALLOCATE(sparsecb(5, lrepmatw5_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error allocating sparsecb', 1)
      !
      io_u(1) = iunepmatcb_merge
      io_u(2) = iunsparseqcb_merge
      io_u(3) = iunsparsekcb_merge
      io_u(4) = iunsparseicb_merge
      io_u(5) = iunsparsejcb_merge
      io_u(6) = iunsparsetcb_merge
      !
      dirname(1) = 'Fepmatkqcb1'
      dirname(2) = 'Fsparsecb'
      !
      filename(1) = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkqcb1' 
      filename(2) = 'sparseqcb'
      filename(3) = 'sparsekcb'
      filename(4) = 'sparseicb'
      filename(5) = 'sparsejcb'
      filename(6) = 'sparsetcb'
      !
      path_to_files(1) = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkqcb1' // '_'
      path_to_files(2) = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparsecb' // '_'
      !
      lrepmatw5_tot = 0
      lrepmatw5_tot(my_pool_id + 1) = lrepmatw5_merge
      CALL mp_sum(lrepmatw5_tot, world_comm)
      DO ich = 1, 6
        CALL mp_barrier(world_comm)
        CALL MPI_FILE_OPEN(world_comm, filename(ich), MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, io_u(ich), ierr)
      ENDDO
      !
      DO ich = 1, 2
        ! Read files per processor
        WRITE(my_pool_id_ch, "(I0)") my_pool_id
        filint = TRIM(path_to_files(ich)) // TRIM(my_pool_id_ch)
        OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
        IF (ich == 1) THEN
          DO i2 = 1, lrepmatw5_merge
            READ(iunepmatcb) trans_probcb(i2)
          ENDDO
        ELSE
          DO i2 = 1, lrepmatw5_merge
            DO ifil = 1 ,5
              READ(iunepmatcb) sparsecb(ifil, i2)
            ENDDO
          ENDDO
        ENDIF
        CLOSE(iunepmatcb, STATUS = 'delete')
        IF (ich == 1) THEN
          lrepmatw = INT(SUM(lrepmatw5_tot(1:my_pool_id + 1)) - lrepmatw5_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
          & 8_MPI_OFFSET_KIND 
          lsize = INT(lrepmatw5_merge, KIND = MPI_OFFSET_KIND) 
          CALL MPI_FILE_WRITE_AT(io_u(1), lrepmatw, trans_probcb(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        ELSE
          DO ifil = 1, 5
            lrepmatw = INT(SUM(lrepmatw5_tot(1:my_pool_id + 1)) - lrepmatw5_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * &
            & 4_MPI_OFFSET_KIND 
            lsize = INT(lrepmatw5_merge, KIND = MPI_OFFSET_KIND) 
            CALL MPI_FILE_WRITE_AT(io_u(ifil + 1), lrepmatw, sparsecb(ifil, :), lsize, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          ENDDO
        ENDIF
      ENDDO
      !
      DO ich = 1, 6
        CALL MPI_FILE_CLOSE(io_u(ich), ierr)
      ENDDO
      ! 
      DEALLOCATE(trans_probcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating trans_probcb', 1)
      DEALLOCATE(sparsecb, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge_parallel', 'Error deallocating sparsecb', 1)
      !
    ENDIF ! in all other cases it is still to decide which files to open
    ! 
#endif
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_merge_parallel
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_open(ind_tot, ind_totcb, lrepmatw2_restart, lrepmatw5_restart)
    !----------------------------------------------------------------------------
    ! 
    ! This routine opens all the files needed to save scattering rates for the IBTE.
    ! 
    USE kinds,            ONLY : DP
    USE io_files,         ONLY : tmp_dir, prefix, create_directory, delete_if_present
    USE io_var,           ONLY : iunepmat, iunsparseq, iunsparsek,                 &
                                 iunsparsei, iunsparsej, iunsparset, iunsparseqcb, &
                                 iunsparsekcb, iunsparseicb, iunsparsejcb,         &
                                 iunsparsetcb, iunepmatcb, iunrestart
    USE mp_global,        ONLY : world_comm, my_pool_id, npool
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE elph2,            ONLY : lrepmatw2_merge, lrepmatw5_merge
    USE epwcom,           ONLY : int_mob, carrier, ncarrier
    USE io_global,        ONLY : ionode_id
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND
#endif
    ! 
    IMPLICIT NONE
    !  
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! To restart opening files
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! To restart opening files
#if defined(__MPI)
    INTEGER (KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER (KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of component for the conduction band
#else
    INTEGER, INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER, INTENT(inout) :: ind_totcb
    !! Total number of component for conduction band
#endif  
    ! 
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! my_pool_id in character type
    CHARACTER(LEN = 256) :: dirname(2), dirnamecb(2)
    !! Name of the directory to hold files
    LOGICAL :: exst
    !! Logical for existence of files
    LOGICAL :: exst2
    !! Logical for existence of files
    INTEGER :: ierr
    !! Error index
    INTEGER :: ilrep
    !! index to loop over the reading elements
    INTEGER :: dummy_int
    !! Dummy INTEGER for reading
    INTEGER :: ipool
    !! Pool index
    INTEGER (KIND = 8) :: position_byte
    !! Position in the file in byte
    REAL (KIND = DP) :: dummy_real
    !! Dummy variable for reading
    !
    WRITE(my_pool_id_ch, "(I0)") my_pool_id
    !
    dirname(1)   = 'Fepmatkq1'
    dirname(2)   = 'Fsparse'
    dirnamecb(1) = 'Fepmatkqcb1'
    dirnamecb(2) = 'Fsparsecb'
    ! 
    INQUIRE(FILE = 'restart_ibte.fmt', EXIST = exst)
    !
    IF (my_pool_id == ionode_id) THEN
      IF (exst) THEN
        OPEN(UNIT = iunrestart, FILE = 'restart_ibte.fmt', STATUS = 'old')
        READ (iunrestart,*) 
        READ (iunrestart,*) 
        READ (iunrestart,*) 
        READ (iunrestart,*) 
        DO ipool = 1, npool
          READ (iunrestart,*) lrepmatw2_restart(ipool)
        ENDDO
        DO ipool = 1, npool
          READ (iunrestart,*) lrepmatw5_restart(ipool)
        ENDDO
        CLOSE(iunrestart)
      ENDIF
    ENDIF
    CALL mp_bcast(exst, ionode_id, world_comm )
    CALL mp_bcast(lrepmatw2_restart, ionode_id, world_comm )
    CALL mp_bcast(lrepmatw5_restart, ionode_id, world_comm )
    !
    ! The restart_ibte.fmt exist - we try to restart
    IF (exst) THEN
      ! Hole
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 1E5))) THEN
        !
        filint = './' // ADJUSTL(TRIM(dirname(1))) // '/'//TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw2_restart(my_pool_id + 1) - 1) * 8 + 1  
            READ(iunepmat, POS=position_byte) dummy_real 
          ENDIF
        ELSE 
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fepmatkq1 folder', 1) 
        ENDIF
        ! 
        filint = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN 
            position_byte = (5 * lrepmatw2_restart(my_pool_id + 1) - 1) * 4 + 1
            READ(iunsparseq, POS = position_byte) dummy_int
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fsparse folder', 1) 
        ENDIF
        !
      ENDIF ! Hole
      ! Electron
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 1E5))) THEN
        !
        filint = './' // ADJUSTL(TRIM(dirnamecb(1))) // '/' // TRIM(prefix) // '.epmatkqcb1'//'_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw5_restart(my_pool_id + 1) - 1) * 8 + 1  
            READ(iunepmatcb, POS = position_byte) dummy_real
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fepmatkqcb1 folder', 1)
        ENDIF
        ! 
        filint = './' // ADJUSTL(TRIM(dirnamecb(2))) // '/' // 'sparsecb' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (5 * lrepmatw5_restart(my_pool_id + 1) - 1) * 4 + 1 
            READ(iunsparseqcb, POS = position_byte) dummy_int
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fsparse folder', 1)
        ENDIF
        !
      ENDIF ! electron
      lrepmatw2_merge = lrepmatw2_restart(my_pool_id + 1)
      lrepmatw5_merge = lrepmatw5_restart(my_pool_id + 1)
      !  
    ELSE ! no restart file present
      ! Hole
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 1E5))) THEN
        ! 
        CALL create_directory(ADJUSTL(TRIM(dirname(1))))
        CALL create_directory(ADJUSTL(TRIM(dirname(2))))
        ! 
        filint = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF 
        ! 
        filint = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF 
        ! 
      ENDIF ! Hole
      ! Electron
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 1E5))) THEN
        ! 
        CALL create_directory(ADJUSTL(TRIM(dirnamecb(1))))
        CALL create_directory(ADJUSTL(TRIM(dirnamecb(2))))        
        ! 
        filint = './' // ADJUSTL(TRIM(dirnamecb(1))) // '/' // TRIM(prefix) // '.epmatkqcb1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !  
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        ! 
        filint = './' // ADJUSTL(TRIM(dirnamecb(2))) // '/' // 'sparsecb' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint)
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        !  
      ENDIF !electron 
      lrepmatw2_merge = 0
      lrepmatw5_merge = 0
      ! 
    ENDIF ! restart 
    !
    ind_tot   = 0
    ind_totcb = 0
    lrepmatw2_restart(:) = 0
    lrepmatw5_restart(:) = 0
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_open
    !----------------------------------------------------------------------------    
    !
    !----------------------------------------------------------------------------
    SUBROUTINE scattering_write(itemp, etemp, ef0, etf_all)
    !----------------------------------------------------------------------------
    !! 
    !! Write scattering rates
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilscatt_rate
    USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, inv_tau_all, nbndfst, nktotf
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
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndsub, nkqtotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file.
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: temp
    !! Temporary file name used to write scattering rate to file. 
    !
    WRITE(stdout, '(/5x,"Writing scattering rate to file"/)')
    !
    IF (mpime == ionode_id) THEN
      !
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1,'(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1,'(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1,'(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate, FILE = name1, FORM = 'formatted')
      WRITE(iufilscatt_rate, '(a)') '# Inverse scattering time (ps)'
      WRITE(iufilscatt_rate, '(a)') '#      ik       ibnd                 E(ibnd)    scattering rate(1/ps)'
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        DO ibnd = 1, nbndfst
          !
          ! note that ekk does not depend on q
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
          !
          WRITE(iufilscatt_rate, '(i9,2x)', ADVANCE = 'no') ik
          WRITE(iufilscatt_rate, '(i9,2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
          WRITE(iufilscatt_rate, '(E22.14)', ADVANCE = 'no') ryd2ev * ekk
          WRITE(iufilscatt_rate, '(E26.16E3)') ryd2mev * meV2invps * inv_tau_all(itemp, ibnd, ik)
          !
        ENDDO
        !
      ENDDO
      !
      CLOSE(iufilscatt_rate)
    ENDIF
    !CALL mp_barrier(inter_pool_comm)
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE scattering_write
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE scattering_read(etemp, ef0, etf_all, inv_tau_all)
    !----------------------------------------------------------------------------
    !!
    !! Read scattering files
    !!  
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilscatt_rate
    USE elph2,     ONLY : ibndmax, ibndmin, nkqtotf, nktotf, nbndfst
    USE epwcom,    ONLY : nbndsub, nstemp
    USE constants_epw, ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                              meV2invps, eps4
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(out) :: etf_all(nbndsub, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), INTENT(out) :: inv_tau_all(nstemp, nbndfst, nktotf)
    !! Inverse scattering rates
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file. 
    CHARACTER(LEN = 256) :: dummy1
    !! Dummy variable to store the text of the scattering_rate file 
    INTEGER :: ik
    !! K-point index
    INTEGER :: ik_tmp
    !! K-point index read from file
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ibnd_tmp
    !! Local band index read from file
    INTEGER :: ios
    !! Status of reading file
    REAL(KIND = DP) :: temp
    !! Temporary file name used to write scattering rate to file. 
    ! 
    WRITE(stdout,'(/5x,"Reading scattering rate from file"/)')
    !
    IF (mpime == ionode_id) THEN
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1, '(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1, '(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1, '(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate, FILE = name1, STATUS = 'old', IOSTAT = ios)
      WRITE(stdout,'(a16,a22)') '     Open file: ',name1   
      ! There are two comment line at the beginning of the file
      READ(iufilscatt_rate, *) dummy1
      READ(iufilscatt_rate, *) dummy1
      !
      DO ik = 1, nktotf
        !
        DO ibnd = 1, nbndfst
          !
          READ(iufilscatt_rate, *) ik_tmp, ibnd_tmp, etf_all(ibndmin - 1 + ibnd, ik), inv_tau_all(1, ibnd, ik) 
          inv_tau_all(1, ibnd, ik) = inv_tau_all(1, ibnd, ik) / (ryd2mev * meV2invps) 
          !
          ! Check that the file corresponds to the run we are making
          IF (ABS(ibnd_tmp - ibndmin - ibnd + 1) > 0)  CALL errore('scattering_read', &
            'Band read from the scattering_rate file do not match current calculation ', 1)
          ! 
        ENDDO
        ! Check that the file corresponds to the run we are making
        IF (ABS(ik_tmp - ik) > 0)  CALL errore('scattering_read', &
          'k-point read from the scattering_rate file do not match current calculation ', 1)
        !
      ENDDO
      !
      etf_all = etf_all / ryd2ev
      etf_all = etf_all + ef0
      !
      CLOSE(iufilscatt_rate)
    ENDIF
    CALL mp_bcast(etf_all, ionode_id, world_comm)
    CALL mp_bcast(inv_tau_all, ionode_id, world_comm)
    ! 
    WRITE(stdout,'(/5x,"Scattering rate read from file"/)')
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE scattering_read
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE electron_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !! 
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : ibndmax, ibndmin, lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufilsigma_all
    USE io_files,  ONLY : diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: sigmar_all(nbndfst, nktotf)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: sigmai_all(nbndfst, nktotf)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: zi_all(nbndfst, nktotf)
    !! Z parameter of electron-phonon self-energy accross all pools
    ! 
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf + 2)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      !
      lsigma_all = 3 * nbndfst * nktotf + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          i = i + 1
          aux(i) = sigmar_all(ibnd, ik)
        ENDDO
      ENDDO
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          i = i + 1
          aux(i) = sigmai_all(ibnd, ik)
        ENDDO
      ENDDO
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          i = i + 1
          aux(i) = zi_all(ibnd, ik)
        ENDDO
      ENDDO
      CALL diropn(iufilsigma_all, 'sigma_restart', lsigma_all, exst)
      CALL davcio(aux, lsigma_all, iufilsigma_all, 1, +1)
      CLOSE(iufilsigma_all)
    ENDIF
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN 
      sigmar_all(:, 1:lower_bnd - 1) = zero
      sigmai_all(:, 1:lower_bnd - 1) = zero
      zi_all(:, 1:lower_bnd - 1) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      sigmar_all(:, upper_bnd + 1:nktotf) = zero
      sigmai_all(:, upper_bnd + 1:nktotf) = zero
      zi_all(:, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE electron_write
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE electron_read(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !! 
    !! Self-energy reading
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : ibndmax, ibndmin, lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufilsigma_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE constants_epw, ONLY :  zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: sigmar_all(nbndfst, nktotf)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: sigmai_all(nbndfst, nktotf)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: zi_all(nbndfst, nktotf)
    !! Z parameter of electron-phonon self-energy accross all pools
    ! 
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
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
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf + 2)
    !! Vector to store the array
    ! 
    CHARACTER(LEN = 256) :: name1
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.sigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.sigma_restart'
#endif    
      INQUIRE(FILE = name1, EXIST = exst)
      ! 
      IF (exst) THEN ! read the file
        !
        lsigma_all = 3 * nbndfst * nktotf + 2
        CALL diropn(iufilsigma_all, 'sigma_restart', lsigma_all, exst)
        CALL davcio(aux, lsigma_all, iufilsigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('electron_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        ! 
        i = 2
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            sigmar_all(ibnd, ik) = aux(i)
          ENDDO
        ENDDO
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            sigmai_all(ibnd, ik) = aux(i)
          ENDDO
        ENDDO
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            zi_all(ibnd, ik) = aux(i)
          ENDDO
        ENDDO
        CLOSE(iufilsigma_all)
      ENDIF
    ENDIF
    ! 
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(sigmar_all, ionode_id, world_comm)
      CALL mp_bcast(sigmai_all, ionode_id, world_comm)
      CALL mp_bcast(zi_all, ionode_id, world_comm)
      ! 
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        sigmar_all(:, 1:lower_bnd - 1) = zero
        sigmai_all(:, 1:lower_bnd - 1) = zero
        zi_all(:, 1:lower_bnd - 1) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        sigmar_all(:, upper_bnd + 1:nktotf) = zero
        sigmai_all(:, upper_bnd + 1:nktotf) = zero
        zi_all(:, upper_bnd + 1:nktotf) = zero
      ENDIF
      ! 
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE electron_read
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE tau_write(iqq, totq, nktotf, second)
    !----------------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nstemp
    USE io_global, ONLY : meta_ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : diropn
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE constants_epw, ONLY : zero
    !! 
    !! Write scattering rates
    !! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! q-point from the selected ones within the fstick window. 
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    LOGICAL, INTENT(in) :: second
    !! IF we have two Fermi level
    ! 
    ! Local variable
    LOGICAL :: exst
    !! Does the file exists
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
    REAL(KIND = DP) :: aux(2 * nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array inv_tau_all and zi_all
    !
    IF (mpime == meta_ionode_id) THEN
      !
      ltau_all = 2 * nstemp * (nbndfst) * nktotf + 2
      ! First element is the iteration number
      aux(1) = REAL(iqq - 1, KIND = DP)   ! -1 because we will start at the next one. 
      aux(2) = REAL(totq, KIND = DP)
      i = 2
      ! 
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = inv_tau_all(itemp, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i +1
            aux(i) = zi_allvb(itemp, ibnd, ik) 
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iufiltau_all, 'tau_restart', ltau_all, exst)
      CALL davcio(aux, ltau_all, iufiltau_all, 1, +1 )
      CLOSE(iufiltau_all)
      ! 
      IF (second) THEN
        ! First element is the iteration number
        aux(1) = iqq - 1   ! -1 because we will start at the next one. 
        aux(2) = totq
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              aux(i) = inv_tau_allcb(itemp, ibnd, ik)
            ENDDO
          ENDDO
        ENDDO
        !
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              aux(i) = zi_allcb(itemp, ibnd, ik)     
            ENDDO
          ENDDO
        ENDDO
        ! 
        CALL diropn(iufiltau_all, 'tau_restart_CB', ltau_all, exst)
        CALL davcio(aux, ltau_all, iufiltau_all, 1, +1)
        CLOSE(iufiltau_all)   
      ENDIF
      ! 
    ENDIF
    ! 
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) inv_tau_all(:, :, 1:lower_bnd - 1) = zero
    IF (upper_bnd < nktotf) inv_tau_all(:, :, upper_bnd + 1:nktotf) = zero
    IF (second) THEN
      IF (lower_bnd > 1) inv_tau_allcb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) inv_tau_allcb(:, :, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! Same for the Znk factor
    IF (lower_bnd > 1) zi_allvb(:, :, 1:lower_bnd - 1) = zero
    IF (upper_bnd < nktotf) zi_allvb(:, :, upper_bnd + 1:nktotf) = zero
    IF (second) THEN
      IF (lower_bnd > 1) zi_allcb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) zi_allcb(:, :, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE tau_write
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE tau_read(iqq, totq, nktotf, second)
    !----------------------------------------------------------------------------
    !!
    !! Scattering read
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, meta_ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE epwcom,    ONLY : nstemp
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : world_comm
    USE mp_world,  ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point from selecq.fmt
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    LOGICAL, INTENT(in) :: second
    !! IF we have two Fermi level
    ! 
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    CHARACTER(LEN = 256) :: name1
    !! Name of the file
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
    REAL(KIND = DP) :: aux(2 * nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array
    ! 
    IF (mpime == meta_ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart'
#endif 
      INQUIRE(FILE = name1, EXIST = exst)
      ! 
      IF (exst) THEN ! read the file
        !
        ltau_all = 2 * nstemp * nbndfst * nktotf + 2
        CALL diropn(iufiltau_all, 'tau_restart', ltau_all, exst)
        CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('tau_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        ! 
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              inv_tau_all(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        ! 
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              zi_allvb(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO 
        CLOSE(iufiltau_all)
      ENDIF
      ! 
      IF (second) THEN
        ! First inquire if the file exists
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart_CB1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart_CB'
#endif 
        INQUIRE(FILE = name1, EXIST = exst)
        ! 
        IF (exst) THEN ! read the file
          !
          ltau_all = nstemp * nbndfst * nktotf + 2
          CALL diropn(iufiltau_all, 'tau_restart_CB', ltau_all, exst)
          CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
          !
          ! First element is the iteration number
          iqq = INT(aux(1))
          iqq = iqq + 1 ! we need to start at the next q
          nqtotf_read = INT(aux(2))
          IF (nqtotf_read /= totq) CALL errore('tau_read',&
            &'Error: The current total number of q-point is not the same as the read one. ', 1)
          ! 
          i = 2
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                i = i + 1
                inv_tau_allcb(itemp, ibnd, ik) = aux(i)
              ENDDO
            ENDDO
          ENDDO
          ! 
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                i = i + 1
                zi_allcb(itemp, ibnd, ik) = aux(i)
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufiltau_all)
          WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau_CB: ', iqq, '/', totq 
        ENDIF
      ENDIF ! second
    ENDIF
    ! 
    CALL mp_bcast(exst, meta_ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, meta_ionode_id, world_comm)
      CALL mp_bcast(inv_tau_all, meta_ionode_id, world_comm)
      CALL mp_bcast(zi_allvb, meta_ionode_id, world_comm)
      IF (second) CALL mp_bcast(inv_tau_allcb, meta_ionode_id, world_comm)
      IF (second) CALL mp_bcast(zi_allcb, meta_ionode_id, world_comm)
      ! 
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1)      inv_tau_all(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) inv_tau_all(:, :, upper_bnd + 1:nktotf) = zero
      IF (lower_bnd > 1)      zi_allvb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) zi_allvb(:, :, upper_bnd + 1:nktotf) = zero
      !  
      IF (second) THEN
        ! Make everythin 0 except the range of k-points we are working on
        IF (lower_bnd > 1)      inv_tau_allcb(:, :, 1:lower_bnd - 1) = zero
        IF (upper_bnd < nktotf) inv_tau_allcb(:, :, upper_bnd + 1:nktotf) = zero
        IF (lower_bnd > 1)      zi_allcb(:, :, 1:lower_bnd - 1) = zero
        IF (upper_bnd < nktotf) zi_allcb(:, :, upper_bnd + 1:nktotf) = zero
      ENDIF 
      ! 
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau: ', iqq, '/', totq
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE tau_read
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE merge_read(nktotf, nqtotf_new, inv_tau_all_new)
    !----------------------------------------------------------------------------
    !!
    !! File merging
    !! 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : ibndmax, ibndmin, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : tmp_dir, diropn
    USE epwcom,    ONLY : nstemp, restart_filq
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    INTEGER, INTENT(out) :: nqtotf_new
    !! Total number of q-points
    REAL(KIND = DP), INTENT(inout) :: inv_tau_all_new(nstemp, nbndfst, nktotf)
    !! Scattering rate read from file restart_filq
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name of the file 
    LOGICAL :: exst
    !! Does the variable exist
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
    INTEGER(KIND = 8) :: unf_recl
    !! Record length unit
    REAL(KIND = DP) :: aux(nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array 
    REAL(KIND = DP) :: dummy
    !! Test what the record length is
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
      name1 = TRIM(tmp_dir) // TRIM(restart_filq)
      INQUIRE(FILE = name1, EXIST = exst)
      ! 
      IF (exst) THEN ! read the file
        !
        ltau_all = nstemp * nbndfst * nktotf + 2
        !CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
        ! 
        INQUIRE(IOLENGTH = unf_recl) dummy  
        unf_recl = unf_recl * INT(ltau_all, KIND = KIND(unf_recl))
        OPEN(UNIT = iufiltau_all, FILE = restart_filq, IOSTAT = ios, FORM ='unformatted', &
             STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
        !  
        CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
        !
        ! First element is the iteration number
        iq = INT(aux(1))
        iq = iq + 1 ! we need to start at the next q
        nqtotf_new = INT(aux(2))
        ! 
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i +1
              inv_tau_all_new(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufiltau_all)
      ENDIF
    ENDIF
    ! 
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(nqtotf_new, ionode_id, world_comm)
      CALL mp_bcast(inv_tau_all_new, ionode_id, world_comm)
      ! 
      WRITE(stdout, '(a,a)' ) '     Correctly read file ', restart_filq
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE merge_read
    !----------------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------
    SUBROUTINE rwepmatw(epmatw, nbnd, np, nmodes, nrec, iun, iop)
    !-----------------------------------------------------------------
    !!
    !! A simple wrapper to the davcio routine to read/write arrays
    !! instead of vectors 
    !!
    !-----------------------------------------------------------------
    USE kinds, ONLY : DP
    USE mp,    ONLY : mp_barrier
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nbnd
    !! Total number of bands
    INTEGER, INTENT(in) :: np
    !! np is either nrr_k or nq (epmatwe and epmatwp have the same structure)
    INTEGER, INTENT(in) :: nmodes
    !! Number of modes
    INTEGER, INTENT(in) :: nrec
    !! Place where to start reading/writing
    INTEGER, INTENT(in) :: iun
    !! Record number
    INTEGER, INTENT(in) :: iop
    !! If -1, read and if +1 write the matrix
    COMPLEX(KIND = DP), INTENT(inout) :: epmatw(nbnd, nbnd, np, nmodes)
    !! El-ph matrix to read or write
    !
    ! Local variables
    INTEGER :: lrec
    !! Record length
    INTEGER :: i
    !! Index number 
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: imode
    !! Mode index
    INTEGER :: ip
    !! REal space index (either nrr_k or nq)
    COMPLEX(KIND = DP):: aux(nbnd * nbnd * np * nmodes)
    !! 1-D vector to store the matrix elements. 
    !
    lrec = 2 * nbnd * nbnd * np * nmodes
    !
    IF (iop == -1) then
      !
      !  read matrix
      !
      CALL davcio(aux, lrec, iun, nrec, -1)
      !
      i = 0
      DO imode = 1, nmodes
        DO ip = 1, np
          DO jbnd = 1, nbnd
            DO ibnd = 1, nbnd
              i = i + 1
              epmatw(ibnd, jbnd, ip, imode) = aux(i)
              ! 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
    ELSEIF (iop == 1) THEN 
      !
      !  write matrix
      !
      i = 0
      DO imode = 1, nmodes
        DO ip = 1, np
          DO jbnd = 1, nbnd
            DO ibnd = 1, nbnd
              i = i + 1
              aux(i) = epmatw(ibnd, jbnd, ip, imode)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      CALL davcio(aux, lrec, iun, nrec, +1)
      !
    ELSE
      !
      CALL errore('rwepmatw', 'iop not permitted', 1)
      !
    ENDIF
    !
    !----------------------------------------------------------------------
    END SUBROUTINE rwepmatw
    !----------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------
    SUBROUTINE epw_write(nrr_k, nrr_q, nrr_g, w_centers)
    !----------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem
    USE pwcom,     ONLY : ef, nelec
    USE elph2,     ONLY : chw, rdw, cdmew, cvmew, chw_ks, &
                          zstar, epsi, epmatwp
    USE ions_base, ONLY : amass, ityp, nat, tau
    USE cell_base, ONLY : at, bg, omega, alat
    USE phcom,     ONLY : nmodes  
    USE io_var,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp, &
                          crystal
    USE noncollin_module, ONLY : noncolin              
    USE io_files,  ONLY : prefix, diropn
    USE mp,        ONLY : mp_barrier
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id, stdout
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS vectors for the electrons
    INTEGER, INTENT(in) :: nrr_q
    !! Number of WS vectors for the phonons
    INTEGER, INTENT(in) :: nrr_g
    !! Number of WS vectors for the electron-phonons
    REAL(KIND = DP), INTENT(in) :: w_centers(3, nbndsub)
    !! Wannier center
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: filint
    !! Name of the file
    LOGICAL             :: exst
    !! The file exists
    INTEGER :: ibnd, jbnd
    !! Band index
    INTEGER :: jmode, imode
    !! Mode index        
    INTEGER :: irk, irq, irg
    !! WS vector looping index on electron, phonons and el-ph
    INTEGER :: ipol
    !! Cartesian direction (polarison direction)
    INTEGER :: lrepmatw
    !! Record length
    !
    WRITE(stdout,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwdata, FILE = 'epwdata.fmt')
      OPEN(UNIT = crystal, FILE = 'crystal.fmt')
      IF (vme) THEN 
        OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt')
      ELSE
        OPEN(UNIT = iundmedata, FILE = 'dmedata.fmt')
      ENDIF
      IF (eig_read) OPEN(UNIT = iunksdata, FILE = 'ksdata.fmt')
      WRITE(crystal,*) nat
      WRITE(crystal,*) nmodes
      WRITE(crystal,*) nelec
      WRITE(crystal,*) at
      WRITE(crystal,*) bg
      WRITE(crystal,*) omega
      WRITE(crystal,*) alat
      WRITE(crystal,*) tau
      WRITE(crystal,*) amass
      WRITE(crystal,*) ityp
      WRITE(crystal,*) noncolin
      WRITE(crystal,*) w_centers
      !
      WRITE(epwdata,*) ef
      WRITE(epwdata,*) nbndsub, nrr_k, nmodes, nrr_q, nrr_g
      WRITE(epwdata,*) zstar, epsi
      !
      DO ibnd = 1, nbndsub
        DO jbnd = 1, nbndsub
          DO irk = 1, nrr_k
            WRITE (epwdata,*) chw(ibnd, jbnd, irk)
            IF (eig_read) WRITE (iunksdata,*) chw_ks(ibnd, jbnd, irk)
            DO ipol = 1, 3
              IF (vme) THEN 
                WRITE(iunvmedata,*) cvmew(ipol, ibnd, jbnd, irk)
              ELSE
                WRITE(iundmedata,*) cdmew(ipol, ibnd, jbnd, irk)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      DO imode = 1, nmodes
        DO jmode = 1, nmodes
          DO irq = 1, nrr_q
            WRITE(epwdata,*) rdw(imode, jmode, irq) 
          ENDDO
        ENDDO
      ENDDO
      !
      IF (etf_mem == 0) THEN
        ! SP: The call to epmatwp is now inside the loop
        !     This is important as otherwise the lrepmatw INTEGER 
        !     could become too large for integer(kind=4).
        !     Note that in Fortran the record length has to be a integer
        !     of kind 4. 
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        filint   = TRIM(prefix)//'.epmatwp'
        CALL diropn(iunepmatwp, 'epmatwp', lrepmatw, exst)
        DO irg = 1, nrr_g
          CALL davcio(epmatwp(:, :, :, :, irg), lrepmatw, iunepmatwp, irg, +1)
        ENDDO
        ! 
        CLOSE(iunepmatwp)
      ENDIF 
      !
      CLOSE(epwdata)
      CLOSE(crystal)
      IF (vme) THEN 
        CLOSE(iunvmedata)
      ELSE
        CLOSE(iundmedata)
      ENDIF
      IF (eig_read) CLOSE(iunksdata)
      !
    ENDIF
    !--------------------------------------------------------------------------------
    END SUBROUTINE epw_write
    !--------------------------------------------------------------------------------
    ! 
    !--------------------------------------------------------------------------------
    SUBROUTINE epw_read(nrr_k, nrr_q, nrr_g)
    !--------------------------------------------------------------------------------
    USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem, lifc, nqc1, nqc2, nqc3
    USE pwcom,     ONLY : ef
    USE elph2,     ONLY : chw, rdw, epmatwp, cdmew, cvmew, chw_ks, zstar, epsi
    USE ions_base, ONLY : nat
    USE phcom,     ONLY : nmodes  
    USE io_global, ONLY : stdout
    USE io_files,  ONLY : prefix, diropn
    USE io_var,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp
    USE constants_epw, ONLY : czero, zero
#if defined(__NAG)
    USE f90_unix_io,ONLY : flush
#endif
    USE io_global, ONLY : ionode_id
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : inter_pool_comm, world_comm
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: nrr_k
    !! Number of WS vectors for the electrons 
    INTEGER, INTENT(out) :: nrr_q
    !! Number of WS vectors for the phonons
    INTEGER, INTENT(out) :: nrr_g
    !! Number of WS vectors for the electron-phonons
    ! 
    ! Local variables
    ! 
    CHARACTER(LEN = 256) :: filint                                                                                   
    !! Name of the file
    LOGICAL :: exst                           
    !! The file exists
    INTEGER :: ibnd, jbnd
    !! Band index
    INTEGER :: jmode, imode
    !! Mode index        
    INTEGER :: irk, irq, irg
    !! WS vector looping index on electron, phonons and el-ph
    INTEGER :: ipol
    !! Cartesian direction (polarison direction)
    INTEGER :: lrepmatw
    !! Record length
    INTEGER :: ios 
    !! Status of files
    INTEGER :: ierr
    !! Error status
    ! 
    WRITE(stdout,'(/5x,"Reading Hamiltonian, Dynamical matrix and EP vertex in Wann rep from file"/)')
    FLUSH(stdout)
    ! 
    ! This is important in restart mode as zstar etc has not been allocated
    ALLOCATE(zstar(3, 3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating zstar', 1)
    ALLOCATE(epsi(3, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating epsi', 1)
    ! 
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwdata, FILE = 'epwdata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('ephwann_shuffle', 'error opening epwdata.fmt', epwdata)
      IF (eig_read) OPEN(UNIT = iunksdata, FILE = 'ksdata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (eig_read .AND. ios /= 0) CALL errore ('ephwann_shuffle', 'error opening ksdata.fmt', iunksdata)
      IF (vme) THEN 
        OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt', STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) CALL errore ('ephwann_shuffle', 'error opening vmedata.fmt', iunvmedata)
      ELSE
        OPEN(UNIT = iundmedata, FILE = 'dmedata.fmt', STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) CALL errore ('ephwann_shuffle', 'error opening dmedata.fmt', iundmedata)
      ENDIF
      READ(epwdata,*) ef
      READ(epwdata,*) nbndsub, nrr_k, nmodes, nrr_q, nrr_g
      READ(epwdata,*) zstar, epsi
      ! 
    ENDIF
    CALL mp_bcast(ef,      ionode_id, world_comm)
    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr_k,   ionode_id, world_comm)
    CALL mp_bcast(nmodes,  ionode_id, world_comm)
    CALL mp_bcast(nrr_q,   ionode_id, world_comm)
    CALL mp_bcast(nrr_g,   ionode_id, world_comm)
    CALL mp_bcast(zstar,   ionode_id, world_comm)
    CALL mp_bcast(epsi,    ionode_id, world_comm)
    !
    ALLOCATE(chw(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating chw', 1)
    ALLOCATE(chw_ks(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating chw_ks', 1)
    ALLOCATE(rdw(nmodes, nmodes,  nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating rdw', 1)
    IF (vme) THEN 
      ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cvmew', 1)
    ELSE
      ALLOCATE(cdmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cdmew', 1)
    ENDIF
    !
    IF (mpime == ionode_id) THEN
     !
     DO ibnd = 1, nbndsub
       DO jbnd = 1, nbndsub
         DO irk = 1, nrr_k
           READ(epwdata,*) chw(ibnd, jbnd, irk)
           IF (eig_read) READ(iunksdata,*) chw_ks(ibnd, jbnd, irk)
           DO ipol = 1,3
             IF (vme) THEN 
               READ(iunvmedata,*) cvmew(ipol, ibnd, jbnd, irk)
             ELSE
               READ(iundmedata,*) cdmew(ipol, ibnd, jbnd, irk)
             ENDIF
           ENDDO
         ENDDO
       ENDDO
     ENDDO
     !
     IF (.NOT. lifc) THEN
       DO imode = 1, nmodes
         DO jmode = 1, nmodes
           DO irq = 1, nrr_q
             READ(epwdata,*) rdw(imode, jmode, irq)
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     !
    ENDIF
    !
    CALL mp_bcast(chw, ionode_id, world_comm)
    !
    IF (eig_read) CALL mp_bcast(chw_ks, ionode_id, world_comm)
    IF (.NOT. lifc) CALL mp_bcast(rdw, ionode_id, world_comm)
    !
    IF (vme) THEN 
      CALL mp_bcast(cvmew, ionode_id, world_comm)
    ELSE
      CALL mp_bcast(cdmew, ionode_id, world_comm)
    ENDIF
    !
    IF (lifc) THEN
      CALL read_ifc
    ENDIF
    !
    IF (etf_mem == 0) THEN
      ALLOCATE(epmatwp(nbndsub, nbndsub, nrr_k, nmodes, nrr_g), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating epmatwp', 1)
      epmatwp = czero
      IF (mpime == ionode_id) THEN
        ! SP: The call to epmatwp is now inside the loop
        !     This is important as otherwise the lrepmatw INTEGER 
        !     could become too large for integer(kind=4).
        !     Note that in Fortran the record length has to be a integer
        !     of kind 4.      
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        filint   = TRIM(prefix)//'.epmatwp'
        CALL diropn(iunepmatwp, 'epmatwp', lrepmatw, exst)
        DO irg = 1, nrr_g
          CALL davcio(epmatwp(:, :, :, :, irg), lrepmatw, iunepmatwp, irg, -1)
        ENDDO
        !  
        CLOSE(iunepmatwp)
      ENDIF
      !
      CALL mp_bcast(epmatwp, ionode_id, world_comm)
      !
    ENDIF
    !
    !CALL mp_barrier(inter_pool_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(epwdata)
      IF (vme) THEN 
        CLOSE(iunvmedata)
      ELSE
        CLOSE(iundmedata)
      ENDIF
    ENDIF
    !
    WRITE(stdout, '(/5x,"Finished reading Wann rep data from file"/)')
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epw_read
    !------------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat_param(fildyn, ntyp, nat)
    !----------------------------------------------------------------------------
    !! 
    !! Read paramters from the dynamical matrix
    !! 
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256), INTENT(in) :: fildyn
    !! Name of the file to read
    INTEGER, INTENT(out) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(out) :: nat
    !! Number of atoms
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    ! 
    IF (meta_ionode) THEN
      !
      CALL iotk_free_unit(iudyn, ierr)
      !
    ENDIF
    !
    CALL errore('read_dyn_mat_param', 'no free units to write ', ierr)
    IF (meta_ionode) THEN
      !
      ! Open XML descriptor
      ierr = 0
      CALL iotk_open_read(iudyn, FILE = TRIM(fildyn) // '.xml', BINARY = .FALSE., IERR = ierr)
    ENDIF
    !
    CALL errore('read_dyn_mat_param', 'error opening the dyn mat file ', ierr)
    !
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "GEOMETRY_INFO")
      CALL iotk_scan_dat(iudyn, "NUMBER_OF_TYPES", ntyp)
      CALL iotk_scan_dat(iudyn, "NUMBER_OF_ATOMS", nat)
      CALL iotk_scan_end(iudyn, "GEOMETRY_INFO")
    ENDIF
    !  
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat_param
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag,     &
               celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
               nqs, lrigid, epsil, zstareu, lraman, ramtns)
    !----------------------------------------------------------------------------
    !!   
    !! Read the dynamical matrix
    !!
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 3), INTENT(out) :: atm(ntyp)
    !! Atom 
    LOGICAL, INTENT(out), OPTIONAL :: lrigid
    !! 
    LOGICAL, INTENT(out), OPTIONAL :: lraman
    !! Raman 
    INTEGER, INTENT(in) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(out) :: ibrav
    !! Bravais lattice
    INTEGER, INTENT(out) :: nspin_mag
    !! 
    INTEGER, INTENT(out) :: nqs 
    !! 
    INTEGER,  INTENT(out) :: ityp(nat)
    !! Atom type
    REAL(KIND = DP), INTENT(out) :: celldm(6)
    !! Celldm
    REAL(KIND = DP), INTENT(out) :: at(3, 3)
    !! Real-space lattice
    REAL(KIND = DP), INTENT(out) :: bg(3, 3)
    !! Reciprocal-space latrice
    REAL(KIND = DP), INTENT(out) :: omega
    !! Volume of primitive cell
    REAL(KIND = DP), INTENT(out) :: amass(ntyp)
    !! Atom mass
    REAL(KIND = DP), INTENT(out) :: tau(3, nat)
    !! Atom position
    REAL(KIND = DP), INTENT(out) :: m_loc(3, nat)
    !! 
    REAL(KIND = DP), INTENT(out), OPTIONAL :: epsil(3, 3)
    !! Dielectric cst
    REAL(KIND = DP), INTENT(out), OPTIONAL :: zstareu(3, 3, nat)
    !! 
    REAL(KIND = DP), INTENT(out), OPTIONAL :: ramtns(3, 3, 3, nat)
    !! 
    ! 
    ! Local work
    CHARACTER(iotk_attlenx) :: attr
    !! Attribute
    LOGICAL :: found_z
    !! 
    LOGICAL :: lrigid_
    !!  
    INTEGER :: nt
    !! Type of atoms
    INTEGER :: na
    !! Number of atoms
    INTEGER :: kc
    !! Cartesian direction
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxillary
    !
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "GEOMETRY_INFO")
      CALL iotk_scan_dat(iudyn, "BRAVAIS_LATTICE_INDEX", ibrav)
      CALL iotk_scan_dat(iudyn, "SPIN_COMPONENTS", nspin_mag)
      CALL iotk_scan_dat(iudyn, "CELL_DIMENSIONS", celldm)
      CALL iotk_scan_dat(iudyn, "AT", at)
      CALL iotk_scan_dat(iudyn, "BG", bg)
      CALL iotk_scan_dat(iudyn, "UNIT_CELL_VOLUME_AU", omega)
      ! 
      DO nt = 1, ntyp
        CALL iotk_scan_dat(iudyn, "TYPE_NAME"//TRIM(iotk_index(nt)), atm(nt))
        CALL iotk_scan_dat(iudyn, "MASS" // TRIM(iotk_index(nt)), amass(nt))
      ENDDO
      DO na = 1, nat
        CALL iotk_scan_empty(iudyn,"ATOM" // TRIM(iotk_index(na)), attr)
        CALL iotk_scan_attr(attr, "INDEX",  ityp(na))
        CALL iotk_scan_attr(attr, "TAU", tau(:, na))
        IF (nspin_mag == 4) THEN
          CALL iotk_scan_dat(iudyn, "STARTING_MAG_" // TRIM(iotk_index(na)), m_loc(:, na))
        ENDIF 
      ENDDO
      CALL iotk_scan_dat(iudyn, "NUMBER_OF_Q", nqs)

      CALL iotk_scan_end(iudyn, "GEOMETRY_INFO")
      IF (PRESENT(lrigid)) lrigid = .FALSE.
      IF (PRESENT(epsil)) THEN
        CALL iotk_scan_begin(iudyn, "DIELECTRIC_PROPERTIES", FOUND = lrigid_)
        IF (PRESENT(lrigid)) lrigid = lrigid_
        IF (lrigid_) THEN
          CALL iotk_scan_dat(iudyn, "EPSILON", epsil)
          CALL iotk_scan_begin(iudyn, "ZSTAR", FOUND = found_z)
          IF (found_z) THEN
            DO na = 1, nat
              CALL iotk_scan_dat(iudyn, "Z_AT_" // TRIM(iotk_index(na)), aux(:, :))
              IF (PRESENT(zstareu)) zstareu(:, :, na) = aux
            ENDDO
            CALL iotk_scan_end(iudyn, "ZSTAR")
          ELSE
            IF (PRESENT(zstareu)) zstareu = 0.0_DP
          ENDIF
          IF (PRESENT(lraman)) THEN
            CALL iotk_scan_begin(iudyn, "RAMAN_TENSOR_A2", found = lraman)
            IF (lraman) THEN
              DO na = 1, nat
                DO kc = 1, 3
                  CALL iotk_scan_dat(iudyn, "RAMAN_S_ALPHA" // TRIM(iotk_index(na)) &
                      // TRIM(iotk_index(kc)), aux)
                  IF (PRESENT(ramtns)) ramtns(:, :, kc, na) = aux(:, :)
                ENDDO
              ENDDO
              CALL iotk_scan_END(iudyn, "RAMAN_TENSOR_A2")
            ELSE
              IF (PRESENT(ramtns)) ramtns = 0.0_DP
            ENDIF
          ENDIF
          CALL iotk_scan_end(iudyn, "DIELECTRIC_PROPERTIES")
        ELSE
          IF (PRESENT(epsil)) epsil = 0.0_DP
          IF (PRESENT(zstareu)) zstareu = 0.0_DP
          IF (PRESENT(ramtns)) ramtns = 0.0_DP
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat_header
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat(nat, iq, xq, dyn)
    !----------------------------------------------------------------------------
    !!
    !! This routine reads the dynamical matrix file. The file is assumed to
    !! be already opened. iq is the number of the dynamical matrix to read.
    !!
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    REAL(KIND = DP), INTENT(out) :: xq(3)
    !! Q-point value
    COMPLEX(KIND = DP), INTENT(out) :: dyn(3, 3, nat, nat)
    !! Dynamical matrix
    ! 
    ! Local variables
    INTEGER :: na, nb
    !! Number of atoms
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "DYNAMICAL_MAT_" // TRIM(iotk_index(iq)))
      CALL iotk_scan_dat(iudyn, "Q_POINT", xq)
      ! 
      DO na = 1, nat
        DO nb = 1, nat
          CALL iotk_scan_dat(iudyn, "PHI"//TRIM(iotk_index(na)) // TRIM(iotk_index(nb)), dyn(:, :, na, nb))
        ENDDO
      ENDDO
      !  
      CALL iotk_scan_end(iudyn, "DYNAMICAL_MAT_" // TRIM(iotk_index(iq)))
    ENDIF
    !  
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_ifc_param(nr1, nr2, nr3)
    !----------------------------------------------------------------------------
    !! 
    !! Read IFC parameters 
    !! 
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end
    USE kinds,       ONLY : DP
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: nr1, nr2, nr3
    !! Grid size
    ! Local varialbes
    INTEGER :: meshfft(3)
    !! Mesh
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
      CALL iotk_scan_dat(iudyn, "MESH_NQ1_NQ2_NQ3", meshfft)
      nr1 = meshfft(1)
      nr2 = meshfft(2)
      nr3 = meshfft(3)
      CALL iotk_scan_end(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
    ENDIF
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc_param
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_ifc_xml(nr1, nr2, nr3, nat, phid)
    !----------------------------------------------------------------------------
    !! 
    !! Read IFC in XML format
    !!  
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Grid size
    INTEGER, INTENT(in) :: nat
    !! Number of atoms  
    REAL(KIND = DP), INTENT(out) :: phid(nr1 * nr2 * nr3, 3, 3, nat, nat)
    !! 
    ! Local variables
    INTEGER :: na, nb
    !! Atoms
    INTEGER :: nn
    !! 
    INTEGER :: m1, m2, m3
    !! nr dimension 
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxillary
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
      DO na = 1, nat
        DO nb = 1, nat
          nn = 0
          DO m3 = 1, nr3
            DO m2 = 1, nr2
              DO m1 = 1, nr1
                nn = nn + 1
                CALL iotk_scan_begin(iudyn, "s_s1_m1_m2_m3" //     &
                    TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) // &
                    TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) // &
                    TRIM(iotk_index(m3)))
                CALL iotk_scan_dat(iudyn, 'IFC', aux)
                phid(nn, :, :, na, nb) = aux(:, :)
                CALL iotk_scan_end(iudyn, "s_s1_m1_m2_m3" //        &
                     TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) // &
                     TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) // &
                     TRIM(iotk_index(m3)))
              ENDDO ! m1
            ENDDO ! m2
          ENDDO ! m3
        ENDDO ! nb
      ENDDO ! na
      CALL iotk_scan_end(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
      CALL iotk_close_read(iudyn)
    ENDIF ! meta_ionode
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc_xml
    !----------------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------------------
    SUBROUTINE read_ifc()
    !---------------------------------------------------------------------------------
    !!
    !! Read IFC in real space from the file generated by q2r. 
    !! Adapted from PH/matdyn.x by C. Verdi and S. Ponce
    !! 
    !
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : ifc, zstar, epsi
    USE epwcom,    ONLY : asr_typ, dvscf_dir, nqc1, nqc2, nqc3
    USE ions_base, ONLY : nat
    USE cell_base, ONLY : ibrav, omega, at, bg, celldm, alat
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iunifc
    USE noncollin_module, ONLY : nspin_mag
    USE io_global, ONLY : ionode_id
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, root_pool
    USE mp_world,  ONLY : mpime, world_comm
#if defined(__NAG)
    USE f90_unix_io, ONLY : flush
#endif
    !
    IMPLICIT NONE
    !
    ! Local variables 
    LOGICAL :: lpolar_
    !! Polar flag
    LOGICAL :: has_zstar
    !! Does it has Born effective charges
    LOGICAL :: is_plain_text_file
    !! Is the file txt
    LOGICAL :: is_xml_file
    !! Is the file XML
    CHARACTER(LEN = 80) :: line
    !! 
    CHARACTER(LEN = 256) :: tempfile
    !! 
    CHARACTER(LEN = 3), ALLOCATABLE :: atm(:)
    !! 
    INTEGER :: ios
    !! 
    INTEGER :: i, j
    !! 
    INTEGER :: m1, m2, m3
    !! 
    INTEGER :: na, nb
    !! 
    INTEGER :: idum
    !! 
    INTEGER :: ibid, jbid
    !! 
    INTEGER :: nabid
    !! 
    INTEGER :: nbbid
    !! 
    INTEGER :: m1bid, m2bid, m3bid
    !! 
    INTEGER :: ntyp_
    !! 
    INTEGER :: nat_
    !! 
    INTEGER :: ibrav_
    !! 
    INTEGER :: ityp_(nat)
    !! 
    INTEGER :: nqs
    !! 
    INTEGER :: ierr
    !! Error status
    INTEGER, PARAMETER :: ntypx = 10
    !! 
    REAL(KIND = DP):: tau_(3, nat)
    !! 
    REAL(KIND = DP):: amass2(ntypx)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: m_loc(:, :)
    !! 
    !
    WRITE(stdout, '(/5x,"Reading interatomic force constants"/)')
    FLUSH(stdout)
    ! 
    ! Generic name for the ifc.q2r file. If it is xml, the file will be named ifc.q2r.xml instead
    tempfile = TRIM(dvscf_dir) // 'ifc.q2r'
    ! The following function will check if the file exists in xml format
    CALL check_is_xml_file(tempfile, is_xml_file)
    ! 
    IF (mpime == ionode_id) THEN
      ! 
      IF (is_xml_file) THEN
        ! pass the 'tempfile' as the '.xml' extension is added in the next routine
        CALL read_dyn_mat_param(tempfile, ntyp_, nat_)
        ALLOCATE(m_loc(3, nat_), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error allocating m_loc', 1)
        ALLOCATE(atm(ntyp_), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error allocating atm', 1)
        CALL read_dyn_mat_header(ntyp_, nat_, ibrav, nspin_mag, &
                 celldm, at, bg, omega, atm, amass2, &
                 tau_, ityp_,  m_loc, nqs, has_zstar, epsi, zstar)
        call volume(alat, at(1, 1), at(1, 2), at(1, 3), omega)
        CALL read_ifc_param(nqc1, nqc2, nqc3)
        CALL read_ifc_xml(nqc1, nqc2, nqc3, nat_, ifc)
        DEALLOCATE(m_loc, STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error deallocating m_loc', 1)
        DEALLOCATE(atm, STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error deallocating atm', 1)
        ! 
      ELSE
        !
        OPEN(UNIT = iunifc, FILE = tempfile, STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) call errore ('read_ifc', 'error opening ifc.q2r', iunifc)
        !
        !  read real-space interatomic force constants
        !
        READ(iunifc,'(3i4)') ntyp_ , nat_ , ibrav_
        IF (ibrav_ == 0) THEN
          DO i = 1, 3
            READ(iunifc, *) line
          ENDDO
        ENDIF
        DO i = 1, ntyp_
          READ(iunifc, '(a)') line
        ENDDO
        DO na = 1, nat
          READ(iunifc, *) idum, idum, (tau_(j, na), j = 1, 3)
        ENDDO
        READ(iunifc, *) lpolar_
        !
        IF (lpolar_) THEN
          READ (iunifc,*) ((epsi(i, j), j = 1, 3), i = 1, 3)
          DO na = 1, nat
             READ(iunifc, *) idum
             READ(iunifc, *) ((zstar(i, j, na), j = 1, 3), i = 1, 3)
          ENDDO
          WRITE(stdout, '(5x,a)') "Read Z* and epsilon"
        ENDIF
        !
        READ (iunifc,*) idum
        !
        ifc = 0.d0
        DO i = 1, 3
          DO j = 1, 3
            DO na = 1, nat
              DO nb = 1, nat
                READ(iunifc, *) ibid, jbid, nabid, nbbid
                IF (i /= ibid .OR. j /= jbid .OR. na /= nabid .OR. nb /= nbbid)  &
                  CALL errore('read_epw', 'error in reading ifc', 1)
                READ(iunifc, *) (((m1bid, m2bid, m3bid, ifc(m1, m2, m3, i, j, na, nb), &
                           m1 = 1, nqc1), m2 = 1, nqc2), m3 = 1, nqc3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF ! noncol
    ENDIF
    !
    ! It has to be casted like this because mpi cannot cast 7 indices
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            CALL mp_bcast(ifc(:, :, :, i, j, na, nb), ionode_id, inter_pool_comm)
            CALL mp_bcast(ifc(:, :, :, i, j, na, nb), root_pool, intra_pool_comm)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_bcast(zstar, ionode_id, world_comm)
    CALL mp_bcast(epsi, ionode_id, world_comm)
    CALL mp_bcast(tau_, ionode_id, world_comm)
    CALL mp_bcast(ibrav_, ionode_id, world_comm)
    !
    WRITE(stdout,'(5x,"IFC last ", 1f12.7)') ifc(nqc1, nqc2, nqc3, 3, 3, nat, nat)
    !
    CALL set_asr2(asr_typ, nqc1, nqc2, nqc3, ifc, zstar, nat, ibrav_, tau_)
    !
    !CALL mp_barrier(inter_pool_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(iunifc)
    ENDIF
    !
    WRITE(stdout, '(/5x,"Finished reading ifcs"/)')
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE read_ifc
    !-------------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE set_asr2(asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
    !-----------------------------------------------------------------------
    !!
    !! Set the acoustic sum rule. 
    !! Taken directly from PHonon/PH/q2trans.f90
    !! It would be better to take the set_asr for /Modules/. 
    !! However they are different (frc) and to be consitent with q2r.x we take this one.
    !
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    !
    IMPLICIT NONE
    ! 
    CHARACTER(LEN = 10), INTENT(in) :: asr
    !! Acoustic sum rule
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! 
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(in) :: ibrav
    !! Bravais lattice
    REAL(KIND = DP), INTENT(in) :: tau(3, nat)
    !! Atomic position
    REAL(KIND = DP), INTENT(inout) :: frc(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Real-space IFC
    REAL(KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! Z*
    !
    ! Local variables
    TYPE vector
      REAL(KIND = DP), pointer :: vec(:, :, :, :, :, :, :)
    END TYPE vector
    !
    TYPE (vector) u(6 * 3 * nat)
    ! These are the "vectors" associated with the sum rules on force-constants
    INTEGER :: axis
    !! 
    INTEGER :: n
    !! 
    INTEGER :: i, j
    !! 
    INTEGER :: na, nb
    !! 
    INTEGER :: n1, n2, n3
    !! 
    INTEGER :: m, p, k, l, q, r
    !! 
    INTEGER :: i1, j1 
    !! 
    INTEGER :: na1
    !! 
    INTEGER :: ierr
    !! Error status
    INTEGER :: u_less(6 * 3 * nat)
    !! indices of the vectors u that are not independent to the preceding ones
    INTEGER :: n_less
    !! Number of index that are no independent
    INTEGER :: i_less
    !! temporary parameter
    INTEGER :: zeu_less(6 * 3)
    !! indices of the vectors zeu_u that are not independent to the preceding ones,
    INTEGER :: nzeu_less
    !! number of such vectors
    INTEGER :: izeu_less 
    !!  temporary parameter
    INTEGER, ALLOCATABLE :: ind_v(:, :, :)
    !! 
    REAL(KIND = DP) :: zeu_new(3, 3, nat)
    !! 
    REAL(KIND = DP) :: scal
    !! 
    REAL(KIND = DP) :: norm2
    !! 
    REAL(KIND = DP) :: sum
    !! 
    REAL(KIND = DP) :: zeu_u(6 * 3, 3, 3, nat)
    !! These are the "vectors" associated with the sum rules on effective charges
    REAL(KIND = DP) :: zeu_w(3, 3, nat)
    !! 
    REAL(KIND = DP) :: zeu_x(3, 3, nat)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: w(:, :, :, :, :, :, :)
    !! temporary vectors and parameters
    REAL(KIND = DP), ALLOCATABLE :: x(:, :, :, :, :, :, :)
    !! temporary vectors and parameters
    REAL(KIND = DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
    !
    REAL(KIND = DP), ALLOCATABLE :: v(:, :)
    !! These are the "vectors" associated with symmetry conditions, coded by
    !! indicating the positions (i.e. the seven indices) of the non-zero elements (there
    !! should be only 2 of them) and the value of that element. We do so in order
    !! to limit the amount of memory used.
    !
    ! Initialization. n is the number of sum rules to be considered (if
    ! asr/='simple')
    ! and 'axis' is the rotation axis in the case of a 1D system
    ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if
    ! axis='3')
    !
    IF ((asr /= 'simple') .AND. (asr /= 'crystal') .AND. (asr /= 'one-dim') .AND. (asr /= 'zero-dim')) THEN
      CALL errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
    ENDIF
    !
    IF (asr == 'simple') THEN
      !
      ! Simple Acoustic Sum Rule on effective charges
      !
      DO i = 1, 3
        DO j = 1, 3
          sum = 0.0d0
          DO na = 1, nat
            sum = sum + zeu(i, j, na)
          ENDDO
          DO na = 1, nat
            zeu(i, j, na) = zeu(i, j, na) - sum / nat
          ENDDO
        ENDDO
      ENDDO
      !
      ! Simple Acoustic Sum Rule on force constants in real space
      !
      DO i = 1, 3
        DO j = 1, 3
          DO na = 1, nat
            sum = 0.0d0
            DO nb = 1, nat
              DO n1 = 1, nr1
                DO n2 = 1, nr2
                  DO n3 = 1, nr3
                    sum = sum + frc(n1, n2, n3, i, j, na, nb)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            frc(1, 1, 1, i, j, na, na) = frc(1, 1, 1, i, j, na, na) - sum
          ENDDO
        ENDDO
      ENDDO
      WRITE (stdout,'(5x,a)') " Imposed simple ASR"
      RETURN
      !
    ENDIF ! simple
    ! 
    IF (asr == 'crystal') n = 3
    IF (asr == 'one-dim') THEN
      ! the direction of periodicity is the rotation axis
      ! It will work only if the crystal axis considered is one of
      ! the cartesian axis (typically, ibrav = 1, 6 or 8, or 4 along the z-direction)
      IF (nr1 * nr2 * nr3 == 1) axis = 3
      IF ((nr1 /= 1) .AND. (nr2 * nr3 == 1)) axis = 1
      IF ((nr2 /= 1) .AND. (nr1 * nr3 == 1)) axis = 2
      IF ((nr3 /= 1) .AND. (nr1 * nr2 == 1)) axis = 3
      IF (((nr1 /= 1) .AND. (nr2 /=1 )) .OR. ((nr2 /= 1) .AND. (nr3 /= 1)) .OR. ((nr1 /= 1) .AND. (nr3 /=1 ))) THEN
        CALL errore('set_asr', 'too many directions of periodicity in 1D system', axis)
      ENDIF
      IF ((ibrav /= 1) .AND. (ibrav /= 6) .AND. (ibrav /= 8) .AND. ((ibrav /= 4) .OR. (axis /= 3))) THEN
        WRITE(stdout, *) 'asr: rotational axis may be wrong'
      ENDIF
      WRITE(stdout, '("asr rotation axis in 1D system= ", I4)') axis
      n = 4
    ENDIF
    IF(asr == 'zero-dim') n = 6
    !
    ! Acoustic Sum Rule on effective charges
    !
    ! generating the vectors of the orthogonal of the subspace to project the effective charges matrix on
    !
    zeu_u(:, :, :, :) = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          zeu_new(i, j, na) = zeu(i, j, na)
        ENDDO
      ENDDO
    ENDDO
    !
    p = 0
    DO i = 1, 3
      DO j = 1, 3
         ! These are the 3*3 vectors associated with the
         ! translational acoustic sum rules
         p = p + 1
         zeu_u(p, i, j, :) = 1.0d0
         !
      ENDDO
    ENDDO
    !
    IF (n == 4) THEN
      DO i = 1, 3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p = p + 1
        DO na = 1, nat
          zeu_u(p, i, MOD(axis, 3) + 1, na) = -tau(MOD(axis + 1, 3) + 1,na)
          zeu_u(p, i, MOD(axis + 1, 3) + 1, na) = tau(MOD(axis, 3) + 1, na)
        ENDDO
      ENDDO
    ENDIF
    !
    IF (n == 6) THEN
      DO i = 1, 3
        DO j = 1, 3
          ! These are the 3*3 vectors associated with the
          ! three rotational sum rules (0D system - typ. molecule)
          p = p + 1
          DO na = 1, nat
            zeu_u(p, i, MOD(j, 3) + 1, na) = -tau(MOD(j + 1, 3) + 1, na)
            zeu_u(p, i, MOD(j + 1, 3) + 1, na) = tau(MOD(j, 3) + 1, na)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Gram-Schmidt orthonormalization of the set of vectors created.
    !
    nzeu_less = 0
    DO k = 1, p
      zeu_w(:, :, :) = zeu_u(k, :, :, :)
      zeu_x(:, :, :) = zeu_u(k, :, :, :)
      DO q = 1, k - 1
        r = 1
        DO izeu_less = 1, nzeu_less
          IF (zeu_less(izeu_less) == q) r = 0
        ENDDO
        IF (r /= 0) THEN
          CALL sp_zeu(zeu_x, zeu_u(q, :, :, :), nat,scal)
          zeu_w(:, :, :) = zeu_w(:, :, :) - scal * zeu_u(q, :, :, :)
        ENDIF
      ENDDO
      CALL sp_zeu(zeu_w, zeu_w, nat, norm2)
      IF (norm2 > 1.0d-16) THEN
        zeu_u(k,:,:,:) = zeu_w(:, :, :) / DSQRT(norm2)
      ELSE
        nzeu_less = nzeu_less + 1
        zeu_less(nzeu_less) = k
      ENDIF
    ENDDO
    !
    ! Projection of the effective charge "vector" on the orthogonal of the
    ! subspace of the vectors verifying the sum rules
    !
    zeu_w(:, :, :) = 0.0d0
    DO k = 1, p
      r = 1
      DO izeu_less = 1,nzeu_less
        IF (zeu_less(izeu_less) == k) r = 0
      ENDDO
      IF (r /= 0) THEN
        zeu_x(:, :, :) = zeu_u(k, :, :, :)
        CALL sp_zeu(zeu_x, zeu_new, nat, scal)
        zeu_w(:, :, :) = zeu_w(:, :, :) + scal * zeu_u(k, :, :, :)
      ENDIF
    ENDDO
    !
    ! Final substraction of the former projection to the initial zeu, to get
    ! the new "projected" zeu
    !
    zeu_new(:, :, :) = zeu_new(:, :, :) - zeu_w(:, :, :)
    CALL sp_zeu(zeu_w, zeu_w, nat, norm2)
    WRITE(stdout, '(5x,"Norm of the difference between old and new effective charges: ", 1f12.7)') SQRT(norm2)
    !
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          zeu(i, j, na) = zeu_new(i, j, na)
        ENDDO
      ENDDO
    ENDDO
    !
    ! Acoustic Sum Rule on force constants
    !
    ! generating the vectors of the orthogonal of the subspace to project
    ! the force-constants matrix on
    !
    DO k = 1, 18 * nat
      ALLOCATE(u(k) % vec(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating u(k) % vec', 1)
      u(k) % vec (:, :, :, :, :, :, :) = 0.0d0
    ENDDO
    ALLOCATE(frc_new(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating frc_new', 1)
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  frc_new(n1, n2, n3, i, j, na, nb) = frc(n1, n2, n3, i, j, na, nb)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    p = 0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          ! These are the 3*3*nat vectors associated with the
          ! translational acoustic sum rules
          p = p + 1
          u(p) % vec(:, :, :, i, j, na, :) = 1.0d0
          !
        ENDDO
      ENDDO
    ENDDO
    !
    IF (n == 4) THEN
      DO i = 1, 3
        DO na = 1, nat
          ! These are the 3*nat vectors associated with the
          ! single rotational sum rule (1D system)
          p = p + 1
          DO nb = 1, nat
            u(p) % vec(: ,:, :, i, MOD(axis, 3) + 1, na, nb) = -tau(MOD(axis + 1, 3) + 1, nb)
            u(p) % vec(:, :, :, i, MOD(axis + 1, 3) + 1, na, nb) = tau(MOD(axis, 3) + 1, nb)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    IF (n == 6) THEN
      DO i = 1, 3
        DO j = 1, 3
          DO na = 1, nat
            ! These are the 3*3*nat vectors associated with the
            ! three rotational sum rules (0D system - typ. molecule)
            p = p + 1
            DO nb = 1, nat
              u(p) % vec(:, :, :, i, MOD(j, 3) + 1, na, nb) = -tau(MOD(j + 1, 3) + 1, nb)
              u(p) % vec(:, :, :, i, MOD(j + 1, 3) + 1, na, nb) = tau(MOD(j, 3) + 1, nb)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ALLOCATE(ind_v(9 * nat * nat * nr1 * nr2 * nr3, 2, 7), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating ind_v', 1)
    ALLOCATE(v(9 * nat * nat * nr1 * nr2 * nr3, 2), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating v', 1)
    m = 0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  ! These are the vectors associated with the symmetry
                  ! constraints
                  q = 1
                  l = 1
                  DO WHILE((l <= m) .AND. (q /= 0))
                    IF ((ind_v(l, 1, 1) == n1) .AND. (ind_v(l, 1, 2) == n2) .AND. &
                        (ind_v(l, 1, 3) == n3) .AND. (ind_v(l, 1, 4) == i) .AND. &
                        (ind_v(l, 1, 5) == j) .AND. (ind_v(l, 1, 6) == na) .AND. &
                        (ind_v(l, 1, 7) == nb)) q = 0
                    IF ((ind_v(l, 2, 1) == n1) .AND. (ind_v(l, 2, 2) == n2) .AND. &
                        (ind_v(l, 2, 3) == n3) .AND. (ind_v(l, 2, 4) == i) .AND. &
                        (ind_v(l, 2, 5) == j) .AND. (ind_v(l, 2, 6) == na) .AND. &
                        (ind_v(l, 2, 7) == nb)) q = 0
                    l = l + 1
                  ENDDO
                  IF ((n1 == MOD(nr1 + 1 - n1, nr1) + 1) .AND. (n2 == MOD(nr2 + 1 - n2, nr2) + 1) &
                       .AND. (n3 == MOD(nr3 + 1 - n3, nr3) + 1) .AND. (i == j) .AND. (na == nb)) q = 0
                  IF (q /= 0) THEN
                    m = m + 1
                    ind_v(m, 1, 1) = n1
                    ind_v(m, 1, 2) = n2
                    ind_v(m, 1, 3) = n3
                    ind_v(m, 1, 4) = i
                    ind_v(m, 1, 5) = j
                    ind_v(m, 1, 6) = na
                    ind_v(m, 1, 7) = nb
                    v(m, 1) = 1.0d0 / DSQRT(2.0d0)
                    ind_v(m, 2, 1) = MOD(nr1 + 1 - n1, nr1) + 1
                    ind_v(m, 2, 2) = MOD(nr2 + 1 - n2, nr2) + 1
                    ind_v(m, 2, 3) = MOD(nr3 + 1 - n3, nr3) + 1
                    ind_v(m, 2, 4) = j
                    ind_v(m, 2, 5) = i
                    ind_v(m, 2, 6) = nb
                    ind_v(m, 2, 7) = na
                    v(m, 2) = -1.0d0 / DSQRT(2.0d0)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    ! Gram-Schmidt orthonormalization of the set of vectors created.
    ! Note that the vectors corresponding to symmetry constraints are already
    ! orthonormalized by construction.
    !
    n_less = 0
    ALLOCATE(w(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating w', 1)
    ALLOCATE(x(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating x', 1)
    DO k = 1, p
      w(:, :, :, :, :, :, :) = u(k) % vec(:, :, :, :, :, :, :)
      x(:, :, :, :, :, :, :) = u(k) % vec(:, :, :, :, :, :, :)
      DO l = 1, m
        !
        CALL sp2(x, v(l, :), ind_v(l, :, :), nr1, nr2, nr3, nat, scal)
        DO r = 1, 2
          n1 = ind_v(l, r, 1)
          n2 = ind_v(l, r, 2)
          n3 = ind_v(l, r, 3)
          i = ind_v(l, r, 4)
          j = ind_v(l, r, 5)
          na = ind_v(l, r, 6)
          nb = ind_v(l, r, 7)
          w(n1, n2, n3, i, j, na, nb) = w(n1, n2, n3, i, j, na, nb) - scal * v(l, r)
        ENDDO
      ENDDO
      IF (k <= (9 * nat)) THEN
        na1 = MOD(k, nat)
        IF (na1 == 0) na1 = nat
        j1 = MOD((k - na1) / nat, 3) + 1
        i1 = MOD((((k - na1) / nat) - j1 + 1) / 3, 3) + 1
      ELSE
        q = k - 9 * nat
        IF (n == 4) THEN
          na1 = MOD(q, nat)
          IF (na1 == 0) na1 = nat
          i1 = MOD((q - na1) / nat, 3) + 1
        ELSE
          na1 = MOD(q, nat)
          IF (na1 == 0) na1 = nat
          j1 = MOD((q - na1) / nat, 3) + 1
          i1 = MOD((((q - na1) / nat) - j1 + 1) / 3, 3) + 1
        ENDIF
      ENDIF
      DO q = 1, k - 1
        r = 1
        DO i_less = 1, n_less
          IF (u_less(i_less) == q) r = 0
        ENDDO
        IF (r /= 0) THEN
          CALL sp3(x, u(q) % vec(:, :, :, :, :, :, :), i1, na1, nr1, nr2, nr3, nat, scal)
          w(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) - scal * u(q) % vec(:, :, :, :, :, :, :)
        ENDIF
      ENDDO
      CALL sp1(w, w, nr1, nr2, nr3, nat, norm2)
      IF (norm2 > 1.0d-16) THEN
        u(k) % vec(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) / DSQRT(norm2)
      ELSE
        n_less = n_less + 1
        u_less(n_less) = k
      ENDIF
    ENDDO
    !
    ! Projection of the force-constants "vector" on the orthogonal of the
    ! subspace of the vectors verifying the sum rules and symmetry contraints
    !
    w(:, :, :, :, :, :, :) = 0.0d0
    DO l = 1, m
      CALL sp2(frc_new, v(l, :), ind_v(l, :, :), nr1, nr2, nr3, nat, scal)
      DO r = 1, 2
        n1 = ind_v(l, r, 1)
        n2 = ind_v(l, r, 2)
        n3 = ind_v(l, r, 3)
        i = ind_v(l, r, 4)
        j = ind_v(l, r, 5)
        na = ind_v(l, r, 6)
        nb = ind_v(l, r, 7)
        w(n1, n2, n3, i, j, na, nb) = w(n1, n2, n3, i, j, na, nb) + scal * v(l, r)
      ENDDO
    ENDDO
    DO k = 1, p
      r = 1
      DO i_less = 1,n_less
        IF (u_less(i_less) == k) r = 0
      ENDDO
      IF (r /= 0) THEN
        x(:, :, :, :, :, :, :) = u(k) % vec(:, :, :, :, :, :, :)
        call sp1(x, frc_new, nr1, nr2, nr3, nat, scal)
        w(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) + scal * u(k)%vec(:, :, :, :, :, :, :)
      ENDIF
      DEALLOCATE(u(k)%vec, STAT = ierr)
      IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating u(k)%vec', 1)
    ENDDO
    !
    ! Final substraction of the former projection to the initial frc, to get
    ! the new "projected" frc
    !
    frc_new(:, :, :, :, :, :, :) = frc_new(:, :, :, :, :, :, :) - w(:, :, :, :, :, :, :)
    CALL sp1(w, w, nr1, nr2, nr3, nat, norm2)
    ! 
    WRITE(stdout, '(5x,"Norm of the difference between old and new force-constants: ", 1f12.7)') SQRT(norm2)
    !
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  frc(n1, n2, n3, i, j, na, nb) = frc_new(n1, n2, n3, i, j, na, nb)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(x, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating x', 1)
    DEALLOCATE(w, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating w', 1)
    DEALLOCATE(v, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating v', 1)
    DEALLOCATE(ind_v, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating ind_v', 1)
    DEALLOCATE(frc_new, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating frc_new', 1)
    WRITE(stdout, '(5x,a)') "Imposed crystal ASR"
    !
    RETURN
    !----------------------------------------------------------------------
    END SUBROUTINE set_asr2
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp_zeu(zeu_u, zeu_v, nat, scal)
    !-----------------------------------------------------------------------
    !!
    !! Does the scalar product of two effective charges matrices zeu_u and zeu_v
    !! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
    !!
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE
    !  
    INTEGER, INTENT(in) :: nat
    !! 
    REAL(KIND = DP), INTENT(in) :: zeu_u(3, 3, nat)
    !! 
    REAL(KIND = DP), INTENT(in) :: zeu_v(3, 3, nat)
    !! 
    REAL(KIND = DP), INTENT(inout) :: scal
    !! 
    ! Local variables
    INTEGER :: i
    !! 
    INTEGER :: j
    !! 
    INTEGER :: na
    !! 
    !
    scal = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          scal = scal + zeu_u(i, j, na) * zeu_v(i, j, na)
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sp_zeu
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp1(u, v, nr1, nr2, nr3, nat, scal)
    !-----------------------------------------------------------------------
    !!
    !! Does the scalar product of two force-constants matrices u and v (considered as
    !! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
    !!
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE 
    ! 
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Supercell dims
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(in) :: u(nr1, nr2, nr3, 3, 3, nat, nat)
    !! First force-constent matrix 
    REAL(KIND = DP), INTENT(in) :: v(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Second force-constent matrix
    REAL(KIND = DP), INTENT(out) :: scal
    !! Scalar product
    ! 
    ! Local variables
    INTEGER ::  i, j 
    !! Cartesian direction
    INTEGER :: na, nb
    !! Atoms index
    INTEGER :: n1, n2, n3
    !! Supercell index
    !
    scal = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  scal = scal + u(n1, n2, n3, i, j, na, nb) * v(n1, n2, n3, i, j, na, nb)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
    !
    !----------------------------------------------------------------------
    END SUBROUTINE sp1
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp2(u, v, ind_v, nr1, nr2, nr3, nat, scal)
    !-----------------------------------------------------------------------
    !!
    !! Does the scalar product of two force-constants matrices u and v (considered as
    !! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
    !! but v is coded as explained when defining the vectors corresponding to the symmetry constraints
    !!
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Supercell dims
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(in) :: ind_v(2, 7)
    !! Index vector
    REAL(KIND = DP), INTENT(in) :: u(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Input vector
    REAL(KIND = DP), INTENT(in) :: v(2)
    !! input vector
    REAL(KIND = DP), INTENT(out) :: scal
    !! Output vector
    ! 
    ! Local variables
    INTEGER ::  i
    !! Index
    !
    scal = 0.0d0
    DO i = 1, 2
      scal = scal + u(ind_v(i, 1), ind_v(i, 2), ind_v(i, 3), ind_v(i, 4), ind_v(i, 5), ind_v(i, 6), &
           ind_v(i, 7)) * v(i)
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sp2
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp3(u, v, i, na, nr1, nr2, nr3, nat, scal)
    !-----------------------------------------------------------------------
    !
    ! like sp1, but in the particular case when u is one of the u(k)%vec
    ! defined in set_asr (before orthonormalization). In this case most of the
    ! terms are zero (the ones that are not are characterized by i and na), so
    ! that a lot of computer time can be saved (during Gram-Schmidt).
    !
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Supercell dims
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(in) :: u(nr1, nr2, nr3, 3, 3, nat, nat)
    !! First force-constent matrix 
    REAL(KIND = DP), INTENT(in) :: v(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Second force-constent matrix
    REAL(KIND = DP), INTENT(out) :: scal
    !! Scalar product
    !
    ! Local variables
    INTEGER ::  i, j
    !! Cartesian direction
    INTEGER :: na, nb
    !! Atoms index
    INTEGER :: n1, n2, n3
    !! Supercell index
    !
    scal = 0.0d0
    DO j = 1,3
      DO nb = 1, nat
        DO n1 = 1, nr1
          DO n2 = 1, nr2
            DO n3 = 1, nr3
              scal = scal + u(n1, n2, n3, i, j, na, nb) * v(n1, n2, n3, i, j, na, nb)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    return
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE sp3
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE check_is_xml_file(filename, is_xml_file)
    !-------------------------------------------------------------------------------
    !!
    !! This SUBROUTINE checks if a file is formatted in XML. It does so by
    !! checking if the file exists and if the file + '.xml' in its name exists.
    !! If both of them or none of them exists, an error is raised. If only one of
    !! them exists, it sets the 'is_xml_file' to .TRUE. of .FALSE. depending of
    !! the file found.
    !!
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256), INTENT(in) :: filename
    !! The name of the file to check if formatted in XML format
    !! This string is assumed to be trimmed
    LOGICAL, INTENT(out) :: is_xml_file
    !! Is .TRUE. if the file is in xml format. .FALSE. otherwise.
    !
    ! Local variables
    CHARACTER(LEN = 256) :: filename_xml
    !! File name
    CHARACTER(LEN = 256) :: errmsg
    !! Error message 
    LOGICAL :: is_plain_text_file
    !! Plain tex  t
    ! 
    filename_xml = TRIM(filename) // '.xml'
    filename_xml = TRIM(filename_xml)
    INQUIRE(FILE = filename, EXIST = is_plain_text_file)
    INQUIRE(FILE = filename_xml, EXIST = is_xml_file)
    ! Tell user if any inconsistencies
    IF (is_xml_file .AND. is_plain_text_file) THEN
      ! 2 different type of files exist => warn user
      errmsg = "Detected both: '" // filename // "' and '" // filename_xml // &
              &"' which one to choose?"
      CALL errore('check_is_xml_file', errmsg, 1)
    ELSEIF (.NOT. is_xml_file .AND. .NOT. is_plain_text_file) THEN
      errmsg = "Expected a file named either '" // filename //"' or '"&
              &// filename_xml // "' but none was found."
      CALL errore('check_is_xml_file', errmsg, 1)
    ENDIF
    ! else one of the file in exists
    !------------------------------------------------------------------------------
    END SUBROUTINE check_is_xml_file
    !------------------------------------------------------------------------------
    ! 
    !-------------------------------------------------------------
    SUBROUTINE readdvscf(dvscf, recn, iq, nqc)
    !-------------------------------------------------------------
    !!
    !! Open dvscf files as direct access, read, and close again
    !!
    !! SP - Nov 2017
    !! Replaced fstat by Fortran instric function inquire. 
    !! 
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !-------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_files,  ONLY : prefix
    USE units_ph,  ONLY : lrdrho
    USE fft_base,  ONLY : dfftp
    USE epwcom,    ONLY : dvscf_dir
    USE io_var,    ONLY : iudvscf
    USE low_lvl,   ONLY : set_ndnmbr 
    USE noncollin_module, ONLY : nspin_mag
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: recn
    !! perturbation number
    INTEGER, INTENT(in) :: iq
    !! the current q-point
    INTEGER, INTENT(in) :: nqc
    !! the total number of q-points in the list
    COMPLEX(KIND = DP), INTENT(inout) :: dvscf(dfftp%nnr, nspin_mag) 
    !! dVscf potential is read from file
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: tempfile
    !! Temp file 
    CHARACTER(LEN = 3) :: filelab
    !! File number
    INTEGER :: unf_recl
    !! Rcl unit
    INTEGER :: ios
    !! Error number
    INTEGER(KIND = 8) :: mult_unit
    !! Record length
    INTEGER(KIND = 8) :: file_size
    !! File size
    REAL(KIND = DP) :: dummy
    !! Dummy variable 
    !
    !  the call to set_ndnmbr is just a trick to get quickly
    !  a file label by exploiting an existing subroutine
    !  (if you look at the sub you will find that the original 
    !  purpose was for pools and nodes)
    !   
    CALL set_ndnmbr(0, iq, 1, nqc, filelab)
    tempfile = TRIM(dvscf_dir) // TRIM(prefix) // '.dvscf_q' // filelab
    INQUIRE(IOLENGTH = unf_recl) dummy 
    unf_recl = unf_recl  * lrdrho
    mult_unit = unf_recl
    mult_unit = recn * mult_unit
    !
    !  open the dvscf file, read and close
    !
    OPEN(iudvscf, FILE = tempfile, FORM = 'unformatted', &
         ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl, STATUS = 'old')
    IF (ios /= 0) CALL errore('readdvscf', 'error opening ' // tempfile, iudvscf)
    !
    ! check that the binary file is long enough
    INQUIRE(FILE = tempfile, SIZE = file_size)
    IF (mult_unit > file_size) CALL errore('readdvscf', &
         TRIM(tempfile) //' too short, check ecut', iudvscf)
    !
    READ(iudvscf, REC = recn) dvscf
    CLOSE(iudvscf, STATUS = 'keep')
    !
    RETURN
    !
    !-------------------------------------------------------------
    END SUBROUTINE readdvscf
    !-------------------------------------------------------------
    ! 
    !------------------------------------------------------------
    SUBROUTINE readwfc(ipool, recn, evc0)
    !------------------------------------------------------------
    !!
    !! Open wfc files as direct access, read, and close again
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !
    USE kinds,    ONLY : DP
    USE io_files, ONLY : prefix, tmp_dir
    USE units_lr, ONLY : lrwfc, iuwfc
    USE wvfct,    ONLY : npwx
    USE pwcom,    ONLY : nbnd
    USE low_lvl,  ONLY : set_ndnmbr
    USE noncollin_module, ONLY : npol
    USE mp_global,        ONLY : nproc_pool, me_pool, npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: recn
    !! kpoint number
    INTEGER, INTENT(in) :: ipool
    !! poolfile number to be read (not used in serial case)
    COMPLEX(KIND = DP), INTENT(out) :: evc0(npwx * npol, nbnd)
    !! wavefunction is read from file
    !
    ! Local variables
    CHARACTER(LEN = 256) :: tempfile
    !! Temp file 
    CHARACTER(LEN = 3) :: nd_nmbr0
    !! File number
    INTEGER :: unf_recl
    !! Rcl unit
    INTEGER :: ios
    !! Error number
    REAL(KIND = DP) :: dummy
    !! Dummy variable 
    !
    ! Open the wfc file, read and close
    CALL set_ndnmbr(ipool, me_pool, nproc_pool, npool, nd_nmbr0)
    !
#if defined(__MPI)
    tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc' // nd_nmbr0
# else
    tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc'
#endif
    INQUIRE(IOLENGTH = unf_recl) dummy
    unf_recl = unf_recl * lrwfc
    !
    OPEN(iuwfc, FILE = tempfile, FORM = 'unformatted', ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl)
    IF (ios /= 0) CALL errore('readwfc', 'error opening wfc file', iuwfc)
    READ(iuwfc, REC = recn) evc0
    CLOSE(iuwfc, STATUS = 'keep')
    !
    RETURN
    !
    !------------------------------------------------------------
    END SUBROUTINE readwfc
    !------------------------------------------------------------
    ! 
    !--------------------------------------------------------------
    SUBROUTINE readgmap(nkstot, ngxx, ng0vec, g0vec_all_r, lower_bnd) 
    !--------------------------------------------------------------
    !!
    !! read map of G vectors G -> G-G_0 for a given q point
    !! (this is used for the folding of k+q into the first BZ) 
    !!
    ! 
    USE kinds,    ONLY : DP
    USE mp_global,ONLY : inter_pool_comm, world_comm
    USE mp,       ONLY : mp_bcast, mp_max
    use io_global,ONLY : meta_ionode, meta_ionode_id
    use io_var,   ONLY : iukgmap, iukmap
    use pwcom,    ONLY : nks
    use elph2,    ONLY : shift, gmap, igk_k_all, ngk_all
    USE io_files, ONLY : prefix
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nkstot
    !! Total number of k-points
    INTEGER, INTENT(out) :: ngxx
    !! Maximum number of G-vectors over all pools
    INTEGER, INTENT(out) :: ng0vec
    !! Number of G_0 vectors
    INTEGER, INTENT(in) :: lower_bnd
    !! Lower bound for the k-parallellization
    REAL(KIND = DP), INTENT(out) :: g0vec_all_r(3, 125)
    !! G_0 vectors needed to fold the k+q grid into the k grid, cartesian coord.
    !
    ! Lork variables
    INTEGER :: ik
    !! Counter on k-points 
    INTEGER :: ik1, itmp
    !! Temporary indeces when reading kmap and kgmap files
    INTEGER :: ig0
    !! Counter on G_0 vectors
    INTEGER :: ishift
    !! Counter on G_0 vectors
    INTEGER :: ig
    !! Counter on G vectors
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: tmp
    !! Temporary variable
    !
    ! OBSOLETE: now we read directly the igkq to get the proper ngxx
    ! read only a piece of the map to save time 
    ! the proper allocation bound would be ngxx = max(max(igkq))
    ! where the max is taken over the ig and the ik
    ! Here I use a simpler estimate: take the sphere npwx + two
    ! extra shells. This may not work for strange shapes of the
    ! reciproc latt. In this case just set ngxx = ngm_g
    !
    ! ngxx = NINT(4./3.*3.14*(2+(3.0/4.0/3.14*DBLE(npwx))**(1./3.))**3.)
    !
    ! Note that the k+q point below does not correspond to the actual (true) 
    ! k+q, but since we only need to take the max over k and k+q this
    ! does not matter
    !
    ngxx = 0
    DO ik = 1, nks
      !
      IF (MAXVAL(igk_k_all(1:ngk_all(ik + lower_bnd - 1), ik + lower_bnd - 1)) > ngxx) THEN
        ngxx = MAXVAL(igk_k_all(1:ngk_all(ik + lower_bnd - 1), ik + lower_bnd - 1))
      ENDIF
      !
    ENDDO
    !
#if defined(__MPI)
    tmp = DBLE(ngxx)
    CALL mp_max(tmp, inter_pool_comm)  
    ngxx = NINT(tmp)
#endif
    !
    IF (meta_ionode) THEN
      !
      OPEN(iukgmap, FILE = TRIM(prefix)//'.kgmap', FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
      IF (ios /=0) CALL errore('readgmap', 'error opening kgmap file', iukgmap)
      !
      DO ik = 1, nkstot
        READ(iukgmap, *) ik1, shift(ik1)
      ENDDO
      READ(iukgmap, *) ng0vec
      !
      !  the following seems crazy but I make it for compatibility
      !  with versions up to 2.1.5:
      !
      !  iukgmap has been created by ../PW/set_kplusq.f90 and has
      !  the correct gmap(), but the wrong shift() (actually the
      !  shift for a specific q-point)
      !
      !  since createkmap.f90 has regenerated the shifts for the
      !  present k-point I read them again in kmap.dat. The above 
      !  'fake' reading is because the gmap appears *after* the
      !  wrong kmap.
      !
      OPEN(iukmap, FILE = TRIM(prefix)//'.kmap', FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('readgmap', 'error opening kmap file', iukmap)
      DO ik = 1, nkstot
        READ(iukmap,*) ik1, itmp, shift(ik1)
      ENDDO
      CLOSE(iukmap) 
      !
    ENDIF
    !
    ! first node broadcasts ng0vec to all nodes for allocation of gmap
    !
    CALL mp_bcast(ng0vec, meta_ionode_id, world_comm)
    !
    ALLOCATE(gmap(ngxx * ng0vec), STAT = ierr)
    IF (ierr /= 0) CALL errore('readgmap', 'Error allocating gmap', 1)
    !
    IF (meta_ionode) THEN
       !
      DO ig0 = 1, ng0vec
        READ(iukgmap,*) g0vec_all_r(:,ig0)
      ENDDO
      DO ig = 1, ngxx
        ! 
        ! at variance with the nscf calculation, here gmap is read as a vector,
        ! 
        READ(iukgmap,*) (gmap(ng0vec * ( ig - 1 ) + ishift), ishift = 1, ng0vec)
      ENDDO
      !
      CLOSE(iukgmap)
      !
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    !
    CALL mp_bcast(g0vec_all_r, meta_ionode_id, world_comm)
    CALL mp_bcast(shift, meta_ionode_id, world_comm)
    CALL mp_bcast(gmap, meta_ionode_id, world_comm)
    !
    !--------------------------------------------------------------
    END SUBROUTINE readgmap
    !--------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io_epw
  !------------------------------------------------------------------------------
