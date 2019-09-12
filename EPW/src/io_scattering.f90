  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !                                                                            
  !----------------------------------------------------------------------
  MODULE io_scattering
  !----------------------------------------------------------------------
  !! 
  !! This module contains various printing routines
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
    USE io_epw,        ONLY : iufilFi_all
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
    INTEGER :: itemp
    !! Temperature index
    ! 
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
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
    USE io_epw,    ONLY : iufilFi_all
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
    INTEGER :: itemp
    !! Temperature index
    ! 
    CHARACTER(LEN = 256) :: name1
 
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
      ! 
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
    USE io_epw,           ONLY : iunepmat_merge, iunepmat, iunepmatcb_merge,&
                                 iunepmatcb, iunsparseq_merge, iunsparsek_merge,iunsparsei_merge, &
                                 iunsparsej_merge,iunsparset_merge, iunepmatcb_merge,&
                                 iunsparseqcb_merge, iunsparsekcb_merge,&
                                 iunsparseicb_merge, iunsparsejcb_merge, iunsparsetcb_merge
    USE mp_global,        ONLY : my_pool_id, npool, world_comm 
    USE io_files,         ONLY : tmp_dir, prefix 
    USE mp,               ONLY : mp_sum, mp_barrier
    USE io_global,        ONLY : stdout
    USE elph2,            ONLY : lrepmatw2_merge, lrepmatw5_merge
    USE epwcom,           ONLY : int_mob, carrier, ncarrier
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE,MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_INTEGER
#endif
    !
    IMPLICIT NONE
    !
    INTEGER :: i2, i4, i5, i6
    !! Indexes to loop over file sizes
    INTEGER :: ipool
    !! Process index
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! Pool number, character
    REAL(KIND = DP), ALLOCATABLE :: trans_prob(:)
    !! Variable for reading and writing trans_prob
    REAL(KIND = DP), ALLOCATABLE :: trans_probcb(:)
    !! Variable for reading and writing trans_prob
    INTEGER :: lrepmatw2_tot(npool)
    !! Lenght of each file
    INTEGER :: lrepmatw5_tot(npool)
    !! Lenght of each file
    CHARACTER(LEN = 256) :: dirname(2)
    !! Name of the directory to hold files
    CHARACTER(LEN = 256) :: filename(6)
    !! Name of the files to merge files
    CHARACTER(LEN = 256) :: path_to_files(2)
    !! Name of the path to files
    INTEGER :: ich
    !! Loop over directories
    INTEGER, ALLOCATABLE :: sparse(:, :)
    !! Vaariable for reading and writing the files
    INTEGER, ALLOCATABLE :: sparsecb(:, :)
    !! Vaariable for reading and writing the files
    INTEGER :: ierr
    !! Error variable for MPI
    INTEGER :: ifil
    !! Index over the files
    INTEGER :: io_u(6)
    !! Input output units
#if defined(__MPI)
    INTEGER (KIND = MPI_OFFSET_KIND) :: lsize
    !! Size of what we write
    INTEGER (KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset while writing scattering to files
    !
    IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 0.0))) THEN
      !
      ALLOCATE(trans_prob(lrepmatw2_merge)) !There may be a problem if lrepmatw2_merge == 0
      ALLOCATE(sparse(5, lrepmatw2_merge))
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
      path_to_files(1)='./'//ADJUSTL(TRIM(dirname(1)))//'/'//TRIM(prefix)//'.epmatkq1'//'_'
      path_to_files(2)='./'//ADJUSTL(TRIM(dirname(2)))//'/'//'sparse'//'_'
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
        WRITE(my_pool_id_ch,"(I0)") my_pool_id
        filint = TRIM(path_to_files(ich))//TRIM(my_pool_id_ch)
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
      DEALLOCATE(trans_prob)
      DEALLOCATE(sparse)
      !
    ENDIF
    IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 0.0))) THEN
      !
      ALLOCATE(trans_probcb(lrepmatw5_merge))
      ALLOCATE(sparsecb(5, lrepmatw5_merge))
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
      DEALLOCATE(trans_probcb)
      DEALLOCATE(sparsecb)
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
    ! This SUBROUTINE opens all the files needed to save scattering rates for the IBTE.
    ! 
    USE kinds,            ONLY : DP
    USE io_files,         ONLY : tmp_dir, prefix, create_directory, delete_if_present
    USE io_epw,           ONLY : iunepmat, iunsparseq, iunsparsek, &
                                 iunsparsei, iunsparsej, iunsparset, iunsparseqcb, &
                                 iunsparsekcb, iunsparseicb, iunsparsejcb,&
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
    INQUIRE(FILE = 'restart_ibte.fmt',EXIST = exst)
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
        filint = './'//ADJUSTL(TRIM(dirname(1)))//'/'//TRIM(prefix)//'.epmatkq1'//'_'//TRIM(my_pool_id_ch)
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
        filint = './'//ADJUSTL(TRIM(dirname(2)))//'/'//'sparse'//'_'//TRIM(my_pool_id_ch)
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
        filint = './'//ADJUSTL(TRIM(dirnamecb(1)))//'/'//TRIM(prefix)//'.epmatkqcb1'//'_'//TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        ! 
        IF (exst2) THEN
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw5_restart(my_pool_id + 1) - 1) * 8 + 1  
            READ(iunepmatcb, POS=position_byte) dummy_real
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart_ibte.fmt is present but not the Fepmatkqcb1 folder', 1)
        ENDIF
        ! 
        filint = './'//ADJUSTL(TRIM(dirnamecb(2)))//'/'//'sparsecb'//'_'//TRIM(my_pool_id_ch)
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
        filint = './'//ADJUSTL(TRIM(dirname(1)))//'/'//TRIM(prefix)//'.epmatkq1'//'_'//TRIM(my_pool_id_ch)
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
        filint = './'//ADJUSTL(TRIM(dirname(2)))//'/'//'sparse'//'_'//TRIM(my_pool_id_ch)
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
        filint = './'//ADJUSTL(TRIM(dirnamecb(1)))//'/'//TRIM(prefix)//'.epmatkqcb1'//'_'//TRIM(my_pool_id_ch)
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
        filint = './'//ADJUSTL(TRIM(dirnamecb(2)))//'/'//'sparsecb'//'_'//TRIM(my_pool_id_ch)
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
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_epw,    ONLY : iufilscatt_rate
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
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file. 
    CHARACTER(LEN = 256) :: dummy1
    !! Dummy variable to store the text of the scattering_rate file 
    ! 
    WRITE(stdout,'(/5x,"Reading scattering rate from file"/)')
    !
    IF (mpime == ionode_id) THEN
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1,'(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1,'(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1,'(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate,FILE = name1, status='old',iostat=ios)
      WRITE(stdout,'(a16,a22)') '     Open file: ',name1   
      ! There are two comment line at the beginning of the file
      READ(iufilscatt_rate,*) dummy1
      READ(iufilscatt_rate,*) dummy1
      !
      DO ik = 1, nktotf
        !
        DO ibnd = 1, nbndfst
          !
          READ(iufilscatt_rate,*) ik_tmp, ibnd_tmp, &
               etf_all(ibndmin - 1 + ibnd, ik), inv_tau_all(1, ibnd, ik) 
          inv_tau_all(1, ibnd, ik) = inv_tau_all(1, ibnd, ik) / (ryd2mev * meV2invps) 
       
          !
          ! Check that the file corresponds to the run we are making
          IF (ABS(ibnd_tmp - ibndmin - ibnd +1 ) > 0)  CALL errore('io_scattering', &
            'Band read from the scattering_rate file do not match current calculation ', 1)
          ! 
        ENDDO
        ! Check that the file corresponds to the run we are making
        IF (ABS(ik_tmp - ik) > 0)  CALL errore('io_scattering', &
          'k-point read from the scattering_rate file do not match current calculation ', 1)
        !
      ENDDO
      !
      etf_all = etf_all / ryd2ev
      etf_all = etf_all + ef0
      !
      CLOSE(iufilscatt_rate)
      ! 
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
    !
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : ibndmax, ibndmin, lower_bnd, upper_bnd, nbndfst
    USE io_epw,    ONLY : iufilsigma_all
    USE io_files,  ONLY : diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    ! Local variable
    LOGICAL :: exst
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
      ! 
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
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : ibndmax, ibndmin, lower_bnd, upper_bnd, nbndfst
    USE io_epw,    ONLY : iufilsigma_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE constants_epw, ONLY :  zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    ! Local variable
    LOGICAL :: exst
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
        IF (nqtotf_read /= totq) CALL errore('io_scattering',&
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
    USE io_epw,    ONLY : iufiltau_all
    USE io_files,  ONLY : diropn
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE constants_epw, ONLY : zero
    !
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
          DO ibnd = 1, (nbndfst)
            i = i + 1
            aux(i) = inv_tau_all(itemp, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, (nbndfst)
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
        ! 
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, (nbndfst)
              i = i + 1
              aux(i) = inv_tau_allcb(itemp, ibnd, ik)
            ENDDO
          ENDDO
        ENDDO
        !
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, (nbndfst)
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
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, meta_ionode_id
    USE elph2,     ONLY : ibndmax, ibndmin, inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_epw,    ONLY : iufiltau_all
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
        IF (nqtotf_read /= totq) CALL errore('io_scattering',&
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
          IF (nqtotf_read /= totq) CALL errore('io_scattering',&
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
    !----------------------------------------------------------------------------
    SUBROUTINE merge_read(nktotf, nqtotf_new, inv_tau_all_new)
    !----------------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : ibndmax, ibndmin, nbndfst
    USE io_epw,    ONLY : iufiltau_all
    USE io_files,  ONLY : tmp_dir, diropn
    USE epwcom,    ONLY : nstemp, restart_filq
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    ! Local variable
    LOGICAL :: exst
    !
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    INTEGER, INTENT(OUT) :: nqtotf_new
    !! Total number of q-points
    REAL(KIND = DP), INTENT(inout) :: inv_tau_all_new(nstemp, nbndfst, nktotf)
    !! Scattering rate read from file restart_filq
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name of the file 
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
    !! 
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
        INQUIRE (IOLENGTH = unf_recl) dummy  
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
      CALL errore('rwepmatw','iop not permitted', 1)
      !
    ENDIF
    !
    !----------------------------------------------------------------------
    END SUBROUTINE rwepmatw
    !----------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io_scattering
  !------------------------------------------------------------------------------
