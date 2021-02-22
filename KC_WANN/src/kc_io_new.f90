!----------------------------------------------------------------------------
MODULE io_kcwann
  !----------------------------------------------------------------------------
  !
  !! this module contains some common subroutines used to read and write
  !! the data produced by KC_WANN. 
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, iunpun
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
  !
  IMPLICIT NONE
  PRIVATE
  !
  LOGICAL,       SAVE      :: rho_binary = .TRUE.
  !
  PUBLIC :: rho_binary
  PUBLIC :: write_rhowann, read_rhowann
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rhowann( file_base, rho, fft_desc, ionode, inter_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho in real-space, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE mp,        ONLY : mp_get, mp_sum, mp_rank, mp_size
#if defined __HDF5
      USE hdf5_qe,  ONLY  : write_rho_hdf5, h5fclose_f, &
                            prepare_for_writing_final, add_attributes_hdf5, rho_hdf5_write  
#endif
      USE fft_types
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN) :: file_base
      COMPLEX(DP),       INTENT(IN) :: rho(:)
      TYPE(fft_type_descriptor),INTENT(IN) :: fft_desc
      INTEGER,           INTENT(IN) :: inter_group_comm
      LOGICAL,           INTENT(IN) :: ionode
      !
      INTEGER               :: nr1,nr2,nr3, nr1x, nr2x,nr3x
      INTEGER               :: rhounit, ierr, i, j, jj, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: rho_file
      CHARACTER(LEN=256)    :: rho_file_hdf5
      CHARACTER(LEN=10)     :: rho_extension
      COMPLEX(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: my_group_id, me_group, me_group2, me_group3, &
                                            nproc_group, nproc_group2, nproc_group3, &
                                            io_group_id, io_group2, io_group3
      INTEGER,  EXTERNAL    :: find_free_unit
      !

      my_group_id = mp_rank( inter_group_comm )

      me_group = fft_desc%mype ; me_group2 = fft_desc%mype2 ; me_group3 = fft_desc%mype3
      nproc_group = fft_desc%nproc ; nproc_group2 = fft_desc%nproc2 ; nproc_group3 = fft_desc%nproc3
      !
      nr1  = fft_desc%nr1  ; nr2  = fft_desc%nr2  ; nr3  = fft_desc%nr3
      nr1x = fft_desc%nr1x ; nr2x = fft_desc%nr2x ; nr3x = fft_desc%nr3x
      !
      rho_extension = '.dat'
      IF ( .NOT. rho_binary ) rho_extension = '.xml'
      !
      rho_file = TRIM( file_base ) // TRIM( rho_extension )
      rhounit = find_free_unit ()
      !
      IF ( ionode ) THEN 
#if defined  __HDF5
         rho_file_hdf5 = TRIM( file_base ) // '.hdf5'
         CALL prepare_for_writing_final(rho_hdf5_write, 0 ,rho_file_hdf5)
         CALL add_attributes_hdf5(rho_hdf5_write,nr1,"nr1")
         CALL add_attributes_hdf5(rho_hdf5_write,nr2,"nr2")
         CALL add_attributes_hdf5(rho_hdf5_write,nr3,"nr3")
#else
         IF (rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='unformatted')
         IF (.NOT. rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='formatted')
         CALL errore( 'write_rhowann', 'cannot open ' // TRIM( rho_file ) // ' file for writing', ierr )
#endif
      END IF 
      !
#if !defined __HDF5
      IF ( ionode ) THEN
         !
         IF (rho_binary) THEN 
           WRITE(rhounit) nr1, nr2, nr3
         ELSE
           WRITE(rhounit, '("mesh ", 3i10)') nr1, nr2, nr3
         ENDIF
         !
      END IF
#endif
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the group (pool) that will write rho
      !
      io_group_id = 0
      !
      IF ( ionode ) io_group_id = my_group_id
      !
      CALL mp_sum( io_group_id, fft_desc%comm )
      CALL mp_sum( io_group_id, inter_group_comm ) ! io_group_id is the (pool) group that contains the ionode
      !
      ! ... find the index of the ionode within Y and Z  groups
      !
      io_group2 = 0 ; IF ( ionode ) io_group2 = me_group2
      CALL mp_sum( io_group2, fft_desc%comm )  ! io_group2 is the group index of the ionode in the Y group (nproc2)
      io_group3 = 0 ; IF ( ionode ) io_group3 = me_group3
      CALL mp_sum( io_group3, fft_desc%comm )  ! io_group3 is the group index of the ionode in the Z group (nproc3)
      !
      ! ... find out the owner of each "z" plane
      !
      DO ip = 1, nproc_group3
         !
         kowner( (fft_desc%i0r3p(ip)+1):(fft_desc%i0r3p(ip)+fft_desc%nr3p(ip)) ) = ip - 1
         !
      END DO
      !
      ldr = nr1x*fft_desc%my_nr2p
      !
      IF ( ( my_group_id == io_group_id ) ) THEN ! only the group of ionode collects and writes the data
         !
         DO k = 1, nr3
            !
            !  Only one subgroup write the charge density
            ! 
            rho_plane = (0.D0, 0.D0)
            IF( ( kowner(k) == me_group3 ) ) THEN
               !
               kk = k - fft_desc%my_i0r3p
               ! 
               DO jj = 1, fft_desc%my_nr2p
                  !
                  j = jj + fft_desc%my_i0r2p
                  DO i = 1, nr1
                     !
                     rho_plane(i+(j-1)*nr1) = rho(i+(jj-1)*nr1x+(kk-1)*ldr)
                     !
                  END DO
                  !
               END DO
               call mp_sum(rho_plane, fft_desc%comm2 ) ! collect the data over the Y group (nproc2)
               !
            END IF
            !
            ! if this processor is in the same comm3 group as ionode (me_group2==io_group2) 
            IF ( kowner(k) /= io_group3 .and. me_group2==io_group2) & 
               CALL mp_get( rho_plane, rho_plane, me_group3, io_group3, kowner(k), k, fft_desc%comm3 )
            !
            IF ( ionode ) THEN
#if defined __HDF5
              CALL write_rho_hdf5(rho_hdf5_write,k,rho_plane)
#else
              !
              IF (rho_binary) THEN 
                 WRITE(rhounit) k
                 WRITE(rhounit) rho_plane
              ELSE
                 WRITE(rhounit, '("z= ", i5)') k
                 WRITE(rhounit, '(E23.16)') rho_plane
              ENDIF
              !
#endif
            ENDIF
            !
         END DO
         !
      END IF
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
#if defined __HDF5
         CALL h5fclose_f(rho_hdf5_write%file_id,ierr)
#else
         CLOSE (rhounit, STATUS='keep')
#endif       
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rhowann
    !
   !------------------------------------------------------------------------
   SUBROUTINE read_rhowann( rho_file_base, fft_desc, rho )
     !------------------------------------------------------------------------
     !
     ! ... Reads charge density rho, one plane at a time, to avoid 
     ! ... collecting the entire charge density on a single processor
     !
     USE io_global, ONLY : ionode, ionode_id
     USE mp_images, ONLY : intra_image_comm
     USE mp,        ONLY : mp_put, mp_sum, mp_rank, mp_size
#if defined __HDF5
     USE hdf5_qe,   ONLY : read_rho_hdf5, read_attributes_hdf5, &
          prepare_for_reading_final, h5fclose_f, rho_hdf5_write, hdf5_type
#endif
     USE fft_types
     USE io_files,  ONLY : check_file_exist
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=*),  INTENT(IN)  :: rho_file_base
     TYPE(fft_type_descriptor),INTENT(IN) :: fft_desc
     COMPLEX(DP),          INTENT(OUT) :: rho(:)
     !
     INTEGER               :: rhounit, ierr, i, j, jj, k, kk, ldr, ip, k_
     INTEGER               :: nr( 3 ), nr1_, nr2_, nr3_, nr1, nr2, nr3, nr1x, nr2x, nr3x
     INTEGER               :: me_group, me_group2, me_group3, &
                              nproc_group, nproc_group2, nproc_group3
     CHARACTER(LEN=256)    :: rho_file
     CHARACTER(LEN=256)    :: rho_file_hdf5
     COMPLEX(DP), ALLOCATABLE :: rho_plane(:)
     INTEGER,  ALLOCATABLE :: kowner(:)
     LOGICAL               :: exst
     INTEGER,  EXTERNAL    :: find_free_unit
#if defined(__HDF5)
     TYPE(hdf5_type),ALLOCATABLE   :: h5desc
#endif
     CHARACTER(LEN=256) :: string
     CHARACTER(LEN=10)     :: rho_extension
     !
     me_group = fft_desc%mype ; me_group2 = fft_desc%mype2 ; me_group3 = fft_desc%mype3
     nproc_group = fft_desc%nproc ; nproc_group2 = fft_desc%nproc2 ; nproc_group3 = fft_desc%nproc3
     !
     nr1  = fft_desc%nr1  ; nr2  = fft_desc%nr2  ; nr3  = fft_desc%nr3
     nr1x = fft_desc%nr1x ; nr2x = fft_desc%nr2x ; nr3x = fft_desc%nr3x
     !
#if defined(__HDF5)
     rho_file_hdf5 = TRIM( rho_file_base ) // '.hdf5'
     exst = check_file_exist(TRIM(rho_file_hdf5))
     IF ( .NOT. exst ) CALL errore ('read_rhowann', 'searching for '// TRIM(rho_file_hdf5),10)
#else 
     rho_extension = '.dat'
     IF ( .NOT. rho_binary ) rho_extension = '.xml'
     rhounit = find_free_unit ( )
     rho_file = TRIM( rho_file_base ) // TRIM( rho_extension )
     exst = check_file_exist( TRIM(rho_file) ) 
     !
     IF ( .NOT. exst ) CALL errore('read_rhowann', 'searching for '//TRIM(rho_file), 10)
#endif
     !
     IF ( ionode ) THEN
#if defined (__HDF5)
        ALLOCATE ( h5desc)
        CALL prepare_for_reading_final(h5desc, 0 ,rho_file_hdf5)
        CALL read_attributes_hdf5(h5desc, nr1_,"nr1")
        CALL read_attributes_hdf5(h5desc, nr2_,"nr2")
        CALL read_attributes_hdf5(h5desc, nr3_,"nr3")
        nr = [nr1_,nr2_,nr3_]
#else
        IF (rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='unformatted', STATUS='old')
        IF (.NOT. rho_binary) OPEN (rhounit, FILE=rho_file, IOSTAT=ierr, FORM='formatted', STATUS='old')
        CALL errore( 'read_rho_xml', 'cannot open ' // TRIM( rho_file ) // ' file for reading', ierr )
        IF (rho_binary) READ(rhounit) nr(1), nr(2), nr(3)
        IF (.NOT. rho_binary) READ(rhounit,*) string,nr(1), nr(2), nr(3)
#endif
        !
        IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
           CALL errore( 'read_rhowann', 'dimensions do not match', 1 )
        !
     END IF
     !
     ALLOCATE( rho_plane( nr1*nr2 ) )
     ALLOCATE( kowner( nr3 ) )
     !
     DO ip = 1, nproc_group3
        !
        kowner( (fft_desc%i0r3p(ip)+1):(fft_desc%i0r3p(ip)+fft_desc%nr3p(ip)) ) = ip - 1
        !
     END DO
     !
     ldr = nr1x*fft_desc%my_nr2p
     !
     ! ... explicit initialization to zero is needed because the physical
     ! ... dimensions of rho may exceed the true size of the FFT grid 
     !
     rho(:) = (0.0_DP, 0.0_DP)
     !
     DO k = 1, nr3
        !
        ! ... only ionode reads the charge planes
        !
        IF ( ionode ) THEN
#if defined __HDF5
           CALL  read_rho_hdf5(h5desc , k,rho_plane)
#else
           IF (rho_binary) THEN 
              READ(rhounit) k_
              READ(rhounit) rho_plane
           ELSE
              READ(rhounit,*) string, k_
              READ(rhounit, '(E23.16)') rho_plane
           ENDIF
#endif
        ENDIF
        !
        ! ... planes are sent to the destination processor (all processors in this image)
        !
        CALL mp_bcast( rho_plane, ionode_id, intra_image_comm )
        !
        IF( kowner(k) == me_group3 ) THEN
           !
           kk = k - fft_desc%my_i0r3p
           DO jj = 1, fft_desc%my_nr2p
              j = jj + fft_desc%my_i0r2p
              DO i = 1, nr1
                 rho(i+(jj-1)*nr1x+(kk-1)*ldr) = rho_plane(i+(j-1)*nr1)
              END DO
           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( rho_plane )
     DEALLOCATE( kowner )
     !
     IF ( ionode ) THEN
        !
#if defined __HDF5
        CALL h5fclose_f(h5desc%file_id,ierr)
        DEALLOCATE ( h5desc)
#else
        CLOSE (rhounit)
#endif    
     END IF
     !
     RETURN
     !
   END SUBROUTINE read_rhowann



END MODULE io_kcwann
