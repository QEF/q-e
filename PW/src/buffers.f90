!
!        ~~~ BUffer Input/Output Library. ~~~
! Copyright Lorenzo Paulatto <paulatz@gmail.com> 2013
!
! Contains a few changes by PG wrt the original implementation:
! - data is complex, not real
! - most routines are functions that return error status instead of stopping
! - added possibility to store file name info in the linked list
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt 
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE buiol
  USE kinds, ONLY : DP
  !
  PUBLIC :: init_buiol          ! init the linked chain of i/o units
  PUBLIC :: stop_buiol          ! destroy the linked chain, dealloc everything
  PUBLIC :: report_buiol        ! report on total number of units and memory usage
  PUBLIC :: buiol_open_unit     ! (unit, recl, ext, dir) open a new unit
  PUBLIC :: buiol_close_unit    ! (unit) close the unit, dealloc the space
  PUBLIC :: buiol_check_unit    ! (unit) returns recl, if opened, -1 if closed
  PUBLIC :: buiol_get_ext       ! (unit) returns file extension
  PUBLIC :: buiol_get_dir       ! (unit) returns dir where file is opened
  PUBLIC :: buiol_report_unit   ! (unit, mem?) report about unit status (on stdout)
  PUBLIC :: buiol_write_record  ! (unit, recl, nrec, DATA) write DATA(recl) in record nrec of unit
  PUBLIC :: buiol_read_record   ! (unit, recl, nrec, DATA) read DATA(recl) from record nrec of unit
  !
  PRIVATE
  ! initial number of records in the buffer (each record will only be allocated on write!)
  INTEGER,PARAMETER :: nrec0 = 1024
  ! when writing beyond the last available record increase the index by AT LEAST this factor..
  REAL(DP),PARAMETER :: fact0 = 1.5_dp
  ! .. furthermore, allocate up to AT LEAST this factor times the required overflowing nrec
  REAL(DP),PARAMETER :: fact1 = 1.2_dp
  ! NOTE: the new buffer size will be determined with both methods, taking the MAX of the two
  !
  ! Size of the single item of the record (for memory usage report only)
  INTEGER,PARAMETER :: size0 = DP ! 8 bytes
  !
  ! base element of the linked chain of buffers
  TYPE index_of_list
    TYPE(data_in_the_list),POINTER :: index(:)
    INTEGER :: nrec, unit, recl
    CHARACTER(LEN=256) :: extension, save_dir
    TYPE(index_of_list),POINTER :: next => null()
  END TYPE
  !
  ! sub-structure containing the data buffer
  TYPE data_in_the_list 
    COMPLEX(DP), POINTER :: data(:) => null()
  END TYPE
  !
  ! beginning of the linked chain, statically allocated (for implementation simplicity)
  TYPE(index_of_list),SAVE,POINTER :: ENTRY => null()
  !
  ! set to true when the library has been initialized
  LOGICAL,SAVE :: is_init_buiol = .false.
  !
  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  !
  SUBROUTINE init_buiol
    IMPLICIT NONE
    ! avoid initializing twice, or we will loose the head of the list!
    IF (is_init_buiol) THEN 
#if defined(__DEBUG)
       CALL infomsg('buiol', 'already initialized')
#endif
       RETURN
    ENDIF
    !
    ALLOCATE(ENTRY)
    ALLOCATE(ENTRY%index(0))
    ENTRY%nrec =  0
    ENTRY%unit = -1
    ENTRY%recl = -1
    ENTRY%extension= ' '
    ENTRY%save_dir = ' '
    NULLIFY(ENTRY%next)
    is_init_buiol = .true.
    !
    RETURN
  END SUBROUTINE init_buiol
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE stop_buiol
    IMPLICIT NONE
    TYPE(index_of_list),POINTER :: CURSOR, AUX
    IF (.not.is_init_buiol) RETURN
    IF (.not.associated(ENTRY) ) CALL errore('stop_buiol', 'ENTRY was lost.',1)
    !
    CURSOR => ENTRY
    DO WHILE (associated(CURSOR%NEXT))
      AUX => CURSOR
      CURSOR => CURSOR%NEXT
      CALL dealloc_buffer(AUX)
    ENDDO
    CALL dealloc_buffer(CURSOR)
    !
    is_init_buiol=.false.
    RETURN
  END SUBROUTINE stop_buiol
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE report_buiol
    IMPLICIT NONE
    TYPE(index_of_list),POINTER :: CURSOR
    INTEGER :: mem
    !
    IF (.not.is_init_buiol) THEN
      WRITE(*,'(2x,a,3i14)') "[BUIOL] not even initialized"
      RETURN
    ENDIF
    !
    WRITE(*,'(2x,106("-") )')
    mem = 0
    CURSOR => ENTRY
    DO WHILE (associated(CURSOR%NEXT))
      CALL buiol_report_buffer(CURSOR, mem)
      CURSOR => CURSOR%NEXT
    ENDDO
    CALL buiol_report_buffer(CURSOR, mem)
    WRITE(*,'(2x,106("-"))')
    WRITE(*,'(2x,a,3i14)') "[BUIOL] total memory used B/KB/MB", mem, mem/1024, mem/1024**2
    WRITE(*,'(2x,106("-"))')

    RETURN
  END SUBROUTINE report_buiol
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_open_unit(unit, recl, extension, save_dir) RESULT (ierr)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit, recl
    CHARACTER(LEN=*), INTENT(in) :: extension, save_dir
    INTEGER :: ierr
    TYPE(index_of_list),POINTER :: CURSOR
    !
    IF (.not.is_init_buiol) CALL errore('buiol_open_unit', 'You must init before open',1)
    IF(recl<0) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_open_unit', 'wrong recl')
#endif
       ierr = 1
       RETURN
    END IF
    !
    ! check if the unit is already opened
    CURSOR => find_unit(unit)
    IF(associated(CURSOR)) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_open_unit', 'unit already opened')
#endif
       ierr = -1
       RETURN
    END IF
    !
    ! all is fine, allocate a new unit with standard size
    CURSOR => alloc_buffer(unit, recl, nrec0, extension, save_dir)
    !
    ! place it at the beginning of the chain
    CURSOR%next => ENTRY%next
    ENTRY%next  => CURSOR
    ierr = 0
    !
    RETURN
    !
  END FUNCTION buiol_open_unit
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_close_unit(unit) RESULT (ierr)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    INTEGER :: ierr
    TYPE(index_of_list),POINTER :: CURSOR, AUX
    !
    ! find the unit to close
    CURSOR => find_prev_unit(unit)
    IF(.not.associated(CURSOR))  THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_close_unit', 'cannot close this unit')
#endif
       ierr = 1
    END IF
    IF(.not.associated(CURSOR%next)) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_close_unit', 'cannot find unit to close',1)
#endif
       ierr = 2
    END IF
    !
    ! replace this unit with the next, but keep track of it
    AUX => CURSOR%next
    CURSOR%next => AUX%next ! <--- works even if %next is null()
    !
    ! destroy the closed unit
    CALL dealloc_buffer(AUX)
    ierr = 0
    !
    RETURN
    !
  END FUNCTION buiol_close_unit
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_check_unit(unit) RESULT(recl)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    INTEGER :: recl
    TYPE(index_of_list),POINTER :: CURSOR
    !
    ! find the unit
    CURSOR => find_unit(unit)
    IF(.not.associated(CURSOR)) THEN
      recl = -1
    ELSE
      recl = CURSOR%recl
    ENDIF
    !
    RETURN
    !
  END FUNCTION buiol_check_unit
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_get_ext(unit) RESULT(extension)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    CHARACTER(LEN=256) :: extension
    TYPE(index_of_list),POINTER :: CURSOR
    !
    ! find the unit
    CURSOR => find_unit(unit)
    IF(.not.associated(CURSOR)) THEN
      extension = ' '
    ELSE
      extension = CURSOR%extension
    ENDIF
    !
    RETURN
    !
  END FUNCTION buiol_get_ext
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_get_dir(unit) RESULT(save_dir)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    CHARACTER(LEN=256) :: save_dir
    TYPE(index_of_list),POINTER :: CURSOR
    !
    ! find the unit
    CURSOR => find_unit(unit)
    IF(.not.associated(CURSOR)) THEN
      save_dir = ' '
    ELSE
      save_dir = CURSOR%save_dir
    ENDIF
    !
    RETURN
    !
  END FUNCTION buiol_get_dir
  ! \/o\______\\_______________________________________/^>
  SUBROUTINE increase_nrec(nrec_new, CURSOR)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nrec_new
    TYPE(index_of_list),POINTER,INTENT(inout) :: CURSOR
    !
    INTEGER :: i
    TYPE(data_in_the_list),POINTER :: new(:), old(:)
    !
    IF(nrec_new < CURSOR%nrec) CALL errore('increase_nrec', 'wrong new nrec',1)
    !
    ! create a new index with more space
    ALLOCATE(new(nrec_new))
    !
    ! associate the data to the new unit
    old => CURSOR%index
    DO i = 1, CURSOR%nrec
      new(i)%data => old(i)%data ! <-- also the null() are copied
    ENDDO
    CURSOR%index => new
    !
    ! clean the old index
    CURSOR%nrec = nrec_new
    DEALLOCATE(old)
    !
    RETURN
    !
  END SUBROUTINE increase_nrec
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_write_record(unit, recl, nrec, DATA) RESULT (ierr)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit, recl, nrec
    COMPLEX(dp),INTENT(in) :: DATA(recl)
    INTEGER :: ierr
    !
    TYPE(index_of_list),POINTER :: CURSOR
    INTEGER :: nrec_new
    !
    ! find the unit, if it exists
    CURSOR => find_unit(unit)
    IF(.not.associated(CURSOR)) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_write_record', 'cannot write: unit not opened')
#endif
       ierr = 1
       RETURN
    END IF
    IF(CURSOR%recl/=recl) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_write_record', 'cannot write: wrong recl')
#endif
       ierr = 2
       RETURN
    END IF
    !
    ! increase size of index, if necessary
    IF(CURSOR%nrec<nrec) THEN
      nrec_new = NINT(MAX(fact0*DBLE(CURSOR%nrec),fact1*DBLE(nrec)))
      CALL increase_nrec(nrec_new, CURSOR ) 
    ENDIF
    !
    IF(.not.associated(CURSOR%index(nrec)%data)) &
      ALLOCATE( CURSOR%index(nrec)%data(recl) )
    !
    ! copy the data
    CURSOR%index(nrec)%data = DATA
    ierr = 0
    RETURN
    !
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  FUNCTION buiol_read_record(unit, recl, nrec, DATA) RESULT (ierr)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit, recl, nrec
    COMPLEX(dp),INTENT(out) :: DATA(recl)
    INTEGER :: ierr
    !
    TYPE(index_of_list),POINTER :: CURSOR
    !
    ! sanity checks
    CURSOR => find_unit(unit)
    IF(.not.associated(CURSOR)) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_read_record', 'cannot read: unit not opened')
#endif
       ierr = 1
       RETURN
    END IF
    IF(CURSOR%recl/=recl) THEN
#if defined(__DEBUG)
        CALL infomsg('buiol_read_record', 'cannot read: wrong recl')
#endif
       ierr = 1
       RETURN
    END IF
    IF(CURSOR%nrec<nrec) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_read_record', 'cannot read: wrong nrec')
#endif
       ierr =-1
       RETURN
    END IF
    IF(.not.associated(CURSOR%index(nrec)%data)) THEN
#if defined(__DEBUG)
       CALL infomsg('buiol_read_record', 'cannot read: virgin nrec')
#endif
       ierr =-1
       RETURN
    END IF
    !
    DATA = CURSOR%index(nrec)%data
    ierr = 0
    RETURN
    !
  END FUNCTION buiol_read_record
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE buiol_report_unit(unit)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    !
    TYPE(index_of_list),POINTER :: CURSOR
    ! sanity checks
    CURSOR => find_unit(unit)
#if defined(__DEBUG)
    IF(.not.associated(CURSOR)) CALL errore('buiol_report_unit', 'cannot report: unit not opened',1)
#endif
    CALL buiol_report_buffer(CURSOR)
    RETURN
    !
  END SUBROUTINE buiol_report_unit
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE buiol_report_buffer(CURSOR, mem)
    IMPLICIT NONE
    TYPE(index_of_list),INTENT(in) :: CURSOR
    INTEGER,OPTIONAL,INTENT(inout) :: mem
    !
    INTEGER :: i, ndata, bytes
    !
    ndata = 0
    DO i = 1,CURSOR%nrec
      IF(associated(CURSOR%index(i)%data)) ndata=ndata+1
    ENDDO
    !
    bytes = ndata*CURSOR%recl*size0
    WRITE(*,'(2x,a,2(a,i8),(a,2i8),(a,i12))') "[BUIOL] ", &
             "unit:", CURSOR%unit, &
        "   | recl:", CURSOR%recl, &
        "   | nrec (idx/alloc):", CURSOR%nrec, ndata, &
        "   | memory used:", bytes
    IF(present(mem)) mem = mem+bytes
    RETURN
    !
  END SUBROUTINE buiol_report_buffer
  ! \/o\________\\\_________________________________________/^>
  FUNCTION find_unit(unit) RESULT(CURSOR)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    TYPE(index_of_list),POINTER :: CURSOR
    !
    IF (.not.is_init_buiol) CALL errore('find_unit', 'You must init before find_unit',1)
    !
    CURSOR => ENTRY
    DO WHILE (associated(CURSOR%NEXT))
      CURSOR => CURSOR%NEXT
      IF(CURSOR%unit == unit) RETURN ! <-- found
    ENDDO
    CURSOR => null() ! <------------------ not found 
    RETURN
  END FUNCTION find_unit
  ! \/o\________\\\_________________________________________/^>
  FUNCTION find_prev_unit(unit) RESULT(CURSOR)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit
    TYPE(index_of_list),POINTER :: CURSOR
    !
    IF (.not.is_init_buiol) CALL errore('find_prev_unit', 'You must init before find_prev_unit',1)
    !
    CURSOR => ENTRY
    DO WHILE (associated(CURSOR%NEXT))
      IF(CURSOR%next%unit == unit) RETURN ! <-- found
      CURSOR => CURSOR%NEXT
    ENDDO
    CURSOR => null() ! <------------------ not found 
    RETURN
  END FUNCTION find_prev_unit
  ! \/o\________\\\_________________________________________/^>
  FUNCTION alloc_buffer(unit, recl, nrec, extension, save_dir)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: unit, recl, nrec
    CHARACTER(LEN=*), INTENT(in) :: extension, save_dir
    TYPE(index_of_list),POINTER :: alloc_buffer
    TYPE(index_of_list),POINTER :: CURSOR
    !
    ALLOCATE(CURSOR)
    CURSOR%unit = unit
    CURSOR%recl = recl
    CURSOR%nrec = nrec0
    CURSOR%extension = extension
    CURSOR%save_dir  = save_dir
    NULLIFY(CURSOR%next)
    ALLOCATE(CURSOR%index(CURSOR%nrec))
    !
    alloc_buffer => CURSOR
    RETURN
  END FUNCTION alloc_buffer
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE dealloc_buffer(CURSOR)
    IMPLICIT NONE
    TYPE(index_of_list),POINTER,INTENT(inout) :: CURSOR
    !
    INTEGER :: i
    DO i = 1,CURSOR%nrec
      IF(associated(CURSOR%index(i)%data)) THEN
         DEALLOCATE(CURSOR%index(i)%data)
         NULLIFY(CURSOR%index(i)%data)
      ENDIF
    ENDDO
    DEALLOCATE(CURSOR%index)
    CURSOR%unit = -1
    CURSOR%recl = -1
    CURSOR%nrec = -1
    DEALLOCATE(CURSOR)
    NULLIFY(CURSOR)
    !
  END SUBROUTINE dealloc_buffer
  ! \/o\________\\\_________________________________________/^>
END MODULE buiol
! <<^V^\\=========================================//-//-//========//O\\//

Module buffers

  use kinds, only: dp
  use buiol, only: init_buiol, buiol_open_unit, buiol_close_unit, &
                   buiol_check_unit, buiol_get_ext, buiol_get_dir, &
                   buiol_read_record, buiol_write_record
  implicit none
  !
  ! QE interfaces to BUIOL module
  !
  PUBLIC :: open_buffer, get_buffer, save_buffer, close_buffer
  !
  PRIVATE
  INTEGER:: nunits = 0
  !
contains

  !----------------------------------------------------------------------------
  SUBROUTINE open_buffer (unit, extension, nword, io_level, exst, exst_file, direc)
    !---------------------------------------------------------------------------
    !
    !   io_level>0: connect unit "unit" to file "wfc_dir"/"prefix"."extension"
    !   (or "direc"/"prefix"."extension" if optional variable direc specified)
    !   for direct I/O access, with record length = nword complex numbers;
    !   on output, exst=T(F) if the file (does not) exists
    !
    !   io_level=0: open a buffer for storing records of length nword complex
    !   numbers; store in memory file-related variables for later usage.
    !   on output, exst=T(F) if the buffer is already allocated
    !
    !   on output, optional variable exst_file=T(F) if file is present (absent)
    !
    USE io_files,  ONLY : diropn, wfc_dir
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: extension
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: direc
    INTEGER, INTENT(IN) :: unit, nword, io_level
    LOGICAL, INTENT(OUT) :: exst
    LOGICAL, INTENT(OUT), OPTIONAL :: exst_file
    CHARACTER(LEN=256) :: save_dir
    !
    INTEGER :: ierr
    !
    !   not-so-elegant way to initialize the linked chain with units
    !
    IF ( nunits == 0 ) CALL init_buiol( )
    !
    IF (extension == ' ') &
       CALL errore ('open_buffer','filename extension not given',1)
    !
    IF (present(direc)) THEN
       save_dir=TRIM(direc)
    ELSE
       save_dir=TRIM(wfc_dir)
    ENDIF
    !
    IF ( io_level <= 0 ) THEN
       CALL diropn ( unit, extension, -1, exst, save_dir )      
       IF (present(exst_file)) exst_file=exst
       ierr = buiol_open_unit ( unit, nword, extension, save_dir )
       IF ( ierr > 0 ) CALL errore ('open_buffer', ' cannot open unit', 2)
       exst = ( ierr == -1 )
       IF (exst) THEN
          CALL infomsg ('open_buffer', 'unit already opened')
          nunits = nunits - 1
       END IF
    ELSE
       CALL diropn ( unit, extension, 2*nword, exst, save_dir )      
       IF (present(exst_file)) exst_file=exst
    ENDIF
    nunits = nunits + 1
    !
    RETURN
    !
  END SUBROUTINE open_buffer
  !----------------------------------------------------------------------------
  SUBROUTINE save_buffer( vect, nword, unit, nrec )
    !---------------------------------------------------------------------------
    !
    ! ... copy vect(1:nword) into the "nrec"-th record of a previously
    ! ... allocated buffer / opened direct-access file, depending upon
    ! ... how "open_buffer" was called
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nword, unit, nrec
    COMPLEX(DP), INTENT(IN) :: vect(nword)
    INTEGER :: ierr
    !
    ierr = buiol_check_unit (unit)
    IF( ierr > 0 ) THEN
       ierr = buiol_write_record ( unit, nword, nrec, vect )
       if ( ierr > 0 ) &
           CALL errore ('save_buffer', 'cannot write record', unit)
#if defined(__DEBUG)
       print *, 'save_buffer: record', nrec, ' written to unit', unit
#endif
    ELSE 
       CALL davcio ( vect, 2*nword, unit, nrec, +1 )
    END IF
    !
  END SUBROUTINE save_buffer
  !
  !----------------------------------------------------------------------------
  SUBROUTINE get_buffer( vect, nword, unit, nrec )
    !---------------------------------------------------------------------------
    !
    ! ... copy vect(1:nword) from the "nrec"-th record of a previously
    ! ... allocated buffer / opened direct-access file, depending upon
    ! ... how "open_buffer" was called. If buffer access was chosen 
    ! ... but buffer is not allocated, open the file, read from file
    !
    USE io_files, ONLY : diropn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nword, unit, nrec
    COMPLEX(DP), INTENT(OUT) :: vect(nword)
    CHARACTER(LEN=256) :: extension, save_dir
    INTEGER :: ierr
    LOGICAL :: opnd
    !
    ierr = buiol_check_unit (unit)
    IF( ierr > 0 ) THEN
       ierr = buiol_read_record ( unit, nword, nrec, vect )
#if defined(__DEBUG)
       print *, 'get_buffer: record', nrec, ' read from unit', unit
#endif
       if ( ierr < 0 ) then
          ! record not found: open file if not opened, read from it...
          INQUIRE( UNIT = unit, OPENED = opnd )
          IF ( .NOT. opnd ) THEN
             extension = buiol_get_ext (unit)
             save_dir  = buiol_get_dir (unit)
             CALL diropn ( unit, extension, 2*nword, opnd, save_dir )      
          END IF
          CALL davcio ( vect, 2*nword, unit, nrec, -1 )
          ! ... and save to memory
          ierr =  buiol_write_record ( unit, nword, nrec, vect )
          if ( ierr /= 0 ) CALL errore ('get_buffer', &
                                  'cannot store record in memory', unit)
#if defined(__DEBUG)
          print *, 'get_buffer: record', nrec, ' read from file', unit
#endif
       end if
#if defined(__DEBUG)
       print *, 'get_buffer: record', nrec, ' read from unit', unit
#endif
    ELSE
       CALL davcio ( vect, 2*nword, unit, nrec, -1 )
    END IF
    !
  END SUBROUTINE get_buffer

  SUBROUTINE close_buffer ( unit, status )
    !
    !     close unit with status "status" ('keep' or 'delete')
    !     deallocate related buffer if any; if "status='keep'"
    !     save it to file (opening it if not already opened).
    !     Does not complain if closing an already closed unit
    !
    USE io_files, ONLY : diropn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=*), INTENT(IN) :: status
    !
    COMPLEX(dp), ALLOCATABLE :: vect(:)
    CHARACTER(LEN=256) :: extension, save_dir
    INTEGER :: n, ierr, nrec, nword
    LOGICAL :: opnd
    !
    nword = buiol_check_unit (unit)
    !
    IF( nword > 0 ) THEN
       ! data is in memory buffer
       IF ( status == 'keep' .or. status == 'KEEP' ) then
          ! open file if not previously opened
          INQUIRE( UNIT = unit, OPENED = opnd )
          IF ( .NOT. opnd ) THEN
             extension = buiol_get_ext (unit)
             save_dir  = buiol_get_dir (unit)
             CALL diropn ( unit, extension, 2*nword, opnd, save_dir )      
          END IF
          allocate (vect(nword))
          n = 1
  10      continue
             ierr = buiol_read_record ( unit, nword, n, vect )
             IF ( ierr /= 0 ) go to 20
             CALL davcio ( vect, 2*nword, unit, n, +1 )
             n = n+1
          go to 10
  20      deallocate (vect)
       end if
       ierr = buiol_close_unit ( unit )
       if ( ierr < 0 ) &
            CALL errore ('close_buffer', 'error closing', ABS(unit))
#if defined(__DEBUG)
       print *, 'close_buffer: unit ',unit, 'closed'
#endif
    END IF
    INQUIRE( UNIT = unit, OPENED = opnd )
    IF ( opnd ) CLOSE( UNIT = unit, STATUS = status )
    nunits = nunits - 1
    !
  END SUBROUTINE close_buffer

  ! end interface for old "buffers" module

end module buffers
