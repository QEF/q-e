!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!==-----------------------------------------------------------------------==!
      MODULE environment
!==-----------------------------------------------------------------------==!

        USE kinds
        USE io_files, ONLY: crash_file, crashunit, &
                            stop_file, stopunit

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(dbl)          :: start_seconds

        PUBLIC :: environment_start
        PUBLIC :: environment_end
        PUBLIC :: opening_date_and_time, closing_date_and_time
      
!==-----------------------------------------------------------------------==!
      CONTAINS
!==-----------------------------------------------------------------------==!
    

        SUBROUTINE environment_start( ionode, mpime, nproc, cp_version )

          USE io_global, ONLY: stdout
          USE parser, ONLY: int_to_char

          LOGICAL, INTENT(IN) :: ionode
          INTEGER, INTENT(IN) :: mpime, nproc
          CHARACTER(LEN=*), INTENT(IN) :: cp_version

          LOGICAL    :: texst
          REAL(dbl)  :: elapsed_seconds
          EXTERNAL      elapsed_seconds
          INTEGER    :: nchar
          CHARACTER(LEN=80) :: uname


          CALL init_clocks( .TRUE. )
          CALL start_clock( 'FPMD' )

          start_seconds = elapsed_seconds()

          ! ...  search for file CRASH and delete it

          IF( ionode ) THEN
            INQUIRE( FILE=TRIM(crash_file), EXIST=texst )
            IF( texst ) THEN
              OPEN(  UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD' )
              CLOSE( UNIT=crashunit, STATUS='DELETE' )
            END IF
          END IF

! ...       use IBM SP Hardware performance monitor      

! ...       each processor other than me=1 opens its own standard output file,
! ...       this is mainly for debugging, usually only ionode writes to stdout

          IF( .NOT. ionode ) THEN

            uname = 'fort_6.' // int_to_char( mpime )
            nchar = INDEX(uname,' ') - 1
            OPEN( unit = stdout, file = uname(1:nchar), status = 'unknown', form = 'formatted')

          END IF

          ! ...  turn on stream buffering (Cray only)

#if defined __STREAMS_BUF
          CALL set_d_stream(1)
#endif

          CALL opening_date_and_time(ionode, cp_version)

          WRITE( stdout,100)  nproc, mpime
100       FORMAT(3X,'Tasks =',I5,'  This task id =',I5)

          RETURN
        END SUBROUTINE

!==-----------------------------------------------------------------------==!

        SUBROUTINE environment_end( ionode, mpime )

          USE io_global, ONLY: stdout

          LOGICAL, INTENT(IN) :: ionode
          INTEGER, INTENT(IN) :: mpime

          REAL(dbl)  :: total_seconds

          REAL(dbl)  :: elapsed_seconds
          EXTERNAL      elapsed_seconds

! ...     turn off stream buffering (Cray only)

#if defined __STREAMS_BUF
          CALL set_d_stream(0)                                                             
#endif

          IF(ionode) THEN
            WRITE( stdout,*)
          END IF

          CALL print_clock( 'FPMD' )
          CALL stop_clock( 'FPMD' )

          CALL closing_date_and_time(ionode)

          total_seconds = elapsed_seconds() - start_seconds

          IF(ionode) THEN
            WRITE( stdout,'(A,F7.1)') '   ELAPSED SECONDS: ', total_seconds
            WRITE( stdout,'(A)')      '   JOB DONE.'
            WRITE( stdout,3335)
          END IF
 3335     FORMAT('=',78('-'),'=')

          RETURN
        END SUBROUTINE

!==-----------------------------------------------------------------------==!

        SUBROUTINE opening_date_and_time( ionode, cp_version )

          USE io_global, ONLY: stdout

          LOGICAL, INTENT(IN) :: ionode
          CHARACTER(LEN=*), INTENT(IN) :: cp_version
          CHARACTER(LEN=9) :: cdate, ctime

          CALL date_and_tim(cdate, ctime)

! ...     write program heading
          IF(ionode) THEN
            WRITE( stdout,3333) cp_version
            WRITE( stdout,3334) 'THIS RUN WAS STARTED ON:  ' // ctime // ' ' // cdate
          END IF

 3333     FORMAT('=',78('-'),'=',/,&
      18X,'AB-INITIO COSTANT PRESSURE MOLECULAR DYNAMICS',/,&
      25X,'A Car-Parrinello Parallel Code',/,&
          '=',78('-'),'=',//, &
          '  Version: ',A60,/,&
          '  Authors: Carlo Cavazzoni, Guido Chiarotti, Sandro Scandolo,',/,&
          '    Paolo Focher, Gerardo Ballabio',//, &
          '=',78('-'),'=',/)

 3334     FORMAT(3X,A60,/)
          RETURN
        END SUBROUTINE

!==-----------------------------------------------------------------------==!

        SUBROUTINE closing_date_and_time( ionode )

          USE io_global, ONLY: stdout

          LOGICAL, INTENT(IN) :: ionode
          CHARACTER(LEN=9) :: cdate, ctime

          CALL date_and_tim(cdate, ctime)

          IF( ionode ) THEN
            WRITE( stdout,*)
            WRITE( stdout,3334) 'THIS RUN WAS TERMINATED ON:  ' // ctime // ' ' // cdate
            WRITE( stdout,3335)
          END IF

 3334     FORMAT(3X,A60,/)
 3335     FORMAT('=',78('-'),'=')

          RETURN
        END SUBROUTINE

!==-----------------------------------------------------------------------==!
      END MODULE environment
!==-----------------------------------------------------------------------==!
