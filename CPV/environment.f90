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
    

        SUBROUTINE environment_start( )

          USE io_global, ONLY: stdout, ionode
          USE mp_global, ONLY: mpime, nproc
          USE control_flags, ONLY: program_name
          USE parser, ONLY: int_to_char
          use para_mod, ONLY: me, node
          use mp, only: mp_env
          USE version

          LOGICAL    :: texst
          REAL(dbl)  :: elapsed_seconds
          EXTERNAL      elapsed_seconds
          INTEGER    :: nchar
          CHARACTER(LEN=80) :: uname
          CHARACTER(LEN=80) :: cp_version


          CALL init_clocks( .TRUE. )
          CALL start_clock( program_name )

          start_seconds = elapsed_seconds()

          cp_version = TRIM (version_number) // " - " // TRIM (version_date)

          ! ...  Temporary for para_mod

          me = mpime + 1

          ! ...  search for file CRASH and delete it

          IF( ionode ) THEN
            INQUIRE( FILE=TRIM(crash_file), EXIST=texst )
            IF( texst ) THEN
              OPEN(  UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD' )
              CLOSE( UNIT=crashunit, STATUS='DELETE' )
            END IF
          END IF

          ! ...       each processor other than me=1 opens its own standard output file,
          ! ...       this is mainly for debugging, usually only ionode writes to stdout

          IF( .NOT. ionode ) THEN

            uname = 'out.' // int_to_char( mpime )
            nchar = INDEX(uname,' ') - 1

            IF( program_name == 'CP90' ) THEN
              !
              ! useful for debugging purposes 
              !     open(6,file=file = uname(1:nchar),status='unknown')
              open( unit = stdout, file='/dev/null', status='unknown' )

            ELSE IF( program_name == 'FPMD' ) THEN

              OPEN( unit = stdout, file = uname(1:nchar), status = 'unknown' )

            END IF

          END IF

          ! ...  Temporary for para_mod

          if (me < 10) then
             write(node,'(i1,2x)') me
          else if (me < 100) then
             write(node,'(i2,1x)') me
          else if (me < 1000) then
             write(node,'(i3)') me
          else
             call errore('startup','wow, >1000 nodes !!',nproc)
          end if

          IF( program_name == 'CP90' ) THEN

            if ( ionode ) then
              WRITE( stdout,'(72("*"))')
              WRITE( stdout,'(4("*"),64x,4("*"))')
              WRITE( stdout,'(4("*"),"  CPV: variable-cell Car-Parrinello ", &
                   &  "molecular dynamics          ",4("*"))')
              WRITE( stdout,'(4("*"),"  using ultrasoft Vanderbilt ", &
                   &  "pseudopotentials - v.",a6,8x,4("*"))') version_number
              WRITE( stdout,'(4("*"),64x,4("*"))')
              WRITE( stdout,'(72("*"))')
              WRITE( stdout,'(/5x,''Parallel version (MPI)'')')
            end if

          ELSE IF( program_name == 'FPMD' ) THEN

            CALL opening_date_and_time( cp_version )

          END IF

          WRITE( stdout,100)  nproc, mpime
100       FORMAT(3X,'Tasks =',I5,'  This task id =',I5)

          RETURN
        END SUBROUTINE

!==-----------------------------------------------------------------------==!

        SUBROUTINE environment_end( )

          USE io_global, ONLY: stdout, ionode
          USE control_flags, ONLY: program_name

          REAL(dbl)  :: total_seconds

          REAL(dbl)  :: elapsed_seconds
          EXTERNAL      elapsed_seconds

          IF(ionode) THEN
            WRITE( stdout,*)
          END IF

          CALL print_clock( program_name )
          CALL stop_clock( program_name )

          CALL closing_date_and_time( )

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

        SUBROUTINE opening_date_and_time( cp_version )

          USE io_global, ONLY: stdout, ionode

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

        SUBROUTINE closing_date_and_time( )

          USE io_global, ONLY: stdout, ionode

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
