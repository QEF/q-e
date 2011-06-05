
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE read_pseudo_module_fpmd
!=----------------------------------------------------------------------------=!

        USE kinds
        USE io_files,     ONLY: pseudounit, psfile, pseudo_dir
        USE pseudo_types, ONLY: pseudo_upf
        USE pseudo_types, ONLY: nullify_pseudo_upf, deallocate_pseudo_upf
        USE uspp_param,   ONLY: upf

        IMPLICIT NONE

        SAVE

        PRIVATE

        REAL(DP) :: TOLMESH = 1.d-5

        PUBLIC :: readpp
        PUBLIC :: check_file_type

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!
INTEGER FUNCTION check_file_type( is )
  !
  ! ...   This subroutine guesses the pseudopotential type
  ! on return:
  ! -1   file is nonexistent
  !  0   file is unknown (guess: old CPV norm-conserving format) 
  !  1   file is *.vdb or *.van  Vanderbilt US pseudopotential
  !  2   file is *.RRKJ3         Andrea's   US new code 
  ! 11   file is NUMERIC (FPMD only) no more supported use UPF
  ! 12   file is ANALYTIC (FPMD only) no more supported use UPF
  ! 20   file is UPF
  !
  INTEGER, INTENT(IN) :: is
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=80) :: dummy
  LOGICAL, EXTERNAL :: matches
  INTEGER :: ios, info, l
  LOGICAL :: exst
  !
  info = 0  
  ios  = 0
  filename = TRIM( pseudo_dir ) // TRIM( psfile(is) )
  !
  INQUIRE ( FILE = TRIM(filename), EXIST=exst )
  IF ( .NOT. exst) THEN
     check_file_type = -1
     return
  END IF
  OPEN( UNIT = pseudounit, FILE = TRIM(filename), ACTION = 'READ', &
        STATUS = 'OLD' )
  header_loop: do while (ios == 0)
    read ( pseudounit, *, iostat = ios, err = 200) dummy  
    if (matches ("<PP_HEADER", dummy) ) then
      info = 20
      exit header_loop
    endif
  enddo header_loop
  200 continue

  IF( info == 0 ) THEN
    REWIND( pseudounit )
    READ ( pseudounit, *, iostat = ios, err = 300)
    dummy = ' ' 
    READ ( pseudounit, *, iostat = ios, err = 300) dummy  
    IF( matches( "NUMERIC", dummy ) ) THEN
      info = 11
    ELSE IF( matches( "ANALYTIC", dummy ) ) THEN
      info = 12
    END IF
  END IF
  300 continue

  CLOSE( pseudounit )

  IF( info == 0 ) THEN
    l = len_trim ( filename )
    if (filename (l - 3:l) .eq.'.vdb'.or.filename (l - 3:l) .eq.'.van') &
      info = 1
    if (l > 5) then
      if (filename (l - 5:l) .eq.'.RRKJ3') info = 2
    end if
  END IF

  check_file_type = info

  RETURN
END FUNCTION check_file_type

!=----------------------------------------------------------------------------=!

   SUBROUTINE readpp( xc_type )

     !  this subroutine reads pseudopotential parameters from file
     !
     !  See check_file_type for Allowed format
     !  
     !  
     !  ----------------------------------------------

      USE mp, ONLY: mp_bcast, mp_sum
      USE io_global, ONLY: stdout, ionode, ionode_id
      USE uspp, ONLY : okvan
      USE core, ONLY : nlcc_any
      USE uspp_param, ONLY : oldvan, nvb
      use ions_base, only: zv, nsp
      use upf_module, only: read_upf
      use read_uspp_module, only: readvan, readrrkj
      use funct, only: get_iexch, get_icorr, get_igcx, get_igcc, set_dft_from_name, dft_is_hybrid
      USE upf_to_internal, ONLY: set_pseudo_upf
      USE atom,            ONLY :  msh, rgrid
      use radial_grids, ONLY : deallocate_radial_grid, nullify_radial_grid

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: xc_type

! ... declare other variables
      CHARACTER(LEN=20)  :: dft_name
      CHARACTER(LEN=20)  :: pottyp
      CHARACTER(LEN=80)  :: error_msg
      CHARACTER(LEN=256) :: filename
      INTEGER            :: is, ierr, info
      INTEGER            :: iexch_, icorr_, igcx_, igcc_

!  end of declarations
!  ----------------------------------------------

      nvb     = 0  ! number of Vanderbilt pseudo
      !
      nlcc_any = .false. ! core corrections

      IF( nsp < 1 ) THEN
        CALL errore(' READPOT ',' nsp less than one! ', 1 )
      END IF

      IF( ALLOCATED( rgrid ) ) THEN
         DO is = 1, SIZE( rgrid )
            CALL deallocate_radial_grid( rgrid( is ) )
            CALL nullify_radial_grid( rgrid( is ) )
         END DO
         DEALLOCATE( rgrid )
         DEALLOCATE( msh )
      END IF

      ALLOCATE( rgrid( nsp ), msh(nsp ) )

      DO is = 1, nsp
         CALL nullify_radial_grid( rgrid( is ) )
      END DO

      IF( ALLOCATED( upf ) ) THEN
        DO is = 1, SIZE( upf )
          CALL deallocate_pseudo_upf( upf( is ) )
          CALL nullify_pseudo_upf( upf( is ) )
        END DO
        DEALLOCATE( upf )
      END IF

      ALLOCATE( upf( nsp ) )

      !  nullify upf objects as soon as they are instantiated

      DO is = 1, nsp
         CALL nullify_pseudo_upf( upf( is ) )
      END DO

      ierr = 0
      info = 0
      error_msg = 'none'
     
      IF( ionode ) THEN
        WRITE( stdout,4)
    4   FORMAT(//,3X,'Atomic Pseudopotentials Parameters',/, &
                  3X,'----------------------------------' )
      END IF

      DO is = 1, nsp

        filename = TRIM( pseudo_dir ) // TRIM( psfile(is) )
        !
        upf(is)%nlcc  = .FALSE.
        upf(is)%nbeta = 0
        upf(is)%tvanp = .FALSE.
        !
        IF( ionode ) THEN
          WRITE( stdout,6) is, TRIM(filename)
    6     FORMAT( /,3X,'Reading pseudopotential for specie # ',I2,' from file :',/,3X,A)
        END IF

        IF( ionode ) THEN
          info = check_file_type( is )
          SELECT CASE (info)
          CASE (0)
             WRITE( stdout,"(3X,'file type is ',I2,': Old CPV NC PP')") info
          CASE (1)
             WRITE( stdout,"(3X,'file type is ',I2,': Vanderbilt US PP')") info
          CASE (2)
             WRITE( stdout,"(3X,'file type is ',I2,': RRKJ3')") info
          CASE (11)
             WRITE( stdout,"(3X,'file type is ',I2,': Old FPMD Numeric')") info
          CASE (12)
             WRITE( stdout,"(3X,'file type is ',I2,': Old FPMD Analytic')") info
          CASE (20)
             WRITE( stdout, "(3X,'file type is ',I2,': UPF')") info
          END SELECT
        END IF
        CALL mp_bcast( info, ionode_id )
        IF (info == -1) CALL errore ('readpp', &
                            'file '//TRIM(filename)//' not found',is)

        !  Now each processor read the pseudopotential file
  
        ierr = 0

        OPEN( UNIT = pseudounit, FILE = TRIM(filename), ACTION = 'READ', &
              STATUS = 'OLD' )
        !
        ! used only by obsolete Vanderbilt format with Herman-Skillman grid
        !
        oldvan(is)  = .false.
        !
        IF( info == 20 ) THEN
           !
           !  ...      Pseudopotential form is UPF
           !
           call read_upf(upf(is), rgrid(is), ierr, unit=pseudounit)
           !
           IF ( ierr /= 0 ) THEN
             CALL deallocate_pseudo_upf( upf(is) )
           ELSE
             call set_pseudo_upf( is, upf( is ) )
           END IF

        ELSE IF( info == 1 ) THEN

           CALL readvan( pseudounit, is, upf(is) )
           CALL set_pseudo_upf( is, upf( is ), rgrid( is ) )

        ELSE IF( info == 2 ) THEN

           CALL readrrkj( pseudounit, is, upf(is) )
           CALL set_pseudo_upf( is, upf( is ), rgrid( is ) )

        ELSE IF( info == 11 ) THEN

          error_msg = ' type no more supported, convert to UPF using fpmd2upf '
          ierr = info

        ELSE IF( info == 12 ) THEN

          error_msg = ' type no more supported, convert to UPF using fpmd2upf '
          ierr = info

        ELSE IF( info == 0 ) THEN

            CALL errore(' readpp ', ' file format no longer supported ', 2 )

        END IF

        CLOSE( pseudounit )

        CALL mp_sum( ierr )
        IF( ierr /= 0 ) THEN
          CALL errore(' readpseudo ', error_msg, ABS(ierr) )
        END IF

        ! ... Zv = valence charge of the (pseudo-)atom, read from PP files,
        ! ... is set equal to Zp = pseudo-charge of the pseudopotential
        !     (should be moved out from here)
 
        zv(is) = upf(is)%zp
        !
        !     Ultrasoft formats: UPF, AdC, Vanderbilt ("old" and new)
        !     norm-conserving formats: UPF
        !
        !     check on input ordering: US first, NC later 
        !
        if(is > 1) then
          if ( (.NOT. upf(is-1)%tvanp) .AND. upf(is)%tvanp ) then
             call errore ('readpp', 'ultrasoft PPs must precede norm-conserving',is)
          endif
        endif
        !
        !     count u-s vanderbilt species 
        !
        if (upf(is)%tvanp) nvb=nvb+1
        !
        !     check for core corrections
        !
        nlcc_any = nlcc_any .OR. upf(is)%nlcc
        !
        if ( xc_type /= 'none' ) then
          ! 
          !  DFT xc functional, given from input
          !
          dft_name = TRIM( xc_type )
          CALL set_dft_from_name( dft_name )

          WRITE( stdout, fmt="(/,3X,'Warning XC functionals forced to be: ',A)" ) dft_name
          !
        else
          !
          ! check for consistency of DFT
          !
          if (is == 1) then
            iexch_ = get_iexch()
            icorr_ = get_icorr()
            igcx_ =  get_igcx()
            igcc_ =  get_igcc()
          else
            if ( iexch_ /= get_iexch() .or. icorr_ /= get_icorr() .or. &
                 igcx_  /= get_igcx()  .or. igcc_ /= get_igcc() ) then
               CALL errore( 'readpp','inconsistent DFT read',is)
            end if
          end if
        end if
 
        IF ( dft_is_hybrid() ) &
            CALL errore( 'readpp', 'HYBRID XC not implemented in CPV', 1 )

      END DO

      okvan = ( nvb > 0 )
      !
      RETURN
      END SUBROUTINE readpp

!=----------------------------------------------------------------------------=!
   END MODULE read_pseudo_module_fpmd
!=----------------------------------------------------------------------------=!
!
