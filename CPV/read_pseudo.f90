
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
   MODULE read_pseudo_module_fpmd
!=----------------------------------------------------------------------------=!

!  this module handles the reading of pseudopotential data
!  ----------------------------------------------


! ...   declare modules

        USE kinds, ONLY: DP
        USE io_files, ONLY: pseudounit
        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE pseudo_types, ONLY: nullify_pseudo_upf, deallocate_pseudo_upf

        IMPLICIT NONE

        SAVE

        PRIVATE

        REAL(DP) :: TOLMESH = 1.d-5
        INTEGER   :: nspnl = 0  ! number of non local species

        TYPE (pseudo_ncpp), ALLOCATABLE, TARGET :: ap(:)
        TYPE (pseudo_upf),  ALLOCATABLE, TARGET :: upf(:)

        PUBLIC :: ap, upf, nspnl, readpp
        PUBLIC :: pseudo_filename, check_file_type

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

!
!  -----------
!  subroutines
!  -----------
!

CHARACTER(LEN=256) FUNCTION pseudo_filename( is )
  USE io_files, ONLY: psfile, pseudo_dir
  INTEGER, INTENT(IN) :: is
  IF (TRIM(pseudo_dir) == ' ' ) then
     pseudo_filename=TRIM(psfile(is))
  ELSE
     pseudo_filename=TRIM(pseudo_dir)//TRIM(psfile(is))
  END IF
  RETURN
END FUNCTION pseudo_filename

!=----------------------------------------------------------------------------=!

INTEGER FUNCTION check_file_type( is )
  !
  ! ...   This subroutine guesses the pseudopotential type
  ! on return:
  ! -1   file is nonexistent
  !  0   file is unknown (guess: old CPV norm-conserving format) 
  !  1   file is *.vdb or *.van  Vanderbilt US pseudopotential
  !  2   file is *.RRKJ3         Andrea's   US new code 
  ! 11   file is NUMERIC (FPMD only)
  ! 12   file is ANALYTIC (FPMD only)
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
  filename = pseudo_filename( is )
  !
  INQUIRE ( FILE = TRIM(filename), EXIST=exst )
  IF ( .NOT. exst) THEN
     check_file_type = -1
     return
  END IF
  OPEN( UNIT = pseudounit, FILE = TRIM(filename), STATUS = 'OLD' )
  header_loop: do while (ios == 0)
    read ( pseudounit, *, iostat = ios, err = 200) dummy  
    if (matches ("<PP_HEADER>", dummy) ) then
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

SUBROUTINE check_types_order( )
  USE ions_base, ONLY: nsp
  IMPLICIT NONE
  INTEGER :: is, il
  LOGICAL :: tvanp
  !
  !   With Vanderbilt, only UPF are allowed
  !
  IF( ANY( upf(1:nsp)%tvanp  ) ) THEN
    CALL errore( &
           ' check_types_order ', ' vanderbilt pseudo, not yet implemented in FPMD ', 1 )
  END IF
  !
  !   non-local species must be ahead the local one,
  !
  il = 0
  DO is = 1, nsp
    IF ( ap(is)%nbeta == 0 ) THEN
      il = 1
    ELSE IF ( il == 1 ) THEN
      CALL errore( &
           ' check_types_order ', ' Local pseudopotentials should follow non local ones ', 1 )
    END IF
  END DO
  RETURN
END SUBROUTINE check_types_order

!=----------------------------------------------------------------------------=!

REAL(DP) FUNCTION calculate_dx( a, m )
  USE constants, ONLY: eps14
  REAL(DP), INTENT(IN) :: a(:)
  INTEGER, INTENT(IN) :: m 
  INTEGER :: n, nn
  REAL(DP) :: ra, rb 
  n  = MIN( SIZE( a ), m )
  nn = n
  IF( a(1) < eps14 ) THEN
     ra = a(2)
     nn = n - 1
  ELSE
     ra = a(1)
  END IF
  rb = a(n)
  calculate_dx = LOG( rb / ra ) / DBLE( nn - 1 )
  RETURN
END FUNCTION calculate_dx

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
      USE uspp_param, ONLY : zp, tvanp, oldvan
      USE atom, ONLY: numeric, nlcc, oc, lchi, nchi
      USE cvan, ONLY: nvb
      use ions_base, only: zv, nsp
      use read_upf_module, only: read_pseudo_upf
      use read_uspp_module, only: readvan, readrrkj
      use control_flags, only: program_name
      use funct, only: get_iexch, get_icorr, get_igcx, get_igcc, set_dft_from_name, dft_is_hybrid
      USE upf_to_internal, ONLY: set_pseudo_upf

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

      nspnl   = 0  ! number of non local pseudo
      nvb     = 0  ! number of Vanderbilt pseudo
      !
      oc      = 0  ! init atomic wf occupation
      lchi    = 0  ! init atomic wf angular momentum
      nchi    = 0  ! init numbero of atomic wf

      IF( nsp < 1 ) THEN
        CALL errore(' READPOT ',' nsp less than one! ', 1 )
      END IF

      IF( ALLOCATED( ap  ) ) DEALLOCATE( ap )
      IF( ALLOCATED( upf ) ) THEN
        DO is = 1, SIZE( upf )
          CALL deallocate_pseudo_upf( upf( is ) )
          CALL nullify_pseudo_upf( upf( is ) )
        END DO
        DEALLOCATE( upf )
      END IF

      ALLOCATE( ap( nsp )  )
      ALLOCATE( upf( nsp ) )

      ierr = 0
      info = 0
      error_msg = 'none'
     
      IF( ionode ) THEN
        WRITE( stdout,4)
    4   FORMAT(//,3X,'Atomic Pseudopotentials Parameters',/, &
                  3X,'----------------------------------' )
      END IF

      DO is = 1, nsp

        filename = TRIM( pseudo_filename( is ) )
        !
        CALL nullify_pseudo_upf( upf( is ) )
        !
        ap(is)%tnlcc  = .FALSE.
        ap(is)%nbeta    = 0
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

        OPEN( UNIT = pseudounit, FILE = filename, STATUS = 'OLD' )

        numeric(is) = .true.
        !
        ! used only by obsolete "bhs" format of CP
        !
        oldvan(is)  = .false.
        !
        ! used only by obsolete Vanderbilt format with Herman-Skillman grid
        !
        IF( info == 20 ) THEN
           !
           !  ...      Pseudopotential form is UPF
           !
           ap(is)%pottyp = 'UPF'
           !
           call read_pseudo_upf(pseudounit, upf(is), ierr)
           !
           IF ( ierr /= 0 ) THEN
             CALL deallocate_pseudo_upf( upf(is) )
           ELSE
             call set_pseudo_upf( is, upf( is ) )
             !
             IF( .NOT. upf(is)%tvanp ) THEN
               CALL upf2ncpp( upf(is), ap(is) )
             END IF
             !
           END IF

        ELSE IF( info == 1 ) THEN

           call readvan( is, pseudounit )

        ELSE IF( info == 2 ) THEN

           call readrrkj( is, pseudounit )

        ELSE IF( info == 11 ) THEN

          CALL read_head_pp( pseudounit, ap(is), error_msg, ierr)
          CALL read_numeric_pp( pseudounit, ap(is), error_msg, ierr)
          CALL ncpp2internal ( ap(is), is, xc_type, ierr )

        ELSE IF( info == 12 ) THEN

          CALL read_head_pp( pseudounit, ap(is), error_msg, ierr)
          CALL read_analytic_pp( pseudounit, ap(is), error_msg, ierr)
          CALL ncpp2internal ( ap(is), is, xc_type, ierr )

        ELSE IF( info == 0 ) THEN

          IF( program_name == 'FPMD' ) THEN
            CALL errore(' readpp ', ' file format not supported ', 1 )
          ELSE
            call readbhs(is,pseudounit)
          END IF

        END IF

        CLOSE( pseudounit )

        ! ... Zv = valence charge of the (pseudo-)atom, read from PP files,
        ! ... is set equal to Zp = pseudo-charge of the pseudopotential
 
        zv(is) = zp(is)

        CALL mp_sum( ierr )
        IF( ierr /= 0 ) THEN
          CALL errore(' readpseudo ', error_msg, ABS(ierr) )
        END IF
      
        IF( program_name == 'FPMD' ) THEN
          !
          IF( ap(is)%nbeta > 0 ) nspnl = nspnl + 1
          IF( upf(is)%tvanp  ) nvb   = nvb + 1
          IF( ionode ) THEN
            CALL ap_info( ap(is) )
          END IF
          !
        ELSE IF( program_name == 'CP90' ) THEN
          !
          !     Ultrasoft formats: UPF, AdC, Vanderbilt ("old" and new)
          !     norm-conserving formats: hsc, bhs, UPF
          !
          !     check on input ordering: US first, NC later 
          !
          if(is > 1) then
            if ( (.NOT. tvanp(is-1)) .AND. tvanp(is) ) then
               call errore ('readpp', &
                            'ultrasoft PPs must precede norm-conserving',is)
            endif
          endif
          !
          !     count u-s vanderbilt species 
          !
          if (tvanp(is)) nvb=nvb+1
          !
        END IF

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

      IF( program_name == 'FPMD' ) THEN
        CALL check_types_order()
      END IF

      okvan = ( nvb > 0 )

      RETURN
      END SUBROUTINE readpp

!=----------------------------------------------------------------------------=!

      SUBROUTINE analytic_to_numeric(ap)

!       This subroutine converts an Analytic pseudo into a numeric one

        USE constants, ONLY: pi
        USE pseudo_types, ONLY: pseudo_ncpp

        IMPLICIT NONE

        TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
        INTEGER   :: ir, mesh, lmax, l, n, il, ib, ll
        REAL(DP) :: xmin, zmesh, dx, x

! ...   declare external function
        REAL(DP) :: erf, erfc
        EXTERNAL erf, erfc

        IF( ap%mesh == 0 ) THEN
! ...     Local pseudopotential, define a logaritmic grid
          mesh  =  400
          xmin  = -5.0d0
          zmesh =  6.0d0
          dx    =  0.025d0
          DO ir = 1, mesh
            x = xmin + DBLE(ir-1) * dx
            ap%rw(ir)  = EXP(x) / zmesh
          END DO
          ap%mesh = mesh
          ap%dx   = dx
          ap%rab  = ap%dx * ap%rw
        END IF

        ap%vnl  = 0.0d0
        ap%vloc = 0.0d0
        ap%vrps = 0.0d0
        do l = 1, 3
          do ir = 1, ap%mesh
            ap%vnl(ir,l)= - ( ap%wrc(1) * erf( SQRT( ap%rc(1) ) * ap%rw(ir) ) + &
                              ap%wrc(2) * erf( SQRT( ap%rc(2) ) * ap%rw(ir) )   & 
                            ) * ap%zv / ap%rw(ir)
          end do
          do ir = 1, ap%mesh
            do n = 1, ap%igau
              ap%vnl(ir,l) = ap%vnl(ir,l) + &
                             ( ap%al(n,l) + ap%bl(n,l) * ap%rw(ir)**2 ) * &
                             EXP( - ap%rcl(n,l) * ap%rw(ir)**2 )
            end do
          end do
        end do

! ...   Copy local component to a separate array
        ap%vloc( : ) = ap%vnl( :, ( ap%lloc + 1 ) )
        DO l = 1, ap%nbeta
          ll = ap%lll(l) + 1  ! find out the angular momentum (ll-1) of the  
                              ! component stored in position l
          ap%vrps( :, l ) = ( ap%vnl( :, ll ) - ap%vloc( : ) ) * ap%rps( :, ll )
        END DO

        RETURN
      END SUBROUTINE analytic_to_numeric

!=----------------------------------------------------------------------------=!

      SUBROUTINE ap_info( ap )
        USE pseudo_types, ONLY: pseudo_ncpp
        USE io_global, ONLY: stdout

        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        INTEGER   :: in1, in2, in3, in4, m, il, ib, l, i

        WRITE( stdout, * ) 
        IF (ap%nbeta > 0) THEN
          WRITE( stdout,10) ap%pottyp
          IF (ap%tmix) THEN
            WRITE( stdout,107) 
            WRITE( stdout,106)  (ap%lll(l),l=1,ap%nbeta)
            WRITE( stdout,105)  (ap%wgv(l),l=1,ap%nbeta)
          ELSE
            WRITE( stdout,50 )   ap%lloc   
          END IF
          WRITE( stdout,60) (ap%lll(l),l=1,ap%nbeta)
        ELSE
! ...     A local pseudopotential has been read.
          WRITE( stdout,11) ap%pottyp
          WRITE( stdout,50) ap%lloc 
        END IF
        IF( ap%tnlcc ) THEN
          WRITE( stdout,12)
        END IF

   10   FORMAT(   3X,'Type is ',A10,' and NONLOCAL. ')
  107   FORMAT(   3X,'Mixed reference potential:')
  106   FORMAT(   3X,'  L     :',3(9X,i1))
  105   FORMAT(   3X,'  Weight:',3(2X,F8.5))
   50   FORMAT(   3X,'Local component is ..... : ',I3)
   60   FORMAT(   3X,'Non local components are : ',4I3)
   11   FORMAT(   3X,'Type is ',A10,' and LOCAL. ')
   12   FORMAT(   3X,'Using non local core corcorrections for this pseudo')
   20   FORMAT(   3X,'Pseudo charge : ',F8.3)

        WRITE( stdout,20) ap%zv

        IF( ap%pottyp /= 'ANALYTIC' ) THEN

          WRITE( stdout,131) ap%nchan, ap%mesh, ap%dx
          in1=1
          in2=ap%mesh/4
          in3=ap%mesh/2
          in4=ap%mesh
          WRITE( stdout,132)
          WRITE( stdout,120) in1,ap%rw(in1),ap%vloc(in1),(ap%vrps(in1,m),m=1,ap%nbeta)
          WRITE( stdout,120) in2,ap%rw(in2),ap%vloc(in2),(ap%vrps(in2,m),m=1,ap%nbeta)
          WRITE( stdout,120) in3,ap%rw(in3),ap%vloc(in3),(ap%vrps(in3,m),m=1,ap%nbeta)
          WRITE( stdout,120) in4,ap%rw(in4),ap%vloc(in4),(ap%vrps(in4,m),m=1,ap%nbeta)
  131     FORMAT(/, 3X,'Pseudopotentials Grid    : Channels = ',I2,&
                   ', Mesh = ',I5,/,30X,'dx   = ',F16.14)
  132     FORMAT(   3X,'point    radius        vloc         ( vnl - vloc )')
  120     FORMAT(I8,E14.6,5E14.6)

        ELSE

          WRITE( stdout,25) ap%igau
          WRITE( stdout,30)
          WRITE( stdout,104) ap%wrc(1),ap%rc(1),ap%wrc(2),ap%rc(2)
   25     FORMAT(/, 3X,'Gaussians used : ',I2,'. Parameters are : ')
   30     FORMAT(   3X,'C (core), Alfa(core) : ')
  104     FORMAT(4(3X,F8.4))

          WRITE( stdout,40)
          DO il=1,3
            DO ib=1,ap%igau
              WRITE( stdout,103) ap%rcl(ib,il),ap%al(ib,il),ap%bl(ib,il)
            END DO
          END DO
   40     FORMAT(   3X,'Hsc radii and coeff. A and B :')
  103     FORMAT(3X,F8.4,2(3X,F15.7))


        END IF

        IF( ap%nrps > 0 .AND. ap%mesh > 0 ) THEN
          WRITE( stdout,141) ap%nrps, ap%mesh, ap%dx
          in1=1
          in2=ap%mesh/4
          in3=ap%mesh/2
          in4=ap%mesh
          WRITE( stdout,145) (ap%oc(i),i=1,ap%nrps)
          WRITE( stdout,142)
          WRITE( stdout,120) in1,ap%rw(in1),(ap%rps(in1,m),m=1,ap%nrps)
          WRITE( stdout,120) in2,ap%rw(in2),(ap%rps(in2,m),m=1,ap%nrps)
          WRITE( stdout,120) in3,ap%rw(in3),(ap%rps(in3,m),m=1,ap%nrps)
          WRITE( stdout,120) in4,ap%rw(in4),(ap%rps(in4,m),m=1,ap%nrps)
        END IF

  141   FORMAT(/, 3X,'Atomic wavefunction Grid : Channels = ',I2,&
                   ', Mesh = ',I5,/,30X,'dx   = ',F16.14)
  142   FORMAT(   3X,'point      radius        wavefunction')
  145   FORMAT(   3X,'Channels occupation number : ',5F10.4)

        IF( ap%tnlcc ) THEN
          WRITE( stdout,151) ap%mesh, ap%dx
          in1 = 1
          in2 = ap%mesh / 4
          in3 = ap%mesh / 2
          in4 = ap%mesh
          WRITE( stdout,152)
          WRITE( stdout,120) in1,ap%rw(in1),ap%rhoc(in1)
          WRITE( stdout,120) in2,ap%rw(in2),ap%rhoc(in2)
          WRITE( stdout,120) in3,ap%rw(in3),ap%rhoc(in3)
          WRITE( stdout,120) in4,ap%rw(in4),ap%rhoc(in4)
        END IF

  151   FORMAT(/, 3X,'Core correction Grid     : Mesh = ',I5, &
             ', dx   = ',F16.14)
  152   FORMAT(   3X,'point      radius        rho core')

        RETURN
      END SUBROUTINE ap_info

!=----------------------------------------------------------------------------=!

SUBROUTINE read_atomic_wf( iunit, ap, err_msg, ierr)

  USE pseudo_types, ONLY: pseudo_ncpp
  USE parser, ONLY: field_count

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, mesh
  REAL(DP) :: rdum

! ... read atomic wave functions
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic wf '

  ap%rps  = 0.0_DP
  ap%nrps = 0
  ap%oc   = 0.0d0

  ap%lrps = 0

  ! this is for local pseudopotentials
  IF( ap%nbeta == 0 ) RETURN
              
  READ(iunit,'(A80)',end=100) input_line
  CALL field_count(nf, input_line)

  strlen = len_trim(input_line)

  IF( nf == 2 ) THEN
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps
  ELSE
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps, ( ap%oc(j), j=1, MIN(ap%nrps,SIZE(ap%oc)) )
  END IF
  IF( ap%nrps > SIZE(ap%rps,2) ) THEN
    ierr = 2   
    err_msg = ' NCHAN NOT PROGRAMMED '
    GO TO 110
  END IF
  IF( mesh > SIZE(ap%rw) .OR. mesh < 0) THEN
    ierr = 4
    err_msg = ' WAVMESH OUT OF RANGE '
    GO TO 110
  END IF

  DO j = 1, mesh
    READ(iunit,*,IOSTAT=ierr) rdum, (ap%rps(j,m),m=1,ap%nrps)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( ABS(rdum - ap%rw(j))/(rdum+ap%rw(j)) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' radial meshes do not match '
      GO TO 110
    END IF
  END DO

  ! In this format each columns is an atomic functions of 
  ! increasing angular momentum ( starting from L=0 )
  !
  DO m = 1, ap%nrps
    ap%lrps( m ) = m - 1
  END DO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE read_atomic_wf

!=----------------------------------------------------------------------------=!

SUBROUTINE read_numeric_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, l, ll

! ... read numeric atomic pseudopotential
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic numeric pseudo '

  IF(ap%tmix) THEN
    READ(iunit,*) (ap%wgv(l),l=1,ap%nbeta)
  END IF

  READ(iunit,*,IOSTAT=ierr) ap%zv
  READ(iunit,*,IOSTAT=ierr) ap%mesh, ap%nchan

  IF((ap%nchan > SIZE(ap%vnl,2) ) .OR. (ap%nchan < 1)) THEN
    ierr = 1
    err_msg = ' NCHAN NOT PROGRAMMED '
    GO TO 110
  END IF
  IF((ap%mesh > SIZE(ap%rw) ) .OR. (ap%mesh < 0)) THEN
    info = 2
    err_msg = ' NPOTMESH OUT OF RANGE '
    GO TO 110
  END IF

  ap%rw = 0.0d0
  ap%vnl = 0.0d0
  ap%vloc = 0.0d0
  ap%vrps = 0.0d0
  DO j = 1, ap%mesh
    READ(iunit,*,IOSTAT=ierr) ap%rw(j), (ap%vnl(j,l),l=1,ap%nchan)
  END DO

  IF( MINVAL( ap%rw(1:ap%mesh) ) <= 0.0d0 ) THEN
    info = 30
    err_msg = ' ap rw too small '
    GO TO 110
  END IF

! ...  mixed reference potential is in vr(lloc)
  IF(ap%tmix) THEN
    DO j=1,ap%mesh
      ap%vnl( j, ( ap%lloc + 1 ) )= 0.d0
      DO l=1,ap%nchan
        IF( l /= ( ap%lloc + 1 ) ) THEN
          ap%vnl(j, ( ap%lloc + 1 ) )=  ap%vnl(j, ( ap%lloc + 1 ) ) + ap%wgv(l) * ap%vnl(j,l)
        END IF
      END DO
    END DO
  END IF
  ap%vloc(:) = ap%vnl( :, ap%lloc + 1 )
  ap%dx = calculate_dx( ap%rw, ap%mesh )
  ap%rab  = ap%dx * ap%rw

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GO TO 110

  DO l = 1, ap%nbeta
    ll=ap%lll(l) + 1
    ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
  END DO

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap,  err_msg, ierr)
    IF( ierr /= 0 ) GO TO 110
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE read_numeric_pp

!=----------------------------------------------------------------------------=!
!
!
!
!---------------------------------------------------------------------
subroutine ncpp2internal ( ap, is, xc_type, ierr )
  !---------------------------------------------------------------------
  !
  !   convert and copy "is"-th pseudopotential in the Unified Pseudopotential 
  !   Format to internal PWscf variables
  !   return error code in "ierr" (success: ierr=0)
  !
  ! CP90 modules
  !
  use uspp_param, only: zp, qfunc, qfcoef, rinner, qqq, vloc_at, &
                   lll, nbeta, kkbeta,  nqlc, nqf, betar, dion, tvanp
  use atom, only: chi, lchi, nchi, rho_atc, rgrid, nlcc, numeric, oc
  use funct, only: set_dft_from_name, dft_is_hybrid
  !
  use pseudo_types
  !
  implicit none
  !
  integer :: ierr 
  INTEGER, INTENT(IN) :: is
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(len=*) :: xc_type
  !
  !     Local variables
  !
  integer :: il, ir, ic
  real(DP), allocatable :: fint(:)
  !
  zp(is)     = ap%zv
  tvanp(is)  = .FALSE.
  nlcc(is)   = ap%tnlcc
  rgrid(is)%mesh   = ap%mesh
  nchi(is)   = ap%nrps
  nbeta(is)  = ap%nbeta
  kkbeta(is) = rgrid(is)%mesh
  nqlc(is)   = 0 ! upf%nqlc
  nqf (is)   = 0 ! upf%nqf
  if (rgrid(is)%mesh > SIZE(rgrid(is)%r) ) call errore('read_pseudo','increase ndmx',rgrid(is)%mesh)
  !
  call set_dft_from_name( TRIM( xc_type ) )
  IF ( dft_is_hybrid() ) &
     CALL errore( 'read_pseudo', 'HYBRID XC not implemented in CPV', 1 )

  !
  !
  oc  ( 1 : ap%nrps, is ) = ap%oc  ( 1 : ap%nrps )
  lchi( 1 : ap%nrps, is ) = ap%lrps( 1 : ap%nrps )
  chi ( 1 : ap%mesh, 1 : ap%nrps, is ) = ap%rps( 1 : ap%mesh, 1 : ap%nrps )
  !
  betar( 1 : ap%mesh, 1 : ap%nbeta, is ) = 2.0d0 * ap%vrps( 1 : ap%mesh, 1 : ap%nbeta )
  !
  lll  ( 1 : ap%nbeta, is ) = ap%lll( 1 : ap%nbeta )  ! = upf%lll( 1:upf%nbeta )

  rinner(:,is) = 0.0d0
  qqq(:,:,is)  = 0.0d0
  qfunc(:,:,is) = 0.0d0
  qfcoef(:,:,:,:,is) = 0.0d0

  !
  rgrid(is)%r  (1:ap%mesh) = ap%rw( 1:ap%mesh )         ! = upf%r  (1:upf%mesh)
  rgrid(is)%rab(1:ap%mesh) = ap%dx * ap%rw( 1:ap%mesh ) ! = upf%rab(1:upf%mesh)
  !
  rho_atc (:,is) = 0.d0
  if ( ap%tnlcc ) then
     rho_atc (1:ap%mesh, is) = ap%rhoc(1:ap%mesh)        ! = upf%rho_atc(1:upf%mesh)
  end if
  !
  ! rsatom (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  ! lloc(is) = 1
  !
  vloc_at (:, is) = 0.0d0 
  vloc_at (1:ap%mesh, is) = 2.0d0 * ap%vloc( 1:ap%mesh ) ! = upf%vloc(1:upf%mesh)

  dion(:,:,is) = 0.0d0                                   ! upf%dion(1:upf%nbeta, 1:upf%nbeta)
  allocate(fint(rgrid(is)%mesh))
  do il = 1, nbeta(is)
    do ic = 1, nchi(is)
      if( lchi( ic, is ) == lll( il, is ) ) exit
    end do
    do ir = 1, rgrid(is)%mesh
      fint(ir) = chi( ir, ic, is ) * 2.0d0 * ap%vrps( ir, il )
    end do
    call simpson_cp90( rgrid(is)%mesh, fint, rgrid(is)%rab, dion(il,il,is) )
    dion(il,il,is) = 1.0d0/dion(il,il,is)
  end do
  deallocate(fint)
  !
  return
end subroutine ncpp2internal

!=----------------------------------------------------------------------------=!


SUBROUTINE upf2ncpp( upf, ap )

  !
  !   convert and copy upf norm conserving pseudo to internal FPMD variables
  !

  use pseudo_types
  use constants, ONLY: eps14

  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  TYPE (pseudo_upf ), INTENT(INOUT) :: upf

  integer :: l, il, i
  integer :: which_lloc( 0 : upf%nbeta )

  ap%rw    = 0.0d0
  ap%vnl   = 0.0d0
  ap%vrps  = 0.0d0
  ap%rps   = 0.0d0
  ap%oc    = 0.0d0
  ap%tmix  = .FALSE.
  ap%tnlcc = upf%nlcc
  !
  ap%zv   = upf%zp
  ap%nbeta  = upf%nbeta

  ap%lll( 1:upf%nbeta ) = upf%lll( 1:upf%nbeta )

  !  Calculate lloc
  ap%lloc = upf%nbeta
  which_lloc = 0
  DO l = 1, ap%nbeta
    which_lloc( ap%lll( l ) ) = 1
  END DO
  !
  !  the first "l" which is not non-local
  !  is taken as the "l" of the local part of the pseudo
  !
  LLOC: DO l = 0, ap%nbeta
    IF( which_lloc( l ) == 0 ) THEN
      ap%lloc = l
      exit LLOC
    END IF
  END DO LLOC

  ap%nchan = upf%nbeta + 1   ! projectors and local part
  ap%mesh  = upf%mesh
  ap%rw( 1:upf%mesh )     = upf%r( 1:upf%mesh )
  ap%vnl( 1:upf%mesh, 1 ) = upf%vloc( 1:upf%mesh ) / 2.0d0  ! Rydberg to Hartree atomic units
  ap%dx   = calculate_dx( ap%rw, ap%mesh )
  IF( ap%rw( 1 ) < eps14 ) THEN
     ap%rab  = upf%rab 
  ELSE
     ap%rab  = ap%dx * ap%rw
  END IF
  !
  ap%vloc( 1:upf%mesh ) = upf%vloc( 1:upf%mesh ) / 2.0d0
  ap%nrps = upf%nwfc
  ap%lrps( 1:upf%nwfc ) = upf%lchi( 1:upf%nwfc )
  ap%oc  ( 1:upf%nwfc ) = upf%oc( 1:upf%nwfc )
  ap%rps( 1:upf%mesh, 1:upf%nwfc ) = upf%chi( 1:upf%mesh, 1:upf%nwfc )
  
  DO l = 1, ap%nbeta

    !  vrps(i, l) = ( vnl(i, l) - vloc(i) ) * rps(i, l)
     
    ap%vrps( 1:upf%mesh, l ) = upf%beta( 1:upf%mesh, l ) / 2.0d0 

  END DO
 
  IF( ap%tnlcc ) THEN
    ap%rhoc = 0.0d0
    ap%rhoc(1:upf%mesh) = upf%rho_atc(1:upf%mesh)
  END IF

  RETURN
END SUBROUTINE upf2ncpp

!=----------------------------------------------------------------------------=!

SUBROUTINE read_head_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  INTEGER :: i, l

! ... read pseudo header

  ierr = 0
  err_msg = ' error while reading header pseudo '

  ap%lll = 0
  READ(iunit, *) ap%tnlcc, ap%tmix
  READ(iunit, *) ap%pottyp, ap%lloc, ap%nbeta, (ap%lll(l), l = 1, MIN(ap%nbeta, SIZE(ap%lll)) ) 

  ap%lll = ap%lll - 1
  ap%lloc = ap%lloc - 1

  IF( ap%nbeta > SIZE(ap%lll) .OR. ap%nbeta < 0 ) THEN
    ierr = 1
    err_msg = 'nbeta out of range'
    GO TO 110
  END IF
  IF( ( ap%lloc + 1 ) < 1 .OR. ( ap%lloc + 1 ) > SIZE( ap%vnl, 2 ) ) THEN
    ierr = 3
    err_msg = 'LLOC out of range'
    GO TO 110
  END IF
  IF( ap%tmix .AND. ap%pottyp /= 'NUMERIC' ) THEN
    ierr = 4
    err_msg = 'tmix not implemented for pseudo ' // ap%pottyp
    GO TO 110
  END IF
  DO l = 2, ap%nbeta
    IF( ap%lll(l) <= ap%lll(l-1)) THEN
      ierr = 5
      err_msg =' NONLOCAL COMPONENTS MUST BE GIVEN IN ASCENDING ORDER'
      GO TO 110
    END IF
  END DO
  DO l = 1, ap%nbeta
    IF( ap%lll(l) == ap%lloc ) THEN
      ierr = 6
      err_msg = ' LLOC.EQ.L NON LOCAL!!' 
      GO TO 110
    END IF
  END DO

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE read_head_pp

!=----------------------------------------------------------------------------=!

SUBROUTINE read_analytic_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  INTEGER :: i, l

! ... read analytic pseudo gaussians

  ierr = 0
  err_msg = ' error while reading atomic analytic pseudo '

  READ(iunit,*,IOSTAT=ierr) ap%zv, ap%igau

  ap%mesh = 0 
  ap%nchan = 0 
  ap%dx = 0.0d0
  ap%rab  = 0.0d0
  ap%rw   = 0.0d0
  ap%vnl   = 0.0d0
  ap%vloc   = 0.0d0
  ap%vrps   = 0.0d0

  SELECT CASE (ap%igau)
    CASE ( 1 )
      READ(iunit,*,IOSTAT=ierr) ap%rc(1)
      ap%wrc(1) = 1.d0
      ap%wrc(2) = 0.d0
      ap%rc(2)  = 0.d0
    CASE ( 3 )
      READ(iunit,*,IOSTAT=ierr) ap%wrc(1), ap%rc(1), ap%wrc(2), ap%rc(2)
    CASE DEFAULT
      ierr = 1
      err_msg = ' IGAU NOT PROGRAMMED '
      GO TO 110
  END SELECT

  DO l=1,3
    DO i=1,ap%igau
      READ(iunit,*,IOSTAT=ierr) ap%rcl(i,l), ap%al(i,l), ap%bl(i,l)
    END DO
  END DO

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GO TO 110

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap, err_msg, ierr)
    IF( ierr /= 0 ) GO TO 110
  END IF

! ... Analytic pseudo are not supported anymore, conversion
! ... to numeric form is forced
  CALL analytic_to_numeric( ap )

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE read_analytic_pp

!=----------------------------------------------------------------------------=!

SUBROUTINE read_atomic_cc( iunit, ap, err_msg, ierr )

  !  this subroutine reads core correction charge mesh

  USE pseudo_types, ONLY: pseudo_ncpp

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: j, mesh
  REAL(DP) :: rdum

! ... read atomic core

  ierr = 0
  err_msg = ' error while reading atomic core pseudo '

  ap%rhoc = 0.0d0

  READ( iunit, *, IOSTAT = ierr ) mesh
  IF( mesh > SIZE( ap%rw ) .OR. mesh < 0 ) THEN
    ierr = 17
    err_msg = '  CORE CORRECTION MESH OUT OF RANGE '
    GO TO 110
  END IF
  DO j = 1, mesh
    READ( iunit, *, IOSTAT = ierr ) rdum, ap%rhoc(j)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( ABS( rdum - ap%rw(j) ) / ( rdum + ap%rw(j) ) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' core cor. radial mesh does not match '
      GO TO 110
    END IF
  END DO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx   = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE read_atomic_cc

!=----------------------------------------------------------------------------=!
   END MODULE read_pseudo_module_fpmd
!=----------------------------------------------------------------------------=!
!
!
!
!     
!---------------------------------------------------------------------
      subroutine readbhs( is, iunps )
!---------------------------------------------------------------------
!
      use atom, only: rgrid, nlcc, rho_atc, numeric
      use uspp_param, only: zp, betar, dion, vloc_at, lll, nbeta, kkbeta
      use bhs, only: rcl, rc2, bl, al, wrc1, lloc, wrc2, rc1
      use funct, only: set_dft_from_name, dft_is_hybrid
      use io_global, only: stdout

!
      implicit none
!
      integer is, iunps
!
      integer meshp, ir, ib, il, i, j, jj
      real(8), allocatable:: fint(:), vnl(:)
      real(8) rdum, alpha, z, zval, cmesh, cmeshp, exfact
      character(len=20) :: dft_name
!
! nlcc is unfortunately not read from file
!
      numeric(is) = .false.
      nlcc(is)=.false.
      read(iunps,*) z,zp(is),nbeta(is),lloc(is),exfact
      if (zp(is) < 1 .or. zp(is) > 100 ) then
         call errore('readbhs','wrong potential read',15)
      endif

      call dftname_cp (nint(exfact), dft_name)
      call set_dft_from_name( dft_name )
      IF ( dft_is_hybrid() ) &
         CALL errore( 'readbhs', 'HYBRID XC not implemented in CPV', 1 )
!
      if(lloc(is).eq.2)then 
         lll(1,is)=0
         lll(2,is)=1
      else if(lloc(is).ne.2) then
         call errore('readbhs','kb-ization for lloc=2 only',10)
      endif
!     
!     see eqs. (2.21) and (2.22) of bhs, prb 26, 4199 (1982).
!
!     wrc1  =c_core(1)
!     wrc2  =c_core(2)
!     rc1   =alpha_core(1)
!     rc2   =alpha_core(2)
!     al(i) =a(i)          i=1,3
!     bl(i) =a(i+3)        i=1,3
!     rcl(i)=alpha(i)      i=1,3 
!
!     ------------------------------------------------------------------
!     pp parameters are read from file iunps
!     bhs 's coefficients have been turned into lengths
!     ------------------------------------------------------------------
      read(iunps,*) wrc1(is),rc1(is),wrc2(is),rc2(is)  
      rc1(is)=1.0d0/sqrt(rc1(is))
      rc2(is)=1.0d0/sqrt(rc2(is))
      do il=1,3
         do ib=1,3
            read(iunps,*) rcl(ib,is,il),al(ib,is,il),bl(ib,is,il)
            rcl(ib,is,il)=1.0d0/sqrt(rcl(ib,is,il))
         end do
      end do
!
!     ------------------------------------------------------------------
!     wavefunctions are read from file iunps
!     ------------------------------------------------------------------
      do il=1,nbeta(is)
         read(iunps,*) rgrid(is)%mesh,cmesh
!
! kkbeta is for compatibility with Vanderbilt PP
!
         kkbeta(is)=rgrid(is)%mesh
         do j=1,rgrid(is)%mesh
            read(iunps,*) jj,rgrid(is)%r(j),betar(j,il,is)
         end do
      end do
!     
!     ------------------------------------------------------------------
!     core charge is read from unit 15
!     ------------------------------------------------------------------
!
      if(nlcc(is)) then
         read(15,*) meshp,cmeshp
         if ( meshp.ne.rgrid(is)%mesh .or. cmeshp.ne.cmesh ) then
            call errore('readbhs','core charge mesh mismatch',is)
         endif
         do ir=1,rgrid(is)%mesh
            read(15,*) rdum, rho_atc(ir,is)
         end do
      endif
!
!  rab(i) is the derivative of the radial mesh
!
      do ir=1,rgrid(is)%mesh
         rgrid(is)%rab(ir)=rgrid(is)%r(ir) * log(cmesh)
      end do
!
!     ------------------------------------------------------------------
!     local potential 
!     ------------------------------------------------------------------
      lloc(is)=lloc(is)+1
!
! NB: the following is NOT the local potential: the -ze^2/r term is missing
!
      do ir=1,rgrid(is)%mesh
         vloc_at(ir,is)=0.d0
         do i=1,3
            vloc_at(ir,is) = vloc_at(ir,is)                             &
     &            +(al(i,is,lloc(is))+bl(i,is,lloc(is))*rgrid(is)%r(ir)**2)    &
     &            *exp(-(rgrid(is)%r(ir)/rcl(i,is,lloc(is)))**2)
         end do
      end do
!
!     ------------------------------------------------------------------
!     nonlocal potentials: kleinman-bylander form 
!     (1) definition of betar   (2) calculation of dion 
!     ------------------------------------------------------------------
      allocate(fint(rgrid(is)%mesh), vnl(rgrid(is)%mesh))
      do il=1,nbeta(is)
         do ir=1,rgrid(is)%mesh
            vnl(ir)=0.d0
            do i=1,3
               vnl(ir) = vnl(ir) + (al(i,is,il)+bl(i,is,il)*rgrid(is)%r(ir)**2)&
     &                    * exp(-(rgrid(is)%r(ir)/rcl(i,is,il))**2)
            end do
            vnl(ir) = vnl(ir) - vloc_at(ir,is)
            fint(ir)= betar(ir,il,is)**2*vnl(ir)
            betar(ir,il,is)=vnl(ir)*betar(ir,il,is)
         end do
         call simpson_cp90(rgrid(is)%mesh,fint,rgrid(is)%rab,dion(il,il,is))
         dion(il,il,is) = 1.0d0/dion(il,il,is)
      end do
      deallocate(vnl, fint)
!     
!     ------------------------------------------------------------------
!     output: pp info 
!     ------------------------------------------------------------------
      WRITE( stdout,3000) z,zp(is)
3000  format(2x,'bhs pp for z=',f3.0,2x,'zv=',f3.0)

      WRITE( stdout,'(2x,a20)') dft_name
      WRITE( stdout,3002) lloc(is)-1 
3002  format(2x,'   local angular momentum: l=',i3)
      WRITE( stdout,3005) nbeta(is)
3005  format(2x,'number of nl ang. mom. nbeta=',i3)
      do il=1,nbeta(is)
         WRITE( stdout,3010) lll(il,is)
3010     format(2x,'nonlocal angular momentum: l=',i3)
      end do
      WRITE( stdout,3030) 
3030  format(2x,'pseudopotential parameters:')
      WRITE( stdout,3035) wrc1(is),1.0d0/rc1(is)**2
3035  format(2x,'core:',2x,'c1_c=',f7.4,' alpha1_c=',f7.4)
      WRITE( stdout,3036) wrc2(is),1.0d0/rc2(is)**2
3036  format(2x,'     ',2x,'c2_c=',f7.4,' alpha2_c=',f7.4)
      WRITE( stdout,3038)
3038  format(2x,'other table parameters:')
      do il=1,3
         WRITE( stdout,3040) il-1
3040     format(2x,'l=',i3)
         do i =1,3
            alpha=1.0d0/rcl(i,is,il)**2
            WRITE( stdout,3050) i,alpha,i,al(i,is,il),i+3,bl(i,is,il)
         end do
      end do
3050  format(2x,'alpha',i1,'=',f6.2,'  a',i1,'=',f16.7,                 &
     &           '  a',i1,'=',f16.7)
      WRITE( stdout,*)
!     
      return
      end subroutine readbhs
!
!  Description of the Native FPMD pseudopotential format
!
!  The format of the file must be as follows
!  (lowercase text and }'s are comments):
!
!  When POTTYP = 'ANALYTIC' the layout is:
!
!    TCC      TMIX                additional stuff on each line is ignored
!    POTTYP   LLOC LNL ( INDL(i), i = 1, LNL )
!    ( WGV(i), i = 1, LNL )       this line only if tmix(is) is true
!    ZV       IGAU                igau must be 1 or 3     }
!    WRC(1) RC(1) WRC(2) RC(2)    this line if igau = 3   }
!    RC(1)                        this one if igau = 1    }
!    RCL(1,1)    AL(1,1)    BL(1,1)         }             }  this
!     ...         ...        ...            }  l = 0      }  section
!    RCL(IGAU,1) AL(IGAU,1) BL(IGAU,1)      }             }  only if
!    RCL(1,2)    AL(1,2)    BL(1,2)      }                }  pottyp is
!     ...         ...        ...         }     l = 1      }  'ANALYTIC'
!    RCL(IGAU,2) AL(IGAU,2) BL(IGAU,2)   }                }
!    RCL(1,3)    AL(1,3)    BL(1,3)         }             }
!     ...         ...        ...            }  l = 2      }
!    RCL(IGAU,3) AL(IGAU,3) BL(IGAU,3)      }             }
!    NMESH NCHAN                                       }
!    RW( 1 )     ( RPS( 1, j ), j = 1, NCHAN )         }  pseudowave
!     ...         ...              ...                 }
!    RW( NMESH ) ( RPS( NMESH, j ), j = 1, NCHAN )     }
!
!  
!  When POTTYP = 'NUMERIC' the layout is:
!
!    TCC      TMIX             additional stuff on each line is ignored
!    POTTYP   LLOC LNL  ( INDL(i), i = 1, LNL )
!    ( WGV(i), i = 1, LNL )       this line only if tmix(is) is true
!    ZV                                             }
!    NMESH NCHAN                                    }    this if
!    RW( 1 )     ( VR( 1, j ), j = 1, NCHAN )       }    pottyp is
!     ...       ...             ...                 }    'NUMERIC'
!    RW( NMESH ) ( VR( NMESH, j ), j = 1, NCHAN )   }
!    NMESH NCHAN                                       }
!    RW( 1 )     ( RPS( 1, j ), j = 1, NCHAN )         }  pseudowave
!     ...         ...              ...                 }
!    RW( NMESH ) ( RPS( NMESH, j ), j = 1, NCHAN )     }
!
!  DETAILED DESCRIPTION OF INPUT PARAMETERS:
!
!    TCC      (logical)   True if Core Correction are required for this 
!                         pseudo
!
!    TMIX     (logical)   True if we want to mix nonlocal pseudopotential 
!                         components 
!
!    WGV(i)   (real)      wheight of the nonlocal components in the 
!                         pseudopotential mixing scheme 
!                         These parameters are present only if TMIX = .TRUE.
!                         1 <= i <= LNL
!                           
!    POTTYP   (character) pseudopotential type
!                         pottyp = 'ANALYTIC' : use an analytic expression
!                         pottyp = 'NUMERIC'  : read values from a table
!
!    ZV       (integer)   valence for each species
!
!    IGAU     (integer)   number of Gaussians in the pseudopotentials
!                         expression used only if pottyp='ANALYTIC'
!
!  parameters from Bachelet-Hamann-Schluter's table:
!
!    WRC(2)   (real)      c1, c2 (core)  parameters
!    RC(2)    (real)      alpha1, alpha2 parameters
!
!    RCL(i,3) (real)      alpha1, alpha2, alpha3 for each angular momentum
!                         1 <= i <= IGAU
!    AL(i,3)  (real)      parameters for each angular momentum
!                         1 <= i <= IGAU
!    BL(i,3)  (real)      parameters for each angular momentum
!                         1 <= i <= IGAU
!
!  nonlocality
!    IGAU     (integer)   number of Gaussians for analytic pseudopotentials
!    LLOC     (integer)   index of the angular momentum component added to 
!                          the local part  ( s = 1, p = 2, d = 3 )
!    LNL      (integer)   number of non local component
!    INDL(i)  (integer)   indices of non local components
!                         1 <= i <= LNL
!                         ( 1 3 means s and d taken as non local )
!
!  pseudo grids
!    NMESH    (integer)   number of points in the mesh mesh
!    NCHAN    (integer)   numbero of colums, radial components
!    RW(i)    (real)      distance from the core in A.U. (radial mesh)
!                         1 <= i <= NMESH
!    RPS(i,j) (real)      Atomic pseudo - wavefunctions
!                         1 <= i <= NMESH ; 1 <= j <= NCHAN
!    VP(i,j)  (real)      Atomic pseudo - potential
!                         1 <= i <= NMESH ; 1 <= j <= NCHAN
!
!  ----------------------------------------------
!  END manual

