!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
program band_interpolation 
USE io_global,         ONLY : stdout
USE globalmod,         ONLY : print_bands, read_xml_input, method, deallocate_global 
USE idwmod,            ONLY : idw
USE fouriermod,        ONLY : fourier, fourierdiff
USE mp_global,         ONLY : mp_startup
implicit none
  !
  write(stdout,'(A)') 'PROGRAM: band_interpolation '
  !
  CALL mp_startup ( )
  !
  Call set_defaults () 
  !
  Call read_input_file ()
  !
  Call read_xml_input ()
  !
  if(TRIM(method).eq.'idw') then 
    !
    Call idw (1)
    !
  elseif(TRIM(method).eq.'idw-sphere') then 
    !
    Call idw (2)
    !
  elseif(TRIM(method).eq.'fourier') then 
    !
    Call fourier ()
    !
  elseif(TRIM(method).eq.'fourier-diff') then 
    !
    Call fourierdiff ()
    !
  else
    !
    write(stdout, '(A,A)') 'method: ', TRIM(method)
    Call errore('band_interpolation ' , ' wrong method ', 1)
    !
  end if
  !
  Call print_bands (TRIM(method))
  !
  Call deallocate_global ()
  !
  Call stop_pp ()
  !
end program
!----------------------------------------------------------------------------
subroutine set_defaults () 
!
! Set defaults 
!  
USE globalmod,         ONLY : method
USE fouriermod,        ONLY : check_periodicity, miller_max, RoughN, RoughC, NUser 
USE idwmod,            ONLY : p_metric, scale_sphere 
USE input_parameters,  ONLY : k_points 
implicit none
  !
  ! global defauls
  method = 'fourier-diff'
  k_points = 'none'
  !
  ! defaults for IDW methods
  p_metric = 2
  scale_sphere = 4.0d0
  !
  ! defaults for fourier methods
  check_periodicity = .false.
  miller_max = 6  
  RoughN = 1  
  allocate( RoughC(RoughN) )
  RoughC(1) = 1.0d0
  NUser = 0 
  !
  RETURN
  !
end subroutine
!----------------------------------------------------------------------------
subroutine read_input_file ()
!
! Read input file
!
USE io_global,         ONLY : stdout
USE parser,            ONLY : read_line
USE input_parameters,  ONLY : k_points, nkstot 
USE read_cards_module, ONLY : card_kpoints 
USE globalmod,         ONLY : method
USE fouriermod,        ONLY : miller_max, check_periodicity, card_user_stars, card_roughness
USE idwmod,            ONLY : p_metric, scale_sphere
implicit none
  integer, parameter :: iunit = 5
  integer :: ios, i
  CHARACTER(len=256)         :: input_line
  LOGICAL                    :: tend
  CHARACTER(len=80)          :: card
  CHARACTER(len=1), EXTERNAL :: capital
  !
  NAMELIST / interpolation / method, miller_max, check_periodicity, p_metric, scale_sphere
  !
  ios = 0 
  READ( iunit, interpolation, iostat = ios ) 
  !
100   CALL read_line( input_line, end_of_file=tend )
  !
  IF( tend ) GOTO 120
  IF( input_line == ' ' .OR. input_line(1:1) == '#' .OR. &
                             input_line(1:1) == '!' ) GOTO 100
  !
  READ (input_line, *) card
  !
  DO i = 1, len_trim( input_line )
     input_line( i : i ) = capital( input_line( i : i ) )
  ENDDO
  !
  IF ( trim(card) == 'ROUGHNESS' ) THEN
     !
     CALL card_roughness( input_line )
     !
  ELSEIF ( trim(card) == 'USER_STARS' ) THEN
     !
     CALL card_user_stars( input_line )
     !
  ELSEIF ( trim(card) == 'K_POINTS' ) THEN
     !
     CALL card_kpoints( input_line )
     !
  ELSE
     !
     WRITE( *,'(A)') 'Warning: card '//trim(input_line)//' ignored'
     !
  ENDIF
  !
  GOTO 100
  !
120  CONTINUE
  !
  if(k_points.ne.'tpiba') then 
    write(stdout,'(A,A)') 'k_points = ', k_points 
    Call errore('band_interpolation ' , ' K_POINTS card must be specified with tpiba_b ', 1)
  end if
  !
  if(nkstot.le.0) then 
    write(stdout,'(A,I5)') 'nkstot = ', nkstot
    Call errore('band_interpolation ' , ' wrong number of k-points ', 1)
  end if 
  !
  write(stdout,'(A,A)') 'Interpolation method: ', method
  if( TRIM(method).ne.'idw'.and.TRIM(method).ne.'idw-sphere'&
        .and.TRIM(method).ne.'fourier'.and.TRIM(method).ne.'fourier-diff' ) &  
        Call errore('band_interpolation', 'Wrong interpolation method ', 1) 
  !
  RETURN
  !
end subroutine read_input_file
!----------------------------------------------------------------------------
