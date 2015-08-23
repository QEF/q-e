! To be added:
!   atan broadening
!   lorentzian broadening
!   gaussian broadening
!
! Spectra manipulation tool by Oana Bunau and Matteo Calandra, Agust 2015
!
!   This small code allows to
!
!      1) Remove states below a certain energy from the spectrum
!           In large gap insulator it should work as the usual
!           cut_occupied_states option but it is much much faster
!           In metals/semimetals the normal and more time consuming
!           option cut_occupied_states=.true. used in XSpectra 
!           should work better. The procedure use here can sometimes
!           work for metals but care is needed.
!
!      2) For L23 edges, this code generates the full L23 spectrum from
!           the L2 edge only. the L3 edge is obtained multiplying by two the
!           L2 and shifting it by the spin-orbit splitting between 2p1/2
!           and 2p3/2. Thus the spectrum in input must be the L2 edge.
!
!  How to use:
!      
!      1) Use in single processor, no parallelization
!
!      2) Input file description:
!
!              case (i): removal of occupied states from the spectrum
!
!                     &input_manip 
!                       cross_section_file='xanes.dat', 
!                       option='cut_occ_states',
!                       xe0=13.0, 
!                       shift_spectrum=.false.,
!                     &end
!
!                     if shift_spectrum=.true. the spectrum is shifted in xe0
!
!              case (ii): removal of occupied states from the spectrum
!
!                     &input_manip
!                       cross_section_file='xanes.dat',
!                       option='add_L2_L3',
!                       element='Cu',
!                     &end
!
!                     element is used to read the spin orbit splitting 2p1/2, 2p3/2
!
!              The cross_section_file is the xanes.dat file with the cross section.
!                 
!----------------------------------------------------------------------------
Program manip_spectra
  USE kinds, ONLY     : DP
  USE edge_energy, ONLY: getE
  IMPLICIT NONE
! Input
  LOGICAL             :: shift_spectrum
  REAL(kind=dp)       :: xe0
  CHARACTER (LEN=255) :: cross_section_file, option
  CHARACTER (LEN=2)   :: element
  
! Other variables

  LOGICAL                    :: found
  INTEGER                    :: i
  INTEGER                    :: nargs, iiarg, iargc, ierr, ios
  INTEGER                    :: nenergy, istart, i0_l2
  REAL(kind=dp)              :: el2, el3, so_splitting
  REAL(kind=DP), ALLOCATABLE :: cross_section_in(:), energy_in(:)
  REAL(kind=DP), ALLOCATABLE :: cross_section_out(:), energy_out(:)
  REAL(kind=DP), ALLOCATABLE :: cross_section_L3(:), energy_L3(:)

  CHARACTER (LEN=256) :: input_file

  namelist / input_manip / &
       cross_section_file, &
       option, &
       xe0, &
       shift_spectrum, &
       element
  
  !
  ! Defaults
  !
  
  cross_section_file=' '
  option=' '
  xe0=0
  shift_spectrum=.false.
  
  !
  ! Read namelist
  !
  
  nargs = iargc()
  found = .FALSE.
  input_file = ' '
  
  DO iiarg = 1, (nargs-1)
     !
     CALL getarg( iiarg, input_file )
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in'    .OR. &
          TRIM( input_file ) == '-i' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        found = .TRUE.
        EXIT
     ENDIF
     !
  ENDDO
  
  IF (found) THEN
     OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
          STATUS = 'OLD', IOSTAT = ierr )
     IF ( ierr > 0 ) THEN
        CALL errore('iosys', '    input file ' // TRIM( input_file ) // &
             & ' not found' , ierr )
        !        WRITE (stdout, '(/,5x,"*** STOP: input file ",A," not found ***",/)' ) &
        !          TRIM( input_file )
     ENDIF
  ELSE
     ierr = -1
  END IF
  
  
  READ(5, input_manip, err = 200, iostat = ios)
200 CALL errore ('input_xspectra', 'reading input_xspectra namelist', abs (ios) )

  !
  !   Check input
  !

  IF( TRIM(ADJUSTL(option)).ne.'cut_occ_states'.AND.TRIM(ADJUSTL(option)).ne.'add_L2_L3') then
     write(6,*) 'Option not recognized'
     write(6,*) 'Program is stopped'
     stop
  ENDIF

  
  !
  ! Read the cross section
  !

  cross_section_file=TRIM(ADJUSTL(cross_section_file))
  OPEN (unit=277,file=cross_section_file,form='formatted',status='unknown')
  rewind(277)

  call how_many_lines_in_cs_file(277, nenergy)
 
  allocate(cross_section_in(nenergy))
  allocate(energy_in(nenergy))


  read(277,*)
  read(277,*)
  read(277,*)
  read(277,*)

  do i=1,nenergy
     read(277,*) energy_in(i), cross_section_in(i)
  enddo

  close(277)

  

  IF(option.eq.'cut_occ_states') then
     allocate(cross_section_out(nenergy))
     allocate(energy_out(nenergy))
     call removal_of_occupied_states(nenergy, energy_in, cross_section_in, xe0, shift_spectrum, &
          energy_out, cross_section_out)
     OPEN (unit=20,file=trim(adjustl(cross_section_file))//'.cut',form='formatted',status='unknown')
     rewind(20)
     write(20,*)
     write(20,*)
     write(20,*)
     write(20,*)

     do i=1,nenergy
        write(20,*) energy_out(i), cross_section_out(i)
     enddo
     close(20)
     deallocate(cross_section_out)
     deallocate(energy_out)
  ELSEIF(option.eq.'add_L2_L3') then
     !L3=2*L2
     !L3 is at lower energies than L2

     allocate(cross_section_L3(nenergy))
     allocate(energy_L3(nenergy))

     cross_section_L3=2.d0*cross_section_in
     
     el2=getE(element,'L2')
     el3=getE(element,'L3')
     so_splitting=el2-el3


     write(6,*) 'Energy of L2 edge (eV) = ',el2
     write(6,*) 'Energy of L3 edge (eV) = ',el3
     write(6,*) 'SO splitting           = ',so_splitting
     
     energy_L3=energy_in-so_splitting

     if(energy_L3(nenergy).lt.energy_in(1)) then

       OPEN (unit=20,file=trim(adjustl(cross_section_file))//'.L23',&
            form='formatted',status='unknown')

       rewind(20)

       write(20,*)
       write(20,*)
       write(20,*)
       write(20,*)

       do i=1,nenergy
          write(20,*) energy_L3(i), cross_section_L3(i)
       enddo
       write(20,*)
       do i=1,nenergy
          write(20,*) energy_in(i), cross_section_in(i)
       enddo
     
       close(20)

     else
       write(6,*)  'case B'
       allocate(cross_section_out(nenergy))

       cross_section_out=0.d0

        !
        ! determine the first overlap point
        !
  
       do i=1,nenergy-1
          if(energy_L3(i).lt.energy_in(1).and.energy_L3(i+1).gt.energy_in(1)) istart=i
       enddo

       write(6,*) 'istart=',istart

       do i=1,istart
          cross_section_out(i)=cross_section_L3(i)
       enddo

       do i=istart+1, nenergy
          cross_section_out(i)=cross_section_L3(i)+cross_section_in(i-istart)
       enddo 
        
        OPEN (unit=20,file=trim(adjustl(cross_section_file))//'.L23',&
           form='formatted',status='unknown')

        rewind(20)

        write(20,*)
        write(20,*)
        write(20,*)
        write(20,*)

        do i=1,nenergy
           write(20,*) energy_L3(i), cross_section_out(i)
        enddo

        close(20)
 
        deallocate(cross_section_out)
     endif 
     deallocate(cross_section_L3)
     deallocate(energy_L3)
     
  ENDIF


  deallocate(cross_section_in)
  deallocate(energy_in)


end Program manip_spectra

subroutine how_many_lines_in_cs_file( iun, nlines )
  IMPLICIT NONE
  INTEGER :: iun, nlines
  !
  character(256) ::  line
  integer        ::  ios

  rewind(iun)
  !
  !   Read the small header
  !
  read(iun,*)
  read(iun,*)
  read(iun,*)
  read(iun,*)
  
  nlines=0

  !
  !  Count the lines.
  !
  do
     
     read ( iun, '(a)', iostat = ios ) line
     
     if ( ios /= 0 ) then
        exit
     end if
     

     nlines = nlines + 1

  end do

  rewind(iun)

end subroutine how_many_lines_in_cs_file

subroutine removal_of_occupied_states(nenergy, energy_in, cross_section_in, xe0, shift_spectrum, &
     energy_out, cross_section_out)
  USE kinds, ONLY     : DP
  IMPLICIT NONE
  LOGICAL :: shift_spectrum
  INTEGER :: nenergy
  REAL(kind=dp) :: xe0
  REAL(kind=dp) :: energy_in(nenergy), cross_section_in(nenergy)
  REAL(kind=dp) :: energy_out(nenergy), cross_section_out(nenergy)
!
  INTEGER :: i

  energy_out=energy_in
  cross_section_out=0.d0

  do i=1,nenergy
     if(energy_in(i).ge.xe0) cross_section_out(i)=cross_section_in(i) 
  enddo

  if(shift_spectrum) energy_out=energy_out-xe0

  RETURN
end subroutine removal_of_occupied_states
  
