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
!      3) Convolute a given spectrum known on a regular energy grid with
!         a constant lorentz brodening gamma_hole, or with an atan 
!         broadening. In more details (s_in, s_out are cross sections):
!
!        s_out(W)=integral s_in(om)*Gamma(W)/(Gamma(W)**2+(om-W)**2)/pi
!
!        and the integral is over all the space (-infty, infty). 
!
!        If conv_type='lorentz' then Gamma(W)=gamma_hole constant
!
!        If conv_type='lorentz_atan' then 
!              Gamma(W)=gamma_hole   for W < xe0
!                      =Gamma_MAX*(0.5+atan(e-1/(e*e))) otherwise
!          where e=(W-xe0)/(Ectr-xe0)
!          and Ectr is the inflection point of the atan and xe0 is
!          the pre-edge onset.
!          For the definition of e see Eq. 7 and below in 
!           O. Bunau and M. Calandra, PRB 87, 205105 (2013) 
!       
!          Finally emin_conv, emax_conv and nenergy_conv define
!          the energy grid of the convoluted spectrum (energies labeled W above).
!          If omitted, the energy grid is the same as in the input cross section.
!
!          Important, please note that the convolution adds on TOP of the
!          width of your original spectrum. So you should run xspectra with a 
!          tiny xgamma.
!------------------------------------------------------------------------------------
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
!              case (iii) - lorentz
!
!                       &input_manip 
!                         cross_section_file='xanes.dat.L23', 
!                         option='convolution',
!                         shift_spectrum=.false.,
!                         element='Cu',
!                         conv_type='lorentz',
!                         gamma_hole=1.0, 
!                         emin_conv=-50.0,
!                         emax_conv=100.0,
!                         nenergy_conv=1000,
!                       &end
!                        
!                       emin_conv, emax_conv and nenergy_conv can be omitted
!
!
!              case (iii) - lorentz_atan
!
!                    &input_manip 
!                       cross_section_file='xanes.dat.L23', 
!                       option='convolution',
!                       shift_spectrum=.false.,
!                       element='Cu',
!                       conv_type='lorentz_atan',
!                       gamma_hole=0.05, 
!                       gamma_max=6.0, 
!                       ectr=-8.0,
!                       xe0=-19.0,
!                    &end
!
!
!                 
!----------------------------------------------------------------------------
Program manip_spectra
  USE kinds, ONLY     : DP
  USE constants,       ONLY : pi
  USE edge_energy, ONLY: getE
! Input
  LOGICAL             :: shift_spectrum
  REAL(kind=dp)       :: xe0
  CHARACTER (LEN=255) :: cross_section_file, option, conv_type
  CHARACTER (LEN=2)   :: element
  
! Other variables

  LOGICAL                    :: found
  INTEGER                    :: i, j
  INTEGER                    :: nargs, iiarg, ierr, ios
  INTEGER                    :: nenergy, istart, i0_l2, nenergy_conv
  REAL(kind=dp)              :: el2, el3, so_splitting, emin_conv, emax_conv, de
  REAL(kind=dp)              :: Ectr, gamma_hole, gamma_max, ee
  REAL(kind=DP), ALLOCATABLE :: cross_section_in(:), energy_in(:)
  REAL(kind=DP), ALLOCATABLE :: cross_section_out(:), energy_out(:)
  REAL(kind=DP), ALLOCATABLE :: cross_section_L3(:), energy_L3(:)
  REAL(kind=DP), ALLOCATABLE :: gamma_conv(:)

  CHARACTER (LEN=256) :: input_file

  namelist / input_manip / &
       cross_section_file, &
       option, &
       xe0, &
       shift_spectrum, &
       element,        &
       conv_type,      &
       emin_conv, emax_conv,  &
       nenergy_conv, &
       gamma_hole, gamma_max ,&
       Ectr
       
  
  !
  ! Defaults
  !
  
  cross_section_file=' '
  option=' '
  xe0=0.d0
  shift_spectrum=.false.
  conv_type='lorentz'
  emin_conv=0.d0
  emax_conv=0.d0
  gamma_hole=-1.d0
  gamma_max=-1.d0
  Ectr=xe0
  !
  ! Read namelist
  !
  
  nargs = command_argument_count()
  found = .FALSE.
  input_file = ' '
  
  DO iiarg = 1, (nargs-1)
     !
     CALL get_command_argument( iiarg, input_file )
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in'    .OR. &
          TRIM( input_file ) == '-i' ) THEN
        !
        CALL get_command_argument( ( iiarg + 1 ) , input_file )
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

  IF( TRIM(ADJUSTL(option)).ne.'cut_occ_states'.AND.TRIM(ADJUSTL(option)).ne.'add_L2_L3'.AND. &
      TRIM(ADJUSTL(option)).ne.'convolution' ) then
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

  write(6,*) 'CALCULATION TYPE:', trim(adjustl(option))


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
     
  ELSEIF(option.eq.'convolution') then
     
     if(dabs(emin_conv).lt.1.d-8.and.dabs(emax_conv).lt.1.d-8) then
        emin_conv=minval(energy_in)-1.0
        emax_conv=maxval(energy_in)+1.0
        nenergy_conv=nenergy
     endif
     
     de=(emax_conv-emin_conv)/real(nenergy_conv)

     allocate(cross_section_out(nenergy_conv))
     allocate(energy_out(nenergy_conv))
     allocate(gamma_conv(nenergy_conv))
     
     do i=1, nenergy_conv
        energy_out(i)=emin_conv+de*(i-1)
     enddo

     
     if(trim(adjustl(conv_type)).eq.'lorentz') then

        if(gamma_hole.lt.0.d0) then
           write(6,*) 'Please specify core-hole width gamma_hole (eV)'
           stop
        endif

        write(6,*) 'Convolution with constant lorentz broadening gamma_hole=',gamma_hole
        gamma_conv=gamma_hole
     elseif(trim(adjustl(conv_type)).eq.'lorentz_atan') then
        if(gamma_hole.lt.0.d0) then
           write(6,*) 'Please specify core-hole width gamma_hole (eV)'
           stop
        endif
        
        if(gamma_max.lt.0.d0) then
           write(6,*) 'Please specify width gamma_max (eV)'
           stop
        endif
        
        if(dabs(Ectr-xe0).lt.1.d-6) then
           write(6,*) 'Please specify Ectr different from xe0 or use standard lorentz broadening'
           stop
        endif
        
        write(6,*) 'Convolution with variable (atan) lorentz broadening'
        write(6,*) 'gamma_hole=', gamma_hole
        write(6,*) 'gamma_max=', gamma_max
        write(6,*) 'Ectr=', ectr, ' xe0=', xe0
        write(6,*) 'Ectr-xe0=', Ectr-xe0
        
        
        do i=1, nenergy_conv
           if(energy_out(i).gt.xe0) then
              ee=(energy_out(i)-xe0)/(Ectr-xe0)
              gamma_conv(i)=gamma_hole+gamma_max*(0.5+atan(ee-1.d0/(ee*ee))/pi)
           else
              gamma_conv(i)=0.d0
           endif
        enddo
     endif

     cross_section_out=0.d0
     do i=1,nenergy_conv
        if(gamma_conv(i).lt.1.d-6) then 
           cross_section_out(i)=cross_section_in(i)
        else
           de=energy_in(2)-energy_in(1)
           do j=1,nenergy
              cross_section_out(i)=cross_section_out(i)+&
                   cross_section_in(j)*gamma_conv(i)/&
                   (gamma_conv(i)*gamma_conv(i)+&
                   (energy_out(i)-energy_in(j))* (energy_out(i)-energy_in(j)) )&
                   *de
           enddo
        endif
     enddo

     cross_section_out=cross_section_out/pi
     
     if(trim(adjustl(conv_type)).eq.'lorentz') then
        OPEN (unit=20,file=trim(adjustl(cross_section_file))//'.convoluted',&
             form='formatted',status='unknown')
     elseif(trim(adjustl(conv_type)).eq.'lorentz_atan') then
        OPEN (unit=20,file=trim(adjustl(cross_section_file))//'.atan_convoluted',&
             form='formatted',status='unknown')
     endif

     rewind(20)
     
     write(20,*)
     write(20,*)
     write(20,*)
     write(20,*)
     
     do i=1,nenergy_conv
        write(20,*) energy_out(i), cross_section_out(i), gamma_conv(i)
     enddo
     
     close(20)
     
     deallocate(cross_section_out)
     deallocate(energy_out)
     deallocate(gamma_conv)
     
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
  
