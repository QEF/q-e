  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/print_clock_ph - Quantum-ESPRESSO group          
  !-----------------------------------------------------------------------
  subroutine print_clock_epw
  !-----------------------------------------------------------------------
  !  
  !  12/2009 The clock data could be better organized
  !  
  use io_global,     ONLY : stdout
  use uspp,          ONLY : nlcc_any
  implicit none
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    Unfolding on the coarse grid'
  CALL print_clock ('elphon_wrap')
  WRITE( stdout, * )
  CALL print_clock ('ELPHWAN')
  WRITE( stdout,  * ) '    INITIALIZATION: '
  CALL print_clock ('epq_init')
  WRITE( stdout, * )
  CALL print_clock ('epq_init')
  IF (nlcc_any) call print_clock ('set_drhoc')
  CALL print_clock ('init_vloc')
  CALL print_clock ('init_us_1')
  CALL print_clock ('newd')
  CALL print_clock ('dvanqq')
  CALL print_clock ('drho')
  WRITE( stdout, * )
  !
  ! Electron-Phonon interpolation 
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    Electron-Phonon interpolation'   
  CALL print_clock ('ephwann')
  CALL print_clock ('ep-interp')
  CALL print_clock('PH SELF-ENERGY')
  CALL print_clock('ABS SPECTRA')
  CALL print_clock ('crys_cart')
  WRITE( stdout, * )
  CALL print_clock ('load data')
  CALL print_clock ('Ham: step 1')
  CALL print_clock ('Ham: step 2')
  CALL print_clock ('Ham: step 3')
  CALL print_clock ('Ham: step 4')
  CALL print_clock ('ep: step 1')
  CALL print_clock ('ep: step 2')
  CALL print_clock ('ep: step 3')
  CALL print_clock ('ep: step 4')
  CALL print_clock ('DynW2B')
  CALL print_clock ('HamW2B')
  CALL print_clock ('ephW2Bp')
  CALL print_clock ('ephW2Bp1')
  CALL print_clock ('ephW2Bp2') 
  !
  ! Eliashberg equations
  WRITE( stdout, * )
  CALL print_clock( 'ELIASHBERG' )
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    Total program execution'
  CALL print_clock ('EPW') 
  !
  end subroutine print_clock_epw
