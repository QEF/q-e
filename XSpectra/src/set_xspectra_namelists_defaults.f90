subroutine set_xspectra_namelists_defaults()
  use xspectra
  USE io_files,      ONLY : prefix
  USE cut_valence_green, ONLY :&
       cut_ierror, &    ! convergence tolerance for one step in the integral
       cut_stepu , &    ! integration initial step, upper side
       cut_stepl , &    ! integration initial step, lower side
       cut_startt, &    ! integration start value of the t variable
       cut_tinf  , &    ! maximum value of the lower integration boundary
       cut_tsup  , &    ! minimum value of the upper integration boudary
       cut_desmooth,&   ! size of the interval near the fermi energy
                        ! in which cross section is smoothed
       cut_nmemu,&      ! size of the memory of the values of the green function, upper side
       cut_nmeml,&      ! size of the memory of the values of the green function, lower side
       cut_occ_states  ! true if you want tou remove occupied states from the spectrum
  USE gamma_variable_mod, ONLY : gamma_value, gamma_energy, &
                                 gamma_lines, gamma_tab, gamma_points, &
                                 gamma_mode, gamma_file

  IMPLICIT NONE
  
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Set default values for some namelist variables
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  ! ... Namelist input_xspectra
  !calculation='xanes'
  calculation='xanes_dipole'
  prefix=' '
  verbosity='low'
  x_save_file='xanes.sav'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  xniter=50
  xcheck_conv=5
  xonly_plot=.FALSE.
  xread_wf=.FALSE.
  xerror=1.d-2
  xiabs=1               !identify the adsorbing atom
  !<DC>
  !DO i=2,3
  !   xkvec(i)=0.d0
  !   xepsilon(i)=0.d0
  !ENDDO
  xkvec(1)=0.d0       
  xkvec(2)=1.d0
  xkvec(3)=0.d0
  xepsilon(1)=1.d0
  xepsilon(2:3)=0.d0
  !ef_r=0.d0
  xe0=xe0_default 
  !</DC>
  xcoordcrys=.true.
  show_status=.false.
  wf_collect=.false.
  U_projection_type='atomic'
  restart_mode='from_scratch'
  time_limit=1.d8
  edge='K'
  lplus=.false.
  lminus=.false.
  
  
  ! ... Namelist plot
  xnepoint=100
  xemin=0.d0
  xemax=10.d0
  xgamma=0.1d0
  cut_occ_states=.FALSE.
  terminator=.false.
  gamma_mode='constant'
  gamma_file='gamma.dat'
  
  
  ! ... Namelist pseudos 
  filecore='Core.wfc'
  
  ! ... Namelist cut_occ (for cutting the occupied states, paste_fermi function)
  cut_ierror=1.d-7
  cut_stepu=1.d-2
  cut_stepl=1.d-3
  cut_startt=1.d0
  cut_tinf=1.d-6
  cut_tsup=100.d0
  cut_desmooth=1.d-2
  cut_nmemu=100000
  cut_nmeml=100000

end subroutine set_xspectra_namelists_defaults
