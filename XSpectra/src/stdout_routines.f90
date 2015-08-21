
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Writes the relevant parameters of the XSpectra calculation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


subroutine write_sym_param_to_stdout()
  USE kinds, ONLY : DP
  USE io_files,      ONLY : prefix
  USE xspectra
  USE gamma_variable_mod, ONLY : gamma_value, gamma_energy, &
                                 gamma_lines, gamma_tab, gamma_points, &
                                 gamma_mode, gamma_file
  USE io_global,       ONLY : stdout
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

  IMPLICIT NONE

  INTEGER :: i
  
  if(xang_mom.eq.1) then
     WRITE(stdout, '(5x,a,a,/)') 'calculation: ', 'xanes_dipole'
  elseif(xang_mom.eq.2) then
     WRITE(stdout, '(5x,a,a,/)') 'calculation: ', 'xanes_qyadrupole'
  endif

  !<NM> Write xepsilon and (if needed) xkvec
  IF(.NOT.xonly_plot) THEN
     IF ( xcoordcrys ) THEN
        WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
             'xepsilon  [crystallographic coordinates]: ', (xepsilon(i),i=1,3)
     ELSE
        WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
             'xepsilon  [cartesian coordinates]: ', (xepsilon(i),i=1,3)
     ENDIF
     
     IF ( TRIM(ADJUSTL(calculation)).EQ.'xanes_quadrupole' )  THEN
        IF ( xcoordcrys ) THEN
           WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
                'xkvec  [crystallographic coordinates]: ', (xkvec(i),i=1,3)
        ELSE
           WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
                'xkvec [cartesian coordinates]: ', (xkvec(i),i=1,3)
        ENDIF
     ENDIF
     !<NM>
  ENDIF
  
  
  ! ... Writes xonly_plot, its meaning and plot parameters
  
  IF (xonly_plot.EQV..FALSE.) then
     WRITE(stdout,'(5x,a)') 'xonly_plot: FALSE'
     WRITE(stdout,'(8x,a,/)') &
          '=> complete calculation: Lanczos + spectrum plot'
     WRITE(stdout,'(5x,a,a20)') 'filecore (core-wavefunction file): ', &
          filecore
  ELSE
     WRITE(stdout,'(5x,a)') 'xonly_plot: TRUE'
     WRITE(stdout,'(8x,a)') &
          '=> only the spectrum plot'
  ENDIF
  WRITE(stdout,*)
  
  WRITE(stdout,'(5x,a)') 'main plot parameters:'
  IF (cut_occ_states) THEN
     WRITE(stdout,'(8x,a)') 'cut_occ_states: TRUE'
  ELSE
     WRITE(stdout,'(8x,a)') 'cut_occ_states: FALSE'
  ENDIF
  WRITE(stdout,'(8x,a,a8)')   'gamma_mode:  ', gamma_mode
  IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
     WRITE(stdout,'(8x,a,f5.2)') '-> using xgamma [eV]: ', xgamma
  ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
     WRITE(stdout,'(8x,a,a50)')  '-> using gamma_file: ', gamma_file
  ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
     WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)') &
          '-> first, constant up to point (', &
          gamma_energy(1),',',gamma_value(1),') [eV]'
     WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)') &
          '-> then, linear up to point (',gamma_energy(2),',',gamma_value(2),') [eV]'
     WRITE(stdout,'(8x,a)') '-> finally, constant up to xemax'
  ENDIF
  WRITE(stdout,'(8x,a,f6.2)') 'xemin [eV]: ', xemin
  WRITE(stdout,'(8x,a,f6.2)') 'xemax [eV]: ', xemax
  WRITE(stdout,'(8x,a,i4)')   'xnepoint: ', xnepoint
  
  IF (abs(xe0-xe0_default)<1.d-3) THEN
     WRITE(stdout,'(8x,a,/)') &
          'energy zero automatically set to the Fermi level'
     IF (xonly_plot) THEN
        WRITE(stdout,'(5x,a)') 'Fermi level read in x_save_file'
     ELSE
        WRITE(stdout,'(5x,3a)') &
             'Fermi level determined from SCF save directory (', &
             trim(prefix)//'.save',')'
     ENDIF
     WRITE(stdout,'(5x,a)') &
          'NB: For an insulator (SCF calculated with occupations="fixed")'
     WRITE(stdout,'(5x,a)') &
          '    the Fermi level will be placed at the position of HOMO.'
  ELSE
     WRITE(stdout,'(8x,a,f10.6,3a)') 'xe0 [eV]: ', xe0, &
          ' (energy zero read in ','input file',')'
  ENDIF
  WRITE(stdout,*)
  WRITE(stdout,'(5x,a)') 'WARNING: variable ef_r is obsolete'
  !</DC>
  
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Writing the status of the code (working features and things to do)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  IF(show_status) CALL WRITE_status_of_the_code

end subroutine write_sym_param_to_stdout


subroutine calculate_and_write_homo_lumo_to_stdout(ehomo,elumo)
  USE kinds, ONLY : DP
  USE io_files,      ONLY : prefix
  USE xspectra
  USE constants,       ONLY : rytoev
  USE lsda_mod,        ONLY :lsda
  USE ener,            ONLY : ef, ef_up, ef_dw 
  USE io_global,       ONLY : stdout

  IMPLICIT NONE
!
  REAL(DP) :: ehomo, elumo
  ehomo=0.d0
  elumo=0.d0
  CALL get_homo_lumo(ehomo,elumo)
  ehomo = ehomo*rytoev
  elumo = elumo*rytoev


  WRITE(stdout,1000) ! return+line
  WRITE(stdout,'(5x,a)') & 
       '                      Getting the Fermi energy '
  WRITE(stdout,1001) ! line+return
  !
  IF (lsda) THEN
     WRITE(stdout,'(5x,a,a)') 'From SCF save directory',&
          ' (spin polarized work):' 
     !
     IF (abs(ehomo)<1.e+6) THEN ! insulator => HOMO exists
        WRITE(stdout,'(8x,a,f9.4,a)') 'ehomo [eV]: ', ehomo,&
             ' (highest occupied level:max of up and down)'
        ef=ehomo
        IF (abs(elumo)<1.e+6) THEN  ! insulator and LUMO exists 
           WRITE(stdout,'(8x,a,f9.4,a)') 'elumo [eV]: ', elumo,&
                ' (lowest occupied level:min of up and down)'
        ELSE
           WRITE(stdout,'(8x,a)') 'No LUMO values in SCF calculation'
        ENDIF
     ELSE IF (abs(ef)>1.e-4) THEN
        ef=ef*rytoev !ef in eV
     ELSE
        WRITE(stdout,'(8x,a,f9.4)') 'ef_up [eV]: ', ef_up*rytoeV
        WRITE(stdout,'(8x,a,f9.4)') 'ef_dw [eV]: ', ef_dw*rytoeV
        ef=max(ef_dw,ef_up)*rytoeV
        WRITE(stdout,'(8x,a,f9.4)') '-> ef set to the max of ef_up and ef_dw '
     ENDIF
     WRITE(stdout,'(8x,a,f9.4)') 'ef    [eV]: ', ef   
     
     WRITE(stdout,'(/,5x,a)') &
          '-> ef (in eV) will be written in x_save_file'
  ELSE
     WRITE(stdout,'(5x,a)') 'From SCF save directory:'
     !
     ef=ef*rytoev
     IF (abs(ehomo)<1.e+6) THEN ! insulator => HOMO exists
        WRITE(stdout,'(8x,a,f9.4,a)') 'ehomo [eV]: ', ehomo,&
             ' (highest occupied level)'
        ef=ehomo
        IF (abs(elumo)<1.e+6) THEN  ! insulator and LUMO exists 
           WRITE(stdout,'(8x,a,f9.4,a)') 'elumo [eV]: ', elumo,&
                ' (lowest occupied level)'
        ELSE
           WRITE(stdout,'(8x,a)') 'No LUMO value in SCF calculation'
        ENDIF
     ENDIF
     WRITE(stdout,'(8x,a,f9.4)') 'ef    [eV]: ', ef   
     
     WRITE(stdout,'(/,5x,a)') &
          '-> ef (in eV) will be written in x_save_file'
     
  ENDIF
  
  !
  WRITE(stdout,1000) ! return+line
  WRITE(stdout,'(5x,a)') &
       '                      Energy zero of the spectrum '
  WRITE(stdout,1001) ! line+return
  !
  IF (abs(xe0-xe0_default)<1.d-3) THEN ! no xe0 in input
     WRITE(stdout,'(5x,a)') &
          '-> ef will be used as energy zero of the spectrum'
  ELSE
     WRITE(stdout,'(5x,a,/,7x,3a)') &  
          '-> ef will NOT be used as energy zero of the spectrum',&
          '(because xe0 read in ', 'input file',')'  
  ENDIF

1000 FORMAT(/,5x,&
  '-------------------------------------------------------------------------')
1001 FORMAT(5x,&
  '-------------------------------------------------------------------------',&
  /)
  
end subroutine calculate_and_write_homo_lumo_to_stdout

subroutine write_calculation_type(xang_mom, nl_init)
  USE io_global,       ONLY : stdout
  !internal
  implicit none
  integer, intent(in) :: xang_mom
  integer, intent(in), dimension(2) :: nl_init
  
  WRITE(stdout, 1000) ! line 
  WRITE(stdout,'(5x,a)')&
       '                     Starting XANES calculation'
  IF(nl_init(2).eq.0) then
     IF (xang_mom==1) WRITE(stdout,'(5x,a)')&
          '                in the electric dipole approximation'
     IF (xang_mom==2) WRITE(stdout,'(5x,a)')&
          '              in the electric quadrupole approximation'
     WRITE(stdout,1001)  ! line 
  ELSEIF(nl_init(2).eq.1) then
     WRITE(stdout,'(5x,a)')&
          '                in the electric dipole approximation'
  ENDIF
  
  !...  Writes information about the method
  
  WRITE(stdout,'(7(5x,a,/))') &
       "Method of calculation based on the Lanczos recursion algorithm",&
       "--------------------------------------------------------------",&
       "   - STEP 1: Construction of a kpoint-dependent Lanczos basis,",&
       "     in which the Hamiltonian is tridiagonal (each 'iter' ",&
       "     corresponds to the calculation of one more Lanczos vector)",&
       "   - STEP 2: Calculation of the cross-section as a continued fraction",&
       "     averaged over the k-points."
  
  WRITE(stdout,'(5x,"... Begin STEP 1 ...",/)')
  
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Formats 
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
1000 FORMAT(/,5x,&
          '-------------------------------------------------------------------------')
1001 FORMAT(5x,&
          '-------------------------------------------------------------------------',&
          /)
  
end subroutine write_calculation_type


subroutine write_report_cut_occ_states(cut_occ_states, e_core)
  USE kinds, ONLY : DP
  USE io_global,       ONLY : stdout
  USE xspectra,  ONLY: xemin, xemax,xnepoint, xgamma
  USE gamma_variable_mod, ONLY : gamma_value, gamma_energy, &
       gamma_lines, gamma_tab, gamma_points, &
       gamma_mode, gamma_file
  !internal
  implicit none
  LOGICAL, INTENT(in) :: cut_occ_states
  REAL (DP) :: e_core

  IF (cut_occ_states) THEN
     WRITE(stdout,'(8x,a)') 'the occupied states are elimintate from the spectrum'
  ELSE
     WRITE(stdout,'(8x,a)') 'the occupied states are NOT eliminated from the spectrum'
  ENDIF
  WRITE(stdout,'(8x,a,f6.2)') 'xemin [eV]: ', xemin
  WRITE(stdout,'(8x,a,f6.2)') 'xemax [eV]: ', xemax
  WRITE(stdout,'(8x,a,i4)')   'xnepoint: ', xnepoint
  IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
     WRITE(stdout,'(8x,a,f8.3)')'constant broadening parameter [eV]: ', xgamma
  ELSE
     WRITE(stdout,'(8x,a)') 'energy-dependent broadening parameter:'
     IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
        WRITE(stdout,'(8x,a,a30)')' -> using gamma_file: ', gamma_file
     ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
        WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)')                     &
             ' -> first, constant up to point (', gamma_energy(1), &
             ',', gamma_value(1), ') [eV]'
        WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)')                     &
             ' -> then, linear up to point (', gamma_energy(2),   &
             ',', gamma_value(2), ') [eV]'
        WRITE(stdout,'(8x,a)') ' -> finally, constant up to xemax'
     ENDIF
  ENDIF
  WRITE(stdout,'(8x,"Core level energy [eV]:",1x,g11.4)') -e_core
  WRITE(stdout,'(8x,a,/)') &
       ' (from electron binding energy of neutral atoms in X-ray data booklet)'
  
end subroutine write_report_cut_occ_states
