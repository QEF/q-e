!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine do_postproc (nodenumber)  
  !-----------------------------------------------------------------------
  !
  !    This routine performs various postprocessing steps. The action
  !    is controlled through the following variables in namelist inputpp:
  !
  ! prefix      prefix of files saved by program pwscf
  ! tmp_dir     temporary directory where pwscf files resides
  !
  ! filplot     punch file, contains the quantity selected by plot_num
  ! plot_num    selects what is saved in filplot:
  !                0=charge
  !                1=total potential V_bare+V_H + V_xc
  !                2=local ionic potential
  !                3=local density of states at e_fermi
  !                4=local density of electronic entropy
  !                5=STM images
  !                6=spin polarization (rho(up)-rho(down))
  !                7=|psi|^2
  !                8=electron localization function (ELF)
  !                9=planar average of all |psi|^2
  !               10=integrated local density of states from 
  !                  emin to emax (emin, emax in eV)
  !                  if emax is not specified, emax=E_fermi
  !               11=the V_bare + V_H potential
  !
  !             Options for total charge
  !
  ! spin_component 0=total charge (default value),
  !                1=spin up charge,
  !                2=spin down charge.
  !
  !             Options for total potential
  !
  ! spin_component 0=spin averaged potential (default value),
  !                1=spin up potential,
  !                2=spin down potential.
  !
  !             Options for STM images:
  !
  ! sample_bias    the bias of the sample (Ryd) in stm images
  ! stm_wfc_matching     if .t. match the wavefunctions
  ! z           height of matching (in celldm(3) units)
  ! dz          distance of next stm image calculation
  !
  !             Options for |psi|^2:
  !
  ! kpoint      which k-point
  ! kband       which band
  !
  use pwcom  
  use io  

  implicit none
  character :: nodenumber * 3, filband * 14, filplot * 14

  integer :: n_atom_wfc, plot_num, kpoint, kband, spin_component, ios
  logical :: stm_wfc_matching  

  real(kind=DP) :: emin, emax, sample_bias, z, dz  
  ! directory for temporary files

  namelist / inputpp / tmp_dir, prefix, plot_num, stm_wfc_matching, &
       sample_bias, spin_component, z, dz, emin, emax, kpoint, kband,&
       filplot, filband
  !
  nd_nmbr = nodenumber  
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  tmp_dir = './'  
  filplot = ' '  
  filband = ' '  
  plot_num = 0  
  spin_component = 0  
  sample_bias = 0.01d0  
  z = 1.d0  
  dz = 0.05d0  
  stm_wfc_matching = .true.  
  emin = - 999.0d0  
  emax = ef*13.6058d0
  !
  !     reading the namelist inputpp
  !
  read (5, inputpp, err = 200, iostat = ios)  
200 call error ('postproc', 'reading inputpp namelist', abs (ios) )  
  !
  !     Check of namelist variables
  !
  if (filplot.ne.' ') then  

     if (plot_num.lt.0.or.plot_num.gt.11) call error ('postproc', &
          'Wrong plot_num', abs (plot_num) )

     if ( (plot_num.eq.0.or.plot_num.eq.1) .and. ( &
          spin_component.lt.0.or.spin_component.gt.2) ) call error ( &
          'postproc', 'wrong value of spin_component', 1)

     if (plot_num.eq.10) then
        emin = emin / 13.6058d0  
        emax = emax / 13.6058d0  
     end if
  endif

  if (filplot.eq.' '.and.filband.eq.' ') call &
       error ('postproc', 'nothing to do?', 1)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file  
  call openfil  
  call struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
       strf, eigts1, eigts2, eigts3)
  call init_us_1  
  !
  !   Now do whatever you want
  !
  call punch_plot (filplot, plot_num, sample_bias, z, dz, &
       stm_wfc_matching, emin, emax, kpoint, kband, spin_component)
  call punch_band (filband)  

  return  
end subroutine do_postproc
