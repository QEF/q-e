!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE io_files
!=----------------------------------------------------------------------------=!

     USE parameters, ONLY: ntypx

  !
  ! ... The name of the files
  !

     IMPLICIT NONE
     SAVE

     CHARACTER(len=80) :: tmp_dir = './'  ! directory for temporary files
     CHARACTER(len=80) :: prefix  = 'os'  ! prepended to file names
     CHARACTER(len=3)  :: nd_nmbr = '000' ! node number (used only in parallel case)
     CHARACTER(len=80) :: pseudo_dir = './'
     CHARACTER(len=80) :: psfile( ntypx ) = 'UPF'
     CHARACTER(len=80) :: scradir = './'
     CHARACTER(len=80) :: outdir  = './'

     CHARACTER(LEN=80) :: filpun = ' '       ! name of the punch file
     CHARACTER(LEN=80) :: input_drho = ' '   ! name of the file with the input drho
     CHARACTER(LEN=80) :: output_drho = ' '  ! name of the file with the output drho

     CHARACTER(LEN=80) :: band_file = ' '
     CHARACTER(LEN=80) :: tran_file = ' '
     CHARACTER(LEN=80) :: fil_loc = ' '      !  file with 2D eigenvectors and eigenvalues

     CHARACTER(LEN=14), PARAMETER :: rho_name      = 'CHARGE_DENSITY'
     CHARACTER(LEN=17), PARAMETER :: rho_name_up   = 'CHARGE_DENSITY.UP'
     CHARACTER(LEN=19), PARAMETER :: rho_name_down = 'CHARGE_DENSITY.DOWN'
     CHARACTER(LEN=14), PARAMETER :: rho_name_avg  = 'CHARGE_AVERAGE'

     CHARACTER(LEN=15), PARAMETER :: empty_file    = 'EMPTY_STATES.WF'
     CHARACTER(LEN=5 ), PARAMETER :: crash_file    = 'CRASH'
     CHARACTER(LEN=7 ), PARAMETER :: stop_file     = '.cpstop'
     CHARACTER(LEN=2 ), PARAMETER :: ks_file       = 'KS'
     CHARACTER(LEN=6 ), PARAMETER :: ks_emp_file   = 'KS_EMP'
     CHARACTER(LEN=16), PARAMETER :: sfac_file     = 'STRUCTURE_FACTOR'

  !
  ! ... The units where various variables are saved
  !

     INTEGER :: rhounit = 17
     INTEGER :: emptyunit = 19
     INTEGER :: crashunit = 15
     INTEGER :: stopunit = 7
     INTEGER :: ksunit = 18
     INTEGER :: sfacunit = 20
     INTEGER :: pseudounit = 10

     INTEGER :: iunpun           ! unit for saving the final results
     INTEGER :: iunwfc           ! unit with wavefunctions
     INTEGER :: iunat            ! unit for saving orthogonal atomic wfcs
     INTEGER :: iunocc           ! unit for saving the atomic n_{ij}
     INTEGER :: iunoldwfc        ! unit with old wavefunctions (molecular dynamics)
     INTEGER :: iunoldwfc2       ! as above at step -2
     INTEGER :: iunigk           ! unit for saving indices
     INTEGER :: iunres           ! unit for the restart of the run
     INTEGER :: iunneb           ! unit for NEB output ( sdtout or what else )
     INTEGER :: nwordwfc         ! lenght of record in wavefunction file
     INTEGER :: nwordatwfc       ! lenght of record in atomic wfc file

!=----------------------------------------------------------------------------=!
   END MODULE io_files
!=----------------------------------------------------------------------------=!
