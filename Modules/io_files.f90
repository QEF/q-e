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
     CHARACTER (LEN=80) :: &
       dat_file  = 'os.dat',   &! file containing the enegy profile
       int_file  = 'os.int',   &! file containing the interpolated energy profile
       neb_file  = 'os.neb',   &! file containing informations needed to restart a neb simulation
       xyz_file  = 'os.xyz',   &! file containing coordinates of all images in xyz format
       axsf_file = 'os.axsf'    ! file containing coordinates of all images in axsf format
     CHARACTER (LEN=4), PARAMETER :: &
       exit_file = "EXIT"       ! file required for a soft exit  
  !
  ! ... The units where various variables are saved
  !

     INTEGER :: rhounit     = 17
     INTEGER :: emptyunit   = 19
     INTEGER :: crashunit   = 15
     INTEGER :: stopunit    = 7
     INTEGER :: ksunit      = 18
     INTEGER :: sfacunit    = 20
     INTEGER :: pseudounit  = 10

     INTEGER :: iunpun      =  4 ! unit for saving the final results
     INTEGER :: iunwfc      = 10 ! unit with wavefunctions
     INTEGER :: iunat       = 13 ! unit for saving orthogonal atomic wfcs
     INTEGER :: iunocc      = 14 ! unit for saving the atomic n_{ij}
     INTEGER :: iunoldwfc   = 11 ! unit with old wavefunctions
     INTEGER :: iunoldwfc2  = 12 ! as above at step -2
     INTEGER :: iunigk      = 16 ! unit for saving indices
     INTEGER :: iunres      =  1 ! unit for the restart of the run
     INTEGER :: iunbfgs     = 30 ! unit for the bfgs restart file
     !
     INTEGER :: nwordwfc    =  2 ! lenght of record in wavefunction file
     INTEGER :: nwordatwfc  =  2 ! lenght of record in atomic wfc file
     !
     INTEGER :: iunexit     = 26 ! unit for a soft exit  
     INTEGER :: iunupdate   = 27 ! unit for saving old positions (extrapolation)
     INTEGER :: iunpara     = 28 ! unit for parallelization among images
     INTEGER :: iunblock    = 29 ! as above (blocking file)
     !
     ! ... NEB specific
     !
     INTEGER :: iunneb      =  6 ! unit for NEB output ( stdout or what else )
     INTEGER :: iunrestart  = 21 ! unit for saving the restart file ( neb_file )
     INTEGER :: iundat      = 22 ! unit for saving the enegy profile
     INTEGER :: iunint      = 23 ! unit for saving the interpolated energy profile
     INTEGER :: iunxyz      = 24 ! unit for saving coordinates ( xyz format )
     INTEGER :: iunaxsf     = 25 ! unit for saving coordinates ( axsf format )

!=----------------------------------------------------------------------------=!
   END MODULE io_files
!=----------------------------------------------------------------------------=!
