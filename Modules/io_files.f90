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

     IMPLICIT NONE
     SAVE

     CHARACTER(len=80) :: tmp_dir = './'  ! directory for temporary files
     CHARACTER(len=80) :: prefix  = 'os'  ! prepended to file names
     CHARACTER(len=3)  :: nd_nmbr = '000' ! node number (used only in parallel case)
     CHARACTER(len=80) :: pseudo_dir = './'
     CHARACTER(len=80) :: psfile( ntypx ) = 'UPF'
     CHARACTER(len=80) :: scradir = './'
     CHARACTER(len=80) :: outdir  = './'

!=----------------------------------------------------------------------------=!
   END MODULE io_files
!=----------------------------------------------------------------------------=!
