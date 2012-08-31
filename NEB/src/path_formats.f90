!
! Copyright (C) 2003-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE path_formats
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the I/O formats used by all "path"-routines
  !
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
  !
  CHARACTER (LEN=*), PARAMETER :: &
    lattice_vectors  = "(3(2X,F14.10),/,3(2X,F14.10),/,3(2X,F14.10))"
  !
  CHARACTER (LEN=*), PARAMETER :: &
    restart_first  = "(3(2X,F18.12),3(2X,F18.12),3(2X,I1))", &
    restart_others = "(3(2X,F18.12),3(2X,F18.12))"
  !  
  CHARACTER (LEN=*), PARAMETER :: &  
    quick_min = "(9(2X,F18.12))", &
    energy    = "(2X,F18.10)"
  !  
  CHARACTER (LEN=*), PARAMETER :: &  
    dat_fmt  = "(3(2X,F16.10))",    &
    int_fmt  = "(2(2X,F16.10))",    &
    xyz_fmt  = "(A2,3(2X,F14.10))", &
    axsf_fmt = "(A2,6(2X,F14.10))"
  !  
  CHARACTER (LEN=*), PARAMETER :: &
    scf_iter_fmt = "(/,5X,30('-'),(1X,'iteration ',I3,1X),30('-'),/)", &
    scf_fmt      = "(5X,'tcpu = ',F8.1," // &
                 & "'    self-consistency for image ', I3)", &
    scf_fmt_para = "(5X,'cpu = ',I2,'   tcpu = ',F8.1," // &
                 & "'    self-consistency for image ', I3)"
  !
  CHARACTER (LEN=*), PARAMETER :: &
    run_info = "(5X,'image',8X,'energy (eV)',8X,'error (eV/A)',8X,'frozen',/)"
  !
  CHARACTER (LEN=*), PARAMETER :: &
    run_output = "(5X,I5,4X,F15.7,10X,F10.6,12X,L1)"
  !    
  CHARACTER (LEN=*), PARAMETER :: &
    summary_fmt = "(5X,A,T35,' = ',3X,A)"
  !
  CHARACTER (LEN=*), PARAMETER :: &
    final_fmt = "(/,5X,75('-'),/)"
  !
END MODULE path_formats
