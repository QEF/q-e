!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
MODULE formats
  !---------------------------------------------------------------------------
  !
  ! ... this module contains the I/O formats used by all NEB-routines
  !
  CHARACTER (LEN=*), PARAMETER :: &
    lattice_vectors  = "(3(2X,F14.10),/,3(2X,F14.10),/,3(2X,F14.10))"
  !
  CHARACTER (LEN=*), PARAMETER :: &
    restart_first  = "(3(2X,F14.10),3(2X,F14.10),3(2X,I1))", &
    restart_others = "(3(2X,F14.10),3(2X,F14.10))"
  !  
  CHARACTER (LEN=*), PARAMETER :: &  
    velocities = "(3(2X,F14.10))", &
    energy     = "(2X,F16.8)"
  !  
  CHARACTER (LEN=*), PARAMETER :: &  
    dat_fmt  = "(3(2X,F16.10))",    &
    int_fmt  = "(2(2X,F16.10))",    &
    xyz_fmt  = "(A2,3(2X,F14.10))", &
    axsf_fmt = "(A2,6(2X,F14.10))"
  !  
  CHARACTER (LEN=*), PARAMETER :: &  
    scf_fmt            = "(5X,'tcpu = ',F8.2," // &
                       & "';  self-consistency for image ', I3)", &
    scf_fmt_para       = "(5X,'cpu = ',I2,'; tcpu = ',F8.2," // &
                       & "';  self-consistency for image ', I3)", &                   
    run_output         = "(/,5X,'iteration:',I4,4X,'E activation ='," // &
                       & " F6.3,' eV',4X,'error =',F8.4,' eV / bohr'/)", &
    run_output_T_const = "(/,5X,'iteration:',I4,4X,'temperature ='," // &
                       & " F8.2,' K',4X,'forces =',F8.4,' eV / bohr')", &  
    final_output       = "(5X,'image: ',I2,'   E tot = ',F16.8," // &
                       & "' eV   error = ',F8.4,' eV / bohr')"
  !
  CHARACTER (LEN=*), PARAMETER :: &
    stringfmt = "(5X,A,T35,' = ',A)"      
  !
END MODULE formats
