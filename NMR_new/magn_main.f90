!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM magn_main
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the magnetic resposne program. 
  ! ... It controls the initialization routines.
  ! ... Features: NMR chemical shifts
  !...            EPR g-tensor
  ! ... Ported to Espresso by:
  ! ... D. Ceresoli, A. P. Seitsonen, U. Gerstamnn and  F. Mauri
  ! ...
  ! ... References (NMR):
  ! ... F. Mauri and S. G. Louie Phys. Rev. Lett. 76, 4246 (1996)
  ! ... F. Mauri, B. G. Pfrommer, S. G. Louie, Phys. Rev. Lett. 77, 5300 (1996)
  ! ... T. Gregor, F. Mauri, and R. Car, J. Chem. Phys. 111, 1815 (1999)
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. B 63, 245101 (2001)
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. Lett. 91, 196401 (2003)
  ! ...
  ! ... References (g-tensor):
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. Lett. 88, 086403 (2002)
  ! ...
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, ionode, ionode_id
  USE wvfct,           ONLY : gamma_only
  USE io_files,        ONLY : prefix, tmp_dir, nd_nmbr
  USE klist,           ONLY : nks
  USE mp,              ONLY : mp_bcast
  USE cell_base,       ONLY : tpiba
  USE nmr_module
  USE global_version,  ONLY : version_number
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   :: code = 'MAGN'
  !------------------------------------------------------------------------


  CALL init_clocks( .TRUE. )
  CALL start_clock( 'MAGN' )

  ! ... and begin with the initialization part
  CALL startup( nd_nmbr, code, version_number )
  CALL nmr_readin()

  ! read ground state wavefunctions
  call read_file
  call openfil

  if (gamma_only) call errore('MAGN_main', 'gamma_only == .true.', 1)

  CALL nmr_allocate()
  CALL nmr_setup()
  CALL nmr_summary()
  CALL nmr_openfil()
  CALL print_clock( 'MAGN' )

  ! convert q_nmr in units of tpiba
  q_nmr = q_nmr / tpiba

  ! calculation
  select case(trim(job))
    case('nmr')
      call suscept_crystal   

    case('g_tensor')
      call g_tensor_crystal

    case('f-sum')
      call test_f_sum_rule 

    case default
      call errore('magn_main', 'wrong or undefined job in input', 1)
  end select

  ! print timings and stop the code
  CALL print_clock_nmr()
  CALL stop_code( .TRUE. )
  STOP
END PROGRAM magn_main
