! Copyright (C) 2001-2003 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!----------------------------------------------------------------------- 
program initial_state
  !----------------------------------------------------------------------- 
  ! 
  !  compute forces on atoms as a post-process
  ! 
  ! input: namelist "&inputpp", with variables 
  !   prefix      prefix of input files saved by program pwscf 
  !   outdir      temporary directory where files resides 
  ! 
  USE io_global,  ONLY : stdout 
  USE kinds,      ONLY : DP 
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir, iunwfc, nwordwfc
  USE ions_base,  ONLY : nat
  USE klist,      ONLY : nks, xk
  USE wvfct,      ONLY : npw, igk
  USE uspp,       ONLY : nkb, vkb
  USE wavefunctions_module, ONLY : evc
  USE parameters, ONLY : ntypx
#ifdef __PARA 
  use para,       only : me 
  use mp 
#endif 
  implicit none 
  character(len=256) :: outdir 
  integer :: ios, ionode_id = 0 , ik, excite(ntypx)
  namelist / inputpp / outdir, prefix, excite
  ! 
  call start_postproc (nd_nmbr) 
  ! 
  !   set default values for variables in namelist 
  ! 
  excite(:) = 0
  prefix = 'pwscf' 
  outdir = './' 
  ! 
#ifdef __PARA 
  if (me == 1)  then 
#endif 
  read (5, inputpp, err = 200, iostat = ios) 
200 call errore ('postforces', 'reading inputpp namelist', abs (ios) ) 
  ! 
  tmp_dir = trim(outdir) 
  ! 
#ifdef __PARA 
  end if 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast( excite, ionode_id )
#endif 
  ! 
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  call read_file 
  call openfil_pp
  call hinit0 
  call hinit1 
  IF ( nks == 1 ) THEN
     ik = 1
     CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
     IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
  END IF
  call sum_band 
  ! 
  call do_initial_state (excite)
  ! 
  call stop_pp 
  ! 
end program initial_state
 
