!
! Copyright (C) 2001-2005 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE read_input_vib()
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program cpvib.
  !    from prefix.vib.inp input file.
  !
  ! ... modules
  !
  USE constants,        ONLY : AMU_AU
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE ions_base,        ONLY : nat, ityp, pmass, amass
  USE input_parameters, ONLY : atom_label
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast
  USE vibrations,       ONLY : displacement, save_freq,                            &
       trans_inv_max_iter, trans_inv_conv_thr, vib_restart_mode,             &
       trans_inv_flag, trans_rot_inv_flag, animate, isotope
#ifdef DFT_CP
  USE io_files,         ONLY : outdir, prefix
#endif

#ifdef DFT_PW
  USE io_files,         ONLY : prefix, outdir=>tmp_dir
#endif
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER             :: ios, dirlen, i
  CHARACTER (LEN=256) :: restart_label, input_filename
  CHARACTER (LEN=3)   :: tmp_str
  LOGICAL             :: inp_file_exists
  !
  NAMELIST / INPUT_VIB / displacement, save_freq, trans_inv_max_iter,        &
       trans_inv_conv_thr, restart_label, trans_inv_flag,trans_rot_inv_flag, &
       animate
  NAMELIST / ISOTOPES / isotope
  !
  IF ( ionode ) THEN
     !
     ! Set input file name
     !
     dirlen         = INDEX(outdir,' ') - 1
     input_filename = './' // TRIM(prefix) // '.vib.inp'
     input_filename = outdir(1:dirlen) // '/' // input_filename
     WRITE (stdout,*) 'Input information is read from:', input_filename
     !
     ! ... Open input file
     !
     INQUIRE (FILE = input_filename, EXIST = inp_file_exists)
     IF ( .NOT. inp_file_exists) &
          CALL errore( 'read_input_vib', 'input file absent ', 1 )
     OPEN( 99 , FILE = input_filename , STATUS = 'old' )
     !
     ! ... set default values for variables in namelist INPUT_VIB
     !
     displacement       = 0.05
     save_freq          = 1
     trans_inv_flag     = .TRUE.
     trans_rot_inv_flag = .TRUE.
     trans_inv_max_iter = 50
     trans_inv_conv_thr = 1.0D-15
     vib_restart_mode   = 2
     restart_label      = 'auto'
     animate            = .FALSE.
     !
     ! ... set default values for variables in namelist ISOTOPES
     !
#ifdef DFT_PW
     pmass = amass * AMU_AU
#endif
     do i=1,nat
        isotope(i) = pmass(ityp(i)) / AMU_AU
     end do
     !
     ! ...  reading the namelist input_vib
     !
     READ( 99, input_vib, ERR = 200, IOSTAT = ios )
     !
200  CALL errore( 'readin_input_vib', 'reading input_vib namelist',  ABS(ios) )
     !
     ! ...  reading the namelist isotopes
     !
     READ( 99, isotopes, ERR = 201, IOSTAT = ios )
     !
201  CALL errore( 'readin_input_vib', 'reading isotopes namelist',  ABS(ios) )
     !
     ! ... Check all namelist variables
     !
     IF (displacement <= 0.D0) &
          CALL errore (' readin_input_vib', ' Wrong displacement        ', 1)
     !
     IF (trans_inv_max_iter <= 0.D0) &
          CALL errore (' readin_input_vib', ' Wrong trans_inv_max_iter ', 1)
     !
     IF (trans_inv_conv_thr <= 0.D0) &
          CALL errore (' readin_input_vib', ' Wrong trans_inv_conv_thr ', 1)
     !
     SELECT CASE ( TRIM( restart_label ) )
     CASE( 'from_scratch' )
        vib_restart_mode = 0
     CASE ( 'restart'     )
        vib_restart_mode = 1
     CASE ( 'auto'        )
        vib_restart_mode = 2
     CASE DEFAULT
        CALL errore( 'read_input_vib ', &
             'unknown vib_restart_mode ' // TRIM( restart_label ), 1 )
     END SELECT
     !
     do i=1,nat
        if (isotope(i).le.0.0) then
           write (tmp_str,'(I3.3)') i
           CALL errore (' readin_input_vib', &
                ' negative or zero mass for atom '//tmp_str, 1)
        end if
     end do

  END IF
  !
  do i=1,nat
     isotope(i) = isotope(i) * AMU_AU
  end do
  !
  ! ... Broadcasting input variables
  !
  CALL mp_bcast( displacement,       ionode_id )
  CALL mp_bcast( save_freq,          ionode_id )
  CALL mp_bcast( trans_inv_max_iter, ionode_id )
  CALL mp_bcast( trans_inv_conv_thr, ionode_id )
  CALL mp_bcast( trans_inv_flag,     ionode_id )
  CALL mp_bcast( trans_rot_inv_flag, ionode_id )
  CALL mp_bcast( vib_restart_mode,   ionode_id ) 
  CALL mp_bcast( animate,            ionode_id )
  !
  ! ... Printing input data to stdout
  !
  IF (ionode) THEN
     WRITE (stdout,*) '------------------------------------------------------'
     WRITE (stdout,*) ' SUMMERY OF VIBRATIONAL INPUT VARIABLES: '
     WRITE (stdout,*) ''
     WRITE (stdout,101) displacement
     WRITE (stdout,102) save_freq
     WRITE (stdout,103) trans_inv_flag
     WRITE (stdout,104) trans_inv_max_iter
     WRITE (stdout,105) trans_inv_conv_thr
     WRITE (stdout,106) trans_rot_inv_flag
     WRITE (stdout,107) restart_label
     WRITE (stdout,108) animate
     WRITE (stdout,*) ''
     do i=1,nat
        WRITE (stdout,109) i,TRIM(atom_label(ityp(i))),isotope(i) / AMU_AU
     end do


     WRITE (stdout,*) '------------------------------------------------------'
     !
101  FORMAT(3x,'displacement           = ',3x,f10.4)
102  FORMAT(3x,'save_freq              = ',3x,I10)
103  FORMAT(3x,'trans_inv_flag         = ',3x,L10)
104  FORMAT(3x,'trans_inv_max_iter     = ',3x,I10)
105  FORMAT(3x,'trans_inv_conv_thr     = ',3x,D10.4)
106  FORMAT(3x,'trans_rot_inv_flag     = ',3x,L10)
107  FORMAT(3x,'restart_label          = ',3x,A10)
108  FORMAT(3x,'animate                = ',3x,L10)
109  FORMAT(3x,I3.3,3x,A5,F10.4)
     !
     CLOSE(99)
  END IF
  !
  RETURN
  !
END SUBROUTINE read_input_vib
