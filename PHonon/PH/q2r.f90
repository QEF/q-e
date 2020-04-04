!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM q2r
  !----------------------------------------------------------------------------
  !
  !  q2r.x:
  !     reads force constant matrices C(q) produced by the phonon code
  !     for a grid of q-points, calculates the corresponding set of
  !     interatomic force constants (IFC), C(R)
  !
  !  Input data: Namelist "input"
  !     fildyn     :  input file name (character, must be specified)
  !                   "fildyn"0 contains information on the q-point grid
  !                   "fildyn"1-N contain force constants C_n = C(q_n)
  !                   for n=1,...N, where N is the number of q-points
  !                   in the irreducible brillouin zone
  !                   Normally this should be the same as specified
  !                   on input to the phonon code
  !                   In the non collinear/spin-orbit case the files 
  !                   produced by ph.x are in .xml format. In this case 
  !                   fildyn is the same as in the phonon code + the .xml 
  !                   extension.
  !     flfrc      :  output file containing the IFC in real space
  !                   (character, must be specified)
  !     zasr       :  Indicates type of Acoustic Sum Rules used for the Born
  !                   effective charges (character):
  !                   - 'no': no Acoustic Sum Rules imposed (default)
  !                   - 'simple':  previous implementation of the asr used
  !                     (3 translational asr imposed by correction of
  !                     the diagonal elements of the force-constants matrix)
  !                   - 'crystal': 3 translational asr imposed by optimized
  !                      correction of the IFC (projection).
  !                   - 'one-dim': 3 translational asr + 1 rotational asr
  !                     imposed by optimized correction of the IFC (the
  !                     rotation axis is the direction of periodicity; it
  !                     will work only if this axis considered is one of
  !                     the cartesian axis).
  !                   - 'zero-dim': 3 translational asr + 3 rotational asr
  !                     imposed by optimized correction of the IFC.
  !                   Note that in certain cases, not all the rotational asr
  !                   can be applied (e.g. if there are only 2 atoms in a
  !                   molecule or if all the atoms are aligned, etc.).
  !                   In these cases the supplementary asr are cancelled
  !                   during the orthonormalization procedure (see below).
  !     loto_2d    :  set to .true. to activate two-dimensional treatment of LO-TO splitting. 
  !
  !  If a file "fildyn"0 is not found, the code will ignore variable "fildyn"
  !  and will try to read from the following cards the missing information
  !  on the q-point grid and file names:
  !     nr1,nr2,nr3:  dimensions of the FFT grid formed by the q-point grid
  !     nfile      :  number of files containing C(q_n), n=1,nfile
  !  followed by nfile cards:
  !     filin      :  name of file containing C(q_n)
  !  The name and order of files is not important as long as q=0 is the first
  !
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE io_global,   ONLY : ionode_id, ionode, stdout
  USE environment, ONLY : environment_start, environment_end
  USE el_phon,     ONLY : el_ph_nsigma
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: fildyn, filin, flfrc, prefix
  CHARACTER (LEN=10) :: zasr
  LOGICAL            :: la2F, loto_2d
  INTEGER            :: ios
  !
  NAMELIST / input / fildyn, flfrc, prefix, zasr, la2F, loto_2d, el_ph_nsigma
  !
  CALL mp_startup()
  CALL environment_start('Q2R')
  !
  IF (ionode) CALL input_from_file ( )
     !
  fildyn = ' '
  flfrc = ' '
  prefix = ' '
  loto_2d=.false.
  zasr = 'no'
     !
  la2F=.false.
  el_ph_nsigma=10
     !
     !
  IF (ionode)  READ ( 5, input, IOSTAT =ios )

  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('q2r','error reading input namelist', abs(ios))

  CALL mp_bcast(fildyn, ionode_id, world_comm)
  CALL mp_bcast(flfrc, ionode_id, world_comm)
  CALL mp_bcast(prefix, ionode_id, world_comm)
  CALL mp_bcast(zasr, ionode_id, world_comm)
  CALL mp_bcast(loto_2d, ionode_id, world_comm)
  CALL mp_bcast(la2F, ionode_id, world_comm)
  CALL mp_bcast(el_ph_nsigma, ionode_id, world_comm)
  !
  CALL do_q2r(fildyn, flfrc, prefix, zasr, la2F, loto_2d)
  !
  CALL environment_end('Q2R')

  CALL mp_global_end()
  !
END PROGRAM q2r
!
!----------------------------------------------------------------------------
