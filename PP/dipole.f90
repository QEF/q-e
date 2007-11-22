!
! Copyright (C) 2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
PROGRAM dipole
  !-----------------------------------------------------------------------
  !
  !      Calculate dipole in isolated systems, e.g. molecule in supercell
  !      The center of ionic charges will be used as origin
  !      Calculation of Makov-Payne correction for charged supercells
  !      as in : G. Makov and M.C. Payne, PRB 51, 4014 (1995)
  !      (note that Eq. 15 has the wrong sign for the quadrupole term)
  !      Contributed by Giovanni Cantele, Paolo Cazzato, Carlo Sbraccia
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi, rytoev
  USE io_global, ONLY : ionode, ionode_id,  stdout
  USE io_files,  ONLY : prefix, outdir, tmp_dir, trimcheck, nd_nmbr
  USE ener,      ONLY : etot
  USE ions_base, ONLY : nsp, nat, tau, ityp, zv
  USE cell_base, ONLY : at, bg, omega, alat, ibrav
  USE gvect,     ONLY : g, gg, ngm, gstart, igtongl
  USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE lsda_mod,  ONLY : nspin
  USE pfft,      ONLY : npp
  USE mp_global, ONLY : me_pool, intra_pool_comm
  USE vlocal,    ONLY : strf, vloc
  USE mp,        ONLY : mp_sum, mp_bcast
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  CHARACTER(len=256) :: filplot
  INTEGER :: ios
  NAMELIST / inputpp / outdir, prefix
  !
  ! initialise parallel environment
  !
  CALL start_postproc (nd_nmbr)
  IF ( ionode )  CALL input_from_file ( )
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     IF ( ios /= 0 ) &
        CALL infomsg ('dipole', 'reading inputpp namelist')
     !
     tmp_dir = trimcheck ( outdir )
     !
  END IF
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  !
  CALL read_file ( )
  !
  call makov_payne (etot)
  !
  call stop_pp()
  !
END PROGRAM dipole
