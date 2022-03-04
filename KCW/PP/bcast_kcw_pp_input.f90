!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_kcw_pp_input ( )
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the acfdt input to all
  !     the other processors
  !
  !
#ifdef __MPI

  USE run_info,      ONLY : title
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_bcast
  USE mp_global,     ONLY : intra_image_comm
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_kcw
                            
  !
  IMPLICIT NONE
  !
  call mp_bcast ( title,         ionode_id, intra_image_comm ) 
  call mp_bcast ( tmp_dir,       ionode_id, intra_image_comm )
  call mp_bcast ( prefix,        ionode_id, intra_image_comm )
  !
  call mp_bcast ( seedname,            ionode_id, intra_image_comm )
  call mp_bcast ( num_wann,            ionode_id, intra_image_comm )
  call mp_bcast ( kcw_iverbosity,       ionode_id, intra_image_comm )
  call mp_bcast ( mp1,                 ionode_id, intra_image_comm )
  call mp_bcast ( mp2,                 ionode_id, intra_image_comm )
  call mp_bcast ( mp3,                 ionode_id, intra_image_comm )
  call mp_bcast ( use_ws_distance,     ionode_id, intra_image_comm )
  call mp_bcast ( have_empty     ,     ionode_id, intra_image_comm )
   !
#endif
  !
  return
  !
end subroutine bcast_kcw_pp_input
