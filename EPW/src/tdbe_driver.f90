  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE tdbe_driver
  !---------------------------------------------------------------------  
  !!
  !! This module contains subroutines controlling TDBE simulation.
  !! Authored by Yiming Pan (pan@physik.uni-kiel.de) 
  !! and Fabio Caruso (caruso@physik.uni-kiel.de) (Apr. 2025)
  !! based on the equations derived in the following works:
  !! [1] Adv. Phys: X 7, 2095925 (2022) 
  !! [2] Nano Lett. 23, 7463 (2023)
  !
  IMPLICIT NONE
  !
  CONTAINS  
    !
    !--------------------------------------------------------------------- 
    SUBROUTINE tdbe()
    !---------------------------------------------------------------------  
    !!
    !---------------------------------------------------------------------  
    USE input,            ONLY : temp_ph_tdbe, restart_tdbe,  phph_tdbe, mp_mesh_k    
    USE tdbe_mod,         ONLY : map2nkfs, double_grid_map,read_g2_tdbe,            &
                                 deallocate_tdbe, allocate_tdbe,                    &
                                 wanninterp, restart_tdbe_read,                      &
                                 set_initial_dist, write_dist, write_points,        &
                                 double_grid_map, read_g2_tdbe, map2nkfs, check_dg, &
                                 time_propagation, deallocate_tdbe   
    USE phph,             ONLY : psi2_calc,  inv_tau_rta, ph_energy, read_ph_lw,    &
                                 deallocate_phph
    !
    IMPLICIT NONE
    !
    !
    CALL start_clock('TDBE')
    CALL wanninterp()
    !! Allocate space for time-dependent electron and phonon distribution, and collision inegrals 
    CALL allocate_tdbe() 
    !!
    IF (restart_tdbe) THEN
      CALL restart_tdbe_read()  
    ELSE  
      CALL set_initial_dist() 
      CALL write_dist(.TRUE., 0.0d0)
    ENDIF
    !!
    !!
    CALL write_points()
    CALL double_grid_map()
    CALL read_g2_tdbe()
    CALL map2nkfs() 
    IF (mp_mesh_k) THEN
      CALL check_dg() 
    ENDIF
    !!
    IF (phph_tdbe) THEN
      IF (.FALSE.) THEN 
        CALL ph_energy()
        CALL read_ph_lw()
      ELSE  
        CALL psi2_calc()
        CALL inv_tau_rta(temp_ph_tdbe, temp_ph_tdbe + 600.d0, 6)
      ENDIF
    ENDIF
    !! Now start the tdBE calulation
    !!
    CALL time_propagation()
    CALL deallocate_tdbe()
    IF (phph_tdbe) THEN
      CALL deallocate_phph()
    ENDIF
    CALL stop_clock('TDBE')
    !
    !---------------------------------------------------------------------  
    END SUBROUTINE tdbe
    !---------------------------------------------------------------------   
  !-----------------------------------------------------------------------
  END MODULE tdbe_driver 
  !-----------------------------------------------------------------------
