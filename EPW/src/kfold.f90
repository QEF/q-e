  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  module kfold
  !-----------------------------------------------------------------------
  !
  ! ... The variables for folding the k+q gris back into the k grid
  ! 
  USE kinds,      ONLY : DP 
  !
  INTEGER :: &
!       g0vec_all(3,27),&! G-vectors needed to fold the k+q grid into the k grid
       g0vec_all(3,125),&! G-vectors needed to fold the k+q grid into the k grid 
       ng0vec           ! number of inequivalent such translations (27)
  INTEGER, allocatable :: &
       shift (:),      &! for every k+q, index of the G0-vector needed to fold k+q into k+q+G0 
       gmap(:,:)        ! the map G --> G-G_0 in the large (density) set, for every G_0 (27 at most)
  REAL(kind=DP) :: &
!       g0vec_all_r(3,27)! G-vectors needed to fold the k+q grid into the k grid, cartesian coord.
       g0vec_all_r(3,125)! G-vectors needed to fold the k+q grid into the k grid, cartesian coord. - LSCO fix
  !
  end module kfold

