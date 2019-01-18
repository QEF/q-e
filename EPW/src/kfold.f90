  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  MODULE kfold
  !-----------------------------------------------------------------------
  !
  ! ... The variables for folding the k+q gris back into the k grid
  ! 
  USE kinds, ONLY : DP 
  !
  INTEGER :: g0vec_all(3,125)
  !! G-vectors needed to fold the k+q grid into the k grid 
  INTEGER :: ng0vec           
  !! number of inequivalent such translations (125)
  INTEGER, ALLOCATABLE :: shift(:)
  !! for every k+q, index of the G_0-vector needed to fold k+q into k+q+G0 
  INTEGER, ALLOCATABLE :: gmap(:,:)        
  !! the map G --> G-G_0 in the large (density) set, for every G_0 (125 at most)
  !
  REAL(DP) :: g0vec_all_r(3,125)
  !! G-vectors needed to fold the k+q grid into the k grid, cartesian coord. 
  !
  END MODULE kfold
