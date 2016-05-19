  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE gamma_acont( omega, omegap, temp, rgammap, rgammam )
  !-----------------------------------------------------------------------
  !
  ! computes gammam(w,wp)  (notes RM)
  ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
  !
  ! input
  !
  ! omega  - frequency w at point iw on the real-axis
  ! omegap - frequency w' at point iwp on the real-axis
  ! temp   - temperature in eV
  !
  ! output
  !
  ! rgammap = -bose_einstein( w' ) - fermi_dirac(  w + w' )
  ! rgammam =  bose_einstein( w' ) + fermi_dirac( -w + w' )
  !
  USE kinds, ONLY : DP
  ! 
  IMPLICIT NONE
  !
  REAL(DP) :: temp, omega, omegap, rgammap, rgammam
  REAL(DP) :: eps=1.0d-6
  !
  rgammap = 0.d0
  rgammam = 0.d0
  IF ( ABS(temp) < eps ) THEN
     rgammap = 0.d0
     rgammam = 1.d0
  ELSEIF ( omegap .gt. 0.d0 ) THEN 
     rgammap = 0.5d0 * (   tanh( 0.5d0 * ( omega + omegap ) / temp ) &
                         - 1.d0 / tanh( 0.5d0 * omegap / temp ) )
     rgammam = 0.5d0 * (   tanh( 0.5d0 * ( omega - omegap ) / temp ) &
                         + 1.d0 / tanh( 0.5d0 * omegap / temp ) )
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE gamma_acont
  !                         
  !-----------------------------------------------------------------------
