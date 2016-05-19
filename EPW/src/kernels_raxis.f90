  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE kernel_raxis( iw, iwp, itemp, kernelp, kernelm )
  !-----------------------------------------------------------------------
  !
  ! computes kernels K_{+}(w,w',T) and K_{-}(w,w'T)  
  ! reference M. J. Holcomb, PRB 54, 6648 (1996)   
  !
  ! input
  !
  ! iw     - index frequency w : ws(iw)
  ! iwp    - index frequency w' : ws(iwp)
  ! itemp  - index temperature
  !
  ! output
  !
  ! kernelp - phonon kernel K_{+}(w,w',T)
  ! kernelm - phonon kernel K_{-}(w,w',T)
  !
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : pi, ci
  USE epwcom,        ONLY : nqstep
  USE eliashbergcom, ONLY : a2f_iso, bewph, wsph, dwsph, ws, fdwp, estemp
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, iwp, iwph, itemp, ngaussw0
  REAL(DP) :: degaussw0, f1, f2, f3, f4
  COMPLEX(DP) :: e1, e2, e3, e4, kernelp, kernelm
  REAL(DP), EXTERNAL :: wgauss, w0gauss
  REAL(DP) :: eps=1.0d-6
  !
  degaussw0 = 1.d0 * dwsph
  ngaussw0 = 0
  !
  f1 = 0.d0
  f2 = 0.d0
  f3 = 0.d0
  f4 = 0.d0
  kernelp = (0.d0, 0.d0)
  kernelm = (0.d0, 0.d0)
  e1 = (0.d0, 0.d0)
  e2 = (0.d0, 0.d0)
  e3 = (0.d0, 0.d0)
  e4 = (0.d0, 0.d0)
  !
  IF ( .not. ALLOCATED(bewph) ) ALLOCATE( bewph(nqstep) )
  ! Bose-Einstein distribution
  DO iwph = 1, nqstep  ! loop over Omega (integration variable)
     IF ( iw .eq. 1 .AND. iwp .eq. 1 ) THEN
        IF ( ABS(estemp(itemp)) <  eps ) THEN
           bewph(iwph)  = 0.d0
        ELSE
           bewph(iwph) = wgauss( -wsph(iwph) / estemp(itemp), -99 )
           bewph(iwph) = bewph(iwph) / (1.d0 - 2.d0 * bewph(iwph))
        ENDIF
     ENDIF
     !
     ! a small complex number is added to denominator to move the pole away from the real-axis
     !
     ! in order to reduce the numerical noise at very small frequencies coming from the complex number 
     ! added in the denominator, the contribution of the imaginary part is reestimated using 
     ! delta function (RM notes) 
     !
     ! subtract the imaginary part coming from e1 to e4 and add instead the imaginary part 
     ! coming from f1 to f4
     !
     e1 = 1.d0 / ( wsph(iwph) + ws(iwp) + ws(iw) + ci*degaussw0 ) 
     e2 = 1.d0 / ( wsph(iwph) + ws(iwp) - ws(iw) - ci*degaussw0 ) 
     e3 = 1.d0 / ( wsph(iwph) - ws(iwp) + ws(iw) + ci*degaussw0 ) 
     e4 = 1.d0 / ( wsph(iwph) - ws(iwp) - ws(iw) - ci*degaussw0 ) 
     !
     ! estimate of the imaginary part using delta function
     f1 = w0gauss( ( wsph(iwph) + ws(iwp) + ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
     f2 = w0gauss( ( wsph(iwph) + ws(iwp) - ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
     f3 = w0gauss( ( wsph(iwph) - ws(iwp) + ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
     f4 = w0gauss( ( wsph(iwph) - ws(iwp) - ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
     !
     kernelp = kernelp + a2f_iso(iwph) &
             * (  ( 1.d0 - fdwp(iwp) + bewph(iwph) ) * ( e1 - ci*aimag(e1) - ci*pi*f1 + e2 - ci*aimag(e2) + ci*pi*f2 ) & 
                - (        fdwp(iwp) + bewph(iwph) ) * ( e3 - ci*aimag(e3) - ci*pi*f3 + e4 - ci*aimag(e4) + ci*pi*f4 ) )
     kernelm = kernelm + a2f_iso(iwph) &
             * (  ( 1.d0 - fdwp(iwp) + bewph(iwph) ) * ( e1 - ci*aimag(e1) - ci*pi*f1 - e2 + ci*aimag(e2) - ci*pi*f2 ) &
                + (        fdwp(iwp) + bewph(iwph) ) * ( e3 - ci*aimag(e3) - ci*pi*f3 - e4 + ci*aimag(e4) - ci*pi*f4 ) )
  ENDDO ! iwph
  kernelp = kernelp * dwsph 
  kernelm = kernelm * dwsph 
  !
  RETURN
  !
  END SUBROUTINE kernel_raxis
  !
  !-----------------------------------------------------------------------
