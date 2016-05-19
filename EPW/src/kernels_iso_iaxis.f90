  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE kernel_iso_iaxis( itemp )
  !-----------------------------------------------------------------------
  !  
  ! computes kernels K_{+}(n,n',T) and K_{-}(n,n'T)
  ! reference W. E. Pickett, PRB 26, 1186 (1982)
  !
  !
  USE kinds, ONLY : DP
  USE constants_epw, ONLY : pi
  USE eliashbergcom, ONLY : nsiw, estemp, Keri
  ! 
  IMPLICIT NONE
  !
  INTEGER  :: iw, itemp, n
  REAL(DP) :: omega, lambda_eph
  !
  IF ( .not. ALLOCATED(Keri) ) ALLOCATE( Keri(2*nsiw(itemp)) )
  Keri(:) = 0.d0
  !
  DO iw = 1, 2*nsiw(itemp)
     n = iw - 1
     omega = dble(2*n) * pi * estemp(itemp)
     CALL lambdar_iso( omega, lambda_eph )
     Keri(iw) = lambda_eph
  ENDDO 
  !
  RETURN
  !
  END SUBROUTINE kernel_iso_iaxis                                                       
  !
  !-----------------------------------------------------------------------
  SUBROUTINE lambdar_iso( omega, lambda_eph )
  !-----------------------------------------------------------------------
  !
  ! computes lambda(n-n')   
  ! reference W. E. Pickett, PRB 26, 1186 (1982)
  !
  ! input
  !
  ! omega  - frequency 
  !
  ! output
  !
  ! lampda_eph - electron-phonon coupling lambda(n-n')
  !
  !
  USE kinds, ONLY : DP
  USE epwcom, ONLY : nqstep
  USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iwph
  REAL(DP) :: omega, lambda_eph
  !
  lambda_eph = 0.d0
  DO iwph = 1, nqstep  ! loop over Omega (integration variable)
     lambda_eph = lambda_eph + wsph(iwph) * a2f_iso(iwph) & 
                / ( wsph(iwph)**2.d0 + omega**2.d0 )
  ENDDO ! iwph
  lambda_eph = 2.d0 * lambda_eph * dwsph 
  !
  RETURN
  !
  END SUBROUTINE lambdar_iso

  !-----------------------------------------------------------------------
  SUBROUTINE kernel_iso_iaxis_analytic_cont( itemp )
  !-----------------------------------------------------------------------
  !  
  ! computes kernels K_{+}(w,iw_n,T) and K_{-}(w,iw_n,T)
  ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
  !
  !
  USE kinds,         ONLY : DP
  USE epwcom,        ONLY : muc
  USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, Deltai, Dsumi, Zsumi
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, iwp, itemp
  REAL(DP) :: esqrt, kernelp, kernelm
  REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
  COMPLEX(DP) :: lambda_eph
  !
  IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsiw(itemp)) )
  IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsiw(itemp)) )
  IF ( .not. ALLOCATED(Dsumi) )  ALLOCATE( Dsumi(nsw) )
  IF ( .not. ALLOCATED(Zsumi) )  ALLOCATE( Zsumi(nsw) )
  Dsumi(:) = 0.d0
  Zsumi(:) = 0.d0
  !
  DO iw = 1, nsw ! loop over omega
     DO iwp = 1, nsiw(itemp) ! loop over iw_n
        CALL lambdai_iso( ws(iw), wsi(iwp), lambda_eph )
        kernelp = 2.d0 * real(lambda_eph)
        kernelm = 2.d0 * aimag(lambda_eph) 
        IF ( iw .eq. 1 ) THEN
           esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + Deltai(iwp)**2.d0 )
           wesqrt(iwp) =  wsi(iwp) * esqrt
           desqrt(iwp) =  Deltai(iwp) * esqrt
        ENDIF
        Zsumi(iw) = Zsumi(iw) + kernelm * wesqrt(iwp)
        Dsumi(iw) = Dsumi(iw) + ( kernelp - 2.d0 * muc ) * desqrt(iwp)
     ENDDO
  ENDDO
  !
  IF( ALLOCATED(wesqrt) ) DEALLOCATE (wesqrt)
  IF( ALLOCATED(desqrt) ) DEALLOCATE (desqrt)
  !   
  RETURN
  !
  END SUBROUTINE kernel_iso_iaxis_analytic_cont      
  !                                                
  !-----------------------------------------------------------------------
  SUBROUTINE lambdai_iso( omega, omegap, lambda_eph )
  !-----------------------------------------------------------------------
  !
  ! computes lambda(w-iw_n)   
  ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
  !
  ! input
  !
  ! iw     - index frequency w on the real-axis
  ! iwp    - index frequency iw_n on the imaginary-axis
  ! omega  - frequency w at point iw
  ! omegap - frequency w_n at point iwp
  !
  ! output
  !
  ! lampda_eph - electron-phonon coupling lambda(w-iw_n)
  !
  !
  USE kinds, ONLY : DP
  USE epwcom,        ONLY : nqstep
  USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
  USE constants_epw, ONLY : ci
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iwph
  REAL(DP) :: omega, omegap
  COMPLEX(DP) :: lambda_eph
  !
  lambda_eph = (0.d0,0.d0)
  DO iwph = 1, nqstep  ! loop over Omega (integration variable)
     lambda_eph = lambda_eph & 
                + wsph(iwph) * a2f_iso(iwph) / ( wsph(iwph)**2.d0 - (omega - ci*omegap)**2.d0 )
  ENDDO ! iwph
  lambda_eph = lambda_eph * 2.d0 * dwsph 
  !
  RETURN
  !
  END SUBROUTINE lambdai_iso
  !                                        
  !-----------------------------------------------------------------------               
