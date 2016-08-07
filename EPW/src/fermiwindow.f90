  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE fermiwindow 
  !-----------------------------------------------------------------------
  !
  !  find the band indices of the first
  !  and last state falling within the window e_fermi+-efermithickness
  ! 
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf
  USE epwcom,        ONLY : fsthick, nbndsub
  USE pwcom,         ONLY : ef
  USE mp,            ONLY : mp_max, mp_min
  USE mp_global,     ONLY : inter_pool_comm
  implicit none
  integer :: ik, ibnd
  real(kind=DP) :: ebnd, ebndmin, ebndmax
  real(kind=DP) :: tmp
  !
  !
  ibndmin = 100000
  ibndmax = 0
  ebndmin =  1.d8
  ebndmax = -1.d8
  !
  DO ik = 1, nkqf
    DO ibnd = 1, nbndsub
      ebnd = etf (ibnd, ik)
      !
      IF  ( abs(ebnd - ef) .lt. fsthick ) THEN
        ibndmin = min(ibnd,ibndmin)
        ibndmax = max(ibnd,ibndmax)
        ebndmin = min(ebnd,ebndmin)
        ebndmax = max(ebnd,ebndmax)
      ENDIF
      !
    ENDDO
  ENDDO
  !
  tmp = dble (ibndmin)
  CALL mp_min(tmp,inter_pool_comm)
  ibndmin = nint (tmp)
  CALL mp_min(ebndmin,inter_pool_comm)
  !
  tmp = dble (ibndmax)
  CALL mp_max(tmp, inter_pool_comm)
  ibndmax = nint (tmp)
  CALL mp_max(ebndmax,inter_pool_comm)
  !
  WRITE(stdout,'(/14x,a,i5,2x,a,f9.3)') 'ibndmin = ', ibndmin, 'ebndmin = ', ebndmin
  WRITE(stdout,'(14x,a,i5,2x,a,f9.3/)') 'ibndmax = ', ibndmax, 'ebndmax = ', ebndmax
  !
  END SUBROUTINE fermiwindow

