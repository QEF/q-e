  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine, Samuel Ponce
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE plot_band
  !-----------------------------------------------------------------------
  !
  !  This subroutine writes output files for phonon dispersion and band structure 
  !  RM : this subroutine should be tested
  !  SP : Modified so that it works with the current plotband.x of QE 5
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE phcom,      ONLY : nmodes
  USE control_lr, ONLY : lgamma
  USE epwcom,     ONLY : nbndsub, etf_mem 
  USE elph2,      ONLY : etf, nkf, nqtotf, wf, xkf, xqf, nkqtotf
  USE constants_epw, ONLY : ryd2mev, ryd2ev
  USE io_epw,     ONLY : iufilfreq, iufileig
#ifdef __PARA
  USE elph2,      ONLY : nkqf
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_barrier, mp_sum
  USE mp_global,  ONLY : me_pool, inter_pool_comm, my_pool_id, npool
  USE mp_world,   ONLY : mpime
#endif
  !
  IMPLICIT NONE
  !
  real(kind=DP) :: dist, dprev, dcurr
  INTEGER :: ik, ikk, ikq, ibnd, imode, iq
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:), etf_all(:,:)
  !
  INTEGER :: nksqtotf
  !
  IF ( .not. etf_mem ) CALL errore ('plot_band', 'etf_mem should be true', 1)
  !
  nksqtotf =  nkqtotf/2
  !
#ifdef __PARA
  IF ( my_pool_id .eq. ionode_id ) THEN
#endif
  !
  OPEN(iufilfreq, file = "phband.freq", form = 'formatted')
  WRITE(iufilfreq, '(" &plot nbnd=",i4,", nks=",i6," /")') nmodes, nqtotf
  !
  ! crystal to cartesian coordinates
  CALL cryst_to_cart( nqtotf, xqf, bg, 1 )
  !
  dist = 0.d0
  dprev = 0.d0
  dcurr = 0.d0
  DO iq = 1, nqtotf
     !
     IF ( iq .ne. 1 ) THEN  
        dist = sqrt(   ( xqf(1,iq) - xqf(1,iq-1) ) * ( xqf(1,iq) - xqf(1,iq-1) ) & 
                     + ( xqf(2,iq) - xqf(2,iq-1) ) * ( xqf(2,iq) - xqf(2,iq-1) ) & 
                     + ( xqf(3,iq) - xqf(3,iq-1) ) * ( xqf(3,iq) - xqf(3,iq-1) ))
     ELSE 
        dist = 0.d0
     ENDIF
     dcurr = dprev + dist
     dprev = dcurr
     WRITE(iufilfreq,'(10x,3f10.6)') xqf(:,iq)
     WRITE(iufilfreq,'(1000f10.4)') (wf(imode,iq)*ryd2mev, imode=1,nmodes)
     !WRITE(iufilfreq,'(1000f10.6)') dcurr, (wf(imode,iq)*ryd2mev, imode=1,nmodes)
     !
  ENDDO
  CLOSE(iufilfreq)
  !
  ! back from cartesian to crystal coordinates
  CALL cryst_to_cart( nqtotf, xqf, at, -1 )
  !
#ifdef __PARA
  ENDIF
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  DO ik = 1, nkf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
  ENDDO
  !
  ALLOCATE ( xkf_all( 3, nkqtotf) , etf_all( nbndsub, nkqtotf) )
  !
#ifdef __PARA
  !
  CALL poolgather2( 3,       nkqtotf, nkqf, xkf, xkf_all )
  CALL poolgather2( nbndsub, nkqtotf, nkqf, etf, etf_all )
  CALL mp_barrier(inter_pool_comm)
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  !
#endif
  !
#ifdef __PARA
  IF ( my_pool_id .eq. ionode_id ) THEN
#endif
  !
  OPEN(iufileig, file = "band.eig", form = 'formatted')
  WRITE(iufileig, '(" &plot nbnd=",i4,", nks=",i6," /")') nbndsub, nksqtotf
  !
  ! crystal to cartesian coordinates
  CALL cryst_to_cart( nkqtotf, xkf_all, bg, 1 )
  !
  dist = 0.d0
  dprev = 0.d0
  dcurr = 0.d0
  DO ik = 1, nksqtotf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     IF ( ikk .ne. 1 ) THEN
        dist = sqrt(   ( xkf_all(1,ikk) - xkf_all(1,ikk-2) ) * ( xkf_all(1,ikk) - xkf_all(1,ikk-2) ) &
                     + ( xkf_all(2,ikk) - xkf_all(2,ikk-2) ) * ( xkf_all(2,ikk) - xkf_all(2,ikk-2) ) &
                     + ( xkf_all(3,ikk) - xkf_all(3,ikk-2) ) * ( xkf_all(3,ikk) - xkf_all(3,ikk-2) ) )
     ELSE
        dist = 0.d0
     ENDIF
     dcurr = dprev + dist
     dprev = dcurr
     WRITE(iufileig,'(10x,3f10.6)') xkf_all(:,ikk)
     WRITE(iufileig,'(1000f10.4)') (etf_all(ibnd,ikk)*ryd2ev, ibnd=1,nbndsub)
     !WRITE(iufileig,'(1000f10.6)') dcurr, (etf_all(ibnd,ikk)*ryd2ev, ibnd=1,nbndsub)
     !
  ENDDO
  CLOSE(iufileig)
  !
  ! back from cartesian to crystal coordinates
  CALL cryst_to_cart( nkqtotf, xkf_all, at, -1 )
  !
#ifdef __PARA
  ENDIF
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  DEALLOCATE( xkf_all )
  DEALLOCATE( etf_all )
  !
  RETURN
  !
  END SUBROUTINE plot_band
  !
  !-----------------------------------------------------------------------
