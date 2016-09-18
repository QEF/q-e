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
  !!
  !!  This subroutine writes output files for phonon dispersion and band structure 
  !!  RM : this subroutine should be tested
  !!  SP : Modified so that it works with the current plotband.x of QE 5
  !!
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE phcom,      ONLY : nmodes
  USE epwcom,     ONLY : nbndsub, filqf, filkf
  USE elph2,      ONLY : etf, nkf, nqtotf, wf, xkf, xqf, nkqtotf
  USE constants_epw, ONLY : ryd2mev, ryd2ev
  USE io_epw,     ONLY : iufilfreq, iufileig
  USE elph2,      ONLY : nkqf
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_barrier, mp_sum
  USE mp_global,  ONLY : inter_pool_comm, my_pool_id
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !! Global k-point index
  INTEGER :: ikk
  !! Index for the k-point
  INTEGER :: ikq
  !! Index for the q-point
  INTEGER :: ibnd
  !! Band index
  INTEGER :: imode
  !! Mode index
  INTEGER :: iq
  !! Global q-point index
  INTEGER :: nksqtotf
  !! Sum of total number of k+q points
  REAL(kind=DP) :: dist
  !! Distance from Gamma
  REAL(kind=DP) :: dprev
  !! Previous distance
  REAL(kind=DP) :: dcurr
  !! Current distance
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:)
  !! K-points on the full k grid (all pools)
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  !! Eigenenergies on the full k grid (all pools)
  !
  nksqtotf =  nkqtotf/2
  !
  IF (filqf /= ' ') THEN
    ! 
    IF ( my_pool_id .eq. ionode_id ) THEN
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
    ENDIF
    CALL mp_barrier(inter_pool_comm)
  ENDIF ! filqf
  ! 
  IF (filkf /= ' ') THEN
    !
    DO ik = 1, nkf
       !
       ikk = 2 * ik - 1
       ikq = ikk + 1
       !
    ENDDO
    !
    IF ( .not. ALLOCATED(xkf_all) ) ALLOCATE ( xkf_all( 3, nkqtotf)) 
    IF ( .not. ALLOCATED(etf_all) ) ALLOCATE (etf_all( nbndsub, nkqtotf))
    !
#if defined(__MPI)
    CALL poolgather2( 3,       nkqtotf, nkqf, xkf, xkf_all )
    CALL poolgather2( nbndsub, nkqtotf, nkqf, etf, etf_all )
    CALL mp_barrier(inter_pool_comm)
#else    
    !
    xkf_all = xkf
    etf_all = etf
#endif
    !
    IF ( my_pool_id .eq. ionode_id ) THEN
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
         ikk = 2 * ik - 1
         ikq = ikk + 1
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
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(xkf_all)) DEALLOCATE( xkf_all )
    IF ( ALLOCATED(etf_all)) DEALLOCATE( etf_all )
    !
  ENDIF ! filkf
  !
  RETURN
  !
  END SUBROUTINE plot_band
  !
  !-----------------------------------------------------------------------
