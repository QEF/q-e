  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE kernel_aniso_iaxis( itemp )
  !-----------------------------------------------------------------------
  !!  
  !! Compute kernels K_{+}(ik,iq,ibnd,jbnd;n,n',T) and K_{-}(ik,iq,ibnd,jbnd;n,n',T)
  !! and store them in memory
  !!
  USE kinds,         ONLY : DP
  USE epwcom,        ONLY : fsthick
  USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, estemp, AKeri, ekfs, ef0, ixkqf, ixqfs, nqfs
  USE constants_epw, ONLY : pi 
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: itemp
  !! Counter on temperature
  !
  ! Local variables
  INTEGER  :: ik, iq, iq0, ibnd, jbnd, iw, n, lower_bnd, upper_bnd 
  REAL(DP) :: omega, lambda_eph
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  IF ( .not. ALLOCATED(AKeri) ) ALLOCATE( AKeri(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,2*nsiw(itemp)) )
  AKeri(:,:,:,:,:) = 0.d0
  !
  ! RM - if lambdar_aniso_ver2 is used then one needs to CALL evaluate_a2fij
  !
  DO ik = lower_bnd, upper_bnd
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iq = 1, nqfs(ik)
              ! iq0 - index of q-point on the full q-mesh
              iq0 = ixqfs(ik,iq)
              DO jbnd = 1, nbndfs
                 IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                    DO iw = 1, 2*nsiw(itemp)
                       n = iw - 1
                       omega = dble(2*n) * pi * estemp(itemp)
                       CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, omega, lambda_eph )
                       !CALL lambdar_aniso_ver2( ik, iq, ibnd, jbnd, omega, lambda_eph )
                       AKeri(ik,iq,ibnd,jbnd,iw) = lambda_eph
                    ENDDO ! iw
                 ENDIF
              ENDDO ! jbnd
           ENDDO ! iq
        ENDIF
     ENDDO ! ibnd
  ENDDO ! ik
  !
  RETURN
  !
  END SUBROUTINE kernel_aniso_iaxis
  !
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_memlt_aniso_iaxis( itemp )
  !-----------------------------------------------------------------------
  !!  
  !! Estimate the memory requirements for anisotropic Eliashberg equations 
  !! on imaginary axis
  !!
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE epwcom,        ONLY : max_memlt
  USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, nqfs, limag_fly, memlt_pool
  USE mp_global, ONLY : inter_pool_comm, my_pool_id
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  INTEGER  :: itemp, lower_bnd, upper_bnd, imelt
  REAL(DP) :: rmelt
  !
  ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
  my_pool_id = 0
#endif  
  !
  limag_fly = .false.
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  !
  ! get the size of the AKeri kernels that need to stored in each pool
  ! imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * ( 2 * nsiw(itemp) )
  ! rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
  ! RM - avoid problems when imelt is greater than (2^31)-1 (singed long integer) 
  !
  imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:))
  rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
  imelt = nbndfs**2 * ( 2 * nsiw(itemp) )
  rmelt = dble(imelt) * rmelt 
  rmelt = rmelt + memlt_pool(my_pool_id+1)
  !
  memlt_pool(:) = 0.d0
  memlt_pool(my_pool_id+1) = rmelt
  !
  ! collect contributions from all pools
  CALL mp_sum( memlt_pool, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  IF ( maxval(memlt_pool(:)) .gt. max_memlt ) THEN
     WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of required memory per pool :", &
          " ~= ", maxval(memlt_pool(:)), " Gb"
     limag_fly = .true.
     !
     ! remove memory required for AKeri 
     ! imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * ( 2 * nsiw(itemp) )
     ! RM - avoid problems when imelt is greater than (2^31)-1 (singed long integer)
     !
     imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:))
     rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
     imelt = nbndfs**2 * ( 2 * nsiw(itemp) )
     rmelt = dble(imelt) * rmelt
     rmelt = - rmelt + memlt_pool(my_pool_id+1)
     !
     memlt_pool(:) = 0.d0
     memlt_pool(my_pool_id+1) = rmelt
     !
  ! collect contributions from all pools
  CALL mp_sum( memlt_pool, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)

  ENDIF
  !
  IF ( limag_fly ) THEN
     WRITE(stdout,'(/,5x,a/)') "AKeri is calculated on the fly since its size exceedes max_memlt"
  ELSE
     WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of allocated memory per pool :", &
          " ~= ", maxval(memlt_pool(:)), " Gb"
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE eliashberg_memlt_aniso_iaxis
  !  
  !-----------------------------------------------------------------------
  SUBROUTINE lambdar_aniso_ver1( ik, iq, ibnd, jbnd, omega, lambda_eph )
  !-----------------------------------------------------------------------
  !
  ! computes lambda(ik,iq,ibnd,jbnd;n-n')   
  ! reference H. Choi et. al, Physica C 385, 66 (2003)
  !
  ! input
  !
  ! ik - index k-point 
  ! iq - index q-point 
  ! ibnd - index band ibnd at k-point
  ! jbnd - index band jbnd at k+q-point
  ! omega  - frequency 
  !
  ! output
  !
  ! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;n-n')
  !
  !
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE elph2,         ONLY : wf
  USE epwcom,        ONLY : eps_acustic
  USE eliashbergcom, ONLY : ixqfs, g2, dosef
  ! 
  IMPLICIT NONE
  !
  INTEGER :: ik, iq, iq0, ibnd, jbnd, imode
  REAL(DP) :: omega, lambda_eph
  !
  ! iq0 - index of q-point on the full q-mesh
  iq0 = ixqfs(ik,iq)
  lambda_eph = 0.d0
  DO imode = 1, nmodes  ! loop over frequency modes
     IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
        lambda_eph = lambda_eph + g2(ik,iq,ibnd,jbnd,imode) * wf(imode,iq0) & 
                   / ( wf(imode,iq0)**2.d0 + omega**2.d0 )
     ENDIF
  ENDDO 
  lambda_eph = 2.d0 * lambda_eph * dosef
  !
  RETURN
  !
  END SUBROUTINE lambdar_aniso_ver1
  !
  !-----------------------------------------------------------------------
  SUBROUTINE lambdar_aniso_ver2( ik, iq, ibnd, jbnd, omega, lambda_eph )
  !-----------------------------------------------------------------------
  !
  ! computes lambda(ik,iq,ibnd,jbnd;n-n')   
  ! reference H. Choi et. al, Physica C 385, 66 (2003)
  !
  ! input
  !
  ! ik - index k-point 
  ! iq - index q-point 
  ! ibnd - index band ibnd at k-point
  ! jbnd - index band jbnd at k+q-point
  ! omega  - frequency 
  !
  ! output
  !
  ! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;n-n')
  !
  USE kinds,         ONLY : DP
  USE epwcom,        ONLY : nqstep
  USE eliashbergcom, ONLY : a2fij, wsph, dwsph
  !                 
  IMPLICIT NONE        
  !                    
  INTEGER :: ik, iq, ibnd, jbnd, iwph
  REAL(DP) :: omega, lambda_eph 
  !                          
  lambda_eph = 0.d0
  DO iwph = 1, nqstep
     lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik,iq,ibnd,jbnd,iwph) & 
                / ( wsph(iwph)**2.d0 + omega**2.d0 )
  ENDDO 
  lambda_eph = 2.d0 * lambda_eph * dwsph
  !
  RETURN
  !
  END SUBROUTINE lambdar_aniso_ver2
  !
  !-----------------------------------------------------------------------
  SUBROUTINE kernel_aniso_iaxis_analytic_cont( itemp )
  !-----------------------------------------------------------------------
  !!  
  !! computes kernels K_{+}(w,iw_n,T) and K_{-}(w,iw_n,T)
  !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
  !!
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : wqf
  USE epwcom,        ONLY : muc, fsthick
  USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, ADeltai, nkfs, nbndfs, dosef, ixkqf, ixqfs, nqfs, & 
                            w0g, ekfs, ef0, ADsumi, AZsumi
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
  !      
  IMPLICIT NONE
  !
  INTEGER :: iw, iwp, itemp, ik, iq, iq0, ibnd, jbnd, lower_bnd, upper_bnd, imelt
  REAL(DP) :: esqrt, kernelp, kernelm, weight
  REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:)
  COMPLEX(DP) :: lambda_eph
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  !
  ! get memory size required for wesqrt, desqrt, ADsumi, AZsumi
  imelt = 2 * nbndfs * nkfs * nsiw(itemp) + 2 * ( upper_bnd - lower_bnd + 1 ) * nbndfs * nsw 
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nbndfs,nkfs,nsiw(itemp)) )
  IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nbndfs,nkfs,nsiw(itemp)) )
  !
  DO ik = lower_bnd, upper_bnd
     IF ( .not. ALLOCATED(ADsumi) ) ALLOCATE( ADsumi(nbndfs,lower_bnd:upper_bnd,nsw) )
     IF ( .not. ALLOCATED(AZsumi) ) ALLOCATE( AZsumi(nbndfs,lower_bnd:upper_bnd,nsw) )
  ENDDO
  ADsumi(:,:,:) = 0.d0
  AZsumi(:,:,:) = 0.d0
  !
  ! RM - if lambdai_aniso_ver2 is used then one needs to CALL evaluate_a2fij
  !
  DO ik = lower_bnd, upper_bnd
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iq = 1, nqfs(ik)
              ! iq0 - index of q-point on the full q-mesh
              iq0 = ixqfs(ik,iq)
              DO jbnd = 1, nbndfs
                 IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                    weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                    DO iw = 1, nsw ! loop over omega 
                       DO iwp = 1, nsiw(itemp) ! loop over iw_n
                          CALL lambdai_aniso_ver1( ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph )
                          !CALL lambdai_aniso_ver2( ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph )
                          kernelp = 2.d0 * real(lambda_eph)
                          kernelm = 2.d0 * aimag(lambda_eph)
                          IF ( iw .eq. 1 ) THEN
                             esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + ADeltai(jbnd,ixkqf(ik,iq0),iwp)**2.d0 )
                             wesqrt(jbnd,ixkqf(ik,iq0),iwp) =  wsi(iwp) * esqrt
                             desqrt(jbnd,ixkqf(ik,iq0),iwp) =  ADeltai(jbnd,ixkqf(ik,iq0),iwp) * esqrt
                          ENDIF
                          AZsumi(ibnd,ik,iw) = AZsumi(ibnd,ik,iw) & 
                                             + weight * wesqrt(jbnd,ixkqf(ik,iq0),iwp) * kernelm 
                          ADsumi(ibnd,ik,iw) = ADsumi(ibnd,ik,iw) & 
                                             + weight * desqrt(jbnd,ixkqf(ik,iq0),iwp) * ( kernelp - 2.d0 * muc ) 
                       ENDDO ! iwp
                    ENDDO ! iw
                 ENDIF
              ENDDO ! jbnd
           ENDDO ! iq
        ENDIF
     ENDDO ! ibnd
  ENDDO ! ik
  !
  IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
  IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
  !  
  ! remove memory allocated for wesqrt, desqrt
  imelt = 2 * nbndfs * nkfs * nsiw(itemp) 
  CALL mem_size_eliashberg ( -imelt )
  !
  RETURN
  !
  END SUBROUTINE kernel_aniso_iaxis_analytic_cont
  !                                                
  !-----------------------------------------------------------------------
  SUBROUTINE lambdai_aniso_ver1( ik, iq, ibnd, jbnd, omega, omegap, lambda_eph )
  !-----------------------------------------------------------------------
  !!
  !! computes lambda(w-iw_n)   
  !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
  !!
  !! input
  !!  
  !! ik - index k-point
  !! iq - index q-point 
  !! ibnd - index band ibnd at k-point
  !! jbnd - index band jbnd at k+q-point
  !! iw     - index frequency w on the real-axis
  !! iwp    - index frequency iw_n on the imaginary-axis
  !! omega  - frequency w at point iw
  !! omegap - frequency w_n at point iwp
  !!     
  !! output 
  !!        
  !! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
  !!        
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE elph2,         ONLY : wf
  USE epwcom,        ONLY : eps_acustic
  USE eliashbergcom, ONLY : ixqfs, g2, dosef
  USE constants_epw, ONLY : ci
  !     
  IMPLICIT NONE
  !  
  INTEGER :: ik, iq, iq0, ibnd, jbnd, imode
  REAL(DP) :: omega, omegap
  COMPLEX(DP) :: lambda_eph
  !
  ! iq0 - index of q-point on the full q-mesh
  iq0 = ixqfs(ik,iq)
  lambda_eph = (0.d0,0.d0)
  DO imode = 1, nmodes  ! loop over frequency modes
     IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
        lambda_eph = lambda_eph +  g2(ik,iq,ibnd,jbnd,imode) * wf(imode,iq0) & 
                   / ( wf(imode,iq0)**2.d0 - (omega - ci*omegap)**2.d0 )
     ENDIF
  ENDDO ! iwph
  lambda_eph = 2.d0 * lambda_eph * dosef
  !
  RETURN
  !
  END SUBROUTINE lambdai_aniso_ver1
  !
  !-----------------------------------------------------------------------               
  SUBROUTINE lambdai_aniso_ver2( ik, iq, ibnd, jbnd, omega, omegap, lambda_eph )
  !-----------------------------------------------------------------------
  !!
  !! computes lambda(w-iw_n)   
  !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
  !!
  !! input
  !!  
  !! ik - index k-point
  !! iq - index q-point 
  !! ibnd - index band ibnd at k-point
  !! jbnd - index band jbnd at k+q-point
  !! iw     - index frequency w on the real-axis
  !! iwp    - index frequency iw_n on the imaginary-axis
  !! omega  - frequency w at point iw
  !! omegap - frequency w_n at point iwp
  !!     
  !! output 
  !!        
  !! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
  !!        
  USE kinds,         ONLY : DP
  USE epwcom,        ONLY : nqstep
  USE eliashbergcom, ONLY : a2fij, dwsph, wsph
  USE constants_epw, ONLY : ci
  !     
  IMPLICIT NONE
  !  
  INTEGER :: ik, iq, ibnd, jbnd, iwph
  REAL(DP) :: omega, omegap
  COMPLEX(DP) :: lambda_eph
  !
  lambda_eph = (0.d0,0.d0)
  DO iwph = 1, nqstep
     lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik,iq,ibnd,jbnd,iwph) &
                / ( wsph(iwph)**2.d0 - (omega - ci*omegap)**2.d0 )
  ENDDO ! iwph
  lambda_eph = 2.d0 * lambda_eph * dwsph
  !
  RETURN
  !
  END SUBROUTINE lambdai_aniso_ver2
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_a2fij
  !-----------------------------------------------------------------------
  !!
  !! computes the anisotropic spectral function a2F(k,k',w) 
  !!
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE elph2,         ONLY : wf
  USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq
  USE eliashbergcom, ONLY : nkfs, nbndfs, g2, a2fij, ixkqf, ixqfs, nqfs, ekfs, ef0, & 
                            dosef, wsph
  ! 
  IMPLICIT NONE
  !
  INTEGER :: ik, iq, iq0, iwph, ibnd, jbnd, imode, lower_bnd, upper_bnd
  REAL(DP) :: weight
  REAL(DP), EXTERNAL :: w0gauss
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  IF ( .not. ALLOCATED(a2fij) ) ALLOCATE(a2fij(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,nqstep))
  a2fij(:,:,:,:,:) = 0.d0
  !
  DO ik = lower_bnd, upper_bnd 
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iq = 1, nqfs(ik)
              ! iq0 - index of q-point on the full q-mesh
              iq0 = ixqfs(ik,iq)
              DO jbnd = 1, nbndfs
                 IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                    DO imode = 1, nmodes
                       IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
                          DO iwph = 1, nqstep
                             weight  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / degaussq, 0 ) / degaussq
                             a2fij(ik,iq,ibnd,jbnd,iwph) = a2fij(ik,iq,ibnd,jbnd,iwph) &
                                                         + weight * dosef * g2(ik,iq,ibnd,jbnd,imode)
                          ENDDO ! iwph
                       ENDIF ! wf
                    ENDDO ! imode
                 ENDIF
              ENDDO ! jbnd
           ENDDO ! iq
        ENDIF
     ENDDO ! ibnd
  ENDDO ! ik
  !
  RETURN
  !
  END SUBROUTINE evaluate_a2fij
  !
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_memlt_aniso_acon
  !-----------------------------------------------------------------------
  !!  
  !! Estimate the memory requirements for the anisotropic Eliashberg funtion
  !! used for analytic continuation from imaginary to real axis
  !!
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE epwcom,        ONLY : nqstep, max_memlt
  USE eliashbergcom, ONLY : nkfs, nbndfs, nqfs, lacon_fly, memlt_pool
  USE mp_global, ONLY : inter_pool_comm, my_pool_id
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  INTEGER  :: lower_bnd, upper_bnd, imelt
  REAL(DP) :: rmelt
  !
  ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
  my_pool_id = 0
#endif  
  !
  lacon_fly = .false.
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  !
  ! get the size of a2fij that need to stored in each pool
  imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nqstep
  rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
  rmelt = rmelt + memlt_pool(my_pool_id+1)
  !
  memlt_pool(:) = 0.d0
  memlt_pool(my_pool_id+1) = rmelt
  !
  ! collect contributions from all pools
  CALL mp_sum( memlt_pool, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  IF ( maxval(memlt_pool(:)) .gt. max_memlt ) THEN
     WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of required memory per pool :", &
          " ~= ", maxval(memlt_pool(:)), " Gb"
     lacon_fly = .true.
     !
     ! remove memory required for a2fij
     imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nqstep
     CALL mem_size_eliashberg( -imelt )
  ENDIF
  !
  IF ( lacon_fly ) THEN
     WRITE(stdout,'(/,5x,a/)') "a2fij is calculated on the fly since its size exceedes max_memlt"
  ELSE
     WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of allocated memory per pool :", &
          " ~= ", maxval(memlt_pool(:)), " Gb"
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE eliashberg_memlt_aniso_acon
  !   
  !-----------------------------------------------------------------------
