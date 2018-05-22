  !                                                                            
  ! Copyright (C) 2010-2017 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE print_gkk ( iq )
  !-----------------------------------------------------------------------
  !! 
  !! Print the |g| vertex for all n,n' and modes in meV and do average
  !! on degenerate states. 
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            nkf, epf17, xkf, nkqtotf, wf
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iq
  !! Current q-point index 
  !
  ! Local variables 
  ! 
  INTEGER :: lower_bnd
  !! Lower bounds index after k or q paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k or q paral
  INTEGER :: ik
  !! K-point index
  INTEGER :: ikk
  !! K-point index
  INTEGER :: ikq
  !! K+q-point index
  INTEGER :: ibnd
  !! Band index
  INTEGER :: jbnd
  !! Band index
  INTEGER :: pbnd
  !! Band index
  INTEGER :: nu
  !! Mode index
  INTEGER :: mu
  !! Mode index
  INTEGER :: n
  !! Number of modes
  !
  REAL(kind=DP) :: xkf_all(3, nkqtotf)
  !! Collect k-point coordinate from all pools in parallel case
  REAL(kind=DP) :: etf_all(nbndsub,nkqtotf)
  !! Collect eigenenergies from all pools in parallel case
  REAL(kind=DP) :: wq
  !! Phonon frequency 
  REAL(kind=DP) :: w_1
  !! Temporary phonon freq. 1
  REAL(kind=DP) :: w_2
  !! Temporary phonon freq. 2
  REAL(kind=DP) :: gamma
  !! Temporary electron-phonon matrix element
  REAL(kind=DP) :: ekk
  !! Eigenenergies at k
  REAL(kind=DP) :: ekq
  !! Eigenenergies at k+q
  REAL(kind=DP) :: g2
  !! Temporary electron-phonon matrix element square
  REAL(kind=DP), ALLOCATABLE :: epc(:,:,:,:)
  !! g vectex accross all pools 
  REAL(kind=DP), ALLOCATABLE :: epc_sym(:,:,:)
  !! Temporary g-vertex for each pool 
  REAL(kind=DP), PARAMETER :: eps = 0.01/ryd2mev
  !! Tolerence to be degenerate
  ! 
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds( nkqtotf/2, lower_bnd, upper_bnd )
  ! 
  ALLOCATE ( epc (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkqtotf/2) )
  ALLOCATE ( epc_sym (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes) )
  ! 
  epc(:,:,:,:)   = zero
  epc_sym(:,:,:) = zero
  !
  ! First do the average over bands and modes for each pool
  DO ik = 1, nkf
    ikk = 2 * ik - 1
    ikq = ikk + 1
    ! 
    DO nu = 1, nmodes
      wq = wf (nu, iq)
      !DO ibnd = ibndmin, ibndmax
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
          gamma = ( ABS(epf17 (jbnd, ibnd, nu, ik)) )**two
          IF (wq > 0.d0) then
              gamma = gamma / (two * wq)
          ELSE
              gamma = 0.d0
          ENDIF
          gamma = sqrt(gamma)
          ! gamma = |g| [Ry]
          epc(ibnd,jbnd,nu,ik+lower_bnd-1) = gamma
        ENDDO ! jbnd
      ENDDO   ! ibnd        
    ENDDO ! loop on modes
    !
    !  Here we "SYMMETRIZE": actually we simply take the averages over
    !  degenerate states, it is only a convention because g is gauge-dependent!
    !
    ! first the phonons
    DO ibnd = 1, ibndmax-ibndmin+1
      DO jbnd = 1, ibndmax-ibndmin+1
        DO nu = 1, nmodes
          w_1 = wf(nu,iq)
          g2 = 0.d0
          n  = 0
          DO mu = 1, nmodes
            w_2 = wf(mu,iq)
            IF ( abs(w_2-w_1).lt.eps ) THEN
              n = n + 1
              g2 = g2 + epc(ibnd,jbnd,mu,ik+lower_bnd-1)*epc(ibnd,jbnd,mu,ik+lower_bnd-1)
            ENDIF
          ENDDO
          g2 = g2 / float(n)
          epc_sym (ibnd, jbnd, nu) = sqrt (g2)
        ENDDO
      ENDDO
    ENDDO
    epc(:,:,:,ik+lower_bnd-1) = epc_sym
    ! Then the k electrons
    DO nu = 1, nmodes
      DO jbnd = 1, ibndmax-ibndmin+1
        DO ibnd = 1, ibndmax-ibndmin+1
          w_1 = etf (ibndmin-1+ibnd, ikk)
          g2 = 0.d0
          n  = 0
          DO pbnd = 1, ibndmax-ibndmin+1
            w_2 = etf (ibndmin-1+pbnd, ikk)
            IF ( abs(w_2-w_1).lt.eps ) THEN
              n = n + 1
              g2 = g2 + epc(pbnd,jbnd,nu,ik+lower_bnd-1)*epc(pbnd,jbnd,nu,ik+lower_bnd-1)
            ENDIF
          ENDDO
          g2 = g2 / float(n)
          epc_sym (ibnd, jbnd, nu) = sqrt (g2)
        ENDDO
      ENDDO
    ENDDO
    epc(:,:,:,ik+lower_bnd-1) = epc_sym
    !
    ! and finally the k+q electrons
    DO nu = 1, nmodes
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
          w_1 = etf (ibndmin-1+jbnd, ikq)
          g2 = 0.d0
          n  = 0
          DO pbnd = 1, ibndmax-ibndmin+1
            w_2 = etf(ibndmin-1+pbnd, ikq)
            IF ( abs(w_2-w_1).lt.eps ) then
              n = n + 1
              g2 = g2 + epc(ibnd,pbnd,nu,ik+lower_bnd-1)*epc(ibnd,pbnd,nu,ik+lower_bnd-1)
            ENDIF
          ENDDO
          g2 = g2 / float(n)
          epc_sym (ibnd, jbnd, nu) = sqrt (g2)
        ENDDO
      ENDDO
    ENDDO
    epc(:,:,:,ik+lower_bnd-1) = epc_sym
    ! 
  ENDDO ! k-points
  ! 
  ! We need quantity from all the pools
  xkf_all(:,:) = zero
  etf_all(:,:) = zero
  !
#if defined(__MPI)
  !
  ! Note that poolgather2 works with the doubled grid (k and k+q)
  !
  CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
  CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
  CALL mp_sum( epc, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  !
#endif
  !
  ! Only master writes
  IF (mpime.eq.ionode_id) THEN
    !
    WRITE(stdout, '(5x,a)') ' Electron-phonon vertex |g| (meV)'
    !
    WRITE(stdout, '(/5x,"iq = ",i7," coord.: ", 3f12.7)') iq, xqf (:, iq)
    DO ik = 1, nkqtotf/2
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      WRITE(stdout,'(5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:,ikk)
      WRITE(stdout, '(5x,a)') ' ibnd     jbnd     imode   enk[eV]    enk+q[eV]  omega(q)[meV]   |g|[meV]'
      WRITE(stdout,'(5x,a)') repeat('-',78)
      !
      DO ibnd = 1, ibndmax-ibndmin+1
        ekk = etf_all (ibndmin-1+ibnd, ikk) 
        DO jbnd = 1, ibndmax-ibndmin+1
          ekq = etf_all (ibndmin-1+jbnd, ikq) 
          DO nu = 1, nmodes
            WRITE(stdout,'(3i9,3f12.4,1e20.10)') ibndmin-1+ibnd, ibndmin-1+jbnd, nu, ryd2ev * ekk, ryd2ev * ekq, &
                                         ryd2mev * wf(nu,iq), ryd2mev * epc(ibnd,jbnd,nu,ik)
          ENDDO
        ENDDO  
        !
        !
      ENDDO
      WRITE(stdout,'(5x,a/)') repeat('-',78)
      !
    ENDDO
  ENDIF ! master node
  ! 
  DEALLOCATE ( epc )
  DEALLOCATE ( epc_sym )
  !
  END SUBROUTINE print_gkk
