  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_func_pl_q ( iqq, iq, totq )
  !-----------------------------------------------------------------------
  !
  !  Compute the electron spectral function including the  electron-
  !  phonon interaction in the Migdal approximation. 
  !  
  !  We take the trace of the spectral function to simulate the photoemission
  !  intensity. I do not consider the c-axis average for the time being.
  !  The main approximation is constant dipole matrix element and diagonal
  !  selfenergy. The diagonality can be checked numerically. 
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iospectral_sup, iospectral
  USE epwcom,        ONLY : nbndsub, &
                            fsthick, eptemp, ngaussw, degaussw, wmin_specfun,&
                            wmax_specfun, nw_specfun, &
                            efermi_read, fermi_energy,&
                            nel, meff, epsiHEG 
  USE pwcom,         ONLY : nelec, ef, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, &
                            wkf, nkf, wqf, xkf, nkqtotf,&
                            esigmar_all, esigmai_all, a_all,&
                            xqf, dmef  
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm 
  USE cell_base,     ONLY : omega, alat, bg
  USE division,      ONLY : fkbounds
  ! 
  implicit none
  ! 
  INTEGER, INTENT(IN) :: iqq
  !! Q-point index in selecq
  INTEGER, INTENT(IN) :: iq 
  !! Q-point index
  INTEGER, INTENT(IN) :: totq
  !! Total number of q-points in fsthick window
 
  !
  ! variables for collecting data from all pools in parallel case 
  !
  INTEGER :: iw, ik, ikk, ikq, ibnd, jbnd, fermicount
  INTEGER :: nksqtotf, lower_bnd, upper_bnd
  REAL(kind=DP), external :: efermig, dos_ef, wgauss
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, ww, dw, weight
  REAL(kind=DP) :: specfun_sum, esigmar0, tpiba_new
  REAL(kind=DP) :: fermi(nw_specfun)
  REAL(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:)
  REAL(kind=DP) :: kF, vF, fermiHEG, qin, wpl0, eps0, deltaeps, qcut, &
                   qsquared, qTF, dipole, rs, ekk1, degen
  REAL(kind=DP) :: q(3)
  !! The q-point in cartesian unit. 
  !
  ! loop over temperatures can be introduced
  !
  ! energy range and spacing for spectral function
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  !
  IF ( iqq == 1 ) THEN
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Electron Spectral Function in the Migdal Approximation")')
     WRITE(stdout,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick .lt. 1.d3 ) &
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
           'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
     !
  ENDIF
  !
  ! Fermi level and corresponding DOS
  !
  IF ( efermi_read ) THEN
    !
    ef0 = fermi_energy
    !
  ELSE
    !
    ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw,ngaussw,0,isk)
    ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
    ! in ephwann_shuffle
    !
  ENDIF
  !
  IF ( iqq == 1 ) THEN 
     WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
     WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! The total number of k points
  !
  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds( nksqtotf, lower_bnd, upper_bnd )
  !
  IF ( iqq == 1 ) THEN 
    IF ( .not. ALLOCATED(esigmar_all) ) ALLOCATE( esigmar_all(ibndmax-ibndmin+1, nksqtotf, nw_specfun) )
    IF ( .not. ALLOCATED(esigmai_all) ) ALLOCATE( esigmai_all(ibndmax-ibndmin+1, nksqtotf, nw_specfun) )
    esigmar_all(:,:,:) = zero
    esigmai_all(:,:,:) = zero
  ENDIF 
  !
  ! SP: Sum rule added to conserve the number of electron. 
  IF ( iqq == 1 ) THEN
    WRITE (stdout,'(5x,a)') 'The sum rule to conserve the number of electron is enforced.'
    WRITE (stdout,'(5x,a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
    WRITE (stdout,'(5x,a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
    WRITE (stdout,'(a)') ' '
  ENDIF
  !
!
  !nel      =  0.01    ! this should be read from input - # of doping electrons 
  !epsiHEG  =  12.d0   ! this should be read from input - # dielectric constant at zero doping  
  !meff     =  0.25    ! this should be read from input - effective mass 
  tpiba_new=  2.0d0 * pi / alat
  degen    =  1.0d0
  rs       =  (3.d0/(4.d0*pi*nel/omega/degen))**(1.d0/3.d0)*meff*degen ! omega is the unit cell volume in Bohr^3
!  rs       = (3.d0/(4.d0*pi*nel/omega/degen))**(1.d0/3.d0)*meff*degen/epsiHEG ! omega is the unit cell volume in Bohr^3
  kF       =  (3.d0*pi**2*nel/omega/degen )**(1.d0/3.d0)
  vF       =  1.d0/meff * (3.d0*pi**2*nel/omega/degen)**(1.d0/3.d0)
  fermiHEG =  1.d0/(2.d0*meff) * (3.d0*pi**2*nel/omega/degen)**(2.d0/3.d0) * 2.d0 ! [Ryd] multiplication by 2 converts from Ha to Ry
  qTF      =  (6.d0*pi*nel/omega/degen/(fermiHEG/2.d0))**(1.d0/2.d0)    ! [a.u.]
  wpl0     =  sqrt(4.d0*pi*nel/omega/meff/epsiHEG) * 2.d0         ! [Ryd] multiplication by 2 converts from Ha to Ryd
  wq       =  wpl0 ! [Ryd] 
  q(:)     =  xqf(:,iq) 
  CALL cryst_to_cart (1, q, bg, 1)
  qsquared =  (q(1)**2 + q(2)**2 + q(3)**2)
  qin      =  sqrt(qsquared)*tpiba_new
  qcut     =  wpl0 / vF  / tpiba_new / 2.d0 ! 1/2 converts from Ryd to Ha
  !qcut = qcut / 2.d0 ! phenomenological Landau damping
  ! 
  ! qcut2  = kF * ( sqrt( 1.d0 + wpl0 / fermiHEG) - 1.d0 ) / tpiba_new
  CALL get_eps_mahan (qin,rs,kF,eps0) ! qin should be in atomic units for Mahan formula
  !call get_eps_mahan (qin,qTF,kF,eps0) ! qin should be in atomic units for Mahan formula
  deltaeps = -(1.d0/(epsiHEG+eps0-1.d0)-1.d0/epsiHEG)
  !
  IF (iq .EQ. 1) THEN
    WRITE(stdout,'(12x," nel       = ", E15.10)') nel
    WRITE(stdout,'(12x," meff      = ", E15.10)') meff
    WRITE(stdout,'(12x," rs        = ", E15.10)') rs
    WRITE(stdout,'(12x," kF        = ", E15.10)') kF
    WRITE(stdout,'(12x," vF        = ", E15.10)') vF
    WRITE(stdout,'(12x," fermi_en  = ", E15.10)') fermiHEG
    WRITE(stdout,'(12x," qTF       = ", E15.10)') qTF
    WRITE(stdout,'(12x," wpl       = ", E15.10)') wpl0
    WRITE(stdout,'(12x," qcut      = ", E15.10)') qcut
    WRITE(stdout,'(12x," eps0      = ", E15.10)') eps0
    WRITE(stdout,'(12x," epsiHEG   = ", E15.10)') epsiHEG
    WRITE(stdout,'(12x," deltaeps  = ", E15.10)') deltaeps
  ENDIF
  !
  IF (sqrt(qsquared) < qcut) THEN
    !
    ! loop over all k points of the fine mesh
    !
    fermicount = 0 
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      ! (but in this case they are the same)
      !
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        !
        ! Bose occupation
        wgq = wgauss( -wq/eptemp, -99)
        wgq = wgq / ( one - two * wgq )
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          !  the energy of the electron at k (relative to Ef)
          ekk = etf (ibndmin-1+ibnd, ikk) - ef0
          !  
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            !  the fermi occupation for k+q
            ekk1 = etf (ibndmin-1+jbnd, ikk) - ef0
            ekq = etf (ibndmin-1+jbnd, ikq) - ef0
            wgkq = wgauss( -ekq/eptemp, -99)  
            !
            !computation of the dipole
            if (ibnd==jbnd) then
              if(sqrt(qsquared) .gt. 1d-6)then
                dipole = 1./(qsquared * tpiba_new * tpiba_new)
              else 
                dipole = 0.d0 
              endif
            else
              if (abs(ekq-ekk1) > 1d-6) then
                dipole = REAL( dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk) * &
                             conjg(dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk))/((ekk1-ekk)**2 + degaussw**2) )
              else 
                dipole = 0.d0
              endif
            endif 
            !
            g2 = dipole*4.d0*pi * (wq*deltaeps/2.d0)/omega * 2.d0 ! The q^-2 is cancelled by the q->0 limit of the dipole. See e.g., pg. 258 of Grosso Parravicini. 
            !
            DO iw = 1, nw_specfun
              !
              ww = wmin_specfun + dble (iw-1) * dw
              !
              weight = wqf(iq) * real (                                            &
                ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                  ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
              !
              esigmar_all(ibnd,ik+lower_bnd-1,iw) = esigmar_all(ibnd,ik+lower_bnd-1,iw) + g2 * weight 
              ! 
              ! SP : Application of the sum rule
              esigmar0 =  g2 *  wqf(iq) * real (                                   &
                ( (       wgkq + wgq ) / ( -( ekq - wq ) - ci * degaussw )  +  &
                  ( one - wgkq + wgq ) / ( -( ekq + wq ) - ci * degaussw ) ) )
              esigmar_all(ibnd,ik+lower_bnd-1,iw)=esigmar_all(ibnd,ik+lower_bnd-1,iw)-esigmar0
              !
              weight = wqf(iq) * aimag (                                           &
                ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                  ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
              !
              esigmai_all(ibnd,ik+lower_bnd-1,iw) = esigmai_all(ibnd,ik+lower_bnd-1,iw) + g2 * weight
              !
            ENDDO
            !
          ENDDO !jbnd
          !
        ENDDO !ibnd
        !
        !
      ENDIF ! endif  fsthick
      !
    ENDDO ! end loop on k
    !
  ENDIF
  !
  ! The k points are distributed among pools: here we collect them
  !
  IF ( iqq == totq ) THEN
    !
    ALLOCATE ( xkf_all      ( 3,       nkqtotf ), &
               etf_all      ( nbndsub, nkqtotf ) )
    xkf_all(:,:) = zero
    etf_all(:,:) = zero
    !
#if defined(__MPI)
    !
    ! note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
    CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
    CALL mp_sum( esigmar_all, inter_pool_comm )
    CALL mp_sum( esigmai_all, inter_pool_comm )
    CALL mp_sum( fermicount, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    !
#endif
    !
    ! Output electron spectral function here after looping over all q-points 
    ! (with their contributions summed in a etc.)
    !
    WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
    !
    ! construct the trace of the spectral function (assume diagonal selfenergy
    ! and constant matrix elements for dipole transitions)
    !
    IF (.not. ALLOCATED (a_all)) ALLOCATE ( a_all(nw_specfun, nksqtotf) )
    a_all(:,:) = zero
    !
    IF (me_pool == 0) then
      OPEN(unit=iospectral,file='specfun.plself') 
      OPEN(unit=iospectral_sup,file='specfun_sup.plself') 
    ENDIF
    IF (me_pool == 0) then
      WRITE(iospectral, '(/2x,a/)') '#Electron-plasmon spectral function (meV)'
      WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electron-plasmon self-energy (meV)' 
    ENDIF
    IF (me_pool == 0) then
      WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
      WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]   &
&         Real Sigma[meV]  Im Sigma[meV]'
    ENDIF
    !
    DO ik = 1, nksqtotf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f12.7)') ik, xkf_all (:,ikk)
      WRITE(stdout,'(5x,a)') repeat('-',67)
      !
      DO iw = 1, nw_specfun
        !
        ww = wmin_specfun + dble (iw-1) * dw
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          !  the energy of the electron at k
          ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
          !
          a_all(iw,ik) = a_all(iw,ik) + abs( esigmai_all(ibnd,ik,iw) ) / pi / &
             ( ( ww - ekk - esigmar_all(ibnd,ik,iw) )**two + (esigmai_all(ibnd,ik,iw) )**two )
          !
        ENDDO
        !
        WRITE(stdout, 103) ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
        !
      ENDDO
      !
      WRITE(stdout,'(5x,a/)') repeat('-',67)
      !
    ENDDO
    !
    DO ik = 1, nksqtotf
      !
      ! The spectral function should integrate to 1 for each k-point
      specfun_sum = 0.0
      ! 
      DO iw = 1, nw_specfun
         !
         ww = wmin_specfun + dble (iw-1) * dw
         fermi(iw) = wgauss(-ww/eptemp, -99) 
        ! WRITE(stdout,'(2x,i7,2x,f12.4,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
         !
         specfun_sum = specfun_sum + a_all(iw,ik)* fermi(iw) * dw !/ ryd2mev
         !
       IF (me_pool == 0) &
         WRITE(iospectral,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
          !
      ENDDO
      !
      IF (me_pool == 0) &
        WRITE(iospectral,'(a)') ' '
      IF (me_pool == 0) &
        WRITE(iospectral,'(2x,a,2x,e12.5)') '# Integrated spectral function ',specfun_sum
      !
    ENDDO
    !
    IF (me_pool == 0) CLOSE(iospectral)
    !
    DO ibnd = 1, ibndmax-ibndmin+1
      !
      DO ik = 1, nksqtotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        !  the energy of the electron at k
        ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
        !
        DO iw = 1, nw_specfun
          !
          ww = wmin_specfun + dble (iw-1) * dw
          WRITE(stdout,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
            ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd,ik,iw),&
            ryd2mev * esigmai_all(ibnd,ik,iw)
          ! 
          IF (me_pool == 0) &
          WRITE(iospectral_sup,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
            ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd,ik,iw),&
            ryd2mev * esigmai_all(ibnd,ik,iw)
          !
        ENDDO
        !
      ENDDO
      !
      WRITE(stdout,*) ' '
      !
    ENDDO
    !
    IF (me_pool == 0) CLOSE(iospectral_sup)
    !
    IF ( ALLOCATED(xkf_all) )      DEALLOCATE( xkf_all )
    IF ( ALLOCATED(etf_all) )      DEALLOCATE( etf_all )
    IF ( ALLOCATED(esigmar_all) )  DEALLOCATE( esigmar_all )
    IF ( ALLOCATED(esigmai_all) )  DEALLOCATE( esigmai_all )
    IF ( ALLOCATED(a_all) )        DEALLOCATE( a_all )
    !
  ENDIF
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  103 FORMAT(5x,'ik = ',i7,'  w = ',f9.4,' eV   A(k,w) = ',e12.5,' meV^-1')
  !
  RETURN
  !
  END SUBROUTINE spectral_func_pl_q
