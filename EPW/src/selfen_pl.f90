  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_pl_q ( iq )
  !-----------------------------------------------------------------------
  !
  !  Compute the imaginary part of the electron self energy due to electron-
  !  plasmon interaction. 
  ! 
  !  The coupling coefficients have been evaluated analytically employing a 
  !  Lindhard function model for the dielectric function contribution due to 
  !  the extrinsic carriers.
  ! 
  !  There are 3 parameters that the users should provide in the input: 
  !    - DOS effective mass; 
  !    - carrier concentration (Only for doped semiconductors, it shouldnt be used for insulators);
  !    - epsilon_infinity (e.g, from exp. or from RPA).
  !
  !  F. Caruso and S. Ponce - 2017 
  !
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : linewidth_elself
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, lrepmatf, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy, & 
                            nel, meff, epsiHEG 
  USE pwcom,         ONLY : ef
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, &
                            epf17, wkf, nkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew, nqf, & 
                            xqf, dmef  
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm, my_pool_id
  use cell_base,     ONLY : omega, alat
  !USE mp_world,      ONLY : mpime
  implicit none
  !
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  INTEGER :: nksqtotf, lower_bnd, upper_bnd
  INTEGER :: n
  !! Integer for the degenerate average over eigenstates
  REAL(kind=DP) :: tmp
  !! Temporary variable to store real part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp2
  !! Temporary variable to store imag part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp3
  !! Temporary variable to store Z for the degenerate average
  REAL(kind=DP) :: ekk2
  !! Temporary variable to the eigenenergies for the degenerate average
  REAL(kind=DP) :: sigmar_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the real-part of Sigma 
  REAL(kind=DP) :: sigmai_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the imag-part of Sigma 
  REAL(kind=DP) :: zi_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the Z
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, weight,  &
                   w0g1, w0g2, inv_wq, inv_eptemp0, g2_tmp,&
                   inv_degaussw, tpiba_new
  REAL(kind=DP), external :: efermig, dos_ef, wgauss, w0gauss
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:), etf_all(:,:)
  REAL(kind=DP) :: kF, vF, fermiHEG, qin, wpl0, eps0, deltaeps, qcut, qcut2
  REAL(kind=DP) :: qsquared, qTF, dipole, rs, ekk1, degen
  !REAL(kind=DP) :: Nel, epsiHEG, meff, kF, vF, fermiHEG, qin, wpl0, eps0, deltaeps, qcut, qcut2, qsquared, qTF, dipole, rs, ekk1
  ! loop over temperatures can be introduced
  !
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0  = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
  !
  IF ( iq .eq. 1 ) THEN
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Electron-plasmon Self-Energy in the Migdal Approximation")')
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
    ef0 = efnew ! Fermi energy is recalculated on the fine mesh!! 
    !ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw,ngaussw,0,isk)
    ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
    ! in ephwann_shuffle
    !
  ENDIF
  !
  IF ( iq .eq. 1 ) THEN 
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
  IF ( iq .eq. 1 ) THEN 
    IF ( .not. ALLOCATED (sigmar_all) ) ALLOCATE( sigmar_all(ibndmax-ibndmin+1, nksqtotf) )
    IF ( .not. ALLOCATED (sigmai_all) ) ALLOCATE( sigmai_all(ibndmax-ibndmin+1, nksqtotf) )
    IF ( .not. ALLOCATED (zi_all) )     ALLOCATE( zi_all(ibndmax-ibndmin+1, nksqtotf) )
    sigmar_all(:,:) = zero
    sigmai_all(:,:) = zero
    zi_all(:,:) = zero
  ENDIF
  !
  !
  !nel      =  0.01    ! this should be read from input - # of doping electrons 
  !epsiHEG  =  12.d0   ! this should be read from input - # dielectric constant at zero doping  
  !meff     =  0.25    ! this should be read from input - effective mass 
  !
  tpiba_new = 2.0d0 * pi / alat
  degen    = 1.0d0
  rs       = (3.d0/(4.d0*pi*nel/omega/degen))**(1.d0/3.d0)*meff*degen ! omega is the unit cell volume in Bohr^3
  kF       = (3.d0*pi**2*nel/omega/degen )**(1.d0/3.d0) 
  vF       =  1.d0/meff * (3.d0*pi**2*nel/omega/degen)**(1.d0/3.d0) 
  fermiHEG =  1.d0/(2.d0*meff) * (3.d0*pi**2*nel/omega/degen)**(2.d0/3.d0) * 2.d0 ! [Ryd] multiplication by 2 converts from Ha to Ry
  qTF      =  (6.d0*pi*nel/omega/degen/(fermiHEG/2.d0))**(1.d0/2.d0)    ! [a.u.]
  wpl0     =  sqrt(4.d0*pi*nel/omega/meff/epsiHEG) * 2.d0         ! [Ryd] multiplication by 2 converts from Ha to Ryd
  wq       =  wpl0 ! [Ryd] 
  qsquared = (xqf(1,iq)**2 + xqf(2,iq)**2 + xqf(3,iq)**2)
  qin      =  sqrt(qsquared)*tpiba_new
  qcut     = wpl0 / vF / tpiba_new / 2.d0    ! 1/2 converts from Ryd to Ha
  !
  !if (.true.) qcut = qcut / 2.d0 ! renorm to account for Landau damping 
  !
  CALL get_eps_mahan (qin,rs,kF,eps0) ! qin should be in atomic units for Mahan formula
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
  !if (sqrt(qsquared) > qcut) wqf(iq) = 0.d0
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
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
       ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        wgq = wgauss( -wq*inv_eptemp0, -99)
        wgq = wgq / ( one - two * wgq )
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          !  the energy of the electron at k (relative to Ef)
          ekk = etf (ibndmin-1+ibnd, ikk) - ef0
          !
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            ekk1 = etf (ibndmin-1+jbnd, ikk) - ef0
            ekq  = etf (ibndmin-1+jbnd, ikq) - ef0
            wgkq = wgauss( -ekq*inv_eptemp0, -99)  
            !
!            if ( abs(ekq-ekk1) > 1d-8 ) then
!              !dipole = (dmef(1, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)*conjg(dmef(1, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)) +  &
!              !          dmef(2, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)*conjg(dmef(2, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)) +  &
!              !          dmef(3, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)*conjg(dmef(3, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)))   &
!              !          /3.d0/(ekq-ekk)**2 
!
!              !dipole = dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk)*conjg(dmef(1, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk))/(ekq-ekk)**2 
!              dipole = dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk)/((ekk1-ekk)**2 + (degaussw/10)**2  )
!              !dipole = dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk)*conjg(dmef(1, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk))/((ekk1-ekk)**2 + (degaussw/10)**2  )
!                         !abs (dmef(3, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk))**2 )/3.d0/((ekq-ekk)**2 + degaussw**2  )
!              ! THIS EXPRESSION NEEDS TO BE DIVIDED BY THE TRANSITION ENERGIES^2 (SEE GROSSO PARRAVICINI PG 258 FOR THE CORRECT EXPRESSION IN THE Q->0 LIMIT!!!!!!!!!1)
!            !else
!            !  dipole = 1.d0
!            endif 
              
            !computation of the dipole
            IF ( ibnd==jbnd ) THEN
              IF (sqrt(qsquared) .gt. 1d-8) THEN
                dipole = 1./(qsquared * tpiba_new * tpiba_new)
              ELSE
                dipole = 0.d0 
              ENDIF
            ELSE
              IF (abs(ekk-ekk1) > 1d-8) THEN
                dipole = dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk)*&
                             conjg(dmef(1,ibndmin-1+jbnd,ibndmin-1+ibnd,ikk))/((ekk1-ekk)**2 + degaussw**2)
              ELSE 
                dipole = 0.d0
              ENDIF
            ENDIF
            !
            ! this approximates the dipoles as delta_ij
            !if (.false.) then
            !  if (ibnd==jbnd .and. sqrt(qsquared) .gt. 1d-8 ) then
            !    dipole = 1./(qsquared * tpiba_new * tpiba_new)
            !  else 
            !    dipole = 0.d0
            !  endif
            !endif
            !if (.true.) then
            IF ( abs(dipole * (qsquared * tpiba_new * tpiba_new)) .gt.1 ) THEN
              dipole = 1./(qsquared*tpiba_new*tpiba_new)
            ENDIF
            !endif
            !
            !if (ik .eq. 1) WRITE(stdout,'(/12x," dipole    =", f12.7)') dipole
            g2 = dipole*4.d0*pi * (wq*deltaeps/2.d0)/omega * 2.d0 ! The q^-2 is cancelled by the q->0 limit of the dipole. See e.g., pg. 258 of Grosso Parravicini. 
            !if (ik .eq.10) WRITE(stdout,'(/5x," g2 =", f12.7)') g2
              !g2 = dipole * 4.d0*pi / (qsquared*tpiba_new*tpiba_new) * ( wq * deltaeps / 2.d0 ) / omega * 2.d0 ! The last is the spin! 
!            else
            ! g2 = (dmef(1, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)**2 +  &
            !       dmef(2, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)**2 +  &
            !       dmef(3, ibndmin-1+jbnd, ibndmin-1+ibnd, ikk)**2 )/3 & 
            !      * 4.d0*pi * ( wq * deltaeps / 2.d0 ) * 2.d0 ! The last is the spin! 
!              g2 = 0.d0 ! replace here the q -> 0 limit 
!            endif
            !g2 = 1.d0
            !g2 = dipole 
            !
            !
            !  the fermi occupation for k+q
            !
            ! here we take into account the zero-point sqrt(hbar/2M\omega)
            ! with hbar = 1 and M already contained in the eigenmodes
            ! g2 is Ry^2, wkf must already account for the spin factor
            !
            !g2 = (abs(epf (jbnd, ibnd))**two)*inv_wq*g2_tmp
            !g2 = 1.d0 
            !
            ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
            ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
            ! (Many-Particle Physics, 3rd edition)
            ! 
            weight = wqf(iq) * real (                                                   &
                    ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                      ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
            !
!            if (sqrt(qsquared) .lt. 1d-6) then
!            !if (iq.eq.1) then
!              !
!              weight = ( 3.d0 * wqf(iq) / 4.d0 / pi / omega )**(1.d0/3.d0)/pi  *          &
!               real ( ( (     wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
!                      ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!              !
!            endif
            !
!              ecutse needs to be defined if it's used 
!@             if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
            !
            sigmar_all(ibnd,ik+lower_bnd-1) = sigmar_all(ibnd,ik+lower_bnd-1) + g2 * weight
            !
            ! Logical implementation
!            weight = wqf(iq) * aimag (                                                  &
!                    ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
!                      ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!@            if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
            !
            ! Delta implementation 
            w0g1=w0gauss( (ekk-ekq+wq)/degaussw, 0) /degaussw
            w0g2=w0gauss( (ekk-ekq-wq)/degaussw, 0) /degaussw
            weight = pi * wqf(iq) * ( (wgkq+wgq)*w0g1 + (one-wgkq+wgq)*w0g2 )
            !
!            if (sqrt(qsquared) .lt. 1d-6) then
!              !
!              !weight = ( 3.d0 * wqf(iq) / 4.d0 / pi / omega )**(1.d0/3.d0)/pi * &
!              ! aimag ( ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
!              !           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!              weight = ( 3.d0 * wqf(iq) / 4.d0 / pi / omega )**(1.d0/3.d0)/pi  * &   
!                       pi * ( (wgkq+wgq)*w0g1 + (one-wgkq+wgq)*w0g2 )
!              !
!            endif
            !
            sigmai_all(ibnd,ik+lower_bnd-1) = sigmai_all(ibnd,ik+lower_bnd-1) + g2 * weight
            !
            ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
            !
            weight = wqf(iq) * &
                    ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                             ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                      ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                             ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )  
!@            if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
            !
            zi_all(ibnd,ik+lower_bnd-1) = zi_all(ibnd,ik+lower_bnd-1) + g2 * weight
            !
          ENDDO !jbnd
          !
        ENDDO !ibnd
        !
      ENDIF ! endif  fsthick
      !
    ENDDO ! end loop on k
  endif 
  !
  ! The k points are distributed among pools: here we collect them
  !
  IF ( iq .eq. nqtotf ) THEN
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
     CALL mp_sum( sigmar_all, inter_pool_comm )
     CALL mp_sum( sigmai_all, inter_pool_comm )
     CALL mp_sum( zi_all, inter_pool_comm )
     CALL mp_sum(fermicount, inter_pool_comm)
     CALL mp_barrier(inter_pool_comm)
     !
#else
     !
     xkf_all = xkf
     etf_all = etf
     !
#endif
     !
     ! Average over degenerate eigenstates:
     WRITE(stdout,'(5x,"Average over degenerate eigenstates is performed")')
     ! 
     DO ik = 1, nksqtotf
       ikk = 2 * ik - 1
       ikq = ikk + 1
       ! 
       DO ibnd = 1, ibndmax-ibndmin+1
         ekk = etf_all (ibndmin-1+ibnd, ikk)
         n = 0
         tmp = 0.0_DP
         tmp2 = 0.0_DP
         tmp3 = 0.0_DP
         !sigmar_tmp(:) = zero
         DO jbnd = 1, ibndmax-ibndmin+1
           ekk2 = etf_all (ibndmin-1+jbnd, ikk)
           IF ( ABS(ekk2-ekk) < eps6 ) THEN
             n = n + 1
             tmp =  tmp + sigmar_all (jbnd,ik)
             tmp2 =  tmp2 + sigmai_all (jbnd,ik)
             tmp3 =  tmp3 + zi_all (jbnd,ik)
           ENDIF
           ! 
         ENDDO ! jbnd
         sigmar_tmp(ibnd) = tmp / float(n)
         sigmai_tmp(ibnd) = tmp2 / float(n)
         zi_tmp(ibnd) = tmp3 / float(n)
         !
       ENDDO ! ibnd
       sigmar_all (:,ik) = sigmar_tmp(:)
       sigmai_all (:,ik) = sigmai_tmp(:)
       zi_all (:,ik)  = zi_tmp(:)
       ! 
     ENDDO ! nksqtotf
     !
     ! Output electron SE here after looping over all q-points (with their contributions 
     ! summed in sigmar_all, etc.)
     !
     WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     !
     ! Write to file
     OPEN(unit=linewidth_elself,file='linewidth.plself')
     WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
     WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sgima)(meV)'
     ! 
     DO ik = 1, nksqtotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout,'(/5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:,ikk)
        WRITE(stdout,'(5x,a)') repeat('-',67)
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          ! note that ekk does not depend on q 
          ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
          !
          ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
          zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
          !
          WRITE(stdout, 102) ibndmin-1+ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
                             ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik), one/zi_all(ibnd,ik)-one
          WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
          WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
          WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
          WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_all(ibnd,ik)
          !
        ENDDO
        WRITE(stdout,'(5x,a/)') repeat('-',67)
        !
     ENDDO
     !
     DO ibnd = 1, ibndmax-ibndmin+1
        !
        DO ik = 1, nksqtotf
           !
           ikk = 2 * ik - 1
           ikq = ikk + 1
           !
           ! note that ekk does not depend on q 
           ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
           !
           ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
           !zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
           !
           !WRITE(stdout,'(2i9,5f12.4)') ik, ibndmin-1+ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd,ik), &
           !                             ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik),  one/zi_all(ibnd,ik)-one
           ! 
        ENDDO
        !
        !WRITE(stdout,'(a)') '  '
        !
     ENDDO
     !
     CLOSE(linewidth_elself)
     !
     IF ( ALLOCATED(xkf_all) )      DEALLOCATE( xkf_all )
     IF ( ALLOCATED(etf_all) )      DEALLOCATE( etf_all )
     IF ( ALLOCATED(sigmar_all) )   DEALLOCATE( sigmar_all )
     IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
     IF ( ALLOCATED(zi_all) )       DEALLOCATE( zi_all )
     IF ( ALLOCATED(sigmai_mode) )   DEALLOCATE( sigmai_mode )
     !
  ENDIF 
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  103 FORMAT(5x,'k( ',i7,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6)
  !
  RETURN
  !
  END SUBROUTINE selfen_pl_q
  !
  ! 
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  SUBROUTINE get_eps_mahan (q,rs,kF,eps0)
  !subroutine get_eps_mahan (q,qTF,kF,eps0)
  !
  ! Based on Eq. 5.166 of Mahan 2000. 
  !
  USE kinds,        ONLY : DP
  USE constants_epw, ONLY : pi
  implicit none
  ! 
  REAL(kind=DP),    intent (OUT) :: eps0
  REAL(kind=DP),    intent (IN)  :: q
  REAL(kind=DP),    intent (IN)  :: rs
  REAL(kind=DP),    intent (IN)  :: kF
  ! 
! internal
  REAL(kind=DP) :: eta,x,alpha 
  ! 
  eta   = 1.d-6
  alpha = (4.d0/9.d0/pi)**(1.d0/3.d0)
  ! 
  if ( abs(q).gt. 1.d-10) then
    x    = q / 2.d0 / kF 
    eps0 = 1.d0 + (1.d0-x**2)/(2.d0*x)*log (abs((1.d0+x)/(1.d0-x))) 
    eps0 = 1.d0 + alpha * rs / 2.d0 /pi / x**2 * eps0
  else 
    x    = (q + eta) / 2.d0 / kF 
    eps0 = 1.d0 + (1.d0-x**2)/(2.d0*x)*log ( abs ((1.d0+x)/(1.d0-x)) ) 
    eps0 = 1.d0 + alpha * rs / 2.d0 /pi / x**2 * eps0
  endif

!  if ( abs(q).gt. 1.d-10) then
!    x    = q / 2.d0 / kF 
!    eps0 = 1.d0 + (1.d0-x**2)/(2.d0*x)*log (abs((1.d0+x)/(1.d0-x))) 
!    eps0 = 1.d0 + qTF**2/(2.d0*q**2) * eps0
!  else 
!    x    = (q + eta) / 2.d0 / kF 
!    eps0 = 1.d0 + (1.d0-x**2)/(2.d0*x)*log ( abs ((1.d0+x)/(1.d0-x)) ) 
!    eps0 = 1.d0 + qTF**2/(2.d0*(q+eta)**2) * eps0
!  endif

  END SUBROUTINE get_eps_mahan
!--------------------------------------------------------------------------
