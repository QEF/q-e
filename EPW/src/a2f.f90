  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino 
  !                                                                           
  ! This file is distributed under the terms of the GNU General Public       
  ! License. See the file `LICENSE' in the root directory of the          
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .        
  !
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_a2f
  !-----------------------------------------------------------------------
  !!
  !!  Compute the Eliasberg spectral function
  !!  in the Migdal approximation. 
  !!  
  !!  If the q-points are not on a uniform grid (i.e. a line)
  !!  the function will not be correct
  !!  
  !!  02/2009 works in serial on ionode at the moment.  can be parallelized
  !!  03/2009 added transport spectral function -- this involves a v_k dot v_kq term 
  !!          in the quantities coming from selfen_phon.f90.  Not fully implemented  
  !!  10/2009 the code is transitioning to 'on-the-fly' phonon selfenergies
  !!          and this routine is not currently functional
  !!  10/2015 RM: added calcution of Tc based on Allen-Dynes formula 
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes
  USE cell_base, ONLY : omega
  USE epwcom,    ONLY : degaussq, delta_qsmear, nqsmear, nqstep, nsmear, eps_acustic, & 
                        delta_smear, degaussw, fsthick
  USE elph2,     ONLY : nqtotf, wf, wqf, lambda_all, lambda_v_all
  USE constants_epw, ONLY : ryd2mev, ryd2ev, kelvin2eV, two, zero, kelvin2Ry, pi
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_world,  ONLY : mpime, world_comm
  USE io_global, ONLY : ionode_id
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iua2ffil, iudosfil, iua2ftrfil, iures
  USE io_files,  ONLY : prefix
  implicit none
  !
  integer       :: imode, iq, iw, ismear, isig, i, itemp
  real(kind=DP) :: weight
  REAL(KIND=DP) :: temp
  REAL(KIND=DP) :: n
  REAL(KIND=DP) :: be
  REAL(KIND=DP) :: prefact
  real(kind=DP) :: lambda_tot, lambda_tr_tot
  real(kind=DP) :: iomega, sigma, a2F_tmp, a2F_tr_tmp, om_max, dw, w0, l, l_tr, tc, mu
  real(kind=DP), allocatable :: a2F(:,:), a2F_tr(:,:), l_a2F(:), l_a2F_tr(:), dosph(:,:), logavg(:), rho(:,:)
  real(kind=DP), external :: w0gauss
  CHARACTER (len=256) :: fila2f_suffix, fila2ftr, fildos, filres
  !
  !
  CALL start_clock('a2F')
  IF (mpime .eq. ionode_id) THEN
  !
  DO isig = 1, nsmear
    !
    IF ( isig .lt. 10 ) THEN
       WRITE(fila2f_suffix,'(a,a6,i1)') TRIM(prefix),'.a2f.0', isig
    ELSE 
       WRITE(fila2f_suffix,'(a,a5,i2)') TRIM(prefix),'.a2f.', isig
    ENDIF
    OPEN (unit = iua2ffil, file = fila2f_suffix, form = 'formatted')
    !
    IF ( isig .lt. 10 ) THEN
       WRITE(fila2ftr,'(a,a9,i1)') TRIM(prefix),'.a2f_tr.0', isig
    ELSE
       WRITE(fila2ftr,'(a,a8,i2)') TRIM(prefix),'.a2f_tr.', isig
    ENDIF
    OPEN (unit = iua2ftrfil, file = fila2ftr, form = 'formatted')
    !
    IF ( isig .lt. 10 ) THEN
       WRITE(filres,'(a,a6,i1)') TRIM(prefix),'.res.0', isig
    ELSE
       WRITE(filres,'(a,a5,i2)') TRIM(prefix),'.res.', isig
    ENDIF
    OPEN (unit = iures, file = filres, form = 'formatted')
    !
    IF ( isig .lt. 10 ) THEN
       WRITE(fildos,'(a,a8,i1)') TRIM(prefix),'.phdos.0', isig
    ELSE
       WRITE(fildos,'(a,a7,i2)') TRIM(prefix),'.phdos.', isig
    ENDIF
    OPEN (unit = iudosfil, file = fildos, form = 'formatted')
    !
    WRITE(stdout,'(/5x,a)') REPEAT('=',67)
    WRITE(stdout,'(5x,"Eliashberg Spectral Function in the Migdal Approximation")') 
    WRITE(stdout,'(5x,a/)') REPEAT('=',67)
    !
    IF ( .not. ALLOCATED(a2F) )    ALLOCATE( a2F(nqstep, nqsmear) )
    IF ( .not. ALLOCATED(a2F_tr) ) ALLOCATE( a2F_tr(nqstep, nqsmear) )
    IF ( .not. ALLOCATED(dosph) )  ALLOCATE( dosph(nqstep, nqsmear) )
    IF ( .not. ALLOCATED(l_a2F) )  ALLOCATE( l_a2F(nqsmear) )
    IF ( .not. ALLOCATED(l_a2F_tr) )  ALLOCATE( l_a2F_tr(nqsmear) )
    IF ( .not. ALLOCATED(logavg) ) ALLOCATE( logavg(nqsmear) )
    ! 
    ! The resitivity is computed for temperature between 0K-1000K by step of 10
    ! This is hardcoded and needs to be changed here if one wants to modify it
    IF ( .not. ALLOCATED(rho) ) ALLOCATE( rho(100, nqsmear) )
    !
    !om_max = ( MAXVAL( wf(:,:) ) - MINVAL( wf(:,:) ) ) + 5.d0/ryd2mev
    !om_max = MAXVAL( wf(:,:) ) + 1.d0 / ryd2mev
    !dw = om_max / dble(nqstep-1)
    om_max = 1.1d0 * MAXVAL( wf(:,:) ) ! increase by 10%
    dw = om_max / dble(nqstep)
    !
    lambda_tot = zero
    l_a2F(:) = zero
    a2F(:,:) = zero
    lambda_tr_tot = zero
    l_a2F_tr(:) = zero
    a2F_tr(:,:) = zero
    dosph(:,:) = zero
    logavg(:) = zero
    !
    DO ismear = 1, nqsmear
       !
       sigma = degaussq + (ismear-1) * delta_qsmear
       !
       DO iw = 1, nqstep  ! loop over points on the a2F(w) graph
          !
          !iomega = dble(iw-1) * dw ! step through the frequncies we wish to plot
          iomega = dble(iw) * dw ! step through the frequncies we wish to plot
          !
          DO iq = 1, nqtotf ! loop over q-points 
             !
             DO imode = 1, nmodes ! loop over modes
                !
                w0 = wf(imode,iq)
                !
                IF ( w0 .gt. eps_acustic ) THEN 
                   !
                   l  = lambda_all(imode,iq,isig)
                   IF (lambda_all(imode,iq,isig) .lt. 0.d0)  l = 0.d0 ! sanity check
                   ! 
                   a2F_tmp    = wqf(iq) * w0 * l / two
                   !
                   weight = w0gauss( (iomega - w0)/sigma, 0 ) / sigma
                   a2F(iw,ismear) = a2F(iw,ismear) + a2F_tmp * weight
                   dosph(iw,ismear)  = dosph(iw,ismear) + wqf(iq) * weight
                   !
                   l_tr = lambda_v_all(imode,iq,isig)
                   IF (lambda_v_all(imode,iq,isig) .lt. 0.d0)  l_tr = 0.d0 !sanity check
                   ! 
                   a2F_tr_tmp = wqf(iq) * w0 * l_tr / two
                   !
                   a2F_tr(iw,ismear) = a2F_tr(iw,ismear) + a2F_tr_tmp * weight
                   !
                ENDIF
                !
             ENDDO
             !
          ENDDO
          !
          !  output a2F
          !
          IF (ismear .eq. nqsmear) WRITE (iua2ffil,   '(f12.7, 15f12.7)') iomega*ryd2mev, a2F(iw,:)
          IF (ismear .eq. nqsmear) WRITE (iua2ftrfil, '(f12.7, 15f12.7)') iomega*ryd2mev, a2F_tr(iw,:)
          IF (ismear .eq. nqsmear) WRITE (iudosfil,   '(f12.7, 15f12.7)') iomega*ryd2mev, dosph(iw,:)/ryd2mev
          !
          ! do the integral 2 int (a2F(w)/w dw)
          !
          !IF (iomega .gt. eps_acustic) & 
          l_a2F(ismear) = l_a2F(ismear) + two * a2F(iw,ismear) / iomega * dw
          l_a2F_tr(ismear) = l_a2F_tr(ismear) + two * a2F_tr(iw,ismear) / iomega * dw
          logavg(ismear) = logavg(ismear) + two *  a2F(iw,ismear) * log(iomega) / iomega * dw
          !
       ENDDO
       !
       logavg(ismear) = exp( logavg(ismear) / l_a2F(ismear) )
       !
    ENDDO
    !
    DO iq = 1, nqtotf ! loop over q-points 
       DO imode = 1, nmodes ! loop over modes
          IF (lambda_all(imode,iq,isig) .gt. 0.d0 .and. wf(imode,iq) .gt. eps_acustic ) & 
             lambda_tot = lambda_tot + wqf(iq) * lambda_all(imode,iq,isig)
          IF (lambda_v_all(imode,iq,isig) .gt. 0.d0 .and. wf(imode,iq) .gt. eps_acustic ) &
             lambda_tr_tot = lambda_tr_tot + wqf(iq) * lambda_v_all(imode,iq,isig)
       ENDDO
    ENDDO
    WRITE (stdout,'(5x,a,f12.7)') "lambda : ", lambda_tot
    WRITE (stdout,'(5x,a,f12.7)') "lambda_tr : ", lambda_tr_tot
    WRITE (stdout,'(a)') " "
    !
    !
    ! Allen-Dynes estimate of Tc for ismear = 1
    !
    WRITE(stdout,'(5x,a,f12.7,a)') "Estimated Allen-Dynes Tc"
    WRITE (stdout,'(a)') " "
    WRITE(stdout,'(5x,a,f12.7,a,f12.7)') "logavg = ", logavg(1), " l_a2F = ", l_a2F(1)
    DO i = 1, 6
       !
       mu = 0.1d0 + 0.02d0 * dble(i-1)
       tc = logavg(1) / 1.2d0 * exp( - 1.04d0 * ( 1.d0 + l_a2F(1) ) &
                             / ( l_a2F(1) - mu * ( 1.d0 + 0.62d0 * l_a2F(1) ) ))
       ! tc in K
       !
       tc = tc * ryd2ev / kelvin2eV
       !SP: IF Tc is too big, it is not physical
       IF (tc < 1000.0 ) THEN
         WRITE(stdout,'(5x,a,f6.2,a,f22.12,a)') "mu = ", mu, " Tc = ", tc, " K"
       ENDIF 
       !
    ENDDO
    ! 
    rho(:,:) = zero
    ! Now compute the Resistivity of Metal using the Ziman formula
    ! rho(T,smearing) = 4 * pi * me/(n * e**2 * kb * T) int dw hbar w a2F_tr(w,smearing) n(w,T)(1+n(w,T))
    ! n is the number of electron per unit volume and n(w,T) is the Bose-Einstein distribution
    ! Usually this means "the number of electrons that contribute to the mobility" and so it is typically 8 (full shell)
    ! but not always. You might want to check this. 
    ! 
    n = 8.0 / omega
    !print*,'omega ',omega
    WRITE (iures, '(a)') '# Temperature [K]                Resistivity [micro Ohm cm] for different Phonon smearing (meV)        '  
    WRITE (iures, '("#     ", 15f12.7)') ( (degaussq+(ismear-1)*delta_qsmear)*ryd2mev,ismear=1,nqsmear )
    DO ismear = 1, nqsmear
      DO itemp = 1, 100 ! Per step of 10K
        temp = itemp * 10 * kelvin2Ry
        ! omega is the volume of the primitive cell in a.u.  
        ! 
        prefact = 4.0 * pi / ( temp * n )
        DO iw = 1, nqstep  ! loop over points on the a2F(w)
          ! 
          iomega = dble(iw) * dw
          be = 1.0/(exp(iomega/temp)-1); 
          ! Perform the integral with rectangle. 
          rho(itemp,ismear) = rho(itemp,ismear) + prefact * iomega * a2F_tr(iw,ismear) * be * (1.0 + be) * dw  
          ! 
        ENDDO
        ! From a.u. to micro Ohm cm
        ! Conductivity 1 a.u. = 2.2999241E6 S/m
        ! Now to go from Ohm*m to micro Ohm cm we need to multiply by 1E8 
        rho(itemp,ismear) = rho(itemp,ismear) * 1E8 / 2.2999241E6
        IF (ismear .eq. nqsmear) WRITE (iures, '(i8, 15f12.7)') itemp * 10, rho(itemp,:)
      ENDDO
    ENDDO 
    CLOSE(iures)
    !
    WRITE(iua2ffil,*) "Integrated el-ph coupling"
    WRITE(iua2ffil,'("  #         ", 15f12.7)') l_a2F(:)
    WRITE(iua2ffil,*) "Phonon smearing (meV)"
    WRITE(iua2ffil,'("  #         ", 15f12.7)') ( (degaussq+(ismear-1)*delta_qsmear)*ryd2mev,ismear=1,nqsmear )
    WRITE(iua2ffil,'(" Electron smearing (eV)", f12.7)') ((isig-1)*delta_smear+degaussw)*ryd2ev
    WRITE(iua2ffil,'(" Fermi window (eV)", f12.7)') fsthick*ryd2ev
    WRITE(iua2ffil,'(" Summed el-ph coupling ", f12.7)') lambda_tot
    CLOSE(iua2ffil)
    !
    WRITE(iua2ftrfil,*) "Integrated el-ph coupling"
    WRITE(iua2ftrfil,'("  #         ", 15f12.7)') l_a2F_tr(:)
    WRITE(iua2ftrfil,*) "Phonon smearing (meV)"
    WRITE(iua2ftrfil,'("  #         ", 15f12.7)') ( (degaussq+(ismear-1)*delta_qsmear)*ryd2mev,ismear=1,nqsmear )
    WRITE(iua2ftrfil,'(" Electron smearing (eV)", f12.7)') ((isig-1)*delta_smear+degaussw)*ryd2ev
    WRITE(iua2ftrfil,'(" Fermi window (eV)", f12.7)') fsthick*ryd2ev
    WRITE(iua2ftrfil,'(" Summed el-ph coupling ", f12.7)') lambda_tot
    CLOSE(iua2ftrfil)
    !
    CLOSE(iudosfil)
    !
    IF ( ALLOCATED(l_a2F) )     DEALLOCATE(l_a2F)
    IF ( ALLOCATED(l_a2F_tr) )  DEALLOCATE(l_a2F_tr)
    IF ( ALLOCATED(a2F) )       DEALLOCATE(a2F)
    IF ( ALLOCATED(a2F_tr) )    DEALLOCATE(a2F_tr)
    IF ( ALLOCATED(rho) )    DEALLOCATE(rho)
    IF ( ALLOCATED(dosph) )     DEALLOCATE(dosph)
    IF ( ALLOCATED(logavg) )    DEALLOCATE(logavg)
    !
  ENDDO ! isig
  !
  ENDIF
  CALL mp_barrier(world_comm)
  !
  CALL stop_clock('a2F')
  CALL print_clock('a2F')
  !
  RETURN
  !
  END SUBROUTINE eliashberg_a2f
  !-----------------------------------------------------------------------
