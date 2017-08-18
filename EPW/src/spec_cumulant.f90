  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spec_cumulant 
  !-----------------------------------------------------------------------
  !!
  !!  Compute the electron spectral function including the electron-
  !!  phonon interaction using the cumulant expansion method.
  !!  
  !!  Use the frequency-dependent self-energy from the one-shot Migdal approximation.
  !!  See e.g. PRL 77, 2268 (1996) and NCOMMS 8, 15769 (2017) (and S.I.)
  !! 
  !!  CV: implementation based on the original python script from Fabio Caruso (FC)
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP, i4b
  USE constants_epw, ONLY : pi, two, zero, ryd2ev
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iospectral_sup, iospectral_cum
  USE io_files,      ONLY : prefix, tmp_dir, diropn
  USE epwcom,        ONLY : degaussw, eptemp, wmin_specfun, wmax_specfun, nw_specfun, &
                            bnd_cum
  USE elph2,         ONLY : ibndmin, ibndmax
  !
  implicit none
  !
  CHARACTER(len=64) :: line
  !! Auxiliary string
  CHARACTER(len=64) :: filespec
  !! Spectral function file
  INTEGER (kind=i4b) :: maxrecs=1.0E9
  !! Maximum number of lines in the specfun_sup.elself 
  INTEGER :: ios 
  !! Status of opening file
  INTEGER :: im
  !! Dummy counter when reading
  INTEGER :: i1, i2
  !! Auxiliary indices
  INTEGER :: ik
  !! K-point counter
  INTEGER :: iw
  !! Frequency counter
  INTEGER :: iw2
  !! Second freq. counter for the convolution
  INTEGER :: indw
  !! Auxillary index for convolution
  INTEGER :: ibnd
  !! Band index
  INTEGER :: isat
  !! Index to number the satellite
  INTEGER :: nk
  !! Total number of k-points 
  INTEGER :: iqp
  !! Energy index of renormalized quasiparticle
  INTEGER :: iks
  !! Energy index for the KS quasiparticle
  INTEGER :: i0
  !! Energy index of Fermi level
  ! 
  REAL (kind=DP) :: dw
  !! Freq. increment
  REAL (kind=DP) :: tol_e 
  !! Do the cumulant only for states below the Fermi level (or below tol_e)
  REAL (kind=DP) :: a1, a2, a3, a4 
  !! Auxiliary variables
  REAL (kind=DP) :: eptempc
  !! Temperature in eV
  REAL (kind=DP) :: eta
  !! Small broadening for denominator
  REAL (kind=DP) :: eqp
  !! Renpormalized quasiparticle energy
  REAL (kind=DP) :: ekk
  !! K-S energy
  REAL (kind=DP) :: si_ks
  !! ImSigma at KS energy
  REAL (kind=DP) :: si_qp
  !! ImSigma at renormalized energy
  REAL (kind=DP) :: deriv
  !! Derivative of ImSigma at KS energy
  REAL (kind=DP) :: den
  !! Denominator for the part for the satelite
  REAL (kind=DP), ALLOCATABLE :: ww(:)
  !! Frequency variable
  REAL (kind=DP), ALLOCATABLE :: ek(:)
  !! Eigenvalues for the band bnd_cum
  REAL (kind=DP), ALLOCATABLE :: sigmar(:,:) ! kpt, omega
  !! Real self-energy for bnd_cum 
  REAL (kind=DP), ALLOCATABLE :: sigmai(:,:)
  !! Imaginary self-energy for bnd_cum
  REAL (kind=DP), ALLOCATABLE :: a_mig(:,:)
  !! Migdal spectral function (same as in spectral_func.f90 )
  REAL (kind=DP), ALLOCATABLE :: a_qp(:)
  !! Quasiparticle contribution to the spectral function 
  REAL (kind=DP), ALLOCATABLE :: a_s(:)
  !! Temporary quantity needed to compute the satelite 
  REAL (kind=DP), ALLOCATABLE :: a_tmp(:)
  !! Stores the convolution
  REAL (kind=DP), ALLOCATABLE :: conv(:)
  !! Temporary quantity to compute the convolution 
  REAL (kind=DP), ALLOCATABLE :: a_s1(:)
  !! 1st satellite 
  REAL (kind=DP), ALLOCATABLE :: a_s2(:)
  !! 2nd satellite 
  !
  tol_e = 0.d0 ! 0.01 eV above Fermi level
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1) * ryd2ev
  eptempc = eptemp * ryd2ev
  degaussw = degaussw * ryd2ev
  degaussw = 0.005 !as FC
  eta = 0.0001 !as FC (for denominator)
  !
  WRITE(stdout,'(/5x,a)') repeat('=',100)
  WRITE(stdout,'(5x,a)') 'Performing the CUMULANT EXPANSION to obtain the electron spectral function'
  WRITE(stdout,'(5x,a)') 'There is no spin degeneracy factor, as in spectral_func.f90'
  WRITE(stdout,'(5x,a)') 'Warning: the time-ordered cumulant expansion is well-defined only for states below the Fermi energy'
  WRITE(stdout,'(5x,a)') 'Warning: the routine is sequential but very fast.'
  WRITE(stdout,'(5x,a/)') repeat('=',100)
  !
  OPEN (unit=iospectral_sup, file='specfun_sup.elself', status='old', iostat=ios)
  IF (ios /= 0) CALL errore ('spec_cumulant', 'opening file specfun_sup.elself', abs(ios) )
  !
  ! determine number of k points, ibndmin, ibndmax
  DO im=1,6
    READ (iospectral_sup, '(a)') line
  ENDDO
  DO im=1,maxrecs
    READ (iospectral_sup,*,iostat=ios) i1, i2
    IF (im.eq.1) ibndmin = i2
    IF (ios /= 0) EXIT
    IF (im.eq.maxrecs) CALL errore ('spec_cumulant', 'increase maxrecs', 1)
  ENDDO
  !
  REWIND (iospectral_sup)
  !
  nk = i1
  ibndmax = i2
  WRITE(stdout,'(5x,a,i4,a,i4,a,i4,a,f12.6/)') "Check: nk = ", nk, &
         ", ibndmin = ", ibndmin, ", ibndmax = ", ibndmax, " kbT (eV) = ", eptempc
  !eptempc=0.0025 ! as in FC script
  !
  ALLOCATE ( ww(nw_specfun), ek(nk), sigmar(nk,nw_specfun), sigmai(nk,nw_specfun), &
             a_mig(nw_specfun,nk), a_qp(nw_specfun), a_s(nw_specfun), &
             a_tmp(nw_specfun), conv(nw_specfun), a_s1(nw_specfun), a_s2(nw_specfun) )
  !
  ! read and store unperturbed energy, energy grid, real and im sigma for designated band
  DO im=1,6
    READ (iospectral_sup, '(a)') line
  ENDDO
  DO ibnd = 1, ibndmax-ibndmin+1
    DO ik = 1, nk
      DO iw = 1, nw_specfun
        READ (iospectral_sup,*) i1, i2, a1, a2, a3, a4
        IF (i2==bnd_cum) THEN
          ek(ik)=a1
          ww(iw)=a2
          ! meV->eV and multiply by Fermi factor
          sigmar(ik,iw)=a3*0.001d0 ! / ( exp( ww(iw)/eptempc )+1.d0 )
          sigmai(ik,iw)=a4*0.001d0 ! / ( exp( ww(iw)/eptempc )+1.d0 )
          ! spec func as in spectral_func.f90
          a_mig(iw,ik) = abs( sigmai(ik,iw) ) / pi / &
               ( ( ww(iw) - ek(ik) - sigmar(ik,iw) )**two + (sigmai(ik,iw) )**two )
        ENDIF
        ! ek (eV), ww (eV), sigmar (eV), sigmai (eV)
      ENDDO
    ENDDO
  ENDDO
  ! 
  CLOSE (iospectral_sup)
  !
  ! open file for cumulant spectral function
  IF (bnd_cum.lt.10) THEN
    WRITE(filespec,'(a,i1,a)') 'specfun_cum',bnd_cum,'.elself'
  ELSE IF (bnd_cum.gt.9 .and. bnd_cum.lt.100) THEN
    WRITE(filespec,'(a,i2,a)') 'specfun_cum',bnd_cum,'.elself'
  ELSE
    WRITE(filespec,'(a,i3,a)') 'specfun_cum',bnd_cum,'.elself'
  ENDIF
  OPEN (unit=iospectral_cum,file=filespec)
  !
  ! uncomment if needed
  !OPEN (unit=76,file='specfun_k0.dat')
  !
  ! define index corresponding to omega=0
  DO iw = 1, nw_specfun
    IF ( abs(ww(iw)-0.d0).le.(dw/1.8d0) ) THEN
      i0 = iw
      EXIT
    ENDIF
  ENDDO
  ! 
  DO ik = 1, nk
    !
    IF ( ek(ik) .lt. tol_e ) THEN
      !
      ekk = ek(ik)
      ! define index and energy of renormalized qp
      iqp = MAXLOC( a_mig(:,ik), dim=1 )
      eqp = ww(iqp)
      !
      ! define index corresponding to unrenormalized energy
      DO iw = 1, nw_specfun
        IF ( abs(ww(iw)-ekk).le.(dw/1.8d0) ) THEN
          iks = iw
          EXIT
        ENDIF 
      ENDDO
      !
      si_qp = sigmai(ik,iqp)
      si_ks = sigmai(ik,iks)
      !deriv = ( abs(si_qp) - abs(sigmai(ik,iqp-5)) ) / (5.d0*dw) ! as in FC script
      deriv = ( abs(si_ks) - abs(sigmai(ik,iks-5)) ) / (5.d0*dw)
      !
      ! subtract broadening from EPW calculation
      !IF (si_qp.gt.0.01) si_qp = si_qp - degaussw
      !
      ! calculate Aqp and As1
      DO iw = 1, nw_specfun
        !
        !IF ( abs(ww(iw)-ekk) .lt. 0.01 ) THEN ! FC: this is necessary to avoid divergences at the Fermi energy
        !   den = 0.01**two + eta
        !ELSE 
           den = (ww(iw)-ekk)**two + eta
        !ENDIF
        !
        a_qp(iw) = abs(si_qp) / pi / ( (ww(iw)-eqp)**two + (si_qp)**two )
        a_s (iw) = ( abs(sigmai(ik,iw)) - abs(si_ks) - (ww(iw) - ekk) * deriv ) / den / pi 
        conv(iw) = a_qp(iw)
        !
      ENDDO
      !
      DO isat = 1, 2
        !
        a_tmp = zero
        !
        DO iw = 1, nw_specfun
          !
          DO iw2 = 1, nw_specfun
            !
            IF ( ww(iw2) .lt. ekk ) THEN
              indw = i0+iw-iw2
              IF ( indw.le.nw_specfun .and. indw.gt.0 ) THEN
                 a_tmp(iw) = a_tmp(iw) + abs( a_s(iw2) * conv(indw) ) * dw
              ENDIF
            ENDIF
            !
          ENDDO
          !
        ENDDO
        !
        IF (isat.eq.1) a_s1 = a_tmp
        IF (isat.eq.2) a_s2 = a_tmp / 2.d0
        conv = a_tmp
        WRITE(stdout,'(a,i4,a,i4)') ' Done ik ', ik, ' sat ', isat
        !
      ENDDO ! isat
      !
      DO iw = 1, nw_specfun
        !
        a_tmp(iw) = a_qp(iw) + a_s1(iw) + a_s2(iw)
        !
        WRITE(iospectral_cum,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ww(iw), a_tmp(iw) / ( exp( ww(iw)/eptempc )+1.d0 )
        !
        ! uncomment and change if needed
        !IF (ik==114 .or. ik==185) WRITE(76,'(2x,f10.5,2x,f10.5)') ww(iw), a_tmp(iw) ! / ( exp( ww(iw)/eptempc )+1.d0 )
        !
      ENDDO
      !
      WRITE(iospectral_cum,'(a)') ' '
      ! uncomment and change if needed
      !IF (ik==114 .or. ik==185) WRITE (76,'(a)') ' '
      !IF (ik==114 .or. ik==185) WRITE (76,'(a)') ' '
      !
    ENDIF ! only states below Fermi level
    !
  ENDDO ! main loop k
  !
  CLOSE (iospectral_cum)
  !CLOSE (76)
  !
  DEALLOCATE ( ww, ek, sigmar, sigmai )
  DEALLOCATE ( a_mig, a_qp, a_s, a_tmp, conv, a_s1, a_s2 )
  !
  END SUBROUTINE spec_cumulant
  !
