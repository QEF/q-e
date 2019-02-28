  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_cumulant 
  !-----------------------------------------------------------------------
  !!
  !!  Compute the electron spectral function including the electron-
  !!  phonon interaction using the (retarded) cumulant expansion method.
  !!  Use the frequency-dependent self-energy from the one-shot Migdal approximation 
  !!  read from specfun_sup.elself.
  !!
  !!  See e.g. PRL 77, 2268 (1996) and PRB 90, 085112 (2014):
  !!  the cumulant spectral function can be calculated using a series of convolutions 
  !!  in frequency space, or using FFT.
  !!  To converge the FFT with a small dt, we use zero padding of ImSigma outside the 
  !!  calculated frequency range. If the convergence is not satisfactory, one can 
  !!  set 'fact' to a larger value (6 is used as default, see the subroutine cumulant_time). 
  !! 
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP, i4b
  USE constants_epw, ONLY : pi, two, zero, ryd2ev, ryd2mev, ci
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iospectral_sup, iospectral_cum
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
  INTEGER :: ibnd
  !! Band index
  INTEGER :: nk
  !! Total number of k-points 
  INTEGER :: i0
  !! Energy index of Fermi level (w=0)
  ! 
  REAL (kind=DP) :: dw
  !! Freq. increment
  REAL (kind=DP) :: e_thresh
  !! Do the cumulant only for states with energy below e_thresh
  REAL (kind=DP) :: a1, a2, a3, a4 
  !! Auxiliary variables
  REAL (kind=DP) :: ekk
  !! K-S energy
  REAL (kind=DP) :: zeta
  !! Z factor
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
  REAL (kind=DP), ALLOCATABLE :: a_cw(:)
  !! Cumulant spectral function (convolutions)
  REAL (kind=DP), ALLOCATABLE :: a_ct(:), a_tmp(:)
  !! Cumulant spectral function (FFT) 
  !
  ! e_thresh can be changed if one needs e.g. only states below the Fermi level
  e_thresh = 10.d0 / ryd2ev ! referred to the Fermi level
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  !
  WRITE(stdout,'(/5x,a)') repeat('=',75)
  WRITE(stdout,'(5x,a)') 'Performing the CUMULANT EXPANSION of the retarded Greens function to obtain'
  WRITE(stdout,'(5x,a)') 'the electron spectral function.'
  WRITE(stdout,'(5x,a)') 'There is no spin degeneracy factor, as in spectral_func.f90'
  WRITE(stdout,'(5x,a)') 'Warning: the routine is sequential but very fast.'
  WRITE(stdout,'(5x,a/)') repeat('=',75)
  !
  OPEN (unit=iospectral_sup, file='specfun_sup.elself', status='old', iostat=ios)
  IF (ios /= 0) CALL errore ('spectral_cumulant', 'opening file specfun_sup.elself', abs(ios) )
  !
  ! determine number of k points, ibndmin, ibndmax
  DO im=1,6
    READ (iospectral_sup, '(a)') line
  ENDDO
  DO im=1,maxrecs
    READ (iospectral_sup,*,iostat=ios) i1, i2
    IF (im.eq.1) ibndmin = i2
    IF (ios /= 0) EXIT
    IF (im.eq.maxrecs) CALL errore ('spectral_cumulant', 'increase maxrecs', 1)
  ENDDO
  !
  REWIND (iospectral_sup)
  !
  nk = i1
  ibndmax = i2
  WRITE(stdout,'(5x,a/)') "Read self-energy from file specfun_sup.elself"
  WRITE(stdout,'(5x,a,i4,a,i4,a,i4,a,f12.6/)') "Check: nk = ", nk, &
         ", ibndmin = ", ibndmin, ", ibndmax = ", ibndmax, " kbT (eV) = ", eptemp*ryd2ev
  !
  ALLOCATE ( ww(nw_specfun), ek(nk), sigmar(nk,nw_specfun), sigmai(nk,nw_specfun), &
             a_mig(nw_specfun,nk), a_cw(nw_specfun), a_ct(nw_specfun), a_tmp(nw_specfun) )
  !
  ! read and store Kohn-Sham energy, energy grid, real and im sigma for designated band
  DO im=1,6
    READ (iospectral_sup, '(a)') line
  ENDDO
  DO ibnd = 1, ibndmax-ibndmin+1
    DO ik = 1, nk
      DO iw = 1, nw_specfun
        READ (iospectral_sup,*) i1, i2, a1, a2, a3, a4
        IF (i2==bnd_cum) THEN
          ! ek, w read in eV; Sigma read in meV
          ek(ik)=a1/ryd2ev
          ww(iw)=a2/ryd2ev
          sigmar(ik,iw)=a3/ryd2mev ! / ( exp( ww(iw)/eptemp )+1.d0 )
          sigmai(ik,iw)=a4/ryd2mev ! / ( exp( ww(iw)/eptemp )+1.d0 )
          ! spec func as in spectral_func.f90
          a_mig(iw,ik) = abs( sigmai(ik,iw) ) / pi / &
               ( ( ww(iw) - ek(ik) - sigmar(ik,iw) )**two + (sigmai(ik,iw) )**two )
        ENDIF
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
  WRITE(iospectral_cum,'(a)') '#       k   Energy [eV]            A(k,w) [meV^-1]            Z-factor    '
  WRITE(iospectral_cum,'(a)') '#                       with convolutions  |    using FFT    '
  WRITE(stdout,'(8x,a)') 'k   Energy [eV]            A(k,w) [meV^-1]            Z-factor    '
  WRITE(stdout,'(8x,a)') '                with convolutions  |    using FFT '
  !
  ! define index corresponding to omega=0 (Fermi level)
  i0 = MINLOC( abs(ww(:)), dim=1 )
  IF (abs(ww(i0)).gt.dw) CALL errore & 
     ('spectral_cumulant', 'w=0 needs to be included in [wmin:wmax]', 1 )
  !WRITE(stdout,'(5x,a,i4)') "      Check: ind(0) = ", i0
  a_cw = zero
  a_ct = zero
  a_tmp = zero
  ! 
  DO ik = 1, nk
    !
    IF ( ek(ik) .lt. e_thresh ) THEN
      !
      ekk = ek(ik)
      !
      ! cumulant from convolutions in frequency space (ImSigma>0 in EPW)
      CALL cumulant_conv( ekk, ww, sigmar(ik,:), -sigmai(ik,:), a_mig(:,ik), a_cw, zeta )
      !
      ! cumulant calculated in time domain + FFT (ImSigma>0 in EPW)
      IF (.true.) THEN
        CALL cumulant_time( ekk, ww, sigmar(ik,:), -sigmai(ik,:), a_mig(:,ik), a_tmp )
      ENDIF
      !
      DO iw = 1, nw_specfun
        !
        ! map the indices of the FFT frequency grid onto the original one
        IF ( iw.ge.i0) THEN
          a_ct(iw) = a_tmp(iw-i0+1)
        ELSE
          a_ct(iw) = a_tmp(iw+nw_specfun-i0+1)
        ENDIF
        !
        ! write cumulant spectral function on file (in meV^-1, as in spectral_func.f90)
        ! 3rd column: A_cum using convolutions; 4th column: A_cum using FFT
        IF (iw == 1) THEN
          WRITE (iospectral_cum,'(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,f8.4)') &
                 ik, ww(iw)*ryd2ev, a_cw(iw)/ryd2mev, a_ct(iw)/ryd2mev, zeta 
          WRITE (stdout,'(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,f8.4)') &
                 ik, ww(iw)*ryd2ev, a_cw(iw)/ryd2mev, a_ct(iw)/ryd2mev, zeta 
        ELSE 
          WRITE (iospectral_cum,'(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7)') &
                 ik, ww(iw)*ryd2ev, a_cw(iw)/ryd2mev, a_ct(iw)/ryd2mev !/ ( exp( ww(iw)/eptemp )+1.d0 )
          WRITE (stdout,'(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7)') &
                 ik, ww(iw)*ryd2ev, a_cw(iw)/ryd2mev, a_ct(iw)/ryd2mev !/ ( exp( ww(iw)/eptemp )+1.d0 )
          ! uncomment to multiply by Fermi occupation factor
        ENDIF
        !
      ENDDO
      !
      WRITE(iospectral_cum,'(a)') ' '
      WRITE(stdout,'(a)') ' '
      !
    ENDIF ! only states below energy e_thresh
    !
  ENDDO ! main loop k
  ! 
  WRITE(stdout,'(5x,a)') 'The file specfun_cum[BND].elself has been correctly written'
  !
  CLOSE (iospectral_cum)
  !
  DEALLOCATE ( ww, ek, sigmar, sigmai, a_mig, a_cw, a_ct, a_tmp )
  !
  END SUBROUTINE spectral_cumulant
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cumulant_conv(ek, ww, sigmar, sigmai, a_mig, a_cum, zeta )
  !-----------------------------------------------------------------------
  !!
  !! Convolution for the cumulant expansion. 
  !!
  !
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : pi, two, zero, ci, ryd2ev
  USE io_global,     ONLY : stdout
  USE epwcom,        ONLY : degaussw, wmin_specfun, wmax_specfun, nw_specfun
  !
  implicit none
  !
  REAL (kind=DP), INTENT (in) :: ek
  !! K-S energy
  REAL (kind=DP), INTENT (in) :: ww(nw_specfun)
  !! Frequency variable
  REAL (kind=DP), INTENT (in) :: sigmar(nw_specfun) 
  !! Real self-energy 
  REAL (kind=DP), INTENT (in) :: sigmai(nw_specfun)
  !! Imaginary self-energy
  REAL (kind=DP), INTENT (in) :: a_mig(nw_specfun)
  !! Migdal spectral function 
  REAL (kind=DP), INTENT (out) :: a_cum(nw_specfun)
  !! Cumulant spectral function
  REAL (kind=DP), INTENT (out) :: zeta
  !! Z factor
  !
  ! local variables
  !
  INTEGER :: iw
  !! Frequency counter
  INTEGER :: iw2
  !! Second freq. counter for the convolution
  INTEGER :: indw, ind1
  !! Auxillary indices for convolution
  INTEGER :: isat
  !! Index to number the satellite
  INTEGER :: iqp
  !! Energy index of renormalized quasiparticle
  INTEGER :: iks
  !! Energy index for the KS quasiparticle
  INTEGER :: i0
  !! Energy index of Fermi level
  ! 
  REAL (kind=DP) :: dw
  !! Freq. increment
  REAL (kind=DP) :: eqp
  !! Renpormalized quasiparticle energy
  REAL (kind=DP) :: si_ks
  !! ImSigma at KS energy
  REAL (kind=DP) :: si_qp
  !! ImSigma at renormalized energy
  REAL (kind=DP) :: diS
  !! Derivative of ImSigma at KS energy
  REAL (kind=DP) :: drS
  !! Derivative of ReSigma at qp energy
  REAL (kind=DP) :: a_qp(nw_specfun)
  !! Quasiparticle contribution to the spectral function 
  REAL (kind=DP) :: a_s(nw_specfun)
  !! Temporary quantity needed to compute the satelite 
  REAL (kind=DP) :: conv(nw_specfun)
  !! Temporary quantity to compute the convolution 
  REAL (kind=DP) :: a_s1(nw_specfun), a_s2(nw_specfun), a_s3(nw_specfun)
  !! satellites
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  i0 = MINLOC( abs(ww(:)), dim=1 )
  ! define index and energy of renormalized qp (note: qp energy needs to be the absolute max of a_mig)
  iqp = MAXLOC( a_mig(:), dim=1 )
  eqp = ww(iqp)
  !
  ! define index corresponding to unrenormalized (Kohn-Sham) energy
  iks = MINLOC( abs(ww(:)-ek), dim=1 )
  !WRITE(stdout,'(5x,a,i4,a,i4)') "      Check: ind(eks) = ", iks, " ind(eqp) = ", iqp
  !
  si_qp = abs(sigmai(iqp))
  si_ks = abs(sigmai(iks))
  ! finite difference derivatives of Im and Re Sigma
  diS = ( abs(sigmai(iks+1)) - abs(sigmai(iks-1)) ) / (2.d0*dw) 
  drS = ( sigmar(iqp+1) - sigmar(iqp-1) ) / (2.d0*dw)
  zeta = exp( drS )
  !
  ! calculate Aqp and As1
  DO iw = 1, nw_specfun
    !
    a_qp(iw) = zeta * si_qp / pi / ( (ww(iw)-eqp)**two + (si_qp)**two )
    conv(iw) = a_qp(iw)
    !
    ind1=iks+iw-i0
    IF (ind1>0 .and. ind1<nw_specfun) THEN !RC
      a_s (iw) = ( abs(sigmai(ind1)) - si_ks - (ww(iw)) * diS ) * real (1.d0 / (ww(iw)-ci*degaussw)**2.d0 ) / pi
    ELSE
      a_s (iw) = ( - si_ks - (ww(iw)) * diS ) * real (1.d0 / (ww(iw)-ci*degaussw)**2.d0 ) / pi
    ENDIF
    !
  ENDDO
  !
  DO isat = 1, 3
    !
    a_cum = zero
    !
    DO iw = 1, nw_specfun
      !
      DO iw2 = 1, nw_specfun
        !
        indw = i0+iw-iw2
        IF ( indw.le.nw_specfun .and. indw.gt.0 ) THEN
           a_cum(iw) = a_cum(iw) + abs( a_s(iw2) * conv(indw) ) * dw
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    IF (isat.eq.1) a_s1 = a_cum
    IF (isat.eq.2) a_s2 = a_cum / 2.d0
    IF (isat.eq.3) a_s3 = a_cum / 6.d0
    conv = a_cum
    !
  ENDDO ! isat
  !
  DO iw = 1, nw_specfun
    !
    a_cum(iw) = a_qp(iw) + a_s1(iw) + a_s2(iw) + a_s3(iw)
    !
  ENDDO
  !
  END SUBROUTINE cumulant_conv
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cumulant_time( ek, ww, sigmar, sigmai, a_mig, a_cum )
  !-----------------------------------------------------------------------
  !!
  !! Cumulant calculated in the time domain
  !!
  !
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : pi, two, zero, czero, ci, ryd2ev
  USE io_global,     ONLY : stdout
  USE fft_scalar,    ONLY : cfft3d
  USE epwcom,        ONLY : degaussw, wmin_specfun, wmax_specfun, nw_specfun
  !
  implicit none
  !
  REAL (kind=DP), INTENT (in) :: ek
  !! K-S energy
  REAL (kind=DP), INTENT (in) :: ww(nw_specfun)
  !! Frequency variable
  REAL (kind=DP), INTENT (in) :: sigmar(nw_specfun) 
  !! Real self-energy 
  REAL (kind=DP), INTENT (in) :: sigmai(nw_specfun)
  !! Imaginary self-energy
  REAL (kind=DP), INTENT (in) :: a_mig(nw_specfun)
  !! Migdal spectral function 
  REAL (kind=DP), INTENT (out) :: a_cum(nw_specfun)
  !! Cumulant spectral function
  !
  ! local variables
  !
  INTEGER :: iw, it
  !! Frequency, time counters
  INTEGER :: iqp
  !! Energy index of renormalized quasiparticle
  INTEGER :: iks
  !! Energy index for the KS quasiparticle
  INTEGER :: i0
  !! Energy index of Fermi level
  INTEGER :: ind1, ind2
  !! aux indices
  INTEGER :: nw_new
  !! number of points used for FFT 
  !! Zero padding is used for ImSigma outside the energy range computed.
  ! 
  REAL (kind=DP) :: fact
  !! factor used to increase the number of points for the FFT;
  !! fact=4 gives good convergence but it can be increased if needed.
  REAL (kind=DP) :: dw
  !! Freq. increment
  REAL (kind=DP) :: eqp
  !! Renormalized quasiparticle energy
  REAL (kind=DP) :: drS, zeta
  !! Derivative of ReSigma at qp energy, Z factor
  REAL (kind=DP) :: smeart
  !! small broadening
  REAL (kind=DP) :: dt
  !! time increment
  REAL (kind=DP) :: tmin, tmax
  !! min and max of time interval
  REAL (kind=DP) :: tt 
  !! time variable
  COMPLEX (kind=DP), ALLOCATABLE :: cumS(:)
  !! satellite part of the cumulant function
  COMPLEX (kind=DP), ALLOCATABLE :: cum(:)
  !! complete cumulant function
  COMPLEX (kind=DP) :: qpfac
  !! quasiparticle factor
  !
  fact = 6.d0
  ! fact=1 corresponds to FFT with the original no. of points nw_specfun
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  smeart = -degaussw ! correct sign for RC
  !
  tmin = 0.d0
  tmax = 2.d0*pi / dw
  dt = 2.d0*pi / ( (wmax_specfun - wmin_specfun) * fact )
  nw_new = int( fact * (nw_specfun-1) + 1 ) ! to be consistent with dt above
  !
  ALLOCATE( cumS(nw_new), cum(nw_new) )
  !
  i0 = MINLOC( abs(ww(:)), dim=1 )
  !
  ! define index and energy of renormalized qp (note: qp energy needs to be the absolute max of a_mig)
  iqp = MAXLOC( a_mig(:), dim=1 )
  eqp = ww(iqp)
  !
  ! define index corresponding to unrenormalized (Kohn-Sham) energy
  iks = MINLOC( abs(ww(:)-ek), dim=1 )
  !
  drS = ( sigmar(iqp+1) - sigmar(iqp-1) ) / (2.d0*dw)
  qpfac = exp( -ci*( ek + sigmar(iqp) ) + smeart*0.5d0 )
  zeta = exp( drS )
  !
  cumS = czero
  DO iw = 1, nw_new
    !
    ! the w shift is needed because FFT uses positive w \in [0:Omega], Omega=wmax-wmin
    IF ( iw.le.(nw_specfun-i0+1) ) THEN
      ind1 = iw+i0-1
      cumS(iw) = dw * abs(sigmai(ind1))/pi * real(1.d0 / (ek-ww(ind1)-ci*smeart)**2.d0 )
    ELSE IF ( iw.gt.(nw_new-i0+1) ) THEN
      ind2 = iw-nw_new+i0-1
      cumS(iw) = dw * abs(sigmai(ind2))/pi * real(1.d0 / (ek-ww(ind2)-ci*smeart)**2.d0 )
    ENDIF
    !
  ENDDO
  !
  CALL cfft3d ( cumS(:), nw_new,1,1, nw_new,1,1, 1, -1 )
  !this is needed because cfft3d(...,-1) carries a renomalization factor 1/nw_new
  cumS = cumS * nw_new
  !
  DO it = 1, nw_new
    !
    tt = tmin + dble(it-1)*dt
    cumS(it) = cumS(it) * exp(ci*ek*tt)
    cum(it)  = zeta * qpfac**tt * exp(cumS(it))
    !
  ENDDO
  !
  CALL cfft3d ( cum(:), nw_new,1,1, nw_new,1,1, 1, 1 )
  cum = cum *dt / pi
  !
  ! extract the spectral function a_cum on the original w FFT grid (nw_specfun points)
  DO iw = 1, nw_specfun
    IF ( iw.le.(nw_specfun-i0+1) ) THEN
      a_cum(iw) = real(real(cum(iw)))
    ELSE 
      ind1 = iw+nw_new-nw_specfun
      a_cum(iw) = real(real(cum(ind1)))
    ENDIF
  ENDDO
  !
  DEALLOCATE( cumS, cum )
  !
  END SUBROUTINE cumulant_time
