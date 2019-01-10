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
  !!  phonon interaction using the cumulant expansion method.
  !!  
  !!  Use the frequency-dependent self-energy from the one-shot Migdal approximation.
  !!  See e.g. PRL 77, 2268 (1996) and NCOMMS 8, 15769 (2017) (and S.I.):
  !!  the cumulant spectral function is calculated using a series of convolutions in frequency space.
  !! 
  !!  CV: implementation based on the original python script from Fabio Caruso (FC)
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP, i4b
  USE constants_epw, ONLY : pi, two, zero, ryd2ev, ci
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
  INTEGER :: iw2
  !! Second freq. counter for the convolution
  INTEGER :: indw, ind1
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
  REAL (kind=DP) :: e_thresh
  !! Do the cumulant only for states below the Fermi level (or below e_thresh)
  REAL (kind=DP) :: a1, a2, a3, a4 
  !! Auxiliary variables
  REAL (kind=DP) :: eptempc
  !! Temperature in eV
  REAL (kind=DP) :: eqp
  !! Renpormalized quasiparticle energy
  REAL (kind=DP) :: ekk
  !! K-S energy
  REAL (kind=DP) :: si_ks
  !! ImSigma at KS energy
  REAL (kind=DP) :: si_qp
  !! ImSigma at renormalized energy
  REAL (kind=DP) :: diS
  !! Derivative of ImSigma at KS energy
  REAL (kind=DP) :: drS, zeta
  !! Derivative of ReSigma at qp energy, Z factor
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
  REAL (kind=DP), ALLOCATABLE :: a_s1(:), a_s2(:), a_s3(:)
  !! satellites
  !
  ! e_thresh can be changed if one wants e.g. only states below the Fermi level
  e_thresh = 10.d0 ! (eV) referred to the Fermi level
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1) * ryd2ev
  eptempc = eptemp * ryd2ev
  degaussw = degaussw * ryd2ev
  iks = 0
  i0 = 0
  !
  WRITE(stdout,'(/5x,a)') repeat('=',106)
  WRITE(stdout,'(5x,a)') 'Performing the CUMULANT EXPANSION for the retarded Greens function'
  WRITE(stdout,'(5x,a)') 'to obtain the electron spectral function'
  WRITE(stdout,'(5x,a)') 'There is no spin degeneracy factor, as in spectral_func.f90'
  !WRITE(stdout,'(5x,a)') 'Warning: only the hole part of the time-ordered Greens function is calculated (states below the Fermi energy)'
  WRITE(stdout,'(5x,a)') 'Warning: the routine is sequential but very fast.'
  WRITE(stdout,'(5x,a/)') repeat('=',106)
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
  WRITE(stdout,'(5x,a,i4,a,i4,a,i4,a,f12.6/)') "Check: nk = ", nk, &
         ", ibndmin = ", ibndmin, ", ibndmax = ", ibndmax, " kbT (eV) = ", eptempc
  !
  ALLOCATE ( ww(nw_specfun), ek(nk), sigmar(nk,nw_specfun), sigmai(nk,nw_specfun), &
             a_mig(nw_specfun,nk), a_qp(nw_specfun), a_s(nw_specfun), a_tmp(nw_specfun), &
             conv(nw_specfun), a_s1(nw_specfun), a_s2(nw_specfun), a_s3(nw_specfun) )
  a_s = 0.d0
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
          ek(ik)=a1
          ww(iw)=a2
          ! Sigma in meV->eV 
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
  WRITE(iospectral_cum,'(a)') '#       k  Energy [eV]  A(k,w)[meV^-1]'
  WRITE(stdout,'(2x,a)') '      k  Energy [eV]  A(k,w)[meV^-1]'
  ! define index corresponding to omega=0 (Fermi level)
  i0 = MINLOC( abs(ww(:)), dim=1 )
  IF (abs(ww(i0)).gt.dw) CALL errore & 
     ('spectral_cumulant', 'w=0 needs to be included in [wmin:wmax]', 1 )
  !WRITE(stdout,'(5x,a,i4)') "      Check: ind(0) = ", i0
  ! 
  DO ik = 1, nk
    !
    IF ( ek(ik) .lt. e_thresh ) THEN
      !
      ekk = ek(ik)
      ! define index and energy of renormalized qp (note: qp energy needs to be the absolute max of a_mig)
      iqp = MAXLOC( a_mig(:,ik), dim=1 )
      eqp = ww(iqp)
      !
      ! define index corresponding to unrenormalized (Kohn-Sham) energy
      iks = MINLOC( abs(ww(:)-ekk), dim=1 )
      !WRITE(stdout,'(5x,a,i4,a,i4)') "      Check: ind(eks) = ", iks, " ind(eqp) = ", iqp
      !
      si_qp = sigmai(ik,iqp)
      si_ks = sigmai(ik,iks)
      ! finite difference derivatives of Im and Re Sigma
      diS = ( abs(sigmai(ik,iks+1)) - abs(sigmai(ik,iks-1)) ) / (2.d0*dw) 
      drS = ( sigmar(ik,iqp+1) - sigmar(ik,iqp-1) ) / (2.d0*dw)
      zeta = exp( drS )
      WRITE(*,'(a,f12.4)') "      Z-factor = ", zeta
      !
      ! calculate Aqp and As1
      DO iw = 1, nw_specfun
        !
        a_qp(iw) = zeta * abs(si_qp) / pi / ( (ww(iw)-eqp)**two + (si_qp)**two )
        conv(iw) = a_qp(iw)
        !
        ind1=iks+iw-i0
        IF (ind1>0 .and. ind1<nw_specfun) THEN !RC
          a_s (iw) = ( abs(sigmai(ik,ind1)) - abs(si_ks) - (ww(iw)) * diS ) * real (1.d0 / (ww(iw)-ci*degaussw)**2.d0 ) / pi
        ELSE
          a_s (iw) = ( - abs(si_ks) - (ww(iw)) * diS ) * real (1.d0 / (ww(iw)-ci*degaussw)**2.d0 ) / pi
        ENDIF
        !
      ENDDO
      !
      DO isat = 1, 3
        !
        a_tmp = zero
        !
        DO iw = 1, nw_specfun
          !
          DO iw2 = 1, nw_specfun
            !
            indw = i0+iw-iw2
            IF ( indw.le.nw_specfun .and. indw.gt.0 ) THEN
               a_tmp(iw) = a_tmp(iw) + abs( a_s(iw2) * conv(indw) ) * dw
            ENDIF
            !
          ENDDO
          !
        ENDDO
        !
        IF (isat.eq.1) a_s1 = a_tmp
        IF (isat.eq.2) a_s2 = a_tmp / 2.d0
        IF (isat.eq.3) a_s3 = a_tmp / 6.d0
        conv = a_tmp
        !WRITE(stdout,'(a,i4,a,i4)') ' Done ik ', ik, ' sat ', isat
        !
      ENDDO ! isat
      !
      DO iw = 1, nw_specfun
        !
       a_tmp(iw) = a_qp(iw) + a_s1(iw) + a_s2(iw) + a_s3(iw)
        !
        ! write cumulant spectral function on file (in meV^-1, as in spectral_func.f90)
        WRITE(iospectral_cum,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ww(iw), a_tmp(iw)*1.0e-3 !/ ( exp( ww(iw)/eptempc )+1.d0 )
        WRITE(stdout,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ww(iw), a_tmp(iw)*1.0e-3 !/ ( exp( ww(iw)/eptempc )+1.d0 )
        ! uncomment to multiply by Fermi occupation faction
        !
      ENDDO
      !
      WRITE(iospectral_cum,'(a)') ' '
      WRITE(stdout,'(a)') ' '
      !WRITE(stdout,'(5x,a,i4)') ' Done ik ', ik
      !
    ENDIF ! only states below energy e_thresh
    !
  ENDDO ! main loop k
  !
  CLOSE (iospectral_cum)
  !
  DEALLOCATE ( ww, ek, sigmar, sigmai )
  DEALLOCATE ( a_mig, a_qp, a_s, a_tmp, conv, a_s1, a_s2, a_s3 )
  !
  END SUBROUTINE spectral_cumulant
  !
