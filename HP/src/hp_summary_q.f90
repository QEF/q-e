!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_summary_q
  !-----------------------------------------------------------------------
  !
  ! This routine writes on output the quantities which have been read
  ! from the punch file, and the quantities computed in hp_setup_q.
  ! If iverbosity = 1 only a partial summary is done.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : at
  USE klist,         ONLY : lgauss, smearing, degauss, nkstot, xk, wk
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : gcutm, ngm
  USE gvecs,         ONLY : doublegrid, dual, gcutms, ngms
  USE gvecw,         ONLY : ecutwfc
  USE fft_base,      ONLY : dffts
  USE symm_base,     ONLY : s, sr, ftau, sname
  USE funct,         ONLY : write_dft_name
  USE control_flags, ONLY : iverbosity
  USE lr_symm_base,  ONLY : irotmq, minus_q, nsymq
  USE ldaU_hp,       ONLY : conv_thr_chi

  IMPLICIT NONE
  !
  INTEGER :: i, ipol, apol, na, nt, isymq, isym, ik, nsymtot
  ! generic counter
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on atomic types
  ! counter on symmetries
  ! counter on symmetries
  ! counter on k points
  !
  REAL(DP) :: ft1, ft2, ft3, xkg(3)
  ! fractionary translations
  ! k point in crystal coordinates
  !
  WRITE( stdout, '(/,19x,"WRITING LINEAR-RESPONSE SUMMARY:",/)')
  !
  ! Now print the information specific for every q point
  !
  ! Description of symmetries for a given q point
  !
  IF (nsymq.le.1.and..not.minus_q) THEN
     WRITE( stdout, '(5x,"No symmetry (except the identity)!")')
  ELSE
     WRITE( stdout, '(/5x,"Number of symmetries in the small group of q, nsymq = ",i2)') nsymq
     IF (minus_q) WRITE( stdout, '(5x," + the symmetry q -> -q+G ")')
  ENDIF
  !
  ! Description of the symmetry matrices (and vectors of fractional 
  ! translations if f/=0) of the small group of q
  !
  IF (iverbosity > 1) THEN
     !
     WRITE( stdout, '(/5x,"Symmetry matrices (and vectors of fractional translations if f/=0):")')
     !
     IF (minus_q) THEN
        nsymtot = nsymq + 1
     ELSE
        nsymtot = nsymq
     ENDIF
     !
     DO isymq = 1, nsymtot
        !
        IF (isymq.GT.nsymq) THEN
           isym = irotmq
           WRITE( stdout, '(/5x,"This transformation sends q -> -q+G")')
        ELSE
           isym = isymq
        ENDIF
        !
        WRITE( stdout, '(/5x,"isym = ",i2,5x,a45/)') isymq, sname (isym)
        !
        IF (ftau(1,isym).NE.0 .OR. ftau(2,isym).NE.0 .OR. ftau(3,isym).NE.0) THEN
           !
           ft1 = at (1, 1) * ftau (1, isym) / dfftp%nr1 + &
                 at (1, 2) * ftau (2, isym) / dfftp%nr2 + &
                 at (1, 3) * ftau (3, isym) / dfftp%nr3
           ft2 = at (2, 1) * ftau (1, isym) / dfftp%nr1 + &
                 at (2, 2) * ftau (2, isym) / dfftp%nr2 + &
                 at (2, 3) * ftau (3, isym) / dfftp%nr3
           ft3 = at (3, 1) * ftau (1, isym) / dfftp%nr1 + &
                 at (3, 2) * ftau (2, isym) / dfftp%nr2 + &
                 at (3, 3) * ftau (3, isym) / dfftp%nr3
           !
           WRITE( stdout, '(5x,"cryst.",3x,"s(",i2,") = (",3(i6,5x)," )    f =( ",f10.7," )")') &
                & isymq, (s(1,ipol,isym), ipol=1,3), DBLE(ftau(1,isym))  / DBLE(dfftp%nr1)
           WRITE( stdout, '(21x," (",3(i6,5x), " )       ( ",f10.7," )")')  &
                &        (s(2,ipol,isym), ipol=1,3), DBLE(ftau(2,isym))  / DBLE(dfftp%nr2)
           WRITE( stdout, '(21x," (",3(i6,5x)," )       ( ",f10.7," )"/)')  &
                &        (s(3,ipol,isym), ipol=1,3), DBLE(ftau(3,isym))  / DBLE(dfftp%nr3)
           WRITE( stdout, '(5x,"cart.",4x,"s(",i2,") = (",3f11.7, " )    f =( ",f10.7," )")') &
                & isymq, (sr(1,ipol,isym), ipol=1,3), ft1
           WRITE( stdout, '(21x," (",3f11.7, " )       ( ",f10.7," )")')    &
                &        (sr(2,ipol,isym), ipol=1,3), ft2
           WRITE( stdout, '(21x," (",3f11.7, " )       ( ",f10.7," )"/)')   &
                &        (sr(3,ipol,isym), ipol=1,3), ft3
           !
        ELSE
           !
           WRITE( stdout, '(5x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                & isymq,  (s(1,ipol,isym), ipol=1,3)
           WRITE( stdout, '(21x," (",3(i6,5x)," )")')  &
                &         (s (2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(21x," (",3(i6,5x)," )"/)') &
                &         (s(3,ipol,isym), ipol=1,3)
           WRITE( stdout, '(5x,"cart.",4x,"s(",i2,") = (",3f11.7, " )")')    &
               & isymq,   (sr(1,ipol,isym), ipol=1,3)
           WRITE( stdout, '(21x," (",3f11.7," )")')    &
               &          (sr(2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(21x," (",3f11.7," )"/)')   &
               &          (sr(3,ipol,isym), ipol=1,3)
           !
        ENDIF
        !
     ENDDO ! isymq
     !
  ENDIF
  !
  ! Description of the G cutoff and the FFT grid
  !
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (",i7," G-vectors)", &
                 & "     FFT grid: (",i3,",",i3,",",i3,")")')     &
                 & gcutm, ngm, dfftp%nr1, dfftp%nr2, dfftp%nr3
  !
  IF (doublegrid) &
  WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", i7," G-vectors)", &
                & "  smooth grid: (",i3, ",",i3,",",i3,")")')     &
                & gcutms, ngms, dffts%nr1, dffts%nr2, dffts%nr3
  !
  ! Number of k (and k+q if q/=0) points
  !
  IF (.NOT.lgauss) THEN
     WRITE( stdout, '(/5x,"Number of k (and k+q if q/=0) points =",i6,/)') nkstot
  ELSE
     WRITE( stdout, '(/5x,"Number of k (and k+q if q/=0) points =", i6, 2x, &
                   & a," smearing, width (Ry) =",f8.4,/)') &
                   & nkstot, TRIM(smearing), degauss
  ENDIF
  !
  ! Coordinates of the k (and k+q if q/=0) points
  !
  IF ( iverbosity > 1 .OR. (nkstot<100) ) THEN
     !
     ! cartesian coordinates
     !
     WRITE( stdout, '(23x,"cart. coord. (in units 2pi/alat)")')
     DO ik = 1, nkstot
        WRITE( stdout, '(8x,"k (",i5,") = (",3f12.7,"), wk =",f10.7)') &
             & ik, (xk(ipol,ik), ipol=1,3), wk(ik)
     ENDDO
     !
     ! crystal coordinates
     !
     WRITE( stdout, '(/23x,"cryst. coord.")')
     DO ik = 1, nkstot
        DO ipol = 1, 3
           xkg (ipol) = at (1, ipol) * xk (1, ik) + &
                        at (2, ipol) * xk (2, ik) + &
                        at (3, ipol) * xk (3, ik)
        ENDDO
        WRITE( stdout, '(8x,"k (",i5,") = (",3f12.7,"), wk =",f10.7)') &
             & ik, (xkg(ipol), ipol=1,3), wk(ik)
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE hp_summary_q
