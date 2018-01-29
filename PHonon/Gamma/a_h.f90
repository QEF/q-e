!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE A_h(npw,e,h,ah)
  !-----------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE cell_base,ONLY : alat, omega, tpiba2
  USE uspp,     ONLY : vkb, nkb
  USE lsda_mod, ONLY : current_spin, nspin
  USE wvfct, ONLY: nbnd, npwx, g2kin
  USE wavefunctions_module,  ONLY: evc, psic
  USE scf,      ONLY : vrs, rho
  USE fft_base, ONLY : dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect,    ONLY : gstart, g, gg
  USE constants,  ONLY: degspin, e2, fpi
  USE becmod, ONLY: bec_type, becp, calbec
  USE cgcom
  USE funct, ONLY: dft_is_gradient
  !
  IMPLICIT NONE
  INTEGER :: npw, j, jkb, ibnd, na,nt,ih
  real(DP) :: e(nbnd)
  COMPLEX(DP) :: h(npwx,nbnd), ah(npwx,nbnd)
  !
  COMPLEX(DP) :: fp, fm
  COMPLEX(DP), POINTER :: dpsic(:), drhoc(:)
  REAL(dp), allocatable :: dv(:)
  real(DP), POINTER  :: drho(:)
  !
  CALL start_clock('a_h')
  !
  drho  => auxr
  dpsic => aux2
  drhoc => aux3
  !
  drho(:) = 0.d0
  !
  ! [(k+G)^2 - e ]psi
  DO ibnd = 1,nbnd
     ! set to zero the imaginary part of h at G=0
     ! needed for numerical stability
     IF (gstart==2) h(1,ibnd) = cmplx( dble(h(1,ibnd)),0.d0,kind=DP)
     DO j = 1,npw
        ah(j,ibnd) = (g2kin(j)-e(ibnd)) * h(j,ibnd)
     ENDDO
  ENDDO
  !     V_Loc psi
  DO ibnd = 1,nbnd, 2
     dpsic(:)= (0.d0, 0.d0)
     psic(:) = (0.d0, 0.d0)
     IF (ibnd<nbnd) THEN
        ! two ffts at the same time
        DO j = 1,npw
           psic (dfftp%nl (j)) = evc(j,ibnd) + (0.0d0,1.d0)* evc(j,ibnd+1)
           dpsic(dfftp%nl (j)) =   h(j,ibnd) + (0.0d0,1.d0)*   h(j,ibnd+1)
           psic (dfftp%nlm(j))= conjg(evc(j,ibnd)-(0.0d0,1.d0)* evc(j,ibnd+1))
           dpsic(dfftp%nlm(j))= conjg(  h(j,ibnd)-(0.0d0,1.d0)*   h(j,ibnd+1))
        ENDDO
     ELSE
        DO j = 1,npw
           psic (dfftp%nl (j)) = evc(j,ibnd)
           dpsic(dfftp%nl (j)) =   h(j,ibnd)
           psic (dfftp%nlm(j)) = conjg( evc(j,ibnd))
           dpsic(dfftp%nlm(j)) = conjg(   h(j,ibnd))
        ENDDO
     ENDIF
     CALL invfft ('Wave', psic, dffts)
     CALL invfft ('Wave',dpsic, dffts)
     DO j = 1,dfftp%nnr
        drho(j) = drho(j) - 2.0d0*degspin/omega *   &
              dble(psic(j)*conjg(dpsic(j)))
        dpsic(j) = dpsic(j) * vrs(j,current_spin)
     ENDDO
     CALL fwfft ('Wave',dpsic, dffts)
     IF (ibnd<nbnd) THEN
        ! two ffts at the same time
        DO j = 1,npw
           fp = (dpsic (dfftp%nl(j)) + dpsic (dfftp%nlm(j)))*0.5d0
           fm = (dpsic (dfftp%nl(j)) - dpsic (dfftp%nlm(j)))*0.5d0
           ah(j,ibnd  ) = ah(j,ibnd)  +cmplx( dble(fp), aimag(fm),kind=DP)
           ah(j,ibnd+1) = ah(j,ibnd+1)+cmplx(aimag(fp),- dble(fm),kind=DP)
        ENDDO
     ELSE
        DO j = 1,npw
           ah(j,ibnd) = ah(j,ibnd)  + dpsic (dfftp%nl(j))
        ENDDO
     ENDIF
  ENDDO
  !
  ! V_NL psi
  CALL calbec ( npw, vkb, h, becp)
  IF (nkb > 0) CALL add_vuspsi (npwx, npw, nbnd, ah)
  !
  DO j = 1,dfftp%nnr
     drhoc(j) = cmplx(drho(j),0.d0,kind=DP)
  ENDDO
  CALL fwfft ('Rho', drhoc, dfftp)
  DO j = 1,dfftp%ngm
     dpsic(j) = drhoc(dfftp%nl(j))
  ENDDO
  !
  ! drho is deltarho(r)
  ! drhoc is deltarho(g) on the FFT grid
  ! dpsic is deltarho(g) on the G-vector grid
  !
  !  mu'(n(r)) psi(r) delta psi(r)
  !
  ALLOCATE (dv(dfftp%nnr))
  DO j = 1,dfftp%nnr
     dv(j) = drho(j)*dmuxc(j)
  ENDDO
  !
  !  add gradient correction contribution (if any)
  !
  CALL start_clock('dgradcorr')
  IF (dft_is_gradient() ) CALL dgradcor1  &
       (dfftp, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s,            &
        drho, dpsic, nspin, g, dv)
  CALL stop_clock('dgradcorr')
  NULLIFY(dpsic)
  NULLIFY (drho)
  !
  !  1/|r-r'| * psi(r') delta psi(r')
  !
  ! gstart is the first nonzero G vector (needed for parallel execution)
  !
  IF (gstart==2) drhoc(dfftp%nl(1)) = 0.d0
  !
  DO j = gstart,dfftp%ngm
     drhoc(dfftp%nl (j)) = e2*fpi*drhoc(dfftp%nl(j))/ (tpiba2*gg(j))
     drhoc(dfftp%nlm(j)) = conjg(drhoc(dfftp%nl (j)))
  ENDDO
  CALL invfft ('Rho', drhoc, dfftp)
  !
  ! drhoc now contains deltaV_hartree
  !
  DO j = 1,dfftp%nnr
     dv(j) = - dv(j) - dble(drhoc(j))
  ENDDO
  !
  CALL vloc_psi_gamma(npwx, npw, nbnd, evc, dv, ah)
  !
  NULLIFY(drhoc)
  DEALLOCATE (dv)
  !
  ! set to zero the imaginary part of ah at G=0
  ! needed for numerical stability
  IF (gstart==2) THEN
     DO ibnd = 1, nbnd
        ah(1,ibnd) = cmplx( dble(ah(1,ibnd)),0.d0,kind=DP)
     ENDDO
  ENDIF
  !
  CALL stop_clock('a_h')
  !
  RETURN
END SUBROUTINE A_h
