!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef TESTING
MODULE martyna_tuckerman
  !
  ! ... The variables needed to the Martyna-Tuckerman method for isolated
  !     systems
  !
  USE kinds, ONLY: dp
  USE constants, ONLY : e2, pi, tpi, fpi
  USE ws_base
  !
  IMPLICIT NONE
  !
  TYPE (ws_type) :: ws
  REAL (DP) :: alpha, beta
  REAL (DP), ALLOCATABLE :: wg_corr(:)
  LOGICAL :: wg_corr_is_updated = .FALSE.
  LOGICAL :: do_comp_mt = .FALSE.
  LOGICAL :: gamma_only = .FALSE.
  integer :: gstart = 1
    !
  SAVE

  PRIVATE

  PUBLIC :: tag_wg_corr_as_obsolete, do_comp_mt, &
            wg_corr_ewald, wg_corr_loc, wg_corr_h, wg_corr_force

CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE tag_wg_corr_as_obsolete
!----------------------------------------------------------------------------
     wg_corr_is_updated = .FALSE.
  END SUBROUTINE tag_wg_corr_as_obsolete
!----------------------------------------------------------------------------
  SUBROUTINE wg_corr_h( omega, ngm, rho, v, eh_corr )
!----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ngm
  REAL(DP), INTENT(IN) :: omega
  COMPLEX(DP), INTENT(IN)  :: rho(ngm)
  COMPLEX(DP), INTENT(OUT) :: v(ngm)
  REAL(DP), INTENT(OUT) :: eh_corr

  INTEGER :: ig

  IF (.NOT.wg_corr_is_updated) CALL init_wg_corr
!
  v(:) = (0._dp,0._dp)

  eh_corr =  0._dp
  DO ig = 1,ngm
     v(ig) = e2 * wg_corr(ig) * rho(ig) 
     eh_corr = eh_corr + ABS(rho(ig))**2 * wg_corr(ig)
  END DO
  iF (gamma_only) v(gstart:ngm) = 0.5_dp * v(gstart:ngm)

  eh_corr = 0.5_dp * e2 * eh_corr * omega

  RETURN
  END SUBROUTINE wg_corr_h
!----------------------------------------------------------------------------
  SUBROUTINE wg_corr_loc( omega, ntyp, ngm, zv, strf, v )
!----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ntyp, ngm
  REAL(DP), INTENT(IN) :: omega, zv(ntyp)
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  COMPLEX(DP), INTENT(OUT) :: v(ngm)
  INTEGER :: ig

  IF (.NOT.wg_corr_is_updated) CALL init_wg_corr
!
  do ig=1,ngm
     v(ig) = - e2 * wg_corr(ig) * SUM(zv(1:ntyp)*strf(ig,1:ntyp)) / omega
  end do
  iF (gamma_only) v(gstart:ngm) = 0.5_dp * v(gstart:ngm)

  RETURN
  END SUBROUTINE wg_corr_loc
!----------------------------------------------------------------------------
  SUBROUTINE wg_corr_force( lnuclei, omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, nspin, &
                            rho, force )
!----------------------------------------------------------------------------
  USE cell_base, ONLY : tpiba
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  INTEGER, INTENT(IN) :: nat, ntyp, ityp(nat), ngm, nspin
  REAL(DP), INTENT(IN) :: omega, zv(ntyp), tau(3,nat), g(3,ngm)
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp), rho(ngm,nspin)
  LOGICAL, INTENT(IN) :: lnuclei
  ! this variable is used in wg_corr_force to select if
  ! corr should be done on rho and nuclei or only on rho
  REAL(DP), INTENT(OUT) :: force(3,nat)
  INTEGER :: ig, na
  REAL (DP) :: arg
  COMPLEX(DP), ALLOCATABLE :: v(:)
  COMPLEX(DP) :: rho_tot
  !
  IF (.NOT.wg_corr_is_updated) CALL init_wg_corr
  !
  allocate ( v(ngm) )
  do ig=1,ngm
     rho_tot = rho(ig,1)
     if(lnuclei) rho_tot = rho_tot - SUM(zv(1:ntyp)*strf(ig,1:ntyp)) / omega
     if (nspin==2) rho_tot = rho_tot + rho(ig,2)
     v(ig) = e2 * wg_corr(ig) * rho_tot
  end do
  force(:,:) = 0._dp
  do na=1,nat
     do ig=1,ngm
        arg = tpi * SUM ( g(:,ig)*tau(:, na) ) 
        force(:,na) = force(:,na) + g(:,ig) * CMPLX(SIN(arg),-COS(ARG)) * v(ig)
     end do
     force(:,na) = - force(:,na) * zv(ityp(na))  * tpiba
  end do
  deallocate ( v )
  !
  call mp_sum(  force, intra_bgrp_comm )
  !
  RETURN
  END SUBROUTINE wg_corr_force
!----------------------------------------------------------------------------
  SUBROUTINE init_wg_corr
!----------------------------------------------------------------------------
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE control_flags, ONLY : gamma_only_ => gamma_only
  USE gvect,         ONLY : ngm, gg, gstart_ => gstart, nl, nlm, ecutrho
  USE cell_base,     ONLY : at, alat, tpiba2, omega

  INTEGER :: idx0, idx, ir, i,j,k, ig, nt
  REAL(DP) :: r(3), rws, upperbound, rws2
  COMPLEX (DP), ALLOCATABLE :: aux(:)
  REAL(DP), EXTERNAL :: qe_erfc
#ifdef TESTING
  REAL(DP), ALLOCATABLE :: plot(:)
  CHARACTER (LEN=25) :: filplot
  LOGICAL, SAVE :: first = .TRUE.
#endif

  IF ( ALLOCATED(wg_corr) ) DEALLOCATE(wg_corr)
  ALLOCATE(wg_corr(ngm))
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  !
  alpha = 2.9d0
  upperbound = 1._dp
  DO WHILE ( upperbound > 1.e-7_dp) 
     alpha = alpha - 0.1_dp  
     if (alpha<=0._dp) call errore('init_wg_corr','optimal alpha not found',1)
     upperbound = e2 * sqrt (2.d0 * alpha / tpi) * &
                       qe_erfc ( sqrt ( ecutrho / 4.d0 / alpha) )
  END DO
  beta = 0.5_dp/alpha ! 1._dp/alpha
  ! write (*,*) " alpha, beta MT = ", alpha, beta
  !
  call ws_init(at,ws)
  !
  gstart = gstart_
  gamma_only = gamma_only_
  !
  ! idx0 = starting index of real-space FFT arrays for this processor
  !
  idx0 = dfftp%nr1x*dfftp%nr2x * dfftp%ipp(me_bgrp+1)
  !
  ALLOCATE (aux(dfftp%nnr))
  aux = CMPLX(0._dp,0._dp)
  DO ir = 1, dfftp%nr1x*dfftp%nr2x * dfftp%npl
     !
     ! ... three dimensional indices
     !
     idx = idx0 + ir - 1
     k   = idx / (dfftp%nr1x*dfftp%nr2x)
     idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
     j   = idx / dfftp%nr1x
     idx = idx - dfftp%nr1x*j
     i   = idx

     r(:) = ( at(:,1)/dfftp%nr1*i + at(:,2)/dfftp%nr2*j + at(:,3)/dfftp%nr3*k )

     rws = ws_dist(r,ws)
#ifdef TESTING
     rws2 = ws_dist_stupid(r,ws)
     if (abs (rws-rws2) > 1.e-5 ) then
        write (*,'(4i8)') ir, i,j,k
        write (*,'(5f14.8)') r(:), rws, rws2
        stop
     end if
#endif

     aux(ir) = smooth_coulomb_r( rws*alat )

  END DO

  CALL fwfft ('Dense', aux, dfftp)

  do ig =1, ngm
     wg_corr(ig) = omega * REAL(aux(nl(ig))) - smooth_coulomb_g( tpiba2*gg(ig))
  end do
  wg_corr(:) =  wg_corr(:) * exp(-tpiba2*gg(:)*beta/4._dp)**2
  !
  if (gamma_only) wg_corr(gstart:ngm) = 2.d0 * wg_corr(gstart:ngm)
!
  wg_corr_is_updated = .true.
 
#ifdef TESTING
  if (first) then
     ALLOCATE(plot(dfftp%nnr))

     filplot = 'wg_corr_r'
     CALL invfft ('Dense', aux, dfftp)
     plot(:) = REAL(aux(:))
     call  write_wg_on_file(filplot, plot)

     filplot = 'wg_corr_g'
     aux(:) = CMPLX(0._dp,0._dp)
     do ig =1, ngm
        aux(nl(ig))  = smooth_coulomb_g( tpiba2*gg(ig))/omega
     end do
     if (gamma_only) aux(nlm(1:ngm)) = CONJG( aux(nl(1:ngm)) )

     CALL invfft ('Dense', aux, dfftp)
     plot(:) = REAL(aux(:))
     call  write_wg_on_file(filplot, plot)

     filplot = 'wg_corr_diff'
     aux(:) = CMPLX(0._dp,0._dp)
     aux(nl(1:ngm)) = wg_corr(1:ngm) / omega
     if (gamma_only) then
        aux(:) = 0.5_dp * aux(:) 
        aux(nlm(1:ngm)) = aux(nlm(1:ngm)) + CONJG( aux(nl(1:ngm)) )
     end if
     CALL invfft ('Dense', aux, dfftp)
     plot(:) = REAL(aux(:))
     call  write_wg_on_file(filplot, plot)

     DEALLOCATE (plot)

     first = .false.
  end if
#endif

  DEALLOCATE (aux)

  RETURN

  END SUBROUTINE init_wg_corr 
!----------------------------------------------------------------------------
  SUBROUTINE write_wg_on_file(filplot, plot)
!----------------------------------------------------------------------------
  USE fft_base,        ONLY : dfftp
  USE gvect,           ONLY : gcutm
  USE gvecw,           ONLY : ecutwfc
  USE gvecs,           ONLY : dual
  USE cell_base,       ONLY : at, alat, tpiba2, omega, ibrav, celldm
  USE ions_base,       ONLY : zv, ntyp => nsp, nat, ityp, atm, tau
  CHARACTER (LEN=25), INTENT(IN) :: filplot
  REAL(DP) :: plot(dfftp%nnr)
  CHARACTER (LEN=25) :: title
  INTEGER :: plot_num=0, iflag=+1

  CALL plot_io (filplot, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
     dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
     gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, plot, iflag)
  RETURN
  END SUBROUTINE write_wg_on_file
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION wg_corr_ewald ( omega, ntyp, ngm, zv, strf )
!----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ntyp, ngm
  REAL(DP), INTENT(IN) :: omega, zv(ntyp)
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  INTEGER :: ig
  COMPLEX(DP)  :: rhoion

  IF (.NOT.wg_corr_is_updated) CALL init_wg_corr
!
  wg_corr_ewald = 0._dp
  DO ig=1,ngm
     rhoion = SUM (zv(1:ntyp)* strf(ig,1:ntyp) ) / omega
     wg_corr_ewald = wg_corr_ewald + ABS(rhoion)**2 * wg_corr(ig) 
  END DO
  wg_corr_ewald = 0.5_dp * e2 * wg_corr_ewald * omega
!  write(*,*) "ewald correction   = ", wg_corr_ewald

  END FUNCTION wg_corr_ewald
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION smooth_coulomb_r(r)
!----------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: r
  REAL(DP), EXTERNAL :: qe_erf
!  smooth_coulomb_r = sqrt(2._dp*alpha/tpi)**3 * exp(-alpha*r*r) ! to be modified
  IF (r>1.e-6_dp) THEN
     smooth_coulomb_r = qe_erf(sqrt(alpha)*r)/r
  ELSE
     smooth_coulomb_r = 2._dp/sqrt(pi) * sqrt(alpha)
  END IF

  END FUNCTION smooth_coulomb_r
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION smooth_coulomb_g(q2)
!----------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: q2
!  smooth_coulomb_g = exp(-q2/4._dp/alpha) ! to be modified
  IF (q2>1.e-6_dp) THEN
     smooth_coulomb_g = fpi * exp(-q2/4._dp/alpha)/q2 ! to be modified
  ELSE 
     smooth_coulomb_g = - 1._dp * fpi * (1._dp/4._dp/alpha + 2._dp*beta/4._dp)
  END IF
  END FUNCTION smooth_coulomb_g
!----------------------------------------------------------------------------

END MODULE martyna_tuckerman
