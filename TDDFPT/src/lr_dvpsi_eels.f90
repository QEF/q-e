!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_dvpsi_eels (ik, dvpsi)
  !--------------------------------------------------------------------------
  !
  ! This subroutine calculates the starting Lanczos vectors.
  ! EELS: The external perturbation is an electron (plane wave) 
  ! with a wavevector q.
  !
  ! Writtent by Iurii Timrov (2013)
  !
  USE kinds,                 ONLY : DP
  USE wvfct,                 ONLY : npwx, nbnd
  USE io_global,             ONLY : stdout
  USE fft_base,              ONLY : dffts, dfftp
  USE wvfct,                 ONLY : npw, igk, ecutwfc, g2kin
  USE realus,                ONLY : npw_k, igk_k   
  USE qpoint,                ONLY : npwq, igkq, ikks, ikqs, nksq, xq 
  USE lr_variables,          ONLY : evc0, lr_periodic, q1, q2, q3
  USE eqv,                   ONLY : evq, dpsi 
  USE wavefunctions_module,  ONLY : evc
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
  USE units_ph,              ONLY : lrwfc, iuwfc
  use klist,                 only : xk
  use gvect,                 only : ngm, g
  use cell_base,             only : tpiba2
  USE control_ph,            ONLY : nbnd_occ
  use lsda_mod,              only : current_spin
  USE io_files,              ONLY : iunigk, nwordwfc, iunwfc, prefix, diropn
  use uspp,                  only : vkb, okvan
  USE phus,                  ONLY : becp1
  USE becmod,                ONLY : calbec
  USE mp_bands,              ONLY : ntask_groups
 
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ik 
  COMPLEX(DP), INTENT(out) :: dvpsi(npwx*npol,nbnd)
  INTEGER :: ikk, ikq, ibnd, ig, ios, incr, v_siz
  COMPLEX(DP), ALLOCATABLE :: revc(:,:), &  ! wavefunctions in R-space
                            & dpot(:,:), &  ! the external potential on thick  grid
                            & dpots(:),  &  ! the external potential on smooth grid
                            & spsi(:,:), &  ! temporarry array (needed for US-PP)
                            & tg_psic(:,:)  ! wavefunctions in R-space (for task groups) 
  LOGICAL :: exst
  !
  CALL start_clock ('lr_dvpsi_eels') 
  !
  IF ( ntask_groups > 1 ) dffts%have_task_groups = .TRUE.
  !
  ALLOCATE (revc(dffts%nnr,npol))
  ALLOCATE (dpot(dfftp%nnr,nspin_mag))
  ALLOCATE (dpots(dffts%nnr))
  revc(:,:)  = (0.d0, 0.d0)
  dpot(:,:)  = (0.d0, 0.d0) 
  dpots(:)   = (0.d0, 0.d0)  
  dvpsi(:,:) = (0.d0, 0.d0)
  dpsi(:,:)  = (0.d0, 0.d0)
  !
  incr = 1
  !
  IF ( dffts%have_task_groups ) THEN
     !
     v_siz =  dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE( tg_psic(v_siz,npol) )
     !
     incr = dffts%nogrp
     !
  ENDIF 
  !
  IF (lr_periodic) THEN
     ikk = ik
     ikq = ik
  ELSE
     ikk = ikks(ik)
     ikq = ikqs(ik)
  ENDIF
  !
  ! Determination of npw, igk, and npwq, igkq;
  ! g2kin is used here as a work space.
  !
  CALL gk_sort( xk(1,ikk), ngm, g, ( ecutwfc / tpiba2 ), npw,  igk,  g2kin )  
  CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ), npwq, igkq, g2kin )
  !
  ! Read npw, igk, npwq, igkq.
  !
!  if (nksq.gt.1) then
!      read (iunigk, err = 100, iostat = ios) npw, igk
!100   call errore ('lr_apply_liouvillian', 'reading igk', abs (ios) )
!  endif
  !
!  if (nksq.gt.1) then
!      read (iunigk, err = 200, iostat = ios) npwq, igkq
!200   call errore ('lr_apply_liouvillian', 'reading igkq', abs (ios) )
!  endif
  !
  ! Read unperturbed wavefuctions psi(k) and psi(k+q)
  !
  IF (lr_periodic) THEN
     evc(:,:) = evc0(:,:,ik)
     evq(:,:) = evc0(:,:,ik)
  ELSE
    IF (nksq.GT.1) THEN
       CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
       CALL davcio (evq, lrwfc, iuwfc, ikq, - 1)
    ENDIF
  ENDIF
  !
  ! Re-ordering of the G vectors.
  ! If lr_periodic=.true. there is a multiplication with
  ! the external potential exp(i q*r) in R-space.
  !
  DO ibnd = 1, nbnd_occ(ikk), incr
     !
     IF ( dffts%have_task_groups ) THEN
        !
        ! FFT to R-space
        !
        CALL cft_wave_tg(evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
        !
        IF (lr_periodic) CALL errore('lr_dvpsi_eels', 'Task groups + &
                       & lr_periodic=.true. is not implemented', 1)
        !
        ! back-FFT to G-space
        !
        CALL cft_wave_tg(dvpsi, tg_psic, -1, v_siz, ibnd, nbnd_occ(ikk))
        !
     ELSE
        !
        ! FFT to R-space
        !
        CALL cft_wave(evc(1,ibnd), revc, +1)
        !
        ! Apply the external potential: exp(i q*r) * revc(r)
        !
        IF (lr_periodic) THEN
           !
           ! Define the external potential.
           !
           CALL defextpot_eels(dpot,q1,q2,q3)
           !
           ! Interpolate the potential from the thick grid
           ! to a smooth grid if doublegrid=.true.,
           ! otherwise just copy. 
           ! 
           CALL cinterpolate (dpot(:,1), dpots, -1) 
           !
           DEALLOCATE (dpot)
           ALLOCATE (dpot(dffts%nnr,nspin_mag))
           !
           dpot(:,1) = dpots(:)
           !
           ! Warning: the noncollinear and spin-pollarized cases are
           ! not supported when lr_periodic=.true.
           !
           ! Multiply the external potential with the unperturbed wvct's.
           !
           CALL apply_dpot(dffts%nnr, revc, dpot, current_spin)
           !
        ENDIF
        !
        ! back-FFT to G-space
        !
        CALL cft_wave(dvpsi(1,ibnd), revc, -1) 
        !
     ENDIF
     !
  ENDDO
  !
  dffts%have_task_groups = .FALSE.
  !
  ! Ortogonalize dvpsi to valence states.
  ! Apply -P_c, and then change the sign, because we need +P_c.
  !
  ! Ultrasoft case: calculate beta-functions vkb.
  !
  IF (okvan) CALL init_us_2 (npwq, igkq, xk(1,ikq), vkb)
  !
  CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq) 
  dvpsi = -dvpsi
  !
  ! In the ultrasoft case apply the S^{-1} operator.
  !
  IF (okvan) THEN
     !
     ALLOCATE (spsi(npwx*npol, nbnd))
     spsi(:,:) = (0.0d0, 0.0d0)
     !
     CALL sm1_psi(.TRUE.,ik, npwx, npwq, nbnd_occ(ikk), dvpsi, spsi)
     dvpsi(:,:) = spsi(:,:)
     !
     DEALLOCATE (spsi)
     !
  ENDIF
  !
  DEALLOCATE (revc)
  DEALLOCATE (dpot)
  DEALLOCATE (dpots)
  !
  IF ( ntask_groups > 1) dffts%have_task_groups = .TRUE.
  !
  IF ( dffts%have_task_groups ) THEN
     DEALLOCATE( tg_psic )
  ENDIF
  !
  dffts%have_task_groups = .FALSE.
  !
  CALL stop_clock ('lr_dvpsi_eels')
  !
  ! Calculation of becp1 = <vkb|evc>
  !
  IF (okvan) THEN
     CALL init_us_2(npw, igk, xk(1,ikk), vkb)
     CALL calbec (npw, vkb, evc, becp1(ik))
  ENDIF
  !
  RETURN
  !
CONTAINS

SUBROUTINE defextpot_eels (dvext, q1_, q2_, q3_)
  !
  ! This subroutine defines the external potential on the smooth mesh. 
  ! Inspired by Modules/compute_dipole.f90
  ! Potential is computed on the thick grid dfftp%nnr.
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : pi
  USE cell_base,     ONLY : at, bg
  USE fft_base,      ONLY : dfftp
  USE mp_global,     ONLY : me_bgrp
  USE spin_orb,      ONLY : domag
  
  IMPLICIT NONE

  COMPLEX(DP), INTENT(OUT) :: dvext(dfftp%nnr,nspin_mag)
  ! the external potential
  INTEGER :: i, j, k, ir, ip, nnr, ir_end, index0
  REAL(DP) :: r(3), q1_, q2_, q3_, inv_nr1, inv_nr2, inv_nr3
  !
  dvext(:,:) = (0.0d0, 0.0d0)
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
  nnr = dfftp%nnr
  !
#if defined (__MPI)
  index0 = dfftp%nr1x*dfftp%nr2x*SUM(dfftp%npp(1:me_bgrp))
  ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  index0 = 0
  ir_end = nnr
#endif
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
     ENDDO
     !
     ! ... minimum image convention
     !
     CALL cryst_to_cart( 1, r, bg, -1 )
     !
     r(:) = r(:) - ANINT( r(:) )
     !
     CALL cryst_to_cart( 1, r, at, 1 )
     !
     dvext(ir,1) =  cmplx( cos( 2.0d0*pi*(q1_*r(1) + q2_*r(2) + q3_*r(3)) ), &
                           sin( 2.0d0*pi*(q1_*r(1) + q2_*r(2) + q3_*r(3)) ), kind=dp)
     !
     ! IT: I am not sure that the lines below are correct.
     !
     IF (noncolin .AND. domag) THEN
        !
        dvext(ir,2) = dvext(ir,1)
        dvext(ir,3) = dvext(ir,1)
        dvext(ir,4) = dvext(ir,1)
        !
     ENDIF
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE defextpot_eels

END SUBROUTINE lr_dvpsi_eels
