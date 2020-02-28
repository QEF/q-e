!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_phq
  !-----------------------------------------------------------------------
  !
  ! Dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE klist,         ONLY : nks, nkstot
  USE wvfct,         ONLY : nbnd, npwx
  USE gvect,         ONLY : ngm
  USE lsda_mod,      ONLY : nspin
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE fft_base,      ONLY : dfftp
  USE wavefunctions, ONLY : evc
  USE spin_orb,      ONLY : lspinorb, domag
  USE nc_mag_aux,    ONLY : int1_nc_save, deeq_nc_save
  USE becmod,        ONLY : bec_type, becp, allocate_bec_type
  USE uspp,          ONLY : okvan, nkb, vkb
  USE paw_variables, ONLY : okpaw
  USE uspp_param,    ONLY : nhm
  USE ramanm,        ONLY : ramtns, lraman
  USE phus,          ONLY : int1, int1_nc, int2, int2_so, &
                            int4, int4_nc, int5, int5_so, becsumort, &
                            alphasum, alphasum_nc, becsum_nc, alphap
  USE efield_mod,    ONLY : zstareu, zstareu0, zstarue0, zstarue0_rec, zstarue
  USE units_ph,      ONLY : this_pcxpsi_is_on_file, this_dvkb3_is_on_file
  USE dynmat,        ONLY : dyn00, dyn, dyn_rec, w2
  USE modes,         ONLY : u, npert, name_rap_mode, num_rap_mode
  USE el_phon,       ONLY : el_ph_mat, elph
  USE freq_ph,       ONLY : polar, nfs
  USE lrus,          ONLY : becp1, dpqq, dpqq_so
  USE qpoint,        ONLY : nksq, eigqts, xk_col
  USE eqv,           ONLY : dpsi, evq, vlocq, dmuxc, dvpsi
  USE lr_symm_base,  ONLY : rtau
  USE control_lr,    ONLY : lgamma, ofsbeta
  USE ldaU,          ONLY : lda_plus_u, Hubbard_lmax, nwfcU
  USE ldaU_ph,       ONLY : dnsbare, dnsorth, dnsbare_all_modes, wfcatomk, &
                            dwfcatomk, sdwfcatomk, wfcatomkpq, dwfcatomkpq,  &
                            swfcatomk, swfcatomkpq, sdwfcatomkpq, dvkb, vkbkpq, &
                            dvkbkpq
  USE qpoint_aux,    ONLY : becpt, alphapt

  IMPLICIT NONE
  INTEGER :: ik, ipol, ldim
  !
  !  allocate space for the quantities needed in the phonon program
  !
  if (lgamma) then
     !
     !  q=0  : evq is a pointer to evc
     !
     evq  => evc
  else
     !
     !  q!=0 : evq is allocated and calculated at point k+q
     !
     allocate (evq ( npwx*npol , nbnd))
  endif
  !
  allocate (dvpsi ( npwx*npol , nbnd))
  allocate ( dpsi ( npwx*npol , nbnd))
  !
  allocate (vlocq ( ngm , ntyp))
  allocate (dmuxc ( dfftp%nnr , nspin_mag , nspin_mag))
  !
  allocate (eigqts ( nat))
  allocate (rtau ( 3, 48, nat))
  allocate (u ( 3 * nat, 3 * nat))
  allocate (dyn ( 3 * nat, 3 * nat))
  allocate (dyn_rec ( 3 * nat, 3 * nat))
  allocate (dyn00 ( 3 * nat, 3 * nat))
  allocate (w2 ( 3 * nat))
  allocate (name_rap_mode( 3 * nat))
  allocate (num_rap_mode( 3 * nat ))
  allocate (npert ( 3 * nat))
  ALLOCATE (xk_col(3, nkstot))
  allocate (polar (3, 3,  nfs))
  allocate (zstareu (3, 3,  nat))
  allocate (zstareu0 (3, 3 * nat))
  allocate (zstarue (3 , nat, 3))
  allocate (zstarue0 (3 * nat, 3))
  allocate (zstarue0_rec (3 * nat, 3))
  name_rap_mode=' '
  zstarue=0.0_DP
  zstareu0=(0.0_DP,0.0_DP)
  zstarue0=(0.0_DP,0.0_DP)
  zstarue0_rec=(0.0_DP,0.0_DP)
  IF (noncolin.AND.domag) THEN
     ALLOCATE (becpt(nksq))
     ALLOCATE (alphapt(3,nksq))
     DO ik=1,nksq
        CALL allocate_bec_type ( nkb, nbnd, becpt(ik) )
        DO ipol=1,3
           CALL allocate_bec_type ( nkb, nbnd, alphapt(ipol,ik) )
        ENDDO
     ENDDO
     IF (okvan) THEN
        ALLOCATE(int1_nc_save( nhm, nhm, 3, nat, nspin, 2))
        ALLOCATE (deeq_nc_save( nhm, nhm, nat, nspin, 2))
     ENDIF
     ! AAA Nota: da definire this_pcxpsi_is_on_file_tpw
     !ALLOCATE (this_pcxpsi_is_on_file_tpw(nksq,3,2))
     ! ELSE
     !ALLOCATE (this_pcxpsi_is_on_file_tpw(nksq,3,1))
  ENDIF

  if (okvan) then
     allocate (int1 ( nhm, nhm, 3, nat, nspin_mag))
     allocate (int2 ( nhm , nhm , 3 , nat , nat))
     if (okpaw) then
        allocate (becsumort ( nhm*(nhm+1)/2 , nat , nspin, 3*nat))
     endif
     allocate (int4 ( nhm * (nhm + 1)/2,  3 , 3 , nat, nspin_mag))
     allocate (int5 ( nhm * (nhm + 1)/2 , 3 , 3 , nat , nat))
     allocate (dpqq( nhm, nhm, 3, ntyp))
     IF (noncolin) THEN
        ALLOCATE(int1_nc( nhm, nhm, 3, nat, nspin))
        ALLOCATE(int4_nc( nhm, nhm, 3, 3, nat, nspin))
        ALLOCATE(becsum_nc( nhm*(nhm+1)/2, nat, npol, npol))
        ALLOCATE(alphasum_nc( nhm*(nhm+1)/2, 3, nat, npol, npol))
        IF (lspinorb) THEN
           ALLOCATE(int2_so( nhm, nhm, 3, nat , nat, nspin))
           ALLOCATE(int5_so( nhm, nhm, 3, 3, nat , nat, nspin))
           allocate(dpqq_so( nhm, nhm, nspin, 3, ntyp))
        END IF
     END IF


     allocate (alphasum ( nhm * (nhm + 1)/2 , 3 , nat , nspin_mag))
     allocate (this_dvkb3_is_on_file(nksq))
     this_dvkb3_is_on_file(:)=.false.
  endif
  allocate (this_pcxpsi_is_on_file(nksq,3))
  this_pcxpsi_is_on_file(:,:)=.false.

  ALLOCATE (becp1(nksq))
  ALLOCATE (alphap(3,nksq))
  DO ik=1,nksq
     call allocate_bec_type ( nkb, nbnd, becp1(ik) )
     DO ipol=1,3
        call allocate_bec_type ( nkb, nbnd, alphap(ipol,ik) )
     ENDDO
  END DO
  CALL allocate_bec_type ( nkb, nbnd, becp )

  if (elph) then
    allocate (el_ph_mat( nbnd, nbnd, nksq, 3*nat))
  endif
  allocate ( ramtns (3, 3, 3, nat) )

  ! DFPT+U
  IF (lda_plus_u) THEN
     !
     ldim = 2 * Hubbard_lmax + 1
     !
     ALLOCATE (wfcatomk(npwx,nwfcU))
     ALLOCATE (swfcatomk(npwx,nwfcU))
     ALLOCATE (dwfcatomk(npwx,nwfcU,3))
     ALLOCATE (sdwfcatomk(npwx,nwfcU))
     ALLOCATE (dvkb(npwx,nkb,3))
     ALLOCATE (ofsbeta(nat))
     !
     ALLOCATE (dnsbare(ldim,ldim,nspin,nat,3,nat))
     ALLOCATE (dnsbare_all_modes(ldim,ldim,nspin,nat,3*nat))
     ALLOCATE (dnsorth(ldim,ldim,nspin,nat,3*nat))
     !
     IF (lgamma) THEN
        wfcatomkpq   => wfcatomk
        swfcatomkpq  => swfcatomk
        dwfcatomkpq  => dwfcatomk
        sdwfcatomkpq => sdwfcatomk
        vkbkpq       => vkb
        dvkbkpq      => dvkb
     ELSE
        ALLOCATE (wfcatomkpq(npwx,nwfcU))
        ALLOCATE (swfcatomkpq(npwx,nwfcU))
        ALLOCATE (dwfcatomkpq(npwx,nwfcU,3))
        ALLOCATE (sdwfcatomkpq(npwx,nwfcU))
        ALLOCATE (vkbkpq(npwx,nkb))
        ALLOCATE (dvkbkpq(npwx,nkb,3))
     ENDIF
     !
  ENDIF 

  RETURN

end subroutine allocate_phq
