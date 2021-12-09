!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvpsi_e (ik, ipol)
  !----------------------------------------------------------------------
  !! On output: \(\text{dvpsi}\) contains \(P_c^+ x | \psi_{ik}\rangle\)
  !! in crystal axis (projected on \(\text{at}(*,\text{ipol})\)).
  !
  !! \(\text{dvpsi}\) is read from file if
  !! \(\text{this_pcxpsi_is_on_file}(\text{ik},\text{ipol})=\text{TRUE}\),
  !! otherwise \(\text{dvpsi}\) is COMPUTED and WRITTEN on file 
  !! (\(\text{vkb}\) and \(\text{evc}\) must be set).
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba2, at
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : xk, ngk, igk_k
  USE gvect,           ONLY : g
  USE wvfct,           ONLY : npwx, nbnd, et
  USE wavefunctions, ONLY: evc
  USE buffers,         ONLY : save_buffer, get_buffer
  USE noncollin_module,ONLY : noncolin, npol
  USE becmod,          ONLY : bec_type, becp, calbec, &
                              allocate_bec_type, deallocate_bec_type
  USE mp_bands,        ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE xc_lib,          ONLY : exx_is_active
  USE uspp,            ONLY : okvan, nkb, vkb
  USE uspp_param,      ONLY : nh, nhm
  USE ramanm,          ONLY : eth_rps
  USE units_ph,        ONLY : this_pcxpsi_is_on_file, lrcom, iucom, &
                              lrebar, iuebar

  USE lrus,            ONLY : becp1
  USE qpoint,          ONLY : nksq, ikks
  USE eqv,             ONLY : dpsi, dvpsi
  USE control_lr,      ONLY : nbnd_occ

  implicit none
  !
  integer, intent(IN) :: ipol, ik
  !
  ! ... Local variables
  !
  integer :: npw
  integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0,  &
             nrec, is, js, ijs, ikk
  ! counters

  real(DP), allocatable  :: h_diag (:,:)
  ! the diagonal part of h_scf
  type(bec_type) :: becp2 ! the scalar products
  complex(DP), allocatable :: spsi(:,:)
  real(DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter
  logical :: conv_root
  ! true if convergence has been achieved
  INTEGER :: n_start, n_end

  external ch_psi_all, cg_psi
  !
  call start_clock ('dvpsi_e')
  dpsi=(0.d0, 0.d0)
  dvpsi=(0.d0, 0.d0)
  if (this_pcxpsi_is_on_file(ik,ipol)) then
     nrec = (ipol - 1)*nksq + ik
     call get_buffer(dvpsi, lrebar, iuebar, nrec)
     call stop_clock ('dvpsi_e')
     return
  end if
  !
  ikk=ikks(ik)
  call allocate_bec_type ( nkb, nbnd, becp2)

  ! calculate the commutator [H,x_ipol]  psi > and store it in dpsi

  call commutator_Hx_psi (ikk, nbnd_occ(ikk), at(:, ipol), becp1(ik), becp2, dpsi(:, 1:nbnd_occ(ikk)) )
  !
  !    orthogonalize dpsi to the valence subspace: ps = <evc|dpsi>
  !    Apply -P^+_c
  !    NB it uses dvpsi as workspace
  !
  npw = ngk(ikk)
  CALL orthogonalize(dpsi, evc, ikk, ikk, dvpsi, npw, .false.)
  dpsi=-dpsi
  !
  !   dpsi contains P^+_c [H-eS,x] psi_v for the three crystal polarizations
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !
  CALL g2_kin(ikk)
  !
  ! compute preconditioning matrix h_diag used by cgsolve_all
  !
  allocate (h_diag( npwx*npol, nbnd))
  CALL h_prec (ik, evc, h_diag)
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  !
  thresh = eth_rps
  call cgsolve_all (ch_psi_all, cg_psi, et (1, ikk), dpsi, dvpsi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ (ikk), npol)

  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",es10.3)') &
       ik, ibnd, anorm
  !
  FLUSH( stdout )
  deallocate (h_diag)
  !
  ! we have now obtained P_c x |psi>.
  ! In the case of USPP this quantity is needed for the Born
  ! effective charges, so we save it to disc
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  !
  if (okvan) then
     !
     ! for effective charges
     !
     nrec = (ipol - 1) * nksq + ik
     call save_buffer(dvpsi, lrcom, iucom, nrec)
     !
     allocate (spsi ( npwx*npol, nbnd))
     if (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .and. nbnd>1 ) then
        call divide(inter_bgrp_comm,nbnd,n_start,n_end)
        if (n_end >= n_start) CALL calbec (npw, vkb, dvpsi(:,n_start:n_end), becp , n_end-n_start+1 )
     else
        CALL calbec (npw, vkb, dvpsi, becp )
     end if
     CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
     call dcopy(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
     deallocate (spsi)
     CALL adddvepsi_us(becp1(ik),becp2,ipol,ikk,dvpsi)
  endif

  IF (nkb > 0) call deallocate_bec_type (becp2)

  nrec = (ipol - 1)*nksq + ik
  call save_buffer(dvpsi, lrebar, iuebar, nrec)
  this_pcxpsi_is_on_file(ik,ipol) = .true.
  call stop_clock ('dvpsi_e')
  return
end subroutine dvpsi_e
