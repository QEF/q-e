!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine lr_h_psiq (lda, n, m, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !
  ! Modified by Osman Baris Malcioglu (2009)
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !

  USE kinds,  ONLY : DP
  USE wavefunctions_module,  ONLY : psic, psic_nc
  USE noncollin_module, ONLY : noncolin, npol
  USE lsda_mod, ONLY : current_spin
  use fft_base,             only : dffts
  use fft_interfaces,       only : fwfft, invfft
  USE gvecs,  ONLY : nls
  USE spin_orb, ONLY : domag
  USE scf,    ONLY : vrs
  USE uspp,   ONLY : vkb
  USE wvfct,  ONLY : g2kin,igk
  USE lr_variables,   ONLY : lr_verbosity
  use control_flags,         only : gamma_only
  use io_global,            only : stdout
  !USE qpoint, ONLY : igkq
  implicit none
  !
  !     Here the local variables
  !
  integer :: ibnd
  ! counter on bands

  integer :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  integer :: j
  ! do loop index

  complex(DP) :: psi (lda*npol, m), hpsi (lda*npol, m), spsi (lda*npol, m)
  complex(DP) :: sup, sdwn
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)
  !OBM debug
  !real(DP) :: obm_debug
  !complex(kind=dp), external :: ZDOTC


  call start_clock ('h_psiq')
  If (lr_verbosity > 5) WRITE(stdout,'("<lr_h_psiq>")')
  if (gamma_only) then
   call lr_h_psiq_gamma()
  else
   call lr_h_psiq_k()
  endif
  call stop_clock ('h_psiq')
  return
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!k point part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lr_h_psiq_k()

    USE becmod, ONLY : bec_type, becp, calbec

    IMPLICIT NONE
  call start_clock ('init')

  call calbec ( n, vkb, psi, becp, m)
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  hpsi=(0.d0,0.d0)
  do ibnd = 1, m
     do j = 1, n
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO j = 1, n
           hpsi (j+lda, ibnd) = g2kin (j) * psi (j+lda, ibnd)
        ENDDO
     ENDDO
  ENDIF
  call stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

  do ibnd = 1, m
     call start_clock ('firstfft')
     IF (noncolin) THEN
        psic_nc = (0.d0, 0.d0)
        do j = 1, n
           psic_nc(nls(igk(j)),1) = psi (j, ibnd)
           psic_nc(nls(igk(j)),2) = psi (j+lda, ibnd)
        enddo
        CALL invfft ('Wave', psic_nc(:,1), dffts)
        CALL invfft ('Wave', psic_nc(:,2), dffts)
     ELSE
        psic(:) = (0.d0, 0.d0)
        do j = 1, n
           psic (nls(igk(j))) = psi (j, ibnd)
        enddo
        CALL invfft ('Wave', psic, dffts)
     END IF
     call stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
     if (noncolin) then
        if (domag) then
           do j=1, dffts%nnr
              sup = psic_nc(j,1) * (vrs(j,1)+vrs(j,4)) + &
                    psic_nc(j,2) * (vrs(j,2)-(0.d0,1.d0)*vrs(j,3))
              sdwn = psic_nc(j,2) * (vrs(j,1)-vrs(j,4)) + &
                    psic_nc(j,1) * (vrs(j,2)+(0.d0,1.d0)*vrs(j,3))
              psic_nc(j,1)=sup
              psic_nc(j,2)=sdwn
           end do
        else
           do j=1, dffts%nnr
              psic_nc(j,1)=psic_nc(j,1) * vrs(j,1)
              psic_nc(j,2)=psic_nc(j,2) * vrs(j,1)
           enddo
        endif
     else
        do j = 1, dffts%nnr
           psic (j) = psic (j) * vrs (j, current_spin)
        enddo
     endif
     !
     !   back to reciprocal space
     !
     call start_clock ('secondfft')
     IF (noncolin) THEN
        CALL fwfft ('Wave', psic_nc(:,1), dffts)
        CALL fwfft ('Wave', psic_nc(:,2), dffts)
     !
     !   addition to the total product
     !
        do j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic_nc (nls(igk(j)), 1)
           hpsi (j+lda, ibnd) = hpsi (j+lda, ibnd) + psic_nc (nls(igk(j)), 2)
        enddo
     ELSE
        CALL fwfft ('Wave', psic, dffts)
     !
     !   addition to the total product
     !
        do j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic (nls(igk(j)))
        enddo
     END IF
     call stop_clock ('secondfft')
  enddo
  !
  !  Here the product with the non local potential V_NL psi
  !

  call add_vuspsi (lda, n, m, hpsi)

  call s_psi (lda, n, m, psi, spsi)

    end subroutine lr_h_psiq_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gamma point part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lr_h_psiq_gamma()

    USE becmod, ONLY : becp, calbec
    USE gvect,  ONLY : gstart
    USE realus, ONLY : real_space, fft_orbital_gamma, &
                       bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                       v_loc_psir, s_psir_gamma, real_space_debug
    use uspp,                  only : nkb

    IMPLICIT NONE

    call start_clock ('init')
    !
    ! Here we apply the kinetic energy (k+G)^2 psi
    !
    if(gstart==2) psi(1,:)=cmplx(real(psi(1,:),dp),0.0d0,dp)
    !
    !!OBM debug
    !  obm_debug=0
    !  do ibnd=1,m
    !     !
    !     obm_debug=obm_debug+ZDOTC(lda*npol,psi(:,ibnd),1,psi(:,ibnd),1)
    !     !
    !  enddo
    !  print *, "lr_h_psiq psi", obm_debug
    !!obm_debug

    do ibnd=1,m
       do j=1,n
          hpsi(j,ibnd)=g2kin(j)*psi(j,ibnd)
       enddo
    enddo
    !!OBM debug
    !  obm_debug=0
    !  do ibnd=1,m
    !     !
    !     obm_debug=obm_debug+ZDOTC(lda*npol,hpsi(:,ibnd),1,hpsi(:,ibnd),1)
    !     !
    !  enddo
    !  print *, "lr_h_psiq hpsi (just after kinetic operator)", obm_debug
    !!obm_debug

    call stop_clock ('init')
      if (nkb > 0 .and. real_space_debug>2) then
        do ibnd=1,m,2
          !call check_fft_orbital_gamma(psi,ibnd,m)
          call fft_orbital_gamma(psi,ibnd,m,.true.) !transform the psi real space, saved in temporary memory
          call calbec_rs_gamma(ibnd,m,becp%r) !rbecp on psi
          call s_psir_gamma(ibnd,m) !psi -> spsi
          call bfft_orbital_gamma(spsi,ibnd,m) !return back to real space
          call fft_orbital_gamma(hpsi,ibnd,m) ! spsi above is now replaced by hpsi
          call v_loc_psir(ibnd,m) ! hpsi -> hpsi + psi*vrs  (psi read from temporary memory)
          call add_vuspsir_gamma(ibnd,m) ! hpsi -> hpsi + vusp
          call bfft_orbital_gamma(hpsi,ibnd,m,.true.) !transform back hpsi, clear psi in temporary memory
        enddo
     ELSE
    call vloc_psi_gamma(lda,n,m,psi,vrs(1,current_spin),hpsi)
    !!OBM debug
    !  obm_debug=0
    !  do ibnd=1,m
    !     !
    !     obm_debug=obm_debug+ZDOTC(lda*npol,hpsi(:,ibnd),1,hpsi(:,ibnd),1)
    !     !
    !  enddo
    !  print *, "lr_h_psiq hpsi (after vloc_psi)", obm_debug
    !!obm_debug

     IF (noncolin) THEN
       call errore ("lr_h_psiq","gamma and noncolin not implemented yet",1)
     ELSE
        call calbec ( n, vkb, psi, becp, m)
     END IF
     !!OBM debug
     ! obm_debug=0
     ! do ibnd=1,m
     !    !
     !    obm_debug=obm_debug+ZDOTC(lda*npol,becp%r(:,ibnd),1,becp%r(:,ibnd),1)
     !    !
     ! enddo
     ! print *, "lr_h_psiq becp", obm_debug
     !!obm_debug


     call add_vuspsi (lda, n, m, hpsi)
     !END IF
     !!OBM debug
     !  obm_debug=0
     !  do ibnd=1,m
     !     !
     !     obm_debug=obm_debug+ZDOTC(lda*npol,hpsi(:,ibnd),1,hpsi(:,ibnd),1)
     !     !
     !  enddo
     !  print *, "lr_h_psiq hpsi (after add_vuspsi)", obm_debug
     !!obm_debug

     call s_psi (lda, n, m, psi, spsi)
     !!OBM debug
     !  obm_debug=0
     !  do ibnd=1,m
     !     !
     !     obm_debug=obm_debug+ZDOTC(lda*npol,spsi(:,ibnd),1,spsi(:,ibnd),1)
     !     !
     !  enddo
     !  print *, "lr_h_psiq spsi ", obm_debug
    !!obm_debug
    ENDIF
    end subroutine lr_h_psiq_gamma

end subroutine lr_h_psiq
