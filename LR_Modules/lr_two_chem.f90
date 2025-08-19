!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE  lr_two_chem
  USE kinds, ONLY : DP
  COMPLEX(DP),SAVE,PUBLIC :: def_val(3)
  COMPLEX(DP),SAVE,PUBLIC :: def_cond(3)
  ! the change of the Fermi energy for each pert. for valence and conduction manifold in the twochem case.
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: drhos_cond(:,:,:)
  !! output: the change of the scf charge
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: drhop_cond(:,:,:)
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: dbecsum_cond(:,:,:,:),dbecsum_cond_nc(:,:,:,:,:,:)
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: ldos_cond(:,:),ldoss_cond(:,:)
  REAL(DP), SAVE, PUBLIC               :: dos_ef_cond
  REAL(DP), ALLOCATABLE, SAVE,PUBLIC    :: becsum1_cond(:,:,:)
  REAL(DP), ALLOCATABLE, SAVE, PUBLIC   :: becsum_cond(:,:,:) ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>, i in the conduction manifold
  COMPLEX (DP), ALLOCATABLE :: becsum_cond_nc(:,:,:,:)     !conduction manifold
  COMPLEX (DP), ALLOCATABLE :: becsumort_cond(:,:,:,:)
  !! it contains \(\text{alphasum}+\sum_i \langle\psi_i | \beta_n\rangle\langle\beta_m| \delta \psi_i \rangle\) (conduction
  !manifold)
  COMPLEX (DP), ALLOCATABLE :: alphasum_cond_nc(:,:,:,:,:)   ! nhm*(nhm+1)/2,3,nat,npol,npol)
  REAL (DP), ALLOCATABLE :: alphasum_cond(:,:,:,:) ! (nhm*(nhm+1)/2,3,nat,nspin)
  !! used to compute modes. It contains \(\sum_i \langle \psi_i| d/du
  !! (|\beta_n><beta_m|) | \psi_i\rangle + (m-n)\) for the conduction manifold
  CONTAINS
!
!-----------------------------------------------------------------------
subroutine ef_shift_twochem (npert, dos_ef,dos_ef_cond,ldos,ldos_cond,&
           drhop,drhop_cond,dbecsum,dbecsum_cond,becsum1,becsum1_cond)
  !-----------------------------------------------------------------------
  !! This routine takes care of the effects of a shift of the two chemical potentials, due to the
  !! perturbation, that can take place in a metal at q=0, in the twochem case
  !! Optionally, update dbecsum using becsum1.
  !
  USE kinds,                ONLY : DP
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg
  USE buffers,              ONLY : get_buffer, save_buffer
  USE uspp_param,           ONLY : nhm
  USE noncollin_module,     ONLY : nspin_mag, nspin_lsda
  !
  IMPLICIT NONE
  !
  ! input/output variables
  !
  INTEGER, INTENT(IN) :: npert
  !! the number of perturbation
  REAL(DP), INTENT(IN) :: dos_ef,dos_ef_cond
  !! density of states at the two chemical potentials
  COMPLEX(DP), INTENT(IN) :: ldos(dfftp%nnr, nspin_mag),ldos_cond(dfftp%nnr, nspin_mag)
  !! local DOS at Ef (with augmentation)
  COMPLEX(DP), INTENT(INOUT) :: drhop(dfftp%nnr, nspin_mag, npert)
  COMPLEX(DP), INTENT(INOUT) :: drhop_cond(dfftp%nnr, nspin_mag,npert)
  !! the change of the charge (with augmentation)
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum((nhm*(nhm+1))/2, nat, nspin_mag, npert)
  !! input:  dbecsum = 2 <psi|beta> <beta|dpsi>
  !! output: dbecsum = 2 <psi|beta> <beta|dpsi> + def * becsum1
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_cond((nhm*(nhm+1))/2, nat, nspin_mag, npert)
  !same as above, but restricted to conduction states
  REAL(DP), INTENT(IN), OPTIONAL :: becsum1((nhm*(nhm+1))/2, nat, nspin_mag)
  !! becsum1 = wdelta * <psi|beta> <beta|psi>
  REAL(DP), INTENT(IN), OPTIONAL :: becsum1_cond((nhm*(nhm+1))/2, nat, nspin_mag)
  !! becsum1_cond = wdelta_cond * <psi|beta> <beta|psi>
  !! (where wdelta is a Dirac-delta-like function)
  !
  ! local variables
  !
  INTEGER :: is
  !! counter on spin polarizations
  INTEGER :: ipert
  !! counter on perturbations
  COMPLEX(DP) :: drhop_val(dfftp%nnr, nspin_mag,npert)
  COMPLEX(DP) :: delta_nv,delta_nc
  !! the change in electron number
  !! This may be complex since perturbation may be complex
  !
  REAL(DP), external :: w0gauss
  !! the smeared delta function
  !
  call start_clock ('ef_shift_twochem')
  !
  ! This routine is used only at q=Gamma where the dimension of irrep never exceeds 3
  IF (npert > 3) CALL errore("ef_shift_twochem", "npert exceeds 3", 1)
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  WRITE( stdout, * )
  drhop_val(:,:,:) = drhop(:,:,:)-drhop_cond(:,:,:)
  do ipert = 1, npert
     delta_nv = (0.d0, 0.d0)
     delta_nc = (0.d0, 0.d0)
     do is = 1, nspin_lsda
        CALL fwfft ('Rho', drhop_val(:,is,ipert), dfftp)
        if (gg(1) < 1.0d-8) delta_nv = delta_nv + omega*drhop_val(dfftp%nl(1),is,ipert)
        CALL invfft ('Rho', drhop_val(:,is,ipert), dfftp)
        !valence states
        CALL fwfft ('Rho', drhop_cond(:,is,ipert), dfftp)
        if (gg(1) < 1.0d-8) delta_nc = delta_nc + omega*drhop_cond(dfftp%nl(1),is,ipert)
        CALL invfft ('Rho', drhop_cond(:,is,ipert), dfftp)
     enddo
     call mp_sum ( delta_nv, intra_bgrp_comm )
     call mp_sum ( delta_nc, intra_bgrp_comm )
     IF ( ABS(dos_ef+dos_ef_cond) > 1.d-18 ) THEN
        def_val (ipert) = - delta_nv / dos_ef
        def_cond (ipert) = - delta_nc / dos_ef_cond
     ELSE
        def_val (ipert) = 0.0_dp
        def_cond (ipert) = 0.0_dp
     ENDIF
  enddo
  !
  ! symmetrizes the Fermi energy shift
  !
  CALL sym_def(def_val)
  WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift valence (Ry) =",2es15.4)')&
       (ipert, def_val (ipert) , ipert = 1, npert )
  CALL sym_def(def_cond)
  WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift conduction (Ry) =",2es15.4)')&
       (ipert, def_cond (ipert) , ipert = 1, npert )
  !
  ! corrects the density response accordingly...
  !
  do ipert = 1, npert
     call zaxpy (dfftp%nnr*nspin_mag, def_val(ipert), ldos, 1, drhop(1,1,ipert), 1)
     call zaxpy (dfftp%nnr*nspin_mag, def_cond(ipert), ldos_cond, 1, drhop(1,1,ipert), 1)
  enddo
  !
  ! In the PAW case there is also a metallic term
  !
  IF (PRESENT(dbecsum) .AND. PRESENT(becsum1)) THEN
     DO ipert = 1, npert
        dbecsum(:,:,:,ipert) = dbecsum(:,:,:,ipert) &
           + def_val(ipert) * CMPLX(becsum1(:,:,:), 0.0_DP, KIND=DP)
     ENDDO
     DO ipert = 1, npert
        dbecsum(:,:,:,ipert) = dbecsum(:,:,:,ipert) &
           + def_cond(ipert) * CMPLX(becsum1_cond(:,:,:), 0.0_DP, KIND=DP)
     ENDDO

  ENDIF
  !
  CALL stop_clock ('ef_shift_twochem')
  !
  end subroutine ef_shift_twochem
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
subroutine ef_shift_wfc_twochem(npert, ldoss,ldoss_cond, drhos)
  !-----------------------------------------------------------------------
  !! This routine takes care of the effects of a shift of both chemical potentials, due to the
  !! perturbation, that can take place in a metal at q=0, on the wavefunctions.
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_sum
  USE wavefunctions,        ONLY : evc
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE buffers,              ONLY : get_buffer, save_buffer
  USE wvfct,                ONLY : npwx, et, nbnd, nbnd_cond
  USE klist,                ONLY : degauss, ngauss, ngk, ltetra,degauss_cond
  USE ener,                 ONLY : ef,ef_cond
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : nbnd_occ
  USE units_lr,             ONLY : iuwfc, lrwfc, lrdwf, iudwf
  USE eqv,                  ONLY : dpsi
  USE dfpt_tetra_mod,       ONLY : dfpt_tetra_delta
  !
  IMPLICIT NONE
  !
  ! input/output variables
  !
  INTEGER, INTENT(IN) :: npert
  !! the number of perturbation
  COMPLEX(DP), INTENT(IN) :: ldoss(dffts%nnr, nspin_mag)
  !! local DOS at Ef without augmentation
  COMPLEX(DP), INTENT(IN) :: ldoss_cond(dffts%nnr, nspin_mag)
  !! local DOS at ef_cond without augmentation, for the conduction chemical potential
  COMPLEX(DP), INTENT(INOUT) :: drhos(dffts%nnr, nspin_mag, npert)
  !! the change of the charge (with augmentation)
  !
  ! local variables
  !
  INTEGER :: npw, ibnd, ik, is, ipert, nrec, ikrec
  ! counter on occupied bands
  ! counter on k-point
  ! counter on spin polarizations
  ! counter on perturbations
  ! record number
  ! record position of wfc at k
  ! auxiliary for spin
  COMPLEX(DP) :: wfshift
  !! the shift coefficient for the wavefunction
  !! This may be complex since perturbation may be complex
  !
  REAL(DP), external :: w0gauss
  ! the smeared delta function
  !
  call start_clock ('ef_shift_wfc_twochem')
  !
  ! This routine is used only at q=Gamma where the dimension of irrep never exceeds 3
  IF (npert > 3) CALL errore("ef_shift_wfc_twochem", "npert exceeds 3", 1)
  !
  ! Update the perturbed wavefunctions according to the Fermi energy shift
  !
  do ik = 1, nksq
     npw = ngk (ik)
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
     ikrec = ik
     if (nksq > 1) call get_buffer (evc, lrwfc, iuwfc, ikrec)
     !
     ! reads delta_psi from iunit iudwf, k=kpoint
     !
     do ipert = 1, npert
        nrec = (ipert - 1) * nksq + ik
        IF (nksq > 1 .OR. npert > 1) CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
        do ibnd = 1, nbnd_occ (ik)
                        if (ibnd.le.(nbnd-nbnd_cond)) then
                 wfshift = 0.5d0 * def_val(ipert) * &
                 w0gauss( (ef-et(ibnd,ik))/degauss, ngauss) / degauss
                 else
                 wfshift = 0.5d0 * def_cond(ipert) * &
                 w0gauss((ef_cond-et(ibnd,ik))/degauss_cond, ngauss) / degauss_cond
           end if
           !
           IF (noncolin) THEN
              call zaxpy (npwx*npol,wfshift,evc(1,ibnd),1,dpsi(1,ibnd),1)
           ELSE
              call zaxpy (npw, wfshift, evc(1,ibnd), 1, dpsi(1,ibnd), 1)
           ENDIF
        enddo
        !
        ! writes corrected delta_psi to iunit iudwf, k=kpoint,
        !
        IF (nksq > 1 .OR. npert > 1) CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
     enddo
  enddo
  !
  do ipert = 1, npert
     do is = 1, nspin_mag
        call zaxpy (dffts%nnr, def_val(ipert), ldoss(1,is), 1, drhos(1,is,ipert), 1)
        call zaxpy (dffts%nnr, def_cond(ipert), ldoss_cond(1,is), 1, drhos(1,is,ipert), 1)
     enddo
  enddo
  !
  CALL stop_clock ('ef_shift_wfc_twochem')
  !
  end subroutine ef_shift_wfc_twochem
  !
!-----------------------------------------------------------------------
subroutine localdos_cond (ldos, ldoss, becsum1, dos_ef)
  !-----------------------------------------------------------------------
  !
  !    This routine compute the local and total density of states
  !    for the conduction chemical potential in the twochem case
  !
  !    Note: this routine use psic as auxiliary variable. it should alread
  !          be defined
  !
  !    NB: this routine works only with gamma
  !
  !
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : omega
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE ener,             ONLY : ef_cond
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : invfft, fft_interpolate
  USE buffers,          ONLY : get_buffer
  USE gvecs,            ONLY : doublegrid
  USE klist,            ONLY : xk, wk, ngk, igk_k, degauss, ngauss, ltetra
  USE lsda_mod,         ONLY : nspin, lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wvfct,            ONLY : nbnd, npwx, et, nbnd_cond
  USE becmod,           ONLY : calbec, bec_type, allocate_bec_type_acc, deallocate_bec_type_acc
  USE wavefunctions,    ONLY : evc, psic, psic_nc
  USE uspp,             ONLY : okvan, nkb, vkb
  USE uspp_param,       ONLY : upf, nh, nhm
  USE qpoint,           ONLY : nksq, ikks
  USE control_lr,       ONLY : nbnd_occ
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE mp_pools,         ONLY : inter_pool_comm
  USE mp,               ONLY : mp_sum
  USE dfpt_tetra_mod,   ONLY : dfpt_tetra_delta
  USE uspp_init,        ONLY : init_us_2
  USE control_flags,    ONLY : offload_type

  implicit none

  complex(DP) :: ldos (dfftp%nnr, nspin_mag), ldoss (dffts%nnr, nspin_mag)
  ! output: the local density of states at Ef
  ! output: the local density of states at Ef without augmentation
  REAL(DP) :: becsum1 ((nhm * (nhm + 1))/2, nat, nspin_mag)
  ! output: the local becsum at ef_cond, for conduction fermi level
  real(DP) :: dos_ef, check
  ! output: the density of states at Ef
  !
  !    local variables for Ultrasoft PP's
  !
  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, nt
  ! counters
  complex(DP), allocatable :: becsum1_nc(:,:,:,:)
  TYPE(bec_type) :: becp
  !
  ! local variables
  !
  real(DP) :: weight, w1, wdelta
  ! weights
  real(DP), external :: w0gauss
  !
  integer :: npw, ik, is, ig, ibnd, j, is1, is2, v_siz
  ! counters
  integer :: ios
  ! status flag for i/o
  !
  !  initialize ldos and dos_ef
  !
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
  !$acc enter data copyin(evc, nl_d)
  v_siz = dffts%nnr

  call start_clock ('localdos_cond')
  IF (noncolin) THEN
     allocate (becsum1_nc( (nhm * (nhm + 1)) / 2, nat, npol, npol))
     becsum1_nc=(0.d0,0.d0)
  ENDIF

  call allocate_bec_type_acc (nkb, nbnd, becp)

  becsum1 (:,:,:) = 0.d0
  ldos (:,:) = (0d0, 0.0d0)
  ldoss(:,:) = (0d0, 0.0d0)
  dos_ef = 0.d0
  !
  !  loop over kpoints
  !
  !$acc data create(psic, psic_nc) copy(ldoss)
  do ik = 1, nksq
     if (lsda) current_spin = isk (ikks(ik))
     npw = ngk(ikks(ik))
     weight = wk (ikks(ik))
     !
     ! unperturbed wfs in reciprocal space read from unit iuwfc
     !
     if (nksq > 1) then
             call get_buffer (evc, lrwfc, iuwfc, ikks(ik))
             !$acc update device(evc)
     endif
     call init_us_2 (npw, igk_k(1,ikks(ik)), xk (1, ikks(ik)), vkb, .true.)
     !
     !$acc data present(vkb, becp, evc)
     call calbec ( offload_type, npw, vkb, evc, becp)
     !$acc end data
     !
     do ibnd = 1+(nbnd-nbnd_cond), nbnd_occ (ikks(ik))
        !
        if(ltetra) then
           wdelta = dfpt_tetra_delta(ibnd,ikks(ik))
        else
           wdelta = w0gauss ( (ef_cond-et(ibnd,ikks(ik))) / degauss, ngauss) / degauss
        end if
        !
        w1 = weight * wdelta / omega
        !
        ! unperturbed wf from reciprocal to real space
        !
        IF (noncolin) THEN
           !$acc kernels
           psic_nc = (0.d0, 0.d0)
           !$acc end kernels
           !$acc parallel loop present(igk_k, psic_nc, nl_d, evc)
           do ig = 1, npw
              psic_nc (nl_d (igk_k(ig,ikks(ik))), 1 ) = evc (ig, ibnd)
              psic_nc (nl_d (igk_k(ig,ikks(ik))), 2 ) = evc (ig+npwx, ibnd)
           enddo
           !$acc end parallel loop
           !$acc host_data use_device(psic_nc)
           CALL invfft ('Rho', psic_nc(:,1), dffts)
           CALL invfft ('Rho', psic_nc(:,2), dffts)
           !$acc end host_data
           !$acc parallel loop present(psic_nc, ldoss)
           do j = 1, v_siz
              ldoss (j, 1) = ldoss (j, 1) + &
                    w1 * ( DBLE(psic_nc(j,1))**2+AIMAG(psic_nc(j,1))**2 + &
                           DBLE(psic_nc(j,2))**2+AIMAG(psic_nc(j,2))**2)
           enddo
           !$acc end parallel loop
           IF (nspin_mag==4) THEN
              !$acc parallel loop present(psic_nc, ldoss)
              DO j = 1, v_siz
              !
                 ldoss(j,2) = ldoss(j,2) + w1*2.0_DP* &
                             (DBLE(psic_nc(j,1))* DBLE(psic_nc(j,2)) + &
                             AIMAG(psic_nc(j,1))*AIMAG(psic_nc(j,2)))

                 ldoss(j,3) = ldoss(j,3) + w1*2.0_DP* &
                             (DBLE(psic_nc(j,1))*AIMAG(psic_nc(j,2)) - &
                              DBLE(psic_nc(j,2))*AIMAG(psic_nc(j,1)))

                 ldoss(j,4) = ldoss(j,4) + w1* &
                             (DBLE(psic_nc(j,1))**2+AIMAG(psic_nc(j,1))**2 &
                             -DBLE(psic_nc(j,2))**2-AIMAG(psic_nc(j,2))**2)
              !
              END DO
              !$acc end parallel loop
           END IF
        ELSE
           !$acc kernels
           psic (:) = (0.d0, 0.d0)
           !$acc end kernels
           !$acc parallel loop present(psic,nl_d, evc)
           do ig = 1, npw
              psic (nl_d (igk_k(ig,ikks(ik)) ) ) = evc (ig, ibnd)
           enddo
           !$acc end parallel loop
           !$acc host_data use_device(psic)
           CALL invfft ('Rho', psic, dffts)
           !$acc end host_data
           !$acc parallel loop present(ldoss, psic)
           do j = 1, v_siz
              ldoss (j, current_spin) = ldoss (j, current_spin) + &
                    w1 * ( DBLE ( psic (j) ) **2 + AIMAG (psic (j) ) **2)
           enddo
           !$acc end parallel loop
        END IF
        !
        !    If we have a US pseudopotential we compute here the becsum term
        !
        if(noncolin) then
        !$acc update self(becp%nc)
        else
        !$acc update self(becp%k)
        endif
        !
        w1 = weight * wdelta
        ijkb0 = 0
        do nt = 1, ntyp
           if (upf(nt)%tvanp ) then
              do na = 1, nat
                 if (ityp (na) == nt) then
                    ijh = 1
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       IF (noncolin) THEN
                          DO is1=1,npol
                             DO is2=1,npol
                                becsum1_nc (ijh, na, is1, is2) = &
                                becsum1_nc (ijh, na, is1, is2) + w1 * &
                                 (CONJG(becp%nc(ikb,is1,ibnd))* &
                                        becp%nc(ikb,is2,ibnd))
                             END DO
                          END DO
                       ELSE
                          becsum1 (ijh, na, current_spin) = &
                            becsum1 (ijh, na, current_spin) + w1 * &
                             DBLE (CONJG(becp%k(ikb,ibnd))*becp%k(ikb,ibnd) )
                       ENDIF
                       ijh = ijh + 1
                       do jh = ih + 1, nh (nt)
                          jkb = ijkb0 + jh
                          IF (noncolin) THEN
                             DO is1=1,npol
                                DO is2=1,npol
                                   becsum1_nc(ijh,na,is1,is2) = &
                                      becsum1_nc(ijh,na,is1,is2) + w1* &
                                      (CONJG(becp%nc(ikb,is1,ibnd))* &
                                             becp%nc(jkb,is2,ibnd) )
                                END DO
                             END DO
                          ELSE
                             becsum1 (ijh, na, current_spin) = &
                               becsum1 (ijh, na, current_spin) + w1 * 2.d0 * &
                                DBLE(CONJG(becp%k(ikb,ibnd))*becp%k(jkb,ibnd) )
                          END IF
                          ijh = ijh + 1
                       enddo
                    enddo
                    ijkb0 = ijkb0 + nh (nt)
                 endif
              enddo
           else
              do na = 1, nat
                 if (ityp (na) == nt) ijkb0 = ijkb0 + nh (nt)
              enddo
           endif
        enddo
        dos_ef = dos_ef + weight * wdelta
     enddo

  enddo
  !$acc end data
  if (doublegrid) then
     do is = 1, nspin_mag
        call fft_interpolate (dffts, ldoss (:, is), dfftp, ldos (:, is))
     enddo
  else
     ldos (:,:) = ldoss (:,:)
  endif

  IF (noncolin.and.okvan) THEN
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==nt) THEN
                 IF (upf(nt)%has_so) THEN
                    CALL transform_becsum_so(becsum1_nc,becsum1,na)
                 ELSE
                    CALL transform_becsum_nc(becsum1_nc,becsum1,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF

  call addusldos (ldos, becsum1)
  !
  !    Collects partial sums on k-points from all pools
  !
  call mp_sum ( ldoss, inter_pool_comm )
  call mp_sum ( ldos, inter_pool_comm )
  call mp_sum ( dos_ef, inter_pool_comm )
  call mp_sum ( becsum1, inter_pool_comm )
  !check
  !      check =0.d0
  !      do is=1,nspin_mag
  !         call fwfft('Rho',ldos(:,is),dfftp)
  !         check = check + omega* DBLE(ldos(nl(1),is))
  !         call invfft('Rho',ldos(:,is),dfftp)
  !      end do
  !      WRITE( stdout,*) ' check ', check, dos_ef
  !check
  !
  IF (noncolin) deallocate(becsum1_nc)
  call deallocate_bec_type_acc(becp)
  !$acc exit data detach(evc) delete(nl_d)
  call stop_clock ('localdos_cond')
  return
end subroutine localdos_cond
!-----------------------------------------------------------------------
subroutine incdrhoscf_cond (drhos, weight, ik, dbecsum, dpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation for the conduction states only in the twochem case
  !     . It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE ions_base,            ONLY : nat
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE wvfct,                ONLY : npwx, nbnd, nbnd_cond
  USE uspp_param,           ONLY : nhm
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : ngk,igk_k
  USE qpoint,               ONLY : ikks, ikqs
  USE control_lr,           ONLY : nbnd_occ
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum
  USE fft_helper_subroutines

  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER, INTENT (IN) :: ik
  ! input: the k point
  REAL(DP), INTENT (IN) :: weight
  ! input: the weight of the k point
  COMPLEX(DP), INTENT (IN) :: dpsi (npwx,nbnd)
  ! input: the perturbed wfc for the given k point
  COMPLEX(DP), INTENT (INOUT) :: drhos (dffts%nnr), dbecsum (nhm*(nhm+1)/2,nat)
  ! input/output: the accumulated change to the charge density and dbecsum
  !
  !   here the local variables
  !
  REAL(DP) :: wgt
  ! the effective weight of the k point

  COMPLEX(DP), ALLOCATABLE :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_dpsi(:), tg_drho(:)

  INTEGER :: npw, npwq, ikk, ikq, itmp
  INTEGER :: ibnd, ir, ir3, ig, incr, v_siz, idx, ioff, ioff_tg, nxyp
  INTEGER :: right_inc, ntgrp
  ! counters

  ! For device buffer
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
  !$acc enter data copyin(nl_d, evc)


  CALL start_clock_gpu ('incdrhoscf_cond')
  !
  ALLOCATE(dpsic(dffts%nnr))
  ALLOCATE(psi(dffts%nnr))
  !
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  incr = 1
  !
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz = dffts%nnr_tg
     !
     ALLOCATE( tg_psi( v_siz ) )
     ALLOCATE( tg_dpsi( v_siz ) )
     ALLOCATE( tg_drho( v_siz ) )
     !
     incr = fftx_ntgrp(dffts)
     !
  ELSE
     v_siz = dffts%nnr
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  !$acc data copyin(dpsi(1:npwx,1:nbnd)) copy(drhos(1:v_siz)) create(psi(1:v_siz),dpsic(1:v_siz))&
  !$acc present(igk_k, evc, nl_d)
  do ibnd = 1+(nbnd-nbnd_cond), nbnd_occ(ikk), incr
  !only conduction states
     !
     IF ( dffts%has_task_groups ) THEN
        !
        tg_drho=(0.0_DP, 0.0_DP)
        tg_psi=(0.0_DP, 0.0_DP)
        tg_dpsi=(0.0_DP, 0.0_DP)
        !
        ioff   = 0
        CALL tg_get_recip_inc( dffts, right_inc )
        ntgrp = fftx_ntgrp( dffts )
        !
        DO idx = 1, ntgrp
           !
           ! ... dtgs%nogrp ffts at the same time. We prepare both
           ! evc (at k) and dpsi (at k+q)
           !
           IF( idx + ibnd - 1 <= nbnd_occ(ikk) ) THEN
              !
              DO ig = 1, npw
                 tg_psi( dffts%nl( igk_k( ig,ikk ) ) + ioff ) = evc( ig, idx+ibnd-1 )
              END DO
              DO ig = 1, npwq
                 tg_dpsi( dffts%nl( igk_k( ig,ikq ) ) + ioff ) = dpsi( ig, idx+ibnd-1 )
              END DO
              !
           END IF
           !
           ioff = ioff + right_inc
           !
        END DO
        CALL invfft ('tgWave', tg_psi, dffts)
        CALL invfft ('tgWave', tg_dpsi, dffts)

        do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
           tg_drho (ir) = tg_drho (ir) + wgt * CONJG(tg_psi (ir) ) *  tg_dpsi (ir)
        enddo
        !
        ! reduce the group charge (equivalent to sum over bands of
        ! orbital group)
        !
        CALL tg_reduce_rho( drhos, tg_drho, dffts )
        !
     ELSE
        !
        ! Normal case: no task groups
        !
        ! Initialize psi and dpsic
        !
        !$acc kernels
        psi (:) = (0.d0, 0.d0)
        dpsic(:) = (0.d0, 0.d0)
        !$acc end kernels
        !
        !$acc parallel loop
        do ig = 1, npw
           itmp = nl_d (igk_k(ig,ikk) )
           psi (itmp ) = evc (ig, ibnd)
        enddo
        !$acc parallel loop
        do ig = 1, npwq
           itmp = nl_d (igk_k(ig,ikq) )
           dpsic ( itmp ) = dpsi (ig, ibnd)
        enddo
        !
        ! FFT to R-space of the unperturbed/perturbed wfcts psi/dpsi
        !
        !$acc host_data use_device(psi)
        CALL invfft ('Wave', psi, dffts)
        !$acc end host_data
        !$acc host_data use_device(dpsic)
        CALL invfft ('Wave', dpsic, dffts)
        !$acc end host_data
        !
        ! Calculation of the response charge-density
        !
        !$acc parallel loop
        do ir = 1, v_siz
           drhos (ir) = drhos (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
        enddo
        !
     ENDIF
     !
  enddo ! loop on bands
  !$acc end data
  !
  ! Ultrasoft contribution
  ! Calculate dbecsum = <evc|vkb><vkb|dpsi>
  !
  CALL addusdbec_cond (ik, weight, dpsi, dbecsum)
  !
  DEALLOCATE(psi)
  DEALLOCATE(dpsic)
  !
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  ENDIF
  !
  !$acc exit data detach(evc) delete(nl_d)
  CALL stop_clock_gpu ('incdrhoscf_cond')
  !
  RETURN
  !
end subroutine incdrhoscf_cond
!
!-----------------------------------------------------------------------
subroutine incdrhoscf_cond_nc (drhos, weight, ik, dbecsum, dpsi, rsign)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation for the conduction states only in the twochem case
  !     . It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, domag, nspin_mag
  USE uspp_param,           ONLY : nhm
  USE wvfct,                ONLY : npwx, nbnd, nbnd_cond
  USE wavefunctions, ONLY : evc
  USE klist,                ONLY : ngk,igk_k
  USE qpoint,               ONLY : ikks, ikqs
  USE control_lr,           ONLY : nbnd_occ
  USE qpoint_aux,           ONLY : becpt
  USE lrus,                 ONLY : becp1
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum
  USE fft_helper_subroutines

  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER, INTENT(IN) :: ik
  ! input: the k point
  REAL(DP), INTENT(IN) :: weight
  REAL(DP), INTENT(IN) :: rsign
  ! input: the weight of the k point
  ! the sign in front of the response of the magnetization density
  COMPLEX(DP), INTENT(IN) :: dpsi(npwx*npol,nbnd)
  ! input: the perturbed wfcs at the given k point
  COMPLEX(DP), INTENT(INOUT) :: drhos (dffts%nnr,nspin_mag), dbecsum (nhm,nhm,nat,nspin)
  ! input/output: the accumulated change of the charge density and dbecsum
  !
  !   here the local variable
  !
  REAL(DP) :: wgt
  ! the effective weight of the k point
  !
  COMPLEX(DP), ALLOCATABLE :: psi (:,:), dpsic (:,:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  !
  COMPLEX(DP), ALLOCATABLE :: tg_psi (:,:), tg_dpsi (:,:), tg_drho(:,:)
  !
  INTEGER :: npw, npwq, ikk, ikq, itmp
  INTEGER :: ibnd, jbnd, ir, ir3, ig, incr, v_siz, idx, ioff, ioff_tg, nxyp
  INTEGER :: ntgrp, right_inc
  ! counters
  !
  ! For device buffer
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
  !$acc enter data copyin(nl_d, evc)
  !
  !
  CALL start_clock_gpu ('incdrhoscf_cond')
  !
  ALLOCATE (dpsic(dffts%nnr, npol))
  ALLOCATE (psi  (dffts%nnr, npol))
  !
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  incr = 1
  !
  IF (dffts%has_task_groups) THEN
     !
     v_siz = dffts%nnr_tg
     !
     ALLOCATE( tg_psi( v_siz, npol ) )
     ALLOCATE( tg_dpsi( v_siz, npol ) )
     ALLOCATE( tg_drho( v_siz, nspin_mag ) )
     !
     incr  = fftx_ntgrp(dffts)
     !
  ELSE
     v_siz = dffts%nnr
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  !$acc data copyin(dpsi(1:npwx*npol,1:nbnd)) copy(drhos(1:v_siz,1:nspin_mag)) &
  !$acc create(psi(1:v_siz,1:npol),dpsic(1:v_siz,1:npol)) present(igk_k, evc, nl_d)
  do ibnd = 1+(nbnd-nbnd_cond), nbnd_occ(ikk), incr
  !only conduction states

     IF (dffts%has_task_groups) THEN
#if defined(__CUDA)
        CALL errore( ' incdrhoscf_cond_nc ', ' taskgroup par not implement with GPU offload', 1 )
#endif
        !
        tg_drho=(0.0_DP, 0.0_DP)
        tg_psi=(0.0_DP, 0.0_DP)
        tg_dpsi=(0.0_DP, 0.0_DP)
        !
        ioff   = 0
        CALL tg_get_recip_inc( dffts, right_inc )
        ntgrp = fftx_ntgrp( dffts )
        !
        DO idx = 1, ntgrp
           !
           ! ... dtgs%nogrp ffts at the same time. We prepare both
           ! evc (at k) and dpsi (at k+q)
           !
           IF( idx + ibnd - 1 <= nbnd_occ(ikk) ) THEN
              !
              DO ig = 1, npw
                 tg_psi( dffts%nl( igk_k( ig,ikk ) ) + ioff, 1 ) = evc( ig, idx+ibnd-1 )
                 tg_psi( dffts%nl( igk_k( ig,ikk ) ) + ioff, 2 ) = evc( npwx+ig, idx+ibnd-1 )
              END DO
              DO ig = 1, npwq
                 tg_dpsi( dffts%nl( igk_k( ig,ikq ) ) + ioff, 1 ) = dpsi( ig, idx+ibnd-1 )
                 tg_dpsi( dffts%nl( igk_k( ig,ikq ) ) + ioff, 2 ) = dpsi( npwx+ig, idx+ibnd-1 )
              END DO
              !
           END IF
           !
           ioff = ioff + right_inc
           !
        END DO
        CALL invfft ('tgWave', tg_psi(:,1), dffts)
        CALL invfft ('tgWave', tg_psi(:,2), dffts)
        CALL invfft ('tgWave', tg_dpsi(:,1), dffts)
        CALL invfft ('tgWave', tg_dpsi(:,2), dffts)

        do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
           tg_drho (ir,1) = tg_drho (ir,1) + wgt * (CONJG(tg_psi (ir,1) )*tg_dpsi (ir,1) &
                                                  + CONJG(tg_psi (ir,2) )*tg_dpsi (ir,2) )
        enddo
        IF (domag) THEN
           do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
              tg_drho(ir,2)= tg_drho(ir,2) + (rsign) *wgt * (CONJG(tg_psi(ir,1))*tg_dpsi(ir,2) &
                                                  + CONJG(tg_psi(ir,2))*tg_dpsi(ir,1) )
              tg_drho(ir,3)= tg_drho(ir,3) + (rsign) *wgt * (CONJG(tg_psi(ir,1))*tg_dpsi(ir,2) &
                                                  - CONJG(tg_psi(ir,2))*tg_dpsi(ir,1) ) * (0.d0,-1.d0)
              tg_drho(ir,4)= tg_drho(ir,4) + (rsign) *wgt * (CONJG(tg_psi(ir,1))*tg_dpsi(ir,1) &
                                                  - CONJG(tg_psi(ir,2))*tg_dpsi(ir,2) )
           enddo
        ENDIF
        !
        ! reduce the group charge (equivalent to sum over the bands of the
        ! orbital group)
        !
        CALL tg_reduce_rho( drhos, tg_drho, dffts )
        !
     ELSE
        !
        ! Normal case: no task groups
        !
        ! Initialize psi and dpsic
        !
        !$acc kernels
        psi (:,:) = (0.d0, 0.d0)
        dpsic (:,:) = (0.d0, 0.d0)
        !$acc end kernels
        !
        !$acc parallel loop
        do ig = 1, npw
           itmp = nl_d ( igk_k(ig,ikk) )
           psi (itmp, 1) = evc (ig, ibnd)
           psi (itmp, 2) = evc (ig+npwx, ibnd)
        enddo
        !$acc parallel loop
        do ig = 1, npwq
           itmp = nl_d (igk_k(ig,ikq))
           dpsic (itmp, 1 ) = dpsi (ig, ibnd)
           dpsic (itmp, 2 ) = dpsi (ig+npwx, ibnd)
        enddo
        !
        ! FFT to R-space of the unperturbed/perturbed wfcts psi/dpsi
        !
        !$acc host_data use_device(psi)
        CALL invfft ('Wave', psi(:,1), dffts)
        CALL invfft ('Wave', psi(:,2), dffts)
        !$acc end host_data
        !$acc host_data use_device(dpsic)
        CALL invfft ('Wave', dpsic(:,1), dffts)
        CALL invfft ('Wave', dpsic(:,2), dffts)
        !$acc end host_data
        !
        ! Calculation of the response charge density
        !
        !$acc parallel loop
        do ir = 1, v_siz
           drhos(ir,1)=drhos(ir,1)+wgt*(CONJG(psi(ir,1))*dpsic(ir,1)  +  &
                                            CONJG(psi(ir,2))*dpsic(ir,2) )
        enddo
        IF (domag) THEN
           !$acc parallel loop
           do ir = 1, v_siz
              drhos(ir,2)=drhos (ir,2) + (rsign) *wgt * (CONJG(psi(ir,1))*dpsic(ir,2) &
                                                  + CONJG(psi(ir,2))*dpsic(ir,1) )
              drhos(ir,3)=drhos (ir,3) + (rsign) *wgt * (CONJG(psi(ir,1))*dpsic(ir,2) &
                                                  - CONJG(psi(ir,2))*dpsic(ir,1) ) * (0.d0,-1.d0)
              drhos(ir,4)=drhos (ir,4) + (rsign) *wgt * (CONJG(psi(ir,1))*dpsic(ir,1) &
                                                  - CONJG(psi(ir,2))*dpsic(ir,2) )
           enddo
        END IF
        !
     END IF
     !
  enddo
  !$acc end data
  !
  ! Ultrasoft contribution
  ! Calculate dbecsum_nc = <evc|vkb><vkb|dpsi>
  !
  IF (rsign==1.0d0) THEN
     CALL addusdbec_nc (ik, weight, dpsi, dbecsum, becp1)
  ELSE
     CALL addusdbec_nc (ik, weight, dpsi, dbecsum, becpt)
  ENDIF
  !
  DEALLOCATE(psi)
  DEALLOCATE(dpsic)
  !
  IF (dffts%has_task_groups) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  END IF
  !
  !$acc exit data detach(evc) delete(nl_d)
  CALL stop_clock_gpu ('incdrhoscf_cond')
  !
  RETURN
  !
end subroutine incdrhoscf_cond_nc
!----------------------------------------------------------------------
subroutine addusdbec_cond (ik, wgt, dpsi, dbecsum)
  !----------------------------------------------------------------------
  !
  !  This routine adds to dbecsum the contribution of this
  !  k point for the conduction states only. It implements Eq. B15 of PRB 64, 235118 (2001).
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE becmod,     ONLY : calbec
  USE wvfct,      ONLY : npwx, nbnd,nbnd_cond
  USE uspp,       ONLY : nkb, vkb, okvan, ijtoh
  USE uspp_param, ONLY : upf, nh, nhm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE klist,      ONLY : ngk
  USE lrus,       ONLY : becp1
  USE qpoint,     ONLY : ikks, ikqs
  USE control_lr, ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  !   the dummy variables
  !
  COMPLEX(DP) :: dbecsum (nhm*(nhm+1)/2, nat), dpsi(npwx,nbnd)
  ! inp/out: the sum kv of bec *
  ! input  : contains delta psi
  INTEGER :: ik
  ! input: the k point
  REAL(DP) :: wgt
  ! input: the weight of this k point
  !
  !     here the local variables
  !
  INTEGER :: na, nt, ih, jh, ibnd, ikb, jkb, ijh, startb, &
       lastb, ijkb0
  ! counter on atoms
  ! counter on atomic type
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on the bands
  ! the real k point
  ! counter on solid becp
  ! counter on solid becp
  ! composite index for dbecsum
  ! divide among processors the sum
  ! auxiliary variable for counting
  INTEGER :: ikk, ikq, npwq
  ! index of the point k
  ! index of the point k+q
  ! number of the plane-waves at point k+q

  COMPLEX(DP), ALLOCATABLE :: dbecq (:,:)
  ! the change of becq

  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('addusdbec_cond')
  !
  ALLOCATE (dbecq( nkb, nbnd))
  !
  ikk  = ikks(ik)
  ikq  = ikqs(ik)
  npwq = ngk(ikq)
  !
  ! First compute the product of dpsi and vkb
  !
  CALL calbec (npwq, vkb, dpsi, dbecq)
  !
  !  And then we add the product to becsum
  !
  !  Band parallelization: each processor takes care of its slice of bands
  !
  CALL divide (intra_bgrp_comm, nbnd_occ (ikk)-(nbnd-nbnd_cond), startb, lastb)
  !conduction states only
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              !
              !  And qgmq and becp and dbecq
              !
              do ih = 1, nh(nt)
                 ikb = ijkb0 + ih
                 ijh=ijtoh(ih,ih,nt)
                 do ibnd = startb, lastb
                 !conduction states only
                    dbecsum (ijh, na) = dbecsum (ijh, na) + &
                         wgt * ( CONJG(becp1(ik)%k(ikb,ibnd)) * dbecq(ikb,ibnd) )
                 enddo
                 do jh = ih + 1, nh (nt)
                    ijh=ijtoh(ih,jh,nt)
                    jkb = ijkb0 + jh
                    do ibnd = startb, lastb
                       dbecsum (ijh, na) = dbecsum (ijh, na) + &
                         wgt*( CONJG(becp1(ik)%k(ikb,ibnd))*dbecq(jkb,ibnd) + &
                               CONJG(becp1(ik)%k(jkb,ibnd))*dbecq(ikb,ibnd) )
                    enddo
                    ijh = ijh + 1
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
        enddo
     endif
  enddo
  !
  DEALLOCATE (dbecq)
  !
  CALL stop_clock ('addusdbec_cond')
  !
  RETURN
  !
end subroutine addusdbec_cond
!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine addusdbec_cond_nc (ik, wgt, dpsi, dbecsum_nc, becp1)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the dbecsum the term which correspond to this
  !  k point. After the accumulation the additional part of the charge
  !  is computed in addusddens.
  !
  USE kinds,            ONLY : DP
  USE lsda_mod,         ONLY : nspin
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE becmod,           ONLY : calbec, bec_type
  USE wvfct,            ONLY : npwx, nbnd,nbnd_cond
  USE uspp,             ONLY : nkb, vkb, okvan
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_param,       ONLY : upf, nh, nhm
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE klist,            ONLY : ngk
  USE qpoint,           ONLY : nksq, ikks, ikqs
  USE control_lr,       ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  !   the dummy variables
  !
  COMPLEX(DP) :: dbecsum_nc (nhm,nhm,nat,nspin), dpsi(npwx*npol,nbnd)
  ! inp/out: the sum kv of bec *
  ! input  : contains delta psi
  TYPE(bec_type) :: becp1(nksq)
  !
  INTEGER :: ik
  ! input: the k point
  REAL(DP) :: wgt
  ! input: the weight of this k point
  !
  !     here the local variables
  !
  INTEGER :: na, nt, ih, jh, ibnd, ikb, jkb, startb, &
       lastb, ijkb0, is1, is2, ijs
  ! counter on atoms
  ! counter on atomic type
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on the bands
  ! the real k point
  ! counter on solid becp
  ! counter on solid becp
  ! composite index for dbecsum
  ! divide among processors the sum
  ! auxiliary variable for counting
  INTEGER :: ikk, ikq, npwq
  ! index of the point k
  ! index of the point k+q
  ! number of the plane-waves at point k+q

  COMPLEX(DP), ALLOCATABLE :: dbecq_nc(:,:,:)
  ! the change of becq

  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('addusdbec_cond_nc')
  !
  ALLOCATE (dbecq_nc( nkb,npol, nbnd))
  !
  ikk  = ikks(ik)
  ikq  = ikqs(ik)
  npwq = ngk(ikq)
  !
  ! First compute the product of dpsi and vkb
  !
  CALL calbec (npwq, vkb, dpsi, dbecq_nc)
  !
  !  And then we add the product to becsum
  !
  !  Band parallelization: each processor takes care of its slice of bands
  !
  CALL divide (intra_bgrp_comm, nbnd_occ (ikk)-(nbnd-nbnd_cond), startb, lastb)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd = startb, lastb
                       ijs=0
                       DO is1=1,npol
                          DO is2=1,npol
                             ijs=ijs+1
                             dbecsum_nc(ih,jh,na,ijs)=dbecsum_nc(ih,jh,na,ijs)+&
                                wgt*CONJG(becp1(ik)%nc(ikb,is1,ibnd))          &
                                        *dbecq_nc(jkb,is2,ibnd)
                          ENDDO
                       ENDDO
                    ENDDO
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
        enddo
     endif
  enddo
  !
  DEALLOCATE (dbecq_nc)
  !
  CALL stop_clock ('addusdbec_cond_nc')
  !
  RETURN
  !
end subroutine addusdbec_cond_nc
!
SUBROUTINE sternheimer_kernel_twochem(first_iter, time_reversed, npert, lrdvpsi, iudvpsi, &
         thresh, dvscfins, all_conv, avg_iter, drhoout, dbecsum, dbecsum_nc, &
         drhoout_cond,dbecsum_cond,dbecsum_cond_nc,exclude_hubbard)
   !----------------------------------------------------------------------------
   !This is a copy of the sternheimer kernel to be used when twochem and lmetq0.
   !In addition to the usual density response, it also calculates the conduction bands only
   !density response, needed to determine the conduction and valence manifold chemical potential
   !shift that can be present in the q=0 case
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE io_global,             ONLY : stdout
   USE mp,                    ONLY : mp_sum
   USE mp_pools,              ONLY : inter_pool_comm
   USE buffers,               ONLY : get_buffer, save_buffer
   USE fft_base,              ONLY : dffts
   USE ions_base,             ONLY : nat
   USE klist,                 ONLY : xk, wk, ngk, igk_k
   USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                 ONLY : nbnd, npwx, et
   USE wavefunctions,         ONLY : evc
   USE noncollin_module,      ONLY : noncolin, domag, npol, nspin_mag
   USE uspp,                  ONLY : vkb
   USE uspp_param,            ONLY : nhm
   USE uspp_init,             ONLY : init_us_2
   USE ldaU,                  ONLY : lda_plus_u
   USE units_lr,              ONLY : iuwfc, lrwfc, lrdwf, iudwf
   USE control_lr,            ONLY : nbnd_occ, lgamma
   USE qpoint,                ONLY : nksq, ikks, ikqs
   USE qpoint_aux,            ONLY : ikmks, ikmkmqs, becpt
   USE eqv,                   ONLY : dpsi, dvpsi, evq
   USE apply_dpot_mod,        ONLY : apply_dpot_bands
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! true if the first iteration.
   LOGICAL, INTENT(IN) :: time_reversed
   !! true if solving for time reversed wave functions
   LOGICAL, INTENT(IN), OPTIONAL :: exclude_hubbard
   !! true if ignoring the Hubbard response term
   INTEGER, INTENT(IN) :: npert
   !! number of perturbations
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(IN) :: thresh
   !! threshold for solving the linear equation
   LOGICAL, INTENT(OUT) :: all_conv
   !! True if converged at all k points and perturbations
   REAL(DP), INTENT(OUT) :: avg_iter
   !! average number of iterations for the linear equation solver
   COMPLEX(DP), POINTER, INTENT(IN) :: dvscfins(:, :, :)
   !! dV_ind calculated in the previous iteration
   COMPLEX(DP), INTENT(INOUT) :: drhoout(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: drhoout_cond(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT) :: dbecsum_cond(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_cond_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   !
   LOGICAL :: conv_root
   !! true if linear system is converged
   LOGICAL :: exclude_hubbard_
   !! Local variable to set the default of exclude_hubbard to false
   INTEGER :: ikk, ikq, npw, npwq, ipert, num_iter, ik, nrec, ikmk, ikmkmq
   !! counters
   INTEGER :: tot_num_iter
   !! total number of iterations in cgsolve_all
   INTEGER :: tot_cg_calls
   !! total number of cgsolve_all calls
   REAL(DP) :: anorm
   !! the norm of the error of the conjugate gradient solution
   REAL(DP) :: rsign
   !! sign of the term in the magnetization
   REAL(DP), ALLOCATABLE :: h_diag(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning
   COMPLEX(DP) , ALLOCATABLE :: aux2(:, :)
   !! temporary storage used in apply_dpot_bands
   !
   EXTERNAL ch_psi_all, cg_psi
   !! functions passed to cgsolve_all
   !
   ! Initialization
   !
   CALL start_clock("sth_kernel")
   !
   exclude_hubbard_ = .FALSE.
   IF (PRESENT(exclude_hubbard)) exclude_hubbard_ = exclude_hubbard
   !
   ALLOCATE(h_diag(npwx*npol, nbnd))
   ALLOCATE(aux2(npwx*npol, nbnd))
   !
   !$acc enter data create(aux2(1:npwx*npol, 1:nbnd))
   !
   all_conv = .TRUE.
   tot_num_iter = 0
   tot_cg_calls = 0
   !
   DO ik = 1, nksq
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      !
      ! Set time-reversed k and k+q points
      !
      IF (time_reversed) THEN
         ikmk = ikmks(ik)
         ikmkmq = ikmkmqs(ik)
         rsign = -1.0_DP
      ELSE
         ikmk = ikk
         ikmkmq = ikq
         rsign = 1.0_DP
      ENDIF
      !
      IF (lsda) current_spin = isk(ikk)
      !
      ! reads unperturbed wavefunctions psi_k in G_space, for all bands
      ! if q=0, evq is a pointer to evc
      !
      IF (nksq > 1 .OR. (noncolin .AND. domag)) THEN
         IF (lgamma) THEN
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
         ELSE
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            CALL get_buffer(evq, lrwfc, iuwfc, ikmkmq)
         ENDIF
      ENDIF
      !
      ! compute beta functions and kinetic energy for k-point ik
      ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
      !
      CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
      !$acc update host(vkb)
      CALL g2_kin(ikq)
      !
      ! compute preconditioning matrix h_diag used by cgsolve_all
      !
      CALL h_prec(ik, evq, h_diag)
      !
      DO ipert = 1, npert
         !
         ! read P_c^+ x psi_kpoint into dvpsi.
         !
         nrec = (ipert - 1) * nksq + ik
         IF (time_reversed) nrec = nrec + npert * nksq
         !
         CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
         !
         IF (.NOT. first_iter) THEN
            !
            ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
            ! dvscf_q from previous iteration (mix_potential)
            !
            CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvscfins(:, :, ipert), evc, aux2)
            dvpsi = dvpsi + aux2
            !
            !  In the case of US pseudopotentials there is an additional
            !  selfconsist term which comes from the dependence of D on
            !  V_{eff} on the bare change of the potential
            !
            IF (time_reversed) THEN
               CALL adddvscf_ph_mag(ipert, ik)
            ELSE
               CALL adddvscf(ipert, ik)
            ENDIF
            !
            ! DFPT+U: add to dvpsi the scf part of the response
            ! Hubbard potential dV_hub
            !
            IF (lda_plus_u .AND. (.NOT. exclude_hubbard_)) CALL adddvhubscf(ipert, ik)
            !
         ENDIF
         !
         ! Orthogonalize dvpsi to valence states
         !
         CALL orthogonalize(dvpsi, evq, ikmk, ikmkmq, dpsi, npwq, .FALSE.)
         !
         ! Initial guess for dpsi
         !
         IF (first_iter) THEN
            !
            !  At the first iteration dpsi is set to zero
            !
            dpsi(:, :) = (0.d0,0.d0)
         ELSE
            !
            ! starting value for delta_psi is read from iudwf
            !
            CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
         ENDIF
         !
         ! iterative solution of the linear system (H-e)*dpsi=dvpsi
         ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
         !
         conv_root = .TRUE.
         !
         ! TODO: should nbnd_occ(ikk) be nbnd_occ(ikmk)?
         CALL cgsolve_all(ch_psi_all, cg_psi, et(1, ikmk), dvpsi, dpsi, h_diag, &
            npwx, npwq, thresh, ik, num_iter, conv_root, anorm, nbnd_occ(ikk), npol)
         !
         tot_num_iter = tot_num_iter + num_iter
         tot_cg_calls = tot_cg_calls + 1
         !
         IF (.NOT. conv_root) THEN
            all_conv = .FALSE.
            WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
               &root not converged, thresh < ', es10.3)") ik, anorm
         ENDIF
         !
         ! writes delta_psi on iunit iudwf, k=kpoint,
         !
         CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
         !
         ! calculates dvscf, sum over k => dvscf_q_ipert
         !
         IF (noncolin) THEN
            CALL incdrhoscf_nc(drhoout(1,1,ipert), wk(ikk), ik, &
                               dbecsum_nc(1,1,1,1,ipert), dpsi, rsign)
            CALL incdrhoscf_cond_nc(drhoout_cond(1,1,ipert), wk(ikk), ik, &
                               dbecsum_cond_nc(1,1,1,1,ipert), dpsi, rsign)
         ELSE
            CALL incdrhoscf(drhoout(1,current_spin,ipert), wk(ikk), &
                            ik, dbecsum(1,1,current_spin,ipert), dpsi)
            CALL incdrhoscf_cond(drhoout_cond(1,current_spin,ipert), wk(ikk), &
                            ik, dbecsum_cond(1,1,current_spin,ipert), dpsi)
         ENDIF
      ENDDO ! ipert
   ENDDO ! ik
   !
   CALL mp_sum(tot_num_iter, inter_pool_comm)
   CALL mp_sum(tot_cg_calls, inter_pool_comm)
   avg_iter = REAL(tot_num_iter, DP) / REAL(tot_cg_calls, DP)
   !
   !$acc exit data delete(aux2)
   !
   DEALLOCATE(aux2)
   DEALLOCATE(h_diag)
   !
   CALL stop_clock("sth_kernel")
   !
!----------------------------------------------------------------------------
END SUBROUTINE sternheimer_kernel_twochem
!
SUBROUTINE allocate_twochem(npe, nsolv)
   USE fft_base,             ONLY : dfftp, dffts
   USE ions_base,            ONLY : nat
   USE uspp_param,           ONLY : nhm
   USE lsda_mod,             ONLY : nspin
   USE paw_variables,        ONLY : okpaw
   USE noncollin_module,     ONLY : noncolin, nspin_mag
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: npe, nsolv
   !
   IF (noncolin) allocate (dbecsum_cond_nc (nhm,nhm, nat , nspin , npe, nsolv))
   allocate (drhos_cond ( dffts%nnr, nspin_mag , npe))
   allocate (drhop_cond ( dfftp%nnr, nspin_mag , npe))
   allocate (dbecsum_cond ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
   allocate ( ldos_cond( dfftp%nnr  , nspin_mag) )
   allocate ( ldoss_cond( dffts%nnr , nspin_mag) )
   allocate (becsum1_cond ( (nhm * (nhm + 1))/2 , nat , nspin_mag))
   call localdos_cond ( ldos_cond , ldoss_cond , becsum1_cond, dos_ef_cond )
   IF (.NOT.okpaw) deallocate(becsum1_cond)
END SUBROUTINE allocate_twochem
!
SUBROUTINE deallocate_twochem
   USE paw_variables,        ONLY : okpaw
   USE noncollin_module,     ONLY : noncolin
   !
   IMPLICIT NONE
   !
   !deallocate for twochem calculation at gamma
   if (allocated(ldoss_cond)) deallocate (ldoss_cond)
   if (allocated(ldos_cond)) deallocate (ldos_cond)
   deallocate (dbecsum_cond)
   IF (okpaw) THEN
            if (allocated(becsum1_cond)) deallocate (becsum1_cond)
   ENDIF
   IF (noncolin) deallocate (dbecsum_cond_nc)
   deallocate (drhop_cond)
   deallocate (drhos_cond)
END SUBROUTINE deallocate_twochem
!
END MODULE lr_two_chem
