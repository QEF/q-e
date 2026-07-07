!
! Copyright (C) 2001-2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE localdos_mod
  !
CONTAINS
  !-----------------------------------------------------------------------
  subroutine localdos_new (ldos_data, ef, firstband )
  !-----------------------------------------------------------------------
  !
  !!    This routine computes the local and total density of states
  !!    at Fermi energy ef (maybe it should be contained in ldos_data?)
  !!    If a two-chemical-potential calculations is required,
  !!    specify variable "firstband" wit self-explanatory meaning
  !!
  !!    Workspace psic must be allocated before this routine is called
  !
  !
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : omega
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : invfft, fft_interpolate
  USE buffers,          ONLY : get_buffer
  USE gvecs,            ONLY : doublegrid
  USE klist,            ONLY : xk, wk, ngk, igk_k, degauss, ngauss, ltetra
  USE lsda_mod,         ONLY : nspin, lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wvfct,            ONLY : nbnd, npwx, et
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
  USE dfpt_type,        ONLY : dfpt_ldos_type
  USE paw_variables,    ONLY : okpaw

  implicit none

  TYPE(dfpt_ldos_type), INTENT(INOUT) :: ldos_data
  !! Local density of states at Ef for conduction states
  !! Contains: dos_ef, ldos, ldoss, becsum_dos
  REAL(DP), INTENT(IN) :: Ef
  !! Fermi energy - maybe it should be copied into ldos_data?
  !! For the twochem case: Fermi energy of the conduction bands
  
  INTEGER, INTENT(IN), OPTIONAL :: firstband
  !! For the twochem case: first conduction band
  !! (the code loops only over conduction bands)
  !
  !    local variables for Ultrasoft PP's
  !
  integer :: nb1, ikb, jkb, ijkb0, ih, jh, na, ijh, nt
  ! counters
  complex(DP), allocatable :: becsum1_nc(:,:,:,:)
  TYPE(bec_type) :: becp
  !
  ! local variables
  !
  real(DP) :: weight, w1, wdelta, dos_ef
  ! weights, DOS(E_F) 
  real(DP), external :: w0gauss
  !
  integer :: npw, ik, ikk, is, ig, ibnd, j, is1, is2, v_siz
  ! counters
  integer :: ios
  ! status flag for i/o
  !
  !  initialize ldos and dos_ef
  !
  v_siz = dffts%nnr
  dos_ef = 0.0_dp
  !
  nb1 = 1
  IF ( PRESENT(firstband) ) THEN
     nb1 = firstband
     IF (nb1 /= 1 .AND. okvan ) &
          CALL errore('localdos','twochem+USPP not working',1)
  END IF
  call start_clock ('localdos_cond')
  IF (noncolin) THEN
     allocate (becsum1_nc( (nhm * (nhm + 1)) / 2, nat, npol, npol))
     becsum1_nc=(0.d0,0.d0)
  ENDIF

  call allocate_bec_type_acc (nkb, nbnd, becp)

  ldos_data%becsum_dos (:,:,:) = 0.d0
  ldos_data%ldos (:,:) = (0d0, 0.0d0)
  ldos_data%ldoss(:,:) = (0d0, 0.0d0)
  ldos_data%dos_ef = 0.0_dp
  !
  ! NB: auxiliary dos_ef variable is used to compute DOS(E_F) on CPU
  ! NB: instead of ldos_data%dos_ef because !$acc copy(ldos_dat) 
  ! NB: overwrites it with GPU value at the end of the $acc region
  !
  !  loop over kpoints
  !
  !$acc data create(psic, psic_nc) copy(ldos_data,ldos_data%ldoss) &
  !$acc      present(vkb, evc, becp)
  do ik = 1, nksq
     ikk = ikks(ik)
     if (lsda) current_spin = isk (ikk)
     npw = ngk(ikk)
     weight = wk (ikk)
     !
     ! unperturbed wfs in reciprocal space read from unit iuwfc
     !
     if (nksq > 1) then
        call get_buffer (evc, lrwfc, iuwfc, ikk)
        !$acc update device(evc)
     endif
     !
     call init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb, .true.)
     call calbec ( offload_type, npw, vkb, evc, becp)
     !
     do ibnd = nb1, nbnd_occ (ikk)
        !
        if(ltetra) then
           wdelta = dfpt_tetra_delta(ibnd,ikk)
        else
           wdelta = w0gauss ( (ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
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
           !$acc parallel loop present(igk_k, psic_nc, dffts, dffts%nl)
           do ig = 1, npw
              psic_nc (dffts%nl (igk_k(ig,ikk)), 1 ) = evc (ig, ibnd)
              psic_nc (dffts%nl (igk_k(ig,ikk)), 2 ) = evc (ig+npwx, ibnd)
           enddo
           !$acc end parallel loop
           !$acc host_data use_device(psic_nc)
           CALL invfft ('Rho', psic_nc(:,1), dffts)
           CALL invfft ('Rho', psic_nc(:,2), dffts)
           !$acc end host_data
           !$acc parallel loop present(psic_nc, ldos_data, ldos_data%ldoss)
           do j = 1, v_siz
              ldos_data%ldoss (j, 1) = ldos_data%ldoss (j, 1) + &
                    w1 * ( DBLE(psic_nc(j,1))**2+AIMAG(psic_nc(j,1))**2 + &
                           DBLE(psic_nc(j,2))**2+AIMAG(psic_nc(j,2))**2)
           enddo
           !$acc end parallel loop
           IF (nspin_mag==4) THEN
              !$acc parallel loop present(psic_nc, ldos_data, ldos_data%ldoss)
              DO j = 1, v_siz
              !
                 ldos_data%ldoss(j,2) = ldos_data%ldoss(j,2) + w1*2.0_DP* &
                             (DBLE(psic_nc(j,1))* DBLE(psic_nc(j,2)) + &
                             AIMAG(psic_nc(j,1))*AIMAG(psic_nc(j,2)))

                 ldos_data%ldoss(j,3) = ldos_data%ldoss(j,3) + w1*2.0_DP* &
                             (DBLE(psic_nc(j,1))*AIMAG(psic_nc(j,2)) - &
                              DBLE(psic_nc(j,2))*AIMAG(psic_nc(j,1)))

                 ldos_data%ldoss(j,4) = ldos_data%ldoss(j,4) + w1* &
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
           !$acc parallel loop present(psic,igk_k,dffts,dffts%nl)
           do ig = 1, npw
              psic (dffts%nl (igk_k(ig,ikk) ) ) = evc (ig, ibnd)
           enddo
           !$acc end parallel loop
           !$acc host_data use_device(psic)
           CALL invfft ('Rho', psic, dffts)
           !$acc end host_data
           !$acc parallel loop present(ldos_data, ldos_data%ldoss, psic)
           do j = 1, v_siz
              ldos_data%ldoss (j, current_spin) = ldos_data%ldoss (j, current_spin) + &
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
                          ldos_data%becsum_dos (ijh, na, current_spin) = &
                            ldos_data%becsum_dos (ijh, na, current_spin) + w1 * &
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
                             ldos_data%becsum_dos (ijh, na, current_spin) = &
                               ldos_data%becsum_dos (ijh, na, current_spin) + w1 * 2.d0 * &
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
        call fft_interpolate (dffts, ldos_data%ldoss (:, is), dfftp, ldos_data%ldos (:, is))
     enddo
  else
     ldos_data%ldos (:,:) = ldos_data%ldoss (:,:)
  endif

  IF (noncolin.and.okvan) THEN
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==nt) THEN
                 IF (upf(nt)%has_so) THEN
                    CALL transform_becsum_so(becsum1_nc,ldos_data%becsum_dos,na)
                 ELSE
                    CALL transform_becsum_nc(becsum1_nc,ldos_data%becsum_dos,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF

  call addusldos (ldos_data%ldos, ldos_data%becsum_dos)
  !
  !    Collects partial sums on k-points from all pools
  !
  call mp_sum ( ldos_data%ldoss, inter_pool_comm )
  call mp_sum ( ldos_data%ldos, inter_pool_comm )
  call mp_sum ( ldos_data%becsum_dos, inter_pool_comm )
  call mp_sum ( dos_ef, inter_pool_comm )
  ldos_data%dos_ef = dos_ef
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
  !
  ! For non-PAW calculations, becsum_dos is not needed anymore so deallocate it.
  ! In the non-PAW but USPP case, becsum_dos still needs to be allocated and 
  ! computed in ! localdos since it is used to update ldos.
  ! In the NCPP case, becsum_dos is not needed at all but is being computed.
  ! This can be ! optimized in the future.
  !
  IF (.NOT. okpaw) DEALLOCATE(ldos_data%becsum_dos)
  !
  call stop_clock ('localdos_cond')
  return
end subroutine localdos_new
!
END MODULE localdos_mod
