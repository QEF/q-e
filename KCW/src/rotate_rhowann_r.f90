!
!Giovanni Cistaro
!
!Based on the function from pw2wannier. WARNING! the rotated wfs is only the periodic part,
! not the full wf as it may seem from the header.
   SUBROUTINE rotate_rhowann_r(isym, psi, gpsi)
      !-----------------------------------------------------------------------
      ! g u_k = u_k(rS) WARNING it might need a factor with fractional translation 
      !-----------------------------------------------------------------------
      USE kinds,           ONLY : DP
      USE wvfct,           ONLY : npwx
      USE wavefunctions,   ONLY : evc, psic, psic_nc
      USE fft_base,        ONLY : dffts, dfftp
      USE fft_interfaces,  ONLY : fwfft, invfft
      USE cell_base,       ONLY : bg
      USE constants,       ONLY : tpi
      USE gvect,           ONLY : g, ngm
      USE klist,           ONLY : igk_k, ngk
      USE mp,              ONLY : mp_sum
      USE mp_pools,        ONLY : intra_pool_comm
      USE fft_interfaces,  ONLY : invfft
      USE scatter_mod,     ONLY : gather_grid, scatter_grid
      USE control_kcw,     ONLY : num_wann, rir
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN):: isym
      !(Non colin) COMPLEX(DP), INTENT(IN):: psi(npwx*npol,nbnd)
      !(Non colin) COMPLEX(DP), INTENT(OUT):: gpsi(npwx*npol,nbnd)
      COMPLEX(DP), INTENT(IN):: psi(dffts%nnr)
      COMPLEX(DP)            :: gpsi(dffts%nnr)
      !
      !INTEGER:: ig, igk_local, igk_global, npw1, npw2, n, nxxs, ipol, istart, isym0
      !REAL(DP)                 :: kdiff_cart(3), srt(3,3)
      !REAL(DP)                 :: phase_arg
      !COMPLEX(DP)              :: phase_factor, u_spin(2,2), u_spin2(2,2)
      !COMPLEX(DP), ALLOCATABLE :: phase(:), gpsi_tmp(:,:)
      COMPLEX(DP), ALLOCATABLE :: temppsi_all(:)
      COMPLEX(DP), ALLOCATABLE :: psi_all(:)
      INTEGER :: nxxs
      nxxs = dffts%nr1x *dffts%nr2x *dffts%nr3x !Giovanni Cistaro ->  global number of r points
      !IF(.not. ALLOCATED(gpsi)) ALLOCATE(gpsi(dffts%nnr))
      ALLOCATE( psi_all(nxxs), temppsi_all(nxxs))!, gpsi_tmp(npwx,2))!, stat=ierr)
      !IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_all/temppsic_all/gpsi_tmp', 1)
      !ALLOCATE( phase(dffts%nnr))!, stat=ierr)
      !IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase', 1)
      !
      ! for spin space rotation (only used for noncollinear case)
      !
      !Noncolin -> srt = transpose(sr2(:,:,isym))
      !Noncolin -> CALL find_u(srt, u_spin)
      !
      gpsi = CMPLX(0.0, 0.0, kind=DP)
      !
      !npw1 = ngk(ik1)
      !npw2 = ngk(ik2)
      !DO n=1, nbnd
         !
         ! real space rotation
         !
         !gpsi_tmp = 0
         !DO ipol=1,npol
            !istart = npwx*(ipol-1)
            !psic(:) = (0.d0, 0.d0)
            !psic(dffts%nl (igk_k (1:npw1,ik1) ) ) = psi (istart+1:istart+npw1, n)
            !
            !IF (igk_global > 0 .or. isym > 1) THEN
               !CALL invfft ('Wave', psic, dffts)
               IF (isym > 1) THEN
#if defined(__MPI)
                 ! gather among all the CPUs
                 CALL gather_grid(dffts, psi, temppsi_all) 
                 ! gathers all the piecies of psi to temppsi_all
                 ! (Giovanni Cistaro)
                 !
                 ! apply rotation
                 psi_all(1:nxxs) = temppsi_all(rir(1:(nxxs),isym)) 
                 !rotates the space of temppsic_all and builds 
                 !psic_all with it Giovanni Cistaro)
                 !
                 ! scatter back a piece to each CPU
                 CALL scatter_grid(dffts, psi_all, gpsi)
#else
                 gpsi(1:nnxs) = psi(rir(1:(nxxs),isym))
#endif
               ELSE
                 gpsi = psi 
               END IF
               ! Non colin -> IF(t_rev2(isym) == 1) psic = conjg(psic)
               ! apply phase e^{-iGr}
               !IF(igk_global > 0) psic(1:dffts%nnr) = psic(1:dffts%nnr) * conjg(phase(1:dffts%nnr))
            !END IF
            !
         !END DO
         !(Non collinear)
         ! spin space rotation
         !
         !DO ipol=1,npol
         !   istart = npwx*(ipol-1)
         !   IF (noncolin) THEN
         !      u_spin2 = u_spin
         !      IF (t_rev2(isym) == 1) u_spin2 = matmul(t_rev_spin, conjg(u_spin))
         !      gpsi(istart+1:istart+npw2,n) = matmul(gpsi_tmp(1:npw2,:), u_spin2(ipol,:))
         !   ELSE
         !      gpsi(istart+1:istart+npw2,n) = gpsi_tmp(1:npw2,ipol)
         !   END IF
         !END DO
         !gpsi(istart+1:istart+npw2,n) = gpsi_tmp(1:npw2,ipol)
      !END DO
      !
      !phase_arg = -tpi * dot_product(xkc(:,ik1), ft2(:,isym))
      !IF (t_rev2(isym) == 1) phase_arg = -phase_arg
      !phase_factor = CMPLX(COS(phase_arg), SIN(phase_arg), KIND=DP)
      !gpsi = gpsi * phase_factor
      !
      !DEALLOCATE( phase, psic_all, temppsic_all, gpsi_tmp )
      DEALLOCATE( psi_all, temppsi_all )
   END SUBROUTINE rotate_rhowann_r



