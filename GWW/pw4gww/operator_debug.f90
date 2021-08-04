SUBROUTINE operator_debug(npw, e, x, u)

   USE wannier_gw
   USE gvect
   USE becmod,               ONLY : becp,allocate_bec_type,deallocate_bec_type
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE g_psi_mod,            ONLY : h_diag, s_diag
   USE gvecs,                ONLY : doublegrid
   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,                ONLY : DP 
   USE klist,                ONLY : xk,igk_k
   USE mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world,             ONLY : world_comm, mpime, nproc
   USE uspp,                 ONLY : vkb, nkb, okvan
   USE wvfct,                ONLY : g2kin, npwx, nbnd, et

   IMPLICIT NONE

! Dummy Variables 
   REAL(KIND=DP) :: e ! eigenvalue (needed for call not used)
   INTEGER :: npw
   COMPLEX(KIND=DP) :: x(npw), u(npw,1) ! Upon call should contain v_c|x>

   WRITE(*,*) 'Inside operator'
   u(1:npw,1) = x(1:npw) 
   WRITE(*,*) 'Leaving operator'
   RETURN

END SUBROUTINE operator_debug
