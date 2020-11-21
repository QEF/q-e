!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Joshua Elliott JDE
!
! This subroutine computes the operator (1-vP) for overlapping 
! wannier functions
! 



SUBROUTINE operator_1_vp(npw, e, x, u)

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
   USE wavefunctions,        ONLY : evc

   IMPLICIT NONE

! Dummy Variables 
   INTEGER :: npw
   REAL(KIND=DP) :: e ! eigenvalue (needed for call not used)
   COMPLEX(KIND=DP) :: x(npw), u(npw,1) ! Upon call should contain v_c|x>

! variables for computing bare Coulomb interaction
   REAL(KIND=DP), ALLOCATABLE :: fac(:)
   REAL(KIND=DP) :: qq
   INTEGER :: ig

! variables for computing P|x> = |op_x>
   COMPLEX(KIND=DP), ALLOCATABLE :: op_x(:)
   COMPLEX(KIND=DP), ALLOCATABLE :: ovp_x(:)
   COMPLEX(KIND=DP), ALLOCATABLE :: fcw_state(:,:) ! Only needed for call, not used
   REAL(KIND=DP), ALLOCATABLE :: hdiag(:) ! Only needed for call, not used (pmat_type=2)
   REAL(KIND=DP), ALLOCATABLE :: fcw_mat(:,:) ! Only needed for call, not used
   REAL(KIND=DP), ALLOCATABLE :: v_states(:,:) ! The valence states
   INTEGER :: fcw_number ! Only needed for call, not used

   INTEGER :: loop, lp_v

   ALLOCATE(fac(1:npw))
   ALLOCATE(op_x(1:npw))
   ALLOCATE(ovp_x(1:npw))

   !v_c
   IF (l_truncated_coulomb) THEN
     DO ig=1,npw
          qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
        IF (qq > 1.d-8) THEN
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-DCOS(DSQRT(qq)*truncation_radius*tpiba))
        ELSE
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
        END IF
     END DO
     fac(:)=fac(:)/omega
   ELSE
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   END IF   
   ! Allocation and assignment of "fake conduction" states
   ! NOTE: not used during call since pmat_type=0
   fcw_number=1
   ALLOCATE(fcw_state(npw,fcw_number))
   ALLOCATE(fcw_mat(fcw_number,fcw_number))
   ALLOCATE(hdiag(1))

   ! Have to reobtain the valence_states to avoid having them in call
   ALLOCATE(v_states(dfftp%nnr,num_nbndv(1)))
   CALL evc_to_real(num_nbndv(1), v_states) ! gives v_states in real space

   ! P operator onto state
   CALL o_1psi_gamma(num_nbndv(1), v_states, x, op_x, .FALSE., hdiag, &
                     0, fcw_number, fcw_state, fcw_mat, pmat_ethr) 

   ! -1*(v_c) operator onto state
   ovp_x(1:npw) = op_x(1:npw) * fac(1:npw)
   ovp_x(1:npw) = -4.d0 * ovp_x(1:npw)

   ! calcualte (1-vP)|x> = |x> + |ovp_x>
   DO loop = 1,npw
      u(loop,1) = x(loop) + ovp_x(loop)
   END DO

   DEALLOCATE(v_states)
   DEALLOCATE(hdiag)
   DEALLOCATE(fcw_state)
   DEALLOCATE(fcw_mat)
   DEALLOCATE(fac)
   DEALLOCATE(op_x)
   DEALLOCATE(ovp_x)

   RETURN

END SUBROUTINE operator_1_vp 

