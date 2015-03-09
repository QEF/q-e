!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE addusdens(rho)
  !----------------------------------------------------------------------
  !
  USE realus,               ONLY : addusdens_r
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  REAL(kind=dp), intent(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  IF ( tqr ) THEN
     CALL addusdens_r(rho,.true.)
  ELSE
#if defined(__CUDA) && !defined(__DISABLE_CUDA_ADDUSDENS)
     CALL addusdens_g_gpu(rho)
#else
     CALL addusdens_g(rho)
#endif
  END IF
  !
  RETURN
  !
END SUBROUTINE addusdens
!
!----------------------------------------------------------------------
subroutine addusdens_g(rho)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, nl, nlm, gg, g, &
                                   eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : becsum, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  implicit none
  !
  REAL(kind=dp), intent(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ijh, is, nab, nij
  ! counters

  real(DP), allocatable :: tbecsum(:,:)
  real(DP), allocatable :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics

  complex(DP), allocatable :: skk(:), aux2(:,:)
  complex(DP), allocatable ::  aux (:,:), qgm(:,:)
  ! work space for rho(G,nspin)
  ! Fourier transform of q

  if (.not.okvan) return

  call start_clock ('addusdens')

  allocate (aux ( ngm, nspin_mag))    
  allocate (qmod( ngm))    
  allocate (ylmk0( ngm, lmaxq * lmaxq))    
  ALLOCATE ( skk(ngm), aux2(ngm,nspin_mag) )

  aux (:,:) = (0.d0, 0.d0)
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
  !
  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        allocate (qgm(ngm,nij), tbecsum(nij,nspin_mag) )
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              ijh = ijh + 1
              call qvan2 (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0)
           end do
        end do
        !
        do na = 1, nat
           IF ( ityp(na) == nt ) THEN
              !
              tbecsum(:,:) = becsum(1:nij,na,1:nspin_mag)
              !
!$omp parallel do default(shared) private(ig)
              do ig = 1, ngm
                 skk(ig) = eigts1 (mill (1,ig), na) * &
                           eigts2 (mill (2,ig), na) * &
                           eigts3 (mill (3,ig), na)
              end do
!$omp end parallel do
              CALL dgemm( 'N', 'N', 2*ngm, nspin_mag, nij, 1.0_dp, qgm, 2*ngm,&
                          tbecsum, nij, 0.0_dp, aux2, 2*ngm )
              do is = 1, nspin_mag
!$omp parallel do default(shared) private(ig)
                 do ig = 1, ngm
                    aux(ig,is)=aux(ig,is) + aux2(ig,is)*skk(ig)
                 enddo
!$omp end parallel do
              enddo
           endif
        enddo
        deallocate (tbecsum, qgm)
     endif
  enddo
  !
  deallocate (aux2, skk)
  deallocate (ylmk0)
  deallocate (qmod)
  !
  !     convert aux to real space and add to the charge density
  !
#ifdef DEBUG_ADDUSDENS
  call start_clock ('addus:fft')
#endif
  do is = 1, nspin_mag
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     if (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is))
     CALL invfft ('Dense', psic, dfftp)
     rho(:, is) = rho(:, is) +  DBLE (psic (:) )
  enddo
#ifdef DEBUG_ADDUSDENS
  call stop_clock ('addus:fft')
#endif
  deallocate (aux)

  call stop_clock ('addusdens')
  return
end subroutine addusdens_g

