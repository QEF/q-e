!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_tab_qrad (omega, intra_bgrp_comm)
  !----------------------------------------------------------------------
  !
  ! Compute interpolation table qrad(i,nm,l+1,nt) = Q^{(L)}_{nm,nt}(q_i)
  ! of angular momentum L, for atom of type nt, on grid q_i, where
  ! nm = combined (n,m) index; n,m = 1,...,nbeta (number of beta functions)
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi
  USE atom,         ONLY : rgrid
  USE uspp_param,   ONLY : upf, lmaxq, nbetam, nsp
  USE uspp_data,    ONLY : nqxq, dq, qrad, qrad_d
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  real(DP), intent(in) :: omega
  integer,  intent(in) :: intra_bgrp_comm
  !
  INTEGER :: ndm, startq, lastq, nt, l, nb, mb, ijv, iq, ir
  ! various indices
  REAL(dp) :: prefr
  ! the prefactor of the Q functions
  REAL(dp) :: q
  REAL(dp), ALLOCATABLE :: aux (:), besr (:)
  ! various work space
  !
  prefr = fpi / omega
  ndm = MAXVAL ( upf(:)%kkbeta )
  ALLOCATE (aux ( ndm))
  ALLOCATE (besr( ndm))
  !
  CALL divide (intra_bgrp_comm, nqxq, startq, lastq)
  !
  qrad(:,:,:,:)= 0.d0
  DO nt = 1, nsp
     if ( upf(nt)%tvanp ) then
        DO l = 0, upf(nt)%nqlc -1
           !
           !     l is the true (combined) angular momentum
           !     Note that the index of array qfuncl runs from 0 to l,
           !     while the same index for qrad runs from 1 to l+1
           !     FIXME: qrad has "holes" if USPP/PAW do not precede NCPP
           !
           DO iq = startq, lastq
              !
              q = (iq - 1) * dq
              !
              !     here we compute the spherical bessel function for each q_i
              !
              CALL sph_bes ( upf(nt)%kkbeta, rgrid(nt)%r, q, l, besr)
              !
              DO nb = 1, upf(nt)%nbeta
                 !
                 !    the Q are symmetric with respect to nb,nm indices
                 !
                 DO mb = nb, upf(nt)%nbeta
                    ijv = mb * (mb - 1) / 2 + nb
                    IF ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .AND. &
                         ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .AND. &
                         (mod(l+upf(nt)%lll(nb)+upf(nt)%lll(mb),2)==0) ) THEN
                       DO ir = 1, upf(nt)%kkbeta
                          aux  (ir) = besr (ir) * upf(nt)%qfuncl(ir,ijv,l)
                       ENDDO
                       !
                       !   and then we integrate with all the Q functions
                       !
                       CALL simpson ( upf(nt)%kkbeta, aux, rgrid(nt)%rab, &
                                     qrad(iq,ijv,l+1, nt) )
                    ENDIF
                 ENDDO
              ENDDO
              ! igl
           ENDDO
           ! l
        ENDDO
        qrad (:, :, :, nt) = qrad (:, :, :, nt)*prefr
        CALL mp_sum ( qrad (:, :, :, nt), intra_bgrp_comm )
     ENDIF
     ! nsp
  ENDDO
  !
  DEALLOCATE (besr)
  DEALLOCATE (aux)
  !
  ! update GPU memory (taking care of zero-dim allocations)
  !
#if defined __CUDA
  if ( nbetam > 0 .and. lmaxq > 0 ) qrad_d=qrad
#endif
  !
END SUBROUTINE init_tab_qrad
