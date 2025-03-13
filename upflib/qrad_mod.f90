!
! Copyright (C) 2021-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qrad_mod
  !
  !! Variables and routines for augmentation charges in numerical form
  !! Contains generation of interpolation tables in reciprocal space,
  !! interpolation routines and other utility routines
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: dq, tab_qrad
  PUBLIC :: init_tab_qrad
  PUBLIC :: scale_tab_qrad
  PUBLIC :: deallocate_tab_qrad
  !
  SAVE
  !
  INTEGER :: nqx = 0
  !! size of interpolation table
  REAL(DP), PARAMETER:: dq = 0.01_dp
  !! grid step for interpolation table
  REAL(DP) :: qmax = 0.0_dp 
  !! max q covered by the interpolation table
  REAL(DP), ALLOCATABLE :: tab_qrad(:,:,:,:)
  !! interpolation table for numerical pseudopotentials
  !
CONTAINS
!----------------------------------------------------------------------
  SUBROUTINE init_tab_qrad (qmax_, omega, comm, ierr)
  !----------------------------------------------------------------------
  !
  !! Allocate and fill interpolation table tab_qrad:
  !!    tab_qrad(i,nm,l+1,nt) = Q^{(L)}_{nm,nt}(q_i)
  !! for angular momentum L, for atom of type nt, on grid q_i, where
  !! nm = combined (n,m) index; n,m = 1,...,nbeta (number of beta functions)
  !
  USE atom,         ONLY : rgrid
  USE uspp_param,   ONLY : upf, lmaxq, nbetam, nsp
  USE mp,           ONLY : mp_sum
  !
  INTEGER, INTENT(IN)  :: comm
  !! MPI communicator, to split the workload
  INTEGER, INTENT(OUT) :: ierr
  !! error code: ierr = 0 if interpolation table (IT) was allocated
  !!             ierr =-1 if IT had insufficient dimension and was re-allocated
  !!             ierr =-2 if IT was already present and nothing is done
  !!             ierr =-3 if IT not needed and nothing is done
  REAL(dp), INTENT(IN) :: omega
  !! Unit-cell volume
  REAL(dp), INTENT(IN) :: qmax_
  !! Interpolate q up to qmax_ (sqrt(Ry), q^2 is an energy)
  !
  INTEGER :: ndm, startq, lastq, nt, l, nb, mb, ijv, iq, ir
  ! various indices
  REAL(dp) :: q
  REAL(dp), ALLOCATABLE :: aux (:), besr (:)
  ! various work space
  !
  ierr = -3
  IF ( lmaxq <= 0 .OR. ALL ( .NOT. upf(1:nsp)%tvanp ) ) RETURN
  IF ( .NOT. ALLOCATED(tab_qrad) ) THEN
     !! table not yet allocated
     qmax = qmax_
     ierr = 0
  ELSE IF ( qmax_ > qmax ) THEN
     !! table ìs allocated but dimension insufficient: re-allocate
     !! (with some margin so that this does not happen too often)
     !$acc exit data delete(tab_qrad)
     DEALLOCATE ( tab_qrad )
     qmax = qmax_ + MAX(dq*100,qmax_-qmax)
     ierr =-1
  ELSE
     !! table already computed: exit
     ierr =-2
     RETURN
  END IF
  nqx = INT( qmax/dq + 4)
  ALLOCATE (tab_qrad(nqx,nbetam*(nbetam+1)/2, lmaxq, nsp))
  !$acc enter data create(tab_qrad)
  !
  ndm = MAXVAL ( upf(:)%kkbeta )
  ALLOCATE (aux ( ndm))
  ALLOCATE (besr( ndm))
  !
  CALL divide (comm, nqx, startq, lastq)
  !
  tab_qrad(:,:,:,:)= 0.0_dp
  DO nt = 1, nsp
     if ( upf(nt)%tvanp ) then
        DO l = 0, upf(nt)%nqlc -1
           !
           !     l is the true (combined) angular momentum
           !     Note that the index of array qfuncl runs from 0 to l,
           !     while the same index for tab_qrad runs from 1 to l+1
           !     FIXME: tab_qrad has "holes" if USPP/PAW do not precede NCPP
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
                                     tab_qrad(iq,ijv,l+1, nt) )
                    ENDIF
                 ENDDO
              ENDDO
              ! igl
           ENDDO
           ! l
        ENDDO
        tab_qrad (:, :, :, nt) = tab_qrad (:, :, :, nt) * fpi / omega

        CALL mp_sum ( tab_qrad (:, :, :, nt), comm )
     ENDIF
     ! nsp
  ENDDO
  !
  DEALLOCATE (besr)
  DEALLOCATE (aux)
  !
  !$acc update device(tab_qrad)
  !
END SUBROUTINE init_tab_qrad
!
SUBROUTINE scale_tab_qrad ( vol_ratio_m1 )
  !
  REAL(DP), INTENT(in) :: vol_ratio_m1
  !! vol_ratio_m1 = omega_old / omega
  !
  IF ( ALLOCATED ( tab_qrad ) ) THEN
     tab_qrad(:,:,:,:) = tab_qrad(:,:,:,:) * vol_ratio_m1
     !$acc update device (tab_qrad)
  END IF
  !
END SUBROUTINE scale_tab_qrad
!
SUBROUTINE deallocate_tab_qrad ( )
  !
  IF ( ALLOCATED ( tab_qrad ) ) THEN
     !$acc exit data delete(tab_qrad)
     DEALLOCATE (tab_qrad)
  END IF
  qmax = 0.0_dp
  !
END SUBROUTINE deallocate_tab_qrad
!
END MODULE qrad_mod
