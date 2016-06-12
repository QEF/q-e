!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
SUBROUTINE lr_addus_dvpsi ( ik, lda, n, m, psi, dpsi )
  !---------------------------------------------------------------------------
  !
  ! ... Calculate the ultrasoft term of the perturbation exp(iq*r),
  ! ... and then sum up the input wavefunction and the ultrasoft term.
  !
  ! ... input:
  !
  ! ...    ik    given k point
  ! ...    lda   leading dimension of the array psi
  ! ...    n     true dimension of psi
  ! ...    m     number of bands of psi
  !
  ! Written by Iurii Timrov (2015)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE paw_variables,        ONLY : okpaw
  USE noncollin_module,     ONLY : npol, noncolin
  USE uspp,                 ONLY : vkb, nkb, indv_ijkb0
  USE cell_base,            ONLY : omega
  USE lrus,                 ONLY : becp1
  USE qpoint,               ONLY : xq, eigqts
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ik, lda, n, m
  COMPLEX(DP), INTENT(in) :: psi(lda*npol,m)
  ! input: wavefunction u_n,k
  COMPLEX(DP), INTENT(out) :: dpsi(lda*npol,m)
  ! output: sum of the input wavefunction and the USPP term
  !
  ! the local variables
  !
  INTEGER :: na, nt, ih, jh
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  !
  REAL(DP) :: qmod                              ! the modulus of q
  REAL(DP), ALLOCATABLE :: ylmk0(:)             ! the spherical harmonics
  COMPLEX(DP) :: qgm(1)                         ! FFT of Q(r)
  COMPLEX(DP), ALLOCATABLE :: ps(:,:), qqc(:,:) ! auxiliary arrays
  !
  dpsi = psi
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('lr_addus_dvpsi')
  !
  IF (noncolin) CALL errore( 'lr_addus_dvpsi', 'Noncollinear case is not supported', 1 )
  IF (okpaw)    CALL errore( 'lr_addus_dvpsi', 'PAW is not supported', 1 )
  !
  ALLOCATE (ylmk0(lmaxq*lmaxq))
  ALLOCATE (ps(nkb,m))
  ps(:,:) = (0.d0, 0.d0)
  !
  qmod = xq(1)**2 + xq(2)**2 + xq(3)**2 
  !
  ! Calculate sphrecial harmonics
  !
  CALL ylmr2 (lmaxq*lmaxq, 1, xq, qmod, ylmk0)
  !
  qmod = sqrt(qmod)
  !
  DO nt = 1, ntyp
    IF (upf(nt)%tvanp) THEN
      !
      ! Calculate the Fourier transform of the Q functions.
      !
      ALLOCATE (qqc(nh(nt),nh(nt)))
      qqc(:,:) = (0.d0, 0.d0)
      !
      DO ih = 1, nh(nt)
        DO jh = ih, nh(nt) 
           !
           CALL qvan2 (1, ih, jh, nt, qmod, qgm(1), ylmk0)
           qqc(ih,jh) = omega * qgm(1)
           qqc(jh,ih) = qqc(ih,jh)
           !
        ENDDO
      ENDDO
      !
      ! Calculate a product of Q^*(q) with <beta|evc>
      !
      DO na = 1, nat
         IF (ityp (na).eq.nt) THEN
           !
           qqc = qqc * eigqts(na)
           !
           CALL ZGEMM('C','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                    & qqc, nh(nt), becp1(ik)%k(indv_ijkb0(na)+1,1), nkb, &
                    & (0.0_dp,0.0_dp), ps(indv_ijkb0(na)+1,1), nkb )
           !
         ENDIF
      ENDDO
      !
      DEALLOCATE (qqc)
      !
    ENDIF
  ENDDO
  !
  ! Sum up the normal and ultrasoft terms.
  !
  IF ( m == 1 ) THEN
     !
     CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb, lda, &
               & ps, 1, ( 1.D0, 0.D0 ), dpsi, 1 )
     !
  ELSE
     !
     CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb, lda, &
              & ps, nkb, ( 1.D0, 0.D0 ), dpsi, lda )
     !
  ENDIF
  !
  DEALLOCATE (ylmk0)
  DEALLOCATE (ps)
  !
  CALL stop_clock ('lr_addus_dvpsi')
  !
  RETURN
  !
END SUBROUTINE lr_addus_dvpsi
