!
! Copyright (C) 2001-2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dy_base &
     ( npw, npwx, igk, xk, nat, tau, ityp, ntyp, tpiba, &
       omega, nr1, nr2, nr3, eigts1, eigts2, eigts3, mill, g, u, dvkb )
  !----------------------------------------------------------------------
  !! Calculates the beta functions of the pseudopotential with the
  !! derivative of the spherical harmonics projected on vector u
  !
  USE upf_kinds,   ONLY: dp
  USE upf_const,   ONLY: tpi
  USE uspp,        ONLY: nkb, indv, nhtol, nhtolm
  USE uspp_data,   ONLY: nqx, tab, dq
  USE uspp_param,  ONLY: upf, lmaxkb, nbetam, nh
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw
  !! number ok plane waves 
  INTEGER, INTENT(IN) :: npwx
  !! max number ok plane waves across k-points
  INTEGER, INTENT(IN) :: igk(npw)
  !! indices of plane waves k+G
  REAL(dp), INTENT(IN) :: xk(3)
  !! k-point
  INTEGER, INTENT(IN) :: nat
  !! number of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! index of type per atom
  INTEGER, INTENT(IN) :: ntyp
  !! number of atomic types
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions (cc alat units)
  REAL(DP), INTENT(IN) :: tpiba
  !! rec.lattice units 2pi/a
  REAL(DP), INTENT(IN) :: omega
  !! cell volume
  INTEGER, INTENT(IN) :: nr1,nr2,nr3
  !! fft dims (dense grid)
  COMPLEX(DP), INTENT(IN) :: eigts1(-nr1:nr1,nat)
  !! structure factor 1
  COMPLEX(DP), INTENT(IN) :: eigts2(-nr2:nr2,nat)
  !! structure factor 2
  COMPLEX(DP), INTENT(IN) :: eigts3(-nr3:nr3,nat)
  !! structure factor 3
  INTEGER, INTENT(IN) :: mill(3,*)
  !! miller index map
  REAL(DP), INTENT(IN) :: g(3,*)
  !! g vectors (2pi/a units)
  REAL(DP), INTENT(IN) :: u(3)
  !! input: projection vector
  COMPLEX(DP), INTENT(OUT) :: dvkb(npwx,nkb)
  !! output: kleinman-bylander pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: na, nt, nb, ih, l, lm, ikb, iig, ipol, i0, i1, i2, &
             i3, ig
  REAL(DP), ALLOCATABLE :: gk(:,:), q(:)
  REAL(DP) :: px, ux, vx, wx, arg
  !
  REAL(DP), ALLOCATABLE :: vkb0(:,:,:), dylm(:,:), dylm_u(:,:)
  ! dylm = d Y_lm/dr_i in cartesian axes
  ! dylm_u as above projected on u
  !
  COMPLEX(DP), ALLOCATABLE :: sk(:)
  COMPLEX(DP) :: phase, pref
  !
  INTEGER :: iq
  !
  dvkb(:,:) = (0.d0, 0.d0)
  IF (lmaxkb <= 0) RETURN
  !
  ALLOCATE( vkb0(npw,nbetam,ntyp), dylm_u(npw,(lmaxkb+1)**2), gk(3,npw) )
  ALLOCATE( q(npw) )
  !
  DO ig = 1, npw
     iig = igk(ig)
     gk(1, ig) = xk(1) + g(1, iig)
     gk(2, ig) = xk(2) + g(2, iig)
     gk(3, ig) = xk(3) + g(3, iig)
     q(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  ENDDO
  !
  ALLOCATE( dylm(npw,(lmaxkb+1)**2) )
  dylm_u(:,:) = 0.d0
  DO ipol = 1, 3
     CALL dylmr2( (lmaxkb+1)**2, npw, gk, q, dylm, ipol )
     CALL daxpy( npw * (lmaxkb + 1) **2, u(ipol), dylm, 1, dylm_u, 1 )
  ENDDO
  DEALLOCATE( dylm )
  !
  DO ig = 1, npw
     q(ig) = SQRT(q(ig)) * tpiba
  ENDDO
  !
  DO nt = 1, ntyp
     ! calculate beta in G-space using an interpolation table
     DO nb = 1, upf(nt)%nbeta
        DO ig = 1, npw
           px = q(ig)/dq - INT(q(ig)/dq)
           ux = 1.d0 - px
           vx = 2.d0 - px
           wx = 3.d0 - px
           i0 = q(ig)/dq + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           vkb0(ig, nb, nt) = tab(i0, nb, nt) * ux * vx * wx / 6.d0 + &
                              tab(i1, nb, nt) * px * vx * wx / 2.d0 - &
                              tab(i2, nb, nt) * px * ux * wx / 2.d0 + &
                              tab(i3, nb, nt) * px * ux * vx / 6.d0
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE( q )
  !
  ALLOCATE( sk(npw) )
  !
  ikb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           arg = (xk(1) * tau(1, na) + xk(2) * tau(2, na) &
                + xk(3) * tau(3, na) ) * tpi
           phase = CMPLX( COS(arg), -SIN(arg), KIND=DP )
           DO ig = 1, npw
              iig = igk(ig)
              sk(ig) = eigts1(mill (1,iig), na) * &
                       eigts2(mill (2,iig), na) * &
                       eigts3(mill (3,iig), na) * phase
           ENDDO
           !
           DO ih = 1, nh(nt)
              nb = indv(ih, nt)
              l = nhtol(ih, nt)
              lm = nhtolm(ih, nt)
              ikb = ikb + 1
              pref = (0.d0, -1.d0)**l
              !
              DO ig = 1, npw
                 dvkb(ig, ikb) = vkb0(ig, nb, nt) * sk(ig) * dylm_u(ig, lm) &
                      * pref / tpiba
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  IF (ikb /= nkb) CALL upf_error( 'gen_us_dy', 'unexpected error', 1 )
  !
  DEALLOCATE( sk )
  DEALLOCATE( vkb0, dylm_u, gk )
  !
  RETURN
  !
END SUBROUTINE gen_us_dy_base
!
