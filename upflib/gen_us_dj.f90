!
! Copyright (C) 2001-2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dj_base &
     ( npw, npwx, igk, xk, nat, tau, ityp, ntyp, tpiba, &
       omega, nr1, nr2, nr3, eigts1, eigts2, eigts3, mill, g, dvkb )
  !----------------------------------------------------------------------
  !! Calculates the beta function pseudopotentials with
  !! the derivative of the Bessel functions.
  !
  USE upf_kinds,  ONLY: dp
  USE upf_const,  ONLY: tpi
  USE uspp,       ONLY: nkb, indv, nhtol, nhtolm
  USE uspp_data,  ONLY: nqx, tab, dq
  USE uspp_param, ONLY: upf, lmaxkb, nbetam, nh
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
  COMPLEX(DP), INTENT(OUT) :: dvkb(npwx, nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: ikb, nb, ih, ig, i0, i1, i2, i3, nt
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! counter on G vectors
  ! index of the first nonzero point in the r
  ! counter on atomic type
  !
  REAL(DP) :: arg, px, ux, vx, wx
  ! argument of the atomic phase factor
  !
  COMPLEX(DP) :: phase, pref
  ! atomic phase factor
  ! prefactor
  !
  INTEGER :: na, l, iig, lm, iq
  REAL(DP), ALLOCATABLE :: djl(:,:,:), ylm(:,:), q(:), gk(:,:)
  REAL(DP) :: qt
  COMPLEX(DP), ALLOCATABLE :: sk(:)
  !
  IF (nkb == 0) RETURN
  !
  CALL start_clock( 'stres_us31' )
  !
  ALLOCATE( djl(npw,nbetam,ntyp)   )    
  ALLOCATE( ylm(npw,(lmaxkb+1)**2) )    
  ALLOCATE( gk(3,npw) )    
  ALLOCATE( q(npw)    )    
  !
  DO ig = 1, npw
     iig = igk(ig)
     gk(1, ig) = xk(1) + g(1, iig)
     gk(2, ig) = xk(2) + g(2, iig)
     gk(3, ig) = xk(3) + g(3, iig)
     q(ig) = gk(1, ig)**2 + gk(2, ig)**2 + gk(3, ig)**2
  ENDDO
  !
  CALL stop_clock( 'stres_us31' )
  CALL start_clock( 'stres_us32' )
  CALL ylmr2( (lmaxkb+1)**2, npw, gk, q, ylm )
  CALL stop_clock( 'stres_us32' )
  CALL start_clock( 'stres_us33' )
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nbeta
        !
        DO ig = 1, npw
           qt = SQRT(q (ig)) * tpiba
           px = qt / dq - INT(qt/dq)
           ux = 1.d0 - px
           vx = 2.d0 - px
           wx = 3.d0 - px
           i0 = qt / dq + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           djl(ig,nb,nt) = ( tab(i0, nb, nt) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                             tab(i1, nb, nt) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                             tab(i2, nb, nt) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                             tab(i3, nb, nt) * (+ux*vx-px*vx-px*ux)/6.d0 )/dq
        ENDDO
        !
     ENDDO
  ENDDO
  !
  CALL stop_clock( 'stres_us33' )
  CALL start_clock( 'stres_us34' )
  !
  DEALLOCATE ( q  )
  DEALLOCATE ( gk )
  !
  ALLOCATE ( sk(npw) )
  ikb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        !
        IF (ityp(na) == nt) THEN
           arg = ( xk(1) * tau(1,na) + &
                   xk(2) * tau(2,na) + &
                   xk(3) * tau(3,na) ) * tpi
           phase = CMPLX( COS(arg), -SIN(arg) ,KIND=DP )
           DO ig = 1, npw
              iig = igk(ig)
              sk (ig) = eigts1(mill (1,iig), na) * &
                        eigts2(mill (2,iig), na) * &
                        eigts3(mill (3,iig), na) * phase
           ENDDO
           DO ih = 1, nh(nt)
              nb = indv(ih, nt)
              l = nhtol(ih, nt)
              lm= nhtolm(ih, nt)
              ikb = ikb + 1
              pref = (0.d0, -1.d0) **l
              !
              DO ig = 1, npw
                 dvkb(ig, ikb) = djl(ig, nb, nt) * sk(ig) * ylm(ig, lm) &
                      * pref
              ENDDO
           ENDDO
        ENDIF
        !
     ENDDO
  ENDDO
  !
  CALL stop_clock('stres_us34')
  !
  IF (ikb /= nkb) CALL upf_error('gen_us_dj', 'unexpected error', 1)
  DEALLOCATE( sk  )
  DEALLOCATE( ylm )
  DEALLOCATE( djl )
  !
  RETURN
  !
END SUBROUTINE gen_us_dj_base

