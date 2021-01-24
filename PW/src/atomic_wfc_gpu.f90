!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_gpu( ik, wfcatom_d )
  !-----------------------------------------------------------------------
  !! This routine computes the superposition of atomic wavefunctions
  !! for k-point "ik" - output in "wfcatom".
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi, fpi, pi
  USE cell_base,        ONLY : omega, tpiba
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,            ONLY : natomwfc
  !USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,            ONLY : xk, ngk, igk_k_d !, igk_k
  USE wvfct,            ONLY : npwx
  USE us,               ONLY : tab_at, dq
  USE uspp_param,       ONLY : upf
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE spin_orb,         ONLY : lspinorb, rot_ylm, fcoef, lmaxx, domag, &
                               starting_spin_angle
  USE mp_bands,         ONLY : inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  USE gvect_gpum,       ONLY : mill_d, eigts1_d, eigts2_d, eigts3_d, g_d
  USE us_gpum,          ONLY : using_tab_at, using_tab_at_d, tab_at_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom_d(npwx,npol,natomwfc)
  !! Superposition of atomic wavefunctions
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm, npw
  COMPLEX(DP) :: kphase, lphase
  REAL(DP)    :: arg, px, ux, vx, wx
  INTEGER     :: ig_start, ig_end, mil1, mil2, mil3
  !
  REAL(DP) :: xk1, xk2, xk3, qgr
  REAL(DP), ALLOCATABLE :: chiq_d(:,:,:)
  REAL(DP), ALLOCATABLE :: ylm_d(:,:), gk_d(:,:), qg_d(:)
  COMPLEX(DP), ALLOCATABLE :: sk_d(:), aux_d(:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfcatom_d, ylm_d, gk_d, qg_d, sk_d, chiq_d, aux_d
#endif  
  !
  CALL start_clock( 'atomic_wfc' )

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, ntyp
     lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  END DO
  !
  nwfcm = MAXVAL( upf(1:ntyp)%nwfc )
  npw = ngk(ik)
  !
  ALLOCATE( ylm_d(npw,(lmax_wfc+1)**2), gk_d(3,npw), qg_d(npw) )
  ALLOCATE( chiq_d(npw,nwfcm,ntyp), sk_d(npw) )
  !
  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, npw
     iig = igk_k_d(ig,ik)
     gk_d(1,ig) = xk1 + g_d(1,iig)
     gk_d(2,ig) = xk2 + g_d(2,iig)
     gk_d(3,ig) = xk3 + g_d(3,iig)
     qg_d(ig) = gk_d(1,ig)**2 +  gk_d(2,ig)**2 + gk_d(3,ig)**2
  END DO
  !
  !  ylm = spherical harmonics
  !
  CALL ylmr2_gpu( (lmax_wfc+1)**2, npw, gk_d, qg_d, ylm_d )
  !
  ! from now to the end of the routine the ig loops are distributed across bgrp
  CALL divide( inter_bgrp_comm, npw, ig_start, ig_end )
  !
  ! set now q=|k+G| in atomic units
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ig = ig_start, ig_end
     qg_d(ig) = SQRT( qg_d(ig) )*tpiba
  END DO
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  CALL using_tab_at(0)
  CALL using_tab_at_d(0)
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           !
           !$cuf kernel do (1) <<<*,*>>>
           DO ig = ig_start, ig_end
              qgr = qg_d(ig)
              px = qgr / dq - DBLE(INT(qgr/dq))
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT(qgr/dq) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq_d(ig,nb,nt) = &
                     tab_at_d(i0,nb,nt) * ux * vx * wx / 6.d0 + &
                     tab_at_d(i1,nb,nt) * px * vx * wx / 2.d0 - &
                     tab_at_d(i2,nb,nt) * px * ux * wx / 2.d0 + &
                     tab_at_d(i3,nb,nt) * px * ux * vx / 6.d0
           END DO
           !
        END IF
     END DO
  END DO

  DEALLOCATE( qg_d, gk_d )
  !
  IF (noncolin) ALLOCATE( aux_d(npw) )
  !
  wfcatom_d(:,:,:) = (0.0_dp, 0.0_dp)
  !
  DO na = 1, nat
     arg = (xk1*tau(1,na) + xk2*tau(2,na) + xk3*tau(3,na)) * tpi
     kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
     !
     !     sk is the structure factor
     !
     !$cuf kernel do (1) <<<*,*>>>
     DO ig = ig_start, ig_end
        iig = igk_k_d(ig,ik)
        mil1 = mill_d(1,iig)
        mil2 = mill_d(2,iig)
        mil3 = mill_d(3,iig)
        sk_d(ig) = kphase * eigts1_d(mil1,na) * &
                            eigts2_d(mil2,na) * &
                            eigts3_d(mil3,na)
     END DO
     !
     nt = ityp(na)
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           l = upf(nt)%lchi(nb)
           lphase = (0.d0,1.d0)**l
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 IF (starting_spin_angle.OR..NOT.domag) THEN
                    CALL atomic_wfc_so_gpu( )
                 ELSE
                    CALL atomic_wfc_so_mag_gpu( )
                 END IF
                 !
              ELSE
                 !
                 CALL atomic_wfc_nc_gpu( )
                 !
              END IF
              !
           ELSE
              !
              CALL atomic_wfc___gpu( )
              !
           END IF
           !
        END IF
        !
     END DO
     !
  END DO

  IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
       'internal error: some wfcs were lost ', 1 )

  IF (noncolin) DEALLOCATE( aux_d )
  DEALLOCATE( sk_d, chiq_d, ylm_d )

  ! collect results across bgrp
  CALL mp_sum( wfcatom_d, inter_bgrp_comm )

  CALL stop_clock( 'atomic_wfc' )
  RETURN

CONTAINS
!----------------------------------------------------------------
  SUBROUTINE atomic_wfc_so_gpu( )
   !------------------------------------------------------------
   !! Spin-orbit case.
   !
   REAL(DP) :: fact(2), fact_is, j
   COMPLEX(DP) :: rot_ylm_in1
   REAL(DP), EXTERNAL :: spinor
   INTEGER :: ind, ind1, n1, is, sph_ind
   !
   j = upf(nt)%jchi(nb)
   DO m = -l-1, l
      fact(1) = spinor(l,j,m,1)
      fact(2) = spinor(l,j,m,2)
      IF ( ABS(fact(1)) > 1.d-8 .OR. ABS(fact(2)) > 1.d-8 ) THEN
         n_starting_wfc = n_starting_wfc + 1
         IF (n_starting_wfc > natomwfc) CALL errore &
              ('atomic_wfc_so', 'internal error: too many wfcs', 1)
         !     
         DO is = 1, 2
            fact_is = fact(is)
            IF (ABS(fact(is)) > 1.d-8) THEN
               ind = lmaxx + 1 + sph_ind(l,j,m,is)
               aux_d = (0.d0,0.d0)
               DO n1 = 1, 2*l+1
                  ind1 = l**2+n1
                  rot_ylm_in1 = rot_ylm(ind,n1)
                  IF (ABS(rot_ylm_in1) > 1.d-8) THEN
                    !$cuf kernel do (1) <<<*,*>>>
                    DO ig = ig_start, ig_end
                      aux_d(ig) = aux_d(ig) + rot_ylm_in1 * &
                                  CMPLX(ylm_d(ig,ind1))
                    ENDDO
                  ENDIF
               ENDDO
               !$cuf kernel do (1) <<<*,*>>>
               DO ig = ig_start, ig_end
                  wfcatom_d(ig,is,n_starting_wfc) = lphase * &
                                sk_d(ig)*aux_d(ig)*CMPLX(fact_is* &
                                chiq_d(ig,nb,nt))
               END DO
            ELSE
                wfcatom_d(:,is,n_starting_wfc) = (0.d0,0.d0)
            END IF
         END DO
         !
      END IF
   END DO
   !
  END SUBROUTINE atomic_wfc_so_gpu
  ! 
  SUBROUTINE atomic_wfc_so_mag_gpu( )
   !
   !! Spin-orbit case, magnetization along "angle1" and "angle2"
   !! In the magnetic case we always assume that magnetism is much larger
   !! than spin-orbit and average the wavefunctions at l+1/2 and l-1/2
   !! filling then the up and down spinors with the average wavefunctions,
   !! according to the direction of the magnetization, following what is
   !! done in the noncollinear case.
   !
   REAL(DP) :: alpha, gamman, j
   COMPLEX(DP) :: fup, fdown  
   REAL(DP), ALLOCATABLE :: chiaux_d(:)
   INTEGER :: nc, ib, ig
   !
#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: chiaux_d
#endif
   !
   j = upf(nt)%jchi(nb)
   !
   !  This routine creates two functions only in the case j=l+1/2 or exit in the
   !  other case 
   !    
   IF (ABS(j-l+0.5_DP)<1.d-4) RETURN

   ALLOCATE( chiaux_d(npw) )
   !
   !  Find the functions j=l-1/2
   !
   IF (l == 0)  THEN
      !$cuf kernel do
      DO ig = 1, npw
         chiaux_d(ig) = chiq_d(ig,nb,nt)
      END DO
   ELSE
      DO ib = 1, upf(nt)%nwfc
         IF ((upf(nt)%lchi(ib) == l).AND. &
                      (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
            nc = ib
            EXIT
         ENDIF
      ENDDO
      !
      !  Average the two functions
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = 1, npw
        chiaux_d(ig) = (chiq_d(ig,nb,nt)*DBLE(l+1)+chiq_d(ig,nc,nt)*l)/ &
                       DBLE(2*l+1)
      ENDDO
      !
   ENDIF
   !
   !  and construct the starting wavefunctions as in the noncollinear case.
   !
   alpha = angle1(nt)
   gamman = - angle2(nt) + 0.5d0*pi
   !
   DO m = 1, 2*l+1
      lm = l**2+m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc + 2*l+1 > natomwfc ) CALL errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = ig_start, ig_end
        aux_d(ig) = sk_d(ig)* CMPLX(ylm_d(ig,lm)*chiaux_d(ig))
      END DO
      !
      ! now, rotate wfc as needed
      ! first : rotation with angle alpha around (OX)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = ig_start, ig_end
         fup = CMPLX(COS(0.5d0*alpha))*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*alpha))*aux_d(ig)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom_d(ig,1,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman)) &
                        +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fup
         wfcatom_d(ig,2,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman)) &
                        -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = CMPLX(COS(0.5d0*(alpha+pi)))*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*(alpha+pi)))*aux_d(ig)
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom_d(ig,1,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman)) &
                  +(0.d0,1.d0)*CMPLX(SIN(0.5d0 *gamman)))*fup
         wfcatom_d(ig,2,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman)) &
                  -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fdown
      END DO
   END DO
   !
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
   DEALLOCATE( chiaux_d )
   !
  END SUBROUTINE atomic_wfc_so_mag_gpu
  !
  !
  SUBROUTINE atomic_wfc_nc_gpu( )
   !
   !! noncolinear case, magnetization along "angle1" and "angle2"
   !
   REAL(DP) :: alpha, gamman
   COMPLEX(DP) :: fup, fdown
   INTEGER :: m, lm, ig  
   !
   alpha = angle1(nt)
   gamman = - angle2(nt) + 0.5d0*pi
   !
   DO m = 1, 2*l+1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = ig_start, ig_end
         aux_d(ig) = sk_d(ig)*CMPLX(ylm_d(ig,lm)*chiq_d(ig,nb,nt))
      END DO
      !
      ! now, rotate wfc as needed
      ! first : rotation with angle alpha around (OX)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = ig_start, ig_end
         fup = CMPLX(COS(0.5d0*alpha))*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*alpha))*aux_d(ig)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom_d(ig,1,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman)) &
                        +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fup
         wfcatom_d(ig,2,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman)) &
                        -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = CMPLX(COS(0.5d0*(alpha+pi)))*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*(alpha+pi)))*aux_d(ig)
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom_d(ig,1,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman)) &
                  +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fup
         wfcatom_d(ig,2,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman)) &
                  -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman)))*fdown
      END DO
   END DO
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
  END SUBROUTINE atomic_wfc_nc_gpu
  !
  !
  SUBROUTINE atomic_wfc___gpu( )
   !
   ! ... LSDA or nonmagnetic case
   !
   INTEGER :: m, lm, ig
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc > natomwfc) CALL errore &
         ('atomic_wfc___', 'internal error: too many wfcs', 1)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = ig_start, ig_end
         wfcatom_d(ig,1,n_starting_wfc) = lphase * &
            sk_d(ig) * CMPLX(ylm_d(ig,lm) * chiq_d(ig,nb,nt))
      ENDDO
      !
   END DO
   !
  END SUBROUTINE atomic_wfc___gpu
  !
END SUBROUTINE atomic_wfc_gpu
