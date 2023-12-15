!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
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
  USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,            ONLY : xk, ngk, igk_k_d !, igk_k
  USE wvfct,            ONLY : npwx
  USE uspp_param,       ONLY : upf, nwfcm
  USE noncollin_module, ONLY : noncolin, domag, npol, angle1, angle2, &
                               starting_spin_angle
  USE upf_spinorb,      ONLY : rot_ylm, lmaxx
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
             i0, i1, i2, i3, npw
  COMPLEX(DP) :: kphase, lphase
  REAL(DP)    :: arg, px, ux, vx, wx
  !
  REAL(DP) :: xk1, xk2, xk3, qgr
  REAL(DP), ALLOCATABLE :: chiq(:,:,:), qg(:)
  REAL(DP), ALLOCATABLE :: ylm_d(:,:), gk_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk_d(:), aux_d(:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfcatom_d, ylm_d, gk_d, sk_d, aux_d
#endif  
  !
  CALL start_clock( 'atomic_wfc' )

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, ntyp
     lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  END DO
  !
  npw = ngk(ik)
  !
  ALLOCATE( ylm_d(npw,(lmax_wfc+1)**2), gk_d(3,npw), sk_d(npw) ) 
  ALLOCATE( chiq(npw,nwfcm,ntyp), qg(npw) )
  !
  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)
  !
  !$acc parallel loop
  DO ig = 1, npw
     iig = igk_k_d(ig,ik)
     gk_d(1,ig) = xk1 + g(1,iig)
     gk_d(2,ig) = xk2 + g(2,iig)
     gk_d(3,ig) = xk3 + g(3,iig)
  END DO
  !
  !  ylm = spherical harmonics
  !
  !$acc data create(chiq) present(eigts1,eigts2,eigts3,mill)
  !$acc data create(qg)
  !$acc parallel loop
  DO ig = 1, npw
     qg(ig) = gk_d(1,ig)**2 +  gk_d(2,ig)**2 + gk_d(3,ig)**2
  END DO
  !$acc host_data use_device(qg)
  CALL ylmr2_gpu( (lmax_wfc+1)**2, npw, gk_d, qg, ylm_d )
  !$acc end host_data
  !
  ! set now q=|k+G| in atomic units
  !
  !$acc parallel loop
  DO ig = 1, npw
     qg(ig) = SQRT( qg(ig) )*tpiba
  END DO
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  CALL interp_atwfc ( npw, qg, nwfcm, chiq )
  !
  !$acc end data
  DEALLOCATE( qg, gk_d )
  !
  IF (noncolin) ALLOCATE( aux_d(npw) )
  !
  wfcatom_d(:,:,:) = (0.0_dp, 0.0_dp)
  !
  !$acc host_data use_device(chiq)
  DO na = 1, nat
     arg = (xk1*tau(1,na) + xk2*tau(2,na) + xk3*tau(3,na)) * tpi
     kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
     !
     !     sk is the structure factor
     !
     !$acc parallel loop
     DO ig = 1, npw
        iig = igk_k_d(ig,ik)
        sk_d(ig) = kphase * eigts1(mill(1,iig),na) * &
                            eigts2(mill(2,iig),na) * &
                            eigts3(mill(3,iig),na)
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
                    CALL atomic_wfc_so_gpu( chiq )
                 ELSE
                    CALL atomic_wfc_so_mag_gpu( chiq )
                 END IF
                 !
              ELSE
                 !
                 CALL atomic_wfc_nc_gpu( chiq )
                 !
              END IF
              !
           ELSE
              !
              CALL atomic_wfc___gpu( chiq )
              !
           END IF
           !
        END IF
        !
     END DO
     !
  END DO
  !$acc end host_data

  IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
       'internal error: some wfcs were lost ', 1 )

  IF (noncolin) DEALLOCATE( aux_d )
  DEALLOCATE( sk_d, ylm_d )
  !$acc end data
  DEALLOCATE (chiq)
  
  CALL stop_clock( 'atomic_wfc' )
  RETURN

CONTAINS
!----------------------------------------------------------------
  SUBROUTINE atomic_wfc_so_gpu( chiq_d )
   !------------------------------------------------------------
   !! Spin-orbit case.
   !
   REAL(DP) :: fact(2), fact_is, j
   COMPLEX(DP) :: rot_ylm_in1
   REAL(DP), EXTERNAL :: spinor
   INTEGER :: ind, ind1, n1, is, sph_ind
   REAL(DP) :: chiq_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: chiq_d
#endif  
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
                    DO ig = 1, npw
                      aux_d(ig) = aux_d(ig) + rot_ylm_in1 * &
                                  CMPLX(ylm_d(ig,ind1), KIND=DP)
                    ENDDO
                  ENDIF
               ENDDO
               !$cuf kernel do (1) <<<*,*>>>
               DO ig = 1, npw
                  wfcatom_d(ig,is,n_starting_wfc) = lphase * &
                                sk_d(ig)*aux_d(ig)*CMPLX(fact_is* &
                                chiq_d(ig,nb,nt), KIND=DP)
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
  SUBROUTINE atomic_wfc_so_mag_gpu( chiq_d )
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
   REAL(DP) :: chiq_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: chiq_d
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
      DO ig = 1, npw
        aux_d(ig) = sk_d(ig)* CMPLX(ylm_d(ig,lm)*chiaux_d(ig), KIND=DP)
      END DO
      !
      ! now, rotate wfc as needed
      ! first : rotation with angle alpha around (OX)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = 1, npw
         fup = CMPLX(COS(0.5d0*alpha), KIND=DP)*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*alpha), KIND=DP)*aux_d(ig)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom_d(ig,1,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fup
         wfcatom_d(ig,2,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = CMPLX(COS(0.5d0*(alpha+pi)), KIND=DP)*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*(alpha+pi)))*aux_d(ig)
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom_d(ig,1,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  +(0.d0,1.d0)*CMPLX(SIN(0.5d0 *gamman), KIND=DP))*fup
         wfcatom_d(ig,2,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
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
  SUBROUTINE atomic_wfc_nc_gpu( chiq_d )
   !
   !! noncolinear case, magnetization along "angle1" and "angle2"
   !
   REAL(DP) :: alpha, gamman
   COMPLEX(DP) :: fup, fdown
   INTEGER :: m, lm, ig  
   REAL(DP) :: chiq_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: chiq_d
#endif
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
      DO ig = 1, npw
         aux_d(ig) = sk_d(ig)*CMPLX(ylm_d(ig,lm)*chiq_d(ig,nb,nt), KIND=DP)
      END DO
      !
      ! now, rotate wfc as needed
      ! first : rotation with angle alpha around (OX)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = 1, npw
         fup = CMPLX(COS(0.5d0*alpha), KIND=DP)*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*alpha), KIND=DP)*aux_d(ig)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom_d(ig,1,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fup
         wfcatom_d(ig,2,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = CMPLX(COS(0.5d0*(alpha+pi)), KIND=DP)*aux_d(ig)
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*(alpha+pi)), KIND=DP)*aux_d(ig)
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom_d(ig,1,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fup
         wfcatom_d(ig,2,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
      END DO
   END DO
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
  END SUBROUTINE atomic_wfc_nc_gpu
  !
  !
  SUBROUTINE atomic_wfc___gpu( chiq_d )
    !
    ! ... LSDA or nonmagnetic case
    !
    REAL(dp) :: chiq_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: chiq_d
#endif

   INTEGER :: m, lm, ig
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc > natomwfc) CALL errore &
         ('atomic_wfc___', 'internal error: too many wfcs', 1)
      !
      !$cuf kernel do (1) <<<*,*>>>
      DO ig = 1, npw
         wfcatom_d(ig,1,n_starting_wfc) = lphase * &
            sk_d(ig) * CMPLX(ylm_d(ig,lm) * chiq_d(ig,nb,nt), KIND=DP)
      ENDDO
      !
   END DO
   !
  END SUBROUTINE atomic_wfc___gpu
  !
END SUBROUTINE atomic_wfc_gpu
