!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dfunct

CONTAINS
!---------------------------------------

SUBROUTINE newq(vr,deeq,skip_vltot)
  !
  !   This routine computes the integral of the perturbed potential with
  !   the Q function
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, gstart, mill, &
                                   eigts1, eigts2, eigts3, nl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vltot
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE spin_orb,             ONLY : lspinorb, domag
  USE noncollin_module,     ONLY : nspin_mag
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  !
  ! Input: potential , output: contribution to integral
  REAL(kind=dp), intent(in)  :: vr(dfftp%nnr,nspin)
  REAL(kind=dp), intent(out) :: deeq( nhm, nhm, nat, nspin )
  LOGICAL, intent(in) :: skip_vltot !If .false. vltot is added to vr when necessary
  ! INTERNAL
  INTEGER :: ig, nt, ih, jh, na, is, ijh, nij, nb, nab
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:), aux(:,:), qgm(:,:)
    ! work space
  REAL(DP), ALLOCATABLE :: ylmk0(:,:), qmod(:), deeaux(:,:)
    ! spherical harmonics, modulus of G
  REAL(DP) :: fact
  !
  IF ( gamma_only ) THEN
     fact = 2.0_dp
  ELSE
     fact = 1.0_dp
  END IF
  !
  deeq(:,:,:,:) = 0.D0
  !
  ALLOCATE( vaux(ngm,nspin_mag), qmod( ngm ), ylmk0( ngm, lmaxq*lmaxq ) )
  !
  CALL ylmr2( lmaxq * lmaxq, ngm, g, gg, ylmk0 )
  qmod(1:ngm) = SQRT( gg(1:ngm) )
  !
  ! ... fourier transform of the total effective potential
  !
  DO is = 1, nspin_mag
     !
     IF ( (nspin_mag == 4 .AND. is /= 1) .or. skip_vltot ) THEN 
!$omp parallel do default(shared) private(ig)
        do ig=1,dfftp%nnr
           psic(ig) = vr(ig,is)
        end do
!$omp end parallel do
     ELSE
!$omp parallel do default(shared) private(ig)
        do ig=1,dfftp%nnr
           psic(ig) = vltot(ig) + vr(ig,is)
        end do
!$omp end parallel do
     END IF
     CALL fwfft ('Dense', psic, dfftp)
!$omp parallel do default(shared) private(ig)
        do ig=1,ngm
           vaux(ig, is) = psic(nl(ig))
        end do
!$omp end parallel do
     !
  END DO

  DO nt = 1, ntyp
     !
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        ALLOCATE ( qgm(ngm,nij) )
        !
        ! ... Compute and store Q(G) for this atomic species 
        ! ... (without structure factor)
        !
        ijh = 0
        DO ih = 1, nh(nt)
           DO jh = ih, nh(nt)
              ijh = ijh + 1
              CALL qvan2 ( ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
           END DO
        END DO
        !
        ! count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        END DO
        ALLOCATE ( aux (ngm, nab ), deeaux(nij, nab) )
        !
        ! ... Compute and store V(G) times the structure factor e^(-iG*tau)
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF ( ityp(na) == nt ) THEN
                 nb = nb + 1
!$omp parallel do default(shared) private(ig)
                 do ig=1,ngm
                    aux(ig, nb) = vaux(ig,is) * CONJG ( &
                      eigts1(mill(1,ig),na) * &
                      eigts2(mill(2,ig),na) * &
                      eigts3(mill(3,ig),na) )
                 end do
!$omp end parallel do
              END IF
           END DO
           !
           ! ... here we compute the integral Q*V for all atoms of this kind
           !
           CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, aux, &
                    2*ngm, 0.0_dp, deeaux, nij )
           IF ( gamma_only .AND. gstart == 2 ) &
                CALL DGER(nij, nab,-1.0_dp, qgm, 2*ngm, aux, 2*ngm, deeaux, nij)
           !
           nb = 0
           DO na = 1, nat
              IF ( ityp(na) == nt ) THEN
                 nb = nb + 1
                 ijh = 0
                 DO ih = 1, nh(nt)
                    DO jh = ih, nh(nt)
                       ijh = ijh + 1
                       deeq(ih,jh,na,is) = omega * deeaux(ijh,nb)
                       if (jh > ih) deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
                    END DO
                 END DO
              END IF
           END DO
           !
        END DO
        !
        DEALLOCATE ( deeaux, aux, qgm )
        !
     END IF
     !
  END DO
  !
  DEALLOCATE( qmod, ylmk0, vaux )
  CALL mp_sum( deeq( :, :, :, 1:nspin_mag ), intra_bgrp_comm )
  !
END SUBROUTINE newq
  !
!----------------------------------------------------------------------------
SUBROUTINE newd( ) 
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the integral of the effective potential with
  ! ... the Q function and adds it to the bare ionic D term which is used
  ! ... to compute the non-local term in the US scheme.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : deeq, dvan, deeq_nc, dvan_so, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE spin_orb,             ONLY : lspinorb, domag
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : nhtol, nhtolm
  USE scf,                  ONLY : v
  USE realus,        ONLY : newq_r
  USE control_flags, ONLY : tqr
  USE ldaU,          ONLY : lda_plus_U, U_projection
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht, nb, mb
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  !
  !
  IF ( .NOT. okvan ) THEN
     !
     ! ... no ultrasoft potentials: use bare coefficients for projectors
     !
     DO na = 1, nat
        !
        nt  = ityp(na)
        nht = nh(nt)
        !
        IF ( lspinorb ) THEN
           !
           deeq_nc(1:nht,1:nht,na,1:nspin) = dvan_so(1:nht,1:nht,1:nspin,nt)
           !
        ELSE IF ( noncolin ) THEN
           !
           deeq_nc(1:nht,1:nht,na,1) = dvan(1:nht,1:nht,nt)
           deeq_nc(1:nht,1:nht,na,2) = ( 0.D0, 0.D0 )
           deeq_nc(1:nht,1:nht,na,3) = ( 0.D0, 0.D0 )
           deeq_nc(1:nht,1:nht,na,4) = dvan(1:nht,1:nht,nt)
           !
        ELSE
           !
           DO is = 1, nspin
              !
              deeq(1:nht,1:nht,na,is) = dvan(1:nht,1:nht,nt)
              !
           END DO
           !
        END IF
        !
     END DO
     !
     ! ... early return
     !
     RETURN
     !
  END IF
  !
  CALL start_clock( 'newd' )
  !
  IF (tqr) THEN
     CALL newq_r(v%of_r,deeq,.false.)
  ELSE
     CALL newq(v%of_r,deeq,.false.)
  END IF
  !
  IF (noncolin) call add_paw_to_deeq(deeq)
  !
  atoms : &
  DO na = 1, nat
     !
     nt  = ityp(na)
     if_noncolin:&
     IF ( noncolin ) THEN
        !
        IF (upf(nt)%has_so) THEN
           !
           CALL newd_so(na)
           !
        ELSE
           !
           CALL newd_nc(na)
           !
        END IF
        !
     ELSE if_noncolin
        !
        DO is = 1, nspin
           !
           DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                 deeq(ih,jh,na,is) = deeq(ih,jh,na,is) + dvan(ih,jh,nt)
                 deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
              END DO
           END DO
           !
        END DO
        !
     END IF if_noncolin
     !
  END DO atoms
  !
  IF (.NOT.noncolin) CALL add_paw_to_deeq(deeq)
  !
  IF (lda_plus_U .AND. (U_projection == 'pseudo')) CALL add_vhub_to_deeq(deeq)
  !
  CALL stop_clock( 'newd' )
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_so(na)
      !------------------------------------------------------------------------
      !
      USE spin_orb, ONLY : fcoef
      !
      IMPLICIT NONE
      !
      INTEGER :: na

      INTEGER :: ijs, is1, is2, kh, lh
      !
      !
      nt=ityp(na)
      ijs = 0
      !
      DO is1 = 1, 2
         !
         DO is2 =1, 2
            !
            ijs = ijs + 1
            !
            IF (domag) THEN
               DO ih = 1, nh(nt)
                  !
                  DO jh = 1, nh(nt)
                     !
                     deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                     !
                     DO kh = 1, nh(nt)
                        !
                        DO lh = 1, nh(nt)
                           !
                           deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                deeq (kh,lh,na,1)*            &
                             (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  + &
                             fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt)) + &
                             deeq (kh,lh,na,2)*            &
                             (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  + &
                             fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                             (0.D0,-1.D0)*deeq (kh,lh,na,3)*            &
                             (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  - &
                             fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                             deeq (kh,lh,na,4)*            &
                             (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  - &
                             fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt))   
                           !
                        END DO
                        !
                     END DO
                     !
                  END DO
                  !
               END DO
               !
            ELSE
               !
               DO ih = 1, nh(nt)
                  !
                  DO jh = 1, nh(nt)
                     !
                     deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                     !
                     DO kh = 1, nh(nt)
                        !
                        DO lh = 1, nh(nt)
                           !
                           deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                deeq (kh,lh,na,1)*            &
                             (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  + &
                             fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt) ) 
                           !
                        END DO
                        !
                     END DO
                     !
                  END DO
                  !
               END DO
               !
            END IF
            !
         END DO
         !
      END DO
      !
    RETURN
      !
    END SUBROUTINE newd_so
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_nc(na)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER :: na
      !
      nt = ityp(na)
      !
      DO ih = 1, nh(nt)
         !
         DO jh = 1, nh(nt)
            !
            IF (lspinorb) THEN
               deeq_nc(ih,jh,na,1) = dvan_so(ih,jh,1,nt) + &
                                     deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
               !                      
               deeq_nc(ih,jh,na,4) = dvan_so(ih,jh,4,nt) + &
                                     deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
               !
            ELSE
               deeq_nc(ih,jh,na,1) = dvan(ih,jh,nt) + &
                                     deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
               !                      
               deeq_nc(ih,jh,na,4) = dvan(ih,jh,nt) + &
                                     deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
               !
            END IF
            deeq_nc(ih,jh,na,2) = deeq(ih,jh,na,2) - &
                                  ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
            !                      
            deeq_nc(ih,jh,na,3) = deeq(ih,jh,na,2) + &
                                  ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
            !                      
         END DO
         !
      END DO
      !
    RETURN
    END SUBROUTINE newd_nc
    !
END SUBROUTINE newd

END MODULE dfunct
