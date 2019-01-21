  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi,
  ! Feliciano Giustino
  !
  ! Copyright (C) 2001-2003 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from LR_Modules/newdq.f90 (QE)
  !
  !----------------------------------------------------------------------
  SUBROUTINE newdq2( dvscf, npe )
  !----------------------------------------------------------------------
  !!
  !! This routine computes the contribution of the selfconsistent
  !! change of the potential to the known part of the linear
  !! system and adds it to dvpsi.
  !!
  !! Roxana Margine - Jan 2019: Updated based on QE 6.3
  !!
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, mill, eigts1, eigts2, eigts3
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE lrus,                 ONLY : int3
  USE qpoint,               ONLY : xq, eigqts
  USE constants_epw,        ONLY : czero
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npe
  !! Number of perturbations for this irr representation
  COMPLEX(DP), INTENT(in) :: dvscf(dfftp%nnr, nspin_mag, npe)
  !! Change of the selfconsistent potential
  !
  !   Local variables
  !
  INTEGER :: na
  !! counter on atoms
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! counter on atomic types
  INTEGER ::  ir
  !! counter on real mesh
  INTEGER :: ipert
  !! counter on change of Vscf due to perturbations
  INTEGER :: is
  !! counter on spin
  INTEGER :: ih
  !! Counter on beta functions
  INTEGER :: jh
  !! Counter on beta functions
  !
  REAL(DP), ALLOCATABLE :: qmod(:)
  !! the modulus of q+G
  REAL(DP), ALLOCATABLE :: qg(:,:)
  !! the values of q+G
  REAL(DP), ALLOCATABLE :: ylmk0(:,:)
  !! the spherical harmonics
  !
  COMPLEX(DP), EXTERNAL :: zdotc
  !! the scalar product function
  COMPLEX(DP), ALLOCATABLE :: aux1(:), aux2(:,:), veff(:), qgm(:)
  !
  IF (.not.okvan) RETURN
  !
  CALL start_clock('newdq2')
  !
  int3(:,:,:,:,:) = czero
  !
  ALLOCATE( aux1(ngm) )
  ALLOCATE( aux2(ngm,nspin_mag) )
  ALLOCATE( veff(dfftp%nnr) )
  ALLOCATE( ylmk0(ngm, lmaxq * lmaxq) )
  ALLOCATE( qgm(ngm) )
  ALLOCATE( qmod(ngm) )
  ALLOCATE( qg(3,ngm) )
  !
  !    first compute the spherical harmonics
  !
  CALL setqmod( ngm, xq, g, qmod, qg )
  CALL ylmr2( lmaxq * lmaxq, ngm, qg, qmod, ylmk0 )
  DO ig = 1, ngm
    qmod(ig) = sqrt( qmod(ig) )
  ENDDO
  !
  !     and for each perturbation of this irreducible representation
  !     integrate the change of the self consistent potential and
  !     the Q functions
  !
  DO ipert = 1, npe
    !
    DO is = 1, nspin_mag
      DO ir = 1, dfftp%nnr
         veff(ir) = dvscf(ir,is,ipert)
      ENDDO
      CALL fwfft('Rho', veff, dfftp)
      DO ig = 1, ngm
         aux2(ig,is) = veff( dfftp%nl(ig) )
      ENDDO
    ENDDO
    !
    DO nt = 1, ntyp
      IF (upf(nt)%tvanp ) THEN
        DO ih = 1, nh(nt)
          DO jh = ih, nh(nt)
            CALL qvan2( ngm, ih, jh, nt, qmod, qgm, ylmk0 )
            DO na = 1, nat
              IF (ityp(na) == nt) THEN
                DO ig = 1, ngm
                  aux1(ig) = qgm(ig) * eigts1(mill(1,ig),na) * &
                                       eigts2(mill(2,ig),na) * &
                                       eigts3(mill(3,ig),na) * &
                                       eigqts(na)
                ENDDO
                DO is = 1, nspin_mag
                  int3(ih,jh,na,is,ipert) = omega * &
                                     zdotc(ngm,aux1,1,aux2(1,is),1)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DO na = 1, nat
          IF (ityp(na) == nt) THEN
            !
            !    We use the symmetry properties of the ps factor
            !
            DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                DO is = 1, nspin_mag
                  int3(jh,ih,na,is,ipert) = int3(ih,jh,na,is,ipert)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  !
  CALL mp_sum(int3, intra_pool_comm)
  !
  IF (noncolin) CALL set_int3_nc( npe )
  !
  DEALLOCATE (aux1)
  DEALLOCATE (aux2)
  DEALLOCATE (veff)
  DEALLOCATE (ylmk0)
  DEALLOCATE (qgm)
  DEALLOCATE (qmod)
  DEALLOCATE (qg)
  !
  CALL stop_clock('newdq2')
  !
  RETURN
  !
  END SUBROUTINE newdq2
