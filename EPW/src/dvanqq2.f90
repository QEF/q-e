  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from PH/dvanqq (QE)
  ! 
  !----------------------------------------------------------------------
  SUBROUTINE dvanqq2
  !----------------------------------------------------------------------
  !!
  !! New
  !! This routine calculates two integrals of the Q functions and
  !! its derivatives with V_loc and V_eff which are used
  !! to compute term dV_bare/dtau * psi  in addusdvqpsi.
  !! The result is stored in int1, int2, int4, int5. The routine is called
  !! for each q in nqc. 
  !! int1 -> Eq. B20 of Ref.[1]
  !! int2 -> Eq. B21 of Ref.[1]
  !! int4 -> Eq. B23 of Ref.[1]
  !! int5 -> Eq. B24 of Ref.[1]
  !!
  !! [1] PRB 64, 235118 (2001).
  !! 
  !! RM - Nov/Dec 2014 
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !! Roxana Margine - Dec 2018: Updated based on QE 6.3
  !!
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE spin_orb,         ONLY : lspinorb
  USE cell_base,        ONLY : tpiba2, omega, tpiba
  USE gvect,            ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE scf,              ONLY : v, vltot
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE phcom,            ONLY : int1, int2, int4, int4_nc, int5, int5_so, & 
                               vlocq
  USE qpoint,           ONLY : xq, eigqts
  USE uspp_param,       ONLY : upf, lmaxq, nh
  USE uspp,             ONLY : okvan, ijtoh
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE fft_interfaces,   ONLY : fwfft
  USE constants_epw,    ONLY : zero, czero
  !
  IMPLICIT NONE
  !
  !   Local variables
  !
  INTEGER :: na
  !! counter on atoms
  INTEGER :: nb
  !! counter on atoms  
  INTEGER :: ntb
  !! counter on atomic types (species)
  INTEGER :: nta
  !! index of atomic type (specie)
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: ir
  !! counter on FFT mesh points
  INTEGER :: ih
  !! counter on beta functions per atomic type
  INTEGER :: jh
  !! counter on beta functions per atomic type
  INTEGER :: ijh
  !! correspondence beta indexes ih,jh -> composite index ijh
  INTEGER :: ipol
  !! counter on polarizations
  INTEGER :: jpol
  !! counter on polarizations
  INTEGER :: is
  !! counter on spin
  ! 
  REAL(kind=DP), ALLOCATABLE :: qmod(:)
  !! the modulus of q+G
  REAL(kind=DP), ALLOCATABLE :: qmodg(:)
  !! the modulus of G
  REAL(DP), ALLOCATABLE :: qpg(:,:)
  !! the q+G vectors
  REAL(kind=DP), ALLOCATABLE :: ylmkq(:,:)
  !! the spherical harmonics at q+G
  REAL(kind=DP), ALLOCATABLE ::  ylmk0(:,:)
  !! the spherical harmonics at G
  !
  COMPLEX(kind=DP) :: fact
  !! e^{-i q * \tau} * conjg(e^{-i q * \tau}) 
  COMPLEX(kind=DP) :: fact1
  !! -i * omega
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  !! the scalar product function
  COMPLEX(kind=DP), ALLOCATABLE :: aux1(:), aux2(:), &
       aux3(:), aux5(:), sk(:)
  COMPLEX(kind=DP), ALLOCATABLE :: veff(:,:)
  !! effective potential
  COMPLEX(kind=DP), ALLOCATABLE, TARGET :: qgm(:)
  !! the augmentation function at G
  COMPLEX(kind=DP), POINTER :: qgmq(:)
  !! the augmentation function at q+G
  ! 
  IF (.not.okvan) RETURN
  !
  CALL start_clock('dvanqq2')
  ! 
  int1(:,:,:,:,:) = czero
  int2(:,:,:,:,:) = czero
  int4(:,:,:,:,:) = czero
  int5(:,:,:,:,:) = czero
  ALLOCATE( sk(ngm) )    
  ALLOCATE( aux1(ngm) )    
  ALLOCATE( aux2(ngm) )    
  ALLOCATE( aux3(ngm) )    
  ALLOCATE( aux5(ngm) )    
  ALLOCATE( qmodg(ngm) )    
  ALLOCATE( qmod(ngm) )
  ALLOCATE( qgmq(ngm) )
  ALLOCATE( qgm(ngm))
  ALLOCATE( ylmk0(ngm, lmaxq * lmaxq) )    
  ALLOCATE( ylmkq(ngm, lmaxq * lmaxq) )    
  !
  ! compute spherical harmonics
  !
  CALL ylmr2( lmaxq * lmaxq, ngm, g, gg, ylmk0 )
  DO ig = 1, ngm
    qmodg(ig) = sqrt( gg(ig) )
  ENDDO
  ! 
  ALLOCATE( qpg(3, ngm) )    
  CALL setqmod( ngm, xq, g, qmod, qpg )
  CALL ylmr2(lmaxq * lmaxq, ngm, qpg, qmod, ylmkq)
  DEALLOCATE(qpg)
  DO ig = 1, ngm
    qmod(ig) = sqrt( qmod(ig) )
  ENDDO
  !
  !   we start by computing the FT of the effective potential
  !
  ALLOCATE (veff(dfftp%nnr,nspin_mag))    
  DO is = 1, nspin_mag
    IF (nspin_mag.ne.4 .or. is==1) THEN
      DO ir = 1, dfftp%nnr
        veff(ir,is) = CMPLX(vltot(ir) + v%of_r(ir,is), zero, kind=DP)
      ENDDO
    ELSE
      DO ir = 1, dfftp%nnr
        veff(ir,is) = CMPLX(v%of_r(ir,is), zero, kind=DP)
      ENDDO
    ENDIF
    CALL fwfft('Rho', veff(:,is), dfftp)
  ENDDO
  !
  !    We compute here two of the three integrals needed in the phonon
  !
  fact1 = CMPLX(0.d0, - tpiba * omega, kind=DP)
  !
  DO ntb = 1, ntyp
    IF (upf(ntb)%tvanp ) THEN
      !
      DO ih = 1, nh(ntb)
        DO jh = ih, nh(ntb)
          ijh = ijtoh(ih,jh,ntb)
          !
          !    compute the augmentation function
          !
          CALL qvan2( ngm, ih, jh, ntb, qmodg, qgm, ylmk0 )
          CALL qvan2( ngm, ih, jh, ntb, qmod, qgmq, ylmkq )
          !
          !     NB: for this integral the moving atom and the atom of Q
          !     do not necessarily coincide
          !
          DO nb = 1, nat
            IF (ityp(nb) == ntb) THEN
              DO ig = 1, ngm
                aux1(ig) = qgmq(ig) * eigts1(mill(1,ig),nb) &
                                    * eigts2(mill(2,ig),nb) &
                                    * eigts3(mill(3,ig),nb)
              ENDDO
              !
              DO na = 1, nat
                fact = eigqts(na) * conjg( eigqts(nb) )
                !
                !    nb is the atom of the augmentation function
                !
                nta = ityp(na)
                DO ig = 1, ngm
                   sk(ig) = vlocq(ig,nta) * eigts1(mill(1,ig),na) &
                                          * eigts2(mill(2,ig),na) &
                                          * eigts3(mill(3,ig),na) 
                ENDDO
                !
                DO ipol = 1, 3
                   DO ig = 1, ngm
                     aux5(ig) = sk(ig) * ( g(ipol,ig) + xq(ipol) )
                   ENDDO
                   int2(ih,jh,ipol,na,nb) = fact * fact1 * &
                         zdotc(ngm, aux1, 1, aux5, 1)
                   ! 
                   DO jpol = 1, 3
                      IF (jpol >= ipol) THEN
                         DO ig = 1, ngm
                            aux3(ig) = aux5(ig) * &
                                     ( g(jpol,ig) + xq(jpol) )
                         ENDDO
                         int5(ijh,ipol,jpol,na,nb) = &
                             conjg(fact) * tpiba2 * omega * &
                             zdotc(ngm, aux3, 1, aux1, 1)
                      ELSE
                         int5(ijh,ipol,jpol,na,nb) = &
                             int5(ijh,jpol,ipol,na,nb)
                      ENDIF
                   ENDDO
                ENDDO !ipol
                !
              ENDDO !na
              !
              DO ig = 1, ngm
                 aux1(ig) = qgm(ig) * eigts1(mill(1,ig),nb) &
                                    * eigts2(mill(2,ig),nb) &
                                    * eigts3(mill(3,ig),nb)
              ENDDO
              !
              DO is = 1, nspin_mag
                DO ipol = 1, 3
                  DO ig = 1, ngm
                    aux2(ig) = veff(dfftp%nl(ig),is) * g(ipol,ig)
                  ENDDO
                  int1(ih,jh,ipol,nb,is) = - fact1 * &
                     zdotc(ngm, aux1, 1, aux2, 1)
                  DO jpol = 1, 3
                    IF (jpol >= ipol) THEN
                      DO ig = 1, ngm
                         aux3(ig) = aux2(ig) * g(jpol,ig)
                      ENDDO
                      int4(ijh,ipol,jpol,nb,is) = - tpiba2 * &
                          omega * zdotc(ngm, aux3, 1, aux1, 1)
                    ELSE
                      int4(ijh,ipol,jpol,nb,is) = &
                          int4(ijh,jpol,ipol,nb,is)
                    ENDIF
                  ENDDO ! jpol
                ENDDO ! ipol
              ENDDO ! is
              !
            ENDIF ! ityp
          ENDDO ! nb
          !
        ENDDO ! jh
      ENDDO ! ih
      !
      DO ih = 1, nh(ntb)
        DO jh = ih + 1, nh(ntb)
          !
          !    We use the symmetry properties of the integral factor
          !
          DO nb = 1, nat
            IF (ityp(nb) == ntb) THEN
              DO ipol = 1, 3
                DO is = 1, nspin_mag
                  int1(jh,ih,ipol,nb,is) = int1(ih,jh,ipol,nb,is)
                ENDDO
                DO na = 1, nat
                  int2(jh,ih,ipol,na,nb) = int2(ih,jh,ipol,na,nb)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          !
        ENDDO ! jh
      ENDDO ! ih
      !
    ENDIF ! upf
  ENDDO ! ntb
  CALL mp_sum(int1, intra_pool_comm)
  CALL mp_sum(int2, intra_pool_comm)
  CALL mp_sum(int4, intra_pool_comm)
  CALL mp_sum(int5, intra_pool_comm)
  !
  IF (noncolin) THEN
    CALL set_int12_nc(0)
    int4_nc = czero
    IF (lspinorb) int5_so = czero
    DO ntb = 1, ntyp
      IF ( upf(ntb)%tvanp ) THEN
        DO na = 1, nat
          IF (ityp(na) == ntb) THEN
            IF (upf(ntb)%has_so)  THEN
              CALL transform_int4_so(int4,na)
              CALL transform_int5_so(int5,na)
            ELSE
              CALL transform_int4_nc(int4,na)
              IF (lspinorb) CALL transform_int5_nc(int5,na)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !
!DBRM
  !write(*,'(a,e20.12)') 'int1 = ', &
  !SUM((REAL(REAL(int1(:,:,:,:,:))))**2)+SUM((REAL(AIMAG(int1(:,:,:,:,:))))**2)
  !write(*,'(a,e20.12)') 'int2 = ', &
  !SUM((REAL(REAL(int2(:,:,:,:,:))))**2)+SUM((REAL(AIMAG(int2(:,:,:,:,:))))**2)
  !write(*,'(a,e20.12)') 'int4 = ', &
  !SUM((REAL(REAL(int4(:,:,:,:,:))))**2)+SUM((REAL(AIMAG(int4(:,:,:,:,:))))**2)
  !write(*,'(a,e20.12)') 'int5 = ', &
  !SUM((REAL(REAL(int5(:,:,:,:,:))))**2)+SUM((REAL(AIMAG(int5(:,:,:,:,:))))**2)
!END
  !
  DEALLOCATE(sk)
  DEALLOCATE(aux1)
  DEALLOCATE(aux2)
  DEALLOCATE(aux3)
  DEALLOCATE(aux5)
  DEALLOCATE(qmodg)
  DEALLOCATE(qmod)
  DEALLOCATE(qgmq)
  DEALLOCATE(qgm)
  DEALLOCATE(ylmk0)
  DEALLOCATE(ylmkq)
  DEALLOCATE(veff)
  !
  CALL stop_clock ('dvanqq2')
  RETURN
  ! 
  END SUBROUTINE dvanqq2
