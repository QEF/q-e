  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! 
  ! Copyright (C) 2001-2003 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! adapted from PH/dvqpsi_us (QE)
  !
  !----------------------------------------------------------------------
  SUBROUTINE dvqpsi_us3( ik, uact, addnlcc, xxkq, xq0 )
  !----------------------------------------------------------------------
  !!
  !! This routine calculates dV_bare/dtau * psi for one perturbation
  !! with a given q. The displacements are described by a vector u.
  !! The result is stored in dvpsi. The routine is called for each k point
  !! and for each pattern u. It computes simultaneously all the bands.
  !! It implements Eq. B29 of PRB 64, 235118 (2001). The contribution
  !! of the local pseudopotential is calculated here, that of the nonlocal
  !! pseudopotential in dvqpsi_us_only3.
  !!
  !! RM - Nov/Dec 2014 
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !! Roxana Margine - Jan 2019: Updated based on QE 6.3
  !!
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat, ityp
  USE cell_base,             ONLY : tpiba
  USE fft_base,              ONLY : dfftp, dffts
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE gvect,                 ONLY : eigts1, eigts2, eigts3, mill, g, ngm
  USE gvecs,                 ONLY : ngms, doublegrid
  USE lsda_mod,              ONLY : lsda, isk
  USE scf,                   ONLY : rho, rho_core
  USE noncollin_module,      ONLY : nspin_lsda, nspin_gga, npol
  use uspp_param,            ONLY : upf
  USE wvfct,                 ONLY : npwx
  USE wavefunctions,         ONLY : evc
  USE nlcc_ph,               ONLY : drc
  USE uspp,                  ONLY : nlcc_any
  USE eqv,                   ONLY : dvpsi, dmuxc, vlocq
  USE qpoint,                ONLY : eigqts, npwq 
  USE klist,                 ONLY : ngk
  USE gc_lr,                 ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE funct,                 ONLY : dft_is_gradient, dft_is_nonlocc
  USE elph2,                 ONLY : igkq, igk, lower_band, upper_band
  USE constants_epw,         ONLY : czero, eps12
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: addnlcc
  !! True if NLCC is present
  !
  INTEGER, INTENT(in) :: ik
  !! Counter on k-point
  ! 
  REAL(kind=DP), INTENT (in) :: xq0(3)
  !! Current coarse q-point coordinate
  REAL(kind=DP), INTENT (in) :: xxkq(3)
  !! k+q point coordinate 
  ! 
  COMPLEX(kind=DP), INTENT(in) :: uact(3 * nat)
  !! the pattern of displacements
  !
  ! Local variables
  !
  INTEGER :: na
  !! counter on atoms
  INTEGER :: mu
  !! counter on modes
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! counter on atomic types
  INTEGER :: ibnd
  !! counter on bands
  INTEGER ::  ir
  !! counter on real mesh
  INTEGER :: is
  !! counter on spin
  INTEGER :: ip
  !! counter on polarizations
  INTEGER :: npw
  !! Number of k+G-vectors inside 'ecut sphere'
  !
  REAL(kind=DP) :: fac
  !! spin degeneracy factor
  !
  COMPLEX(kind=DP) :: gtau
  !! e^{-i G * \tau}
  COMPLEX(kind=DP) :: u1, u2, u3
  !! components of displacement pattern u 
  COMPLEX(kind=DP) :: gu0
  !! scalar product q * u
  COMPLEX(kind=DP) :: gu
  !! q * u + G * u
  COMPLEX(kind=DP) :: fact
  !! e^{-i q * \tau}
  COMPLEX(kind=DP), ALLOCATABLE, TARGET :: aux(:)
  COMPLEX(kind=DP), ALLOCATABLE :: aux1(:), aux2(:)
  COMPLEX(kind=DP), POINTER :: auxs(:)
  COMPLEX(kind=DP), ALLOCATABLE :: drhoc(:)
  !! response core charge density
  !
  CALL start_clock('dvqpsi_us3')
  !
  npw = ngk(ik)
  !
  IF (nlcc_any .AND. addnlcc) THEN
     ALLOCATE( drhoc(dfftp%nnr) )
     ALLOCATE( aux (dfftp%nnr) )
     ALLOCATE( auxs(dffts%nnr) )
  ENDIF
  ALLOCATE( aux1(dffts%nnr) )
  ALLOCATE( aux2(dffts%nnr) )
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  dvpsi(:,:) = czero
  aux1(:) = czero
  DO na = 1, nat
    fact = tpiba * (0.d0, -1.d0) * eigqts(na)
    mu = 3 * (na - 1)
    u1 = uact(mu+1)
    u2 = uact(mu+2)
    u3 = uact(mu+3)
    IF (abs(u1) + abs(u2) + abs(u3) .gt. eps12) THEN
      nt = ityp(na)
      gu0 = xq0(1) * u1 + xq0(2) * u2 + xq0(3) * u3
      DO ig = 1, ngms
        gtau = eigts1(mill(1,ig),na) * &
               eigts2(mill(2,ig),na) * & 
               eigts3(mill(3,ig),na)
        gu = gu0 + g(1,ig) * u1 + g(2,ig) * u2 + g(3,ig) * u3
        aux1(dffts%nl(ig)) = aux1(dffts%nl(ig)) + vlocq(ig,nt) * gu * fact * gtau
      ENDDO
    ENDIF
  ENDDO
  !
  ! add NLCC when present
  !
  IF (nlcc_any .AND. addnlcc) THEN
    drhoc(:) = czero
    DO na = 1, nat
      fact = tpiba * (0.d0, -1.d0) * eigqts(na)
      mu = 3 * (na - 1)
      u1 = uact(mu+1)
      u2 = uact(mu+2)
      u3 = uact(mu+3)
      IF (abs(u1) + abs(u2) + abs(u3) .gt. eps12) THEN
        nt = ityp(na)
        gu0 = xq0(1) * u1 + xq0(2) * u2 + xq0(3) * u3
        IF (upf(nt)%nlcc) THEN
          DO ig = 1,ngm
            gtau = eigts1(mill(1,ig),na) * &
                   eigts2(mill(2,ig),na) * &
                   eigts3(mill(3,ig),na)
            gu = gu0 + g(1,ig) * u1 + g(2,ig) * u2 + g(3,ig) * u3
            drhoc(dfftp%nl(ig)) = drhoc(dfftp%nl(ig)) + drc(ig,nt) * gu * fact * gtau
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    !
    CALL invfft('Rho', drhoc, dfftp)
    !
    IF (.not.lsda) THEN
      DO ir = 1, dfftp%nnr
        aux(ir) = drhoc(ir) * dmuxc(ir,1,1)
      ENDDO
    ELSE
      is = isk(ik)
      DO ir = 1, dfftp%nnr
        aux(ir) = drhoc(ir) * 0.5d0 * ( dmuxc(ir,is,1) + dmuxc(ir,is,2) )
      ENDDO
    ENDIF
    !
    fac = 1.d0 / dble(nspin_lsda)
    DO is = 1, nspin_lsda
      rho%of_r(:,is) = rho%of_r(:,is) + fac * rho_core
    ENDDO
    !
    IF ( dft_is_gradient() ) &
      CALL dgradcorr( dfftp, rho%of_r, grho, &
             dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq0, drhoc, &
             1, nspin_gga, g, aux )
    !
    IF ( dft_is_nonlocc() ) &
      CALL dnonloccorr( rho%of_r, drhoc, xq0, aux )
    !
    DO is = 1, nspin_lsda
      rho%of_r(:,is) = rho%of_r(:,is) - fac * rho_core
    ENDDO
    !
    CALL fwfft('Rho', aux, dfftp)
    !
    ! This is needed also when the smooth and the thick grids coincide to
    ! cut the potential at the cut-off
    ! 
    auxs(:) = czero
    DO ig = 1, ngms
      auxs(dffts%nl(ig)) = aux(dfftp%nl(ig))
    ENDDO
    aux1(:) = aux1(:) + auxs(:)
  ENDIF
  !
  ! Now we compute dV_loc/dtau in real space
  !
  CALL invfft('Rho', aux1, dffts)
  DO ibnd = lower_band, upper_band
    DO ip = 1, npol
      aux2(:) = czero
      IF ( ip == 1 ) THEN
        DO ig = 1, npw
          aux2( dffts%nl( igk(ig) ) ) = evc(ig,ibnd)
        ENDDO
      ELSE
        DO ig = 1, npw
          aux2( dffts%nl( igk(ig) ) ) = evc(ig+npwx,ibnd)
        ENDDO
      ENDIF
      !
      !  This wavefunction is computed in real space
      !
      CALL invfft('Wave', aux2, dffts)
      DO ir = 1, dffts%nnr
        aux2(ir) = aux2(ir) * aux1(ir)
      ENDDO
      !
      ! and finally dV_loc/dtau * psi is transformed in reciprocal space
      !
      CALL fwfft('Wave', aux2, dffts)
      IF ( ip == 1 ) THEN
        DO ig = 1, npwq
          dvpsi(ig,ibnd) = aux2( dffts%nl( igkq(ig) ) )
        ENDDO
      ELSE
        DO ig = 1, npwq
          dvpsi(ig+npwx,ibnd) = aux2( dffts%nl( igkq(ig) ) )
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  ! 
  IF (nlcc_any .AND. addnlcc) THEN
    DEALLOCATE(drhoc)
    DEALLOCATE(aux)
    DEALLOCATE(auxs)
  ENDIF
  DEALLOCATE(aux1)
  DEALLOCATE(aux2)
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients in the perturbat
  !
  CALL dvqpsi_us_only3( ik, uact, xxkq )
  !
  CALL stop_clock('dvqpsi_us3')
  !
  RETURN
  !
  END SUBROUTINE dvqpsi_us3
