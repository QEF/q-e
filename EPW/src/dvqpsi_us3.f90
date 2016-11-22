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
  SUBROUTINE dvqpsi_us3 (ik, uact, addnlcc, xxk, xq0)
  !----------------------------------------------------------------------
  !!
  !! This routine calculates dV_bare/dtau * psi for one perturbation
  !! with a given q. The displacements are described by a vector u.
  !! The result is stored in dvpsi. The routine is called for each k point
  !! and for each pattern u. It computes simultaneously all the bands.
  !!
  !! RM - Nov/Dec 2014 
  !! Imported the noncolinear case implemented by xlzhang
  !!
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat, ityp
  USE cell_base,             ONLY : tpiba
  USE fft_base,              ONLY : dfftp, dffts
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE gvect,                 ONLY : eigts1, eigts2, eigts3, mill, g, nl, &
                                    ngm
  USE gvecs,                 ONLY : ngms, doublegrid, nls
  USE lsda_mod,              ONLY : lsda, isk
  USE noncollin_module,      ONLY : npol
  use uspp_param,            ONLY : upf
  USE wvfct,                 ONLY : nbnd, npwx
  USE wavefunctions_module,  ONLY : evc
  USE nlcc_ph,               ONLY : drc
  USE uspp,                  ONLY : nlcc_any
  USE eqv,                   ONLY : dvpsi, dmuxc, vlocq
  USE qpoint,                ONLY : eigqts, npwq !, ikks
  USE klist,                 ONLY : ngk
  USE elph2,                 ONLY : igkq, igk, lower_band, upper_band

  implicit none
  !
  LOGICAL, INTENT (in) :: addnlcc
  !! True if NLCC is present
  !
  INTEGER, INTENT (in) :: ik
  !! Counter on k-point
  ! 
  REAL(kind=DP), INTENT (in) :: xq0(3)
  !! Current coarse q-point coordinate
  REAL(kind=DP), INTENT (in) :: xxk(3)
  !! k-point coordinate
  ! 
  COMPLEX(kind=DP), INTENT (in) :: uact (3 * nat)
  !! the pattern of displacements
  !
  ! Local variables
  integer :: na, mu, ig, nt, ibnd, ir, is, ip, npw
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh

  complex(kind=DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(kind=DP) , ALLOCATABLE, TARGET :: aux (:)
  complex(kind=DP) , ALLOCATABLE :: aux1 (:), aux2 (:)
  complex(kind=DP) , POINTER :: auxs (:)
  ! work space
  ! SP : seems to be useless ?
  !logical :: htg


  CALL start_clock ('dvqpsi_us')
  !IF (nlcc_any) THEN
  ! SP - changed according to QE5/PH/dvqpsi_us
  npw  = ngk(ik)

  IF (nlcc_any.AND.addnlcc) THEN
     ALLOCATE (aux( dfftp%nnr))
     IF (doublegrid) then
        ALLOCATE (auxs(dffts%nnr))
     ELSE
        auxs => aux
     ENDIF
  ENDIF
  ALLOCATE (aux1(dffts%nnr))
  ALLOCATE (aux2(dffts%nnr))
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  aux1(:) = (0.d0, 0.d0)
  DO na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     IF (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
          3) ) .gt.1.0d-12) THEN
        nt = ityp (na)
        u1 = uact (mu + 1)
        u2 = uact (mu + 2)
        u3 = uact (mu + 3)
        gu0 = xq0 (1) * u1 + xq0 (2) * u2 + xq0 (3) * u3
        DO ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * eigts3 ( &
                mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           aux1 (nls (ig) ) = aux1 (nls (ig) ) + vlocq (ig, nt) * gu * &
                fact * gtau
        ENDDO
     ENDIF
  ENDDO
  !
  ! add NLCC when present
  !
   IF (nlcc_any.and.addnlcc) THEN
      aux(:) = (0.d0, 0.d0)
      DO na = 1,nat
         fact = tpiba*(0.d0,-1.d0)*eigqts(na)
         mu = 3*(na-1)
         IF (abs(uact(mu+1))+abs(uact(mu+2))  &
                         +abs(uact(mu+3)).gt.1.0d-12) then
            nt=ityp(na)
            u1 = uact(mu+1)
            u2 = uact(mu+2)
            u3 = uact(mu+3)
            gu0 = xq0(1)*u1 +xq0(2)*u2+xq0(3)*u3
            IF (upf(nt)%nlcc) THEN
               DO ig = 1,ngm
                  gtau = eigts1(mill(1,ig),na)*   &
                         eigts2(mill(2,ig),na)*   &
                         eigts3(mill(3,ig),na)
                  gu = gu0+g(1,ig)*u1+g(2,ig)*u2+g(3,ig)*u3
                  aux(nl(ig))=aux(nl(ig))+drc(ig,nt)*gu*fact*gtau
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      CALL invfft ('Dense', aux, dfftp)
      IF (.not.lsda) THEN
         DO ir=1,dfftp%nnr
            aux(ir) = aux(ir) * dmuxc(ir,1,1)
         ENDDO
      ELSE
         ! RM - I think it should be ik to be consistent with dvqpsi_us_only3
         !is=isk(ikk)
         is=isk(ik)
         DO ir=1,dfftp%nnr
            aux(ir) = aux(ir) * 0.5d0 *  &
                 (dmuxc(ir,is,1)+dmuxc(ir,is,2))
         ENDDO
      ENDIF
      CALL fwfft ('Dense', aux, dfftp)
      IF (doublegrid) THEN
         auxs(:) = (0.d0, 0.d0)
         DO ig=1,ngms
            auxs(nls(ig)) = aux(nl(ig))
         ENDDO
      ENDIF
      aux1(:) = aux1(:) + auxs(:)
   ENDIF
  !
  ! Now we compute dV_loc/dtau in real space
  !
  CALL invfft ('Smooth', aux1, dffts)
  DO ibnd = lower_band, upper_band
     DO ip = 1, npol
        aux2(:) = (0.d0, 0.d0)
        IF ( ip == 1 ) THEN
           DO ig = 1, npw
              aux2 (nls (igk (ig) ) ) = evc (ig, ibnd)
           ENDDO
        ELSE
           DO ig = 1, npw
              aux2 (nls (igk (ig) ) ) = evc (ig+npwx, ibnd)
           ENDDO
        ENDIF
        !
        !  This wavefunction is computed in real space
        !
        CALL invfft ('Wave', aux2, dffts)
        DO ir = 1, dffts%nnr
           aux2 (ir) = aux2 (ir) * aux1 (ir)
        ENDDO
        !
        ! and finally dV_loc/dtau * psi is transformed in reciprocal space
        !
        CALL fwfft ('Wave', aux2, dffts)
        IF ( ip == 1 ) THEN
           DO ig = 1, npwq
              dvpsi (ig, ibnd) = aux2 (nls (igkq (ig) ) )
           ENDDO
        ELSE
           DO ig = 1, npwq
              dvpsi (ig+npwx, ibnd) = aux2 (nls (igkq (ig) ) )
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  ! 
  DEALLOCATE (aux2)
  DEALLOCATE (aux1)
  IF (nlcc_any.and.addnlcc) THEN
     DEALLOCATE (aux)
     IF (doublegrid) DEALLOCATE (auxs)
  ENDIF
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients in the perturbat
  !
  CALL dvqpsi_us_only3 (ik, uact, xxk)
  CALL stop_clock ('dvqpsi_us')
  RETURN
END SUBROUTINE dvqpsi_us3
