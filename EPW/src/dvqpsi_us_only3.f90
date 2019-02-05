  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! 
  ! Copyright (C) 2001-2008 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! adapted from PH/dvqpsi_us_only (QE)
  !
  !----------------------------------------------------------------------
  subroutine dvqpsi_us_only3( ik, uact, xxkq )
  !----------------------------------------------------------------------
  !!
  !! This routine calculates dV_bare/dtau * psi for one perturbation
  !! with a given q. The displacements are described by a vector uact.
  !! The result is stored in dvpsi. The routine is called for each k point
  !! and for each pattern u. It computes simultaneously all the bands.
  !! This routine implements Eq. B29 of PRB 64, 235118 (2001).
  !! Only the contribution of the nonlocal potential is calculated here.
  !!
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : tpiba
  USE gvect,      ONLY : g
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,   ONLY : lsda, current_spin, isk, nspin
  USE spin_orb,   ONLY : lspinorb
  USE wvfct,      ONLY : npwx, et
  USE uspp,       ONLY : okvan, nkb, vkb
  USE uspp_param, ONLY : nh, nhm
  USE qpoint,     ONLY : npwq
  USE phus,       ONLY : int1, int1_nc, int2, int2_so, alphap
  USE lrus,       ONLY : becp1
  USE eqv,        ONLY : dvpsi
  USE elph2,      ONLY : igkq, lower_band, upper_band
  USE noncollin_module, ONLY : noncolin, npol
  USE constants_epw,    ONLY : czero, cone, eps12
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ik
  !! the k point
  REAL(kind=DP), INTENT(in) :: xxkq(3) 
  !! the k+q point (cartesian coordinates)
  COMPLEX(kind=DP), INTENT(in) :: uact(3 * nat)
  !! the pattern of displacements
  !
  ! Local variables
  !
  INTEGER :: na
  !! Counter on atoms
  INTEGER :: nb
  !! Counter on atoms
  INTEGER :: mu
  !! Counter on modes
  INTEGER :: nu
  !! Counter on modes
  INTEGER :: ig
  !! Counter on G vectors
  INTEGER :: igg
  !! Auxiliary counter on G vectors
  INTEGER :: nt
  !! Counter on atomic types
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: ijkb0
  !! Auxiliary variable for counting
  INTEGER :: ikb
  !! Counter on becp functions
  INTEGER :: jkb
  !! Counter on becp functions
  INTEGER :: ipol
  !! Counter on polarizations
  INTEGER :: ih
  !! Counter on beta functions
  INTEGER :: jh
  !! Counter on beta functions
  INTEGER :: is
  !! Counter on polarization
  INTEGER :: js
  !! Counter on polarization
  INTEGER ::  ijs
  !! Counter on combined is and js polarization
  !
  REAL(kind=DP), ALLOCATABLE :: deff(:,:,:)
  !
  COMPLEX(kind=DP), ALLOCATABLE :: ps1(:,:), ps2(:,:,:), aux(:), deff_nc(:,:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: ps1_nc(:,:,:), ps2_nc(:,:,:,:)
  !
  LOGICAL :: ok
  ! 
  CALL start_clock('dvqpsi_us_on')
  IF (noncolin) THEN
    ALLOCATE( ps1_nc(nkb, npol, lower_band:upper_band) )
    ALLOCATE( ps2_nc(nkb, npol, lower_band:upper_band, 3) )
    ALLOCATE( deff_nc(nhm, nhm, nat, nspin) )
  ELSE
    ALLOCATE( ps1(nkb, lower_band:upper_band) )
    ALLOCATE( ps2(nkb, lower_band:upper_band, 3) )
    ALLOCATE( deff(nhm, nhm, nat) )
  ENDIF
  ALLOCATE( aux(npwx) )
  IF (lsda) current_spin = isk(ik)
  !
  !   we first compute the coefficients of the vectors
  !
  IF (noncolin) THEN
    ps1_nc(:,:,:)   = czero
    ps2_nc(:,:,:,:) = czero
  ELSE
    ps1(:,:)   = czero
    ps2(:,:,:) = czero
  ENDIF
  !
  DO ibnd = lower_band, upper_band
    IF (noncolin) THEN
      CALL compute_deff_nc( deff_nc, et(ibnd,ik) )
    ELSE
      CALL compute_deff( deff, et(ibnd,ik) )
    ENDIF
    !
    ijkb0 = 0
    DO nt = 1, ntyp
      DO na = 1, nat
        IF (ityp(na) .eq. nt) THEN
          mu = 3 * (na - 1)
          DO ih = 1, nh(nt)
            ikb = ijkb0 + ih
            DO jh = 1, nh(nt)
              jkb = ijkb0 + jh
              DO ipol = 1, 3
                IF ( abs( uact(mu+1) ) + abs( uact(mu+2) ) + abs( uact(mu+3) ) > eps12 ) THEN
                  IF (noncolin) THEN
                    ijs = 0
                    DO is = 1, npol
                      DO js = 1, npol
                        ijs = ijs + 1
                        ps1_nc(ikb,is,ibnd) = ps1_nc(ikb,is,ibnd) + &
                               deff_nc(ih,jh,na,ijs) *           &
                               alphap(ipol,ik)%nc(jkb,js,ibnd) * uact(mu+ipol)
                        ps2_nc(ikb,is,ibnd,ipol) = ps2_nc(ikb,is,ibnd,ipol) + &
                               deff_nc(ih,jh,na,ijs) * becp1(ik)%nc(jkb,js,ibnd) * &
                               (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                      ENDDO
                    ENDDO
                  ELSE
                    ps1(ikb,ibnd) = ps1(ikb,ibnd) + &
                        deff(ih,jh,na) * &
                        alphap(ipol,ik)%k(jkb,ibnd) * uact(mu+ipol)
                    ps2(ikb,ibnd,ipol) = ps2(ikb,ibnd,ipol) + &
                        deff(ih,jh,na) * becp1(ik)%k(jkb,ibnd) * &
                        (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                  ENDIF
!                  IF (okvan) THEN
!                    IF (noncolin) THEN
!                      ijs = 0
!                      DO is = 1, npol
!                        DO js = 1, npol
!                          ijs = ijs + 1
!                          ps1_nc(ikb,is,ibnd) = ps1_nc(ikb,is,ibnd) + &
!                             int1_nc(ih,jh,ipol,na,ijs) * &
!                             becp1(ik)%nc(jkb,js,ibnd) * uact(mu+ipol)
!                        ENDDO
!                      ENDDO
!                    ELSE
!                      ps1(ikb,ibnd) = ps1(ikb, ibnd) + &
!                          int1(ih,jh,ipol,na,current_spin) * &
!                          becp1(ik)%k(jkb,ibnd) * uact(mu+ipol)
!                    ENDIF
!                  ENDIF ! okvan
                ENDIF  ! uact>0
!                IF (okvan) THEN
!                  DO nb = 1, nat
!                    nu = 3 * (nb - 1)
!                    IF (noncolin) THEN
!                      IF (lspinorb) THEN
!                        ijs = 0
!                        DO is = 1, npol
!                          DO js = 1, npol
!                            ijs = ijs + 1
!                            ps1_nc(ikb,is,ibnd) = ps1_nc(ikb,is,ibnd) + &
!                               int2_so(ih,jh,ipol,nb,na,ijs) * &
!                               becp1(ik)%nc(jkb,js,ibnd) * uact(nu+ipol)
!                          ENDDO
!                        ENDDO
!                      ELSE
!                        DO is = 1, npol
!                          ps1_nc(ikb,is,ibnd) = ps1_nc(ikb,is,ibnd) + &
!                             int2(ih,jh,ipol,nb,na) * &
!                             becp1(ik)%nc(jkb,is,ibnd) * uact(nu+ipol)
!                        ENDDO
!                      ENDIF
!                    ELSE
!                      ps1(ikb,ibnd) = ps1(ikb,ibnd) + &
!                          int2(ih,jh,ipol,nb,na) * &
!                          becp1(ik)%k(jkb,ibnd) * uact(nu+ipol)
!                    ENDIF
!                  ENDDO
!                ENDIF  ! okvan
              ENDDO ! ipol
            ENDDO ! jh
          ENDDO ! ih
          ijkb0 = ijkb0 + nh(nt)
        ENDIF
      ENDDO  ! na
    ENDDO ! nt
  ENDDO ! nbnd
  !
  !      This term is proportional to beta(k+q+G)
  !
  IF (nkb.gt.0) THEN
    IF (noncolin) THEN
      CALL zgemm( 'n', 'n', npwq, (upper_band-lower_band+1)*npol, nkb, &
                  cone, vkb, npwx, ps1_nc, nkb, cone, dvpsi, npwx )
    ELSE
      CALL zgemm( 'n', 'n', npwq, (upper_band-lower_band+1), nkb, &
                  cone, vkb, npwx, ps1, nkb, cone, dvpsi, npwx )
    ENDIF
  ENDIF
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  DO ikb = 1, nkb
    DO ipol = 1, 3
      ok = .false.
      IF (noncolin) THEN
        DO ibnd = lower_band, upper_band
           ok = ok .OR. ( abs( ps2_nc(ikb,1,ibnd,ipol) ) .gt. eps12 ) .OR. &
                        ( abs( ps2_nc(ikb,2,ibnd,ipol) ) .gt. eps12 )
        ENDDO
      ELSE
        DO ibnd = lower_band, upper_band
           ok = ok .OR. ( abs( ps2(ikb,ibnd,ipol) ) .gt. eps12 )
        ENDDO
      ENDIF
      IF (ok) THEN
        DO ig = 1, npwq
          igg = igkq(ig)
          !aux(ig) =  vkb(ig,ikb) * ( xk(ipol,ikq) + g(ipol,igg) )
          aux(ig) = vkb(ig,ikb) * ( xxkq(ipol) + g(ipol,igg) )
        ENDDO
        DO ibnd = lower_band, upper_band
          IF (noncolin) THEN
             CALL zaxpy( npwq, ps2_nc(ikb,1,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1 )
             CALL zaxpy( npwq, ps2_nc(ikb,2,ibnd,ipol), aux, 1, dvpsi(1+npwx,ibnd), 1 )
          ELSE
             CALL zaxpy( npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1 )
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  !
  DEALLOCATE(aux)
  IF (noncolin) THEN
    DEALLOCATE(ps2_nc)
    DEALLOCATE(ps1_nc)
    DEALLOCATE(deff_nc)
  ELSE
    DEALLOCATE(ps2)
    DEALLOCATE(ps1)
    DEALLOCATE(deff)
  ENDIF
  !
  CALL stop_clock('dvqpsi_us_on')
  !
  RETURN
  !
  END SUBROUTINE dvqpsi_us_only3
