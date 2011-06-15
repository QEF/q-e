!!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!----------------------------------------------------------------------
!
!
! Modified by Osman Baris Malcioglu (2009)
SUBROUTINE lr_dvpsi_e(ik,ipol,dvpsi)
!----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set) !OBM: This is now handled elesewhere
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba, at
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : xk
  USE gvect,           ONLY : g, gstart
  USE wvfct,           ONLY : npw, npwx, nbnd, igk, g2kin, et
  USE wavefunctions_module, ONLY: evc
  USE lsda_mod,        ONLY : current_spin
  USE spin_orb,        ONLY : lspinorb
  USE noncollin_module,ONLY : noncolin, npol
  USE uspp,            ONLY : okvan, nkb, vkb, qq, qq_so, deeq, deeq_nc
  USE uspp_param,      ONLY : nh
  USE control_flags,   ONLY : gamma_only
  !USE ramanm,          ONLY : eth_rps
  !USE eqv,             ONLY : d0psi, eprec
  !USE phus,            ONLY : becp1, becp1_nc !becp1 is calculated here for every k point
  !USE qpoint,          ONLY : npwq, nksq
  !USE units_ph,        ONLY : this_pcxpsi_is_on_file, lrcom, iucom, &
  !                            lrebar, iuebar !related to removed reading part
  USE control_ph,      ONLY : nbnd_occ

  USE mp_global,       ONLY: intra_pool_comm
  USE mp,              ONLY: mp_sum

  USE realus,                ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma,npw_k, real_space_debug

   USE lr_variables,   ONLY : lr_verbosity, evc0
   USE io_global,      ONLY : stdout
!DEBUG
  USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,check_vector_gamma,check_vector_f
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ipol, ik
  !
  COMPLEX(kind=dp),INTENT(out) :: dvpsi(npwx,nbnd)
  real(kind=dp) :: atnorm

  !

  !
  ! Local variables
  !  !
  COMPLEX(kind=dp),ALLOCATABLE :: d0psi(:,:)
  !

  INTEGER :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0, nrec
  ! counters

  real(DP), ALLOCATABLE  :: gk (:,:), h_diag (:,:), eprec(:)

  ! the derivative of |k+G|
  real(DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter

  LOGICAL :: conv_root
  ! true if convergence has been achieved

  COMPLEX(DP), ALLOCATABLE :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:),  &
       work (:,:), spsi(:,:), psc(:,:,:,:)
  real(kind=dp), EXTERNAL :: ddot
  COMPLEX(DP), EXTERNAL :: ZDOTC
  ! the scalar products
  EXTERNAL ch_psi_all, cg_psi
  !
  !obm debug
  real(DP) ::obm_debug
  !
  CALL start_clock ('lr_dvpsi_e')

  IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dvpsi_e>")')

  ALLOCATE(d0psi(npwx*npol,nbnd))
  d0psi=(0.d0, 0.d0)
  dvpsi=(0.d0, 0.d0)
  !if (this_pcxpsi_is_on_file(ik,ipol)) then
  !   nrec = (ipol - 1)*nksq + ik
  !   call davcio(dvpsi, lrebar, iuebar, nrec, -1)
  !   call stop_clock ('dvpsi_e')
  !   return
  !end if
  !
  ALLOCATE (work ( npwx, max(nkb,1)))

  ALLOCATE (gk ( 3, npwx))
  ALLOCATE (h_diag( npwx*npol, nbnd))
  !OBM!!!! eprec is also calculated on the fly for each k point
  ALLOCATE(eprec(nbnd))
  !OBM!!!!
  evc(:,:)=evc0(:,:,ik)
  IF (nkb > 0) THEN
     ALLOCATE (dvkb (npwx, nkb), dvkb1(npwx, nkb))
     dvkb (:,:) = (0.d0, 0.d0)
     dvkb1(:,:) = (0.d0, 0.d0)
  ENDIF
  DO ig = 1, npw_k(ik)
     gk (1, ig) = (xk (1, ik) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, ik) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, ik) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
  ENDDO

  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  DO ibnd = 1, nbnd_occ (ik)
     DO ig = 1, npw_k(ik)
        d0psi (ig, ibnd) = (at(1, ipol) * gk(1, ig) + &
             at(2, ipol) * gk(2, ig) + &
             at(3, ipol) * gk(3, ig) ) &
             *(0.d0,-2.d0)*evc (ig, ibnd)
     ENDDO
     IF (noncolin) THEN
        DO ig = 1, npw_k(ik)
           d0psi (ig+npwx, ibnd) = (at(1, ipol) * gk(1, ig) + &
                at(2, ipol) * gk(2, ig) + &
                at(3, ipol) * gk(3, ig) ) &
                 *(0.d0,-2.d0)*evc (ig+npwx, ibnd)
        ENDDO
     ENDIF
  ENDDO



!
! Uncomment this goto and the continue below to calculate
! the matrix elements of p without the commutator with the
! nonlocal potential.
!
!  goto 111
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  CALL gen_us_dj (ik, dvkb)
  CALL gen_us_dy (ik, at (1, ipol), dvkb1)


  DO ig = 1, npw_k(ik)
     IF (g2kin (ig) < 1.0d-10) THEN
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     ELSE
        gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
     ENDIF
  ENDDO

  jkb = 0
  work=(0.d0,0.d0)
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (nt == ityp (na)) THEN
           DO ikb = 1, nh (nt)
              jkb = jkb + 1
              DO ig = 1, npw_k(ik)
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (at (1, ipol) * gk (1, ig) + &
                       at (2, ipol) * gk (2, ig) + &
                       at (3, ipol) * gk (3, ig) )
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE (gk)
  !OBM!!!be careful, from bwalker, why?!!!!!
  work(:,:)=(0.0d0,1.0d0)*work(:,:)
  !OBM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IF(gamma_only) THEN
     CALL lr_dvpsi_e_gamma()
  ELSEIF (noncolin) THEN
     CALL lr_dvpsi_e_noncolin()
  ELSE
     CALL lr_dvpsi_e_k()
  ENDIF

  IF (nkb > 0) THEN
     DEALLOCATE (dvkb1, dvkb)
  ENDIF
  !OBM!!!! End of Brent's seperation

  DEALLOCATE (h_diag)
  DEALLOCATE (work)

  DEALLOCATE (eprec)
  DEALLOCATE (d0psi)

  !OBM!!!! Addendum to PH dvpsi
  IF (okvan) THEN
    ALLOCATE (spsi ( npwx*npol, nbnd))
    CALL sm1_psi(.true.,ik,npwx,npw_k(ik),nbnd,dvpsi,spsi)
    dvpsi(:,:) = spsi(:,:)
    DEALLOCATE(spsi)
  ENDIF

  !
  ! For some ibrav the crystal axes are not normalized
  ! Here we include the correct normalization
  ! for Lanczos initial wfcs
  atnorm=dsqrt(at(1,ipol)**2+at(2,ipol)**2+at(3,ipol)**2)
  !
  dvpsi(:,:)=dvpsi(:,:)/atnorm
  !
  !nrec = (ipol - 1)*nksq + ik
  !call davcio(dvpsi, lrebar, iuebar, nrec, 1)
  !this_pcxpsi_is_on_file(ik,ipol) = .true.
  CALL stop_clock ('lr_dvpsi_e')
  !
  RETURN
  !
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!k-point modifications
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE lr_dvpsi_e_k()
  !OBM!! The k point part is essentially the same as the phonon dvpsi_e, with noncolin part removed
  USE becmod,          ONLY : bec_type, becp, calbec

  IMPLICIT NONE

  COMPLEX(kind=dp), ALLOCATABLE :: becp1(:,:),becp2(:,:)

  IF (nkb > 0) THEN
        ALLOCATE (becp1 (nkb, nbnd))
        ALLOCATE (becp2 (nkb, nbnd))
  ENDIF
     becp1=0.0
     becp2=0.0
     becp%k=0.0

     ! OBM : Change calbec to use calbec_type as soon as possible
     !
     CALL calbec (npw_k(ik), vkb, evc, becp1)

     IF(nkb>0) THEN !If there are no beta functions, below loop returns an error
       CALL calbec (npw_k(ik), work, evc, becp2)

     ijkb0 = 0
     ALLOCATE (ps2 ( nkb, nbnd, 2))
     ps2=(0.d0,0.d0)
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (nt == ityp (na)) THEN
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd = 1, nbnd_occ (ik)
                         !OBM!!! Notice that in the original version ikb is 1, I
                         !       think this was a bug of the original code
                          ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                              (deeq(ih,jh,na,current_spin) &
                              -et(ibnd,ik)*qq(ih,jh,nt))
                          ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd) * &
                              (-1.d0,0.d0)*(deeq(ih,jh,na,current_spin)&
                              -et(ibnd,ik)*qq(ih,jh,nt))
                         !original from Phonon
                         !ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         ! ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd,ik) * &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         !OBM!!! But why?
                    ENDDO
                 ENDDO
              ENDDO
              ijkb0=ijkb0+nh(nt)
           ENDIF
        ENDDO
     ENDDO
     IF (ikb /= nkb .or. jkb /= nkb) CALL errore ('lr_dvpsi_e', 'unexpected error',1)
     !OBM!!! In the archaic version Brent used, this part was embedded in the
     !       above loop as ZAXPY
     !

        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik), nkb, &
             (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik), nkb, &
             (1.d0,0.d0),work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        DEALLOCATE (ps2)
     ENDIF
     !print *, "new version k point, before ortho", d0psi(1:3,1)
     !
     !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
     !    Apply -P^+_c
     !OBM!! lr_ortho no longer calculates sevc, do it beforhand

     IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd)
     CALL s_psi (npwx, npw_k(ik), nbnd, evc, dvpsi)
     CALL lr_ortho(d0psi, evc, ik, ik, dvpsi,.false.)
     !d0psi=-d0psi

    !print *, "new version k point, after ortho", d0psi(1:3,1)


     !
     !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
     !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
     !
     !OBM!!! thresh = eth_rps
     thresh = 1.d-5
     h_diag=0.d0

     !OBM!!! on the fly calculation of eprec
     ! I am using the work as a working memory, it seems not to be used anymore
        DEALLOCATE (work)
        ALLOCATE (work(npwx,nbnd))
        DO ibnd=1,nbnd
          !
          work = 0.d0
          !
          conv_root = .true.
          !
          DO ig=1,npw_k(ik)
             work(ig,1)=g2kin(ig)*evc(ig,ibnd)
          ENDDO
          !
          eprec(ibnd)=1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1)
          !
       ENDDO
       !
#ifdef __PARA
       CALL mp_sum(eprec, intra_pool_comm)
#endif
       !OBM!!!
       !print *, "eprec", eprec
       DO ibnd = 1, nbnd_occ (ik)
          DO ig = 1, npw_k(ik)
             h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
          ENDDO
          IF (noncolin) THEN
             DO ig = 1, npw_k(ik)
                h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
             ENDDO
          ENDIF
       ENDDO
       !
       !print *, h_diag(1:10,1)
       !OBM!!! upto here, the dvpsi was used as a scratch
       !
       dvpsi(:,:) = (0.d0, 0.d0)
       !
       !OBM!!!! Original was
       !call lr_cgsolve_all(et(1,ik),d0psi,dvpsi,h_diag,npwx,npw_k(ik),&
       !            thresh,ik,lter,conv_root,anorm,nbnd_occ (ik),lr_alpha_pv)
       CALL cgsolve_all (ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
            h_diag, npwx, npw_k(ik), thresh, ik, lter, conv_root, anorm, &
            nbnd_occ(ik), 1)
       !
       IF (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
            & " linter: root not converged ",e10.3)') &
            ik, ibnd, anorm
       !
       CALL flush_unit( stdout )
       !print *, "after cg_solve_all"
       !CALL lr_normalise( dvpsi(:,:), anorm)
       !
       !
       ! we have now obtained P_c x |psi>.
       ! In the case of USPP this quantity is needed for the Born
       ! effective charges, so we save it to disc
       !
       ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
       ! therefore we apply S again, and then subtract the additional term
       ! furthermore we add the term due to dipole of the augmentation charges.
       !
       IF (okvan) THEN
          !
          ! for effective charges
          !
          ! nrec = (ipol - 1) * nksq + ik
          ! call davcio (dvpsi, lrcom, iucom, nrec, 1)
          !
          ALLOCATE (spsi ( npwx*npol, nbnd))
          CALL calbec (npw_k(ik), vkb, dvpsi, becp )
          CALL s_psi(npwx,npw_k(ik),nbnd,dvpsi,spsi)
          CALL DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
          DEALLOCATE (spsi)
          CALL lr_adddvepsi_us_k(becp1,becp2,ipol,ik,dvpsi)
          !call adddvepsi_us(becp2,ipol,ik)
       ENDIF

       !OBM!!!! Addendum to PH/dvpsi_e
       !
       work = 0.d0 !Reset working space in any case
       ! orthogonalize dvpsi to the valence subspace
       !
       !OBM!!!! due to modifications, lr_ortho requires sevc as input, putting it into work
       IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd)
       CALL s_psi (npwx, npw_k(ik), nbnd, evc, work)
       CALL lr_ortho(dvpsi, evc, ik, ik, work,.false.) !very strange, check: wothout the last false, normconserving pps fail
       !dvpsi=-dvpsi
       !OBM!!!! end of orthogonolization
       IF (nkb > 0) THEN
             DEALLOCATE(becp1)
             DEALLOCATE(becp2)
       ENDIF
       !
       !OBM!!!!
  END SUBROUTINE lr_dvpsi_e_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!noncolin modifications
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE lr_dvpsi_e_noncolin()
  !OBM!! essentially the same as the phonon dvpsi_e, only with non-collinear part
  USE becmod,          ONLY : bec_type, becp, calbec

  IMPLICIT NONE

  COMPLEX(kind=dp), ALLOCATABLE :: becp1_nc(:,:,:), becp2_nc(:,:,:)

  CALL errore ('lr_dvpsi_e', 'non collinear not implemented', 1) !well, everything is here, but have to understand the formulation, do not have the time
  IF (nkb > 0) THEN
        ALLOCATE (becp1_nc (nkb, npol, nbnd))
        ALLOCATE (becp2_nc (nkb, npol, nbnd))
  ENDIF

     !OBM!!! This part is seperated as gamma_only and k_point in Brent's version
        !OBM!!becp1 is also calculated
        CALL calbec (npw_k(ik), vkb, evc, becp1_nc)
        IF(nkb>0) CALL calbec (npw_k(ik), work, evc, becp2_nc)

     ijkb0 = 0
        ALLOCATE (psc ( nkb, npol, nbnd, 2))
        psc=(0.d0,0.d0)
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (nt == ityp (na)) THEN
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd = 1, nbnd_occ (ik)
                          IF (lspinorb) THEN
                             psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+(0.d0,-1.d0)* &
                                (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+       &
                                 becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)- &
                                          et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                             psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+(0.d0,-1.d0)*  &
                                (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+       &
                                 becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)- &
                                          et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                             psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+(0.d0,-1.d0)* &
                                (becp1_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+      &
                                becp1_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)-  &
                                          et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                             psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+(0.d0,-1.d0)*  &
                                (becp1_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+      &
                                becp1_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)-  &
                                          et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                          ELSE
                             psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+ (0.d0,-1.d0)* &
                                 ( becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                                -et(ibnd,ik)*qq(ih,jh,nt)) + &
                                   becp2_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                             psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+ (0.d0,-1.d0)* &
                                 ( becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                                -et(ibnd,ik)*qq(ih,jh,nt))+  &
                                   becp2_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                             psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+ (0.d0,-1.d0)* &
                                 ( becp1_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                                -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                   becp1_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                             psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+ (0.d0,-1.d0)* &
                                 ( becp1_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                                -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                   becp1_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                          ENDIF
                    ENDDO
                 ENDDO
              ENDDO
              ijkb0=ijkb0+nh(nt)
           ENDIF
        ENDDO
     ENDDO
     IF (ikb /= nkb .or. jkb /= nkb) CALL errore ('lr_dvpsi_e', 'unexpected error',1)
     !OBM!!! In the archaic version Brent used, this part was embedded in the
     !       above loop as ZAXPY
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik)*npol, nkb, &
             (1.d0,0.d0), vkb(1,1), npwx, psc(1,1,1,1), nkb, (1.d0,0.d0), &
             d0psi, npwx )
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik)*npol, nkb, &
             (1.d0,0.d0),work(1,1), npwx, psc(1,1,1,2), nkb, (1.d0,0.d0), &
             d0psi, npwx )
        DEALLOCATE (psc)

     !
     !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
     !    Apply -P^+_c
     !OBM!! lr_ortho no longer calculates sevc, do it beforhand
     IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd )
     CALL s_psi (npwx, npw_k(ik), nbnd, evc, dvpsi)
     CALL lr_ortho(d0psi, evc, ik, ik, dvpsi, .false.)
     !d0psi=-d0psi
     !
     !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
     !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
     !
     !OBM!!! thresh = eth_rps
     thresh = 1.d-5
     h_diag=0.d0

     !OBM!!! on the fly calculation of eprec
     ! I am using the work as a working memory, it seems not to be used anymore
        DEALLOCATE (work)
        ALLOCATE (work(npwx,nbnd))
        DO ibnd=1,nbnd
          !
          work = 0.d0
          !
          conv_root = .true.
          !
          DO ig=1,npw_k(ik)
             work(ig,1)=g2kin(ig)*evc(ig,ibnd)
          ENDDO
          !
          eprec(ibnd)=1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1)
          !
       ENDDO
       !
#ifdef __PARA
       CALL mp_sum(eprec, intra_pool_comm)
#endif
       !OBM!!!
       !print *, "eprec", eprec
       DO ibnd = 1, nbnd_occ (ik)
          DO ig = 1, npw_k(ik)
             h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
          ENDDO
             DO ig = 1, npw_k(ik)
                h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
             ENDDO
       ENDDO
       !
       !OBM!!! upto here, the dvpsi was used as a scratch
       !
       dvpsi(:,:) = (0.d0, 0.d0)
       !
       !OBM!!!! Original was
       !call lr_cgsolve_all(et(1,ik),d0psi,dvpsi,h_diag,npwx,npw_k(ik),&
       !            thresh,ik,lter,conv_root,anorm,nbnd_occ (ik),lr_alpha_pv)
       CALL cgsolve_all (ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
            h_diag, npwx, npw_k(ik), thresh, ik, lter, conv_root, anorm, &
            nbnd_occ(ik), 1)
       !
       IF (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
            & " linter: root not converged ",e10.3)') &
            ik, ibnd, anorm
       !
       CALL flush_unit( stdout )
       !
       !
       ! we have now obtained P_c x |psi>.
       ! In the case of USPP this quantity is needed for the Born
       ! effective charges, so we save it to disc
       !
       ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
       ! therefore we apply S again, and then subtract the additional term
       ! furthermore we add the term due to dipole of the augmentation charges.
       !
       IF (okvan) THEN
          !
          ! for effective charges
          !
          ! nrec = (ipol - 1) * nksq + ik
          ! call davcio (dvpsi, lrcom, iucom, nrec, 1)
          !
          ALLOCATE (spsi ( npwx*npol, nbnd))
          CALL calbec (npw_k(ik), vkb, dvpsi, becp )
          CALL s_psi(npwx,npw_k(ik),nbnd,dvpsi,spsi)
          CALL DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
          DEALLOCATE (spsi)
          !OBM!!! non collinear in lr_addd_vespi is still not working
          !call adddvepsi_us(becp2_nc,ipol,ik)
       ENDIF
       !OBM!!!! Addendum to PH/dv_psi_e
       !
       work = 0.d0 !Reset working space in any case
       ! orthogonalize dvpsi to the valence subspace
       !
       !OBM!!!! due to modifications, lr_ortho requires sevc as input, putting it into work
       IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd )
       CALL s_psi (npwx, npw_k(ik), nbnd, evc, work)
       CALL lr_ortho(dvpsi, evc, ik, ik, work,.false.)
       !dvpsi=-dvpsi
       !OBM!!!! end of orthogonolization
       IF (nkb > 0) THEN
             DEALLOCATE(becp1_nc)
             DEALLOCATE(becp2_nc)
       ENDIF
       !
       !OBM!!!!
  END SUBROUTINE lr_dvpsi_e_noncolin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!GAMMA point modifications
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lr_dvpsi_e_gamma
  !OBM!! modified for gamma point algorithms
  USE becmod,          ONLY : bec_type,becp, calbec

  IMPLICIT NONE

  real(kind=dp), ALLOCATABLE :: becp1(:,:),becp2(:,:)

  IF (nkb > 0) THEN
        ALLOCATE (becp1 (nkb, nbnd))
        ALLOCATE (becp2 (nkb, nbnd))
  ENDIF
     becp1=0.0
     becp2=0.0
     becp%r=0.0

     CALL calbec (npw, vkb, evc, becp1)

     IF(nkb>0) THEN

      CALL calbec (npw, work, evc, becp2)

     ijkb0 = 0
     ALLOCATE (ps2 ( nkb, nbnd, 2))
     ps2=(0.d0,0.d0)
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (nt == ityp (na)) THEN
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd = 1, nbnd_occ (ik)
                         !OBM!!! Notice that in the original version ikb is 1
                          ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                              (deeq(ih,jh,na,current_spin) &
                              -et(ibnd,ik)*qq(ih,jh,nt))
                          ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd) * &
                              (-1.d0,0.d0)*(deeq(ih,jh,na,current_spin)&
                              -et(ibnd,ik)*qq(ih,jh,nt))
                         !original from Phonon
                         !ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         ! ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd,ik) * &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         !OBM!!! But why?
                    ENDDO
                 ENDDO
              ENDDO
              ijkb0=ijkb0+nh(nt)
           ENDIF
        ENDDO
     ENDDO
     IF (ikb /= nkb .or. jkb /= nkb) CALL errore ('lr_dvpsi_e', 'unexpected error',1)
     !OBM!!! In the archaic version Brent used, this part was embedded in the
     !       above loop as ZAXPY
        CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
             (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
             (1.d0,0.d0),work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        DEALLOCATE (ps2)
     ENDIF
     !
     !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
     !    Apply -P^+_c
     !OBM!! lr_ortho no longer calculates sevc, do it beforhand
     !print *, "new version,gamma before ortho", d0psi(1:3,1)
     IF (okvan) CALL calbec ( npw, vkb, evc, becp, nbnd)
     CALL s_psi (npwx, npw, nbnd, evc, dvpsi)
     CALL lr_ortho(d0psi, evc, ik, ik, dvpsi,.false.) !Strange bug, without last .false., norm conserving pps fail
     !d0psi=-d0psi
     !print *, "new version,gamma after ortho", d0psi(1:3,1)

     !
     !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
     !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
     !
     !OBM!!! thresh = eth_rps
     thresh = 1.d-5
     h_diag=0.d0

     !OBM!!! on the fly calculation of eprec
     ! I am using the work as a working memory, it seems not to be used anymore
        DEALLOCATE (work)
        ALLOCATE (work(npwx,nbnd))
        DO ibnd=1,nbnd
          !
          work = 0.d0
          !
          conv_root = .true.
          !
          DO ig=1,npw
             work(ig,1)=g2kin(ig)*evc(ig,ibnd)
          ENDDO
          !
          eprec(ibnd)=2.0d0*DDOT(2*npw,evc(1,ibnd),1,work,1)
          !
          IF(gstart==2) eprec(ibnd)=eprec(ibnd)-dble(evc(1,ibnd))*dble(work(1,ibnd))
          !
          eprec(ibnd)=1.35d0*eprec(ibnd)
          !OBM!!! was : eprec(ibnd)=1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1)
          !
       ENDDO
       !
#ifdef __PARA
       CALL mp_sum(eprec, intra_pool_comm)
#endif
       !print *, eprec
       !OBM!!!
       !print *, "eprec", eprec
       DO ibnd = 1, nbnd_occ (ik)
          DO ig = 1, npw
             h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
          ENDDO
          IF (noncolin) THEN
             DO ig = 1, npw
                h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
             ENDDO
          ENDIF
       ENDDO
       !print *, h_diag(1:10,1)
       !
       !OBM!!! upto here, the dvpsi was used as a scratch
       !
       dvpsi(:,:) = (0.d0, 0.d0)
       !
       !OBM!!!! Original was
       !call lr_cgsolve_all(et(1,ik),d0psi,dvpsi,h_diag,npwx,npw,&
       !            thresh,ik,lter,conv_root,anorm,nbnd_occ (ik),lr_alpha_pv)
       CALL cgsolve_all (ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
            h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
            nbnd_occ(ik), 1)
       !
       IF (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
            & " linter: root not converged ",e10.3)') &
            ik, ibnd, anorm
       !
       CALL flush_unit( stdout )
       !print *, "after cg_solve_all"
       !CALL lr_normalise( dvpsi(:,:), anorm)
       !
       !
       ! we have now obtained P_c x |psi>.
       ! In the case of USPP this quantity is needed for the Born
       ! effective charges, so we save it to disc
       !
       ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
       ! therefore we apply S again, and then subtract the additional term
       ! furthermore we add the term due to dipole of the augmentation charges.
       !
       IF (okvan) THEN
          !
          ! for effective charges
          !
          ! nrec = (ipol - 1) * nksq + ik
          ! call davcio (dvpsi, lrcom, iucom, nrec, 1)
          !
          ALLOCATE (spsi ( npwx*npol, nbnd))
          CALL calbec (npw, vkb, dvpsi, becp )
          CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
          CALL DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
          DEALLOCATE (spsi)
          CALL lr_adddvepsi_us_gamma(becp1,becp2,ipol,ik,dvpsi)
          !call adddvepsi_us(becp2,ipol,ik)
       ENDIF

       !OBM!!!! Addendum to PH dvpsi_e
       !
       work = 0.d0 !Reset working space in any case
       ! orthogonalize dvpsi to the valence subspace
       !
       !OBM!!!! due to modifications, lr_ortho requires sevc as input, putting it into work
       IF (okvan) CALL calbec ( npw, vkb, evc, becp, nbnd)
       CALL s_psi (npwx, npw_k(ik), nbnd, evc, work)
       CALL lr_ortho(dvpsi, evc, ik, ik, work,.false.)
       !dvpsi=-dvpsi
       !OBM!!!! end of orthogonolization

       IF (nkb > 0) THEN
             DEALLOCATE(becp1)
             DEALLOCATE(becp2)
       ENDIF
       !
       !OBM!!!!
  END SUBROUTINE lr_dvpsi_e_gamma

END SUBROUTINE lr_dvpsi_e


