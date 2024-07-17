!
! Copyright (C) 2005 Paolo Umari
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE h_epsi_her_apply( lda, n, nbande, psi, hpsi, pdir, e_field )
  !-----------------------------------------------------------------------
  !! This subroutine applies w_k+w_k* on psi, (as in Souza et al.
  !! PRB B 69, 085106 (2004)). The output is put into hpsi.
  !
  !! * evcel must contain the wavefunctions from previous iteration;
  !! * spin polarized systems supported only with fixed occupations.
  !
  USE kinds,                ONLY : DP
  USE noncollin_module,     ONLY : noncolin, npol, lspinorb
  USE wvfct,                ONLY : npwx, nbnd, ik => current_k
  USE lsda_mod,             ONLY : current_spin, nspin
  USE gvect
  USE uspp,                 ONLY : okvan, nkb, vkb, qq_so, qq_at
  USE uspp_param,           ONLY : nh, nhm, nbetam
  USE bp
  USE klist
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE ions_base,            ONLY : ityp, tau, nat, ntyp => nsp
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE fixed_occ
  USE io_global,            ONLY : stdout
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: pdir
  !! direction on which the polarization is calculated
  REAL(DP) :: e_field
  !! electric field along pdir    
  INTEGER :: lda
  !! leading dimension
  INTEGER ::  n
  !! total number of wavefunctions 
  INTEGER :: nbande
  !! number of wavefunctions to be calculated 
  COMPLEX(DP) :: psi(lda*npol, nbande)
  !! the wavefunction
  COMPLEX(DP) :: hpsi(lda*npol,nbande)
  !! output: (w_k+w_k*) psi
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE  :: evct(:,:) !temporary array
  COMPLEX(DP) :: ps(nkb,nbnd*npol)
  !
  TYPE(bec_type) :: becp0
  INTEGER :: nkbtona(nkb)
  INTEGER :: nkbtonh(nkb)
  COMPLEX(DP) :: sca, sca1, pref
  INTEGER :: npw, nb,mb, jkb, nhjkb, na, np, nhjkbm,jkb1,i,j,iv
  INTEGER :: jkb_bp,nt,ig, ijkb0,ibnd,jh,ih,ikb
  REAL(dp) :: eps
  COMPLEX(kind=DP), ALLOCATABLE :: sca_mat(:,:),sca_mat1(:,:)
  COMPLEX(kind=DP) :: pref0(4)
  !
  !
  !  --- Define a small number ---
  eps = 0.000001d0
  IF (ABS(e_field)<eps) RETURN
  CALL start_clock( 'h_epsi_apply' )
  !
  ALLOCATE( evct(npwx*npol,nbnd) )
  CALL allocate_bec_type( nkb, nbnd, becp0 )
  npw = ngk(ik) 
  IF (okvan) THEN
     ! --- Initialize arrays ---
     jkb_bp = 0
      DO nt = 1, ntyp
         DO na = 1, nat
            IF (ityp(na) == nt) THEN
               DO i = 1, nh(nt)
                  jkb_bp = jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      CALL calbec( npw, vkb, psi, becp0, nbande )
  ENDIF
  !
  ALLOCATE( sca_mat(nbnd,nbande), sca_mat1(nbnd,nbande) )
  !
  CALL ZGEMM( 'C', 'N', nbnd, nbande, npw, (1.d0,0.d0), evcel, npwx*npol, &
              psi, npwx*npol, (0.d0,0.d0), sca_mat, nbnd )
  IF (noncolin) THEN
     CALL ZGEMM( 'C', 'N', nbnd, nbande, npw, (1.d0,0.d0), evcel(npwx+1,1), &
                 npwx*npol, psi(npwx+1,1), npwx*npol, (1.d0,0.d0), sca_mat, nbnd )
  ENDIF
  !
  CALL mp_sum( sca_mat, intra_bgrp_comm )
  !
  IF (okvan) THEN
     CALL start_clock( 'h_eps_van2' )
     !
!apply w_k 
        !
        !  
        DO nb = 1, nbande
           DO jkb = 1, nkb
              nhjkb = nkbtonh(jkb)
              na = nkbtona(jkb)
              np = ityp(na)
              nhjkbm = nh(np)
              jkb1 = jkb - nhjkb
              pref0 = (0.d0,0.d0)
              DO j = 1,nhjkbm
                 ! bec_evcel is relative to ik
                 IF (lspinorb) THEN
                    pref0(1) = pref0(1)+becp0%nc(jkb1+j,1,nb) &
                               *qq_so(nhjkb,j,1,np)
                    pref0(2) = pref0(2)+becp0%nc(jkb1+j,2,nb) &
                               *qq_so(nhjkb,j,2,np)
                    pref0(3) = pref0(3)+becp0%nc(jkb1+j,1,nb) &
                               *qq_so(nhjkb,j,3,np)
                    pref0(4) = pref0(4)+becp0%nc(jkb1+j,2,nb) &
                               *qq_so(nhjkb,j,4,np)

                 ELSE
                    pref0(1) = pref0(1)+becp0%k(jkb1+j,nb) &
                               *qq_at(nhjkb,j,na)
                 ENDIF
              ENDDO
              !
              DO mb = 1, nbnd
                 IF (lspinorb) THEN
                    pref=(0.d0,0.d0)
                    pref = pref+CONJG(bec_evcel%nc(jkb,1,mb))*pref0(1)
                    pref = pref+CONJG(bec_evcel%nc(jkb,1,mb))*pref0(2)
                    pref = pref+CONJG(bec_evcel%nc(jkb,2,mb))*pref0(3)
                    pref = pref+CONJG(bec_evcel%nc(jkb,2,mb))*pref0(4)

                 ELSE
                    pref = CONJG(bec_evcel%k(jkb,mb))*pref0(1)
                 ENDIF
                 sca_mat(mb,nb) = sca_mat(mb,nb)+pref
              ENDDO
              !
           ENDDO
        ENDDO
        !
        !
     CALL stop_clock( 'h_eps_van2' )
  ENDIF
  !
  CALL ZGEMM( 'N', 'N', npw, nbande, nbnd, fact_hepsi(ik,pdir), evcelm(1,1,pdir), npwx*npol, &
              sca_mat,nbnd, (1.d0,0.d0), hpsi, npwx*npol )
  CALL ZGEMM( 'N', 'N', npw, nbande, nbnd, -fact_hepsi(ik,pdir), evcelp(1,1,pdir), npwx*npol,&
              sca_mat,nbnd, (1.d0,0.d0), hpsi, npwx*npol )
  IF (noncolin) THEN
     CALL ZGEMM( 'N', 'N', npw, nbande, nbnd, fact_hepsi(ik,pdir), evcelm(1+npwx,1,pdir), npwx*npol, &
                 sca_mat, nbnd, (1.d0,0.d0), hpsi(1+npwx,1), npwx*npol )
     CALL ZGEMM( 'N', 'N', npw, nbande, nbnd, -fact_hepsi(ik,pdir), evcelp(1+npwx,1,pdir), npwx*npol,&
                 sca_mat, nbnd, (1.d0,0.d0), hpsi(1+npwx,1), npwx*npol )
  ENDIF
  !
!apply w_k*
  !
  IF (.NOT.okvan) THEN
     !
     DO nb = 1, nbande
        DO mb = 1, nbnd!index on states of evcel        
           sca = dot_product(evcelm(1:npw,mb,pdir),psi(1:npw,nb))
           IF (noncolin) sca = sca+dot_product(evcelm(1+npwx:npw+npwx,mb,pdir),psi(1+npwx:npw+npwx,nb))
           sca1 = dot_product(evcelp(1:npw,mb,pdir),psi(1:npw,nb))
           IF (noncolin) sca1 = sca1 + dot_product(evcelp(1+npwx:npw+npwx,mb,pdir),psi(1+npwx:npw+npwx,nb))
           CALL mp_sum( sca,  intra_bgrp_comm )
           CALL mp_sum( sca1, intra_bgrp_comm )
           !
           DO ig = 1, npw
              !
              hpsi(ig,nb) = hpsi(ig,nb) + &
                            CONJG(fact_hepsi(ik,pdir))*evcel(ig,mb)*(sca-sca1)
              IF (noncolin) hpsi(ig+npwx,nb) = hpsi(ig+npwx,nb) + &
                            CONJG(fact_hepsi(ik,pdir))*evcel(ig+npwx,mb)*(sca-sca1)
           ENDDO
        ENDDO
     ENDDO
     !
  ELSE ! US case
     !
     CALL start_clock( 'h_eps_ap_van' )
! copy evcel into evct
     DO iv = 1, nbnd
        DO ig = 1, npwx*npol
           evct(ig,iv) = evcel(ig,iv)
        ENDDO
     ENDDO
!  calculate S|evct>
     CALL start_clock( 'h_eps_van2' )
     ps(:,:) = (0.d0, 0.d0)
     ijkb0 = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           if (ityp(na) == nt) THEN
              DO ibnd = 1, nbnd
                 DO jh = 1, nh(nt)
                    jkb = ijkb0 + jh
                    DO ih = 1, nh(nt)
                       ikb = ijkb0 + ih
                       IF (lspinorb) THEN
                          ps(ikb,(ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1) + &
                                                    qq_so(ih,jh,1,nt)* bec_evcel%nc(jkb,1,ibnd)
                          ps(ikb,(ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1) + &
                                                    qq_so(ih,jh,2,nt)* bec_evcel%nc(jkb,2,ibnd)
                          ps(ikb,(ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2) + &
                                                    qq_so(ih,jh,3,nt)* bec_evcel%nc(jkb,1,ibnd)
                          ps(ikb,(ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2) + &
                                                    qq_so(ih,jh,4,nt)* bec_evcel%nc(jkb,2,ibnd)
                       ELSE
                          ps(ikb, ibnd) = ps(ikb,ibnd) + &
                                          qq_at(ih,jh,na)* bec_evcel%k(jkb,ibnd)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (nt)
           ENDIF
        ENDDO
     ENDDO
     !
     CALL stop_clock( 'h_eps_van2' )
     !
     CALL ZGEMM( 'N', 'N', npw, nbnd*npol, nkb, (1.d0, 0.d0), vkb, &!vkb is relative to the last ik read
                 npwx, ps, nkb, (1.d0, 0.d0), evct, npwx )
     !
     CALL ZGEMM( 'C', 'N', nbnd, nbande, npw, (1.d0,0.d0), evcelm(1,1,pdir), npwx*npol, &
                 psi, npwx*npol, (0.d0,0.d0), sca_mat, nbnd )
     IF (noncolin) THEN
        CALL ZGEMM( 'C', 'N', nbnd, nbande, npw, (1.d0,0.d0), evcelm(npwx+1,1,pdir), npwx*npol, &
                    psi(npwx+1,1), npwx*npol, (1.d0,0.d0), sca_mat,nbnd )
     ENDIF
     CALL mp_sum( sca_mat, intra_bgrp_comm )
     CALL ZGEMM( 'C', 'N', nbnd, nbande, npw, (1.d0,0.d0), evcelp(1,1,pdir), npwx*npol, &
                 psi, npwx*npol, (0.d0,0.d0), sca_mat1, nbnd )
     IF (noncolin) THEN
        CALL ZGEMM( 'C', 'N', nbnd, nbande, npw, (1.d0,0.d0), evcelp(npwx+1,1,pdir), npwx*npol, &
                    psi(npwx+1,1), npwx*npol, (1.d0,0.d0), sca_mat1, nbnd )
     ENDIF
     CALL mp_sum( sca_mat1, intra_bgrp_comm )
     !
     sca_mat(1:nbnd,1:nbande) = sca_mat(1:nbnd,1:nbande)-sca_mat1(1:nbnd,1:nbande)
     CALL ZGEMM( 'N', 'N', npw, nbande, nbnd, DCONJG(fact_hepsi(ik,pdir)), evct(1,1), npwx*npol, &
                 sca_mat, nbnd, (1.d0,0.d0), hpsi, npwx*npol )
     IF (noncolin) THEN
        CALL ZGEMM( 'N', 'N', npw, nbande, nbnd, DCONJG(fact_hepsi(ik,pdir)), evct(1+npwx,1), &
                    npwx*npol, sca_mat, nbnd, (1.d0,0.d0), hpsi(1+npwx,1), npwx*npol )
     ENDIF
     !
     CALL stop_clock( 'h_eps_ap_van' )
     !
  ENDIF
  !
  DEALLOCATE( evct )
  !
  CALL deallocate_bec_type( becp0 )
  !
  CALL stop_clock( 'h_epsi_apply' )
  !
  DEALLOCATE( sca_mat  )
  DEALLOCATE( sca_mat1 )
  !
  RETURN
  !
 END SUBROUTINE h_epsi_her_apply
