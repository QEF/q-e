!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE vofrho_x( nfi, rhor, drhor, rhog, drhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
!-----------------------------------------------------------------------
!     computes: the one-particle potential v in real space,
!               the total energy etot,
!               the forces fion acting on the ions,
!               the derivative of total energy to cell parameters h
!     rhor input : electronic charge on dense real space grid
!                  (plus core charge if present)
!     rhog input : electronic charge in g space (up to density cutoff)
!     rhos input : electronic charge on smooth real space grid
!     rhor output: total potential on dense real space grid
!     rhos output: total potential on smooth real space grid
!
      USE kinds,            ONLY: dp
      USE control_flags,    ONLY: iprint, iverbosity, thdyn, tpre, tfor, &
                                  tprnfor, iesr, textfor
      USE io_global,        ONLY: stdout
      USE ions_base,        ONLY: nsp, na, nat, rcmax, compute_eextfor
      USE ions_base,        ONLY: ind_srt, ind_bck
      USE cell_base,        ONLY: omega, r_to_s
      USE cell_base,        ONLY: alat, at, tpiba2, h, ainv
      USE cell_base,        ONLY: ibrav, isotropic  !True if volume option is chosen for cell_dofree
      USE cp_main_variables, ONLY: iprint_stdout    !print control
      USE gvect,            ONLY: gstart, gg, g
      USE electrons_base,   ONLY: nspin
      USE constants,        ONLY: pi, fpi, au_gpa, e2
      USE energies,         ONLY: etot, eself, enl, ekin, epseu, esr, eht, &
                                  exc, eextfor, exx 
      USE local_pseudo,     ONLY: vps, dvps, rhops
      USE uspp,             ONLY: nlcc_any
      USE smallbox_gvec
      USE dener,            ONLY: detot, dekin, dps, dh, dsr, dxc, denl, denlc,&
                                  detot6, dekin6, dps6, dh6, dsr6, dxc6, denl6
      USE mp,               ONLY: mp_sum
      USE mp_global,        ONLY: intra_bgrp_comm
      USE funct,            ONLY: dft_is_meta, dft_is_nonlocc, nlc, get_inlc,&
                                  dft_is_hybrid, exx_is_active
      USE vdW_DF,           ONLY: stress_vdW_DF
      use rVV10,            ONLY: stress_rVV10 
      USE pres_ai_mod,      ONLY: abivol, abisur, v_vol, P_ext, volclu,  &
                                  Surf_t, surfclu
      USE fft_interfaces,   ONLY: fwfft, invfft
      USE sic_module,       ONLY: self_interaction, sic_alpha
      USE energies,         ONLY: self_exc, self_ehte
      USE cp_interfaces,    ONLY: pseudo_stress, compute_gagb, stress_hartree, &
                                  add_drhoph, stress_local, force_loc, self_vofhar
      USE fft_base,         ONLY: dfftp, dffts
      USE ldaU_cp,          ONLY: e_hubbard
      USE control_flags,    ONLY: ts_vdw
      USE tsvdw_module,     ONLY: tsvdw_calculate
      USE tsvdw_module,     ONLY: EtsvdW,UtsvdW,FtsvdW,HtsvdW
      USE mp_global,        ONLY: me_image
      USE exx_module,       ONLY: dexx_dh, exxalfa  ! exx_wf related
      USE fft_rho
      USE fft_helper_subroutines
      
      USE plugin_variables, ONLY: plugin_etot

      IMPLICIT NONE
!
      LOGICAL     :: tlast, tfirst
      INTEGER     :: nfi
      REAL(DP)    :: rhor(:,:), drhor(:,:,:,:), rhos(:,:), fion(:,:)
      REAL(DP)    :: rhoc(:), tau0(:,:)
      ! COMPLEX(DP) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat), ei3(-nr3:nr3,nat)
      COMPLEX(DP) :: ei1(:,:), ei2(:,:), ei3(:,:)
      COMPLEX(DP) :: eigrb(:,:)
      COMPLEX(DP) :: rhog(:,:), drhog(:,:,:,:)
      COMPLEX(DP) :: sfac(:,:)
      INTEGER     :: irb(:,:)
      !
      INTEGER iss, isup, isdw, ig, ir, i, j, k, ij, is, ia, inlc
      REAL(DP) :: vtxc, vave, ebac, wz, eh, ehpre, enlc
      COMPLEX(DP)  fp, fm, ci, drhop, zpseu, zh
      COMPLEX(DP), ALLOCATABLE :: rhotmp(:), vtemp(:)
      COMPLEX(DP), ALLOCATABLE :: drhot(:,:)
      REAL(DP), ALLOCATABLE    :: gagb(:,:), rhosave(:,:), rhocsave(:)
      !
      REAL(DP), ALLOCATABLE :: fion1( :, : )
      REAL(DP), ALLOCATABLE :: stmp( :, : )
      !
      COMPLEX(DP), ALLOCATABLE :: self_vloc(:)
      COMPLEX(DP)              :: self_rhoeg
      REAL(DP)                 :: self_ehtet, fpibg
      LOGICAL                  :: ttsic
      REAL(DP)                 :: detmp( 3, 3 ), desr( 6 ), deps( 6 )
      REAL(DP)                 :: detmp2( 3, 3 )
      REAL(DP)                 :: ht( 3, 3 )
      REAL(DP)                 :: deht( 6 )
      COMPLEX(DP)              :: screen_coul( 1 )
      REAL(DP)                 :: dexx(3,3) ! stress tensor from exact exchange exx_wf related
!
      INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)

      ! ...  dalbe(:) = delta( alpha(:), beta(:) )
      REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)


      CALL start_clock( 'vofrho' )

      !
      !     TS-vdW calculation (RAD)
      !
      IF (ts_vdw) THEN
        !
        CALL start_clock( 'ts_vdw' )
        ALLOCATE (stmp(3,nat), rhocsave(dfftp%nnr) )
        stmp(:,:) = tau0(:,ind_bck(:))
        !
        IF ( nspin==2 ) THEN
           rhocsave(:) = rhor(:,1) + rhor(:,2) 
        ELSE
           rhocsave(:) = rhor(:,1)
        ENDIF
        !
        CALL tsvdw_calculate(stmp,rhocsave)
        !
        DEALLOCATE (rhocsave,stmp)
        CALL stop_clock( 'ts_vdw' )
        !
      END IF
      !
      ci = ( 0.0d0, 1.0d0 )
      !
      !     wz = factor for g.neq.0 because of c*(g)=c(-g)
      !
      wz = 2.0d0
      !
      ht = TRANSPOSE( h )
      !
      ALLOCATE( vtemp( dfftp%ngm ) )
      ALLOCATE( rhotmp( dfftp%ngm ) )
      !
      IF ( tpre ) THEN
         ALLOCATE( drhot( dfftp%ngm, 6 ) )
         ALLOCATE( gagb( 6, dfftp%ngm ) )
         CALL compute_gagb( gagb, g, dfftp%ngm, tpiba2 )
      END IF
!
!     ab-initio pressure and surface tension contributions to the potential
!
      if (abivol.or.abisur) call vol_clu(rhor,rhog,nfi)
      !
      !     compute plugin contributions to the potential, add it later
      !
      CALL plugin_get_potential(rhor,nfi)
      !
      !     compute plugin contributions to the energy
      !
      plugin_etot = 0.0_dp
      !
      CALL plugin_energy(rhor,plugin_etot)
      !
      ttsic = ( ABS( self_interaction ) /= 0 )
      !
      IF( ttsic ) ALLOCATE( self_vloc( dfftp%ngm ) )
      !
      !     first routine in which fion is calculated: annihilation
      !
      fion  = 0.d0
      !
      !     forces on ions, ionic term in real space
      !
      IF( tprnfor .OR. tfor .OR. tfirst .OR. tpre ) THEN
         !
         ALLOCATE( stmp( 3, nat ) )
         !
         CALL r_to_s( tau0, stmp, na, nsp, ainv )
         !
         CALL vofesr( iesr, esr, dsr6, fion, stmp, tpre, h )
         !
         call mp_sum( fion, intra_bgrp_comm )
         !
         DEALLOCATE( stmp )
         !
      END IF
!
!$omp parallel default(shared), private(ig,is,ij,i,j,k)
!$omp workshare
      rhotmp( 1:dfftp%ngm ) = rhog( 1:dfftp%ngm, 1 )
!$omp end workshare
      IF( nspin == 2 ) THEN
!$omp workshare
         rhotmp( 1:dfftp%ngm ) = rhotmp( 1:dfftp%ngm ) + rhog( 1:dfftp%ngm, 2 )
!$omp end workshare
      END IF
      !
      IF( tpre ) THEN
!$omp do
         DO ij = 1, 6
            i = alpha( ij )
            j = beta( ij )
            drhot( :, ij ) = 0.0d0
            DO k = 1, 3
               drhot( :, ij ) = drhot( :, ij ) +  drhog( :, 1, i, k ) * ht( k, j )
            END DO
         END DO
!$omp end do
         IF( nspin == 2 ) THEN
!$omp do
            DO ij = 1, 6
               i = alpha( ij )
               j = beta( ij )
               DO k = 1, 3
                  drhot( :, ij ) = drhot( :, ij ) +  drhog( :, 2, i, k ) * ht( k, j )
               END DO
            END DO
!$omp end do
         ENDIF
      END IF
      !
      !     calculation local potential energy
      !
!$omp master
      zpseu = 0.0d0
!$omp end master 
      !
!$omp do
      DO ig = 1, SIZE(vtemp)
         vtemp(ig)=(0.d0,0.d0)
      END DO
      DO is=1,nsp
!$omp do
         DO ig=1,dffts%ngm
            vtemp(ig)=vtemp(ig)+CONJG(rhotmp(ig))*sfac(ig,is)*vps(ig,is)
         END DO
      END DO
!$omp do reduction(+:zpseu)
      DO ig=1,dffts%ngm
         zpseu = zpseu + vtemp(ig)
      END DO
!$omp end parallel

      epseu = wz * DBLE(zpseu)
      !
      IF (gstart == 2) epseu = epseu - DBLE( vtemp(1) )
      !
      CALL mp_sum( epseu, intra_bgrp_comm )
      epseu = epseu * omega
      !
      IF( tpre ) THEN
         !
         CALL stress_local( dps6, epseu, gagb, sfac, rhotmp, drhot, omega )
         !
      END IF
      !
      !     
      !     calculation hartree energy
      !    
      !
      self_ehtet = 0.d0  
      !
      IF( ttsic ) self_vloc = 0.d0 

      zh = 0.0d0

!$omp parallel default(shared), private(ig,is)

      DO is=1,nsp
!$omp do
         DO ig=1,dffts%ngm
            rhotmp(ig)=rhotmp(ig)+sfac(ig,is)*rhops(ig,is)
         END DO
      END DO
      !
!$omp do
      DO ig = gstart, dfftp%ngm
         vtemp(ig) = CONJG( rhotmp( ig ) ) * rhotmp( ig ) / gg( ig )
      END DO

!$omp do reduction(+:zh)
      DO ig = gstart, dfftp%ngm
         zh = zh + vtemp(ig)
      END DO

!$omp end parallel

      eh = DBLE( zh ) * wz * 0.5d0 * fpi / tpiba2
!
      CALL mp_sum( eh, intra_bgrp_comm )
      !
      IF ( ttsic ) THEN
         !
         CALL self_vofhar( .false., self_ehte, self_vloc, rhog, omega, h )
         !
         eh = eh - self_ehte / omega
         !
      END IF
      !
      IF(tpre) THEN
         !
         CALL add_drhoph( drhot, sfac, gagb )
         CALL stress_hartree(dh6, eh*omega, sfac, rhotmp, drhot, gagb, omega )
         !
         DEALLOCATE( gagb )
         DEALLOCATE( drhot )
         !
      END IF
      !    
      !     forces on ions, ionic term in reciprocal space
      !     
      ALLOCATE( fion1( 3, nat ) )
      !
      fion1 = 0.d0
      !
      IF( tprnfor .OR. tfor .OR. tpre) THEN
          vtemp( 1:dfftp%ngm ) = rhog( 1:dfftp%ngm, 1 )
          IF( nspin == 2 ) THEN
             vtemp( 1:dfftp%ngm ) = vtemp(1:dfftp%ngm) + rhog( 1:dfftp%ngm, 2 )
          END IF
          CALL force_loc( .false., vtemp, fion1, rhops, vps, ei1, ei2, ei3, sfac, omega, screen_coul )
      END IF
      !
      !     calculation hartree + local pseudo potential
      !
      !
      IF (gstart == 2) vtemp(1)=(0.d0,0.d0)

!$omp parallel default(shared), private(ig,is)
!$omp do
      DO ig=gstart,dfftp%ngm
         vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*gg(ig))
      END DO
      !
      DO is=1,nsp
!$omp do
         DO ig=1,dffts%ngm
            vtemp(ig)=vtemp(ig)+sfac(ig,is)*vps(ig,is)
         END DO
      END DO
!$omp end parallel

      DEALLOCATE (rhotmp)
!
!     vtemp = v_loc(g) + v_h(g)
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
      !
      ! ... UGLY HACK WARNING: rhor must be saved before exch_corr_h
      ! ... overwrites it and before add_cc adds to it the core charge
      ! ... We also need an allocated rhoc array even in absence of core charge
      !
      IF ( dft_is_nonlocc() ) THEN
         ALLOCATE ( rhosave(dfftp%nnr,nspin),  rhocsave(dfftp%nnr) )
         rhosave(:,:) = rhor(:,:)
         IF ( SIZE(rhoc) == dfftp%nnr ) THEN
            rhocsave(:)= rhoc(:)
         ELSE
            rhocsave(:)= 0.0_dp
         ENDIF
      END IF
      !
      IF ( nlcc_any ) CALL add_cc( rhoc, rhog, rhor )
      CALL exch_corr_h( nspin, rhog, rhor, rhoc, sfac, exc, dxc, self_exc )
      !
      ! ... add non local corrections (if any)
      !
      IF ( dft_is_nonlocc() ) THEN
         !
         ! ... UGLY HACK WARNING: nlc adds nonlocal term IN RYDBERG A.U.
         ! ... to input energy and potential (potential is stored in rhor)
         !
         enlc = 0.0_dp
         vtxc = 0.0_dp
         rhor = rhor*e2
         CALL nlc( rhosave, rhocsave, nspin, enlc, vtxc, rhor )
         rhor = rhor/e2
         CALL mp_sum( enlc, intra_bgrp_comm )
         exc = exc + enlc / e2
         !
         ! ... non-local XC contribution to stress brought to crystal axis,
         ! ... transformed into energy derivative (Ha), added to dxc
         !
         IF ( tpre ) THEN
             CALL mp_sum( vtxc, intra_bgrp_comm )
             DO i=1,3
                DO j=1,3
                   dxc( i, j ) = dxc( i, j ) + (enlc - vtxc)/e2*ainv( j, i )
                END DO
             END DO
             denlc(:,:) = 0.0_dp
             inlc = get_inlc()
             !
             !^^ ... TEMPORARY FIX (newlsda) ...
             IF ( nspin==2 ) THEN
               rhosave(:,1) = rhosave(:,1) + rhosave(:,2) 
               CALL errore('stres_vdW', 'LSDA+stress+vdW-DF not implemented',1)
             END IF
             !^^.......................
             !   
             if (inlc==1 .or. inlc==2) then
               CALL stress_vdW_DF(rhosave(:,1), rhocsave, nspin, denlc )
             elseif (inlc == 3) then
               CALL stress_rVV10(rhosave(:,1), rhocsave, nspin, denlc )
             end if
             !
             dxc(:,:) = dxc(:,:) - omega/e2 * MATMUL(denlc,TRANSPOSE(ainv))
         END IF
         DEALLOCATE ( rhocsave, rhosave )
      ELSE
         denlc(:,:) = 0.0_dp
      END IF
      !
      !     Add TS-vdW wavefunction forces to rhor here... (RAD)
      !
      IF (ts_vdw) THEN
        !
        CALL fftx_add_field( rhor, UtsvdW, dfftp )
        !
      END IF
      !
      !     add plugin contributions to potential here... 
      !
      CALL plugin_add_potential( rhor )
      !
!
!     rhor contains the xc potential in r-space
!
!     ===================================================================
!     fourier transform of xc potential to g-space (dense grid)
!     -------------------------------------------------------------------
!
      IF( abivol .or. abisur ) THEN
         CALL rho_r2g ( dfftp, rhor, rhog, v_vol )
      ELSE
         CALL rho_r2g ( dfftp, rhor, rhog )
      END IF
       
      IF( nspin == 1 ) THEN
         rhog( 1:dfftp%ngm, 1 ) = rhog( 1:dfftp%ngm, 1 ) + vtemp(1:dfftp%ngm) 
      ELSE
         isup=1
         isdw=2
         rhog( 1:dfftp%ngm, isup ) = rhog( 1:dfftp%ngm, isup ) + vtemp(1:dfftp%ngm) 
         rhog( 1:dfftp%ngm, isdw ) = rhog( 1:dfftp%ngm, isdw ) + vtemp(1:dfftp%ngm) 
         IF( ttsic ) THEN
            rhog( 1:dfftp%ngm, isup ) = rhog( 1:dfftp%ngm, isup ) - self_vloc(1:dfftp%ngm) 
            rhog( 1:dfftp%ngm, isdw ) = rhog( 1:dfftp%ngm, isdw ) - self_vloc(1:dfftp%ngm) 
         END IF
      END IF

      DEALLOCATE (vtemp)
      IF( ttsic ) DEALLOCATE( self_vloc )
!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
      IF( tprnfor .OR. tfor ) THEN

         IF ( nlcc_any ) CALL force_cc( irb, eigrb, rhor, fion1 )

         CALL mp_sum( fion1, intra_bgrp_comm )
         !
         !    add g-space ionic and core correction contributions to fion
         !
         fion = fion + fion1
         !
         !    Add TS-vdW ion forces to fion here... (RAD)
         !
         IF (ts_vdw) THEN
            fion1(:,:) = FtsvdW(:,ind_srt(:))
            fion = fion + fion1
            !fion=fion+FtsvdW
         END IF
         !
         !     plugin patches on internal forces
         !
         CALL plugin_int_forces(fion)

      END IF

      DEALLOCATE( fion1 )
!
!
!     ===================================================================
!     fourier transform of total potential to r-space (dense grid)
!     -------------------------------------------------------------------
      
      CALL rho_g2r( dfftp, rhog, rhor )

      IF(nspin.EQ.1) THEN
         vave=SUM(rhor(:,1))/DBLE( dfftp%nr1* dfftp%nr2* dfftp%nr3)
      ELSE
         isup=1
         isdw=2
         vave=(SUM(rhor(:,isup))+SUM(rhor(:,isdw))) / 2.0d0 / DBLE(  dfftp%nr1 *  dfftp%nr2 *  dfftp%nr3 )
      END IF

      CALL mp_sum( vave, intra_bgrp_comm )

      !
      !     fourier transform of total potential to r-space (smooth grid)
      !
      CALL rho_g2r ( dffts, rhog, rhos )

      IF( dft_is_meta() ) CALL vofrho_meta( )

      ebac = 0.0d0
      !
      eht = eh * omega + esr - eself
      !
      eextfor = 0.0_DP
      IF( textfor ) eextfor = compute_eextfor( tau0 )
      !
      !     etot is the total energy ; ekin, enl were calculated in rhoofr
      !
      etot = ekin + eht + epseu + enl + exc + ebac +e_hubbard + eextfor
      !
      !     Add EXX energy to etot here. exx_wf related
      !
      IF(dft_is_hybrid().AND.exx_is_active()) THEN
        !
        etot = etot - exxalfa*exx 
        !
      END IF
      !
      !     Add TS-vdW energy to etot here... (RAD)
      !
      IF (ts_vdw)  etot=etot+EtsvdW
      !
      if (abivol) etot = etot + P_ext*volclu
      if (abisur) etot = etot + Surf_t*surfclu
      !
      !     Add possible external contribution from plugins to the energy
      !
      etot = etot + plugin_etot 
      !
      IF( tpre ) THEN
         !
         detot6 = dekin6 + dh6 + dps6 + dsr6
         !
         call mp_sum( detot6, intra_bgrp_comm )
         !
         DO k = 1, 6
            detmp( alpha(k), beta(k) ) = detot6(k)
            detmp( beta(k), alpha(k) ) = detmp( alpha(k), beta(k) )
         END DO
         !
         detot = MATMUL( detmp(:,:), TRANSPOSE( ainv(:,:) ) )
         !
         detot = detot + denl + dxc
         !
         !     Add TS-vdW cell derivatives to detot here... (RAD)
         !
         IF (ts_vdw) THEN
           !
           detot = detot + HtsvdW
           !
           ! BS / RAD start print ts_vdW pressure ------------- 
           !
           IF(MOD(nfi,iprint_stdout).EQ.0)  THEN
             detmp = HtsvdW 
             detmp = -1.d0 * (MATMUL( detmp(:,:), TRANSPOSE(h) )) / omega
             detmp = detmp * au_gpa   ! GPa 
             WRITE( stdout,9013) (detmp(1,1)+detmp(2,2)+detmp(3,3))/3.0_DP , nfi
             9013  FORMAT (/,5X,'TS-vdW Pressure (GPa)',F15.5,I7)
           END IF
           !
         END IF
         !
         ! BS / RAD / HK
         ! Adding the stress tensor from exact exchange here. exx_wf related
         !
         IF(dft_is_hybrid().AND.exx_is_active()) THEN
           !
           IF (isotropic .and. (ibrav.eq.1)) THEN
             !
             ! BS / RAD
             ! This part is dE/dV; so works only for cubic cells and isotropic change
             ! in simulation cell while doing variable cell calculation .. 
             ! dE/dV = -(1/3) * (-exx * 0.25_DP) / V
             ! dexx(3,3) = dE/dh = (dE/dV) * (dV/dh) = (dE/dV) * V * (Transpose h)^-1
             !
             DO k = 1, 6
               IF(alpha(k) .EQ. beta(k)) THEN
                 detmp( alpha(k), beta(k) ) = - (1.0_DP/3.0_DP) * (-exx * exxalfa) 
               ELSE
                 detmp( alpha(k), beta(k) ) = 0.0_DP
               END IF
               detmp( beta(k), alpha(k) ) = detmp( alpha(k), beta(k) )
             END DO
             !
             dexx = MATMUL( detmp(:,:), TRANSPOSE( ainv(:,:) ) )
             !
           ELSE
             !
             ! HK: general case is computed in exx_gs.f90
             !    (notice that the negative sign comes from the definition of
             !     exchange energy being as positive value)
             !
             dexx = -dexx_dh*exxalfa
             !
           END IF
           !
           detot = detot + dexx 
           !
           IF(MOD(nfi,iprint_stdout).EQ.0) WRITE( stdout,9014) (-1.0_DP/3.0_DP)*&
               (dexx(1,1)+dexx(2,2)+dexx(3,3))*au_gpa , nfi
           9014  FORMAT (5X,'EXX Pressure (GPa)',F15.5,I7)
           !
         END IF  
         ! BS / RAD : stress tensor from exx ends here ... 
         !
         ! BS / RAD start print total electronic pressure ------------- 
         !
         IF(MOD(nfi,iprint_stdout).EQ.0)  THEN
           detmp = detot
           detmp = -1.d0 * (MATMUL( detmp(:,:), TRANSPOSE(h) )) / omega
           detmp = detmp * au_gpa   ! GPa 
           WRITE( stdout,9015) (detmp(1,1)+detmp(2,2)+detmp(3,3))/3.0_DP , nfi
           9015  FORMAT (5X,'Total Electronic Pressure (GPa)',F15.5,I7)
         END IF 
         !
      END IF
      !
      CALL stop_clock( 'vofrho' )
      !
      IF ( tpre ) THEN
         !
         IF( ( iverbosity > 1 ) .AND. ( MOD( nfi - 1, iprint) == 0 ) ) THEN  
            !
            WRITE( stdout,*)
            WRITE( stdout,*) "From vofrho:"
            WRITE( stdout,*) "cell parameters h"
            WRITE( stdout,5555) (at(i,1)*alat, at(i,2)*alat, at(i,3)*alat,i=1,3)
            !
            WRITE( stdout,*)
            WRITE( stdout,*) "derivative of e(tot)"
            WRITE( stdout,5555) ((detot(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( detot, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*)
            WRITE( stdout,*) "derivative of e(kin)"
            WRITE( stdout,5555) ((dekin(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dekin, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(h)"
            WRITE( stdout,5555) ((dh(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dh, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
             !
            WRITE( stdout,*) "derivative of e(sr)"
            WRITE( stdout,5555) ((dsr(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dsr, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(ps)"
            WRITE( stdout,5555) ((dps(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dps, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(nl)"
            WRITE( stdout,5555) ((denl(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( denl, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(xc)"
            WRITE( stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dxc, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            IF (ts_vdw) THEN
               WRITE( stdout,*) "derivative of e(TS-vdW)"
               WRITE( stdout,5555) ((HtsvdW(i,j),j=1,3),i=1,3)
               WRITE( stdout,*) "kbar"
               detmp = -1.0d0 * MATMUL( HtsvdW, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
               WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            END IF
         ENDIF
      ENDIF

      RETURN

5555  FORMAT(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!

END SUBROUTINE vofrho_x
