!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE set_rhoc
  !-----------------------------------------------------------------------
  !
  !    This routine computes the core charge on the real space 3D mesh
  !    (rho_core) and in reciprocal space (rhog_core)
  !
  !
  USE kinds,     ONLY : dp
  USE io_global, ONLY : stdout
  USE atom,      ONLY : msh, rgrid
  USE uspp_param,ONLY : upf
  USE ions_base, ONLY : ntyp => nsp
  USE cell_base, ONLY : omega, tpiba2
  USE fft_base,  ONLY : dfftp
  USE fft_rho,   ONLY : rho_g2r
  USE gvect,     ONLY : ngm, ngl, gl, igtongl, ecutrho
  USE vlocal,    ONLY : strf
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum, mp_max
  USE scf,       ONLY : rho_core, rhog_core
  USE cellmd,    ONLY : cell_factor
  USE rhoc_mod,  ONLY : init_tab_rhc, interp_rhc
  !
  IMPLICIT NONE
  !
  REAL(DP) , ALLOCATABLE ::  rhocg(:)
  ! the radial fourier transform
  REAL(DP) ::  qmax, rhoneg
  ! used to check the core charge
  INTEGER :: ir, nt, ng
  ! counter on mesh points
  ! counter on atomic types
  ! counter on g vectors

  rhog_core(:) = 0.0_DP
  rho_core(:)  = 0.0_DP

  IF ( ANY( upf(1:ntyp)%nlcc ) ) THEN
     !
     qmax = tpiba2 * MAXVAL ( gl )
     CALL mp_max (qmax, intra_bgrp_comm)
     !! this is the actual maximum |G|^2 needed in the interpolation table
     !! for variable-cell calculations. It may exceed ecutrho, so we use
     !! "cell_factor" (1.2 or so) as below, in order to avoid too frequent
     !! re-allocations of the interpolation table
     !
     qmax = MAX (sqrt(qmax), sqrt(ecutrho)*cell_factor)
     CALL init_tab_rhc  ( qmax, omega, intra_bgrp_comm, ir )
     !
     ALLOCATE (rhocg( ngl))
     !$acc data create(rhocg) copyin(gl, strf, rhog_core, rho_core) present(igtongl)
     !
     !    the sum is on atom types
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%nlcc ) THEN
           !
           ! interp_rhc computes the radial fourier transform for each shell of g vec
           !
           CALL interp_rhc (nt, ngl, gl, tpiba2, rhocg)
           !
           !     multiply by the structure factor and sum
           !
           !$acc parallel loop
           DO ng = 1, ngm
              rhog_core(ng) = rhog_core(ng) + strf(ng,nt) * rhocg(igtongl(ng))
           END DO
        ENDIF
     ENDDO
     !
     CALL rho_g2r( dfftp, rhog_core, rho_core )
     !
     !    test on the charge and computation of the core energy
     !
     rhoneg = 0.d0
     !$acc parallel loop reduction(+:rhoneg)
     DO ir = 1, dfftp%nnr
        rhoneg = rhoneg + min (0.d0, rho_core (ir) )
        !
        ! NOTE: Core charge is computed in reciprocal space and brought to real
        ! space by FFT. For non smooth core charges (or insufficient cut-off)
        ! this may result in negative values in some grid points.
        ! Up to October 1999 the core charge was forced to be positive definite.
        ! This induces an error in force, and probably stress, calculation if
        ! the number of grid points where the core charge would be negative is
        ! large. The error disappears for sufficiently high cut-off, but may be
        ! rather large and it is better to leave the core charge as it is.
        ! If you insist to have it positive definite (with the possible problems
        ! mentioned above) uncomment the following lines.  SdG, Oct 15 1999
        !
        !         rho_core(ir) = MAX (0.0_dp, rho_core(ir))
     ENDDO
     !
     rhoneg = rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     CALL mp_sum(  rhoneg, intra_bgrp_comm )
     !
     IF (rhoneg < -1.0d-6) &
        WRITE(stdout,'(/5x,"Check: negative core charge=",2f12.6)') rhoneg
     !
     ! calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
     ! The term was present in previous versions of the code but it shouldn't
     !
     !   call create_scf_type(dum)
     !   dum%of_r(:,:) = 0.0_DP
     !   dum%of_g(:,:) = (0.0_DP, 0.0_DP)
     !   
     !   call v_xc( dum, rho_core, rhog_core, etxcc, vtxcc, aux )
     ! 
     !   call destroy_scf_type(dum)
     !   WRITE( stdout, 9000) etxcc
     ! 9000 format (5x,'core-only xc energy         = ',f15.8,' Ry')
     !   WRITE( stdout,  * ) 'BEWARE it will be subtracted from total energy !'
     !
     !$acc exit data copyout(rho_core, rhog_core)
     !$acc end data
     DEALLOCATE (rhocg)
  END IF
  !
  RETURN

END SUBROUTINE set_rhoc
