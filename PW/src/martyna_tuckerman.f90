!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef TESTING
!---------------------------------------------------------------------------
MODULE martyna_tuckerman
  !-------------------------------------------------------------------------
  !! The variables and the routines needed by the Martyna-Tuckerman method
  !! for isolated systems.
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : e2, pi, tpi, fpi
  USE ws_base
  !
  IMPLICIT NONE
  !
  TYPE(ws_type) :: ws
  !! Wigner-Seitz cell
  REAL(DP) :: alpha
  !! convergence parameter [it corresponds to \(\alpha_\text{ewd}\) of
  !! J.Chem.Phys. 110, 2810 (1999)]
  REAL(DP) :: beta
  !! convergence parameter
  REAL(DP), ALLOCATABLE :: wg_corr(:)
  !! W(G) function [see Eq.(4.3) of J.Chem.Phys. 110, 2810 (1999)]
  LOGICAL :: wg_corr_is_updated = .FALSE.
  !! if .FALSE. initialize wg correction stuff
  LOGICAL :: do_comp_mt = .FALSE.
  !! if .TRUE. apply MT corrections
  LOGICAL :: gamma_only = .FALSE.
  !! .TRUE. if Gamma point only
  INTEGER :: gstart = 1
  !! index of the first G vector whose module is > 0
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: tag_wg_corr_as_obsolete, do_comp_mt, &
            wg_corr_ewald, wg_corr_loc, wg_corr_h, wg_corr_force
  !
CONTAINS
  !----------------------------------------------------------------------------
  SUBROUTINE tag_wg_corr_as_obsolete
    !----------------------------------------------------------------------------
    !! wg_corr_is_updated = .FALSE.
    !
    wg_corr_is_updated = .FALSE.
    !
  END SUBROUTINE tag_wg_corr_as_obsolete
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE wg_corr_h( omega, ngm, rho, v, eh_corr )
    !----------------------------------------------------------------------------
    !! Calculates the Hartree contribution to the energy and potential.
    !
    INTEGER, INTENT(IN) :: ngm
    !! local number of G vectors
    REAL(DP), INTENT(IN) :: omega
    !! cell volume
    COMPLEX(DP), INTENT(IN) :: rho(ngm)
    !! charge density
    COMPLEX(DP), INTENT(OUT) :: v(ngm)
    !! Hartree contribution to the potential
    REAL(DP), INTENT(OUT) :: eh_corr
    !! Hartree contribution to the energy
    !
    ! ... local variables
    !
    INTEGER :: ig
    !
    IF (.NOT. wg_corr_is_updated) CALL init_wg_corr
    !
    v(:) = (0._dp,0._dp)
    !
    eh_corr =  0._dp
    DO ig = 1, ngm
       v(ig) = e2 * wg_corr(ig) * rho(ig) 
       eh_corr = eh_corr + ABS(rho(ig))**2 * wg_corr(ig)
    ENDDO
    IF (gamma_only) v(gstart:ngm) = 0.5_dp * v(gstart:ngm)
    !
    eh_corr = 0.5_dp * e2 * eh_corr * omega
    !
    RETURN
    !
  END SUBROUTINE wg_corr_h
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE wg_corr_loc( omega, ntyp, ngm, zv, strf, v )
    !----------------------------------------------------------------------------
    !! Calculates the local contribution to the potential.
    !
    INTEGER, INTENT(IN) :: ntyp
    !! number of species
    INTEGER, INTENT(IN) :: ngm
    !! local number of G vectors
    REAL(DP), INTENT(IN) :: omega
    !! cell volume
    REAL(DP), INTENT(IN) :: zv(ntyp)
    !! (pseudo-) atomic charge
    COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
    !! the structure factor
    COMPLEX(DP), INTENT(OUT) :: v(ngm)
    !! local contribution to the potential
    !
    ! ... local variables
    !
    INTEGER :: ig
    !
    IF (.NOT. wg_corr_is_updated) CALL init_wg_corr
    !
    DO ig = 1, ngm
       v(ig) = - e2 * wg_corr(ig) * SUM(zv(1:ntyp)*strf(ig,1:ntyp)) / omega
    ENDDO
    IF (gamma_only) v(gstart:ngm) = 0.5_dp * v(gstart:ngm)
    !
    RETURN
    !
  END SUBROUTINE wg_corr_loc
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE wg_corr_force( lnuclei, omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                            rho, force )
    !----------------------------------------------------------------------------
    !! MT correction to the force.
    !
    USE cell_base,   ONLY : tpiba
    USE mp_bands,    ONLY : intra_bgrp_comm
    USE mp,          ONLY : mp_sum
    !
    INTEGER, INTENT(IN) :: nat
    !! number of atoms
    INTEGER, INTENT(IN) :: ntyp
    !! number of species
    INTEGER, INTENT(IN) :: ityp(nat)
    !! indices for atomic species
    INTEGER, INTENT(IN) :: ngm
    !! local number of G vectors
    REAL(DP), INTENT(IN) :: omega
    !! cell volume
    REAL(DP), INTENT(IN) :: zv(ntyp)
    !! (pseudo-) atomic charge
    REAL(DP), INTENT(IN) :: tau(3,nat)
    !! atomic positions read from stdin (in Bohr)
    REAL(DP), INTENT(IN) :: g(3,ngm)
    !! G-vectors
    COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
    !! the structure factor
    COMPLEX(DP), INTENT(IN) :: rho(ngm)
    !! the charge density
    LOGICAL, INTENT(IN) :: lnuclei
    !! this variable is used to select if the correction should be done
    !! on rho and nuclei or only on rho.
    REAL(DP), INTENT(OUT) :: force(3,nat)
    !! the correction to the force
    !
    ! ... local variables
    !
    INTEGER :: ig, na
    REAL(DP) :: arg
    COMPLEX(DP), ALLOCATABLE :: v(:)
    COMPLEX(DP) :: rho_tot
    !
    IF (.NOT. wg_corr_is_updated) CALL init_wg_corr
    !
    ALLOCATE( v(ngm) )
    DO ig = 1, ngm
       rho_tot = rho(ig)
       IF (lnuclei) rho_tot = rho_tot - SUM(zv(1:ntyp)*strf(ig,1:ntyp)) / omega
       v(ig) = e2 * wg_corr(ig) * rho_tot
    ENDDO
    force(:,:) = 0._dp
    DO na = 1, nat
       DO ig = 1, ngm
          arg = tpi * SUM ( g(:,ig)*tau(:, na) ) 
          force(:,na) = force(:,na) + g(:,ig) * CMPLX(SIN(arg),-COS(ARG), KIND=dp) * v(ig)
       ENDDO
       force(:,na) = - force(:,na) * zv(ityp(na))  * tpiba
    ENDDO
    DEALLOCATE( v )
    !
    CALL mp_sum( force, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE wg_corr_force
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE init_wg_corr
    !----------------------------------------------------------------------------
    !! Calculates W(g) function [see Eq.(4.3) of J.Chem.Phys. 110, 2810 (1999)].
    !
    USE mp_bands,         ONLY : me_bgrp
    USE fft_base,         ONLY : dfftp
    USE fft_interfaces,   ONLY : invfft
    USE fft_rho,          ONLY : rho_r2g
    USE fft_types,        ONLY : fft_index_to_3d
    USE control_flags,    ONLY : gamma_only_ => gamma_only
    USE gvect,            ONLY : ngm, gg, gstart_ => gstart, ecutrho
    USE cell_base,        ONLY : at, alat, tpiba2, omega
    !
    ! ... local variables
    !
    INTEGER :: ir, i,j,k, ig, nt
    LOGICAL :: offrange
    REAL(DP) :: r(3), rws, upperbound, rws2
    REAL(DP), ALLOCATABLE :: auxr(:)
    COMPLEX(DP), ALLOCATABLE :: aux(:), auxg(:,:)
    !
#ifdef TESTING
    REAL(DP), ALLOCATABLE :: plot(:)
    CHARACTER (LEN=25) :: filplot
    LOGICAL, SAVE :: first = .TRUE.
#endif
    !
    IF ( ALLOCATED(wg_corr) ) DEALLOCATE( wg_corr )
    ALLOCATE( wg_corr(ngm) )
    !
    ! ... choose alpha in order to have convergence in the sum over G.
    ! upperbound is a safe upper bound for the error in the sum over G.
    !
    alpha = 2.9d0
    upperbound = 1._dp
    DO WHILE ( upperbound > 1.e-7_dp) 
       alpha = alpha - 0.1_dp  
       IF (alpha<=0._dp) CALL errore( 'init_wg_corr', 'optimal alpha not found', 1 )
       upperbound = e2 * SQRT(2.d0 * alpha / tpi) * &
                         erfc( SQRT( ecutrho / 4.d0 / alpha) )
    ENDDO
    beta = 0.5_dp/alpha ! 1._dp/alpha
    ! write (*,*) " alpha, beta MT = ", alpha, beta
    !
    CALL ws_init( at, ws )
    !
    gstart = gstart_
    gamma_only = gamma_only_
    !
    ALLOCATE( auxr(dfftp%nnr), auxg(dfftp%nnr,1) )
    auxr = (0._dp,0._dp)
    !
    DO ir = 1, dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
       !
       ! ... three dimensional indexes
       !
       CALL fft_index_to_3d (ir, dfftp, i,j,k, offrange)
       IF ( offrange ) CYCLE
       !
       r(:) = ( at(:,1)/dfftp%nr1*i + at(:,2)/dfftp%nr2*j + at(:,3)/dfftp%nr3*k )
       !
       rws = ws_dist(r,ws)
#ifdef TESTING
       rws2 = ws_dist_stupid(r,ws)
       IF (ABS(rws-rws2) > 1.e-5 ) THEN
          WRITE(*,'(4i8)') ir, i,j,k
          WRITE(*,'(5f14.8)') r(:), rws, rws2
          STOP
       ENDIF
#endif
       !
       auxr(ir) = smooth_coulomb_r( rws*alat )
       !
    ENDDO
    !
    CALL rho_r2g( dfftp, auxr, auxg )
    !
    DO ig = 1, ngm
       wg_corr(ig) = omega * REAL(auxg(ig,1)) - smooth_coulomb_g( tpiba2*gg(ig))
    ENDDO
    wg_corr(:) =  wg_corr(:) * EXP(-tpiba2*gg(:)*beta/4._dp)**2
    !
    IF (gamma_only) wg_corr(gstart:ngm) = 2.d0 * wg_corr(gstart:ngm)
    !
    wg_corr_is_updated = .TRUE.
    !
#ifdef TESTING
    IF (first) THEN
       ALLOCATE( aux(dfftp%nnr), plot(dfftp%nnr) )
       !
       filplot = 'wg_corr_r'
       CALL invfft( 'Rho', aux, dfftp )
       plot(:) = REAL(aux(:))
       CALL write_wg_on_file( filplot, plot )
       !
       filplot = 'wg_corr_g'
       aux(:) = (0._dp,0._dp)
       DO ig = 1, ngm
          aux(dfftp%nl(ig)) = smooth_coulomb_g( tpiba2*gg(ig))/omega
       ENDDO
       IF (gamma_only) aux(dfftp%nlm(1:ngm)) = CONJG( aux(dfftp%nl(1:ngm)) )
       !
       CALL invfft( 'Rho', aux, dfftp )
       plot(:) = REAL(aux(:))
       CALL write_wg_on_file( filplot, plot )
       !
       filplot = 'wg_corr_diff'
       aux(:) = (0._dp,0._dp)
       aux(dfftp%nl(1:ngm)) = wg_corr(1:ngm) / omega
       IF (gamma_only) THEN
          aux(:) = 0.5_dp * aux(:) 
          aux(dfftp%nlm(1:ngm)) = aux(dfftp%nlm(1:ngm)) + CONJG( aux(dfftp%nl(1:ngm)) )
       ENDIF
       CALL invfft( 'Rho', aux, dfftp )
       plot(:) = REAL(aux(:))
       CALL write_wg_on_file( filplot, plot )
       !
       DEALLOCATE( aux, plot )
       !
       first = .FALSE.
    ENDIF
#endif
    !
    DEALLOCATE( auxr, auxg )
    !
    RETURN
    !
  END SUBROUTINE init_wg_corr 
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE write_wg_on_file( filplot, plot )
    !----------------------------------------------------------------------------
    !! Write W(g) function on file.
    !
    USE fft_base,        ONLY : dfftp
    USE gvect,           ONLY : gcutm
    USE gvecw,           ONLY : ecutwfc
    USE gvecs,           ONLY : dual
    USE cell_base,       ONLY : at, alat, tpiba2, omega, ibrav, celldm
    USE ions_base,       ONLY : zv, ntyp => nsp, nat, ityp, atm, tau
    !
    CHARACTER (LEN=25), INTENT(IN) :: filplot
    !! name of the plot file
    REAL(DP) :: plot(dfftp%nnr)
    !! plot data
    !
    ! ... local variables
    !
    CHARACTER (LEN=75) :: title
    INTEGER :: plot_num=0, iflag=+1
    !
    CALL plot_io( filplot, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
       dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
       gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, plot, iflag )
    !
    RETURN
    !
  END SUBROUTINE write_wg_on_file
  !
  !
  !----------------------------------------------------------------------------
  REAL(DP) FUNCTION wg_corr_ewald( omega, ntyp, ngm, zv, strf )
    !----------------------------------------------------------------------------
    !! MT Ewald correction.
    !
    INTEGER, INTENT(IN) :: ntyp, ngm
    REAL(DP), INTENT(IN) :: omega, zv(ntyp)
    COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
    INTEGER :: ig
    COMPLEX(DP)  :: rhoion
    !
    IF (.NOT. wg_corr_is_updated) CALL init_wg_corr
    !
    wg_corr_ewald = 0._dp
    DO ig = 1, ngm
       rhoion = SUM(zv(1:ntyp)* strf(ig,1:ntyp) ) / omega
       wg_corr_ewald = wg_corr_ewald + ABS(rhoion)**2 * wg_corr(ig) 
    ENDDO
    wg_corr_ewald = 0.5_dp * e2 * wg_corr_ewald * omega
    ! write(*,*) "ewald correction   = ", wg_corr_ewald
    !
  END FUNCTION wg_corr_ewald
  !
  !
  !----------------------------------------------------------------------------
  REAL(DP) FUNCTION smooth_coulomb_r( r )
    !----------------------------------------------------------------------------
    !! Smooth factor for Coulomb potential in r-space [see Eq.(4.1) of J.Chem.
    !! Phys. 110, 2810 (1999)].
    !
    REAL(DP), INTENT(IN) :: r
    !
    !  smooth_coulomb_r = SQRT(2._dp*alpha/tpi)**3 * exp(-alpha*r*r) ! to be modified
    IF (r > 1.e-6_dp) THEN
       smooth_coulomb_r = erf(SQRT(alpha)*r)/r
    ELSE
       smooth_coulomb_r = 2._dp/SQRT(pi) * SQRT(alpha)
    ENDIF
    !
  END FUNCTION smooth_coulomb_r
  !
  !
  !----------------------------------------------------------------------------
  REAL(DP) FUNCTION smooth_coulomb_g( q2 )
    !----------------------------------------------------------------------------
    !! Smooth factor for Coulomb potential in g-space [see Eq.(4.3) of J.Chem.
    !! Phys. 110, 2810 (1999)].
    !
    REAL(DP), INTENT(IN) :: q2
    !
    !  smooth_coulomb_g = exp(-q2/4._dp/alpha) ! to be modified
    IF (q2 > 1.e-6_dp) THEN
       smooth_coulomb_g = fpi * EXP(-q2/4._dp/alpha)/q2 ! to be modified
    ELSE 
       smooth_coulomb_g = - 1._dp * fpi * (1._dp/4._dp/alpha + 2._dp*beta/4._dp)
    ENDIF
    !
  END FUNCTION smooth_coulomb_g
  !
  !
END MODULE martyna_tuckerman
