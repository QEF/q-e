!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE xc_lda_lsda
!
USE kinds,     ONLY: DP 
!
IMPLICIT NONE
!
PRIVATE
SAVE
!
!  LDA and LSDA exchange-correlation drivers
PUBLIC :: xc, xc_spin, select_functionals
!
PUBLIC :: iexch_l, icorr_l
PUBLIC :: exx_started_l, exx_fraction_l
PUBLIC :: is_there_finite_size_corr, finite_size_cell_volume_l
!
!  indexes defining xc functionals
INTEGER  :: iexch_l, icorr_l
!
!  variables for hybrid exchange and finite_size_cell_volume correction
LOGICAL  :: exx_started_l, is_there_finite_size_corr
REAL(DP) :: exx_fraction_l, finite_size_cell_volume_l
!
!
 CONTAINS
!
!
!----------------------------------------------------------------------------
!----- Select functionals by the corresponding indexes ----------------------
!----------------------------------------------------------------------------
SUBROUTINE select_functionals( iexch, icorr, exx_fraction, finite_size_cell_volume )
   !
   IMPLICIT NONE
   !
   INTEGER,  INTENT(IN) :: iexch, icorr
   REAL(DP), INTENT(IN), OPTIONAL :: exx_fraction, finite_size_cell_volume
   !
   ! exchange-correlation indexes
   iexch_l = iexch
   icorr_l = icorr
   !
   ! hybrid exchange vars
   exx_started_l  = .false.
   exx_fraction_l = 0._DP
   IF ( present(exx_fraction) ) THEN
      exx_started_l  = .true.
      exx_fraction_l = exx_fraction
   ENDIF
   !
   ! finite size correction vars
   is_there_finite_size_corr = .false.
   finite_size_cell_volume_l = -1.0_DP
   IF ( present(finite_size_cell_volume) ) THEN
      is_there_finite_size_corr = .true.
      finite_size_cell_volume_l = finite_size_cell_volume
   ENDIF
   !
   !
   RETURN
   !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
!-------  LDA-LSDA DRIVERS --------------------------------------------------
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE xc( length, rho, ex, ec, vx, vc )
  !--------------------------------------------------------------------------
  !     lda exchange and correlation functionals - Hartree a.u.
  !
  !     exchange   :  Slater, relativistic Slater
  !     correlation:  Ceperley-Alder (Perdew-Zunger parameters)
  !                   Vosko-Wilk-Nusair
  !                   Lee-Yang-Parr
  !                   Perdew-Wang
  !                   Wigner
  !                   Hedin-Lundqvist
  !                   Ortiz-Ballone (Perdew-Zunger formula)
  !                   Ortiz-Ballone (Perdew-Wang formula)
  !                   Gunnarsson-Lundqvist
  !
  !     input : rho=rho(r)
  !     definitions: E_x = \int E_x(rho) dr, E_x(rho) = rho\epsilon_c(rho)
  !                  same for correlation
  !     output: ex = \epsilon_x(rho) ( NOT E_x(rho) )
  !             vx = dE_x(rho)/drho  ( NOT d\epsilon_x(rho)/drho )
  !             ec, vc as above for correlation
  !
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc, ex, vx
  !
  ! local vars
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: ec_, vc_
  REAL(DP), PARAMETER :: small = 1.E-10_DP,  third = 1.0_DP / 3.0_DP, &
                         pi34 = 0.6203504908994_DP  ! pi34=(3/4pi)^(1/3)
  REAL(DP) :: rs(length)
  !
  !
  WHERE ( rho > small )
     rs = pi34 / rho**third
  ELSEWHERE
     rs = 1.0_DP
  END WHERE
  !
  IF ( icorr_l >= 12 .AND. icorr_l <= 14 ) THEN
     ALLOCATE( ec_(length) )
     ALLOCATE( vc_(length) )
  ENDIF
  !
  !
  ! ... EXCHANGE
  !
  SELECT CASE( iexch_l )
  CASE( 1 )                      ! 'sla'
     !
     CALL slater( length, rs, ex, vx )
     !
  CASE( 2 )                      ! 'sl1'
     !
     CALL slater1( length, rs, ex, vx )
     !
  CASE( 3 )                      ! 'rxc'
     !
     CALL slater_rxc( length, rs, ex, vx )
     !
  CASE( 4, 5 )                   ! 'oep','hf'
     !
     IF ( exx_started_l ) THEN
        ex = 0.0_DP
        vx = 0.0_DP
     ELSE
        CALL slater( length, rs, ex, vx )
     ENDIF
     !
  CASE( 6, 7 )                   ! 'pb0x' or 'DF-cx-0', or 'DF2-0',
     !                           ! 'B3LYP'
     CALL slater( length, rs, ex, vx )
     IF ( exx_started_l ) THEN
        ex = (1.0_DP - exx_fraction_l) * ex
        vx = (1.0_DP - exx_fraction_l) * vx
     ENDIF
     !
  CASE( 8 )                      ! 'sla+kzk'
     !
     IF (.NOT. is_there_finite_size_corr) CALL errore ('XC',&
          'finite size corrected exchange used w/o initialization',1)
     CALL slaterKZK( length, rs, ex, vx, finite_size_cell_volume_l )
     !
  CASE( 9 )                      ! 'X3LYP'
     !
     CALL slater( length, rs, ex, vx )
     IF ( exx_started_l ) THEN
        ex = (1.0_DP - exx_fraction_l) * ex
        vx = (1.0_DP - exx_fraction_l) * vx
     ENDIF
     !
  CASE DEFAULT
     !
     ex = 0.0_DP
     vx = 0.0_DP
     !
  END SELECT
  !
  !
  ! ... CORRELATION
  !
  SELECT CASE( icorr_l )
  CASE( 1 )
     !
     CALL pz( length, rs, 1, ec, vc )
     !
  CASE( 2 )
     !
     CALL vwn( length, rs, ec, vc )
     !
  CASE( 3 )
     !
     CALL lyp( length, rs, ec, vc )
     !
  CASE( 4 )
     !
     CALL pw( length, rs, 1, ec, vc )
     !
  CASE( 5 )
     !
     CALL wignerc( length, rs, ec, vc )
     !
  CASE( 6 )
     !
     CALL hl( length, rs, ec, vc )
     !
  CASE( 7 )
     !
     CALL pz( length, rs, 2, ec, vc )
     ! 
  CASE( 8 )
     !
     CALL pw( length, rs, 2, ec, vc )
     !
  CASE( 9 )
     !
     CALL gl( length, rs, ec, vc )
     !
  CASE( 10 )
     !
     IF (.NOT. is_there_finite_size_corr) CALL errore ('XC',&
          'finite size corrected exchange used w/o initialization',2)
     CALL pzKZK ( length, rs, ec, vc, finite_size_cell_volume_l )
     !
  CASE( 11 )
     !
     CALL vwn1_rpa( length, rs, ec, vc )
     !
  CASE( 12 )                ! 'B3LYP'
     !
     CALL vwn( length, rs, ec, vc )
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     !
     CALL lyp( length, rs, ec_, vc_ )
     ec = ec + 0.81_DP * ec_
     vc = vc + 0.81_DP * vc_
     !
  CASE( 13 )                ! 'B3LYP-V1R'
     !
     CALL vwn1_rpa ( length, rs, ec, vc )
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     !
     CALL lyp( length, rs, ec_, vc_ )
     ec = ec + 0.81_DP * ec_
     vc = vc + 0.81_DP * vc_
     !
  CASE( 14 )                ! 'X3LYP'
     !
     CALL vwn1_rpa( length, rs, ec, vc )
     ec = 0.129_DP * ec
     vc = 0.129_DP * vc
     !
     CALL lyp( length, rs, ec_, vc_ )
     ec = ec + 0.871_DP * ec_
     vc = vc + 0.871_DP * vc_
     !
  CASE DEFAULT
     !
     ec = 0.0_DP
     vc = 0.0_DP
     !
  END SELECT
  !
  !
  WHERE ( rho <= small )
     ex = 0.0_DP
     ec = 0.0_DP
     vx = 0.0_DP
     vc = 0.0_DP
  END WHERE
  !
  IF ( icorr_l >= 12 .AND. icorr_l <= 14 ) THEN
     DEALLOCATE( ec_ )
     DEALLOCATE( vc_ )
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE xc
!
!
!
SUBROUTINE xc_spin( length, rho, zeta, ex, ec, vx, vc )
  !---------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length)   :: rho, zeta
  REAL(DP), INTENT(OUT), DIMENSION(length)   :: ex, ec               !^^^
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vx, vc
  !
  !   local vars
  REAL(DP), DIMENSION(length) :: rs 
  REAL(DP), ALLOCATABLE, DIMENSION(:)   :: ec_
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: vc_
  !
  REAL(DP), PARAMETER :: small= 1.E-10_DP, third = 1.0_DP/3.0_DP, &
                         pi34= 0.6203504908994_DP ! pi34=(3/4pi)^(1/3)
  !
  WHERE ( rho > small )
     rs = pi34 / rho**third
  ELSEWHERE
     rs = 1.0_DP
  END WHERE
  !
  IF ( icorr_l >= 12 .AND. icorr_l <= 14 ) THEN
     ALLOCATE( ec_(length) )
     ALLOCATE( vc_(length,2) )
  ENDIF
  !
  !
  ! ... EXCHANGE
  !
  SELECT CASE( iexch_l )
  CASE( 1 )                                      ! 'sla'
     !
     CALL slater_spin( length, rho, zeta, ex, vx )
     !
  CASE( 2 )                                      ! 'sl1'
     !
     CALL slater1_spin( length, rho, zeta, ex, vx )
     !
  CASE( 3 )                                      ! 'rxc'
     !
     CALL slater_rxc_spin( length, rho, zeta, ex, vx )
     !
  CASE( 4, 5 )                                   ! 'oep','hf'
     !
     IF ( exx_started_l ) THEN
        ex = 0.0_DP
        vx = 0.0_DP
     ELSE
        CALL slater_spin( length, rho, zeta, ex, vx )
     ENDIF
     !
  CASE( 6 )                                      ! 'pb0x'
     !
     CALL slater_spin( length, rho, zeta, ex, vx )
     IF ( exx_started_l ) THEN
        ex = (1.0_DP - exx_fraction_l) * ex
        vx = (1.0_DP - exx_fraction_l) * vx
     ENDIF
     !
  CASE( 7 )                                      ! 'B3LYP'
     !
     CALL slater_spin( length, rho, zeta, ex, vx )
     IF ( exx_started_l ) THEN
        ex = (1.0_DP - exx_fraction_l) * ex
        vx = (1.0_DP - exx_fraction_l) * vx
     ENDIF
     !
  CASE( 9 )                                      ! 'X3LYP'
     !
     CALL slater_spin( length, rho, zeta, ex, vx )
     IF ( exx_started_l ) THEN
        ex = (1.0_DP - exx_fraction_l) * ex
        vx = (1.0_DP - exx_fraction_l) * vx
     ENDIF
     !
  CASE DEFAULT
     !
     ex = 0.0_DP
     vx = 0.0_DP
     !
  END SELECT
  !
  !
  ! ... CORRELATION
  !
  SELECT CASE( icorr_l )
  CASE( 0 )
     !
     ec = 0.0_DP
     vc = 0.0_DP
     !
  CASE( 1 )
     !
     CALL pz_spin( length, rs, zeta, ec, vc )
     !
  CASE( 2 )
     !
     CALL vwn_spin( length, rs, zeta, ec, vc )
     !
  CASE( 3 )
     !
     CALL lsd_lyp( length, rho, zeta, ec, vc )       ! from CP/FPMD (more_functionals)
     !
  CASE( 4 )
     !
     CALL pw_spin( length, rs, zeta, ec, vc )
     !
  CASE( 12 )                                         ! 'B3LYP'
     !
     CALL vwn_spin( length, rs, zeta, ec, vc )
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     !
     CALL lsd_lyp( length, rho, zeta, ec_, vc_ )   ! from CP/FPMD (more_functionals)
     ec = ec + 0.81_DP * ec_
     vc = vc + 0.81_DP * vc_
     !     
  CASE( 13 )                                         ! 'B3LYP-V1R'
     !
     CALL vwn1_rpa_spin( length, rs, zeta, ec, vc )
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     !
     CALL lsd_lyp( length, rho, zeta, ec_, vc_ )   ! from CP/FPMD (more_functionals)
     ec = ec + 0.81_DP * ec_
     vc = vc + 0.81_DP * vc_
     !
  CASE( 14 )                                         ! 'X3LYP
     !
     CALL vwn1_rpa_spin( length, rs, zeta, ec, vc )
     ec = 0.129_DP * ec
     vc = 0.129_DP * vc
     !
     CALL lsd_lyp( length, rho, zeta, ec_, vc_ )   ! from CP/FPMD (more_functionals)
     ec = ec + 0.871_DP * ec_
     vc = vc + 0.871_DP * vc_
     !
  CASE DEFAULT
     !
     CALL errore( 'lsda_functional (xc_spin)', 'not implemented', icorr_l )
     !
  END SELECT
  !
  !
  WHERE ( rho <= small )
     ex = 0.0_DP
     ec = 0.0_DP
     vx(:,1) = 0.0_DP
     vx(:,2) = 0.0_DP
     vc(:,1) = 0.0_DP
     vc(:,2) = 0.0_DP
  END WHERE
  !
  !
  RETURN
  !
END SUBROUTINE xc_spin
!
!
END MODULE xc_lda_lsda
