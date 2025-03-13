!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_rho_g( rhocg, nspina )
  !-----------------------------------------------------------------------
  !! Compute superposition of atomic charges in reciprocal space.
  !
  !! Three cases:
  !
  !! * if \(\text{nspina}=1\) the total atomic charge density is calculated;
  !! * if \(\text{nspina}=2\) collinear case. The total density is calculated
  !!               in the first component and the magnetization in 
  !!               the second;
  !! * if \(\text{nspina}=4\) noncollinear case. Total density in the first
  !!               component and magnetization vector in the
  !!               other three.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps8
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : tpiba2, omega
  USE gvect,                ONLY : ngm, ngl, gl, igtongl, ecutrho
  USE lsda_mod,             ONLY : starting_magnetization
  USE starting_scf,         ONLY : starting_charge
  USE vlocal,               ONLY : strf
  USE noncollin_module,     ONLY : angle1, angle2
  USE uspp_param,           ONLY : upf
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_max
  USE cellmd,               ONLY : cell_factor
  USE rhoat_mod,            ONLY : init_tab_rhoat, interp_rhoat
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspina
  !! number of spin components to be calculated. It may differ from
  !! nspin because in some cases the total charge only is needed, 
  !! even in a LSDA calculation.
  COMPLEX(DP), INTENT(OUT) :: rhocg(ngm,nspina)
  !! contains G-space components of the superposition of atomic charges
  !! contained in the array upf%rho_at (read from pseudopotential files).
  !
  ! ... local variables
  !
  REAL(DP) :: rhoscale, fac
  REAL(DP), ALLOCATABLE :: rhoatg(:)
  REAL(DP) :: angular(nspina)
  REAL(DP) :: qmax
  INTEGER :: ir, is, ig, igl, nt, ierr
  !
  qmax = tpiba2 * MAXVAL ( gl )
  CALL mp_max (qmax, intra_bgrp_comm)
  !! this is the actual maximum |G|^2 needed in the interpolation table
  !! for variable-cell calculations. It may exceed ecutrho, so we use
  !! "cell_factor" (1.2 or so) as below, in order to avoid too frequent
  !! re-allocations of the interpolation table
  !
  qmax = MAX (sqrt(qmax), sqrt(ecutrho)*cell_factor)
  CALL init_tab_rhoat (qmax, omega, intra_bgrp_comm, ierr)
  !! Initialize  interpolation tables (if not already done)
  !
  ALLOCATE (rhoatg( ngl))
  !$acc data create(rhoatg) copyin( gl, strf ) present ( igtongl )
  !
  !$acc kernels
  rhocg(:,1:nspina) = (0.0_dp, 0.0_dp)
  !$acc end kernels
  DO nt = 1, ntyp
     !
     ! interpolate atomic rho(G)
     !
     CALL interp_rhoat( nt, ngl, gl, tpiba2, rhoatg )
     !
     IF (upf(nt)%zp > eps8) THEN
        rhoscale = MAX(0.0_dp, upf(nt)%zp - starting_charge(nt)) / upf(nt)%zp
     ELSE
        rhoscale = 1.0_dp
     ENDIF
     !
     !$acc parallel loop
     DO ig = 1, ngm
        rhocg(ig,1) = rhocg(ig,1) + &
                strf(ig,nt) * rhoscale * rhoatg(igtongl(ig))
     ENDDO
     !
     IF ( nspina >= 2 ) THEN
        !
        angular(1) = 1._dp
        IF ( nspina == 4 ) THEN
           angular(1) = sin(angle1(nt))*cos(angle2(nt))
           angular(2) = sin(angle1(nt))*sin(angle2(nt))
           angular(3) = cos(angle1(nt))
        ENDIF
        !
        DO is = 2, nspina
           fac = starting_magnetization(nt) * angular(is-1) * rhoscale
           !$acc parallel loop
           DO ig = 1, ngm
              rhocg(ig,is) = rhocg(ig,is) + fac * &
                            strf(ig,nt) * rhoatg(igtongl(ig))
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !$acc end data
  DEALLOCATE (rhoatg)

END SUBROUTINE atomic_rho_g
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_rho( rhoa, nspina )
  !-----------------------------------------------------------------------
  !! Same as \(\texttt{atomic_rho_g}\), with real-space output charge
  !! \(\text{rhoa}(:,\text{nspina})\).
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE control_flags,        ONLY : gamma_only
  USE lsda_mod,             ONLY : lsda
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_g2r
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspina
  !! number of spin components to be calculated. It may differ from
  !! nspin because in some cases the total charge only is needed, 
  !! even in a LSDA calculation.
  REAL(DP), INTENT(OUT) :: rhoa(dfftp%nnr,nspina)
  !! contains R-space components of the superposition of atomic charges.
  !
  ! ... local variables
  !
  REAL(DP) :: rhoneg
  COMPLEX(DP), allocatable :: rhocg (:,:)
  INTEGER :: ir, is
  !
  ! allocate work space 
  !
  ALLOCATE (rhocg(dfftp%ngm, nspina))
  !
  CALL atomic_rho_g (rhocg, nspina)
  !
  ! bring to real space
  !
  rhoa(:,:) = 0.d0
  CALL rho_g2r ( dfftp, rhocg, rhoa )
  DEALLOCATE (rhocg)
  !
  DO is = 1, nspina
     !
     ! check on negative charge
     !
     rhoneg = 0.0_dp
     DO ir = 1, dfftp%nnr
        rhoneg = rhoneg + MIN (0.0_dp,  DBLE (rhoa (ir,is)) )
     ENDDO
     rhoneg = omega * rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     !
     CALL mp_sum(  rhoneg, intra_bgrp_comm )
     !
     IF ( (is == 1) .OR. lsda ) THEN
        !
        IF ( (rhoneg < -1.0d-4) ) THEN
           IF ( lsda ) THEN 
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
                   &"(component",i1,"):",f12.6)') is, rhoneg
           ELSE
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
          &          f12.6)') rhoneg
           END IF
        END IF
     END IF
     !
     ! it is useless to set negative terms to zero in real space: 
     ! negative charge will re-appear when Fourier-transformed back and forth
     !
  ENDDO

END SUBROUTINE atomic_rho

