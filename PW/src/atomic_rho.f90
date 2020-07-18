!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
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
  USE atom,                 ONLY : rgrid, msh
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : tpiba, omega
  USE gvect,                ONLY : ngm, ngl, gstart, gl, igtongl
  USE lsda_mod,             ONLY : starting_magnetization
  USE vlocal,               ONLY : starting_charge, strf
  USE noncollin_module,     ONLY : angle1, angle2
  USE uspp_param,           ONLY : upf
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
  REAL(DP) :: rhoneg, rhoima, rhoscale, gx
  REAL(DP), ALLOCATABLE :: rhocgnt(:), aux(:)
  REAL(DP) :: angular(nspina)
  INTEGER :: ir, is, ig, igl, nt, ndm
  !
  ! allocate work space 
  !
  ndm = MAXVAL ( msh(1:ntyp) )
  ALLOCATE (rhocgnt( ngl))
  !
!$omp parallel private(aux, gx, rhoscale, angular)
  !
  call threaded_nowait_memset(rhocg, 0.0_dp, ngm*nspina*2)
  !
  ALLOCATE (aux(ndm))
  !
  DO nt = 1, ntyp
     !
     ! Here we compute the G=0 term
     !
!$omp master
     IF (gstart == 2) then
        DO ir = 1, msh (nt)
           aux (ir) = upf(nt)%rho_at (ir)
        ENDDO
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (1) )
     ENDIF
!$omp end master
     !
     ! Here we compute the G<>0 term
     !
!$omp do
     DO igl = gstart, ngl
        gx = sqrt (gl (igl) ) * tpiba
        DO ir = 1, msh (nt)
           IF (rgrid(nt)%r(ir) < eps8) then
              aux(ir) = upf(nt)%rho_at(ir)
           ELSE
              aux(ir) = upf(nt)%rho_at(ir) * &
                        sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
           ENDIF
        ENDDO
        CALL simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (igl) )
     ENDDO
!$omp end do
     !
     ! we compute the 3D atomic charge in reciprocal space
     !
     IF (upf(nt)%zp > eps8) THEN
        rhoscale = MAX(0.0_dp, upf(nt)%zp - starting_charge(nt)) / upf(nt)%zp
     ELSE
        rhoscale = 1.0_dp
     ENDIF
     !
     !
!$omp do
     DO ig = 1, ngm
        rhocg(ig,1) = rhocg(ig,1) + &
                strf(ig,nt) * rhoscale * rhocgnt(igtongl(ig)) / omega
     ENDDO
!$omp end do nowait
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
!$omp do
           DO ig = 1, ngm
              rhocg(ig,is) = rhocg(ig,is) + &
                            starting_magnetization(nt) * angular(is-1) * &
                            strf(ig,nt) * rhoscale * rhocgnt(igtongl(ig)) / omega
           ENDDO
!$omp end do nowait
        ENDDO
        !
     ENDIF
     ! must complete the computation of rhocg before updating rhocgnt
     ! for the next type
!$omp barrier
  ENDDO

  DEALLOCATE (aux)
!$omp end parallel

  DEALLOCATE (rhocgnt)

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
  USE cell_base,            ONLY : tpiba, omega
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

