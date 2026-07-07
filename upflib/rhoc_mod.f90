!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE rhoc_mod
  !
  !! Variables and routines for nonlinear core correction
  !! Contains generation of interpolation tables in reciprocal space,
  !! interpolation routines and other utility routines
  !! Code moved to upflib and restructured by Paolo Giannozzi, 2023
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: init_tab_rhc
  PUBLIC :: interp_rhc
  PUBLIC :: interp_tac
  PUBLIC :: interp_drhc
  PUBLIC :: interp_dtac
  PUBLIC :: scale_tab_rhc
  PUBLIC :: deallocate_tab_rhc
  PUBLIC :: deallocate_tab_tac
  PRIVATE :: interp_tab, interp_dtab
  !
  SAVE
  !
  INTEGER :: nqx = 0
  !! size of interpolation table
  REAL(DP), PARAMETER:: dq = 0.01_dp
  !! grid step for interpolation table
  REAL(DP) :: qmax = 0.0_dp 
  !! max q covered by the interpolation table
  REAL(DP), ALLOCATABLE :: tab_rhc(:,:)
  !! interpolation table for atomic pseudo-core charge density
  REAL(DP), ALLOCATABLE :: tab_tac(:,:)
  !! interpolation table for atomic pseudo-core kinetic energy density
  !
CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE init_tab_rhc (qmax_, omega, comm, ierr)
  !----------------------------------------------------------------------
   !! Compute interpolation table for atomic core (pseudo-)charge density
   !! and kinetic-energy density:
   !! tab_rhc(i,nt) = rhoc_nt(q_i), tab_tac(i,nt) = tauc_nt(q_i)
   !! for atom of type nt, on grid q_i
  !
  USE atom,         ONLY : rgrid, msh
  USE uspp_param,   ONLY : upf, nsp
  USE mp,           ONLY : mp_sum
  !
  INTEGER, INTENT(IN)  :: comm
  !! MPI communicator, to split the workload
  INTEGER, INTENT(OUT) :: ierr
  !! return code: ierr = 0 if interpolation table (IT) was allocated
  !!              ierr =-1 if IT had insufficient dimension and was re-allocated
  !!              ierr =-2 if IT was already present and nothing is done
  REAL(dp), INTENT(IN) :: qmax_
  !! Interpolate q up to qmax_ (sqrt(Ry), q^2 is an energy)
  REAL(dp), INTENT(IN) :: omega
  !! Unit-cell volume
  !
  INTEGER :: ndm, startq, lastq, nt, iq, ir
  !! Various indices
  REAL(dp) :: q
  REAL(dp), ALLOCATABLE :: aux (:)
  !! Work space
  !
  ierr = 0
  IF ( .NOT. ANY(upf(1:nsp)%nlcc) .AND. .NOT. ANY(upf(1:nsp)%with_metagga_info) ) RETURN
  !
  IF ( .NOT. ALLOCATED(tab_rhc) .OR. .NOT. ALLOCATED(tab_tac) ) THEN
     !! one of the tables not yet allocated
     IF ( ALLOCATED(tab_rhc) ) CALL deallocate_tab_rhc ( )
     IF ( ALLOCATED(tab_tac) ) CALL deallocate_tab_tac ( )
     qmax = qmax_
  ELSE IF ( qmax_ > qmax ) THEN
     !! tables are allocated but dimension insufficient: re-allocate
     !! (with some margin so that this does not happen too often)
     CALL deallocate_tab_rhc ( )
     CALL deallocate_tab_tac ( )
     qmax = qmax_ + MAX(dq*100,qmax_-qmax)
     ierr =-1
  ELSE
     !! tables already computed: exit
     ierr =-2
     RETURN
  END IF
  nqx = INT( qmax/dq + 4)
  ALLOCATE ( tab_rhc(nqx,nsp) )
  ALLOCATE ( tab_tac(nqx,nsp) )
  !$acc enter data create(tab_rhc)
  !$acc enter data create(tab_tac)
  !
  ndm = MAXVAL( msh(1:nsp) )
  ALLOCATE (aux(ndm))
  !
  CALL divide (comm, nqx, startq, lastq)
  !
  DO nt = 1, nsp
     !
     tab_rhc(:,nt)= 0.d0
     DO iq = startq, lastq
        !
        IF ( upf(nt)%nlcc ) THEN
           q = (iq - 1) * dq
           DO ir = 1, msh(nt)
              IF ( iq > 1 .AND. rgrid(nt)%r(ir) > eps8 ) then
                 !! check above prevents divide by zero
                 aux (ir) = upf(nt)%rho_atc(ir) * rgrid(nt)%r2(ir) * &
                      sin(q*rgrid(nt)%r(ir)) / (q*rgrid(nt)%r(ir))
              ELSE
                 aux (ir) = upf(nt)%rho_atc(ir) * rgrid(nt)%r2(ir)
              ENDIF
           ENDDO
           !
           CALL simpson ( msh(nt), aux, rgrid(nt)%rab, tab_rhc(iq,nt) )
           tab_rhc (iq,nt) = fpi * tab_rhc (iq,nt) / omega
        END IF
     ENDDO
     !
  END DO
  !
  CALL mp_sum ( tab_rhc (:,1:nsp), comm )
  !$acc update device(tab_rhc)
  !
  DO nt = 1, nsp
     !
     tab_tac(:,nt)= 0.d0
     IF ( .NOT. upf(nt)%with_metagga_info ) CYCLE
     !
     DO iq = startq, lastq
        !
        q = (iq - 1) * dq
        DO ir = 1, msh(nt)
           IF ( iq > 1 .AND. rgrid(nt)%r(ir) > eps8 ) then
              !! check above prevents divide by zero
              aux (ir) = upf(nt)%tau_core(ir) * rgrid(nt)%r2(ir) * &
                   sin(q*rgrid(nt)%r(ir)) / (q*rgrid(nt)%r(ir))
           ELSE
              aux (ir) = upf(nt)%tau_core(ir) * rgrid(nt)%r2(ir)
           ENDIF
        ENDDO
        !
        CALL simpson ( msh(nt), aux, rgrid(nt)%rab, tab_tac(iq,nt) )
        tab_tac (iq,nt) = fpi * tab_tac (iq,nt) / omega
        !
     ENDDO
     !
  END DO
  !
  CALL mp_sum ( tab_tac (:,1:nsp), comm )
  !$acc update device(tab_tac)
  !
  DEALLOCATE (aux)
  RETURN
  !
  END SUBROUTINE init_tab_rhc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_rhc( nt, ngl, gl, tpiba2, rhocg )
  !-----------------------------------------------------------------------
  !! Calculates the radial Fourier transform of the core charge.
  !
  INTEGER,  INTENT(IN)  :: nt, ngl
  REAL(DP), INTENT(IN)  :: gl(ngl), tpiba2
  REAL(DP), INTENT(OUT) :: rhocg(ngl)
  !
  CALL interp_tab( tab_rhc, nt, ngl, gl, tpiba2, rhocg )
  !
  END SUBROUTINE interp_rhc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_tac( nt, ngl, gl, tpiba2, tacg )
  !-----------------------------------------------------------------------
  !! Calculates the radial Fourier transform of the core kinetic energy density.
  !
  INTEGER,  INTENT(IN)  :: nt, ngl
  REAL(DP), INTENT(IN)  :: gl(ngl), tpiba2
  REAL(DP), INTENT(OUT) :: tacg(ngl)
  !
  CALL interp_tab( tab_tac, nt, ngl, gl, tpiba2, tacg )
  !
  END SUBROUTINE interp_tac
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_drhc( nt, ngl, gl, tpiba2, drhocg )
  !-----------------------------------------------------------------------
  !! Calculates the Fourier transform of d Rho_c / dG.
  !
  INTEGER,  INTENT(IN)  :: nt, ngl
  REAL(DP), INTENT(IN)  :: gl(ngl), tpiba2
  REAL(DP), INTENT(OUT) :: drhocg(ngl)
  !
  CALL interp_dtab( tab_rhc, nt, ngl, gl, tpiba2, drhocg )
  !
  END SUBROUTINE interp_drhc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_dtac( nt, ngl, gl, tpiba2, dtacg )
  !-----------------------------------------------------------------------
  !! Calculates the Fourier transform of d tau_c / dG.
  !
  INTEGER,  INTENT(IN)  :: nt, ngl
  REAL(DP), INTENT(IN)  :: gl(ngl), tpiba2
  REAL(DP), INTENT(OUT) :: dtacg(ngl)
  !
  CALL interp_dtab( tab_tac, nt, ngl, gl, tpiba2, dtacg )
  !
  END SUBROUTINE interp_dtac
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_tab( tab, nt, ngl, gl, tpiba2, fg )
  !-----------------------------------------------------------------------
  !! Cubic interpolation of a tabulated radial Fourier transform.
  !
  REAL(DP), INTENT(IN)  :: tab(:,:)
  INTEGER,  INTENT(IN)  :: nt, ngl
  REAL(DP), INTENT(IN)  :: gl(ngl), tpiba2
  REAL(DP), INTENT(OUT) :: fg(ngl)
  !
  REAL(DP) :: gx, px, ux, vx, wx
  INTEGER  :: igl, i0, i1, i2, i3
  !
  !$acc data copyin(gl, tab) copyout(fg)
  !$acc parallel loop
  DO igl = 1, ngl
     gx = SQRT( gl(igl) * tpiba2 )
     px = gx / dq - INT(gx/dq)
     ux = 1.d0 - px
     vx = 2.d0 - px
     wx = 3.d0 - px
     i0 = INT(gx/dq) + 1
     i1 = i0 + 1
     i2 = i0 + 2
     i3 = i0 + 3
     fg(igl) = tab(i0,nt) * ux * vx * wx / 6.d0 + &
               tab(i1,nt) * px * vx * wx / 2.d0 - &
               tab(i2,nt) * px * ux * wx / 2.d0 + &
               tab(i3,nt) * px * ux * vx / 6.d0
  ENDDO
  !$acc end data
  !
  END SUBROUTINE interp_tab
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_dtab( tab, nt, ngl, gl, tpiba2, dfg )
  !-----------------------------------------------------------------------
  !! Cubic interpolation of d(tabulated FT) / dG.
  !
  REAL(DP), INTENT(IN)  :: tab(:,:)
  INTEGER,  INTENT(IN)  :: nt, ngl
  REAL(DP), INTENT(IN)  :: gl(ngl), tpiba2
  REAL(DP), INTENT(OUT) :: dfg(ngl)
  !
  REAL(DP) :: gx, px, ux, vx, wx
  INTEGER  :: igl, i0, i1, i2, i3
  !
  !$acc data copyin(gl, tab) copyout(dfg)
  !$acc parallel loop
  DO igl = 1, ngl
     gx = SQRT( gl(igl) * tpiba2 )
     px = gx / dq - INT(gx/dq)
     ux = 1.d0 - px
     vx = 2.d0 - px
     wx = 3.d0 - px
     i0 = INT(gx/dq) + 1
     i1 = i0 + 1
     i2 = i0 + 2
     i3 = i0 + 3
     dfg(igl) = (- tab(i0,nt) * (ux*vx + vx*wx + ux*wx) / 6.0_dp &
                 + tab(i1,nt) * (wx*vx - px*wx - px*vx) / 2.0_dp &
                 - tab(i2,nt) * (wx*ux - px*wx - px*ux) / 2.0_dp &
                 + tab(i3,nt) * (ux*vx - px*ux - px*vx) / 6.0_dp ) / dq
  ENDDO
  !$acc end data
  !
  END SUBROUTINE interp_dtab
  !
  subroutine scale_tab_rhc( vol_ratio_m1 )
     ! vol_ratio_m1 = omega_old / omega
     real(DP), intent(in) :: vol_ratio_m1
     !
     if ( allocated(tab_rhc) ) then
        tab_rhc(:,:)  = tab_rhc(:,:) * vol_ratio_m1
        !$acc update device (tab_rhc)
     end if
     if ( allocated(tab_tac) ) then
        tab_tac(:,:)  = tab_tac(:,:) * vol_ratio_m1
        !$acc update device (tab_tac)
     end if
     !
  end subroutine scale_tab_rhc
  !
  subroutine deallocate_tab_rhc(  )
     !
     if ( allocated(tab_rhc) ) then
        !$acc exit data delete(tab_rhc)
        deallocate (tab_rhc)
     end if
     !
  end subroutine deallocate_tab_rhc
  !
  subroutine deallocate_tab_tac(  )
     !
     if ( allocated(tab_tac) ) then
        !$acc exit data delete(tab_tac)
        deallocate (tab_tac)
     end if
     !
  end subroutine deallocate_tab_tac
  !
END MODULE rhoc_mod

