!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE rhoat_mod
  !
  !! Variables and routines for atomic charge density
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
  PUBLIC :: init_tab_rhoat
  PUBLIC :: interp_rhoat
  PUBLIC :: scale_tab_rhoat
  !
  SAVE
  !
  INTEGER :: nqx = 0
  !! size of interpolation table
  REAL(DP), PARAMETER:: dq = 0.01_dp
  !! grid step for interpolation table
  REAL(DP) :: qmax = 0.0_dp 
  !! max q covered by the interpolation table
  REAL(DP), ALLOCATABLE :: tab_rhoat(:,:)
  !! interpolation table for atomic pseudo-valence charge density
  !
CONTAINS
  !
  
  !----------------------------------------------------------------------
  SUBROUTINE init_tab_rhoat (qmax_, omega, comm, ierr)
  !----------------------------------------------------------------------
  !
  !! Compute interpolation table for atomic charge density
  !! tab_rhoat(i,nt) = \rho_nt(q_i) for atom of type nt, on grid q_i
  !! Output in tab_rhoat
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
  IF ( .NOT. ALLOCATED(tab_rhoat) ) THEN
     !! table not yet allocated
     qmax = qmax_
     ierr = 0
  ELSE IF ( qmax_ > qmax ) THEN
     !! table ìs allocated but dimension insufficient: re-allocate
     !! (with some margin so that this does not happen too often)
     CALL deallocate_tab_rhoat ( )
     qmax = qmax_ + MAX(dq*100,qmax_-qmax)
     ierr =-1
  ELSE
     !! table already computed: exit
     ierr =-2
     RETURN
  END IF
  nqx = INT( qmax/dq + 4)
  ALLOCATE ( tab_rhoat(nqx,nsp) )
  !$acc enter data create(tab_rhoat)
  !
  ndm = MAXVAL( msh(1:nsp) )
  ALLOCATE (aux(ndm))
  !
  CALL divide (comm, nqx, startq, lastq)
  !
  DO nt = 1, nsp
     !
     tab_rhoat(:,nt)= 0.d0
     !
     DO iq = startq, lastq
        !
        q = (iq - 1) * dq
        DO ir = 1, msh(nt)
           IF ( iq > 1 .AND. rgrid(nt)%r(ir) > eps8 ) then
              !! check above prevents divide by zero
              aux (ir) = upf(nt)%rho_at(ir) * & 
                        sin(q*rgrid(nt)%r(ir)) / (q*rgrid(nt)%r(ir))
           ELSE
              aux (ir) = upf(nt)%rho_at(ir)
           ENDIF 
        ENDDO
        !
        CALL simpson ( msh(nt), aux, rgrid(nt)%rab, tab_rhoat(iq,nt) )
        tab_rhoat (iq,nt) = tab_rhoat (iq,nt) / omega 
        !
     ENDDO
     !
  END DO
  !
  CALL mp_sum ( tab_rhoat (:,1:nsp), comm )
  !$acc update device (tab_rhoat)
  !
  DEALLOCATE (aux)
  !
END SUBROUTINE init_tab_rhoat
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_rhoat( nt, ngl, gl, tpiba2, rhoag )
  !-----------------------------------------------------------------------
  !! Calculates the radial Fourier transform of the core charge.
  !
  INTEGER :: nt
  !! input: atomic type
  INTEGER :: ngl
  !! input: the number of g shell
  REAL(DP) :: gl(ngl)
  !! input: the number of G shells
  REAL(DP) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP) :: rhoag(ngl)
  !! output: the Fourier transform of the atomic charge
  !
  ! ... local variables
  !
  REAL(DP) :: gx, px, ux, vx, wx
  ! the modulus of g for a given shell
  ! variables used for interpolation
  INTEGER :: igl, i0, i1, i2, i3
  ! counters
  !
  !$acc data present_or_copyin(gl) present_or_copyout(rhoag) present(tab_rhoat)
  !$acc parallel loop
  DO igl = 1, ngl
     gx = SQRT(gl(igl) * tpiba2)
     px = gx / dq - int (gx/dq)
     ux = 1.d0 - px
     vx = 2.d0 - px
     wx = 3.d0 - px
     i0 = INT(gx/dq) + 1
     i1 = i0 + 1
     i2 = i0 + 2
     i3 = i0 + 3
     rhoag (igl) = tab_rhoat(i0, nt) * ux * vx * wx / 6.d0 + &
                   tab_rhoat(i1, nt) * px * vx * wx / 2.d0 - &
                   tab_rhoat(i2, nt) * px * ux * wx / 2.d0 + &
                   tab_rhoat(i3, nt) * px * ux * vx / 6.d0

  ENDDO
  !$acc end data
  !
END SUBROUTINE interp_rhoat
  !
  subroutine scale_tab_rhoat( vol_ratio_m1 )
     ! vol_ratio_m1 = omega_old / omega
     real(DP), intent(in) :: vol_ratio_m1
     !
     if ( allocated(tab_rhoat) ) then
         tab_rhoat(:,:)  = tab_rhoat(:,:) * vol_ratio_m1
         !$acc update device (tab_rhoat)
     end if
     !
  end subroutine scale_tab_rhoat
  !
  subroutine deallocate_tab_rhoat(  )
     !
     if ( allocated(tab_rhoat) ) then
         !$acc exit data delete(tab_rhoat)
         deallocate (tab_rhoat)
     end if
     !
  end subroutine deallocate_tab_rhoat
  !
END MODULE rhoat_mod

