!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!---------------------------------------------------------------------
MODULE plot_wan2odd
  !-------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  !
  !
  IMPLICIT NONE 
  !
  PRIVATE
  !
  PUBLIC :: plot_wann
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE plot_wann( list, nrtot, nwann )
    !-------------------------------------------------------------------
    !
    ! ...  This routine generates a XSF file, in a format readable
    ! ...  by XCrySDen, with the plot of the Wannier functions in list
    !
    USE io_global,           ONLY : ionode, stdout
    USE mp,                  ONLY : mp_sum
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE fft_interfaces,      ONLY : invfft
    USE buffers,             ONLY : get_buffer, close_buffer
    USE cell_base,           ONLY : alat, omega, at
    USE constants,           ONLY : BOHR_RADIUS_ANGS
    USE ions_base,           ONLY : atm
    USE noncollin_module,    ONLY : npol
    USE parameters,          ONLY : ntypx
    USE scell_wfc,           ONLY : bcast_psic
    USE fft_supercell,       ONLY : dfftcp, at_cp, nat_cp, tau_cp, ityp_cp, &
                                    ngmcp, npwxcp, iunwann, nwordwann, gamma_only_x
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: list(:)        ! list of WFs to plot
    INTEGER, INTENT(IN) :: nrtot          ! number of R-vectors
    INTEGER, INTENT(IN) :: nwann          ! number of WFs
    !
    CHARACTER(LEN=30) :: filename
    INTEGER :: fileunit=324
    INTEGER :: ir, ibnd, iw
    INTEGER :: i, j, rr
    INTEGER :: nnrg
    REAL(DP) :: alang
    REAL(DP) :: orig(3), dirs(3,3)
    COMPLEX(DP) :: evc(npwxcp*npol,nwann)  ! Wannier function to plot (plane waves)
    COMPLEX(DP) :: psic(dfftcp%nnr)
    COMPLEX(DP), ALLOCATABLE :: psicg(:)
    !
    !
    CALL start_clock( 'plot_wann' )
    !
    nnrg = dfftcp%nnr
    CALL mp_sum( nnrg, intra_bgrp_comm )
    ALLOCATE( psicg(nnrg) )
    !
    IF ( ionode ) THEN
      WRITE(stdout,*) "Plot of Wannier functions:    iw    ir  ibnd"
    ENDIF
    !
    iw = 0
    !
    DO ir = 1, nrtot
      !
      CALL get_buffer( evc, nwordwann, iunwann, ir )
      !
      DO ibnd = 1, nwann
        !
        iw = iw + 1
        !
        DO j = 1, SIZE( list )
          IF ( list(j) == iw ) GOTO 500
        ENDDO
        !
        CYCLE
        !
        !
500     psic(:) = ( 0.D0, 0.D0 )
        !
        IF ( ionode ) THEN
          WRITE(stdout,'(10x, "Wannier function:", 3I6)') iw, ir, ibnd
        ENDIF
        !
        WRITE( filename, 100 ) ir, ibnd
        evc(:,ibnd) = evc(:,ibnd) / SQRT(DBLE(nrtot))
        psic(dfftcp%nl(1:npwxcp)) = evc(1:npwxcp,ibnd)
        IF ( gamma_only_x ) psic(dfftcp%nlm(1:npwxcp)) = CONJG( evc(1:npwxcp,ibnd) )
        CALL invfft( 'Wave', psic, dfftcp )
        CALL bcast_psic( psic, psicg, dfftcp )
        CALL real_wann_max_mod( psicg, nnrg )
        !
        alang = alat * BOHR_RADIUS_ANGS     ! alat in angstrom
        !
        orig(:) = (/ 0.0, 0.0, 0.0 /)       ! origin of the datagrid
        dirs(:,1) = at_cp(:,1) * alang      ! 1st spanning vector datagrid
        dirs(:,2) = at_cp(:,2) * alang      ! 2nd spanning vector datagrid
        dirs(:,3) = at_cp(:,3) * alang      ! 3rd spanning vector datagrid
        !
        !
        IF ( ionode ) THEN
          !
          OPEN( UNIT=fileunit, FILE=trim(filename), STATUS='unknown', FORM='formatted' ) 
          !
          WRITE( fileunit, 201 ) at_cp(:,1)*alang, at_cp(:,2)*alang, at_cp(:,3)*alang
          WRITE( fileunit, 202 ) at_cp(:,1)*alang, at_cp(:,2)*alang, at_cp(:,3)*alang
          WRITE( fileunit, 203 ) nat_cp
          WRITE( fileunit, 204 ) ( atm(ityp_cp(i)), tau_cp(:,i)*alang, i=1,nat_cp )
          WRITE( fileunit, 205 )
          WRITE( fileunit, 206 ) dfftcp%nr1, dfftcp%nr2, dfftcp%nr3, &
                                 orig, dirs(:,1), dirs(:,2), dirs(:,3) 
          WRITE( fileunit, 207 ) ( REAL(psicg(rr)), rr=1,nnrg )
          WRITE( fileunit, 208 ) 
          !
          CLOSE( fileunit )
          !
        ENDIF
        !
      ENDDO
    ENDDO
    !
    !
    CALL close_buffer( iunwann, 'delete' )
    CALL stop_clock( 'plot_wann' )
    !
    !
100 FORMAT( 'WF_R', I4.4, '_B', I4.4, '.xsf' )    ! ex: WF_R0001_B0002.xsf
    !
201 FORMAT( 'CRYSTAL', /,'PRIMVEC', /, 3F12.7, /, 3F12.7, /, 3F12.7 )
202 FORMAT( 'CONVVEC', /, 3F12.7, /, 3F12.7, /, 3F12.7 )
203 FORMAT( 'PRIMCOORD', /, I6, '  1' )
204 FORMAT( A2, 3X, 3F12.7 )
205 FORMAT( /, 'BEGIN_BLOCK_DATAGRID_3D', /, '3D_field', /, 'BEGIN_DATAGRID_3D_WANNIER' )
206 FORMAT( 3I6, /, 3F12.6, /, 3F12.7, /, 3F12.7, /, 3F12.7 )
207 FORMAT( 6E13.5 )
208 FORMAT( 'END_DATAGRID_3D', /, 'END_BLOCK_DATAGRID_3D' )
    !
    !
  END SUBROUTINE plot_wann
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE real_wann_max_mod( psi, nnr )
    !-------------------------------------------------------------------
    !
    ! ...  This routine normalizes the input Wannier function in 
    ! ...  order to be real at the point where it has the max.
    ! ...  modulus
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: psi(:)
    INTEGER, INTENT(IN) :: nnr
    !
    INTEGER :: ir
    COMPLEX(DP) :: wmod
    REAL(DP) :: tmax, tmax_
    !
    !
    tmax = 0.D0
    wmod = ( 0.D0, 1.D0 )
    !
    DO ir = 1, nnr
      !
      tmax_ = DBLE( psi(ir) )**2 + AIMAG( psi(ir) )**2
      !
      IF ( tmax_ > tmax ) THEN
        !
        tmax = tmax_
        wmod = psi(ir)
        !
      ENDIF
      ! 
    ENDDO
    !
    wmod = wmod / SQRT( DBLE(wmod)**2 + AIMAG(wmod)**2 )
    psi(:) = psi(:) / wmod
    !
    !
  END SUBROUTINE real_wann_max_mod
  !
  !
END MODULE plot_wan2odd
