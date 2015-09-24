!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE o_rinitcgg( npwx, npw, nstart, numwp, psi, o_evc, e, numv, v_states,hdiag,ptype,fcw_number,fcw_state,fcw_mat,ethr)
  !----------------------------------------------------------------------------
  !
  ! ... Operator O  diagonalization in the subspace spanned
  ! ... by nstart states psi (atomic or random wavefunctions).
  ! ... Produces on output numwp eigenvectors (numwp <= nstart) in o_evc.
  ! ... Minimal memory use - o_evc and psi may overlap
  ! ... Calls o_1psi to calculate O|psi>
  !
  USE kinds, ONLY : DP
  USE gvect, ONLY : gstart
  USE io_global, ONLY : stdout
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  USE fft_base,             ONLY : dffts
  !
  IMPLICIT NONE
  !
  INTEGER :: npw, npwx, nstart, numwp
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  COMPLEX (DP) :: psi(npwx,nstart), o_evc(npwx,numwp)
    ! input and output eigenvectors (may overlap) 
  REAL(DP) :: e(numwp)
  
  INTEGER, INTENT(in) ::  numv!number of valence states
  REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv) !valence states in real space
  REAL(kind=DP), INTENT(in) :: hdiag(npw)!inverse of estimation of diagonal part of hamiltonian  
  INTEGER, INTENT(in) :: ptype!type of approximation for O operator
  INTEGER, INTENT(in) :: fcw_number!number of "fake conduction" states for O matrix method
  COMPLEX(kind=DP) :: fcw_state(npw,fcw_number)! "fake conduction" states for O matrix method
  REAL(kind=DP) :: fcw_mat(fcw_number,fcw_number)! "fake conduction" matrix
  REAL(kind=DP), INTENT(in) :: ethr!threshold for o_1psi_gamma 


    ! eigenvalues
  !
  !... local variables
  !
  INTEGER                        :: m, i, j, npw2, npwx2
  REAL (DP)                 :: rtmp(2)
  COMPLEX (DP), ALLOCATABLE :: aux(:,:)
  COMPLEX (DP), ALLOCATABLE :: ctmp(:)
  REAL (DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:)
  REAL (DP),    ALLOCATABLE :: en(:)
  !
  !
  CALL start_clock( 'wfcrot1' )
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx
  !
  ALLOCATE( aux( npwx, 2 ) )
  ALLOCATE( ctmp( numwp ) )
  ALLOCATE( hr( nstart, nstart, 2 ) )
  ALLOCATE( sr( nstart, nstart ) )
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix
  !
  DO m = 1, nstart
     !
     !CALL h_1psi( npwx, npw, psi(1,m), aux(1,1), aux(1,2) )
   
     write(stdout,*) 'Call o_1psi_gamma',m,nstart
     FLUSH(stdout)
     call o_1psi_gamma( numv, v_states, psi(1,m), aux(1,1),.false.,hdiag,ptype,fcw_number,fcw_state,fcw_mat,ethr) 
     write(stdout,*) 'Done'
     FLUSH(stdout)

     !call o_1psi_gamma_real( numv, v_states, psi(1,m), aux(1,1))
     aux(:,2)=psi(:,m)

     !
     CALL DGEMV( 'T', npw2, 2, 2.D0, aux, npwx2, psi(1,m), 1, 0.D0, rtmp, 1 )
     !
     IF ( gstart == 2 ) rtmp(:) = rtmp(:) - psi(1,m) * aux(1,:)
     !
     hr(m,m,1) = rtmp(1)
     sr(m,m)   = rtmp(2)
     !
     DO j = m + 1, nstart
        !
        CALL DGEMV( 'T', npw2, 2, 2.D0, aux, npwx2, psi(1,j), 1, 0.D0, rtmp, 1 )
        !
        IF ( gstart == 2 ) rtmp(:) = rtmp(:) - psi(1,j) * aux(1,:)
        !
        hr(j,m,1) = rtmp(1)
        sr(j,m)   = rtmp(2)
        !
        hr(m,j,1) = rtmp(1)
        sr(m,j)   = rtmp(2)
        !
     END DO
     !
  END DO
  !
  !CALL reduce(  nstart * nstart, hr(1,1,1) )
  call mp_sum(hr(:,:,1),world_comm)
  !CALL reduce(  nstart * nstart, sr(1,1) )
  CALL mp_sum(sr(:,:),world_comm)
  !
  ! ... diagonalize
  !
  write(stdout,*) 'Call rdiaghg'
  FLUSH(stdout)

  CALL rdiaghg( nstart, numwp, hr, sr, nstart, en, hr(1,1,2) )
  write(stdout,*) 'Done'
  FLUSH(stdout)

 !
  e(1:numwp) = en(1:numwp)
  !
  ! ... update the basis set
  !
  DO i = 1, npw
     !
     DO m = 1, numwp
        !
        ctmp(m) = SUM( hr(:,m,2) * psi(i,:) )
        !
     END DO
     !
     o_evc(i,1:numwp) = ctmp(1:numwp)
     !
  END DO
  !
  DEALLOCATE( en )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( ctmp )
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'wfcrot1' )
  !
  RETURN
  !
END SUBROUTINE o_rinitcgg
