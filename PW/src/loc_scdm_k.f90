! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE loc_scdm_k
  !--------------------------------------
  !! Variables and subroutines for localizing molecular orbitals based
  !! on a modified SCDM approach. Original SCDM method:
  !
  !! A. Damle, L. Lin, L. Ying: J. Comput. Phys., 334 (2017) 1-5 
  !
  !! Ivan Carnimeo: SCDM has been implemented with a parallel prescreening algorithm 
  !! in order to save CPU time and memory for large grids.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE exx,                  ONLY : dfftt, exxbuff, exxmat
  USE exx_base,             ONLY : nkqs
  !
  IMPLICIT NONE
  !
  SAVE
  !
  REAL(DP), PARAMETER :: Zero=0.0_DP, One=1.0_DP, Two=2.0_DP, Three=2.0_DP
  INTEGER :: n_scdm = 1
  !
 CONTAINS
  !
!------------------------------------------------------------------------
SUBROUTINE localize_orbitals_k( )
  !----------------------------------------------------------------------
  !! Driver for SCDM_k orbital localization:
  !
  !! * \(\texttt{SCDM_PGG_k}\): performs SCDM localization using the 
  !!                            parallel prescreening algorithm;
  !! * \(\texttt{measure_localization_k}\): computes the absolute overlap
  !!                     integrals required by \(\texttt{vexx_loc_k in}\)
  !!                     exx.f90 and puts them in \(\texttt{exxmat}\).
  !
  !! NOTE: \(\texttt{localize_orbitals_k}\) does not have iscdm/nscdm because 
  !! alignment is not implemented for K-Points.
  !
  USE noncollin_module,     ONLY : npol
  USE exx,                  ONLY : x_occupation
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : nbnd
  USE mp_bands,             ONLY : me_bgrp
  USE control_flags,        ONLY : lmd
  !   
  IMPLICIT NONE
  !
  INTEGER :: NGrid, ikq, jk, k, NBands, nnk
  CHARACTER(LEN=1) :: HowTo
  REAL(DP) :: Tot1, Tot2, Ave1, tot, ave 
  !
  IF ( n_scdm /= 1 ) CALL errore( 'localize_orbitals_k','nscdm for K-points NYI.', 1 )
  IF ( lmd ) CALL errore( 'localize_orbitals_k', 'localization with K-points not tested.', 1 )
  !
  NGrid = dfftt%nnr * npol
  HowTo = 'G'  ! How to compute the absolute overlap integrals (only G for k)
  !
  exxmat = One ! initialized to one because the absolute overlaps between
  !              occupied and virtual orbitals are assumed always large 
  !              and all ov integrals will be computed by vexx_loc_k.
  !
  NBands = INT(SUM(x_occupation(:,1))) 
  !
  WRITE (stdout,*) ' '
  WRITE (stdout,*) 'NBands = ', NBands, ' nks = ', nks, ' nkqs = ', nkqs
  WRITE (stdout,'(5X,A)') 'Canonical Orbitals '
  !
  Tot1 = Zero 
  Ave1 = Zero 
  Tot2 = Zero 
  nnk = 0
  !
  DO ikq = 1, nkqs
     !
     CALL measure_localization_k( NBands, ikq, tot, ave )
     !
     Tot1 = Tot1 + tot
     Ave1 = Ave1 + ave
     !
     DO jk = 1, nks
       !
       nnk = nnk + 1
       !
       CALL AbsOvG_k( NBands, ikq, jk, tot, ave )
       !
       Tot2 = Tot2 + tot
       !
     ENDDO
     !
  ENDDO
  !
  Ave1 = Ave1/DBLE(nkqs)
  WRITE (stdout,'(7X,A,f24.6)') 'Total AbsOv          =', Tot2 
  WRITE (stdout,'(7X,A,f24.6)') 'Aver. AbsOv          =', Tot2/DBLE(nnk)
  WRITE (stdout,'(7X,A,f24.6)') 'Total Spread [A**2]  =', Tot1 
  WRITE (stdout,'(7X,A,f24.6)') 'Aver. Spread [A**2]  =', Ave1 
  !
  ! IF (me_bgrp == 0) THEN
  !   CALL AbsOv_histogram_k( NBands, 'hist1.dat' )
  ! ENDIF 
  !
  WRITE (stdout,'(5X,A)') 'SCDM-PGG_k localization'
  DO ikq = 1, nkqs
     CALL SCDM_PGG_k( NGrid, NBands, ikq )
  ENDDO 
  WRITE (stdout,'(7X,A)') 'SCDM-PGG_k done ' 
  !
  WRITE (stdout,'(5X,A)') 'Localized Orbitals '
  !
  Tot1 = Zero 
  Ave1 = Zero 
  Tot2 = Zero 
  nnk = 0
  !
  DO ikq = 1, nkqs
    !
    CALL measure_localization_k( NBands, ikq, tot, ave )
    !
    Tot1 = Tot1 + tot
    Ave1 = Ave1 + ave
    !
    DO jk = 1, nks
       !
       nnk = nnk + 1
       !
       CALL AbsOvG_k( NBands, ikq, jk, tot, ave )
       !
       Tot2 = Tot2 + tot
       !
    ENDDO 
    !
  ENDDO 
  !
  Ave1 = Ave1 / DBLE(nkqs)
  WRITE (stdout,'(7X,A,f24.6)') 'Total AbsOv         =', Tot2
  WRITE (stdout,'(7X,A,f24.6)') 'Aver. AbsOv         =', Tot2/DBLE(nnk)
  WRITE (stdout,'(7X,A,f24.6)') 'Total Spread [A**2] =', Tot1
  WRITE (stdout,'(7X,A,f24.6)') 'Aver. Spread [A**2] =', Ave1
  !
  ! IF (me_bgrp == 0) THEN
  !   CALL AbsOv_histogram_k( NBands, 'hist2.dat' )
  ! ENDIF
  !
END SUBROUTINE localize_orbitals_k
!
!
!-------------------------------------------------------------------------------
SUBROUTINE measure_localization_k( NBands, IK, TotSpread, AveSpread )
  !-----------------------------------------------------------------------------
  !! Computes the absolute overlap integrals required by \(\texttt{vexx_loc_k}\)
  !! in exx.f90 and puts them in \(\texttt{exxmat}\).
  !
  USE noncollin_module,    ONLY : npol
  USE cell_base,           ONLY : alat, omega, at, bg
  USE exx,                 ONLY : compute_density_k
  USE constants,           ONLY : bohr_radius_angs 
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: TotSpread
  !! Total Spread increment per IK
  REAL(DP), INTENT(OUT) :: AveSpread
  !! Average Spread increment per IK
  INTEGER :: NBands
  !! number of bands (with auxiliary functions)
  INTEGER :: IK
  !! index on k+q points
  !
  ! ... local variables
  !
  INTEGER :: ibnd
  REAL(DP) :: RJunk(4), SpreadPBC(3)
  !
  CALL start_clock( 'measure' )
  !
  TotSpread = Zero
  AveSpread = Zero
  !
  DO ibnd = 1, NBands
    !
    CALL compute_density_k( .FALSE., .FALSE., RJunk(1:3), SpreadPBC, RJunk(4), &
                            exxbuff(1,ibnd,IK), exxbuff( 1,ibnd,IK), &
                            dfftt%nnr*npol, ibnd, ibnd )
    TotSpread = TotSpread + SpreadPBC(1) + SpreadPBC(2) + SpreadPBC(3)
    !
  ENDDO
  !
  TotSpread = TotSpread * bohr_radius_angs**2
  AveSpread = TotSpread / DBLE(NBands) 
  !
  ! WRITE (stdout,'(7X,A,I3)')  'IK = ', IK
  ! WRITE (stdout,'(7X,A,f12.6,I3)') 'Total Spread [A**2]   =', TotSpread 
  ! WRITE (stdout,'(7X,A,f12.6,I3)') 'Aver. Spread [A**2]   =', AveSpread 
  !
  CALL stop_clock( 'measure' )
  !
END SUBROUTINE measure_localization_k
!
!
!---------------------------------------------------------------------
SUBROUTINE AbsOvG_k( NBands, IKQ, JK, loc_diag, loc_off )
  !--------------------------`----------------------------------------
  !! Computes the Absolute Overlap in G-space (cutoff might not be 
  !! accurate for the moduli of the wavefunctions).
  !
  USE noncollin_module,    ONLY : npol
  USE fft_interfaces,      ONLY : fwfft
  USE wvfct,               ONLY : npwx
  USE exx_band,            ONLY : igk_exx
  USE klist,               ONLY : nks, ngk
  !
  IMPLICIT NONE
  !
  INTEGER :: NBands
  !! number of bands (with auxiliary functions)
  INTEGER :: IKQ
  !! index on k+q points
  INTEGER :: JK
  !! index on k-point
  REAL(DP), INTENT(OUT) :: loc_diag
  !! Total Charge
  REAL(DP), INTENT(OUT) :: loc_off
  !! Total Abs. Overlap
  !
  ! ... local variables
  !
  INTEGER :: ibnd, jbnd, ig, kk, npw
  REAL(DP) :: tmp
  COMPLEX(DP), ALLOCATABLE :: buffer(:), GorbtI(:,:), GorbtJ(:,:)
  INTEGER, EXTERNAL  :: global_kpoint_index
  COMPLEX(DP), ALLOCATABLE :: Mat(:,:)
  !
  CALL start_clock( 'measure' )
  !
  ! WRITE (stdout, '(5X,A)') ' ' 
  ! WRITE (stdout, '(5X,A)') 'Absolute Overlap calculated in G-space (K)'
  !
  kk = global_kpoint_index ( nks, JK )
  !
  ALLOCATE( Mat(NBands,NBands) )
  ALLOCATE( buffer(dfftt%nnr*npol), GorbtI(npwx,NBands), GorbtJ(npwx,NBands) )
  !
  GorbtJ = (Zero,Zero) 
  GorbtI = (Zero,Zero) 
  npw = ngk(JK)
  !
  DO jbnd = 1, NBands 
     !
     buffer(:) = ABS(exxbuff(:,jbnd,JK)) 
     ! buffer(:) = exxbuff(:,jbnd,JK)
     !
     CALL fwfft( 'Wave' , buffer, dfftt )
     !
     DO ig = 1, npw
       GorbtJ(ig,jbnd) = buffer(dfftt%nl(igk_exx(ig,kk))) 
     ENDDO
     !
     buffer(:) = ABS(exxbuff(:,jbnd,IKQ)) 
     ! buffer(:) = exxbuff(:,jbnd,IKQ)
     !
     CALL fwfft( 'Wave' , buffer, dfftt )
     !
     DO ig = 1, npw
       GorbtI(ig,jbnd) = buffer(dfftt%nl(igk_exx(ig,kk)))  
     ENDDO
     !
  ENDDO
  !
  ! IF (IKQ == JK) THEN
  !   CALL matcalc_k( 'AbsOv-',.FALSE., 2, 0, npwx, NBands, NBands, GorbtI, GorbtJ, Mat, tmp )
  ! ELSE
  CALL matcalc_k( 'AbsOv-',.FALSE., 0, 0, npwx, NBands, NBands, GorbtI, GorbtJ, Mat, tmp )
  ! END IF 
  !
  loc_diag = Zero 
  loc_off  = Zero
  !
  DO ibnd = 1, NBands
     !
     loc_diag = loc_diag + DBLE(Mat(ibnd,ibnd))
     !
     DO jbnd = 1, ibnd-1 
        !
        loc_off = loc_off + DBLE(Mat(ibnd,jbnd)) + DBLE(Mat(jbnd,ibnd))
        exxmat(ibnd, IKQ, jbnd, JK) = DBLE(Mat(ibnd,jbnd))
        exxmat(jbnd, IKQ, ibnd, JK) = DBLE(Mat(jbnd,ibnd))
        !
     ENDDO
     !
  ENDDO
  !
  WRITE (stdout,'(7X,5(A,I3),2(A,f12.6))')  'IKQ = ', IKQ, & 
        '  JK = ', JK, '  kk = ', kk, ' NBands = ', NBands, ' size = ', SIZE(exxbuff,3), &
        '  Total Charge =', loc_diag, '  Total Abs. Overlap =', loc_off 
  !
  DEALLOCATE( buffer, GorbtI, GorbtJ )
  DEALLOCATE( Mat ) 
  !
  CALL stop_clock( 'measure' )
  !
END SUBROUTINE AbsOvG_k
!
!
!------------------------------------------------------------------------------------
SUBROUTINE AbsOv_histogram_k( n, filename )
  !----------------------------------------------------------------------------------
  !! Prints the distribution of absolute overlap integrals (as an histogram) in the
  !! range [0,1] with NHist intervals.
  !
  !! WARNING: HUGE amount of integrals, use only for VERY small systems with a few
  !! K-Points.
  !
  USE klist,   ONLY: nks
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: filename
  !! file with histogram data
  INTEGER, INTENT(IN) :: n
  !! number of bands
  !
  ! ... local variables
  !
  INTEGER :: i,j,k, iik, jjk, NHist
  INTEGER,  ALLOCATABLE :: histogram(:)
  REAL(DP), ALLOCATABLE :: XHist(:)
  REAL(DP) :: xstart, xstep, integral
  INTEGER :: io_histogram
  INTEGER, EXTERNAL :: find_free_unit
  !
  NHist = 1000
  xstep  = One/FLOAT(NHist)
  xstart = One/FLOAT(NHist)/Two
  !
  WRITE (stdout,'(A,I7,2(A,f12.6))') 'NHist  = ', NHist , ' xstep = ', xstep, ' xstart = ', xstart
  !
  ALLOCATE( histogram(NHist), XHist(NHist) )
  !
  XHist = Zero
  histogram = 0
  !
  DO k = 1, NHist 
    XHist(k) = xstart + FLOAT(k-1)*xstep 
  ENDDO
  !
  DO i = 1, n 
     DO iik = 1, nkqs 
        !
        DO j = 1, n
           DO jjk = 1, nks
              !
              integral = exxmat(i, iik, j, jjk)
              !
              IF (integral < Zero) THEN 
                 CALL errore( 'AbsOv_histogram_k','Abs. Ov. < 0 found.', 1 )
              ELSE
                 DO k = 1, NHist 
                    IF (integral >= (XHist(k)-xstart) .AND. integral < (XHist(k)+xstart)) &
                                                          histogram(k) = histogram(k) + 1
                 ENDDO 
              ENDIF
              !
           ENDDO
        ENDDO
        !
     ENDDO 
  ENDDO
  !
  io_histogram = find_free_unit()
  OPEN( io_histogram, file=filename, status='unknown' )
  !
  DO k = 1, NHist
    WRITE (io_histogram, '(f12.6,2I10)') XHist(k), histogram(k)
  ENDDO 
  !
  CLOSE(io_histogram)
  !
  DEALLOCATE( histogram, XHist )
  !
END SUBROUTINE AbsOv_histogram_k
!
!
!-----------------------------------------------------------------------------
SUBROUTINE SCDM_PGG_k( NGrid, NBands, IKQ )
  !---------------------------------------------------------------------------
  !! Density matrix localization (I/O in psi). K-Points version.
  !
  USE mp_bands,      ONLY: nproc_bgrp
  USE loc_scdm,      ONLY: scdm_thresholds, scdm_points, scdm_prescreening
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NGrid
  !! Number of grid points
  INTEGER, INTENT(IN) :: NBands
  !! Number of bands
  INTEGER, INTENT(IN) :: IKQ
  !! Index on k+q points
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: QRbuff(:,:), mat(:,:)
  INTEGER,  ALLOCATABLE :: pivot(:), list(:), cpu_npt(:)
  REAL(DP), ALLOCATABLE :: den(:), grad_den(:,:)
  INTEGER  :: nptot
  REAL(DP) :: ThrDen, ThrGrd
  !
  CALL start_clock( 'localization' )
  !
  ! WRITE (stdout,'(5X,A)') 'SCDM localization with prescreening'
  !
  ALLOCATE( den(dfftt%nnr), grad_den(3, dfftt%nnr) )
  !
  CALL scdm_thresholds( den, grad_den, ThrDen, ThrGrd )
  !
  ALLOCATE( cpu_npt(0:nproc_bgrp-1) )
  !
  CALL scdm_points( den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot )
  !
  ALLOCATE( list(nptot), pivot(nptot) )
  !
  CALL scdm_prescreening_k( NGrid, NBands, exxbuff(1,1,IKQ), den, grad_den, &
                            ThrDen, ThrGrd, cpu_npt, nptot, list, pivot )
  DEALLOCATE( den, grad_den ) 
  !
  ! Psi(pivot(1:NBands),:) in mat
  !
  ALLOCATE( mat(NBands,NBands) )
  !
  mat = (Zero, Zero)
  !
  CALL scdm_fill_k( .TRUE., nptot, NGrid, NBands, cpu_npt, pivot, list, &
                    exxbuff(1,1,IKQ), Mat )
  !
  ! Pc = Psi * Psi(pivot(1:NBands),:)' in QRbuff
  !
  ALLOCATE( QRbuff(NGrid, NBands) )
  !
  QRbuff = (Zero, Zero)
  !
  CALL ZGEMM( 'N' , 'N' , NGrid, NBands, NBands, (One,Zero), exxbuff(1,1,IKQ), &
              NGrid, mat, NBands, (Zero,Zero), QRbuff, NGrid )
  !
  ! Orthonormalization
  ! Pc(pivot(1:NBands),:) in mat 
  !
  mat = (Zero, Zero)
  !
  CALL scdm_fill_k( .FALSE., nptot, NGrid, NBands, cpu_npt, pivot, list, QRBuff, mat )
  !
  DEALLOCATE( cpu_npt )
  !
  ! Cholesky(psi)^(-1) in mat 
  CALL invchol_k(NBands,mat)
  !
  ! Phi = Pc * Chol^(-1) = QRbuff * mat
  CALL ZGEMM( 'N' , 'T' , NGrid, NBands, NBands, (One,Zero), QRbuff, NGrid, mat, &
              NBands, (Zero,Zero), exxbuff(1,1,IKQ), NGrid )
  !
  DEALLOCATE( QRbuff, mat, pivot, list )
  !
  ! WRITE (stdout,'(7X,A)') 'SCDM-PGG_k done '
  !
  CALL stop_clock( 'localization' )
  !
END SUBROUTINE SCDM_PGG_k
!
!
!-----------------------------------------------------------------------------------------------------
SUBROUTINE scdm_prescreening_k( NGrid, NBands, psi, den, grad_den, ThrDen, ThrGrd, &
                                cpu_npt, nptot, list, pivot ) 
  !---------------------------------------------------------------------------------------------------
  !! Get List from ThrDen and ThrGrd, and Pivot from the QRCP of small.
  !
  USE mp,                ONLY : mp_sum 
  USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: list(nptot)
  INTEGER, INTENT(OUT) :: pivot(nptot)
  INTEGER, INTENT(IN) :: cpu_npt(0:nproc_bgrp-1)
  !! Number of relevant points per processor
  INTEGER, INTENT(IN) :: nptot
  !! Number of relevant points for the allocation
  INTEGER, INTENT(IN) :: NGrid
  !! Number of grid points
  INTEGER, INTENT(IN) :: NBands
  !! Number of bands
  COMPLEX(DP), INTENT(IN) :: psi(NGrid,NBands)
  REAL(DP), INTENT(IN) :: den(dfftt%nnr)
  REAL(DP), INTENT(IN) :: grad_den(3,dfftt%nnr)
  REAL(DP), INTENT(IN) :: ThrDen
  REAL(DP), INTENT(IN) :: ThrGrd
  !
  ! ... local variables
  !
  INTEGER :: ir, ir_end, INFO, lwork
  INTEGER :: n, npt, ncpu_start
  REAL(DP) :: grad
  COMPLEX(DP), ALLOCATABLE :: small(:,:), tau(:), work(:)
  REAL(DP), ALLOCATABLE :: rwork(:)
  !
#if defined (__MPI)
  ir_end = dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
#else
  ir_end = dfftt%nnr
#endif
  !
  ! find the map of the indeces
  ALLOCATE( small(NBands,nptot) )
  small = (Zero, Zero)
  list = 0
  n = 0
  !
  DO ir = 1, ir_end
     grad = One
     IF (den(ir) > ThrDen) THEN
        grad = SQRT( grad_den(1,ir)**2 + grad_den(2,ir)**2 + grad_den(3,ir)**2 ) 
        IF (grad < ThrGrd) THEN
           n = n + 1
           ncpu_start = SUM(cpu_npt(0:me_bgrp-1))
           npt = ncpu_start+n
           small(:,npt) = psi(ir,:)
           list(npt) = ir
        ENDIF 
     ENDIF
  ENDDO
  !
  CALL mp_sum( small, intra_bgrp_comm )
  CALL mp_sum( list,  intra_bgrp_comm )
  !
  ! perform the QRCP on the small matrix and get pivot
  lwork = 4*nptot
  ALLOCATE( tau(nptot), work(lwork), rwork(2*nptot) )
  tau = (Zero, Zero)
  work = (Zero, Zero)
  pivot = 0
  INFO = -1 
  CALL ZGEQP3( NBands, nptot, small, NBands, pivot, tau, work, lwork, rwork, INFO )
  DEALLOCATE( tau, work, rwork, small )
  !
END SUBROUTINE scdm_prescreening_k
!
!
!----------------------------------------------------------------------------------
SUBROUTINE scdm_fill_k( op, nptot, NGrid, NBands, CPUPts, Pivot, List, Vect, Mat )
  !--------------------------------------------------------------------------------
  !! Fills the matrix Mat with the elements of Vect mapped by CPUPts, Pivot
  !! and List.
  !
  USE mp,            ONLY: mp_sum
  USE mp_bands,      ONLY: intra_bgrp_comm, me_bgrp, nproc_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NBands
  INTEGER, INTENT(IN) :: NGrid
  INTEGER, INTENT(IN) :: nptot
  INTEGER, INTENT(IN) :: CPUPts(0:nproc_bgrp-1)
  INTEGER, INTENT(IN) :: Pivot(nptot)
  INTEGER, INTENT(IN) :: List(nptot)
  COMPLEX(DP), INTENT(IN)  :: Vect(NGrid,NBands)
  COMPLEX(DP), INTENT(OUT) :: Mat(NBands,NBands)
  !
  ! ... local variables
  !
  INTEGER :: i, NStart, NEnd
  LOGICAL :: op
  !
  Mat = (Zero, Zero)
  !
  DO i = 1, NBands
     NStart = SUM(CPUPts(0:me_bgrp-1))
     NEnd   = SUM(CPUPts(0:me_bgrp))
     IF (Pivot(i) <= NEnd .AND. Pivot(i) >= NStart+1) THEN
        IF (op) THEN
           Mat(:,i) = CONJG(Vect(List(pivot(i)),:))
        ELSE 
           Mat(:,i) = Vect(List(pivot(i)),:)
        ENDIF
     ENDIF 
  ENDDO 
  !
  CALL mp_sum( Mat, intra_bgrp_comm )
  !
END SUBROUTINE scdm_fill_k
!
!
!
END MODULE loc_scdm_k
!-----------------------------------------------------------------------
