! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
MODULE loc_scdm 
  !------------------------------------------------------------------
  !! Variables and subroutines for localizing molecular orbitals based
  !! on a modified SCDM approach.
  !
  !! Original SCDM (Selected Columns of the Density Matrix) method:  
  !! A. Damle, L. Lin, L. Ying: J. Chem. Theory Comput. 2015, 11, 1463.  
  !! Ivan Carnimeo: SCDM has been implemented with a parallel prescreening
  !! algorithm in order to save CPU time and memory for large grids.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE exx,                  ONLY : dfftt, x_nbnd_occ, locbuff, locmat
  USE exx_base,             ONLY : nkqs

  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: use_scdm = .FALSE.
  !! if .TRUE. enable Lin Lin's SCDM localization. 
  !! Currently implemented only within ACE formalism.
  INTEGER, ALLOCATABLE :: fragments(:)
  !! fragment list for dipole moments.
  REAL(DP) :: scdm_den
  !! scdm density
  REAL(DP) :: scdm_grd
  !! scdm gradient
  INTEGER :: n_scdm
  !! number of scdm
  INTEGER :: iscdm=0
  !! counter for QRCP/SVD 
  COMPLEX(DP), ALLOCATABLE :: locbuff_G(:,:,:)
  !! buffer for localized orbitals in G-space (MD only).
  REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Three=3.0d0
  !
 CONTAINS
  !
  !------------------------------------------------------------------------
SUBROUTINE localize_orbitals()
  !-----------------------------------------------------------------------
  !! Driver for SCDM orbital localization:  
  !! * iscdm/nscdm: control whether the localization is done via SCDM
  !!   (QR decomposition) or with orbital alignment (SVD);
  !! * SCDM_PGG: perform SCDM localization using the parallel
  !!   prescreening algorithm;
  !! * measure_localization: compute the absolute overlap integrals 
  !!   required by \texttt{vexx_loc} in \texttt{exx.f90} and put them
  !!   in locmat.
  !
  USE noncollin_module,     ONLY : npol
  USE wvfct,                ONLY : nbnd, npwx, current_k
  USE exx,                  ONLY : x_occupation
  USE wavefunctions,        ONLY : evc
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers,              ONLY : get_buffer
  USE control_flags,        ONLY : lmd
  USE funct,                ONLY : dft_is_hybrid
  !   
  implicit none
  integer :: NGrid, ikq, NBands, npw
  real(DP) :: tmp
  real(DP), allocatable :: MatQ(:,:)
  complex(DP), allocatable :: evcbuff(:,:), MatC(:,:)
  character(len=1) :: HowTo
  Logical :: QRCP
  
  call start_clock('localization')

  IF( lmd.and.(dft_is_hybrid()).and.(n_scdm.ne.1) ) &
         Call errore('localize_orbitals','MD+exx+nscdm NYI',1)

  QRCP = iscdm.eq.0.or.(mod(iscdm,n_scdm).eq.0) ! if .false. localize with SVD

  write(stdout,'(A,I10,A)') 'QRCP every ',n_scdm, ' steps.'
  write(stdout,'(A,I6,A,L1)')  'localize_orbitals: iscdm=',iscdm,' QRCP=',QRCP

  NGrid = dfftt%nnr * npol

! HowTo = 'R'  ! Compute the absolute overlap integrals in R-space (exact but
               ! EXTREMELY slow)
  HowTo = 'G'  ! How to compute the absolute overlap integrals (approx but fast)

  locmat = One ! initialized to one because the absolute overlaps between                         
               ! occupied and virtual orbitals are assumed always large 
               ! and all ov integrals will be computed by vexx_loc

  IF(.not.allocated(locbuff_G)) allocate( locbuff_G(npwx*npol, nbnd ,nks) )

  DO ikq = 1, nkqs
    npw = ngk (ikq)
    current_k = ikq
    NBands = int(sum(x_occupation(:,ikq)))
    allocate( MatQ(NBands,NBands), MatC(NBands,NBands), evcbuff(npwx*npol, nbnd)   )
    IF ( lsda ) current_spin = isk(ikq) 
    IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ikq)
    locmat(:,:,ikq) = One
    CALL measure_localization(HowTo,NBands,ikq)  ! compute: 
                                                 ! the matrix of absolute overlap integrals
                                                 ! average spread of the localized functions
!   Call AbsOv_histogram(ikq,NBands,'hist1.dat')
    IF(QRCP) THEN 
      CALL SCDM_PGG(locbuff(1,1,ikq), NGrid, NBands)  ! localization with SCDM in R-space
      CALL wave_to_G(locbuff(1,1,ikq), locbuff_G(1,1,ikq), NGrid, NBands) ! bring the localized functions to G-space
    END IF
!   FROM HERE: align the orbitals in evc (step n) to the orbitals in locbuff_G (step n-1)
!   MatQ contains the overlap between orbitals in evc and orbitals in locbuff_G  
    Call matcalc ('<Psi2|W1>',   .false.,   0, npw, NBands, NBands, locbuff_G(1,1,ikq), evc, MatQ, tmp)
    Call MatCheck(MatQ, NBands)  ! Orthonormalization check
    Call PTSVD(MatQ, NBands)  ! SVD of MatQ
    MatC = (One,Zero)*MatQ
!   align the orbitals in evc rotating them with the unitary matrix in MatC 
    Call ZGEMM('N','T',npw,NBands,NBands,One,evc,npw,MatC,NBands,Zero,evcbuff,npw) 
    CALL matcalc('<W(t)|W(t)>-',.false.,0,npw,NBands,NBands,locbuff_G(1,1,ikq),evcbuff,MatQ,tmp)
    Call MatCheck(MatQ, NBands) ! Orthonormalization check 
!   TO HERE: align the orbitals in evc (step n) to the orbitals in locbuff_G (step n-1)

!   IF(QRCP) just delete locbuff_G. 
    IF(.not.QRCP) THEN  ! save locbuff_G for the next iteration  
      CALL wave_to_R(evcbuff, locbuff(1,1,1), NGrid, NBands)
      locbuff_G(1:npwx*npol,1:nbnd,ikq) = evcbuff  ! save the orbitals at step n for later use at step n+1
    END IF
    deallocate( MatQ, MatC, evcbuff )

    locmat(:,:,ikq) = One
    CALL measure_localization(HowTo,NBands,ikq)
!   Call AbsOv_histogram(ikq,NBands,'hist2.dat')

  END DO

  iscdm = iscdm + 1 

  call stop_clock('localization')

END SUBROUTINE localize_orbitals
!
!------------------------------------------------------------------------
SUBROUTINE AbsOv_histogram( iik, n, filename )
  !------------------------------------------------------------------------
  !! Print the distribution of absolute overlap integrals (as an histogram)
  !! in the range [0,1] with NHist intervals.
  !
  !! WARNING: HUGE amount of integrals, use only for VERY small systems.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: filename
  !! file where the distribution is printed
  INTEGER, INTENT(IN) :: n
  !! dimension of the absolute overlap matrix
  INTEGER, INTENT(IN) :: iik
  !! k+q point index (see \texttt{exx_base})
  !
  ! ... local variables
  !
  INTEGER :: i,j,k, NHist
  INTEGER, ALLOCATABLE :: histogram(:)
  REAL(DP), ALLOCATABLE :: XHist(:)
  REAL(DP) :: xstart, xstep, integral
  INTEGER :: io_histogram
  INTEGER, EXTERNAL :: find_free_unit
  !
  NHist = 1000
  xstep = One/float(NHist)
  xstart = One/float(NHist)/Two
  write(stdout,'(A,I7,2(A,f12.6))') 'NHist  = ', NHist , ' xstep = ', xstep, ' xstart = ', xstart
  ALLOCATE(histogram(NHist),XHist(NHist))
  XHist = Zero
  histogram = 0
  do k = 1, NHist 
    XHist(k) = xstart + float(k-1)*xstep 
  end do  
  do i = 1, n 
    do j = 1, n
      integral = locmat(i,j,iik)
      if(integral.lt.Zero) then 
        Call errore('AbsOv_histogram','Abs. Ov. < 0 found.',1)
      else
        do k = 1, NHist 
          IF(integral.ge.(XHist(k)-xstart).and.integral.lt.(XHist(k)+xstart)) &
               histogram(k)=histogram(k)+1
        end do 
      end if 
    end do 
  end do 
  io_histogram = find_free_unit()
  open(io_histogram,file=filename,status='unknown')
  do k = 1, NHist
    write(io_histogram,'(f12.6,2I10)') XHist(k), histogram(k)
  end do 
  close(io_histogram)
  DEALLOCATE(histogram,XHist)
  !
END SUBROUTINE AbsOv_histogram
!
!--------------------------------------------------------------------
SUBROUTINE measure_localization( CFlag, NBands, IKK )
  !---------------------------------------------------------------------
  !! Analyze various localization criteria for the wavefunctions in orbt.
  !
  USE noncollin_module,  ONLY : npol
  USE cell_base,         ONLY : alat, omega, at, bg
  USE exx,               ONLY : compute_density
  USE constants,         ONLY : bohr_radius_angs 
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1) :: CFlag
  !! Calculation of absolute overlap:
  !! * CFlag = 'R' real space integral (exact but slow);
  !! * CFlag = 'G' via FFT (fast but less accurate).
  INTEGER :: NBands
  !! number of bands
  INTEGER :: IKK
  !! k+q point index (see \texttt{exx_base})
  !
  ! ... local variables
  !
  INTEGER :: jbnd, kbnd
  REAL(DP) :: loc_diag, loc_off, tmp, DistMax
  REAL(DP) :: RDist(3),  SpreadPBC(3), TotSpread
  REAL(DP), ALLOCATABLE :: CenterPBC(:,:), Mat(:,:)
  
  REAL(DP), PARAMETER :: epss=0.0010d0

  ALLOCATE( Mat(NBands,NBands), CenterPBC(3,NBands) )

  IF(CFlag.eq.'R') then 
    Call AbsOvR(NBands, IKK, Mat)  ! exact but slow
  ELSEIF(CFlag.eq.'G') then     
    call AbsOvG(NBands, IKK, Mat)  ! approx but fast (default)
  ELSE
    call errore('measure_localization','Wrong CFlag',1)
  END IF 

  loc_diag  = Zero  
  loc_off   = Zero  
  TotSpread = Zero  
  DistMax   = Zero 
  DO jbnd = 1, NBands
    loc_diag = loc_diag + Mat(jbnd,jbnd) 
    call compute_density(.false.,.false.,CenterPBC(1,jbnd), SpreadPBC, tmp, &
                         locbuff(1,jbnd,IKK), locbuff(1,jbnd,IKK), &
                         dfftt%nnr*npol, jbnd, jbnd)
    TotSpread = TotSpread + SpreadPBC(1) + SpreadPBC(2) + SpreadPBC(3) 
    DO kbnd = 1, jbnd - 1 
      loc_off = loc_off + Mat(jbnd,kbnd) 
      RDist(1) = (CenterPBC(1,jbnd) - CenterPBC(1,kbnd))/alat  
      RDist(2) = (CenterPBC(2,jbnd) - CenterPBC(2,kbnd))/alat  
      RDist(3) = (CenterPBC(3,jbnd) - CenterPBC(3,kbnd))/alat  
      CALL cryst_to_cart( 1, RDist, bg, -1 )
      RDist(:) = RDist(:) - ANINT( RDist(:) )
      CALL cryst_to_cart( 1, RDist, at, 1 )
      tmp = alat *bohr_radius_angs* sqrt( RDist(1)**2 + RDist(2)**2 + RDist(3)**2 )
      IF(DistMax.lt.tmp) DistMax = tmp  
    ENDDO
  ENDDO
  write(stdout,'(7X,A,f12.6,A)')  'Max Dist [A]      = ', bohr_radius_angs*alat*sqrt(Three)/Two, ' (sqrt(3)*L/2)'
  write(stdout,'(7X,A,f12.6)')    'Max Dist Found [A] =', DistMax 
  write(stdout,'(7X,A,f12.6,I3)') 'Total Charge =', loc_diag 
  write(stdout,'(7X,A,f12.6,I3)') 'Total Abs. Overlap =', loc_off 
  write(stdout,'(7X,A,f12.6,I3)') 'Total Spread [A**2]   =', TotSpread * bohr_radius_angs**2
  write(stdout,'(7X,A,f12.6,I3)') 'Aver. Spread [A**2]   =', TotSpread * bohr_radius_angs**2/ dble(NBands) 
  locmat(1:NBands,1:NBands,IKK) = Mat 
  DEALLOCATE( CenterPBC, Mat ) 

END SUBROUTINE measure_localization
!
!--------------------------------------------------------------------------
SUBROUTINE AbsOvG( NBands, IKK, Mat )
  !----------------------------------------------------------------------
  !! Compute the Absolute Overlap in G-space (cutoff might not be accurate
  !! for the moduli of the wavefunctions).
  !
  USE noncollin_module,  ONLY : npol
  USE fft_interfaces,    ONLY : fwfft
  USE wvfct,             ONLY : npwx
  !
  IMPLICIT NONE
  !
  INTEGER :: NBands
  !! Number of bands
  INTEGER :: IKK
  !! k+q point index (see \texttt{exx_base})
  REAL(DP) :: Mat(NBands,NBands)
  !! Absolute overlap matrix
  !
  ! ... local variables
  !
  INTEGER :: jbnd, ig
  REAL(DP) :: tmp
  COMPLEX(DP), ALLOCATABLE :: buffer(:), Gorbt(:,:)

  call start_clock('measure')

  write(stdout,'(5X,A)') ' ' 
  write(stdout,'(5X,A)') 'Absolute Overlap calculated in G-space'

! Localized functions to G-space and Overlap matrix onto localized functions
  allocate( buffer(dfftt%nnr * npol), Gorbt(npwx,NBands) )

  Mat = Zero 
  buffer = (Zero,Zero)
  Gorbt = (Zero,Zero) 
  DO jbnd = 1, NBands 
    buffer(:) = abs(dble(locbuff(:,jbnd, IKK))) + (Zero,One)*Zero  
    CALL fwfft( 'Wave' , buffer, dfftt )
    DO ig = 1, npwx
      Gorbt(ig,jbnd) = buffer(dfftt%nl(ig))
    ENDDO
  ENDDO
  CALL matcalc('Coeff-',.false.,0,npwx,NBands,NBands,Gorbt,Gorbt,Mat,tmp)
  deallocate ( buffer, Gorbt )

  call stop_clock('measure')

END SUBROUTINE AbsOvG 
!
!----------------------------------------------------------------------
SUBROUTINE AbsOvR( NBands, IKK, Mat )
  !-------------------------------------------------------------------
  !! Compute the Absolute Overlap in R-space (Exact but slow).
  !
  USE mp,                ONLY : mp_sum
  USE mp_bands,          ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: NBands
  !! Number of bands
  INTEGER :: IKK
  !! k+q point index (see \texttt{exx_base})
  REAL(DP) :: Mat(NBands,NBands)
  !! Absolute overlap matrix
  !
  ! ... local variables
  !
  INTEGER :: nxxs, ir, jbnd, kbnd
  REAL(DP) ::  cost, tmp

  call start_clock('measure')

  write(stdout,'(5X,A)') ' ' 
  write(stdout,'(5X,A)') 'Absolute Overlap calculated in R-space'

  nxxs = dfftt%nr1x *dfftt%nr2x *dfftt%nr3x
  cost = One/dble(nxxs)
  Mat = Zero 
  DO jbnd = 1, NBands
    Mat(jbnd,jbnd) = Mat(jbnd,jbnd) + cost * sum( abs(locbuff(:,jbnd,IKK)) * abs(locbuff(:,jbnd,IKK)))
    DO kbnd = 1, jbnd - 1 
        tmp = cost * sum( abs(locbuff(:,jbnd,IKK)) * abs(locbuff(:,kbnd,IKK)) )
        Mat(jbnd,kbnd) = Mat(jbnd,kbnd) + tmp 
        Mat(kbnd,jbnd) = Mat(kbnd,jbnd) + tmp
    ENDDO
  ENDDO
  call mp_sum(mat,intra_bgrp_comm)

  call stop_clock('measure')

END SUBROUTINE AbsOvR 
!
!-----------------------------------------------------------------
SUBROUTINE SCDM_PGG( psi, NGrid, NBands )
  !-------------------------------------------------------------
  !! Density matrix localization (I/O in psi) - Gamma version.
  !
  USE mp_bands,          ONLY : nproc_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: NGrid
  !! grid points
  INTEGER, INTENT(IN) :: NBands
  !! Number of bands
  REAL(DP), INTENT(INOUT) :: psi(NGrid,NBands)
  !! localized functions as output
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: QRbuff(:,:), mat(:,:)
  INTEGER,  ALLOCATABLE :: pivot(:), list(:), cpu_npt(:)
  REAL(DP), ALLOCATABLE :: den(:), grad_den(:,:)
  INTEGER  :: nptot
  REAL(DP) :: ThrDen, ThrGrd

  write(stdout,'(5X,A)') ' ' 
  write(stdout,'(5X,A)') 'SCDM localization with prescreening'

  allocate( den(dfftt%nnr), grad_den(3, dfftt%nnr) )

  Call scdm_thresholds( den, grad_den, ThrDen, ThrGrd )

  allocate( cpu_npt(0:nproc_bgrp-1) )

  Call scdm_points( den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot )

  allocate( list(nptot), pivot(nptot) )

  Call scdm_prescreening( NGrid, NBands, psi, den, grad_den, ThrDen, ThrGrd, &
                                                     cpu_npt, nptot, list, pivot )

  deallocate( den, grad_den ) 

! Psi(pivot(1:NBands),:) in mat
  allocate( mat(NBands,NBands) )
  Call scdm_fill( nptot, NGrid, NBands, cpu_npt, pivot, list, psi, Mat)

! Pc = Psi * Psi(pivot(1:NBands),:)' in QRbuff
  allocate( QRbuff(NGrid, NBands) )
  QRbuff = Zero 
  CALL DGEMM( 'N' , 'N' , NGrid, NBands, NBands, One, psi, NGrid, mat, NBands, Zero, QRbuff, NGrid) 

! Orthonormalization

! Pc(pivot(1:NBands),:) in mat 
  Call scdm_fill( nptot, NGrid, NBands, cpu_npt, pivot, list, QRBuff, mat)
  deallocate( cpu_npt )

! Cholesky(psi)^(-1) in mat 
! CALL invchol(NBands,mat)
  CALL MatChol(NBands,mat)
  CALL MatInv('L',NBands,mat)
  Call MatSymm('U','L',mat, NBands)

! Phi = Pc * Chol^(-1) = QRbuff * mat
  psi = Zero
  CALL DGEMM( 'N' , 'N' , NGrid, NBands, NBands, One, QRbuff, NGrid, mat, NBands, Zero, psi, NGrid) 
  deallocate( QRbuff, mat, pivot, list )
  write(stdout,'(7X,A)') 'SCDM-PGG done ' 

END SUBROUTINE SCDM_PGG
!
!---------------------------------------------------------------
SUBROUTINE scdm_thresholds( den, grad_den, ThrDen, ThrGrd )
  !-----------------------------------------------------------------
  !! Interpolate density to the exx grid.
  !
  USE cell_base,         ONLY : omega
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fft_interpolate
  USE scf,               ONLY : rho
  USE lsda_mod,          ONLY : nspin
  USE mp,                ONLY : mp_sum, mp_max
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE exx,               ONLY : gt
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: den(dfftt%nnr)
  !! density
  REAL(DP), INTENT(OUT) :: grad_den(3,dfftt%nnr)
  !! gradient of the density
  REAL(DP), INTENT(OUT) :: ThrDen
  !! Density Threshold
  REAL(DP), INTENT(OUT) :: ThrGrd
  !! Gradient Threshold
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: temp(:)
  REAL(DP) :: charge, grad, DenAve, GrdAve, DenMax, GrdMax
  INTEGER :: ir, ir_end, nxxs, nxtot
  !
  ! interpolate density to the exx grid
  allocate( temp(dfftp%nnr))
  temp(:) = rho%of_r(:,1)
  Call fft_interpolate(dfftp, temp, dfftt, den)
  deallocate( temp ) 

#if defined (__MPI)
  ir_end = dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
#else
  ir_end = dfftt%nnr
#endif
  nxtot = dfftt%nr1x *dfftt%nr2x *dfftt%nr3x
  nxxs = dfftt%nnr 

  charge = Zero
  DenAve = Zero
  DenMax = Zero
  do ir = 1, ir_end 
    charge = charge + den(ir) * omega / dble(nxtot) 
    DenAve = DenAve + den(ir)
    IF(DenMax.lt.den(ir)) DenMax=den(ir)
  end do 
  call mp_sum(DenAve,intra_bgrp_comm)
  call mp_sum(charge,intra_bgrp_comm)
  call mp_max(DenMax,intra_bgrp_comm)
  DenAve = DenAve / dble(nxtot)
! write(stdout,'(7x,A,f12.6)') 'Charge  = ', charge
! write(stdout,'(7x,A,f12.6)') 'DenAve  = ', DenAve 
! write(stdout,'(7x,A,f12.6)') 'DenMax  = ', DenMax 

! gradient on the exx grid 
  call fft_gradient_r2r ( dfftt, den, gt, grad_den )
  charge  = Zero
  GrdAve = Zero 
  GrdMax = Zero 
  do ir = 1, ir_end 
    grad  = sqrt( grad_den(1,ir)**2  +  grad_den(2,ir)**2  +  grad_den(3,ir)**2  )
    charge = charge + grad * omega / dble(nxtot)
    GrdAve = GrdAve + grad
    IF(GrdMax.lt.grad) GrdMax=grad
  end do 
  call mp_sum(GrdAve,intra_bgrp_comm)
  call mp_sum(charge,intra_bgrp_comm)
  call mp_max(GrdMax,intra_bgrp_comm)
  GrdAve = GrdAve / dble(nxtot)
! write(stdout,'(7X,A,f12.6)') 'GradTot = ', charge
! write(stdout,'(7X,A,f12.6)') 'GrdAve  = ', GrdAve 
! write(stdout,'(7X,A,f12.6)') 'GrdMax  = ', GrdMax 

  ThrDen = scdm_den * DenAve  
  ThrGrd = scdm_grd * GrdAve  
! write(stdout,'(7x,2(A,f12.6))') 'scdm_den = ', scdm_den, ' scdm_grd = ',scdm_grd
! write(stdout,'(7x,2(A,f12.6))') 'ThrDen   = ', ThrDen,   ' ThrGrd   = ',ThrGrd  

END SUBROUTINE scdm_thresholds

!---------------------------------------------------------------------------
SUBROUTINE scdm_fill( nptot, NGrid, NBands, CPUPts, Pivot, List, Vect, Mat )
  !-------------------------------------------------------------------------
  !! Fill the matrix Mat with the elements of Vect mapped by CPUPts, 
  !! Pivot and List.
  !
  USE mp,                ONLY : mp_sum
  USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NBands
  ! number of bands
  INTEGER, INTENT(IN) :: NGrid
  ! grid points
  INTEGER, INTENT(IN) :: nptot
  ! total number of relevant points
  INTEGER, INTENT(IN) :: CPUPts(0:nproc_bgrp-1)
  ! relevant points per CPU
  INTEGER, INTENT(IN) :: Pivot(nptot)
  INTEGER, INTENT(IN) :: List(nptot)
  REAL(DP), INTENT(IN) :: Vect(NGrid,NBands)
  REAL(DP), INTENT(OUT) :: Mat(NBands,NBands)
  !
  ! ... local variables
  !
  INTEGER :: i, NStart, NEnd

  Mat = Zero 
  do i = 1, NBands
    NStart = sum(CPUPts(0:me_bgrp-1))
    NEnd   = sum(CPUPts(0:me_bgrp))
    if(Pivot(i).le.NEnd.and.Pivot(i).ge.NStart+1) Mat(:,i) = Vect(List(pivot(i)),:)
  end do
  call mp_sum(Mat,intra_bgrp_comm)

END SUBROUTINE scdm_fill
!
!-------------------------------------------------------------------------
SUBROUTINE scdm_points( den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot )
  !------------------------------------------------------------------------
  !! Find the relevant points for the allocation.
  !
  USE mp,                ONLY : mp_sum 
  USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: den(dfftt%nnr)
  !! density
  REAL(DP), INTENT(IN) :: grad_den(3,dfftt%nnr)
  !! density gradient
  REAL(DP), INTENT(IN) :: ThrDen
  !! density threshold
  REAL(DP), INTENT(IN) :: ThrGrd 
  !! gradient threshold
  INTEGER, INTENT(OUT) :: cpu_npt(0:nproc_bgrp-1)
  !! number of relevant points per cpu
  INTEGER, INTENT(OUT) :: nptot
  !! total number of relevant points
  !
  ! ... local variables
  !
  INTEGER :: npt, ir, ir_end
  REAL(DP) :: grad

#if defined (__MPI)
  ir_end = dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
#else
  ir_end = dfftt%nnr
#endif

  npt = 0
  cpu_npt(:) = 0
  do ir = 1, ir_end 
    if(den(ir).gt.ThrDen) then 
      grad = sqrt( grad_den(1,ir)**2 +  grad_den(2,ir)**2 +  grad_den(3,ir)**2 ) 
      if(grad.lt.ThrGrd) then    
        npt = npt + 1
      end if 
    end if 
  end do 
  nptot = npt
  cpu_npt(me_bgrp) = npt
  call mp_sum(nptot,intra_bgrp_comm)
  if(nptot.le.0) call errore('SCDM_PGG', 'No points prescreened. Loose the thresholds', 1) 
  call mp_sum(cpu_npt,intra_bgrp_comm)
! write(stdout,'(7X,2(A,I8))')  'Max npt = ', maxval(cpu_npt(:)), ' Min npt = ', minval(cpu_npt(:))
! write(stdout,'(7X,2(A,I10))') 'Reduced matrix, allocate: ', nptot, ' out of ', dfftt%nr1x *dfftt%nr2x *dfftt%nr3x
!
END SUBROUTINE scdm_points
!
!-----------------------------------------------------------------------------
SUBROUTINE scdm_prescreening( NGrid, NBands, psi, den, grad_den, ThrDen, &
                              ThrGrd, cpu_npt, nptot, list, pivot  ) 
  !--------------------------------------------------------------------------
  !! Get List from ThrDen and ThrGrd, and Pivot from the QRCP of small.
  !
  USE mp,                ONLY : mp_sum 
  USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: cpu_npt(0:nproc_bgrp-1), nptot
  INTEGER, INTENT(IN) :: NGrid, NBands
  INTEGER, INTENT(OUT) :: list(nptot), pivot(nptot)
  REAL(DP), INTENT(IN) :: psi(NGrid,NBands) 
  REAL(DP), INTENT(IN) :: den(dfftt%nnr), grad_den(3, dfftt%nnr) 
  REAL(DP), INTENT(IN) :: ThrDen, ThrGrd 
  !
  INTEGER :: ir, ir_end, INFO, lwork
  INTEGER :: n, npt, ncpu_start
  REAL(DP) :: grad
  REAL(DP), ALLOCATABLE :: small(:,:), tau(:), work(:)
  !
#if defined (__MPI)
  ir_end = dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
#else
  ir_end = dfftt%nnr
#endif
  !
  ! find the map of the indeces
  allocate( small(NBands,nptot) )
  small = Zero
  list = 0
  n = 0
  do ir = 1, ir_end
    grad = One
    if(den(ir).gt.ThrDen) then 
      grad = sqrt( grad_den(1,ir)**2 +  grad_den(2,ir)**2 +  grad_den(3,ir)**2 ) 
      if(grad.lt.ThrGrd) then    
        n = n + 1
        ncpu_start = sum(cpu_npt(0:me_bgrp-1))
        npt = ncpu_start+n
        small(:,npt) = psi(ir,:)
        list(npt) = ir
      end if 
    end if 
  end do 
  call mp_sum(small,intra_bgrp_comm)
  call mp_sum(list,intra_bgrp_comm)

! perform the QRCP on the small matrix and get pivot
  lwork = 4*nptot
  allocate( tau(nptot), work(lwork) )
  tau = Zero
  work = Zero  
  pivot = 0
  INFO = -1 
  CALL DGEQP3( NBands, nptot, small, NBands, pivot, tau, work, lwork, INFO )
  deallocate( tau, work, small )

END SUBROUTINE scdm_prescreening
!
!--------------------------------------------------------------------------
SUBROUTINE wave_to_R( psiG, psiR, NGrid, NBands )
  !----------------------------------------------------------------------
  !! psi functions from G- to R-space.
  !
  USE cell_base,         ONLY : omega 
  USE kinds,             ONLY : DP
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE wvfct,             ONLY : npwx
  USE exx,               ONLY : npwt
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: psiG(npwx,NBands)
  !! G-space psi
  REAL(DP), INTENT(OUT) :: psiR(NGrid,NBands)
  !! R-space psi
  INTEGER :: NGrid
  !! grid points
  INTEGER :: NBands
  !! number of bands
  !
  ! ... local variables
  !
  INTEGER :: ir, jbnd, ig
  REAL(DP) :: tmp
  COMPLEX(DP), ALLOCATABLE :: buffer(:)
  REAL(DP), ALLOCATABLE :: Mat(:,:)
  !
  write(stdout,'(A)') 'Wave to R '
  !
  allocate( buffer(NGrid) )
  psiR = Zero
  DO jbnd = 1, NBands 
    buffer = (Zero,Zero)
    DO ig=1,npwt
       buffer(dfftt%nl (ig)) = psiG(ig,jbnd)
       buffer(dfftt%nlm(ig)) = conjg( psiG(ig,jbnd) )
    ENDDO
    CALL invfft( 'Wave' , buffer, dfftt )
    psiR(1:NGrid,jbnd) = dble(buffer(1:NGrid))
  ENDDO
  deallocate ( buffer )
! allocate( buffer(NGrid), Mat(NBands,NBands) )
! tmp = One/dble(NGrid)
! CALL DGEMM( 'T' , 'N' , NBands, NBands, NGrid, tmp, PsiR, NGrid, PsiR, NGrid, Zero, Mat, NBands) 
! Call MatPrt('check',NBands,NBands,MAt)
! deallocate ( Mat )
!
END SUBROUTINE wave_to_R
!
!-------------------------------------------------------------------------
SUBROUTINE wave_to_G( psiR, psiG, NGrid, NBands )
  !----------------------------------------------------------------------
  !! psi functions from R- to G-space.
  !
  USE noncollin_module,     ONLY : npol
  USE kinds,                ONLY : DP
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE wvfct,                ONLY : npwx
  !
  IMPLICIT NONE
  !
  INTEGER :: NGrid
  !! grid points
  INTEGER :: NBands
  !! number of bands
  REAL(DP), INTENT(IN) :: psiR(NGrid,NBands)
  !! R-space psi
  COMPLEX(DP), INTENT(OUT) :: psiG(npwx*npol,NBands)
  !! G-space psi
  !
  ! ... local variables
  !
  REAL(DP) :: tmp
  INTEGER :: jbnd, ig
  COMPLEX(DP), ALLOCATABLE :: buffer(:)
  REAL(DP), ALLOCATABLE :: Mat(:,:)

  write(stdout,'(A)') 'Wave to G '

  allocate( buffer(NGrid) )
  buffer = (Zero,Zero)
  psiG = (Zero,Zero) 
  DO jbnd = 1, NBands 
    buffer(:) = dble(psiR(:,jbnd)) + (Zero,One)*Zero
    CALL fwfft( 'Wave' , buffer, dfftt )
    DO ig = 1, npwx
      psiG(ig,jbnd) = buffer(dfftt%nl(ig))
    ENDDO
  ENDDO
  deallocate ( buffer )
  allocate( Mat(NBands,NBands) )
  CALL matcalc('Check',.true.,1,npwx,NBands,NBands,psiG,psiG,Mat,tmp)
  deallocate ( Mat )

END SUBROUTINE wave_to_G 
!
!
END MODULE loc_scdm 
