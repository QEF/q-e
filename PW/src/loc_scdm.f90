! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE loc_scdm 
  !--------------------------------------
  !
  ! Variables and subroutines for localizing molecular orbitals based
  ! on a modified SCDM approach. Original SCDM method:
  ! A. Damle, L. Lin, L. Ying: J. Chem. Theory Comput. 2015, 11, 1463
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE exx,                  ONLY : dfftt, locbuff, locmat, nkqs

  IMPLICIT NONE
  SAVE
  LOGICAL ::  use_scdm=.false. ! if .true. enable Lin Lin's SCDM localization
                               ! currently implemented only within ACE formalism
  REAL(DP) :: scdm_den, scdm_grd 
  REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Three=2.0d0

 CONTAINS
  !
  !------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE localize_orbitals( )
  !
  ! Driver for SCDM orbital localization 
  !
  USE noncollin_module,  ONLY : npol
  USE wvfct,             ONLY : nbnd
  USE control_flags,     ONLY : gamma_only
  USE exx,               ONLY : x_occupation
  !   
  implicit none
  integer :: NGrid, ikq, NBands
  character(len=1) :: HowTo
  
  if(.not.gamma_only) CALL errore('localize_orbitals', 'k-points NYI.',1)    

  NGrid = dfftt%nnr * npol
  HowTo = 'G'  ! How to compute the absolute overlap integrals

  DO ikq = 1, nkqs
    NBands = int(sum(x_occupation(:,ikq)))
    locmat(:,:,ikq) = One
    CALL measure_localization(HowTo,NBands,ikq)
    CALL SCDM_PGG(locbuff(1,1,ikq), NGrid, NBands)
    locmat(:,:,ikq) = One
    CALL measure_localization(HowTo,NBands,ikq)
  END DO

END SUBROUTINE localize_orbitals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE measure_localization(CFlag, NBands, IKK) 
USE noncollin_module,  ONLY : npol
USE cell_base,         ONLY : alat, omega, at, bg
USE exx,               ONLY : compute_density
USE constants,         ONLY : bohr_radius_angs 
implicit none
!
! Analyze various localization criteria for the wavefunctions in orbt
! Calculation of absolute overlap: 
!              CFlag = 'R' real space integral (exact but slow) 
!                      'G' via FFT (fast but less accurate) 
!
  INTEGER :: NBands, jbnd, kbnd, IKK
  REAL(DP) :: loc_diag, loc_off, tmp, DistMax
  REAL(DP) :: RDist(3),  SpreadPBC(3), TotSpread
  REAL(DP), ALLOCATABLE :: CenterPBC(:,:), Mat(:,:)
  CHARACTER(LEN=1) :: CFlag
  REAL(DP), PARAMETER :: epss=0.0010d0

  ALLOCATE( Mat(NBands,NBands), CenterPBC(3,NBands) )

  IF(CFlag.eq.'R') then 
    Call AbsOvR(NBands, IKK, Mat) 
  ELSEIF(CFlag.eq.'G') then 
    call AbsOvG(NBands, IKK, Mat) 
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AbsOvG(NBands, IKK, Mat) 
USE noncollin_module,  ONLY : npol
USE fft_interfaces,    ONLY : fwfft
USE wvfct,             ONLY : npwx
implicit none
!
! Compute the Absolute Overlap in G-space 
! (cutoff might not be accurate for the moduli of the wavefunctions)
!
  INTEGER :: NBands, jbnd, ig, IKK
  REAL(DP) :: Mat(NBands,NBands), tmp
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AbsOvR(NBands, IKK, Mat) 
USE mp,                ONLY : mp_sum
USE mp_bands,          ONLY : intra_bgrp_comm
implicit none
!
! Compute the Absolute Overlap in R-space 
! (Exact but slow)
!
  INTEGER :: NBands, nxxs, ir, jbnd, kbnd, IKK
  REAL(DP) :: Mat(NBands,NBands)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCDM_PGG(psi, NGrid, NBands)
USE mp_bands,          ONLY : nproc_bgrp
!
! density matrix localization (I/O in psi)
!
IMPLICIT NONE
  INTEGER,  INTENT(IN)    :: NGrid, NBands
  REAL(DP), INTENT(INOUT) :: psi(NGrid,NBands)

  REAL(DP), ALLOCATABLE :: QRbuff(:,:), mat(:,:)
  INTEGER,  ALLOCATABLE :: pivot(:), list(:), cpu_npt(:)
  REAL(DP), ALLOCATABLE :: den(:), grad_den(:,:)
  INTEGER  :: nptot
  REAL(DP) :: ThrDen, ThrGrd

  call start_clock('localization')

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
  CALL invchol(NBands,mat)
  Call MatSymm('U','L',mat, NBands)

! Phi = Pc * Chol^(-1) = QRbuff * mat
  psi = Zero
  CALL DGEMM( 'N' , 'N' , NGrid, NBands, NBands, One, QRbuff, NGrid, mat, NBands, Zero, psi, NGrid) 
  deallocate( QRbuff, mat, pivot, list )
  write(stdout,'(7X,A)') 'SCDM-PGG done ' 

  call stop_clock('localization')

END SUBROUTINE SCDM_PGG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE scdm_thresholds( den, grad_den, ThrDen, ThrGrd )
USE cell_base,         ONLY : omega
USE fft_base,          ONLY : dfftp
USE fft_interfaces,    ONLY : fft_interpolate
USE scf,               ONLY : rho
USE lsda_mod,          ONLY : nspin
USE mp,                ONLY : mp_sum, mp_max
USE mp_bands,          ONLY : intra_bgrp_comm
USE exx,               ONLY : gt
IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: den(dfftt%nnr), grad_den(3, dfftt%nnr) 
  REAL(DP), INTENT(OUT) :: ThrDen, ThrGrd 

  REAL(DP), ALLOCATABLE :: temp(:) 
  REAL(DP) :: charge, grad, DenAve, GrdAve, DenMax, GrdMax
  INTEGER :: ir, ir_end, nxxs, nxtot

! interpolate density to the exx grid
  allocate( temp(dfftp%nnr))
  temp(:) = rho%of_r(:,1)
  IF ( nspin == 2 ) temp(:) = temp(:) + rho%of_r(:,2) 
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
  write(stdout,'(7x,A,f12.6)') 'Charge  = ', charge
  write(stdout,'(7x,A,f12.6)') 'DenAve  = ', DenAve 
  write(stdout,'(7x,A,f12.6)') 'DenMax  = ', DenMax 

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
  write(stdout,'(7X,A,f12.6)') 'GradTot = ', charge
  write(stdout,'(7X,A,f12.6)') 'GrdAve  = ', GrdAve 
  write(stdout,'(7X,A,f12.6)') 'GrdMax  = ', GrdMax 

  ThrDen = scdm_den * DenAve  
  ThrGrd = scdm_grd * GrdAve  
  write(stdout,'(7x,2(A,f12.6))') 'scdm_den = ', scdm_den, ' scdm_grd = ',scdm_grd
  write(stdout,'(7x,2(A,f12.6))') 'ThrDen   = ', ThrDen,   ' ThrGrd   = ',ThrGrd  

END SUBROUTINE scdm_thresholds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE scdm_fill( nptot, NGrid, NBands, CPUPts, Pivot, List, Vect, Mat)
USE mp,                ONLY : mp_sum
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
!
! Fill the matrix Mat with the elements of Vect
! mapped by CPUPts, Pivot and List 
!
IMPLICIT NONE
  INTEGER,  INTENT(IN)  :: NBands, NGrid, nptot
  INTEGER,  INTENT(IN)  :: CPUPts(0:nproc_bgrp-1), Pivot(nptot), List(nptot)
  REAL(DP), INTENT(IN)  :: Vect(NGrid, NBands)
  REAL(DP), INTENT(OUT) :: Mat(NBands,NBands)
  INTEGER :: i, NStart, NEnd

  Mat = Zero 
  do i = 1, NBands
    NStart = sum(CPUPts(0:me_bgrp-1))
    NEnd   = sum(CPUPts(0:me_bgrp))
    if(Pivot(i).le.NEnd.and.Pivot(i).ge.NStart+1) Mat(:,i) = Vect(List(pivot(i)),:)
  end do 
  call mp_sum(Mat,intra_bgrp_comm)

END SUBROUTINE scdm_fill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE scdm_points ( den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot  ) 
USE mp,                ONLY : mp_sum 
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
!
! find the relevant points for the allocation
!
IMPLICIT NONE
  INTEGER, INTENT(OUT) :: cpu_npt(0:nproc_bgrp-1), nptot
  REAL(DP), INTENT(IN) :: den(dfftt%nnr), grad_den(3, dfftt%nnr) 
  REAL(DP), INTENT(IN) :: ThrDen, ThrGrd 

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
  write(stdout,'(7X,2(A,I8))')  'Max npt = ', maxval(cpu_npt(:)), ' Min npt = ', minval(cpu_npt(:))
  write(stdout,'(7X,2(A,I10))') 'Reduced matrix, allocate: ', nptot, ' out of ', dfftt%nr1x *dfftt%nr2x *dfftt%nr3x

END SUBROUTINE scdm_points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE scdm_prescreening ( NGrid, NBands, psi, den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot, &
                                                                                       list, pivot  ) 
USE mp,                ONLY : mp_sum 
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
!
!  Get List from ThrDen and ThrGrd, and Pivot from the QRCP of small
!
IMPLICIT NONE
  INTEGER, INTENT(OUT) :: list(nptot), pivot(nptot)
  INTEGER, INTENT(IN)  :: cpu_npt(0:nproc_bgrp-1), nptot
  INTEGER, INTENT(IN)  :: NGrid, NBands
  REAL(DP), INTENT(IN) :: psi(NGrid,NBands) 
  REAL(DP), INTENT(IN) :: den(dfftt%nnr), grad_den(3, dfftt%nnr) 
  REAL(DP), INTENT(IN) :: ThrDen, ThrGrd 

  INTEGER :: ir, ir_end, INFO, lwork
  integer :: n, npt, ncpu_start
  REAL(DP) :: grad
  REAL(DP), ALLOCATABLE :: small(:,:), tau(:), work(:)

#if defined (__MPI)
  ir_end = dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
#else
  ir_end = dfftt%nnr
#endif

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE loc_scdm 
!-----------------------------------------------------------------------
