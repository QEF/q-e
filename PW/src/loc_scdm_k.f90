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
  !
  ! Variables and subroutines for localizing molecular orbitals based
  ! on a modified SCDM approach. Original SCDM method:
  ! A. Damle, L. Lin, L. Ying: J. Comput. Phys., 334 (2017) 1-5 
  ! Ivan Carnimeo: SCDM has been implemented with a parallel prescreening algorithm 
  !                in order to save CPU time and memory for large grids 
  ! 
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE exx,                  ONLY : dfftt, exxbuff, exxmat
  USE exx_base,             ONLY : nkqs

  IMPLICIT NONE
  SAVE
  REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Three=2.0d0
  integer :: n_scdm = 1

 CONTAINS
  !
  !------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE localize_orbitals_k( )
  !
  ! Driver for SCDM_k orbital localization 
  !         SCDM_PGG_k ... perform SCDM localization using the parallel
  !                        prescreening algorithm
  !         measure_localization_k ... compute the absolute overlap integrals 
  !                                    required by vexx_loc_k in exx.f90 
  !                                    and put them in exxmat 
  ! NOTE: localize_orbitals_k does not have iscdm/nscdm because alignment is not 
  !       implemented for K-Points
  !
  USE noncollin_module,     ONLY : npol
  USE exx,                  ONLY : x_occupation
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : nbnd
  USE mp_bands,             ONLY : me_bgrp
  USE control_flags,        ONLY : lmd
  !   
  implicit none
  integer :: NGrid, ikq, jk, k, NBands, nnk
  character(len=1) :: HowTo
  real(DP) :: Tot1, Tot2, Ave1, tot, ave 

  IF( n_scdm.ne.1 ) Call errore('localize_orbitals_k','nscdm for K-points NYI.',1)
  IF( lmd ) Call errore('localize_orbitals_k','localization with K-points not tested.',1)

  NGrid = dfftt%nnr * npol
  HowTo = 'G'  ! How to compute the absolute overlap integrals (only G for k)

  exxmat = One ! initialized to one because the absolute overlaps between
               ! occupied and virtual orbitals are assumed always large 
               ! and all ov integrals will be computed by vexx_loc_k

  NBands = int(sum(x_occupation(:,1))) 

  write(stdout,*) ' '
  write(stdout,*) 'NBands = ', NBands, ' nks = ', nks, ' nkqs = ', nkqs
  write(stdout,'(5X,A)') 'Canonical Orbitals '
  Tot1 = Zero 
  Ave1 = Zero 
  Tot2 = Zero 
  nnk = 0 
  DO ikq= 1, nkqs
    CALL measure_localization_k(NBands, ikq, tot, ave)
    Tot1 = Tot1 + tot
    Ave1 = Ave1 + ave
    DO jk = 1, nks 
      nnk = nnk + 1
      Call AbsOvG_k(NBands, ikq, jk, tot, ave) 
      Tot2 = Tot2 + tot
    END DO 
  END DO 
  Ave1 = Ave1 / dble(nkqs)
  write(stdout,'(7X,A,f24.6)') 'Total AbsOv          =', Tot2 
  write(stdout,'(7X,A,f24.6)') 'Aver. AbsOv          =', Tot2/dble(nnk) 
  write(stdout,'(7X,A,f24.6)') 'Total Spread [A**2]  =', Tot1 
  write(stdout,'(7X,A,f24.6)') 'Aver. Spread [A**2]  =', Ave1 

! if(me_bgrp.eq.0) then 
!   Call AbsOv_histogram_k(NBands,'hist1.dat')
! end if 

  write(stdout,'(5X,A)') 'SCDM-PGG_k localization'
  DO ikq= 1, nkqs
    CALL SCDM_PGG_k(NGrid, NBands, ikq)
  END DO 
  write(stdout,'(7X,A)') 'SCDM-PGG_k done ' 

  write(stdout,'(5X,A)') 'Localized Orbitals '
  Tot1 = Zero 
  Ave1 = Zero 
  Tot2 = Zero 
  nnk = 0 
  DO ikq= 1, nkqs
    CALL measure_localization_k(NBands, ikq, tot, ave)
    Tot1 = Tot1 + tot
    Ave1 = Ave1 + ave
    DO jk = 1, nks 
      nnk = nnk + 1
      Call AbsOvG_k(NBands, ikq, jk, tot, ave) 
      Tot2 = Tot2 + tot
    END DO 
  END DO 
  Ave1 = Ave1 / dble(nkqs)
  write(stdout,'(7X,A,f24.6)') 'Total AbsOv         =', Tot2 
  write(stdout,'(7X,A,f24.6)') 'Aver. AbsOv         =', Tot2/dble(nnk) 
  write(stdout,'(7X,A,f24.6)') 'Total Spread [A**2] =', Tot1 
  write(stdout,'(7X,A,f24.6)') 'Aver. Spread [A**2] =', Ave1 

! if(me_bgrp.eq.0) then 
!   Call AbsOv_histogram_k(NBands,'hist2.dat')
! end if 

END SUBROUTINE localize_orbitals_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE measure_localization_k(NBands, IK, TotSpread, AveSpread) 
USE noncollin_module,  ONLY : npol
USE cell_base,         ONLY : alat, omega, at, bg
USE exx,               ONLY : compute_density_k
USE constants,         ONLY : bohr_radius_angs 
implicit none
!
! Compute absolute overlap integrals
!
  REAL(DP), INTENT(OUT) :: TotSpread, AveSpread
  INTEGER :: NBands, ibnd, IK
  REAL(DP) :: RJunk(4), SpreadPBC(3)

  call start_clock('measure')

  TotSpread = Zero   
  AveSpread = Zero 
  DO ibnd = 1, NBands
    call compute_density_k(.false.,.false.,RJunk(1:3), SpreadPBC, RJunk(4), &
                           exxbuff(1,ibnd,IK), exxbuff(1,ibnd,IK), &
                           dfftt%nnr*npol, ibnd, ibnd)
    TotSpread = TotSpread + SpreadPBC(1) + SpreadPBC(2) + SpreadPBC(3) 
  ENDDO
  TotSpread = TotSpread * bohr_radius_angs**2
  AveSpread = TotSpread / dble(NBands) 

! write(stdout,'(7X,A,I3)')  'IK = ', IK
! write(stdout,'(7X,A,f12.6,I3)') 'Total Spread [A**2]   =', TotSpread 
! write(stdout,'(7X,A,f12.6,I3)') 'Aver. Spread [A**2]   =', AveSpread 

  call stop_clock('measure')

END SUBROUTINE measure_localization_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AbsOvG_k(NBands, IKQ, JK, loc_diag, loc_off) 
USE noncollin_module,  ONLY : npol
USE fft_interfaces,    ONLY : fwfft
USE wvfct,             ONLY : npwx
USE exx_band,          ONLY : igk_exx
USE klist,             ONLY : nks, ngk
implicit none
!
! Compute the Absolute Overlap in G-space 
! (cutoff might not be accurate for the moduli of the wavefunctions)
!
  INTEGER :: NBands, ibnd, jbnd, ig, kk, IKQ, JK, npw
  REAL(DP) :: tmp
  REAL(DP), INTENT(OUT) :: loc_diag, loc_off
  COMPLEX(DP), ALLOCATABLE :: buffer(:), GorbtI(:,:), GorbtJ(:,:)
  INTEGER, EXTERNAL  :: global_kpoint_index
  COMPLEX(DP), ALLOCATABLE :: Mat(:,:)

  call start_clock('measure')

! write(stdout,'(5X,A)') ' ' 
! write(stdout,'(5X,A)') 'Absolute Overlap calculated in G-space (K)'

  kk = global_kpoint_index ( nks, JK )
  ALLOCATE( Mat(NBands,NBands) )
  allocate( buffer(dfftt%nnr * npol), GorbtI(npwx,NBands), GorbtJ(npwx,NBands) )
  GorbtJ = (Zero,Zero) 
  GorbtI = (Zero,Zero) 
  npw = ngk (JK)
  DO jbnd = 1, NBands 
    buffer(:) = abs(exxbuff(:,jbnd,JK)) 
!   buffer(:) = exxbuff(:,jbnd,JK)
    CALL fwfft( 'Wave' , buffer, dfftt )
    DO ig = 1, npw
      GorbtJ(ig,jbnd) = buffer(dfftt%nl(igk_exx(ig,kk))) 
    ENDDO
    buffer(:) = abs(exxbuff(:,jbnd,IKQ)) 
!   buffer(:) = exxbuff(:,jbnd,IKQ)
    CALL fwfft( 'Wave' , buffer, dfftt )
    DO ig = 1, npw
      GorbtI(ig,jbnd) = buffer(dfftt%nl(igk_exx(ig,kk)))  
    ENDDO
  ENDDO
! IF(IKQ.eq.JK) THEN
!   CALL matcalc_k('AbsOv-',.false., 2, 0, npwx, NBands, NBands, GorbtI, GorbtJ, Mat, tmp)
! ELSE
    CALL matcalc_k('AbsOv-',.false., 0, 0, npwx, NBands, NBands, GorbtI, GorbtJ, Mat, tmp)
! END IF 
  loc_diag = Zero 
  loc_off  = Zero 
  DO ibnd = 1, NBands 
    loc_diag = loc_diag + dble(Mat(ibnd,ibnd))
    DO jbnd = 1, ibnd-1 
      loc_off = loc_off + dble(Mat(ibnd,jbnd)) + dble(Mat(jbnd,ibnd))
      exxmat(ibnd, IKQ, jbnd, JK) = dble(Mat(ibnd,jbnd))
      exxmat(jbnd, IKQ, ibnd, JK) = dble(Mat(jbnd,ibnd))
    END DO 
  END DO 
  write(stdout,'(7X,5(A,I3),2(A,f12.6))')  'IKQ = ', IKQ, & 
        '  JK = ', JK, '  kk = ', kk, ' NBands = ', NBands, ' size = ', size(exxbuff,3), &
        '  Total Charge =', loc_diag, '  Total Abs. Overlap =', loc_off 
  deallocate ( buffer, GorbtI, GorbtJ )

  DEALLOCATE( Mat ) 
  call stop_clock('measure')

END SUBROUTINE AbsOvG_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AbsOv_histogram_k(n,filename) 
  USE klist,                ONLY : nks
implicit none
!
!  Print the distribution of absolute overlap integrals (as an histogram)
!  in the range [0,1] with NHist intervals
!  WARNING: HUGE amount of integrals, use only for VERY small systems with a few
!           K-Points
!
  character(len=*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: n
  INTEGER :: i,j,k, iik, jjk, NHist
  INTEGER,  ALLOCATABLE :: histogram(:)
  REAL(DP), ALLOCATABLE :: XHist(:)
  REAL(DP) :: xstart, xstep, integral
  integer :: io_histogram
  integer, external :: find_free_unit

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
    do iik = 1, nkqs 
      do j = 1, n
        do jjk = 1, nks
          integral = exxmat(i, iik, j, jjk)
          if(integral.lt.Zero) then 
            Call errore('AbsOv_histogram_k','Abs. Ov. < 0 found.',1)
          else
            do k = 1, NHist 
              IF(integral.ge.(XHist(k)-xstart).and.integral.lt.(XHist(k)+xstart)) &
                   histogram(k)=histogram(k)+1
            end do 
          end if 
        end do 
      end do 
    end do 
  end do 
  io_histogram = find_free_unit()
  open(io_histogram,file=filename,status='unknown')
  do k = 1, NHist
    write(io_histogram,'(f12.6,2I10)') XHist(k), histogram(k)
  end do 
  close(io_histogram)
  DEALLOCATE(histogram,XHist)

END SUBROUTINE AbsOv_histogram_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCDM_PGG_k(NGrid, NBands, IKQ)
USE mp_bands,          ONLY : nproc_bgrp
USE loc_scdm,          ONLY : scdm_thresholds, scdm_points, scdm_prescreening
!
! density matrix localization (I/O in psi) 
! (K-Points version)
!
IMPLICIT NONE
  INTEGER,  INTENT(IN)    :: NGrid, NBands, IKQ
  COMPLEX(DP), ALLOCATABLE :: QRbuff(:,:), mat(:,:)
  INTEGER,  ALLOCATABLE :: pivot(:), list(:), cpu_npt(:)
  REAL(DP), ALLOCATABLE :: den(:), grad_den(:,:)
  INTEGER  :: nptot
  REAL(DP) :: ThrDen, ThrGrd

  call start_clock('localization')
! write(stdout,'(5X,A)') 'SCDM localization with prescreening'
  allocate( den(dfftt%nnr), grad_den(3, dfftt%nnr) )
  Call scdm_thresholds( den, grad_den, ThrDen, ThrGrd )
  allocate( cpu_npt(0:nproc_bgrp-1) )
  Call scdm_points( den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot )
  allocate( list(nptot), pivot(nptot) )
  Call scdm_prescreening_k ( NGrid, NBands, exxbuff(1,1,IKQ), den, grad_den, ThrDen, ThrGrd, &
                                                     cpu_npt, nptot, list, pivot )
  deallocate( den, grad_den ) 
! Psi(pivot(1:NBands),:) in mat
  allocate( mat(NBands,NBands) )
  mat = (Zero, Zero)
  Call scdm_fill_k( .true., nptot, NGrid, NBands, cpu_npt, pivot, list, exxbuff(1,1,IKQ), Mat)
! Pc = Psi * Psi(pivot(1:NBands),:)' in QRbuff
  allocate( QRbuff(NGrid, NBands) )
  QRbuff = (Zero, Zero)
  CALL ZGEMM( 'N' , 'N' , NGrid, NBands, NBands, (One,Zero), exxbuff(1,1,IKQ), NGrid, mat, NBands, (Zero,Zero), QRbuff, NGrid) 
! Orthonormalization
! Pc(pivot(1:NBands),:) in mat 
  mat = (Zero, Zero)
  Call scdm_fill_k( .false., nptot, NGrid, NBands, cpu_npt, pivot, list, QRBuff, mat)
  deallocate( cpu_npt )
! Cholesky(psi)^(-1) in mat 
  CALL invchol_k(NBands,mat)
! Phi = Pc * Chol^(-1) = QRbuff * mat
  CALL ZGEMM( 'N' , 'T' , NGrid, NBands, NBands, (One,Zero), QRbuff, NGrid, mat, NBands, (Zero,Zero), exxbuff(1,1,IKQ), NGrid) 
  deallocate( QRbuff, mat, pivot, list )
! write(stdout,'(7X,A)') 'SCDM-PGG_k done ' 
  call stop_clock('localization')

END SUBROUTINE SCDM_PGG_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE scdm_prescreening_k ( NGrid, NBands, psi, den, grad_den, ThrDen, ThrGrd, cpu_npt, nptot, &
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
  COMPLEX(DP), INTENT(IN) :: psi(NGrid,NBands) 
  REAL(DP), INTENT(IN) :: den(dfftt%nnr), grad_den(3, dfftt%nnr) 
  REAL(DP), INTENT(IN) :: ThrDen, ThrGrd 

  INTEGER :: ir, ir_end, INFO, lwork
  integer :: n, npt, ncpu_start
  REAL(DP) :: grad
  COMPLEX(DP), ALLOCATABLE :: small(:,:), tau(:), work(:)
  REAL(DP), ALLOCATABLE :: rwork(:)

#if defined (__MPI)
  ir_end = dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
#else
  ir_end = dfftt%nnr
#endif

! find the map of the indeces
  allocate( small(NBands,nptot) )
  small = (Zero, Zero)
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
  allocate( tau(nptot), work(lwork), rwork(2*nptot) )
  tau = (Zero, Zero)
  work = (Zero, Zero)
  pivot = 0
  INFO = -1 
  CALL ZGEQP3( NBands, nptot, small, NBands, pivot, tau, work, lwork, rwork, INFO )
  deallocate( tau, work, rwork, small )

END SUBROUTINE scdm_prescreening_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE scdm_fill_k( op, nptot, NGrid, NBands, CPUPts, Pivot, List, Vect, Mat)
USE mp,                ONLY : mp_sum
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
!
! Fill the matrix Mat with the elements of Vect
! mapped by CPUPts, Pivot and List 
!
IMPLICIT NONE
  INTEGER,  INTENT(IN)  :: NBands, NGrid, nptot
  INTEGER,  INTENT(IN)  :: CPUPts(0:nproc_bgrp-1), Pivot(nptot), List(nptot)
  COMPLEX(DP), INTENT(IN)  :: Vect(NGrid, NBands)
  COMPLEX(DP), INTENT(OUT) :: Mat(NBands,NBands)
  INTEGER :: i, NStart, NEnd
  LOGICAL :: op

  Mat = (Zero, Zero)
  do i = 1, NBands
    NStart = sum(CPUPts(0:me_bgrp-1))
    NEnd   = sum(CPUPts(0:me_bgrp))
    if(Pivot(i).le.NEnd.and.Pivot(i).ge.NStart+1) THEN
      IF(op) THEN
        Mat(:,i) = conjg(Vect(List(pivot(i)),:))
      ELSE 
        Mat(:,i) = Vect(List(pivot(i)),:)
      END IF
    END IF 
  end do 
  call mp_sum(Mat,intra_bgrp_comm)

END SUBROUTINE scdm_fill_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE loc_scdm_k
!-----------------------------------------------------------------------
