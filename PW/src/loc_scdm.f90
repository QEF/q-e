! Copyright (C) 2005-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE loc_scdm 
  !--------------------------------------
  !
  ! Variables and subroutines related to the calculation of 
  ! SCDM localization of molecular orbitals 
  ! Implements ACE: Lin Lin, J. Chem. Theory Comput. 2016, 12, 2242
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE exx,                  ONLY : exx_fft, x_nbnd_occ, locbuff, nkqs

  IMPLICIT NONE
  SAVE
  LOGICAL ::  use_scdm=.false. ! if .true. enable Lin Lin's SCDM localization
                               ! currently implemented only within ACE formalism
  LOGICAL ::  scdm_dipole=.false.
  REAL(DP) :: scdm_den, scdm_grd 

 CONTAINS
  !
  !------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCDM_PGG(psi, NQR, nbnd_eff)
USE cell_base,         ONLY : omega
USE mp,                ONLY : mp_stop, mp_barrier, mp_sum
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
USE vdW_DF,            ONLY : numerical_gradient
USE fft_base,          ONLY : dfftp
USE scf,               ONLY : rho
USE lsda_mod,          ONLY : nspin
!
! density matrix decomposition (I/O in psi)
!
IMPLICIT NONE
  integer :: NQR, nbnd_eff, lwork, INFO, i, j, n, ir, jbnd, ir_end, nnr
  real(DP), allocatable :: QRbuff(:,:), tau(:), work(:), mat(:,:),mat2(:,:)
  integer, allocatable :: pivot(:)
  real(DP) :: charge, grad, psi(NQR,nbnd_eff)
  integer :: npt, nptot, icpu, ncpu_start, ncpu_end, nxxs
  integer, allocatable :: list(:), cpu_npt(:)
  real(DP), allocatable :: small(:,:)
  real(DP), allocatable :: den(:), grad_den(:,:)
! real(DP), parameter :: ThDen = 0.10d0, ThGrad=0.200d0 

  call start_clock('localization')
  write(stdout,'(A)') '--------------------------------'
  write(stdout,'(A)') 'Gradient-based SCDM localization'
  write(stdout,'(A)') '--------------------------------'
  write(stdout,'(2(A,f12.6))') '   scdm_den = ', scdm_den, ' scdm_grd = ',scdm_grd

  if(exx_fft%dfftt%nnr.ne.dfftp%nnr) then 
    write(stdout,*) exx_fft%dfftt%nnr, dfftp%nnr
    call errore( 'SCDM_PGG', 'density and orbital grids not identical',1)
  end if 

  nnr = dfftp%nnr

#if defined (__MPI)
!  ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
  ir_end = dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
#else
  ir_end = nnr
#endif

  nxxs = exx_fft%dfftt%nr1x *exx_fft%dfftt%nr2x *exx_fft%dfftt%nr3x

  allocate( den(nnr), grad_den(nnr, 3) )
  charge = 0.0d0
  den(:) = rho%of_r(:,1)
  IF ( nspin == 2 ) den(:) = den(:) + rho%of_r(:,2) 
  do ir = 1, ir_end 
    charge = charge + den(ir) * omega / float(nxxs)
  end do 
  call mp_sum(charge,intra_bgrp_comm)
  write(stdout,'(A,f12.6)') '    charge = ', charge

! find numerical gradient of the density
  call numerical_gradient (den, grad_den)

  charge = 0.0d0
  do ir = 1, ir_end 
    grad = sqrt( grad_den(ir,1)**2 +  grad_den(ir,2)**2 +  grad_den(ir,3)**2 )
    charge = charge + grad * omega / float(nxxs)
  end do 
  call mp_sum(charge,intra_bgrp_comm)
  write(stdout,'(A,f12.6)') '    grad   = ', charge

! find the relevant point for the allocation
  allocate( cpu_npt(0:nproc_bgrp-1) )
  npt = 0
  cpu_npt(:) = 0
  do ir = 1, ir_end 
    grad = 1.0d0
    if(den(ir).gt.scdm_den) then 
      grad = sqrt( grad_den(ir,1)**2 +  grad_den(ir,2)**2 +  grad_den(ir,3)**2 ) 
      if(grad.lt.scdm_grd) then    
        npt = npt + 1
      end if 
    end if 
  end do 
  nptot = npt
  cpu_npt(me_bgrp) = npt
  call mp_sum(nptot,intra_bgrp_comm)
  if(nptot.le.0) call errore('SCDM_PGG', 'No points prescreened. Loose the thresholds', 1) 
  call mp_sum(cpu_npt,intra_bgrp_comm)
  write(*,'(2(A,I8))') '    npt = ', npt, '  procID= ', me_bgrp 
  write(stdout,*) '    reduced matrix, allocate: ', nptot
  write(stdout,*) '    cpu_npt = ', cpu_npt(:)


! find the map of the index 
  allocate( small(nbnd_eff,nptot), list(nptot) )
  small = 0.0d0
  list = 0
  n = 0
  do ir = 1, ir_end
    grad = 1.0d0
    if(den(ir).gt.scdm_den) then 
      grad = sqrt( grad_den(ir,1)**2 +  grad_den(ir,2)**2 +  grad_den(ir,3)**2 ) 
      if(grad.lt.scdm_grd) then    
        n = n + 1
        ncpu_start = sum(cpu_npt(0:me_bgrp-1))
        icpu = ncpu_start+n
        small(:,icpu) = psi(ir,:)
        list(icpu) = ir
      end if 
    end if 
  end do 
! call matprt('small',nbnd_eff,nptot,small)
  call mp_sum(small,intra_bgrp_comm)
! call matprt('small',nbnd_eff,nptot,small)
  call mp_sum(list,intra_bgrp_comm)

  lwork = 4*nptot
  allocate( pivot(nptot), tau(nptot), work(lwork) )
  tau = 0.0d0
  work = 0.0d0  
  pivot = 0
  INFO = -1 
  CALL DGEQP3( nbnd_eff, nptot, small, nbnd_eff, pivot, tau, work, lwork, INFO )
  do i = 1, nbnd_eff
    j = list(pivot(i))
    grad = sqrt( grad_den(j,1)**2 +  grad_den(j,2)**2 +  grad_den(j,3)**2 )
    ncpu_start = 0
    ncpu_end   = 0
    ncpu_start = sum(cpu_npt(0:me_bgrp-1))
    ncpu_end   = sum(cpu_npt(0:me_bgrp))
!   if(pivot(i).le.ncpu_end.and.pivot(i).ge.ncpu_start+1) &
!             write(*,'(A,3I9,2f16.8)') 'pivoting: ', me_bgrp, i, pivot(i), den(j), grad 
  end do 
  deallocate( den, grad_den) 
  deallocate( tau, work )
  deallocate( small )

! Psi(pivot(1:nbnd_eff),:) in mat
  lwork = 3*nbnd_eff
  allocate( mat(nbnd_eff,nbnd_eff), mat2(nbnd_eff,nbnd_eff), tau(nbnd_eff), work(lwork) )
  mat = 0.0d0
  do i = 1, nbnd_eff
!   mat(:,i) = QRbuff(:,list(pivot(i)))
    ncpu_start = 0
    ncpu_end   = 0
    ncpu_start = sum(cpu_npt(0:me_bgrp-1))
    ncpu_end   = sum(cpu_npt(0:me_bgrp))
    if(pivot(i).le.ncpu_end.and.pivot(i).ge.ncpu_start+1) mat(:,i) = psi(list(pivot(i)),:)
  end do 
! call matprt('Q',nbnd_eff,nbnd_eff,mat)
  call mp_sum(mat,intra_bgrp_comm)
! call matprt('Q',nbnd_eff,nbnd_eff,mat)
! deallocate( QRbuff ) 
  deallocate( mat2, tau, work )

! Pc = Psi * Psi(pivot(1:nbnd_eff),:)' in QRbuff
  allocate( QRbuff(NQR, nbnd_eff) )
  QRbuff = 0.0d0
  CALL DGEMM( 'N' , 'N' , NQR, nbnd_eff, nbnd_eff, 1.0d0, psi, NQR, mat, nbnd_eff, 0.0d0, QRbuff, NQR) 

! Orthonormalization
! Pc(pivot(1:nbnd_eff),:) in mat 
  mat = 0.0d0
  do i = 1, nbnd_eff
    ncpu_start = 0
    ncpu_end   = 0
    ncpu_start = sum(cpu_npt(0:me_bgrp-1))
    ncpu_end   = sum(cpu_npt(0:me_bgrp))
    if(pivot(i).le.ncpu_end.and.pivot(i).ge.ncpu_start+1) mat(i,:) = QRBuff(list(pivot(i)),:) 
!   mat(i,:) = QRbuff(list(pivot(i)),:)
  end do 
  deallocate( cpu_npt )
! call matprt('Q2',nbnd_eff,nbnd_eff,mat)
  call mp_sum(mat,intra_bgrp_comm)
! call matprt('Q2',nbnd_eff,nbnd_eff,mat)

! Cholesky(psi)^(-1) in mat 
  CALL invchol(nbnd_eff,mat)
  allocate( mat2(nbnd_eff,nbnd_eff) )
  mat2 = 0.0d0
  do i = 1, nbnd_eff
    do j = i, nbnd_eff
      mat2(i,j) = mat(j,i)
    end do 
  end do 
  deallocate( mat )
! Phi = Pc * Chol^(-1) = QRbuff * mat
  psi = 0.0d0
  CALL DGEMM( 'N' , 'N' , NQR, nbnd_eff, nbnd_eff, 1.0d0, QRbuff, NQR, mat2, nbnd_eff, 0.0d0, psi, NQR) 
  deallocate( QRbuff, mat2, pivot, list )
  write(stdout,'(A)') '    SCDM-PGG done ' 

  call stop_clock('localization')

END SUBROUTINE SCDM_PGG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invchol(n,A)
USE exx, ONLY : errinfo
IMPLICIT NONE
!
! given a matrix A, returns the inverse of the Cholesky decomposition of A
! for real matrices
!
  real(DP) :: A(n,n)
  integer :: n, INFO

  INFO = -1
  CALL DPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('DPOTRF','Cholesky failed in invchol.',INFO)
  INFO = -1
  CALL DTRTRI( 'L', 'N', n, A, n, INFO )
  CALL errinfo('DTRTRI','inversion failed in invchol.',INFO)

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE localize_orbitals( )
  !
  ! Driver for SCDM orbital localization 
  !
  USE noncollin_module,  ONLY : npol
  USE exx,               ONLY : locbuff, measure_localization, measure_localization_G
  !   
  implicit none
  integer :: NQR
 
  NQR = exx_fft%dfftt%nnr * npol

! CALL measure_localization(locbuff(1,1,1), NQR, x_nbnd_occ)
  CALL measure_localization_G(locbuff(1,1,1), NQR, x_nbnd_occ)
  CALL SCDM_PGG(locbuff(1,1,1), NQR, x_nbnd_occ)
! CALL measure_localization(locbuff(1,1,1), NQR, x_nbnd_occ)
  CALL measure_localization_G(locbuff(1,1,1), NQR, x_nbnd_occ)

END SUBROUTINE localize_orbitals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE loc_scdm 
!-----------------------------------------------------------------------
