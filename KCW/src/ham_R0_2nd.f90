! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO CMPLX(0.D0,0.D0,kind=DP)
#define ONE CMPLX(0.D0,1.D0, kind=DP)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE ham_R0_2nd ()
  !---------------------------------------------------------------------
  !
  !! This routine compute the KI correction to second order to the the 
  !! Hamiltonian in real space in the basis if Wannier Fuction. 
  !! ONLY Diagonal CORRECTIONS. Obsolete, use hamilt.f90 instead. 
  !
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE klist,                ONLY : nkstot
  USE mp,                   ONLY : mp_sum
  USE control_kcw,          ONLY : kcw_iverbosity, spin_component, num_wann, iorb_start, x_q, &
                                   iorb_end, alpha_final, num_wann_occ, irvect, seedname
  USE buffers,              ONLY : get_buffer, save_buffer
  USE io_global,            ONLY : stdout, ionode
  USE control_kcw,          ONLY : iurho_wann
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE fft_interfaces,       ONLY : fwfft
  USE control_lr,           ONLY : lgamma
  USE lsda_mod,             ONLY : nspin
  USE gvecs,                ONLY : ngms
  USE solve_linter_koop_mod 
  USE qpoint,               ONLY : xq
  !
  !USE mp_world,             ONLY : mpime
  !
  USE cell_base,            ONLY : omega
  ! 
  !
  IMPLICIT NONE
  ! 
  INTEGER :: iq, nqs
  ! ... Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
  !
  INTEGER :: iwann, jwann, lrrho
  ! ... Band counters, leght of the rho record
  !
  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
  ! The periodic part of the wannier orbital density and potential 
  !
  COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), drhog_scf(:,:), delta_vg_(:,:) 
  !
  COMPLEX(DP) :: pi_q_unrelax, pi_q_relax, pi_q_relax_rs
  !
  COMPLEX(DP), ALLOCATABLE :: deltaHq(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: deltaHR(:,:,:)
  !
  !
  COMPLEX(DP) :: drho_zero
  !
  REAL(DP) :: weight(nkstot)
  !
  CHARACTER (len=256) :: hocc_file, hemp_file
  INTEGER :: ierr, ieff, jeff
  !
  CHARACTER(LEN=9)  :: cdate, ctime
  CHARACTER(LEN=33) :: header
  INTEGER :: i
  !
  nqs = nkstot/nspin
  !
  ALLOCATE( deltaHq(num_wann,num_wann,nqs) )
  ALLOCATE( deltaHR(num_wann,num_wann,nqs) )
  !
#ifdef DEBUG
  write(*,'(/,"DEBUG: The list of G  vectors")')
  do ig = 1, 10
     xq_ = g(:,ig)
     CALL cryst_to_cart(1, xq_, at, -1)
     write(*,'("i = ", i3, 3x, "G(i) = ", 3x, 3f8.4, "  [Cryst]",3x, "|G(i)| = ", f12.6 )') & 
                       ig, (xq_(iq), iq=1,3),  sqrt(sum (xq_(:)*xq_(:)))
  enddo
#endif
  !
  pi_q_relax = ZERO
  pi_q_relax_rs = ZERO
  pi_q_unrelax = ZERO
  !
  drho_zero = ZERO
  deltaHq = ZERO
  deltaHR = ZERO
  !
  WRITE( stdout, '(5X,"INFO: KI[2nd, (R=0,i=j)] CALCULATION ...")')
  !
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    xq(1:3)  = x_q(1:3,iq)
    !
    !
    lrrho=num_wann*dffts%nnr
    CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    !! Retrive the rho_wann_q(r) from buffer in REAL space
    IF (kcw_iverbosity .gt. 1 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED")') 
    !
    ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), drhog_scf (ngms, nspin), delta_vg_(ngms,nspin) )
    !
    IF ( lgamma ) CALL check_density (rhowann) 
    !! CHECK: For q==0 the sum over k and v should give the density. If not something wrong...
    !
    weight(iq) = 1.D0/nqs ! No SYMM 
    !
    !WRITE(stdout, '("weight =", i5, f12.8)') iq, weight(iq)
    !
    !DO iwann = 1, num_wann  ! for each band, that is actually the perturbation
    DO iwann = iorb_start, iorb_end 
       !
       drhog_scf (:,:) = ZERO
       rhog(:)         = ZERO
       delta_vg(:,:)   = ZERO
       vh_rhog(:)      = ZERO
       rhor(:)         =ZERO
       !
       rhor(:) = rhowann(:,iwann)
       !! The periodic part of the orbital desity in real space
       !
       CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ ) 
       !! The periodic part of the perturbation DeltaV_q(G)
       ! 
       deltaHq(iwann, iwann, iq) = sum (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
       !
       !
    ENDDO
    ! 
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, drhog_scf, delta_vg_ )
    !
    !    
  ENDDO ! qpoints
  !
  CALL mp_sum (deltaHq, intra_bgrp_comm)
  ! ... Sum over different processes (G vectors) 
  !
  ! At this point we have sum_G (rho^jwann_q(G) DeltaV^iwann_q(G)) stored in deltaHq
  !
  deltaHR(:,:,:) = CMPLX(0.D0, 0.D0,kind=DP)
  ! ... initialize deltaH(R)
  !
  ! First Build up the H_nm(R=0)
  ! This is the only contribution for occupied state (in the KI approximation),
  ! and the leading one for empty
  !
  DO iwann = 1, num_wann
    !
    DO iq = 1, nqs
      !
      ! ... only the diagonal n=m are different from zero when R=0
      deltaHR(iwann,iwann,1) = deltaHR(iwann,iwann,1) + deltaHq(iwann,iwann,iq)
    ENDDO
    !
    IF (iwann .le. num_wann_occ) THEN 
       deltaHR(iwann,iwann,1) = -0.5D0 * deltaHR(iwann,iwann,1) * alpha_final(iwann)
    ELSE
       deltaHR(iwann,iwann,1) = +0.5D0 * deltaHR(iwann,iwann,1) * alpha_final(iwann)
    ENDIF
    !
  ENDDO
  !
  hocc_file=TRIM( seedname ) // '_KC_ham_occ.dat'
  hemp_file=TRIM( seedname ) // '_KC_ham_emp.dat'
  !
  ! ... Write the result on file 
  IF (ionode ) THEN 
    !
    CALL date_and_tim( cdate, ctime )
    header = 'Written on '//cdate//' at '//ctime
    !
    ! ... Occupied states
    OPEN (UNIT = 1001, FILE = hocc_file, FORM = 'formatted', STATUS = 'unknown', IOSTAT=ierr )
    !
    WRITE(1001, * ) header
    WRITE(1001,'(i5)') num_wann_occ
    ! ... number of wannier 
    WRITE(1001,'(i5)') nqs
    WRITE(1001,'(15i5)') (1, i=1,nqs)
    ! ... number of R points (identical to q points) 
    DO iq = 1, nqs
      DO iwann = 1, num_wann_occ
        DO jwann = 1, num_wann_occ
         WRITE(1001,'(3i4,3x,2i4,2f15.8)') irvect(:,iq), iwann, jwann, & 
                                           REAL(deltaHR(iwann,jwann,iq)), AIMAG(deltaHR(iwann,jwann,iq))
        ENDDO
      ENDDO
    ENDDO
    CLOSE (1001)
    !
    CALL date_and_tim( cdate, ctime )
    header = 'Written on '//cdate//' at '//ctime
    !
    ! ... Empty states
    OPEN (UNIT = 1001, FILE = hemp_file, FORM = 'formatted', STATUS = 'unknown', IOSTAT=ierr )
    !
    WRITE(1001, * ) header
    WRITE(1001,'(i5)') num_wann-num_wann_occ
    WRITE(1001,'(i5)') nqs
    WRITE(1001,'(15i5)') (1, i=1,nqs)
    DO iq = 1, nqs
      DO iwann = 1, num_wann-num_wann_occ
        ieff = iwann + num_wann_occ 
        DO jwann = 1, num_wann-num_wann_occ
          jeff = jwann + num_wann_occ 
         WRITE(1001,'(3i4,3x,2i4,2f15.8)') irvect(:,iq), iwann, jwann, &
                                           REAL(deltaHR(ieff,jeff,iq)), AIMAG(deltaHR(ieff,jeff,iq))
        ENDDO
      ENDDO
    ENDDO
    CLOSE (1001)
    !
  ENDIF
  !
  WRITE(stdout,'(/,5x, "The diagonal H(R=0)")') 
  WRITE(stdout,'("     iwann   alpha       Real(H)       Im(H)")') 
  DO iwann = 1, num_wann
     WRITE(stdout,'(5x, i3, 3x, f8.4, 3x, 2f12.8)') iwann, alpha_final(iwann), &
                                                    REAL(deltaHR(iwann,iwann,1)), AIMAG(deltaHR(iwann,iwann,1))
  ENDDO
  !
  WRITE( stdout, '(5X,"INFO: KI[2nd, (R=0,i=j)] CALCULATION ... DONE")')
  !
END subroutine 

