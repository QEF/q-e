!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_calc_dens( evc1, response_calc )
  !---------------------------------------------------------------------
  !
  ! Optical case
  ! This subroutine calculates the response charge density 
  ! from linear response orbitals and ground state orbitals.
  !
  ! Input : evc1 (qdash etc) 
  ! Output: rho = 2 * \sum_v ( revc0_v(r) . revc1_v(r,w) )
  !
  ! In case of US PP, becsum is also calculated here.
  ! In case of charge response calculation, the rho_tot is calculated here.
  !
  ! Modified by Osman Baris Malcioglu in 2009
  ! Modified by Simone Binnie in 2012 
  !
  USE ions_base,              ONLY : ityp, nat, ntyp=>nsp
  USE cell_base,              ONLY : omega
  USE ener,                   ONLY : ef
  USE gvecs,                  ONLY : nls, nlsm, doublegrid
  USE fft_base,               ONLY : dffts, dfftp
  USE fft_interfaces,         ONLY : invfft
  USE io_global,              ONLY : stdout
  USE kinds,                  ONLY : dp
  USE klist,                  ONLY : nks, xk, wk, ngk, igk_k
  USE lr_variables,           ONLY : evc0,revc0,rho_1,lr_verbosity,&
                                     & charge_response, itermax,&
                                     & cube_save, LR_iteration,&
                                     & LR_polarization, project,&
                                     & evc0_virt, F,nbnd_total,&
                                     & n_ipol, becp1_virt, rho_1c,&
                                     & lr_exx
  USE lsda_mod,               ONLY : current_spin, isk
  USE wavefunctions_module,   ONLY : psic
  USE wvfct,                  ONLY : nbnd, et, wg, npwx
  USE control_flags,          ONLY : gamma_only
  USE uspp,                   ONLY : vkb, nkb, okvan, qq, becsum
  USE uspp_param,             ONLY : upf, nh
  USE io_global,              ONLY : ionode, stdout
  USE io_files,               ONLY : tmp_dir, prefix
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : inter_pool_comm, intra_bgrp_comm,&
                                     inter_bgrp_comm 
  USE realus,                 ONLY : addusdens_r
  USE charg_resp,             ONLY : w_T, lr_dump_rho_tot_cube,&
                                     & lr_dump_rho_tot_xyzd, &
                                     & lr_dump_rho_tot_xcrys,&
                                     & resonance_condition, epsil,&
                                     & rho_1_tot, rho_1_tot_im
  USE noncollin_module,       ONLY : nspin_mag, npol
  USE control_flags,          ONLY : tqr
  USE becmod,                 ONLY : becp
  USE lr_exx_kernel,          ONLY : lr_exx_kernel_int, revc_int,&
                                     & revc_int_c
  USE constants,              ONLY : eps12
  !
  IMPLICIT NONE
  !
  CHARACTER(len=6), EXTERNAL   :: int_to_char
  COMPLEX(kind=dp), INTENT(in) :: evc1(npwx*npol,nbnd,nks)
  LOGICAL, INTENT(in)          :: response_calc
  ! functions
  REAL(kind=dp) :: ddot
  !
  ! Local variables
  !
  INTEGER       :: ir, ik, ibnd, jbnd, ig, ijkb0, np, na
  INTEGER       :: ijh,ih,jh,ikb,jkb ,ispin 
  INTEGER       :: i, j, k, l
  REAL(kind=dp) :: w1, w2, scal, rho_sum
  ! These are temporary buffers for the response 
  REAL(kind=dp), ALLOCATABLE :: rho_sum_resp_x(:), rho_sum_resp_y(:),&
                              & rho_sum_resp_z(:)  
  CHARACTER(len=256) :: tempfile, filename
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_calc_dens>")')
  ENDIF
  !
  CALL start_clock('lr_calc_dens')
  !
  ALLOCATE( psic(dfftp%nnr) )
  psic(:)    = (0.0d0,0.0d0)
  !
  IF (gamma_only) THEN
     rho_1(:,:) =  0.0d0
  ELSE
     rho_1c(:,:) =  0.0d0
  ENDIF
  !
  IF (gamma_only) THEN
     !
     IF (lr_exx) revc_int = 0.0d0
     !
     CALL lr_calc_dens_gamma()
     !
     ! If a double grid is used, interpolate onto the fine grid
     !
     IF ( doublegrid ) CALL interpolate(rho_1,rho_1,1)
     !
#if defined(__MPI)
     CALL mp_sum(rho_1, inter_bgrp_comm)
#endif
     !
  ELSE
     !
     IF (lr_exx) revc_int_c = (0.0d0,0.d0)
     !
     CALL lr_calc_dens_k()
     !
     ! If a double grid is used, interpolate onto the fine grid
     !
     IF ( doublegrid ) CALL cinterpolate(rho_1c,rho_1c,1)
     !
  ENDIF
  !
  ! Here we add the Ultrasoft contribution to the charge density
  ! response. 
  !
  IF (okvan) THEN
     IF (tqr) THEN
        CALL addusdens_r(rho_1,.FALSE.)
     ELSE
        CALL addusdens(rho_1)
     ENDIF
  ENDIF
  !
  ! The psic workspace can present a memory bottleneck
  !
  DEALLOCATE ( psic )
  !
#if defined(__MPI)
  IF(gamma_only) THEN
     CALL mp_sum(rho_1, inter_pool_comm)
  ELSE
     CALL mp_sum(rho_1c, inter_pool_comm)
  ENDIF
#endif
  !
  ! check response charge density sums to 0
  !
  IF (lr_verbosity > 0) THEN
     ! 
     DO ispin = 1, nspin_mag
        !
        rho_sum = 0.0d0
        rho_sum = SUM(rho_1(:,ispin))
        !
#if defined(__MPI)
        CALL mp_sum(rho_sum, intra_bgrp_comm )
#endif
        !
        rho_sum = rho_sum * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
        !
        IF (ABS(rho_sum) > eps12) THEN
           !
           IF (tqr) THEN
              !
              WRITE(stdout,'(5X, "lr_calc_dens: Charge drift due to &
                   &real space implementation = " ,1X,e12.5)') rho_sum
              !
           ELSE
              !
              WRITE(stdout,'(5X,"lr_calc_dens: ****** response &
                   &charge density does not sum to zero")')
              !
              WRITE(stdout,'(5X,"lr_calc_dens: ****** response &
                   &charge density =",1X,e12.5)') rho_sum
              !
              WRITE(stdout,'(5X,"lr_calc_dens: ****** response &
                   &charge density, US part =",1X,e12.5)') scal
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
  ENDIF
  !
  IF (charge_response == 2 .AND. LR_iteration /=0) THEN
     !
     ALLOCATE( rho_sum_resp_x( dfftp%nr1 ) )
     ALLOCATE( rho_sum_resp_y( dfftp%nr2 ) )
     ALLOCATE( rho_sum_resp_z( dfftp%nr3 ) )
     !
     rho_sum_resp_x = 0.D0
     rho_sum_resp_y = 0.D0
     rho_sum_resp_z = 0.D0
     !
     DO ir=1,dfftp%nnr
        !
        i=cube_save(ir,1)+1
        j=cube_save(ir,2)+1
        k=cube_save(ir,3)+1
        !
        rho_sum_resp_x(i)=rho_sum_resp_x(i)+rho_1(ir,1)
        rho_sum_resp_y(j)=rho_sum_resp_y(j)+rho_1(ir,1)
        rho_sum_resp_z(k)=rho_sum_resp_z(k)+rho_1(ir,1)
        !
     ENDDO
     !
#if defined(__MPI)
     CALL mp_sum(rho_sum_resp_x, intra_bgrp_comm)
     CALL mp_sum(rho_sum_resp_y, intra_bgrp_comm)
     CALL mp_sum(rho_sum_resp_z, intra_bgrp_comm)
     IF (ionode) THEN
#endif
        WRITE(stdout,'(5X,"Dumping plane sums of densities for &
             &iteration ",I4)') LR_iteration
        !
        filename = TRIM(prefix) // ".density_x"
        tempfile = TRIM(tmp_dir) // TRIM(filename)
        !
        OPEN (158, file = tempfile, form = 'formatted', status = &
             &'unknown', position = 'append') 
        !
        DO i=1,dfftp%nr1
           WRITE(158,*) rho_sum_resp_x(i)
        ENDDO
        !
        CLOSE(158)
        !
        filename = TRIM(prefix) // ".density_y"
        tempfile = TRIM(tmp_dir) // TRIM(filename)
        !
        OPEN (158, file = tempfile, form = 'formatted', status = &
             &'unknown', position = 'append')
        !
        DO i=1,dfftp%nr2
           WRITE(158,*) rho_sum_resp_y(i)
        ENDDO
        !
        CLOSE(158)
        !
        filename = TRIM(prefix) // ".density_z"
        tempfile = TRIM(tmp_dir) // TRIM(filename)
        !
        OPEN (158, file = tempfile, form = 'formatted', status = &
             &'unknown', position = 'append')
        !
        DO i=1,dfftp%nr3
           WRITE(158,*) rho_sum_resp_z(i)
        ENDDO
        !
        CLOSE(158)
        !
#if defined(__MPI)
     ENDIF
#endif
     !
     DEALLOCATE( rho_sum_resp_x )
     DEALLOCATE( rho_sum_resp_y )
     DEALLOCATE( rho_sum_resp_z )
     !
  ENDIF
  !
  IF (charge_response == 1 .AND. response_calc) THEN
    !
    IF (LR_iteration < itermax) WRITE(stdout,'(5x,"Calculating total &
         &response charge density")')
    ! the charge response, it is actually equivalent to an element of
    ! V^T . phi_v where V^T is the is the transpose of the Krylov
    ! subspace generated by the Lanczos algorithm. The total charge
    ! density can be written as,
    !
    ! \sum_(lanczos iterations) (V^T.phi_v) . w_T
    !
    ! Where w_T is the corresponding eigenvector from the solution of,
    !
    ! (w-L)e_1 = w_T
    !
    ! notice that rho_1 is already reduced across pools above, so no
    ! parallelization is necessary 
    !
    ! the lr_calc_dens corresponds to q of x only in even iterations
    !
    IF (resonance_condition) THEN
       !
       ! Singular matrix, the broadening term dominates, phi' has
       ! strong imaginary component 
       !
       ! Using BLAS here would result in cmplx(rho_1(:,1),0.0d0 ,dp)
       ! being copied into a NEW array due to the call being to an
       !  F77 funtion. 
       !
       rho_1_tot_im(1:dfftp%nnr,:) = rho_1_tot_im(1:dfftp%nnr,:) &
            & +  w_T(LR_iteration) * cmplx(rho_1(1:dfftp%nnr,:),0.0d0,dp) 
       !
    ELSE
       !
       ! Not at resonance.
       ! The imaginary part is neglected, these are the non-absorbing
       !  oscillations
       !
       rho_1_tot(1:dfftp%nnr,:) = rho_1_tot(1:dfftp%nnr,:) &
            & +  dble( w_T(LR_iteration) ) * rho_1(1:dfftp%nnr,:)
       !   
    ENDIF
    !
 ENDIF
 !
 CALL stop_clock('lr_calc_dens')
 !
 RETURN
 !
CONTAINS
  !
  SUBROUTINE lr_calc_dens_gamma
    !
    ! Gamma_only case
    !
    USE becmod,              ONLY : bec_type, becp, calbec
    USE lr_variables,        ONLY : becp_1, tg_revc0
    USE io_global,           ONLY : stdout
    USE realus,              ONLY : real_space, invfft_orbital_gamma,&
                                    & initialisation_level,&
                                    & calbec_rs_gamma,&
                                    & add_vuspsir_gamma, v_loc_psir,&
                                    & real_space_debug 
    USE mp_global,           ONLY : ibnd_start, ibnd_end, inter_bgrp_comm, &
                                    me_bgrp, me_pool
    USE mp,                  ONLY : mp_sum
    USE realus,              ONLY : tg_psic
    USE fft_base,            ONLY : dffts, dtgs

    IMPLICIT NONE
    !
    INTEGER :: ibnd_start_gamma, ibnd_end_gamma
    INTEGER :: v_siz, incr, ioff, idx
    REAL(DP), ALLOCATABLE :: tg_rho(:)
    !
    ibnd_start_gamma = ibnd_start
    IF (MOD(ibnd_start, 2)==0) ibnd_start_gamma = ibnd_start + 1
    ibnd_end_gamma = MAX(ibnd_end, ibnd_start_gamma)
    !
    incr = 2
    !
    IF ( dtgs%have_task_groups ) THEN
       !
       v_siz =  dtgs%tg_nnr * dtgs%nogrp
       !
       incr = 2 * dtgs%nogrp
       !
       ALLOCATE( tg_rho( v_siz ) )
       tg_rho= 0.0_DP
       !
    ENDIF
    !
    DO ibnd = ibnd_start_gamma, ibnd_end_gamma, incr
       !
       ! FFT: evc1 -> psic
       !
       CALL invfft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
       !
       IF (dtgs%have_task_groups) THEN
          !
          ! Now the first proc of the group holds the first two bands
          ! of the 2*dtgs%nogrp bands that we are processing at the same time,
          ! the second proc. holds the third and fourth band
          ! and so on.
          !
          ! Compute the proper factor for each band
          !
          DO idx = 1, dtgs%nogrp
             IF( dtgs%nolist( idx ) == me_bgrp ) EXIT
          ENDDO
          !
          ! Remember two bands are packed in a single array :
          ! proc 0 has bands ibnd   and ibnd+1
          ! proc 1 has bands ibnd+2 and ibnd+3
          ! ....
          !
          idx = 2 * idx - 1
          !
          IF ( idx + ibnd - 1 < nbnd ) THEN
!         IF( idx + ibnd - 1 < ibnd_end_gamma ) THEN
             w1 = wg( idx + ibnd - 1, 1) / omega
             w2 = wg( idx + ibnd    , 1) / omega
          ELSE IF( idx + ibnd - 1 == nbnd ) THEN
!         ELSE IF( idx + ibnd - 1 == ibnd_end_gamma ) THEN
             w1 = wg( idx + ibnd - 1, 1) / omega
             w2 = w1
          ELSE
             w1 = 0.0d0
             w2 = w1
          END IF
          !
          DO ir = 1, dtgs%tg_npp( me_bgrp + 1 ) * dffts%nr1x * dffts%nr2x
             tg_rho(ir) = tg_rho(ir) &
                  + 2.0d0*(w1*real(tg_revc0(ir,ibnd,1),dp)*real(tg_psic(ir),dp)&
                  + w2*aimag(tg_revc0(ir,ibnd,1))*aimag(tg_psic(ir)))
          ENDDO
          !
       ELSE
          !
          ! Set weights of the two real bands now in psic
          !
          w1 = wg(ibnd,1)/omega
          !
          IF (ibnd<nbnd) THEN
             w2 = wg(ibnd+1,1)/omega
          ELSE
             w2 = 0.d0
          ENDIF
          !
          ! (n'(r,w) = 2*sum_v (psi_v(r) . q_v(r,w))
          ! where psi are the ground state valance orbitals
          ! and q_v are the standard batch representation (rotated)
          ! response orbitals.
          ! Here, since the i-th iteration is the best approximation we
          ! have for the most dominant eigenvalues/vectors, an estimate
          ! for the response charge density can be calculated. This is
          ! in no way the final response charge density.  
          ! The loop is over real space points.
          !
          DO ir = 1, dffts%nnr
             rho_1(ir,1) = rho_1(ir,1) &
                  + 2.0d0*(w1*real(revc0(ir,ibnd,1),dp)*real(psic(ir),dp)&
                  + w2*aimag(revc0(ir,ibnd,1))*aimag(psic(ir)))
          ENDDO
          !
          ! OBM - psic now contains the response functions in real space.
          ! Eagerly putting all the real space stuff at this point. 
          !
          ! Notice that betapointlist() is called in lr_readin at the
          ! very beginning.
          !
          IF ( real_space_debug > 6 .AND. okvan) THEN
             ! The rbecp term
             CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
             !
          ENDIF
          !
          ! End of real space stuff
          !
          IF (lr_exx) CALL lr_exx_kernel_int ( evc1(:,:,1), ibnd, nbnd, 1 )
          !
       ENDIF
       !
    ENDDO
    !
    IF (dtgs%have_task_groups) THEN
       !
       ! reduce the group charge
       !
       CALL mp_sum( tg_rho, gid = dtgs%ogrp_comm )
       !
       ioff = 0
       DO idx = 1, dtgs%nogrp
          IF ( me_bgrp == dtgs%nolist( idx ) ) EXIT
          ioff = ioff + dffts%nr1x * dffts%nr2x * dffts%npp( dtgs%nolist( idx ) + 1 )
       END DO
       !
       ! copy the charge back to the processor location
       !
       DO ir = 1, dffts%nnr
          rho_1(ir,1) = rho_1(ir,1) + tg_rho(ir+ioff)
       ENDDO
       !
    ENDIF
    !
    ! If we have a US pseudopotential we compute here the becsum term. 
    ! This corresponds to the right hand side of the formula (36) in
    ! the ultrasoft paper. 
    !
    ! Be careful about calling lr_calc_dens, as it modifies this
    ! globally.
    !
    IF ( okvan ) THEN
       !
       scal = 0.0d0
       becsum(:,:,:) = 0.0d0
       !
       IF ( real_space_debug <= 6) THEN 
          ! In real space, the value is calculated above
          CALL calbec(ngk(1), vkb, evc1(:,:,1), becp)
          !
       ENDIF
       !
       CALL start_clock( 'becsum' )
       !
       DO ibnd = 1, nbnd
          !
          scal = 0.0d0
          w1 = wg(ibnd,1)
          ijkb0 = 0
          !
          DO np = 1, ntyp
             !
             IF ( upf(np)%tvanp ) THEN
                !
                DO na = 1, nat
                   !
                   IF ( ityp(na) == np ) THEN
                      !
                      ijh = 1
                      !
                      DO ih = 1, nh(np)
                         !
                         ikb = ijkb0 + ih
                         !
                         becsum(ijh,na,current_spin) = &
                              becsum(ijh,na,current_spin) + &
                              2.d0 * w1 * becp%r(ikb,ibnd) *&
                              & becp_1(ikb,ibnd)
                         !
                         scal = scal + qq(ih,ih,np) *1.d0 *&
                              &  becp%r(ikb,ibnd) * becp_1(ikb,ibnd)
                         !
                         ijh = ijh + 1
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            !
                            jkb = ijkb0 + jh
                            !
                            becsum(ijh,na,current_spin) = &
                                 becsum(ijh,na,current_spin) + &
                                 w1 * 2.D0 * (becp_1(ikb,ibnd) * &
                                 &becp%r(jkb,ibnd) + & 
                                 becp_1(jkb,ibnd) * becp%r(ikb,ibnd))
                            !
                            scal = scal + qq(ih,jh,np) * 1.d0 *&
                                 & (becp%r(ikb,ibnd) * &
                                 &becp_1(jkb, ibnd) + &
                                 &becp%r(jkb,ibnd) * becp_1(ikb,ibnd))
                            !
                            ijh = ijh + 1
                            !
                         ENDDO
                         !
                      ENDDO
                      !
                      ijkb0 = ijkb0 + nh(np)
                      !
                   ENDIF
                   !
                ENDDO
                !
             ELSE
                !
                DO na = 1, nat
                   !
                   IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                   !
                ENDDO
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
       CALL stop_clock( 'becsum' )
       !
    ENDIF
    !
    IF ( dtgs%have_task_groups ) THEN
       DEALLOCATE( tg_rho )
    END IF
    !   
    RETURN
    !
  END SUBROUTINE lr_calc_dens_gamma
!-------------------------------------------------------------------  
 SUBROUTINE lr_calc_dens_k
    !
    ! Generalised k-points case
    ! I. Timrov: Task groups parallelisation is not implemented here.
    !
    USE becmod,              ONLY : bec_type, becp, calbec
    USE lr_variables,        ONLY : becp1_c
    USE realus,              ONLY : invfft_orbital_k
    !
    IMPLICIT NONE
    !
    DO ik=1,nks
       DO ibnd=1,nbnd
          !
          ! FFT: evc1(:,:,ik) -> psic(:)
          !
          psic(:) = (0.0d0,0.0d0)
          !
          DO ig = 1, ngk(ik)
             psic(nls(igk_k(ig,ik)))=evc1(ig,ibnd,ik)
          ENDDO
          !
          CALL invfft ('Wave', psic, dffts)
          !
          ! I. Timrov: Try to use invfft_orbital_k
          ! CALL invfft_orbital_k(evc1(:,:,ik), ibnd, nbnd, ik)
          !
          w1 = wg(ibnd,ik)/omega
          !
          ! loop over real space points
          !
          DO ir=1,dffts%nnr
             rho_1c(ir,:) = rho_1c(ir,:) &
                  + 2.0d0*w1*conjg(revc0(ir,ibnd,ik))*psic(ir)
          ENDDO
          !
          IF (lr_exx) CALL lr_exx_kernel_int (evc1(:,:,ik), ibnd, nbnd, ik )
          !
       ENDDO
    ENDDO
    !
    ! If we have a US pseudopotential we compute here the becsum term
    !
    IF ( okvan .and. nkb>0 ) THEN
       !
       DO ik =1,nks
          !
          ! Calculate the beta-functions vkb
          !
          CALL init_us_2(ngk(ik),igk_k(1,ik),xk(1,ik),vkb)
          !
          scal = 0.0d0
          becsum(:,:,:) = 0.0d0
          !
          ! Calculate the product of beta-functions vkb with the 
          ! wavefunctions evc1 : becp%k = <vkb|evc1>
          !
          CALL calbec(ngk(ik),vkb,evc1(:,:,ik),becp)
          !
          CALL start_clock( 'becsum' )
          !
          DO ibnd = 1, nbnd
             !
             scal = 0.0d0
             w1 = wg(ibnd,ik)
             ijkb0 = 0
             !
             DO np = 1, ntyp
                !
                IF ( upf(np)%tvanp ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) THEN
                         !
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            becsum(ijh,na,current_spin) = &
                                 becsum(ijh,na,current_spin) + &
                                 2.d0 * w1 * &
                                 &DBLE(CONJG(becp%k(ikb,ibnd)) *&
                                 & becp1_c(ikb,ibnd,ik)) 
                            !
                            scal = scal + qq(ih,ih,np) * 1.d0 *&
                                 &  DBLE(CONJG(becp%k(ikb,ibnd)) *&
                                 & becp1_c(ikb,ibnd,ik))
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               becsum(ijh,na,current_spin) = &
                                    becsum(ijh,na,current_spin) + &
                                    w1 * 2.d0 * DBLE(&
                                    & CONJG(becp1_c(ikb,ibnd,ik)) *&
                                    & becp%k(jkb,ibnd) + &
                                    becp1_c(jkb,ibnd,ik) *&
                                    & CONJG(becp%k(ikb,ibnd)))
                               !
                               scal = scal + qq(ih,jh,np) * 1.d0 * &
                                    & DBLE(CONJG(becp%k(ikb,ibnd)) *&
                                    & becp1_c(jkb,ibnd,ik)+&
                                    & becp%k(jkb,ibnd) * &
                                    & CONJG(becp1_c(ikb,ibnd,ik)))
                               !
                               ijh = ijh + 1
                               !
                            ENDDO
                            !
                         ENDDO
                         !
                         ijkb0 = ijkb0 + nh(np)
                         !
                      ENDIF
                      !
                   ENDDO
                   !
                ELSE
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                      !
                   ENDDO
                   !
                ENDIF
                !
             ENDDO
             !
          ENDDO
          !
          CALL stop_clock( 'becsum' )
          !
       ENDDO
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_calc_dens_k
!--------------------------------------------------------------------
END SUBROUTINE lr_calc_dens
!--------------------------------------------------------------------

