!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_apply_liouvillian( evc1, evc1_new, interaction )
  !---------------------------------------------------------------------
  !
  ! Applies the linear response operator to response wavefunctions
  ! S^{-1} P_c^+(k) { (H - E*S)*psi(k) + V_HXC*psi0(k) } 
  !
  ! Note 1: S^{-1} P_c^+(k) = P_c(k) S^{-1}
  ! Note 2: In the norm-conserving case: S=1, S^{-1}=1.
  ! 
  ! Or to be more exact this routine is responsible for calculating
  ! L.q(i) and (L^T).p(i), where q is evc1.
  ! It returns partial qdash(i+1) and pdash(i+1) in evc1_new.
  !
  ! interaction=.true. corresponds to eq.(32)
  ! interaction=.false. corresponds to eq.(33)
  ! in B. Walker and R. Gebauer, J. Chem. Phys., 127, 164106 (2007)
  !
  ! Modified by Osman Baris Malcioglu in 2009
  ! Modified by Simone Binnie in 2012 (EXX)
  ! Modified by Xiaochuan Ge  in 2013 (Davidson)
  ! Modified by Iurii Timrov  in 2014 (Environ)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : ityp, nat, ntyp=>nsp
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dffts, dfftp, dtgs
  USE fft_interfaces,       ONLY : fwfft
  USE gvecs,                ONLY : nls, nlsm
  USE gvect,                ONLY : nl, ngm, gstart, g, gg
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lr_variables,         ONLY : evc0, sevc0, revc0, rho_1, rho_1c, &
                                 & ltammd, size_evc, no_hxc, lr_exx, &
                                 & scissor, davidson, lr_verbosity
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE uspp_param,           ONLY : nhm, nh
  USE wavefunctions_module, ONLY : psic
  USE wvfct,                ONLY : nbnd, npwx, g2kin, et
  USE control_flags,        ONLY : gamma_only
  USE realus,               ONLY : real_space, invfft_orbital_gamma,&
                                   & initialisation_level,&
                                   & fwfft_orbital_gamma,&
                                   & calbec_rs_gamma, newq_r, &
                                   & add_vuspsir_gamma, v_loc_psir,   &
                                   & s_psir_gamma, real_space_debug,  &
                                   & betasave, box_beta, maxbox_beta
  USE dfunct,               ONLY : newq
  USE control_flags,        ONLY : tqr
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE mp_global,            ONLY : intra_bgrp_comm
  USE noncollin_module,     ONLY : npol
  USE becmod,               ONLY : bec_type, becp, calbec
  USE lr_exx_kernel
  USE dv_of_drho_lr
#if defined(__ENVIRON)
  USE plugin_flags,         ONLY : use_environ
  USE scf,                  ONLY : rho
  USE solvent_tddfpt,       ONLY : calc_vsolvent_tddfpt
#endif
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: evc1(npwx*npol,nbnd,nks)
  COMPLEX(DP), INTENT(OUT) :: evc1_new(npwx*npol,nbnd,nks)
  LOGICAL,     INTENT(IN)  :: interaction
  !
  ! Local variables
  !
  INTEGER :: ir, ibnd, ik, ig, ia, mbia
  INTEGER :: ijkb0, na, nt, ih, jh, ikb, jkb, iqs,jqs
  REAL(DP), ALLOCATABLE :: dvrs(:,:), dvrss(:), d_deeq(:,:,:,:), &
                           & w1(:), w2(:)
  COMPLEX(DP), ALLOCATABLE :: dvrs_temp(:,:), spsi1(:,:), dvrsc(:,:), &
                              & dvrssc(:), sevc1_new(:,:,:)
  !
  ! Environ related arrays
  !
  REAL(DP), ALLOCATABLE :: &
          dv_pol(:), &  ! response polarization potential
          dv_epsilon(:) ! response dielectric potential
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_apply_liouvillian>")')
  ENDIF
  !
  CALL start_clock('lr_apply')
  !
  IF (interaction)      CALL start_clock('lr_apply_int')
  IF (.not.interaction) CALL start_clock('lr_apply_no')
  !
  ALLOCATE( d_deeq(nhm, nhm, nat, nspin) )
  d_deeq(:,:,:,:)=0.0d0
  !
  ALLOCATE( spsi1(npwx, nbnd) )
  spsi1(:,:)=(0.0d0,0.0d0)
  !
  ALLOCATE(sevc1_new(npwx*npol,nbnd,nks))
  sevc1_new(:,:,:) = (0.0d0,0.0d0)
  !
  evc1_new(:,:,:) = (0.0d0,0.0d0)
  !
  IF ( interaction ) THEN 
     !
     ! Calculate the full L
     !
     IF (gamma_only) THEN
        ALLOCATE( dvrs(dfftp%nnr, nspin) )
        ALLOCATE( dvrss(dffts%nnr) )
        dvrs(:,:)=0.0d0
        dvrss(:)=0.0d0
     ELSE
        ALLOCATE( dvrsc(dfftp%nnr, nspin) )
        ALLOCATE( dvrssc(dffts%nnr) )
        dvrsc(:,:)=0.0d0
        dvrssc(:) =0.0d0
     ENDIF
     !
     ! Calculation of the charge density response
     !
     CALL lr_calc_dens( evc1, .false. )
     !
     ! Given the change of the charge density
     ! we calculate the change of the Hartree and XC potential and
     ! put it in dvrss
     !
     IF (no_hxc) THEN
        !
        ! With no_hxc=.true. we recover the independent electron 
        ! approximation, so we zero the interation.
        !
        IF (gamma_only) THEN
           dvrs(:,1) = 0.0d0
           CALL interpolate (dvrs(:,1),dvrss,-1)
        ELSE
           dvrsc(:,1) = 0.0d0
           CALL cinterpolate (dvrsc(:,1),dvrssc,-1)
        ENDIF
        !
     ELSE
        !
        IF (gamma_only) THEN
           !
           dvrs(:,1) = rho_1(:,1)
           !
           ! In the gamma_only case dvrs is real, but dv_of_drho expects
           ! a complex array on input, hence this temporary variable.
           !
           ALLOCATE( dvrs_temp(dfftp%nnr, nspin) )
           !
           dvrs_temp = CMPLX( dvrs, 0.0d0, kind=DP )         
           !
           DEALLOCATE ( dvrs )  ! to save memory
           !
           CALL dv_of_drho(dvrs_temp,.FALSE.)
           !
           ALLOCATE ( dvrs(dfftp%nnr, nspin) )
           ! 
           dvrs = DBLE(dvrs_temp)
           !
           DEALLOCATE(dvrs_temp)
           !
#if defined(__ENVIRON)
           !
           IF ( use_environ ) THEN
              !
              ALLOCATE( dv_pol(dfftp%nnr) )
              ALLOCATE( dv_epsilon(dfftp%nnr) )
              dv_pol(:) = 0.0d0
              dv_epsilon(:) = 0.0d0
              !
              IF (.not.davidson) THEN
                 WRITE( stdout, '(5x,"ENVIRON: Calculate the response &
                             & polarization and dielectric potentials")' )
              ENDIF
              !
              CALL calc_vsolvent_tddfpt(dfftp%nnr, nspin, rho%of_r(:,1), &
                                       & rho_1(:,1), dv_pol, dv_epsilon)
              !
              ! Add the response polarization and dielectric potentials
              ! to the response HXC potential.
              !
              dvrs(:,1) = dvrs(:,1) + dv_pol(:) + dv_epsilon(:)
              !
              DEALLOCATE( dv_pol )
              DEALLOCATE( dv_epsilon )
              !
           ENDIF
#endif           
           !
        ELSE
           !
           dvrsc(:,1) = rho_1c(:,1)
           !
           CALL dv_of_drho(dvrsc,.FALSE.)
           !
        ENDIF
        !
        IF ( okvan )  THEN
           IF ( tqr ) THEN
              CALL newq_r(dvrs,d_deeq,.TRUE.)
           ELSE
              ALLOCATE( psic(dfftp%nnr) )
              psic(:)=(0.0d0,0.0d0)
              !
              CALL newq(dvrs,d_deeq,.TRUE.)
              !
              DEALLOCATE( psic )
           ENDIF
        ENDIF
        !
        CALL add_paw_to_deeq(d_deeq)
        !
        ! Put the interaction on the smooth grid.
        !
        IF (gamma_only) THEN
           CALL interpolate (dvrs(:,1),dvrss,-1)
        ELSE
           CALL cinterpolate (dvrsc(:,1),dvrssc,-1)
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
  ! S. Binnie: Make sure the psic workspace is availible.
  !
  ALLOCATE ( psic (dffts%nnr) )
  !
  IF ( gamma_only ) THEN
     CALL lr_apply_liouvillian_gamma()
  ELSE
     CALL lr_apply_liouvillian_k()
  ENDIF
  !
  DEALLOCATE ( psic )
  !
  IF ( (interaction .or. lr_exx) .and. (.not.ltammd)  ) THEN
     !
     !   Normal interaction
     !
     if (.not. davidson) WRITE(stdout,'(5X,"lr_apply_liouvillian: &
                                 & applying interaction: normal")')
     !
     ! Here we add the two terms: 
     ! [H(k) - E(k)] * evc1(k) + dV_HXC * revc0(k)
     !
     CALL zaxpy(size_evc,CMPLX(1.0d0,0.0d0,kind=DP),&
                       & evc1_new(:,:,:), 1, sevc1_new(:,:,:), 1)
     !
  ELSEIF ( interaction .and. ltammd ) THEN
     !
     !   Tamm-Dancoff approximation
     !
     IF (.not. davidson) WRITE(stdout,'(5X,"lr_apply_liouvillian:&
                          & applying interaction: Tamm-Dancoff")')
     !
     !   Here evc1_new contains the interaction
     !
     CALL zaxpy(size_evc,CMPLX(0.5d0,0.0d0,kind=DP),&
                       & evc1_new(:,:,:) , 1, sevc1_new(:,:,:),1)
     !
  ELSE
     !
     !   Non interacting
     !
     IF (.not. davidson) WRITE(stdout,'(5X,"lr_apply_liouvillian:&
                                    & not applying interaction")') 
     !
  ENDIF
  !
  IF (gstart == 2 .AND. gamma_only ) sevc1_new(1,:,:) = &
                & CMPLX( REAL( sevc1_new(1,:,:), DP ), 0.0d0, DP )
  !
  IF (gstart==2 .and. gamma_only) THEN
     DO ik=1,nks
        DO ibnd=1,nbnd
           IF (abs(aimag(sevc1_new(1,ibnd,ik)))>1.0d-12) THEN
              !
              CALL errore(' lr_apply_liouvillian ',&
                   'Imaginary part of G=0 '// &
                   'component does not equal zero',1)
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  ! Apply the projector on empty states P_c^+.
  ! Note: The projector P_c^+ can be applied only
  ! to the responce HXC term, but in order to increase
  ! the stability of the Lanczos chain we apply it 
  ! also to the [H(k) - E(k)] * evc1(k) term.
  ! However, we can apply the projector P_c^+ to
  ! [H(k) - E(k)] * evc1(k) only in the case of systems
  ! with the energy gap (insulators, molecules, etc.).
  ! In the case of metals, P_c^+ is a more complicated object
  ! (see orthogonalize.f90), and hence we cannot apply it to 
  ! [H(k) - E(k)] * evc1(k). Keep this in mind if you want to
  ! generalize this subroutine to metals. 
  !
  DO ik = 1, nks
     !
     CALL orthogonalize(sevc1_new(:,:,ik), evc0(:,:,ik), ik, ik, &
                                  & sevc0(:,:,ik), ngk(ik), .true.)
     sevc1_new(:,:,ik) = -sevc1_new(:,:,ik)
     !
  ENDDO 
  !
  ! Here we apply the S^{-1} operator.
  ! See equations after Eq.(47) of B. Walker et al., J. Chem. Phys.
  !  127, 164106 (2007).  
  !
  ! interaction=.false.: S^{-1} P_c^+(k) (H(k)-E(k)*S) * evc1(k) 
  ! interaction=.true.:  S^{-1} P_c^+(k) { (H(k)-E(k)*S) * evc1(k) + dV_HXC * revc0(k) } 
  !
  DO ik = 1, nks
     !
     CALL lr_sm1_psi (.FALSE., ik, npwx, ngk(ik), nbnd, &
                       & sevc1_new(1,1,ik), evc1_new(1,1,ik))
     ! 
  ENDDO
  !
  IF (allocated(dvrs)) DEALLOCATE(dvrs)
  IF (allocated(dvrss)) DEALLOCATE(dvrss)
  DEALLOCATE(d_deeq)
  DEALLOCATE(spsi1)
  DEALLOCATE(sevc1_new)
  !
  IF (interaction)      CALL stop_clock('lr_apply_int')
  IF (.not.interaction) CALL stop_clock('lr_apply_no')
  !
  CALL stop_clock('lr_apply')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_apply_liouvillian_gamma()
    !
    ! Gamma only case 
    !
    USE lr_variables,             ONLY : becp_1, tg_revc0
    USE realus,                   ONLY : tg_psic
    USE mp_global,                ONLY : me_bgrp
    USE fft_parallel,             ONLY : tg_gather
    USE mp_global,                ONLY : ibnd_start, ibnd_end, inter_bgrp_comm
    USE mp,                       ONLY : mp_sum
    USE lr_exx_kernel,            ONLY : lr_exx_sum_int
    !
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE :: becp2(:,:)
    REAL(DP), ALLOCATABLE :: tg_dvrss(:)
    INTEGER :: v_siz, incr, ioff
    INTEGER :: ibnd_start_gamma, ibnd_end_gamma
    !
    incr = 2
    !
    IF ( nkb > 0 .and. okvan ) THEN
       !
       ALLOCATE(becp2(nkb,nbnd))
       becp2(:,:) = 0.0d0
       !
    ENDIF
    !
    ! Now apply to the ground state wavefunctions
    ! and convert to real space
    !
    IF ( interaction ) THEN
       !
       CALL start_clock('interaction')
       !
       IF (nkb > 0 .and. okvan) THEN
          ! calculation of becp2
          becp2(:,:) = 0.0d0
          !
          ijkb0 = 0
          !
          DO nt = 1, ntyp
             !
             DO na = 1, nat
                !
                IF ( ityp(na) == nt ) THEN
                   !
                   DO ibnd = 1, nbnd
                      !
                      DO jh = 1, nh(nt)
                         !
                         jkb = ijkb0 + jh
                         !
                         DO ih = 1, nh(nt)
                            !
                            ikb = ijkb0 + ih
                            becp2(ikb, ibnd) = becp2(ikb, ibnd) + &
                                 d_deeq(ih,jh,na,1) * becp_1(jkb,ibnd)
                            !
                         ENDDO
                         !
                      ENDDO
                      !
                   ENDDO
                   !
                   ijkb0 = ijkb0 + nh(nt)
                   !
                ENDIF
                !
             ENDDO
             !
          ENDDO
          !end: calculation of becp2
       ENDIF

      IF ( dtgs%have_task_groups ) THEN
         !
         v_siz =  dtgs%tg_nnr * dtgs%nogrp
         !
         incr = 2 * dtgs%nogrp
         !
         ALLOCATE( tg_dvrss(1:v_siz) )
         tg_dvrss=0.0d0
         !
         CALL tg_gather(dffts, dtgs, dvrss, tg_dvrss)
         !
      ENDIF
       !
       ! evc1_new is used as a container for the interaction
       !
       evc1_new(:,:,:) = (0.0d0,0.0d0)
       !
       ibnd_start_gamma = ibnd_start
       IF (MOD(ibnd_start, 2)==0) ibnd_start_gamma = ibnd_start + 1
       ibnd_end_gamma = MAX(ibnd_end, ibnd_start_gamma)
       !
       IF (lr_exx) CALL lr_exx_sum_int()
       !
       DO ibnd = ibnd_start_gamma ,ibnd_end_gamma, incr
!      DO ibnd = 1,nbnd,2
          !
          ! Product with the potential vrs = (vltot+vr)
          ! revc0 is on smooth grid. psic is used up to smooth grid
          !
          IF (dtgs%have_task_groups) THEN
             !
             DO ir=1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
                !
                tg_psic(ir) = tg_revc0(ir,ibnd,1)*CMPLX(tg_dvrss(ir),0.0d0,DP)
                !
             ENDDO
             !
          ELSE
             !
             DO ir = 1,dffts%nnr
                !
                psic(ir) = revc0(ir,ibnd,1)*CMPLX(dvrss(ir),0.0d0,DP)
                !
             ENDDO
             !
          ENDIF
          !
          IF (lr_exx) THEN
             CALL lr_exx_apply_revc_int(psic, ibnd, nbnd,1)
          ENDIF
          !
          IF (real_space_debug > 7 .and. okvan .and. nkb > 0) THEN
          !THE REAL SPACE PART (modified from s_psi)
                  !fac = sqrt(omega)
                  !
                  ijkb0 = 0
                  iqs = 0
                  jqs = 0
                  !
                  DO nt = 1, ntyp
                  !
                    DO ia = 1, nat
                      !
                      IF ( ityp(ia) == nt ) THEN
                        !
                        mbia = maxbox_beta(ia)
                        ALLOCATE( w1(nh(nt)),  w2(nh(nt)) )
                        w1 = 0.D0
                        w2 = 0.D0
                        !
                        DO ih = 1, nh(nt)
                        !
                          DO jh = 1, nh(nt)
                          !
                          jkb = ijkb0 + jh
                          w1(ih) = w1(ih) + becp2(jkb, ibnd)
                          IF ( ibnd+1 <= nbnd ) w2(ih) = w2(ih) + &
                               &becp2(jkb, ibnd+1)
                          !
                          ENDDO
                        !
                        ENDDO
                        !
                        !w1 = w1 * fac
                        !w2 = w2 * fac
                        ijkb0 = ijkb0 + nh(nt)
                        !
                        DO ih = 1, nh(nt)
                          !
                          DO ir = 1, mbia
                          !
                           iqs = jqs + ir
                           psic( box_beta(ir,ia) ) = &
                                &psic(  box_beta(ir,ia) ) + &
                                &betasave(ir,ih,ia)*&
                                &CMPLX( w1(ih), w2(ih) )
                          !
                          ENDDO
                        !
                        jqs = iqs
                        !
                        ENDDO
                        !
                      DEALLOCATE( w1, w2 )
                        !
                      ENDIF
                    !
                    ENDDO
                  !
                 ENDDO

          ENDIF
          !
          !   Back to reciprocal space 
          !
          CALL fwfft_orbital_gamma (evc1_new(:,:,1), ibnd, nbnd,.false.)
          !
       ENDDO
       !
#if defined(__MPI)
       CALL mp_sum( evc1_new(:,:,1), inter_bgrp_comm )
#endif
       IF (dtgs%have_task_groups) DEALLOCATE (tg_dvrss)
       !
       IF( nkb > 0 .and. okvan .and. real_space_debug <= 7) THEN
          !The non real_space part
          CALL dgemm( 'N', 'N', 2*ngk(1), nbnd, nkb, 1.d0, vkb, &
               2*npwx, becp2, nkb, 1.d0, evc1_new, 2*npwx )
          !
       ENDIF
       !
       CALL stop_clock('interaction')
       !
    ENDIF
    !
    IF (lr_exx .AND. .NOT.interaction) CALL lr_exx_kernel_noint(evc1,evc1_new)
    !
    ! The kinetic energy g2kin was already computed when
    ! calling the routine lr_solve_e.
    !
    ! Compute sevc1_new = H*evc1
    !
    CALL h_psi(npwx,ngk(1),nbnd,evc1(1,1,1),sevc1_new(1,1,1))
    !
    ! Compute spsi1 = S*evc1 
    !
    IF (real_space_debug > 9 ) THEN
        DO ibnd = 1,nbnd,2
           CALL invfft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
           CALL s_psir_gamma(ibnd,nbnd)
           CALL fwfft_orbital_gamma(spsi1,ibnd,nbnd)
        ENDDO
    ELSE
       CALL s_psi(npwx,ngk(1),nbnd,evc1(1,1,1),spsi1)
    ENDIF
    !
    !   Subtract the eigenvalues
    !
    DO ibnd = 1,nbnd
       !
       CALL zaxpy(ngk(1), CMPLX(-(et(ibnd,1)-scissor),0.0d0,DP), &
                      & spsi1(:,ibnd), 1, sevc1_new(:,ibnd,1), 1)
       !
    ENDDO
    !
    IF ( nkb > 0 .and. okvan ) DEALLOCATE(becp2)
    !
    RETURN
    !
END SUBROUTINE lr_apply_liouvillian_gamma
  
SUBROUTINE lr_apply_liouvillian_k()
    !
    ! Generalised k-point case
    !
    USE lr_variables,        ONLY : becp1_c
    USE wvfct,               ONLY : current_k
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), ALLOCATABLE :: becp2(:,:)
    !
    IF ( nkb > 0 .AND. okvan ) THEN
       !
       ALLOCATE(becp2(nkb,nbnd))
       becp2(:,:) = (0.0d0,0.0d0)
       !
    ENDIF
    !
    ! Now apply to the ground state wavefunctions
    ! and convert to real space
    !
    IF ( interaction ) THEN
       !
       CALL start_clock('interaction')
       !
       !   evc1_new is used as a container for the interaction
       !
       evc1_new(:,:,:) = (0.0d0,0.0d0)
       !
       DO ik=1,nks
          !
          DO ibnd = 1,nbnd
             !
             !   Product with the potential vrs = (vltot+vr)
             !
             DO ir = 1,dffts%nnr
                !
                psic(ir) = revc0(ir,ibnd,ik)*dvrssc(ir)
                !
             ENDDO
             !
!            IF (lr_exx) psic(:) = psic(:) + 0.5d0 * revc_int_c(:,ibnd,ik)
             !
             IF (lr_exx) THEN
                CALL lr_exx_apply_revc_int(psic, ibnd, nbnd, ik)
             ENDIF
             !
             !   Back to reciprocal space
             !
             CALL fwfft ('Wave', psic, dffts)
             !
             DO ig = 1,ngk(ik)
                !
                evc1_new(ig,ibnd,ik) = psic(nls(igk_k(ig,ik)))
                !
             ENDDO
             !
          ENDDO
          !
       ENDDO
       !
       CALL stop_clock('interaction')
       !
       IF ( nkb > 0 .and. okvan ) THEN
          !
          DO ik = 1,nks
             !
             CALL init_us_2(ngk(ik),igk_k(1,ik),xk(1,ik),vkb)
             !
             becp2(:,:) = 0.0d0
             !
             ijkb0 = 0
             !
             DO nt = 1, ntyp
                !
                DO na = 1, nat
                   !
                   IF ( ityp(na) == nt ) THEN
                      !
                      DO ibnd = 1, nbnd
                         !
                         DO jh = 1, nh(nt)
                            !
                            jkb = ijkb0 + jh
                            !
                            DO ih = 1, nh(nt)
                               !
                               ikb = ijkb0 + ih
                               becp2(ikb, ibnd) = becp2(ikb, ibnd) + &
                                    d_deeq(ih,jh,na,1) * becp1_c(jkb,ibnd,ik)
                                !
                            ENDDO
                            !
                         ENDDO
                         !
                      ENDDO
                      !
                      ijkb0 = ijkb0 + nh(nt)
                      !
                   ENDIF
                   !
                ENDDO
                !
             ENDDO
             !
             !evc1_new(ik) = evc1_new(ik) + vkb*becp2(ik)
             CALL zgemm( 'N', 'N', ngk(ik), nbnd, nkb, (1.d0,0.d0), vkb, &
                  npwx, becp2, nkb, (1.d0,0.d0), evc1_new(:,:,ik), npwx )
             !
          ENDDO
          !
       ENDIF
       !
    ENDIF
    !
    !   Call h_psi on evc1
    !
    IF (lr_exx .AND. .NOT.interaction) CALL lr_exx_kernel_noint(evc1,evc1_new)
    !
    DO ik = 1, nks
       !
       ! US case: Compute the beta function vkb at point k
       !
       CALL init_us_2(ngk(ik), igk_k(1,ik), xk(1,ik), vkb)
       !
       ! Compute the kinetic energy g2kin: (k+G)^2
       !
       CALL g2_kin(ik)
       !
       current_k = ik
       !
       ! Compute sevc1_new = H*evc1
       !
       CALL h_psi(npwx,ngk(ik),nbnd,evc1(1,1,ik),sevc1_new(1,1,ik))
       !
       ! Compute spsi1 = S*evc1
       !
       CALL s_psi(npwx,ngk(ik),nbnd,evc1(1,1,ik),spsi1)
       !
       ! Subtract the eigenvalues
       !
       DO ibnd = 1, nbnd
          !
          CALL zaxpy(ngk(ik), CMPLX(-(et(ibnd,ik)-scissor),0.0d0,dp), &
                      & spsi1(:,ibnd), 1, sevc1_new(:,ibnd,ik), 1)
          !
       ENDDO
       !
    ENDDO ! end k loop
    !
    IF ( nkb > 0 .and. okvan ) DEALLOCATE(becp2)
    !
    RETURN
    !
  END SUBROUTINE lr_apply_liouvillian_k
  !
END SUBROUTINE lr_apply_liouvillian
!-----------------------------------------------------------------------
