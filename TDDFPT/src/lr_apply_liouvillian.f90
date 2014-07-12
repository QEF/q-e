!-----------------------------------------------------------------------
SUBROUTINE lr_apply_liouvillian( evc1, evc1_new, sevc1_new, interaction )
  !---------------------------------------------------------------------
  ! Applies linear response operator to response wavefunctions
  !   (H - E)*psi(k+q) + HXC   
  ! 
  ! Or to be more exact this routine is responsible for calculating
  !    L.q(i) and (L^T).p(i)
  ! Where q is evc1
  ! It returns partial qdash(i+1) and pdash(i+1) in evc1_new.
  !
  ! interaction=.true. corresponds to eq.(32)
  ! interaction=.false. corresponds to eq.(33)
  ! in Ralph Gebauer, Brent Walker  J. Chem. Phys., 127, 164106 (2007)
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu in 2009
  !
  USE ions_base,            ONLY : ityp, nat, ntyp=>nsp
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft
  USE gvecs,                ONLY : nls, nlsm
  USE gvect,                ONLY : nl, ngm, gstart, g, gg
  USE fft_base,             ONLY : dfftp, dffts
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nks, xk
  USE lr_variables,         ONLY : evc0, revc0, rho_1, lr_verbosity,&
                                 & ltammd, size_evc, no_hxc, lr_exx, &
                                 & scissor, rho_1c,davidson
  USE realus,               ONLY : igk_k,npw_k
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE uspp_param,           ONLY : nhm, nh
  USE wavefunctions_module, ONLY : psic
  USE wvfct,                ONLY : nbnd, npwx, igk, g2kin, et
  USE control_flags,        ONLY : gamma_only
  USE realus,               ONLY : real_space, fft_orbital_gamma,&
                                   & initialisation_level,&
                                   & bfft_orbital_gamma,&
                                   & calbec_rs_gamma,&
                                   & add_vuspsir_gamma, v_loc_psir,&
                                   & s_psir_gamma, real_space_debug,&
                                   & betasave, box_beta, maxbox_beta,&
                                   & newq_r 
  USE lr_variables,         ONLY : lr_verbosity
  USE io_global,            ONLY : stdout
  USE dfunct,               ONLY : newq
  USE control_flags,        ONLY : tqr
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE mp_global,            ONLY : intra_bgrp_comm
  use lr_exx_kernel
#ifdef __ENVIRON
  USE plugin_flags,         ONLY : use_environ
  USE scf,                  ONLY : rho
  USE solvent_tddfpt,       ONLY : calc_vsolvent_tddfpt
#endif
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(kind=dp),INTENT(in)  :: evc1(npwx,nbnd,nks)
  COMPLEX(kind=dp),INTENT(out) :: evc1_new(npwx,nbnd,nks),&
                               & sevc1_new(npwx,nbnd,nks)
  ! output : sevc1_new = S * evc1_new
  !
  LOGICAL, INTENT(in) :: interaction
  !
  !   Local variables
  !
  INTEGER :: ir, ibnd, ik, ig, ia, mbia
  INTEGER :: ijkb0, na, nt, ih, jh, ikb, jkb, iqs,jqs
  !
  REAL(kind=dp), ALLOCATABLE :: dvrs(:,:), dvrss(:)
  REAL(kind=dp), ALLOCATABLE :: d_deeq(:,:,:,:)
  !
  ! Environ related arrays
  !
  REAL(kind=dp), ALLOCATABLE :: &
          dv_pol(:), &  ! response polarization potential
          dv_epsilon(:) ! response dielectric potential
  !
  COMPLEX(kind=dp), ALLOCATABLE :: dvrs_temp(:,:)   
  COMPLEX(kind=dp), ALLOCATABLE :: spsi1(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: dvrsc(:,:), dvrssc(:)
  !
  COMPLEX(kind=dp) :: fp, fm
  !
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: w1, w2
  !
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
  IF( interaction ) THEN 
     !
     ! Calculate the full L
     !
     IF(gamma_only) THEN
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
     ! put it in dvrss()
     !
     IF (no_hxc) THEN
        !
        ! With no_hxc=.true. we recover the independent electrion 
        ! approximation, so we zero the interation.
        !
        IF(gamma_only) THEN
           dvrs(:,1)=0.0d0
           CALL interpolate (dvrs(:,1),dvrss,-1)
        ELSE
           dvrsc(:,1)=0.0d0
           CALL cinterpolate (dvrsc(:,1),dvrssc,-1)
        ENDIF
        !
     ELSE
        !
        IF(gamma_only) THEN
           dvrs(:,1)=rho_1(:,1)
           !
           ! In the gamma_only case dvrs is real, but dv_of_drho expects
           ! a complex array on input, hence this temporary variable.
           !
           ALLOCATE( dvrs_temp(dfftp%nnr, nspin) )
           dvrs_temp = CMPLX( dvrs, 0.0d0, kind=dp )         
           !
           DEALLOCATE ( dvrs )
           !
           CALL dv_of_drho(0,dvrs_temp,.FALSE.)
           !
           ALLOCATE ( dvrs(dfftp%nnr, nspin) ) 
           dvrs=DBLE(dvrs_temp)
           DEALLOCATE(dvrs_temp)
           !
#ifdef __ENVIRON
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
           dvrsc(:,1)=rho_1c(:,1)
           !
           CALL dv_of_drho(0,dvrsc,.FALSE.)
           !
        ENDIF
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
        CALL add_paw_to_deeq(d_deeq)
        !
        ! Put the nteraction on the smooth grid.
        !
        IF (gamma_only) THEN
           CALL interpolate (dvrs(:,1),dvrss,-1)
        ELSE
           CALL cinterpolate (dvrsc(:,1),dvrssc,-1)
        ENDIF
     ENDIF
     !
  ENDIF
  !
  ! Make sure the psic workspace is availible
  !
  ALLOCATE ( psic (dffts%nnr) )
  !
  IF( gamma_only ) THEN
     !
     CALL lr_apply_liouvillian_gamma()
     !
  ELSE
     !
     CALL lr_apply_liouvillian_k()
     !
  ENDIF
  !
  DEALLOCATE ( psic )
  !
  IF ( (interaction .or. lr_exx) .and. (.not.ltammd)  ) THEN
     !
     !   Normal interaction
     !
     if(.not. davidson) WRITE(stdout,'(5X,"lr_apply_liouvillian: &
          &applying interaction: normal")')
     !
     ! Here we add the two terms: 
     ! [H(k) - E(k)] * evc1(k)  +  P_c(k) [dV_{HXC} * revc0(k)]
     !
     CALL zaxpy(size_evc,CMPLX(1.0d0,0.0d0,kind=dp),&
          &evc1_new(:,:,:), 1, sevc1_new(:,:,:), 1)
     !
     !
  ELSEIF ( interaction .and. ltammd ) THEN
     !
     !   Tamm-dancoff interaction
     !
     if(.not. davidson) WRITE(stdout,'(5X,"lr_apply_liouvillian:&
          & applying interaction: tamm-dancoff")')
     !
     !   Here evc1_new contains the interaction
     !
     CALL zaxpy(size_evc,CMPLX(0.5d0,0.0d0,kind=dp),&
          &evc1_new(:,:,:) , 1 , sevc1_new(:,:,:),1)
     !
     !
  ELSE
     !
     !   Non interacting
     !
     if(.not. davidson) WRITE(stdout,'(5X,"lr_apply_liouvillian:&
          & not applying interaction")') 
     !
  ENDIF
  !
  IF (gstart == 2 .AND. gamma_only ) sevc1_new(1,:,:) = &
       &CMPLX( REAL( sevc1_new(1,:,:), dp ), 0.0d0 ,dp )
  !
  IF(gstart==2 .and. gamma_only) THEN
     !
     DO ik=1,nks
        !
        DO ibnd=1,nbnd
           !
           !
           IF (abs(aimag(sevc1_new(1,ibnd,ik)))>1.0d-12) THEN
              !
              CALL errore(' lr_apply_liouvillian ',&
                   'Imaginary part of G=0 '// &
                   'component does not equal zero',1)
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  ! Here we apply the S^{-1} operator.
  ! See equations after Eq.(47) of B. Walker et al., J. Chem. Phys.
  !  127, 164106 (2007).  
  !
  ! S^{-1} ( H(k)*evc1(k) - E(k) * S evc1(k) ) 
  ! or
  ! S^{-1} ( H(k)*evc1(k) - E(k) * S evc1(k)  +  P_c(k) [dV_{HXC} *
  !  revc0(k)] ) 
  !
  DO ik=1,nks
     !
     CALL sm1_psi(.FALSE., ik, npwx, npw_k(ik), nbnd, &
          &sevc1_new(1,1,ik), evc1_new(1,1,ik))
     !
  ENDDO
  !
  IF (allocated(dvrs)) DEALLOCATE(dvrs)
  IF (allocated(dvrss)) DEALLOCATE(dvrss)
  DEALLOCATE(d_deeq)
  DEALLOCATE(spsi1)
  !
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
    USE lr_variables,             ONLY : becp1, tg_revc0
    USE wavefunctions_module,     ONLY : psic
    USE realus,                   ONLY : tg_psic
    USE mp_global,                ONLY : me_bgrp
    USE fft_base,                 ONLY : dffts, tg_gather
    USE mp_global,                ONLY : ibnd_start, ibnd_end, inter_bgrp_comm
    USE mp,                       ONLY : mp_sum
    USE lr_exx_kernel,            ONLY : lr_exx_sum_int
    !
    REAL(kind=dp), ALLOCATABLE :: becp2(:,:)
    REAL(kind=dp), ALLOCATABLE :: tg_dvrss(:)
    LOGICAL :: use_tg
    INTEGER :: v_siz, incr, ioff
    INTEGER   :: ibnd_start_gamma, ibnd_end_gamma
    !
    use_tg=dffts%have_task_groups
    incr = 2
    !
    IF ( nkb > 0 .and. okvan ) THEN
       !
       ALLOCATE(becp2(nkb,nbnd))
       becp2(:,:)=0.0d0
       !
    ENDIF
    !
    !    Now apply to the ground state wavefunctions
    !    and convert to real space
    !
    IF ( interaction ) THEN
       !
       CALL start_clock('interaction')

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
                                 d_deeq(ih,jh,na,1) * becp1(jkb,ibnd)
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

      IF( dffts%have_task_groups ) THEN
         !
         v_siz =  dffts%tg_nnr * dffts%nogrp
         !
         incr = 2 * dffts%nogrp
         !
         ALLOCATE( tg_dvrss(1:v_siz) )
         tg_dvrss=0.0d0
         !
         CALL tg_gather(dffts, dvrss, tg_dvrss)
         !
      ENDIF
       !
       !   evc1_new is used as a container for the interaction
       !
       evc1_new(:,:,:)=(0.0d0,0.0d0)
       !
       ibnd_start_gamma=ibnd_start
       IF (MOD(ibnd_start, 2)==0) ibnd_start_gamma = ibnd_start+1
       ibnd_end_gamma = MAX(ibnd_end, ibnd_start_gamma)
       !
       IF (lr_exx) CALL lr_exx_sum_int()
       !
       DO ibnd=ibnd_start_gamma,ibnd_end_gamma,incr
!       DO ibnd=1,nbnd,2
          !
          ! Product with the potential vrs = (vltot+vr)
          ! revc0 is on smooth grid. psic is used up to smooth grid
          !
          IF(dffts%have_task_groups) THEN
             !
             DO ir=1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
                !
                tg_psic(ir)=tg_revc0(ir,ibnd,1)*CMPLX(tg_dvrss(ir),0.0d0,dp)
                !
             ENDDO
             !
          ELSE
             !
             DO ir=1,dffts%nnr
                !
                psic(ir)=revc0(ir,ibnd,1)*CMPLX(dvrss(ir),0.0d0,dp)
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
                                &betasave(ia,ih,ir)*&
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

          CALL bfft_orbital_gamma (evc1_new(:,:,1), ibnd, nbnd,.false.)
          !
       ENDDO
       !
#ifdef __MPI
       CALL mp_sum( evc1_new(:,:,1), inter_bgrp_comm )
#endif
       IF(dffts%have_task_groups) DEALLOCATE (tg_dvrss)
       !
       !
       IF( nkb > 0 .and. okvan .and. real_space_debug <= 7) THEN
         !The non real_space part
          CALL dgemm( 'N', 'N', 2*npw_k(1), nbnd, nkb, 1.d0, vkb, &
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
    !   Call h_psi on evc1 such that h.evc1 = sevc1_new
    !
    CALL h_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),sevc1_new(1,1,1))
    !
    ! spsi1 = s*evc1
    !
    IF (real_space_debug > 9 ) THEN
        DO ibnd=1,nbnd,2
         CALL fft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
         CALL s_psir_gamma(ibnd,nbnd)
         CALL bfft_orbital_gamma(spsi1,ibnd,nbnd)
        ENDDO
    ELSE
       CALL s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi1)
    ENDIF
    !
    !   Subtract the eigenvalues
    !
    DO ibnd=1,nbnd
       !
       CALL zaxpy(npw_k(1), CMPLX(-(et(ibnd,1)-scissor), &
            &0.0d0,dp), spsi1(:,ibnd), 1, sevc1_new(:,ibnd,1), 1)
       !
    ENDDO
    !
    IF( nkb > 0 .and. okvan ) DEALLOCATE(becp2)
    !
    RETURN
    !
  END SUBROUTINE lr_apply_liouvillian_gamma
  !
  SUBROUTINE lr_apply_liouvillian_k()
    !
    USE lr_variables,        ONLY : becp1_c
    USE wvfct,               ONLY : current_k, npw
    !
    COMPLEX(kind=dp), ALLOCATABLE :: becp2(:,:)
    !
    IF( nkb > 0 .AND. okvan ) THEN
       !
       ALLOCATE(becp2(nkb,nbnd))
       becp2(:,:)=(0.0d0,0.0d0)
       !
    ENDIF
    !
    !   Now apply to the ground state wavefunctions
    !   and  convert to real space
    !
    IF ( interaction ) THEN
       !
       CALL start_clock('interaction')
       !
       !   evc1_new is used as a container for the interaction
       !
       evc1_new(:,:,:)=(0.0d0,0.0d0)
       !
       DO ik=1,nks
          !
          DO ibnd=1,nbnd
             !
             !   Product with the potential vrs = (vltot+vr)
             !
             DO ir=1,dffts%nnr
                !
                psic(ir)=revc0(ir,ibnd,ik)*dvrssc(ir)
                !
             ENDDO
!             IF (lr_exx) psic(:) = psic(:) + 0.5d0 * revc_int_c(:,ibnd,ik)
             IF (lr_exx) THEN
                CALL lr_exx_apply_revc_int(psic, ibnd, nbnd, ik)
             ENDIF

             !
             !   Back to reciprocal space
             !
             CALL fwfft ('Wave', psic, dffts)
             !
             DO ig=1,npw_k(ik)
                !
                evc1_new(ig,ibnd,ik)=psic(nls(igk_k(ig,ik)))
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
          DO ik=1,nks
             !
             CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
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
             CALL zgemm( 'N', 'N', npw_k(ik), nbnd, nkb, (1.d0,0.d0), vkb, &
                  npwx, becp2, nkb, (1.d0,0.d0), evc1_new(:,:,ik), npwx )
             !
          ENDDO
          !
       ENDIF
       !
    ENDIF
    !
    !   Call h_psi on evc1
    !   h_psi uses arrays igk and npw, so restore those
    !
    IF (lr_exx .AND. .NOT.interaction) CALL lr_exx_kernel_noint(evc1,evc1_new)
    DO ik=1,nks
       !
       CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
       !
       DO ig=1,npw_k(ik)
          !
          g2kin(ig)=((xk(1,ik)+g(1,igk_k(ig,ik)))**2 &
               +(xk(2,ik)+g(2,igk_k(ig,ik)))**2 &
               +(xk(3,ik)+g(3,igk_k(ig,ik)))**2)*tpiba2
          !
       ENDDO
       !
       igk(:)=igk_k(:,ik)
       current_k = ik
       npw=npw_k(ik)
       !
       CALL h_psi(npwx,npw_k(ik),nbnd,evc1(1,1,ik),sevc1_new(1,1,ik))
       !
       CALL s_psi(npwx,npw_k(ik),nbnd,evc1(1,1,ik),spsi1)
       !
       !   Subtract the eigenvalues
       !
       DO ibnd=1,nbnd
          !
          DO ig=1,npw_k(ik)
             !
             sevc1_new(ig,ibnd,ik)=sevc1_new(ig,ibnd,ik) &
                  -cmplx(et(ibnd,ik),0.0d0,dp)*spsi1(ig,ibnd)
             !
          ENDDO
          !
       ENDDO
       !
    ENDDO ! end k loop
    !
    IF( nkb > 0 .and. okvan ) DEALLOCATE(becp2)
    !
    RETURN
  END SUBROUTINE lr_apply_liouvillian_k
  !
END SUBROUTINE lr_apply_liouvillian
!-----------------------------------------------------------------------
