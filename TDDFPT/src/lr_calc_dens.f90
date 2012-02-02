!-----------------------------------------------------------------------
SUBROUTINE lr_calc_dens( evc1, response_calc )
  !---------------------------------------------------------------------
  ! ... calculates response charge density from linear response
  ! ... orbitals and ground state orbitals
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu in 2009
  !
  ! Input : evc1 (qdash etc) Output: rho_1 (=2*sum_v (revc0_v(r) . revc1_v(r,w)  v:valance state index, r denotes a transformation t
  ! in case of uspps becsum is also calculated here
  ! in case of charge response calculation, the rho_tot is calculated here

#include "f_defs.h"
  !
  USE ions_base,                    ONLY : ityp,nat,ntyp=>nsp
  USE cell_base,                    ONLY : omega
  USE ener,                     ONLY : ef
  USE gvecs,                  ONLY : nls,nlsm,doublegrid
  USE fft_base,                 ONLY : dffts, dfftp
  USE fft_interfaces,           ONLY : invfft
  USE io_global,                ONLY : stdout
  USE kinds,                    ONLY : dp
  USE klist,                    ONLY : nks,xk,wk
  USE lr_variables,             ONLY : evc0,revc0,rho_1,lr_verbosity, &
                                       charge_response, itermax,&
                                       cube_save, rho_1_tot,rho_1_tot_im, &
                                       LR_iteration, LR_polarization, &
                                       project,evc0_virt,F,nbnd_total,n_ipol, becp1_virt
  USE lsda_mod,                 ONLY : current_spin, isk
  USE wavefunctions_module,     ONLY : psic
  USE wvfct,                    ONLY : nbnd,et,wg,npwx,npw
  USE control_flags,            ONLY : gamma_only
  USE uspp,                     ONLY : vkb,nkb,okvan,qq,becsum
  USE uspp_param,               ONLY : upf, nh
  USE io_global,                ONLY : ionode, stdout
  USE io_files,                 ONLY : tmp_dir, prefix
  USE mp,                       ONLY : mp_sum
  USE mp_global,                ONLY : inter_pool_comm, intra_pool_comm,nproc
  USE realus,                   ONLY : igk_k,npw_k,addusdens_r
  USE charg_resp,               ONLY : w_T, lr_dump_rho_tot_cube,&
                                       lr_dump_rho_tot_xyzd, &
                                       lr_dump_rho_tot_xcrys,&
                                       resonance_condition,epsil
  USE noncollin_module,     ONLY : nspin_mag
  USE control_flags,         ONLY : tqr
  USE becmod,              ONLY : becp
  !
  IMPLICIT NONE
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  COMPLEX(kind=dp), INTENT(in) :: evc1(npwx,nbnd,nks)
  LOGICAL, INTENT(in) :: response_calc
  !
  ! functions
  real(kind=dp) :: ddot
  !
  ! local variables
  INTEGER :: ir,ik,ibnd,jbnd,ig,ijkb0,np,na,ijh,ih,jh,ikb,jkb,ispin
  INTEGER :: i, j, k, l
  real(kind=dp) :: w1,w2,scal
  real(kind=dp) :: rho_sum!,weight
  real(kind=dp), ALLOCATABLE :: rho_sum_resp_x(:),rho_sum_resp_y(:),rho_sum_resp_z(:) ! These are temporary buffers for response cha

  !complex(kind=dp), external :: ZDOTC
  !complex(kind=dp), allocatable :: spsi(:,:)
  !
  CHARACTER(len=256) :: tempfile, filename
  !
  !OBM DEBUG
  COMPLEX(kind=dp),EXTERNAL :: lr_dot

  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_calc_dens>")')
  ENDIF
  !
  CALL start_clock('lr_calc_dens')
  !
  !allocate(spsi(npwx,nbnd))
  !spsi(:,:)=(0.0d0,0.0d0)
  !
  ALLOCATE( psic(dfftp%nnr) )
  psic(:)=(0.0d0,0.0d0)
  rho_1(:,:)=0.0d0
  !
  IF(gamma_only) THEN
     CALL lr_calc_dens_gamma()
  ELSE
     CALL lr_calc_dens_k()
  ENDIF
  !print *, "rho_1 after lr_calc_dens calculates",SUM(rho_1)
  !print *, "norm of evc1 after lr_calc_dens calculates", lr_dot(evc1(1,1,1),evc1(1,1,1))

  !
  ! ... If a double grid is used, interpolate onto the fine grid
  !
  IF ( doublegrid ) CALL interpolate(rho_1,rho_1,1)
  !
  ! ... Here we add the Ultrasoft contribution to the charge
  !
  !IF ( okvan ) CALL lr_addusdens(rho_1)
  !print *, "rho_1 before addusdens",SUM(rho_1)
  !call start_clock('lrcd_usdens') !TQR makes a huge gain here
  IF(okvan) THEN
   IF (tqr) THEN
    CALL addusdens_r(rho_1,.false.)
   ELSE
    CALL addusdens(rho_1)
   ENDIF
  ENDIF
  DEALLOCATE ( psic )
  !call stop_clock('lrcd_usdens')
  !
  !print *, "rho_1 after addusdens",SUM(rho_1)
#ifdef __MPI
  !call poolreduce(nrxx,rho_1)
  CALL mp_sum(rho_1, inter_pool_comm)
#endif
  !
  ! check response charge density sums to 0
  !call start_clock('lrcd_sp') !Minimal lag, no need to improve
IF (lr_verbosity > 0) THEN

  DO ispin = 1, nspin_mag
   rho_sum=0.0d0
   rho_sum=sum(rho_1(:,ispin))
   !
#ifdef __MPI
   CALL mp_sum(rho_sum, intra_pool_comm )
#endif
   !
   rho_sum=rho_sum*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
   !
   IF(abs(rho_sum)>1.0d-12) THEN
     IF (tqr) THEN
      WRITE(stdout,'(5X, "lr_calc_dens: Charge drift due to real space implementation = " ,1X,e12.5)')&
           rho_sum
      !seems useless
      !rho_sum=rho_sum/(1.0D0*nrxxs)
      !rho_1(:,ispin)=rho_1(:,ispin)-rho_sum
     ELSE
      WRITE(stdout,'(5X,"lr_calc_dens: ****** response charge density does not sum to zero")')
      !
      WRITE(stdout,'(5X,"lr_calc_dens: ****** response charge density =",1X,e12.5)')&
           rho_sum
      !
      WRITE(stdout,'(5X,"lr_calc_dens: ****** response charge density, US part =",1X,e12.5)')&
           scal
      !     call errore(' lr_calc_dens ','Linear response charge density '// &
      !          & 'does not sum to zero',1)
     ENDIF
   ENDIF
  ENDDO
  !
ENDIF
   IF (charge_response == 2 .and. LR_iteration /=0) THEN
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
#ifdef __MPI
     CALL mp_sum(rho_sum_resp_x, intra_pool_comm)
     CALL mp_sum(rho_sum_resp_y, intra_pool_comm)
     CALL mp_sum(rho_sum_resp_z, intra_pool_comm)
     IF (ionode) THEN
#endif
     WRITE(stdout,'(5X,"Dumping plane sums of densities for iteration ",I4)') LR_iteration
     !
     filename = trim(prefix) // ".density_x"
     tempfile = trim(tmp_dir) // trim(filename)
     !
     OPEN (158, file = tempfile, form = 'formatted', status = 'unknown', position = 'append')
     !
     DO i=1,dfftp%nr1
        WRITE(158,*) rho_sum_resp_x(i)
     ENDDO
     !
     CLOSE(158)
     !
     filename = trim(prefix) // ".density_y"
     tempfile = trim(tmp_dir) // trim(filename)
     !
     OPEN (158, file = tempfile, form = 'formatted', status = 'unknown', position = 'append')
     !
     DO i=1,dfftp%nr2
        WRITE(158,*) rho_sum_resp_y(i)
     ENDDO
     !
     CLOSE(158)
     !
     filename = trim(prefix) // ".density_z"
     tempfile = trim(tmp_dir) // trim(filename)
     !
     OPEN (158, file = tempfile, form = 'formatted', status = 'unknown', position = 'append')
     !
     DO i=1,dfftp%nr3
        WRITE(158,*) rho_sum_resp_z(i)
     ENDDO
     !
     CLOSE(158)
     !
#ifdef __MPI
     ENDIF
#endif
     !
     DEALLOCATE( rho_sum_resp_x )
     DEALLOCATE( rho_sum_resp_y )
     DEALLOCATE( rho_sum_resp_z )
     !
  ENDIF
  IF (charge_response == 1 .and. response_calc) THEN
    IF (LR_iteration < itermax) WRITE(stdout,'(5x,"Calculating total response charge density")')
    ! the charge response, it is actually equivalent to an element of
    ! V^T . phi_v where V^T is the is the transpose of the Krylov subspace generated
    ! by the Lanczos algorithm. The total charge density can be written
    ! as
    ! \sum_(lanczos iterations) (V^T.phi_v) . w_T
    ! Where w_T is the corresponding eigenvector from the solution of
    ! (w-L)e_1 = w_T
    !
    ! notice that rho_1 is already reduced across pools above, so no parallelization is necessary
    !
    ! the lr_calc_dens corresponds to q of x only in even iterations
    !
    !print *,"1"
    !print *,"weight",(-1.0d0*AIMAG(w_T(LR_iteration)))
    IF (resonance_condition) THEN
    !singular matrix, the broadening term dominates, phi' has strong imaginary component
     !DO ir=1,nrxx
     ! rho_1_tot_im(ir,:)=rho_1_tot_im(ir,:)+cmplx(rho_1(ir,:),0.0d0,dp)*w_T(LR_iteration)
     !enddo
     CALL zaxpy(dfftp%nnr, w_T(LR_iteration),cmplx(rho_1(:,1),0.0d0,dp),1,rho_1_tot_im(:,1),1) !spin not implemented
    ELSE
    !not at resonance, the imaginary part is neglected ,these are the non-absorbing oscillations
     !DO ir=1,nrxx
     ! rho_1_tot(ir,:)=rho_1_tot(ir,:)+rho_1(ir,:)*dble(w_T(LR_iteration))
     !enddo
     CALL daxpy(dfftp%nnr, dble(w_T(LR_iteration)),rho_1(:,1),1,rho_1_tot(:,1),1) !spin not implemented
    ENDIF
  !
  ENDIF
  !
  !deallocate(spsi)
  !
  CALL stop_clock('lr_calc_dens')
  !
  !call stop_clock('lrcd_sp')
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_calc_dens_gamma
    !
    USE becmod,              ONLY : bec_type, becp, calbec
    USE lr_variables,        ONLY : becp1   !,real_space
    !use real_beta,           only : ccalbecr_gamma, fft_orbital_gamma
    USE io_global,           ONLY : stdout
    USE realus,              ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                           bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, v_loc_psir,&
                           real_space_debug

    !

    DO ibnd=1,nbnd,2
       CALL fft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
       !
       w1=wg(ibnd,1)/omega
       !
       IF(ibnd<nbnd) THEN
          w2=wg(ibnd+1,1)/omega
       ELSE
          w2=w1
       ENDIF
       !call start_clock('lrcd-lp1')
       ! OBM:
       ! (n'(r,w)=2*sum_v (psi_v(r) . q_v(r,w))
       ! where psi are the ground state valance orbitals
       ! and q_v are the standart batch representation (rotated)
       ! response orbitals
       ! Here, since the ith iteration is the best approximate we have
       ! for the most dominant eigenvalues/vectors, an estimate for the response
       ! charge density can be calculated. This is in no way the final
       ! response charge density.
       ! the loop is over real space points.
       DO ir=1,dffts%nnr
          rho_1(ir,1)=rho_1(ir,1) &
               +2.0d0*(w1*real(revc0(ir,ibnd,1),dp)*real(psic(ir),dp)&
               +w2*aimag(revc0(ir,ibnd,1))*aimag(psic(ir)))
       ENDDO
       !
       !call stop_clock('lrcd-lp1')
       ! OBM - psic now contains the response functions at
       ! real space, eagerly putting all the real space stuff at this point.
       ! notice that betapointlist() is called in lr_readin at the very start
       IF ( real_space_debug > 6 .and. okvan) THEN
        ! The rbecp term
        CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
       ENDIF
       ! End of real space stuff
    ENDDO
    !
    ! ... If we have a US pseudopotential we compute here the becsum term
    ! This corresponds to the right hand side of the formula (36) in Ultrasoft paper
    ! be careful about calling lr_calc_dens, as it modifies this globally
    !call start_clock('lrcd-us')
    IF ( okvan ) THEN
       !
       scal = 0.0d0
       becsum(:,:,:) = 0.0d0
       !
       IF ( real_space_debug <= 6) THEN !in real space, the value is calculated above
          !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,evc1,npwx,rbecp,nkb)
          CALL calbec(npw_k(1), vkb, evc1(:,:,1), becp)
       ENDIF
       !
       CALL start_clock( 'becsum' )
       !
       DO ibnd = 1, nbnd
          scal = 0.0d0
          !
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
                              2.d0 * w1 * becp%r(ikb,ibnd) * becp1(ikb,ibnd)
                         scal = scal + qq(ih,ih,np) *1.d0 *  becp%r(ikb,ibnd) * becp1(ikb,ibnd)
                         !
                         ijh = ijh + 1
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            !
                            jkb = ijkb0 + jh
                            !
                            becsum(ijh,na,current_spin) = &
                                 becsum(ijh,na,current_spin) + &
                                 w1 * 2.D0 * (becp1(ikb,ibnd) * becp%r(jkb,ibnd) + &
                                 becp1(jkb,ibnd) * becp%r(ikb,ibnd))
                            scal = scal + qq(ih,jh,np) *1.d0  * (becp%r(ikb,ibnd) * becp1(jkb,ibnd)+&
                                 becp%r(jkb,ibnd) * becp1(ikb,ibnd))
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
           ! OBM debug
           !write(stdout,'(5X,"lr_calc_dens: ibnd,scal=",1X,i3,1X,e12.5)')&
           !     ibnd,scal
       ENDDO
       !
       CALL stop_clock( 'becsum' )
       !
    ENDIF
    !call stop_clock('lrcd-us')
    !
    RETURN
    !
  END SUBROUTINE lr_calc_dens_gamma
  !-----------------------------------------------------------------------
  SUBROUTINE lr_calc_dens_k
    !
    USE becmod,              ONLY : bec_type, becp, calbec
    USE lr_variables,        ONLY : becp1_c
    !
    DO ik=1,nks
       DO ibnd=1,nbnd
          psic(:)=(0.0d0,0.0d0)
          DO ig=1,npw_k(ik)
             psic(nls(igk_k(ig,ik)))=evc1(ig,ibnd,ik)
          ENDDO
          !
          CALL invfft ('Wave', psic, dffts)
          !
          w1=wg(ibnd,ik)/omega
          !
          ! loop over real space points
          DO ir=1,dffts%nnr
             rho_1(ir,:)=rho_1(ir,:) &
                  +2.0d0*w1*real(conjg(revc0(ir,ibnd,ik))*psic(ir),dp)
          ENDDO
          !
       ENDDO
    ENDDO
    !
    ! ... If we have a US pseudopotential we compute here the becsum term
    !
    IF ( okvan ) THEN
       !
       DO ik =1,nks
          !
          CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
          !
          scal = 0.0d0
          becsum(:,:,:) = 0.0d0
          !
          IF ( nkb > 0 .and. okvan ) THEN
             ! call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp,vkb,evc1)
             CALL calbec(npw_k(ik),vkb,evc1(:,:,ik),becp)
          ENDIF
          !
          CALL start_clock( 'becsum' )
          !
          DO ibnd = 1, nbnd
             scal = 0.0d0
             !
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
                                 2.d0 * w1 * real(conjg(becp%k(ikb,ibnd)) * becp1_c(ikb,ibnd,ik))
                            scal = scal + qq(ih,ih,np) *1.d0 *  real(conjg(becp%k(ikb,ibnd)) * becp1_c(ikb,ibnd,ik))
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               becsum(ijh,na,current_spin) = &
                                    becsum(ijh,na,current_spin) + &
                                    w1 * 2.d0 * real(conjg(becp1_c(ikb,ibnd,ik)) * becp%k(jkb,ibnd) + &
                                    becp1_c(jkb,ibnd,ik) * conjg(becp%k(ikb,ibnd)))
                               scal = scal + qq(ih,jh,np) *1.d0  * real(conjg(becp%k(ikb,ibnd)) * becp1_c(jkb,ibnd,ik)+&
                                    becp%k(jkb,ibnd) * conjg(becp1_c(ikb,ibnd,ik)))
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
             ! write(stdout,'(5X,"lr_calc_dens: ibnd,scal=",1X,i3,1X,e12.5)')&
             !      ibnd,scal
          ENDDO
          CALL stop_clock( 'becsum' )
          !
       ENDDO
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_calc_dens_k
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE lr_calc_dens
!-----------------------------------------------------------------------
