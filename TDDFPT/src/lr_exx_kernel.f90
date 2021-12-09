!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------

MODULE lr_exx_kernel
  !------------------------------------------------------------------
  !
  !  This module handles the EXX part of the Liouvillian. It requires
  !  all the correct initialization to have been done for the routines
  !  in PW/src/exx.f90 (if you can't call vexx() then this module
  !  won't work). It respects the divergence treatments etc set in
  !  exx.f90 and also uses the custom grid defined there for gamma_only
  !  calculations.
  !
  !  lr_exx_alloc must be called first to allocate the appropriate 
  !  workspaces, lr_exx_dealloc can be used to dealloc them. Also
  !  lr_exx_revc0_init must be called at the beginning to set up the
  !  EXX operator.
  !  
  !------------------------------------------------------------------
  !

  USE kinds,                  ONLY : DP
  USE lr_variables,           ONLY : evc0, revc0
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : invfft, fwfft
  USE lsda_mod,               ONLY : nspin
  USE wvfct,                  ONLY : nbnd, npwx, wg
  USE gvect,                  ONLY : g, ngm
  USE klist,                  ONLY : xk, wk, nks
  USE lr_variables,           ONLY : gamma_only, lr_verbosity
  USE realus,                 ONLY : invfft_orbital_gamma, fwfft_orbital_gamma,&
                                   & invfft_orbital_k, fwfft_orbital_k
  USE wavefunctions,          ONLY : psic
  USE cell_base,              ONLY : omega
  USE exx_base,               ONLY : g2_convolution
  USE exx,                    ONLY : exxalfa, npwt, gt, dfftt 


  REAL(kind=dp),    PUBLIC, ALLOCATABLE :: revc_int(:,:)
  COMPLEX(kind=dp), PUBLIC, ALLOCATABLE :: revc_int_c(:,:,:)

  ! Workspaces used in the k*d_term functions.
  
  ! Used for the density like \phi(r') \psi(r') products
  COMPLEX(DP), PRIVATE, ALLOCATABLE :: pseudo_dens_c(:)
  ! Used to hold the result of \int \phi(r') \psi(r') 1/|r_12| dr'
  COMPLEX(DP), PRIVATE, ALLOCATABLE :: vhart(:,:)

  ! This holds the groundstate orbitals, in real space on the 
  ! custom EXX grid.
  COMPLEX(DP), PRIVATE, ALLOCATABLE :: red_revc0(:,:,:)

  ! This provides the k -> q point correspondence for the
  ! EXX operator.
  INTEGER, PRIVATE, ALLOCATABLE :: k2q(:)

CONTAINS
  !------------------------------------------------------------------------
  SUBROUTINE lr_exx_restart( set_ace )
     !------------------------------------------------------------------------
     !This SUBROUTINE is called when restarting an exx calculation
     USE xc_lib,    ONLY : xclib_get_exx_fraction, start_exx, &
                           exx_is_active, get_screening_parameter
     
     USE cell_base, ONLY : at
     USE exx_base,  ONLY : exxdiv, erfc_scrlen, exx_divergence, exx_grid_init,&
                           exx_div_check
     ! FIXME: are these variable useful?
     USE exx,       ONLY : fock0, exxenergy2, local_thr, use_ace
     USE exx,       ONLY : exxinit, aceinit, exx_gvec_reinit 

     IMPLICIT NONE
     LOGICAL, INTENT(in) :: set_ace
     !
     CALL exx_grid_init( reinit=.true. )
     CALL exx_gvec_reinit( at )
     CALL exx_div_check()
     !
     use_ace = set_ace
     erfc_scrlen = get_screening_parameter()
     
     exxdiv = exx_divergence()
     exxalfa = xclib_get_exx_fraction()
     CALL start_exx()
     CALL weights()
     ! FIXME: is this useful ?
     IF(local_thr.gt.0.0d0) CALL errore('exx_restart','SCDM with restart NYI',1)
     CALL exxinit(.false.)
     IF (use_ace) CALL aceinit ( DOLOC = .FALSE. )
     ! FIXME: are these variable useful?
     fock0 = exxenergy2()
     !
     RETURN
     !------------------------------------------------------------------------
   END SUBROUTINE lr_exx_restart
  !------------------------------------------------------------------------

  SUBROUTINE lr_exx_alloc()

  USE exx_base,    ONLY : nkqs
  USE klist,       ONLY : nks
  
  IMPLICIT NONE
  INTEGER :: nrxxs

  IF (gamma_only) THEN
     nrxxs= dfftt%nnr
  ELSE
     nrxxs= dffts%nnr
  ENDIF
  !
  ALLOCATE (vhart(nrxxs,nspin))
  ALLOCATE (pseudo_dens_c(nrxxs))
  ALLOCATE (red_revc0(nrxxs,nbnd,nkqs))
  red_revc0 = (0.0_dp, 0.0_dp)
  !
  IF (gamma_only) THEN
     ALLOCATE (revc_int(nrxxs,nbnd))
  ELSE
     ALLOCATE(revc_int_c(nrxxs,nbnd,nks))
     ALLOCATE(k2q(nks)) 
     k2q=0
  ENDIF
  !
END SUBROUTINE lr_exx_alloc
!
SUBROUTINE lr_exx_dealloc()

  IMPLICIT NONE

  DEALLOCATE(pseudo_dens_c, vhart, red_revc0)
  !
  IF (gamma_only) THEN
     DEALLOCATE(revc_int)
  ELSE
     DEALLOCATE(revc_int_c, k2q)
  ENDIF
  !
END SUBROUTINE lr_exx_dealloc
!
SUBROUTINE lr_exx_sum_int()
  USE mp_global, ONLY : inter_bgrp_comm
  USE mp,        ONLY : mp_sum

  CALL start_clock('lr_exx_sum')

  CALL mp_sum(revc_int, inter_bgrp_comm)

  CALL stop_clock('lr_exx_sum')

  RETURN

END SUBROUTINE lr_exx_sum_int
!
SUBROUTINE lr_exx_apply_revc_int(psi, ibnd, nbnd, ik)
  !------------------------------------------------------------------
  !
  !  This routine takes the interaction generated by lr_exx_kernel_int
  !  and applies it to psi for a specified band, ibnd.
  !  It handles the conversion to the smooth grid but assumes that
  !  lr_exx_kernel_int has been called for every band and
  !  lr_exx_sum_int has been called.
  !------------------------------------------------------------------
  !

  USE klist,       ONLY : ngk

  IMPLICIT NONE
  
  COMPLEX(DP), INTENT(INOUT) :: psi(:)
  INTEGER, INTENT(IN) :: ibnd, nbnd, ik
  
  COMPLEX(DP), ALLOCATABLE :: tempphic(:,:),temppsic(:,:), psi_t(:)
  INTEGER :: nrxxs, j
  INTEGER :: npw ! number of plane waves at point k
  
  COMPLEX(DP) :: fp, fm

  CALL start_clock('lr_exx_apply')

  nrxxs= dfftt%nnr

  IF (gamma_only) THEN
     !
     npw = ngk(1)
     !
     ALLOCATE ( tempphic(dffts%nnr,2),temppsic(dffts%nnr,2),&
          & psi_t(dffts%nnr) ) 
     tempphic = (0.0_dp, 0.0_dp)
     temppsic = (0.0_dp, 0.0_dp)
     psi_t    = (0.0_dp, 0.0_dp)

     IF (ibnd < nbnd) THEN
        !
        ! We need to do two bands at a time as per gamma tricks.
        !
        tempphic(1:nrxxs,1)= 0.5d0 * CMPLX( revc_int(1:nrxxs,ibnd),&
             & revc_int(1:nrxxs,ibnd+1), kind=dp ) 
        !
        ! To g-space
        !
        CALL fwfft ('Wave', tempphic(:,1), dfftt)
        !
        ! Now separate the two bands and apply the correct nl mapping
        !
        DO j = 1, npw
           fp = (tempphic(dfftt%nl(j),1) + &
                 tempphic(dfftt%nlm(j),1))*0.5d0 
           fm = (tempphic(dfftt%nl(j),1) - &
                 tempphic(dfftt%nlm(j),1))*0.5d0 
           temppsic( j, 1) = CMPLX( DBLE(fp), AIMAG(fm),kind=DP)
           temppsic( j, 2) = CMPLX(AIMAG(fp),- DBLE(fm),kind=DP)
        ENDDO

        psi_t(dffts%nl(1:npw))= temppsic(1:npw,1)+ (0.0_dp,1.0_dp)&
             &*temppsic(1:npw,2) 
        psi_t(dffts%nlm(1:npw))= CONJG(temppsic(1:npw,1)- (0.0_dp,1.0_dp)&
             &*temppsic(1:npw,2))
     ELSE
        tempphic(1:nrxxs,1) = 0.5d0 * CMPLX( revc_int(1:nrxxs,ibnd),&
             & 0.0_dp, kind=dp ) 
        !
        ! To g-space
        !
        CALL fwfft ('Wave', tempphic(:,1), dfftt)
        !
        ! Correct the nl mapping for the two grids.
        !
        temppsic(1:npw,1)=tempphic(dfftt%nl(1:npw),1)
        psi_t(dffts%nl(1:npw))=temppsic(1:npw,1)
        psi_t(dffts%nlm(1:npw))=CONJG(temppsic(1:npw,1))
        !
     ENDIF
     !
     DEALLOCATE ( tempphic, temppsic )
     !
     CALL invfft ('Wave', psi_t, dffts)
     !
     psi(:) = psi(:) + psi_t(:)    
     !
     DEALLOCATE ( psi_t )
     !
  ELSE
     !
     psi(:) = psi(:) + revc_int_c(:,ibnd,ik)
     !
  ENDIF
  CALL stop_clock('lr_exx_apply')
  RETURN
         
END SUBROUTINE lr_exx_apply_revc_int


SUBROUTINE lr_exx_revc0_init(orbital, ik)
  !------------------------------------------------------------------
  !
  !  This routine should be called at the start of a TDDFT EXX calculation
  !  and it take sthe ground state orbitals (passed into 'orbital'), 
  !  interpolates them onto the custom FFT grid used by exx.f90, and stores
  !  the real-space result in red_revc0().
  !
  !  For the k-points version also it performs the necessary rotations to
  !  obtain the q points from their respective k points.
  !------------------------------------------------------------------
  !
  
  USE mp_global,    ONLY : me_bgrp
  USE exx_base,     ONLY : rir, nkqs, index_sym, index_xk
  USE scatter_mod,  ONLY : gather_grid, scatter_grid
  USE symm_base,    ONLY : sname

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN) :: orbital(:,:,:)
  INTEGER, INTENT(IN) :: ik

  INTEGER :: ibnd, nnr_, ikq, isym, nxxs, nrxxs
  COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:), temppsic(:)

  IF (gamma_only) THEN
     !
     nnr_= dfftt%nnr
     !
     DO ibnd=1,nbnd,2
        !
        CALL invfft_orbital_custom_gamma(orbital(:,:,1), ibnd, nbnd, &
             npwt, dfftt)
        red_revc0(1:nnr_,ibnd,1)=psic(1:nnr_)
        !
     ENDDO
  ELSE
     nxxs=dffts%nr1x * dffts%nr2x * dffts%nr3x
     nrxxs=dffts%nnr
     !
     ALLOCATE(temppsic_all(1:nxxs), psic_all(1:nxxs), temppsic(1:nrxxs))
     !
     DO ibnd=1,nbnd
        !
        CALL invfft_orbital_k ( orbital(:,:,ik), ibnd, nbnd, ik)
        !
        DO ikq=1,nkqs
           !
           IF (index_xk(ikq) .NE. ik) CYCLE
           !
           ! Now rotate k to q
           isym = ABS(index_sym(ikq) )
           IF (index_sym(ikq) > 0) THEN
              IF (TRIM(sname(index_sym(ikq))) .EQ. "identity" ) &
                   & k2q(ik) = ikq
           ENDIF
           !
#if defined(__MPI)
           CALL gather_grid (dffts, psic, psic_all)
           IF ( me_bgrp == 0 ) &
                temppsic_all(1:nxxs) = psic_all(rir(1:nxxs, isym))
           CALL scatter_grid(dffts, temppsic_all, temppsic)
#else
           temppsic(1:nrxxs) = psic(rir(1:nrxxs, isym))
#endif
           IF (index_sym(ikq) < 0 ) &
                &temppsic(1:nrxxs) = CONJG(temppsic(1:nrxxs))
           !
           red_revc0(1:nrxxs, ibnd, ikq) = temppsic(:)
           !
        ENDDO
        !
     ENDDO
     !
     DEALLOCATE(temppsic_all, psic_all)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE lr_exx_revc0_init
!
SUBROUTINE lr_exx_kernel_noint ( evc, int_vect )
  !------------------------------------------------------------------
  !
  !  This routine computes the change of the self consistent potential
  !  due to the perturbation. More specifically it combines the K^1d
  !  and K^2d terms from Ref (1) in the 'non-interacting' (lower SBR
  !  vectors) case. 
  !  (In the hybrid/exx CASE non-interacting is a misnomer as Hartree
  !  like terms are applied even to the lower batch). 
  !
  !   (1)  Rocca, Lu and Galli, J. Chem. Phys., 133, 164109 (2010)
  !------------------------------------------------------------------
  !

  USE kinds,                  ONLY : DP
  USE lr_variables,           ONLY : gamma_only, lr_verbosity
  USE exx_base,               ONLY : g2_convolution, nqs, index_sym, &
                                     index_xkq, index_xk, rir, nkqs
  USE exx,                    ONLY : exxalfa
  USE symm_base,              ONLY : s
  USE cell_base,              ONLY : bg, at
  USE xc_lib,                 ONLY : exx_is_active
  USE io_global,              ONLY : stdout
  USE mp_global,              ONLY : inter_bgrp_comm, ibnd_start, ibnd_end,&
                                   & me_bgrp
  USE mp,                     ONLY : mp_sum
  USE scatter_mod,            ONLY : gather_grid, scatter_grid

  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in)  :: evc(npwx,nbnd,nks)
  COMPLEX(DP), INTENT(out) :: int_vect(npwx,nbnd,nks)
  !
  !
  REAL(kind=dp), ALLOCATABLE    :: revc_int(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: revc_int_c(:,:,:)
  ! Container for the interaction
  !
  INTEGER :: ibnd, ik
  ! Counters
  REAL(kind=dp) :: w1, w2, scale
  !
  REAL(kind=DP), ALLOCATABLE :: fac(:)
  ! Varible to hold to actual interaction (1/|r-r'| or similar in g-space).
  INTEGER :: nxxs, nrxxs
  !
  ! q-points related stuff
  INTEGER :: iq, ikq, isym, ikk
  REAL(DP) :: xk_cryst(3), sxk(3), xkq(3)
  COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:), temppsic(:)
  !
  INTEGER   :: ibnd_start_gamma, ibnd_end_gamma

  LOGICAL   :: exst


  ! Offset ibnd_start and ibnd_ed to avoid conflicts with gamma_tricks for
  ! even values
  !
  CALL start_clock('lr_exx_noint')
  !
  IF (lr_verbosity > 5 ) WRITE(stdout,'("<lr_exx_kernel_noint>")')
  !
  ! Setup the variables that describe the FFT grid in use.
  IF(gamma_only) THEN
     ALLOCATE( fac(dfftt%ngm) )
     nrxxs= dfftt%nnr
  ELSE
     ALLOCATE( fac(ngm) )
     nxxs=dffts%nr1x * dffts%nr2x * dffts%nr3x
     nrxxs = dffts%nnr
     ALLOCATE(temppsic_all(1:nxxs), psic_all(1:nxxs), temppsic(1:nrxxs))
  ENDIF
  !
  fac=0.d0
  !
  IF (exx_is_active()) THEN
     scale = exxalfa
  ELSE
     scale=1.d0
  ENDIF
  !
  !
  IF( gamma_only ) THEN
     !
     ! Put the appropriate interaction in fac(). Note g2_convolution respects
     ! the choice of divergence treatment etc set in the initial PWscf run.
     !
     CALL g2_convolution(dfftt%ngm, gt, xk(:,1), xk(:,1), fac)
     !
     ALLOCATE(revc_int(nrxxs,nbnd))
     !
     revc_int=0.d0
     int_vect=(0.d0,0.d0)
     !
     !
     ! If ibnd_start is even offset it so bands aren't skipped/double counted
     ! when we use gamma_tricks.
     !
     ibnd_start_gamma=ibnd_start
     IF (MOD(ibnd_start, 2)==0) ibnd_start_gamma = ibnd_start+1
     ibnd_end_gamma = MAX(ibnd_end, ibnd_start_gamma)
     !
     DO ibnd=ibnd_start_gamma,ibnd_end_gamma,2
        !
        CALL invfft_orbital_custom_gamma(evc(:,:,1), ibnd, nbnd, &
             npwt, dfftt)
        !
        w1=wg(ibnd,1)/omega
        !
        IF (ibnd<nbnd) THEN
           w2=wg(ibnd+1,1)/omega
        ELSE
           w2=0.0d0
        ENDIF
        ! Update the container with the actual interaction for this band(s).
        revc_int(1:dfftt%nnr,:)= revc_int(1:dfftt%nnr,:) &
             & -1.d0 * scale * k1d_term_gamma(w1,w2,psic,fac,ibnd,evc(:,:,1)) & 
             & +1.d0 * scale * k2d_vexx_term_gamma(w1,w2,psic,fac,ibnd,.false.)
        !
     ENDDO
     !
     CALL mp_sum( revc_int(:,:), inter_bgrp_comm )
     !
     ! Put the constructed interaction into the psic workspace and
     ! then on to int_vec. 
     !
     DO ibnd=ibnd_start_gamma,ibnd_end_gamma,2!1,nbnd,2
        !
        psic=(0.0_dp, 0.0_dp)
        !
        IF (ibnd<nbnd) psic(1:nrxxs)=CMPLX(revc_int(1:nrxxs,ibnd),&
             & revc_int(1:nrxxs,ibnd+1),dp)
        IF (ibnd==nbnd) psic(1:nrxxs)=CMPLX(revc_int(1:nrxxs,ibnd)&
             &,0.d0,dp)
        !
        CALL fwfft_orbital_custom_gamma (int_vect(:,:,1), ibnd, nbnd, &
             npwt, dfftt) 
        !
     ENDDO
     !
     CALL mp_sum( int_vect(:,:,1), inter_bgrp_comm )
     !
     DEALLOCATE(revc_int, fac)
     !
     int_vect=int_vect*0.5d0 
     !
  ELSE
     !
     ! In the case of k-points...
     !
     ALLOCATE(revc_int_c(dffts%nnr,nbnd,nks))
     !
     int_vect=(0.d0,0.d0)
     revc_int_c=(0.d0,0.d0)
     !
     ! Loop over all k-points
     !
     DO ikk=1,nks
        !
        DO ibnd=1,nbnd,1
           !
           !
           CALL invfft_orbital_k (evc(:,:,ikk), ibnd, nbnd, ikk)
           !
#if defined(__MPI)
           psic_all=(0.d0,0.d0)
           CALL gather_grid(dffts, psic, psic_all)
#endif
           !
           ! Now psic_all has the collected band ibnd at kpoint ikk.
           ! We now search to see which q points are identical via symmetry
           !
           DO ikq=1,nkqs
              !
              IF (ikk /= index_xk(ikq)) CYCLE
              !
              ! Now rotate k to q, scatter and put in temppsic
              !
              isym = ABS(index_sym(ikq))
#if defined(__MPI)
              IF ( me_bgrp == 0 ) &
                   temppsic_all(1:nxxs) = psic_all(rir(1:nxxs, isym))
              CALL scatter_grid(dffts, temppsic_all, temppsic)
#else

              temppsic(1:nrxxs) = psic(rir(1:nrxxs, isym))
#endif
              IF (index_sym(ikq) < 0 ) &
                   &temppsic(1:nrxxs) = CONJG(temppsic(1:nrxxs))
              !
              ! We now search over all k+q to find which correspond
              ! to the k+q point in temppsic
              !
              DO ik=1,nks
                 DO iq=1,nqs
                    ! 
                    !
                    IF (ikq /= index_xkq(ik,iq)) CYCLE
                    !
                    xk_cryst(:) = at(1,:)*xk(1,ikk) + at(2,:)*xk(2,ikk) &
                         + at(3,:)*xk(3,ikk)
                    !
                    IF (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
                    !
                    sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                         s(:,2,isym)*xk_cryst(2) + &
                         s(:,3,isym)*xk_cryst(3) 
                    xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)
                    !
                    CALL g2_convolution(ngm, g, xk(:,ik), xkq, fac)
                    !
                    !
                    w1 = wg(ibnd,ik)/(wk(ik) * nqs)
                    !
                    ! Update the container with the actual interaction for this band(s).
                    revc_int_c(:,:, ik)= revc_int_c(:,:, ik) &
                         &-1.d0 * scale * &
                         &k1d_term_k(w1,temppsic,fac,ibnd,ik,ikq) &
                         &+1.d0 * scale * &
                         &k2d_term_k(w1,temppsic,fac,ibnd,ik,ikq)
                 END DO
              END DO
           ENDDO
        ENDDO
     ENDDO
     !
     !
     CALL mp_sum( revc_int_c, inter_bgrp_comm )
     !
     DO ik=1,nks
        !
        DO ibnd=1,nbnd,1
           !
           psic(:)=(0.0_dp,0.0_dp)
           !
           ! Put the constructed interaction into the psic workspace
           ! and then on to int_vec. 
           !
           psic(:)=revc_int_c(:,ibnd,ik)
           !
           CALL fwfft_orbital_k (int_vect(:,:,ik), ibnd, nbnd, ik)
           !
        ENDDO
        !
     ENDDO
     !
     DEALLOCATE(revc_int_c, fac)
     DEALLOCATE(temppsic_all, psic_all, temppsic)
     !
  ENDIF
  !
  !
  CALL stop_clock('lr_exx_noint')
  RETURN
  !
END SUBROUTINE lr_exx_kernel_noint
!
SUBROUTINE lr_exx_kernel_int ( orbital, ibnd, nbnd, ikk )
  !------------------------------------------------------------------
  !
  !  This routine computes the change of the self consistent potential
  !  due to the perturbation. More specifically it combines the K^1d
  !  and K^2d terms from Ref (1) in the 'interacting' (upper SBR
  !  vectors) case.
  !
  !  (1)  Rocca, Lu and Galli, J. Chem. Phys., 133, 164109 (2010)
  !
  !  Unlike lr_exx_kernel_noint above this routine should be called
  !  a band at a time. Also it does sum the interaction over the band
  !  group (obviously) and therefore does not place the full interaction
  !  in the revc_int container. See lr_bse_sum and lr_bse_apply_revc_int.
  !------------------------------------------------------------------
  !

  USE kinds,                  ONLY : DP
  USE klist,                  ONLY : xk
  USE io_global,              ONLY : stdout  
  USE exx,                    ONLY : exxalfa
  USE exx_base,               ONLY : g2_convolution, nqs, index_sym, &
                                     index_xkq, index_xk, rir, nkqs
  USE symm_base,              ONLY : s
  USE cell_base,              ONLY : bg, at
  USE xc_lib,                 ONLY : exx_is_active
  USE mp_global,              ONLY : me_bgrp
  USE scatter_mod,            ONLY : gather_grid, scatter_grid
  USE lr_variables,           ONLY : ltammd

  IMPLICIT NONE

  COMPLEX(DP), INTENT(in) :: orbital(:,:)
  INTEGER, INTENT(in)     :: ibnd   ! Current band under consideration
  INTEGER, INTENT(in)     :: ikk    ! K-point of the current band
  INTEGER, INTENT(in)     :: nbnd   ! Total number of bands


  REAL(kind=DP), ALLOCATABLE :: fac(:)
  ! Varible to hold to actual interaction (1/|r-r'| or similar in g-space).

  REAL(kind=DP) :: scale
  ! Variable to hold exxalfa
  ! (or otherwise if interaction is used whilst exx_is_active == .false.)

  INTEGER :: nxxs, nrxxs
  REAL(kind=dp) :: w1, w2
  ! Weights for the band(s)

  ! q-points related stuff
  INTEGER :: iq, ikq, isym, ik
  REAL(DP) :: xk_cryst(3), sxk(3), xkq(3)
  COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:), temppsic(:)

  CALL start_clock('lr_exx_int')
  IF (lr_verbosity > 5 ) WRITE(stdout,'("<lr_exx_kernel_int>")')

  IF(gamma_only) THEN
     ALLOCATE( fac(dfftt%ngm) )
     nrxxs= dfftt%nnr
  ELSE
     ALLOCATE( fac(ngm) )
     nxxs=dffts%nr1x * dffts%nr2x * dffts%nr3x
     nrxxs = dffts%nnr
     ALLOCATE(temppsic_all(1:nxxs), psic_all(1:nxxs), temppsic(1:nrxxs))
  ENDIF
  !
  fac=0.d0
  !
  IF (exx_is_active()) THEN
     scale = exxalfa
  ELSE
     scale=1.d0
  ENDIF
  !
  IF( gamma_only ) THEN
     !
     CALL invfft_orbital_custom_gamma( orbital, ibnd, nbnd, npwt, dfftt )
     !
     w1=wg(ibnd,1)/omega
     !
     IF (ibnd<nbnd) THEN
        w2=wg(ibnd+1,1)/omega
     ELSE
        w2=0.0d0
     ENDIF
     !
     CALL g2_convolution(dfftt%ngm, gt,xk(:,1), xk(:,1), fac)
     !
     IF (.NOT.ltammd) THEN
        revc_int(1:dfftt%nnr,:)= revc_int(1:dfftt%nnr,:)&
             & -1.d0 * scale * k1d_term_gamma(w1,w2,psic,fac,ibnd,orbital) &
             & -1.d0 * scale * k2d_vexx_term_gamma(w1,w2,psic,fac,ibnd,.true.)
     ELSE
        !
        ! Slightly different interaction in the Tamm--Dancoff case.
        ! Note the factor of 2 to account for the fact lr_apply_liovillian
        ! scales the *whole interaction term* by a factor of 0.5 in the TD
        ! case.
        !
        revc_int(1:dfftt%nnr,:)= revc_int(1:dfftt%nnr,:)&
             & -2.d0 * scale * k1d_term_gamma(w1,w2,psic,fac,ibnd,orbital)
     ENDIF
     !
  ELSE
     !
     CALL invfft_orbital_k (orbital(:,:), ibnd, nbnd, ikk)
     !
#if defined(__MPI)
     CALL gather_grid(dffts, psic, psic_all)
#endif
     !
     ! Now psic_all has the collected band ibnd at kpoint ikk.
     ! We now search to see which q points are identical via symmetry
     !
     DO ikq=1,nkqs
        !
        IF (ikk /= index_xk(ikq)) CYCLE
        !
        ! Now rotate k to q, scatter and put in temppsic
        !
        isym = ABS(index_sym(ikq))
#if defined(__MPI)
        IF ( me_bgrp == 0 ) &
             temppsic_all(1:nxxs) = psic_all(rir(1:nxxs, isym))
        CALL scatter_grid(dffts, temppsic_all, temppsic)
#else
        
        temppsic(1:nrxxs) = psic(rir(1:nrxxs, isym))
#endif
        IF (index_sym(ikq) < 0 ) &
                   &temppsic(1:nrxxs) = CONJG(temppsic(1:nrxxs))
        !
        ! We now search over all k+q to find which correspond
        ! to the k+q point in temppsic
        !
        DO ik=1,nks
           DO iq=1,nqs
              ! 
              !
              IF (ikq /= index_xkq(ik,iq)) CYCLE
              !
              xk_cryst(:) = at(1,:)*xk(1,ikk) + at(2,:)*xk(2,ikk) &
                   + at(3,:)*xk(3,ikk)
              !
              IF (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
              !
              sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                   s(:,2,isym)*xk_cryst(2) + &
                   s(:,3,isym)*xk_cryst(3) 
              xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)
              !
              CALL g2_convolution(ngm, g, xk(:,ik), xkq, fac)
              !
              !
              w1=wg(ibnd,ik)/(wk(ik) * nqs)
              !
              IF(.NOT.ltammd) THEN
                 revc_int_c(:,:,ik)= revc_int_c(:,:,ik) & 
                      & -1.d0 * scale * &
                      &k1d_term_k(w1,temppsic,fac,ibnd,ik,ikq) &
                      & -1.d0 * scale * &
                      &k2d_term_k(w1,temppsic,fac,ibnd,ik,ikq)
              ELSE
                 !
                 ! Slightly different interaction in the Tamm--Dancoff case.
                 ! Note the factor of 2 to account for the fact
                 ! lr_apply_liovillian scales the *whole interaction term*
                 ! by a factor of 0.5 in the TD case.
                 !
                 revc_int_c(:,:,ik)= revc_int_c(:,:,ik) & 
                      & -2.d0 * scale * &
                      &k1d_term_k(w1,temppsic,fac,ibnd,ik,ikq) 
              ENDIF
              !
           ENDDO
        ENDDO
     ENDDO
     !
     DEALLOCATE(temppsic_all, psic_all, temppsic)
     !
  ENDIF
  IF (lr_verbosity > 5 ) WRITE(stdout,'("<lr_exx_kernel_int> exit")')
  CALL stop_clock('lr_exx_int')

  ! Note at this point revc_int(_c) may be incomplete
  ! (depending if all the bands have been processed) and still needs to be summed
  ! and converted to the smooth FFT grid.

  RETURN
  
END SUBROUTINE lr_exx_kernel_int
!
FUNCTION k1d_term_gamma(w1, w2, psi, fac_in, ibnd, orbital) RESULT (psi_int)
  !------------------------------------------------------------------
  !
  !   This routine computes the  K^1d term, Eq (21) from Eq Ref (1).
  !   As the integral in this term remains the same throughout the 
  !   chain it can also be calculated once for each combination of
  !   bands and stored in k1d_pot().
  !
  !   psi contains two bands of |a> with w1, w2 the associated weights
  !   fac_in contains the interaction W(G) and ibnd the band index n'.
  !
  !   in gamma, bands are real, so is not necessary computes vhart 
  !   integrals for each couple of bands (integral in eq(21), 
  !   because v, v' is the same of v', v, so only nbnd(nbnd+1)/2 
  !   couple are computed    
  !
  !   (1)  Rocca, Lu and Galli, J. Chem. Phys., 133, 164109 (2010)
  !
  !------------------------------------------------------------------
  !

  IMPLICIT NONE
  
  COMPLEX(KIND=DP),DIMENSION(:), INTENT(IN)  :: psi
  REAL(KIND=DP),DIMENSION(:), INTENT(IN)     :: fac_in
  REAL(kind=DP), INTENT(IN)  :: w1, w2
  REAL(KIND=DP), ALLOCATABLE :: psi_int(:,:)
  INTEGER, INTENT(IN) :: ibnd
  COMPLEX(DP), INTENT(IN) :: orbital(:,:)
  COMPLEX(DP) :: psitemp(dfftt%nnr)
  !
  ! Workspaces
  !
  INTEGER                  :: ibnd2, is, npw_, ngm_, nnr_
  INTEGER                  :: nrec
  !
  npw_=npwt
  ngm_=dfftt%ngm
  nnr_=dfftt%nnr
  !  
  ALLOCATE(psi_int(nnr_, nbnd))
  psi_int = 0.d0
  !
  !
  IF (ibnd < nbnd) then
     !
     vhart(:,:) = (0.d0,0.d0)
     pseudo_dens_c(:) = (0.d0,0.d0)
     !
     ! calculate vhart for couples ibnd,ibnd and ibnd+1,ibnd
     !
     pseudo_dens_c(1:nnr_) = CMPLX( w1*DBLE(red_revc0(1:nnr_,ibnd,1)) *&
          & DBLE(red_revc0(1:nnr_,ibnd,1)), &
          & w2*AIMAG(red_revc0(1:nnr_,ibnd, 1)) *&
          & DBLE(red_revc0(1:nnr_,ibnd,1)), kind=DP )
     !
     CALL fwfft ('Rho', pseudo_dens_c, dfftt)
     !
     DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) =&
                & pseudo_dens_c(dfftt%nl(1:ngm_)) *&
                & fac_in(1:ngm_)
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_)
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
     ENDDO
     !
     ! Finally return the interaction psi_int for terms:
     ! ibnd,   ibnd,   ibnd 
     ! ibnd+1, ibnd+1, ibnd 
     ! ibnd,   ibnd,   ibnd+1
     !
     psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
          & + DBLE(vhart(1:nnr_, 1)) * DBLE(psi(1:nnr_)) &
          & + AIMAG(vhart(1:nnr_,1)) * AIMAG(psi(1:nnr_))
     !
     psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
           & + AIMAG(vhart(1:nnr_,1)) * DBLE(psi(1:nnr_))
     !
     ! calculate vhart for couple ibnd+1,ibnd+1
     !
     vhart(:,:) = (0.d0,0.d0)
     pseudo_dens_c(:) = (0.d0,0.d0)
     !
     pseudo_dens_c(1:nnr_) = CMPLX( w2*AIMAG(red_revc0(1:nnr_,ibnd,1)) *&
          & AIMAG(red_revc0(1:nnr_,ibnd,1)), 0.0d0, kind=DP )
     !
     CALL fwfft ('Rho', pseudo_dens_c, dfftt)
     !
     DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) =&
                & pseudo_dens_c(dfftt%nl(1:ngm_)) *&
                & fac_in(1:ngm_)
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_)
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
     ENDDO
     !
     ! Finally return the interaction psi_int for term:
     ! ibnd+1, ibnd+1, ibnd+1
     ! 
     psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
           & + DBLE(vhart(1:nnr_,1)) * AIMAG(psi(1:nnr_))
     !
     !
     ! start second loop over bands
     !
     DO ibnd2=1,ibnd-1,1  
        !
        ! calculate vhart for couples ibnd,ibnd2 and ibnd+1,ibnd2
        !
        vhart(:,:) = (0.d0,0.d0)
        pseudo_dens_c(:) = (0.d0,0.d0)
        !
        ! The following code is to check if the value of ibnd2 is even or odd
        ! and therefore whether the REAL or IMAGINARY part of red_revc0 is to be
        ! used. This is because red_revc0 is stored using gamma_tricks.
        !
        IF (MOD(ibnd2,2)==1) THEN
           pseudo_dens_c(1:nnr_) = CMPLX( w1*DBLE(red_revc0(1:nnr_,ibnd,1)) *&
                & DBLE(red_revc0(1:nnr_,ibnd2,1)), &
                & w2*AIMAG(red_revc0(1:nnr_,ibnd, 1)) *&
                & DBLE(red_revc0(1:nnr_,ibnd2,1)), kind=DP )
        ELSE
           pseudo_dens_c(1:nnr_) = CMPLX( w1*DBLE(red_revc0(1:nnr_,ibnd,1)) *&
                &AIMAG(red_revc0(1:nnr_,ibnd2-1,1)),&
                &w2*AIMAG(red_revc0(1:nnr_,ibnd,1)) *&
                &AIMAG(red_revc0(1:nnr_,ibnd2-1,1)), kind=DP )
        ENDIF
        !
        CALL fwfft ('Rho', pseudo_dens_c, dfftt)
        !
        ! hartree contribution is computed in reciprocal space
        !
        DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) =&
                & pseudo_dens_c(dfftt%nl(1:ngm_)) *&
                & fac_in(1:ngm_) 
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_) 
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
           !
        ENDDO
        !
        ! Finally return the interaction psi_int for terms:
        ! ibnd,   ibnd,   ibnd2
        ! ibnd+1, ibnd+1, ibnd2
        ! ibnd2,  ibnd2,  ibnd
        ! ibnd2,  ibnd2,  ibnd+1
        !
        psi_int(1:nnr_,ibnd2) = psi_int(1:nnr_,ibnd2) &
             & + DBLE(vhart(1:nnr_, 1)) * DBLE(psi(1:nnr_)) &
             & + AIMAG(vhart(1:nnr_,1)) * AIMAG(psi(1:nnr_))
        !
        CALL invfft_orbital_ibnd2_gamma(orbital(:,:), psitemp, ibnd2, npw_, dfftt)
        !
        !
        psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
             & + DBLE(vhart(1:nnr_, 1)) * DBLE(psitemp(1:nnr_))
        !
        psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
             & + AIMAG(vhart(1:nnr_, 1)) * DBLE(psitemp(1:nnr_))
        !
     ENDDO
     !
  ELSE
     !
     ! calculate vhart for couple ibnd,ibnd 
     !
     vhart(:,:) = (0.d0,0.d0)
     pseudo_dens_c(:) = (0.d0,0.d0)
     !
     pseudo_dens_c(1:nnr_) = CMPLX( w1*DBLE(red_revc0(1:nnr_,ibnd,1)) *&
          & DBLE(red_revc0(1:nnr_,ibnd,1)), 0.0d0, kind=DP )
     !
     CALL fwfft ('Rho', pseudo_dens_c, dfftt)
     !
     DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) =&
                & pseudo_dens_c(dfftt%nl(1:ngm_)) *&
                & fac_in(1:ngm_)
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_)
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
     ENDDO
     !
     ! Finally return the interaction psi_int for term:
     ! ibnd,   ibnd,   ibnd
     ! 
     psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
           & + DBLE(vhart(1:nnr_,1)) * DBLE(psi(1:nnr_))
     !
     ! start second loop over bands
     !
     DO ibnd2=1,ibnd-1,1
        !
        !
        ! Set up the pseudo density for this Hartree like interaction.
        ! calculate vhart for couple ibnd,ibnd2
        !
        vhart(:,:) = (0.d0,0.d0)
        pseudo_dens_c(:) = (0.d0,0.d0)
        !
        ! The following code is to check if the value of ibnd2 is even or odd
        ! and therefore whether the REAL or IMAGINARY part of red_revc0 is to be
        ! used. This is because red_revc0 is stored using gamma_tricks.
        !
        IF (MOD(ibnd2,2)==1) THEN
           pseudo_dens_c(1:nnr_) = CMPLX( w1*DBLE(red_revc0(1:nnr_,ibnd,1)) *&
                & DBLE(red_revc0(1:nnr_,ibnd2,1)), 0.0d0, kind=DP )
        ELSE
           pseudo_dens_c(1:nnr_) = CMPLX( w1*DBLE(red_revc0(1:nnr_,ibnd,1)) *&
                & AIMAG(red_revc0(1:nnr_,ibnd2-1,1)), 0.0d0, kind=DP )
        ENDIF
        !
        CALL fwfft ('Rho', pseudo_dens_c, dfftt)
        !
        ! hartree contribution is computed in reciprocal space
        !
        DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) =&
                & pseudo_dens_c(dfftt%nl(1:ngm_)) *&
                & fac_in(1:ngm_)
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_)
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
           !
        ENDDO
        !
        ! Finally return the interaction psi_int for terms:
        ! ibnd, ibnd,  ibnd2
        ! ibn2, ibnd2, ibnd
        !
        psi_int(1:nnr_,ibnd2) = psi_int(1:nnr_,ibnd2) &
             & + DBLE(vhart(1:nnr_, 1)) * DBLE(psi(1:nnr_))
        !
        CALL invfft_orbital_ibnd2_gamma(orbital(:,:), psitemp, ibnd2, npw_, dfftt)
        !
        !
        psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
             & + DBLE(vhart(1:nnr_, 1)) * DBLE(psitemp(1:nnr_))
        !
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END FUNCTION k1d_term_gamma
!
FUNCTION k1d_term_k(w1, psi, fac_in, ibnd, ik,ikq) RESULT (psi_int)
  !------------------------------------------------------------------
  !
  !   This routine computes the  K^1d term, Eq (21) from Eq Ref (1).
  !   As the integral in this term remains the same throughout the 
  !   chain it can also be calculated once for each combination of
  !   bands and stored in k1d_pot().
  !
  !   psi contains two bands of |a> with w1, w2 the associated weights
  !   fac_in contains the interaction W(G) and ibnd the band index n'.
  !
  !   (1)  Rocca, Lu and Galli, J. Chem. Phys., 133, 164109 (2010)
  !------------------------------------------------------------------
  !
    
  IMPLICIT NONE
  
  COMPLEX(KIND=DP),DIMENSION(:), INTENT(IN)  :: psi
  REAL(KIND=DP),DIMENSION(:), INTENT(IN)     :: fac_in
  REAL(kind=DP), INTENT(IN)  :: w1
  COMPLEX(KIND=DP), ALLOCATABLE :: psi_int(:,:)
  INTEGER, INTENT(IN) :: ibnd, ik,ikq
  !
  ! Workspaces
  !
  INTEGER                  :: ibnd2, is
  
  ALLOCATE(psi_int(dffts%nnr, nbnd))
  psi_int = (0.d0,0.d0)
  !
  DO ibnd2=1,nbnd,1
     !
     !
     ! Set up the pseudo density for this Hartree like interaction.
     !
     vhart(:,:) = (0.d0,0.d0)
     pseudo_dens_c(:) = (0.d0,0.d0)
     !
     pseudo_dens_c(:) = CONJG(red_revc0(:,ibnd,ikq))*&
          &red_revc0(:,ibnd2,k2q(ik))/omega
     !
     CALL fwfft ('Rho', pseudo_dens_c, dffts)
     !
     ! hartree contribution is computed in reciprocal space
     !
     DO is = 1, nspin
        !
        vhart(dffts%nl(1:ngm),is) = w1*pseudo_dens_c(dffts%nl(1:ngm)) *&
             & fac_in(1:ngm) 
        !
        !  and transformed back to real space
        !
        CALL invfft ('Rho', vhart (:, is), dffts)
        !
     ENDDO
     !
     ! Finally return the interaction
     !
     psi_int(:,ibnd2) = psi_int(:,ibnd2) + vhart(:,1) * psi(:)
     !
  ENDDO
  !
END FUNCTION k1d_term_k
!
FUNCTION k2d_vexx_term_gamma(w1, w2, psi, fac_in, ibnd, interaction) RESULT (psi_int)
  !------------------------------------------------------------------
  !
  !   This routine computes the  K^2d term, Eq (22) from Eq Ref (1).
  !   psi contains two bands of |b> with w1, w2 the associated weights
  !   fac_in contains the interaction W(G) and ibnd the band index n'.
  !
  !   Also computes the Vexx term, Eq (15) from Eq Ref (2). 
  !   So in gamma calculations, h_psi soubrite doesn't call anymore 
  !   vexx routine 
  !
  !   (1)  Rocca, Lu and Galli, J. Chem. Phys., 133, 164109 (2010)
  !   (2)  Ge, Binnie, Rocca, Gebauer, Baroni, Comp. Phys. Comm., 
  !        185, 2080 (2014)
  !------------------------------------------------------------------
  !

  IMPLICIT NONE
  
  COMPLEX(KIND=DP),DIMENSION(:), INTENT(IN)  :: psi
  REAL(KIND=DP),DIMENSION(:), INTENT(IN)     :: fac_in
  REAL(kind=DP), INTENT(IN)  :: w1, w2
  INTEGER, INTENT(IN) :: ibnd
  REAL(KIND=DP), ALLOCATABLE :: psi_int(:,:)
  LOGICAL, INTENT(IN) :: interaction
  !
  ! Workspaces
  !
  INTEGER                  :: ibnd2, is, npw_,ngm_, nnr_
  !
  nnr_ = dfftt%nnr
  ngm_ = dfftt%ngm
  npw_ = npwt
  !
  ALLOCATE(psi_int(nnr_, nbnd))
  psi_int = 0.d0
  !
  !
  DO ibnd2=1,nbnd,1
     !
     ! Set up the pseudo density for this Hartree like interaction.
     !
     vhart(:,:) = (0.d0,0.d0)
     pseudo_dens_c(:) = (0.d0,0.d0)
     !
     ! The following code is to check if the value of ibnd2 is even or odd
     ! and therefore whether the REAL or IMAGINARY part of red_revc0 is to be
     ! used. This is because red_revc0 is stored using gamma_tricks.
     !
     IF (MOD(ibnd2,2)==1) THEN
        pseudo_dens_c(1:nnr_) = CMPLX( &
             & w1*DBLE(psi(1:nnr_))*DBLE(red_revc0(1:nnr_,ibnd2,1)),&
             & w2*AIMAG(psi(1:nnr_))*DBLE(red_revc0(1:nnr_,ibnd2,1)), kind=DP )
        !
        CALL fwfft ('Rho', pseudo_dens_c, dfftt)
        !
        ! hartree contribution is computed in reciprocal space
        !
        DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) = pseudo_dens_c(dfftt%nl(1:ngm_)) * fac_in(1:ngm_)
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_)
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
           !
        ENDDO
        !
        !
        ! Finally return the interaction
        !
        psi_int(1:nnr_,ibnd2) = psi_int(1:nnr_,ibnd2) &
             & +DBLE(vhart(1:nnr_,1)) * DBLE(red_revc0(1:nnr_,ibnd,1)) &
             & +AIMAG(vhart(1:nnr_,1)) * AIMAG(red_revc0(1:nnr_,ibnd,1))
        !
        IF (interaction) then 
           psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
                & +DBLE(vhart(1:nnr_,1)) * DBLE(red_revc0(1:nnr_,ibnd2,1))
           IF (ibnd < nbnd) then
              psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
                   & +AIMAG(vhart(1:nnr_,1)) * DBLE(red_revc0(1:nnr_,ibnd2,1))
           ENDIF
        ELSE
           psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
                & -DBLE(vhart(1:nnr_,1)) * DBLE(red_revc0(1:nnr_,ibnd2,1)) 
           IF (ibnd < nbnd) then
              psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
                   & -AIMAG(vhart(1:nnr_,1)) * DBLE(red_revc0(1:nnr_,ibnd2,1))
           ENDIF
        ENDIF
        !
     ELSE
        pseudo_dens_c(1:nnr_) = CMPLX( &
             & w1*DBLE(psi(1:nnr_))*AIMAG(red_revc0(1:nnr_,ibnd2-1,1)),&
             & w2*AIMAG(psi(1:nnr_))*AIMAG(red_revc0(1:nnr_,ibnd2-1,1)), kind=DP )
        !
        CALL fwfft ('Rho', pseudo_dens_c, dfftt)
        !
        ! hartree contribution is computed in reciprocal space
        !
        DO is = 1, nspin
           !
           vhart(dfftt%nl(1:ngm_),is) = pseudo_dens_c(dfftt%nl(1:ngm_)) * fac_in(1:ngm_)
           IF (gamma_only) vhart(dfftt%nlm(1:ngm_),is) = &
                & pseudo_dens_c(dfftt%nlm(1:ngm_)) *&
                & fac_in(1:ngm_)
           !
           !  and transformed back to real space
           !
           CALL invfft ('Rho', vhart (:, is), dfftt)
           !
        ENDDO
        !
        !
        ! Finally return the interaction
        !
        psi_int(1:nnr_,ibnd2) = psi_int(1:nnr_,ibnd2) &
             & +DBLE(vhart(1:nnr_,1)) * DBLE(red_revc0(1:nnr_,ibnd,1)) &
             & +AIMAG(vhart(1:nnr_,1)) * AIMAG(red_revc0(1:nnr_,ibnd,1))
        !
        ! 
        ! and interaction for Vexx term
        !
        IF (interaction) then
           psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
                & +DBLE(vhart(1:nnr_,1)) * AIMAG(red_revc0(1:nnr_,ibnd2-1,1)) 
           IF (ibnd < nbnd) then
              psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
                   & +AIMAG(vhart(1:nnr_,1)) * AIMAG(red_revc0(1:nnr_,ibnd2-1,1))
           ENDIF
        ELSE
           psi_int(1:nnr_,ibnd) = psi_int(1:nnr_,ibnd) &
                & -DBLE(vhart(1:nnr_,1)) * AIMAG(red_revc0(1:nnr_,ibnd2-1,1))
           IF (ibnd < nbnd) then
              psi_int(1:nnr_,ibnd+1) = psi_int(1:nnr_,ibnd+1) &
                   & -AIMAG(vhart(1:nnr_,1)) * AIMAG(red_revc0(1:nnr_,ibnd2-1,1))
           ENDIF
        ENDIF
        !
     ENDIF
     !
     !
  ENDDO
  !
  RETURN
  !
END FUNCTION k2d_vexx_term_gamma
!
FUNCTION k2d_term_k(w1, psi, fac_in, ibnd, ik, ikq) RESULT (psi_int)
  !------------------------------------------------------------------
  !
  !   This routine computes the  K^2d term, Eq (22) from Eq Ref (1).
  !   psi contains two bands of |b> with w1, w2 the associated weights
  !   fac_in contains the interaction W(G) and ibnd the band index n'.
  !
  !   (1)  Rocca, Lu and Galli, J. Chem. Phys., 133, 164109 (2010)
  !------------------------------------------------------------------
  !
  
  IMPLICIT NONE
  
  COMPLEX(KIND=DP),DIMENSION(:), INTENT(IN)  :: psi
  REAL(KIND=DP),DIMENSION(:), INTENT(IN)     :: fac_in
  REAL(kind=DP), INTENT(IN)  :: w1
  INTEGER, INTENT(IN) :: ibnd, ik, ikq
  COMPLEX(KIND=DP), ALLOCATABLE :: psi_int(:,:)
  !
  ! Workspaces
  !
  INTEGER                  :: ibnd2, is
  !
  ALLOCATE(psi_int(dffts%nnr, nbnd))
  psi_int = (0.d0,0.d0)
  !
  DO ibnd2=1,nbnd,1
     !
     ! Set up the pseudo density for this Hartree like interaction.
     !
     vhart(:,:) = (0.d0,0.d0)
     pseudo_dens_c(:) = (0.d0,0.d0)
     !
     pseudo_dens_c(:) =  CONJG(psi(:))*red_revc0(:,ibnd2,k2q(ik))/omega
     !
     CALL fwfft ('Rho', pseudo_dens_c, dffts)
     !
     ! hartree contribution is computed in reciprocal space
     !
     DO is = 1, nspin
        !
        vhart(dffts%nl(1:ngm),is) = w1*pseudo_dens_c(dffts%nl(1:ngm)) *&
             & fac_in(1:ngm) 
        !
        !  and transformed back to real space
        !
        CALL invfft ('Rho', vhart (:, is), dffts)
        !
     ENDDO
     !
     !
     ! Finally return the interaction
     !
     psi_int(:,ibnd2) = psi_int(:,ibnd2)  + vhart(:,1) *&
          &red_revc0(:,ibnd,ikq)
     !
  ENDDO
  !
END FUNCTION k2d_term_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! The following routines are wrapper functions for inv/fw fft handling both
!! the custom FFT grids and gamma tricks. They are the analogues of those found
!! in PW/src/realus.f90. Ideally both these and their counterparts should be
!! moved somewhere else but for now they live here.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invfft_orbital_custom_gamma(orbital, ibnd, nbnd, npwt, dfftt)

  USE kinds,        ONLY : DP
  USE fft_types,    ONLY : fft_type_descriptor

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)    :: orbital(:,:)
  INTEGER, INTENT(IN)        :: ibnd, nbnd, npwt
  TYPE(fft_type_descriptor), INTENT(IN)  :: dfftt
  !
  psic=(0.0_dp, 0.0_dp)
  !
  IF (ibnd < nbnd) THEN
     !
     psic(dfftt%nl(1:npwt)) = orbital(1:npwt,ibnd) + &
          &(0.0_dp, 1.0_dp) * orbital(1:npwt,ibnd+1) 
     psic(dfftt%nlm(1:npwt)) = CONJG(orbital(1:npwt,ibnd) - &   
          &(0.0_dp, 1.0_dp) * orbital(1:npwt,ibnd+1))
     !
  ELSE
     !
     psic(dfftt%nl(1:npwt))  = orbital(1:npwt,ibnd)
     psic(dfftt%nlm(1:npwt)) =CONJG(orbital(1:npwt,ibnd))
     !
  ENDIF
  !
  CALL invfft ('Wave', psic, dfftt)
  !
  RETURN
  !
END SUBROUTINE invfft_orbital_custom_gamma

SUBROUTINE fwfft_orbital_custom_gamma(orbital, ibnd, nbnd, npwt, dfftt)

  USE kinds,        ONLY : DP
  USE fft_types,    ONLY : fft_type_descriptor

  IMPLICIT NONE

  COMPLEX(DP), INTENT(INOUT)    :: orbital(:,:)
  INTEGER, INTENT(IN)        :: ibnd, nbnd, npwt
  TYPE(fft_type_descriptor), INTENT(IN)  :: dfftt

  ! Workspaces
  COMPLEX(DP) :: fp, fm
  ! Counters
  INTEGER :: j
  !
  CALL fwfft ('Wave', psic(:), dfftt)
  !
  IF (ibnd < nbnd) THEN
     !
     ! two ffts at the same time
     DO j = 1, npwt
        fp = (psic(dfftt%nl(j)) + psic(dfftt%nlm(j)))*0.5d0
        fm = (psic(dfftt%nl(j)) - psic(dfftt%nlm(j)))*0.5d0
        orbital( j, ibnd)   = CMPLX( DBLE(fp), AIMAG(fm),kind=DP)
        orbital( j, ibnd+1) = CMPLX(AIMAG(fp),- DBLE(fm),kind=DP)
     ENDDO
     !
  ELSE
     !
     orbital(1:npwt,ibnd)=psic(dfftt%nl(1:npwt))
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE fwfft_orbital_custom_gamma

SUBROUTINE invfft_orbital_ibnd2_gamma(orbital, psitemp, ibnd2, npwt, dfftt)

  USE kinds,        ONLY : DP
  USE fft_types,    ONLY : fft_type_descriptor

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)    :: orbital(:,:)
  INTEGER, INTENT(IN)        :: ibnd2, npwt
  TYPE(fft_type_descriptor), INTENT(IN)  :: dfftt
  COMPLEX(DP), INTENT(OUT)   :: psitemp(:) 
  !
  psitemp=(0.0_dp, 0.0_dp)
  !
  psitemp(dfftt%nl(1:npwt))  = orbital(1:npwt,ibnd2)
  psitemp(dfftt%nlm(1:npwt)) = CONJG(orbital(1:npwt,ibnd2))
  !
  CALL invfft ('Wave', psitemp, dfftt)
  !
  RETURN
  !
END SUBROUTINE invfft_orbital_ibnd2_gamma

END MODULE lr_exx_kernel
