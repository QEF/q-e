
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  ! --- Wrapper routine: performs bgrp parallelization on non-distributed bands
  ! --- if suitable and required, calls old H\psi routine h_psi_
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE funct,            ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, set_bgrp_indices, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)   
  !
  INTEGER     :: m_start, m_end
  !
  CALL start_clock( 'h_psi_bgrp' )

  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. exact exchange is not active (if it is, band parallelization is already
  !    used in exx routines called by Hpsi)
  ! 3. there is more than one band, otherwise there is nothing to parallelize
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     ! use band parallelization here
     hpsi(:,:) = (0.d0,0.d0)
     CALL set_bgrp_indices(m,m_start,m_end)
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL h_psi_( lda, n, m_end-m_start+1, psi(1,m_start), hpsi(1,m_start) )
     CALL mp_sum(hpsi,inter_bgrp_comm)
  ELSE
     ! don't use band parallelization here
     CALL h_psi_( lda, n, m, psi, hpsi )
  END IF

  CALL stop_clock( 'h_psi_bgrp' )
  RETURN
  !
END SUBROUTINE h_psi
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi_( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,    ONLY : DP
  USE wavefunctions_module, ONLY : psic
  USE bp,       ONLY : lelfield,l3dstring,gdir, efield, efield_cry
  USE becmod,   ONLY : bec_type, becp, calbec
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs  
  USE wvfct,    ONLY : g2kin
  USE uspp,     ONLY : vkb, nkb
  USE ldaU,     ONLY : lda_plus_u, U_projection
  USE gvect,    ONLY : gstart
  USE funct,    ONLY : dft_is_meta
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY: npol, noncolin
  USE realus,   ONLY : real_space, &
                       invfft_orbital_gamma, fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, & 
                       invfft_orbital_k, fwfft_orbital_k, calbec_rs_k, add_vuspsir_k, & 
                       v_loc_psir_inplace
  USE fft_base, ONLY : dffts, dtgs
  USE exx,      ONLY : vexx, vexxace_gamma, vexxace_k
  USE funct,    ONLY : exx_is_active
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)   
  !
  INTEGER     :: ipol, ibnd, incr
  REAL(dp)    :: ee
  !
  CALL start_clock( 'h_psi' )

  hpsi (:, 1:m) = (0.0_dp, 0.0_dp)

  CALL start_clock( 'h_psi:pot' )
  !
  ! ... the local potential V_Loc psi
  !
  !
  IF ( gamma_only ) THEN
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dtgs%have_task_groups ) then 
           incr = 2 * dtgs%nogrp
        ELSE
           incr = 2
        ENDIF
        DO ibnd = 1, m, incr
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_gamma(psi,ibnd,m) 
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock( 'h_psi:calbec' )
           CALL calbec_rs_gamma(ibnd,m,becp%r) 
     CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace(ibnd,m) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_gamma(ibnd,m)
           ! ... transform psic back in reciprocal space and assign it to hpsi
           CALL fwfft_orbital_gamma(hpsi,ibnd,m) 
        END DO
        !
     ELSE
        ! ... usual reciprocal-space algorithm
        CALL vloc_psi_gamma ( lda, n, m, psi, vrs(1,current_spin), hpsi ) 
        !
     ENDIF 
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL vloc_psi_nc ( lda, n, m, psi, vrs, hpsi )
     !
  ELSE  
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dtgs%have_task_groups ) then 
           incr = dtgs%nogrp
        ELSE
           incr = 1
        ENDIF
        DO ibnd = 1, m
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_k(psi,ibnd,m) 
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock( 'h_psi:calbec' )
           CALL calbec_rs_k(ibnd,m) 
     CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace(ibnd,m) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_k(ibnd,m)
           ! ... transform psic back in reciprocal space and assign it to hpsi
           CALL fwfft_orbital_k(hpsi,ibnd,m) 
        END DO
        !
     ELSE
        !
        CALL vloc_psi_k ( lda, n, m, psi, vrs(1,current_spin), hpsi )
        !
     END IF  
     !
  END IF  

  IF ( nkb > 0 .AND. .NOT. real_space) THEN
     !
     CALL start_clock( 'h_psi:calbec' )
     CALL calbec ( n, vkb, psi, becp, m )
     CALL stop_clock( 'h_psi:calbec' )
     CALL add_vuspsi( lda, n, m, hpsi )
     !
  END IF
  !  
  CALL stop_clock( 'h_psi:pot' )
  !
  ! ... Here we add the kinetic energy (k+G)^2 psi
  !
  DO ibnd = 1, m
     hpsi (1:n, ibnd) = hpsi(1:n, ibnd) + g2kin (1:n) * psi (1:n, ibnd)
     IF ( noncolin ) THEN
        hpsi (lda+1:lda+n, ibnd) = hpsi(lda+1:lda+n,ibnd) + g2kin (1:n) * psi (lda+1:lda+n, ibnd)
     END IF
  END DO
  !
  if (dft_is_meta()) call h_psi_meta (lda, n, m, psi, hpsi)
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u .AND. U_projection.NE."pseudo" ) THEN
     !
     IF (noncolin) THEN
        CALL vhpsi_nc( lda, n, m, psi, hpsi )
     ELSE
        call vhpsi( lda, n, m, psi, hpsi )
     ENDIF
     !
  ENDIF
  !
  !
  ! ... Here the product with the non local potential V_NL psi
  ! ... (not in the real-space case: it is done together with V_loc)
  !
  IF ( exx_is_active() ) THEN
#if defined(__EXX_ACE) 
     IF (gamma_only) THEN
        CALL vexxace_gamma(lda,m,psi,ee,hpsi)
     ELSE
        CALL vexxace_k(lda,m,psi,ee,hpsi) 
     END IF
#else  
     CALL vexx( lda, n, m, psi, hpsi, becp )
#endif
  END IF
  !
  ! ... electric enthalpy if required
  !
  IF ( lelfield ) THEN
     !
     IF ( .NOT.l3dstring ) THEN
        CALL h_epsi_her_apply( lda, n, m, psi, hpsi,gdir, efield )
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_apply( lda, n, m, psi, hpsi,ipol,efield_cry(ipol) )
        END DO
     END IF
     !
  END IF
  !
  ! ... Gamma-only trick: set to zero the imaginary part of hpsi at G=0
  !
  IF ( gamma_only .AND. gstart == 2 ) &
      hpsi(1,1:m) = CMPLX( DBLE( hpsi(1,1:m) ), 0.D0 ,kind=DP)
  !
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
END SUBROUTINE h_psi_
