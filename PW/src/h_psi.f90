
! Copyright (C) 2002-2009 Quantum ESPRESSO group
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
  ! --- bgrp parallelization allowed 
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE funct,            ONLY : exx_is_active
  USE mp_bands,         ONLY : tbgrp, set_bgrp_indices, inter_bgrp_comm
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

  if (tbgrp .and. .not. exx_is_active() ) then
      hpsi(:,:) = (0.d0,0.d0)
      call set_bgrp_indices(m,m_start,m_end)
      if (m_end >= m_start)  & !! at least one band in this band group
          call h_psi_( lda, n, m_end-m_start+1, psi(1,m_start), hpsi(1,m_start) )
      call mp_sum(hpsi,inter_bgrp_comm)
   else ! no one else to communicate with 
      call h_psi_( lda, n, m, psi, hpsi )
   end if

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
  USE realus,   ONLY : real_space, invfft_orbital_gamma, initialisation_level, &
                       fwfft_orbital_gamma, calbec_rs_gamma, &
                       add_vuspsir_gamma, v_loc_psir
  USE fft_base, ONLY : dffts
  USE exx,      ONLY : vexx
  USE funct,    ONLY : exx_is_active
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)   
  !
  INTEGER     :: ipol, ibnd, incr
  !
  CALL start_clock( 'h_psi' )
  !  
  ! ... Here we apply the kinetic energy (k+G)^2 psi
  !
  DO ibnd = 1, m
     hpsi (1:n, ibnd) = g2kin (1:n) * psi (1:n, ibnd)
     hpsi (n+1:lda,ibnd) = (0.0_dp, 0.0_dp)
     IF ( noncolin ) THEN
        hpsi (lda+1:lda+n, ibnd) = g2kin (1:n) * psi (lda+1:lda+n, ibnd)
        hpsi (lda+n+1:lda*npol, ibnd) = (0.0_dp, 0.0_dp)
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
  ! ... the local potential V_Loc psi
  !
  CALL start_clock( 'h_psi:vloc' )
  !
  IF ( gamma_only ) THEN
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%have_task_groups .AND. ( m >= dffts%nogrp )) then 
           incr = 2 * dffts%nogrp
        ELSE
           incr = 2
        ENDIF
        DO ibnd = 1, m, incr
           ! ... transform psi to real space, saved in temporary memory
           CALL invfft_orbital_gamma(psi,ibnd,m,.true.) 
           ! ... becp%r = < beta|psi> on psi in real space
           CALL calbec_rs_gamma(ibnd,m,becp%r) 
           ! ... psi is now replaced by hpsi ??? WHAT FOR ???
           CALL invfft_orbital_gamma(hpsi,ibnd,m)
           ! ... hpsi -> hpsi + psi*vrs  (psi read from temporary memory)
           CALL v_loc_psir(ibnd,m) 
           ! ... hpsi -> hpsi + vusp
           CALL  add_vuspsir_gamma(ibnd,m)
           ! ... transform back hpsi, clear psi in temporary memory
           CALL fwfft_orbital_gamma(hpsi,ibnd,m,.true.) 
        END DO
        !
     ELSE
        ! ... usual reciprocal-space algorithm
        CALL vloc_psi_gamma ( lda, n, m, psi, vrs(1,current_spin), hpsi ) 
     ENDIF 
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL vloc_psi_nc ( lda, n, m, psi, vrs, hpsi )
     !
  ELSE  
     !
     CALL vloc_psi_k ( lda, n, m, psi, vrs(1,current_spin), hpsi )
     !
  END IF  
  CALL stop_clock( 'h_psi:vloc' )
  !
  ! ... Here the product with the non local potential V_NL psi
  ! ... (not in the real-space case: it is done together with V_loc)
  !
  IF ( nkb > 0 .AND. .NOT. real_space) THEN
     !
     CALL start_clock( 'h_psi:vnl' )
     CALL calbec ( n, vkb, psi, becp, m )
     CALL add_vuspsi( lda, n, m, hpsi )
     CALL stop_clock( 'h_psi:vnl' )
     !
  END IF
  IF ( exx_is_active() ) CALL vexx( lda, n, m, psi, hpsi, becp )
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
