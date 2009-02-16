
! Copyright (C) 2002-2007 Quantum-ESPRESSO group
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
  USE kinds,    ONLY : DP
  USE bp,       ONLY : lelfield,l3dstring,gdir, efield, efield_cry
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY: npol, noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)   
  !
  INTEGER     :: ipol
  !
  CALL start_clock( 'h_psi' )
  !  
  IF ( gamma_only ) THEN
     !
     CALL h_psi_gamma( lda, n, m, psi, hpsi )
     !
  ELSE IF ( noncolin ) THEN
     !
     CALL h_psi_nc( lda, n, m, psi, hpsi )
     !
  ELSE  
     !
     CALL h_psi_k ( lda, n, m, psi, hpsi )
     !
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
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
END SUBROUTINE h_psi
!
!-----------------------------------------------------------------------
SUBROUTINE h_psi_gamma ( lda, n, m, psi, hpsi )
  !-----------------------------------------------------------------------
  ! 
  ! ... gamma version
  !
  USE kinds,    ONLY : DP
  USE becmod,   ONLY : rbecp, calbec
  USE uspp,     ONLY : vkb, nkb
  USE wvfct,    ONLY : g2kin
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs  
  USE gvect,    ONLY : gstart
#ifdef EXX
  USE exx,      ONLY : vexx
  USE funct,    ONLY : exx_is_active
#endif
  USE funct,    ONLY : dft_is_meta
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda,m)   
  !
  INTEGER :: ibnd, j
  !
  !
  CALL start_clock( 'h_psi:init' )
  !
  ! ... Here we apply the kinetic energy (k+G)^2 psi
  !
  DO ibnd = 1, m
     !
     hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
     !
  END DO
  !
  CALL stop_clock( 'h_psi:init' )
  !      
  if (dft_is_meta()) call h_psi_meta (lda, n, m, psi, hpsi)
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u ) CALL vhpsi( lda, n, m, psi, hpsi )
  !
  ! ... the local potential V_Loc psi
  !
  CALL vloc_psi( lda, n, m, psi, vrs(1,current_spin), hpsi )
  !
  ! ... Here the product with the non local potential V_NL psi
  !
  IF ( nkb > 0 ) THEN
     !
     CALL start_clock( 'h_psi:vnl' )
     CALL calbec ( n, vkb, psi, rbecp, m )
     CALL add_vuspsi( lda, n, m, psi, hpsi )
     CALL stop_clock( 'h_psi:vnl' )
     !
  END IF
  !
  ! ... set to zero the imaginary part of hpsi at G=0
  !
  IF ( gstart == 2 ) hpsi(1,1:m) = CMPLX( DBLE( hpsi(1,1:m) ), 0.D0 )
  !
#ifdef EXX
  IF ( exx_is_active() ) CALL vexx( lda, n, m, psi, hpsi )
#endif
  !
  RETURN
  !
END SUBROUTINE h_psi_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE h_psi_k ( lda, n, m, psi, hpsi )
  !-----------------------------------------------------------------------
  !
  ! ... k-points version
  !
  USE kinds,    ONLY : DP
  USE becmod,   ONLY : becp, calbec
  USE uspp,     ONLY : vkb, nkb
  USE wvfct,    ONLY : g2kin
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs  
#ifdef EXX
  USE exx,      ONLY : vexx
  USE funct,    ONLY : exx_is_active
#endif
  USE funct,    ONLY : dft_is_meta
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi(lda,m) 
  COMPLEX(DP), INTENT(OUT):: hpsi(lda,m)   
  !
  INTEGER :: ibnd, j
  !
  CALL start_clock( 'h_psi:init' )
  !
  ! ... Here we apply the kinetic energy (k+G)^2 psi
  !
  DO ibnd = 1, m
     !
     hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
     !
  END DO
  !
  CALL stop_clock( 'h_psi:init' )
  !      
  if (dft_is_meta()) call h_psi_meta (lda, n, m, psi, hpsi)
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u ) CALL vhpsi( lda, n, m, psi, hpsi )
  !
  ! ... the local potential V_Loc psi
  !
  CALL vloc_psi_k( lda, n, m, psi, vrs(1,current_spin), hpsi )
  !
  ! ... Here the product with the non local potential V_NL psi
  !
  IF ( nkb > 0 ) THEN
     !
     CALL start_clock( 'h_psi:vnl' )
     CALL calbec( n, vkb, psi, becp, m )
     CALL add_vuspsi( lda, n, m, psi, hpsi )
     CALL stop_clock( 'h_psi:vnl' )
     !
  END IF
  !
#ifdef EXX
  IF ( exx_is_active() ) CALL vexx( lda, n, m, psi, hpsi )
#endif
  !
  RETURN
  !
END SUBROUTINE h_psi_k     
!
!-----------------------------------------------------------------------
subroutine h_psi_nc(lda, n, m, psi, hpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     matrix with m wavefunctions contained in psi
  ! input:
  !     lda   leading dimension of arrays psi, hpsi
  !     n     true dimension of psi, hpsi
  !     m     number of states psi
  !     psi
  ! output:
  !     hpsi  H*psi
  !
  USE kinds, ONLY : DP
  use uspp, only: vkb, nkb
  use wvfct, only: g2kin
  use ldaU, only : lda_plus_u
  use scf, only: vrs
  use becmod, only: becp_nc, calbec
  use noncollin_module,     only: npol
  !
  implicit none
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT):: hpsi(lda,npol,m)   
  !
  complex(DP) :: sup, sdwn
  integer :: ibnd,j,ipol
  !
  call start_clock ('h_psi:init')
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  do ibnd = 1, m
     do ipol = 1, npol
        do j = 1, n
           hpsi (j, ipol, ibnd) = g2kin (j) * psi (j+(ipol-1)*lda, ibnd)
        enddo
     enddo
  enddo

  call stop_clock ('h_psi:init')
  !
  ! Here we add the Hubbard potential times psi
  !
  !!!if (lda_plus_u) call vhpsi_nc (lda, n, m, psi(1,1), hpsi(1,1,1))
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  CALL vloc_psi_nc ( lda, n, m, psi, vrs, hpsi )
  !
  !  Here the product with the non local potential V_NL psi
  !
  if (nkb > 0) then
     !
     CALL start_clock( 'h_psi:vnl' )
     CALL calbec ( n, vkb, psi, becp_nc, m )
     CALL add_vuspsi_nc (lda, n, m, psi, hpsi(1,1,1))
     CALL stop_clock( 'h_psi:vnl' )
     !
  end if
  return
end subroutine h_psi_nc
