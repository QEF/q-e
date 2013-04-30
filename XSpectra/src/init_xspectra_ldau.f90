!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  this routine initializes some variables for DFT+U
!  calculation, in particular atomic wfc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


subroutine init_xanes_ldau
  USE ldaU,             ONLY : lda_plus_u, Hubbard_U, Hubbard_l, &
                               Hubbard_alpha, Hubbard_lmax, U_projection 
  USE basis,            only : natomwfc  
  USE xspectra,         only : xread_wf, U_projection_type
  use wvfct,            ONLY : npwx,nbndx,nbnd,npw
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE lsda_mod,         ONLY : lsda, nspin
  USE uspp_param,       ONLY : upf
  implicit none

  integer :: ldim, nt
  INTEGER :: set_Hubbard_l

  call init_at_1
  Hubbard_lmax = -1
  !
  DO nt = 1, ntyp
     !
     IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
        !
        Hubbard_l(nt) = set_Hubbard_l( upf(nt)%psd )
        Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
        WRITE( UNIT = stdout, FMT = * ) &
           ' HUBBARD L FOR TYPE ',upf(nt)%psd,' IS ', Hubbard_l(nt)
        !
     END IF
     !
  END DO
  !
  WRITE( UNIT = stdout, FMT = * ) &
       ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
  !
  IF ( Hubbard_lmax == -1 ) CALL errore( 'setup', &
                   & 'lda_plus_u calculation but Hubbard_l not set', 1 )

  ldim = 2 * Hubbard_lmax + 1
  U_projection=U_projection_type

END SUBROUTINE init_xanes_ldau

! Here is a slightly modified version of orthoatwfc, for a 
! selected k-point and without writing the result
! modified by CG

!-----------------------------------------------------------------------
SUBROUTINE init_xanes_ldau_2(ik)
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize all the atomic wfcs. This is
  ! useful when we want to compute the occupation of the atomic orbitals
  ! in order to make lda+U calculations
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, nwordwfcU, diropn
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk
  USE ldaU,       ONLY : U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx, npw, igk
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol

  ! 
  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol
  ! the k point under consideration
  ! counter on bands
  LOGICAL :: orthogonalize_wfc, normalize_only
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)
  LOGICAL :: exst

  ALLOCATE (wfcatom( npwx*npol, natomwfc))    

  IF (U_projection=="file") THEN

     WRITE( stdout,*) 'LDA+U Projector read from file '
     CALL diropn( iunhub, 'hub', 2*nwordwfcU, exst )
     CALL davcio( wfcU, nwordwfcU, iunhub, ik, -1 )
     CLOSE( UNIT = iunhub, STATUS = 'KEEP' )

     RETURN
  END IF

  IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are NOT orthogonalized'
  ELSE IF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are orthogonalized'
     IF (gamma_only) CALL errore('init_xanes_ldau2', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are normalized but NOT orthogonalized'
     IF (gamma_only) CALL errore('init_xanes_ldau2', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE( stdout,*) "U_projection_type =", U_projection
     CALL errore ("init_xanes_ldau2"," this U_projection_type is not valid",1)
  END IF

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  CALL atomic_wfc (ik, wfcatom)
  CALL calbec (npw, vkb, wfcatom, becp)
  CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

  IF (orthogonalize_wfc) &
     CALL ortho_swfc ( normalize_only, natomwfc, wfcatom, swfcatom )
  !
  CALL copy_U_wfc (swfcatom)
  !
  DEALLOCATE (wfcatom)
  CALL deallocate_bec_type (becp )
  !
  RETURN
     
END SUBROUTINE init_xanes_ldau_2


