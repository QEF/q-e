!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
!
MODULE ldaU
  !
  ! ... The quantities needed in lda+U calculations
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lqmax, ntypx
  !
  SAVE
  !
  INTEGER, PARAMETER :: nspinx=2
  COMPLEX(DP), ALLOCATABLE :: &
       swfcatom(:,:),         &! orthogonalized atomic wfcs
       d_spin_ldau(:,:,:)      ! the rotations in spin space for all the symmetries
  REAL(DP) :: &
       eth,                  &! the Hubbard contribution to the energy
       Hubbard_U(ntypx),     &! the Hubbard U
       Hubbard_J0(ntypx),    &! the Hubbard J, in simplified LDA+U
       Hubbard_J(3,ntypx),   &! extra Hubbard parameters for full LDA+U:  
                              !     p: J(1) = J
                              !     d: J(1) = J, J(2) =  B 
                              !     f: J(1) = J, J(2) = E2, J(3) = E3 
       Hubbard_alpha(ntypx), &! the Hubbard alpha (used to calculate U)
       Hubbard_beta(ntypx),  &! the Hubbard beta (used to calculate J0)
       starting_ns(lqmax,nspinx,ntypx) !
  INTEGER :: &
       niter_with_fixed_ns,  &! no. of iterations with fixed ns
       lda_plus_u_kind,      &! 1/0 --> full/simplified(old) LDA+U calculation
       Hubbard_l(ntypx),     &! the angular momentum of Hubbard states
       Hubbard_lmax = 0       ! maximum angular momentum of Hubbard states
  LOGICAL :: &
       is_hubbard(ntypx),    &! .TRUE. if this atom species has U correction
       lda_plus_u,           &! .TRUE. if lda+u calculation is performed
       conv_ns                ! .TRUE. if ns are converged
  CHARACTER(len=30) :: &      ! 'atomic', 'ortho-atomic', 'file'
       U_projection           ! specifies how input coordinates are given
  INTEGER, ALLOCATABLE :: &
       oatwfc(:)              ! offset of atomic wfcs used for projections
  REAL(DP), ALLOCATABLE :: &
       q_ae(:,:,:),          &! coefficients for projecting onto beta functions
       q_ps(:,:,:)            ! (matrix elements on AE and PS atomic wfcs)
  !
CONTAINS
  !
  SUBROUTINE init_lda_plus_u ( psd, noncolin )
    !
    USE ions_base, ONLY: nat, ntyp => nsp
    !
    IMPLICIT NONE
    CHARACTER (LEN=2), INTENT(IN) :: psd(:)
    LOGICAL, INTENT(IN) :: noncolin
    !
    INTEGER, EXTERNAL :: set_Hubbard_l
    INTEGER :: nt
    !
    !
    IF ( .NOT. lda_plus_u ) THEN
       Hubbard_lmax = 0
       is_hubbard(:) = .FALSE.
       RETURN
    END IF
    !
    Hubbard_lmax = -1
    ! Set the default of Hubbard_l for the species which have
    ! Hubbard_U=0 (in that case set_Hubbard_l will not be called)
    Hubbard_l(:) = -1
    !
    if ( lda_plus_u_kind == 0 ) then
       !
       DO nt = 1, ntyp
          !
          is_hubbard(nt) = Hubbard_U(nt)/= 0.0_dp .OR. &
                           Hubbard_alpha(nt) /= 0.0_dp .OR. &
                           Hubbard_J0(nt) /= 0.0_dp .OR. &
                           Hubbard_beta(nt)/= 0.0_dp
          !
          IF ( is_hubbard(nt) ) THEN
             Hubbard_l(nt) = set_Hubbard_l( psd(nt) )
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
          END IF
          !
       END DO
       !
    ELSE IF ( lda_plus_u_kind == 1 ) THEN
       !
       IF ( U_projection == 'pseudo' ) CALL errore( 'setup', &
            & 'full LDA+U not implemented with pseudo projection type', 1 )
       !
       IF (noncolin) THEN
          ALLOCATE( d_spin_ldau(2,2,48) )
          call comp_dspinldau ()
       END IF
       
       DO nt = 1, ntyp
          IF (Hubbard_alpha(nt)/=0.d0 ) CALL errore( 'setup', &
               'full LDA+U does not support Hubbard_alpha calculation', 1 )

          is_hubbard(nt) = Hubbard_U(nt)/= 0.0_dp .OR. &
                         ANY( Hubbard_J(:,nt) /= 0.0_dp )
        
          IF ( is_hubbard(nt) ) THEN
             !
             Hubbard_l(nt) = set_Hubbard_l( psd(nt) )
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
             !
             if (Hubbard_U(nt) == 0.0_dp) Hubbard_U(nt) = 1.d-14

             if ( Hubbard_l(nt) == 2 ) then
                if ( Hubbard_J(2,nt) == 0.d0 ) &
                     Hubbard_J(2,nt) = 0.114774114774d0 * Hubbard_J(1,nt)
             elseif ( Hubbard_l(nt) == 3 ) then
                if ( Hubbard_J(2,nt) == 0.d0 ) &
                     Hubbard_J(2,nt) = 0.002268d0 * Hubbard_J(1,nt)
                if ( Hubbard_J(3,nt)==0.d0 ) &
                     Hubbard_J(3,nt) = 0.0438d0 * Hubbard_J(1,nt)
             endif
          END IF
          !
       END DO
    else
       CALL errore( 'setup', 'lda_plus_u_kind should be 0 or 1', 1 )
    endif
    IF ( Hubbard_lmax == -1 ) CALL errore( 'setup', &
         'lda_plus_u calculation but Hubbard_l not set', 1 )
    IF ( Hubbard_lmax > 3 ) &
         CALL errore( 'setup', 'Hubbard_l should not be > 3 ', 1 )

    ! compute index of atomic wfcs used as projectors
    IF ( .NOT.allocated(oatwfc)) ALLOCATE ( oatwfc(nat) )
    CALL offset_atom_wfc ( nat, oatwfc )
    !
  END SUBROUTINE init_lda_plus_u

END MODULE ldaU
!
