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
  USE basis,      ONLY : natomwfc
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  !
  SAVE
  !
  INTEGER, PARAMETER :: nspinx=2
  COMPLEX(DP), ALLOCATABLE :: &
       wfcU(:,:),             &! atomic wfcs with U term
       d_spin_ldau(:,:,:)      ! the rotations in spin space for all symmetries
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
       nwfcU,                &! total no. of atomic wavefunctions having U term
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
       oatwfc(:), offsetU(:)  ! offset of atomic wfcs used for projections
  REAL(DP), ALLOCATABLE :: &
       q_ae(:,:,:),          &! coefficients for projecting onto beta functions
       q_ps(:,:,:)            ! (matrix elements on AE and PS atomic wfcs)
  !
CONTAINS
  !
  SUBROUTINE init_lda_plus_u ( psd, noncolin )
    !
    IMPLICIT NONE
    CHARACTER (LEN=2), INTENT(IN) :: psd(:)
    LOGICAL, INTENT(IN) :: noncolin
    !
    INTEGER, EXTERNAL :: set_Hubbard_l
    INTEGER :: na, nt
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
       IF ( U_projection == 'pseudo' ) CALL errore( 'init_lda_plus_u', &
            & 'full LDA+U not implemented with pseudo projection type', 1 )
       !
       IF (noncolin) THEN
          ALLOCATE( d_spin_ldau(2,2,48) )
          call comp_dspinldau ()
       END IF
       
       DO nt = 1, ntyp
          IF (Hubbard_alpha(nt)/=0.d0 ) CALL errore( 'init_lda_plus_u', &
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
       CALL errore( 'init_lda_plus_u', 'lda_plus_u_kind should be 0 or 1', 1 )
    endif
    IF ( Hubbard_lmax == -1 ) CALL errore( 'init_lda_plus_u', &
         'lda_plus_u calculation but Hubbard_l not set', 1 )
    IF ( Hubbard_lmax > 3 ) &
         CALL errore( 'init_lda_plus_u', 'Hubbard_l should not be > 3 ', 1 )

    ! compute index of atomic wfcs used as projectors
    IF ( .NOT.allocated(oatwfc)) ALLOCATE ( oatwfc(nat) )
    CALL offset_atom_wfc ( .false., oatwfc, nwfcU )
    ! nwfcU is set to natomwfc by the routine above
    IF ( nwfcU .NE.natomwfc ) &
         CALL errore ('offset_atom_wfc', 'wrong number of wavefunctions', 1)
    ! for each atom, compute index of its projectors (among projectors only)
    IF ( .NOT.allocated(offsetU)) ALLOCATE ( offsetU(nat) )
    CALL offset_atom_wfc ( .true., offsetU, nwfcU )
    !
  END SUBROUTINE init_lda_plus_u
  !
  SUBROUTINE deallocate_ldaU ( flag )
  !
  LOGICAL, INTENT (in) :: flag
  !
  IF ( flag ) THEN
     IF ( ALLOCATED( oatwfc ) )     DEALLOCATE( oatwfc )
     IF ( ALLOCATED( offsetU ) )    DEALLOCATE( offsetU )
     IF ( ALLOCATED( q_ae ) )       DEALLOCATE( q_ae )
     IF ( ALLOCATED( q_ps ) )       DEALLOCATE( q_ps )
  END IF
  IF ( ALLOCATED( wfcU ) )       DEALLOCATE( wfcU )
  !
  END SUBROUTINE deallocate_ldaU
  !
  SUBROUTINE copy_U_wfc ( swfcatom, noncolin )
  !
  !  Copy (orthogonalized) atomic wavefunctions "swfcatom"
  !  having a Hubbard U correction to array "wfcU"
  !
  IMPLICIT NONE
  COMPLEX (KIND=dp), INTENT (IN) :: swfcatom(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: noncolin
  LOGICAL :: twice
  INTEGER :: na, nt, m1, m2

  IF ( PRESENT (noncolin) ) THEN
     twice = noncolin
  ELSE
     twice = .FALSE.
  END IF
  DO na=1,nat
     nt = ityp(na)
     if ( is_hubbard(nt) ) then
        m1 = 1
        m2 = 2*hubbard_l(nt)+1
        IF ( twice ) m2 = 2*m2
        wfcU(:,offsetU(na)+m1:offsetU(na)+m2) = swfcatom(:,oatwfc(na)+m1:oatwfc(na)+m2)
     end if
  END DO

  END SUBROUTINE copy_U_wfc

END MODULE ldaU
!
