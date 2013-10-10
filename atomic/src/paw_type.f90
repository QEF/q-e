!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE paw_type
   USE kinds,       ONLY : DP
   USE radial_grids, ONLY : radial_grid_type, nullify_radial_grid, deallocate_radial_grid
   !
   PUBLIC
   !
   TYPE :: paw_t
      !
      ! Type describing a PAW dataset (temporary).
      ! Functions are defined on a logarithmic radial mesh.
      !
      CHARACTER(LEN=2) :: symbol
      REAL (DP)        :: zval
      REAL (DP)        :: z
      CHARACTER(LEN=80):: dft
      TYPE(radial_grid_type) :: grid
      REAL (DP)        :: rmatch_augfun  ! the matching radius for augmentation charges
      LOGICAL          :: nlcc ! nonlinear core correction
      INTEGER          :: nwfc ! number of wavefunctions/projectors
      INTEGER          :: lmax ! maximum angular momentum of projectors
      INTEGER          :: rel  ! the relativistic level
      INTEGER, POINTER :: l(:) !l(nwfsx) ! angular momentum of projectors
      INTEGER, POINTER :: ikk(:) !ikk(nwfsx) ! cutoff radius for the projectors
      INTEGER          :: irc ! r(irc) = radius of the augmentation sphere
      CHARACTER(LEN=2),POINTER ::  els (:) ! the name of the wavefunction
      REAL (DP), POINTER :: &
         oc(:), &          !(nwfsx) the occupations
         enl(:), &         !(nwfsx) the energy of the wavefunctions
         jj (:), &         ! the total angular momentum
         rcutus (:), &     ! the cutoff
         aewfc(:,:), &     !(ndmx,nwfsx) all-electron wavefunctions
         aewfc_rel(:,:), &     !(ndmx,nwfsx) all-electron wavefunctions
         pswfc(:,:), &     !(ndmx,nwfsx) pseudo wavefunctions
         proj(:,:), &      !(ndmx,nwfsx) projectors
         augfun(:,:,:,:), &!(ndmx,nwfsx,nwfsx,0:2*lmaxx+1),
         augmom(:,:,:), &  !(nwfsx,nwfsx,0:2*lmaxx) moments of the augmentation functions
         aeccharge(:), &   !(ndmx) AE core charge * 4PI r^2
         psccharge(:), &   !(ndmx) PS core charge * 4PI r^2
         pscharge(:), &    !(ndmx) PS charge * 4PI r^2
         aeloc(:), &       !(ndmx) descreened AE potential: v_AE-v_H[n1]-v_XC[n1+nc]
         psloc(:), &       !(ndmx) descreened local PS potential: v_PS-v_H[n~+n^]-v_XC[n~+n^+n~c]
         kdiff(:,:), &     !(nwfsx,nwfsx) kinetic energy differences
         dion(:,:)         !(nwfsx,nwfsx) descreened D coeffs
   !!!  Notes about screening:
   !!!       Without nlcc, the local PSpotential is descreened with n~+n^ only.
   !!!       The local AEpotential is descreened ALWAYS with n1+nc. This improves
   !!!       the accuracy, and will not cost in the plane wave code (atomic
   !!!       contribution only).
   END TYPE paw_t

 CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Nullify, allocate and deallocate for paw_t type. Used only
! in atomic code.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE nullify_pseudo_paw( paw )
        TYPE( paw_t ), INTENT(INOUT) :: paw
        CALL nullify_radial_grid( paw%grid ) ! nullify grid object, here grid is not a POINTER!
        NULLIFY( paw%l, paw%ikk )
        NULLIFY( paw%oc, paw%enl, paw%aewfc, paw%aewfc_rel, paw%pswfc, paw%proj)
        NULLIFY( paw%augfun, paw%augmom, paw%aeccharge, paw%psccharge, paw%pscharge )
        NULLIFY( paw%aeloc, paw%psloc, paw%dion )
        NULLIFY( paw%kdiff )
        RETURN
        END SUBROUTINE nullify_pseudo_paw
        
        SUBROUTINE allocate_pseudo_paw( paw, size_mesh, size_nwfc, size_lmax )
        TYPE( paw_t ), INTENT(INOUT) :: paw
        INTEGER, INTENT(IN) :: size_mesh, size_nwfc, size_lmax
        !WRITE(0,"(a,3i5)") "Allocating PAW setup: ",size_mesh, size_nwfc, size_lmax
        ALLOCATE ( paw%l(size_nwfc) )
        ALLOCATE ( paw%jj(size_nwfc) )
        ALLOCATE ( paw%ikk(size_nwfc) )
        ALLOCATE ( paw%oc(size_nwfc) )
        ALLOCATE ( paw%rcutus(size_nwfc) )
        ALLOCATE ( paw%els(size_nwfc) )
        ALLOCATE ( paw%enl(size_nwfc) )
        ALLOCATE ( paw%aewfc(size_mesh,size_nwfc) )
        ALLOCATE ( paw%aewfc_rel(size_mesh,size_nwfc) )
        ALLOCATE ( paw%pswfc(size_mesh,size_nwfc) )
        ALLOCATE ( paw%proj (size_mesh,size_nwfc) )
        ALLOCATE ( paw%augfun(size_mesh,size_nwfc,size_nwfc,0:2*size_lmax) )
        ALLOCATE ( paw%augmom(size_nwfc,size_nwfc,0:2*size_lmax) )
        ALLOCATE ( paw%aeccharge(size_mesh) )
        ALLOCATE ( paw%psccharge(size_mesh) )
        ALLOCATE ( paw%pscharge(size_mesh) )
        ALLOCATE ( paw%aeloc(size_mesh) )
        ALLOCATE ( paw%psloc(size_mesh) )
        ALLOCATE ( paw%kdiff(size_nwfc,size_nwfc) )
        ALLOCATE ( paw%dion (size_nwfc,size_nwfc) )
        END SUBROUTINE allocate_pseudo_paw
        
        SUBROUTINE deallocate_pseudo_paw( paw )
        TYPE( paw_t ), INTENT(INOUT) :: paw
        CALL deallocate_radial_grid( paw%grid ) ! here grid is not a POINTER!
        IF( ASSOCIATED( paw%l ) ) DEALLOCATE( paw%l )
        IF( ASSOCIATED( paw%jj ) ) DEALLOCATE( paw%jj )
        IF( ASSOCIATED( paw%ikk ) ) DEALLOCATE( paw%ikk )
        IF( ASSOCIATED( paw%oc ) ) DEALLOCATE( paw%oc )
        IF( ASSOCIATED( paw%els ) ) DEALLOCATE( paw%els )
        IF( ASSOCIATED( paw%rcutus ) ) DEALLOCATE( paw%rcutus )
        IF( ASSOCIATED( paw%enl ) ) DEALLOCATE( paw%enl )
        IF( ASSOCIATED( paw%aewfc_rel ) ) DEALLOCATE( paw%aewfc_rel ) 
        IF( ASSOCIATED( paw%aewfc ) ) DEALLOCATE( paw%aewfc )
        IF( ASSOCIATED( paw%pswfc ) ) DEALLOCATE( paw%pswfc )
        IF( ASSOCIATED( paw%proj ) ) DEALLOCATE( paw%proj )
        IF( ASSOCIATED( paw%augfun ) ) DEALLOCATE( paw%augfun )
        IF( ASSOCIATED( paw%augmom ) ) DEALLOCATE( paw%augmom )
        IF( ASSOCIATED( paw%aeccharge ) ) DEALLOCATE( paw%aeccharge )
        IF( ASSOCIATED( paw%psccharge ) ) DEALLOCATE( paw%psccharge )
        IF( ASSOCIATED( paw%pscharge ) ) DEALLOCATE( paw%pscharge )
        IF( ASSOCIATED( paw%aeloc ) ) DEALLOCATE( paw%aeloc )
        IF( ASSOCIATED( paw%psloc ) ) DEALLOCATE( paw%psloc )
        IF( ASSOCIATED( paw%kdiff ) ) DEALLOCATE( paw%kdiff )
        IF( ASSOCIATED( paw%dion ) ) DEALLOCATE( paw%dion )
        RETURN
        END SUBROUTINE deallocate_pseudo_paw

END MODULE paw_type
