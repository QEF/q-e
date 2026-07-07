!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE scf
  !--------------------------------------------------------------------------
  !! This module contains variables and auxiliary routines needed for
  !! the self-consistent cycle.  
  !
  USE kinds,           ONLY : DP
  USE lsda_mod,        ONLY : nspin
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax, lda_plus_u_kind, ldmx, &
                              ldmx_b, ldmx_tot, max_num_neighbors, &
                              is_hubbard_back
  USE ions_base,       ONLY : nat
  USE xc_lib,          ONLY : xclib_dft_is
  USE fft_base,        ONLY : dfftp
  USE gvect,           ONLY : ngm
  USE gvecs,           ONLY : ngms
  USE ions_base,       ONLY : ntyp => nsp
  USE paw_variables,   ONLY : okpaw
  USE uspp_param,      ONLY : nhm
  USE control_flags,   ONLY : lxdm, sic
  !
  SAVE
  !
  ! Details of PAW implementation:
  !
  ! NOTE1: scf_type is used for two different quantities: density and potential.
  !        These correspond, for PAW, to becsum and D coefficients.
  !        Due to interference with the ultrasoft routines only the becsum part
  !        is stored in the structure (at the moment).
  !
  ! NOTE2: rho%bec is different from becsum for two reasons:
  !        1. rho%bec is mixed, while becsum is not
  !        2. for npool > 1 rho%bec is collected, becsum is not
  !           ( this is necessary to make the stress work)
  !
  !
  TYPE scf_type
     !! Used for two different quantities: density and potential
     REAL(DP),    ALLOCATABLE :: of_r(:,:)
     !! the charge density in R-space
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)
     !! the charge density in G-space
     REAL(DP),    ALLOCATABLE :: kin_r(:,:)
     !! the kinetic energy density in R-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:)
     !! the kinetic energy density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)
     !! the DFT+U occupation matrix
     REAL(DP),    ALLOCATABLE :: nsb(:,:,:,:)
     !! the DFT+U occupation matrix (background states)
     COMPLEX(DP), ALLOCATABLE :: ns_nc(:,:,:,:)
     !! the DFT+U occupation matrix - noncollinear case
     COMPLEX(DP), ALLOCATABLE :: nsg(:,:,:,:,:)
     !! the DFT+U+V generalized occupation matrix
     !! Matrix nsg(at1,m1,viz,m2,sp) stores the expectation value:
     !! <C^\dagger_{at1,m1,sp}C_{viz,m2,sp}>, where sp = spin and
     !! viz identifies the atom in the neighborhood of at1.
     REAL(DP),    ALLOCATABLE :: bec(:,:,:)
     !! the PAW hamiltonian elements
     REAL(DP),   ALLOCATABLE :: pol_r(:,:) 
     !! the polaron density in R-space
     COMPLEX(DP),ALLOCATABLE :: pol_g(:,:) 
     !! the polaron density in G-space
     REAL(DP) :: el_dipole
     !! electronic dipole, if a dipole field is present
  END TYPE scf_type
  !
  TYPE(scf_type) :: rho
  !! the charge density and its other components
  TYPE(scf_type) :: v
  !! the scf potential
  TYPE(scf_type) :: vnew
  !! used to correct the forces
  !
  REAL(DP) :: v_of_0
  !! vltot(G=0)      
  REAL(DP), ALLOCATABLE :: vltot(:)
  !! the local potential in real space
  REAL(DP), ALLOCATABLE :: vrs(:,:)
  !! the total pot. in real space (smooth grid)
  REAL(DP), ALLOCATABLE :: rho_core(:)
  !! the core charge in real space
  REAL(DP), ALLOCATABLE :: kedtau(:,:)
  !! position dependent kinetic energy enhancement factor
  COMPLEX(DP), ALLOCATABLE :: rhog_core(:)
  !! the core charge in reciprocal space
  REAL(DP), ALLOCATABLE :: tau_core(:)
  !! the kinetic energy density of the core charge in real space
  COMPLEX(DP), ALLOCATABLE :: taug_core(:)
  !! the kinetic energy density of the core charge in reciprocal space
  !
  !! DFT+U, colinear and noncolinear cases
  !! These variables are set every time create_scf_type is called
  !
  LOGICAL, PRIVATE :: lda_plus_u_co  ! collinear case
  LOGICAL, PRIVATE :: lda_plus_u_cob ! collinear case (background states)
  LOGICAL, PRIVATE :: lda_plus_u_nc  ! noncollinear case
  LOGICAL, PRIVATE :: lda_plus_u_v   ! U+V case
  !
CONTAINS
 !
 !----------------------------------------------------------
 SUBROUTINE create_scf_type( rho, do_not_allocate_becsum )
   !----------------------------------------------------------
   !! Creates an \(\text{scf_type}\) object by allocating all the 
   !! different terms.
   !
   IMPLICIT NONE
   !
   TYPE(scf_type) :: rho
   !! the object to create
   LOGICAL, INTENT(IN), OPTIONAL :: do_not_allocate_becsum ! PAW hack
   !! if true, the PAW part is ignored.
   !
   ! ... local variable
   !
   LOGICAL :: allocate_becsum ! PAW hack
   !
   ALLOCATE( rho%of_r(dfftp%nnr,nspin) )
   ALLOCATE( rho%of_g(ngm,nspin) )
   IF (xclib_dft_is('meta') .OR. lxdm) THEN
      ALLOCATE( rho%kin_r(dfftp%nnr,nspin) )
      ALLOCATE( rho%kin_g(ngm,nspin) )
   ELSE
      ALLOCATE( rho%kin_r(1,1) )
      ALLOCATE( rho%kin_g(1,1) )
   ENDIF
   !
   lda_plus_u_co  = lda_plus_u .AND. .NOT. ( nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_nc  = lda_plus_u .AND.       ( nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_cob = lda_plus_u_co .AND. ANY( is_hubbard_back(1:ntyp) )
   lda_plus_u_v   =  ( lda_plus_u_kind == 2 )
   !
   IF (lda_plus_u_co)  ALLOCATE( rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
   IF (lda_plus_u_cob) ALLOCATE( rho%nsb(ldmx_b,ldmx_b,nspin,nat) )
   IF (lda_plus_u_nc)  ALLOCATE( rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
   IF (lda_plus_u_v ) ALLOCATE ( rho%nsg(ldmx_tot,ldmx_tot,max_num_neighbors,nat,nspin) )
   !
   IF (okpaw) THEN ! See the top of the file for clarification
      IF ( PRESENT(do_not_allocate_becsum) ) THEN
         allocate_becsum = .NOT. do_not_allocate_becsum
      ELSE
         allocate_becsum = .TRUE.
      ENDIF
      IF (allocate_becsum) ALLOCATE( rho%bec(nhm*(nhm+1)/2,nat,nspin) )
   ENDIF
   !
   rho%el_dipole = 0._dp
   IF (sic) THEN
      IF(.NOT. ALLOCATED(rho%pol_r)) ALLOCATE(rho%pol_r(dfftp%nnr,nspin)) 
      IF(.NOT. ALLOCATED(rho%pol_g)) ALLOCATE(rho%pol_g(ngm,nspin)) 
   END IF
   !
   RETURN
   !
 END SUBROUTINE create_scf_type
 !
 !
 !----------------------------------------------------
 SUBROUTINE destroy_scf_type( rho )
   !---------------------------------------------------
   !! Deallocates an scf_type object
   !
   IMPLICIT NONE
   !
   TYPE(scf_type) :: rho
   !! the object to deallocate
   !
   IF (ALLOCATED(rho%of_r) )  DEALLOCATE( rho%of_r  )
   IF (ALLOCATED(rho%of_g) )  DEALLOCATE( rho%of_g  )
   IF (ALLOCATED(rho%kin_r))  DEALLOCATE( rho%kin_r )
   IF (ALLOCATED(rho%kin_g))  DEALLOCATE( rho%kin_g )
   IF (ALLOCATED(rho%ns)   )  DEALLOCATE( rho%ns    )
   IF (ALLOCATED(rho%nsb)  )  DEALLOCATE( rho%nsb   )
   IF (ALLOCATED(rho%ns_nc))  DEALLOCATE( rho%ns_nc )
   IF (ALLOCATED(rho%nsg  ))  DEALLOCATE( rho%nsg   )
   IF (ALLOCATED(rho%bec)  )  DEALLOCATE( rho%bec   )
   IF (ALLOCATED(rho%pol_r))  DEALLOCATE( rho%pol_r )
   IF (ALLOCATED(rho%pol_g))  DEALLOCATE( rho%pol_g )
   !
   RETURN
   !
 END SUBROUTINE destroy_scf_type
 !
 !----------------------------------------------------------------------------
 SUBROUTINE scf_type_COPY( X, Y )
  !----------------------------------------------------------------------------
  !! Works like DCOPY for \(\text{scf_type}\) copy variables: \(Y = X\).
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(IN)    :: X
  TYPE(scf_type), INTENT(INOUT) :: Y
  !
  Y%of_r = X%of_r
  Y%of_g = X%of_g
  !
  IF (xclib_dft_is('meta') .OR. lxdm) THEN
     Y%kin_r = X%kin_r
     Y%kin_g = X%kin_g
  ENDIF
  !
  CALL scf_ns_copy (X, Y)
  !
  IF (okpaw)          Y%bec   = X%bec
  Y%el_dipole = X%el_dipole
  IF (sic) THEN
     Y%pol_r = X%pol_r
     Y%pol_g = X%pol_g
  END IF
  !
  RETURN
  !
 END SUBROUTINE scf_type_COPY
 !
 !-----------------------------------------------------------------------
 SUBROUTINE scf_ns_copy ( rho1, rho2 )
  !-----------------------------------------------------------------------
  !! Copy Hubbard ns from rho1 into rho2
  !
  IMPLICIT NONE
  TYPE(scf_type), INTENT(IN)    :: rho1
  TYPE(scf_type), INTENT(INOUT) :: rho2
  !  
  IF (lda_plus_u_co)  rho2%ns(:,:,:,:)   = rho1%ns(:,:,:,:)
  IF (lda_plus_u_cob) rho2%nsb(:,:,:,:)  = rho1%nsb(:,:,:,:)
  IF (lda_plus_u_nc)  rho2%ns_nc(:,:,:,:)= rho1%ns_nc(:,:,:,:)
  IF (lda_plus_u_v)   rho2%nsg(:,:,:,:,:)= rho1%nsg(:,:,:,:,:)
  !
END SUBROUTINE scf_ns_copy
!-------------------------------------------------------------------------------
SUBROUTINE bcast_scf_type( rho, root, comm )
  !----------------------------------------------------------------------------
  !! Broadcast all mixed quantities from first pool to all others.
  !! Needed to prevent divergencies in k-point parallelization.
  !
  USE mp,   ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  INTEGER, INTENT(IN) :: root, comm
  !
  CALL mp_bcast( rho%of_g, root, comm )
  CALL mp_bcast( rho%of_r, root, comm )
  IF ( xclib_dft_is('meta') .OR. lxdm) THEN
     CALL mp_bcast( rho%kin_g, root, comm )
     CALL mp_bcast( rho%kin_r, root, comm )
  END IF
  IF (lda_plus_u_co)  CALL mp_bcast( rho%ns,    root, comm )
  IF (lda_plus_u_cob) CALL mp_bcast( rho%nsb,   root, comm )
  IF (lda_plus_u_nc)  CALL mp_bcast( rho%ns_nc, root, comm )
  IF (lda_plus_u_v)   CALL mp_bcast( rho%nsg,   root, comm )
  IF (okpaw)          CALL mp_bcast( rho%bec,   root, comm )
  IF (sic) THEN
     CALL mp_bcast ( rho%pol_r, root, comm )
     CALL mp_bcast ( rho%pol_g, root, comm )
  END IF
  !
END SUBROUTINE bcast_scf_type
!
!---------------------------------------------------------------------------
SUBROUTINE rhoz_or_updw( rho, sp, dir )
  !--------------------------------------------------------------------------
  !! Converts rho(up,dw) into rho(up+dw,up-dw) if dir='->rhoz' and
  !! vice versa if dir='->updw'.
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  !! the charge density
  CHARACTER(LEN=*), INTENT(IN) :: dir
  !! direction of the conversion
  CHARACTER(LEN=*), INTENT(IN) :: sp
  !! g-space ('only_g') or r-space ('only_r') or both
  !
  ! ... local variables
  !
  INTEGER :: ir, dfftp_nnr
  REAL(DP) :: vi
  !
  IF ( nspin /= 2 ) RETURN
  !
 !$acc data present_or_copy(rho)
  !
  vi = 0._dp
  IF (dir == '->updw')  vi = 0.5_dp
  IF (dir == '->rhoz')  vi = 1.0_dp
  IF (vi  == 0._dp)  CALL errore( 'rhoz_or_updw', 'wrong input', 1 )
  !
  IF ( sp /= 'only_g' ) THEN
     !
     dfftp_nnr = dfftp%nnr
    !$acc parallel loop present_or_copy(rho%of_r)
     DO ir = 1, dfftp_nnr
        rho%of_r(ir,1) = ( rho%of_r(ir,1) + rho%of_r(ir,nspin) ) * vi
        rho%of_r(ir,nspin) = rho%of_r(ir,1) - rho%of_r(ir,nspin) * vi * 2._dp
     ENDDO
     !
  ENDIF
  IF ( sp /= 'only_r' ) THEN
     !
    !$acc parallel loop present_or_copy(rho%of_g)
     DO ir = 1, ngm
        rho%of_g(ir,1) = ( rho%of_g(ir,1) + rho%of_g(ir,nspin) ) * vi
        rho%of_g(ir,nspin) = rho%of_g(ir,1) - rho%of_g(ir,nspin) * vi * 2._dp
     ENDDO
     !
  ENDIF
  !
 !$acc end data
  !
  RETURN
  !
  END SUBROUTINE rhoz_or_updw
  !
  !
END MODULE scf
