!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE scf
  !--------------------------------------------------------------------------
  !! This module contains variables and auxiliary routines needed for
  !! the self-consistent cycle.  
  !
  USE kinds,           ONLY : DP
  USE lsda_mod,        ONLY : nspin
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax, lda_plus_u_kind, ldmx, &
                              ldmx_b, is_hubbard_back
  USE ions_base,       ONLY : nat
  USE buffers,         ONLY : open_buffer, close_buffer, get_buffer, save_buffer
  USE funct,           ONLY : dft_is_meta
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : invfft
  USE gvect,           ONLY : ngm
  USE gvecs,           ONLY : ngms
  USE ions_base,       ONLY : ntyp => nsp
  USE paw_variables,   ONLY : okpaw
  USE uspp_param,      ONLY : nhm
  USE extfield,        ONLY : dipfield, emaxpos, eopreg, edir
  USE control_flags,   ONLY : lxdm
  !
  SAVE
  !
  ! Details of PAW implementation:
  !
  ! NOTE1: scf_type is used for two different quantities: density and potential.
  !        These correspond, for PAW, to becsum and D coefficients.
  !        Due to interference with the ultrasoft routines only the becsum part
  !        is stored in the structure (at the moment).
  !        This only holds for scf_type; mix_type is not affected.
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
     REAL(DP),    ALLOCATABLE :: bec(:,:,:)
     !! the PAW hamiltonian elements
  END TYPE scf_type
  !
  TYPE mix_type
     !! For quantities directly involved in the mixing
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)
     !! the charge density in G-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:)
     !! the charge density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)
     !! the DFT+U occupation matrix
     REAL(DP),    ALLOCATABLE :: nsb(:,:,:,:)
     !! the DFT+U occupation matrix (background states)
     COMPLEX(DP), ALLOCATABLE :: ns_nc(:,:,:,:)
     !! the DFT+U occupation matrix noncollinear case 
     REAL(DP),    ALLOCATABLE :: bec(:,:,:)
     !! PAW corrections to hamiltonian
     REAL(DP) :: el_dipole
     !! electrons dipole
  END TYPE mix_type
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
  !
  INTEGER, PRIVATE  :: record_length, &
                       rlen_rho=0,  rlen_kin=0,  rlen_ldaU=0,  rlen_bec=0,&
                       rlen_dip=0, rlen_ldaUb=0, &
                       start_rho=0, start_kin=0, start_ldaU=0, start_bec=0, &
                       start_dipole=0, start_ldaUb=0
  INTEGER :: nt
  ! DFT+U, colinear and noncolinear cases
  LOGICAL, PRIVATE :: lda_plus_u_co  ! collinear case
  LOGICAL, PRIVATE :: lda_plus_u_cob ! collinear case (background states)
  LOGICAL, PRIVATE :: lda_plus_u_nc  ! noncollinear case
  COMPLEX(DP), PRIVATE, ALLOCATABLE:: io_buffer(:)
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
   IF (dft_is_meta() .OR. lxdm) THEN
      ALLOCATE( rho%kin_r(dfftp%nnr,nspin) )
      ALLOCATE( rho%kin_g(ngm,nspin) )
   ELSE
      ALLOCATE( rho%kin_r(1,1) )
      ALLOCATE( rho%kin_g(1,1) )
   ENDIF
   !
   lda_plus_u_co  = lda_plus_u .AND. .NOT. ( nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_nc  = lda_plus_u .AND.       ( nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_cob = .FALSE.
   IF (lda_plus_u_co) THEN
      DO nt = 1, ntyp
         IF (is_hubbard_back(nt)) lda_plus_u_cob = .TRUE.
      ENDDO
   ENDIF
   !
   IF (lda_plus_u_co)  ALLOCATE( rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
   IF (lda_plus_u_cob) ALLOCATE( rho%nsb(ldmx_b,ldmx_b,nspin,nat) )
   IF (lda_plus_u_nc)  ALLOCATE( rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
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
   IF (ALLOCATED(rho%bec)  )  DEALLOCATE( rho%bec   )
   !
   RETURN
   !
 END SUBROUTINE destroy_scf_type
 !
 !
 !----------------------------------------------------
 SUBROUTINE create_mix_type( rho )
   !--------------------------------------------------
   !
   IMPLICIT NONE
   !
   TYPE(mix_type) :: rho
   !
   ALLOCATE( rho%of_g(ngms,nspin) )
   !
   rho%of_g = 0._dp
   !
   IF (dft_is_meta() .OR. lxdm) THEN
      ALLOCATE( rho%kin_g(ngms,nspin) )
      rho%kin_g = 0._dp
   ENDIF
   !
   lda_plus_u_co = lda_plus_u .AND. .NOT. (nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_nc = lda_plus_u .AND.       (nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_cob = .FALSE.
   IF (lda_plus_u_co) THEN
      DO nt = 1, ntyp
         IF (is_hubbard_back(nt)) lda_plus_u_cob = .TRUE.
      ENDDO
   ENDIF
   !
   IF (lda_plus_u_nc) THEN
      ALLOCATE( rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
      rho%ns_nc = 0._dp
   ENDIF
   !
   IF (lda_plus_u_co) THEN
      ALLOCATE( rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
      rho%ns = 0._dp
   ENDIF
   !
   IF (lda_plus_u_cob) THEN
      ALLOCATE( rho%nsb(ldmx_b,ldmx_b,nspin,nat) )
      rho%nsb = 0._dp
   ENDIF
   !
   IF (okpaw) THEN
      ALLOCATE( rho%bec(nhm*(nhm+1)/2,nat,nspin) )
      rho%bec = 0._dp
   ENDIF
   !
   rho%el_dipole = 0._dp
   !
   RETURN
   !
 END SUBROUTINE create_mix_type
 !
 !
 !------------------------------------------------------
 SUBROUTINE destroy_mix_type( rho )
   !----------------------------------------------------
   !! Deallocates a \(\text{mix_type}\) object.
   !
   IMPLICIT NONE
   !
   TYPE(mix_type) :: rho
   !
   IF (ALLOCATED(rho%of_g) )  DEALLOCATE( rho%of_g  )
   IF (ALLOCATED(rho%kin_g))  DEALLOCATE( rho%kin_g )
   IF (ALLOCATED(rho%ns)   )  DEALLOCATE( rho%ns    )
   IF (ALLOCATED(rho%nsb)  )  DEALLOCATE( rho%nsb   )
   IF (ALLOCATED(rho%ns_nc))  DEALLOCATE( rho%ns_nc )
   IF (ALLOCATED(rho%bec)  )  DEALLOCATE( rho%bec   )
   !
   RETURN
   !
 END SUBROUTINE destroy_mix_type
 !
 !
 !-----------------------------------------------------
 SUBROUTINE assign_scf_to_mix_type( rho_s, rho_m )
   !----------------------------------------------------
   !! It fills a \(\text{mix_type}\) object starting from a
   !! \(\text{scf_type}\) one.
   !
   IMPLICIT NONE
   !
   TYPE(scf_type), INTENT(IN)  :: rho_s
   TYPE(mix_type), INTENT(INOUT) :: rho_m
   !
   REAL(DP) :: e_dipole
   !
   rho_m%of_g(1:ngms,:) = rho_s%of_g(1:ngms,:)
   !
   IF (dft_is_meta() .OR. lxdm) rho_m%kin_g(1:ngms,:) = rho_s%kin_g(1:ngms,:)
   IF (lda_plus_u_nc)  rho_m%ns_nc  = rho_s%ns_nc
   IF (lda_plus_u_co)  rho_m%ns     = rho_s%ns
   IF (lda_plus_u_cob) rho_m%nsb    = rho_s%nsb
   IF (okpaw)          rho_m%bec    = rho_s%bec
   !
   IF (dipfield) THEN
      CALL compute_el_dip( emaxpos, eopreg, edir, rho_s%of_r(:,1), e_dipole )
      rho_m%el_dipole = e_dipole
   ENDIF
   !
   RETURN
   !
 END SUBROUTINE assign_scf_to_mix_type
 !
 !
 !-----------------------------------------------------------------
 SUBROUTINE assign_mix_to_scf_type( rho_m, rho_s )
   !----------------------------------------------------------------
   !! It fills a \(\text{scf_type}\) object starting from a 
   !! \(\text{mix_type}\) one.
   !
   USE wavefunctions,        ONLY : psic
   USE control_flags,        ONLY : gamma_only
   !
   IMPLICIT NONE
   !
   TYPE(mix_type), INTENT(IN) :: rho_m
   TYPE(scf_type), INTENT(INOUT) :: rho_s
   !
   INTEGER :: is
   !   
   rho_s%of_g(1:ngms,:) = rho_m%of_g(1:ngms,:)
   ! define rho_s%of_r 
   !
   DO is = 1, nspin
      psic(:) = ( 0.D0, 0.D0 )
      psic(dfftp%nl(:)) = rho_s%of_g(:,is)
      IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho_s%of_g(:,is) )
      CALL invfft( 'Rho', psic, dfftp )
      rho_s%of_r(:,is) = psic(:)
   ENDDO
   !
   IF (dft_is_meta() .OR. lxdm) THEN
      rho_s%kin_g(1:ngms,:) = rho_m%kin_g(:,:)
      ! define rho_s%kin_r 
      DO is = 1, nspin
         psic(:) = ( 0.D0, 0.D0 )
         psic(dfftp%nl(:)) = rho_s%kin_g(:,is)
         IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho_s%kin_g(:,is) )
         CALL invfft( 'Rho', psic, dfftp )
         rho_s%kin_r(:,is) = psic(:)
      ENDDO
   ENDIF
   !
   IF (lda_plus_u_nc)  rho_s%ns_nc(:,:,:,:) = rho_m%ns_nc(:,:,:,:)
   IF (lda_plus_u_co)  rho_s%ns(:,:,:,:)    = rho_m%ns(:,:,:,:)
   IF (lda_plus_u_cob) rho_s%nsb(:,:,:,:)   = rho_m%nsb(:,:,:,:)
   IF (okpaw)          rho_s%bec(:,:,:)     = rho_m%bec(:,:,:)
   !
   RETURN
   !
 END SUBROUTINE assign_mix_to_scf_type
 !
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
  IF (dft_is_meta() .OR. lxdm) THEN
     Y%kin_r = X%kin_r
     Y%kin_g = X%kin_g
  ENDIF
  !
  IF (lda_plus_u_nc)  Y%ns_nc = X%ns_nc
  IF (lda_plus_u_co)  Y%ns    = X%ns
  IF (lda_plus_u_cob) Y%nsb   = X%nsb
  IF (okpaw)          Y%bec   = X%bec
  !
  RETURN
  !
 END SUBROUTINE scf_type_COPY
 !
 !
 !----------------------------------------------------------------------------
 SUBROUTINE mix_type_AXPY( A, X, Y )
  !----------------------------------------------------------------------------
  !! Works like daxpy for \(\text{scf_type}\) variables: \(Y = A\cdot X + Y\)
  ! NB: A is a REAL(DP) number
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: A
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  !
  Y%of_g = Y%of_g  + A * X%of_g
  !
  IF (dft_is_meta() .OR. lxdm) Y%kin_g     = Y%kin_g     + A * X%kin_g
  IF (lda_plus_u_nc)           Y%ns_nc     = Y%ns_nc     + A * X%ns_nc
  IF (lda_plus_u_co)           Y%ns        = Y%ns        + A * X%ns
  IF (lda_plus_u_cob)          Y%nsb       = Y%nsb       + A * X%nsb
  IF (okpaw)                   Y%bec       = Y%bec       + A * X%bec
  IF (dipfield)                Y%el_dipole = Y%el_dipole + A * X%el_dipole
  !
  RETURN
  !
 END SUBROUTINE mix_type_AXPY
 !
 !
 !----------------------------------------------------------------------------
 SUBROUTINE mix_type_COPY( X, Y )
  !----------------------------------------------------------------------------
  !! Works like DCOPY for \(\text{mix_type}\) copy variables: \(Y = X\).
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  !
  Y%of_g  = X%of_g
  !
  IF (dft_is_meta() .OR. lxdm) Y%kin_g     = X%kin_g
  IF (lda_plus_u_nc)           Y%ns_nc     = X%ns_nc
  IF (lda_plus_u_co)           Y%ns        = X%ns
  IF (lda_plus_u_cob)          Y%nsb       = X%nsb
  IF (okpaw)                   Y%bec       = X%bec
  IF (dipfield)                Y%el_dipole = X%el_dipole
  !
  RETURN
  !
 END SUBROUTINE mix_type_COPY
 !
 !
 !----------------------------------------------------------------------------
 SUBROUTINE mix_type_SCAL( A, X )
  !----------------------------------------------------------------------------
  !! Works like DSCAL for \(\text{mix_type}\) copy variables: \(X = A \cdot X\)  
  !! NB: A is a REAL(DP) number
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP),       INTENT(IN)    :: A
  TYPE(mix_type), INTENT(INOUT) :: X
  !
  X%of_g(:,:) = A * X%of_g(:,:)
  !
  IF (dft_is_meta() .OR. lxdm) X%kin_g     = A * X%kin_g
  IF (lda_plus_u_nc)           X%ns_nc     = A * X%ns_nc
  IF (lda_plus_u_co)           X%ns        = A * X%ns
  IF (lda_plus_u_cob)          X%nsb       = A * X%nsb
  IF (okpaw)                   X%bec       = A * X%bec
  IF (dipfield)                X%el_dipole = A * X%el_dipole
  !
  RETURN
  !
 END SUBROUTINE mix_type_SCAL
 !
 !
 !---------------------------------------------------------------------
 SUBROUTINE high_frequency_mixing( rhoin, input_rhout, alphamix )
   !-------------------------------------------------------------------
   !
   USE wavefunctions,    ONLY : psic
   USE control_flags,    ONLY : gamma_only
   !
   IMPLICIT NONE
   !
   TYPE (scf_type), INTENT(INOUT) :: rhoin
   TYPE (scf_type), INTENT(IN) :: input_rhout
   REAL(DP), INTENT(IN) :: alphamix
   !
   ! ... local variable
   !
   INTEGER :: is
   !
   IF (ngms < ngm ) THEN
      !
      rhoin%of_g = rhoin%of_g + alphamix * (input_rhout%of_g-rhoin%of_g)
      rhoin%of_g(1:ngms,1:nspin) = (0.d0,0.d0)
      ! define rho_s%of_r 
      DO is = 1, nspin
         psic(:) = ( 0.D0, 0.D0 )
         psic(dfftp%nl(:)) = rhoin%of_g(:,is)
         IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rhoin%of_g(:,is) )
         CALL invfft( 'Rho', psic, dfftp )
         rhoin%of_r(:,is) = psic(:)
      ENDDO
      !
      IF (dft_is_meta() .OR. lxdm) THEN
         rhoin%kin_g = rhoin%kin_g + alphamix * ( input_rhout%kin_g-rhoin%kin_g)
         rhoin%kin_g(1:ngms,1:nspin) = (0.d0,0.d0)
         ! define rho_s%of_r 
         DO is = 1, nspin
            psic(:) = ( 0.D0, 0.D0 )
            psic(dfftp%nl(:)) = rhoin%kin_g(:,is)
            IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rhoin%kin_g(:,is) )
            CALL invfft( 'Rho', psic, dfftp )
            rhoin%kin_r(:,is) = psic(:)
         ENDDO
      ENDIF
      !
   ELSE
      !
      rhoin%of_g(:,:)= (0.d0,0.d0)
      rhoin%of_r(:,:)= 0.d0
      IF (dft_is_meta() .OR. lxdm) THEN
         rhoin%kin_g(:,:)= (0.d0,0.d0)
         rhoin%kin_r(:,:)= 0.d0
      ENDIF
      !
   ENDIF
   !
   IF (lda_plus_u_nc)  rhoin%ns_nc(:,:,:,:) = 0.d0
   IF (lda_plus_u_co)  rhoin%ns(:,:,:,:)    = 0.d0
   IF (lda_plus_u_cob) rhoin%nsb(:,:,:,:)   = 0.d0
   !
   RETURN
   !
 END SUBROUTINE high_frequency_mixing 
 !
 !
 !------------------------------------------------------------------------
 SUBROUTINE open_mix_file( iunit, extension, exst )
   !------------------------------------------------------------------------
   !
   USE control_flags,  ONLY : io_level
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: extension
   INTEGER, INTENT(IN) :: iunit
   LOGICAL :: exst
   !
   ! define lengths (in real numbers) of different record chunks
   rlen_rho = 2 * ngms * nspin
   IF (dft_is_meta() .OR. lxdm) rlen_kin  = 2 * ngms * nspin
   IF (lda_plus_u_co)           rlen_ldaU = (2*Hubbard_lmax+1)**2 *nspin*nat
   IF (lda_plus_u_cob)          rlen_ldaUb = (ldmx_b)**2 *nspin*nat
   IF (lda_plus_u_nc)           rlen_ldaU = 2 * (2*Hubbard_lmax+1)**2 *nspin*nat
   IF (okpaw)                   rlen_bec  = (nhm*(nhm+1)/2) * nat * nspin
   IF (dipfield)                rlen_dip  = 1
   !
   ! define the starting point of the different chunks. Beware: each starting point
   ! is the index of a COMPLEX array. When real arrays with odd dimension are copied
   ! to/from the complex array io_buffer, the last complex number will be half-filled
   ! but must still be counted as one!
   start_rho    = 1
   start_kin    = start_rho  + rlen_rho / 2
   start_ldaU   = start_kin  + rlen_kin / 2
   IF (lda_plus_u_cob) THEN
      start_ldaUb = start_ldaU + ( rlen_ldaU + 1 ) / 2
      start_bec = start_ldaUb + ( rlen_ldaUb + 1 ) / 2
   ELSE
      start_bec = start_ldaU + ( rlen_ldaU + 1 ) / 2
   ENDIF
   start_dipole = start_bec  + ( rlen_bec + 1 ) / 2
   !
   ! define total record length, in complex numbers
   record_length = start_dipole + rlen_dip - 1
   !
   ! open file and allocate io_buffer
   CALL open_buffer( iunit, extension, record_length, io_level, exst )
   !
   ALLOCATE( io_buffer(record_length) )
   ! setting to zero -prevents trouble with "holes" due to odd dimensions of real
   ! arrays
   io_buffer(:) = (0.0_dp, 0.0_dp)
   !
   RETURN
   !
 END SUBROUTINE open_mix_file
 !
 !
 !------------------------------------------------------------------------
 SUBROUTINE close_mix_file( iunit, stat )
   !---------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iunit
   CHARACTER(LEN=*), INTENT(IN) :: stat
   !
   DEALLOCATE( io_buffer )
   !
   CALL close_buffer( iunit, TRIM(stat) ) 
   !
   RETURN
   !
 END SUBROUTINE close_mix_file
 !
 !
 !------------------------------------------------------------
 SUBROUTINE davcio_mix_type( rho, iunit, record, iflag )
   !----------------------------------------------------------
   !
   IMPLICIT NONE
   !
   TYPE(mix_type) :: rho
   INTEGER, INTENT(IN) :: iunit, record, iflag
   !
   IF (iflag > 0) THEN
      !
      CALL DCOPY(rlen_rho,rho%of_g,1,io_buffer(start_rho),1)
      !
      IF (dft_is_meta() .OR. lxdm) CALL DCOPY(rlen_kin, rho%kin_g,1,io_buffer(start_kin), 1)
      IF (lda_plus_u_nc)           CALL DCOPY(rlen_ldaU,rho%ns_nc,1,io_buffer(start_ldaU),1)
      IF (lda_plus_u_co)           CALL DCOPY(rlen_ldaU,rho%ns,   1,io_buffer(start_ldaU),1)
      IF (lda_plus_u_cob)          CALL DCOPY(rlen_ldaUb,rho%nsb, 1,io_buffer(start_ldaUb),1)
      IF (okpaw)                   CALL DCOPY(rlen_bec, rho%bec,  1,io_buffer(start_bec), 1)
      !
      IF (dipfield) io_buffer(start_dipole) = CMPLX( rho%el_dipole, 0.0_dp )
      !
      CALL save_buffer( io_buffer, record_length, iunit, record )   
      !
   ELSEIF (iflag < 0 ) THEN
      !
      CALL get_buffer( io_buffer, record_length, iunit, record )
      !
      CALL DCOPY(rlen_rho,io_buffer(start_rho),1,rho%of_g,1)
      !
      IF (dft_is_meta() .OR. lxdm) CALL DCOPY(rlen_kin, io_buffer(start_kin), 1,rho%kin_g,1)
      IF (lda_plus_u_co)           CALL DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns,   1)
      IF (lda_plus_u_cob)          CALL DCOPY(rlen_ldaUb,io_buffer(start_ldaUb),1,rho%nsb,1)
      IF (lda_plus_u_nc)           CALL DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns_nc,1)
      IF (okpaw)                   CALL DCOPY(rlen_bec, io_buffer(start_bec), 1,rho%bec,  1)
      !
      IF (dipfield) rho%el_dipole = DBLE( io_buffer(start_dipole) )
      !
   ENDIF
   !
 END SUBROUTINE davcio_mix_type
 !
 !
 !-----------------------------------------------------------------------------------
FUNCTION rho_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------------
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy.
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : e2, tpi, fpi
  USE cell_base,       ONLY : omega, tpiba2
  USE gvect,           ONLY : gg, gstart
  USE control_flags,   ONLY : gamma_only
  USE paw_onecenter,   ONLY : paw_ddot
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first density matrix
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second density matrix
  INTEGER, INTENT(IN) :: gf
  !! points delimiter
  REAL(DP) :: rho_ddot
  !! output: see function comments
  !
  ! ... local variables
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  fac = e2 * fpi / tpiba2
  !
  rho_ddot = 0.D0
  !
  DO ig = gstart, gf
     rho_ddot = rho_ddot + REAL(CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1), DP) / gg(ig)
  ENDDO
  !
  rho_ddot = fac*rho_ddot
  !
  IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
  !
  IF ( nspin >= 2 )  THEN
     fac = e2*fpi / tpi**2  ! lambda=1 a.u.
     IF ( gstart == 2 ) THEN
        rho_ddot = rho_ddot + &
                fac * SUM(REAL(CONJG( rho1%of_g(1,2:nspin))*(rho2%of_g(1,2:nspin) ), DP))
     ENDIF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        rho_ddot = rho_ddot + &
              fac * SUM(REAL(CONJG( rho1%of_g(ig,2:nspin))*(rho2%of_g(ig,2:nspin) ), DP))
     ENDDO
  ENDIF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL mp_sum( rho_ddot, intra_bgrp_comm )
  !
  IF (dft_is_meta()) rho_ddot = rho_ddot + tauk_ddot( rho1, rho2, gf )
  IF (lda_plus_u )   rho_ddot = rho_ddot + ns_ddot( rho1, rho2 )
  ! 
  ! Beware: paw_ddot has a hidden parallelization on all processors
  !         it must be called on all processors or else it will hang
  ! Beware: commented out because it yields too often negative values
  ! IF (okpaw)         rho_ddot = rho_ddot + paw_ddot(rho1%bec, rho2%bec)
  !
  IF (dipfield) rho_ddot = rho_ddot + (e2/2.0_DP)*(rho1%el_dipole * rho2%el_dipole)*omega/fpi
  !
  RETURN
  !
END FUNCTION rho_ddot
!
!
!----------------------------------------------------------------------------
FUNCTION tauk_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy - kinetic density
  !! version.
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE control_flags, ONLY : gamma_only
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first kinetic density
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second kinetic density
  INTEGER, INTENT(IN) :: gf
  !! point delimiter
  REAL(DP) :: tauk_ddot
  !! output: see function comments
  !
  ! ... local variables
  !
  REAL(DP) :: fac
  INTEGER :: ig
  !
  tauk_ddot = 0.D0
  !
  !  write (*,*) rho1%kin_g(1:4,1)
  !  if (.true. ) stop
  !
  DO ig = gstart, gf
     tauk_ddot = tauk_ddot + DBLE( CONJG( rho1%kin_g(ig,1) )*rho2%kin_g(ig,1) ) 
  ENDDO
  !
  IF ( nspin==1 .AND. gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
  !
  ! ... G=0 term
  !
  IF ( gstart == 2 ) THEN
     tauk_ddot = tauk_ddot + DBLE( CONJG( rho1%kin_g(1,1) ) * rho2%kin_g(1,1) )
  ENDIF
  !
  IF ( nspin >= 2 ) THEN
     DO ig = gstart, gf
        tauk_ddot = tauk_ddot + &
          SUM( REAL( CONJG( rho1%kin_g(1,2:nspin))*(rho2%kin_g(1,2:nspin) ), DP))
     ENDDO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     IF ( gstart == 2 ) THEN
        tauk_ddot = tauk_ddot + &
          SUM(REAL(CONJG( rho1%kin_g(1,1:nspin))*(rho2%kin_g(1,1:nspin) ), DP))
     ENDIF
     !
     IF ( nspin == 2 ) tauk_ddot = 0.5D0 *  tauk_ddot 
  ENDIF
  !
  fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
  !
  tauk_ddot = fac * tauk_ddot * omega * 0.5D0
  !
  CALL mp_sum( tauk_ddot, intra_bgrp_comm )
  !
  RETURN
  !
END FUNCTION tauk_ddot
!
!
!----------------------------------------------------------------------------
FUNCTION ns_ddot( rho1, rho2 )
  !---------------------------------------------------------------------------
  !! Calculates \(U/2 \sum_i \text{ns1}(i)\ \text{ns2}(i)\) used as an estimate
  !! of the self-consistency error on the DFT+U correction to the energy.
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_l, Hubbard_U, Hubbard_U_back, ldim_back, &
                        lda_plus_u_kind, is_hubbard, is_hubbard_back
  USE ions_base, ONLY : nat, ityp
  !
  IMPLICIT NONE  
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first Hubbard ns
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second Hubbard ns
  REAL(DP) :: ns_ddot
  !! output: see function comments
  !
  ! ... local variables
  !
  INTEGER :: na, nt, m1, m2
  !
  ns_ddot = 0.D0
  !
  IF (lda_plus_u_kind.EQ.2) RETURN
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( is_hubbard(nt) ) THEN
        !
        m1 = 2 * Hubbard_l(nt) + 1
        m2 = 2 * Hubbard_l(nt) + 1
        !
        IF (nspin == 4) THEN
          ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( CONJG(rho1%ns_nc(:m1,:m2,:nspin,na))*rho2%ns_nc(:m1,:m2,:nspin,na) )
        ELSE
          ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( rho1%ns(:m1,:m2,:nspin,na)*rho2%ns(:m1,:m2,:nspin,na) )
        ENDIF
        !
     ENDIF
     !
     ! Background part
     !
     IF ( is_hubbard_back(nt) .AND. lda_plus_u_kind.EQ.0 ) THEN
        !
        m1 = ldim_back(nt)
        m2 = ldim_back(nt)
        !
        ns_ddot = ns_ddot + 0.5D0 * Hubbard_U_back(nt) * &
                SUM( rho1%nsb(:m1,:m2,:nspin,na)*rho2%nsb(:m1,:m2,:nspin,na) )
        !
     ENDIF
     !
  ENDDO
  !
  IF ( nspin == 1 ) ns_ddot = 2.D0*ns_ddot
  !
  RETURN
  !
END FUNCTION ns_ddot
!
!----------------------------------------------------------------------------
FUNCTION nsg_ddot( nsg1, nsg2, nspin )
  !----------------------------------------------------------------------------
  !! Calculates \(U/2 \sum_i \text{nsg1}(i) \text{nsg2}(i)\)
  !! used as an estimate of the self-consistency error on the 
  !! DFT+U+V correction to the energy
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_V, max_num_neighbors, is_hubbard, &
                        is_hubbard_back, ldim_u, neighood, ldmx,  &
                        ldmx_tot, at_sc
  USE ions_base, ONLY : nat, ityp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: nsg1(ldmx_tot,ldmx_tot,max_num_neighbors,nat,nspin), &
                             nsg2(ldmx_tot,ldmx_tot,max_num_neighbors,nat,nspin)
  INTEGER,  INTENT(IN) :: nspin
  REAL(DP)             :: nsg_ddot
  !
  INTEGER :: na1, nt1, m1, m2, na2, nt2, viz, equiv_na2, i_type
  INTEGER, EXTERNAL :: find_viz, type_interaction
  !
  nsg_ddot = 0.D0
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( (is_hubbard(nt1) .OR. is_hubbard_back(nt1)) .AND.ldim_u(nt1).GT.0 ) THEN
        DO viz = 1, neighood(na1)%num_neigh
           na2 = neighood(na1)%neigh(viz)
           equiv_na2 = at_sc(na2)%at
           nt2 = ityp(equiv_na2)
           IF ((Hubbard_V(na1,na2,1) .NE. 0.d0) .OR. &
               (Hubbard_V(na1,na2,2) .NE. 0.d0) .OR. &
               (Hubbard_V(na1,na2,3) .NE. 0.d0) .OR. &
               (Hubbard_V(na1,na2,4) .NE. 0.d0) ) THEN
              DO m1=1,ldim_u(nt1)
                 DO m2=1,ldim_u(nt2)
                    i_type = type_interaction(na1,m1,equiv_na2,m2)
                    nsg_ddot = nsg_ddot + 0.5D0 * ABS(Hubbard_V(na1,na2,i_type)) * &
                      REAL(SUM( nsg1(m2,m1,viz, na1,:nspin)*  &
                                CONJG(nsg2(m2,m1,viz, na1,:nspin) ) ) )
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  IF ( nspin == 1 ) nsg_ddot = 2.D0*nsg_ddot
  !
  RETURN
  !
END FUNCTION nsg_ddot
!
!----------------------------------------------------------------------------
FUNCTION local_tf_ddot( rho1, rho2, ngm0 )
  !----------------------------------------------------------------------------
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy - version 
  !! for the case with local-density dependent TF preconditioning to drho.
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : e2, fpi
  USE cell_base,       ONLY : omega, tpiba2
  USE gvect,           ONLY : gg, gstart
  USE control_flags,   ONLY : gamma_only
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngm0
  !! input length (local number of smooth vectors)
  COMPLEX(DP), INTENT(IN) :: rho1(ngm0)
  !! see main comment
  COMPLEX(DP), INTENT(IN) :: rho2(ngm0)
  !! see main comment
  REAL(DP) :: local_tf_ddot
  !! see main comment
  !
  ! ... local variables
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  local_tf_ddot = 0.D0
  !
  fac = e2 * fpi / tpiba2
  !
  !$omp parallel do reduction(+:local_tf_ddot)
  DO ig = gstart, ngm0
     local_tf_ddot = local_tf_ddot + DBLE( CONJG(rho1(ig))*rho2(ig) ) / gg(ig)
  END DO
  !$omp end parallel do
  !
  local_tf_ddot = fac * local_tf_ddot * omega * 0.5D0
  !
  IF (gamma_only) local_tf_ddot = 2.D0 * local_tf_ddot
  !
  CALL mp_sum( local_tf_ddot, intra_bgrp_comm )
  !
  RETURN
  !
END FUNCTION local_tf_ddot
!
!
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
  IF ( dft_is_meta() .OR. lxdm) THEN
     CALL mp_bcast( rho%kin_g, root, comm )
     CALL mp_bcast( rho%kin_r, root, comm )
  END IF
  IF (lda_plus_u_co)  CALL mp_bcast( rho%ns,    root, comm )
  IF (lda_plus_u_cob) CALL mp_bcast( rho%nsb,   root, comm )
  IF (lda_plus_u_nc)  CALL mp_bcast( rho%ns_nc, root, comm )
  IF (okpaw)          CALL mp_bcast( rho%bec,   root, comm )
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE rhoz_or_updw( rho, sp, dir )
  !--------------------------------------------------------------------------
  !! Converts rho(up,dw) into rho(up+dw,up-dw) if dir='->rhoz' and
  !! vice versa if dir='->updw'.
  !
  USE gvect,  ONLY : ngm
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
  INTEGER :: ir
  REAL(DP) :: vi
  !
  IF ( nspin /= 2 ) RETURN
  !
  vi = 0._dp
  IF (dir == '->updw')  vi = 0.5_dp
  IF (dir == '->rhoz')  vi = 1.0_dp
  IF (vi  == 0._dp)  CALL errore( 'rhoz_or_updw', 'wrong input', 1 )
  !
  IF ( sp /= 'only_g' ) THEN
     !
     DO ir = 1, dfftp%nnr  
        rho%of_r(ir,1) = ( rho%of_r(ir,1) + rho%of_r(ir,nspin) ) * vi
        rho%of_r(ir,nspin) = rho%of_r(ir,1) - rho%of_r(ir,nspin) * vi * 2._dp
     ENDDO
     !
  ENDIF
  IF ( sp /= 'only_r' ) THEN
     !
     DO ir = 1, ngm
        rho%of_g(ir,1) = ( rho%of_g(ir,1) + rho%of_g(ir,nspin) ) * vi
        rho%of_g(ir,nspin) = rho%of_g(ir,1) - rho%of_g(ir,nspin) * vi * 2._dp
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE
  !
  !
END MODULE scf
