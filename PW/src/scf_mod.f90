!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE scf
  !  
  !  This module contains variables and auxiliary routines needed for
  !  the self-consistent cycle
  !
  !  ROUTINES: allocate_scf_type
  !
  USE kinds,      ONLY : DP
  !
  USE lsda_mod,     ONLY : nspin
  USE ldaU,         ONLY : lda_plus_u, Hubbard_lmax
  USE ions_base,    ONLY : nat
  USE buffers,      ONLY : open_buffer, close_buffer, get_buffer, save_buffer
  USE funct,        ONLY : dft_is_meta
  USE fft_base,     ONLY : dfftp
  USE fft_interfaces,ONLY: invfft
  USE gvect,        ONLY : ngm
  USE gvecs,        ONLY : ngms
  USE paw_variables,ONLY : okpaw
  USE uspp_param,   ONLY : nhm
  USE extfield,     ONLY : dipfield, emaxpos, eopreg, edir
  USE control_flags,ONLY : lxdm
  !
  SAVE
  !
! Details of PAW implementation:
! NOTE: scf_type is used for two different quantities: density and potential.
!       These correspond, for PAW, to becsum and D coefficients.
!       Due to interference with the ultrasoft routines only the becsum part
!       is stored in the structure (at the moment).
!       This only holds for scf_type; mix_type is not affected.
! NOTE: rho%bec is different from becsum for two reasons:
!       1. rho%bec is mixed, while becsum is not
!       2. for npool > 1 rho%bec is collected, becsum is not
!          ( this is necessary to make the stress work)

  TYPE scf_type
     REAL(DP),   ALLOCATABLE :: of_r(:,:)  ! the charge density in R-space
     COMPLEX(DP),ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     REAL(DP),   ALLOCATABLE :: kin_r(:,:) ! the kinetic energy density in R-space
     COMPLEX(DP),ALLOCATABLE :: kin_g(:,:) ! the kinetic energy density in G-space
     REAL(DP),   ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix
     COMPLEX(DP),ALLOCATABLE :: ns_nc(:,:,:,:)!     ---       noncollinear case
     REAL(DP),   ALLOCATABLE :: bec(:,:,:) ! the PAW hamiltonian elements
  END TYPE scf_type
  !
  TYPE mix_type
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:) ! the charge density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix 
     COMPLEX(DP), ALLOCATABLE :: ns_nc(:,:,:,:)!     ---     noncollinear case 
     REAL(DP),    ALLOCATABLE :: bec(:,:,:) ! PAW corrections to hamiltonian
     REAL(DP)                 :: el_dipole  ! electrons dipole
  END TYPE mix_type

  type (scf_type) :: rho  ! the charge density and its other components
  type (scf_type) :: v    ! the scf potential
  type (scf_type) :: vnew ! used to correct the forces

  REAL(DP) :: v_of_0    ! vltot(G=0)      
  REAL(DP), ALLOCATABLE :: &
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grid)
       rho_core(:),    &! the core charge in real space
       kedtau(:,:)      ! position dependent kinetic energy enhancement factor
  COMPLEX(DP), ALLOCATABLE :: &
       rhog_core(:)     ! the core charge in reciprocal space

  INTEGER, PRIVATE  :: record_length, &
                       rlen_rho=0,  rlen_kin=0,  rlen_ldaU=0,  rlen_bec=0,&
                       rlen_dip=0, &
                       start_rho=0, start_kin=0, start_ldaU=0, start_bec=0, &
                       start_dipole=0
  ! DFT+U, colinear and noncolinear cases
  LOGICAL, PRIVATE :: lda_plus_u_co, lda_plus_u_nc
  COMPLEX(DP), PRIVATE, ALLOCATABLE:: io_buffer(:)
CONTAINS

 SUBROUTINE create_scf_type ( rho, do_not_allocate_becsum )
   IMPLICIT NONE
   TYPE (scf_type) :: rho
   LOGICAL,INTENT(IN),OPTIONAL :: do_not_allocate_becsum ! PAW hack
   LOGICAL                     :: allocate_becsum        ! PAW hack
   allocate ( rho%of_r( dfftp%nnr, nspin) )
   allocate ( rho%of_g( ngm, nspin ) )
   if (dft_is_meta() .or. lxdm) then
      allocate ( rho%kin_r( dfftp%nnr, nspin) )
      allocate ( rho%kin_g( ngm, nspin ) )
   else
      allocate ( rho%kin_r(1,1) )
      allocate ( rho%kin_g(1,1) )
   endif

   lda_plus_u_co = lda_plus_u .and. .not. (nspin == 4 )
   lda_plus_u_nc = lda_plus_u .and.       (nspin == 4 )
   if (lda_plus_u_co) allocate (rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
   if (lda_plus_u_nc) allocate (rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))

   if (okpaw) then ! See the top of the file for clarification
      if(present(do_not_allocate_becsum)) then
         allocate_becsum = .not. do_not_allocate_becsum
      else
         allocate_becsum = .true.
      endif
      if(allocate_becsum) allocate (rho%bec(nhm*(nhm+1)/2,nat,nspin))
   endif
   
 return
 END SUBROUTINE create_scf_type

 SUBROUTINE destroy_scf_type ( rho )
   IMPLICIT NONE
   TYPE (scf_type) :: rho

   if (ALLOCATED(rho%of_r))  deallocate(rho%of_r)
   if (ALLOCATED(rho%of_g))  deallocate(rho%of_g)
   if (ALLOCATED(rho%kin_r)) deallocate(rho%kin_r)
   if (ALLOCATED(rho%kin_g)) deallocate(rho%kin_g)
   if (ALLOCATED(rho%ns))    deallocate(rho%ns)
   if (ALLOCATED(rho%ns_nc))    deallocate(rho%ns_nc)
   if (ALLOCATED(rho%bec))   deallocate(rho%bec)

   return
 END SUBROUTINE destroy_scf_type
 !
 
 SUBROUTINE create_mix_type ( rho )
   IMPLICIT NONE
   TYPE (mix_type) :: rho
   allocate ( rho%of_g( ngms, nspin ) )
   rho%of_g = 0._dp
   if (dft_is_meta() .or. lxdm) then
      allocate (rho%kin_g( ngms, nspin ) )
      rho%kin_g = 0._dp
   end if
   lda_plus_u_co = lda_plus_u .and. .not. (nspin == 4 )
   lda_plus_u_nc = lda_plus_u .and.       (nspin == 4 )
   if (lda_plus_u_nc) then
      allocate (rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
      rho%ns_nc    = 0._dp
   endif
   if (lda_plus_u_co) then
      allocate (rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
      rho%ns       = 0._dp
   endif
   if (okpaw) then
      allocate (rho%bec(nhm*(nhm+1)/2,nat,nspin))
      rho%bec   = 0._dp
   end if
   rho%el_dipole =  0._dp
   
 return
 END SUBROUTINE create_mix_type

 SUBROUTINE destroy_mix_type ( rho )
   IMPLICIT NONE
   TYPE (mix_type) :: rho

   if (ALLOCATED(rho%of_g))  deallocate(rho%of_g)
   if (ALLOCATED(rho%kin_g)) deallocate(rho%kin_g)
   if (ALLOCATED(rho%ns))    deallocate(rho%ns)
   if (ALLOCATED(rho%ns_nc))    deallocate(rho%ns_nc)
   if (ALLOCATED(rho%bec))   deallocate(rho%bec)

   return
 END SUBROUTINE destroy_mix_type
 !
 subroutine assign_scf_to_mix_type(rho_s, rho_m)
   IMPLICIT NONE
   TYPE (scf_type), INTENT(IN)  :: rho_s
   TYPE (mix_type), INTENT(INOUT) :: rho_m
   REAL(DP) :: e_dipole
      
   rho_m%of_g(1:ngms,:) = rho_s%of_g(1:ngms,:)
   
   if (dft_is_meta() .or. lxdm) rho_m%kin_g(1:ngms,:) = rho_s%kin_g(1:ngms,:)
   if (lda_plus_u_nc) rho_m%ns_nc  = rho_s%ns_nc
   if (lda_plus_u_co) rho_m%ns     = rho_s%ns
   if (okpaw)         rho_m%bec = rho_s%bec
   
   if (dipfield) then
      CALL compute_el_dip(emaxpos, eopreg, edir, rho_s%of_r,e_dipole)
      rho_m%el_dipole = e_dipole
   endif
   
 return
 end subroutine assign_scf_to_mix_type
 !
 subroutine assign_mix_to_scf_type(rho_m, rho_s)
   USE wavefunctions_module, ONLY : psic
   USE control_flags,        ONLY : gamma_only
   USE gvect,                ONLY : nl, nlm
   IMPLICIT NONE
   TYPE (mix_type), INTENT(IN) :: rho_m
   TYPE (scf_type), INTENT(INOUT) :: rho_s
   INTEGER :: is
   
   rho_s%of_g(1:ngms,:) = rho_m%of_g(1:ngms,:)
   ! define rho_s%of_r 

   DO is = 1, nspin
      psic(:) = ( 0.D0, 0.D0 )
      psic(nl(:)) = rho_s%of_g(:,is)
      IF ( gamma_only ) psic(nlm(:)) = CONJG( rho_s%of_g(:,is) )
      CALL invfft ('Dense', psic, dfftp)
      rho_s%of_r(:,is) = psic(:)
   END DO

   if (dft_is_meta() .or. lxdm) then
      rho_s%kin_g(1:ngms,:) = rho_m%kin_g(:,:)
      ! define rho_s%kin_r 
      DO is = 1, nspin
         psic(:) = ( 0.D0, 0.D0 )
         psic(nl(:)) = rho_s%kin_g(:,is)
         IF ( gamma_only ) psic(nlm(:)) = CONJG( rho_s%kin_g(:,is) )
         CALL invfft ('Dense', psic, dfftp)
         rho_s%kin_r(:,is) = psic(:)
      END DO
   end if

   if (lda_plus_u_nc) rho_s%ns_nc(:,:,:,:) = rho_m%ns_nc(:,:,:,:)
   if (lda_plus_u_co) rho_s%ns(:,:,:,:)    = rho_m%ns(:,:,:,:)
   if (okpaw)      rho_s%bec(:,:,:)  = rho_m%bec(:,:,:)
       
   return
 end subroutine assign_mix_to_scf_type
 !
 !----------------------------------------------------------------------------
 subroutine scf_type_COPY (X,Y)
  !----------------------------------------------------------------------------
  ! works like DCOPY for scf_type copy variables :  Y = X 
  USE kinds, ONLY : DP
  IMPLICIT NONE
  TYPE(scf_type), INTENT(IN)    :: X
  TYPE(scf_type), INTENT(INOUT) :: Y
  Y%of_r  = X%of_r
  Y%of_g  = X%of_g
  if (dft_is_meta() .or. lxdm) then
     Y%kin_r = X%kin_r
     Y%kin_g = X%kin_g
  end if
  if (lda_plus_u_nc) Y%ns_nc = X%ns_nc
  if (lda_plus_u_co) Y%ns    = X%ns
  if (okpaw)      Y%bec = X%bec
  !
  RETURN
 end subroutine scf_type_COPY
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_AXPY (A,X,Y)
  !----------------------------------------------------------------------------
  ! works like daxpy for scf_type variables :  Y = A * X + Y
  ! NB: A is a REAL(DP) number
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP)                      :: A
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  Y%of_g  = Y%of_g  + A * X%of_g
  if (dft_is_meta() .or. lxdm) Y%kin_g = Y%kin_g + A * X%kin_g
  if (lda_plus_u_nc) Y%ns_nc = Y%ns_nc + A * X%ns_nc
  if (lda_plus_u_co) Y%ns = Y%ns + A * X%ns
  if (okpaw)     Y%bec = Y%bec + A * X%bec
  if (dipfield)  Y%el_dipole =  Y%el_dipole + A * X%el_dipole
  !
  RETURN
 END SUBROUTINE mix_type_AXPY
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_COPY (X,Y)
  !----------------------------------------------------------------------------
  ! works like DCOPY for mix_type copy variables :  Y = X 
  USE kinds, ONLY : DP
  IMPLICIT NONE
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  Y%of_g  = X%of_g
  if (dft_is_meta() .or. lxdm) Y%kin_g = X%kin_g
  if (lda_plus_u_nc) Y%ns_nc  = X%ns_nc
  if (lda_plus_u_co) Y%ns  = X%ns
  if (okpaw)      Y%bec = X%bec
  if (dipfield)   Y%el_dipole =  X%el_dipole
  !
  RETURN
 end subroutine mix_type_COPY
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_SCAL (A,X)
  !----------------------------------------------------------------------------
  ! works like DSCAL for mix_type copy variables :  X = A * X 
  ! NB: A is a REAL(DP) number
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP),       INTENT(IN)    :: A
  TYPE(mix_type), INTENT(INOUT) :: X
  X%of_g(:,:)  = A * X%of_g(:,:)
  if (dft_is_meta() .or. lxdm) X%kin_g = A * X%kin_g
  if (lda_plus_u_nc) X%ns_nc = A * X%ns_nc
  if (lda_plus_u_co) X%ns    = A * X%ns
  if (okpaw)      X%bec= A * X%bec
  if (dipfield)   X%el_dipole =  A * X%el_dipole
  !
  RETURN
 end subroutine mix_type_SCAL
 !
 subroutine high_frequency_mixing ( rhoin, input_rhout, alphamix )
   USE wavefunctions_module, ONLY : psic
   USE control_flags,        ONLY : gamma_only
   USE gvect,                ONLY : nl, nlm
 IMPLICIT NONE
   TYPE (scf_type), INTENT(INOUT)     :: rhoin
   TYPE (scf_type), INTENT(IN)  :: input_rhout
   REAL(DP), INTENT(IN) :: alphamix
   INTEGER :: is
   if (ngms < ngm ) then
      rhoin%of_g = rhoin%of_g + alphamix * ( input_rhout%of_g-rhoin%of_g)
      rhoin%of_g(1:ngms,1:nspin) = (0.d0,0.d0)
      ! define rho_s%of_r 
      DO is = 1, nspin
         psic(:) = ( 0.D0, 0.D0 )
         psic(nl(:)) = rhoin%of_g(:,is)
         IF ( gamma_only ) psic(nlm(:)) = CONJG( rhoin%of_g(:,is) )
         CALL invfft ('Dense', psic, dfftp)
         rhoin%of_r(:,is) = psic(:)
      END DO
      !
      if (dft_is_meta() .or. lxdm) then
         rhoin%kin_g = rhoin%kin_g + alphamix * ( input_rhout%kin_g-rhoin%kin_g)
         rhoin%kin_g(1:ngms,1:nspin) = (0.d0,0.d0)
         ! define rho_s%of_r 
         DO is = 1, nspin
            psic(:) = ( 0.D0, 0.D0 )
            psic(nl(:)) = rhoin%kin_g(:,is)
            IF ( gamma_only ) psic(nlm(:)) = CONJG( rhoin%kin_g(:,is) )
            CALL invfft ('Dense', psic, dfftp)
            rhoin%kin_r(:,is) = psic(:)
         END DO
      end if
   else
      rhoin%of_g(:,:)= (0.d0,0.d0)
      rhoin%of_r(:,:)= 0.d0
      if (dft_is_meta() .or. lxdm) then
         rhoin%kin_g(:,:)= (0.d0,0.d0)
         rhoin%kin_r(:,:)= 0.d0
      endif
   endif
   if (lda_plus_u_nc) rhoin%ns_nc(:,:,:,:) = 0.d0
   if (lda_plus_u_co) rhoin%ns(:,:,:,:)    = 0.d0

   return
 end subroutine high_frequency_mixing 


 subroutine open_mix_file( iunit, extension, exst )
   USE control_flags,        ONLY : io_level
   implicit none
   character(len=*), intent(in) :: extension
   integer, intent(in) :: iunit
   logical :: exst
   ! define lengths (in real numbers) of different record chunks
   rlen_rho = 2 * ngms * nspin
   if (dft_is_meta() .or. lxdm) rlen_kin =  2 * ngms * nspin
   if (lda_plus_u_co)  rlen_ldaU = (2*Hubbard_lmax+1)**2 *nspin*nat
   if (lda_plus_u_nc)  rlen_ldaU = 2 * (2*Hubbard_lmax+1)**2 *nspin*nat
   if (okpaw)          rlen_bec = (nhm*(nhm+1)/2) * nat * nspin
   if (dipfield)       rlen_dip = 1
   ! define the starting point of the different chunks. Beware: each starting point
   ! is the index of a COMPLEX array. When real arrays with odd dimension are copied
   ! to/from the complex array io_buffer, the last complex number will be half-filled
   ! but must still be counted as one!
   start_rho = 1
   start_kin = start_rho + rlen_rho / 2
   start_ldaU = start_kin + rlen_kin / 2
   start_bec = start_ldaU + ( rlen_ldaU + 1 ) / 2 
   start_dipole = start_bec + ( rlen_bec + 1 ) / 2
   ! define total record length, in complex numbers
   record_length = start_dipole + rlen_dip - 1
   ! open file and allocate io_buffer
   call open_buffer ( iunit, extension, record_length, io_level, exst)
   allocate (io_buffer(record_length))
   ! setting to zero -prevents trouble with "holes" due to odd dimensions of real arrays
   io_buffer (:) = (0.0_dp, 0.0_dp)
   !
 return
 end subroutine open_mix_file
 !
 subroutine close_mix_file( iunit, stat )
   implicit none
   integer, intent(in) :: iunit
   character(len=*), intent(in) :: stat
   deallocate (io_buffer)
   call close_buffer ( iunit, trim(stat) ) 
   return
 end subroutine close_mix_file

 subroutine davcio_mix_type( rho, iunit, record, iflag )
   implicit none
   type (mix_type) :: rho
   integer, intent(in) :: iunit, record, iflag

   if (iflag > 0) then      
      call DCOPY(rlen_rho,rho%of_g,1,io_buffer(start_rho),1)
      if (dft_is_meta() .or. lxdm)  call DCOPY(rlen_kin, rho%kin_g,1,io_buffer(start_kin),1)
      if (lda_plus_u_nc)  call DCOPY(rlen_ldaU,rho%ns_nc, 1,io_buffer(start_ldaU),1)
      if (lda_plus_u_co)  call DCOPY(rlen_ldaU,rho%ns,   1,io_buffer(start_ldaU),1)
      if (okpaw)          call DCOPY(rlen_bec, rho%bec,  1,io_buffer(start_bec),1)
      if (dipfield)       io_buffer(start_dipole) = CMPLX ( rho%el_dipole, 0.0_dp )
      CALL save_buffer( io_buffer, record_length, iunit, record )   
   else if (iflag < 0 ) then 
      CALL get_buffer( io_buffer, record_length, iunit, record )
      call DCOPY(rlen_rho,io_buffer(start_rho),1,rho%of_g,1)
      if (dft_is_meta() .or. lxdm) call DCOPY(rlen_kin,io_buffer(start_kin), 1,rho%kin_g,1)
      if (lda_plus_u_co) call DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns,1)
      if (lda_plus_u_nc) call DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns_nc,1)
      if (okpaw)         call DCOPY(rlen_bec, io_buffer(start_bec), 1,rho%bec,1)
      if (dipfield)      rho%el_dipole = REAL ( io_buffer(start_dipole) )
   end if
 end subroutine davcio_mix_type
 !
 !----------------------------------------------------------------------------
 FUNCTION rho_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE spin_orb,      ONLY : domag
  USE control_flags, ONLY : gamma_only
  USE paw_onecenter, ONLY : paw_ddot
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  type(mix_type), INTENT(IN) :: rho1, rho2
  INTEGER,        INTENT(IN) :: gf
  REAL(DP)                :: rho_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  fac = e2 * fpi / tpiba2
  !
  rho_ddot = 0.D0
  
  IF ( nspin == 1 ) THEN
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1), DP ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... first the charge
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1)+rho1%of_g(ig,2) ) * &
                              ( rho2%of_g(ig,1)+rho2%of_g(ig,2) ), DP ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     ! ... then the magnetization
     !
     fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + &
                   fac * REAL( CONJG( rho1%of_g(1,1) - rho1%of_g(1,2) ) * &
                                    ( rho2%of_g(1,1) - rho2%of_g(1,2) ), DP )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   fac * REAL( CONJG( rho1%of_g(ig,1) - rho1%of_g(ig,2) ) * &
                                    ( rho2%of_g(ig,1) - rho2%of_g(ig,2) ), DP )
        !
     END DO
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1), DP ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     IF (domag) THEN
        fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
        !
        IF ( gstart == 2 ) THEN
           !
           rho_ddot = rho_ddot + &
                   fac * ( REAL( CONJG( rho1%of_g(1,2))*(rho2%of_g(1,2) ),DP ) + &
                           REAL( CONJG( rho1%of_g(1,3))*(rho2%of_g(1,3) ),DP ) + &
                           REAL( CONJG( rho1%of_g(1,4))*(rho2%of_g(1,4) ),DP ) )
           !
        END IF
        !
        IF ( gamma_only ) fac = 2.D0 * fac
        !
        DO ig = gstart, gf
           !
           rho_ddot = rho_ddot + &
                   fac *( REAL( CONJG( rho1%of_g(ig,2))*(rho2%of_g(ig,2) ), DP ) + &
                          REAL( CONJG( rho1%of_g(ig,3))*(rho2%of_g(ig,3) ), DP ) + &
                          REAL( CONJG( rho1%of_g(ig,4))*(rho2%of_g(ig,4) ), DP ) )
           !
        END DO
        !
     END IF
     !
  END IF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL mp_sum(  rho_ddot , intra_bgrp_comm )
  !
  IF (dft_is_meta()) rho_ddot = rho_ddot + tauk_ddot( rho1, rho2, gf )
  IF (lda_plus_u )   rho_ddot = rho_ddot + ns_ddot(rho1,rho2)
  ! 
  ! Beware: paw_ddot has a hidden parallelization on all processors
  !         it must be called on all processors or else it will hang
  ! Beware: commented out because it yields too often negative values
  ! IF (okpaw)         rho_ddot = rho_ddot + paw_ddot(rho1%bec, rho2%bec)
  IF (dipfield)      rho_ddot = rho_ddot + (e2/2.0_DP)* &
                                    (rho1%el_dipole * rho2%el_dipole)*omega/fpi

  RETURN
  !
 END FUNCTION rho_ddot
 !
!----------------------------------------------------------------------------
FUNCTION tauk_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
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
  type(mix_type), INTENT(IN) :: rho1, rho2
  INTEGER,     INTENT(IN) :: gf
  REAL(DP)                :: tauk_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  tauk_ddot = 0.D0
  !
!  write (*,*) rho1%kin_g(1:4,1)
!  if (.true. ) stop
  IF ( nspin == 1 ) THEN
     !
     DO ig = gstart, gf
        tauk_ddot = tauk_ddot + &
                    REAL( CONJG( rho1%kin_g(ig,1) )*rho2%kin_g(ig,1) ) 
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + &
                    REAL( CONJG( rho1%kin_g(1,1) ) * rho2%kin_g(1,1) )
        !
     END IF
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     DO ig = gstart, gf
        !
        tauk_ddot = tauk_ddot + &
                         ( REAL( CONJG(rho1%kin_g(ig,1))*rho2%kin_g(ig,1) ) + &
                           REAL( CONJG(rho1%kin_g(ig,2))*rho2%kin_g(ig,2) ) )
        !
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + &
                         ( REAL( CONJG( rho1%kin_g(1,1))*rho2%kin_g(1,1) ) + &
                           REAL( CONJG( rho1%kin_g(1,2))*rho2%kin_g(1,2) ) )
        !
     END IF
     tauk_ddot = 0.5D0 *  tauk_ddot 
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     DO ig = gstart, gf
        !
        tauk_ddot = tauk_ddot + &
                         ( REAL( CONJG(rho1%kin_g(ig,1))*rho2%kin_g(ig,1) ) + &
                           REAL( CONJG(rho1%kin_g(ig,2))*rho2%kin_g(ig,2) ) + &
                           REAL( CONJG(rho1%kin_g(ig,3))*rho2%kin_g(ig,3) ) + &
                           REAL( CONJG(rho1%kin_g(ig,4))*rho2%kin_g(ig,4) ) )
        !
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + &
                         ( REAL( CONJG( rho1%kin_g(1,1))*rho2%kin_g(1,1) ) + &
                           REAL( CONJG( rho1%kin_g(1,2))*rho2%kin_g(1,2) ) + &
                           REAL( CONJG( rho1%kin_g(1,3))*rho2%kin_g(1,3) ) + &
                           REAL( CONJG( rho1%kin_g(1,4))*rho2%kin_g(1,4) ) )
        !
     END IF
     !
  END IF
  !
  fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
  !
  tauk_ddot = fac * tauk_ddot * omega * 0.5D0
  !
  CALL mp_sum(  tauk_ddot , intra_bgrp_comm )
  !
  RETURN
  !
END FUNCTION tauk_ddot
!----------------------------------------------------------------------------
FUNCTION ns_ddot( rho1, rho2 )
  !----------------------------------------------------------------------------
  !
  ! ... calculates U/2 \sum_i ns1(i)*ns2(i)
  ! ... used as an estimate of the self-consistency error on the 
  ! ... LDA+U correction to the energy
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_l, Hubbard_U, Hubbard_alpha
  USE ions_base, ONLY : nat, ityp
  !
  IMPLICIT NONE  
  !
  type(mix_type), INTENT(IN) :: rho1, rho2
  REAL(DP)             :: ns_ddot
  !
  INTEGER :: na, nt, m1, m2
  !
  ns_ddot = 0.D0
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
        m1 = 2 * Hubbard_l(nt) + 1
        m2 = 2 * Hubbard_l(nt) + 1

        if (nspin.eq.4) then
          ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( CONJG(rho1%ns_nc(:m1,:m2,:nspin,na))*rho2%ns_nc(:m1,:m2,:nspin,na) )
        else
          ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( rho1%ns(:m1,:m2,:nspin,na)*rho2%ns(:m1,:m2,:nspin,na) )
        endif

     END IF
  END DO
  !
  IF ( nspin == 1 ) ns_ddot = 2.D0*ns_ddot
  !
  RETURN
  !
END FUNCTION ns_ddot
 !----------------------------------------------------------------------------
 FUNCTION local_tf_ddot( rho1, rho2, ngm0 )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE control_flags, ONLY : gamma_only
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: ngm0
  COMPLEX(DP), INTENT(IN) :: rho1(ngm0), rho2(ngm0)
  REAL(DP)                :: local_tf_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  local_tf_ddot = 0.D0
  !
  fac = e2 * fpi / tpiba2
  !
  DO ig = gstart, ngm0
     local_tf_ddot = local_tf_ddot + REAL( CONJG(rho1(ig))*rho2(ig) ) / gg(ig)
  END DO
  !
  local_tf_ddot = fac * local_tf_ddot * omega * 0.5D0
  !
  IF ( gamma_only ) local_tf_ddot = 2.D0 * local_tf_ddot
  !
  CALL mp_sum(  local_tf_ddot , intra_bgrp_comm )
  !
  RETURN
  !
 END FUNCTION local_tf_ddot
 !
 SUBROUTINE bcast_scf_type ( rho, root, comm )
  !----------------------------------------------------------------------------
  ! ... Broadcast all mixed quantities from first pool to all others
  ! ... Needed to prevent divergencies in k-point parallization
  !
  USE mp,            ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  type(scf_type), INTENT(INOUT) :: rho
  INTEGER, INTENT(IN) :: root, comm
  !
  CALL mp_bcast ( rho%of_g, root, comm )
  CALL mp_bcast ( rho%of_r, root, comm )
  IF ( dft_is_meta() .or. lxdm) THEN
     CALL mp_bcast ( rho%kin_g, root, comm )
     CALL mp_bcast ( rho%kin_r, root, comm )
  END IF
  IF ( lda_plus_u_co) CALL mp_bcast ( rho%ns,    root, comm )
  IF ( lda_plus_u_nc) CALL mp_bcast ( rho%ns_nc, root, comm )
  IF ( okpaw )        CALL mp_bcast ( rho%bec,   root, comm )
  !
  END SUBROUTINE
 !
END MODULE scf
