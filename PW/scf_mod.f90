!
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
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
  USE lsda_mod,   ONLY : nspin
  USE ldaU,       ONLY : lda_plus_u, Hubbard_lmax
  USE ions_base,  ONLY : nat
  USE funct,      ONLY : dft_is_meta
  USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm
  USE gsmooth,    ONLY : ngms
  !
  SAVE
  !
  TYPE scf_type
     REAL(DP),    ALLOCATABLE :: of_r(:,:)  ! the charge density in R-space
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     REAL(DP),    ALLOCATABLE :: kin_r(:,:) ! the kinetic energy density in R-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:) ! the kinetic energy density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix 
  END TYPE scf_type
  !
  TYPE mix_type
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:) ! the charge density in G-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:) ! the charge density in G-space
     REAL(DP),   ALLOCATABLE  :: ns(:,:,:,:)! the LDA+U occupation matrix 
  END TYPE mix_type

  type (scf_type) :: rho

  REAL(DP) :: v_of_0    ! vltot(G=0)      
  REAL(DP), ALLOCATABLE :: &
       vr(:,:),        &! the Hartree + xc potential in real space
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grig)
       rho_core(:),    &! the core charge in real space
!       tauk(:,:),      &! kinetic energy density in real space (dense grid)
       kedtau(:,:),    &! position dependent kinetic energy enhancement factor
                        ! used in META-GGA in real space (smooth grid)
       kedtaur(:,:)     ! position dependent kinetic energy enhancement factor
                        ! used in META-GGA in real space (dense grid)
  COMPLEX(DP), ALLOCATABLE :: &
       rhog_core(:)     ! the core charge in reciprocal space

  INTEGER, PRIVATE  :: record_length, &
                       rlen_rho=0, rlen_kin=0, rlen_ldaU=0, &
                       start_rho=0, start_kin=0, start_ldaU=0
  REAL(DP), PRIVATE, ALLOCATABLE:: io_buffer(:)
CONTAINS

 SUBROUTINE create_scf_type ( rho )
 IMPLICIT NONE
 TYPE (scf_type) :: rho
 allocate ( rho%of_r( nrxx, nspin) )
 allocate ( rho%of_g( ngm, nspin ) )
 if (dft_is_meta()) then
    allocate ( rho%kin_r( nrxx, nspin) )
    allocate ( rho%kin_g( ngm, nspin ) )
 endif
 if (lda_plus_u) allocate (rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
 return
 END SUBROUTINE create_scf_type

 SUBROUTINE destroy_scf_type ( rho )
 IMPLICIT NONE
 TYPE (scf_type) :: rho
 if (allocated(rho%of_r)) deallocate(rho%of_r)
 if (allocated(rho%of_g)) deallocate(rho%of_g)
 if (allocated(rho%kin_r)) deallocate(rho%kin_r)
 if (allocated(rho%kin_g)) deallocate(rho%kin_g)
 if (allocated(rho%ns)) deallocate(rho%ns)
 return
 END SUBROUTINE destroy_scf_type
 !
 SUBROUTINE create_mix_type ( rho )
 TYPE (mix_type) :: rho
 allocate ( rho%of_g( ngms, nspin ) )
 if (dft_is_meta()) allocate ( rho%kin_g( ngms, nspin ) )
 if (lda_plus_u) allocate (rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
 return
 END SUBROUTINE create_mix_type

 SUBROUTINE destroy_mix_type ( rho )
 IMPLICIT NONE
 TYPE (mix_type) :: rho
 if (allocated(rho%of_g)) deallocate(rho%of_g)
 if (allocated(rho%kin_g)) deallocate(rho%kin_g)
 if (allocated(rho%ns)) deallocate(rho%ns)
 return
 END SUBROUTINE destroy_mix_type
 !
 subroutine assign_scf_to_mix_type(rho_s, rho_m)
 IMPLICIT NONE
 TYPE (scf_type), INTENT(IN)  :: rho_s
 TYPE (mix_type), INTENT(INOUT) :: rho_m
 rho_m%of_g(:,:) = rho_s%of_g(1:ngms,:)
 if (dft_is_meta()) rho_m%kin_g(:,:) = rho_s%kin_g(1:ngms,:)
 if (lda_plus_u) rho_m%ns(:,:,:,:) = rho_s%ns(:,:,:,:)
 return
 end subroutine assign_scf_to_mix_type
 !
 subroutine assign_mix_to_scf_type(rho_m, rho_s)
   USE wavefunctions_module, ONLY : psic
   USE wvfct,                ONLY : gamma_only
   USE gvect,                ONLY : nl, nlm
   IMPLICIT NONE
   TYPE (mix_type), INTENT(IN) :: rho_m
   TYPE (scf_type), INTENT(INOUT) :: rho_s
   INTEGER :: is
   rho_s%of_g(1:ngms,:) = rho_m%of_g(:,:)
   ! define rho_s%of_r 
   DO is = 1, nspin
      psic(:) = ( 0.D0, 0.D0 )
      psic(nl(:)) = rho_s%of_g(:,is)
      IF ( gamma_only ) psic(nlm(:)) = CONJG( rho_s%of_g(:,is) )
      CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
      rho_s%of_r(:,is) = psic(:)
   END DO
   if (dft_is_meta()) then
      rho_s%kin_g(1:ngms,:) = rho_m%kin_g(:,:)
      ! define rho_s%kin_r 
      DO is = 1, nspin
         psic(:) = ( 0.D0, 0.D0 )
         psic(nl(:)) = rho_s%kin_g(:,is)
         IF ( gamma_only ) psic(nlm(:)) = CONJG( rho_s%kin_g(:,is) )
         CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
         rho_s%kin_r(:,is) = psic(:)
      END DO
   end if
   if (lda_plus_u) rho_s%ns(:,:,:,:) = rho_m%ns(:,:,:,:)
   return
 end subroutine assign_mix_to_scf_type
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_AXPY (A,X,Y)
  !----------------------------------------------------------------------------
  ! works like DAXPY for scf_type variables :  Y = A * X + Y
  ! NB: A is a REAL(DP) number
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP)                      :: A
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  Y%of_g  = Y%of_g  + A * X%of_g
  if (dft_is_meta()) Y%kin_g = Y%kin_g + A * X%kin_g
  if (lda_plus_u) Y%ns = Y%ns + A * X%ns
  !
  RETURN
 END SUBROUTINE mix_type_AXPY
 !
 subroutine high_frequency_mixing ( rhoin, input_rhout, alphamix )
   USE wavefunctions_module, ONLY : psic
   USE wvfct,                ONLY : gamma_only
   USE gvect,                ONLY : nl, nlm
 IMPLICIT NONE
   TYPE (scf_type), INTENT(INOUT)     :: rhoin
   TYPE (scf_type), INTENT(IN)  :: input_rhout
   REAL(DP), INTENT(IN) :: alphamix
   INTEGER :: is
   if (ngms < ngm ) then
      rhoin%of_g = rhoin%of_g + alphamix * ( input_rhout%of_g-rhoin%of_g)
      rhoin%of_g(1:ngms,nspin) = (0.d0,0.d0)
      ! define rho_s%of_r 
      DO is = 1, nspin
         psic(:) = ( 0.D0, 0.D0 )
         psic(nl(:)) = rhoin%of_g(:,is)
         IF ( gamma_only ) psic(nlm(:)) = CONJG( rhoin%of_g(:,is) )
         CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
         rhoin%of_r(:,is) = psic(:)
      END DO
      !
      if (dft_is_meta()) then
         rhoin%kin_g = rhoin%kin_g + alphamix * ( input_rhout%kin_g-rhoin%kin_g)
         rhoin%kin_g(1:ngms,nspin) = (0.d0,0.d0)
         ! define rho_s%of_r 
         DO is = 1, nspin
            psic(:) = ( 0.D0, 0.D0 )
            psic(nl(:)) = rhoin%kin_g(:,is)
            IF ( gamma_only ) psic(nlm(:)) = CONJG( rhoin%kin_g(:,is) )
            CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
            rhoin%kin_r(:,is) = psic(:)
         END DO
      end if
   else
      rhoin%of_g(:,:)= (0.d0,0.d0)
      rhoin%of_r(:,:)= 0.d0
      if (dft_is_meta()) then
         rhoin%kin_g(:,:)= (0.d0,0.d0)
         rhoin%kin_r(:,:)= 0.d0
      endif
   endif
   if (lda_plus_u) rhoin%ns(:,:,:,:) = 0.d0
   return
 end subroutine high_frequency_mixing 

 subroutine diropn_mix_file( iunit, extension, exst )
 implicit none
 character(len=*), intent(in) :: extension
 integer, intent(in) :: iunit
 logical :: exst
 ! define lengths of different record chunks
 rlen_rho = 2 * ngms * nspin
 if (dft_is_meta() ) rlen_kin = 2 * ngms * nspin
 if (lda_plus_u) rlen_ldaU = (2*Hubbard_lmax+1)**2*nspin*nat
 ! define total record length
 record_length = rlen_rho + rlen_kin + rlen_ldaU
 ! and the starting point of different chunks
 start_rho = 1
 start_kin = start_rho + rlen_rho
 start_ldaU = start_kin + rlen_kin
 ! open file and allocate io_buffer
 call diropn ( iunit, extension, record_length, exst)
 allocate (io_buffer(record_length))
 !
 return
 end subroutine diropn_mix_file
 !
 subroutine close_mix_file( iunit )
 implicit none
 integer, intent(in) :: iunit
 deallocate (io_buffer)
 close(iunit,status='keep')
 return
 end subroutine close_mix_file

 subroutine davcio_mix_type( rho, iunit, record, iflag )
 implicit none
 type (mix_type) :: rho
 integer, intent(in) :: iunit, record, iflag
 if (iflag > 0) then
    call DCOPY(rlen_rho,rho%of_g,1,io_buffer(start_rho),1)
    if (dft_is_meta()) call DCOPY(rlen_kin,rho%kin_g,1,io_buffer(start_kin),1)
    if (lda_plus_u) call DCOPY(rlen_ldaU,rho%ns,1,io_buffer(start_ldaU),1)
 end if
 CALL davcio( io_buffer, record_length, iunit, record, iflag )
 if (iflag < 0) then
    call DCOPY(rlen_rho,io_buffer(start_rho),1,rho%of_g,1)
    if (dft_is_meta()) call DCOPY(start_kin,io_buffer(start_kin),1,rho%kin_g,1)
    if (lda_plus_u) call DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns,1)
 end if
 return
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
  USE wvfct,         ONLY : gamma_only
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
  !
  IF ( nspin == 1 ) THEN
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1) ) / gg(ig)
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
                              ( rho2%of_g(ig,1)+rho2%of_g(ig,2) ) ) / gg(ig)
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
                                    ( rho2%of_g(1,1) - rho2%of_g(1,2) ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   fac * REAL( CONJG( rho1%of_g(ig,1) - rho1%of_g(ig,2) ) * &
                                    ( rho2%of_g(ig,1) - rho2%of_g(ig,2) ) )
        !
     END DO
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1) ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + &
                   fac * ( REAL( CONJG( rho1%of_g(1,2))*(rho2%of_g(1,2) ) ) + &
                           REAL( CONJG( rho1%of_g(1,3))*(rho2%of_g(1,3) ) ) + &
                           REAL( CONJG( rho1%of_g(1,4))*(rho2%of_g(1,4) ) ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   fac *( REAL( CONJG( rho1%of_g(ig,2))*(rho2%of_g(ig,2) ) ) + &
                          REAL( CONJG( rho1%of_g(ig,3))*(rho2%of_g(ig,3) ) ) + &
                          REAL( CONJG( rho1%of_g(ig,4))*(rho2%of_g(ig,4) ) ) )
        !
     END DO
     !
  END IF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL reduce( 1, rho_ddot )
  !
  if (dft_is_meta()) rho_ddot = rho_ddot + tauk_ddot( rho1, rho2, gf )
  IF (lda_plus_u ) rho_ddot = rho_ddot + ns_ddot(rho1,rho2)

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
  USE wvfct,         ONLY : gamma_only
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
  CALL reduce( 1, tauk_ddot )
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
        ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( rho1%ns(:m1,:m2,:nspin,na)*rho2%ns(:m1,:m2,:nspin,na) )
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
  USE wvfct,         ONLY : gamma_only
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
  CALL reduce( 1, local_tf_ddot )
  !
  RETURN
  !
 END FUNCTION local_tf_ddot
 !
END MODULE scf
