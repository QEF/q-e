!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE bp
  !
  ! ... The variables needed for the Berry phase polarization calculation
  !
  USE kinds, ONLY: DP
  USE becmod, ONLY : bec_type
  !
  SAVE
  PRIVATE
  PUBLIC:: lberry, lelfield, lorbm, gdir, nppstr, nberrycyc, evcel, evcelp, evcelm, &
           fact_hepsi, bec_evcel, mapgp_global, mapgm_global, nppstr_3d, &
           ion_pol, el_pol, fc_pol, l_el_pol_old, el_pol_old, el_pol_acc, &
           nx_el, l3dstring, efield, efield_cart, efield_cry, transform_el,&
           mapg_owner, phase_control
  PUBLIC :: allocate_bp_efield, deallocate_bp_efield, bp_global_map
  PUBLIC :: pdl_tot
  !
  LOGICAL :: &
       lberry  =.false., & ! if .TRUE. calculate polarization using Berry phase
       lelfield=.false., & ! if .TRUE. finite electric field using Berry phase
       lorbm=.false.       ! if .TRUE. calculate orbital magnetization (Kubo terms)
  INTEGER :: &
       gdir,        &! G-vector for polarization calculation
       nppstr,      &! number of k-points (parallel vector)
       nberrycyc     ! number of cycles for convergence in electric field 
                     ! without changing the selfconsistent charge
  REAL(DP) :: efield ! electric field intensity in a.u.
  COMPLEX(DP), ALLOCATABLE , TARGET :: evcel(:,:) 
                     ! wavefunctions for calculating the electric field operator
  COMPLEX(DP), ALLOCATABLE , TARGET :: evcelm(:,:,:) 
                     ! wavefunctions for  storing projectors for  electric field operator
  COMPLEX(DP), ALLOCATABLE , TARGET :: evcelp(:,:,:) 
                     ! wavefunctions for  storing projectors for  electric field operator
  COMPLEX(DP), ALLOCATABLE, TARGET :: fact_hepsi(:,:)
                     ! factors for hermitean electric field operators
  !COMPLEX(DP), ALLOCATABLE, TARGET :: bec_evcel(:,:) 
  !                   !for storing bec's factors with evcel
  TYPE(bec_type) :: bec_evcel
  INTEGER, ALLOCATABLE, TARGET :: mapgp_global(:,:)
                     ! map for G'= G+1 correspondence
  INTEGER, ALLOCATABLE, TARGET :: mapgm_global(:,:)
                     ! map for G'= G-1 correspondence
  REAL(DP) :: ion_pol(3) ! the ionic polarization
  REAL(DP) :: el_pol(3)  ! the electronic polarization
  REAL(DP) :: fc_pol(3)  ! the prefactor for the electronic polarization
  LOGICAL  :: l_el_pol_old! if true there is already stored a n older value for the polarization
                          ! neeeded for having correct polarization during MD
  REAL(DP) :: el_pol_old(3)! the old  electronic polarization
  REAL(DP) :: el_pol_acc(3)! accumulator for the electronic polarization

  INTEGER :: nppstr_3d(3)  ! number of element of strings along the reciprocal directions
  INTEGER, ALLOCATABLE :: nx_el(:,:) ! index for string to k-point map, (nks*nspin,dir=3)
  LOGICAL :: l3dstring         ! if true strings are on the 3 three directions
  REAL(DP) :: efield_cart(3)   ! electric field vector in cartesian units
  REAL(DP) :: efield_cry(3)    ! electric field vector in crystal units
  REAL(DP) :: transform_el(3,3)! transformation matrix from cartesian coordinates to normed reciprocal space
  INTEGER, ALLOCATABLE :: mapg_owner(:,:)
  REAL(DP) :: pdl_tot         ! the total phase calculated from bp_c_phase
  INTEGER  :: phase_control! 0 no control, 1 write, 2 read
!
CONTAINS
 
  SUBROUTINE allocate_bp_efield ( ) 

    USE gvect,  ONLY : ngm_g
   ! allocate memory for the Berry's phase electric field
   ! NOTICE: should be allocated ONLY in parallel case, for gdir=1 or 2

   IMPLICIT NONE

   IF ( lberry .OR. lelfield .OR. lorbm ) THEN
      ALLOCATE(mapgp_global(ngm_g,3))
      ALLOCATE(mapgm_global(ngm_g,3))
      ALLOCATE(mapg_owner(2,ngm_g))
   ENDIF

   l_el_pol_old=.false.
   el_pol_acc=0.d0

   RETURN
 END SUBROUTINE allocate_bp_efield

 SUBROUTINE deallocate_bp_efield

   ! deallocate memory used in Berry's phase electric field calculation

   IMPLICIT NONE

   IF ( lberry .OR. lelfield .OR. lorbm ) THEN
      IF ( ALLOCATED(mapgp_global) ) DEALLOCATE(mapgp_global)
      IF ( ALLOCATED(mapgm_global) ) DEALLOCATE(mapgm_global)
      IF ( ALLOCATED(nx_el) ) DEALLOCATE(nx_el)
      IF ( ALLOCATED(mapg_owner) ) DEALLOCATE (mapg_owner)
   ENDIF

   RETURN
 END SUBROUTINE deallocate_bp_efield

 SUBROUTINE bp_global_map

    !this subroutine sets up the global correspondence map G+1 and G-1

    USE mp,                   ONLY : mp_sum
    USE mp_images,            ONLY : me_image, intra_image_comm
    USE gvect,                ONLY : ngm_g, g, ngm, ig_l2g
    USE fft_base,             ONLY : dfftp
    USE cell_base,            ONLY : at

    IMPLICIT NONE

    INTEGER :: ig, mk1,mk2,mk3, idir, imk(3)
    INTEGER, ALLOCATABLE :: ln_g(:,:,:)
    INTEGER, ALLOCATABLE :: g_ln(:,:)

    IF ( .NOT.lberry .AND. .NOT. lelfield .AND. .NOT. lorbm ) RETURN
    ! set up correspondence ln_g ix,iy,iz ---> global g index in
    ! (for now...) coarse grid
    ! and inverse realtion global g (coarse) to ix,iy,iz

    ALLOCATE(ln_g(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3))
    ALLOCATE(g_ln(3,ngm_g))

    ln_g(:,:,:)=0!it means also not found
    DO ig=1,ngm
       mk1=nint(g(1,ig)*at(1,1)+g(2,ig)*at(2,1)+g(3,ig)*at(3,1))
       mk2=nint(g(1,ig)*at(1,2)+g(2,ig)*at(2,2)+g(3,ig)*at(3,2))
       mk3=nint(g(1,ig)*at(1,3)+g(2,ig)*at(2,3)+g(3,ig)*at(3,3))
       ln_g(mk1,mk2,mk3)=ig_l2g(ig)
    ENDDO
    CALL mp_sum(ln_g(:,:,:),intra_image_comm)


    g_ln(:,:)= 0!it means also not found
    DO ig=1,ngm
       mk1=nint(g(1,ig)*at(1,1)+g(2,ig)*at(2,1)+g(3,ig)*at(3,1))
       mk2=nint(g(1,ig)*at(1,2)+g(2,ig)*at(2,2)+g(3,ig)*at(3,2))
       mk3=nint(g(1,ig)*at(1,3)+g(2,ig)*at(2,3)+g(3,ig)*at(3,3))
       g_ln(1,ig_l2g(ig))=mk1
       g_ln(2,ig_l2g(ig))=mk2
       g_ln(3,ig_l2g(ig))=mk3
    ENDDO
    CALL mp_sum(g_ln(:,:),intra_image_comm)

!loop on direction
    DO idir=1,3
!for every g on global array find G+1 and G-1 and put on
       DO ig=1,ngm_g
          imk(:)=g_ln(:,ig)
          imk(idir)=imk(idir)+1
!table array
          mapgp_global(ig,idir)=ln_g(imk(1),imk(2),imk(3))
          imk(idir)=imk(idir)-2
          mapgm_global(ig,idir)=ln_g(imk(1),imk(2),imk(3))
       ENDDO
    ENDDO

    mapg_owner=0
    DO ig=1,ngm
       mapg_owner(1,ig_l2g(ig))=me_image+1
       mapg_owner(2,ig_l2g(ig))=ig
    END DO
    call mp_sum(mapg_owner, intra_image_comm)

    DEALLOCATE(ln_g,g_ln)

    RETURN

  END SUBROUTINE bp_global_map

END MODULE bp
