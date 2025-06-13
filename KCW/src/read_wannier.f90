!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
subroutine read_wannier ()
  !-----------------------------------------------------------------------
  !
  !! This routine read the relevant infos from a previous Wannier90 run
  !! The unitary matrix U are stored in unimatrx, unimatrx, unimatrx_opt.
  !! U(i,j) = <KS_i|wann_j> 
  !! Two cases are addressed:
  !! 1) A standard wannierization (unique manifold)
  !! 2) A separate wannierization where occ and emp states are treated independently. 
  !! In both cases I define a unique U whose dimension is num_wann x num_wann 
  !! and eventually a U_opt whose dimension in num_bands x num_wann
  !
  USE control_kcw,          ONLY : l_unique_manifold
  !
  IMPLICIT NONE 
  !
  IF (l_unique_manifold) THEN 
    !
    ! ... standard wannierization (occ + emp all together) 
    CALL read_wannier_unique_manifold ( )
    !
  ELSE
   !
   ! ... Separate wannierization (No occ-emp mixing)
   CALL read_wannier_two_manifold ( )
   !
  ENDIF
  !
END subroutine read_wannier
  !
  !-----------------------------------------------------------------------
  subroutine read_wannier_unique_manifold ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : unimatrx, seedname, has_disentangle, &
                                   unimatrx_opt, num_wann, kcw_iverbosity
  USE mp_global,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : ionode, ionode_id
  USE klist,                ONLY : nkstot, xk
  USE cell_base,            ONLY : bg
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,      ONLY : npol, nspin_lsda, nspin_gga, nspin_mag
  USE wvfct,                ONLY : nbnd
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE 
  !
  ! Local Variable
  !
  INTEGER :: ik
  ! ... the kpoint index
  !
  INTEGER :: num_wann_file, num_kpoint_file, num_ks_file
  ! ... the manifold space
  !
  REAL (DP) xk_file(3,nkstot), check
  ! ... The unitary matrix elemnts 
  !
  CHARACTER (len=256) :: u_file, u_dis_file, dum
  ! ... the name of the file containing the U matrix
  !
  INTEGER i, j, nkstot_eff
  ! 
  INTEGER ierr
  !
  !! Here and in W90 We treat always one spin component at the time.
  !! While in PWscf nkstot contain k point for spin up (the first nkstot/nspin 
  !! and spin down (the last nkstot.nspin). See also rotate_ks
  !
  num_wann = 0
  !
  IF (nspin == 4) THEN
      nkstot_eff = nkstot
  ELSE IF (nspin ==2)  THEN
      nkstot_eff = nkstot/nspin
  ELSE
      nkstot_eff = nkstot/nspin
  ENDIF
  !
  u_file=TRIM( seedname ) // '_u.mat'
  !
  ! ... The U matrix for empty states 
  !
  IF (has_disentangle) THEN
     !
     u_dis_file=TRIM( seedname ) // '_u_dis.mat'
     !! The Optimal subspace Matrix 
     !
     IF (ionode) THEN
        !
        OPEN (UNIT = 1002, FILE = u_dis_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
        ! 
        IF (ierr /= 0 ) call errore('rotate_orbitals', 'Error while reading Optimal unitary matrix', abs (ierr) )
        !
        READ (1002,*) dum
        READ (1002,*) num_kpoint_file, num_wann_file, num_ks_file
        !
        IF (num_kpoint_file /= nkstot_eff) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Optimal matrix', nkstot_eff)
        !
        IF (num_wann_file /= num_wann .AND. num_wann /= 0) &
              CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Optimal matrix', 1)
        !
        IF (num_ks_file /=  nbnd ) &
             CALL errore ('read_wannier', 'Mismatch between num KS state from PW and Wann90', 1)
        !
     ENDIF
     !
  ENDIF
  !
  IF (ionode) THEN 
     !
     OPEN (UNIT = 1001, FILE =u_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
     ! 
     IF (ierr /= 0 ) call errore('rotate_orbitals', 'reading Empty states unitary matrix', abs (ierr) )
     !
     READ (1001,*) dum 
     READ (1001,*) num_kpoint_file, num_wann_file, num_wann_file
     !
     IF (num_kpoint_file /= nkstot_eff) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Empty matrix', nkstot_eff)
        !
     IF (num_wann_file /= num_wann .AND. num_wann /= 0) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Empty matrix', 1)
     ! 
     num_wann = num_wann_file
     !! Store the numebr of wannier in a global variable
     !
  ENDIF
  !
  CALL mp_bcast( num_wann, ionode_id, intra_image_comm )
  CALL mp_bcast( num_ks_file, ionode_id, intra_image_comm )
  !
  ALLOCATE ( unimatrx( num_wann, num_wann, nkstot_eff) )
  IF ( .NOT. has_disentangle) num_ks_file = num_wann
  ALLOCATE (unimatrx_opt(num_ks_file,num_wann,nkstot_eff))
  !
  unimatrx = CMPLX(0.D0, 0.D0, kind=DP)
  unimatrx_opt = CMPLX(0.D0, 0.D0, kind=DP)
  !
  IF (ionode) THEN
    !
    IF (has_disentangle) THEN 
      !
      DO ik = 1, nkstot_eff
        ! 
        READ (1002, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
        READ (1002,'(f15.10,sp,f15.10)') ((unimatrx_opt(i,j,ik),i=1,num_ks_file),j=1,num_wann)
        !
      ENDDO 
      !
      ! ... transform the kpoints read from U file in cartesian coordinates 
      CALL cryst_to_cart(nkstot_eff, xk_file, bg, 1)
      !
      check = 0.D0
      DO ik=1, nkstot_eff
        check=check+(sum(xk_file(:,ik)-xk(:,ik)))
      ENDDO
      IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
      !
    ELSE
      !
      ! ... If no disentangle than the U-opt is a square matrix of dimension num_wan x num_wann 
      ! ... set it to the identity (in apply_u_mat we anyway apply both Ui_opt and U)
      !
      DO ik = 1, nkstot_eff
        DO i=1, num_wann; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind = DP); ENDDO
      ENDDO
    ENDIF
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_opt, ionode_id, intra_image_comm )
  if ( kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Optimal Matrix READ")') 
  !
  IF (ionode) THEN   
    !
    DO ik = 1, nkstot_eff
      !
      READ (1001, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
      READ (1001,'(f15.10,sp,f15.10)') ((unimatrx(i,j,ik),i=1,num_wann),j=1,num_wann)
      !
    ENDDO
    !
    ! ... transform the kpoints read from U file in cartesian coordinates 
    CALL cryst_to_cart(nkstot_eff, xk_file, bg, 1)
    !
    check = 0.D0
    DO ik=1, nkstot_eff
      check=check+(sum(xk_file(:,ik)-xk(:,ik)))
    ENDDO
    IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx, ionode_id, intra_image_comm )
  !
  CLOSE (1001)
  CLOSE (1002)
  !
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
#if defined (DEBUG)
  ! WRITE
  DO ik =1, nkstot_eff
    WRITE (*, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
    DO i=1,num_ks_file
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx_opt(i,j,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
  DO ik =1, nkstot_eff
    WRITE (*, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
    DO i = 1, num_wann
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx(i,j,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
#endif
  !
  RETURN
  !
END subroutine read_wannier_unique_manifold 
  !
  !
  !-----------------------------------------------------------------------
  subroutine read_wannier_two_manifold ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : unimatrx, have_empty, num_wann_occ, seedname, &
                                   has_disentangle, have_empty, num_wann_emp, & 
                                   unimatrx_opt, num_wann, kcw_iverbosity, spin_component
  USE mp_global,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : ionode, ionode_id
  USE klist,                ONLY : nkstot, xk, nelec, nelup, neldw
  USE cell_base,            ONLY : bg
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_lsda, nspin_gga, nspin_mag
  USE wvfct,                ONLY : nbnd
  USE io_global,            ONLY : stdout
  !
  ! Local Variable
  !
  IMPLICIT NONE
  !
  INTEGER :: ik 
  ! ... the kpoint index
  !
  INTEGER :: num_wann_file, num_kpoint_file, num_ks_file
  ! ... the manifold space
  !
  REAL (DP) U_re, U_im, xk_file(3,nkstot), check
  ! ... The unitary matrix elemnts 
  !
  CHARACTER (len=256) :: u_file, u_dis_file, dum
  ! ... the name of the file containing the U matrix
  !
  INTEGER i, j, norb_occ, nkstot_eff, norb_emp
  ! 
  INTEGER ierr
  !
  COMPLEX(DP), ALLOCATABLE :: unimatrx_occ(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: unimatrx_emp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: unimatrx_emp_opt(:,:,:)
  INTEGER ieff, jeff, nbnd_pw
  !

  !OLD nkstot_eff = nkstot/nspin
  !! Here and in W90 We treat always one spin component at the time.
  !! While in PWscf nkstot contain k point for spin up (the first nkstot/nspin 
  !! and spin down (the last nkstot.nspin). See also rotate_ks
  !
  u_file=TRIM( seedname ) // '_u.mat'
  !
  IF (nspin == 4) THEN
      nkstot_eff = nkstot
      nbnd_pw = nint(nelec)
  ELSE IF (nspin ==2)  THEN
      nkstot_eff = nkstot/nspin
      nbnd_pw = nint(nelup)
      IF (spin_component == 2) nbnd_pw = nint(neldw)
  ELSE
      nkstot_eff = nkstot/nspin
      nbnd_pw= nint(nelec) / 2
  ENDIF 
  !
  IF (ionode) THEN
     !
     OPEN (UNIT = 1001, FILE = u_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
     ! 
     IF (ierr /= 0 ) call errore('rotate_orbitals', 'Error while reading unitary matrix', abs (ierr) )
     !
     READ (1001,*) dum
     READ (1001,*) num_kpoint_file, num_wann_file, num_wann_file
     !
     IF (num_kpoint_file /= nkstot_eff) & 
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U matrix', nkstot_eff)
     !
     IF (num_wann_file /= num_wann_occ .AND. num_wann_occ /= 0) &
              CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U matrix', 1)
     !
     num_wann_occ = num_wann_file
     !! Store the number of occupied wannier in a  global variable
     !
     IF (num_wann_occ /= nbnd_pw) & 
          CALL errore ('read_wannier', 'Mismatch between  num occ bands and num wann', 1)
     !
  ENDIF
  !
  CALL mp_bcast( num_wann_occ, ionode_id, intra_image_comm )
  norb_occ = num_wann_occ ! Store the total number of occupied states
  ALLOCATE (unimatrx_occ(num_wann_occ, num_wann_occ, nkstot_eff))
  !  
  IF ( ionode ) THEN
    !
    check = 0.D0 
    DO ik = 1, nkstot_eff
    !
    !#### In the case occupied and empty manifolds are treated separately
    !#### we assume an insulating systems and NO disentanglement for occupied manifold
      ! 
      READ (1001, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
      !
      DO i = 1, num_wann_occ
        DO j = 1, num_wann_occ
          READ (1001,*) U_re, U_im
          unimatrx_occ(j,i,ik)=CMPLX(U_re, U_im, kind=DP)
        ENDDO
      ENDDO
      !
    ENDDO 
    !
    ! ... transform the kpoints read from U file in cartesian coordinates 
    CALL cryst_to_cart(nkstot_eff, xk_file, bg, 1)
    !
    check = 0.D0
    DO ik=1, nkstot_eff
      check=check+(sum(xk_file(:,ik)-xk(:,ik)))
    ENDDO 
    IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_occ, ionode_id, intra_image_comm )
  !
  CLOSE (1001) 
  !
  num_wann = num_wann_occ
  if (.NOT. have_empty .AND. kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  IF (.NOT. have_empty) THEN
    !
    ! Store the unitary matrix in a global vabiable
    ALLOCATE (unimatrx (num_wann_occ, num_wann_occ, nkstot_eff)) 
    unimatrx = unimatrx_occ
    !
    ! The optimal subspace Matrix
    ALLOCATE (unimatrx_opt (num_wann_occ, num_wann_occ, nkstot_eff)) 
    unimatrx_opt=CMPLX(0.D0,0.D0, kind=DP)
    !
    ! The optimal matrix is simply the identity
    DO i = 1, num_wann_occ; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
    DEALLOCATE (unimatrx_occ)
    !
    RETURN   ! Nothing more to do
    !
  ENDIF
  !
  !####################
  !   EMPTY STATES 
  !###################
  !
  ! ... read unitary matrix empty states ...
  !
  u_file=TRIM( seedname ) // '_emp_u.mat'
  !
  ! ... The U matrix for empty states 
  IF (has_disentangle) THEN
     !
     u_dis_file=TRIM( seedname ) // '_emp_u_dis.mat'
     !! The Optimal subspace Matrix 
     !
     IF (ionode) THEN
        !
        OPEN (UNIT = 1002, FILE = u_dis_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
        ! 
        IF (ierr /= 0 ) call errore('rotate_orbitals', 'Error while reading Optimal unitary matrix', abs (ierr) )
        !
        READ (1002,*) dum
        READ (1002,*) num_kpoint_file, num_wann_file, num_ks_file
        !
        IF (num_kpoint_file /= nkstot_eff) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Optimal matrix', nkstot_eff)
        !
        IF (num_wann_file /= num_wann_emp .AND. num_wann_emp /= 0) &
              CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Optimal matrix', 1)
        !
        IF (num_ks_file /= ( nbnd - num_wann_occ ) ) &
             CALL errore ('read_wannier', 'Mismatch between num empty bands and num empty wann', 1)
        !
     ENDIF
     !
  ENDIF
  !
  IF (ionode) THEN 
     !
     OPEN (UNIT = 1001, FILE =u_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
     ! 
     IF (ierr /= 0 ) call errore('rotate_orbitals', 'reading Empty states unitary matrix', abs (ierr) )
     !
     READ (1001,*) dum 
     READ (1001,*) num_kpoint_file, num_wann_file, num_wann_file
     !
     IF (num_kpoint_file /= nkstot_eff) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Empty matrix', nkstot_eff)
        !
     IF (num_wann_file /= num_wann_emp .AND. num_wann_emp /= 0) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Empty matrix', 1)
     ! 
     num_wann_emp = num_wann_file
     !! Store the numebr of empty wannier in a global variable
     !
  ENDIF
  !
  CALL mp_bcast( num_wann_emp, ionode_id, intra_image_comm )
  CALL mp_bcast( num_ks_file, ionode_id, intra_image_comm )
  norb_emp = num_wann_emp ! store the total numebr of empty variational orbitals
  !
  ALLOCATE ( unimatrx_emp( num_wann_emp, num_wann_emp, nkstot_eff) )
  IF (has_disentangle) ALLOCATE (unimatrx_emp_opt(num_ks_file,num_wann_emp,nkstot_eff))
  !
  IF (has_disentangle) THEN 
    !
    IF (ionode) THEN
      !
      DO ik = 1, nkstot_eff
        ! 
        READ (1002, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
        READ (1002,'(f15.10,sp,f15.10)') ((unimatrx_emp_opt(i,j,ik),i=1,num_ks_file),j=1,num_wann_emp)
        !
      ENDDO 
      !
      ! ... transform the kpoints read from U file in cartesian coordinates 
      CALL cryst_to_cart(nkstot_eff, xk_file, bg, 1)
      !
      check = 0.D0
      DO ik=1, nkstot_eff
        check=check+(sum(xk_file(:,ik)-xk(:,ik)))
      ENDDO
      IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
      !
    ENDIF
    !
    !WRITE(*,*) unimatrx_emp_opt
    CALL mp_bcast( unimatrx_emp_opt, ionode_id, intra_image_comm )
    !
  ENDIF
  !
  if ( kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Optimal Matrix READ")') 
  !
  IF (ionode) THEN   
    !
    DO ik = 1, nkstot_eff
      !
      READ (1001, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
      READ (1001,'(f15.10,sp,f15.10)') ((unimatrx_emp(i,j,ik),i=1,num_wann_emp),j=1,num_wann_emp)
      !
    ENDDO
    !
    ! ... transform the kpoints read from U file in cartesian coordinates 
    CALL cryst_to_cart(nkstot_eff, xk_file, bg, 1)
    !
    check = 0.D0
    DO ik=1, nkstot_eff
      check=check+(sum(xk_file(:,ik)-xk(:,ik)))
    ENDDO
    IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_emp, ionode_id, intra_image_comm )
  !
  CLOSE (1001)
  CLOSE (1002)
  !
  num_wann = num_wann_occ + num_wann_emp
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  ! Store the result in a unique matrix
  ALLOCATE (unimatrx (num_wann, num_wann, nkstot_eff))
  ALLOCATE (unimatrx_opt (nbnd, num_wann, nkstot_eff))
  !
  unimatrx=CMPLX(0.D0,0.D0, kind=DP)
  unimatrx_opt=CMPLX(0.D0,0.D0, kind=DP)
  !
  ! Build the unique disentangle matrix (identity for the occupied state and 
  ! unimatrx_emp_opt for empty states
  !
  ! 1)occ
  DO i = 1, num_wann_occ; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
  !
  ! 2)emp
  IF (.not. has_disentangle) THEN  
    !
    ! ...In this case the optimal matrix is the ideintity (nbnd=num_wann)
    DO i = 1, num_wann; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
    !
  ELSE
    !
    DO i = 1,num_ks_file
      !
      ieff=i+num_wann_occ
      DO j = 1, num_wann_emp
        jeff=j+num_wann_occ
        unimatrx_opt(ieff,jeff,:) = unimatrx_emp_opt(i,j,:)
      ENDDO
      !
    ENDDO
    !
  ENDIF 
  !
  ! Build the unique unitary matrix
  !
  ! 1)occ
  unimatrx(1:num_wann_occ, 1:num_wann_occ, :) = unimatrx_occ(1:num_wann_occ, 1:num_wann_occ, :)
  !
  ! 2)emp
  DO i = 1, num_wann_emp
    !
    ieff=i+num_wann_occ
    DO j = 1, num_wann_emp
      ! 
      jeff=j+num_wann_occ
      unimatrx(ieff,jeff,:) = unimatrx_emp(i,j,:)
      !
    ENDDO
    !
  ENDDO
  !
  DEALLOCATE (unimatrx_emp)
  IF (has_disentangle) DEALLOCATE (unimatrx_emp_opt )
  !

#if defined (DEBUG)
  ! WRITE
  DO ik =1, nkstot_eff
    WRITE (*, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
    DO i=1,nbnd
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx_opt(j,i,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
  DO ik =1, nkstot_eff
    WRITE (*, *) xk_file(1,ik), xk_file(2,ik), xk_file(3,ik)
    DO i = 1, num_wann
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx(j,i,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
#endif
  !
  RETURN
  !
END subroutine read_wannier_two_manifold
