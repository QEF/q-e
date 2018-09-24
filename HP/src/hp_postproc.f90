!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_postproc
!-----------------------------------------------------------------------
  !
  ! This is a post-processing routine for the calculation of Hubbard parameters.  
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, tau, ityp, atm, ntyp => nsp
  USE cell_base,      ONLY : alat, at, bg
  USE io_files,       ONLY : prefix, tmp_dir
  USE io_global,      ONLY : stdout
  USE control_flags,  ONLY : iverbosity
  USE lsda_mod,       ONLY : nspin
  USE matrix_inversion
  USE ldaU,           ONLY : Hubbard_l, Hubbard_lmax, is_hubbard,    &
                             lda_plus_u_kind ! dist_s, ityp_s  DFT+U+V
  USE ldaU_hp,        ONLY : nath, nath_sc, todo_atom, background,   &
                             skip_type, merge_type, skip_atom,       &
                             tmp_dir_save, at_equiv_criterium, magn, &
                             nath_pert, ityp_new, ntyp_new, atm_new, &
                             num_neigh, lmin, rmax, nq1, nq2, nq3
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: dist_sc(:,:),    & ! interatomic distances in the virtual supercell
                           tau_sc(:,:),     & ! coordinates of virtual atoms
                           tau_bohr(:,:),   & ! Atomic coordinates in Cartesian framework (in Bohr)
                           at_bohr(:,:),    & ! Lattice vectors in in Cartesian framework (in Bohr)
                           chi0(:,:),       & ! chi0 including the background terms
                           chi(:,:),        & ! chi including the background terms
                           chi0bg(:,:),     & ! chi0 including the background terms
                           chibg(:,:),      & ! chi including the background terms
                           inv_chi0bg(:,:), & ! inverse of chi0bg
                           inv_chibg(:,:),  & ! inverse of chibg
                           U(:,:)             ! Hubbard U matrix
  !
  INTEGER :: na, nb, nc, nd, & ! dummy indices running over atoms
             ipol,           & ! dummy index running over 3 cartesian coordinates
             nt,             & ! dummy index running over atomic types
             iunitU,         & ! unit number
             na_sc,          & ! dummy index
             nath_scbg         ! the size of chi matrix including the background term
  !
  INTEGER, ALLOCATABLE :: ityp_sc(:),   & ! type of atoms in the supercell
                          ityp_sc0(:),  & ! type of atoms in the supercell
                          spin_sc(:),   & ! spin of atoms in the supercell
                          spin(:),      & ! spin of atoms in the unit cell
                          auxindex(:,:)   ! auxiliary index (needed for mapping
                                          ! between two supercells)
  !
  REAL(DP), PARAMETER :: eps1 = 9.d-3, & ! various thresholds
                         eps2 = 4.d-1, & ! threshold for the magnetization
                         eps3 = 1.d-4    ! the same threshold for the comparison of distances
                                         ! as in PW/src/inter_V.f90 DFT+U+V
  !
  CHARACTER(len=50) :: filenameU
  INTEGER, EXTERNAL :: find_free_unit
  !
  CALL start_clock('hp_calc_U')
  !
  WRITE( stdout, '(/5x,"Post-processing calculation of Hubbard parameters ...",/)')
  !
  ! Allocate various arrays
  CALL alloc_pp()
  !
  ! Read chi0 and chi from file
  CALL read_chi()
  !
  ! If we merge types of atoms (e.g. Ni_up and Ni_down) 
  ! then we have to keep track of their total magnetization
  ! in order to be able to distinguish them
  !
  CALL merge_types_and_determine_spin()
  !
  ! Generate the virtual atoms in order 
  ! to mimic a virtual supercell
  CALL gen_virt_atoms()
  !
  ! Compute interatomic distances
  ! between virtual atoms
  CALL atomic_dist()
  !
  ! Average similar elements in chi0 and chi
  CALL average_similar_elements(chi0) 
  CALL average_similar_elements(chi)
  !
  ! Reconstruct full chi0 and chi using symmetry
  CALL reconstruct_full_chi(chi0)
  CALL reconstruct_full_chi(chi)
  !
  ! Add a background correction if needed
  CALL background_correction(chi0, chi0bg)
  CALL background_correction(chi, chibg)
  !
  ! Invert the matrices chi0 and chi
  CALL invmat (nath_scbg, chi0bg, inv_chi0bg) 
  CALL invmat (nath_scbg, chibg, inv_chibg)
  !
  ! Calculate U and write it to file
  CALL calculate_U()
  !
  ! Deallocate various arrays
  CALL dealloc_pp()
  !
  CALL stop_clock('hp_calc_U')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE alloc_pp
  !
  ! This routine allocates various arrays
  !
  IMPLICIT NONE
  !
  ALLOCATE ( dist_sc(nath_sc, nath_sc) )
  ALLOCATE ( tau_sc(3, nath_sc) )
  ALLOCATE ( ityp_sc(nath_sc) )
  ALLOCATE ( ityp_sc0(nath_sc) )
  ALLOCATE ( tau_bohr(3, nath_sc) )
  ALLOCATE ( at_bohr(3, 3) )
  ALLOCATE ( spin(nat) )
  ALLOCATE ( spin_sc(nath_sc) )
  ALLOCATE ( auxindex(nath,nath_sc) )
  !
  IF (trim(background)=='no') THEN
     nath_scbg = nath_sc
  ELSEIF (trim(background)=='neutral') THEN
     nath_scbg = nath_sc + 1
  ELSE
     CALL errore ('alloc_pp', 'Wrong background', 1)
  ENDIF
  !
  ALLOCATE ( chi(nath_sc, nath_sc) )
  ALLOCATE ( chi0(nath_sc, nath_sc) )
  ALLOCATE ( chibg(nath_scbg, nath_scbg) )
  ALLOCATE ( chi0bg(nath_scbg, nath_scbg) )
  ALLOCATE ( inv_chibg(nath_scbg, nath_scbg) )
  ALLOCATE ( inv_chi0bg(nath_scbg, nath_scbg) )
  ALLOCATE ( U(nath_scbg, nath_scbg) )
  !
  RETURN 
  !
END SUBROUTINE alloc_pp

SUBROUTINE dealloc_pp
  !
  ! This routine deallocates various arrays
  !
  IMPLICIT NONE
  !
  DEALLOCATE (dist_sc)
  DEALLOCATE (tau_sc)
  DEALLOCATE (ityp_sc)
  DEALLOCATE (ityp_sc0)
  DEALLOCATE (tau_bohr)
  DEALLOCATE (at_bohr)
  DEALLOCATE (spin)
  DEALLOCATE (spin_sc)
  DEALLOCATE (auxindex)
  DEALLOCATE (chi)
  DEALLOCATE (chi0)
  DEALLOCATE (chibg)
  DEALLOCATE (chi0bg)
  DEALLOCATE (inv_chibg)
  DEALLOCATE (inv_chi0bg)
  DEALLOCATE (U)
  !
  RETURN
  !
END SUBROUTINE dealloc_pp

SUBROUTINE read_chi
  !
  ! Read chi0 and chi from file 
  !
  IMPLICIT NONE
  !
  INTEGER :: iunitchi ! unit number
  CHARACTER(len=50)  :: filenamechi
  CHARACTER(len=256) :: tempfile
  LOGICAL :: exst
  !
  chi0(:,:) = (0.0d0, 0.0d0)
  chi(:,:)  = (0.0d0, 0.0d0)
  !
  ! Find and open unit
  !  
  iunitchi = find_free_unit()
  filenamechi = TRIM(prefix) // ".chi.dat"
  tempfile = TRIM(tmp_dir) // TRIM(filenamechi)
  !
  INQUIRE (file = tempfile, exist = exst)
  IF (.NOT.exst) THEN
     WRITE( stdout, * ) "    WARNING: " // TRIM(tempfile) // " does not exist !!!"
     CALL errore('hp_calc_U','Missing file with the chi matrix',1)
  ENDIF
  !
  OPEN(iunitchi, file = tempfile, form = 'formatted', status = 'unknown')
  !
  READ(iunitchi,*)
  DO na = 1, nath_sc
     READ(iunitchi,*) (chi0(na,nb), nb=1,nath)
  ENDDO
  !
  READ(iunitchi,*)
  READ(iunitchi,*)
  DO na = 1, nath_sc
     READ(iunitchi,*) (chi(na,nb), nb=1,nath)
  ENDDO
  !
  CLOSE(iunitchi)
  !
  RETURN
  !
END SUBROUTINE read_chi

SUBROUTINE merge_types_and_determine_spin
  !
  ! If we merge types of atoms (e.g. Ni_up and Ni_down) 
  ! then we have to keep track of their total magnetization
  ! (orientation of spin) in order to be able to distinguish them.
  !
  IMPLICIT NONE
  !
  IF ( at_equiv_criterium == 2 .OR. at_equiv_criterium == 3 ) THEN
     !
     ! If it was requested in the input, re-assign the specific type
     ! to another type (e.g., Ni_down to Ni_up)
     !
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt) .AND. skip_type(nt)) THEN
           ityp_new(na) = merge_type(nt)
        ENDIF
     ENDDO
     !
  ENDIF
  !
  ! Determine the spin attribute for every atom
  ! based on the magnetization
  !
  IF ( nspin == 2 ) THEN
     !
     DO na = 1, nat
        IF ( magn(na) > eps2 ) THEN
           spin(na) = 1
        ELSEIF ( magn(na) < -eps2 ) THEN
           spin(na) = -1
        ELSE
           spin(na) = 0
        ENDIF
     ENDDO
     !
  ELSE
     !
     spin(:) = 1
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE merge_types_and_determine_spin

SUBROUTINE gen_virt_atoms
  !
  ! This routine generates virtual atoms in order to mimic 
  ! the supercell and in order to reconstruct the full chi0 and chi
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k
  !
  ! Convert atomic positions in the real cell
  ! from cartesian to crystal framework
  !
  CALL cryst_to_cart( nat, tau, bg, -1 )
  !
  ! We generate virtual atoms in the same order,
  ! as were generated virtual cells (R_points.f90), 
  ! such that virtual atoms are in corresponding virtual cells.
  !
  na_sc = 0
  !
  DO i = 1, nq1
    DO j = 1, nq2
      DO k = 1, nq3
         !
         DO na = 1, nath
            !
            na_sc = na_sc + 1
            !
            tau_sc(1,na_sc) = (tau(1,na) + DBLE(i-1)) / DBLE(nq1)
            tau_sc(2,na_sc) = (tau(2,na) + DBLE(j-1)) / DBLE(nq2)
            tau_sc(3,na_sc) = (tau(3,na) + DBLE(k-1)) / DBLE(nq3)
            !
            ityp_sc(na_sc)  = ityp_new(na)
            ityp_sc0(na_sc) = ityp(na)
            !
            spin_sc(na_sc) = spin(na)
            !
         ENDDO
         !
      ENDDO
    ENDDO
  ENDDO
  !
  ! Convert atomic positions in the real cell and in the virtual 
  ! supercell from crystal to cartesian framework
  !
  CALL cryst_to_cart( nat, tau, at, +1 )
  !
  ! Rescale the lattice vectors in order to mimic
  ! the nq1 x nq2 x nq3 supercell
  !
  at(:,1) = at(:,1) * DBLE(nq1)  ! 1st lattice vector a1
  at(:,2) = at(:,2) * DBLE(nq2)  ! 2nd lattice vector a2
  at(:,3) = at(:,3) * DBLE(nq3)  ! 3rd lattice vector a3
  !
  CALL cryst_to_cart( nath_sc, tau_sc, at, +1 )
  !
  ! Transform the lattice vectors and atomic coordinates 
  ! from alat units to Bohr units
  !
  at_bohr(:,:)  = at(:,:) * alat
  tau_bohr(:,:) = tau_sc(:,:) * alat 
  !
  RETURN
  !
END SUBROUTINE gen_virt_atoms

SUBROUTINE atomic_dist()
  !
  ! This routine computes the interatomic distances
  ! between atoms in the virtual supercell of size
  ! nq1 x nq2 x nq3
  !
  IMPLICIT NONE
  !
  REAL(DP) :: daux  ! auxiliary variable
  INTEGER :: i, j, k, dimn
  LOGICAL, ALLOCATABLE :: found(:)
  !
  DO na = 1, nath_sc
     !
     dist_sc(na, na) = 0.d0
     !
     DO nb = na+1, nath_sc
        !
        dist_sc(na, nb) = 1000.d0
        !
        DO i = -1, 1
          DO j = -1, 1
            DO k = -1, 1
               !
               daux = 0.d0
               !
               DO ipol = 1, 3
                  !
                  daux = daux + ( tau_bohr(ipol,na) &
                                - tau_bohr(ipol,nb) &
                        - DBLE(i) * at_bohr(ipol,1) &
                        - DBLE(j) * at_bohr(ipol,2) &
                        - DBLE(k) * at_bohr(ipol,3) )**2
                  !
               ENDDO
               !
               daux = DSQRT(daux)
               !
               dist_sc(na, nb) = MIN( dist_sc(na, nb), daux )
               !
            ENDDO
          ENDDO
        ENDDO
        !
        dist_sc(nb, na) = dist_sc(na, nb)
        !
     ENDDO
     !
  ENDDO
  !
  IF (lda_plus_u_kind.EQ.2) THEN
     !
     CALL errore("hp_postproc", 'The HP code does not support lda_plus_u_kind=2',1)
     ! 
     ! Number of atoms in the 3x3x3 supercell
     !dimn = nat * (3**3.0d0)
     !
     !ALLOCATE(found(dimn))
     !
     ! Perform a mapping between the nb index of V(na,nb) 
     ! of the nq1xnq2xnq3 supercell with the nc index of V(na,nc)
     ! of the 3x3x3 supercell used in DFT+U+V of PWscf.
     !
     !DO na = 1, nath
     !   found(:) = .FALSE.
     !   DO nb = 1, nath_sc
     !      DO nc = 1, dimn
     !         IF ( ityp_sc0(nb).EQ.ityp_s(nc)                .AND. &
     !              ABS(dist_sc(na,nb)-dist_s(na,nc)).LT.eps3 .AND. &
     !              .NOT.found(nc) ) THEN
     !              !
     !              ! Mapping of index nb to nc
     !              auxindex(na,nb) = nc
     !              found(nc) = .TRUE.
     !              GO TO 9
     !              !
     !         ENDIF 
     !      ENDDO
     !      ! In the case when the mapping is not reached nullify the index
     !      auxindex(na,nb) = 0
!9          CONTINUE         
     !   ENDDO
     !ENDDO   
     !      
     !DEALLOCATE(found)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE atomic_dist

SUBROUTINE average_similar_elements(chi_)
  !
  ! Average similar elements in chi_
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: chi_(nath_sc, nath_sc)
  !
  LOGICAL :: condition
  INTEGER, ALLOCATABLE :: aux(:) ! auxiliary array for indices
  INTEGER  :: nat0
  REAL(DP) :: chi_aux
  !
  ALLOCATE (aux(nath_sc))
  !
  DO na = 1, nath_sc
     !
     DO nb = 1, nath_sc
        !
        IF ( chi_(na,nb).NE.0.d0 ) THEN
           !
           nat0 = 1
           chi_aux = chi_(na,nb)
           aux(nat0) = na
           !
           ! For a given column of chi_, check if there are
           ! similar by value elements.
           !
           DO nc = 1, nath_sc
              !
              condition = nc.NE.na                              .AND. &
                          ityp_sc(nc).EQ.ityp_sc(na)            .AND. &
                          spin_sc(nc).EQ.spin_sc(na)            .AND. & 
                          chi_(nc,nb).NE.0.d0                   .AND. &
                          dist_sc(na,nb).GT.0.d0                .AND. &
                          dist_sc(nc,nb).GT.0.d0                .AND. &
                          ABS(dist_sc(nc,nb)-dist_sc(na,nb)).LE.eps1 
              !
              IF (condition) THEN
                 !
                 nat0 = nat0 + 1
                 aux(nat0) = nc
                 chi_aux = chi_aux + chi_(nc,nb)
                 !
               ENDIF
               !
           ENDDO
           !
           IF ( nat0 > 1 ) THEN
              !
              ! Divide by the number of similar chi_, which were found
              !
              DO nc = 1, nat0
                 chi_(aux(nc), nb) = chi_aux / DBLE(nat0)
              ENDDO
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE (aux) 
  !
  RETURN
  !
END SUBROUTINE average_similar_elements

SUBROUTINE reconstruct_full_chi(chi_)
  !
  ! This routine reconstructs the full matrix of 
  ! the susceptibility using symmetry.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: chi_(nath_sc, nath_sc)
  LOGICAL :: condition
  INTEGER :: nat0
  !
  DO na = 1, nath_sc
     !
     DO nb = 1, nath_sc
        !
        IF ( chi_(na,nb).EQ.0.d0 ) THEN
           !
           nat0 = 0
           !
           DO nc = 1, nath_sc
              !
              IF ( ityp_sc(nc).EQ.ityp_sc(nb) ) THEN
                 !
                 DO nd = 1, nath_sc
                    !
                    IF ( chi_(nd,nc).NE.0.d0 ) THEN
                       !
                       nat0 = nc
                       go to 35
                       !
                    ENDIF
                    !
                 ENDDO
                 !
              ENDIF
              !
           ENDDO
           !
35 continue
           !
           DO nd = 1, nath_sc
              !
              condition = ityp_sc(nd).EQ.ityp_sc(na)                     .AND. &
                          ityp_sc(nat0).EQ.ityp_sc(nb)                   .AND. &
                          ABS(dist_sc(nd,nat0)-dist_sc(na,nb)).LE.eps1   .AND. & 
                          spin_sc(na)*spin_sc(nb).EQ.spin_sc(nat0)*spin_sc(nd)
              !
              IF (condition) THEN
                 !
                 chi_(na,nb) = chi_(nd,nat0)
                 !
              ENDIF
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ! Symmetrization
  !
  DO na = 1, nath_sc
     DO nb = 1, nath_sc
        !
        chi_(na,nb) = 0.5d0 * ( chi_(na,nb) + chi_(nb,na) )
        chi_(nb,na) = chi_(na,nb)
        !
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE reconstruct_full_chi

SUBROUTINE background_correction(chi_, chibg_)
  !
  ! This subroutine adds a background correction (if needed)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: chi_(nath_sc, nath_sc)
  REAL(DP), INTENT(OUT) :: chibg_(nath_scbg, nath_scbg)
  !
  REAL(DP) :: sum_last_column
  !
  chibg_(:,:) = 0.0d0
  sum_last_column = 0.0d0
  !
  ! Copying 
  !
  DO na = 1, nath_sc
     DO nb = 1, nath_sc
        chibg_(na,nb) = chi_(na,nb)
     ENDDO
  ENDDO 
  !
  IF (trim(background)=='neutral') THEN
     !
     ! Compute the last additional column and 
     ! the last additional row in the matrix chibg_
     ! due to the background correction 
     ! 
     DO na = 1, nath_sc
        !
        ! Sum up all the element in a given row with the opposite sign
        ! in order to compute the additional element in the same row
        ! (neutrality condition)
        !
        DO nb = 1, nath_sc
           chibg_(na,nath_scbg) = chibg_(na,nath_scbg) - chi_(na,nb)
        ENDDO
        !
        ! Use the symmetry for the off-diagonal elements
        !
        chibg_(nath_scbg,na) = chibg_(na,nath_scbg)
        !
        ! Sum up all the elements in the last additional column 
        ! in chibg_ with the opposite sign (neutrality condition)
        !
        sum_last_column = sum_last_column - chibg_(na,nath_scbg)
        !
     ENDDO
     !
     ! Assign the last element in the matrix
     !
     chibg_(nath_scbg,nath_scbg) = sum_last_column
     !
     ! In order to be able to perform an inversion of the chibg_ matrix 
     ! (which is singular by construction), we add a small number
     ! to each element of the chibg_ matrix. We can do so because 
     ! we are interested in the difference of the two inverse matrices, 
     ! hence this extra small number will cancel out and will not effect 
     ! the final value of U.
     !
     DO na = 1, nath_scbg
        DO nb = 1, nath_scbg
           chibg_(na,nb) = chibg_(na,nb) + 0.01d0
        ENDDO
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE background_correction

SUBROUTINE calculate_U()
  !
  ! This routine calculates Hubbard U from the inverse matrices
  ! of the susceptibility and writes it to file.
  !
  IMPLICIT NONE
  INTEGER :: nt1, nt2
  !
  ! Find and open unit to write info
  iunitU = find_free_unit()
  filenameU = trim(prefix) // ".Hubbard_parameters.dat"
  OPEN(iunitU, file = filenameU, form = 'formatted', status = 'unknown')
  !
  ! Calculate U = CHI0^{-1} - CHI^{-1}
  !
  DO na = 1, nath_scbg
     DO nb = 1, nath_scbg
        U(na,nb) = inv_chi0bg(na,nb) - inv_chibg(na,nb)
     ENDDO
  ENDDO
  !
  IF ( at_equiv_criterium == 1 .AND. ntyp_new > ntyp ) THEN
     WRITE(iunitU,'(2x,"Warning: The Hubbard parameters listed below were computed by treating some")')
     WRITE(iunitU,'(2x,"         equivalent (by type) Hubbard atoms as if they were non-equivalent")')
     WRITE(iunitU,'(2x,"         (when doing the averaging and reconstruction of the response matrices)")')
     WRITE(iunitU,'(2x,"         because their unperturbed occupations differ. If you want to disable this")')
     WRITE(iunitU,'(2x,"         option then set disable_type_analysis=.true. and redo the post-processing.")')
  ENDIF
  !
  ! Write Hubbard U 
  !
  WRITE(iunitU,'(/2x,"=-------------------------------------------------------------------=",/)')
  WRITE(iunitU,'(27x,"Hubbard U parameters:",/)')
  WRITE(iunitU, '(7x,"site n.",2x,"type",2x,"label",2x,"spin",2x,"new_type",2x,"new_label",2x,"Hubbard U (eV)")')
  DO na = 1, nat
     nt1 = ityp(na)
     nt2 = ityp_new(na)
     IF ( is_hubbard(nt1) ) THEN
        WRITE(iunitU,'(7x,i3,6x,i3,3x,a4,4x,i2,4x,i3,9x,a4,3x,f10.4)') &
                & na, nt1, atm(nt1), spin(na), nt2, atm_new(nt2), U(na,na)
     ENDIF
  ENDDO
  WRITE(iunitU,'(/2x,"=-------------------------------------------------------------------=",/)')
  !
  ! Write Hubbard V
  !
  !IF ( lda_plus_u_kind.EQ.2 ) CALL write_Hubbard_V()
  !
  ! Write the information about the response matrices chi0 and chi,
  ! about their inverse matrices, and about the entire matrix of 
  ! Hubbard parameters.
  !
  IF ( iverbosity > 1 ) THEN
     !
     IF (trim(background)=='neutral') THEN
        ! Subtract the previously added constant for
        ! the inversion of the matrices
        DO na = 1, nath_scbg
           DO nb = 1, nath_scbg
              chi0bg(na,nb) = chi0bg(na,nb) - 0.01d0
              chibg(na,nb)  = chibg(na,nb)  - 0.01d0
           ENDDO
        ENDDO
     ENDIF
     !
     ! Write to file chi0
     WRITE(iunitU,'(/10x,"chi0 matrix :")')
     DO na = 1, nath_scbg
        WRITE(iunitU,'(8(x,f11.6))') (chi0bg(na,nb), nb=1,nath_scbg)
        WRITE(iunitU,*)
     ENDDO
     !
     ! Write to file chi
     WRITE(iunitU,'(/10x,"chi matrix :")')
     DO na = 1, nath_scbg
        WRITE(iunitU,'(8(x,f11.6))') (chibg(na,nb), nb=1,nath_scbg)
        WRITE(iunitU,*)
     ENDDO
     !
     ! Write to file inv_chi0
     WRITE(iunitU,'(/8x,"chi0^{-1} matrix :")')
     DO na = 1, nath_scbg
        WRITE(iunitU,'(8(x,f11.6))') (inv_chi0bg(na,nb), nb=1,nath_scbg)
        WRITE(iunitU,*)
     ENDDO
     !
     ! Write to file inv_chi
     WRITE(iunitU,'(/8x,"chi^{-1} matrix :")')
     DO na = 1, nath_scbg
        WRITE(iunitU,'(8(x,f11.6))') (inv_chibg(na,nb), nb=1,nath_scbg)
        WRITE(iunitU,*)
     ENDDO
     !
     ! Writing U to file
     WRITE(iunitU,'(/9x,"Hubbard matrix :")')
     DO na = 1, nath_scbg
        WRITE(iunitU,'(8(x,f11.6))') (U(na,nb), nb=1,nath_scbg)
        WRITE(iunitU,*)
     ENDDO
     !
  ENDIF
  !
  CLOSE(iunitU)
  !
  RETURN
  !
END SUBROUTINE calculate_U

SUBROUTINE write_Hubbard_V()
  !
  ! Write information about the Hubbard_V parameters
  ! in the order of increase of interatomic distances.
  !
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: dist(:), distord(:)
  INTEGER,  ALLOCATABLE :: typeord(:), indexord(:)
  INTEGER :: ne, nc_min, ipol, tempunit, counter
  CHARACTER(len=80) :: tempfile
  REAL(DP) :: auxdist
  !
  ! Find and open unit to write Hubbard_V parameters
  tempunit = find_free_unit()
  tempfile = TRIM("parameters.out")
  OPEN(tempunit, file = tempfile, form = 'formatted', status = 'unknown')
  !
  WRITE(iunitU,'(/27x,"Hubbard V parameters:")')
  WRITE(iunitU,'(22x,"(adapted for a 3x3x3 supercell)",/)')
  WRITE(iunitU,*) '           Atom 1 ', '    Atom 2 ', '   Distance (Bohr) ', ' Hubbard V (eV)'
  WRITE(iunitU,*)
  !
  WRITE(tempunit,*) '# Atom 1 ', ' Atom 2 ', ' Hubbard V (eV)'
  !
  ALLOCATE(dist(nath_sc))
  ALLOCATE(distord(nath_sc))
  ALLOCATE(typeord(nath_sc))
  ALLOCATE(indexord(nath_sc))
  !
  DO na = 1, nath
     !
     ! All distances from the atom na 
     dist(:) = dist_sc(na,:) 
     !
     ! Order distances
     !
     nc = 1
     DO nb = 1, nath_sc
        IF ( dist(nb).GT.0.d0 ) THEN
            distord(1)  = dist(nb)
            typeord(1)  = ityp_sc(nb)
            indexord(1) = nb
            GO TO 10
        ENDIF
     ENDDO
     !
10   CONTINUE
     !
     DO nb = 1, nath_sc
        IF ( dist(nb).GT.distord(1) ) THEN
           nc = nc + 1
           distord(nc)  = dist(nb)
           typeord(nc)  = ityp_sc(nb)
           indexord(nc) = nb
           GO TO 11
        ENDIF
     ENDDO
     !
     nc = nc + 1
     distord(nc)  = distord(1)
     typeord(nc)  = typeord(1)
     indexord(nc) = indexord(1)
     !
     DO nb = 1, nath_sc
        IF ( dist(nb).GT.0.d0 .AND. dist(nb).LT.distord(nc) ) THEN
            distord(1)  = dist(nb)
            typeord(1)  = ityp_sc(nb)
            indexord(1) = nb
            GO TO 11
        ENDIF
     ENDDO
     !
11   CONTINUE
     !
     DO nb = 1, nath_sc
        IF ( dist(nb).GT.distord(nc) ) THEN
            nc = nc + 1
            distord(nc)  = dist(nb)
            typeord(nc)  = ityp_sc(nb)
            indexord(nc) = nb
        ELSEIF ( dist(nb).LE.distord(1) .AND. nb.NE.indexord(1) ) THEN
            nc_min = 2
            nc = nc + 1
            DO nd = nc, nc_min, -1
               distord(nd)  = distord(nd-1)
               typeord(nd)  = typeord(nd-1)
               indexord(nd) = indexord(nd-1)
            ENDDO
            distord(1)  = dist(nb)
            typeord(1)  = ityp_sc(nb)
            indexord(1) = nb
        ELSE
            DO nd = nc-1, 1, -1
               IF ( dist(nb).GT.distord(nd) .AND. dist(nb).LE.distord(nd+1) ) THEN
                  IF ( nb.NE.indexord(nd+1) ) THEN
                     nc_min = nd + 2
                     nc = nc + 1
                     DO ne = nc, nc_min,-1
                        distord(ne)  = distord(ne-1)
                        typeord(ne)  = typeord(ne-1)
                        indexord(ne) = indexord(ne-1)
                     ENDDO
                     distord(nd+1)  = dist(nb)
                     typeord(nd+1)  = ityp_sc(nb)
                     indexord(nd+1) = nb 
                     GO TO 12
                  ENDIF
               ENDIF
            ENDDO
12          CONTINUE
        ENDIF
     ENDDO
     !
     counter = 0
     !
     DO nb = 1, nath_sc
        IF ( auxindex(na,indexord(nb)) > 0 ) THEN
           counter = counter + 1
           WRITE(iunitU,'(11x,i3,x,a4,x,i5,x,a4,2x,f12.6,4x,f10.4)') &
                 na, atm_new(ityp_sc(na)), auxindex(na,indexord(nb)), atm_new(typeord(nb)), &
                 dist_sc(na,indexord(nb)), U(na,indexord(nb))
           IF ( nb.LE.(num_neigh+1)              .AND. &     
                dist_sc(na,indexord(nb)).LE.rmax .AND. & 
                Hubbard_l(ityp(na)).GE.lmin ) THEN
              WRITE(tempunit,'(3x,i3,4x,i5,3x,f10.4)') &
                    na, auxindex(na,indexord(nb)), U(na,indexord(nb))
           ENDIF
        ENDIF
     ENDDO 
     WRITE(iunitU,'()')
     !
  ENDDO ! na
  !
  WRITE(iunitU,'(/2x,"=-------------------------------------------------------------------=",/)')
  !
  IF ((counter-1).GE.num_neigh) &
    WRITE(iunitU,'(2x,i3,2x,"nearest neighbors of every atom are written to the file parameters.dat")') num_neigh
  !
  DEALLOCATE(dist)
  DEALLOCATE(distord)
  DEALLOCATE(typeord)
  DEALLOCATE(indexord)
  !
  CLOSE(tempunit)
  !
  RETURN
  !
END SUBROUTINE write_Hubbard_V

END SUBROUTINE hp_postproc
