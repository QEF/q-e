!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
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
  USE ldaU,           ONLY : Hubbard_n, Hubbard_l, Hubbard_lmax, is_hubbard, &
                             lda_plus_u_kind, dist_s, ityp_s, &
                             Hubbard_projectors, num_uc
  USE ldaU_hp,        ONLY : nath, nath_sc, todo_atom, background,   &
                             skip_type, equiv_type, skip_atom,       &
                             tmp_dir_save, find_atpert, magn,        &
                             nath_pert, ityp_new, ntyp_new, atm_new, &
                             num_neigh, lmin, rmax, nq1, nq2, nq3,   &
                             determine_num_pert_only, dist_thr
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
                           Hubbard_matrix(:,:)! Matrix of Hubbard parameters
  !
  INTEGER :: na, nb, nc, nd, & ! dummy indices running over atoms
             ipol,           & ! dummy index running over 3 cartesian coordinates
             nt,             & ! dummy index running over atomic types
             unithub,        & ! unit number
             unithub2,       & ! unit number
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
  REAL(DP), PARAMETER :: eps1 = 4.d-1     ! threshold for the magnetization
  !
  CHARACTER(len=256) :: filenamehub, filenamehub2
  INTEGER, EXTERNAL :: find_free_unit
  LOGICAL :: determine_indices_only
  !
  CALL start_clock('hp_postproc')
  !
  IF (determine_num_pert_only) THEN
     IF (lda_plus_u_kind==2) THEN 
        ! DFT+U+V: determine indices of couples for Hubbard V, without computing U and V.
        ! This is useful when DFT+U+V is used for large supercells. So one can determine
        ! indices for a supercell and use V computed for a primitive cell.
        determine_indices_only = .true.
        WRITE( stdout, '(/5x,"Determination of the indices of inter-site couples ...",/)')
     ELSE
        ! DFT+U: exit from this routine
        RETURN
     ENDIF
  ELSE
     determine_indices_only = .false.
     WRITE( stdout, '(/5x,"Post-processing calculation of Hubbard parameters ...",/)')
  ENDIF
  !
  ! Allocate various arrays
  CALL alloc_pp()
  !
  ! Read chi0 and chi from file
  If (.NOT.determine_indices_only) CALL read_chi()
  !
  ! If we merge types of atoms (e.g. Ni_up and Ni_down) 
  ! then we have to keep track of their total magnetization
  ! in order to be able to distinguish them
  !
  CALL equiv_types_and_determine_spin()
  !
  ! Generate the virtual atoms in order 
  ! to mimic a virtual supercell
  CALL gen_virt_atoms()
  !
  ! Compute interatomic distances
  ! between virtual atoms
  CALL atomic_dist()
  !
  IF (determine_indices_only) THEN
     CALL write_uv(.false.)   
     GO TO 15
  ENDIF
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
  ! Calculate Hubbard parameters and write them to file
  CALL calculate_Hubbard_parameters()
  !
15 CONTINUE
  !
  ! Deallocate various arrays
  CALL dealloc_pp()
  !
  CALL stop_clock('hp_postproc')
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
  ALLOCATE ( Hubbard_matrix(nath_scbg, nath_scbg) )
  !
  ! Find and open unit to write info
  unithub = find_free_unit()
  filenamehub = trim(prefix) // ".Hubbard_parameters.dat"
  OPEN(unithub, file = filenamehub, form = 'formatted', status = 'unknown')
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
  DEALLOCATE (Hubbard_matrix)
  !
  CLOSE(unithub)
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

SUBROUTINE equiv_types_and_determine_spin
  !
  ! If we merge types of atoms (e.g. Ni_up and Ni_down) 
  ! then we have to keep track of their total magnetization
  ! (orientation of spin) in order to be able to distinguish them.
  !
  IMPLICIT NONE
  !
  IF ( find_atpert /= 1 ) THEN
     !
     ! If it was requested in the input, re-assign the specific type
     ! to another type (e.g., Ni_down to Ni_up)
     !
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt) .AND. skip_type(nt)) THEN
           ityp_new(na) = equiv_type(nt)
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
        IF ( magn(na) > eps1 ) THEN
           spin(na) = 1
        ELSEIF ( magn(na) < -eps1 ) THEN
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
END SUBROUTINE equiv_types_and_determine_spin

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
     ! Number of atoms in the supercell 
     ! (2*sc_size+1) x (2*sc_size+1) x (2*sc_size+1)
     dimn = num_uc * nat
     !
     ALLOCATE(found(dimn))
     !
     ! Perform a mapping between the nb index of V(na,nb) 
     ! of the nq1xnq2xnq3 supercell with the nc index of V(na,nc)
     ! of the 3x3x3 supercell used in DFT+U+V of the PW code.
     !
     DO na = 1, nath
        found(:) = .FALSE.
        DO nb = 1, nath_sc
           DO nc = 1, dimn
              IF ( ityp_sc0(nb).EQ.ityp_s(nc)                    .AND. &
                   ABS(dist_sc(na,nb)-dist_s(na,nc)).LT.dist_thr .AND. &
                   .NOT.found(nc) ) THEN
                   !
                   ! Mapping of index nb to nc
                   auxindex(na,nb) = nc
                   found(nc) = .TRUE.
                   GO TO 9
                   !
              ENDIF 
           ENDDO
           ! In the case when the mapping is not reached nullify the index
           auxindex(na,nb) = 0
9          CONTINUE         
        ENDDO
     ENDDO   
     !      
     DEALLOCATE(found)
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
                          ABS(dist_sc(nc,nb)-dist_sc(na,nb)).LE.dist_thr 
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
  !
  DO na = 1, nath_sc
     DO nb = 1, nath_sc
        !
        IF ( chi_(na,nb).EQ.0.d0 ) THEN
           !
           DO nc = 1, nath_sc
              IF ( ityp_sc(nc).EQ.ityp_sc(nb) ) THEN
                 !
                 DO nd = 1, nath_sc
                    IF ( ityp_sc(nd).EQ.ityp_sc(na) ) THEN
                       !
                       IF (chi_(nd,nc).NE.0.0d0 .AND. &
                           ABS(dist_sc(nd,nc)-dist_sc(na,nb)).LE.dist_thr  .AND. &
                           spin_sc(na)*spin_sc(nb).EQ.spin_sc(nd)*spin_sc(nc)) THEN
                          !
                          chi_(na,nb) = chi_(nd,nc)
                          GO TO 35
                          !
                       ENDIF
                       !
                       IF (chi_(nc,nd).NE.0.0d0 .AND. &
                           ABS(dist_sc(nc,nd)-dist_sc(na,nb)).LE.dist_thr  .AND. &
                           spin_sc(na)*spin_sc(nb).EQ.spin_sc(nc)*spin_sc(nd)) THEN
                          !
                          chi_(na,nb) = chi_(nc,nd)
                          GO TO 35
                          !
                       ENDIF
                       !
                    ENDIF
                 ENDDO
                 !
              ENDIF
           ENDDO
           !
        ENDIF
        !
35 CONTINUE
        !
     ENDDO
  ENDDO
  !
  ! Check that all elements were found
  !
  IF (ANY(chi_(:,:).EQ.0.0d0)) THEN
     !
     WRITE( stdout, '(/5x,"Existing distances between couples of atoms:"/)')
     DO na = 1, nath_sc
        DO nb = 1, nath_sc
           WRITE( stdout, '(5x,"na=",2x,i4,2x,"nb=",2x,i4,2x,"dist= ",f10.6)') &
              na, nb, dist_sc(na,nb)
        ENDDO
     ENDDO
     !
     DO na = 1, nath_sc
        DO nb = 1, nath_sc
           IF (chi_(na,nb).EQ.0.0d0) WRITE( stdout, '(/5x,"Missing chi element for: na=", &
                2x,i4,2x,"nb=",2x,i4,2x,"dist= ",f10.6/)') na, nb, dist_sc(na,nb)
        ENDDO
     ENDDO
     !
     WRITE( stdout, '(/5x,"Possible solutions:")')
     WRITE( stdout, '(5x, "1. Relax better the structure (in order to have more accurate inter-atomic distances)")')
     WRITE( stdout, '(5x, "2. Increase the value of the parameter dist_thr in the HP input,")')
     WRITE( stdout, '(5x, "   and re-run the postprocessing step by setting compute_hp=.true. in the HP input.")')
     !
     CALL errore ('reconstruct_full_chi', &
            'Reconstruction problem: some chi were not found', 1)
     !
  ENDIF
  !
  ! Symmetrization
  !
  DO na = 1, nath_sc
     DO nb = 1, nath_sc
        chi_(na,nb) = 0.5d0 * ( chi_(na,nb) + chi_(nb,na) )
        chi_(nb,na) = chi_(na,nb)
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

SUBROUTINE calculate_Hubbard_parameters()
  !
  ! This routine calculates Hubbard parameters (U and V) from the inverse matrices
  ! of the susceptibility and writes them to file.
  !
  IMPLICIT NONE
  INTEGER :: nt1, nt2
  CHARACTER(len=2) :: Hubbard_manifold
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=1), EXTERNAL :: l_to_spdf
  !
  ! Calculate the matrix of Hubbard parametres: CHI0^{-1} - CHI^{-1}
  !
  DO na = 1, nath_scbg
     DO nb = 1, nath_scbg
        Hubbard_matrix(na,nb) = inv_chi0bg(na,nb) - inv_chibg(na,nb)
     ENDDO
  ENDDO
  !
  IF ( find_atpert == 1 .AND. ntyp_new > ntyp ) THEN
     WRITE(unithub,'(2x,"Warning: The Hubbard parameters listed below were computed by treating some")')
     WRITE(unithub,'(2x,"         equivalent (by type) Hubbard atoms as if they were non-equivalent")')
     WRITE(unithub,'(2x,"         (when doing the averaging and reconstruction of the response matrices)")')
     WRITE(unithub,'(2x,"         because their unperturbed occupations differ. If you want to disable this")')
     WRITE(unithub,'(2x,"         option then set disable_type_analysis=.true. and redo the post-processing.")')
  ENDIF
  !
  ! Write Hubbard U (i.e. diagonal elements of the Hubbard matrix)
  !
  WRITE(unithub,'(/2x,"=-------------------------------------------------------------------------------=",/)')
  WRITE(unithub,'(33x,"Hubbard U parameters:",/)')
  WRITE(unithub, '(7x,"site n.",2x,"type",2x,"label",2x,"spin",2x,"new_type",2x,"new_label",2x,"manifold",2x,"Hubbard U (eV)")')
  DO na = 1, nat
     nt1 = ityp(na)
     nt2 = ityp_new(na)
     IF ( is_hubbard(nt1) ) THEN
        Hubbard_manifold = TRIM(int_to_char(Hubbard_n(nt1))) // &
                                              l_to_spdf(Hubbard_l(nt1),.FALSE.)
        WRITE(unithub,'(7x,i3,6x,i3,3x,a4,4x,i2,4x,i3,9x,a4,7x,a2,3x,f10.4)') &
                & na, nt1, atm(nt1), spin(na), nt2, atm_new(nt2), Hubbard_manifold, Hubbard_matrix(na,na)
     ENDIF
  ENDDO
  WRITE(unithub,'(/2x,"=-------------------------------------------------------------------------------=",/)')
  !
  ! Write Hubbard parameters to file that can be directly used in the pw.x calculation
  !
  ! Write V
  IF (lda_plus_u_kind==2) CALL write_uv(.true.)
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
     WRITE(unithub,'(/10x,"chi0 matrix :")')
     DO na = 1, nath_scbg
        WRITE(unithub,'(8(x,f11.6))') (chi0bg(na,nb), nb=1,nath_scbg)
        WRITE(unithub,*)
     ENDDO
     !
     ! Write to file chi
     WRITE(unithub,'(/10x,"chi matrix :")')
     DO na = 1, nath_scbg
        WRITE(unithub,'(8(x,f11.6))') (chibg(na,nb), nb=1,nath_scbg)
        WRITE(unithub,*)
     ENDDO
     !
     ! Write to file inv_chi0
     WRITE(unithub,'(/8x,"chi0^{-1} matrix :")')
     DO na = 1, nath_scbg
        WRITE(unithub,'(8(x,f11.6))') (inv_chi0bg(na,nb), nb=1,nath_scbg)
        WRITE(unithub,*)
     ENDDO
     !
     ! Write to file inv_chi
     WRITE(unithub,'(/8x,"chi^{-1} matrix :")')
     DO na = 1, nath_scbg
        WRITE(unithub,'(8(x,f11.6))') (inv_chibg(na,nb), nb=1,nath_scbg)
        WRITE(unithub,*)
     ENDDO
     !
     ! Write the full Hubbard matrix to file
     WRITE(unithub,'(/9x,"Hubbard matrix :")')
     DO na = 1, nath_scbg
        WRITE(unithub,'(8(x,f11.6))') (Hubbard_matrix(na,nb), nb=1,nath_scbg)
        WRITE(unithub,*)
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE calculate_Hubbard_parameters

SUBROUTINE write_uv (lflag)
  !
  ! Write the Hubbard V parameters to file.
  ! V's are written in the order of increasing interatomic distances.
  !
  USE parameters,   ONLY : sc_size
  !
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: dist(:), distord(:)
  INTEGER,  ALLOCATABLE :: typeord0(:), typeord(:), indexord(:)
  INTEGER :: ne, nc_min, ipol, counter, nt, nt2
  REAL(DP) :: auxdist
  LOGICAL :: lflag  ! if .true.  then write V to file
                    ! if .false. then do not write V to file
  CHARACTER(len=2) :: Hubbard_manifold, Hubbard_manifold2
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=1), EXTERNAL :: l_to_spdf
  !
  ! Find and open unit to write info
  unithub2 = find_free_unit()
  filenamehub2 = TRIM("HUBBARD.dat")
  OPEN(unithub2, file = filenamehub2, form = 'formatted', status = 'unknown')
  WRITE(unithub2,'("# Copy this data in the pw.x input file for DFT+Hubbard calculations")')
  WRITE(unithub2,'("HUBBARD {",a,"}")') TRIM(Hubbard_projectors)
  !
  IF (lflag) THEN
     WRITE(unithub,'(/27x,"Hubbard V parameters:")')
  ELSE
     WRITE(unithub,'(/17x,"Indices and distances for inter-site couples:")')
  ENDIF
  WRITE(unithub,'(22x,"(adapted for a supercell",1x,i1,"x",i1,"x",i1,")",/)') &
                       2*sc_size+1, 2*sc_size+1, 2*sc_size+1 
  !
  IF (lflag) THEN
     WRITE(unithub,*) '           Atom 1 ', '    Atom 2 ', '   Distance (Bohr) ', ' Hubbard V (eV)'
     WRITE(unithub,*)
  ELSE
     WRITE(unithub,*) '           Atom 1 ', '    Atom 2 ', '   Distance (Bohr) '
     WRITE(unithub,*)
  ENDIF
  !
  ALLOCATE(dist(nath_sc))
  ALLOCATE(distord(nath_sc))
  ALLOCATE(typeord0(nath_sc))
  ALLOCATE(typeord(nath_sc))
  ALLOCATE(indexord(nath_sc))
  !
  DO na = 1, nath
     !
     nt = ityp_sc0(na)
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
            typeord0(1) = ityp_sc0(nb)
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
           typeord0(nc) = ityp_sc0(nb)
           typeord(nc)  = ityp_sc(nb)
           indexord(nc) = nb
           GO TO 11
        ENDIF
     ENDDO
     !
     nc = nc + 1
     distord(nc)  = distord(1)
     typeord0(nc) = typeord0(1)
     typeord(nc)  = typeord(1)
     indexord(nc) = indexord(1)
     !
     DO nb = 1, nath_sc
        IF ( dist(nb).GT.0.d0 .AND. dist(nb).LT.distord(nc) ) THEN
            distord(1)  = dist(nb)
            typeord0(1) = ityp_sc0(nb)
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
            typeord0(nc) = ityp_sc0(nb)
            typeord(nc)  = ityp_sc(nb)
            indexord(nc) = nb
        ELSEIF ( dist(nb).LE.distord(1) .AND. nb.NE.indexord(1) ) THEN
            nc_min = 2
            nc = nc + 1
            DO nd = nc, nc_min, -1
               distord(nd)  = distord(nd-1)
               typeord0(nd) = typeord0(nd-1)
               typeord(nd)  = typeord(nd-1)
               indexord(nd) = indexord(nd-1)
            ENDDO
            distord(1)  = dist(nb)
            typeord0(1) = ityp_sc0(nb)
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
                        typeord0(ne) = typeord0(ne-1)
                        typeord(ne)  = typeord(ne-1)
                        indexord(ne) = indexord(ne-1)
                     ENDDO
                     distord(nd+1)  = dist(nb)
                     typeord0(nd+1) = ityp_sc0(nb)
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
           IF (lflag) THEN
              ! Print Hubbard U and Hubbard V
              WRITE(unithub,'(11x,i3,x,a4,x,i5,x,a4,2x,f12.6,4x,f10.4)') &
                    na, atm_new(ityp_sc(na)), auxindex(na,indexord(nb)), atm_new(typeord(nb)), &
                    dist_sc(na,indexord(nb)), Hubbard_matrix(na,indexord(nb))
              IF ( nb.LE.(num_neigh+1)              .AND. &     
                   dist_sc(na,indexord(nb)).LE.rmax .AND. & 
                   Hubbard_l(ityp(na)).GE.lmin ) THEN
                   ! write V
                   nt2 = typeord0(nb)
                   Hubbard_manifold = TRIM(int_to_char(Hubbard_n(nt))) // &
                                              l_to_spdf(Hubbard_l(nt),.FALSE.)
                   Hubbard_manifold2 = TRIM(int_to_char(Hubbard_n(nt2))) // &
                                              l_to_spdf(Hubbard_l(nt2),.FALSE.) 
                   WRITE(unithub2,'("V",2x,a5,"-",a2,2x,a5,"-",a2,1x,i4,1x,i5,1x,f8.4)') &
                            ADJUSTR(atm_new(nt)),  Hubbard_manifold,  &
                            ADJUSTR(atm_new(nt2)), Hubbard_manifold2, &
                            na, auxindex(na,indexord(nb)),         &
                            Hubbard_matrix(na,indexord(nb))
              ENDIF
           ELSE
              ! Print couples only
              WRITE(unithub,'(11x,i3,x,a4,x,i5,x,a4,2x,f12.6)') &
                    na, atm_new(ityp_sc(na)), auxindex(na,indexord(nb)), atm_new(typeord(nb)), &
                    dist_sc(na,indexord(nb))
              IF ( nb.LE.(num_neigh+1)              .AND. &
                   dist_sc(na,indexord(nb)).LE.rmax .AND. &
                   Hubbard_l(ityp(na)).GE.lmin) THEN
                   !
                   nt2 = typeord0(nb)
                   Hubbard_manifold = TRIM(int_to_char(Hubbard_n(nt))) // &
                                              l_to_spdf(Hubbard_l(nt),.FALSE.)
                   Hubbard_manifold2 = TRIM(int_to_char(Hubbard_n(nt2))) // &
                                              l_to_spdf(Hubbard_l(nt2),.FALSE.)
                   WRITE(unithub2,'("V",2x,a,"-",a,2x,a,"-",a,1x,i4,1x,i5,3x,"V-not-set")') &
                            TRIM(atm_new(nt)),  Hubbard_manifold,  &
                            TRIM(atm_new(nt2)), Hubbard_manifold2, &
                            na, auxindex(na,indexord(nb))
              ENDIF
           ENDIF
        ENDIF
     ENDDO 
     WRITE(unithub,'()')
     !
  ENDDO ! na
  !
  WRITE(unithub,'(/2x,"=-------------------------------------------------------------------=",/)')
  !
  IF ((counter-1).GE.num_neigh) &
    WRITE(unithub,'(2x,i3,2x,"nearest neighbors of every atom are written to HUBBARD.dat")') num_neigh
  !
  DEALLOCATE(dist)
  DEALLOCATE(distord)
  DEALLOCATE(typeord0)
  DEALLOCATE(typeord)
  DEALLOCATE(indexord)
  CLOSE(unithub2)
  !
  RETURN
  !
END SUBROUTINE write_uv

END SUBROUTINE hp_postproc
