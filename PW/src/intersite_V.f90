!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE alloc_neighborhood()
  !---------------------------------------------------------------------
  !
  ! This routine allocates and assigns values to the neighborhood, at_sc and sc_at 
  !
  USE symm_base,       ONLY : nsym
  USE io_global,       ONLY : stdout
  USE ions_base,       ONLY : nat, tau, ityp
  USE cell_base,       ONLY : at, alat
  USE kinds,           ONLY : DP
  USE constants,       ONLY : rytoev
  USE parameters,      ONLY : sc_size
  USE control_flags,   ONLY : dfpt_hub
  USE ldaU,            ONLY : num_uc, max_num_neighbors, neighood, &
                              at_sc, sc_at, Hubbard_V, is_hubbard, is_hubbard_back, &
                              dist_s, ityp_s, deallocate_at_center_type, eps_dist
  !
  IMPLICIT NONE 
  !
  ! Local variables
  !
  REAL(DP), ALLOCATABLE :: V(:,:,:), tau_sc(:,:)
  REAL(DP) :: daux, distm(7)
  REAL(DP), PARAMETER :: eps1 = 1.d-20  ! threshold for Hubbard_V
  INTEGER :: i, j, k, l, ii, jj, nx, ny, nz, dimn, &
             viz, atom, nb1, nb2, isym, l1, l2, l3, na, nb
  !
  CALL start_clock( 'alloc_neigh' )
  !
  ! Number of atoms in the supercell
  dimn = num_uc * nat       
  !
  ALLOCATE (V(nat,dimn,3))
  ALLOCATE (tau_sc(3,dimn))
  IF (.NOT.ALLOCATED(ityp_s))   ALLOCATE (ityp_s(dimn))
  IF (.NOT.ALLOCATED(dist_s))   ALLOCATE (dist_s(nat,dimn))
  !
13 CONTINUE
  !
  V(:,:,:) = 0.d0
  !
  ! Build V from Hubbard_V (input plus equivalent neighborhood)
  !
  DO na = 1, nat
     DO nb = 1, dimn
        IF (ABS(Hubbard_V(na,nb,1)) /= 0.0_dp .OR. &
            ABS(Hubbard_V(na,nb,2)) /= 0.0_dp .OR. &
            ABS(Hubbard_V(na,nb,3)) /= 0.0_dp ) THEN
            !
            V(na,nb,1:3) = Hubbard_V(na,nb,1:3)
            ! 
        ENDIF
     ENDDO
  ENDDO
  !
  IF (.NOT.ALLOCATED(neighood)) ALLOCATE (neighood(nat))
  IF (.NOT.ALLOCATED(at_sc))    ALLOCATE (at_sc(dimn))
  IF (.NOT.ALLOCATED(sc_at))    ALLOCATE (sc_at(1:nat,-sc_size:sc_size,-sc_size:sc_size,-sc_size:sc_size))
  ! 
  ! Initialization of various quantities in the supercell
  ! First initialization in the real unit cell
  !
  atom = 0
  !
  DO na = 1, nat
     atom = atom + 1
     at_sc(atom)%n(1) = 0
     at_sc(atom)%n(2) = 0
     at_sc(atom)%n(3) = 0
     at_sc(atom)%at = na        ! it is equivalent to itself
     sc_at(na,0,0,0) = atom
     tau_sc(:,atom) = tau(:,na)
     ityp_s(atom) = ityp(na)
  ENDDO
  !
  ! Now initialization in the virtual cells of the supercell
  ! (in the Cartesian framework)
  !
  DO nx = -sc_size, sc_size
     DO ny = -sc_size, sc_size
        DO nz = -sc_size, sc_size
           !
           IF ( nx.NE.0 .OR. ny.NE.0 .OR. nz.NE.0 ) THEN
              !  
              DO na = 1, nat
                 !
                 atom = atom + 1
                 at_sc(atom)%n(1) = nx
                 at_sc(atom)%n(2) = ny
                 at_sc(atom)%n(3) = nz
                 at_sc(atom)%at = na
                 sc_at(na,nx,ny,nz) = atom
                 tau_sc(:,atom) = tau(:,na) + DBLE(nx)*at(:,1) &
                                            + DBLE(ny)*at(:,2) &
                                            + DBLE(nz)*at(:,3) 
                 ityp_s(atom) = ityp(na)
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
     ENDDO
  ENDDO
  !
  ! Initialize is_hubbard and is_hubbard_back
  !
  DO na = 1, nat
     ! 
     DO nb = 1, dimn
        !
        IF (ABS(V(na,nb,1)) /= 0.0_dp .OR. &
            ABS(V(na,nb,2)) /= 0.0_dp .OR. &
            ABS(V(na,nb,3)) /= 0.0_dp ) THEN
            !
            IF (ABS(V(na,nb,1)) /= 0.0_dp) THEN
               is_hubbard(ityp(na))   = .TRUE.
               is_hubbard(ityp_s(nb)) = .TRUE.
            ENDIF
            !
            IF (ABS(V(na,nb,2)) /= 0.0_dp) THEN
               is_hubbard(ityp(na))        = .TRUE.
               is_hubbard_back(ityp_s(nb)) = .TRUE.
            ENDIF
            !
            IF (ABS(V(na,nb,3)) /= 0.0_dp) THEN
               is_hubbard_back(ityp(na))   = .TRUE.
               is_hubbard_back(ityp_s(nb)) = .TRUE.
            ENDIF
            !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ! Compute the distances between atom na and all other atoms 
  ! in the supercell
  !
  dist_s(:,:) = 0.d0
  !
  DO na = 1, nat
     !
     DO nb = 1, dimn
        !
        DO k = 1, 3
           dist_s(na,nb) = dist_s(na,nb) + (tau_sc(k,na) - tau_sc(k,nb))**2
        ENDDO
        dist_s(na,nb) = DSQRT(dist_s(na,nb)) * alat
        !
        ! Uncomment the following lines if you want that the
        ! distances are printed in the output file
        !
        !IF (.NOT.dfpt_hub) WRITE(stdout,*) na, nb, ityp_s(na), ityp_s(nb), dist_s(na,nb)
        !
     ENDDO
     !
  ENDDO
  !
  ! Order distances
  !
  IF (.NOT.dfpt_hub) WRITE(stdout,'(5x,"First shells distances (in Bohr):")')
  !
  ! Here we determine 7 smallest distances from the first atom
  ! to its neighbors (in the increasing order). If larger distances
  ! are needed, then 7 has to be changed to a larger number.  
  !
  l = 0
  distm(:) = 10000.d0
  !
  DO na = 1, 7
     DO nb = 1, dimn
        DO k = 1, l
           IF ( ABS(dist_s(1,nb)-distm(k)).LE.1.d-6 ) GO TO 49
        ENDDO
        distm(l+1) = MIN(distm(l+1), dist_s(1,nb))
49 CONTINUE
     ENDDO
     l = l + 1
     IF (.NOT.dfpt_hub) WRITE(stdout,'(5x,"shell:",2x,i2,2x,f10.6)') l, distm(l)
  ENDDO
  !
  ! Find the missing Hubbard_V(i,j,2) i.e. standard-background interactions. 
  ! This loop is needed only when the background is included.
  !
  IF ( ANY(V(:,:,2).GE.eps1) ) THEN 
     !
     DO i = 1, nat
        DO j = 1, dimn
           DO l = 1, nat
              DO k = 1, dimn
                 !
                 IF ( ityp_s(k).EQ.ityp_s(j)                   .AND. &
                      ityp(l).EQ.ityp(i)                       .AND. &
                      ABS(dist_s(l,k)-dist_s(i,j)).LE.eps_dist .AND. &
                      ABS(V(l,k,2)).GE.eps1 )                  THEN 
                      !
                      IF (Hubbard_V(i,j,2).EQ.0.d0) Hubbard_V(i,j,2) = V(l,k,2)
                      GO TO 22
                      !
                 ENDIF
                 !
              ENDDO ! k
           ENDDO ! l    
22         CONTINUE
        ENDDO ! j
     ENDDO ! i
     !
     DO i = 1, nat
        DO j = 1, dimn
           DO l = 1, nat
              DO k = 1, dimn
                 !
                 IF ( ityp_s(k).EQ.ityp(i)                      .AND. &
                      ityp(l).EQ.ityp_s(j)                      .AND. &
                      ABS(dist_s(l,k)-dist_s(i,j)).LE.eps_dist  .AND. &
                      ABS(V(l,k,2)).GE.eps1 )                   THEN 
                      !
                      IF (Hubbard_V(i,j,4).EQ.0.d0) Hubbard_V(i,j,4) = V(l,k,2)
                      GO TO 23
                      !
                 ENDIF
                 !
              ENDDO ! k
           ENDDO ! l   
23         CONTINUE
        ENDDO ! j
     ENDDO ! i
  ENDIF
  !
  ! Find the missing Hubbard_V(i,j,1) i.e. standard-standard interactions.
  !
  DO i = 1, nat
     DO j = 1, dimn
        DO l = 1, nat
           DO k = 1, dimn
              !
              IF ( ityp_s(k).EQ.ityp_s(j)                    .AND. &
                   ityp(l).EQ.ityp(i)                        .AND. &
                   ABS(dist_s(l,k)-dist_s(i,j)).LE.eps_dist  .AND. &
                  ( ABS(V(l,k,1)).GE.eps1 ) )                THEN
                  !
                  IF (Hubbard_V(i,j,1).EQ.0.d0) Hubbard_V(i,j,1) = V(l,k,1)
                  GO TO 32
                  !
              ENDIF
              !
              IF ( ityp_s(k).EQ.ityp(i)                      .AND. &
                   ityp(l).EQ.ityp_s(j)                      .AND. &
                   ABS(dist_s(l,k)-dist_s(i,j)).LE.eps_dist  .AND. &
                  ( ABS(V(l,k,1)).GE.eps1 ) )                THEN
                  !
                  IF (Hubbard_V(i,j,1).EQ.0.d0) Hubbard_V(i,j,1) = V(l,k,1)
                  GO TO 32
                  !
              ENDIF
              !
           ENDDO
        ENDDO      
32      CONTINUE
     ENDDO
  ENDDO
  !
  ! Find the missing Hubbard_V(i,j,3) i.e. background-background interactions.
  ! This loop is needed only when the background is included.
  !
  IF ( ANY(V(:,:,3).GE.eps1) ) THEN 
     DO i = 1, nat
        DO j = 1, dimn
           DO l = 1, nat
              DO k = 1, dimn
                 !
                 IF ( ityp_s(k).EQ.ityp_s(j)                   .AND. &
                      ityp(l).EQ.ityp(i)                       .AND. &
                      ABS(dist_s(l,k)-dist_s(i,j)).LE.eps_dist .AND. &
                     ( ABS(V(l,k,3)).GE.eps1 ) )               THEN
                     !
                     IF (Hubbard_V(i,j,3).EQ.0.d0) Hubbard_V(i,j,3) = V(l,k,3)
                     GO TO 33
                     !
                 ENDIF
                 !
                 IF ( ityp_s(k).EQ.ityp(i)                     .AND. &
                      ityp(l).EQ.ityp_s(j)                     .AND. &
                      ABS(dist_s(l,k)-dist_s(i,j)).LE.eps_dist .AND. &
                     ( ABS(V(l,k,3)).GE.eps1) )                THEN
                     !
                     IF (Hubbard_V(i,j,3).EQ.0.d0) Hubbard_V(i,j,3) = V(l,k,3)
                     GO TO 33
                     !
                 ENDIF
                 !
              ENDDO ! k
           ENDDO ! l     
33         CONTINUE
        ENDDO ! j
     ENDDO ! i
  ENDIF
  !
  IF (.NOT.dfpt_hub) WRITE(stdout,'(/5x,"i",4x,"j",2x,"dist (Bohr)", &
                                & 7x,"stan-stan",x,"stan-bac",x,"bac-bac",x,"bac-stan")')
  ! stan-stan = standard-standard
  ! stan-bac  = standard-background
  ! bac-bac   = background-background
  ! bac-stan  = background-standard
  !
  ! Determine how many neighbors there are and what are their indices.
  !
  max_num_neighbors = 0
  !
  DO i = 1, nat
     !
     neighood(i)%num_neigh = 0
     !
     DO j = 1, dimn
        IF ( ABS(Hubbard_V(i,j,1)).GE.eps1 .OR. &
             ABS(Hubbard_V(i,j,2)).GE.eps1 .OR. &
             ABS(Hubbard_V(i,j,3)).GE.eps1 .OR. &
             i==j ) THEN
             ! 
             ! Diagonal term i=j needed in init_nsg, nsg_adj and mix_rho
             !
             neighood(i)%num_neigh = neighood(i)%num_neigh + 1
             !
        ENDIF
     ENDDO
     !
     IF (max_num_neighbors .LT. neighood(i)%num_neigh) THEN
        max_num_neighbors = neighood(i)%num_neigh
     ENDIF
     !
     ALLOCATE (neighood(i)%neigh(1:neighood(i)%num_neigh))
     ! 
     viz = 0
     !
     DO j = 1, dimn
        IF ( ABS(Hubbard_V(i,j,1)).GE.eps1 .OR. &
             ABS(Hubbard_V(i,j,2)).GE.eps1 .OR. &
             ABS(Hubbard_V(i,j,3)).GE.eps1 .OR. &
             i==j ) THEN
             !
             ! Diagonal term i=j needed in init_nsg, nsg_adj and mix_rho
             !
             IF (.NOT.dfpt_hub) WRITE(stdout,'(2x,i4,x,i4,x,f12.8,x,a5,4(x,f8.4))') &
                        i, j, dist_s(i,j), ' V = ', (Hubbard_V(i,j,k)*rytoev, k=1,4)
             !
             viz = viz + 1
             !
             neighood(i)%neigh(viz) = j
             ! 
        ENDIF
     ENDDO
     !
  ENDDO
  !
  nb2 = 0
  !
  DO i = 1, nat
     !
     DO nb1 = 1, neighood(i)%num_neigh
        !
        j = neighood(i)%neigh(nb1)
        !
        DO isym = 1, nsym
           !
           IF ( ABS(Hubbard_V(i,j,1)) /= 0.0_dp .OR. &
                ABS(Hubbard_V(i,j,2)) /= 0.0_dp .OR. &
                ABS(Hubbard_V(i,j,3)) /= 0.0_dp ) THEN
                !
                CALL symonpair(i,j,isym,ii,jj)
                !
                IF (ABS(dist_s(i,j)-dist_s(ii,jj)).GT.eps_dist) THEN
                   WRITE(stdout,'(/2x,"Different distances between couples!")')
                   WRITE(stdout,'(2x,"Original couple:      ",2x,i4,2x,i4,2x,"dist =",1x,f8.4,1x,"(Bohr)")') &
                   i, j, dist_s(i,j)
                   WRITE(stdout,'(2x,"New additional couple:",2x,i4,2x,i4,2x,"dist =",1x,f8.4,1x,"(Bohr)")') &
                   ii, jj, dist_s(ii,jj)
                   CALL errore('alloc_neighborhood', 'Probably a larger sc_size is needed',1) 
                ENDIF
                !
                IF ( ABS(Hubbard_V(ii,jj,1)).LT.eps1 .AND. &
                     ABS(Hubbard_V(ii,jj,2)).LT.eps1 .AND. &
                     ABS(Hubbard_V(ii,jj,3)).LT.eps1 ) THEN
                     !
                     Hubbard_V(ii,jj,:) = Hubbard_V(i,j,:)
                     WRITE(stdout,'(2x,"Found a new Hubbard_V element from the symmetry analysis!")')
                     WRITE(stdout,'(2x,"Original couple:      ",2x,i4,2x,i4,2x,"dist =",1x,f8.4,1x,"(Bohr)")') &
                     i, j, dist_s(i,j)
                     WRITE(stdout,'(2x,"New additional couple:",2x,i4,2x,i4,2x,"dist =",1x,f8.4,1x,"(Bohr)")') &
                     ii, jj, dist_s(ii,jj)
                     !
                     nb2 = nb2 + 1
                     !
                ELSEIF ( ABS(Hubbard_V(ii,jj,1)-Hubbard_V(i,j,1)) + &
                         ABS(Hubbard_V(ii,jj,2)-Hubbard_V(i,j,2)) + &
                         ABS(Hubbard_V(ii,jj,3)-Hubbard_V(i,j,3)).GT.eps1 ) THEN
                     !
                     ! This V element has been set to something different before
                     ! possibly, wrong assignement/symm ops...
                     !
                     ! Some more output for the standard-standard case
                     !
                     IF (ABS(Hubbard_V(ii,jj,1)-Hubbard_V(i,j,1)).GT.eps1) THEN
                        ! 
                        WRITE(stdout,'(/4x,"WARNING! Equivalent couples with some small variations in V :")')
                        WRITE(stdout,'(4x,"V (",i3,",",i5,") =",1x,f10.6)') &
                                          i, j, Hubbard_V(i,j,1)*rytoev
                        WRITE(stdout,'(4x,"V (",i3,",",i5,") =",1x,f10.6)') &
                                          ii, jj, Hubbard_V(ii,jj,1)*rytoev
                        WRITE(stdout,'(4x,"delta(V) =",1x,f10.6)') &
                                          ABS(Hubbard_V(ii,jj,1)-Hubbard_V(i,j,1))*rytoev
                        !
                        IF ( ABS(Hubbard_V(ii,jj,1)-Hubbard_V(i,j,1))*rytoev .GE. 5.0d-3 ) THEN
                           WRITE(stdout,*) ' i,  j, V:',  i,  j, Hubbard_V(i,j,1)*rytoev
                           WRITE(stdout,*) 'ii, jj, V:', ii, jj, Hubbard_V(ii,jj,1)*rytoev
                           CALL errore('alloc_neighborhood', &
                               & 'Problem with Hubbard_V after a symmetry analysis',1)
                        ELSE
                           WRITE(stdout,'(4x,"Setting them to be equal...")')
                           Hubbard_V(ii,jj,1) = Hubbard_V(i,j,1)   
                        ENDIF      
                        !          
                     ENDIF
                     !
                     IF ( ABS(Hubbard_V(ii,jj,2)-Hubbard_V(i,j,2))*rytoev .GE. 5.0d-3   .OR. &  
                          ABS(Hubbard_V(ii,jj,3)-Hubbard_V(i,j,3))*rytoev .GE. 5.0d-3 ) THEN
                          CALL errore('alloc_neighborhood', &
                             & 'Problem with Hubbard_V after a symmetry analysis',1)
                     ENDIF
                     !
                ENDIF
                ! 
           ENDIF
           !
        ENDDO ! isym
        !
     ENDDO ! nb1
     !
  ENDDO ! i
  !
  IF ( nb2.GT.0 ) THEN
     DO na = 1, nat
        CALL deallocate_at_center_type ( neighood(na) )
     ENDDO
     WRITE(stdout,'(2x,"Check all Hubbard_V again based on the distances...")')
     GO TO 13
  ENDIF
  !
  DEALLOCATE (V)
  DEALLOCATE (tau_sc)
  !
  CALL stop_clock( 'alloc_neigh' )
  !
  RETURN
  !
END SUBROUTINE alloc_neighborhood
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION type_interaction (na1, m1, na2, m2)
  !-----------------------------------------------------------------------
  !
  ! This function determines type of the interaction
  !
  USE ions_base,   ONLY : ityp
  USE ldaU,        ONLY : Hubbard_l         
  !
  IMPLICIT NONE
  INTEGER :: na1, na2, & ! atom number (in the unit cell)
             m1, m2      ! magnetic quantum number
  INTEGER :: type_interaction, nt1, nt2
  !
  nt1 = ityp(na1)
  nt2 = ityp(na2)
  !
  IF ( m1 .LE. 2*Hubbard_l(nt1)+1 .AND. &
       m2 .LE. 2*Hubbard_l(nt2)+1) THEN
       !
       ! standard-standard term
       type_interaction = 1
       !
  ELSEIF ( m1 .LE. 2*Hubbard_l(nt1)+1  .AND. &
           m2 .GT. 2*Hubbard_l(nt2)+1 ) THEN
       !    
       ! standard-background term
       type_interaction = 2
       !
  ELSEIF ( m1 .GT. 2*Hubbard_l(nt1)+1   .AND. &
           m2 .GT. 2*Hubbard_l(nt2)+1 ) THEN
       !
       ! background-background term
       type_interaction = 3
       !
  ELSE
       !
       ! background-standard term
       type_interaction = 4 
       !
  ENDIF
  !
  RETURN
  !
END FUNCTION type_interaction
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE symonpair (at1, at2, p_sym, rat1, rat2)
  !-----------------------------------------------------------------------
  !
  !  This subroutine acts the symmetry operatorion number p_sym over the 
  !  atomic position of atom at1 in the reference unit cell and over the
  !  atomic position of atom at2 in the supercell. 
  !  Then it identifies, in the reference unit cell,
  !  the equivalent atom to the first 'rotation' above, storing the result
  !  in rat1 and the equivalent atom to the second 'rotation', storing the
  !  result in rat2 (temporarily).
  !  The difference dx = O(p_sym)at1 - rat1 is also subtracted from 
  !  O(p_sym)at2 to find the final position of the second atom.
  !  To determine in which unit cell the second atom ended up, we
  !  compute the difference between its final position and the position
  !  of its equivalent in the reference unit cell.
  !  Finally, we store in rat2 the position, in the full list of atoms in
  !  the supercell, the rotated and translated second atom corresponds to.
  !     
  USE symm_base,       ONLY : s,ft,nsym,irt
  USE ions_base,       ONLY : nat,ityp
  USE cell_base,       ONLY : bg
  USE fft_base,        ONLY : dfftp
  USE ldaU,            ONLY : atom_pos, at_sc, sc_at, num_uc
  USE kinds
  USE io_global,       ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: at1, at2, p_sym
  INTEGER, INTENT(OUT) :: rat1, rat2
  ! rat1 = atom equiv to at1 after 'rotation' and eventual integer translation
  ! rat2 = atom equiv to at2 after 'rotation' and the same translation above.
  !
  ! Local variables
  !
  INTEGER :: i, j, at, dr(3), equiv_2, dimn
  ! dr(1), dr(2), dr(3) = location of the unit cell where at2 goes after sym. operation
  !
  REAL(DP) :: diff, x2(3), r1(3), r2(3), dx(3), ss(3,3)
  REAL(DP), PARAMETER :: eps = 5.d-6
  !
  ! Number of atoms in the supercell
  dimn = num_uc * nat
  !
  ! Convert the symmetry matrix from integer to real type
  ! for a given symmetry operation p_sym
  !
  DO i = 1, 3
     DO j = 1, 3
        ss(i,j) = DBLE(s(i,j,p_sym))
     ENDDO
  ENDDO
  !
  ! Determine the coordinates of at2 in the supercell
  !
  DO i = 1, 3
     x2(i) = atom_pos(at_sc(at2)%at,i) + DBLE(at_sc(at2)%n(i))
  ENDDO
  !
  ! Apply the symmetry operation: S*at - f
  !
  DO i = 1, 3
     !
     r1(i) = 0.d0
     r2(i) = 0.d0
     !
     ! Rotate at1 and at2
     ! 
     DO j = 1, 3
        r1(i) = r1(i) + atom_pos(at1,j) * ss(j,i)
        r2(i) = r2(i) + x2(j) * ss(j,i) 
     ENDDO
     !
     ! Subtract vectors of fractional translations
     ! 
     r1(i) = r1(i) - ft(i, p_sym)
     r2(i) = r2(i) - ft(i, p_sym)
     !
  ENDDO
  !
  ! Atom equivalent to r2 : irt(p_sym,at2)
  !
  equiv_2 = at_sc(at2)%at
  !
  diff = 1.d0
  at = 1
  !
  DO WHILE ( at.LE.nat .AND. diff > eps )
     !
     IF ( ityp(at).EQ.ityp(equiv_2) ) THEN
        diff = 0.d0
        DO i = 1, 3
           dx(i) = r2(i) - atom_pos(at,i)
           diff = diff + ABS(dx(i)-NINT(dx(i)))
        ENDDO
     ELSE
        diff = 1.d0
     ENDIF
     !
     at = at + 1
     ! 
  ENDDO
  !
  IF ( diff > eps ) THEN
     WRITE(stdout,*) "diff > 0, diff= ", diff, "at1= ", at1, "at2= ", at2
     CALL errore('symonpair', 'No atom equivalent to r2', 1)
  ENDIF
  !
  rat2 = at - 1
  !
  ! Atom equivalent to r1 : irt(p_sym,at1)
  !
  diff = 1.d0
  at = 1
  !
  DO WHILE ( at.LE.nat .AND. diff > eps )
     !
     IF ( ityp(at).EQ.ityp(at1) ) THEN
        diff = 0.d0
        DO i = 1, 3
           dx(i) = r1(i) - atom_pos(at,i)
           diff = diff + ABS(dx(i)-NINT(dx(i)))
        ENDDO
     ELSE
        diff = 1.d0
     ENDIF
     !
     at = at + 1
     !
  ENDDO
  !    
  IF ( diff > eps ) THEN
    WRITE(stdout,*) "diff > 0, diff= ", diff, "at1= ", at1, "at2= ", at2 
    CALL errore('symonpair', 'No atom equivalent to r1', 1)
  ENDIF
  !
  rat1 = at - 1
  !
  IF (rat1 > nat .OR. rat1 < 1) THEN
     WRITE(stdout,*) "Index of the first rotated atom=", rat1
     WRITE(stdout,*) "Number of atoms in the original unit cell=", nat
     CALL errore('symonpair', 'Out of bounds', 1)
  ENDIF
  !
  DO i = 1, 3 
     r2(i) = r2(i) - dx(i)
     dr(i) = NINT(r2(i) - atom_pos(rat2,i))
  ENDDO
  !
  rat2 = sc_at(rat2, dr(1), dr(2), dr(3))
  !
  IF (rat2 > dimn) THEN
     WRITE(stdout,*) "Index of the second rotated atom=", rat2
     WRITE(stdout,*) "Number of atoms in the supercell=", dimn
     WRITE(stdout,*) "Probably a larger sc_size is needed"
     CALL errore('symonpair', 'Out of bounds', 1) 
  ELSEIF (rat2 < 1) THEN
     WRITE(stdout,*) "Index of the second rotated atom=", rat2
     CALL errore('symonpair', 'Out of bounds', 1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE symonpair
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION find_viz (center, atom)
  !------------------------------------------------------------------------
  !
  ! This function finds by which number the atom 'atom' (in the full supercell 
  ! list) is indexed amongst the neighbour atoms of the atom 'center' (in the 
  ! reference unit cell)
  !
  USE ldaU,   ONLY : neighood
  !
  IMPLICIT NONE
  INTEGER :: center, atom
  INTEGER :: i, find_viz
  !
  i = 1
  !
  DO WHILE ( i.LE.neighood(center)%num_neigh .AND. &
             neighood(center)%neigh(i).NE.atom )
     i = i + 1
  ENDDO
  !
  IF ( i.LE.neighood(center)%num_neigh ) THEN
     find_viz = i
  ELSE
     find_viz = -1
     WRITE(*,*) "find_viz(",center,atom,")", neighood(center)%num_neigh, i
     CALL errore('find_viz', 'atom is not neighbour of center', 1)
  ENDIF
  ! 
  RETURN
  !
END FUNCTION find_viz
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE phase_factor (ik)
  !------------------------------------------------------------------------
  !
  ! This routine computes the phases for each atom at a given k point 'ik'
  !
  USE kinds,           ONLY : DP
  USE klist,           ONLY : xk
  USE ions_base,       ONLY : nat, ityp
  USE cell_base,       ONLY : at, tpiba
  USE ldaU,            ONLY : at_sc, ldim_u, neighood, phase_fac, num_uc
  USE constants,       ONLY : tpi
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik
  !
  ! Local variables
  ! 
  INTEGER :: na1, viz, na2, nt1, nt2, i, j, dimn
  REAL(DP) :: angle, sum_j
  !
  ! Number of atoms in the supercell 
  dimn = num_uc * nat
  !
  IF (.NOT.ALLOCATED(phase_fac)) ALLOCATE(phase_fac(dimn))
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( ldim_u(nt1).GT.0 ) THEN
        !
        DO viz = 1, neighood(na1)%num_neigh
           !
           na2 = neighood(na1)%neigh(viz)
           angle = 0.d0
           !  
           DO i = 1, 3
              sum_j = 0.d0
              DO j = 1, 3
                 sum_j = sum_j + DBLE(at_sc(na2)%n(j)) * at(i,j)
              ENDDO
              angle = angle + sum_j * xk(i,ik)
           ENDDO
           !
           angle = tpi * angle
           phase_fac(na2) = CMPLX(COS(angle),SIN(angle),kind=DP)
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE phase_factor
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE alloc_atom_pos()
  !-----------------------------------------------------------------------
  !
  ! This routine computes the atomic position in the primitive basis coordinates
  !
  USE ions_base,        ONLY : nat, tau
  USE cell_base,        ONLY : bg
  USE ldaU,             ONLY : atom_pos
  ! 
  IMPLICIT NONE
  INTEGER :: na, ipol
  ! counters on atoms and coordinates
  !
  ALLOCATE(atom_pos(nat,3))
  !
  DO na = 1, nat
     DO ipol = 1, 3
        atom_pos(na,ipol) = bg(1,ipol)*tau(1,na) + &
                            bg(2,ipol)*tau(2,na) + &
                            bg(3,ipol)*tau(3,na)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE alloc_atom_pos
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE write_V
  !-----------------------------------------------------------------------
  !
  ! This routine writes Hubbard_V to file
  !
  USE io_files,       ONLY : seqopn
  USE io_global,      ONLY : ionode
  USE ldaU,           ONLY : Hubbard_V
  !
  IMPLICIT NONE
  INTEGER :: i, j, k, iunit
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit
  !
  iunit = find_free_unit()
  ! 
  IF (ionode) THEN
     CALL seqopn( iunit, 'HubbardV.txt', 'FORMATTED', exst )
     DO i = 1, SIZE(Hubbard_V,1)
        DO j = 1, SIZE(Hubbard_V,2)
           DO k = 1, SIZE(Hubbard_V,3)
              IF (Hubbard_V(i,j,k) > 1.d-20) WRITE(iunit,*) i, j, k, Hubbard_V(i,j,k)
           ENDDO
        ENDDO
     ENDDO
     CLOSE(UNIT=iunit, STATUS='KEEP')
  ENDIF
  !
  RETURN
  ! 
END SUBROUTINE write_V
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE read_V
  !-----------------------------------------------------------------------
  !
  ! This routine reads Hubbard_V from file
  !
  USE kinds,          ONLY : DP
  USE io_files,       ONLY : seqopn
  USE io_global,      ONLY : ionode, ionode_id
  USE ldaU,           ONLY : Hubbard_V
  USE mp_images,      ONLY : intra_image_comm
  USE mp,             ONLY : mp_bcast
  !
  IMPLICIT NONE
  REAL(DP) :: V
  INTEGER :: i, j, k, iunit, ierr
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit
  !
  iunit = find_free_unit()
  ! 
  Hubbard_V(:,:,:) = 0.0d0
  !
  IF (ionode) THEN
     CALL seqopn( iunit, 'HubbardV.txt', 'FORMATTED', exst )
     IF (exst) THEN 
10      READ(iunit,*,END=11,IOSTAT=ierr) i, j, k, V
        IF ( ierr/=0 ) THEN
           CALL mp_bcast( ierr, ionode_id, intra_image_comm )
           CALL errore('read_V', 'Reading Hubbard_V', 1)
        ENDIF
        Hubbard_V(i,j,k) = V
        GO TO 10
     ELSE
        CALL errore('read_V','File HubbardV.txt was not found...',1)
     ENDIF
11   CLOSE( UNIT=iunit, STATUS='KEEP' )
  ENDIF
  ! Broadcast Hubbard_V across all processors
  CALL mp_bcast( Hubbard_V, ionode_id, intra_image_comm )
  !
  RETURN
  ! 
END SUBROUTINE read_V
!-------------------------------------------------------------------------
