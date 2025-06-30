!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE init_nsg
   !-----------------------------------------------------------------------
   !
   ! This routine computes the starting ns (for DFT+U+V calculation) filling
   ! up the Hubbard manifold (we are only interested in the on-site potential 
   ! for the moment) according to the Hund's rule (valid for the isolated atoms 
   ! on which starting potential is built), and to the starting_magnetization:
   ! majority spin levels are populated first, then the remaining electrons
   ! are equally distributed among the minority spin states
   !
   USE kinds,       ONLY : DP
   USE ions_base,   ONLY : nat, ityp
   USE uspp_param,  ONLY : upf
   USE lsda_mod,    ONLY : nspin, starting_magnetization
   USE ldaU,        ONLY : Hubbard_l, Hubbard_l2, Hubbard_l3, hubbard_occ, &
                           nsg, nsgnew, ldim_u, backall, is_hubbard, is_hubbard_back
   USE noncollin_module, ONLY : angle1, angle2, noncolin
#if defined (__OSCDFT)
   USE plugin_flags,     ONLY : use_oscdft
   USE oscdft_base,      ONLY : oscdft_ctx
#endif   
   !
   IMPLICIT NONE
   REAL(DP) :: totoc, totoc_b, cosin 
   COMPLEX(DP) :: esin, n, m, ns(4)
   INTEGER :: ldim, ldim2, na, nt, is, m1, majs, mins, viz
   LOGICAL :: nm        ! true if the atom is non-magnetic
   INTEGER, EXTERNAL :: find_viz
   !
   nsg(:,:,:,:,:) = (0.d0, 0.d0)
   !
#if defined (__OSCDFT)
   IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2) .AND. &
      .NOT.oscdft_ctx%inp%constraint_diag) THEN
      CALL oscdft_nsg(1)
      nsg = nsgnew
   ENDIF
#endif
   ! 
   DO na = 1, nat
      ! 
      ! The index of atom 'na' in the neighbors list
      !
      viz = find_viz(na,na) 
      !
      nt = ityp(na)
      !
      IF ( is_hubbard(nt) ) THEN 
         !
#if defined (__OSCDFT)
      IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2) .AND. &
         .NOT.oscdft_ctx%inp%constraint_diag) THEN
         IF (ANY(DBLE(nsg(:,:,:,na,:))/=0.d0)) GO TO 7
      ENDIF
#endif
         !
         ldim = 2*Hubbard_l(nt)+1
         nm = .TRUE.
         !
         totoc = hubbard_occ(nt,1)
         !
         IF (nspin .GE. 2) THEN
            IF (starting_magnetization(nt).GT.0.d0) THEN 
               nm = .FALSE.
               majs = 1  
               mins = 2  
            ELSEIF (starting_magnetization(nt).LT.0.d0) THEN 
               nm = .FALSE.
               majs = 2  
               mins = 1  
            ENDIF  
         ENDIF
         !
         IF (.NOT.nm) THEN
            ! Atom is magnetic
            ! 
            IF (noncolin) THEN
               !-- parameters for rotating occ. matrix
               cosin   = COS(angle1(nt)) 
               esin    = ( COS(angle2(nt)) + (0.d0,1.d0)*SIN(angle2(nt)) ) * SIN(angle1(nt))
               !--
               !-- occ. matrix in quantiz. axis  
               IF (totoc>ldim) THEN                 
                  ns(majs) = 1.d0
                  ns(mins) = (totoc -ldim ) / ldim                           
               ELSE
                  ns(majs) = totoc / ldim
                  ns(mins) = 0.d0
               ENDIF
               !--
               !-- charge and moment
               n =  ns(1) + ns(2) 
               m =  ns(1) - ns(2)  
               !--
               !-- rotating occ. matrix
               ns(1) = ( n + m*cosin ) / 2.d0 
               ns(2) = m * esin / 2.d0
               ns(3) = m * CONJG( esin ) / 2.d0 
               ns(4) = ( n - m*cosin ) / 2.d0 
               DO m1 = 1, ldim
                  nsg(m1,m1,viz,na,:) = ns(:)
               ENDDO
            ELSE
               !
               IF (totoc.GT.ldim) THEN
                  DO m1 = 1, ldim
                     nsg (m1,m1,viz,na,majs) = 1.d0
                     nsg (m1,m1,viz,na,mins) = (totoc - ldim) / ldim
                  ENDDO
               ELSE
                  DO m1 = 1, ldim
                     nsg (m1,m1,viz,na,majs) = totoc / ldim
                  ENDDO
               ENDIF
               !   
            ENDIF
         ELSE  
            ! Atom is non-magnetic
            IF (noncolin) THEN
               DO m1 = 1, ldim
                  nsg (m1,m1,viz,na,1) = totoc / 2.d0 / ldim
                  nsg (m1,m1,viz,na,4) = totoc / 2.d0 / ldim
               ENDDO
            ELSE
               DO is = 1, nspin
                  DO m1 = 1, ldim  
                     nsg (m1,m1,viz,na,is) = totoc /  2.d0 / ldim
                  ENDDO  
               ENDDO  
            ENDIF
         ENDIF  
         !
         ! Background part
         !
         IF ( is_hubbard_back(nt) ) THEN
            !     
            IF (.NOT.backall(nt)) THEN
               ! Fill in the second Hubbard manifold
               ldim2 = 2*Hubbard_l2(nt)+1
               totoc_b = hubbard_occ(nt,2)
               DO is = 1, nspin
                  DO m1 = ldim+1, ldim_u(nt)
                     nsg (m1,m1,viz,na,is) = totoc_b / 2.d0 / ldim2
                  ENDDO
               ENDDO
            ELSE
               ! Fill in the second Hubbard manifold
               ldim2 = 2*Hubbard_l2(nt)+1
               totoc_b = hubbard_occ(nt,2)
               DO is = 1, nspin
                  DO m1 = ldim+1, ldim+2*Hubbard_l2(nt)+1
                     nsg (m1,m1,viz,na,is) = totoc_b / 2.d0 / ldim2
                  ENDDO
               ENDDO
               ! Fill in the third Hubbard manifold
               ldim2 = 2*(Hubbard_l2(nt) + Hubbard_l3(nt)) + 2
               totoc_b = hubbard_occ(nt,3)
               DO is = 1, nspin
                  DO m1 = ldim+2*Hubbard_l2(nt)+2, ldim_u(nt)
                     nsg (m1,m1,viz,na,is) = totoc_b / 2.d0 / ldim2
                  ENDDO
               ENDDO
            ENDIF
            !
         ENDIF
         !
#if defined (__OSCDFT)
      IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2) .AND. &
         .NOT.oscdft_ctx%inp%constraint_diag) THEN
7        CONTINUE
      ENDIF
#endif
         !
      ENDIF
      ! 
   ENDDO ! na 
   !
   RETURN
   !  
END SUBROUTINE init_nsg
!-----------------------------------------------------------------------

