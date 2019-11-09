!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE init_ns
   !-----------------------------------------------------------------------
   !! This routine computes the starting ns (for lda+U calculation) filling
   !! up the d states (the only interested by the on-site potential for the
   !! moment) according to the Hund's rule (valid for the isolated atoms on
   !! which starting potential is built), and to the starting_magnetization:
   !! majority spin levels are populated first, THEN the remaining electrons
   !! are equally distributed among the minority spin states.
   !
   USE kinds,        ONLY : DP
   USE ions_base,    ONLY : nat, ityp
   USE lsda_mod,     ONLY : nspin, starting_magnetization
   USE ldaU,         ONLY : hubbard_u, hubbard_alpha, hubbard_l
   USE scf,          ONLY : rho
   USE uspp_param,   ONLY : upf
   !
   IMPLICIT NONE
   !
   REAL(DP) :: totoc
   REAL(DP), EXTERNAL :: hubbard_occ
   !
   INTEGER :: ldim, na, nt, is, m1, majs, mins
   LOGICAL :: nm        ! true if the atom is non magnetic
   !
   rho%ns(:,:,:,:) = 0.d0
   !
   DO na = 1, nat
      nt = ityp(na)
      IF (Hubbard_U(nt)/=0.d0 .OR. Hubbard_alpha(nt)/=0.d0) THEN
         ldim = 2*Hubbard_l(nt)+1
         totoc = hubbard_occ ( upf(nt)%psd )
         nm = .TRUE.
         IF (nspin==2) THEN
            IF (starting_magnetization(nt) > 0.d0) THEN  
               nm = .FALSE.
               majs = 1  
               mins = 2  
            ELSEIF (starting_magnetization(nt) < 0.d0) THEN  
               nm = .FALSE.
               majs = 2  
               mins = 1  
            ENDIF  
         ENDIF
         IF (.NOT.nm) THEN  
            IF (totoc>ldim) THEN  
               DO m1 = 1, ldim  
                  rho%ns(m1,m1,majs,na) = 1.d0  
                  rho%ns(m1,m1,mins,na) = (totoc - ldim) / ldim
               ENDDO  
            ELSE  
               DO m1 = 1, ldim  
                  rho%ns(m1,m1,majs,na) = totoc / ldim
               ENDDO  
            ENDIF  
         ELSE  
            DO is = 1,nspin
               DO m1 = 1, ldim  
                  rho%ns(m1,m1,is,na) = totoc /  2.d0 / ldim
               ENDDO  
            ENDDO  
         ENDIF  
      ENDIF  
   ENDDO
   !
   RETURN
   !
END SUBROUTINE init_ns
!
!
!-----------------------------------------------------------------------
SUBROUTINE init_ns_nc
   !--------------------------------------------------------------------
   !! Noncollinear version of \(\textrm{init}\_\textrm{ns}\) (A. Smogunov).
   !
   USE kinds,            ONLY : DP
   USE ions_base,        ONLY : nat, ityp
   USE lsda_mod,         ONLY : nspin, starting_magnetization
   USE ldaU,             ONLY : hubbard_u, hubbard_l
   USE noncollin_module, ONLY : angle1, angle2 
   USE scf,              ONLY : rho
   USE uspp_param,       ONLY : upf
   !
   IMPLICIT NONE
   !
   REAL(DP) :: totoc, cosin 
   REAL(DP), EXTERNAL :: hubbard_occ
   COMPLEX(DP) :: esin, n, m, ns(4)  
   !
   INTEGER :: ldim, na, nt, is, m1, m2, majs, isym, mins
   LOGICAL :: nm        ! true if the atom is non magnetic
   !
   rho%ns_nc(:,:,:,:) = 0.d0
   !
   DO na = 1, nat
      nt = ityp (na)
      IF (Hubbard_U(nt) /= 0.d0) THEN
         ldim = 2*Hubbard_l(nt) + 1
         totoc = hubbard_occ( upf(nt)%psd )
         nm=.TRUE.
         IF (starting_magnetization(nt) > 0.d0) THEN  
            nm = .FALSE.
            majs = 1  
            mins = 2  
         ELSEIF (starting_magnetization(nt) < 0.d0) THEN  
            nm = .FALSE.
            majs = 2  
            mins = 1  
         ENDIF
         !
         IF (.NOT.nm) THEN  
            !
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
              rho%ns_nc(m1,m1,:,na) = ns(:)
            ENDDO
            !--
         ELSE
            !
            DO m1 = 1, ldim  
              rho%ns_nc(m1,m1,1,na) = totoc /  2.d0 / ldim
              rho%ns_nc(m1,m1,4,na) = totoc /  2.d0 / ldim
            ENDDO
            !  
         ENDIF
         !
      ENDIF  
   ENDDO  
   !
   RETURN
   !
END SUBROUTINE init_ns_nc


