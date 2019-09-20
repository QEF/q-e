  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE rotate
  !----------------------------------------------------------------------
  !! 
  !! This module contains the routines to perform symmetry rotations.
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE rotate_eigenm(iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2) 
    !--------------------------------------------------------------------------
    !!
    !!  Here:
    !! 
    !!  1). determine the unitary matrix which gives the transformed eigenmodes
    !!     from the first q in the star to the present q
    !! 
    !!  2). rotate eigenmodes from the first q in the star (cz1) to the current 
    !!      q (cz2)
    !!  
    !!  In rotate_epmat.f90:
    !! 
    !!  3). rotate the electron-phonon matrix from the cartesian representation
    !!     of the first qpoint of the star to the eigenmode representation 
    !!     (using cz1).
    !!  
    !!  4). rotate the electron-phonon matrix from the eigenmode representation
    !!     to the cartesian representation of the qpoint iq in the star (with cz2).
    !! 
    !!     Step 3 requires using the rotated eigenmodes to set the phases 
    !!     (the diagonalizers of the rotated dyn mat will not
    !!     work because of the gages in deg. subspaces and the phases).
    !! 
    !!  W.r.t. the standard method of q2qstar_ph.f90, here I construct the
    !!  unitary matrix Gamma defined in Maradudin & Vosko, RMP 1968.
    !!  In this way all rotations are just zgemm operations and we throw away  
    !!  all nasty indices.
    !!
    !!  SP - Sept 2019: Cleaning of the routine.  
    !
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE elph2,         ONLY : dynq
    USE phcom,         ONLY : nmodes
    USE constants_epw, ONLY : cone, czero, twopi, rydcm1, eps10, cmm12meV 
    USE control_flags, ONLY : iverbosity
    USE cell_base,     ONLY : at, bg
    USE ions_base,     ONLY : nat, amass, nat, ityp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq_first
    !! Originating q point of the star
    INTEGER, INTENT(in) :: nqc 
    !! Total number of qpoints
    INTEGER, INTENT(in) :: isym
    !! The symmetry under consideration
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! The symmetry operation for the eigenmodes
    INTEGER, INTENT(in) :: irt(48, nat)
    !! The rotated of each atom
    INTEGER, INTENT(in) :: invs(48)
    !! The index of the inverse operation
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! The rotated q vector
    REAL(KIND = DP), INTENT(in) :: rtau(3, 48, nat)
    !! The relative position of the rotated atom to the original one
    COMPLEX(KIND = DP), INTENT(inout) :: cz1(nmodes, nmodes)
    !! The eigenvectors for the first q in the star
    COMPLEX(KIND = DP), INTENT(inout) :: cz2(nmodes, nmodes)
    !! The rotated eigenvectors, for the current q in the star
    ! 
    ! Local variables
    INTEGER :: neig
    !! Lapack nb of eig. 
    INTEGER :: info
    !! 
    INTEGER :: ifail(nmodes)
    !! 
    INTEGER :: iwork(5 * nmodes)
    !! 
    INTEGER :: imode
    !! Mode index
    INTEGER :: jmode
    !! Mode index
    INTEGER :: nu
    !! Mode index
    INTEGER :: mu
    !! Mode index
    INTEGER :: sna
    !! 
    INTEGER :: ism1
    !! 
    INTEGER :: na
    !! Atom index
    INTEGER :: nb
    !! Atom index
    REAL(KIND = DP) :: arg
    !! 
    REAL(KIND = DP) :: massfac
    !! 
    REAL(KIND = DP) :: w1(nmodes)
    !! Phonon freq
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon freq
    REAL(KIND = DP) :: rwork(7 * nmodes)
    !! 
    REAL(KIND = DP) :: wtmp(nmodes)
    !! 
    REAL(KIND = DP) :: scart(3, 3)
    !! 
    COMPLEX(KIND = DP) :: gamma(nmodes, nmodes) 
    !! The Gamma matrix for the symmetry operation on the dyn mat 
    COMPLEX(KIND = DP) :: cwork(2 * nmodes) 
    !! 
    COMPLEX(KIND = DP) :: dynp(nmodes * (nmodes + 1) / 2) 
    !! Complex dynmat packed (upper triangular part for zhpevx)
    COMPLEX(KIND = DP) :: cfac 
    !! 
    COMPLEX(KIND = DP) :: dyn1(nmodes, nmodes) 
    !! Dynamical matrix
    COMPLEX(KIND = DP) :: dyn2(nmodes, nmodes) 
    !! Dynamical matrix
    !
    ! ------------------------------------------------------------------
    ! diagonalize dynq(iq_first) --> w1, cz1
    ! ------------------------------------------------------------------
    !
    ! first divide by the square root of masses 
    !
    DO na = 1, nat
      DO nb = 1, nat
        massfac = 1.d0 / SQRT(amass(ityp(na)) * amass(ityp(nb)))
        dyn1(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
        dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, iq_first) * massfac
      ENDDO
    ENDDO
    !
    DO jmode = 1, nmodes
      DO imode = 1, jmode
        dynp(imode + (jmode - 1) * jmode / 2) = (dyn1(imode, jmode) + CONJG(dyn1(jmode, imode))) / 2.d0
      ENDDO
    ENDDO
    !
    CALL ZHPEVX('V', 'A', 'U', nmodes, dynp , 0.0, 0.0, &
                 0, 0, -1.0, neig, w1, cz1, nmodes, cwork, rwork, iwork, ifail, info)
    !
    IF (iverbosity == 1) THEN
      !
      ! check the frequencies
      DO nu = 1, nmodes
        IF (w1(nu) > 0.d0) THEN
          wtmp(nu) =  SQRT(ABS(w1(nu)))
        ELSE
          wtmp(nu) = -SQRT(ABS(w1(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the first q in the star (cm^{-1})")') 
      WRITE(stdout, '(6(2x,f10.5))') (wtmp(nu) * rydcm1, nu = 1, nmodes)
      !
    ENDIF
    !
    ! here dyn1 is dynq(:,:,iq_first) after dividing by the masses
    ! (the true dyn matrix D_q)
    !
    ! -----------------------------------------------------------------------
    ! the matrix gamma (Maradudin & Vosko, RMP, eq. 2.37)   
    ! -----------------------------------------------------------------------
    !
    ! I have built the matrix by following step-by-step q2star_ph.f90 and
    ! rotate_and_add_dyn.f90
    !
    ism1 = invs(isym)
    !
    ! the symmetry matrix in cartesian coordinates 
    ! (so that we avoid going back and forth with the dynmat)  
    ! note the presence of both at and bg in the transform!
    !
    scart = DBLE(s(:, :, ism1))
    scart = MATMUL(MATMUL(bg, scart), TRANSPOSE(at))
    ! 
    gamma = czero
    DO na = 1, nat
      !
      ! the corresponding of na in {S|v}
      sna = irt(isym, na)
      !
      ! cfac = exp[iSq*(tau_K - {S|v} tau_k)]   (Maradudin&Vosko RMP Eq. 2.33)
      ! [v can be ignored since it cancels out, see endnotes. xq is really Sq]
      ! rtau(:,isym,na) = s(:,:,invs(isym)) * tau(:, na) - tau(:,irt(isym,na))) (cartesian)
      !
      arg = twopi * DOT_PRODUCT(xq, rtau (:, isym, na)) 
      cfac = dcmplx(COS(arg), -SIN(arg))
      !
      ! the submatrix (sna,na) contains the rotation scart 
      !
      gamma(3 * (sna - 1) + 1:3 * sna, 3 * (na - 1) + 1:3 * na) = cfac * scart 
      !
    ENDDO
    !
    !  possibly run some consistency checks
    !
    IF (iverbosity == 1) THEN
      !
      !  D_{Sq} = gamma * D_q * gamma^\dagger (Maradudin & Vosko, RMP, eq. 3.5) 
      ! 
      CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
             nmodes, dyn1 , nmodes, czero , dyn2, nmodes)
      CALL ZGEMM('n', 'c', nmodes, nmodes, nmodes, cone, dyn2, &
             nmodes, gamma, nmodes, czero, dyn1, nmodes)
      !
      DO jmode = 1, nmodes
        DO imode = 1, jmode
          dynp(imode + (jmode - 1) * jmode / 2) = (dyn1(imode, jmode) + CONJG(dyn1(jmode, imode))) / 2.d0
        ENDDO
      ENDDO
      !
      CALL ZHPEVX('V', 'A', 'U', nmodes, dynp, 0.0, 0.0, &
                  0, 0, -1.0, neig, w2, cz2, nmodes, cwork, rwork, iwork, ifail, info)
      !
      ! Check the frequencies
      !
      DO nu = 1, nmodes
        IF (w2(nu) > 0.d0) THEN
          wtmp(nu) =  SQRT(ABS(w2(nu)))
        ELSE
          wtmp(nu) = -SQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the first q in the star (meV)")') 
      WRITE(stdout, '(6(2x,f10.5))') (wtmp(nu) * rydcm1 * cmm12meV, nu = 1, nmodes)
      !
    ENDIF
    !
    !
    ! -----------------------------------------------------------------------
    ! transformation of the eigenvectors: e_{Sq} = gamma * e_q  (MV eq. 2.36)
    ! -----------------------------------------------------------------------
    !
    CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
         nmodes, cz1 , nmodes, czero , cz2, nmodes)
    !
    !   now check that cz2 diagonalizes dyn2 (after dividing by the masses)
    !
    DO na = 1, nat
      DO nb = 1, nat
        massfac = 1.d0 / SQRT(amass(ityp(na)) * amass(ityp(nb)))
        dyn2(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
          dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, nqc) * massfac
      ENDDO
    ENDDO
    !
    ! the rotated matrix and the one read from file
    IF (iverbosity == 1) WRITE(stdout, '(2f15.10)') dyn2 - dyn1 
    !
    ! here I have checked that the matrix rotated with gamma 
    ! is perfectly equal to the one read from file for this q in the star
    !
    CALL ZGEMM('c', 'n', nmodes, nmodes, nmodes, cone, cz2, nmodes, dyn2, nmodes, czero, dyn1, nmodes)
    CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, dyn1, nmodes, cz2, nmodes, czero, dyn2, nmodes)
    !
    DO nu = 1, nmodes
      w2(nu) = ABS(dyn2(nu, nu))
      DO mu = 1, nmodes
        IF (mu /= nu .AND. ABS(dyn2(mu, nu)) > eps10) CALL errore('rotate_eigenm', 'problem with rotated eigenmodes', 0)
      ENDDO
    ENDDO
    !
    IF (iverbosity == 1) THEN
      !
      ! A simple check on the frequencies
      !
      DO nu = 1, nmodes
        IF (w2(nu) > 0.d0 ) then
          wtmp(nu) =  SQRT(ABS(w2(nu)))
        ELSE
          wtmp(nu) = -SQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the current q in the star (cm^{-1})")') 
      WRITE(stdout, '(6(2x,f10.5))') (wtmp(nu) * rydcm1, nu = 1, nmodes)
      !
    ENDIF
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE rotate_eigenm
    !--------------------------------------------------------------------------
    !                                                                            
    !---------------------------------------------------------------------------
    SUBROUTINE rotate_epmat(cz1, cz2, xq, iq, lwin, lwinq, exband)
    !---------------------------------------------------------------------------
    !!
    !! 1). rotate the electron-phonon matrix from the cartesian representation
    !!    of the first qpoint of the star to the eigenmode representation 
    !!    (using cz1).
    !! 
    !! 2). rotate the electron-phonon matrix from the eigenmode representation
    !!     to the cartesian representation of the qpoint iq (with cz2).
    !!
    !! SP - Sep. 2019: Cleaning. 
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : epmatq, zstar, epsi, bmat
    USE epwcom,        ONLY : lpolar, nqc1, nqc2, nqc3
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : cone, czero, one, ryd2mev, eps8
    USE pwcom,         ONLY : nbnd, nks
    USE ions_base,     ONLY : amass, ityp
    USE rigid_epw,     ONLY : rgd_blk_epw
    ! 
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lwin(nbnd, nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: lwinq(nbnd, nks)
    !! Bands at k+q within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    INTEGER, INTENT(in) :: iq
    !!  Current qpoint
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !  Rotated q vector
    COMPLEX(KIND = DP), INTENT(inout) :: cz1(nmodes, nmodes)
    !! eigenvectors for the first q in the star
    COMPLEX(KIND = DP), INTENT(inout) :: cz2(nmodes, nmodes)
    !!  Rotated eigenvectors for the current q in the star
    !
    ! Local variables 
    INTEGER :: mu
    !! Counter on phonon branches
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: i
    !! Counter on band index
    INTEGER :: j
    !! Counter on band index
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    REAL(KIND = DP) :: massfac
    !! square root of mass 
    COMPLEX(KIND = DP) :: eptmp(nmodes)
    !! temporary e-p matrix elements
    COMPLEX(KIND = DP) :: epmatq_opt(nbnd, nbnd, nks, nmodes)
    !! e-p matrix elements in the outer window
    COMPLEX(KIND = DP) :: epmatq_tmp(nbnd, nbnd, nks, nmodes)
    !! temporary e-p matrix 
    COMPLEX(KIND = DP) :: cz_tmp(nmodes, nmodes)
    !! temporary variables
    COMPLEX(KIND = DP) :: cz2t(nmodes, nmodes)
    !! temporary variables
    !
    ! the mass factors: 
    !  1/sqrt(M) for the  direct transform
    !  SQRT(M)   for the inverse transform 
    !
    ! if we set cz1 = cz2 here and we calculate below
    ! cz1 * cz2 we find the identity
    !
    cz2t = cz2
    !
    DO mu = 1, nmodes
      na = (mu - 1) / 3 + 1
      massfac = SQRT(amass(ityp(na)))
      cz1(mu, :) = cz1(mu, :) / massfac
      cz2(mu, :) = cz2(mu, :) * massfac
      cz2t(mu, :) = cz2t(mu, :) / massfac
    ENDDO
    !
    ! the inverse transform also requires the hermitian conjugate
    !
    cz_tmp = CONJG(TRANSPOSE(cz2))
    cz2 = cz_tmp
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF
    ENDDO
    ! 
    ! slim down to the first ndimwin(ikq), ndimwin(ik) states within the outer window
    !
    epmatq_opt = czero
    epmatq_tmp = czero
    IF (nexband_tmp > 0) THEN
      DO ik = 1, nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          IF (lwin(j, ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              IF (lwinq(i, ik)) THEN
                ibnd = ibnd + 1
                epmatq_tmp(ibnd, jbnd, ik, :) = epmatq(i, j, ik, :, iq)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
                ibnd = ibnd + 1
                epmatq_opt(i, j, ik, :) = epmatq_tmp(ibnd, jbnd, ik, :)
            ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        jbnd = 0
        DO j = 1, nbnd
          IF (lwin(j, ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (lwinq(i, ik)) THEN
                ibnd = ibnd + 1
                epmatq_opt(ibnd, jbnd, ik, :) = epmatq(i, j, ik, :, iq)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    ! 
    ! ep_mode(j) = cfac * sum_i ep_cart(i) * u(i,j)
    !
    epmatq(:, :, :, :, iq) = czero
    DO ik = 1, nks
      DO jbnd = 1, nbnd
        DO ibnd = 1, nbnd
          !
          ! bring e-p matrix from the cartesian representation of the
          ! first q in the star to the corresponding eigenmode representation
          !
          CALL ZGEMV('t', nmodes, nmodes, cone, cz1, nmodes,  &
                     epmatq_opt(ibnd, jbnd, ik, :), 1, czero, eptmp, 1)
          !
          IF (lpolar) THEN
            IF ((ABS(xq(1)) > eps8) .OR. (ABS(xq(2)) > eps8) .OR. (ABS(xq(3)) > eps8)) THEN
              CALL rgd_blk_epw(nqc1, nqc2, nqc3, xq, cz2t, eptmp, &
                       nmodes, epsi, zstar, bmat(ibnd, jbnd, ik, iq), -one)
            ENDIF
          ENDIF
          !
          ! rotate epmat in the cartesian representation for this q in the star
          !
          CALL ZGEMV('t', nmodes, nmodes, cone, cz2, nmodes, &
                    eptmp, 1, czero, epmatq(ibnd, jbnd, ik, :, iq), 1)
        ENDDO
      ENDDO
    ENDDO
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE rotate_epmat
    !---------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  END MODULE rotate
  !-----------------------------------------------------------------------------

