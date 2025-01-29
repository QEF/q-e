!!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------
SUBROUTINE write_epsilon(npe, drhop)
  !----------------------------------------------
  !! F. Macheda (2024)
  !! This routine takes the density (drhop) as input, computes the associated potential in real space,
  !! transforms in reciprocal space and takes its macroscopic component. Various response quantities
  !! are printed to file (fildvscf). The most important quantity is the inverse of the macroscopic dielectric function,
  !! corresponding to Eqs. (43) and (45) PRB 110, 094306 (2024)
  ! ---------------------------------------------
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : ionode
  USE fft_base,             ONLY : dfftp 
  USE units_ph,             ONLY : iurhoun 
  USE gvect,                ONLY : gg
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_interfaces,       ONLY : fwfft
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE constants,            ONLY : fpi,e2
  USE cell_base,            ONLY : tpiba
  USE qpoint,               ONLY : xq
  use mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE output,               ONLY : fildvscf
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npe
  !! input: the number of perturbation
  COMPLEX(DP), INTENT(IN) :: drhop(dfftp%nnr, nspin_mag, npe)
  !! input: change of the charge density (smooth and hard parts, dfftp)
  INTEGER, EXTERNAL :: find_free_unit
  !! to find available printing units 
  REAL(DP) :: vq
  !! Coulombian
  COMPLEX(DP) :: epsm1, chi, chi0, barchi
  !! auxiliary variables to enable printing
  COMPLEX(DP), ALLOCATABLE :: dvscf_toprint(:, :, :), drho_toprint(:, :, :)
  !! auxiliary variables to enable printing
  COMPLEX(DP) :: drhoaux, dvaux
  !! auxiliary variables to enable printing
  INTEGER :: ierr
  !! error status
  !
  ALLOCATE(dvscf_toprint(dfftp%nnr, nspin_mag , npe), STAT = ierr)
  IF (ierr /= 0) CALL errore('write_epsilon', 'Error allocating dvscf_toprint', 1)
  ALLOCATE(drho_toprint(dfftp%nnr, nspin_mag , npe), STAT = ierr)
  IF (ierr /= 0) CALL errore('write_epsilon', 'Error allocating drho_toprint', 1)
  !
  CALL zcopy(dfftp%nnr*nspin_mag, drhop, 1, dvscf_toprint, 1)
  CALL zcopy(dfftp%nnr*nspin_mag, drhop, 1, drho_toprint, 1)
  !
  CALL dv_of_drho(dvscf_toprint(1, 1, 1)) !nlcc cannot be applied since they depend on the perturbation, that we are taking as a scalar here
  !
  CALL fwfft('Rho', dvscf_toprint(:, 1, 1), dfftp)
  CALL fwfft('Rho', drho_toprint(:, 1, 1), dfftp)
  !
  drhoaux=0d0
  dvaux=0d0
  IF (gg(1) < 1d-8) THEN
    !      
    drhoaux = drho_toprint(dfftp%nl(1), 1, 1)
    dvaux   = dvscf_toprint(dfftp%nl(1), 1, 1)
    !
  ENDIF
  CALL mp_sum(drhoaux, intra_bgrp_comm)
  CALL mp_sum(dvaux, intra_bgrp_comm)
  !
  IF (ionode) THEN
    !      
    iurhoun=find_free_unit()
    OPEN(unit=iurhoun, file=fildvscf)
    !
    WRITE(iurhoun, '(a)') "#  Re(bar rho^{el}_q),Im(bar rho^{el}_q),Re(bar V^{tot}_q),Im(bar V^{tot}_q)"
    !
    WRITE(iurhoun,'(4f18.12)') REAL(drhoaux), AIMAG(drhoaux),&
                               REAL(dvaux)  , AIMAG(dvaux)
    !
    vq = fpi*e2/tpiba**2/DOT_PRODUCT(xq, xq)
    !
    chi0 = drhoaux/(1d0+dvaux)
    WRITE(iurhoun, '(a)') "#  $bar chi^0_q$"
    WRITE(iurhoun,'(2f18.12)') REAL(chi0), AIMAG(chi0)
    !
    WRITE(iurhoun, '(a)') "#  $bar chi_q$   Eq. (32) PRB 110, 094306 (2024)"
    !
    barchi = drhoaux
    WRITE(iurhoun,'(2f18.12)') REAL(barchi), AIMAG(barchi)
    !
    WRITE(iurhoun, '(a)') "# $1/epsilon^{-1}_{L}(q)$   Eqs. (43) and (45) PRB 110, 094306 (2024)"
    WRITE(iurhoun,'(2f18.12)') REAL(1d0 - vq * drhoaux), AIMAG(1d0 - vq * drhoaux)
    !
  ENDIF
  !
  CLOSE(iurhoun)
  DEALLOCATE (dvscf_toprint)
  DEALLOCATE (drho_toprint)
  !
!----------------------------------------------------------
END SUBROUTINE write_epsilon
!----------------------------------------------------------
!
!----------------------------------------------------------
SUBROUTINE write_drhoun
  !----------------------------------------------
  !! F. Macheda (2024)
  !! This routine takes the computed induced density by an atomic displacement, contained in dyn, and
  !! performs its symmetrization with respect to the small group of q. Then, it prints the results to file. 
  ! ---------------------------------------------
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout, ionode
  USE units_ph,       ONLY : iudumpdrho 
  USE modes,          ONLY : u
  USE ions_base,      ONLY : tau, nat, ityp
  USE gvect,          ONLY : mill, ig_l2g
  USE qpoint,         ONLY : xq
  USE cell_base,      ONLY : tpiba, omega
  USE uspp_param,     ONLY : upf
  USE constants,      ONLY : tpi, e2 
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE eqv,            ONLY : vlocq
  USE gvecs,          ONLY : ngms, ngms_g
  USE dynmat,         ONLY : dyn
  use mp,             ONLY : mp_sum
  USE output,         ONLY : fildrho
  !
  IMPLICIT NONE
  !
  INTEGER  :: ig, i, j, ipol, jpol, na, nb 
  !! loop indexes
  INTEGER, ALLOCATABLE :: itmp_mill(:, :)
  !! Miller indexes
  COMPLEX(DP) :: Im_i=(0._dp, 1._dp)
  !! Imaginary unit
  COMPLEX(DP), ALLOCATABLE :: phase(:)
  !! Phase factor
  REAL(DP)  :: zval
  !! Atomic charge
  REAL(DP)  :: arg 
  !! Argument of the phase
  INTEGER, EXTERNAL :: find_free_unit
  !! to find available printing units 
  COMPLEX(DP) :: phi(3, nat) 
  !! working vector to reconstruct the charge
  INTEGER :: ierr
  !! error status
  !
  ALLOCATE(phase(nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('write_drhoun', 'Error allocating phase', 1)
  !
  ALLOCATE(itmp_mill(3,ngms_g), STAT = ierr)
  IF (ierr /= 0) CALL errore('write_drhoun', 'Error allocating itmp_mill', 1)
  !
  itmp_mill =0
  !
  DO ig=1,ngms
    !
    itmp_mill(1, ig_l2g(ig)) = mill(1, ig)
    itmp_mill(2, ig_l2g(ig)) = mill(2, ig)
    itmp_mill(3, ig_l2g(ig)) = mill(3, ig)
    !
  ENDDO
  !
  CALL mp_sum(itmp_mill , intra_bgrp_comm)
  !
  !
  phi = (0.0d0, 0.0d0)
  !
  DO na = 1, nat
    !
    DO ipol = 1, 3
      !
      i = (na-1)*3 + ipol
      DO nb=1,nat
        !
        DO jpol = 1, 3
          !
          j = (nb-1)*3 + jpol
          phi(ipol, na) = phi(ipol,na) + dyn (i, j) * CONJG(u (i, j)) 
          !
        ENDDO
        !
      ENDDO
      !
    ENDDO
    !
  ENDDO
  !
  CALL symvectorq(nat,phi)
  !
  IF ( ionode ) THEN
    !
    iudumpdrho = find_free_unit()
    !
    OPEN(unit=iudumpdrho,file=fildrho)
    WRITE(iudumpdrho,'(a)') '#     Re(bar rho^{tot}_{qsx}),Im(bar rho^{tot}_{qsx}),Re(bar rho^{tot}_{qsy}),&
                         Im(bar rho^{tot}_{qsy}),Re(bar rho^{tot}_{qsz}),Im(bar rho^{tot}_{qsz}) '
    !
    DO ig = 1, ngms_g 
      !
      IF(ABS(itmp_mill(1, ig)) < 1d-8.and.&
         ABS(itmp_mill(2, ig)) < 1d-8.and.&
         ABS(itmp_mill(3, ig)) < 1d-8)THEN 
        !
        DO na = 1, nat
          !
          arg = (xq(1) * tau(1, na) + &
                 xq(2) * tau(2, na) + &
                 xq(3) * tau(3, na))*tpi
          !
          phase(na)= CMPLX(COS(arg), SIN(arg), KIND=DP)
          !
        ENDDO
        !
        DO na = 1, nat
          !
          zval = upf(ityp(na))%zp
          !
          write(iudumpdrho, '(6f18.12)') (-sqrt(e2)/omega*(phi(ipol, na)*phase(na)+&
                                           Im_i*xq(ipol)*tpiba*zval), ipol = 1, 3)
          !
        ENDDO 
        !
      ENDIF
      !
    ENDDO
    !
    CLOSE(iudumpdrho)
    !
  ENDIF
  !
  DEALLOCATE(phase)
  DEALLOCATE(itmp_mill)
  !
  WRITE(stdout, *) "------------------------------------------------------------------------"
  WRITE(stdout, *) " The code is printing the induced charges"
  WRITE(stdout, *) " Please refer to:"
  WRITE(stdout, *) " Macheda F., Barone P. & Mauri, F. (2024), "
  WRITE(stdout, *) " First-principles calculations of dynamical Born effective charges, quadrupoles,"
  WRITE(stdout, *) " and higher order terms from the charge response in large semiconducting and metallic systems"
  WRITE(stdout, *) " Physical Review B, 110, 094306. https://doi.org/10.1103/PhysRevB.110.094306"
  WRITE(stdout, *) "------------------------------------------------------------------------"
  !
  CONTAINS
    !--------------------------------------------------------------------------
    SUBROUTINE symvectorq(nat, vect)
      !-----------------------------------------------------------------------
      !! Symmetrize a function \(f(i,na)\) (e.g. the forces in cartesian axis),
      !! where \(i\) is the cartesian component, \(na\) the atom index.
      USE kinds,          ONLY : DP
      USE cell_base,      ONLY : at, bg
      USE symm_base,      ONLY : s, nsym, t_rev, irt, invs 
      USE lr_symm_base,   ONLY : minus_q, irotmq, nsymq, rtau
      USE qpoint,         ONLY : xq
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat
      !! number of atoms
      COMPLEX(DP), INTENT(INOUT) :: vect(3, nat)
      !! vector function to symmetrize
      !
      ! ... local variables
      !
      INTEGER :: na 
      COMPLEX(DP), ALLOCATABLE :: phi(:, :), work(:), phip(:, :)
      REAL(DP) :: arg, fase
      COMPLEX(DP) :: faseq(48)
      INTEGER :: iflb(nat), isymq, kpol, sna, irot
      INTEGER :: ierr
      !! error status
      !
      IF (nsym == 1) RETURN
      !
      ALLOCATE (phip(3,nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('symvectorq', 'Error allocating phip', 1)
      ALLOCATE (phi(3,nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('symvectorq', 'Error allocating phi', 1)
      ALLOCATE (work(3), STAT = ierr)
      IF (ierr /= 0) CALL errore('symvectorq', 'Error allocating work', 1)
      !
      ! bring vector to crystal axis
      !
      DO na = 1, nat
        !
        phi(:,na) =  vect(1,na)*at(1,:) + &
                     vect(2,na)*at(2,:) + &
                     vect(3,na)*at(3,:)
        !
      ENDDO
      !
      !    If no other symmetry is present we quit here
      !
      IF ((nsymq == 1) .AND. (.NOT. minus_q)) RETURN
      !
      !    Then we impose the symmetry q -> -q+G IF present
      !
      IF (minus_q) THEN
        !      
        DO na = 1, nat
          !
          DO ipol = 1, 3
            !
            work(:) = (0.d0, 0.d0)
            sna = irt(irotmq, na)
            arg = 0.d0
            !
            DO kpol = 1, 3
              !
              arg = arg + (xq(kpol) * (-rtau (kpol, irotmq, na))) 
              !
            ENDDO
            !
            arg = arg * tpi
            fase = CMPLX(cos(arg), sin (arg), kind=DP)
            !
            DO kpol = 1, 3
              !
              work(ipol) = work(ipol) + &
                           s(ipol, kpol, irotmq) &
                           * phi (kpol, sna) * fase
              !
            ENDDO
            !
            phip(ipol, na) = (phi(ipol, na) + &
                              CONJG(work(ipol))) * 0.5d0
            !
          ENDDO
          !
        ENDDO
        !
        phi = phip
        !
      ENDIF
      !
      !    Here we symmetrize with respect to the small group of q
      !
      IF (nsymq == 1) return
      !
      iflb(:) = 0
      !
      DO na = 1, nat
        !
        IF (iflb (na) == 0) THEN
          !
          work(:) = (0.d0, 0.d0)
          !
          DO isymq = 1, nsymq
            !
            irot = isymq
            sna = irt(irot, na)
            arg = 0.d0
            !
            DO ipol = 1, 3
              !
              arg = arg + (xq(ipol) * (-rtau (ipol, irot, na)))
              !
            ENDDO
            !
            arg = arg * tpi
            faseq(isymq) = CMPLX(cos(arg), sin(arg) , kind=DP)
            !
            DO ipol = 1, 3
              !
              DO kpol = 1, 3
                !
                IF (t_rev(isymq)==1) THEN
                  !       
                  work(ipol) = work(ipol) + &
                                s(ipol, kpol, irot) &
                                * CONJG(phi(kpol, sna) * faseq(isymq))
                  !
                ELSE
                  !
                  work(ipol) = work(ipol) + &
                                s(ipol, kpol, irot) &
                                * phi(kpol, sna) * faseq(isymq)
                  !
                ENDIF
                !
              ENDDO
              !
            ENDDO
            !
          ENDDO
          !
          DO isymq = 1, nsymq
            !
            irot = isymq
            sna = irt(irot, na)
            !
            DO ipol = 1, 3
              !
              phi(ipol, sna) = (0.d0, 0.d0)
              !
              DO kpol = 1, 3
                !
                IF (t_rev(isymq)==1) THEN
                  !
                  phi(ipol, sna) = phi(ipol, sna) &
                                + s(ipol, kpol, invs(irot))&
                                * CONJG(work(kpol)*faseq(isymq))
                  !
                ELSE
                  !
                  phi(ipol,sna) = phi(ipol, sna) &
                                + s(ipol, kpol, invs(irot))&
                                * work(kpol) * CONJG(faseq(isymq))
                  !
                ENDIF
                !
              ENDDO
              !
            ENDDO
            !
            iflb(sna) = 1
            !
          ENDDO
          !
        ENDIF
        !
      ENDDO
      !
      phi(:, :) = phi(:, :)/DBLE(nsymq)
      !
      ! bring vector back to cartesian axis
      !
      DO na = 1, nat
        !
        vect(:, na) = phi(1, na)*bg(:, 1) + &
                     phi(2, na)*bg(:, 2) + &
                     phi(3, na)*bg(:, 3)
        !
      END DO
      !
      DEALLOCATE (phi)
      DEALLOCATE (work)
      !
    !----------------------------------------------------------
    END SUBROUTINE symvectorq
    !----------------------------------------------------------
    !
!----------------------------------------------------------
END SUBROUTINE write_drhoun
!----------------------------------------------------------
!
!----------------------------------------------------------
COMPLEX(DP) FUNCTION Vaeps_dvloc(uact, ind_ig)
  !----------------------------------------------------------
  !! F. Macheda (2024)
  !! This routine takes the mode-th irreducible component of the potential, pot, and refers its macroscopic components to the value
  !!of the Coulomb potential (Eq. (B33) of Ref. PRB 110, 094306 (2024)). The index of the macroscopic component is given by ind_ig
  ! ---------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE fft_base,       ONLY : dffts
  USE modes,          ONLY : nmodes
  USE ions_base,      ONLY : nat, ityp, ntyp => nsp
  USE gvect,          ONLY : g, mill, eigts1, eigts2, eigts3, ngm
  USE qpoint,         ONLY : xq, eigqts
  USE cell_base,      ONLY : tpiba, tpiba2, omega
  USE uspp_param,     ONLY : upf
  USE gvecs,          ONLY : ngms
  !
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN) :: uact(nmodes)
  !! pattern of the representation
  INTEGER, INTENT(IN) :: ind_ig
  !! Index to be passed, that identifies the macroscopic component
  !
  ! ... local variables
  !
  REAL (DP), ALLOCATABLE :: vlocq(:, :)  ! ngm, ntyp)
  INTEGER :: na, mu, ig, nt
  INTEGER, ALLOCATABLE :: nl_d(:)
  COMPLEX(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  COMPLEX(DP) , ALLOCATABLE :: aux1 (:)
  REAL(DP) :: zval
  INTEGER :: ierr
  ! 
  ALLOCATE(nl_d(dffts%ngm), STAT = ierr )
  IF (ierr /= 0) CALL errore('Vaeps_dvloc', 'Error allocating nl_d', 1)
  nl_d  = dffts%nl
  ALLOCATE(aux1(dffts%nnr), STAT = ierr) !aux1 is dvlocin
  IF (ierr /= 0) CALL errore('Vaeps_dvloc', 'Error allocating aux1', 1)
  ALLOCATE(vlocq(ngm,ntyp), STAT = ierr)
  IF (ierr /= 0) CALL errore('Vaeps_dvloc', 'Error allocating vlocq', 1)
  !
  DO nt = 1, ntyp
    !
    zval=upf(nt)%zp
    CALL setlocq_coul(xq, zval, tpiba2, ngm, g, omega, vlocq(:,nt))
    !
  ENDDO
  !
  aux1 = 0.0d0
  !
  DO na = 1, nat
    !
    fact = tpiba * (0.d0, -1.d0) * eigqts(na)
    mu = 3 * (na - 1)
    !
    IF (abs(uact(mu + 1)) + abs(uact(mu + 2)) + &
        abs(uact(mu + 3)) > 1.0d-12) THEN
       !
      nt = ityp(na)
      u1 = uact(mu + 1)
      u2 = uact(mu + 2)
      u3 = uact(mu + 3)
      gu0 = xq(1) * u1 + xq(2) * u2 + xq(3) * u3
      !
      DO ig = 1, ngms
        !
        gtau = eigts1(mill(1,ig), na) * eigts2(mill(2,ig), na) * &
               eigts3(mill(3,ig), na)
        gu = gu0 + g(1, ig) * u1 + g(2, ig) * u2 + g(3, ig) * u3
        aux1 (dffts%nl(ig)) = aux1 (dffts%nl(ig)) + vlocq(ig, nt) &
                              * gu * fact * gtau
        !
      ENDDO
      !
    ENDIF
    !
  ENDDO
  !
  Vaeps_dvloc = - aux1(ind_ig)
  DEALLOCATE(aux1)
  !
  CONTAINS
    !----------------------------------------------------------------------
    SUBROUTINE setlocq_coul (xq, zp, tpiba2, ngm, g, omega, vloc)
     !----------------------------------------------------------------------
     !! Fourier transform of the Coulomb potential - For all-electron
     !! calculations, in specific cases only, for testing purposes.
     !
     USE kinds, ONLY: DP
     USE constants, ONLY : fpi, e2, eps8
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: ngm
     REAL(DP) :: xq (3), zp, tpiba2, omega, g(3, ngm)
     REAL(DP), INTENT (OUT) :: vloc(ngm)
     !
     REAL(DP) :: g2a
     INTEGER :: ig
     !
     DO ig = 1, ngm
       !
       g2a = (xq(1) + g(1, ig)) **2 + (xq(2) + g(2, ig)) **2 + &
             (xq(3) + g(3, ig)) **2
       !
       IF (g2a < eps8) THEN
         !
         vloc (ig) = 0.d0
         !
       ELSE
         !
         vloc (ig) = - fpi * zp *e2 / omega / tpiba2 / g2a
         !
       ENDIF
       !
     ENDDO
     !
    !----------------------------------------------------------------------
    END SUBROUTINE setlocq_coul
    !----------------------------------------------------------------------
    !
!----------------------------------------------------------
END FUNCTION Vaeps_dvloc
!----------------------------------------------------------------------
!
!----------------------------------------------------------
SUBROUTINE init_rho(npe, drhos, drhop, iq_dummy)
  ! ---------------------------------------------
  !! F. Macheda (2024)
  !! This routine simply reads the density drhop from file, and associate the correct drhos
  ! ---------------------------------------------
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : ionode
  USE fft_base,             ONLY : dfftp, dffts 
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_interfaces,       ONLY : fft_interpolate
  USE save_ph,              ONLY : tmp_dir_save
  USE output,               ONLY : fildrho
  USE units_ph,             ONLY : iudrho, lrdrho
  USE io_files,             ONLY : prefix, diropn
  USE qpoint,               ONLY : xq
  USE cell_base,            ONLY : at
  USE gvecs,                ONLY : doublegrid
  USE dfile_autoname,       ONLY : dfile_name
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npe
  !! input: the number of perturbation
  COMPLEX(DP), INTENT(INOUT) :: drhos(dffts%nnr, nspin_mag, npe)
  !! input/output: change of the charge density (smooth and hard parts, dfftp)
  COMPLEX(DP), INTENT(INOUT) :: drhop(dfftp%nnr, nspin_mag, npe)
  !! input/output: change of the charge density (smooth parts, dffts)
  INTEGER, INTENT(IN) :: iq_dummy
  !! needed for the filename that is being read
  !
  ! ...local variables
  !
  INTEGER :: ipert, is
  CHARACTER(LEN=256) :: filename
  LOGICAL :: exst
  !
  DO ipert = 1, npe
    !
    IF (fildrho .NE. ' ') THEN
      !
      IF (ionode) THEN
        ! 
        INQUIRE(UNIT = iudrho, OPENED = exst)
        IF (exst) CLOSE (UNIT = iudrho, STATUS='KEEP')
        filename = dfile_name(xq, at, fildrho, TRIM(tmp_dir_save)//prefix, generate=.TRUE., index_q=iq_dummy)
        CALL diropn (iudrho, filename, lrdrho, exst)
        !
      ENDIF ! ionode
      !     
      CALL davcio_drho (drhop(1,1,ipert), lrdrho, iudrho, 1, -1)
      !
    ENDIF
    !
  ENDDO
  CLOSE (UNIT = iudrho, STATUS='KEEP')
  !
  IF (doublegrid) THEN
    !
    DO is = 1, nspin_mag
      !
      DO ipert = 1, npe
        !
        CALL fft_interpolate (dfftp, drhop(:, is, ipert), dffts, drhos(:, is, ipert))
        !
      ENDDO
      !
    ENDDO
    !
  ELSE
    !
    CALL zcopy (npe*nspin_mag*dfftp%nnr, drhop, 1, drhos, 1)
    !
  ENDIF
  !
!----------------------------------------------------------
END SUBROUTINE init_rho
!----------------------------------------------------------

