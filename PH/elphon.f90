!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE elphon()
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  !
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau
  USE pwcom
  USE kinds, ONLY : DP
  USE phcom
  USE el_phon
  !
  IMPLICIT NONE
  !
  INTEGER :: irr, imode0, ipert, is
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(kind=DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)

  CALL start_clock ('elphon')

  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  !!! rewind (iudvscf)
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin ( nrxx , nspin , npert(irr)) )
     !!! read (iudvscf) dvscfin
     DO ipert = 1, npert (irr)
        CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
                           imode0 + ipert, -1 )
     END DO
     IF (doublegrid) THEN
        ALLOCATE (dvscfins ( nrxxs , nspin , npert(irr)) )
        DO is = 1, nspin
           DO ipert = 1, npert(irr)
              CALL cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           ENDDO 
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     CALL newdq (dvscfin, npert(irr))
     CALL elphel (npert (irr), imode0, dvscfins)
     !
     imode0 = imode0 + npert (irr)
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
  ENDDO
  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
  !
  IF (.NOT.trans) CALL readmat (iudyn, ibrav, celldm, nat, ntyp, &
       ityp, omega, amass, tau, xq, w2, dyn)
  !
  CALL stop_clock ('elphon')
  RETURN
END SUBROUTINE elphon
!
!-----------------------------------------------------------------------
SUBROUTINE readmat (iudyn, ibrav, celldm, nat, ntyp, ityp, omega, &
     amass, tau, q, w2, dyn)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  ! Input
  INTEGER :: iudyn, ibrav, nat, ntyp, ityp (nat)
  REAL(kind=DP) :: celldm (6), amass (ntyp), tau (3, nat), q (3), &
       omega
  ! output
  REAL(kind=DP) :: w2 (3 * nat)
  COMPLEX(kind=DP) :: dyn (3 * nat, 3 * nat)
  ! local (control variables)
  INTEGER :: ntyp_, nat_, ibrav_, ityp_
  REAL(kind=DP) :: celldm_ (6), amass_, tau_ (3), q_ (3)
  ! local
  REAL(kind=DP) :: dynr (2, 3, nat, 3, nat)
  CHARACTER(len=80) :: line
  CHARACTER(len=3)  :: atm
  INTEGER :: nt, na, nb, naa, nbb, nu, mu, i, j
  !
  !
  REWIND (iudyn)
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
  IF (ntyp.NE.ntyp_.OR.nat.NE.nat_.OR.ibrav_.NE.ibrav.OR.ABS ( &
       celldm_ (1) - celldm (1) ) .GT.1.0d-5) CALL errore ('readmat', &
       'inconsistent data', 1)
  DO nt = 1, ntyp
     READ (iudyn, * ) i, atm, amass_
     IF (nt.NE.i.OR.ABS (amass_ - amass (nt) ) .GT.1.0d-5) CALL errore ( &
          'readmat', 'inconsistent data', 1 + nt)
  ENDDO
  DO na = 1, nat
     READ (iudyn, * ) i, ityp_, tau_
     IF (na.NE.i.OR.ityp_.NE.ityp (na) ) CALL errore ('readmat', &
          'inconsistent data', 10 + na)
  ENDDO
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (line (11:80), * ) (q_ (i), i = 1, 3)
  READ (iudyn, '(a)') line
  DO na = 1, nat
     DO nb = 1, nat
        READ (iudyn, * ) naa, nbb
        IF (na.NE.naa.OR.nb.NE.nbb) CALL errore ('readmat', 'error reading &
             &file', nb)
        READ (iudyn, * ) ( (dynr (1, i, na, j, nb), dynr (2, i, na, j, nb) &
             , j = 1, 3), i = 1, 3)
     ENDDO
  ENDDO
  !
  ! divide the dynamical matrix by the masses
  !
  DO nb = 1, nat
     DO j = 1, 3
        DO na = 1, nat
           DO i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! solve the eigenvalue problem.
  ! NOTA BENE: eigenvectors are overwritten on dyn
  !
  CALL cdiagh (3 * nat, dynr, 3 * nat, w2, dyn)
  !
  ! divide by sqrt(mass) to get displacements
  !
  DO nu = 1, 3 * nat
     DO mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu) = dyn (mu, nu) / SQRT (amass (ityp (na) ) )
     ENDDO
  ENDDO
  !
  !
  RETURN
END SUBROUTINE readmat
!
!-----------------------------------------------------------------------
SUBROUTINE elphel (npe, imode0, dvscfins)
  !-----------------------------------------------------------------------
  !
  !      Calculation of the electron-phonon matrix elements el_ph_mat
  !         <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !      Original routine written by Francesco Mauri
  !
  USE pwcom
  USE wavefunctions_module,  ONLY: evc
  USE kinds, ONLY : DP
  USE io_files, ONLY: iunigk
  USE phcom
  USE el_phon
  IMPLICIT NONE
  !
  INTEGER :: npe, imode0
  COMPLEX(kind=DP) :: dvscfins (nrxxs, nspin, npe)
  ! LOCAL variables
  INTEGER :: ik, ikk, ikq, ipert, mode, nrec, ibnd, jbnd, ir, ig, &
       ios
  COMPLEX(kind=DP) , ALLOCATABLE :: aux1 (:), elphmat (:,:,:)
  COMPLEX(kind=DP) :: ZDOTC
  !
  ALLOCATE (aux1    ( nrxxs))    
  ALLOCATE (elphmat ( nbnd , nbnd , npe))    
  !
  !  Start the loops over the k-points
  !
  IF (nksq.GT.1) REWIND (unit = iunigk)
  DO ik = 1, nksq
     IF (nksq.GT.1) THEN
        READ (iunigk, err = 100, iostat = ios) npw, igk
100     CALL errore ('elphel', 'reading igk', ABS (ios) )
     ENDIF
     !
     !  ik = counter of k-points with vector k
     !  ikk= index of k-point with vector k
     !  ikq= index of k-point with vector k+q
     !       k and k+q are alternated if q!=0, are the same if q=0
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
        npwq = npw
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     IF (lsda) current_spin = isk (ikk)
     IF (.NOT.lgamma.AND.nksq.GT.1) THEN
        READ (iunigk, err = 200, iostat = ios) npwq, igkq
200     CALL errore ('elphel', 'reading igkq', ABS (ios) )
     ENDIF
     !
     CALL init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     !
     ! read unperturbed wavefuctions psi(k) and psi(k+q)
     !
     IF (nksq.GT.1) THEN
        IF (lgamma) THEN
           CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
        ELSE
           CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
           CALL davcio (evq, lrwfc, iuwfc, ikq, - 1)
        ENDIF
     ENDIF
     !
     DO ipert = 1, npe
        nrec = (ipert - 1) * nksq + ik
        !
        !  dvbare_q*psi_kpoint is read from file (if available) or recalculated
        !
        IF (trans) THEN
           CALL davcio (dvpsi, lrbar, iubar, nrec, - 1)
        ELSE
           mode = imode0 + ipert
           ! TODO : .false. or .true. ???
           CALL dvqpsi_us (ik, mode, u (1, mode), .FALSE. )
        ENDIF
        !
        ! calculate dvscf_q*psi_k
        !
        DO ibnd = 1, nbnd
           aux1(:) = (0.d0, 0.d0)
           DO ig = 1, npw
              aux1 (nls (igk (ig) ) ) = evc (ig, ibnd)
           ENDDO
           CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
           DO ir = 1, nrxxs
              aux1 (ir) = aux1 (ir) * dvscfins (ir, current_spin, ipert)
           ENDDO
           CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
           DO ig = 1, npwq
              dvpsi (ig, ibnd) = dvpsi (ig, ibnd) + aux1 (nls (igkq (ig) ) )
           ENDDO
        END DO
        CALL adddvscf (ipert, ik)

        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        DO ibnd =1, nbnd
           DO jbnd = 1, nbnd
              elphmat (jbnd, ibnd, ipert) = ZDOTC (npwq, evq (1, jbnd), 1, &
                   dvpsi (1, ibnd), 1)
           ENDDO
           !
        ENDDO
     ENDDO
     !
     CALL reduce (2 * nbnd * nbnd * npe, elphmat)
     !
     !  save all e-ph matrix elements into el_ph_mat
     !
     DO ipert = 1, npe
        DO jbnd = 1, nbnd
           DO ibnd = 1, nbnd
              el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = elphmat (ibnd, jbnd, ipert)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE (elphmat)
  DEALLOCATE (aux1)
  !
  RETURN
END SUBROUTINE elphel
!
!-----------------------------------------------------------------------
SUBROUTINE elphsum
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri
  !
  USE ions_base, ONLY : nat
  USE pwcom
  USE kinds, ONLY : DP
  USE phcom
  USE el_phon
  USE mp_global, ONLY : me_pool, root_pool
  !
  IMPLICIT NONE
  ! eps = 20 cm^-1, in Ry
  REAL(kind=DP) :: eps
  PARAMETER (eps = 20.d0 / 13.6058d0 / 8065.5d0)
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, nsig, iuelph, ios
  REAL(kind=DP) :: weight, w0g1, w0g2, w0gauss, degauss1, dosef, dos_ef, &
       ef1, phase_space, lambda, gamma
  EXTERNAL dos_ef
  !
  COMPLEX(kind=DP) :: el_ph_sum (3*nat,3*nat)

  !
  WRITE (6, '(5x,"electron-phonon interaction  ..."/)')
  ngauss1 = 1
  nsig = 10
  IF (filelph.NE.' ') THEN

     ! parallel case: only first node writes
     IF ( me_pool /= root_pool ) THEN
        iuelph = 0
     ELSE
        !
        iuelph = 4
        OPEN (unit = iuelph, file = filelph, status = 'unknown', err = &
             100, iostat = ios)
100     CALL errore ('elphon', 'opening file'//filelph, ABS (ios) )
        REWIND (iuelph)
        WRITE (iuelph, '(3f15.8,2i8)') xq, nsig, 3 * nat
        WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
        !
     END IF
     !
  ELSE
     iuelph = 0
  ENDIF
  !
  !
  DO isig = 1, nsig
     degauss1 = 0.01 * isig
     el_ph_sum(:,:) = (0.d0, 0.d0)
     phase_space = 0.d0
     !
     ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
     ! for this gaussian broadening
     !
     ! Note that the weights of k+q points must be set to zero for the
     ! following call to yield correct results
     !
     CALL efermig (et, nbnd, nks, nelec, wk, degauss1, ngauss1, ef1, 0, isk)
     dosef = dos_ef (ngauss1, degauss1, ef1, et, wk, nks, nbnd)
     ! N(Ef) is the DOS per spin, not summed over spin
     dosef = dosef / 2.d0
     !
     ! Sum over bands with gaussian weights
     !
     DO ik = 1, nksq

        !
        ! see subroutine elphel for the logic of indices
        !
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        DO ibnd = 1, nbnd
           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
                / degauss1
           DO jbnd = 1, nbnd
              w0g2 = w0gauss ( (ef1 - et (jbnd, ikq) ) / degauss1, ngauss1) &
                   / degauss1
              ! note that wk(ikq)=wk(ikk)
              weight = wk (ikk) * w0g1 * w0g2
              DO jpert = 1, 3 * nat
                 DO ipert = 1, 3 * nat
                    el_ph_sum (ipert, jpert) = el_ph_sum (ipert, jpert)  +  weight * &
                                        CONJG (el_ph_mat (jbnd, ibnd, ik, ipert) ) * &
                                               el_ph_mat (jbnd, ibnd, ik, jpert)
                 ENDDO
              ENDDO
              phase_space = phase_space+weight
           ENDDO
        ENDDO

     ENDDO
     !
     ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
     !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
     !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
     !
     ! collect contributions from all pools (sum over k-points)
     !
     CALL poolreduce (2 * 3 * nat * 3 * nat, el_ph_sum)
     CALL poolreduce (1, phase_space)
     !
     ! symmetrize el_ph_sum(mu,nu) : it transforms as the dynamical matrix
     !
     CALL symdyn_munu (el_ph_sum, u, xq, s, invs, rtau, irt, irgq, at, &
          bg, nsymq, nat, irotmq, minus_q)
     !
     WRITE (6, 9000) degauss1, ngauss1
     WRITE (6, 9005) dosef, ef1 * 13.6058
     WRITE (6, 9006) phase_space
     IF (iuelph.NE.0) THEN
        WRITE (iuelph, 9000) degauss1, ngauss1
        WRITE (iuelph, 9005) dosef, ef1 * 13.6058
     ENDIF
     !
     DO nu = 1, nmodes
        gamma = 0.0
        DO mu = 1, 3 * nat
           DO vu = 1, 3 * nat
              gamma = gamma + REAL (CONJG (dyn (mu, nu) ) * el_ph_sum (mu, vu) &
                   * dyn (vu, nu) )
           ENDDO
        ENDDO
        gamma = 3.1415926 * gamma / 2.d0
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        IF (SQRT (ABS (w2 (nu) ) ) .GT.eps) THEN
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lambda = gamma / 3.1415926 / w2 (nu) / dosef
        ELSE
           lambda = 0.0
        ENDIF
        ! 3.289828x10^6 is the conversion factor from Ry to GHz
        WRITE (6, 9010) nu, lambda, gamma * 3.289828d6
        IF (iuelph.NE.0) WRITE (iuelph, 9010) nu, lambda, gamma * &
             3.289828d6
     ENDDO
  ENDDO
9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
       &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i2,')=',f10.6,'   gamma=',f10.6,' GHz')
  !
  !
  IF (iuelph.NE.0) CLOSE (unit = iuelph)
  RETURN
END SUBROUTINE elphsum
!
!-----------------------------------------------------------------------
FUNCTION dos_ef (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(kind=DP) :: dos_ef
  INTEGER :: ngauss, nbnd, nks
  REAL(kind=DP) :: et (nbnd, nks), wk (nks), ef, degauss
  !
  INTEGER :: ik, ibnd
  REAL(kind=DP) :: w0gauss
  !
  !     Compute DOS at E_F (states per Ry per unit cell)
  !
  dos_ef = 0.0
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        dos_ef = dos_ef + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     ENDDO
  ENDDO
  !
  !    Collects partial sums on k-points from all pools
  !
  CALL poolreduce (1, dos_ef)
  !
  RETURN
END FUNCTION dos_ef
