!
! Copyright (C) 2016  Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Author: Mitsuaki Kawamura, U. Tokyo
!----------------------------------------------------------------------------
MODULE elph_tetra_mod
  !--------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  LOGICAL,SAVE :: &
  &  lshift_q = .false., &
  &  in_alpha2f = .FALSE.
  !
  INTEGER,SAVE :: elph_tetra = 0 ! switch to output electron-phonon matrix
  !
  PUBLIC elph_tetra, lshift_q, in_alpha2f, &
  &      elph_tetra_lambda, elph_tetra_gamma
  !
  CONTAINS
 !
!----------------------------------------------------------------------------
SUBROUTINE elph_tetra_lambda()
  !--------------------------------------------------------------------------
  !
  ! This routine computes the electron-phonon matrix
  ! in the irreducible Brillouin zone and
  ! expand that to whole BZ.
  !
  USE ener, ONLY : ef
  USE constants, ONLY : pi, ry_to_cmm1, ry_to_ghz, rytoev
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum, mp_bcast
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_images, ONLY : me_image, nproc_image, intra_image_comm
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE cell_base, ONLY : at, bg
  USE ions_base, ONLY : nat
  USE symm_base, ONLY : s, irt, invs
  USE klist,  ONLY: nks, nkstot
  USE wvfct, ONLY: et, nbnd
  USE qpoint, ONLY : xq, nksq, ikks
  USE dynmat, ONLY : dyn, w2
  USE el_phon, ONLY : el_ph_mat, elph_nbnd_min, elph_nbnd_max, done_elph, gamma_disp
  USE control_lr,  ONLY : lgamma
  USE control_ph, ONLY : current_iq, qplot
  USE modes, ONLY : u, nirr
  USE lr_symm_base, ONLY : minus_q, nsymq, rtau, irotmq
  USE lsda_mod,   ONLY : nspin
  USE ktetra,     ONLY : ntetra, tetra, opt_tetra_dos_t
  USE output, ONLY : fildyn
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: eps = 20_dp/ry_to_cmm1 ! eps = 20 cm^-1, in Ry
  INTEGER :: iq, ipert, jpert, nu, iuelph, nbnd_fs, ntpp, rest, &
  &          tfst, tlst, ios, irr
  REAL(dp) :: dosef(2), lambda(3 * nat), gamma, phase_space
  COMPLEX(dp) :: el_ph_sum (3*nat,3*nat)
  !
  REAL(dp),ALLOCATABLE :: wght(:,:,:), et_col(:,:)
  !
  character(len=80) :: filelph
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  !
  DO irr=1,nirr
     IF (.NOT.done_elph(irr)) RETURN
  ENDDO
  !
  ! Work Sharing
  !
  ntpp = ntetra / nproc_image
  rest = MOD(ntetra, nproc_image)
  IF(me_image < rest) THEN
     tfst = (ntpp + 1) *  me_image + 1
     tlst = (ntpp + 1) * (me_image + 1)
  ELSE
     tfst = ntpp *  me_image + 1  + rest
     tlst = ntpp * (me_image + 1) + rest
  END IF
  !
  IF(lgamma) THEN
     iq = 0
  ELSE
     iq = 1
  END IF
  !
  nbnd_fs = elph_nbnd_max - elph_nbnd_min + 1
  !
  WRITE(stdout,'(a)') ""
  WRITE(stdout,'(a,i10)') "     [elph_tetra]   Lowest band which contains FS : ", elph_nbnd_min
  WRITE(stdout,'(a,i10)') "     [elph_tetra]  Highest band which contains FS : ", elph_nbnd_max
  WRITE(stdout,'(a,i10)') "     [elph_tetra]    # of bands which contains FS : ", nbnd_fs
  !
  ! Collect eigval
  !
  ALLOCATE(wght(nbnd_fs,nbnd_fs,nkstot), et_col(nbnd, nkstot))
  !
  CALL poolcollect(nbnd, nks, et, nkstot, et_col)
  !
  IF(lgamma) THEN
     !
     wght(1:nbnd_fs,1:nbnd_fs,1:nkstot) = 0.0_dp
     phase_space = 0.0_dp
     !
     WRITE(stdout,'(a)') "     ##############  Please Chack  ##############"
     WRITE(stdout,'(a)') "       \lambda_{q \nu} is singular at q = \Gamma"
     WRITE(stdout,'(a)') "       You should use shifted q grid            "
     WRITE(stdout,'(a)') "     ############################################"
     !
  ELSE
     !
     CALL elph_tetra_delta1(nbnd_fs,iq,tfst,tlst,et_col,wght)
     phase_space = SUM(wght(1:nbnd_fs,1:nbnd_fs,1:nkstot))
     CALL poolscatter( nbnd_fs*nbnd_fs, nkstot, wght, nks, wght)
     CALL elph_tetra_average_weight(1,nbnd_fs,wght)
     !
  END IF
  !
  DO jpert = 1, 3 * nat
     DO ipert = 1, 3 * nat
        !
        el_ph_sum (ipert, jpert) = SUM(wght( 1:nbnd_fs, 1:nbnd_fs, ikks(1:nksq)) &
        & *       el_ph_mat(elph_nbnd_min:elph_nbnd_max, elph_nbnd_min:elph_nbnd_max, 1:nksq, jpert) &
        & * CONJG(el_ph_mat(elph_nbnd_min:elph_nbnd_max, elph_nbnd_min:elph_nbnd_max, 1:nksq, ipert)) )
     END DO
  END DO
  !
  CALL mp_sum(el_ph_sum, inter_pool_comm)
  !
  CALL symdyn_munu_new (el_ph_sum, u, xq, s, invs, rtau, irt,  at, bg, nsymq, nat, irotmq, minus_q)
  !
  dosef(1:2) = 0.0_dp
  CALL opt_tetra_dos_t (et_col, nspin, nbnd, nkstot, ef, dosef)
  !
  dosef(1:2) = 0.5_dp * dosef(1:2)
  !
  DO nu = 1, 3 * nat
     !
     lambda(nu) = 0.0_dp
     DO ipert = 1, 3 * nat
        DO jpert = 1, 3 * nat
           lambda(nu) = lambda(nu) &
           & + REAL(CONJG(dyn(jpert, nu)) * el_ph_sum(jpert, ipert) * dyn(ipert, nu), dp)
        END DO
     END DO
     !
     IF(w2(nu) >= 0.0_dp) THEN
        lambda(nu) = lambda(nu) / (2.0_dp * w2(nu) * SUM(dosef(1:2)))
     ELSE
        lambda(nu) = 0.0_dp
     END IF
     !
  END DO
  !
  ! Output
  !
  filelph=TRIM(fildyn)//'.elph.'//TRIM(int_to_char(current_iq))
  !
  ! parallel case: only first node writes
  IF ( ionode ) THEN
     !
     iuelph = find_free_unit()
     OPEN (unit = iuelph, file = TRIM(filelph), status = 'unknown', err = &
          100, iostat = ios)
     REWIND (iuelph)
  ELSE
     iuelph = 0
     !
  END IF
100 CONTINUE
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  CALL errore ('elph_tetra_lambda', 'opening file '//filelph, ABS (ios) )
  !
  IF (ionode) THEN
     WRITE (iuelph, '(3f15.8,i8)') xq, 3 * nat
     WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, 3 * nat)
  ENDIF
  WRITE (stdout, *)
  WRITE (stdout, 9000)
  WRITE (stdout, 9005) SUM(dosef(1:2)), ef * rytoev
  WRITE (stdout, 9006) phase_space
  !
  IF (ionode) THEN
     WRITE (iuelph, 9000)
     WRITE (iuelph, 9005) SUM(dosef(1:2)), ef * rytoev
  ENDIF
  !
  DO nu = 1, 3 * nat
     !
     gamma = lambda(nu) * pi * w2(nu) * SUM(dosef(1:2))
     !
     WRITE (stdout, 9010) nu, lambda(nu), gamma * ry_to_gHz
     IF (ionode) WRITE (iuelph, 9010) nu, lambda(nu), gamma * ry_to_gHz
     IF (qplot) gamma_disp(nu,1,current_iq) = gamma * ry_to_gHz
     !
  END DO
  !
9000 FORMAT(5x,'Tetrahedron method')
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i5,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  IF (ionode) CLOSE (unit = iuelph)
  !
  DEALLOCATE(wght, et_col)
  !
END SUBROUTINE elph_tetra_lambda
!
!-----------------------------------------------------------------------
SUBROUTINE elph_tetra_delta1(nbnd_fs,iq,tfst,tlst,et_col,wght)
  !---------------------------------------------------------------------
  !
  ! This routine computed the weight for the double-delta function.
  !
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum
  USE mp_images, ONLY : intra_image_comm
  USE el_phon, ONLY : elph_nbnd_min
  USE ener,   ONLY : ef
  USE wvfct, ONLY: nbnd
  USE klist,  ONLY: nkstot
  USE lsda_mod,   ONLY : nspin
  USE ktetra, ONLY : tetra, ntetra, nntetra, wlsm
  !
  INTEGER,INTENT(IN) :: nbnd_fs, iq, tfst, tlst
  REAL(dp),INTENT(IN) :: et_col(nbnd, nkstot)
  REAL(dp),INTENT(OUT) :: wght(nbnd_fs, nbnd_fs, nkstot)
  !
  INTEGER :: nspin_lsda, ns, nk, nt, ii, jj, ibnd, itetra(4), ik
  REAL(dp) :: ei0(4,nbnd_fs), ej0(4,nbnd_fs), e(4), a(4,4), V, tsmall(4,4), &
  &           ej1(3,nbnd_fs), w0(nbnd_fs,nbnd_fs,4), w1(nbnd_fs,3)
  !
  wght(1:nbnd_fs,1:nbnd_fs,1:nkstot) = 0.0_dp
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  DO ns = 1, nspin_lsda
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     END IF
     !
     DO nt = tfst, tlst
        !
        ei0(1:4, 1:nbnd_fs) = 0.0_dp
        ej0(1:4, 1:nbnd_fs) = 0.0_dp
        DO ii = 1, nntetra
           !
           ik = tetra(ii, nt) + nk
           DO ibnd = 1, nbnd_fs
              ei0(1:4, ibnd) = ei0(1:4, ibnd) &
              &             + wlsm(1:4,ii) * (et_col(ibnd + elph_nbnd_min - 1, ik)    - ef)
              ej0(1:4, ibnd) = ej0(1:4, ibnd) &
              &             + wlsm(1:4,ii) * (et_col(ibnd + elph_nbnd_min - 1, ik+iq) - ef)
           END DO
           !
        END DO
        !
        w0(1:nbnd_fs,1:nbnd_fs,1:4) = 0.0_dp
        !
        DO ibnd = 1, nbnd_fs
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4, ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(1:4,ii) = (0.0_dp - e(ii)) / (e(1:4) - e(ii))
           END DO
           !
           IF(e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) THEN
              !
              ! A
              !
              !V = 3.0_dp * a(2,1) * a(3,1) * a(4,1) / (0.0_dp - e(1))
              V = 3.0_dp * a(2,1) * a(3,1)           / (e(4) - e(1))
              !
              tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              !
              ej1(1:3,1:nbnd_fs) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs))
              !
              CALL elph_tetra_delta2(nbnd_fs,ej1,w1)
              !
              w0(1:nbnd_fs,ibnd,itetra(1:4)) = w0(1:nbnd_fs,ibnd,itetra(1:4)) &
              &           + V * MATMUL(w1(1:nbnd_fs,1:3), tsmall(1:3,1:4))
              !
           ELSE IF( e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) THEN
              !
              ! B - 1
              !
              !V = 3.0_dp * a(3,1) * a(4,1) * a(2,4) / (0.0_dp - e(1))
              V = 3.0_dp           * a(4,1) * a(2,4) / (e(3) - e(1))
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              !
              ej1(1:3,1:nbnd_fs) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs))
              !
              CALL elph_tetra_delta2(nbnd_fs,ej1,w1)
              !
              w0(1:nbnd_fs,ibnd,itetra(1:4)) = w0(1:nbnd_fs,ibnd,itetra(1:4)) &
              &           + V * MATMUL(w1(1:nbnd_fs,1:3), tsmall(1:3,1:4))
              !
              ! B - 2
              !
              !V = 3.0_dp * a(2,3) * a(3,1) * a(4,2) / (0.0_dp - e(1))
              V = 3.0_dp * a(2,3)           * a(4,2) / (e(3) - e(1))
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
              tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              !
              ej1(1:3,1:nbnd_fs) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs))
              !
              CALL elph_tetra_delta2(nbnd_fs,ej1,w1)
              !
              w0(1:nbnd_fs,ibnd,itetra(1:4)) = w0(1:nbnd_fs,ibnd,itetra(1:4)) &
              &           + V * MATMUL(w1(1:nbnd_fs,1:3), tsmall(1:3,1:4))
              !
           ELSE IF(e(3) < 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C
              !
              !V = 3.0_dp * a(1,4) * a(2,4) * a(3,4) / (e(4) - 0.0_dp)
              V = 3.0_dp * a(1,4) * a(2,4)           / (e(4) - e(3))
              !
              tsmall(1, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
              !
              ej1(1:3,1:nbnd_fs) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs))
              !
              CALL elph_tetra_delta2(nbnd_fs,ej1,w1)
              !
              w0(1:nbnd_fs,ibnd,itetra(1:4)) = w0(1:nbnd_fs,ibnd,itetra(1:4)) &
              &           + V * MATMUL(w1(1:nbnd_fs,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        END DO
        !
        DO ii = 1, nntetra
           !
           ik = tetra(ii, nt) + nk
           DO jj = 1, 4
              wght(1:nbnd_fs,1:nbnd_fs,ik) = wght(1:nbnd_fs,1:nbnd_fs,ik) &
              &             + wlsm(jj,ii) * w0(1:nbnd_fs, 1:nbnd_fs, jj)
           END DO
           !
        END DO ! ii
        !
     END DO ! nt
     !
  END DO ! ns
  !
  wght(1:nbnd_fs,1:nbnd_fs,1:nkstot) = wght(1:nbnd_fs,1:nbnd_fs,1:nkstot) &
       / REAL(ntetra, dp)
  IF(nspin == 1) wght(1:nbnd_fs,1:nbnd_fs,1:nkstot) = 2.0_dp * &
                     wght(1:nbnd_fs,1:nbnd_fs,1:nkstot)
  !
  CALL mp_sum(wght, intra_image_comm)
  !
END SUBROUTINE elph_tetra_delta1
!
!-------------------------------------------------------------------------------
SUBROUTINE elph_tetra_delta2(nbnd_fs,ej0,w)
  !-----------------------------------------------------------------------------
  !
  ! 2nd step of tetrahedra method.
  !
  USE kinds, ONLY : dp
  !
  INTEGER,INTENT(IN) :: nbnd_fs
  REAL(dp),INTENT(IN) :: ej0(3,nbnd_fs)
  REAL(dp),INTENT(OUT) :: w(nbnd_fs,3)
  !
  INTEGER :: ibnd, itetra(3), ii
  REAL(dp) :: e(3), a(3,3), V
  !
  w(1:nbnd_fs, 1:3) = 0.0_dp
  !
  DO ibnd = 1, nbnd_fs
     !
     IF(MAXVAL(ABS(ej0(ibnd,1:3))) < 1e-10_dp) &
     & CALL errore("elph_tetra_delta2", "Nesting occurs.", ibnd)
     !
     itetra(1) = 0
     e(1:3) = ej0(1:3,ibnd)
     call hpsort (3, e, itetra)
     !
     DO ii = 1, 3
        a(1:3,ii) = (0.0_dp - e(ii)) / (e(1:3) - e(ii))
     END DO
     !
     IF((e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) .OR. (e(1) <= 0.0_dp .AND. 0.0_dp < e(2))) THEN
        !
        !V = a(2,1) * a(3,1) / (0.0_dp - e(1)) 
        V = a(2,1)           / (e(3) - e(1)) 
        !
        w(ibnd,itetra(1)) = V * (a(1,2) + a(1,3))
        w(ibnd,itetra(2)) = V * a(2,1)
        w(ibnd,itetra(3)) = V * a(3,1)
        !
     ELSE IF((e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) .OR. (e(2) < 0.0_dp .AND. 0.0_dp <= e(3))) THEN
        !
        !V = a(1,3) * a(2,3) / (e(3) - 0.0_dp) 
        V = a(1,3)           / (e(3) - e(2)) 
        !
        w(ibnd,itetra(1)) = V * a(1,3)
        w(ibnd,itetra(2)) = V * a(2,3)
        w(ibnd,itetra(3)) = V * (a(3,1) + a(3,2))
        !
     END IF
     !
  END DO ! ib
  !
END SUBROUTINE elph_tetra_delta2
!
!----------------------------------------------------------------------------
SUBROUTINE elph_tetra_gamma()
  !--------------------------------------------------------------------------
  !
  ! This routine computes the electron-phonon matrix
  ! in the irreducible Brillouin zone and
  ! expand that to whole BZ.
  !
  USE ener, ONLY : ef
  USE constants, ONLY : pi, ry_to_cmm1, ry_to_ghz, rytoev
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum, mp_bcast
  USE mp_pools, ONLY : inter_pool_comm
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE cell_base, ONLY : at, bg
  USE ions_base, ONLY : nat
  USE symm_base, ONLY : s, irt, invs
  USE klist,  ONLY: nks, nkstot
  USE wvfct, ONLY: et, nbnd
  USE qpoint, ONLY : xq, nksq, ikks
  USE dynmat, ONLY : dyn, w2
  USE el_phon, ONLY : el_ph_mat, elph_nbnd_min, elph_nbnd_max, done_elph, gamma_disp
  USE control_lr,  ONLY : lgamma
  USE control_ph, ONLY : current_iq, qplot
  USE modes, ONLY : u, nirr
  USE lr_symm_base, ONLY : minus_q, nsymq, rtau, irotmq
  USE lsda_mod,   ONLY : nspin
  USE ktetra,     ONLY : ntetra, tetra, opt_tetra_dos_t
  USE mp_images, ONLY : me_image, nproc_image, intra_image_comm
  USE output, ONLY : fildyn
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: eps = 20_dp/ry_to_cmm1 ! eps = 20 cm^-1, in Ry
  INTEGER :: ipert, jpert, iuelph, nu, nbnd_fs, iq, ntpp, rest, &
  &          tfst, tlst, ios, irr
  REAL(dp) :: dosef(2), lambda(3 * nat), gamma, phase_space(3*nat)
  COMPLEX(dp) :: el_ph_sum (3*nat,3*nat,3*nat)
  !
  REAL(dp),ALLOCATABLE :: wght(:,:,:,:), et_col(:,:)
  !
  character(len=80) :: filelph
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  !
  DO irr=1,nirr
     IF (.NOT.done_elph(irr)) RETURN
  ENDDO
  !
  ! Work Sharing
  !
  ntpp = ntetra / nproc_image
  rest = MOD(ntetra, nproc_image)
  IF(me_image < rest) THEN
     tfst = (ntpp + 1) *  me_image + 1
     tlst = (ntpp + 1) * (me_image + 1)
  ELSE
     tfst = ntpp *  me_image + 1  + rest
     tlst = ntpp * (me_image + 1) + rest
  END IF
  !
  IF(lgamma) THEN
     iq = 0
  ELSE
     iq = 1
  END IF
  !
  nbnd_fs = elph_nbnd_max - elph_nbnd_min + 1
  !
  WRITE(stdout,'(a)') ""
  WRITE(stdout,'(a,i10)') "     [elph_tetra]   Lowest band which contains FS : ", elph_nbnd_min
  WRITE(stdout,'(a,i10)') "     [elph_tetra]  Highest band which contains FS : ", elph_nbnd_max
  WRITE(stdout,'(a,i10)') "     [elph_tetra]    # of bands which contains FS : ", nbnd_fs
  !
  ! Collect eigval
  !
  ALLOCATE(wght(3 * nat,nbnd_fs,nbnd_fs,nkstot), et_col(nbnd, nkstot))
  !
  CALL poolcollect(nbnd, nks, et, nkstot, et_col)
  !
  CALL elph_tetra_step1(nbnd_fs,iq,tfst,tlst,et_col,wght)
  DO nu = 1, 3* nat
     IF(w2(nu) > 0.0_dp) THEN
        phase_space(nu) = SUM(wght(nu,1:nbnd_fs,1:nbnd_fs,1:nkstot)) / (REAL(3 * nat, dp) * SQRT(w2(nu)))
     ELSE
        phase_space(nu) = 0.0_dp
        wght(nu,1:nbnd_fs,1:nbnd_fs,1:nkstot) = 0.0_dp
     END IF
  END DO
  CALL poolscatter( 3*nat*nbnd_fs*nbnd_fs, nkstot, wght, nks, wght)
  CALL elph_tetra_average_weight(3 * nat,nbnd_fs,wght)
  !
  DO jpert = 1, 3 * nat
     DO ipert = 1, 3 * nat
        DO nu = 1, 3 * nat
           el_ph_sum(ipert, jpert, nu) = SUM(wght(nu, 1:nbnd_fs, 1:nbnd_fs, ikks(1:nksq)) &
           &  * CONJG( el_ph_mat(elph_nbnd_min:elph_nbnd_max, elph_nbnd_min:elph_nbnd_max, 1:nksq, ipert)) &
           &         * el_ph_mat(elph_nbnd_min:elph_nbnd_max, elph_nbnd_min:elph_nbnd_max, 1:nksq, jpert)  )
        END DO
     END DO
  END DO
  !
  CALL mp_sum(el_ph_sum, inter_pool_comm)
  !
  DO nu = 1, 3 * nat
     CALL symdyn_munu_new (el_ph_sum(1:3 * nat,1:3 * nat, nu), u, xq, s, invs, &
     &                     rtau, irt,  at, bg, nsymq, nat, irotmq, minus_q)
  END DO
  !
  dosef(1:2) = 0.0_dp
  CALL opt_tetra_dos_t (et_col, nspin, nbnd, nkstot, ef, dosef)
  !
  dosef(1:2) = 0.5_dp * dosef(1:2)
  !
  DO nu = 1, 3 * nat
     !
     lambda(nu) = 0.0_dp
     !
     DO ipert = 1, 3 * nat
        DO jpert = 1, 3 * nat
           lambda(nu) = lambda(nu) &
           &  + REAL(CONJG(dyn(jpert,nu)) * el_ph_sum(jpert,ipert,nu) * dyn(ipert,nu), dp)
        END DO
     END DO
     !
     IF(w2(nu) >= 0.0_dp) THEN
        lambda(nu) = lambda(nu) / (2.0_dp * w2(nu) * SQRT(w2(nu)) * SUM(dosef(1:2)))
     ELSE
        lambda(nu) = 0.0_dp
     END IF
     !
  END DO
  !
  ! Output
  !
  filelph=TRIM(fildyn)//'.elph.'//TRIM(int_to_char(current_iq))
  !
  ! parallel case: only first node writes
  IF ( ionode ) THEN
     !
     iuelph = find_free_unit()
     OPEN (unit = iuelph, file = TRIM(filelph), status = 'unknown', err = &
          100, iostat = ios)
     REWIND (iuelph)
  ELSE
     iuelph = 0
     !
  END IF
100 CONTINUE
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  CALL errore ('elph_tetra_gamma', 'opening file '//filelph, ABS (ios) )
  !
  IF (ionode) THEN
     WRITE (iuelph, '(3f15.8,i8)') xq, 3 * nat
     WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, 3 * nat)
  ENDIF
  WRITE (stdout, *)
  WRITE (stdout, 9000)
  WRITE (stdout, 9005) SUM(dosef(1:2)), ef * rytoev
  WRITE (stdout, 9006) phase_space
  !
  IF (ionode) THEN
     WRITE (iuelph, 9000)
     WRITE (iuelph, 9005) SUM(dosef(1:2)), ef * rytoev
  ENDIF
  !
  DO nu = 1, 3 * nat
     !
     gamma = lambda(nu) * pi * w2(nu) * SUM(dosef(1:2))
     !
     WRITE (stdout, 9010) nu, lambda(nu), gamma * ry_to_gHz
     IF (ionode) WRITE (iuelph, 9010) nu, lambda(nu), gamma * ry_to_gHz
     IF (qplot) gamma_disp(nu,1,current_iq) = gamma * ry_to_gHz
     !
  END DO
  !
9000 FORMAT(5x,'Tetrahedron method')
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i5,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  IF (ionode) CLOSE (unit = iuelph)
  !
  DEALLOCATE(wght, et_col)
  !
END SUBROUTINE elph_tetra_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE elph_tetra_step1(nbnd_fs,iq,tfst,tlst,et_col,wght)
  !---------------------------------------------------------------------
  !
  ! This routine computed the weight for the double-delta function.
  !
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum
  USE mp_images, ONLY : intra_image_comm
  USE ions_base, ONLY : nat
  USE el_phon, ONLY : elph_nbnd_min
  USE ener,   ONLY : ef
  USE wvfct, ONLY: nbnd
  USE klist,  ONLY: nkstot
  USE lsda_mod,   ONLY : nspin
  USE ktetra, ONLY : ntetra, tetra, nntetra, wlsm
  !
  INTEGER,INTENT(IN) :: nbnd_fs, iq, tfst, tlst
  REAL(dp),INTENT(IN) :: et_col(nbnd, nkstot)
  REAL(dp),INTENT(OUT) :: wght(3 * nat, nbnd_fs, nbnd_fs, nkstot)
  !
  INTEGER :: nt, nspin_lsda, ns, nk, ibnd, ii, jj, itetra(4), ik
  REAL(dp) :: e(4), a(4,4), tsmall(4,4), V, thr = 1e-10_dp, ei0(4,nbnd_fs), ej0(4,nbnd_fs), &
  &           ei1(4), ej1(4,nbnd_fs), w0(3*nat,nbnd_fs,nbnd_fs,4), w1(3*nat,nbnd_fs,4)
  !
  wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,1:nkstot) = 0.0_dp
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  DO ns = 1, nspin_lsda
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     END IF
     !
     DO nt = tfst, tlst
        !
        ei0(1:4, 1:nbnd_fs) = 0.0_dp
        ej0(1:4, 1:nbnd_fs) = 0.0_dp
        DO ii = 1, nntetra
           !
           ik = tetra(ii, nt) + nk
           DO ibnd = 1, nbnd_fs
              ei0(1:4, ibnd) = ei0(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd + elph_nbnd_min - 1, ik)    - ef)
              ej0(1:4, ibnd) = ej0(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd + elph_nbnd_min - 1, ik+iq) - ef)
           END DO
           !
        END DO
        !
        w0(1:3 * nat,1:nbnd_fs,1:nbnd_fs,1:4) = 0.0_dp
        !
        DO ibnd = 1, nbnd_fs
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4, ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(1:4,ii) = (0.0_dp - e(ii) ) / (e(1:4) - e(ii) )
           END DO
           !
           IF( e(1) <= 0.0_dp .AND. 0.0_dp < e(2) ) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(4, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
           ELSE IF( e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
           ELSE IF( e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
                 !
                 CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
                 !
                 DO ii = 1, 3 * nat
                    w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                    &           + V * MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
                 END DO
                 !
              END IF
              !
           ELSE IF( e(4) <= 0.0_dp ) THEN
              !
              ! D - 1
              !
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ej1(1:4,1:nbnd_fs) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs))
              !
              CALL elph_tetra_step2(nbnd_fs,ei1,ej1,w1)
              !
              DO ii = 1, 3 * nat
                 w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) = w0(ii,1:nbnd_fs,ibnd,itetra(1:4)) &
                 &           + MATMUL(w1(ii,1:nbnd_fs,1:4), tsmall(1:4,1:4))
              END DO
              !
           END IF
           !
        END DO ! ibnd
        !
        DO ii = 1, nntetra
           !
           ik = tetra(ii, nt) + nk
           !
           DO jj = 1, 4
              wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,ik) = wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,ik) &
              &           + wlsm(jj,ii) * w0(1:3 * nat,1:nbnd_fs,1:nbnd_fs,jj)
           END DO
           !
        END DO ! ii
        !
     END DO ! nt
     !
  END DO ! ns
  !
  wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,1:nkstot) = wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,1:nkstot) / REAL(ntetra, dp)
  IF(nspin == 1) wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,1:nkstot) = wght(1:3 * nat,1:nbnd_fs,1:nbnd_fs,1:nkstot) * 2.0_dp
  !
  CALL mp_sum(wght, intra_image_comm)
  !
END SUBROUTINE elph_tetra_step1
!
!------------------------------------------------------------------------------------------------
SUBROUTINE elph_tetra_step2(nbnd_fs,ei0,ej0,w)
  !------------------------------------------------------------------------------------------------
  !
  ! This routine compute the second step function in the gamma
  !
  USE kinds, ONLY : dp
  USE ions_base, ONLY : nat
  !
  INTEGER,INTENT(IN) :: nbnd_fs
  REAL(dp),INTENT(IN) :: ei0(4), ej0(4,nbnd_fs)
  REAL(dp),INTENT(OUT) :: w(3*nat,nbnd_fs,4)
  !
  INTEGER :: ibnd, itetra(4), ii
  REAL(dp) :: e(4), a(4,4), de(4), w1(3*nat,4), tsmall(4,4), V, thr = 1e-8_dp
  !
  w(1:3*nat,1:nbnd_fs,1:4) = 0.0_dp
  !
  DO ibnd = 1, nbnd_fs
     !
     e(1:4) = ej0(1:4,ibnd)
     itetra(1) = 0
     call hpsort (4, e, itetra)
     !
     DO ii = 1, 4
        a(1:4,ii) = (0.0_dp - e(ii) ) / (e(1:4) - e(ii))
     END DO
     !
     IF(0.0_dp <= e(1)) THEN
        !
        ! A - 1
        !
        de(1:4) = e(1:4) - ei0(itetra(1:4))
        !
        CALL elph_tetra_delta3(nbnd_fs,de,w1)
        !
        w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) + w1(1:3*nat,1:4)
        !
     ELSE IF((e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) .OR. (e(1) <= 0.0_dp .AND. 0.0_dp < e(2))) THEN
        !
        ! B - 1
        !
        V = a(1,2)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
           tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !       
        END IF
        !
        ! B - 2
        !
        V = a(1,3) * a(2,1)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
           tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        ! B - 3
        !
        V = a(1,4) * a(2,1) * a(3,1)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
           tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !       
        END IF
        !          
     ELSE IF((e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) .OR. (e(2) <= 0.0_dp .AND. 0.0_dp < e(3))) THEN
        !          
        ! C - 1
        !
        V = a(2,4) * a(1,4) * a(3,1)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
           tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !      
        END IF
        !
        ! C - 2
        !
        V = a(1,3) * a(2,3)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        ! C - 3
        ! 
        V = a(1,3) * a(2,4) * a(3,2)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !
        END IF
        !          
     ELSE IF((e(3) < 0.0_dp .AND. 0.0_dp <= e(4)) .OR. (e(3) <= 0.0_dp .AND. 0.0_dp < e(4))) THEN
        !
        ! D - 1
        !
        V = a(3,4) * a(2,4) * a(1,4) 
        !          
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
           tsmall(2, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4) - ei0(itetra(1:4)))
           !
           CALL elph_tetra_delta3(nbnd_fs,de,w1)
           !
           w(1:3*nat,ibnd,itetra(1:4)) = w(1:3*nat,ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:3*nat,1:4), tsmall(1:4,1:4))
           !        
        END IF
        !
     END IF
     !
  END DO
  !
END SUBROUTINE elph_tetra_step2
!
!-------------------------------------------------------------------------------------------
SUBROUTINE elph_tetra_delta3(nbnd_fs,de,w)
  !-------------------------------------------------------------------------------------------
  !
  USE ions_base, ONLY : nat
  USE dynmat, ONLY : w2
  USE kinds, ONLY : dp
  !
  INTEGER,INTENT(IN) :: nbnd_fs
  REAL(dp),INTENT(IN) :: de(4)
  REAL(dp),INTENT(INOUT) :: w(3 * nat,4)
  !
  INTEGER :: nu, ii, itetra(4)
  REAL(dp) :: a(4,4), e(4), V, e0(3 * nat)
  !
  e0(1:3 * nat) = SQRT(ABS(w2(1:3 * nat)))
  !
  e(1:4) = de(1:4)
  itetra(1) = 0
  call hpsort (4, e, itetra)
  !
  w(1:3*nat,1:4) = 0.0_dp
  !
  DO nu = 1, 3 * nat
     !
     DO ii = 1, 4
        a(1:4,ii) = (e0(nu) - e(ii)) / (e(1:4) - e(ii))
     END DO
     !
     IF(e(1) < e0(nu) .AND. e0(nu) <= e(2)) THEN
        !
        V = a(2,1) * a(3,1) * a(4,1) / (e0(nu) - e(1) )
        w(nu,itetra(1)) = a(1,2) + a(1,3) + a(1,4)
        w(nu,itetra(2:4)) = a(2:4,1)
        w(nu,1:4) = w(nu,1:4) * V
        !
     ELSE IF(e(2) < e0(nu) .AND. e0(nu) <= e(3)) THEN
        !
        V = a(2,3) * a(3,1) + a(3,2) * a(2,4)
        w(nu,itetra(1)) = a(1,4) * V + a(1,3) * a(3,1) * a(2,3)
        w(nu,itetra(2)) = a(2,3) * V + a(2,4) * a(2,4) * a(3,2)
        w(nu,itetra(3)) = a(3,2) * V + a(3,1) * a(3,1) * a(2,3)
        w(nu,itetra(4)) = a(4,1) * V + a(4,2) * a(2,4) * a(3,2)
        V  = 1.0_dp / ( e(4) - e(1) )
        w(nu,1:4) = w(nu,1:4) * V
        !
     ELSE IF(e(3) < e0(nu) .AND. e0(nu) < e(4)) THEN
        !
        V = a(1,4) * a(2,4) * a(3,4) / ( e(4) - e0(nu) )
        w(nu,itetra(1:3))  = a(1:3,4)
        w(nu,itetra(4))  = a(4,1) + a(4,2) + a(4,3)
        w(nu,1:4) = w(nu,1:4) * V
        !
     END IF
     !
  END DO ! nu
  !
END SUBROUTINE elph_tetra_delta3
!
!--------------------------------------------------------------------------
SUBROUTINE elph_tetra_average_weight(nmode,nbnd_fs,wght)
  !--------------------------------------------------------------------------
  !
  ! Average weights of degenerated states
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : et
  USE klist, ONLY : nks
  USE qpoint,ONLY : nksq, ikks, ikqs
  USE el_phon, ONLY : elph_nbnd_min
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nmode, nbnd_fs
  REAL(dp),INTENT(INOUT) :: wght(nmode,nbnd_fs,nbnd_fs,nks)
  !
  INTEGER :: ibnd, jbnd, kbnd, ik
  REAL(dp) :: wght2(nmode,nbnd_fs)
  !
  DO ik = 1, nksq
     !
     DO ibnd = 1, nbnd_fs
        !
        wght2(1:nmode,1:nbnd_fs) = wght(1:nmode,1:nbnd_fs,ibnd,ikks(ik))
        !
        DO jbnd = ibnd + 1, nbnd_fs
           !
           IF(ABS(et(ibnd - 1 + elph_nbnd_min,ikks(ik)) - et(jbnd - 1 + elph_nbnd_min,ikks(ik))) < 1e-6_dp) THEN
              wght2(1:nmode,1:nbnd_fs) = wght2(1:nmode,1:nbnd_fs) + wght(1:nmode,1:nbnd_fs,jbnd,ikks(ik))
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 wght(1:nmode,1:nbnd_fs,kbnd,ikks(ik)) = wght2(1:nmode,1:nbnd_fs) / real(jbnd - ibnd, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
        wght2(1:nmode,1:nbnd_fs) = wght(1:nmode,ibnd,1:nbnd_fs,ikks(ik))
        !
        DO jbnd = ibnd + 1, nbnd_fs
           !
           IF(ABS(et(ibnd - 1 + elph_nbnd_min,ikqs(ik)) - et(jbnd - 1 + elph_nbnd_min,ikqs(ik))) < 1e-6_dp) THEN
              wght2(1:nmode,1:nbnd_fs) = wght2(1:nmode,1:nbnd_fs) + wght(1:nmode,jbnd,1:nbnd_fs,ikks(ik))
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 wght(1:nmode,kbnd,1:nbnd_fs,ikks(ik)) = wght2(1:nmode,1:nbnd_fs) / real(jbnd - ibnd, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
END SUBROUTINE elph_tetra_average_weight
!
END MODULE elph_tetra_mod
