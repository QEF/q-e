!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Author: Mitsuaki Kawamura, U. Tokyo
!----------------------------------------------------------------------------
MODULE elph_scdft_mod
  !--------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC elph_scdft
  !
  CONTAINS
 !
!----------------------------------------------------------------------------
SUBROUTINE elph_scdft()
  !--------------------------------------------------------------------------
  !
  ! This routine computes the electron-phonon matrix
  ! in the irreducible Brillouin zone and
  ! expand that to whole BZ.
  !
  USE kinds, ONLY : dp
  USE parameters,  ONLY : npk
  USE mp, ONLY : mp_sum
  USE mp_images, ONLY : nproc_image, intra_image_comm, me_image
  USE mp_pools, ONLY : npool, inter_pool_comm, my_pool_id, &
  &                    nproc_pool, intra_pool_comm, me_pool
  USE io_global,   ONLY : stdout
  USE cell_base, ONLY : at, bg
  USE ions_base, ONLY : nat
  USE symm_base, ONLY : s, irt
  USE klist,  ONLY: nks, nkstot, xk
  USE noncollin_module, ONLY: nspin_lsda
  USE start_k, ONLY: nk1, nk2, nk3
  USE wvfct, ONLY: nbnd, et
  USE qpoint, ONLY : xq, nksq
  USE dynmat, ONLY : dyn, w2
  USE el_phon, ONLY : el_ph_mat, elph_nbnd_min, elph_nbnd_max
  USE control_lr,  ONLY : lgamma
  USE control_ph, ONLY : current_iq
  USE modes, ONLY : u
  USE lr_symm_base, ONLY : minus_q, nsymq, rtau
  !
  INTEGER :: ik, ik2, ikk, ikq, ib, jb, ii, jj, isym, nu, jsym, elph_unit, &
  &          nn, nn2, nks_real, ikv2(3), nsymqbz, nksq2, nrcv, ipe
  INTEGER :: &
  & nbnd_fs,  & ! = elph_nbnd_max - elph_nbnd_min + 1
  & nkBZ        ! = nk1 * nk2 * nk3
  !
  REAL(dp) :: kv1(3), kv2(3)
  REAL(dp) :: xk_col(3,nkstot)
  COMPLEX(dp) :: epmat(3 * nat, 3 * nat)
  !
  INTEGER,ALLOCATABLE :: symq(:,:,:), eqBZ(:), sBZ(:), indx(:), cnt(:), dsp(:)
  REAL(dp),ALLOCATABLE :: gg_col(:,:), gg(:,:,:,:)
  !
  CHARACTER(100) :: chiq, elphname
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  WRITE(stdout,*) "[elph_scdft]  write elph.dat with symmetries (only on Fermi surfaces)"
  !
  nbnd_fs = elph_nbnd_max - elph_nbnd_min + 1
  !
  WRITE(stdout,*) "[elph_scdft]   Lowest band which contains FS : ", elph_nbnd_min
  WRITE(stdout,*) "[elph_scdft]  Highest band which contains FS : ", elph_nbnd_max
  WRITE(stdout,*) "[elph_scdft]    # of bands which contains FS : ", nbnd_fs
  !
  ! Symmetries of this q
  !
  nsymqbz = nsymq
  IF(minus_q) nsymqbz = nsymq * 2
  !
  ALLOCATE(symq(3,3,nsymqbz))
  !
  symq(1:3,1:3,1:nsymq) = s(1:3,1:3,1:nsymq)
  IF(minus_q) symq(1:3,1:3,nsymq+1:nsymq+nsymq) = - s(1:3,1:3,1:nsymq)
  !
  WRITE(stdout,*) "[elph_scdft]  # of symmetries of this q : ", nsymqbz
  !
  ! Find which k-points of a uniform grid are in the IBZ
  !
  CALL elph_scdft_cnt_and_dsp(inter_pool_comm,npool,my_pool_id,nks,cnt,dsp)
  !
  xk_col(1:3,1:nkstot) = 0_dp
  xk_col(1:3,dsp(my_pool_id) + 1: dsp(my_pool_id) + nks) = xk(1:3,1:nks)
  CALL mp_sum(xk_col, inter_pool_comm)
  !
  nkBZ = nk1 * nk2 * nk3
  nks_real = nkstot / nspin_lsda
  ALLOCATE(eqBZ(nkBZ), sBZ(nkBZ), indx(nkBZ))
  !
  IF(lgamma)THEN
     CALL lint (nsymq, s, minus_q, at, bg, npk, 0, 0, 0, &
     &          nk1, nk2, nk3, nks_real, xk_col, 1, nkBZ, eqBZ, sBZ)
  ELSE
     CALL lint (nsymq, s, minus_q, at, bg, npk, 0, 0, 0, &
     &          nk1, nk2, nk3, nks_real, xk_col, 2, nkBZ, eqBZ, sBZ)
  END IF
  !
  ! Querry of El-Ph matrix
  !
  nksq2 = nksq / nproc_pool
  IF(me_pool < MOD(nksq, nproc_pool)) nksq2 = nksq2 + 1
  !  
  CALL elph_scdft_cnt_and_dsp(intra_pool_comm,nproc_pool,me_pool,nksq2,cnt,dsp)
  !
  indx(1:nkBZ) = 0
  nksq2 = 0
  DO ik = dsp(me_pool) + 1, dsp(me_pool) + cnt(me_pool)
     !
     IF(lgamma) THEN
        ikk = ik
     ELSE
        ikk = 2 * ik - 1
     END IF
     !
     kv1(1:3) = MATMUL(xk(1:3,ikk), at(1:3, 1:3))
     !
     DO isym = 1, nsymqbz
        !
        kv2(1:3) = MATMUL(REAL(symq(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nk1, nk2, nk3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1e-5_dp)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nk1, nk2, nk3/))
        nn2 = 1 + ikv2(1) + nk1 * ikv2(2) + nk1 * nk2 * ikv2(3)
        !
        IF(indx(nn2) /= 0) CYCLE
        !
        nksq2 = nksq2 + 1
        indx(nn2) = nksq2
        !
     END DO
     !
  END DO
  !
  ! Collect El-Ph matrix -> PE0
  !
  ALLOCATE(gg(3 * nat, nbnd_fs, nbnd_fs, nksq2))
  !
  gg(1:3 * nat, 1:nbnd_fs, 1:nbnd_fs, 1:nksq2) = 0_dp
  !
  nksq2 = 0
  indx(1:nkBZ) = 0
  DO ik = dsp(me_pool) + 1, dsp(me_pool) + cnt(me_pool)
     !
     IF(lgamma) THEN
        ikk = ik
     ELSE
        ikk = 2 * ik - 1
     END IF
     !
     kv1(1:3) = MATMUL(xk(1:3,ikk), at(1:3, 1:3))
     !
     DO isym = 1, nsymqbz
        !
        kv2(1:3) = MATMUL(REAL(symq(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nk1, nk2, nk3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1e-5_dp)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nk1, nk2, nk3/))
        nn  = ikv2(3) + ikv2(2) * nk3 + ikv2(1) * nk2 * nk3 + 1
        nn2 = ikv2(1) + ikv2(2) * nk1 + ikv2(3) * nk2 * nk1 + 1
        !
        IF(indx(nn2) /= 0) CYCLE
        !
        jsym = sBZ(nn)
        nksq2 = nksq2 + 1
        indx(nn2) = nksq2
        !
        DO ib = 1, nbnd_fs
           DO jb = 1, nbnd_fs
              !
              DO ii = 1, 3 * nat
                 DO jj = 1, 3 * nat
                    epmat(ii,jj) = CONJG(el_ph_mat(elph_nbnd_min - 1 + jb,elph_nbnd_min - 1 + ib,ik,ii)) &
                    &                  * el_ph_mat(elph_nbnd_min - 1 + jb,elph_nbnd_min - 1 + ib,ik,jj)
                 END DO
              END DO
              !
              CALL symm(epmat, u, xq, s, jsym, rtau, irt, at, bg, nat)
              !
              epmat(1:3 * nat,1:3 * nat) = MATMUL(epmat(1:3 * nat,1:3 * nat), dyn(1:3 * nat,1:3 * nat))
              !
              DO nu = 1, 3 * nat
                 !
                 gg(nu,jb,ib,nksq2) = REAL(DOT_PRODUCT(dyn(1:3 * nat, nu), epmat(1:3 * nat, nu)), dp)
                 !
              END DO ! nu
              !
           END DO ! jb
        END DO ! ib
        !
     END DO ! isym
     !
  END DO ! ik
  !
  DO nu = 1, 3 * nat
     IF(w2(nu) <= 0_dp) THEN
        gg(nu, 1:nbnd_fs, 1:nbnd_fs, 1:nksq2) = 0_dp
     ELSE
        gg(nu, 1:nbnd_fs, 1:nbnd_fs, 1:nksq2) = gg(nu,1:nbnd_fs,1:nbnd_fs,1:nksq2) &
        &                                  / SQRT(w2(nu)) * 0.5_dp
     END IF
  END DO
  !
  ! indices of k; ik is at indx(ik)
  !
  CALL elph_scdft_cnt_and_dsp(intra_image_comm,nproc_image,me_image,nksq2,cnt,dsp)
  !
  DO ik = 1, nkBZ
     IF(indx(ik) /= 0) indx(ik) = indx(ik) + dsp(me_image)
  END DO
  !
  CALL mp_sum(indx, intra_image_comm)
  !
  ! Gather El-Ph matrix inter pool
  !
  IF(me_image == 0) THEN
     nrcv = 3 * nat * nbnd_fs * nbnd_fs * nkBZ
     ALLOCATE(gg_col(3 * nat * nbnd_fs * nbnd_fs, nkBZ))
  ELSE
     nrcv = 1
     ALLOCATE(gg_col(1,1))
  END IF
  !
  CALL elph_scdft_gather_r(gg,3 * nat * nbnd_fs * nbnd_fs * nksq2,gg_col,nrcv, &
  &                me_image, nproc_image, intra_image_comm)
  !
  CALL mp_sum(nksq2, intra_image_comm)
  !
  WRITE(stdout,*) "Verify SUM(nksq2), nkBZ : ", nksq2, nkBZ
  !
  ! Write el-ph to file elph.dat
  !
  IF(me_image == 0) THEN
     !
     elph_unit = find_free_unit()
     !
     WRITE(chiq,*) current_iq
     WRITE(elphname,'(a,a,a)') "elph", TRIM(ADJUSTL(chiq)), ".dat"
     !
     OPEN(elph_unit,file = TRIM(elphname))
     !
     !# of Monkhost-Pack grid
     !
     WRITE(elph_unit,*) nk1, nk2, nk3
     ! 
     !# of band
     !
     WRITE(elph_unit,*) nbnd_fs, elph_nbnd_min
     !
     ! q-vector(Crystal cordinate)
     !
     WRITE(elph_unit,*) MATMUL(xq(1:3), at(1:3, 1:3))
     !
     !# of mode
     !
     WRITE(elph_unit,*) 3 * nat
     !
     ! Frequences[Ryd]
     !
     DO ii = 1, 3 * nat
        WRITE(elph_unit,*) SQRT(ABS(w2(ii)))
     END DO
     WRITE(elph_unit,*) ""
     !
     ! Electron-Phonon matrix
     !
     DO ik = 1, nkBZ
        WRITE(elph_unit,'(3e25.15)') gg_col(1:3 * nat * nbnd_fs * nbnd_fs, indx(ik))
     END DO
     !
     CLOSE(elph_unit)
     !
  END IF ! IF(ionode)
  !
  DEALLOCATE(sBZ, eqBZ, gg, gg_col)
  !
  DEALLOCATE(symq, cnt, dsp)
  !
  CALL elph_scdft_dmuxc()
  !
END SUBROUTINE elph_scdft
!
!----------------------------------------------------------
SUBROUTINE elph_scdft_cnt_and_dsp(comm,npe,mype,nn,cnt,dsp)
  !--------------------------------------------------------
  !
  ! This routine compute coutnt and disples for MPI 
  ! Communicater : comm
  !
  USE mp, ONLY : mp_sum
  !
  INTEGER,INTENT(IN) :: npe, mype, comm, nn
  INTEGER,INTENT(OUT),ALLOCATABLE :: cnt(:), dsp(:)
  !
  INTEGER :: ipe
  !
  IF(ALLOCATED(cnt)) DEALLOCATE(cnt)
  IF(ALLOCATED(dsp)) DEALLOCATE(dsp)
  !
  ALLOCATE(cnt(0:npe-1), dsp(0:npe-1))
  !
  cnt(0:npe - 1) = 0
  cnt(mype) = nn
  CALL mp_sum(cnt, comm)
  !
  dsp(0) = 0
  DO ipe = 1, npe - 1
     dsp(ipe) = dsp(ipe - 1) + cnt(ipe - 1)
  END DO
  !
END SUBROUTINE elph_scdft_cnt_and_dsp
!
!------------------------------------------------------------------------
SUBROUTINE elph_scdft_gather_r(snd,nsnd,rcv,nrcv,mype,npe,comm)
  !----------------------------------------------------------------------
  !
  ! This routine gathers a real matrix to PE 0.
  !
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum, mp_gather
  !
  INTEGER,INTENT(IN) :: nsnd, nrcv, mype, npe, comm
  REAL(dp),INTENT(IN) :: snd(nsnd)
  REAL(dp),INTENT(OUT) :: rcv(nrcv)
  !
  INTEGER :: cnt(0:npe - 1), dsp(0:npe - 1), ipe
  !
  cnt(0:npe - 1) = 0
  cnt(mype) = nsnd
  !
  CALL mp_sum(cnt, comm)
  !
  dsp(0) = 0
  DO ipe = 1, npe - 1
     dsp(ipe) = dsp(ipe - 1) + cnt(ipe - 1)
  END DO
  !
  CALL mp_gather(snd(1:nsnd), rcv(1:nrcv), cnt, dsp, 0, comm)
  !
END SUBROUTINE elph_scdft_gather_r
!
!--------------------------------------------------------
SUBROUTINE elph_scdft_dmuxc()
  !------------------------------------------------------
  !
  ! This routine output the f_{XC} for LDA in G space to a file.
  !
  USE mp,        ONLY : mp_sum, mp_max, mp_min
  USE kinds,     ONLY : dp
  USE gvect,     ONLY : mill
  USE wvfct,     ONLY : npwx
  USE klist,     ONLY : ngk, igk_k
  USE fft_base,  ONLY : dffts
  USE mp_images, ONLY : me_image
  USE mp_pools,  ONLY : intra_pool_comm
  USE eqv,       ONLY : dmuxc
  !
  INTEGER :: gmin(3), gmax(3), ig, dmxc_unit
  COMPLEX(dp) :: dmxc1(dffts%nnr), dmxc2(npwx)
  COMPLEX(dp),ALLOCATABLE :: dmxc3(:,:,:)
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  ! Define whole G-grid
  !
  gmin(1:3) =   100
  gmax(1:3) = - 100
  DO ig = 1, ngk(1)
     gmax(1:3) = max(gmax(1:3), mill(1:3, igk_k(ig,1)))
     gmin(1:3) = min(gmin(1:3), mill(1:3, igk_K(ig,1)))
  END DO
  !
  ! FW_FFT dmuxc_r -> dmuxc_G
  !
  dmxc1(1:dffts%nnr) = cmplx(dmuxc(1:dffts%nnr,1,1), 0_dp)
  CALL elph_scdft_fft(dmxc2(1:npwx), dmxc1(1:dffts%nnr),  ngk(1),  igk_k(:,1),  -1)
  !
  CALL mp_max(gmax, intra_pool_comm)
  CALL mp_min(gmin, intra_pool_comm)
  !
  ALLOCATE(dmxc3(gmin(1):gmax(1), gmin(2):gmax(2), gmin(3):gmax(3)))
  !
  dmxc3(gmin(1):gmax(1),gmin(2):gmax(2),gmin(3):gmax(3)) = cmplx(0_dp,0_dp)
  DO ig = 1, ngk(1)
     dmxc3(mill(1,igk_k(ig,1)), mill(2,igk_k(ig,1)), mill(3,igk_k(ig,1))) = dmxc2(ig)
  END DO
  !
  CALL mp_sum(dmxc3, intra_pool_comm)
  !
  IF(me_image == 0) THEN
     !
     dmxc_unit = find_free_unit()
     !
     OPEN(dmxc_unit, file = "dmuxc.dat")
     !
     WRITE(dmxc_unit,*) gmin(1:3)
     WRITE(dmxc_unit,*) gmax(1:3)
     !
     WRITE(dmxc_unit,*) ""
     WRITE(dmxc_unit,'(2e25.15)') dmxc3(gmin(1):gmax(1),gmin(2):gmax(2),gmin(3):gmax(3))
     !
     CLOSE(dmxc_unit)
     !
  END IF
  !
  DEALLOCATE(dmxc3)
  !
END SUBROUTINE elph_scdft_dmuxc
!
!--------------------------------------------------------------------
SUBROUTINE elph_scdft_fft(evc_g, evc_r, npw, igk, isw)
  !------------------------------------------------------------------
  !
  ! This routine perform inverse-FFT for f_{XC}
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : npwx
  USE fft_base,   ONLY: dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvecs, ONLY : nls
  !
  INTEGER,INTENT(IN) :: isw
  INTEGER,INTENT(IN) :: npw, igk(npw)
  COMPLEX(dp),INTENT(INOUT) :: evc_g(npwx), evc_r(dffts%nnr)
  !
  INTEGER :: ig
  !
  IF (isw.eq.1) THEN
     !
     evc_r = (0_dp, 0_dp)
     !     
     DO ig = 1, npw
        evc_r(nls(igk(ig))) = evc_g(ig)
     END DO
     !
     CALL invfft ('Wave', evc_r(1:dffts%nnr), dffts)
     !
  ELSE IF(isw.eq.-1) THEN
     !
     CALL fwfft ('Wave', evc_r(1:dffts%nnr), dffts)
     !
     DO ig = 1, npw
        evc_g(ig) = evc_r(nls(igk(ig)))
     END DO
     !
  END IF
  !
END SUBROUTINE elph_scdft_fft
!
END MODULE elph_scdft_mod
