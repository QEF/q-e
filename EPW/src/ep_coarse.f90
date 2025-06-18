  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !----------------------------------------------------------------------
  MODULE ep_coarse
  !----------------------------------------------------------------------
  !!
  !! This module contains routines for the EPC computed on the coarse grid.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE dvscf_read(imode0, iq_irr, nqc_irr, npe, xq0, dvscfin, timerev)
    !-----------------------------------------------------------------------
    !!
    !! Electron-phonon calculation from data saved in fildvscf
    !! Shuffle2 mode (shuffle on electrons + load all phonon q's)
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !! Roxana Margine - Jan 2019: Updated based on QE 6.3 for US
    !! HL - Mar 2020: PAW added based on QE v6.5
    !! SP - This routine cannot be placed in io/ due to cyclic dependence.
    !!
    !-----------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE mp_pools,         ONLY : my_pool_id, inter_pool_comm, root_pool
    USE ions_base,        ONLY : nat
    USE gvecs,            ONLY : doublegrid
    USE lrus,             ONLY : int3, int3_nc, int3_paw
    USE uspp,             ONLY : okvan
    USE paw_variables,    ONLY : okpaw
    USE lsda_mod,         ONLY : nspin
    USE fft_base,         ONLY : dfftp, dffts
    USE io,               ONLY : readdvscf, readint3paw
    USE uspp_param,       ONLY : nhm
    USE ep_constants,     ONLY : czero, cone
    USE fft_interfaces,   ONLY : fft_interpolate
    USE noncollin_module, ONLY : nspin_mag, noncolin
    USE dvqpsi,           ONLY : newdq2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: imode0
    !! Current mode number
    INTEGER, INTENT(in) :: iq_irr
    !! Current ireducible q-point
    INTEGER, INTENT(in) :: nqc_irr
    !! Total number of irreducible q-points in the list
    INTEGER, INTENT(in) :: npe
    !! Number of perturbations for this irr representation
    REAL(KIND = DP), INTENT(in) :: xq0(3)
    !! The first q-point in the star (cartesian coords.)
    COMPLEX(KIND = DP), INTENT(inout) :: dvscfin(dfftp%nnr, nspin_mag, npe)
    !! Change of the scf potential
    LOGICAL, INTENT(in) :: timerev
    !!  true if we are using time reversal
    !
    ! Local variables
    INTEGER :: ipert
    !! Change of Vscf due to perturbations
    INTEGER :: is
    !! Counter on spin
    INTEGER :: ierr
    !! Error status
    COMPLEX(KIND = DP), ALLOCATABLE :: dvscfin_tmp(:, :, :)
    !! Change of the scf potential
    !
    CALL start_clock('dvscf_read')
    !
    IF (okvan) THEN
      ALLOCATE(int3(nhm, nhm, nat, nspin_mag, npe), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvscf_read', 'Error allocating int3', 1)
      IF (okpaw) THEN
        ALLOCATE(int3_paw(nhm, nhm, nat, nspin_mag, npe), STAT = ierr)
        IF (ierr /= 0) CALL errore('dvscf_read', 'Error allocating int3_paw', 1)
      ENDIF
      IF (noncolin) THEN
        ALLOCATE(int3_nc(nhm, nhm, nat, nspin, npe), STAT = ierr)
        IF (ierr /= 0) CALL errore('dvscf_read', 'Error allocating int3_nc', 1)
      ENDIF
    ENDIF
    !
    ! read the <prefix>.dvscf_q[iq] files
    !
    dvscfin = czero
    IF (okpaw) int3_paw = czero
    !
    !! 2020.03 HL
    !! The part below should be changed accordingly when parallelization over G vectors is implemented.
    !
    IF (my_pool_id == 0) THEN
      DO ipert = 1, npe
        CALL readdvscf(dvscfin(:, :, ipert), imode0 + ipert, iq_irr, nqc_irr)
        IF (okpaw) CALL readint3paw(int3_paw(:, :, :, :, ipert), imode0 + ipert, iq_irr, nqc_irr)
      ENDDO
    ENDIF
    !
    !! 2020.03 HL
    !! Below root_pool should be interpreted as the index of root pool group.
    !! Currently, there is no root_pool_id like root_bgrp_id in bands groups.
    !
    CALL mp_bcast(dvscfin, root_pool, inter_pool_comm)
    IF (okpaw) CALL mp_bcast(int3_paw, root_pool, inter_pool_comm)
    !
    IF (doublegrid) THEN
      ALLOCATE(dvscfin_tmp(dffts%nnr, nspin_mag, npe) , STAT = ierr)
      IF (ierr /= 0) CALL errore('dvscf_read', 'Error allocating dvscfin_tmp', 1)
      DO is = 1, nspin_mag
        DO ipert = 1, npe
          CALL fft_interpolate(dfftp, dvscfin(:, is, ipert), dffts, dvscfin_tmp(:, is, ipert))
        ENDDO
      ENDDO
      dvscfin(:, :, :) = dvscfin_tmp(:, :, :)
    ENDIF
    !
    CALL newdq2(dvscfin, npe, xq0, timerev)
    !
    IF (doublegrid) THEN
      DEALLOCATE(dvscfin_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvscf_read', 'Error deallocating dvscfin_tmp', 1)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE dvscf_read
    !-----------------------------------------------------------------------
    !---------------------------------------------------------------------
    SUBROUTINE ep_coarse_compute(irr, npe, imode0, dvscfins, gmapsym, eigv, &
                                 isym, xq0, nqc, el_ph_mat, epmatq, timerev)
    !---------------------------------------------------------------------
    !!
    !! Calculation of the electron-phonon matrix elements el_ph_mat
    !! <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
    !!
    !! Written by Feliciano Giustino based on the routine PH/elphon.f90/elphel.
    !! Main difference w.r.t. to original routine is gauge fixing,
    !! shuffle (umklapp) mode and all-q implementation.
    !!
    !! Shuffle mode implemented on may 7 2006
    !!
    !! Nota Bene: this SUBROUTINE is intended only for one proc per pool,
    !! i.e. with no G-vector parallelization (some work on the igkq is
    !! required for that in the g-mapping)
    !!
    !! In order to allow a pool reading the wfc file of another
    !! pool, I had to modify the bound npwx in PW/n_plane_waves.f90
    !! which is now the max across all pools. In this way lrwfc is
    !! the same for all pools.
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !! SP - Nov 2015
    !! We want g(k,Sq) = < k+S(q)(r) | V_S(q)(r) | k(r) >
    !!                 = < k+S(q)(r) | V_q({S|v}^-1 r) | k(r) >
    !!                 = < k+S(q)({S|v}r) | V_q (r) | k({S|v}r) >
    !!
    !! It is important to note that the KB projectors that are applied to the V need
    !! to be computed at (r) and not ({S|v}r). Therefore, for the KB proj (computed in
    !! init_us_2, we need to provide the < Sk+q (r)| and |Sk (r)>.
    !! See Eq. 11.40 and 11.41 of the R. Martin Electronic Structure book.
    !!
    !! Note that in QE Sq is defined as S^-1(q)
    !!
    !! In case of time-reversal
    !! ------------------------
    !!  g(k,-Sq) = < k-S(q)({S|v}r) | V^loc_-q (r) + (V^nloc_q)* | k({S|v}r) >
    !!  where V^loc_{-q} is obtained with setlocq and V^nloc_q = CONGJ(u_pattern)*dvscfins*u_pattern.
    !!  We have to do this splitting because we do not have V^nloc_-q and
    !!  V^loc has to be computed at -q to be mappable with the vkb of the wavefunctions
    !!  computed in init_us_2.
    !!
    !! Roxana Margine - Jan 2019: Updated based on QE 6.3 for US
    !! SP - Sept. 2019: Cleaning
    !!
    !! JL - Feb 2020
    !! Roation of the <k| and <k+q| with SU(2) in the case of noncolin wavefunctions
    !!
    !---------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE mp_global,        ONLY : my_pool_id, nproc_pool, intra_pool_comm, &
                                 inter_pool_comm, inter_image_comm, world_comm
    USE mp,               ONLY : mp_barrier, mp_bcast, mp_sum
    USE wavefunctions,    ONLY : evc
    USE io_files,         ONLY : diropn
    USE wvfct,            ONLY : npwx
    USE pwcom,            ONLY : current_spin, lsda, nbnd, nks
    USE cell_base,        ONLY : at, bg
    USE input,            ONLY : xk_loc, xk_all, isk_loc, et_all
    USE cell_base,        ONLY : tpiba
    USE gvect,            ONLY : g
    USE uspp,             ONLY : vkb
    USE symm_base,        ONLY : s, ft
    USE modes,            ONLY : u, nirr, nmodes
    USE ions_base,        ONLY : nat
    USE qpoint,           ONLY : xq, npwq
    USE eqv,              ONLY : dvpsi
    USE phcom,            ONLY : evq
    USE units_lr,         ONLY : lrwfc, iuwfc
    USE phus,             ONLY : alphap
    USE lrus,             ONLY : becp1
    USE becmod,           ONLY : calbec
    USE global_var,       ONLY : igk_k_all, xkq, etq, ngk_all, lower_band,  &
                                 upper_band, ibndkept, nbndep, ngxx, ngxxf, &
                                 ng0vec, shift, gmap, g0vec_all_r
    USE fft_base,         ONLY : dffts
    USE ep_constants,     ONLY : czero, cone, ci, zero, eps8
    USE klist,            ONLY : nkstot
    USE parallelism,      ONLY : kpointdivision, fkbounds, fkbounds_bnd
    USE kfold,            ONLY : ktokpmq
    USE low_lvl,          ONLY : fractrasl, rotate_cart, s_crystocart
    USE io,               ONLY : readwfc, readkmap
    USE noncollin_module, ONLY : noncolin, npol, nspin_mag
    USE dvqpsi,           ONLY : dvqpsi_us3, adddvscf2
    USE uspp_init,        ONLY : init_us_2
    USE input,            ONLY : nqc1, nqc2, nqc3
    USE lrus,             ONLY : int3, int3_nc, int3_paw
    USE uspp,             ONLY : okvan
    USE paw_variables,    ONLY : okpaw
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: irr
    !! Current irr perturbation
    INTEGER, INTENT(in) :: npe
    !! Number of perturbations for this irr representation
    INTEGER, INTENT(in) :: imode0
    !! Current mode number
    INTEGER, INTENT(in) :: gmapsym(ngxxf)
    !! Correspondence  G->S(G)
    INTEGER, INTENT(in) :: isym
    !! The symmetry which generates the current q in the star
    INTEGER, INTENT(in) :: nqc
    !! Current q-point index
    REAL(KIND = DP), INTENT(in) :: xq0(3)
    !! The first q-point in the star (cartesian coords.)
    COMPLEX(KIND = DP), INTENT(in) :: dvscfins(dffts%nnr, nspin_mag, npe)
    !! Delta scf potential
    COMPLEX(KIND = DP), INTENT(in) :: eigv(ngxxf)
    !! $e^{iGv}$ for 1...nsym (v the fractional translation)
    COMPLEX(KIND = DP), INTENT(inout) :: el_ph_mat(nbndep, nbndep, nks, 3 * nat)
    !! ep matrix element for this q-point
    COMPLEX(KIND = DP), INTENT(inout) :: epmatq(nbndep, nbndep, nks, nmodes, nqc1 * nqc2 * nqc3)
    !! ep-mat on full BZ
    LOGICAL, INTENT(in) :: timerev
    !!  true if we are using time reversal
    !
    ! Local variables
    LOGICAL :: exst
    !! logical variable to check file exists
    INTEGER :: ik
    !! Counter on k-points in the pool
    INTEGER :: ik0
    !! Index of the first k-point block in this pool - 1
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: ipert
    !! Counter on change of Vscf due to perturbations
    INTEGER :: mode
    !! Counter on modes plus pertubations
    INTEGER :: ig
    !! Counter on G-vectors
    INTEGER :: ipooltmp
    !! Index of pool for k
    INTEGER :: ipool
    !! Index of pool for k+q
    INTEGER :: igkq_tmp(npwx)
    !! Correspondence k+q+G <-> G
    INTEGER :: imap
    !! Index in gmap for G-sphere igkq translation with G_0
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: npw
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: nbnd_loc
    !! Local number of bands for band parallelization
    INTEGER :: nkq
    !! Index of k+q-point in the pool
    INTEGER :: nkq_abs
    !! Absolute index of k+q-point
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: igk(:)
    !! Index for k+G
    INTEGER, ALLOCATABLE :: igkq(:)
    !! Index for k+q+G
    ! Local variables for rotating the wavefunctions (in order to use q in the irr wedge)
    REAL(KIND = DP) :: xkqtmp(3)
    !! Temporary k+q vector for KB projectors
    REAL(KIND = DP) :: sxk(3)
    !! Rotated k-point xk
    REAL(KIND = DP) :: zero_vect(3)
    !! Temporary zero vector
    REAL(KIND = DP) :: s_cart(3, 3)
    !! Symmetry matrix in Cartesian basis.
    COMPLEX(KIND = DP) :: su2(2, 2)
    !! SU2 rotation matrix to act WFs
    COMPLEX(KIND = DP) :: su2_no_dagger(2, 2)
    !! SU2 that comes out of find_u, we need the dagger of this, which we assign to su2
    COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:, :)
    !! Auxillary wavefunction
    COMPLEX(KIND = DP), ALLOCATABLE :: aux2(:, :)
    !! Auxillary wavefunction
    COMPLEX(KIND = DP), ALLOCATABLE :: aux4(:, :)
    !! Rotated psi_m,k+q WF by SU(2)
    COMPLEX(KIND = DP), ALLOCATABLE :: aux5(:, :)
    !! Rotated psi_nk WF by SU(2)
    COMPLEX(KIND = DP), ALLOCATABLE :: elphmat(:, :, :)
    !! arrays for e-ph matrix elements
    COMPLEX(KIND = DP), EXTERNAL :: zdotc
    !! Important for NAG compiler
    !
    ALLOCATE(elphmat(nbndep, nbndep, npe), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating elphmat', 1)
    ALLOCATE(aux1(dffts%nnr, npol), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating aux1', 1)
    elphmat(:, :, :) = czero
    aux1(:, :) = czero
    zero_vect = zero
    xkq(:, :) = zero
    !
    IF (nproc_pool > 1) THEN
      CALL errore('ep_coarse_compute', 'only one proc per pool in shuffle mode', 1)
    ENDIF
    !
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)
    !
    ! SP: Bound for band parallelism
    CALL fkbounds_bnd(nbndep, lower_band, upper_band)
    nbnd_loc = upper_band - lower_band + 1
    !
    ALLOCATE(aux2(npwx * npol, nbndep), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating aux2', 1)
    ALLOCATE(aux4(npwx * npol, nbndep), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating aux4', 1)
    ALLOCATE(aux5(npwx * npol, nbndep), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating aux5', 1)
    ALLOCATE(dvpsi(npwx * npol, lower_band:upper_band), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating dvpsi', 1)
    dvpsi(:, :) = czero
    !
    ! setup for k+q folding
    !
    CALL kpointdivision(ik0)
    !
    CALL readkmap(nkstot)
    !
    ! close all sequential files in order to re-open them as direct access
    ! close all .wfc files in order to prepare shuffled read
    !
    CLOSE(iuwfc, STATUS = 'keep')
    ! never remove this barrier
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(etq(nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating etq', 1)
    DO ik = 1, nks
      !
      etq(:, :) = zero
      elphmat(:, :, :) = czero
      IF (lsda) THEN
        current_spin = isk_loc(ik)
      ENDIF
      !DBSP
      !b = zero
      !c = zero
      !d = zero
      !
      ! find index, and possibly pool, of k+q
      ! the index nkq (nkq_abs) takes into account the even/odd ordering
      ! of the nscf calc
      ! we also redefine the ikq points and the corresponding energies
      ! (we need to make sure that xk(:,ikq) is really k+q for the KB projectors
      ! below and also that the eigenvalues are taken correctly in ephwann)
      !
      CALL ktokpmq(xk_loc(:, ik), xq, +1, ipool, nkq, nkq_abs)
      !
      !   we define xkq(:,ik) and etq(:,ik) for the current xq
      !
      xkq(:, ik) = xk_all(:, nkq_abs)
      etq(:, ik) = et_all(:, nkq_abs)
      !
      ipooltmp = my_pool_id + 1
      !
      ! in serial execution ipool is not used in the called subroutines,
      ! in parallel ipooltmp is for k and ipool is for k+q
      !
      ! read unperturbed wavefunctions psi(k) and psi(k+q)
      !
      CALL readwfc(ipooltmp, ik, evc)
      CALL readwfc(ipool, nkq, evq)
      !
      ! --------------------------------------------------
      !  Calculate SO(3) in cartesian from s_axis_to_cart
      !  then calculate SU(2) transformation matrix from
      !  find_u.
      ! --------------------------------------------------
      s_cart(:, :) = zero
      su2_no_dagger(:, :) = czero
      su2(:, :) = czero
      CALL s_crystocart(s(:,:,isym), s_cart, at, bg)
      CALL find_u(s_cart, su2_no_dagger)
      su2 = TRANSPOSE(CONJG(su2_no_dagger))
      !
      ! Now we define the igk and igkq from the global igk_k_all
      !
      npw  = ngk_all(ik + lower_bnd - 1)
      npwq = ngk_all(nkq_abs)
      !
      ALLOCATE(igk(npw), STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating igk', 1)
      ALLOCATE(igkq(npwq), STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error allocating igkq', 1)
      !
      igk = igk_k_all(1:npw, ik + lower_bnd - 1)
      igkq = igk_k_all(1:npwq, nkq_abs)
      !
      IF (nks > 1 .AND. MAXVAL(igkq(1:npwq)) > ngxx) THEN
        CALL errore('ep_coarse_compute', 'ngxx too small', 1)
      ENDIF
      !
      !
      !------------------------------------------------------------
      ! Rotate evc --> axu5 and evq --> aux4 by SU2^{dagger}
      !------------------------------------------------------------
      aux4(:, :) = czero
      aux5(:, :) = czero
      DO ibnd = 1, nbndep
        jbnd = ibndkept(ibnd)
        IF (noncolin) THEN
          DO ig = 1, npw
            aux5(ig, ibnd) = su2(1, 1) * evc(ig, jbnd) + su2(1, 2) * evc(ig + npwx, jbnd)
            aux5(ig + npwx, ibnd) = su2(2, 1) * evc(ig, jbnd) + su2(2, 2) * evc(ig + npwx, jbnd)
          ENDDO
          DO ig = 1, npwq
            aux4(ig, ibnd) = su2(1, 1) * evq(ig, jbnd) + su2(1, 2) * evq(ig + npwx, jbnd)
            aux4(ig + npwx, ibnd) = su2(2, 1) * evq(ig, jbnd) + su2(2, 2) * evq(ig + npwx, jbnd)
          ENDDO
        ELSE
          aux5(1:npw,  ibnd) = evc(1:npw,  jbnd)
          aux4(1:npwq, ibnd) = evq(1:npwq, jbnd)
        ENDIF
      ENDDO
      !
      ! ----------------------------------------------------------------
      ! Set the gauge for the eigenstates: unitary transform and phases
      ! ----------------------------------------------------------------
      !
      ! With this option, different compilers and different machines
      ! should always give the same wavefunctions.
      !
!      CALL ktokpmq(xk_loc(:, ik),  zero_vect, +1, ipool, nkk, nkk_abs)
!      CALL ktokpmq(xkq(:, ik), zero_vect, +1, ipool, nkk, nkq_abs)
      !
!      umat(:, :, ik)  = umat_all(:, :, nkk_abs)
!      umatq(:, :, ik) = umat_all(:, :, nkq_abs)
      !
      ! the k-vector needed for the KB projectors
      xkqtmp = xkq(:, ik)
      !
      ! --------------------------------------------------
      !   Fourier translation of the G-sphere igkq
      ! --------------------------------------------------
      !
      !  Translate by G_0 the G-sphere where evq is defined,
      !  none of the G-points are lost.
      !
      IF (ANY( ABS(g0vec_all_r(:, shift(ik + ik0))) > eps8 )) THEN
        DO ig = 1, npwq
          imap = ng0vec * (igkq(ig) - 1) + shift(ik + ik0)
          igkq_tmp(ig) = gmap(imap)
          !  the old matrix version...
          !  igkq_tmp(ig) = gmap( igkq(ig), shift(ik+ik0) )
        ENDDO
        igkq = igkq_tmp
      ENDIF
      !
      !  find k+q from k+q+G_0
      !  (this is needed in the calculation of the KB terms
      !  for nonlocal pseudos)
      !
      xkqtmp = xkq(:, ik) - g0vec_all_r(:, shift(ik + ik0))
      !
      ! ---------------------------------------------------------------------
      ! phase factor arising from fractional traslations
      ! ---------------------------------------------------------------------
      !
      !  u_{k+q+G_0} carries an additional factor e^{i G_0 v}
      !
      IF (ANY( ABS(ft(:, isym)) > eps8 )) THEN
        CALL fractrasl(nbndep, npw,  igk,  aux5, eigv, cone)
        CALL fractrasl(nbndep, npwq, igkq, aux4, eigv, cone)
      ENDIF
      !
      ! ---------------------------------------------------------------------
      ! wave function rotation to generate matrix elements for the star of q
      ! ---------------------------------------------------------------------
      !
      ! ps. don't use npwx instead of npw, npwq since the unused elements
      ! may be large and blow up gmapsym (personal experience)
      !
      IF (isym /= 1) THEN
         igk(1:npw) = gmapsym(igk(1:npw))
         igkq(1:npwq) = gmapsym(igkq(1:npwq))
      ENDIF
      !
      ! In dvqpsi_us_only3 we need becp1 and alphap for the rotated wfs.
      ! The other quantities (deeq and qq) do not depend on the wfs, in
      ! particular in the KB case (not ultrasoft), the deeq's are the
      ! unscreened coefficients, and the qq's are zero.
      !
      ! For the KB part, remember dV_NL[q_0] ~ |S^-1(k)+q_0> <S^-1(k)|
      ! the total momentum transfer must be q_0 and the rotation
      ! tranforms k+Sq_0 into S^-1(k)+q_0, k into S^-1(k)
      ! [see Eqs. (A9),(A14) Baroni et al. RMP]
      !
      ! Since in QE a normal rotation s is defined as S^-1 we have here
      ! sxk = S(k).
      !
      CALL rotate_cart(xk_loc(:, ik), s(:, :, isym), sxk)
      !
      ! here we generate vkb on the igk() set and for k ...
      CALL init_us_2(npw, igk, sxk, vkb)
      !
      ! ... and we recompute the becp terms with the wfs (rotated through igk)
      !
      CALL calbec(npw, vkb, aux5, becp1(ik), nbndep)
      !
      ! we also recompute the derivative of the becp terms with the (rotated) wfs
      !
      DO ipol = 1, 3
        aux2(:, :) = czero
        DO ibnd = lower_band, upper_band
          DO ig = 1, npw
            aux2(ig, ibnd) = aux5(ig, ibnd) * tpiba * ci * (sxk(ipol) + g(ipol,igk(ig)))
          END DO
          IF (noncolin) THEN
            DO ig = 1, npw
              aux2(ig + npwx, ibnd) = aux5(ig + npwx, ibnd) * tpiba * ci * (sxk(ipol) + g(ipol, igk(ig)))
            ENDDO
          ENDIF
        ENDDO
        CALL calbec(npw, vkb, aux2(:, lower_band:upper_band), alphap(ipol, ik), upper_band - lower_band + 1)
      ENDDO
      !
      ! now we generate vkb on the igkq() set because dvpsi is needed on that set
      ! we need S(k)+q_0 in the KB projector: total momentum transfer must be q_0
      !
      xkqtmp = sxk + xq0
      CALL init_us_2(npwq, igkq, xkqtmp, vkb)
      !
      ! --------------------------------------------------
      !   Calculation of the matrix element
      ! --------------------------------------------------
      !
      DO ipert = 1, npe
        !
        !  recalculate dvbare_q*psi_k
        !  the call to dvqpsi_us3 differs from the old one to dvqpsi_us
        !  only the xkqtmp passed.
        !
        !  The bare potential from nonlinear core correction is not compute in dvqpsi_us3
        !  (because of the .FALSE. parameter) because it is part of the induced potential
        !  dvscfins written by ph.x. This is done so because the bare potential is spin
        !  independent but the potential due to core correction is not.
        !
        !  we have to use the first q in the star in the dvqpsi_us3 call below (xq0)
        !
        mode = imode0 + ipert
        IF (timerev) THEN
          CALL dvqpsi_us3(ik, CONJG(u(:, mode)), .FALSE., xkqtmp, xq0, igk, igkq, npw, npwq, &
                          aux5, CONJG(dvscfins(:, :, ipert)))
        ELSE
          CALL dvqpsi_us3(ik, u(:, mode), .FALSE., xkqtmp, xq0, igk, igkq, npw, npwq, &
                          aux5, dvscfins(:, :, ipert))
        ENDIF
        !DBSP
        ! b = b+SUM((REAL(REAL(dvpsi(:, :))))**2)+SUM((REAL(AIMAG(dvpsi(:, :))))**2)
        !
        CALL adddvscf2(ipert, ik)
        ! DBPS
        ! c = b+SUM((REAL(REAL(dvpsi(:, :))))**2)+SUM((REAL(AIMAG(dvpsi(:, :))))**2)
        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        CALL ZGEMM('C', 'N', nbndep, nbnd_loc, npwx * npol, &
                  (1.d0, 0.d0), aux4(1, 1), npwx * npol, dvpsi(1, lower_band), npwx * npol, &
                  (0.d0, 0.d0), elphmat(1, 1, ipert), nbndep)
        !
      ENDDO
      !
      CALL mp_sum(elphmat, intra_pool_comm)
      CALL mp_sum(elphmat, inter_image_comm)
      !
      !DBSP
      !IF (ik==2) THEN
      !  write(*,*)'SUM dvpsi b ', b
      !  write(*,*)'SUM dvpsi c ', c
      !  write(*,*)'elphmat(:, :, :)**2', SUM((REAL(REAL(elphmat(:, :, :))))**2)+SUM((REAL(AIMAG(elphmat(:, :, :))))**2)
      !ENDIF
      !
      ! 02/2020
      ! The part below is commented since umat and umatq originating from umat_all in setphases_wrap
      ! are identity matrices.
      !
      ! Rotate elphmat with the gauge matrices (this should be equivalent
      ! to calculate elphmat with the truely rotated eigenstates)
      !
      ! DO ipert = 1, npe
      !   !
      !   ! the two zgemm call perform the following ops:
      !   !  elphmat = umat(k+q)^\dagger * [ elphmat * umat(k) ]
      !   !
      !   CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, elphmat(:, :, ipert), &
      !              nbnd, umat(:, :, ik), nbnd, czero, eptmp, nbnd)
      !   CALL ZGEMM('c', 'n', nbnd, nbnd, nbnd, cone, umatq(:, :, ik), &
      !              nbnd, eptmp, nbnd, czero, elphmat(:, :, ipert), nbnd)
      !   !
      ! ENDDO
      !
      !  save eph matrix elements into el_ph_mat
      !
      DO ipert = 1, npe
        DO jbnd = 1, nbndep
          DO ibnd = 1, nbndep
            el_ph_mat(ibnd, jbnd, ik, ipert + imode0) = elphmat(ibnd, jbnd, ipert)
          ENDDO
        ENDDO
      ENDDO
      !
      DEALLOCATE(igk, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating igk', 1)
      DEALLOCATE(igkq, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating igkq', 1)
      !
    ENDDO ! ik
    !
    ! When all irr perturbations have been done, rotate from pattern to Cartesian basis
    IF (irr == nirr) THEN
      ! epmat_{CART} = conjg ( U ) * epmat_{PATTERN}
      ! note it is not U^\dagger but u_pattern!
      ! Have a look to symdyn_munu.f90 for comparison
      !
      DO ibnd = 1, nbndep
        DO jbnd = 1, nbndep
          DO ik = 1, nks
            !
            ! Here is where we calculate epmatq, it appears to be
            ! epmatq = cone * conjug(u) * el_ph_mat + czero
            IF (timerev) THEN
              CALL ZGEMV('n', nmodes, nmodes, cone, u, nmodes, &
                el_ph_mat(ibnd, jbnd, ik, :), 1, czero, epmatq(ibnd, jbnd, ik, :, nqc), 1)
            ELSE
              CALL ZGEMV('n', nmodes, nmodes, cone, CONJG(u), nmodes, &
                el_ph_mat(ibnd, jbnd, ik, :), 1, czero, epmatq(ibnd, jbnd, ik, :, nqc), 1)
            ENDIF
            !
          ENDDO ! ik
        ENDDO ! jbnd
      ENDDO ! ibnd
      !
    ENDIF ! irr == nirr
    !
    !  restore original configuration of files
    !
    CALL diropn(iuwfc, 'wfc', lrwfc, exst)
    ! never remove this barrier - > insures that wfcs are restored to each pool before moving on
    CALL mp_barrier(world_comm)
    !
    IF (okvan) THEN
      DEALLOCATE(int3, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating int3', 1)
      IF (okpaw) THEN
        DEALLOCATE(int3_paw, STAT = ierr)
        IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating int3_paw', 1)
      ENDIF
      IF (noncolin) THEN
        DEALLOCATE(int3_nc, STAT = ierr)
        IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating int3_nc', 1)
      ENDIF
    ENDIF
    !
    DEALLOCATE(elphmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating elphmat', 1)
    DEALLOCATE(aux1, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating aux1', 1)
    DEALLOCATE(aux2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating aux2', 1)
    DEALLOCATE(aux4, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating aux4', 1)
    DEALLOCATE(aux5, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating aux5', 1)
    DEALLOCATE(dvpsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating dvpsi', 1)
    DEALLOCATE(etq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_compute', 'Error deallocating etq', 1)
    !
    !------------------------------------------------------------
    END SUBROUTINE ep_coarse_compute
    !------------------------------------------------------------
  !-----------------------------------------------------------------------------
  END MODULE ep_coarse
  !-----------------------------------------------------------------------------
