  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE elphel2_shuffle( npe, imode0, dvscfins, gmapsym, eigv, isym, xq0, timerev )
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
  USE io_global,        ONLY : stdout
  USE wavefunctions,    ONLY : evc
  USE io_files,         ONLY : diropn
  USE wvfct,            ONLY : npwx
  USE pwcom,            ONLY : current_spin, lsda, nbnd, nks
  USE cell_base,        ONLY : at, bg
  USE klist_epw,        ONLY : xk_loc, xk_all, isk_loc, et_all
  USE cell_base,        ONLY : tpiba
  USE gvect,            ONLY : g
  USE uspp,             ONLY : vkb
  USE symm_base,        ONLY : s, ft
  USE modes,            ONLY : u
  USE qpoint,           ONLY : xq, npwq
  USE eqv,              ONLY : dvpsi
  USE phcom,            ONLY : evq
  USE units_lr,         ONLY : lrwfc, iuwfc
  USE phus,             ONLY : alphap
  USE lrus,             ONLY : becp1
  USE becmod,           ONLY : calbec
  USE elph2,            ONLY : el_ph_mat, igk_k_all, xkq, etq, ngk_all, lower_band, &
                               upper_band, ibndstart, ibndend, nbndep, ngxx, ngxxf, &
                               ng0vec, shift, gmap, g0vec_all_r
  USE fft_base,         ONLY : dffts
  USE constants_epw,    ONLY : czero, cone, ci, zero, eps8
  USE control_flags,    ONLY : iverbosity
  USE klist,            ONLY : nkstot
  USE division,         ONLY : kpointdivision, fkbounds, fkbounds_bnd
  USE kfold,            ONLY : ktokpmq
  USE low_lvl,          ONLY : fractrasl, rotate_cart, s_crystocart
  USE io_epw,           ONLY : readwfc, readkmap
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE dvqpsi,           ONLY : dvqpsi_us3, adddvscf2
  USE uspp_init,        ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npe
  !! Number of perturbations for this irr representation
  INTEGER, INTENT(in) :: imode0
  !! Current mode number
  INTEGER, INTENT(in) :: gmapsym(ngxxf)
  !! Correspondence  G->S(G)
  INTEGER, INTENT(in) :: isym
  !! The symmetry which generates the current q in the star
  REAL(KIND = DP), INTENT(in) :: xq0(3)
  !! The first q-point in the star (cartesian coords.)
  COMPLEX(KIND = DP), INTENT(in) :: dvscfins(dffts%nnr, nspin_mag, npe)
  !! Delta scf potential
  COMPLEX(KIND = DP), INTENT(in) :: eigv (ngxxf)
  !! $e^{iGv}$ for 1...nsym (v the fractional translation)
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
  INTEGER :: nkk
  !! Index of k-point in the pool
  INTEGER :: nkk_abs
  !! Absolute index of k-point
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
  COMPLEX(KIND = DP) :: aux4(npwx * npol, nbnd)
  !! Rotated psi_m,k+q WF by SU(2)
  COMPLEX(KIND = DP) :: aux5(npwx * npol, nbnd)
  !! Rotated psi_nk WF by SU(2)
  COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:, :)
  !! Auxillary wavefunction
  COMPLEX(KIND = DP), ALLOCATABLE :: aux2(:, :)
  !! Auxillary wavefunction
  COMPLEX(KIND = DP), ALLOCATABLE :: aux3(:, :)
  !! Auxillary wavefunction
  COMPLEX(KIND = DP), ALLOCATABLE :: elphmat(:, :, :)
  !! arrays for e-ph matrix elements
  COMPLEX(KIND = DP), EXTERNAL :: zdotc
  !! Important for NAG compiler
  !
  ALLOCATE(elphmat(nbndep, nbndep, npe), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating elphmat', 1)
  ALLOCATE(aux1(dffts%nnr, npol), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating aux1', 1)
  ALLOCATE(aux2(npwx * npol, nbnd), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating aux2', 1)
  elphmat(:, :, :) = czero
  aux1(:, :) = czero
  aux2(:, :) = czero
  zero_vect = zero
  xkq(:, :) = zero
  !
  IF (nproc_pool > 1) THEN
    CALL errore('elphel2_shuffle', 'only one proc per pool in shuffle mode', 1)
  ENDIF
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds(nkstot, lower_bnd, upper_bnd)
  !
  ! SP: Bound for band parallelism
  CALL fkbounds_bnd(nbndep, lower_band, upper_band)
  !
  ALLOCATE(aux3(npwx * npol, lower_band:upper_band), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating aux3', 1)
  ALLOCATE(dvpsi(npwx * npol, lower_band:upper_band), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating dvpsi', 1)
  aux3(:, :) = czero
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
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating etq', 1)
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
    aux4(:, :) = czero
    aux5(:, :) = czero
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
    IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating igk', 1)
    ALLOCATE(igkq(npwq), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error allocating igkq', 1)
    !
    igk = igk_k_all(1:npw, ik + lower_bnd - 1)
    igkq = igk_k_all(1:npwq, nkq_abs)
    !
    IF (nks > 1 .AND. MAXVAL(igkq(1:npwq)) > ngxx) THEN
      CALL errore('elphel2_shuffle', 'ngxx too small', 1)
    ENDIF
    !
    !
    !------------------------------------------------------------
    ! Rotate evc --> axu5 and evq --> aux4 by SU2^{dagger}
    !------------------------------------------------------------
    IF (noncolin) THEN
      DO ibnd = ibndstart, ibndend
        DO ig = 1, npw
          aux5(ig, ibnd) = su2(1, 1) * evc(ig, ibnd) + su2(1, 2) * evc(ig + npwx, ibnd)
          aux5(ig + npwx, ibnd) = su2(2, 1) * evc(ig, ibnd) + su2(2, 2) * evc(ig + npwx, ibnd)
        ENDDO
        DO ig = 1, npwq
          aux4(ig, ibnd) = su2(1, 1) * evq(ig, ibnd) + su2(1, 2) * evq(ig + npwx, ibnd)
          aux4(ig + npwx, ibnd) = su2(2, 1) * evq(ig, ibnd) + su2(2, 2) * evq(ig + npwx, ibnd)
        ENDDO
      ENDDO
    ELSE
      aux5 = evc
      aux4 = evq
    ENDIF
    !
    ! ----------------------------------------------------------------
    ! Set the gauge for the eigenstates: unitary transform and phases
    ! ----------------------------------------------------------------
    !
    ! With this option, different compilers and different machines
    ! should always give the same wavefunctions.
    !
!    CALL ktokpmq(xk_loc(:, ik),  zero_vect, +1, ipool, nkk, nkk_abs)
!    CALL ktokpmq(xkq(:, ik), zero_vect, +1, ipool, nkk, nkq_abs)
    !
!    umat(:, :, ik)  = umat_all(:, :, nkk_abs)
!    umatq(:, :, ik) = umat_all(:, :, nkq_abs)
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
      CALL fractrasl(npw,  igk,  aux5, eigv, cone)
      CALL fractrasl(npwq, igkq, aux4, eigv, cone)
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
    CALL calbec(npw, vkb, aux5, becp1(ik))
    !
    ! we also recompute the derivative of the becp terms with the (rotated) wfs
    !
    DO ipol = 1, 3
      aux2 = czero
      DO ibnd = ibndstart, ibndend
        DO ig = 1, npw
          aux2(ig, ibnd) = aux5(ig,ibnd) * tpiba * ci * (sxk(ipol) + g(ipol,igk(ig)))
        END DO
        IF (noncolin) THEN
          DO ig = 1, npw
            aux2(ig + npwx, ibnd) = aux5(ig + npwx, ibnd) * tpiba * ci * (sxk(ipol) + g(ipol, igk(ig)))
          ENDDO
        ENDIF
      ENDDO
      CALL calbec(npw, vkb, aux2, alphap(ipol, ik))
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
      !  we have to use the first q in the star in the dvqpsi_us3 call below (xq0)
      !
      mode = imode0 + ipert
      IF (timerev) THEN
        CALL dvqpsi_us3(ik, CONJG(u(:, mode)), .FALSE., xkqtmp, xq0, igk, igkq, npw, npwq, aux5)
      ELSE
        CALL dvqpsi_us3(ik, u(:, mode), .FALSE., xkqtmp, xq0, igk, igkq, npw, npwq, aux5)
      ENDIF
      !DBSP
      ! b = b+SUM((REAL(REAL(dvpsi(:, :))))**2)+SUM((REAL(AIMAG(dvpsi(:, :))))**2)
      !
      !  calculate dvscf_q*psi_k
      !
      CALL start_clock('dvscf_q*psi_k')
      !
      aux3 = czero
      DO ibnd = lower_band, upper_band
        CALL invfft_wave(npw, igk, aux5(:,ibnd + ibndstart - 1), aux1)
        IF (timerev) THEN
          CALL apply_dpot(dffts%nnr, aux1, CONJG(dvscfins(:, :, ipert)), current_spin)
        ELSE
          CALL apply_dpot(dffts%nnr, aux1, dvscfins(:, :, ipert), current_spin)
        ENDIF
        CALL fwfft_wave(npwq, igkq, aux3(:, ibnd), aux1)
      ENDDO
      dvpsi = dvpsi + aux3
      !
      !DBSP
      !c = c+SUM((REAL(REAL(dvpsi(:, :))))**2)+SUM((REAL(AIMAG(dvpsi(:, :))))**2)
      !
      CALL adddvscf2(ipert, ik)
      ! DBPS
      !d = c+SUM((REAL(REAL(dvpsi(:, :))))**2)+SUM((REAL(AIMAG(dvpsi(:, :))))**2)
      !
      ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
      !
      !
      DO ibnd =lower_band, upper_band
        DO jbnd = 1, nbndep
          elphmat(jbnd, ibnd, ipert) = ZDOTC(npwq, aux4(1, jbnd + ibndstart - 1), 1, dvpsi(1, ibnd), 1)
          IF (noncolin) THEN
            elphmat(jbnd, ibnd, ipert) = elphmat(jbnd, ibnd, ipert) + &
               ZDOTC(npwq, aux4(npwx + 1, jbnd + ibndstart - 1), 1, dvpsi(npwx + 1, ibnd), 1)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum(elphmat, intra_pool_comm)
    CALL mp_sum(elphmat, inter_image_comm)
    !
!DBSP
!     IF (ik==2) THEN
!       write(*,*)'SUM dvpsi b ', b
!       write(*,*)'SUM dvpsi c ', c
!       write(*,*)'SUM dvpsi d ', d
!       write(*,*)'elphmat(:, :, :)**2', SUM((REAL(REAL(elphmat(:, :, :))))**2)+SUM((REAL(AIMAG(elphmat(:, :, :))))**2)
!     ENDIF
!END
    !! 02/2020
    !! The part below is commented since umat and umatq originating from umat_all in setphases_wrap
    !! are identity matrices.
    !
    !  Rotate elphmat with the gauge matrices (this should be equivalent
    !  to calculate elphmat with the truely rotated eigenstates)
    !
!    DO ipert = 1, npe
!      !
!      ! the two zgemm call perform the following ops:
!      !  elphmat = umat(k+q)^\dagger * [ elphmat * umat(k) ]
!      !
!      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, elphmat(:, :, ipert), &
!                 nbnd, umat(:, :, ik), nbnd, czero, eptmp, nbnd)
!      CALL ZGEMM('c', 'n', nbnd, nbnd, nbnd, cone, umatq(:, :, ik), &
!                 nbnd, eptmp, nbnd, czero, elphmat(:, :, ipert), nbnd)
!      !
!    ENDDO
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
    IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating igk', 1)
    DEALLOCATE(igkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating igkq', 1)
    !
  ENDDO ! ik
  !
  !  restore original configuration of files
  !
  CALL diropn(iuwfc, 'wfc', lrwfc, exst)
  ! never remove this barrier - > insures that wfcs are restored to each pool before moving on
  CALL mp_barrier(world_comm)
  !
  DEALLOCATE(elphmat, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating elphmat', 1)
  DEALLOCATE(aux1, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating aux1', 1)
  DEALLOCATE(aux2, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating aux2', 1)
  DEALLOCATE(aux3, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating aux3', 1)
  DEALLOCATE(dvpsi, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating dvpsi', 1)
  DEALLOCATE(etq, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphel2_shuffle', 'Error deallocating etq', 1)
  !
  !------------------------------------------------------------
  END SUBROUTINE elphel2_shuffle
  !------------------------------------------------------------
