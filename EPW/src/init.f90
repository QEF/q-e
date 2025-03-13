  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from the code PH/phq_init - Quantum-ESPRESSO group
  !----------------------------------------------------------------------------
  SUBROUTINE init(first_run)
  !----------------------------------------------------------------------------
  !
  !! This initialization is done nqc_irr times from elphon_shuffle_wrap
  !! not all of the following code is necessary.  More adaptation from
  !! phq_init is needed
  !!
  !! Roxana Margine - Dec 2018: Updated based on QE 6.3
  !!
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ntyp => nsp, tau
  USE becmod,           ONLY : allocate_bec_type
  USE lrus,             ONLY : becp1
  USE uspp,             ONLY : nlcc_any, okvan, nkb
  USE pwcom,            ONLY : npwx, nbnd, nks
  USE constants,        ONLY : tpi
  USE ep_constants,     ONLY : zero, czero
  USE klist,            ONLY : ngk, igk_k, nkstot
  USE gvect,            ONLY : ngm
  USE noncollin_module, ONLY : noncolin, nspin_mag, lspinorb
  USE uspp_param,       ONLY : nhm
  USE phcom,            ONLY : vlocq
  USE qpoint,           ONLY : xq, eigqts
  USE nlcc_ph,          ONLY : drc
  USE global_var,       ONLY : igk_k_all, ngk_all, ngxx, veff, ig_s, ig_e
  USE mp,               ONLY : mp_barrier
  USE lsda_mod,         ONLY : nspin
  USE phus,             ONLY : int1, int1_nc, int2, int2_so, alphap
  USE parallelism,      ONLY : poolgather_int, poolgather_int1
  USE dvqpsi,           ONLY : dvanqq2
  USE scf,              ONLY : v, vltot
  USE fft_base,         ONLY : dfftp
  USE fft_interfaces,   ONLY : fwfft
  USE mp_images,        ONLY : nproc_image, me_image
  ! --------------------------------------------------------------------------------
  ! Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
  ! Shell implementation for future use.
  !USE bzgrid,           ONLY : loadqmesh_serial, loadkmesh_para
  ! --------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: first_run
  !
  ! Local variables
  INTEGER :: ik
  !! counter on k points
  INTEGER :: ipol
  !! counter on polarizations
  INTEGER :: na
  !! counter on atoms
  INTEGER :: ir
  !! counter on FFT mesh points
  INTEGER :: is
  !! counter on spin
  INTEGER :: ierr
  !! Error status
  INTEGER :: itmp
  !! Temporary variable
  INTEGER, EXTERNAL :: ldim_block, gind_block
  REAL(KIND = DP) :: arg
  !! the argument of the phase
  !
  CALL start_clock('init')
  !
  IF (first_run) THEN
    ALLOCATE(vlocq(ngm, ntyp), STAT = ierr)
    IF (ierr /= 0) CALL errore('init', 'Error allocating vlocq', 1)
    ALLOCATE(eigqts(nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('init', 'Error allocating eigqts', 1)
    IF (okvan) THEN
      ALLOCATE(veff(dfftp%nnr, nspin_mag), STAT = ierr)
      IF (ierr /= 0) CALL errore('init', 'Error allocating veff', 1)
      !
      !   we start by computing the FT of the effective potential
      !
      veff(:, :) = czero
      !
      DO is = 1, nspin_mag
        IF (nspin_mag /= 4 .OR. is == 1) THEN
          DO ir = 1, dfftp%nnr
            veff(ir, is) = CMPLX(vltot(ir) + v%of_r(ir, is), zero, KIND = DP)
          ENDDO
        ELSE
          DO ir = 1, dfftp%nnr
            veff(ir, is) = CMPLX(v%of_r(ir, is), zero, KIND = DP)
          ENDDO
        ENDIF
        CALL fwfft('Rho', veff(:, is), dfftp)
      ENDDO
      ALLOCATE(int1(nhm, nhm, 3, nat, nspin_mag), STAT = ierr)
      IF (ierr /= 0) CALL errore('init', 'Error allocating int1', 1)
      ALLOCATE(int2(nhm, nhm, 3, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('init', 'Error allocating int2', 1)
      IF (noncolin) THEN
        ALLOCATE(int1_nc(nhm, nhm, 3, nat, nspin), STAT = ierr)
        IF (ierr /= 0) CALL errore('init', 'Error allocating int1_nc', 1)
        IF (lspinorb) THEN
          ALLOCATE(int2_so(nhm, nhm, 3, nat, nat, nspin), STAT = ierr)
          IF (ierr /= 0) CALL errore('init', 'Error allocating int2_so', 1)
        ENDIF
      ENDIF ! noncolin
    ENDIF ! okvan
    !
    ALLOCATE(becp1(nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('init', 'Error allocating becp1', 1)
    ALLOCATE(alphap(3, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('init', 'Error allocating alphap', 1)
    !
    DO ik = 1, nks
      CALL allocate_bec_type(nkb, nbnd, becp1(ik))
      DO ipol = 1, 3
        CALL allocate_bec_type(nkb, nbnd, alphap(ipol,ik))
      ENDDO
    ENDDO
  ENDIF
  !
  DO na = 1, nat
    !
    ! xq here is the first q of the star
    arg = (xq(1) * tau(1, na) + &
           xq(2) * tau(2, na) + &
           xq(3) * tau(3, na)) * tpi
    !
    eigqts(na) = CMPLX(COS(arg), - SIN(arg), KIND = DP)
    !
  END DO
  !
  ! ... a0) compute rhocore for each atomic-type if needed for nlcc
  !
  IF (nlcc_any) CALL set_drhoc(xq, drc)
  !
  ! ... b) the fourier components of the local potential at q+G
  !        From PHonon/PH/phq_init.f90
  !
  CALL init_vlocq ( xq )
  !
  IF (first_run) THEN
    ALLOCATE(igk_k_all(npwx, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('init', 'Error allocating igk_k_all', 1)
    ALLOCATE(ngk_all(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('init', 'Error allocating ngk_all', 1)
    !
#if defined(__MPI)
    !
    CALL poolgather_int(npwx, nkstot, nks, igk_k(:, 1:nks), igk_k_all)
    CALL poolgather_int1(nkstot, nks, ngk(1:nks), ngk_all)
    !CALL mp_barrier(inter_pool_comm)
    !
#else
    !
    igk_k_all = igk_k
    ngk_all = ngk
    !
#endif
    !
    ! OBSOLETE: now we read directly the igkq to get the proper ngxx
    ! read only a piece of the map to save time
    ! the proper allocation bound would be ngxx = max(max(igkq))
    ! where the max is taken over the ig and the ik
    ! Here I use a simpler estimate: take the sphere npwx + two
    ! extra shells. This may not work for strange shapes of the
    ! reciproc latt. In this case just set ngxx = ngm_g
    !
    ! ngxx = NINT(4./3.*3.14*(2+(3.0/4.0/3.14*DBLE(npwx))**(1./3.))**3.)
    !
    ! Note that the k+q point below does not correspond to the actual (true)
    ! k+q, but since we only need to take the max over k and k+q this
    ! does not matter
    !
    ngxx = 0
    DO ik = 1, nkstot
      !
      itmp = MAXVAL(igk_k_all(1:ngk_all(ik), ik))
      IF (itmp > ngxx) THEN
        ngxx = itmp
      ENDIF
      !
    ENDDO
    !
    IF (okvan) THEN
#if defined(__MPI)
      !
      itmp = ldim_block(ngm, nproc_image, me_image)
      ig_s = gind_block(1, ngm, nproc_image, me_image)
      ig_e = ig_s + itmp - 1
      !
#else
      !
      ig_s = 1
      ig_e = ngm
      !
#endif
    ENDIF
    !
  ENDIF
  !
  IF (.NOT. first_run) THEN
    CALL dvanqq2()
  ENDIF
  !
  !
  CALL stop_clock('init')
  !
  !----------------------------------------------------------------------------
  END SUBROUTINE init
  !----------------------------------------------------------------------------
