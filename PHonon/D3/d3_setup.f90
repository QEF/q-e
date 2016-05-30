!
! Copyright (C) 2001-2008 Quantm-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3_setup()
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  d3toten program:
  !  1) computes the total local potential (external+scf) on the smoot
  !     grid to be used in h_psi and similia
  !  2) computes dmuxc 3.1) with GC if needed
  !  3) for metals sets the occupated bands
  !  4) computes alpha_pv
  !  5.1) computes the variables needed to pass to the pattern representat
  !       of the small group of q
  !     u      the patterns
  !     t      the matrices of the small group of q on the pattern basis
  !     tmq    the matrix of the symmetry which sends q -> -q + G
  !     gi     the G associated to each symmetry operation
  !     gimq   the G of the q -> -q+G symmetry
  !     irgq   the small group indices
  !     nsymq  the order of the small group of q
  !     irotmq the index of the q->-q+G symmetry
  !     nirr   the number of irreducible representation
  !     npert  the dimension of each irreducible representation
  !     nmodes the number of modes
  !     minus_q true if there is a symmetry sending q -> -q+G
  !  5.2) computes the variables needed to pass to the pattern representat
  !       of the group of the crystal
  !     ug0     the patterns
  !     tg0     the matrices of the group on the pattern basis
  !     nsymg0  the order of the group of the crystal
  !     nirrg0  the number of irreducible representation
  !     npertg0 the dimension of each irreducible representation
  !  6) set the variables needed to deal with nlcc
  !  7) set the variables needed to distribute one loop between pools
  !  8) set the variables needed to calculate only selected q=0 modes
  !
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp, tau
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE io_files,      ONLY : tmp_dir
  USE kinds,         ONLY : DP
  USE pwcom
  USE fft_base,      ONLY : dfftp
  USE scf, only : rho, rho_core, v, vltot, vrs, kedtau
  USE symm_base,     ONLY : nrot, nsym, s, ftau, irt, invs, inverse_s, &
                            s_axis_to_cart, find_sym, copy_sym, s_axis_to_cart
  USE uspp_param,    ONLY : upf
  USE uspp,          ONLY : nlcc_any
  USE control_flags, ONLY : iverbosity, modenum
  USE constants,     ONLY : degspin
  USE qpoint,        ONLY : xq, ikks, ikqs, nksq
  USE phcom
  USE d3com,         ONLY : q0mode, wrmode, nsymg0, npertg0, nirrg0, &
                            npert_i, npert_f, q0mode_todo, allmodes, ug0, &
                            fild0rho
  USE mp_global,     ONLY : npool, my_pool_id, inter_pool_comm, intra_image_comm
  USE mp,            ONLY : mp_max, mp_min, mp_bcast
  USE funct,         ONLY : dmxc, dmxc_spin

  USE lr_symm_base, ONLY : nsymq, irotmq, irgq, gi, gimq, minus_q, rtau
  USE control_lr,   ONLY : alpha_pv, nbnd_occ, lgamma
  !
  IMPLICIT NONE
  !
  REAL (DP) :: rhotot, rhoup, rhodw, TARGET, small, fac, xmax, emin, &
       emax, wrk, xqck(3)
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy
  ! working array

  INTEGER :: ir, isym, jsym, iinv, irot, jrot, ik, &
       ibnd, ipol, mu, nu, imode0, irr, ipert, nt, ii, nu_i
  ! counters
  LOGICAL :: sym (48), magnetic_sym
  ! the symmetry operations
  REAL (DP) :: mdum(3)
  CHARACTER(LEN=256) :: tmp_dir_save
#ifdef __MPI
  INTEGER :: nlength_w, nlength (npool), nresto
#endif
  CALL start_clock ('d3_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smoot grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Computes the derivative of the xc potential
  !
  dmuxc (:,:,:) = 0.d0
  IF (lsda) THEN
     DO ir = 1, dfftp%nnr
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        CALL dmxc_spin (rhoup, rhodw, dmuxc (ir, 1, 1), &
             dmuxc (ir, 2, 1), dmuxc (ir, 1, 2), dmuxc (ir, 2, 2) )
     ENDDO
  ELSE
     DO ir = 1, dfftp%nnr
        rhotot = rho%of_r (ir, nspin) + rho_core (ir)
        IF (rhotot > 1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
        IF (rhotot < - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
     ENDDO
  ENDIF
  !
  ! 3) Computes the number of occupated bands for each k point
  !
  call setup_nbnd_occ()
  !
  ! 4) Computes alpha_pv
  !
  emin = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emin = MIN (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the minimum across pools

  CALL mp_min( emin, inter_pool_comm )
  emax = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emax = MAX (emax, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the maximum across pools

  CALL mp_max( emax, inter_pool_comm )
  alpha_pv = 2.d0 * (emax - emin)
  ! avoid zero value for alpha_pv
  alpha_pv = MAX (alpha_pv, 1.0d-2)
  !
  ! 5) set all the variables needed to use the pattern representation
  !
  ! 5.0) Computes the inverse of each matrix
  !
  ! TEMP TEMP TEMP TEMP: this should not be needed any longer
  !
  modenum = 0
  magnetic_sym = .false.
  CALL find_sym ( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
               magnetic_sym, mdum )
  sym(:)       =.false.
  sym(1:nsym)=.true.
  !
  ! Here we re-order all rotations in such a way that true sym.ops.
  ! are the first nsymq; rotations that are not sym.ops. follow
  !
  call smallg_q (xq, modenum, at, bg, nsym, s, ftau, sym, minus_q)
  nsymq  = copy_sym ( nsym, sym )
  !
  nsymg0 = nsym
  CALL inverse_s ( )
  CALL s_axis_to_cart ( )
  nsym = nsymq
  !
  !  the first nsymq matrices are symmetries of the small group of q
  !
  ! 5.1) Finds the variables needeed for the pattern representation
  !      of the small group of q
  !
  sym(1:nsymg0)=.true.
  CALL sgam_ph (at, bg, nsymg0, s, irt, tau, rtau, nat, sym)
  nmodes = 3 * nat
  ! if minus_q=.t. set_irr will search for
  ! Sq=-q+G symmetry. On output minus_q=.t.
  ! if such a symmetry has been found
  minus_q = (modenum .eq. 0)
  !
  ! BEWARE: In set_irr, smallgq is called
  !
  ! FIXME: workaround for filename mess - needed to find where
  !        the patterns are
  tmp_dir_save=tmp_dir
  if ( lgamma ) tmp_dir=TRIM(tmp_dir)//'_ph0/'
  ! FIXME END
  IF (modenum .ne. 0) THEN
     npertx=1
     CALL allocate_pert_d3()
     CALL set_irr_mode (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u,    &
          npert, nirr, gi, gimq, iverbosity, modenum)
  ELSE
      IF(ionode) CALL io_pattern ( nat, fildrho, nirr, npert, u, xqck, tmp_dir, -1 )
      call mp_bcast(u,     ionode_id, intra_image_comm)
      call mp_bcast(nirr,  ionode_id, intra_image_comm)
      call mp_bcast(npert, ionode_id, intra_image_comm)
      call mp_bcast(xqck,  ionode_id, intra_image_comm)
         IF(SUM(ABS(xqck(:)-xq(:))) > 1.d-4) CALL errore('d3_setup', 'Wrong drho for q', 1)
      npertx = 0
      DO irr = 1, nirr
         npertx = max (npertx, npert (irr) )
      ENDDO
      IF (.not.lgamma) THEN
         IF(ionode) call io_pattern ( nat, fild0rho, nirrg0, npertg0, ug0, xqck, tmp_dir, -1 )
         call mp_bcast(ug0,     ionode_id, intra_image_comm)
         call mp_bcast(nirrg0,  ionode_id, intra_image_comm)
         call mp_bcast(npertg0, ionode_id, intra_image_comm)
         call mp_bcast(xqck,    ionode_id, intra_image_comm)
         IF(SUM(ABS(xqck(:))) > 1.d-4) CALL errore('d3_setup', 'Wrong drho for Gamma', 2)
         DO irr = 1, nirrg0
            npertx = max (npertx, npertg0 (irr) )
         ENDDO
      ENDIF
      CALL allocate_pert_d3()
      CALL set_sym_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
           irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u,   &
           npert, nirr, gi, gimq, iverbosity)
  ENDIF

  IF ( lgamma ) THEN
     !
     nksq = nks
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik=1,nksq
        ikks(ik) = ik
        ikqs(ik) = ik
     ENDDO
     !
  ELSE
     !
     nksq = nks / 2
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik=1,nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  END IF

  !
  ! 5.2) Finds the variables needeed for the pattern representation
  !      of the small group of the crystal
  !
  IF (lgamma) THEN
     nirrg0 = nirr
  ELSE
     !
     ! Calculates the variables need for the pattern representation
     ! for the q=0 symmetries
     !
     CALL set_d3irr ( )
     !
  ENDIF
  !
  ! FIXME: workaround for filename mess - needed to find where
  !        the patterns are
  tmp_dir=tmp_dir_save
  ! FIXME END
  npertx = 0
  do irr = 1, nirr
      npertx = max (npertx, npert (irr) )
  enddo
  do irr = 1, nirrg0
      npertx = max (npertx, npertg0 (irr) )
  enddo
  !
  ! 6) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  IF (nlcc_any) ALLOCATE (drc( ngm, ntyp))
  !
  ! 7) Sets up variables needed to distribute one loop between pools
  !
  npert_i = 1
  npert_f = 3 * nat
#ifdef __MPI
  nlength_w = (3 * nat) / npool
  nresto = 3 * nat - nlength_w * npool
  DO ii = 1, npool
     IF (ii <= nresto) THEN
        nlength (ii) = nlength_w + 1
     ELSE
        nlength (ii) = nlength_w
     ENDIF
  ENDDO
  npert_i = 1
  DO ii = 1, my_pool_id
     npert_i = npert_i + nlength (ii)
  ENDDO

  npert_f = npert_i - 1 + nlength (my_pool_id+1)
#endif
  !
  ! 8) Sets up variables needed to calculate only selected
  !    modes at q=0 --the first index of the third order matrix--
  !
  IF (q0mode_todo (1) <= 0) THEN
     DO ii = 1, 3 * nat
        q0mode (ii) = .TRUE.
     ENDDO
  ELSE
     DO ii = 1, 3 * nat
        q0mode (ii) = .FALSE.
     ENDDO
     ii = 1
     DO WHILE (q0mode_todo (ii) > 0)
        q0mode (q0mode_todo (ii) ) = .TRUE.
        ii = ii + 1
     ENDDO
  ENDIF
  !
  ! if you want to compute all the modes; and lgamma=.true.
  ! the calculation can be simplyfied, in this case allmodes
  ! is set .true.
  !
  allmodes = lgamma .AND. (q0mode_todo (1) <= 0)
  !
  ! Sets up variables needed to write only selected
  ! modes at q=0 --the first index of the third order matrix--
  !
  DO ii = 1, 3 * nat
     wrk = 0.d0
     DO nu_i = 1, 3 * nat
        IF (q0mode (nu_i) ) THEN
           wrk = wrk + ug0 (ii, nu_i) * CONJG (ug0 (ii, nu_i) )
        ENDIF
     ENDDO
     wrmode (ii) = .FALSE.
     IF (wrk > 1.d-8) wrmode (ii) = .TRUE.
  ENDDO
  CALL stop_clock ('d3_setup')
  RETURN
END SUBROUTINE d3_setup
