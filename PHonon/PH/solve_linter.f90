!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter (irr, imode0, dfpt_data)
  !-----------------------------------------------------------------------
  !! Driver routine for the solution of the linear system that
  !! defines the change of the wavefunction due to a lattice distorsion.
  !! It performs the following tasks:
  !! a) computes the bare potential term \(\Delta V | \psi \rangle\)
  !!    and an additional term in the case of US pseudopotentials.
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) compute also the bare
  !!    potential term Delta \(V_\text{hub} | \psi \rangle\);
  !! b) adds to it the screening term \(\Delta V_\text{SCF} | \psi \rangle\).
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) computes also the SCF part
  !!    of the response Hubbard potential;
  !! c) applies \(P_c^+\) (orthogonalization to valence states);
  !! d) calls \(\text{cgsolve_all}\) to solve the linear system;
  !! e) computes \(\Delta\rho\), \(\Delta V_\text{SCF}\) and symmetrizes
  !!    them;
  !! f) If \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
  !!    occupation matrices \(\text{dnsscf}\);
  !! g) --Introduced in February 2020-- If \(\text{noncolin}=\text{TRUE}\)
  !!    and \(\text{domag}=\text{TRUE}\), the linear system is solved twice
  !!    (\(\text{nsolv}=2\), the case \(\text{isolv}=2\) needs the time-reversed
  !!    wave functions). For the theoretical background, please refer to:
  !!    Phys. Rev. B 100, 045115 (2019).
  !! Step b, c, d are done inside sternheimer_kernel.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : ionode
  USE io_files,             ONLY : prefix, diropn
  USE wavefunctions,        ONLY : evc
  USE cell_base,            ONLY : at
  USE ions_base,            ONLY : nat
  USE klist,                ONLY : xk, ngk, igk_k
  USE fft_base,             ONLY : dfftp, dffts
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, vkb, deeq_nc
  USE uspp_param,           ONLY : nhm
  USE noncollin_module,     ONLY : noncolin, domag, nspin_mag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential
  USE buffers,              ONLY : save_buffer, get_buffer
  USE ldaU,                 ONLY : lda_plus_u
  USE control_ph,           ONLY : ext_recover, lmultipole
  USE el_phon,              ONLY : elph
  USE uspp,                 ONLY : nlcc_any
  USE units_ph,             ONLY : iudrho, lrdrho, iubar, lrbar, iudvscf, iuint3paw, &
                                   lint3paw, iudrhous, lrdrhous
  USE units_lr,             ONLY : iuwfc, lrwfc
  USE output,               ONLY : fildrho, fildvscf
  USE phus,                 ONLY : alphap, int1_nc
  USE modes,                ONLY : u
  USE recover_mod,          ONLY : read_rec
  USE dfile_autoname,       ONLY : dfile_name
  USE save_ph,              ONLY : tmp_dir_save
  USE lrus,                 ONLY : int3_paw, becp1
  USE eqv,                  ONLY : dvpsi
  USE qpoint,               ONLY : xq, nksq, ikks, ikqs
  USE qpoint_aux,           ONLY : ikmks, becpt, alphapt
  USE control_lr,           ONLY : convt, rec_code, rec_code_read, where_rec
  USE uspp_init,            ONLY : init_us_2
  USE lr_nc_mag,            ONLY : int1_nc_save, deeq_nc_save
  USE dfpt_type,            ONLY : dfpt_data_type
  USE dfpt_kernels,         ONLY : dfpt_kernel
  USE phus,                 ONLY : becsumort
  USE recover_mod,          ONLY : write_rec
  !
  IMPLICIT NONE
  !
  EXTERNAL :: stop_smoothly_ph
  !
  integer :: irr
  !! input: the irreducible representation
  integer :: imode0
  !! input: the position of the modes
  TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
  !! Output: Data that describes linear response quantities
  !
  ! ... local variables
  !
  real(DP) :: dr2
  ! dr2   : self-consistency error

  logical :: exst
  !! used to open the recover file

  integer :: iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             nrec,       & ! the record number for dvpsi and dpsi
             mode,       & ! mode index
             isolv,      & ! counter on linear systems
             nsolv,      & ! number of linear systems
             ikmk          ! index of mk
  integer  :: npe
  integer  :: npw, npwq
  integer  :: iq_dummy
  character(len=256) :: filename

  integer :: nnr
  !
  IF (rec_code_read > 20 ) RETURN

  call start_clock ('solve_linter')
!
!  This routine is task group aware
!
  nsolv=1
  IF (noncolin.AND.domag) nsolv=2
  !
  npe = dfpt_data%npert
  !
  nnr = dffts%nnr
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnr, 1:nspin_mag, 1:npe))
  !
  if (rec_code_read == 10.AND.ext_recover) then
     ! restart from Phonon calculation
     CALL read_rec(dr2, iter0, dfpt_data)
     IF (okpaw .AND. convt) THEN
        CALL PAW_dpotential(dfpt_data%dbecsum,rho%bec,int3_paw,npe)
     ENDIF
     rec_code=0
  else
    iter0 = 0
    convt =.FALSE.
    where_rec='no_recover'
    dr2 = 0.d0
  endif
  !
  IF (lmultipole) THEN
    !
    IF (irr == 1) THEN
      iter0 = 0
      convt =.FALSE.
      where_rec='no_recover'
    ELSE
      convt=.TRUE.
      CALL init_rho(npe, dfpt_data%drhos, dfpt_data%drhop, iq_dummy)
    ENDIF
    !
  ENDIF
  !
  IF (convt) GOTO 155
  !
  !
  ! In this case it has recovered after computing the contribution
  ! to the dynamical matrix. This is a new iteration that has to
  ! start from the beginning.
  !
  IF (iter0==-1000) iter0=0
  !
  ! Compute dV_bare * psi and write to buffer iubar
  !
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npw  = ngk(ikk)
     npwq = ngk(ikq)
     !
     IF (lsda) current_spin = isk(ikk)
     !
     ! compute beta functions for k-point ikq
     !
     CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb)
     !
     DO isolv = 1, nsolv
        IF (isolv == 1) THEN
           ikmk = ikks(ik)
        ELSE
           ikmk = ikmks(ik)
        ENDIF
        !
        ! read unperturbed wavefunctions psi(k) and psi(k+q)
        !
        IF (nksq > 1 .OR. nsolv == 2) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
           !$acc update device(evc)
        ENDIF
        !
        DO ipert = 1, npe
           mode = imode0 + ipert
           nrec = (isolv-1) * npe * nksq + (ipert - 1) * nksq + ik
           !
           IF (isolv==1) THEN
              !
              CALL dvqpsi_us(ik, u(1, mode), .FAlSE., becp1, alphap)
              !
           ELSE
              IF (okvan) THEN
                 deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,2)
                 !$acc update device(deeq_nc)
                 int1_nc(:,:,:,:,:) = int1_nc_save(:,:,:,:,:,2)
              ENDIF
              CALL dvqpsi_us(ik, u(1, mode), .FAlSE., becpt, alphapt)
              IF (okvan) THEN
                 deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,1)
                 !$acc update device(deeq_nc)
                 int1_nc(:,:,:,:,:) = int1_nc_save(:,:,:,:,:,1)
              ENDIF
           ENDIF
           !
           CALL save_buffer(dvpsi, lrbar, iubar, nrec)
           !
        ENDDO ! ipert
     ENDDO ! isolv
  ENDDO ! ik
  !
  ! Compute the change of core charge due to atomic displacement
  ! drhoc is computed only once, stored in drhoc, passed to dfpt_kernel.
  !
  ! If lmultipole is .TRUE., one is actually considering a exp(iqr) perturbation, not
  ! a phonon perturbation. So we do not compute the core charge perturbation.
  !
  IF (.NOT. lmultipole) THEN
     ALLOCATE(dfpt_data%drhoc(dfftp%nnr, npe))
     DO ipert = 1, npe
        CALL addcore(u(1, imode0+ipert), dfpt_data%drhoc(1, ipert))
     ENDDO
     !
     ! Set Pulay correction terms
     !
     IF (okvan) THEN
        ! Pulay correction to density due to augmentation charge
        ALLOCATE(dfpt_data%drhop_pulay(dfftp%nnr, nspin_mag, npe))
        DO ipert = 1, npe
           CALL get_buffer(dfpt_data%drhop_pulay(1, 1, ipert), lrdrhous, iudrhous, imode0 + ipert)
        ENDDO
     ENDIF
     !
     IF (okpaw) THEN
        ! Pulay correction to dbecsum due to augmentation charge
        ALLOCATE(dfpt_data%dbecsum_pulay((nhm * (nhm + 1))/2, nat, nspin_mag, npe))
        DO ipert = 1, npe
           dfpt_data%dbecsum_pulay(:,:,:,ipert) = becsumort(:,:,:,imode0+ipert)
        ENDDO
     ENDIF
  ENDIF
  !
  ! Set records for restart
  !
  rec_code = 10
  where_rec = 'solve_lint'
  !
  !    Solve DFPT fixed-point equation
  !
  CALL dfpt_kernel('PHONON', npe, iter0, lrbar, iubar, dr2, dfpt_data, irr, imode0, &
                   write_rec_callback = write_rec, stop_callback = stop_smoothly_ph)
  !
155 CONTINUE
  !
  !    A part of the dynamical matrix requires the integral of
  !    the self consistent change of the potential and the variation of
  !    the charge due to the displacement of the atoms.
  !    We compute it here.
  !
  IF (convt) THEN
     !
     IF (lmultipole) CALL write_epsilon(npe, dfpt_data%drhop)
     !
     CALL drhodvus (irr, imode0, dfpt_data%dvscfp, npe)
     IF (nlcc_any) CALL dynmat_nlcc (imode0, dfpt_data%drhop, npe)
     !
     ! Write charge density to file
     !
     IF (fildrho /= ' ') THEN
        IF (ionode) THEN
           INQUIRE(UNIT = iudrho, OPENED = exst)
           IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
           filename = dfile_name(xq, at, fildrho, TRIM(tmp_dir_save)//prefix, generate=.true., index_q=iq_dummy)
           CALL diropn (iudrho, filename, lrdrho, exst)
        ENDIF ! ionode
        !
        DO ipert = 1, npe
           CALL davcio_drho(dfpt_data%drhop(1,1,ipert), lrdrho, iudrho, imode0+ipert, +1)
        ENDDO
        CLOSE (UNIT = iudrho, STATUS='keep')
        !
     ENDIF ! fildrho
     !
     if (fildvscf.ne.' ') then
        do ipert = 1, npe
           call davcio_drho ( dfpt_data%dvscfp(1,1,ipert),  lrdrho, iudvscf, imode0 + ipert, +1 )
           IF (okpaw.AND.ionode) CALL davcio( int3_paw(:,:,:,:,ipert), lint3paw, &
                                              iuint3paw, imode0+ipert, + 1 )
        end do
        if (elph) call elphel (irr, npe, imode0, dfpt_data%dvscfs)
     ENDIF ! fildvscf
     !
     IF (lda_plus_u) CALL dnsq_store(npe, imode0)
     !
  ENDIF ! convt
  !
  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)

  call stop_clock ('solve_linter')

  RETURN

END SUBROUTINE solve_linter
