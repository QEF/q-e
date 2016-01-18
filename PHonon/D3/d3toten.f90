!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program d3toten
  !-----------------------------------------------------------------------
  !
  use pwcom
  use qpoint,        ONLY : xq
  use phcom
  use d3com
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp, zv, tau
  USE io_global,     ONLY : stdout
  use io_files,      ONLY : prefix
  use control_flags, ONLY : gamma_only
  USE mp_global,     ONLY : mp_startup
  USE environment,   ONLY : environment_start

  implicit none
  character(len=9) :: cdate, ctime, code = 'D3TOTEN'
  integer :: nu_i, nu_i0, irecv
  real (DP) :: t0, t1, get_clock
  !
  !
  gamma_only = .false.
  all_done=.false.
  !
  ! Initialize MPI, clocks, print initial messages
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( code )
  !
  ! Initialization routines
  !
  call d3_readin
  call allocate_d3
  call d3_setup
  call d3_summary
  call openfild3
  call d3_init
  call print_clock ('D3TOTEN')
  !
  ! Used for testing purposes: if wraux=.true. it writes
  ! different terms of the third derivative matrix in different files.
  !
  if (wraux) call write_aux (1)

  d3dyn(:,:,:) = (0.d0, 0.d0)
  !
  nu_i0 = 1
  if (recv) then
     !
     ! If recv.eq.true. this is a recover run
     !
     call d3_recover (irecv, - 1)

     WRITE( stdout,  * ) ' Recover Run             index:', irecv
     if (irecv.ge.401.and.irecv.lt.499) then
        nu_i0 = irecv - 400
        goto 304
     else
        goto (301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, &
             312, 313) irecv
     endif
  endif
  !
  !  Non-selfconsistent calculation of the wavefunctions
  !
  write( stdout, '(/,5x,"Nscf calculating of the perturbed wavefunctions")')
  !
  ! It calculates the variation of wavefunctions | d/du(q) psi(k) >
  !
  t0 = get_clock ('D3TOTEN')
  if (.not.lgamma) then
!     WRITE( stdout, '(/,5x,"calling gen_dwfc(1)")')
     write( stdout, '(/,5x,"Calculating for the wavevector q")')

     call gen_dwfc (1)
     call d3_recover (1, + 1)
     t1 = get_clock ('D3TOTEN') - t0
     t0 = get_clock ('D3TOTEN')
     WRITE( stdout, '(5x,"gen_dwfc(1)   cpu time:",f9.2, &
          &         " sec    Total time:",f12.2," sec")') t1, t0
  endif
  if (istop.eq.1) stop
  !
  ! It calculates the variation of wavefunctions | d/du(q=0) psi(k) >
  !
301 continue
!  WRITE( stdout, '(/,5x,"calling gen_dwfc(3)")')
  write( stdout, '(/,5x,"Calculating for the wavevector q=0 at the original k-points")')

  call gen_dwfc (3)
  call d3_recover (2, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"gen_dwfc(3)   cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.2) stop
  !
  ! It calculates the variation of wavefunctions | d/du(q=0) psi(k+q) >
  !  to be used for the terms < dpsi | dpsi ><psi| dH |psi>
  !
302 continue
  if (.not.lgamma) then
     write( stdout, '(/,5x,"Calculating for the wavevector q=0 at the (k+q)-points")')
     WRITE( stdout, '(/,5x,"calling gen_dwfc(2)")')

     call gen_dwfc (2)
     call d3_recover (3, + 1)
     t1 = get_clock ('D3TOTEN') - t0
     t0 = get_clock ('D3TOTEN')
     WRITE( stdout, '(5x,"gen_dwfc(2)   cpu time:",f9.2, &
          &          " sec    Total time:",f12.2," sec")') t1, t0
  endif

  write( stdout, '(/,5x,"Finished the ncf calculation of the perturbed wavefunctions")')

  if (istop.eq.3) stop
  !
  ! It writes on files terms of the type: <dpsi| dH | psi>, that
  ! will be used for the metallic case
  !
303 continue
  WRITE( stdout, '(/,5x,"calling gen_dpdvp")')
  call gen_dpdvp
  call d3_recover (4, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"gen_dpdvp     cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.4) stop
  !
  ! It calculates the term < dpsi| dH | dpsi >
  !
304 continue

  WRITE( stdout, '(/,5x,"Calculating the matrix elements <dpsi |dH |dpsi>")')
  do nu_i = nu_i0, 3 * nat
     if (q0mode (nu_i) ) then
        WRITE( stdout, '(/,5x,"calling dpsidvdpsi:",i3)') nu_i
        call dpsidvdpsi (nu_i)
        call d3_recover (401 + nu_i, + 1)
        t1 = get_clock ('D3TOTEN') - t0
        t0 = get_clock ('D3TOTEN')

        WRITE( stdout, '(5x,"dpsidvdpsi",i3," cpu time:",f9.2, &
             &   " sec    Total time:",f12.2," sec")') nu_i, t1, t0

        if (istop.gt.400.and.nu_i.eq. (istop - 400) ) stop
     endif
  enddo
  call d3_recover (5, + 1)
  if (istop.eq.5) stop
  !
  ! It calculates the term < dpsi| dpsi > < psi | dH | psi>
  !
305 continue
  WRITE( stdout, '(/,5x,"Calculating the matrix elements <dpsi|dpsi>< psi|dH|psi> ")')
  WRITE( stdout, '(/,5x,"calling dpsidpsidv")')
  call dpsidpsidv
  call d3_recover (6, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"dpsidpsidv    cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.6) stop
  !
  ! It calculates the term   drho * d2V
  !
306 continue
  WRITE( stdout, '(/,5x,"Calculating the matrix elements <psi |d^2 v |dpsi>")')
  WRITE( stdout, '(/,5x,"calling drhod2v")')
  call drhod2v
  call d3_recover (7, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"drhod2v       cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.7) stop
  !
  ! It calculates the term   rho * d3V
  !
307 continue
  WRITE( stdout, '(/,5x,"Calculating the matrix elements <psi |d^3v |psi>")')
  WRITE( stdout, '(/,5x,"calling d3vrho")')
  call d3vrho
  call d3_recover (8, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"d3vrho        cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.8) stop
  !
  ! It calculates the contribution due to ionic term
  !
308 continue
  WRITE( stdout, '(/,5x,"Calculating the Ewald contribution")')
  WRITE( stdout, '(/,5x,"calling d3ionq")')
  call d3ionq (nat, ntyp, ityp, zv, tau, alat, omega, xq, at, bg, g, &
       gg, ngm, gcutm, nmodes, u, ug0, npert_i, npert_f, q0mode, d3dyn)
  call d3_recover (9, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"d3ionq        cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.9) stop
  !
  ! In the metallic case some additional terms are needed
  !
309 continue
  WRITE( stdout, '(/,5x,"Calculating the valence contribution")')
  WRITE( stdout, '(/,5x,"calling d3_valence")')
  call d3_valence
  call d3_recover (10, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"d3_valence    cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.10) stop
  !
  ! drho_cc(+1) adds to the variation or the charge -written on a file-
  ! the variation of the core charge. The variation of the charge,
  ! modified this way is used by the routines d3_exc and d3dyn_cc.
  ! drho_cc(-1) restores drho as it was before (useless)
  !
310 continue
  WRITE( stdout, '(/,5x,"calling drho_cc(+1)")')
  call drho_cc ( + 1)
  call d3_recover (11, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"drho_cc(+1)   cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  !
  ! It calculates d3Ei * drho * drho * drho, where drho is the variation
  ! of the charge and d3Ei is the third derivative of the
  ! Kohn-Sham-Energy term depending on the charge density.
  !
311 continue
  WRITE( stdout, '(/,5x,"Calculating the exchange-correlation contribution")')
  WRITE( stdout, '(/,5x,"calling d3_exc")')
  call d3_exc
  call d3_recover (12, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  WRITE( stdout, '(5x,"d3_exc        cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  !
  ! It calculates additional terms due to non_linear-core-corrections
  !
312 continue
  WRITE( stdout, '(/,5x,"Calculating the core-correction contribution")')
  WRITE( stdout, '(/,5x,"calling d3dyn_cc")')
  call d3dyn_cc
  call d3_recover (13, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')

  WRITE( stdout, '(5x,"d3dyn_cc      cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  !
  ! drho is restored as it was before
  !
  !      WRITE( stdout,'(/,5x,"calling drho_cc(-1)")')
  !      call drho_cc(-1)
  !      t1 = get_clock('D3TOTEN') - t0
  !      t0 = get_clock('D3TOTEN')
  !      WRITE( stdout,'(5x,"drho_cc(-1)   time: ",f12.2,
  !     +       " sec    Total time:",f12.2," sec")') t1,t0
  if (wraux) call write_aux (2)
  !
  ! Symmetrizes d3dyn, calculates the q in the star and writes the result
  ! for every q on a file.
  !
313 continue
  WRITE( stdout, '(/,5x,"Symmetrizing and writing the tensor to disc")')
  WRITE( stdout, '(/,5x,"calling d3matrix")')
  call d3matrix
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')

  WRITE( stdout, '(5x,"d3matrix      cpu time:",f9.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (wraux) call write_aux (3)

  call stop_d3 (.true.)
end program d3toten
