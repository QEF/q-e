!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=======================================================================
!
module cg_sub
   implicit none
contains
   subroutine runcg_uspp(nfi, tfirst, tlast, eigr, bec, irb, eigrb, rhor, rhog, rhos, &
                         rhoc, ei1, ei2, ei3, sfac, fion, ema0bg, becdr, lambdap, lambda, &
                         nlam, vpot, c0, c0_d, cm, phi, dbec, l_cprestart)

!! please see https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.79.1337 (ensemble DFT)
!! and        https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.64.1045 (conjugate gradient)

      use kinds, only: dp
      use control_flags, only: tpre, iverbosity, tfor, tprnfor

!---ensemble-DFT
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin,          &
     &                    atot, entropy, egrand
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                nbspx, nbsp, ispin

      use ensemble_dft, only: tens, ef, z0t, c0diag, &
                              becdiag, fmat0, e0, id_matrix_init
!---
      use smallbox_gvec, only: ngb
      use gvecw, only: ngw, g2kin
      use gvect, only: gstart
      use ions_base, only: na, nat, nax, nsp, rcmax, ityp
      use cell_base, only: omega, alat, tpiba2
      use local_pseudo, only: vps, rhops
      use io_global, ONLY: stdout, ionode, ionode_id
      use mp_global, ONLY: intra_bgrp_comm
      use dener
      use constants, only: pi, au_gpa
      USE io_files, ONLY: tmp_dir, prefix
      use uspp, only: nkb, nkbus, &
                      betae => vkb, rhovan => becsum, &
                      deeq, qq_nt, nlcc_any, ofsbeta
      use uspp_param, only: nh, upf
      use cg_module, only: ene_ok, maxiter, niter_cg_restart, &
                           conv_thr, passop, enever, itercg, c0old, pre_state
      use ions_positions, only: tau0
      use efield_module, only: tefield, evalue, ctable, qmat, detq, ipolp, &
                               berry_energy, ctabin, gqq, gqqm, df, pberryel, &
                               tefield2, evalue2, ctable2, qmat2, detq2, ipolp2, &
                               berry_energy2, ctabin2, gqq2, gqqm2, pberryel2
      use mp, only: mp_sum, mp_bcast
      use cp_electronic_mass, ONLY: emass_cutoff
      use orthogonalize_base, ONLY: calphi_bgrp
      use cp_interfaces, ONLY: rhoofr, dforce, compute_stress, vofrho, nlfl_bgrp, prefor
      use cp_interfaces, ONLY: nlsm2_bgrp, calbec, caldbec_bgrp, nlfq_bgrp, runcp_uspp
      USE cp_main_variables, ONLY: idesc, drhor, drhog
      USE mp_global, ONLY: me_image, my_image_id, nbgrp
      USE fft_base, ONLY: dffts, dfftp
      use wave_gauge, only: project_parallel_gauge_2
      use wave_base, only: print_norm_square_difference ! for debugging

!
      implicit none
!

      LOGICAL, INTENT(in) :: l_cprestart !if true prepares a CG->CP restart

      CHARACTER(LEN=80) :: uname
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      integer :: nfi, nlam
      logical :: tfirst, tlast
      complex(dp) :: eigr(ngw, nat)
      real(dp) :: bec(nkb, nbspx)
      real(dp) :: becdr(nkb, nbspx, 3)
      integer irb(3, nat)
      complex(dp) :: eigrb(ngb, nat)
      real(dp) :: rhor(dfftp%nnr, nspin)
      real(dp) :: vpot(dfftp%nnr, nspin)
      complex(dp) :: rhog(dfftp%ngm, nspin)
      real(dp) :: rhos(dffts%nnr, nspin)
      real(dp) :: rhoc(dfftp%nnr)
      complex(dp) :: ei1(-dfftp%nr1:dfftp%nr1, nat)
      complex(dp) :: ei2(-dfftp%nr2:dfftp%nr2, nat)
      complex(dp) :: ei3(-dfftp%nr3:dfftp%nr3, nat)
      complex(dp) :: sfac(dffts%ngm, nsp)
      real(dp) :: fion(3, nat)
      real(dp) :: ema0bg(ngw)
      real(dp) :: lambdap(nlam, nlam, nspin)
      real(dp) :: lambda(nlam, nlam, nspin)
      complex(dp) :: c0(ngw, nbspx)
      complex(dp) :: c0_d(:, :)
      complex(dp) :: cm(ngw, nbspx)
      complex(dp) :: phi(ngw, nbspx)
      real(dp) :: dbec(:,:,:,:)
!
      include 'laxlib.fh'
!
      integer :: i, j, ig, k, is, iss, ia, iv, jv, il, ii, jj, kk, ip, nrlx
      integer :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot, comm
      real(dp) :: enb, enbi, x
      complex(dp), allocatable :: c2(:)
      complex(dp), allocatable :: c3(:)
      real(dp) :: gamma, entmp, sta
      complex(dp), allocatable :: hpsi(:, :), hpsi0(:, :), gi(:, :), hi(:, :)
      complex(dp), allocatable :: hpsi_d(:, :), gi_d(:, :)
#if defined(__CUDA)
      attributes(device) :: hpsi_d, c0_d, gi_d, phi
#endif
      real(DP), allocatable::               s_minus1(:, :)!factors for inverting US S matrix
      real(DP), allocatable::               k_minus1(:, :)!factors for inverting US preconditioning matrix
      real(DP), allocatable :: lambda_dist(:, :) ! replicated copy of lambda
      real(dp) :: sca, dumm(1)
      logical  :: newscheme, firstiter, filee
      integer :: maxiter3
!
!
      real(kind=DP), allocatable :: bec0(:, :), becm(:, :), becdrdiag(:, :, :)
      real(kind=DP), allocatable :: ave_ene(:)!average kinetic energy for preconditioning
      real(kind=DP), allocatable :: fmat_(:, :)!average kinetic energy for preconditioning

      !logical :: pre_state!if .true. does preconditioning state by state

      real(DP) esse, essenew !factors in c.g.
      logical ltresh!flag for convergence on energy
      real(DP) passo!step to minimum
      real(DP) etotnew, etotold!energies
      real(DP) spasso!sign of small step
      logical restartcg!if .true. restart again the CG algorithm, performing a SD step
      integer numok!counter on converged iterations
      integer iter3
      real(DP) passof, passov !step to minimum: effective, estimated
      real(DP) ene0, ene1, dene0, enesti !energy terms for linear minimization along hi

      INTEGER :: i_max
      REAL(kind=DP) :: max_sca

      nrlx = MAXVAL(idesc(LAX_DESC_NRLX, :))

      allocate (bec0(nkb, nbspx), becm(nkb, nbspx), becdrdiag(nkb, nbspx, 3))
      allocate (ave_ene(nbspx))
      allocate (c2(ngw), c3(ngw))

      call start_clock('runcg_uspp')

      if (nbgrp > 1) &
         call errore(' runcg_uspp ', ' parallelization over bands not yet implemented ', 1)
#if defined(__CUDA)
      if (nkbus > 0 ) &
         call errore('runcg_uspp','USPP for GPU not present in this version', 1)
      if (tens) &
         call errore('runcg_uspp','Ensemble DFT for GPU not present in this version', 1)
#endif
      if (pre_state .and. nkbus > 0) &
         call errore(' runcg_uspp ', ' preconditioning with kinetic energy not implemented for ultrasoft pseudopotentials')
      newscheme = .false.
      firstiter = .true.

      !pre_state = .true.!normally is disabled

      maxiter3 = 250

      if (ionode) then
         uname = TRIM(tmp_dir)//trim(prefix)//'.' &
                 //trim(int_to_char(my_image_id))//'_'//trim(int_to_char(me_image))
         !open(37,file='convergence.dat',status='unknown')!for debug and tuning purposes
         open (37, file=uname, status='unknown')!for debug and tuning purposes
      endif
      if (tfirst .and. ionode) &
         write (stdout, *) 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'

!set tpa preconditioning
!eq. 5.16 of https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.64.1045

      call emass_precond_tpa(ema0bg, tpiba2, emass_cutoff)

      call prefor(eigr, betae)

      ltresh = .false.
      itercg = 1
      etotold = 1.d8
      restartcg = .true.
      passof = passop
      ene_ok = .false.

      !orthonormalize c0

#if defined (__CUDA)
!TODO: is this necessary?
!NOTE (R.B.) : at the moment the array c0_d is used as a buffer. Most of the
!operation are done on the cpu, in particular for the ensemble dft case
      c0 = c0_d
#endif

!$acc data copy(betae, bec, c0)
      call calbec(nbsp, betae, c0, bec)
      CALL gram_bgrp(betae, bec, nkb, c0, ngw)

      !$acc update self(c0,betae,bec)

      !calculates phi for pcdaga

!bec stays on cpu
!$acc update self(bec)
      CALL calphi_bgrp(c0, SIZE(c0, 1), bec, nkb, betae, phi, nbsp)
!$acc end data

      !calculates the factors for S and K inversion in US case
      if (nkbus > 0) then
         allocate (s_minus1(nkb, nkb))
         allocate (k_minus1(nkb, nkb))
         call set_x_minus1(betae, s_minus1, dumm, .false.)
         call set_x_minus1(betae, k_minus1, ema0bg, .true.)
      else
         allocate (s_minus1(1, 1))
         allocate (k_minus1(1, 1))
      endif

      !set index on number of converged iterations

      numok = 0

!initialize  z0t
      call id_matrix_init(idesc, nspin)

      allocate (hpsi(ngw, nbspx), hpsi0(ngw, nbspx), gi(ngw, nbspx), hi(ngw, nbspx))
      !loop on cg iterations
      do while (itercg .lt. maxiter .and. (.not. ltresh))

         ENERGY_CHECK: if (.not. ene_ok) then
!$acc data copy(betae, bec) copyin(c0)
            call calbec(nbsp, betae, c0, bec)
!$acc end data
            call compute_energy(c0, bec, ens_goto_diagonal_repr=.true.) !result energy in etot
            if (.not. tens) then
               etotnew = etot
            else
               etotnew = etot + entropy
            end if

         else

            etot = enever
            if (.not. tens) then
               etotnew = etot
            else
               etotnew = etot + entropy
            endif
            ene_ok = .false.

         end if ENERGY_CHECK
         if (ionode) write (37, *) itercg, etotnew, pberryel, pberryel2!for debug and tuning purposes

         if (abs(etotnew - etotold) .lt. conv_thr) then
            numok = numok + 1
         else
            numok = 0
         endif

         if (numok .ge. 4) then
            ltresh = .true.
         endif

         etotold = etotnew
         ene0 = etot
         if (tens .and. newscheme) then
            ene0 = ene0 + entropy
         endif

         !update d

         call newd(vpot, rhovan, fion, .true.)

         call prefor(eigr, betae)
         ! this puts the gradient inside the array hpsi
#if defined (__CUDA)
         c0_d = c0
         !vkb_d = betae ! runcp uses a global variable!!!
#endif
         call runcp_uspp(0, 0.d0, 0.d0, ema0bg, 0.d0, rhos, bec, &
                         c0, c0_d, hpsi, hpsi_d, .false., .false., .true.)
         if (pre_state) call ave_kin(c0, SIZE(c0, 1), nbsp, ave_ene)

!$acc data copy(c0,hpsi)
         call pcdaga2(c0, phi, hpsi)
!$acc end data
         hpsi0 = hpsi
         gi = hpsi
! becm = <betae | hpsi >
         !$acc data copy(betae,hpsi,becm,c0,bec,gi,bec0)
         call calbec(nbsp, betae, hpsi, becm)
!        (1+S) |hpsi>
         call xminus1(hpsi,betae,dumm,becm,s_minus1,.false., .false., tpiba2, ave_ene)
! becm = <betae| (1+S) |hpsi>
         call calbec(nbsp, betae, hpsi, becm)
! project (1+S) |hpsi>
         call pc2(c0, bec, hpsi, becm)

         !$acc update self(hpsi)
         call xminus1(gi,betae,ema0bg,becm,k_minus1,.true., pre_state, tpiba2, ave_ene )
         call calbec(nbsp, betae, gi, becm)
         call pc2(c0, bec, gi, becm)

         if (tens) call calcmt(nrlx, f, z0t, fmat0)
         call calbec(nbsp, betae, hpsi, bec0)
         !$acc end data

!  calculates gamma
         gamma = 0.d0

         if (.not. tens) then
            !$acc data copyin(gi,hpsi,gstart)
            !$acc parallel loop reduction(+:gamma)
            do i = 1, nbsp
               !$acc loop reduction(+:gamma)
               do ig = 1, ngw
                  gamma = gamma + 2.d0*DBLE(CONJG(gi(ig, i))*hpsi(ig, i))
               enddo
               if (gstart == 2) then
                  gamma = gamma - DBLE(CONJG(gi(1, i))*hpsi(1, i))
               endif
            enddo
            !$acc end data
            call mp_sum(gamma, intra_bgrp_comm)

            if (nkbus .gt. 0) then
               do i = 1, nbsp
                  do ia = 1, nat
                     is = ityp(ia)
                     IF (upf(is)%tvanp) THEN
                        do iv = 1, nh(is)
                           do jv = 1, nh(is)
                              inl = ofsbeta(ia) + iv
                              jnl = ofsbeta(ia) + jv
                              gamma = gamma + qq_nt(iv, jv, is)*becm(inl, i)*bec0(jnl, i)
                           end do
                        end do
                     END IF
                  end do
               enddo
            endif

         else

            do iss = 1, nspin
               nss = nupdwn(iss)
               istart = iupdwn(iss)
               me_rot = idesc(LAX_DESC_MYPE, iss)
               np_rot = idesc(LAX_DESC_NPC, iss)*idesc(LAX_DESC_NPR, iss)
               allocate (fmat_(nrlx, nudx))
               do ip = 1, np_rot
                  if (me_rot == (ip - 1)) then
                     fmat_ = fmat0(:, :, iss)
                  end if
                  nrl = ldim_cyclic(nss, np_rot, ip - 1)
                  CALL mp_bcast(fmat_, ip - 1, intra_bgrp_comm)
                  do i = 1, nss
                     jj = ip
                     do j = 1, nrl
                        do ig = 1, ngw
                           gamma = gamma + 2.d0*DBLE(CONJG(gi(ig, i + istart - 1))*hpsi(ig, jj + istart - 1))*fmat_(j, i)
                        enddo
                        if (gstart == 2) then
                           gamma = gamma - DBLE(CONJG(gi(1, i + istart - 1))*hpsi(1, jj + istart - 1))*fmat_(j, i)
                        endif
                        jj = jj + np_rot
                     enddo
                  enddo
               enddo
               deallocate (fmat_)
            enddo
            if (nkbus .gt. 0) then
               do iss = 1, nspin
                  nss = nupdwn(iss)
                  istart = iupdwn(iss)
                  me_rot = idesc(LAX_DESC_MYPE, iss)
                  np_rot = idesc(LAX_DESC_NPC, iss)*idesc(LAX_DESC_NPR, iss)
                  allocate (fmat_(nrlx, nudx))
                  do ip = 1, np_rot
                     if (me_rot == (ip - 1)) then
                        fmat_ = fmat0(:, :, iss)
                     end if
                     nrl = ldim_cyclic(nss, np_rot, ip - 1)
                     CALL mp_bcast(fmat_, ip - 1, intra_bgrp_comm)

                     do i = 1, nss
                        jj = ip
                        do j = 1, nrl
                           do ia = 1, nat
                              is = ityp(ia)
                              IF (upf(is)%tvanp) THEN
                                 do iv = 1, nh(is)
                                    do jv = 1, nh(is)
                                       inl = ofsbeta(ia) + iv
                                       jnl = ofsbeta(ia) + jv
                                       gamma = gamma + qq_nt(iv, jv, is)*becm(inl, i + istart - 1) &
                                               *bec0(jnl, jj + istart - 1)*fmat_(j, i)
                                    end do
                                 end do
                              END IF
                           enddo
                           jj = jj + np_rot
                        enddo
                     enddo
                  end do
                  deallocate (fmat_)
               enddo
            endif
            call mp_sum(gamma, intra_bgrp_comm)
         endif

         !case of first iteration

         if (itercg == 1 .or. (mod(itercg, niter_cg_restart) .eq. 1) .or. restartcg) then

            restartcg = .false.
            passof = passop
            hi = gi!hi is the search direction
            esse = gamma

         else

            !find direction hi for general case
            !calculates gamma for general case, not using Polak Ribiere

            essenew = gamma
            gamma = gamma/esse
            esse = essenew

            hi = gi + gamma*hi

         endif
!note that hi, is saved  on gi, because we need it before projection on conduction states

         !find minimum along direction hi:

         !project hi on conduction sub-space
         !$acc data copy(betae, hi, bec0,c0,bec)
         call calbec(nbsp, betae, hi, bec0)
         call pc2(c0, bec, hi, bec0)
         !$acc end data

         !do quadratic minimization
         !
         !calculate derivative with respect to  lambda along direction hi

         dene0 = 0.
         if (.not. tens) then
            !$acc data copyin(hi,hpsi0,gstart)
            !$acc parallel loop reduction(+:dene0)
            do i = 1, nbsp
               !$acc loop reduction(+:dene0)
               do ig = 1, ngw
                  dene0 = dene0 - 4.d0*DBLE(CONJG(hi(ig, i))*hpsi0(ig, i))
               enddo
               if (gstart == 2) then
                  dene0 = dene0 + 2.d0*DBLE(CONJG(hi(1, i))*hpsi0(1, i))
               endif
            end do
            !$acc end data
         else
            !in the ensemble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
            !     calculation of the kinetic energy x=xmin
            call calcmt(nrlx, f, z0t, fmat0)
            do iss = 1, nspin
               nss = nupdwn(iss)
               istart = iupdwn(iss)
               me_rot = idesc(LAX_DESC_MYPE, iss)
               np_rot = idesc(LAX_DESC_NPC, iss)*idesc(LAX_DESC_NPR, iss)
               allocate (fmat_(nrlx, nudx))
               do ip = 1, np_rot
                  if (me_rot == (ip - 1)) then
                     fmat_ = fmat0(:, :, iss)
                  end if
                  nrl = ldim_cyclic(nss, np_rot, ip - 1)
                  CALL mp_bcast(fmat_, ip - 1, intra_bgrp_comm)
                  do i = 1, nss
                     jj = ip
                     do j = 1, nrl
                        do ig = 1, ngw
                           dene0 = dene0 - 2.d0*DBLE(CONJG(hi(ig, i + istart - 1))*hpsi0(ig, jj + istart - 1))*fmat_(j, i)
                           dene0 = dene0 - 2.d0*DBLE(CONJG(hpsi0(ig, i + istart - 1))*hi(ig, jj + istart - 1))*fmat_(j, i)
                        enddo
                        if (gstart == 2) then
                           dene0 = dene0 + DBLE(CONJG(hi(1, i + istart - 1))*hpsi0(1, jj + istart - 1))*fmat_(j, i)
                           dene0 = dene0 + DBLE(CONJG(hpsi0(1, i + istart - 1))*hi(1, jj + istart - 1))*fmat_(j, i)
                        end if
                        jj = jj + np_rot
                     enddo
                  enddo
               end do
               deallocate (fmat_)
            enddo
         endif

         call mp_sum(dene0, intra_bgrp_comm)

         !if the derivative is positive, search along opposite direction
         if (dene0 .gt. 0.d0) then
            spasso = -1.D0
         else
            spasso = 1.d0
         endif

         !calculates wave-functions on a point on direction hi

         cm = c0 + spasso*passof*hi

         !orthonormalize
         !$acc data copy(betae, cm, becm)
         call calbec(nbsp, betae, cm, becm)
         CALL gram_bgrp(betae, becm, nkb, cm, ngw)
         !$acc end data

         !calculate energy
         call compute_energy(cm, becm, ens_goto_diagonal_repr=.false.) ! result in etot
         ene1 = etot
         if (tens .and. newscheme) then
            ene1 = ene1 + entropy
         endif

         !find the minimum

         call minparabola(ene0, spasso*dene0, ene1, passof, passo, enesti)

         if (iverbosity > 0) write (6, *) ene0, dene0, ene1, passo, gamma, esse

         !set new step

         passov = passof
         passof = 2.d0*passo

         !calculates wave-functions at minimum

         cm = c0 + spasso*passo*hi
         if (gstart == 2) then
            cm(1, :) = 0.5d0*(cm(1, :) + CONJG(cm(1, :)))
         endif

         !$acc data copy(betae, cm, becm)
         call calbec(nbsp, betae, cm, becm)
         CALL gram_bgrp(betae, becm, nkb, cm, ngw)
         !$acc end data

         !test on energy: check the energy has really diminished
         call compute_energy(cm, becm, ens_goto_diagonal_repr=.false.) ! result in etot

         enever = etot
         if (tens .and. newscheme) then
            enever = enever + entropy
         endif
         if (tens .and. newscheme) then
            if (ionode) write (37, '(a3,4f20.10)') 'CG1', ene0, ene1, enesti, enever
            if (ionode) write (37, '(a3,4f10.7)') 'CG2', spasso, passov, passo, (enever - ene0)/passo/dene0
         else
            if (ionode) write (37, '(a3,4f20.10)') 'CG1', ene0 + entropy, ene1 + entropy, enesti + entropy, enever + entropy
            if (ionode) write (37, '(a3,4f10.7)') 'CG2', spasso, passov, passo, (enever - ene0)/passo/dene0
         endif
         !check with  what supposed

         if (ionode) then
            if (iverbosity > 1) then
               write (stdout, *) 'cg_sub: estimate :', (enesti - enever)/(ene0 - enever)
               write (stdout, *) 'cg_sub: minmum   :', enever, passo, passov
            endif
         endif

         !if the energy has diminished with respect to  ene0 and ene1 , everything ok
         if (((enever .lt. ene0) .and. (enever .lt. ene1)) .or. (tefield .or. tefield2)) then
            c0(:, :) = cm(:, :)
            bec(:, :) = becm(:, :)
            ene_ok = .true.
         elseif ((enever .ge. ene1) .and. (enever .lt. ene0)) then
            if (ionode) then
               write (stdout, *) 'cg_sub: missed minimum, case 1, iteration', itercg
            endif
            c0 = c0 + spasso*passov*hi
            restartcg = .true.
            !$acc data copy(betae, c0, bec)
            call calbec(nbsp, betae, c0, bec)
            CALL gram_bgrp(betae, bec, nkb, c0, ngw)
            !$acc end data
            ene_ok = .false.
            !if  ene1 << energy <  ene0; go to  ene1
         else if ((enever .ge. ene0) .and. (ene0 .gt. ene1)) then
            if (ionode) then
               write (stdout, *) 'cg_sub: missed minimum, case 2, iteration', itercg
            endif
            c0 = c0 + spasso*passov*hi
            restartcg = .true.!ATTENZIONE
            !$acc data copy(betae, c0, bec)
            call calbec(nbsp, betae, c0, bec)
            CALL gram_bgrp(betae, bec, nkb, c0, ngw)
            !$acc end data
            !if ene > ene0,en1 do a steepest descent step
            ene_ok = .false.
         else if ((enever .ge. ene0) .and. (ene0 .le. ene1)) then
            if (ionode) then
               write (stdout, *) 'cg_sub: missed minimum, case 3, iteration', itercg
            endif

            iter3 = 0
            do while (enever .gt. ene0 .and. iter3 .lt. maxiter3)
               iter3 = iter3 + 1
               passov = passov*0.5d0
               cm = c0 + spasso*passov*hi
               ! chenge the searching direction
               spasso = spasso*(-1.d0)
               !$acc data copy(betae, cm, becm)
               call calbec(nbsp, betae, cm, becm)
               CALL gram_bgrp(betae, bec, nkb, cm, ngw)
               call calbec(nbsp, betae, cm, becm)
               !$acc end data
               call compute_energy(cm, becm, ens_goto_diagonal_repr=.false.)

               enever = etot
               if (tens .and. newscheme) then
                  enever = enever + entropy
               endif

            end do
            if (iter3 == maxiter3) write (stdout, *) 'missed minimun: iter3 = maxiter3'
            c0(:, :) = cm(:, :)
            restartcg = .true.
            ene_ok = .false.
         end if

         if (tens .and. newscheme) enever = enever - entropy

         if (.not. ene_ok) then
            !$acc data copy(bec,betae) copyin(c0)
            call calbec(nbsp, betae, c0, bec)
            !$acc end data
         end if

!$acc data copyin(betae)
#if defined (__CUDA)
         c0_d = c0
         CALL calphi_bgrp(c0_d, SIZE(c0, 1), bec, nkb, betae, phi, nbsp)
         c0 = c0_d
#else
         !calculates phi for pc_daga
         CALL calphi_bgrp(c0, SIZE(c0, 1), bec, nkb, betae, phi, nbsp)
#endif
!$acc end data
         !=======================================================================
         !
         !                 start of the inner loop
         !                 (Uij degrees of freedom)
         !
         !=======================================================================
         if (tens .and. .not. newscheme) then
            call inner_loop_cold(nfi, tfirst, tlast, eigr, irb, eigrb, &
                                 rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, c0, bec, dbec, firstiter, vpot)
!the following sets up the new energy
            enever = etot
         endif

         !=======================================================================
         !                 end of the inner loop
         !=======================================================================

         itercg = itercg + 1

!   restore hi
!        hi(:,:)=gi(:,:)

      end do !on conjugate gradiient iterations

      !================================================================
      !calculates atomic forces and lambda

      if (tpre) then!if pressure is need the following is written because of caldbec
         if (.not. tens) then
            call caldbec_bgrp(eigr, c0, dbec, idesc)
            call rhoofr(nfi, c0(:, :), irb, eigrb, bec, dbec, rhovan, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6)
         else

            !     calculation of the rotated quantities
            call rotate(nrlx, z0t, c0(:, :), bec, c0diag, becdiag)
            !     calculation of rho corresponding to the rotated wavefunctions
            call caldbec_bgrp(eigr, c0diag, dbec, idesc)
            call rhoofr(nfi, c0diag, irb, eigrb, becdiag, dbec, rhovan, rhor, &
                        drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6)
         endif

         !calculates the potential
         !
         !     put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(rhoc)

         !
         !---ensemble-DFT

         vpot = rhor

         call vofrho(nfi, vpot, drhor, rhog, drhog, rhos, rhoc, tfirst, tlast,             &
                &        ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)

      endif

      call calcmt(nrlx, f, z0t, fmat0)


      call newd(vpot, rhovan, fion, .true.)
!$acc data copy(betae)
#if defined (__CUDA)
      if (.not. tens) then
         c0_d = c0
         if (tfor .or. tprnfor) call nlfq_bgrp(c0_d, betae, bec, becdr, fion)
      else
         c0_d = c0diag
         if (tfor .or. tprnfor) call nlfq_bgrp(c0_d, betae, becdiag, becdrdiag, fion)
      endif
#else
      if (.not. tens) then
         if (tfor .or. tprnfor) call nlfq_bgrp(c0, betae, bec, becdr, fion)
      else
         if (tfor .or. tprnfor) call nlfq_bgrp(c0diag, betae, becdiag, becdrdiag, fion)
      endif
#endif
!$acc end data
      call prefor(eigr, betae)
      ! this puts the gradient inside the array gi
#if defined (__CUDA)
      c0_d = c0
      !vkb_d = betae ! runcp uses a global variable!!!
#endif
      call runcp_uspp(0, 0.d0, 0.d0, ema0bg, 0.d0, rhos, bec, &
                      c0, c0_d, gi, gi_d, .false., .false., .true.)

      !$acc data copyin(c0, gi)
      call compute_lambda(c0,gi,lambda,nupdwn, iupdwn, nudx, nspin, ngw, intra_bgrp_comm, gstart)

      if (l_cprestart .and. .not. tens .and. nspin == 1 .and. nkbus < 1) then

!if required project c0 on previous manifold of occupied states
!NOT IMPLEMENTED YET FOR ENSEMBLE DFT AND NSPIN==2
!NOT IMPLEMENTED FOR US PSEUDOPOTENTIALS
         cm(:, :) = c0(:, :)
         !$acc data copyin(cm,c0old)
         nss = nupdwn(1)
         call project_parallel_gauge_2(c0old, cm, c0, &
                                       nss, ngw, ngw, gstart)
         !$acc end data
         !$acc update host(c0)
         !$acc data copy(betae, bec)
         call calbec(nbsp, betae, c0, bec)
         CALL gram_bgrp(betae, bec, nkb, c0, ngw)
         call calbec(nbsp, betae, c0, bec)
         !$acc end data
         !$acc update host(c0)
#if defined (__CUDA)
         c0_d = c0
         !vkb_d = betae ! runcp uses a global variable!!!
#endif
         call runcp_uspp(0, 0.d0, 0.d0, ema0bg, 0.d0, rhos, bec, &
                         c0, c0_d, gi, gi_d, .false., .false., .true.)
         !$acc update device(gi)
         !$acc update device(c0)

         call compute_lambda(c0,gi,lambda,nupdwn, iupdwn, nudx, nspin, ngw, intra_bgrp_comm, gstart)

         cm(:, :) = c0(:, :)

         !$acc data copy(betae, becm) copyin(cm)
         call calbec(nbsp, betae, cm, becm)
         !$acc end data

      endif
      !$acc end data

#if defined (__CUDA)
!TODO: is this necessary?
!NOTE (R.B.) : at the moment the array c0_d is used as a buffer. Most of the
!operation are done on the cpu, in particular for the ensemble dft case
      c0_d = c0
#endif

      if (tens) then
         !
         ! in the ensemble case matrix labda must be multiplied with f

         ALLOCATE (lambda_dist(nlam, nlam))

         do iss = 1, nspin
            !
            nss = nupdwn(iss)
            !
            lambdap(:, :, iss) = 0.0d0
            !
            CALL cyc2blk_redist(nss, fmat0(:, :, iss), nrlx, SIZE(fmat0, 2), &
                                lambda_dist, nlam, nlam, idesc(:, iss))
            !
            ! Perform lambdap = lambda * fmat0
            !
            CALL sqr_mm_cannon('N', 'N', nss, 1.0d0, lambda(:, :, iss), nlam, lambda_dist, nlam, &
                               0.0d0, lambdap(:, :, iss), nlam, idesc(:, iss))
            !
            lambda_dist = lambda(:, :, iss)
            lambda(:, :, iss) = lambdap(:, :, iss)
            lambdap(:, :, iss) = lambda_dist
            !
         end do
         !
         DEALLOCATE (lambda_dist)
         !
         call nlsm2_bgrp(ngw, nkb, betae, c0, becdr, nbspx, nbsp)
         !
      endif
      !

      !
      CALL nlfl_bgrp(bec, becdr, lambda, idesc, fion)

      ! bforceion adds the force term due to electronic berry phase
      ! only in US-case

      if (tefield .and. (evalue .ne. 0.d0)) then
         call bforceion(fion, tfor .or. tprnfor, ipolp, qmat, bec, becdr, gqq, evalue)

      endif
      if (tefield2 .and. (evalue2 .ne. 0.d0)) then
         call bforceion(fion, tfor .or. tprnfor, ipolp2, qmat2, bec, becdr, gqq2, evalue2)
      endif
      deallocate (hpsi0, hpsi, gi, hi)
      deallocate (s_minus1, k_minus1)
      if (ionode) close (37)!for debug and tuning purposes
      call stop_clock('runcg_uspp')

      deallocate (bec0, becm, becdrdiag)
      deallocate (ave_ene)
      deallocate (c2, c3)

      return

   CONTAINS

      SUBROUTINE compute_energy(cx, becx, ens_goto_diagonal_repr)

         COMPLEX(dp), intent(inout) :: cx(:, :) !input wf, it can be resetted to diagonal rep for ens dft
         REAL(dp), intent(inout) :: becx(:, :) !input bec, it can be resetted to diagonal rep for ens dft
         logical, intent(in) :: ens_goto_diagonal_repr
         !note: every declaration of the parent subroutine is valid here
         call start_clock('compute_energy')
         if (.not. tens) then
#if defined(__CUDA)
            c0_d = cx
#endif
            call rhoofr(nfi, cx, c0_d(:, :), becx, dbec, rhovan, rhor, &
                        drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6)
         else

            if (newscheme .or. firstiter) then
               call inner_loop_cold(nfi, tfirst, tlast, eigr, irb, eigrb, &
                                    rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, cx, becx, dbec, firstiter, vpot)
               firstiter = .false.
            endif
            !     calculation of the rotated quantities

            call rotate(nrlx, z0t, cx(:, :), becx, c0diag, becdiag)
            !     calculation of rho corresponding to the rotated wavefunctions
#if defined(__CUDA)
            c0_d = c0diag
#endif
            call rhoofr(nfi, cx, c0_d, becdiag, dbec, rhovan, rhor, &
                        drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6)
            call rhoofr(nfi, c0diag, irb, eigrb, becdiag, dbec, &
                        rhovan, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6)
         endif

!when cycle is restarted go to diagonal representation

         if (ens_goto_diagonal_repr .and. mod(itercg, niter_cg_restart) == 1 .and. itercg >= 2) then

            call rotate(nrlx, z0t, cx(:, :), becx, c0diag, becdiag)
            cx(:, :) = c0diag(:, :)
            becx(:, :) = becdiag(:, :)
            call id_matrix_init(idesc, nspin)
         endif

         !calculates the potential
         !
         !     put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(rhoc)

         !
         !---ensemble-DFT

         vpot = rhor

         call vofrho(nfi, vpot, drhor, rhog, drhog, rhos, rhoc, tfirst, tlast,             &
                &        ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
         if (.not. tens) then
            etotnew = etot
         else
            etotnew = etot + entropy
         end if

         if (tefield) then!just in this case calculates elfield stuff at zeo field-->to be bettered

            call berry_energy(enb, enbi, becx, cx(:, :), fion)
            etot = etot + enb + enbi
         endif
         if (tefield2) then!just in this case calculates elfield stuff at zeo field-->to be bettered

            call berry_energy2(enb, enbi, becx, cx(:, :), fion)
            etot = etot + enb + enbi
         endif
         call stop_clock('compute_energy')

      END SUBROUTINE
  

   END SUBROUTINE runcg_uspp

   subroutine compute_lambda(c0_, gi_, lambda, nupdwn, iupdwn, nudx, nspin, ngw, intra_bgrp_comm, gstart)
      use kinds, only : dp
      use mp, only : mp_sum
      USE cp_main_variables, ONLY: idesc
      implicit none
      include 'laxlib.fh'
      complex(kind=dp), intent(in) :: c0_(:,:),gi_(:,:)
      !$acc declare present(c0_,gi_)
      real(kind=dp), intent(inout) :: lambda(:,:,:)
      real(DP), allocatable :: lambda_repl(:, :)
      integer, intent(in) :: nupdwn(:), iupdwn(:), nudx, ngw, nspin, intra_bgrp_comm,  gstart

      real(kind=dp) ::  entmp
      integer :: is, nss, istart, i,j, iv, jv,ii,jj,ig

      

      ALLOCATE (lambda_repl(nudx, nudx))
      !$acc data create(lambda_repl)
      do is = 1, nspin
         !
         nss = nupdwn(is)
         istart = iupdwn(is)

         !
         !$acc parallel loop private(i,j,jj,ii,entmp,ig) present(c0_,gi_,lambda_repl)
         do jv = 0, nss*(nss + 1)/2 - 1
            i = jv/nss 
            j = mod(jv, nss)
            if (i > j) then
               j = nss - j - 1
               i = nss - i
            end if
            j = j + 1 ! equivalent to the indexes of a double triangular loop
            i = i + 1 ! do i = 1, nss do j = i, nss .
            !Now every cycle of the outer loop (the only one) has the same amount of work
            
            ii = i + istart - 1
            jj = j + istart - 1
            entmp = 0.0_dp
            !$acc loop vector reduction(+:entmp)
            do ig = 1, ngw
               entmp = entmp - 2.d0*DBLE(CONJG(c0_(ig, ii))*gi_(ig, jj))
            enddo
            if (gstart == 2) then
               entmp = entmp + DBLE(CONJG(c0_(1, ii))*gi_(1, jj))
            endif
            lambda_repl(j, i) = entmp
            lambda_repl(i, j) = entmp
         enddo
         !$acc end parallel
         !$acc host_data use_device(lambda_repl)
         CALL mp_sum(lambda_repl, intra_bgrp_comm)
         !$acc end host_data
         !$acc update host(lambda_repl)
         !
         CALL distribute_lambda(lambda_repl, lambda(:, :, is), idesc(:, is))
         !
      enddo
      !$acc end data
      deallocate(lambda_repl)         
      end subroutine

   !-----------------------------------------------------------------------
   subroutine calcmt(nrlx, fdiag, zmat, fmat)
      !-----------------------------------------------------------------------
      !! Constructs \(\text{fmat}=\text{z0}^t\cdot\text{fdiag}\cdot\text{z0},
      !!               \qquad \text{zmat}=\text{z0}^t\)
      !
      ! Constructs fmat=z0^t.fdiag.z0    zmat = z0^t
      !
      USE kinds, ONLY: DP
      use electrons_base, ONLY: nudx, nspin, nupdwn, iupdwn, nx => nbspx
      USE cp_main_variables, ONLY: idesc
      USE mp, ONLY: mp_sum, mp_bcast

      implicit none

      include 'laxlib.fh'

      integer, intent(in)  :: nrlx
      real(DP) :: zmat(nrlx, nudx, nspin), fmat(nrlx, nudx, nspin)
      !  NOTE: zmat and fmat are distributed by row across processors
      real(dp) :: fdiag(nx) !  fdiag is replicated

      integer  :: iss, nss, istart, i, j, k, ii, jj, kk
      integer  :: np_rot, me_rot, nrl, comm_rot, ip, nrl_ip

      real(DP), ALLOCATABLE :: mtmp(:, :)
      real(DP) :: f_z0t

      call start_clock('calcmt')

      fmat = 0.0d0

      DO iss = 1, nspin

         nss = nupdwn(iss)
         istart = iupdwn(iss)
         np_rot = idesc(LAX_DESC_NPR, iss)*idesc(LAX_DESC_NPC, iss)
         me_rot = idesc(LAX_DESC_MYPE, iss)
         nrl = idesc(LAX_DESC_NRL, iss)
         comm_rot = idesc(LAX_DESC_COMM, iss)

         IF (idesc(LAX_DESC_ACTIVE_NODE, iss) > 0) THEN

            ALLOCATE (mtmp(MAXVAL(idesc(LAX_DESC_NRLX, :)), nudx))

            DO ip = 1, np_rot

               IF (me_rot == (ip - 1)) THEN
                  mtmp = zmat(:, :, iss)
               END IF
               nrl_ip = ldim_cyclic(nss, np_rot, ip - 1)
               CALL mp_bcast(mtmp, ip - 1, comm_rot)

               DO j = 1, nss
                  ii = ip
                  DO i = 1, nrl_ip
                     f_z0t = fdiag(j + istart - 1)*mtmp(i, j)
                     DO k = 1, nrl
                        fmat(k, ii, iss) = fmat(k, ii, iss) + zmat(k, j, iss)*f_z0t
                     END DO
                     ii = ii + np_rot
                  END DO
               END DO

            END DO

            DEALLOCATE (mtmp)

         END IF

      END DO

      call stop_clock('calcmt')

      RETURN
   END SUBROUTINE calcmt

!-----------------------------------------------------------------------
   subroutine rotate(nrlx, z0, c0, bec, c0diag, becdiag)
!-----------------------------------------------------------------------
      use kinds, only: dp
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx, n => nbsp
      use uspp_param, only: nh, upf
      use uspp, only: nkb, qq_nt, ofsbeta
      use gvecw, only: ngw
      use ions_base, only: nat, ityp
      USE cp_main_variables, ONLY: idesc
      USE cp_interfaces, ONLY: protate

      implicit none
      include 'laxlib.fh'
      integer, intent(in) :: nrlx
      real(kind=DP) z0(nrlx, nudx, nspin)
      real(kind=DP) bec(nkb, n), becdiag(nkb, n)
      complex(kind=DP) c0(ngw, nx), c0diag(ngw, nx)
      integer :: np_rot, me_rot, nrl, comm_rot
      integer iss, nss, istart

      CALL start_clock('rotate')

      DO iss = 1, nspin
         istart = iupdwn(iss)
         nss = nupdwn(iss)
         np_rot = idesc(LAX_DESC_NPR, iss)*idesc(LAX_DESC_NPC, iss)
         me_rot = idesc(LAX_DESC_MYPE, iss)
         nrl = idesc(LAX_DESC_NRL, iss)
         comm_rot = idesc(LAX_DESC_COMM, iss)
         CALL protate(c0, bec, c0diag, becdiag, ngw, nss, istart, z0(:, :, iss), nrl, &
                      ityp, nat, ofsbeta, nh, np_rot, me_rot, comm_rot)
      END DO

      CALL stop_clock('rotate')
      return
   end subroutine rotate

   subroutine minparabola(ene0, dene0, ene1, passop, passo, stima)
!this subroutines finds the minimum of a quadratic real function

      use kinds, only: dp

      implicit none
      real(dp) ene0, dene0, ene1, passop, passo, stima
      real(dp) a, b, c!a*x^2+b*x+c

      c = ene0
      b = dene0
      a = (ene1 - b*passop - c)/(passop**2.d0)

      passo = -b/(2.d0*a)
      if (a .lt. 0.d0) then
         if (ene1 .lt. ene0) then
            passo = passop
         else
            passo = 0.5d0*passop
         endif
      endif

      stima = a*passo**2.d0 + b*passo + c

      return
   end subroutine minparabola

   subroutine pc2(a, beca, b, becb)
      !
      !! This subroutine applies the Pc operator:  
      !! a input: unperturbed wavefunctions;  
      !! b input: first order wavefunctions;  
      !! b output: \(b_i =|b_i-a_j\rangle\langle a_j|S|b_i\rangle\).
      !
      use kinds, only: dp
      use ions_base, only: na, nsp, nat, ityp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin, nupdwn, iupdwn, nspin
      use uspp_param, only: nh, upf
      use uspp, only: nkb, nkbus, ofsbeta
      use uspp, only: qq_nt

      implicit none

      complex(kind=DP) a(ngw, n), b(ngw, n)

      real(kind=DP) beca(nkb, n), becb(nkb, n)
      !$acc declare present(a,b,beca,becb)

      ! ... local variables

      integer is, iv, jv, ia, inl, jnl, i, j, ig
      real(kind=DP) sca
      real(DP), allocatable :: bectmp(:, :), qq_tmp(:, :), qqb_tmp(:, :)
      complex(DP), allocatable :: zbectmp(:, :)
      integer :: nl_max
      integer :: nss, iss, istart

      !logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock('pc2')

      do iss = 1, nspin
         nss = nupdwn(iss)
         istart = iupdwn(iss)
         allocate (bectmp(nss, nss))
         bectmp(:, :) = 0.d0
         allocate (zbectmp(nss, nss))
         !$acc data copyin(gstart,zbectmp,bectmp)
         !$acc host_data use_device(a,b,zbectmp)
         call myzgemm('C', 'N', nss, nss, ngw, (1.d0, 0.d0), a(:, istart), ngw, b(:, istart), ngw, (0.d0, 0.d0), zbectmp, nss)
         !$acc end host_data
         !$acc kernels loop
         do j = 1, nss
            do i = 1, nss
               bectmp(i, j) = 2.d0*dble(zbectmp(i, j))
               if (gstart == 2) bectmp(i, j) = bectmp(i, j) - DBLE(a(1, j))*DBLE(b(1, i))

            enddo
         enddo
         !$acc update self(bectmp)
         call mp_sum(bectmp(:, :), intra_bgrp_comm)
         !$acc update host(bectmp)
         if (nkbus >= 0) then

            nl_max = 0
            do is = 1, nsp
               nl_max = nl_max + nh(is)*na(is)
            enddo
            allocate (qq_tmp(nl_max, nl_max))
            allocate (qqb_tmp(nl_max, nss))
            qq_tmp(:, :) = 0.d0
            do ia = 1, nat
               is = ityp(ia)
               IF (upf(is)%tvanp) THEN
                  do iv = 1, nh(is)
                     do jv = 1, nh(is)
                        inl = ofsbeta(ia) + iv
                        jnl = ofsbeta(ia) + jv
                        qq_tmp(inl, jnl) = qq_nt(iv, jv, is)
                     enddo
                  enddo
               ENDIF
            enddo
            !$acc data copyin(qq_tmp,qqb_tmp)
            !$acc host_data use_device(qq_tmp,qqb_tmp,becb,beca,bectmp)
            call mydgemm('N', 'N', nl_max, nss, nl_max, 1.d0, qq_tmp, nl_max, becb(:, istart), nkb, 0.d0, qqb_tmp, nl_max)
            call mydgemm('T', 'N', nss, nss, nl_max, 1.d0, beca(:, istart), nkb, qqb_tmp, nl_max, 1.d0, bectmp, nss)
            !$acc end host_data
            !$acc end data
            deallocate (qq_tmp, qqb_tmp)
         endif
         !$acc kernels loop
         do i = 1, nss
            do j = 1, nss
               zbectmp(i, j) = CMPLX(bectmp(i, j), 0.d0, kind=dp)
            enddo
         enddo
         !$acc host_data use_device(a,b,zbectmp,bectmp,beca,becb)
         call myzgemm('N', 'N', ngw, nss, nss, (-1.d0, 0.d0), a(:, istart), ngw, zbectmp, nss, (1.d0, 0.d0), b(:, istart), ngw)
         call mydgemm('N', 'N', nkb, nss, nss, 1.0d0, beca(:, istart), nkb, bectmp, nss, 1.0d0, becb(:, istart), nkb)
         !$acc end host_data
         !$acc end data
         deallocate (zbectmp)
         deallocate (bectmp)
      enddo!on spin
      CALL stop_clock('pc2')
      return
   end subroutine pc2

   subroutine pcdaga2(a, as, b)
      !
      !! This subroutine applies the Pc operator:  
      !! a input: unperturbed wavefunctions;  
      !! b input: first order wavefunctions;  
      !! b output: \(b_i =b_i-S|a_j\rangle\langle a_j|b_i\rangle\).
      !
      use kinds
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw, n), b(ngw, n), as(ngw, n)
      !$acc declare present(a,b)
#if defined(__CUDA)
      attributes(device) :: as
#endif
      ! local variables
      integer is, iv, jv, ia, inl, jnl, i, j, ig
      real(dp) sca
      real(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2')
      allocate (scar(n))
      !$acc data copyin(ispin, gstart) create(scar)
      do j = 1, n
         !$acc kernels loop private(sca)
         do i = 1, n
            sca = 0.0d0
            if (ispin(i) == ispin(j)) then
               if (gstart == 2) b(1, i) = CMPLX(dble(b(1, i)), 0.0d0, kind=dp)
               !$acc loop vector reduction(+:sca)
               do ig = 1, ngw           !loop on g vectors
                  sca = sca + DBLE(CONJG(a(ig, j))*b(ig, i))
               enddo
               sca = sca*2.0d0  !2. for real weavefunctions
               if (gstart == 2) sca = sca - dble(a(1, j))*dble(b(1, i))
            endif
            scar(i) = sca
         enddo

         !$acc update self(scar)
         call mp_sum(scar, intra_bgrp_comm)
         !$acc update host(scar)

         !$acc kernels loop private(sca)
         do i = 1, n
            if (ispin(i) == ispin(j)) then
               sca = scar(i)
               !$acc loop vector
               do ig = 1, ngw
                  b(ig, i) = b(ig, i) - sca*as(ig, j)
               enddo
               ! this to prevent numerical errors
               if (gstart == 2) b(1, i) = CMPLX(dble(b(1, i)), 0.0d0, kind=dp)
            endif
         enddo
      enddo
      !$acc end data
      deallocate (scar)
      call stop_clock('pcdaga2')
      return
   end subroutine pcdaga2

   subroutine set_x_minus1(betae, m_minus1, ema0bg, use_ema)
      !
      !! This function calculates the factors for the inverse of the US K matrix.
      !! It takes care of the preconditioning.
      !
      use kinds, only: DP
      use ions_base, only: nsp, nat, ityp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum, mp_bcast
      use electrons_base, only: n => nbsp, ispin
      use uspp_param, only: nh, upf
      use uspp, only: nkb, qq_nt, ofsbeta
      use io_global, ONLY: ionode, ionode_id

      implicit none

      complex(DP) :: betae(ngw, nkb)
      real(DP)    :: m_minus1(nkb, nkb)
      real(DP)    :: ema0bg(ngw)
      logical     :: use_ema

! local variables
      real(DP), allocatable :: q_matrix(:, :), b_matrix(:, :), c_matrix(:, :)
      integer is, iv, jv, ia, inl, jnl, i, j, k, ig, js, ja
      real(DP) sca
      integer info, lwork
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: work(:)

      call start_clock('set_x_minus1')
      allocate (ipiv(nkb))
      allocate (work(nkb))

      lwork = nkb

      allocate (q_matrix(nkb, nkb), c_matrix(nkb, nkb))
!construct q matrix
      q_matrix(:, :) = 0.d0

      do ia = 1, nat
         is = ityp(ia)
         IF (upf(is)%tvanp) THEN
            do iv = 1, nh(is)
               do jv = 1, nh(is)
                  inl = ofsbeta(ia) + iv
                  jnl = ofsbeta(ia) + jv
                  q_matrix(inl, jnl) = qq_nt(iv, jv, is)
               enddo
            enddo
         END IF
      enddo

!construct b matrix
! m_minus1 used to be b matrix
      m_minus1(:, :) = 0.d0
      do ia = 1, nat
         is = ityp(ia)
         IF (upf(is)%tvanp) THEN
            do iv = 1, nh(is)
               do ja = 1, nat
                  js = ityp(ja)
                  IF (upf(js)%tvanp) THEN
                     do jv = 1, nh(js)
                        inl = ofsbeta(ia) + iv
                        jnl = ofsbeta(ja) + jv
                        sca = 0.d0
                        if (use_ema) then
                           ! k_minus case
                           do ig = 1, ngw           !loop on g vectors
                              sca = sca + ema0bg(ig)*DBLE(CONJG(betae(ig, inl))*betae(ig, jnl))
                           enddo
                           sca = sca*2.0d0  !2. for real wavefunctions
                           if (gstart == 2) sca = sca - ema0bg(1)*DBLE(CONJG(betae(1, inl))*betae(1, jnl))
                        else
                           ! s_minus case
                           do ig = 1, ngw           !loop on g vectors
                              sca = sca + DBLE(CONJG(betae(ig, inl))*betae(ig, jnl))
                           enddo
                           sca = sca*2.0d0  !2. for real weavefunctions
                           if (gstart == 2) sca = sca - DBLE(CONJG(betae(1, inl))*betae(1, jnl))
                        endif
                        m_minus1(inl, jnl) = sca
                     enddo
                  END IF
               enddo
            enddo
         END IF
      enddo
      call mp_sum(m_minus1, intra_bgrp_comm)

!calculate -(1+QB)**(-1) * Q
      CALL dgemm('N', 'N', nkb, nkb, nkb, 1.0d0, q_matrix, nkb, m_minus1, nkb, 0.0d0, c_matrix, nkb)

      do i = 1, nkb
         c_matrix(i, i) = c_matrix(i, i) + 1.d0
      enddo

      if (ionode) then
         call dgetrf(nkb, nkb, c_matrix, nkb, ipiv, info)
         if (info .ne. 0) write (stdout, *) 'set_k_minus1 Problem with dgetrf :', info
         call dgetri(nkb, c_matrix, nkb, ipiv, work, lwork, info)
         if (info .ne. 0) write (stdout, *) 'set_k_minus1 Problem with dgetri :', info
      endif
      call mp_bcast(c_matrix, ionode_id, intra_bgrp_comm)

      CALL dgemm('N', 'N', nkb, nkb, nkb, -1.0d0, c_matrix, nkb, q_matrix, nkb, 0.0d0, m_minus1, nkb)

      deallocate (q_matrix, c_matrix)
      deallocate (ipiv, work)
      call stop_clock('set_x_minus1')
      return
   end subroutine set_x_minus1
!

   SUBROUTINE emass_precond_tpa(ema0bg, tpiba2, emaec)
      use kinds, ONLY: dp
      use gvecw, ONLY: g2kin, ngw
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: ema0bg(ngw)
      REAL(DP), INTENT(IN) ::  tpiba2, emaec
      INTEGER :: i

      real(DP) :: x

      call start_clock('emass_p_tpa')
      do i = 1, ngw

         x = 0.5d0*tpiba2*g2kin(i)/emaec
         ema0bg(i) = 1.d0/(1.d0 + (16.d0*x**4)/(27.d0 + 18.d0*x + 12.d0*x**2 + 8.d0*x**3))
      end do
      call stop_clock('emass_p_tpa')
      RETURN
   END SUBROUTINE emass_precond_tpa

      subroutine xminus1(c0,betae,ema0bg,beck,m_minus1,do_k, pre_state, tpiba2, ave_ene)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0> 
! endif
      use kinds, only: dp
      use ions_base, only: nsp, nat, ityp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use uspp_param, only: nh, upf
      use uspp, only :nkb, nkbus, qq_nt, ofsbeta
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw, g2kin
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use gvect, only: gstart
!
      implicit none
      complex(dp) c0(ngw,n), betae(ngw,nkb)
      real(dp)    beck(nkb,n), ema0bg(ngw), tpiba2 
      real(DP)    :: m_minus1(nkb,nkb), ave_ene(n)
      logical :: do_k, pre_state
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      integer is, iv, jv, ia, inl, jnl, i, j, js, ja,ig
      real(dp) becktmp,x


      logical :: mat_par=.true.!if true uses parallel routines      

      call start_clock('xminus1')

      if (nkbus.gt.0) then
!calculates beck
         if (do_k) then
            beck(:,:) = 0.d0
            do ia=1,nat
               is=ityp(ia)
               IF( upf(is)%tvanp ) THEN
                  do iv=1,nh(is)
                     inl = ofsbeta(ia) + iv
                     do i=1,n
                        becktmp = 0.0d0
                        do ig=1,ngw
                           becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                        enddo
                        becktmp = becktmp*2.0d0
                        if (gstart==2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i))
                        beck(inl,i) = beck(inl,i) + becktmp
                     enddo
                  enddo
               END IF
            enddo
            call mp_sum( beck, intra_bgrp_comm )
         endif
!
!
      allocate(phi(ngw,n))
      allocate(qtemp(nkb,n))
      phi(1:ngw,1:n) = 0.0d0
      qtemp(:,:) = 0.0d0
      if(.not.mat_par) then
         call dgemm( 'N', 'N', nkb, n, nkb, 1.0d0, m_minus1,nkb ,    &
                    beck, nkb, 0.0d0, qtemp,nkb )
      else
         call para_dgemm( 'N', 'N', nkb, n, nkb, 1.0d0, m_minus1,nkb ,    &
                    beck, nkb, 0.0d0, qtemp,nkb,intra_bgrp_comm )
      endif

!NB  nkb is the total number of US projectors
!    it works because the first pseudos are the vanderbilt's ones

         CALL dgemm( 'N', 'N', 2*ngw, n, nkb, 1.0d0, betae, 2*ngw,    &
                    qtemp, nkb, 0.0d0, phi, 2*ngw )
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
               end do
            end do
         else
            do j=1,n
               do i=1,ngw
                  c0(i,j)=(phi(i,j)+c0(i,j))
               end do
            end do
         endif
      deallocate(qtemp,phi)

      else 
         if (do_k) then
            if (pre_state) then
               !$acc data present(g2kin) copyin(ave_ene,tpiba2) copy(c0)
               !$acc parallel loop collapse(2) private(x)
               do i = 1, n
                  do ig = 1, ngw
                     ! eq. 5.16 of https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.64.1045
                     x = tpiba2*g2kin(ig)/ave_ene(i)
                     c0(ig, i) = c0(ig, i)* &
                                 (27.0_dp + 18.0_dp*x + 12.0_dp*x**2 + 8.0_dp*x**3)/ &
                                 (27.0_dp + 18.0_dp*x + 12.0_dp*x**2 + 8.0_dp*x**3 + 16.0_dp*x**4)
                  end do
               end do
               !$acc end data
            else
               !$acc parallel loop collapse(2) copyin(ema0bg) copy(c0)
               do j=1,n
                  do ig=1,ngw
                     c0(ig,j)=c0(ig,j)*ema0bg(ig)
                  end do
               end do
            endif
         endif
      endif
      call stop_clock('xminus1')
      return
     end subroutine xminus1



   subroutine ave_kin(c, ngwx, n, ene_ave)
      !
      !! This subroutine calculates the average kinetic energy of
      !! each state, to be used for preconditioning.
      !
      USE kinds, ONLY: DP
      USE constants, ONLY: pi, fpi
      USE gvecw, ONLY: ngw
      USE gvect, ONLY: gstart
      USE gvecw, ONLY: g2kin
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: intra_bgrp_comm
      USE cell_base, ONLY: tpiba2

      IMPLICIT NONE

      ! input

      INTEGER, INTENT(IN) :: ngwx, n
      COMPLEX(kind=DP), INTENT(IN) :: c(ngwx, n)
      REAL(kind=DP), INTENT(out) :: ene_ave(n)!average kinetic energy to be calculated
      !
      ! local

      INTEGER  :: ig, i
      real(kind=dp) :: tmp

      !
      !$acc parallel loop private(tmp) present(g2kin) copyin(c,gstart,ngw) copyout(ene_ave)
      DO i = 1, n
         tmp = 0.d0
         !$acc loop vector reduction(+:tmp)
         DO ig = gstart, ngw
            tmp = tmp + DBLE(CONJG(c(ig, i))*c(ig, i))*g2kin(ig)
         END DO
         ene_ave(i) = tmp
      END DO

      CALL mp_sum(ene_ave(1:n), intra_bgrp_comm)
      ene_ave(:) = ene_ave(:)*tpiba2

      RETURN
   END subroutine ave_kin

!
! ... some simple routines for parallel linear algebra (the matrices are
! ... always replicated on all the cpus)
!
! ... written by carlo sbraccia ( 2006 )
!
!----------------------------------------------------------------------------
   SUBROUTINE para_dgemm( transa, transb, m, n, k, &
                          alpha, a, lda, b, ldb, beta, c, ldc, comm )
      !----------------------------------------------------------------------------
      !! Trivial parallelization (splitting matrix B by columns) of \(\text{dgemm}\).
      !
      USE kinds, ONLY: DP
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=1), INTENT(IN)    :: transa, transb
      INTEGER, INTENT(IN)    :: m, n, k
      REAL(DP), INTENT(IN)    :: alpha, beta
      INTEGER, INTENT(IN)    :: lda, ldb, ldc
      REAL(DP), INTENT(INOUT) :: a(lda, *), b(ldb, *), c(ldc, *)
      INTEGER, INTENT(IN)    :: comm
      !
      ! ... quick return if possible
      !
      IF (m == 0 .OR. n == 0 .OR. &
          ((alpha == 0.0_DP .OR. k == 0) .AND. beta == 1.0_DP)) RETURN
      !
!write(*,*) 'DEBUG: para_dgemm'
      !
      CALL rep_matmul_drv(transa, transb, m, n, k, &
                          alpha, a, lda, b, ldb, beta, c, ldc, comm)
      RETURN
      !
   END SUBROUTINE para_dgemm

end module
