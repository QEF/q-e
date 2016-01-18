!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_linter_d3 (irr, imode0, npe, isw_sl)
  !-----------------------------------------------------------------------
  !    This routine is a driver for the solution of the linear system whic
  !    defines the change of the wavefunction due to the perturbation.
  !    It reads from a file the charge variation due to perturbation
  !    and calculates variation of the wavefunctions.
  !
  ! 1) It writes on file the proiection on conduction band of the variation
  !    of the wavefunction with respect to the perturbation
  !
  !    Several cases are possible:
  !  isw_sl = 1  : calculates | Pc d/du(q) psi_k >    and writes on:  iudqwf
  !  isw_sl = 2  : calculates | Pc d/du(0) psi_k+q >  and writes on:  iud0qwf
  !  isw_sl = 3  : calculates | Pc d/du(0) psi_k >    and writes on:  iudwf
  !
  ! 2) It writes on a file the scalar product of the wavefunctions with the
  !    K-S Hamiltonian
  !  isw_sl = 1  : calculates <psi_k+q|dH/du(q)|psi_k > and writes on: iupdqvp
  !  isw_sl = 3  : calculates <psi_k  |dH/du(0)|psi_k > and writes on: iupd0vp
  !
  USE ions_base,  ONLY : nat
  USE cell_base,  ONLY : tpiba2
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunigk
  USE gvect,      ONLY : g
  USE fft_base,   ONLY : dfftp 
  USE ener,       ONLY : ef
  USE klist,      ONLY : xk, wk, degauss, ngauss
  USE wvfct,      ONLY : nbnd, npwx, npw, igk, g2kin, et
  USE kinds, only : DP
  USE uspp, ONLY : vkb
  USE wavefunctions_module,  ONLY : evc
  use qpoint,     ONLY : xq, igkq, npwq, nksq
  use phcom
  use d3com
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum

  implicit none
  integer :: irr, npe, imode0, isw_sl
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes
  ! input: a switch

  real (DP) :: thresh, wg1, wg2, wwg, deltae, theta, anorm, averlt, &
       eprec1, aux_avg (2), tcpu, xq_ (3)
  ! the convergence threshold
  ! weight for metals
  ! weight for metals
  ! weight for metals
  ! difference of energy
  ! the theta function
  ! the norm of the error
  ! average number of iterations
  ! cut-off for preconditioning
  ! auxiliary variable for avg. iter. coun

  real (DP), external :: w0gauss, wgauss, get_clock
  ! function computing the delta function
  ! function computing the theta function
  ! cpu time

  complex (DP) ::  ps (nbnd), dbecsum, psidvpsi
  ! the scalar products
  ! dummy variable
  ! auxiliary dpsi dV matrix element between k+q  and  k wavefunctions
  complex (DP), external ::  zdotc

  real (DP), allocatable :: h_diag (:,:)
  ! the diagonal part of the Hamiltonian
  complex (DP), allocatable :: drhoscf (:,:), dvloc (:,:),  &
       spsi (:), auxg (:), dpsiaux (:,:)
  ! the variation of the charge
  ! variation of local part of the potential
  ! the function spsi
  logical :: q0mode_f, conv_root, lmetq0
  ! if .true. it is useless to compute this
  ! true if linter is converged
  ! true if xq=(0,0,0) in a metal

  integer :: ipert, ibnd, jbnd, lter, ltaver, lintercall, ik, ikk, &
       ikq, ig, ir, nrec, ios, mode, iuaux
  ! counters
  !
  external ch_psi_all2, cg_psi
  !
  call start_clock ('solve_linter')
  allocate  (drhoscf( dfftp%nnr, npe))
  allocate  (dvloc( dfftp%nnr, npe))
  allocate  (spsi( npwx))
  allocate  (auxg( npwx))
  if (degauss /= 0.d0) allocate  (dpsiaux( npwx, nbnd))
  allocate  (h_diag( npwx, nbnd))
  ltaver = 0
  lintercall = 0
  lmetq0 = (degauss /= 0.d0) .and. (isw_sl >= 3)
  thresh = ethr_ph
  if (isw_sl == 1) then
     xq_ = xq
  else
     xq_ = 0.d0
  endif
  !
  ! calculates the variation of the local part of the K-S potential
  !
  do ipert = 1, npe
     mode = imode0 + ipert
     call dvscf (mode, dvloc (1, ipert), xq_)
  enddo
  drhoscf (:,:) = (0.d0, 0.d0)
  rewind (unit = iunigk)

  do ik = 1, nksq
     read (iunigk, err = 100, iostat = ios) npw, igk
100  call errore ('solve_linter_d3', 'reading igk', abs (ios) )
     if (lgamma) then
        ikk = ik
        ikq = ik
        npwq = npw
     else
        read (iunigk, err = 200, iostat = ios) npwq, igkq
200     call errore ('solve_linter_d3', 'reading igkq', abs (ios) )
        if (isw_sl == 1) then
           ikk = 2 * ik - 1
           ikq = 2 * ik
        elseif (isw_sl == 2) then
           ikk = 2 * ik
           ikq = 2 * ik
           npw = npwq
           do ig = 1, npwx
              igk (ig) = igkq (ig)
           enddo
        elseif (isw_sl == 3) then
           ikk = 2 * ik - 1
           ikq = 2 * ik - 1
           npwq = npw
           do ig = 1, npwx
              igkq (ig) = igk (ig)
           enddo
        endif
     endif
     call init_us_2 (npw , igk , xk (1, ikk), vkb0)
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb )
     !
     ! reads unperturbed wavefuctions psi(k) and psi(k+q)
     !
     call davcio (evc, lrwfc, iuwfc, ikk, - 1)
     if (.not.lgamma) call davcio (evq, lrwfc, iuwfc, ikq, - 1)
     !
     ! compute the kinetic energy
     !
     do ig = 1, npwq
        g2kin (ig) = ( (xk (1, ikq) + g (1, igkq (ig) ) ) **2 + &
                       (xk (2, ikq) + g (2, igkq (ig) ) ) **2 + &
                       (xk (3, ikq) + g (3, igkq (ig) ) ) **2) * tpiba2
     enddo
     !
     do ipert = 1, npe
        q0mode_f = (.not.q0mode (imode0 + ipert) ) .and. (.not.lgamma) &
             .and. (isw_sl /= 1)
        if (q0mode_f) then
           psidqvpsi(:,:) = (0.d0, 0.d0)
           dpsi(:,:) = (0.d0, 0.d0)
           lintercall = 1
           goto 120
        endif
        !
        ! calculates dvscf_q*psi_k in G_space, for all bands
        !
        mode = imode0 + ipert
        call dvdpsi (mode, xq_, dvloc (1, ipert), vkb0, vkb, evc, dvpsi)
        !
        ! calculates matrix element of dvscf between k+q and k  wavefunctions,
        ! that will be written on a file
        !
        if (degauss /= 0.d0) then
           dpsiaux(:,:) = (0.d0, 0.d0)
        end if
        do ibnd = 1, nbnd
           if (isw_sl /= 2) then
              do jbnd = 1, nbnd
                 psidvpsi = zdotc(npwq, evq (1, jbnd), 1, dvpsi (1, ibnd),1)
#ifdef __MPI
                 call mp_sum ( psidvpsi, intra_pool_comm )
#endif
                 psidqvpsi (jbnd, ibnd) = psidvpsi
                 if (degauss /= 0.d0) then
                    deltae = et (ibnd, ikk) - et (jbnd, ikq)
                    !            theta = 2.0d0*wgauss(deltae/degauss,0)
                    theta = 1.0d0
                    if (abs (deltae) > 1.0d-5) then
                       wg1 = wgauss ( (ef-et (ibnd, ikk) ) / degauss, ngauss)
                       wg2 = wgauss ( (ef-et (jbnd, ikq) ) / degauss, ngauss)
                       wwg = (wg1 - wg2) / deltae
                    else
                       wwg = - w0gauss ( (ef - et (ibnd, ikk) ) / degauss, &
                            ngauss) / degauss
                    endif
                    psidvpsi = 0.5d0 * wwg * psidvpsi * theta
                    call zaxpy(npwq,psidvpsi,evq(1,jbnd),1,dpsiaux(1,ibnd),1)
                 endif
              enddo
           endif
        enddo
        !
        ! Ortogonalize dvpsi
        !
        call start_clock ('ortho')
        wwg = 1.0d0
        do ibnd = 1, nbnd_occ (ikk)
           auxg (:) = (0.d0, 0.d0)
           do jbnd = 1, nbnd
              ps (jbnd) = - wwg * zdotc(npwq, evq(1,jbnd), 1, dvpsi(1,ibnd), 1)
           enddo
           call mp_sum ( ps, intra_pool_comm )
           do jbnd = 1, nbnd
              call zaxpy (npwq, ps (jbnd), evq (1, jbnd), 1, auxg, 1)
           enddo
           call zcopy (npwq, auxg, 1, spsi, 1)
           call daxpy (2 * npwq, 1.0d0, spsi, 1, dvpsi (1, ibnd), 1)
        enddo
        call stop_clock ('ortho')
        call dscal (2 * npwx * nbnd, - 1.d0, dvpsi, 1)
        !
        ! solution of the linear system (H-eS)*dpsi=dvpsi,
        ! dvpsi=-P_c^+ (dvscf)*psi
        !
        dpsi (:,:) = (0.d0, 0.d0)
        do ibnd = 1, nbnd_occ (ikk)
           conv_root = .true.
           do ig = 1, npwq
              auxg (ig) = g2kin (ig) * evq (ig, ibnd)
           enddo
           eprec1 = zdotc (npwq, evq (1, ibnd), 1, auxg, 1)
           call mp_sum ( eprec1, intra_pool_comm )
           do ig = 1, npwq
              h_diag (ig, ibnd) = 1.d0/ max (1.0d0, g2kin (ig) / eprec1)
           enddo
        enddo

        call cgsolve_all (ch_psi_all2, cg_psi, et (1, ikk), dvpsi, dpsi, &
             h_diag, npwx, npwq, thresh, ik, lter, conv_root, anorm, &
             nbnd_occ (ikk), 1 )

        ltaver = ltaver + lter
        lintercall = lintercall + 1
        if (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4, &
             & " linter: root not converged ",e10.3)') ikk, ibnd, anorm
120     continue
        !
        ! writes psidqvpsi on iupdqvp
        !
        nrec = imode0 + ipert + (ik - 1) * 3 * nat
        if (isw_sl == 1) then
           call davcio (psidqvpsi, lrpdqvp, iupdqvp, nrec, + 1)
        elseif (isw_sl >= 3) then
           call davcio (psidqvpsi, lrpdqvp, iupd0vp, nrec, + 1)
        endif
        !
        ! writes delta_psi on iunit iudwf, k=kpoint,
        !
        if (isw_sl == 1) then
           iuaux = iudqwf
        elseif (isw_sl >= 3) then
           iuaux = iudwf
        elseif (isw_sl == 2) then
           iuaux = iud0qwf
        endif
        nrec = (imode0 + ipert - 1) * nksq + ik

        call davcio (dpsi, lrdwf, iuaux, nrec, + 1)
        if (q0mode_f) goto 110
        if (isw_sl /= 2) then
           if (degauss /= 0.d0) then
              do ibnd = 1, nbnd
                 wg1 = wgauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss)
                 call dscal (2 * npwq, wg1, dpsi (1, ibnd), 1)
              enddo
              call daxpy (2 * npwx * nbnd, 1.0d0, dpsiaux, 1, dpsi, 1)
           endif
        endif
110     continue
        !
        ! This is used to calculate Fermi energy shift at q=0 in metals
        !
        if (lmetq0) call incdrhoscf2 (drhoscf (1, ipert), wk (ikk), &
             ik, dbecsum, 1, 1)
     enddo

  enddo
  if (lmetq0) then
     do ipert = 1, npe
        call cinterpolate (drhoscf (1, ipert), drhoscf (1, ipert), 1)
     enddo
  endif
#ifdef __MPI
  call mp_sum( drhoscf, inter_pool_comm )
#endif

  if (lmetq0) call set_efsh (drhoscf, imode0, irr, npe)
  aux_avg (1) = DBLE (ltaver)
  aux_avg (2) = DBLE (lintercall)
  call mp_sum( aux_avg, inter_pool_comm )

  averlt = aux_avg (1) / aux_avg (2)
  tcpu = get_clock ('D3TOTEN')

  WRITE( stdout, '(//,5x," thresh=",e10.3," total cpu time : ",f8.1, &
       &      " s   av.# it.: ",f5.1)') thresh, tcpu, averlt
  !
  FLUSH( stdout )
  !
  deallocate (h_diag)
  if (degauss /= 0.d0) deallocate (dpsiaux)
  deallocate (auxg)
  deallocate (spsi)
  deallocate (dvloc)
  deallocate (drhoscf)

  call stop_clock ('solve_linter')
  return
end subroutine solve_linter_d3
