!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine solve_linter (irr, imode0, npe, drhoscf)
  !-----------------------------------------------------------------------
  !
  !    Driver routine for the solution of the linear system which
  !    defines the change of the wavefunction due to the perturbation.
  !    It performs the following tasks:
  !     a) It computes the kinetic energy
  !     b) It adds the term Delta V_{SCF} | psi > and the additional one
  !        in the case of US pseudopotentials
  !     c) It applies P_c^+ to the known part
  !     d) It calls linter to solve the linear system
  !     e) It computes Delta rho, Delta V_{SCF} and symmetrize them
  !
#include "machine.h"

  use pwcom
  use becmod
  use parameters, only : DP
  use phcom
  implicit none


  integer :: irr, npe, imode0
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  complex(kind=DP) :: drhoscf (nrxx, nspin, npe)
  ! output: the change of the scf charge

  real(kind=DP) , allocatable :: h_diag (:,:),eprec (:)
  real(kind=DP) :: dos_ef,  thresh, wg1, w0g, wgp, &
       wwg, weight, deltae, theta, anorm, averlt, aux_avg (2), &
       dr2, w0gauss, wgauss
  ! density of states at Ef
  ! the diagonal part of the Hamiltonia
  ! the convergence threshold
  ! weight for metals
  ! weight for metals
  ! weight for metals
  ! weight for metals
  ! used for summation over k points
  ! difference of energy
  ! the theta function
  ! the norm of the error
  ! average number of iterations
  ! cut-off for preconditioning
  ! auxiliary variable for avg. iter. count
  ! convergence limit
  ! function computing the delta function
  ! function computing the theta function

  complex(kind=DP), pointer :: dvscfin(:,:,:), dvscfins (:,:,:)
  complex(kind=DP), allocatable :: ldos (:,:), ldoss (:,:),&
       drhoscfh (:,:,:), dvscfout (:,:,:),  &
       dbecsum (:,:,:,:), spsi (:), auxg (:), aux1 (:), ps (:)
  ! local density of states af Ef
  ! local density of states af Ef (not a
  ! change of scf potential (input
  ! change of scf potential (input
  ! change of scf potential (output
  ! change of scf potential (smooth
  ! the derivative of becsum
  ! the function spsi
  ! the function spsi
  ! an auxiliary smooth mesh
  ! the scalar products
  complex(kind=DP) :: ZDOTC
  ! the scalar product function

  logical :: conv_root, exst, lmetq0
  ! true if linter is converged
  ! used to open the recover file
  ! true if xq=(0,0,0) in a metal

  integer :: kter, ipert, ibnd, jbnd, iter, lter, ltaver, &
       lintercall, ik, ikk, ikq, ig, ir, is, nrec, nrec1, ios, mode
  ! counter on iterations
  ! counter on perturbations
  ! counter on bands
  ! counter on bands
  ! counter on iterations
  ! counter on iterations of linter
  ! average counter
  ! average number of call to linter
  ! counter on k points
  ! counter on k points
  ! counter on k+q points
  ! counter on G vectors
  ! counter on mesh points
  ! counter on spin polarizations
  ! the record number
  ! the record number for dpsi
  ! integer variable for I/O control
  ! mode index

  real(kind=DP) :: tcpu, get_clock
  ! timing variables

  character (len=42) :: flmixdpot
  ! the name of the file with the mixing potential

  external ch_psi_all, cg_psi
  !
  call start_clock ('solve_linter')
  allocate (ps ( nbnd))    
  allocate (dvscfin ( nrxx , nspin , npe))    
  if (doublegrid) then
     allocate (dvscfins ( nrxxs , nspin , npe))    
  else
     dvscfins => dvscfin
  endif
  allocate (drhoscfh ( nrxx , nspin , npe))    
  allocate (dvscfout ( nrxx , nspin , npe))    
  allocate (auxg (npwx))    
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin , npe))    
  allocate (aux1 ( nrxxs))    
  allocate (spsi ( npwx))    
  allocate (h_diag ( npwx , nbnd))    
  allocate (eprec ( nbnd))    
  !
  !    if this is a recover run
  !
  if (iter0.ne.0) then
     if (okvan) read (iunrec) int3
     read (iunrec) dr2, dvscfin
     close (unit = iunrec, status = 'keep')
     if (doublegrid) then
        do is = 1, nspin
           do ipert = 1, npe
              call cinterpolate (dvscfin (1, is, ipert), dvscfins (1, is, &
                   ipert), - 1)
           enddo
        enddo
     endif
  endif
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !

  lmetq0 = degauss.ne.0.d0.and.lgamma
  if (lmetq0) then
     allocate (ldos ( nrxx   , nspin))    
     allocate (ldoss(  nrxxs , nspin))    
     call localdos (ldos, ldoss   , dos_ef)
  endif
  !
  !   The outside loop is over the iterations
  !
  if (reduce_io) then
     flmixdpot = ' '
  else
     flmixdpot = 'flmixdpot'
  endif
  do kter = 1, niter_ph
     iter = kter + iter0

     convt = .true.
     ltaver = 0

     lintercall = 0
     call setv (2 * nrxx * nspin * npe, 0.d0, drhoscf, 1)
     call setv (npe * nhm * (nhm + 1) * nat * nspin, 0.d0, dbecsum, 1)
     !
     if (nksq.gt.1) rewind (unit = iunigk)
     do ik = 1, nksq
        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_linter', 'reading igk', abs (ios) )
        endif
        if (lgamma) then
           ikk = ik
           ikq = ik
           npwq = npw
        else
           ikk = 2 * ik - 1
           ikq = ikk + 1
        endif
        if (lsda) current_spin = isk (ikk)
        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('solve_linter', 'reading igkq', abs (ios) )

        endif
        call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
        !
        ! reads unperturbed wavefuctions psi(k) and psi(k+q)
        !
        if (nksq.gt.1) then
           if (lgamma) then
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
           else
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
              call davcio (evq, lrwfc, iuwfc, ikq, - 1)
           endif

        endif
        !
        ! compute the kinetic energy
        !
        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                          (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                          (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
        enddo
        !
        ! diagonal elements of the unperturbed hamiltonian
        !
        do ipert = 1, npert (irr)
           mode = imode0 + ipert
           nrec = (ipert - 1) * nksq + ik
           !
           !  and now adds the contribution of the self consistent term
           !
           if (iter.eq.1) then
              !
              !  At the first iteration dpsi and dvscfin are set to zero,
              !  dvbare_q*psi_kpoint is calculated and written to file
              !
              call setv (2 * nbnd * npwx, 0.d0, dpsi, 1)
              call setv (2 * nrxx * nspin, 0.d0, dvscfin (1, 1, ipert), 1)
              call dvqpsi_us (ik, mode, u (1, mode),.false. )
              call davcio (dvpsi, lrbar, iubar, nrec, 1)
              !
              ! starting threshold for the iterative solution of the linear sistem (li
              !
              thresh = 1.0d-2
           else
              !
              !  After the first iteration dvbare_q*psi_kpoint is read from file
              !
              call davcio (dvpsi, lrbar, iubar, nrec, - 1)
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              call start_clock ('vpsifft')
              do ibnd = 1, nbnd_occ (ikk)
                 call setv (2 * nrxxs, 0.d0, aux1, 1)
                 do ig = 1, npw
                    aux1 (nls (igk (ig) ) ) = evc (ig, ibnd)
                 enddo
                 call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
                 do ir = 1, nrxxs
                    aux1 (ir) = aux1 (ir) * dvscfins (ir, current_spin, ipert)
                 enddo
                 call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
                 do ig = 1, npwq
                    dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + aux1(nls(igkq(ig)))
                 enddo
              enddo
              call stop_clock ('vpsifft')
              !
              !  In the case of US pseudopotentials there is an additional 
              !  selfconsist term which comes from the dependence of D on 
              !  V_{eff} on the bare change of the potential
              !
              call adddvscf (ipert, ik)
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
              !
              ! starting value for delta_psi is read from iudwf
              !
              nrec1 = (ipert - 1) * nksq + ik
              if (nksq.gt.1.or.npert (irr) .gt.1.or.kter.eq.1) call davcio ( &
                   dpsi, lrdwf, iudwf, nrec1, - 1)
           endif
           !
           ! Ortogonalize dvpsi
           !
           call start_clock ('ortho')
           do ibnd = 1, nbnd_occ (ikk)
              if (degauss.ne.0.d0) then
                 wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
                 w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
              endif
              call setv (2 * npwx, 0.d0, auxg, 1)
              do jbnd = 1, nbnd
                 if (degauss.ne.0.d0) then
!  metals
                    wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
                    deltae = et (jbnd, ikq) - et (ibnd, ikk)
                    theta = wgauss (deltae / degauss, 0)
                    wwg = wg1 * (1.d0 - theta) + wgp * theta
                    if (jbnd.le.nbnd_occ (ikq) ) then
                       if (abs (deltae) .gt.1.0d-5) then
                          wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
                       else
                          !
                          !  if the two energies are too close takes the limit
                          !  of the 0/0 ratio
                          !
                          wwg = wwg - alpha_pv * theta * w0g
                       endif
                    endif
                 else
!  insulators
                    if (jbnd.le.nint (nelec) / degspin) then
                       wwg = 1.0d0
                    else
                       wwg = 0.0d0
                    endif
                 endif
                 ps(jbnd) = - wwg * ZDOTC(npwq,evq(1,jbnd),1,dvpsi(1,ibnd),1)
              enddo
#ifdef __PARA
              call reduce (2 * nbnd, ps)
#endif
              do jbnd = 1, nbnd
                 call ZAXPY (npwq, ps (jbnd), evq (1, jbnd), 1, auxg, 1)
              enddo
              if (degauss.ne.0.d0) call DSCAL (2*npwq, wg1, dvpsi(1,ibnd), 1)
              call ZCOPY (npwq, auxg, 1, spsi, 1)
              !
              !   In the US case at the end we have to apply the S matrix
              !
              call ccalbec (nkb, npwx, npwq, 1, becp, vkb, auxg)
              call s_psi (npwx, npwq, 1, auxg, spsi)
              call DAXPY (2 * npwq, 1.0d0, spsi, 1, dvpsi (1, ibnd), 1)
           enddo
           call stop_clock ('ortho')
           !
           !    Here we change the sign of the known term
           !
           call DSCAL (2 * npwx * nbnd, - 1.d0, dvpsi, 1)
           !
           ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
           ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           do ibnd = 1, nbnd_occ (ikk)
              conv_root = .true.
              do ig = 1, npwq
                 auxg (ig) = g2kin (ig) * evq (ig, ibnd)
              enddo
              eprec (ibnd) = 1.35d0 * ZDOTC (npwq, evq (1, ibnd), 1, auxg, 1)
           enddo
#ifdef __PARA
           call reduce (nbnd_occ (ikk), eprec)
#endif
           do ibnd = 1, nbnd_occ (ikk)
              do ig = 1, npwq
                 h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
              enddo
           enddo
           conv_root = .true.

           call cgsolve_all (ch_psi_all,cg_psi,et(1,ikk),dvpsi,dpsi, &
                h_diag,npwx,npwq,thresh,ik,lter,conv_root,anorm,nbnd_occ(ikk))

           ltaver = ltaver + lter
           lintercall = lintercall + 1

           if (.not.conv_root) write (6, '(5x,"kpoint",i4," ibnd",i4,  &
                &              " linter: root not converged ",e10.3)') &
                &              ik , ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec1 = (ipert - 1) * nksq + ik
           !               if (nksq.gt.1 .or. npert(irr).gt.1)
           call davcio (dpsi, lrdwf, iudwf, nrec1, + 1)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ikk)
           call incdrhoscf (drhoscf (1, current_spin, ipert), weight, ik, &
                dbecsum (1, 1, current_spin, ipert), mode, 1)
           ! on perturbations
        enddo
        ! on k-points
     enddo
#ifdef __PARA
     call reduce (nhm * (nhm + 1) * nat * nspin * npe, dbecsum)
#endif
     if (doublegrid) then
        do is = 1, nspin
           do ipert = 1, npert (irr)
              call cinterpolate (drhoscfh (1, is, ipert), drhoscf (1, is, &
                   ipert), 1)
           enddo
        enddo
     else
        call ZCOPY (npe * nspin * nrxx, drhoscf, 1, drhoscfh, 1)
     endif
     !
     !    Now we compute for all perturbations the total charge and potential
     !

     call addusddens (drhoscfh, dbecsum, irr, imode0, npe, 0)
#ifdef __PARA
     !
     !   Reduce the delta rho across pools
     !
     call poolreduce (2 * npe * nspin * nrxx, drhoscf)
     call poolreduce (2 * npe * nspin * nrxx, drhoscfh)
#endif
     !
     ! q=0 in metallic case deserve special care (e_Fermi can shift)
     !

     if (lmetq0) call ef_shift(drhoscfh, ldos, ldoss, dos_ef, irr, npe, .false.)
     do ipert = 1, npert (irr)
        if (fildrho.ne.' ') call davcio_drho (drhoscfh (1, 1, ipert) , &
             lrdrho, iudrho, imode0 + ipert, + 1)
        call ZCOPY (nrxx*nspin, drhoscfh(1,1,ipert), 1, dvscfout(1,1,ipert), 1)
        call dv_of_drho (imode0 + ipert, dvscfout (1, 1, ipert), .true.)

     enddo
     !
     !   After the loop over the perturbations we have the change of the pote
     !   for all the modes of this representation. We symmetrize this potenti
     !
#ifdef __PARA
     call psymdvscf (npert (irr), irr, dvscfout)
#else
     call symdvscf (npert (irr), irr, dvscfout)
#endif
     !
     !   And we mix with the old potential
     !

     call mix_potential (2 * npert (irr) * nrxx * nspin, dvscfout, &
          dvscfin, alpha_mix (kter), dr2, npert (irr) * tr2_ph, iter, &
          nmix_ph, flmixdpot, convt)

     if (lmetq0.and.convt) call ef_shift (drhoscf, ldos, ldoss, dos_ef, &
          irr, npe, .true.)
     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin
              call cinterpolate (dvscfin(1,is,ipert), dvscfins(1,is,ipert), -1)
           enddo
        enddo
     endif
     !
     !     with the new change of the potential we compute the integrals
     !     of the change of potential and Q
     !

     call newdq (dvscfin, npe)
#ifdef __PARA
     aux_avg (1) = dfloat (ltaver)
     aux_avg (2) = dfloat (lintercall)
     call poolreduce (2, aux_avg)
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = dfloat (ltaver) / lintercall
#endif
     tcpu = get_clock ('PHONON')

     write (6, '(//,5x," iter # ",i3," total cpu time : ",f7.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / npert (irr)
     write (6, '(5x," thresh=",e10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",e10.3 )') thresh, alpha_mix (kter) , dr2
     !
     !    Here we save the informations for recovering the run from this poin
     !
#ifdef FLUSH
     call flush (6)
#endif
     call start_clock ('write_rec')
     call seqopn (iunrec, 'recover', 'unformatted', exst)
     if (okvan) write (iunrec) int1, int2

     write (iunrec) dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, zstarue0
     if (reduce_io) then
        write (iunrec) irr, 0, convt, done_irr, comp_irr, ifat
     else
        !
        ! save recover information for current iteration, if available
        !
        write (iunrec) irr, iter, convt, done_irr, comp_irr, ifat
        if (okvan) write (iunrec) int3
        write (iunrec) dr2, dvscfin
     endif
     close (unit = iunrec, status = 'keep')

     call stop_clock ('write_rec')
     if (convt.or.tcpu.gt.time_max) goto 155

  enddo
155 continue
  if (tcpu.gt.time_max.and..not.convt) then
     write (6, '(/,5x,"Stopping for time limit ",2f10.0)') tcpu, time_max
     call stop_ph (.false.)
  endif
  !
  !    There is a part of the dynamical matrix which requires the integral
  !    self consistent change of the potential and the variation of the ch
  !    due to the displacement of the atoms. We compute it here because ou
  !    this routine the change of the self-consistent potential is lost.
  !
  if (convt) then
     call drhodvus (irr, imode0, dvscfin, npe)
     if (fildvscf.ne.' ') write (iudvscf) dvscfins
     if (elph) call elphon (npe, imode0, dvscfins)
  endif
  if (convt.and.nlcc_any) call addnlcc (imode0, drhoscfh, npe)
  if (lmetq0) deallocate (ldoss)
  if (lmetq0) deallocate (ldos)
  deallocate (eprec)
  deallocate (h_diag)
  deallocate (spsi)
  deallocate (aux1)
  deallocate (dbecsum)
  deallocate (auxg)
  deallocate (dvscfout)
  deallocate (drhoscfh)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  deallocate (ps)

  call stop_clock ('solve_linter')
  return
end subroutine solve_linter
