!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine solve_e2
  !-----------------------------------------------------------------------
  !
  !   Self consistent cycle to compute the second order derivatives
  !   of the wavefunctions with respect to electric fields
  !
  use kinds, only : DP
  USE io_global, ONLY : stdout
  use pwcom
  use becmod
  USE control_flags, ONLY: reduce_io
  USE io_files, ONLY: prefix, iunigk
  USE ions_base, ONLY: nat
  USE uspp_param, ONLY : nhm
  USE wavefunctions_module,  ONLY: evc
  USE phcom
  USE ramanm
  USE check_stop, ONLY: check_stop_now
  implicit none

  real(DP) ::  thresh, weight, avg_iter, dr2
  ! convergence threshold for the solution of the
  ! linear system
  ! used for summation over k points
  ! average number of iterations
  ! convergence limit

  complex(DP) , pointer :: dvscfin (:,:,:), dvscfins (:,:,:)
  ! change of the scf potential (input)
  ! change of the scf potential (smooth)

  complex(DP) , allocatable :: dvscfout (:,:,:), dbecsum (:,:), &
                                    aux1 (:)
  ! change of the scf potential (output)
  ! auxiliary space
  ! auxiliary space

  complex(DP) :: ZDOTC
  ! the scalar product function

  logical :: exst
  ! used to open the recover file

  integer :: kter, iter0, ipol, ibnd, iter, ik, is, ig, iig, irr, ir, nrec, ios
  ! counter on iterations
  ! counter on perturbations
  ! counter on bands
  ! counter on iterations
  ! counter on k points
  ! counter on G vectors
  ! counter on g vectors
  ! counter on mesh points
  ! the record number
  ! integer variable for I/O control

  character (len=256) :: flmixdpot
  ! the name of the file with the
  ! mixing potential

  external ch_psi_all, cg_psi

  if (lsda) call errore ('solve_e2', ' LSDA not implemented', 1)
  if (okvan) call errore ('solve_e2', ' Ultrasoft PP not implemented', 1)

  call start_clock('solve_e2')
  allocate (dvscfin( nrxx, nspin, 6))
  if (doublegrid) then
     allocate (dvscfins(  nrxxs, nspin, 6))
  else
     dvscfins => dvscfin
  endif
  allocate (dvscfout( nrxx , nspin, 6))
  allocate (dbecsum( nhm*(nhm+1)/2, nat))
  allocate (aux1(  nrxxs))
  if (irr0 == -10) then
     ! restarting in Raman
     read (iunrec) iter0, dr2
     read (iunrec) dvscfin
     if (okvan) read (iunrec) int1, int2, int3 
     close (unit = iunrec, status = 'keep')
     if (doublegrid) then
        do is = 1, nspin
           do ipol = 1, 6
              call cinterpolate (dvscfin (1, is, ipol),      &
                                 dvscfins (1, is, ipol), -1)
           enddo
        enddo
     end if
  else
     iter0 = 0
  end if
  !
  if (degauss.ne.0.d0 .or..not.lgamma) &
        call errore ('solve_e2', 'called in the wrong case', 1)
  !
  !   The outside loop is over the iterations 
  !
  if (reduce_io) then
     flmixdpot = ' '
  else
     flmixdpot = 'mixd'
  endif

  do kter = 1, niter_ph

     iter = kter + iter0
     avg_iter = 0.d0

     dvscfout (:,:,:) = (0.d0, 0.d0)
     dbecsum (:,:) = (0.d0, 0.d0)
     if (nksq.gt.1) rewind (unit = iunigk)

     do ik = 1, nksq
        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_e', 'reading igk', abs (ios) )
        endif
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, -1)
        npwq = npw
        call init_us_2 (npw, igk, xk (1, ik), vkb)
        !
        ! compute the kinetic energy
        !
        do ig = 1, npwq
           iig = igkq (ig)
           g2kin (ig) = ( (xk (1, ik) + g (1, iig) ) **2 + &
                          (xk (2, ik) + g (2, iig) ) **2 + &
                          (xk (3, ik) + g (3, iig) ) **2 ) * tpiba2
        enddo
        !
        ! The counter on the polarizations runs only on the 6 inequivalent
        ! indexes --see the comment on raman.F--
        !
        do ipol = 1, 6
           nrec = (ipol - 1) * nksq + ik

           if (kter.eq.1) then
              dpsi (:,:) = (0.d0, 0.d0)
           else
              call davcio (dpsi, lrd2w, iud2w, nrec, -1)
           endif

           if (iter.eq.1) then
              dvscfin (:,:,:) = (0.d0, 0.d0)
              call davcio (dvpsi, lrba2, iuba2, nrec, -1)
              thresh = 1.0d-2
           else
              call davcio (dvpsi, lrba2, iuba2, nrec, -1)
              do ibnd = 1, nbnd_occ (ik)
                 call cft_wave (evc (1, ibnd), aux1, +1)
                 do ir = 1, nrxxs
                    aux1 (ir) = aux1 (ir) * dvscfins (ir, 1, ipol)
                 enddo
                 call cft_wave (dvpsi (1, ibnd), aux1, -1)
              enddo
              thresh = min (0.1d0 * sqrt(dr2), 1.0d-2)
           endif

           call pcgreen (avg_iter, thresh, ik, et (1, ik) )
           call davcio ( dpsi, lrd2w, iud2w, nrec, +1)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ik)
           call incdrhoscf (dvscfout (1,1,ipol), weight, ik, &
                            dbecsum (1, 1), 1)
           enddo   ! on perturbations
        enddo      ! on k points
#ifdef __PARA
     call reduce (nhm * (nhm + 1) * nat, dbecsum)
#endif
     if (doublegrid) then
        do is = 1, nspin
           do ipol = 1, 6
              call cinterpolate (dvscfout (1, is, ipol),     &
                                 dvscfout (1, is, ipol), 1)
           enddo
        enddo
     endif

!     call addusddense (dvscfout, dbecsum)

     !
     !   After the loop over the perturbations we have the change of the pote
     !   for all the modes, and we symmetrize this potential
     !
#ifdef __PARA
     call poolreduce (2 * 6 * nrxx * nspin, dvscfout)
#endif
     do ipol = 1, 6
        call dv_of_drho (0, dvscfout (1, 1, ipol), .false.)
     enddo

#ifdef __PARA
     call psyme2(dvscfout)
#else
     call syme2(dvscfout)
#endif
     !
     ! Mixing with the old potential
     !
     call mix_potential (2 * 6 * nrxx * nspin, dvscfout, dvscfin,  &
                         alpha_mix (kter), dr2, 6 * tr2_ph, kter,  &
                         nmix_ph, flmixdpot, convt)

     if (doublegrid) then
        do is = 1, nspin
           do ipol = 1, 6
              call cinterpolate (dvscfin (1, is, ipol), &
                                 dvscfins (1, is, ipol), -1)
           enddo
        enddo
     end if

     write (6, "(//,5x,' iter # ',i3, &
          &      '   av.it.: ',f5.1)") iter, avg_iter / (6.d0 * nksq)
     dr2 = dr2 / 6
     write (6, "(5x,' thresh=',e10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',e10.3 )") thresh, alpha_mix (kter), dr2
     !
     CALL flush_unit( stdout )
     !
     call seqopn (iunrec, 'recover', 'unformatted', exst)
     !
     ! irr: state of the calculation
     ! irr=-10  to -19 Raman
     irr = -10
     !
     write (iunrec) irr
     !
     ! partially calculated results
     !
     write (iunrec) dyn, dyn00
     write (iunrec)  epsilon, zstareu, zstarue, zstareu0, zstarue0
     !
     ! info on current iteration (iter=0 potential mixing not available)
     !
     if (reduce_io.or.convt) then
        write (iunrec) 0, dr2
     else
        write (iunrec) iter, dr2
     end if
     write (iunrec) dvscfin
     if (okvan) write (iunrec) int1, int2, int3
     close (unit = iunrec, status = 'keep')

     if ( check_stop_now() ) then
        call stop_ph (.false.)
        goto 155
     endif
     if ( convt ) goto 155

  enddo
155 continue
  deallocate (dvscfin )
  if (doublegrid) deallocate (dvscfins )
  deallocate (dvscfout )
  deallocate (dbecsum )
  deallocate (aux1 )

  call stop_clock('solve_e2')

  return 
end subroutine solve_e2
