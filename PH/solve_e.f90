!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e  
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver for the solution of the linear system whic
  !    defines the change of the wavefunction due to the electric field.
  !    It performs the following tasks:
  !     a) It computes the kinetic energy
  !     b) It adds the term Delta V_{SCF} | psi >
  !     c) It applies P_c to the known part
  !     d) It calls linter to solve the linear system
  !     e) It computes Delta rho, Delta V_{SCF} and symmetrize them
  !
#include "machine.h"
  use pwcom 
  use allocate 
  use parameters, only : DP 
  use phcom
  implicit none 

  real(kind=DP) ::  thresh, weight, anorm, averlt, dr2
  real(kind=DP),pointer :: h_diag (:,:), eprec(:)
  ! the diagonal part of the Hamiltonia
  ! the convergence threshold
  ! used for summation over k points
  ! the norm of the error
  ! average number of iterations
  ! cut-off for preconditioning
  ! convergence limit

  complex(kind=DP) , pointer :: dvscfin (:,:), dvscfins (:,:),&
       dvscfout (:,:), auxg (:), aux1 (:), ps (:)
  ! change of the scf potential (input
  ! change of the scf potential (smoot
  ! change of the scf potential (outpu
  ! auxiliary space
  ! the psi function
  ! the scalar product
  complex(kind=DP) :: dbecsum, ZDOTC
  ! dummy
  ! the scalar product function
  logical :: conv_root, exst  
  ! true if linter is converged
  ! used to open the recover file

  integer :: kter, ipol, ibnd, jbnd, iter, lter, ltaver, lintercall, &
       ik, ig, irr, ir, nrec, nrec1, ios
  ! counter on iterations
  ! counter on perturbations
  ! counter on bands
  ! counter on bands
  ! counter on iterations
  ! counter on iterations of linter
  ! average counter
  ! average number of call to linter
  ! counter on k points
  ! counter on G vectors
  ! the irreducible representation
  ! counter on g vectors
  ! counter on mesh points
  ! the record number
  ! the record number for dpsi
  ! integer variable for I/O control

  real(kind=DP) :: tcpu, get_clock  
  ! timing variables

  character (len=42) :: flmixdpot  
  ! the name of the file with the
  ! mixing potential

  external ch_psi_all, cg_psi  

  if (lsda) call error ('solve_e', ' LSDA not implemented', 1)  

  call start_clock ('solve_e')  
  call mallocate(dvscfin, nrxx , 3)  
  if (doublegrid) then  
     call mallocate(dvscfins,  nrxxs , 3)  
  else  
     dvscfins => dvscfin  
  endif
  call mallocate(dvscfout, nrxx , 3)  
  call mallocate(auxg ,  npwx)  
  call mallocate(aux1 ,  nrxxs)  
  call mallocate(ps ,  nbnd)  
  call mallocate(h_diag , npwx , nbnd)  
  call mallocate(eprec ,  nbnd)  
  if (iter0.ne.0) then  
     read (iunrec) dr2, dvscfin  
     close (unit = iunrec, status = 'keep')  
     if (doublegrid) then  
        do ipol = 1, 3  
           call cinterpolate (dvscfin (1, ipol), dvscfins (1, ipol), &
                - 1)
        enddo
     endif

  endif
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if (degauss.ne.0.d0.or..not.lgamma) call error ('solve_e', &
       'called in the wrong case', 1)
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

     call setv (2 * nrxx * 3, 0.d0, dvscfout, 1)  

     if (nksq.gt.1) rewind (unit = iunigk)  
     do ik = 1, nksq  
        if (nksq.gt.1) then  
           read (iunigk, err = 100, iostat = ios) npw, igk  
100        call error ('solve_e', 'reading igk', abs (ios) )  
        endif
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)  
        npwq = npw  
        call init_us_2 (npw, igk, xk (1, ik), vkb)  
        !
        ! compute the kinetic energy
        !
        do ig = 1, npwq 
           g2kin (ig) = ( (xk (1,ik ) + g (1,igkq (ig)) ) **2 + &
                          (xk (2,ik ) + g (2,igkq (ig)) ) **2 + &
                          (xk (3,ik ) + g (3,igkq (ig)) ) **2 ) * tpiba2
        enddo
        !
        do ipol = 1, 3  
           nrec = (ipol - 1) * nksq + ik  
           !
           !  and now adds the contribution of the self consistent term
           !
           if (iter.eq.1) then  
              !
              !  At the first iteration dpsi and dvscfin are set to zero,
              !  [H,x]*psi_kpoint is calculated and written to file
              !
              call setv (2 * nbnd * npwx, 0.d0, dpsi, 1)  
              call setv (2 * nrxx, 0.d0, dvscfin (1, ipol), 1)  
              call dvpsi_e (ik, ipol)  
              call davcio (dvpsi, lrbar, iubar, nrec, 1)  
              !
              ! starting threshold for the iterative solution of the linear sistem (li
              !
              thresh = 1.d-2  
           else  
              !
              !  After the first iteration [H,x]*psi_kpoint is read from file
              !
              call davcio (dvpsi, lrbar, iubar, nrec, - 1)  
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              do ibnd = 1, nbnd_occ (ik)  
                 call setv (2 * nrxxs, 0.d0, aux1, 1)  
                 do ig = 1, npw  
                    aux1 (nls (igk (ig) ) ) = evc (ig, ibnd)  
                 enddo
                 call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, &
                      + 2)
                 do ir = 1, nrxxs  
                    aux1 (ir) = aux1 (ir) * dvscfins (ir, ipol)  
                 enddo
                 call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, &
                      - 2)
                 do ig = 1, npwq  
                    dvpsi (ig, ibnd) = dvpsi (ig, ibnd) - aux1 (nls (igkq (ig) ) )  
                 enddo
              enddo
              !
              ! starting value for  delta_psi is read from iudwf
              !
              nrec1 = (ipol - 1) * nksq + ik  
              call davcio (dpsi, lrdwf, iudwf, nrec1, - 1)  
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (0.1d0 * sqrt (dr2), 1.0d-2)  
           endif
           !
           ! Orthogonalize dvpsi
           !
           do ibnd = 1, nbnd_occ (ik)  
              call setv (2 * npwx, 0.d0, auxg, 1)  
              do jbnd = 1, nbnd_occ (ik)  
                 ps (jbnd) = - ZDOTC (npwq, evc (1, jbnd), 1, dvpsi (1, ibnd), &
                      1)
              enddo
#ifdef PARA
              call reduce (2 * nbnd, ps)  
#endif
              do jbnd = 1, nbnd_occ (ik)  
                 call ZAXPY (npwq, ps (jbnd), evc (1, jbnd), 1, auxg, 1)  
              enddo
              call DAXPY (2 * npwq, 1.0d0, auxg, 1, dvpsi (1, ibnd), 1)  
           enddo
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           do ibnd = 1, nbnd_occ (ik)  
              do ig = 1, npwq  
                 auxg (ig) = g2kin (ig) * evc (ig, ibnd)  
              enddo
              eprec (ibnd) = 1.35d0 * ZDOTC (npwq, evc (1, ibnd), 1, auxg, 1)  
           enddo
#ifdef PARA
           call reduce (nbnd_occ (ik), eprec)  
#endif
           do ibnd = 1, nbnd_occ (ik)  
              do ig = 1, npwq  
                 h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )  
              enddo
           enddo

           conv_root = .true.  

           call cgsolve_all (ch_psi_all, cg_psi, et (1, ik), dvpsi, dpsi, &
                h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, nbnd_occ ( &
                ik) )
           ltaver = ltaver + lter  
           lintercall = lintercall + 1  
           if (.not.conv_root) write (6, '(5x,"kpoint",i4," ibnd",i4, &
                &                      " linter: root not converged ",e10.3)') ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec1 = (ipol - 1) * nksq + ik  

           call davcio (dpsi, lrdwf, iudwf, nrec1, + 1)  
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ik)  
           call incdrhoscf (dvscfout (1, ipol), weight, ik, dbecsum, 1, 1)  
        enddo

     enddo
     do ipol = 1, 3  
        if (doublegrid) call cinterpolate (dvscfout (1, ipol), dvscfout ( &
             1, ipol), 1)
        if (fildrho.ne.' ') call davcio_drho (dvscfout (1, ipol) , lrdrho, &
             iudrho, ipol, + 1)
        call dv_of_drho (0, dvscfout (1, ipol), .false.)  

     enddo
     !
     !   After the loop over the perturbations we have the change of the pote
     !   for all the modes of this representation. We symmetrize this potenti
     !
#ifdef PARA
     call poolreduce (2 * 3 * nrxx, dvscfout)  
     call psyme (dvscfout)  
#else
     call syme (dvscfout)  
#endif
     !
     !   And we mix with the old potential
     !

     call mix_potential (2 * 3 * nrxx, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph, iter, nmix_ph, flmixdpot, convt)
     if (doublegrid) then  
        do ipol = 1, 3  
           call cinterpolate (dvscfin (1, ipol), dvscfins (1, ipol), &
                - 1)
        enddo

     endif
     averlt = dfloat (ltaver) / dfloat (lintercall)  
     write (6, '(//,5x," iter # ",i3, &
          &      "   av.it.: ",f5.1)') iter, averlt
     dr2 = dr2 / 3  
     write (6, '(5x," thresh=",e10.3, " alpha_mix = ",f6.3, &
          &               " |ddv_scf|^2 = ",e10.3 )') thresh, alpha_mix (kter) , dr2
#ifdef FLUSH
     call flush (6)  
#endif

     call seqopn (iunrec, 'recover', 'unformatted', exst)  
     irr = - 2  

     write (iunrec) dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, &
          zstarue0
     if (reduce_io) then  
        write (iunrec) irr, 0, convt, done_irr, comp_irr, ifat  
     else  
        write (iunrec) irr, iter, convt, done_irr, comp_irr, ifat  
        write (iunrec) dr2, dvscfin  

     endif

     close (unit = iunrec, status = 'keep')  
     tcpu = get_clock ('PHONON')  
     if (convt.or.tcpu.gt.time_max) goto 155  

  enddo
155 continue  
  if (tcpu.gt.time_max) then  
     write (6, '(/,5x,"Stopping for time limit ",2f10.0)') tcpu, &
          time_max
     call stop_ph (.false.)  

  endif
  call mfree (eprec)  
  call mfree (h_diag)  
  call mfree (ps)  
  call mfree (aux1)  
  call mfree (auxg)  
  call mfree (dvscfout)  
  if (doublegrid) call mfree (dvscfins)  
  call mfree (dvscfin)  

  call stop_clock ('solve_e')  
  return  
end subroutine solve_e
