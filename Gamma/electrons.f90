!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine electrons  
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver of the self-consistent cycle.
  !    It uses the routine c_bands for computing the bands at fixed
  !    Hamiltonian, the routine sum_bands to compute the charge
  !    density, the routine v_of_rho to compute the new potential
  !    and the routine mix_potential to mix input and output
  !    potentials.
  !
  !    It prints on output the total energy and its decomposition in
  !    the separate contributions.
  !
#include "machine.h"
  use pwcom  
  use gamma  
  use io, only: prefix
  !
  !     a few local variables
  !
#ifdef PARA
  use para  
#endif
implicit none
#ifdef PARA
  ! number of plane waves summed on all nodes
  integer :: ngkp (npk)  
#define NRXX ncplane*npp(me)
  ! This is needed in mix_pot whenever nproc is not a divisor of nr3.
#else
#define NRXX nrxx
#endif
  character :: flmix * 42  

  real(kind=DP) :: de, dr2, charge, mag, magtot, absmag, tcpu
  ! the correction energy
  ! the norm of the diffence between potential
  ! the total charge
  ! local magnetization
  ! total magnetization
  ! total absolute magnetization

  integer :: i, ir, ig, ik, ibnd, idum, iter, ik_  
  ! counter on polarization
  ! counter on the mesh points
  ! counter on k points
  ! counter on bands
  ! dummy counter on iterations
  ! counter on iterations
  ! used to read ik   from restart file

  real (kind=DP) :: ehart_new,etxc_new,vtxc_new, charge_new
  real (kind=DP), external :: ewald, get_clock
  logical :: exst  

  call start_clock ('electrons')  
  !
  iter = 0  
  ik_ = 0  
  !
  if (restart) then  
     call restart_in_electrons (iter, ik_, dr2)  
     if (ik_.eq. - 1000) then  
        conv_elec = .true.  
        ! jump to the end
        goto 999  
     endif
  endif
  !
  if (lscf) then  
     !   calculates the ewald contribution to total energy
     ewld = ewald (alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
          gg, ngm, gcutm, gstart, gamma_only, strf)
     if (reduce_io) then  
        flmix = ' '  
     else  
        flmix = 'flmix'  
     endif
  endif
10 continue  
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  if (imix.ge.0) then
     do ig=1,ngm
        if (gg(ig).lt. ecutwfc/tpiba2) ngm0=ig
     end do
     ngm0 = ngm
  end if

  do idum = 1, niter  
     tcpu = get_clock ('PWSCF')  
     write (6, 9000) tcpu  
     if (imix.ge.0) call DCOPY(nspin*nrxx,rho,1,rho_save,1)
     iter = iter + 1  
     if (lscf) then  
        write (6, 9010) iter, ecutwfc, mixing_beta
     else  
        write (6, 9009)  
     endif
#ifdef FLUSH
     call flush (6)  
#endif
     ! Convergence threshold for iterative diagonalization
     if (lscf.and.iter.ne.1.and.ik_.eq.0) then
        if (imix.ge.0) then
           ethr=max(min(ethr,dr2/nelec/10.0),tr2/nelec/100.0)
        else
           ethr = max (min (ethr/2.0,sqrt(dr2)/1000.0), 1.d-12)
        end if
     end if
     !
     call c_bands (iter, ik_, dr2)  
     !
     !! skip all the rest if not lscf
     if (.not.lscf) then  
        conv_elec=.true.
#ifdef PARA
        call poolrecover (et, nbndx, nkstot, nks)  
#endif

        do ik = 1, nkstot  
           if (lsda) then  
              if (ik.eq.1) write (6, 9015)  
              if (ik.eq.1 + nkstot / 2) write (6, 9016)  

           endif
           write (6, 9020) (xk (i, ik), i = 1, 3)  
           write (6, 9030) (et (ibnd, ik) * 13.6058, ibnd = 1, nbnd)  
        enddo
        ! jump to the end
        goto 999  

     endif
     tcpu = get_clock ('PWSCF')  
     if (tcpu.gt.time_max) then  
        write (6, '(5x,"Maximum CPU time exceeded",2f15.2)') tcpu, &
             time_max
        call stop_pw (.false.)  
     endif
     !
     call sum_band  
     !
     ! calculate total and absolute magnetization
     !
     if (lsda) then  
        magtot = 0.0d0  
        absmag = 0.0d0  
        do ir = 1, nrxx  
           mag = rho (ir, 1) - rho (ir, 2)  
           magtot = magtot + mag  
           absmag = absmag + ABS (mag)  
        enddo
        magtot = magtot * omega / (nr1 * nr2 * nr3)  
        absmag = absmag * omega / (nr1 * nr2 * nr3)  
#ifdef PARA
        call reduce (1, magtot)  
        call reduce (1, absmag)  
#endif
     endif
     !
     call v_of_rho (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
          nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
          ehart, etxc, vtxc, charge, vnew)

     call delta_e (nr1, nr2, nr3, nrxx, rho, vr, vnew, omega, de, &
          deband, nspin)
     if (imix.ge.0) then

        call mix_rho (rho, rho_save, nsnew, ns, mixing_beta, dr2, iter, &
                      nmix, flmix, conv_elec)

        call DAXPY(nspin*nrxx,-1.d0,vr,1,vnew,1)

        call v_of_rho &
             (rho_save, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
              nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
              ehart_new, etxc_new, vtxc_new, charge_new, vr )

     else ! old style potential mixing


        call vpack (NRXX, nrxx, nspin, vnew, vr, + 1)  

        call mix_potential (nspin * NRXX, vnew, vr, mixing_beta, dr2, tr2, &
             iter, nmix, flmix, conv_elec)

        call vpack (NRXX, nrxx, nspin, vnew, vr, - 1)  
     end if
     !
     ! On output vnew contains V(out)-V(in). Used to correct the forces
     !
     ! define the total local potential (external + scf)
     !
     call set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)  
     !
     !
     !   In the US case we need to recompute the self consistent term in
     !   the nonlocal potential.
     !
     if (lda_plus_u) then  
        if (iter.gt.niter_with_fixed_ns .and. imix.lt.0) &
            call DCOPY(nat*nspin*25,nsnew,1,ns,1)
#ifdef PARA
        if (me.eq.1.and.mypool.eq.1) then  
#endif
           call seqopn (iunocc, trim(prefix)//'.occup', 'formatted', exst)  
              write (iunocc, * ) ns  
           close (unit = iunocc, status = 'keep')  
#ifdef PARA
        endif
#endif
     endif

     call newd  
     !
     !  write the potential (and rho) on file
     !
     if (imix.ge.0) call io_pot(+1,trim(prefix)//'.rho',rho_save,nspin)
     call io_pot(+1,trim(prefix)//'.pot',vr,nspin)
     !
     !  save converged wfc if they have not been written previously
     !
     if (nks.eq.1.and.reduce_io) call davcio(evc,nwordwfc,iunwfc,nks,1)

     !
     !  write recover file
     !

     call save_in_electrons (iter, dr2)  

     if ( (conv_elec.or.mod(iter,iprint).eq.0).and.iswitch.le.2) then

        if (lda_plus_u) call write_ns  
#ifdef PARA
        do ik = 1, nks  
           ngkp (ik) = ngk (ik)  
        enddo
        call ireduce (nks, ngkp)  
        call ipoolrecover (ngkp, 1, nkstot, nks)  
        call poolrecover (et, nbndx, nkstot, nks)  
#endif

        do ik = 1, nkstot  
           if (lsda) then  
              if (ik.eq.1) write (6, 9015)  
              if (ik.eq.1 + nkstot / 2) write (6, 9016)  

           endif
           if (conv_elec) then  
#ifdef PARA
              write (6, 9021) (xk (i, ik), i = 1, 3), ngkp (ik)  
#else
              write (6, 9021) (xk (i, ik), i = 1, 3), ngk (ik)  
#endif
           else  
              write (6, 9020) (xk (i, ik), i = 1, 3)  
           endif
           write (6, 9030) (et (ibnd, ik) * 13.6058, ibnd = 1, nbnd)  
        enddo
        if (lgauss.or.ltetra) write (6, 9040) ef * 13.6058  
     endif
     if (abs (charge-nelec) / charge.gt.1.0e-7) write (6, 9050) charge  
     etot = eband+ (etxc - etxcc) + ewld+ehart + deband+demet + eth  
     if ( (conv_elec.or.mod (iter, iprint) .eq.0) .and.iswitch.le.2) &
          then
        if (imix.ge.0) then
           write (6, 9081) etot, dr2 
        else
           write (6, 9086) etot, dr2 
        end if
        write (6, 9060) eband, eband+deband, ehart, etxc-etxcc, ewld, de
        if (lda_plus_u) write (6, 9065) eth  
        if (degauss.ne.0.0) write (6, 9070) demet  
     elseif (conv_elec.and.iswitch.gt.2) then  
        if (imix.ge.0) then
           write (6, 9081) etot, dr2 
        else
           write (6, 9086) etot, dr2 
        end if
     else  
        if (imix.ge.0) then
           write (6, 9080) etot, dr2 
        else
           write (6, 9085) etot, dr2 
        end if
     endif
     if (lsda) write (6, 9017) magtot, absmag  
     !
#ifdef FLUSH
     call flush (6)  
#endif
     if (conv_elec) then  
        write (6, 9110)  
        ! jump to the end
        goto 999  
     endif

!
! uncomment the following line if you wish to monitor the evolution of the
! force calculation during self-consistency
!
!     call forces

     if (imix.ge.0) call DCOPY(nspin*nrxx,rho_save,1,rho,1)

  enddo

  write (6, 9120)  
  ! <------- jump here if not scf

999 continue  

  if (output_drho.ne.' ') call remove_atomic_rho  

  call stop_clock ('electrons')  

  return  
9000 format (/'     total cpu time spent up to now is ',f9.2,' secs')  
9009 format (/'     Band Structure Calculation')  
9010 format (/'     iteration #',i3,'     ecut=',f9.2,' ryd',5x, &
       &         'beta=',f4.2)
9015 format (/' ------ SPIN UP ------------'/)  
9016 format (/' ------ SPIN DOWN ----------'/)  
9017 format (/'     total magnetization       =',f9.2,' Bohr mag/cell', &
       &        /'     absolute magnetization    =',f9.2,' Bohr mag/cell')
9020 format (/'          k =',3f7.4,'     band energies (ev):'/)  
9021 format (/'          k =',3f7.4,' (',i5,' PWs)   bands (ev):'/)  
9030 format ( '  ',8f9.4)  
9040 format (/'     the Fermi energy is ',f10.4,' ev')  
9050 format (/'     integrated charge         =',f12.5)  
9060 format (/'     band energy sum           =',  f15.8,' ryd' &
       &     /'     one-electron contribution =',  f15.8,' ryd' &
       &     /'     hartree contribution      =',  f15.8,' ryd' &
       &     /'     xc contribution           =',  f15.8,' ryd' &
       &     /'     ewald contribution        =',  f15.8,' ryd' &
       &     /'     scf in/out correction     =',  f15.8,' ryd' )
9065 format ( '     Hubbard energy            =',f15.8,' ryd')  
9070 format ( '     correction for metals     =',f15.8,' ryd')  
9080 format (/'     total energy              =',0pf15.8,' ryd' &
             /'     estimated scf accuracy    <',0pf15.8,' ryd')
9081 format (/'!    total energy              =',0pf15.8,' ryd' &
             /'     estimated scf accuracy    <',0pf15.8,' ryd')
9085 format (/'     total energy              =',0pf15.8,' ryd' &
             /'     potential mean squ. error =',1pe15.1,' ryd^2')
9086 format (/'!    total energy              =',0pf15.8,' ryd' &
             /'     potential mean squ. error =',1pe15.1,' ryd^2')
9090 format (/'     the final potential is written on file ',a14)  
9100 format (/'     this iteration took ',f9.2,' cpu secs')  
9110 format (/'     convergence has been achieved')  
9120 format (/'     convergence NOT achieved. stopping ...')  
end subroutine electrons

