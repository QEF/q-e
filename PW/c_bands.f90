!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine c_bands (iter, ik_, dr2)
  !-----------------------------------------------------------------------
  !
  !   This routine is a driver for the diagonalization routines of the
  !   total Hamiltonian at each k-point.
  !   It reads the Hamiltonian and an initial guess of the wavefunctions
  !   from a file and computes initialization quantities for the
  !   diagonalization routines.
  !   There are three types of iterative diagonalization:
  !   a) Davidson algorithm (all-band)
  !   b) Conjugate Gradient (band-by-band)
  !   c) DIIS algorithm
  !
#include "machine.h"
  use pwcom, ONLY: g, g2kin, tpiba2, ecfixed, qcutz, lda_plus_u, &
 iunwfc, swfcatom, iunat, nwordatwfc, q2sigma, diis_ndim, wg, nbndx, nkstot, &
 okvan, et, istep, ethr, lscf, max_cg_iter, vltot, nrxx, nr1, nr3, nr2, nbnd, &
 nks, npwx, diis_start_cg, isolve, iunigk, vkb, xk, reduce_io, nwordwfc, evc, &
 current_spin, lsda, isk, igk, npw
  use g_psi_mod
  implicit none
  !
  !     First the I/O variables
  !
  integer :: ik_, iter
  ! k-point already done
  ! current iterations
  real(kind=DP) :: dr2
  ! current accuracy of self-consistency
  !
  !     here the local variables
  !

  real(kind=DP) :: avg_iter, cg_iter, diis_iter, v_of_0, dsum, erf
  ! average number of iterations
  ! number of iteration in CG
  ! number of iteration in DIIS
  ! the average of the potential
  ! summation function
  ! error function

  integer :: ik, ig, ibnd, dav_iter, ntry, notconv
  ! counter on k points
  ! counter on G vectors
  ! counter on bands
  ! number of iterations in Davidson
  ! number or repeated call to diagonalization in case of non convergence
  ! number of notconverged elements
  integer, allocatable :: btype(:)
  ! type of band: conduction (1) or valence (0)


  if (ik_ == nks) then
     ik_ = 0
     return
  endif

  call start_clock ('c_bands')
  !
  !   allocate arrays
  !
  allocate (btype(  nbnd))    
  allocate (h_diag( npwx))
  allocate (s_diag( npwx))    

  if (isolve == 0) then
     write (6, '("     Davidson diagonalization (with overlap)")')
  elseif (isolve == 1) then
     write (6, '("     Conjugate-gradient style diagonalization")')
  elseif (isolve == 2) then
     write (6, '("     DIIS style diagonalization")')
     if (iter < diis_start_cg) &
          write (6, '(6x,i3," of ",i3," CG iterations")') iter,diis_start_cg
  else
     call errore ('c_bands', 'isolve not implemented', 1)

  endif
  avg_iter = 0.d0
  !
  ! v_of_0 is (Vloc)(G=0)
  !
  v_of_0 = dsum (nrxx, vltot, 1) / float (nr1 * nr2 * nr3)
#ifdef __PARA
  call reduce (1, v_of_0)
#endif
  !

  if (nks > 1) rewind (iunigk)
  !
  !    For each k point diagonalizes the hamiltonian
  !
  do ik = 1, nks
     if (lsda) current_spin = isk (ik)
     !
     !   Reads the Hamiltonian and the list k+G <-> G of this k point
     !
     if (nks > 1) read (iunigk) npw, igk
     !
     !   do not recalculate k-points if restored from a previous run
     !
     if (ik <= ik_) goto 20
     !
     !   various initializations
     !
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     !   read in wavefunctions from the previous iteration
     !

     if (nks > 1 .or..not.reduce_io) call davcio(evc,nwordwfc,iunwfc,ik,-1)
     ! Needed for LDA+U
     if (lda_plus_u) call davcio (swfcatom, nwordatwfc, iunat, ik,- 1)
     !
     !    sets the kinetic energy
     !
     do ig = 1, npw
        g2kin (ig) =((xk (1, ik) + g (1, igk (ig) ) ) **2 + &
                     (xk (2, ik) + g (2, igk (ig) ) ) **2 + &
                     (xk (3, ik) + g (3, igk (ig) ) ) **2 ) * tpiba2
     enddo
     !
     if (qcutz > 0.d0) then
        do ig = 1, npw
           g2kin (ig) = g2kin (ig) + qcutz * (1.d0 + erf ( (g2kin (ig) &
                - ecfixed) / q2sigma) )
        enddo
     endif
     !
     if (isolve == 1 .or. &
             (isolve == 2 .and. iter <= diis_start_cg)) then
        !
        ! Conjugate-Gradient diagonalization
        ! and first "diis_start_cg" steps of RMM-DIIS diagonalization
        !
        ! h_diag is the precondition matrix
        !
        do ig = 1, npw
           h_diag (ig) = max (1.d0, g2kin (ig) )
        enddo
        ntry = 0
10      continue
        if (iter.ne.1.or.istep.ne.1.or.ntry.gt.0) then
           call cinitcgg (npwx, npw, nbnd, nbnd, evc, evc, et (1, ik) )
           avg_iter = avg_iter + 1.d0
        endif
        call ccgdiagg (npwx, npw, nbnd, evc, et (1, ik), h_diag, ethr, &
             max_cg_iter, .not.lscf, notconv, cg_iter)
        avg_iter = avg_iter + cg_iter
        !
        !   save wave-functions to be used as input for the iterative
        !   diagonalization of the next scf iteration and for rho calculation
        !
        if (nks > 1 .or..not.reduce_io) call davcio(evc,nwordwfc,iunwfc,ik,1)
        ntry = ntry + 1
        if (ntry.le.5.and. ( &
             .not.lscf.and.notconv.gt.0.or.lscf.and.notconv.gt.5) ) goto 10
     elseif (isolve == 2) then
        !
        !  after "diis_start_cg" steps of CG, start the RMM-DIIS method
        !
        do ig = 1, npw
           h_diag (ig) = g2kin (ig) + v_of_0
        enddo
        call usnldiag (h_diag, s_diag)
        ntry = 0
        diis_iter = 0.d0

        btype (:) = 0
        if (iter > 1) then
           do ibnd = 1, nbnd
              if ( wg(ibnd, ik) < 1.0d-4 ) btype (ibnd) = 1
           end do
           !             write(*,'(5f12.6)')(et(ibnd,ik),ibnd=1,nbnd)
           !             write(*,'(20i3)')(btype(ibnd),ibnd=1,nbnd)
           !
        end if
12      continue
        call cdiisg(npw, npwx, nbnd, diis_ndim, evc, et (1, ik), ethr, &
             btype, notconv, diis_iter, iter)
        avg_iter = avg_iter + diis_iter
        ntry = ntry + 1
        !
        !   save wave-functions to be used as input for the iterative
        !   diagonalization of the next scf iteration and for rho calculation
        !
        if (nks > 1 .or..not.reduce_io) call davcio(evc,nwordwfc,iunwfc,ik,1)
        if (ntry.le.5.and. ( &
             .not.lscf.and.notconv.gt.0.or.lscf.and.notconv.gt.5) ) goto 12

     else
        !
        !   Davidson
        !
        !   h_diag are the diagonal matrix elements of the hamiltonian
        !   used in g_psi to evaluate the correction to the trial eigenvectors
        !
        do ig = 1, npw
           h_diag (ig) = g2kin (ig) + v_of_0
        enddo
        call usnldiag (h_diag, s_diag)
        ntry = 0

15      continue

        call cegterg (npw, npwx, nbnd, nbndx, evc, ethr, okvan, &
                et (1, ik), notconv, dav_iter)
        avg_iter = avg_iter + dav_iter
        !
        !   save wave-functions to be used as input for the iterative
        !   diagonalization of the next scf iteration and for rho calculation
        !
        if (nks > 1 .or..not.reduce_io) call davcio(evc,nwordwfc,iunwfc,ik,1)
        ntry = ntry + 1
        if (ntry.le.5.and. ( &
             .not.lscf.and.notconv.gt.0.or.lscf.and.notconv.gt.5) ) goto 15
     endif
     if (notconv.ne.0) write (6, '(" warning : ",i3," eigenvectors not",&
          &" converged after ",i3," attemps")') notconv, ntry
     if (notconv.gt.max (5, nbnd / 4) ) stop
20   continue
     !
     ! save restart information
     !

     call save_in_cbands (iter, ik, dr2)

  enddo
  ik_ = 0
#ifdef __PARA
  call poolreduce (1, avg_iter)
#endif
  avg_iter = avg_iter / nkstot
  write (6, 9000) ethr, avg_iter
  !
  ! deallocate work space
  !
  deallocate (s_diag)
  deallocate (h_diag)
  deallocate (btype)

  call stop_clock ('c_bands')
  return
9000 format(5x,'ethr = ',1pe9.2,',  avg # of iterations =',0pf5.1 )
end subroutine c_bands

