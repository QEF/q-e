!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine work_function (wf)
  !
  ! Print out the workfunction, calculated as the difference between the
  ! potential energy and the fermi energy.
  ! Written for supercells with the main axis along z.
  !
  USE constants, ONLY: rytoev, e2
  USE io_global, ONLY : stdout
  USE ener,      ONLY: ef
  USE lsda_mod,  ONLY: nspin, current_spin
  USE scf,       ONLY: rho, vltot, vr, rho_core
  USE gvect
  USE cell_base, ONLY : omega, alat
#ifdef __PARA
  use para
  use mp
#endif
  implicit none

  integer :: ionode_id = 0
  real(kind=DP) :: wmean1, wmean2, meancharge, wx1, wx2, wxm, vx, vc, ex, &
       ec, rhox, rs, vcca, wf, etxc, vtxc
  integer :: n1, n2, ni, nmean
  logical :: exst
  real(kind=DP), allocatable :: raux1 (:), vaux1 (:), vaux2(:), aux (:)
  real(kind=DP), allocatable :: vxc(:,:)
  ! auxiliary vectors for charge and potential

  allocate (raux1( nrx1 * nrx2 * nrx3))    
  allocate (vaux1( nrx1 * nrx2 * nrx3))    
  allocate (vaux2( nrx1 * nrx2 * nrx3))    

  if (nspin .eq. 4) &
     call errore ('work_function','spin-orbit not implemented',1)

  allocate (vxc(nrxx,nspin))
  call v_xc (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
       nl, ngm, g, nspin, alat, omega, etxc, vtxc, vxc)

#ifdef __PARA
  if (me.eq.1.and.mypool.eq.1) then
#endif
  call seqopn (17, 'workf', 'formatted', exst)
  call seqopn (19, 'charge', 'formatted', exst)
#ifdef __PARA
  endif
#endif

  wf = 0.d0

  do current_spin=1,nspin

#ifdef __PARA
     allocate (aux  ( nrxx))    
     aux(:) = rho(:,current_spin) + rho_core(:)/nspin
     call gather (aux, raux1)
#else
     raux1(1:nrxx) = rho(1:nrxx,current_spin) + rho_core(1:nrxx)/nspin
#endif
     !
     !
#ifdef __PARA
     aux(:) = vltot(:) + vr(:,current_spin)
     call gather (aux, vaux1)
     aux(:) = aux(:) - vxc(:,current_spin)
     call gather (aux, vaux2)
#else
     vaux1(1:nrxx) = vltot(1:nrxx) + vr(1:nrxx,current_spin)
     vaux2(1:nrxx) = vaux1(1:nrxx) -vxc(1:nrxx,current_spin)
#endif
     !
#ifdef __PARA
     deallocate(aux)
     if (me.eq.1.and.mypool.eq.1) then
#endif
     if (nspin.eq.2) then
        if (current_spin.eq.1) then
           WRITE(17,*) " SPIN UP "
           WRITE(19,*) " SPIN UP "
        else
           WRITE(17,*) " SPIN DOWN "
           WRITE(19,*) " SPIN DOWN "
         end if
     endif
     do nmean = 1, nr3
        wmean1 = 0.d0
        wmean2 = 0.d0
        meancharge = 0.d0
        wx1 = 0.d0
        wx2 = 0.d0
        wxm = 0.d0
        do n2 = 1, nr2
           do n1 = 1, nr1
              ni = n1 + (n2 - 1) * nrx1 + (nmean - 1) * nrx1 * nrx2
              meancharge = meancharge+raux1 (ni)
              wxm = wxm + raux1 (ni) **2
              wmean1 = wmean1 + vaux1 (ni)
              wx1 = wx1 + vaux1 (ni) **2
              wmean2 = wmean2 + vaux2 (ni) 
              wx2 = wx2 + vaux2 (ni) **2
           enddo
        enddo
        wmean1 = wmean1 / dble (nr1 * nr2)
        wmean2 = wmean2 / dble (nr1 * nr2)
        meancharge = meancharge / dble (nr1 * nr2)
        wx1 = dsqrt (wx1 / dble (nr1 * nr2) - wmean1 * wmean1)
        wx2 = dsqrt (wx2 / dble (nr1 * nr2) - wmean2 * wmean2)
        wxm = dsqrt (wxm / dble (nr1 * nr2) - meancharge**2)
        if (nmean.eq. (nr3 + 1) / 2) then
           wf = wf + (wmean2 - ef)
           if (nspin.eq.2) then
              if (current_spin.eq.1) then
                 WRITE( stdout,*) " SPIN UP "
              else
                 WRITE( stdout,*) " SPIN DOWN "
              end if
           endif
           WRITE( stdout, 9130) rytoev * (wmean1 - ef), wx1 * rytoev, &
                                rytoev * (wmean2 - ef), wx2 * rytoev
        end if
        write (17, * ) nmean, (wmean1 - ef) * rytoev, wx1 * rytoev, &
                              (wmean2 - ef) * rytoev, wx2 * rytoev
        write (19, * ) nmean, meancharge, wxm
     enddo
#ifdef __PARA
     endif
#endif
  end do
  wf = wf / nspin

#ifdef __PARA
     CALL mp_bcast( wf, ionode_id )
#endif

  WRITE( stdout, '(/5x,"Work function written on file workf")')
  WRITE( stdout, '( 5x,"Planar mean charge written on file charge")')

9130 format (/'     workfunction     = ',f10.4,' +- ',f10.4,' eV', &
    &        /'     without exchcorr = ',f10.4,' +- ',f10.4,' eV')

  close (17)
  close (19)

  deallocate(raux1)
  deallocate(vaux1)
  deallocate(vaux2)
  deallocate(vxc)

  return

end subroutine work_function
