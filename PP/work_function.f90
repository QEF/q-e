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
  use pwcom
#ifdef __PARA
  use para
  use mp
#endif
  implicit none

  integer :: ionode_id = 0
  real(kind=DP) :: wmean1, wmean2, meancharge, wx1, wx2, wxm, vx, vc, ex, &
       ec, rhox, rs, vcca, wf
  integer :: n1, n2, ni, nmean
  logical :: exst
  real(kind=DP), allocatable :: raux1 (:), vaux1 (:), aux (:)
  ! auxiliary vectors for charge and potential

  allocate (raux1( nrx1 * nrx2 * nrx3))    
  allocate (vaux1( nrx1 * nrx2 * nrx3))    

  ! if (.not.lscf) call sum_band
  ! TEMP: lscf no longer read, igk indices not written to disk
  if (nspin .ne. 1) &
     call errore ('work_function','spin polarization not implemented',1)
  current_spin = 1
#ifdef __PARA
  allocate (aux  ( nrxx))    
  aux(:) = rho(:,current_spin) + rho_core(:)
  call gather (aux, raux1)
#else
  raux1(1:nrxx) = rho(1:nrxx,current_spin) + rho_core(1:nrxx)
#endif
  !
#ifdef __PARA
  aux(:) = vltot(:) + vr(:,current_spin)
  call gather (aux, vaux1)
#else
  vaux1(1:nrxx) = vltot(1:nrxx) + vr(1:nrxx,current_spin)
#endif
  !
#ifdef __PARA
  deallocate(aux)
  if (me.eq.1.and.mypool.eq.1) then
#endif
     call seqopn (17, 'workf', 'formatted', exst)
     call seqopn (19, 'charge', 'formatted', exst)
     nmean = (nr3 + 1) / 2
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
              rhox = abs (raux1 (ni) )
              vx = 0.d0
              if (rhox.gt.1.0e-30) then
                 call xc (rhox, ex, ec, vx, vc)
                 vx = e2 * (vx + vc)
              endif
              wmean2 = wmean2 + vaux1 (ni) - vx
              wx2 = wx2 + (vaux1 (ni) - vx) **2
           enddo
        enddo
        wmean1 = wmean1 / dfloat (nr1 * nr2)
        wmean2 = wmean2 / dfloat (nr1 * nr2)
        meancharge = meancharge / dfloat (nr1 * nr2)
        wx1 = dsqrt (wx1 / dfloat (nr1 * nr2) - wmean1 * wmean1)
        wx2 = dsqrt (wx2 / dfloat (nr1 * nr2) - wmean2 * wmean2)
        wxm = dsqrt (wxm / dfloat (nr1 * nr2) - meancharge**2)
        if (nmean.eq. (nr3 + 1) / 2) wf = wmean2 - ef
        write (17, * ) nmean, (wmean1 - ef) * rytoev, wx1 * rytoev, &
             (wmean2 - ef) * rytoev, wx2 * rytoev
        write (19, * ) nmean, meancharge, wxm
        if (nmean.eq. (nr3 + 1) / 2) then
           write (6, 9130) rytoev * (wmean1 - ef), wx1 * rytoev, &
                rytoev * (wmean2 - ef), wx2 * rytoev
        endif
     enddo
#ifdef __PARA
  endif
  CALL mp_bcast( wf, ionode_id )
#endif
  write (6, '(/5x,"Work function written on file workf")')
  write (6, '( 5x,"Planar mean charge written on file charge")')

9130 format (/'     workfunction     = ',f10.4,' +- ',f6.4,' eV', &
    &        /'     without exchcorr = ',f10.4,' +- ',f6.4,' eV')
  close (17)
  close (19)
  deallocate(raux1)
  deallocate(vaux1)
  return

end subroutine work_function
