!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dielec(do_zstar)
  !-----------------------------------------------------------------------
  !
  !      calculates the dielectric tensor and effective charges
  !
#include "f_defs.h"
  USE ions_base, ONLY : nat, zv, ityp
  use pwcom
  use cgcom
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

  implicit none
  logical :: do_zstar
  !
  integer :: ibnd,ipol,jpol,na,nu,kpoint
  character(len=7) :: filbar, fildwf
  real(DP) ::  w, weight
  real(DP), allocatable ::  work(:,:)
  complex(DP), allocatable :: dpsi2(:,:), dpsi3(:,:)
  logical :: done
  !
  call start_clock('dielec')
  !
  allocate (dpsi2( npwx, nbnd))    
  allocate (dpsi3( npwx, nbnd))    
  allocate (work( nbnd, 3))    
  !
  epsilon0(:,:) = 0.d0
  if (do_zstar) zstar (:,:,:) = 0.d0
  !  do kpoint=1,nks
  kpoint=1
  weight = wk(kpoint)
  w = fpi/omega * weight
  !
  !** calculate Effective Charges (<DeltaV*psi(ion)|DeltaPsi(E)>)
  !
  ! read DeltaPsi(E)
  ! pol. 1
  ipol=1
  iudwf=10+ipol
  write(fildwf,'("fildwx",i1)') ipol
  call  seqopn (iudwf,fildwf,'unformatted',done)
  read (iudwf) dpsi
  close(unit=iudwf)
  ! pol. 2
  ipol=2
  iudwf=10+ipol
  write(fildwf,'("fildwx",i1)') ipol
  call  seqopn (iudwf,fildwf,'unformatted',done)
  read (iudwf) dpsi2
  close(unit=iudwf)
  ! pol. 3
  ipol=3
  iudwf=10+ipol
  write(fildwf,'("fildwx",i1)') ipol
  call  seqopn (iudwf,fildwf,'unformatted',done)
  read (iudwf) dpsi3
  close(unit=iudwf)
  !
  if (.not.do_zstar) go to 10
  !
  do nu = 1,nmodes
     na  = (nu-1)/3+1
     if (has_equivalent(na).eq.0) then
        !     DeltaV*psi(ion) for mode nu is recalculated
        call dvpsi_kb(kpoint,nu)
        !
        jpol= mod(nu-1,3)+1
        ! work is the real part of <DeltaV*Psi(ion)|DeltaPsi(E)>
        call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi ,npwx,work(1,1))
        call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi2,npwx,work(1,2))
        call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi3,npwx,work(1,3))
        do ipol = 1,3
           do ibnd = 1,nbnd
              zstar(ipol,jpol,na) = zstar(ipol,jpol,na) + 2.0d0*weight*work(ibnd,ipol)
           end do
        end do
     end if
  end do
10 continue
  !** calculate Dielectric Tensor (<DeltaV*psi(E)\DeltaPsi(E)>)
  !
  do jpol=1,3
     ! read DeltaV*Psi(elec) for polarization jpol
     iubar=jpol
     write(filbar,'("filbar",i1)') iubar
     call  seqopn (iubar,filbar,'unformatted',done)
     read (iubar) dvpsi
     close(iubar)
     ! now work is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
     call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi ,npwx,work(1,1))
     call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi2,npwx,work(1,2))
     call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi3,npwx,work(1,3))
     do ipol = 1,3
        do ibnd = 1,nbnd
           epsilon0(ipol,jpol) = epsilon0(ipol,jpol) + 4.0d0*w*work(ibnd,ipol)
        end do
     end do
  end do
  !     end do
#ifdef __PARA
  if (do_zstar) call mp_sum( zstar, intra_pool_comm )
  call mp_sum( epsilon0, intra_pool_comm )
#endif
  deallocate(work)
  deallocate(dpsi3)
  deallocate(dpsi2)
  !
  ! add the diagonal part
  !
  do ipol=1,3
     epsilon0(ipol,ipol) = epsilon0(ipol,ipol) + 1.0d0
     if (do_zstar) then
        do na=1,nat
           zstar(ipol,ipol,na) = zstar(ipol,ipol,na) + zv(ityp(na))
        end do
     end if
  end do
  !
  call stop_clock('dielec')
  !
  return
end subroutine dielec
