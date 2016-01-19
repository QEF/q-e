!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dielec()
  !-----------------------------------------------------------------------
  !
  !      calculates the dielectric tensor
  !
  USE kinds, only : DP
  USE io_global,  ONLY : stdout
  USE constants, ONLY: fpi
  USE cell_base, ONLY: omega
  USE klist, ONLY: wk
  USE symme, ONLY: symmatrix, crys_to_cart
  USE wvfct, ONLY: npwx
  USE buffers, ONLY : get_buffer
  USE noncollin_module, ONLY : npol
  USE efield_mod, ONLY : epsilon
  USE units_ph, ONLY : lrdwf, iudwf, lrebar, iuebar
  USE eqv, ONLY : dpsi, dvpsi
  USE qpoint, ONLY : nksq
  USE ph_restart, ONLY : ph_writefile
  USE control_lr, ONLY : nbnd_occ
  USE control_ph, ONLY : done_epsil, epsil
  USE mp_pools,   ONLY : inter_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum

  implicit none

  integer :: ibnd, ipol, jpol, nrec, ik, ierr
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  real(DP) :: w, weight

  complex(DP), external :: zdotc

  IF (.NOT.epsil.OR.done_epsil) RETURN

  call start_clock ('dielec')
  epsilon(:,:) = 0.d0
  do ik = 1, nksq
     weight = wk (ik)
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call get_buffer(dvpsi, lrebar, iuebar, nrec)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call get_buffer (dpsi, lrdwf, iudwf, nrec)
           do ibnd = 1, nbnd_occ (ik)
              !
              !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
              !
              epsilon(ipol,jpol)=epsilon(ipol,jpol)-4.d0*w* DBLE( &
                   zdotc(npwx*npol, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1))
           enddo
        enddo
     enddo
  enddo
  call mp_sum ( epsilon, intra_bgrp_comm )
  call mp_sum ( epsilon, inter_pool_comm )
  !
  !      symmetrize
  !
  !       WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
  !       WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)

  call crys_to_cart (epsilon)
  call symmatrix ( epsilon )
  !
  !    pass to cartesian axis
  !
  !      WRITE( stdout,'(/,10x,"Symmetrized in cartesian axis ",/)')
  !      WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     epsilon (ipol, ipol) = epsilon (ipol, ipol) + 1.d0
  enddo
  !
  !  and print the result
  !
  done_epsil=.TRUE.
  CALL summarize_epsilon()
  CALL ph_writefile('tensors',0,0,ierr)

  call stop_clock ('dielec')

  return
end subroutine dielec
