!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine polariz ( iw )
  !-----------------------------------------------------------------------
  !
  !      calculates the frequency dependent polarizability
  !

  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : iunigk
  USE constants,    ONLY : fpi
  USE cell_base,    ONLY : at, bg, omega
  USE klist,        ONLY : wk
  USE symme,        ONLY : symmatrix, crys_to_cart
  USE wvfct,        ONLY : npw, npwx, igk
  USE kinds,        ONLY : DP
  USE efield_mod,   ONLY : epsilon
  USE control_ph,   ONLY : nbnd_occ
  USE units_ph,     ONLY : lrdwf, iudwf, lrebar, iuebar
  USE eqv,          ONLY : dpsi, dvpsi
  USE qpoint,       ONLY : nksq
  USE cell_base,    ONLY : omega
  USE mp_global,    ONLY : inter_pool_comm, intra_pool_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  real(kind=DP) :: iw
  !
  ! local variables
  !
  integer :: ibnd, ipol, jpol, nrec, ik
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  real(kind=DP) :: w, weight

  complex(kind=DP), EXTERNAL :: zdotc

  call start_clock ('polariz')
  epsilon(:,:) = 0.d0
  if (nksq > 1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     weight = wk (ik)
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call davcio (dvpsi, lrebar, iuebar, nrec, - 1)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
           do ibnd = 1, nbnd_occ (ik)
              !
              !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
              !
              epsilon(ipol,jpol)=epsilon(ipol,jpol)-4.d0*w*REAL( &
                   zdotc (npw, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1) )
           enddo
        enddo
     enddo
  enddo
#ifdef __PARA
  call mp_sum ( epsilon, intra_pool_comm )
  call mp_sum ( epsilon, inter_pool_comm )
#endif
  !
  !      symmetrize
  !
  !       WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
  !       WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  call crys_to_cart ( epsilon )
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
  ! compute the polarization
  !
  do ipol = 1, 3
     do jpol = 1, 3
        if ( epsilon (ipol, jpol) .gt. 1.d-4 ) &
        epsilon (ipol, jpol) = (3.d0*omega/fpi) * ( epsilon (ipol, jpol) - 1.d0 ) / &
                                                  ( epsilon (ipol, jpol) + 2.d0 )
     enddo
  enddo
  !
  !  and print the result
  !
  WRITE( stdout, '(/,10x,"Polarizability in cartesian axis at frequency ",f5.2,/)') iw

  WRITE( stdout, '(10x,"(",3f18.9," )")') ((epsilon(ipol,jpol), ipol=1,3), jpol=1,3)

  call stop_clock ('polariz')

  return
end subroutine polariz
