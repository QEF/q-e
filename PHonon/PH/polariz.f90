!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine polariz ( iw, iu )
  !-----------------------------------------------------------------------
  !
  !      calculates the frequency dependent polarizability
  !

  USE io_global,    ONLY : stdout
  USE constants,    ONLY : fpi
  USE cell_base,    ONLY : at, bg, omega
  USE klist,        ONLY : wk, ngk
  USE symme,        ONLY : symmatrix, crys_to_cart
  USE wvfct,        ONLY : npwx
  USE kinds,        ONLY : DP
  USE control_lr,   ONLY : nbnd_occ
  USE units_ph,     ONLY : lrdwf, iudwf, lrebar, iuebar
  USE buffers,      ONLY : get_buffer
  USE freq_ph,      ONLY : polar, done_iu, comp_iu
  USE eqv,          ONLY : dpsi, dvpsi
  USE qpoint,       ONLY : nksq
  USE ph_restart,   ONLY : ph_writefile
  USE cell_base,    ONLY : omega
  USE mp_pools,     ONLY : inter_pool_comm
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(kind=DP) :: iw
  !
  INTEGER, INTENT(IN) :: iu
  !
  ! local variables
  !
  integer :: ibnd, ipol, jpol, nrec, ik, ierr
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  real(kind=DP) :: w, weight, repsilon(3,3)

  complex(kind=DP), EXTERNAL :: zdotc

  call start_clock ('polariz')
  repsilon(:,:) = 0.d0
  do ik = 1, nksq
     weight = wk (ik)
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call get_buffer (dvpsi, lrebar, iuebar, nrec)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call get_buffer(dpsi, lrdwf, iudwf, nrec)
           do ibnd = 1, nbnd_occ (ik)
              !
              !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
              !
              repsilon(ipol,jpol)=repsilon(ipol,jpol)-4.d0*w*REAL( &
                   zdotc ( ngk(ik), dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1) )
           enddo
        enddo
     enddo
  enddo
  call mp_sum ( repsilon, intra_bgrp_comm )
  call mp_sum ( repsilon, inter_pool_comm )
  !
  !      symmetrize
  !
  !       WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
  !       WRITE( stdout,'(10x,"(",3f15.5," )")') ((repsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  call crys_to_cart ( repsilon )
  call symmatrix ( repsilon )
  !
  !    pass to cartesian axis
  !
  !      WRITE( stdout,'(/,10x,"Symmetrized in cartesian axis ",/)')
  !      WRITE( stdout,'(10x,"(",3f15.5," )")') ((repsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     repsilon (ipol, ipol) = repsilon (ipol, ipol) + 1.d0
  enddo
  !
  ! compute the polarization
  !
  do ipol = 1, 3
     do jpol = 1, 3
        if ( repsilon (ipol, jpol) .gt. 1.d-4 ) &
        repsilon (ipol, jpol) = (3.d0*omega/fpi) * ( repsilon (ipol, jpol) - 1.d0 ) / &
                                                  ( repsilon (ipol, jpol) + 2.d0 )
     enddo
  enddo
  !
  !  and print the result
  !
  WRITE( stdout, '(/,10x,"Polarizability in cartesian axis at frequency ",f5.2,/)') iw

  WRITE( stdout, '(10x,"(",3f18.9," )")') ((repsilon(ipol,jpol), ipol=1,3), jpol=1,3)

  polar(:,:,iu)=repsilon(:,:)
  CALL write_polariz(iu)
  done_iu(iu)=.TRUE.
  call ph_writefile('polarization',0,iu,ierr)
  !
  call stop_clock ('polariz')

  return
end subroutine polariz

  SUBROUTINE write_polariz(iu)
!
!  This routine write on output the
!
  USE io_global, ONLY : stdout
  USE constants,    ONLY : BOHR_RADIUS_ANGS
  USE freq_ph, ONLY : fiu, polar
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iu
  INTEGER :: ipol, jpol

  WRITE(stdout,'(2(/),30x," Frequency ",f10.5, "i Ry" )') fiu(iu)
  WRITE(stdout,'(2(/),30x," Cartesian axis " )')
  WRITE(stdout,'(/,5x,"Polarizability (a.u.)^3",20x,"Polarizability (A^3)")')
  WRITE(stdout,'(3f10.2,5x,3f14.4)') ( (polar(ipol,jpol,iu), jpol=1,3), &
                (polar(ipol,jpol,iu)*BOHR_RADIUS_ANGS**3, jpol=1,3), ipol=1,3)
  RETURN
  END SUBROUTINE write_polariz

