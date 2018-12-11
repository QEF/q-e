!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE sym_dns_wrapper (ldim, dns_cart, dns_pattern)
  !-------------------------------------------------------------------
  !
  ! This routine symmetrizes dns_cart. This is done in three steps.
  !
  ! Written by I. Timrov using the code by S. de Gironcoli 
  ! and A. Floris (01.10.2018)
  !  
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat
  USE modes,         ONLY : u, nmodes, nirr, npert
  USE lsda_mod,      ONLY : nspin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim
  COMPLEX(DP), INTENT(INOUT) ::  dns_cart(ldim,ldim,nspin,nat,3,nat)
  COMPLEX(DP), INTENT(OUT) :: dns_pattern(ldim,ldim,nspin,nat,3*nat)
  ! in/out : dns_cart is the dns matrix in the cartesian coordinates;
  !          on the input it is unsymmetrized, on the output it is symmetrized
  ! out    : dns_pattern is the symmetrized dns matrix in the pattern basis
  !
  ! Local variables
  !
  INTEGER :: imode, imode0, na, icart, na_icart, irr, npe
  COMPLEX(DP), ALLOCATABLE :: dns_aux(:,:,:,:,:)
  !
  dns_pattern = (0.d0, 0.d0)
  !
  ! 1 -  transform in the pattern basis
  !
  DO imode = 1, nmodes
     DO na = 1, nat
        DO icart = 1, 3
           na_icart = 3*(na-1) + icart
           dns_pattern(:,:,:,:,imode) = dns_pattern(:,:,:,:,imode) + &
                            dns_cart(:,:,:,:,icart,na) * u(na_icart,imode)
        ENDDO
     ENDDO
  ENDDO
  !
  ! 2 - symmetrize in the pattern basis
  !
  imode0 = 1
  DO irr = 1, nirr
     npe = npert(irr)
     ! allocate
     ALLOCATE (dns_aux(ldim,ldim,nspin,nat,npe))
     ! pack
     dns_aux(:,:,:,:,1:npe) = dns_pattern(:,:,:,:,imode0:imode0-1+npe)
     ! symmetrize
     CALL sym_dns (ldim, npe, irr, dns_aux)
     ! unpack
     dns_pattern(:,:,:,:,imode0:imode0-1+npe) = dns_aux(:,:,:,:,1:npe)
     ! deallocate
     DEALLOCATE (dns_aux)
     ! advance the counter
     imode0 = imode0 + npe
  ENDDO
  !
  ! 3 - back to the cartesian basis
  ! 
  dns_cart = (0.d0, 0.d0)
  !
  DO imode = 1, nmodes
     DO na = 1, nat
        DO icart = 1, 3
           na_icart = 3*(na-1) + icart
           dns_cart(:,:,:,:,icart,na) = dns_cart(:,:,:,:,icart,na) + &
               dns_pattern(:,:,:,:,imode) * CONJG(u(na_icart,imode))
        ENDDO
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE sym_dns_wrapper
!----------------------------------------------------------------------------
