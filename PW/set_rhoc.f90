!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_rhoc
  !-----------------------------------------------------------------------
  !
  !    This routine compute the core charge on the real space 3D mesh
  !
  !
#include "machine.h"
  use pwcom
  implicit none
  !
  !     One parameter
  !
  real(kind=DP) :: eps
  ! a small number

  parameter (eps = 1.d-10)

  complex(kind=DP) , allocatable :: aux (:)
  ! used for the fft of the core ch

  real(kind=DP) , allocatable ::  rhocg(:)
  real(kind=DP) ::  rhoima, rhoneg, rhorea
  real(kind=DP) , allocatable ::  dum(:,:)
  real(kind=DP) ::  vtxcc
  ! the radial fourier trasform
  ! used to check the core charge
  ! used to check the core charge
  ! the real core charge
  ! dummy array containing rho=0
  ! dummy xc energy term

  integer :: ir, nt, ng
  ! counter on mesh points
  ! counter on atomic types
  ! counter on g vectors

  logical :: no_core_only
  ! if .f. subtract etxcc from etot

  etxcc = 0.d0
  do nt = 1, ntyp
     if (nlcc (nt) ) goto 10
  enddo
  call setv (nrxx, 0.d0, rho_core, 1)
  return

10 continue
  allocate (aux( nrxx))    
  allocate (rhocg(  ngl))    
  call setv (2 * nrxx, 0.d0, aux, 1)
  !
  !    the sum is on atom types
  !
  do nt = 1, ntyp
     if (nlcc (nt) ) then
        !
        !     drhoc compute the radial fourier transform for each shell of g vec
        !
        call drhoc (ngl, gl, omega, tpiba2, numeric (nt), a_nlcc (nt), &
             b_nlcc (nt), alpha_nlcc (nt), msh (nt), r (1, nt), rab (1, nt), &
             rho_atc (1, nt), rhocg)
        !
        !     multiply by the structure factor and sum
        !
        do ng = 1, ngm
           aux (nl (ng) ) = aux (nl (ng) ) + strf (ng, nt) * rhocg (igtongl (ng) )
        enddo
     endif
  enddo
  !
  !   the core charge in real space
  !
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  !    test on the charge and computation of the core energy
  !
  rhoneg = 0.d0
  rhoima = 0.d0
  do ir = 1, nrxx
     rhoneg = rhoneg + min (0.d0, DREAL (aux (ir) ) )
     rhoima = rhoima + abs (DIMAG (aux (ir) ) )
     rhorea = max (DREAL (aux (ir) ), eps)
     rho_core(ir) = DREAL (aux(ir))
     !
     ! NOTE: Core charge is computed in reciprocal space and brought to real
     ! space by FFT. For non smooth core charges (or insufficient cut-off)
     ! this may result in negative values in some grid points.
     ! Up to October 1999 the core charge was forced to be positive definite.
     ! This induces an error in the force, and probably stress, calculation i
     ! the number of grid points where the core charge would be otherwise neg
     ! is large. The error disappears for sufficiently high cut-off, but may
     ! rather large and it is better to leave the core charge as it is.
     ! If you insist to have it positive definite (with the possible problems
     ! mentioned above) uncomment the following line.  SdG, Oct 15 1999
     !
     !         rho_core(ir) = rhorea
     !
  enddo
  rhoneg = rhoneg / (nr1 * nr2 * nr3)
  rhoima = rhoima / (nr1 * nr2 * nr3)
#ifdef PARA
  call reduce (1, rhoneg)
  call reduce (1, rhoima)
#endif
  if (rhoneg.lt. - 1.0d-6.or.rhoima.gt.1.0d-6) &
       write (6, '("  warning: negative or imaginary core charge ",2f12.6)')&
       rhoneg, rhoima
  !
  no_core_only = .true.
  if (no_core_only) then
     etxcc = 0.d0
  else
     !
     ! calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
     ! This term is present only for compatibility with previous versions
     !
     allocate (dum(nrxx , nspin))    
     call setv (nspin * nrxx, 0.d0, dum, 1)
     call v_xc (dum, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
          nrxx, nl, ngm, g, nspin, alat, omega, etxcc, vtxcc, aux)
     deallocate(dum)
     write (6, 9000) etxcc
     write (6,  * ) 'BEWARE it will be subtracted from total energy !'

  endif
  deallocate (rhocg)
  deallocate (aux)
  !
  return

9000 format (5x,'core-only xc energy         = ',f15.8,' ryd')

end subroutine set_rhoc

