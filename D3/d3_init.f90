!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include"machine.h"
!
!-----------------------------------------------------------------------
subroutine d3_init
!-----------------------------------------------------------------------

  use pwcom
  use atom, only: numeric, mesh, msh, rab, r
  use phcom
  use d3com
  use para
  USE mp,    ONLY : mp_barrier
  
  implicit none

  integer :: nt, irr, irr1, ipert, imode0, errcode
! counter on atom types
#ifdef DEBUG
  integer :: i, j, k
#endif

  real (kind = dp) :: work (3)                 ! working area

  complex (kind = dp), allocatable :: drhoscf (:,:)
  complex (kind = dp), allocatable :: drhoscf2 (:,:,:)

  allocate (drhoscf( nrxx, 3))    
! allocate (drhoscf2( nrxx, 1, 1))
!
!  the fourier trasform of the core charge both for q=0 and q.ne.0
!
  if (nlcc_any) then
!
!  drc is allocated in phq_setup
!
     if (.not.lgamma) then

        allocate (d0rc( ngm, ntyp))    
        work = 0.d0
        call set_drhoc (work)
        d0rc (:,:) = drc (:,:)
     else
        d0rc => drc
     endif
!
!  drc is calculated in phq_init
!         call set_drhoc(xq)
  endif
!
! uses the same initialization routines as the phonon program
!
  call phq_init
  call write_igk
!
!  the fourier components of the local potential at q+G for q=0
!
  if (.not.lgamma) then
     vlocg0 (:,:) = 0.d0
     work = 0.d0
     do nt = 1, ntyp
        call setlocq (work, lloc(nt), lmax(nt), numeric(nt), &
             mesh(nt), msh(nt), rab(1,nt), r(1,nt), vloc_at(1,nt), &
             cc(1,nt), alpc(1,nt), nlc(nt), nnl(nt), zp(nt), aps(1,0,nt), &
             alps(1,0,nt), tpiba2, ngm, g, omega, vlocg0(1,nt) )
     enddo
  endif
!
! Reads the q=0 variation of the charge --d0rho-- and symmetrizes it
!
#ifdef __PARA
!  if (mypool /= 1) goto 100
#endif
  do irr = 1, nirrg0
     imode0 = 0
     do irr1 = 1, irr - 1
        imode0 = imode0 + npertg0 (irr1)
     enddo
     do ipert = 1, npertg0 (irr)
        call davcio_drho2 (drhoscf(1,ipert), lrdrho, iud0rho, &
                           imode0+ipert, - 1)
     enddo
#ifdef __PARA
     call psymd0rho (npertg0(irr), irr, drhoscf)
#else
     call symd0rho (max_irr_dim, npertg0(irr), irr, drhoscf, s, ftau, nsymg0, &
          irgq, tg0, nat, nr1, nr2, nr3, nrx1, nrx2, nrx3)
#endif
     do ipert = 1, npertg0 (irr)
        call davcio_drho2 (drhoscf(1,ipert), lrdrho, iud0rho, &
                           imode0+ipert, +1)
     enddo
  enddo
!
! Reads the variation of the charge --drho-- and symmetrizes it
!
  if (.not.lgamma) then
     imode0 = 0
     do irr = 1, nirr
        imode0 = 0
        do irr1 = 1, irr - 1
           imode0 = imode0 + npert (irr1)
        enddo

        allocate (drhoscf2( nrxx, nspin,npert(irr) ))

        do ipert = 1, npert (irr)
           call davcio_drho (drhoscf2(1,1,ipert), lrdrho, iudrho, &
                              imode0+ipert, -1)
        enddo
#ifdef __PARA
        call psymdvscf (npert(irr), irr, drhoscf2)
#else
        call symdvscf (npert(irr), irr, drhoscf2)
#endif
        do ipert = 1, npert(irr)
           call davcio_drho (drhoscf2(1,1,ipert), lrdrho, iudrho, &
                              imode0+ipert, +1)
        enddo
        deallocate (drhoscf2)

!        imode0 = imode0 + npert(irr)
     enddo
  endif
#ifdef __PARA
100 continue
  call mp_barrier()
#endif

  deallocate(drhoscf)
!  deallocate(drhoscf2)
  return

end subroutine d3_init
