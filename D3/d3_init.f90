!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3_init  
!-----------------------------------------------------------------------
#include"machine.h"
  use pwcom
  use phcom
  use d3com
  use allocate
#ifdef PARA
  use para
#endif
  implicit none
#ifdef PARA
  include 'mpif.h'  
#endif
  integer :: nt, irr, irr1, ipert, imode0, errcode  
! counter on atom types
  real (8) :: work (3)                 ! working area

  complex (8), pointer :: drhoscf (:,:)  

  call mallocate(drhoscf, nrxx, 3)  
!
!  the fourier trasform of the core charge both for q=0 and q.ne.0
!
  if (nlcc_any) then  
!
!  drc is allocated in phq_setup
!
     if (.not.lgamma) then  

        call mallocate(d0rc, ngm, ntyp)  
        call setv (3, 0.d0, work, 1)  
        call set_drhoc (work)  
        call ZCOPY (ngm * ntyp, drc, 1, d0rc, 1)  
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
     call setv (ngm * ntyp, 0.d0, vlocg0, 1)  
     call setv (3, 0.d0, work, 1)  
     do nt = 1, ntyp  
        call setlocq (work, lloc(nt), lmax(nt), numeric(nt), &
             mesh(nt), msh(nt), rab(1,nt), r(1,nt), vnl(1,lloc(nt),nt), &
             cc(1,nt), alpc(1,nt), nlc(nt), nnl(nt), zp(nt), aps(1,0,nt), &
             alps(1,0,nt), tpiba2, ngm, g, omega, vlocg0(1,nt) )
     enddo
  endif
!
! Reads the q=0 variation of the charge --d0rho-- and symmetrizes it
!
#ifdef PARA
!  if (mypool.ne.1) goto 100  
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
#ifdef PARA
     call psymd0rho (npertg0(irr), irr, drhoscf)  
#else
     call symd0rho (npertg0(irr), irr, drhoscf, s, ftau, nsymg0, &
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
     do irr = 1, nirr  
        imode0 = 0  
        do irr1 = 1, irr - 1  
           imode0 = imode0 + npert (irr1)  
        enddo
        do ipert = 1, npert (irr)  
           call davcio_drho2 (drhoscf(1,ipert), lrdrho, iudrho, &
                              imode0+ipert, -1)
        enddo
#ifdef PARA
        call psymdvscf (npert(irr), irr, drhoscf)  
#else
        call symdvscf (npert(irr), irr, drhoscf)  
#endif
        do ipert = 1, npert(irr)  
           call davcio_drho2 (drhoscf(1,ipert), lrdrho, iudrho, &
                              imode0+ipert, +1)
        enddo
     enddo
  endif
#ifdef PARA
100 continue  
  call MPI_barrier (MPI_COMM_WORLD, errcode)  

  call error ('d3_init', 'at barrier', errcode)  
#endif

  call mfree(drhoscf)  
  return  

end subroutine d3_init
