!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE d3_init
!-----------------------------------------------------------------------
  !
  USE ions_base,  ONLY : nat, ntyp => nsp
  USE pwcom
  USE uspp_param, ONLY : vloc_at, zp
  USE atom,       ONLY: numeric, msh, rgrid
  USE phcom
  USE d3com
  USE mp,         ONLY : mp_barrier
  
  IMPLICIT NONE

  INTEGER :: nt, irr, irr1, ipert, imode0, errcode
  REAL (DP) :: work (3)

  COMPLEX (DP), ALLOCATABLE :: drhoscf (:,:)
  COMPLEX (DP), ALLOCATABLE :: drhoscf2 (:,:,:)

  ALLOCATE (drhoscf( nrxx, 3))    

!
!  the fourier trasform of the core charge both for q=0 and q.ne.0
!
  IF (nlcc_any) THEN
!
!  drc is allocated in phq_setup
!
     IF (.NOT.lgamma) THEN

        ALLOCATE (d0rc( ngm, ntyp))    
        work = 0.d0
        CALL set_drhoc (work)
        d0rc (:,:) = drc (:,:)
     ELSE
        d0rc => drc
     ENDIF
!
!  drc is calculated in phq_init
!         call set_drhoc(xq)
  ENDIF
!
! uses the same initialization routines as the phonon program
!
  CALL phq_init
  CALL write_igk
!
!  the fourier components of the local potential at q+G for q=0
!
  IF (.NOT.lgamma) THEN
     vlocg0 (:,:) = 0.d0
     work = 0.d0
     DO nt = 1, ntyp
        CALL setlocq (work, lloc(nt), lmax(nt), numeric(nt), &
             rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r, vloc_at(1,nt), &
             cc(1,nt), alpc(1,nt), nlc(nt), nnl(nt), zp(nt), aps(1,0,nt), &
             alps(1,0,nt), tpiba2, ngm, g, omega, vlocg0(1,nt) )
     ENDDO
  ENDIF
!
! Reads the q=0 variation of the charge --d0rho-- and symmetrizes it
!

  DO irr = 1, nirrg0
     imode0 = 0
     DO irr1 = 1, irr - 1
        imode0 = imode0 + npertg0 (irr1)
     ENDDO
     DO ipert = 1, npertg0 (irr)
        CALL davcio_drho2 (drhoscf(1,ipert), lrdrho, iud0rho, &
                           imode0+ipert, - 1)
     ENDDO
#ifdef __PARA
     CALL psymd0rho (npertg0(irr), irr, drhoscf)
#else
     CALL symd0rho (max_irr_dim, npertg0(irr), irr, drhoscf, s, ftau, nsymg0, &
          irgq, tg0, nat, nr1, nr2, nr3, nrx1, nrx2, nrx3)
#endif
     DO ipert = 1, npertg0 (irr)
        CALL davcio_drho2 (drhoscf(1,ipert), lrdrho, iud0rho, &
                           imode0+ipert, +1)
     ENDDO
  ENDDO
!
! Reads the variation of the charge --drho-- and symmetrizes it
!
  IF (.NOT.lgamma) THEN
     imode0 = 0
     DO irr = 1, nirr
        imode0 = 0
        DO irr1 = 1, irr - 1
           imode0 = imode0 + npert (irr1)
        ENDDO

        ALLOCATE (drhoscf2( nrxx, nspin,npert(irr) ))

        DO ipert = 1, npert (irr)
           CALL davcio_drho (drhoscf2(1,1,ipert), lrdrho, iudrho, &
                              imode0+ipert, -1)
        ENDDO
#ifdef __PARA
        CALL psymdvscf (npert(irr), irr, drhoscf2)
#else
        CALL symdvscf (npert(irr), irr, drhoscf2)
#endif
        DO ipert = 1, npert(irr)
           CALL davcio_drho (drhoscf2(1,1,ipert), lrdrho, iudrho, &
                              imode0+ipert, +1)
        ENDDO
        DEALLOCATE (drhoscf2)

     ENDDO
  ENDIF

  CALL mp_barrier()

  DEALLOCATE(drhoscf)

  RETURN

END SUBROUTINE d3_init
