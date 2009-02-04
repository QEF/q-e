!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine drhodvus (irr, imode0, dvscfin, npe)
  !-----------------------------------------------------------------------
  !
  !    This subroutine calculates the term of the dynamical matrix
  !    which comes from the interaction of the change of the self consiste
  !    potential with the static change of the charge density.
  !    This term is non zero only if the charge is augmented.
  !
  !
#include "f_defs.h"
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp
  USE gvect,     ONLY : nrxx, nr1, nr2, nr3
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat
  USE spin_orb,  ONLY : domag
  USE uspp,      ONLY : okvan
  USE io_global, ONLY : stdout
  USE uspp_param, ONLY : upf, nh
  USE paw_variables, ONLY : okpaw
  USE noncollin_module, ONLY : nspin_mag

  USE modes,     ONLY : npert, npertx, nirr
  USE dynmat,    ONLY : dyn, dyn_rec
  USE phus,      ONLY : becsumort, int3_paw
  USE units_ph,  ONLY : iudrhous, lrdrhous

  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  implicit none

  integer :: irr, imode0, npe
  ! input: the irreducible representation
  ! input: starting position of this represe
  ! input: the number of perturbations

  complex(DP) :: dvscfin (nrxx, nspin, npe)
  ! input: the change of V_Hxc

  integer :: ipert, irr1, mode0, mu, is, nu_i, nu_j, nrtot, &
             ih, jh, ijh, na, nb, nt
  ! counters
  ! mode0: starting position of the represention
  ! nrtot: the total number of mesh points

  complex(DP) :: dyn1 (3 * nat, 3 * nat)
  ! the dynamical matrix
  complex(DP), allocatable ::  drhous (:,:,:)
  ! the change of the charge
  complex(DP), external :: ZDOTC

  if (.not.okvan) then
     dyn_rec=(0.0_DP,0.0_DP) 
     return
  endif
  call start_clock ('drhodvus')
  allocate (drhous ( nrxx , nspin, npertx))    
  dyn1 (:,:) = (0.d0, 0.d0)
  nrtot = nr1 * nr2 * nr3
  mode0 = 0
  do irr1 = 1, nirr
     do ipert = 1, npert (irr1)
        nu_j = mode0 + ipert
        call davcio (drhous (1, 1, ipert), lrdrhous, iudrhous, nu_j, - 1)
     enddo
     do ipert = 1, npert (irr1)
        nu_j = mode0 + ipert
        do mu = 1, npert (irr)
           nu_i = imode0 + mu
           do is = 1, nspin_mag
              dyn1 (nu_i, nu_j) = dyn1 (nu_i, nu_j) + &
                   ZDOTC (nrxx, dvscfin (1,is,mu), 1, drhous (1,is,ipert), 1) &
                   * omega / DBLE (nrtot)
           enddo
        enddo
     enddo
     mode0 = mode0 + npert (irr1)
  enddo
  deallocate (drhous)
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  !
  call mp_sum ( dyn1, inter_pool_comm )
  call mp_sum ( dyn1, intra_pool_comm )
#endif
!
!  PAW contribution: this part of the dynamical matrix is present only 
!  with PAW. PAW and US dynamical matrices differ only at this point.
!
  IF (okpaw) THEN
     mode0 = 0
     do irr1 = 1, nirr
        do ipert = 1, npert (irr1)
           nu_j = mode0 + ipert
           do mu = 1, npert (irr)
              nu_i = imode0 + mu
              do nt=1,ntyp
                 if (upf(nt)%tpawp) then
                    ijh=0
                    do ih=1,nh(nt)
                       do jh=ih,nh(nt)
                          ijh=ijh+1
                          do na=1,nat
                             if (ityp(na)==nt) then
                                do is = 1, nspin_mag
                                   dyn1(nu_i,nu_j)=dyn1(nu_i,nu_j)+ &
                                       CONJG(int3_paw(ih,jh,mu,na,is))* &
                                       becsumort(ijh,na,is,nu_j)
                                enddo
                             endif
                          enddo
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
        mode0 = mode0 + npert (irr1)
      enddo
   endif
  !       WRITE( stdout,*) 'drhodvus dyn1, dyn'
  !       call tra_write_matrix('drhodvus dyn1',dyn1,u,nat)
  !       call tra_write_matrix('drhodvus dyn',dyn,u,nat)
  !       call stop_ph(.true.)
  dyn (:,:) = dyn (:,:) + dyn1 (:,:)
  dyn_rec(:,:) = dyn1(:,:)
  call stop_clock ('drhodvus')
  return

end subroutine drhodvus
