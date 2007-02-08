!
! Copyright (C) 2001 PWSCF group
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
  USE ions_base, ONLY : nat
  use pwcom
  USE io_global, ONLY : stdout
  USE kinds, only : DP
  use phcom
  implicit none

  integer :: irr, imode0, npe
  ! input: the irreducible representation
  ! input: starting position of this represe
  ! input: the number of perturbations

  complex(DP) :: dvscfin (nrxx, nspin, npe)
  ! input: the change of V_Hxc

  integer :: ipert, irr1, mode0, mu, is, nu_i, nu_j, nrtot, nspin0
  ! counters
  ! mode0: starting position of the represention
  ! nrtot: the total number of mesh points

  complex(DP) :: dyn1 (3 * nat, 3 * nat)
  ! the dynamical matrix
  complex(DP), allocatable ::  drhous (:,:,:)
  ! the change of the charge
  complex(DP), external :: ZDOTC

  if (.not.okvan) return
  nspin0=nspin
  if (nspin==4.and..not.domag) nspin0=1
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
           do is = 1, nspin0
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
  call poolreduce (18 * nat * nat, dyn1)
  call reduce (18 * nat * nat, dyn1)
#endif
  !       WRITE( stdout,*) 'drhodvus dyn1, dyn'
  !       call tra_write_matrix('drhodvus dyn1',dyn1,u,nat)
  !       call tra_write_matrix('drhodvus dyn',dyn,u,nat)
  !       call stop_ph(.true.)
  dyn (:,:) = dyn (:,:) + dyn1 (:,:)
  call stop_clock ('drhodvus')
  return

end subroutine drhodvus
