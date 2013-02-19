!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE el_phon
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  LOGICAL :: elph, elph_mat, elph_simple
  INTEGER :: elph_nbnd_min, elph_nbnd_max
  INTEGER :: el_ph_ngauss, el_ph_nsigma
  INTEGER :: iunwfcwann, lrwfcr
  INTEGER :: npwq_refolded, ikqg
  INTEGER, allocatable :: wan_index_dyn(:)
  INTEGER, allocatable :: kpq(:), g_kpq(:,:),igqg(:)
  REAL(DP) :: el_ph_sigma
  REAL(DP), allocatable :: xk_gamma(:,:)
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       el_ph_mat(:,:,:,:)    !  nbnd, nbnd, nks, 3*nat
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       el_ph_mat_rec(:,:,:,:)    !  nbnd, nbnd, nksq, npe
  COMPLEX(DP), POINTER :: &
       el_ph_mat_rec_col(:,:,:,:)    !  nbnd, nbnd, nksqtot, npe
  CHARACTER (LEN=256) :: auxdvscf
  LOGICAL, ALLOCATABLE :: comp_elph(:), done_elph(:)
  REAL(DP), ALLOCATABLE :: gamma_disp(:,:,:)
  !
END MODULE el_phon
