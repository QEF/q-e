!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module pseudodata
  USE kinds, ONLY: DP
  USE parameters , ONLY: nchix
!  c
!  c The redundant data of KKRJ3 files to be written in UPF header
!  c

  integer          :: &
        nns(nchix),&       ! The pseudo n of a state 
        pseudotype,&       ! the type of pseudopotential
        nwfs,&             ! Number of pseudo wf
         ikk,&             ! the kkbeta for each beta
         ntwfc             ! Number of occupied states
  character(len=2) :: &
        els(nchix)         ! The name of the state(sp notation) 
  real(kind=DP)    :: &
         rcut(nchix),&     ! Pseudization cutoff radius
         rcloc,&           ! Local potential cutoff radius
         rcutus(nchix),&   ! US Pseudization cutoff radii
         etotps,&          ! total energy of the pseudoatom
         rmax              ! Max R of the Mesh
  logical          :: &
         rel               ! if true the atomic calculation is relativistic

  integer          :: &
        ikk2(nchix)        ! Global ikk subst

end module pseudodata
