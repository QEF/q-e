!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! packs the potential and d coefficients into one array
! TEMPORARY: should be merged with vpack
!
subroutine vdpack (ndim, ndimx, nwf, nwfx, nspin, v, d, vd, sflag)
  use kinds, ONLY: DP
  implicit none
  integer :: ndim, ndimx, nwf, nwfx, nspin, is, n, i, ns, ns1
  character(len=4) :: sflag
  real(DP) :: v(ndimx,2), d(nwfx,nwfx,2), vd(ndimx*2+nwfx*nwfx*2)
  select case (sflag)
  case ("PACK")
     i=1
     do is=1, nspin
        do n=1, ndim
           vd(i)=v(n,is)
           i=i+1
        end do
        do ns=1, nwf
           do ns1=1, nwf
              vd(i)=d(ns,ns1,is)
              i=i+1
           end do
        end do
     end do
  case ("UNDO")
     i=1
     do is=1, nspin
        do n=1, ndim
           v(n,is)=vd(i)
           i=i+1
        end do
        do ns=1, nwf
           do ns1=1, nwf
              d(ns,ns1,is)=vd(i)
              i=i+1
           end do
        end do
     end do
  case default
     call errore ('vdpack', ' wrong flag ', 1)
  end select
  return
end subroutine vdpack
