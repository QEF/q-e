!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine convert_to_num (npseu, numeric, ndm, mesh, r, lmaxx, &
     lmax, lloc, nnl, aps, alps, vnl)
  !----------------------------------------------------------------------
  !
  !   Convert PPs from analytical to numerical form
  !
  use parameters, only : DP
  implicit none  
  !
  !    first the dummy variables
  !
  integer :: npseu, ndm, mesh (npseu), lmaxx, lmax (npseu), lloc ( &
       npseu), nnl (npseu)
  ! input: the number of pseudopote
  ! input: maximum dimension of the
  ! input: number of mesh points
  ! input: maximum number of ang. m
  ! input: maximum angular momentum
  ! input: the local part of each p
  ! input: the number of gaussians
  real(kind=DP) :: r (0:ndm, npseu), aps (6, 0:3, npseu), alps (3, 0:3, &
       npseu), vnl (0:ndm, 0:lmaxx, npseu)
  ! input: the logaritmic mesh poin
  ! input: the a and b coefficients
  ! input: the alpha coefficients
  ! output: the potential
  logical :: numeric (npseu)  
  ! input: if true the pseudo is nu
  !
  !    here the parameters
  !
  real(kind=DP) :: rcut, e2  
  ! the cut-off radius avoids und
  ! the square of the electron ch
  parameter (rcut = 10.d0, e2 = 2.0d0)  
  !
  !   and the local variables
  !

  integer :: np, l, ir, i  
  ! counter on pseudopotentials
  ! counter on angular momenta
  ! counter on mesh points
  ! counter on alpha functions
  do np = 1, npseu  
     if (.not.numeric (np) ) then  
        !
        ! PP's in analytical form: bring to numerical form!
        !
        do l = 0, lmax (np)  
           do ir = 1, mesh (np)  
              vnl (ir, l, np) = 0.d0  
              if (r (ir, np) .le.rcut) then  
                 do i = 1, nnl (np)  
                    vnl (ir, l, np) = vnl (ir, l, np) + e2 * (aps (i, l, np) &
                  + aps (i + 3, l, np) * r (ir, np) **2) * exp ( - r (ir, np) &
                         **2 * alps (i, l, np) )
                 enddo
              endif
           enddo
        enddo
        !
        ! assume l=lloc as local part and subtract from the other pseudi
        !
        if (lloc (np) .le.lmax (np) ) then  
           do l = 0, lmax (np)  
              if (l.ne.lloc (np) ) then  
                 do ir = 1, mesh (np)  
                    vnl (ir, l, np) = vnl (ir, l, np) - vnl (ir, lloc (np), &
                         np)
                 enddo
              endif
           enddo
           !cc            numeric(np) =.true.
        endif
     endif

  enddo
  return  
end subroutine convert_to_num
