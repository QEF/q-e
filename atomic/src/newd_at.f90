!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
subroutine newd_at ( )
  !-------------------------------------------------------------------------
  !
  !     this routine computes the new D coefficients
  !
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use ld1inc, only : ddd, bmat, nbeta, nspin, lls, jjs, ikk, qvan, vpstot, &
                     grid, pseudotype, lpaw, which_augfun, qvanl
  implicit none

  integer :: &
       ib,jb,n,is,nst

  real(DP) :: &
       int_0_inf_dr, &   ! the integral function
       gi(ndmx)          ! the gi function

  !
  !    screening the D coefficients
  !
  if (pseudotype == 3) then
     do ib=1,nbeta
        do jb=1,ib
           if (lls(ib).eq.lls(jb).and.abs(jjs(ib)-jjs(jb)).lt.1.0e-7_dp) then
              nst=(lls(ib)+1)*2
              do is=1,nspin
                 IF (which_augfun=='PSQ') then
                    do n=1,ikk(ib)
                       gi(n)=qvanl(n,ib,jb,0)*vpstot(n,is)
                    enddo
                 ELSE
                    do n=1,ikk(ib)
                       gi(n)=qvan(n,ib,jb)*vpstot(n,is)
                    enddo
                 ENDIF
                 ddd(ib,jb,is)= bmat(ib,jb) &
                      + int_0_inf_dr(gi,grid,ikk(ib),nst)
                 ddd(jb,ib,is)=ddd(ib,jb,is)
              enddo
           endif
        enddo
     enddo
  else if (pseudotype == 2) then
     !
     !    non-US separable PP case: just copy unscreened D coeffs
     !
     do is=1,nspin
        ddd(:,:,is)= bmat(:,:)
     enddo
  endif

  return
end subroutine newd_at
