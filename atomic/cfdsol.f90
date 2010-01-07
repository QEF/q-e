!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!   This routine is a f90 translation of the routine present in 
!   Vanderbilt code.
!
!-----------------------------------------------------------------------
!
subroutine cfdsol(zz,yy,jj1,jj2,idim1)
  !
  !-----------------------------------------------------------------------
  !
  !     routine for solving coupled first order differential equations
  !
  !      d yy(x,1)
  !      ---------   =  zz(x,1,1) * yy(x,1) + zz(x,1,2) * yy(2,1)
  !         dx
  !
  !      d yy(x,2)
  !      ---------   =  zz(x,2,1) * yy(x,1) + zz(x,2,2) * yy(2,1)
  !         dx
  !     
  !
  !     using fifth order predictor corrector algorithm
  !
  !     routine integrates from jj1 to jj2 and can cope with both cases
  !     jj1 < jj2 and jj1 > jj2.  first five starting values of yy must 
  !     be provided by the calling program.
  !
  !-----------------------------------------------------------------------
  !
  !
  use kinds, only : DP
  implicit none
  integer :: idim1, jj1, jj2, ip
  real(DP):: zz(idim1,2,2),yy(idim1,2)
  real(DP):: fa(0:5),fb(0:5),abp(1:5),amc(0:4)
  real(DP):: arp, brp

  integer :: isgn, i, j
  !
  !
  !-----------------------------------------------------------------------
  !
  !                   i n i t i a l i s a t i o n
  !
  !     decide whether integrating from:
  !            left to right ---> isgn = + 1
  !        or  right to left ---> isgn = - 1
  !
  isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
  !
  !     run some test just to be conservative
  if ( isgn .eq. + 1 ) then
     if ( jj1 .le. 5 .or. jj2 .gt. idim1 ) then
        write(6,10) isgn,jj1,jj2,idim1
        call errore('cfdsol','stopping jj1 to small or jj2 to large',1)
     endif
  elseif ( isgn .eq. - 1 ) then
     if ( jj1 .ge. ( idim1 - 4 ) .or. jj2 .lt. 1 ) then
        write(6,10) isgn,jj1,jj2,idim1
        call errore('cfdsol','stopping jj1 to large or jj2 too small',1)
     endif
  else
     write(6,10) isgn,jj1,jj2,idim1
  endif

10 format(' ***error in subroutine cfdsol',/,   &
       &  ' isgn =',i2,' jj1 =',i5,' jj2 =',i5,' idim1 =',i5, &
       &  ' are not allowed')
  !
  !     integration coefficients
  !
    abp(1) = 1901.0_dp / 720.0_dp
    abp(2) = -1387.0_dp / 360.0_dp
    abp(3) = 109.0_dp / 30.0_dp
    abp(4) = -637.0_dp / 360.0_dp
    abp(5) = 251.0_dp / 720.0_dp
    amc(0) = 251.0_dp / 720.0_dp
    amc(1) = 323.0_dp / 360.0_dp
    amc(2) = -11.0_dp / 30.0_dp
    amc(3) = 53.0_dp / 360.0_dp
    amc(4) = -19.0_dp / 720.0_dp
    !
    !     set up the arrays of derivatives
    do j = 1,5
       ip = jj1 - isgn * j
       fa(j) = zz(ip,1,1) * yy(ip,1) + zz(ip,1,2) * yy(ip,2)
       fb(j) = zz(ip,2,1) * yy(ip,1) + zz(ip,2,2) * yy(ip,2)
    enddo
    !
    !-----------------------------------------------------------------------
    !
    !                i n t e g r a t i o n  l o o p
    !
    do j = jj1,jj2,isgn
       !
       !       predictor (adams-bashforth)
       !
       arp = yy(j-isgn,1)
       brp = yy(j-isgn,2)
       do i = 1,5
          arp = arp + DBLE(isgn) * abp(i) * fa(i)
          brp = brp + DBLE(isgn) * abp(i) * fb(i)
       enddo

       fa(0) = zz(j,1,1) * arp + zz(j,1,2) * brp
       fb(0) = zz(j,2,1) * arp + zz(j,2,2) * brp
       !
       !       corrector (adams-moulton)
       !
       yy(j,1) = yy(j-isgn,1)
       yy(j,2) = yy(j-isgn,2)
       do i = 0,4,1
          yy(j,1) = yy(j,1) + DBLE(isgn) * amc(i) * fa(i)
          yy(j,2) = yy(j,2) + DBLE(isgn) * amc(i) * fb(i)
       enddo
       !
       !       book keeping
       !
       do i = 5,2,-1
          fa(i) = fa(i-1)
          fb(i) = fb(i-1)
       enddo
       fa(1) = zz(j,1,1) * yy(j,1) + zz(j,1,2) * yy(j,2)
       fb(1) = zz(j,2,1) * yy(j,1) + zz(j,2,2) * yy(j,2)
    enddo

    return
end subroutine cfdsol
