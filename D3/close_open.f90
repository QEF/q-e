!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine close_open (isw)
  !-----------------------------------------------------------------------
  !
  ! Close and open some units. It is useful in case of interrupted run
  !
  !
#include"machine.h"
  use pwcom, only: filpun, degauss
  use phcom, only: iudwf, lrdwf, lgamma
  use d3com
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: isw
  character (len=20) :: filint * 42
  ! the name of the file
  logical :: exst
  ! logical variable to check file existenc

  if (len_trim(filpun).eq.0) call errore ('recv', 'wrong filpun name', 1)
  if (isw.eq.3) then
     !
     ! This is to be used after gen_dwf(3)
     !
#ifdef __PARA
     if (me.ne.1.or.mypool.ne.1) goto 210
#endif
     if (degauss.ne.0.d0) then
        close (unit = iuef, status = 'keep')
        filint = trim(filpun) //'.efs'
        call seqopn (iuef, filint, 'unformatted', exst)
     endif
#ifdef __PARA

210  continue
#endif
     close (unit = iupd0vp, status = 'keep')
     filint = trim(filpun) //'.p0p'
     if (lgamma) filint = trim(filpun) //'.pdp'

     call diropn (iupd0vp, filint, lrpdqvp, exst)
     close (unit = iudwf, status = 'keep')
     filint = trim(filpun) //'.dwf'

     call diropn (iudwf, filint, lrdwf, exst)
  elseif (isw.eq.1) then
     !
     ! This is to be used after gen_dwf(1)
     !

     if (lgamma) call errore (' close_open ', ' isw=1 ; lgamma', 1)
     close (unit = iupdqvp, status = 'keep')
     filint = trim(filpun) //'.pdp'

     call diropn (iupdqvp, filint, lrpdqvp, exst)
     close (unit = iudqwf, status = 'keep')
     filint = trim(filpun) //'.dqwf'

     call diropn (iudqwf, filint, lrdwf, exst)
  elseif (isw.eq.2) then
     !
     ! This is to be used after gen_dwf(2)
     !
     if (lgamma) call errore (' close_open ', ' isw=2 ; lgamma', 1)
     close (unit = iud0qwf, status = 'keep')
     filint = trim(filpun) //'.d0wf'
     call diropn (iud0qwf, filint, lrdwf, exst)
  elseif (isw.eq.4) then
     !
     ! This is to be used after gen_dpdvp
     !

     if (degauss.eq.0.d0) return
     close (unit = iudpdvp_1, status = 'keep')
     filint = trim(filpun) //'.pv1'

     call diropn (iudpdvp_1, filint, lrdpdvp, exst)
     if (.not.lgamma) then
        close (unit = iudpdvp_2, status = 'keep')
        filint = trim(filpun) //'.pv2'

        call diropn (iudpdvp_2, filint, lrdpdvp, exst)
        close (unit = iudpdvp_3, status = 'keep')
        filint = trim(filpun) //'.pv3'
        call diropn (iudpdvp_3, filint, lrdpdvp, exst)
     endif

  endif
  return
end subroutine close_open
