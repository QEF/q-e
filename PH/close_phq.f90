!---------------------------------------
subroutine close_phq(flag)
!----------=========--------------------
!
  USE io_files, ONLY: iunigk
  use control_flags, ONLY : twfcollect
  use phcom
  use us, only : okvan

  implicit none

  logical :: flag

  logical :: exst

  if ( twfcollect ) then
     close (unit = iuwfc, status = 'delete')
  else
     close (unit = iuwfc, status = 'keep')
  end if
  close (unit = iudwf, status = 'keep')
  close (unit = iubar, status = 'keep')
  if(okvan) close(unit = iudrhous, status = 'keep')
  if(epsil.or.zue) close (unit = iuebar, status = 'keep')
#ifdef __PARA
  if (me.ne.1) goto 100
#endif
  if (fildrho.ne.' ') close (unit = iudrho, status = 'keep')
#ifdef __PARA
100 continue
#endif

  if (flag) then
     call seqopn (iunrec, 'recover','unformatted',exst)
     close (unit=iunrec,status='delete')
  end if
  close (unit = iunigk, status = 'delete')

  return
end subroutine close_phq
