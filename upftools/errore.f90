
subroutine errore(a,b,n)
  character(len=*) :: a,b

  if (n.ne.0) then
     write(6,'(//'' program '',a,'':'',a,''.'',8x,i8,8x,''stop'')') a,b,n
     stop
  end if
end subroutine errore
