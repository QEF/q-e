subroutine openfil
use ld1inc
implicit none

integer :: ios

if (file_tests.ne.' '.and.iswitch.gt.1) then
    open(unit=13,file=file_tests,status='unknown', &
              err=1110, iostat=ios,form='formatted')
1110     call errore('openfil','opening file_tests',abs(ios))
else
     open( unit = 13, file = '/dev/null', status = 'unknown' )
endif
return
end
