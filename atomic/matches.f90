!-----------------------------------------------------------------------
logical function matches (string1, string2)
!-----------------------------------------------------------------------
!
! .true. if string 1 is contained in string2, .false. otherwise
!
implicit none
character (len=*) :: string1, string2
integer :: len1, len2, l

len1 = len_trim(string1)
len2 = len_trim(string2)
do l = 1, len2 - len1 + 1
   if (string1 (1:len1) .eq.string2 (l:l + len1 - 1) ) then
      matches = .true.
      return
   endif
enddo

matches = .false.
return
end function matches
