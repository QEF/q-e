subroutine test_input_xml(myunit,lxml)
!
implicit none
!
integer, intent(in) :: myunit
logical, intent(out) :: lxml
!
character(len=256) :: dummy
character :: dummy2(1:256)
integer :: i, j
!
lxml = .false.
dummy = ""
dummy2(:) = ""
!
do while (LEN_TRIM(dummy)<1)
  read(myunit,'(A256)',END=10) dummy
  do i=1,LEN_TRIM(dummy)
    dummy2(i) = dummy(i:i)
  enddo
  if(ANY(dummy2(:)=="<")) lxml=.true.
end do
!
RETURN
!
10 write(0,*) "from test_input_xml: Empty input file .. stopping"
STOP
!
end subroutine test_input_xml
