module averages
   use kinds, only: dp
   implicit none
   type online_average
       logical :: is_vector = .true., initialized=.false.
       integer :: counter = 0
       real(kind=dp) :: average_v(3), average_s
       real(kind=dp) :: var_v(3), var_s
   end type

   interface online_average_do
       module procedure online_average_do_scalar, online_average_do_vector
   end interface

contains
   subroutine online_average_init(t,vector)
       implicit none
       type(online_average), intent(inout) :: t
       logical, intent(in) :: vector

       t%counter = 0
       t%initialized=.true.
       t%is_vector=vector
       t%average_v=0.0_dp
       t%average_s=0.0_dp
       t%var_v=0.0_dp
       t%var_s=0.0_dp
         
   end subroutine


   subroutine online_average_do_scalar(t,v)
       implicit none
       type(online_average), intent(inout) :: t
       real(kind=dp), intent(in) :: v
       real(kind=dp) :: delta

       if (.not. t%initialized) &
           call online_average_init(t,.false.)

       if (.not.t%is_vector) &
           t%counter = t%counter + 1
       delta = v - t%average_s
       t%average_s = t%average_s + delta/real(t%counter, kind=dp)
       t%var_s = t%var_s + delta*(v - t%average_s)
   end subroutine

   subroutine online_average_do_vector(t,v)
       implicit none
       type(online_average), intent(inout) :: t
       real(kind=dp), intent(in) :: v(3)
       real(kind=dp) :: delta(3)

       if (.not. t%initialized) &
           call online_average_init(t,.true.)

       if (t%is_vector) &
           t%counter = t%counter + 1
       delta = v - t%average_v
       t%average_v = t%average_v + delta/real(t%counter, kind=dp)
       t%var_v = t%var_v + delta*(v - t%average_v)
   end subroutine

   subroutine online_average_print(t, iun)
       implicit none
       type(online_average), intent(in) :: t
       integer,intent(in) :: iun
       if (t%is_vector) then
           write (iun, '(6E20.12)', advance='no') t%average_v(:), t%var_v(:)
       else
           write (iun, '(2E20.12)', advance='no') t%average_s, t%var_s
       end if


   end subroutine



end module
