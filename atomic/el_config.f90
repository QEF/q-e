!
!---------------------------------------------------------------
subroutine el_config(config,all_elec,nwf,el,nn,ll,oc,isw)
  !---------------------------------------------------------------
  !
  use kinds, only: dp
  use ld1_parameters
  implicit none
  ! input
  character(len=*):: config
  logical :: all_elec
  ! output: atomic states
  character(len=2) :: el(nwfx)
  integer :: nwf, nn(nwfx), ll(nwfx), isw(nwfx)
  real(kind=dp) :: oc(nwfx)
  ! local variables
  integer ::  i, n, l, len, n0, first, start(18), finish(18)
  character ::  occup*10, core*2, prev*1, curr*1, capital*1
  external :: capital
  ! core states
  character(len=2) :: elc(15)
  integer :: nwfc, nnc(15), llc(15)
  real(kind=dp) :: occ(15)
  data elc/'1S','2S','2P','3S','3P','4S','4P','3D','5S','5P','4D', &
       &   '6S','6P','5D','4F'/
  data nnc/1 ,2 ,2 ,3 ,3 ,4 ,4 ,3 ,5 ,5 ,4 ,6 ,6 ,5 ,4 /
  data llc/0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 3/
  data occ/2.0,2.0,6.0,2.0,6.0,2.0,6.0,10.0,2.0,6.0,10.0,2.0,6.0,10.0,14.0/
  !
  ! len is the length of the string, excluding trailing blanks
  !
  len=len_trim(config)
  if (len.eq.0) call errore('el_config','empty string',1)
  !
  ! first is the position of the first nonblank character in the string
  !
  do first=1,len
     if (config(first:first).ne.' ') go to 10
  end do
10 continue
  !
  ! find core wavefunctions (if any)
  !
  if (all_elec .and. &
       &    config(first  :first  ).eq.'[' .and. &
       &    config(first+3:first+3).eq.']') then
     core = config(first+1:first+2)
     if (core.eq.'He'.or.core.eq.'he'.or.core.eq.'HE') then
        nwfc =1
     else if (core.eq.'Ne'.or.core.eq.'ne'.or.core.eq.'NE') then
        nwfc =3
     else if (core.eq.'Ar'.or.core.eq.'ar'.or.core.eq.'AR') then
        nwfc =5
     else if (core.eq.'Kr'.or.core.eq.'kr'.or.core.eq.'KR') then
        nwfc =8
     else if (core.eq.'Xe'.or.core.eq.'xe'.or.core.eq.'XE') then
        nwfc =11
     else if (core.eq.'Rn'.or.core.eq.'rn'.or.core.eq.'RN') then
        nwfc =15
     else
        call errore('el_config','wrong core: '//core,1)
     end if
     do n=1,nwfc
        el(n)=elc(n)
        ll(n)=llc(n)
        nn(n)=nnc(n)
        oc(n)=occ(n)
        isw(n)=1
     end do
     ! find first nonblank character after core term
     do i=first+4,len
        if (config(i:i).ne.' ') go to 20
     end do
20   first=i
  else
     nwfc=0
     first=1
  end if
  if (first.ge.len) then
     nwf=nwfc
     return
  endif
  !
  ! now count and find the various terms "NLf" (eg: 1s2.0)
  ! N=main quantum number, L=angular quantum number, f=occupancy
  !
  nwf=nwfc+1
  start(nwf)=first
  do i=first+1,len
     prev=config(i-1:i-1)
     curr=config(i  :i  )
     if ((curr.eq.' '.or. curr.eq.',') .and.  &
          &  prev.ne.' '.and.prev.ne.',') then
        finish(nwf)=i-1
     nwf=nwf+1
  else if (curr.ne.' '.and.curr.ne.',' .and.  &
       &  (prev.eq.' '.or. prev.eq.',')) then
     start(nwf)=i
  endif
enddo
finish(nwf)=len
!
! extract N,L,f from "NLf" terms !
!
do n=nwfc+1,nwf
   prev = config(start(n):start(n))
   read(prev,*) nn(n)
   if (nn(n).le.0 .or. nn(n).gt.7) &
  &    call errore('el_config','wrong main quantum number',n)

   curr = capital(config(start(n)+1:start(n)+1))
   ll(n)=-1
   if (curr.eq.'S') ll(n)=0
   if (curr.eq.'P') ll(n)=1
   if (curr.eq.'D') ll(n)=2
   if (curr.eq.'F') ll(n)=3
   if (ll(n).eq.-1) call errore('el_config','l not found:'//curr,n)
   if (ll(n).ge.nn(n)) call errore('el_config', &
       &              'main/angular quantum number mismatch',n)

   el(n)=prev//curr
   isw(n)=1
   do i=1,n-1
      if (el(i).eq.el(n)) then
         if (isw(i).eq.2) then
            call errore('el_config',  &
        &         'wavefunction '//el(n)//' found too many times',n)
         else
            isw(n)=2
         endif   
      endif
   enddo

   if (start(n)+2.gt.finish(n))  &
   &   call errore('el_config','no occupancy field?',n)
      occup = config(start(n)+2:finish(n))
      read(occup,*) oc(n)
      if (oc(n).gt.2*(2*ll(n)+1)) &
   &     call errore('el_config','wrong occupancy:'//occup,n)
enddo
! pseudopotentials: N quantum number .ne. wavefunction label
! we find the lowest N for each L and reassign N as it
if (.not.all_elec) then
   do l=0,3
      n0=1000
      do n=1,nwf
         if (ll(n).eq.l) n0=min(n0,nn(n))
      end do
      do n=1,nwf
         if (ll(n).eq.l) nn(n)=nn(n)-n0+1+l
      end do
   end do
endif

return
end
