!
! Copyright (C) 2004-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine el_config &
     ( config, rel, lsd, all_elec, nwf, el, nn, ll, oc, isw, jj)
  !---------------------------------------------------------------
  !
  ! read electronic configuration from string "config"
  ! example: config = '[Ar] 3d10 4s2 4p2.5'
  ! Replicated wavefunctions are allowed and correspond to:
  !   spin-up, spin-down for lsd;
  !   j=l-1/2, j=j+1/2 for relativistic calculation and l > 0 ;
  !   error in all other cases
  ! See file INPUT_LD1 for full syntax
  !
  use kinds, only: dp
  use ld1_parameters
  use io_global, only : qestdin
  implicit none
  ! input: electronic configuration
  character(len=*), intent(in):: config
  logical, intent(in) :: all_elec
  integer, intent(in) :: lsd, rel
  ! output: atomic states
  character(len=2), intent(out) :: el(nwfx)
  integer, intent(out) :: nwf, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP), intent(out) :: oc(nwfx), jj(nwfx)
  ! local variables
  logical :: toomany
  integer ::  i, n, l, len, n0, first, start(nwfx), finish(nwfx)
  character ::  occup*10, core*2, prev*1, curr*1
  character(len=1), external :: capital
  ! core states
  character(len=2) :: elc(15)
  integer :: nwfc, nnc(15), llc(15)
  real(DP) :: occ(15)
  data elc/'1S','2S','2P','3S','3P','4S','4P','3D','5S','5P','4D', &
       &   '6S','5D','4F','6P'/
  data nnc/1 ,2 ,2 ,3 ,3 ,4 ,4 ,3 ,5 ,5 ,4 ,6 ,5 ,4 ,6 /
  data llc/0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 0, 2, 3 ,1 /
  data occ/2.0,2.0,6.0,2.0,6.0,2.0,6.0,10.0,2.0,6.0,10.0,2.0,10.0,14.0,6.0/
  !
  ! len is the length of the string, excluding trailing blanks
  !
  len=len_trim(config)
  if (len == 0) call errore('el_config','empty string',1)
  !
  ! first is the position of the first nonblank character in the string
  !
  do first=1,len
     if (config(first:first) /= ' ') go to 10
  end do
10 continue
  !
  ! find core wavefunctions (if any)
  !
  if (all_elec .and. &
       &    config(first  :first  ) == '[' .and. &
       &    config(first+3:first+3) == ']') then
     core = config(first+1:first+2)
     if (core == 'He'.or.core == 'he'.or.core == 'HE') then
        nwfc =1
     else if (core == 'Ne'.or.core == 'ne'.or.core == 'NE') then
        nwfc =3
     else if (core == 'Ar'.or.core == 'ar'.or.core == 'AR') then
        nwfc =5
     else if (core == 'Kr'.or.core == 'kr'.or.core == 'KR') then
        nwfc =8
     else if (core == 'Xe'.or.core == 'xe'.or.core == 'XE') then
        nwfc =11
     else if (core == 'Hg'.or.core == 'hg'.or.core == 'HG') then
        nwfc =14
     else if (core == 'Rn'.or.core == 'rn'.or.core == 'RN') then
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
        if (config(i:i) /= ' ') go to 20
     end do
20   first=i
  else
     nwfc=0
     first=1
  end if
  if (first >= len) then
     nwf=nwfc
     return
  endif
  !
  ! now count and find the various terms "NLf" (eg: 1s2.0)
  ! N=main quantum number, L=angular quantum number, f=occupancy
  !
  nwf=nwfc+1
  if (nwf > nwfx) call errore('config','too many states',nwf)
  start(nwf)=first
  do i=first+1,len
     prev=config(i-1:i-1)
     curr=config(i  :i  )
     if ((curr == ' ' .or.  curr == ',') .and.  &
       &  prev /= ' ' .and. prev /= ',') then
        finish(nwf)=i-1
     nwf=nwf+1
  else if (curr /= ' ' .and. curr /= ',' .and.  &
       &  (prev == ' ' .or.  prev == ',')) then
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
   if (nn(n) <= 0 .or. nn(n) > 7) &
  &    call errore('el_config','wrong main quantum number',n)

   curr = capital(config(start(n)+1:start(n)+1))
   ll(n)=-1
   if (curr == 'S') ll(n)=0
   if (curr == 'P') ll(n)=1
   if (curr == 'D') ll(n)=2
   if (curr == 'F') ll(n)=3
   if (ll(n) == -1) call errore('el_config','l not found:'//curr,n)
   if (ll(n).ge.nn(n)) call errore('el_config', &
       &              'main/angular quantum number mismatch',n)

   el(n)=prev//curr

   if (start(n)+2 > finish(n))  &
   &   call errore('el_config','no occupancy field?',n)
   occup = config(start(n)+2:finish(n))
   read(occup,*) oc(n)
   if (oc(n) > 2*(2*ll(n)+1)) &
        call errore('el_config','wrong occupancy:'//occup,n)

enddo
!
! check same labels corresponding to different spin or j value
!
do n=nwfc+1,nwf
   isw(n)=1
   jj(n) =0.0_dp
   toomany = .false.
   do i=1,n-1
      if (el(i) == el(n)) then
         !
         if (rel == 2) then
            !
            ! relativistic case: j=l-1/2, j=l+1/2
            !
            toomany = (jj(i) > 0.0_dp) .or. ll(n) == 0
            !
            jj(i)=ll(i)-0.5_dp
            jj(n)=ll(n)+0.5_dp
            if ( oc(n) > (2.0_dp*jj(n)+1.0_dp) ) &
                 call errore('el_config','occupation wrong',n)
            if ( oc(i) > (2.0_dp*jj(i)+1.0_dp) ) &
                 call errore('el_config','occupation wrong',i)
         else  if (lsd == 1) then
            !
            ! lsda case: spin up, spin down
            !
            toomany = (isw(i) > 1)
            !
            isw(n)=2
         else
            !
            ! lda case: this shouldn't happen
            !
            toomany = .true.
         end if
      endif
   enddo
   !
   if (toomany) call errore('el_config',  &
        'wavefunction '//el(n)//' found too many times',n)
   !
end do
!
! pseudopotentials: N quantum number .ne. wavefunction label
! we find the lowest N for each L and reassign N as it
!
if (.not.all_elec) then
   do l=0,3
      n0=1000
      do n=1,nwf
         if (ll(n) == l) n0=min(n0,nn(n))
      end do
      do n=1,nwf
         if (ll(n) == l) nn(n)=nn(n)-n0+1+l
      end do
   end do
endif

return
end subroutine el_config
!
!---------------------------------------------------------------
subroutine read_config(rel, lsd, nwf, el, nn, ll, oc, isw, jj)
  !---------------------------------------------------------------
  !
  ! read electronic configuration from input. Syntax:
  !    nwf
  ! followed by nwf cards
  !    label N L occ  (example: 4S 4 0 2.0) for lda
  ! or
  !    label N L occ spin (example: 4S 4 0 1.0 1) for lsda
  ! or
  !    label N L occ j (example: 4P 4 1 1.0 1.5) for full relativistic
  ! See file INPUT_LD1 for full syntax
  !
  use kinds, only: dp
  use ld1_parameters, only: nwfx
  use io_global, only : qestdin
  implicit none
  ! input
  integer :: rel, lsd 
  ! output: atomic states
  character(len=2) :: el(nwfx)
  integer :: nwf, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), jj(nwfx)
  ! local variables
  integer :: ios, n, ncheck
  character (len=2) :: label
  character (len=1), external :: capital
  !
  !
  read(qestdin,*,err=200,iostat=ios) nwf
200 call errore('read_config','reading nwf ',abs(ios))
  if (nwf <= 0) call errore('read_config','nwf is wrong',1)
  if (nwf > nwfx) call errore('read_config','too many wfcs',1)
  !
  !     read the occupation of the states
  !
  do n=1,nwf  
     if (rel < 2) then
        jj(n) = 0.0_dp
        if (lsd == 0) then
           read(qestdin,*,err=20,end=20,iostat=ios) &
                el(n), nn(n), ll(n), oc(n)
           isw(n)=1
20         call errore('read_config','reading orbital (lda)',abs(ios))
        else  
           read(qestdin,*,err=21,end=21,iostat=ios) &
                el(n), nn(n), ll(n), oc(n), isw(n)
21         call errore('read_config','reading orbital (lsd)',abs(ios))
           if(isw(n) > 2 .or. isw(n) < 1) &
                call errore('read_config','spin variable wrong ',n)
        endif
     else
        read(qestdin,*,err=22,end=22,iostat=ios) &
             el(n), nn(n), ll(n), oc(n), jj(n)
        isw(n)=1
        if ((abs(ll(n)+0.5_dp-jj(n)) > 1.e-3_dp) .and. &
            (abs(ll(n)-0.5_dp-jj(n)) > 1.e-3_dp) .and. abs(jj(n)) > 1.e-3_dp) &
            call errore('read_config','jj wrong',n)
        if (oc(n) > (2.0_dp*jj(n)+1.0_dp) .and. abs(jj(n)) > 1e-3_dp) &
             call errore('read_config','occupations wrong',n)
22      call errore('read_config','reading orbital (rel)',abs(ios))
     endif
     !
     ! Check: no two same wavefunctions
     !
     do ncheck=1,n-1
        if ( el(ncheck) == el(n) .and. isw(ncheck) == isw(n) .and. &
             jj(ncheck) == jj(n) ) then
           call errore('read_config', &
                'same wavefunction '//el(n)//' appears twice',n)
        endif
     enddo
     !
     ! More sanity checks
     !
     write(label,'(a2)') el(n)
     read (label,'(i1)') ncheck
     if (ncheck /= nn(n)  .or. &
         capital(label(2:2)) == 'S' .and. ll(n) /= 0 .or. &
         capital(label(2:2)) == 'P' .and. ll(n) /= 1 .or. &
         capital(label(2:2)) == 'D' .and. ll(n) /= 2 .or. &
         capital(label(2:2)) == 'F' .and. ll(n) /= 3 .or. &
         oc(n) > 2.0_dp*(2*ll(n)+1) .or. nn(n) < ll(n)+1  ) &
         call errore('read_config',label//' wrong?',n)
  enddo
  !
  return
end subroutine read_config
!
!---------------------------------------------------------------
subroutine read_psconfig (rel, lsd, nwfs, els, nns, lls, ocs, &
     isws, jjs, enls, rcut, rcutus )
  !---------------------------------------------------------------
  ! read pseudostates and pseudization information from input
  ! See file INPUT_LD1 for full syntax
  !
  use kinds, only: dp
  use ld1_parameters, only: nwfsx
  use io_global, only : qestdin
  implicit none
  ! input
  integer :: rel, lsd 
  ! output: atomic states
  character(len=2) :: els(nwfsx)
  integer :: nwfs, nns(nwfsx), lls(nwfsx), isws(nwfsx)
  real(DP) :: ocs(nwfsx), jjs(nwfsx), enls(nwfsx), &
       rcut(nwfsx), rcutus(nwfsx)
  ! local variables
  integer :: ios, n
  character (len=2) :: label
  character (len=1), external :: capital

  read(qestdin,*,err=600,iostat=ios) nwfs
600 call errore('read_psconfig','reading number of pseudo wavefunctions (nwfs)',abs(ios))

  if (nwfs <= 0 .or. nwfs > nwfsx) &
       call errore('read_psconfig','number of pseudo wavefunctions is wrong',1)

  do n=1,nwfs
     if (rel < 2) then
        if (lsd == 1) then
           read(qestdin,*,err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n), isws(n)
           if (isws(n) > 2 .or. isws(n) < 1) &
                call errore('read_psconfig', 'spin variable wrong',n)
           if (ocs(n) > (2.0_dp*lls(n)+1.0_dp))                 &
             call errore('read_psconfig','occupations (ls) wrong',n)
        else
           read(qestdin,*,err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n)
           isws(n)=1
           if (ocs(n) > 2.0_dp*(2.0_dp*lls(n)+1.0_dp))                 &
             call errore('read_psconfig','occupations (l) wrong',n)
        end if
        jjs(n)=0.0_dp
     else
        read(qestdin,*,err=30,end=30,iostat=ios) &
             els(n), nns(n), lls(n), ocs(n), enls(n),     &
             rcut(n), rcutus(n), jjs(n)
        isws(n)=1
        if ((abs(lls(n)+0.5_dp-jjs(n)) > 1.e-3_dp).and.      &
            (abs(lls(n)-0.5_dp-jjs(n)) > 1.e-3_dp).and. abs(jjs(n)) > 1.e-3_dp) &
             call errore('read_psconfig', 'jjs wrong',n)
        if (ocs(n) > (2.0_dp*jjs(n)+1.0_dp).and. abs(jjs(n)) > 1.e-3_dp) &
             call errore('read_psconfig','occupations (j) wrong',n)
     endif
     write(label,'(a2)') els(n)
     if ( capital(label(2:2)) == 'S'.and.lls(n) /= 0.or.   &
          capital(label(2:2)) == 'P'.and.lls(n) /= 1.or.   &
          capital(label(2:2)) == 'D'.and.lls(n) /= 2.or.   &
          capital(label(2:2)) == 'F'.and.lls(n) /= 3.or.   &
          ocs(n) > 2*(2*lls(n)+1).or.                 &
          nns(n) < lls(n)+1 )                         &
          call errore('read_psconfig','ps-label'//' wrong?',n)
     if (rcut(n) > rcutus(n)) &
          call errore('read_psconfig','rcut or rcutus is wrong',1)
  enddo
30 call errore('read_psconfig','reading pseudo wavefunctions configuration',abs(ios))
  !
  return
end subroutine read_psconfig
!------------------------------------------------------------------------
