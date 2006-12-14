!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!---------------------------------------------------------------------
program cpmd2upf  
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in the CPMD format
  !     (TYPE=NORMCONSERVING NUMERIC only, single radial grid)
  !     to unified pseudopotential format
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open (unit = 1, file = filein, status = 'old', form = 'formatted')
  call read_cpmd(1)
  close (1)

  ! convert variables read from CPMD format into those needed
  ! by the upf format - add missing quantities

  call convert_cpmd

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

stop
20 call errore ('cpmd2upf', 'Reading pseudo file name ', 1)

end program cpmd2upf

module cpmd
  !
  ! All variables read from CPMD file format
  !
  character (len=80) title
  !
  integer :: ixc
  real(8) :: alphaxc
  integer :: z, zv
  !
  integer :: mesh_
  real(8) :: amesh
  real(8), allocatable :: r_(:)
  !
  integer ::lmax_
  real(8), allocatable :: vnl(:,:)
  real(8), allocatable :: chi_(:,:)
  !
  logical :: nlcc_
  real(8), allocatable :: rho_atc_(:)
  !
  integer :: maxinfo_, info_lines_
  parameter (maxinfo_ = 100)
  character (len=80), allocatable :: info_sect_(:)
  !------------------------------

end module cpmd
! 
!     ----------------------------------------------------------
subroutine read_cpmd(iunps)
  !     ----------------------------------------------------------
  ! 
  use cpmd
  implicit none
  integer :: iunps
  !  
  integer :: found = 0, closed = 0, unknown = 0
  integer :: i, l, ios
  character (len=80) line
  character (len=4) token
  real (8) :: vnl0(0:3)
  logical, external :: matches
  integer, external :: locate
  !
  nlcc_ = .false.
  info_lines_ = 0
10 read (iunps,'(A)',end=20,err=20) line
  if (matches ("&ATOM", trim(line)) ) then
     found = found + 1
     ! Z
     read (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     read (line(i+1:l),*) z
     ! ZV
     read (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     read (line(i+1:l),*) zv
     ! XC
     read (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     read (line(i+1:l),*) ixc, alphaxc
     ! TYPE
     read (iunps,'(a)',end=200,err=200) line
     if (.not. matches("NORMCONSERVING",line) .or. &
         .not. matches("NUMERIC",line) ) &
             call errore('read_cpmd','unknown type: '//line,1)
  else if (matches ("&INFO", trim(line)) ) then
     found = found + 1
     ! read (iunps,'(a)') title
     ! store info section for later perusal (FIXME: not yet implemented. 2004/10/12, AK)
     allocate (info_sect_(maxinfo_))
     do i=1,maxinfo_
        read (iunps,'(a)',end=20,err=20) title
        if (matches ("&END", trim(title)) ) then
           closed = closed + 1
           goto 10
        else
           info_sect_(i) = trim(title)
           info_lines_ = i
        end if
     enddo
  else if (matches ("&POTENTIAL", trim(line)) ) then
     found = found + 1
     !read (iunps,*) mesh_, amesh
     read (iunps,'(a)') line
     read (line,*,iostat=ios) mesh_, amesh
     if ( ios /= 0) then
        read (line,*,iostat=ios) mesh_
        amesh = -1.0d0
     end if
     allocate (r_(mesh_))
     !
     ! determine the number of angular momenta
     !
     read (iunps, '(a)') line
     ios = 1
     lmax_=4
     do while (ios /= 0)
        lmax_ = lmax_ - 1
        read(line,*,iostat=ios) r_(1),(vnl0(l),l=0,lmax_)
     end do
     allocate (vnl(mesh_,0:lmax_))
     vnl(1,0:lmax_) = vnl0(0:lmax_)
     do i=2,mesh_
        read(iunps, *) r_(i),(vnl(i,l),l=0,lmax_)
     end do
     ! get amesh if not available directly
     if (amesh < 0.0d0) print  "('amesh set to:',f10.6)", exp (r_(1) - r_(0))
     if (amesh < 0.0d0) amesh = exp (r_(1) - r_(0))
  else if (matches ("&WAVEFUNCTION", trim(line)) ) then
     found = found + 1
     ! read (iunps,*) mesh_, amesh
     read (iunps,'(a)') line
     read (line,*,iostat=ios) mesh_
     allocate(chi_(mesh_,lmax_+1))
     do i=1,mesh_
        read(iunps, *) r_(i),(chi_(i,l+1),l=0,lmax_)
     end do
  else if (matches ("&NLCC", trim(line)) ) then
     found = found + 1
     nlcc_ = .true.
     read (iunps, '(a)') line
     if (.not. matches ("NUMERIC", trim(line)) ) &
          call errore('read_cpmd',' only NUMERIC core-correction supported',1)
     read(iunps, *) mesh_
     allocate (rho_atc_(mesh_))
     read(iunps, * ) (r_(i), rho_atc_(i), i=1,mesh_)
  else if (matches ("&ATDENS", trim(line)) ) then
     ! skip over &ATDENS section, add others here, if there are more.
     do while(.not. matches("&END", trim(line)))
        read (iunps,'(a)') line
     end do
  else if (matches ("&END", trim(line)) ) then
     closed = closed + 1
  else
     print*, 'line ignored: ', line
     unknown = unknown + 1
  end if
  go to 10

20 continue
  if (nlcc_ .and. found /= 5 .or. .not.nlcc_ .and. found /= 4) &
       call errore('read_cpmd','some &FIELD card missing',found)
  if (closed /= found) &
       call errore('read_cpmd','some &END card missing',closed)
  if (unknown /= 0 ) print '("WARNING: ",i3," cards not read")', unknown

  return
200 call errore('read_cpmd','error in reading file',1)

end subroutine read_cpmd

!     ----------------------------------------------------------
subroutine convert_cpmd
  !     ----------------------------------------------------------
  !
  use cpmd
  use upf
  implicit none
  real(8), parameter :: rmax = 10.0d0
  real(8), allocatable :: aux(:)
  real(8) :: vll
  character (len=20):: dft  
  character (len=2), external :: atom_name
  integer :: lloc, kkbeta, my_lmax
  integer :: l, i, ir, iv
  !
  write(generated, '("Generated using unknown code")')
  write(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from CPMD format'

  ! NOTE: many CPMD pseudopotentials created with the 'Hamann' code
  ! from Juerg Hutter's homepage have additional (bogus) entries for
  ! pseudo-potential and wavefunction. In the 'report' they have
  ! the same rc and energy eigenvalue than the previous angular momentum.
  ! we need to be able to ignore that part or the resulting UPF file
  ! will be useless. so we first print the info section and ask
  ! for the LMAX to really use. AK 2005/03/30.
  do i=1,info_lines_
        print '(A)', info_sect_(i)
  enddo
  print '("lmax to use. (max.",I2,") > ",$)', lmax_
  read (5,*) my_lmax
  if ((my_lmax <= lmax_) .and. (my_lmax >= 0)) lmax_ = my_lmax
  print '("l local (max.",I2,") > ",$)', lmax_
  read (5,*) lloc
  ! reasonable assumption
  if (z > 18) then
     rel = 1
  else
     rel = 0
  end if
  rcloc = 0.0d0
  nwfs  = lmax_+1
  allocate( els(nwfs), oc(nwfs), epseu(nwfs))
  allocate(lchi(nwfs), nns(nwfs) )
  allocate(rcut (nwfs), rcutus (nwfs))
  do i=1, nwfs
     print '("Wavefunction # ",i1,": label, occupancy > ",$)', i
     read (5,*) els(i), oc(i)
     nns (i)  = 0
     lchi(i)  = i-1
     rcut(i)  = 0.0d0
     rcutus(i)= 0.0d0
     epseu(i) = 0.0d0
  end do
  psd   = atom_name (z)
  pseudotype = 'NC'
  nlcc = nlcc_
  zp = zv
  etotps =0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  if ( lmax_ == lloc) then
     lmax = lmax_-1
  else
     lmax = lmax_
  end if
  nbeta= lmax_
  mesh = mesh_
  ntwfc= nwfs
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  do i=1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  end do
  iexch = ixc/1000
  icorr = (ixc-1000*iexch)/100
  igcx  = (ixc-1000*iexch-100*icorr)/10
  igcc  = (ixc-1000*iexch-100*icorr-10*igcx)
  !
  ! We have igcc=2 (PW91) and 3 (LYP) exchanged wrt CPMD conventions
  !
  if (igcc.eq.3) then
     igcc=2
  else if (igcc.eq.2) then
     igcc=3
  end if

  allocate(rab(mesh))
  allocate(  r(mesh))
  r = r_
  rab = r * log( amesh )

  allocate (rho_atc(mesh))
  if (nlcc) rho_atc = rho_atc_

  allocate (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0(:) = vnl(:,lloc)*2.d0

  if (nbeta > 0) then

     allocate(ikk2(nbeta), lll(nbeta))
     kkbeta=mesh
     do ir = 1,mesh
        if ( r(ir) > rmax ) then
           kkbeta=ir
           exit
        end if
     end do
     ikk2(:) = kkbeta
     allocate(aux(kkbeta))
     allocate(betar(mesh,nbeta))
     allocate(qfunc(mesh,nbeta,nbeta))
     allocate(dion(nbeta,nbeta))
     allocate(qqq (nbeta,nbeta))
     qfunc(:,:,:)=0.0d0
     dion(:,:) =0.d0
     qqq(:,:)  =0.d0
     iv=0
     do i=1,nwfs
        l=lchi(i)
        if (l.ne.lloc) then
           iv=iv+1
           lll(iv)=l
           do ir=1,kkbeta
              ! the factor 2 converts from Hartree to Rydberg
              betar(ir,iv) = 2.d0 * chi_(ir,l+1) * &
                   ( vnl(ir,l) - vnl(ir,lloc) )
              aux(ir) = chi_(ir,l+1) * betar(ir,iv)
           end do
           call simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        end if
     enddo

  end if

  allocate (rho_at(mesh))
  rho_at = 0.d0
  do i=1,nwfs
     rho_at(:) = rho_at(:) + ocw(i) * chi_(:,i) ** 2
  end do
  
  allocate (chi(mesh,ntwfc))
  chi = chi_
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  return
end subroutine convert_cpmd
!
! ------------------------------------------------------------------
integer function locate(onechar,string)
! ------------------------------------------------------------------
  !
  character(len=1) :: onechar
  character(len=*) :: string
  !
  integer:: i
  !
  do i=1,len_trim(string)
     if (string(i:i) .eq. "=") then
        locate = i
        return
     end if
  end do
  locate = 0
  return
end function locate
