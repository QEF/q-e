!
! Copyright (C) 2001-2006 Quantum-espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!---------------------------------------------------------------------
program fhi2upf  
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential file in Fritz-Haber numerical format
  !     either ".cpi" (fhi88pp) or ".fhi" (abinit)
  !     to unified pseudopotential format
  !     Restrictions:
  !     - no semicore states
  !     Adapted from the converter written by Andrea Ferretti 
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open (unit = 1, file = filein, status = 'old', form = 'formatted')
  call read_fhi(1)
  close (1)

  ! convert variables read from FHI format into those needed
  ! by the upf format - add missing quantities

  call convert_fhi

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

stop
20 write (6,'("fhi2upf: error reading pseudopotential file name")')
   stop
end program fhi2upf

module fhi
  !
  ! All variables read from FHI file format
  !

  type angular_comp
     real(8), pointer     :: pot(:)
     real(8), pointer     :: wfc(:)
     real(8), pointer     :: grid(:)
     real(8)              :: amesh
     integer             :: nmesh
     integer             :: lcomp
  end type angular_comp

  !------------------------------

  real(8) :: Zval           ! valence charge
  integer      :: lmax_          ! max l-component used

  logical      :: nlcc_
  real(8), allocatable :: rho_atc_(:) ! core  charge

  type (angular_comp), pointer :: comp(:)  ! PP numerical info
                                           ! (wfc, grid, potentials...)
  !------------------------------

  ! variables for the abinit header

  real(8) :: Zatom, Zion, r2well, rchrg, fchrg, qchrg
  integer :: pspdat = 0, pspcod = 0 , pspxc = 0, lloc_ = -1, mmax = 0
  character(len=256) :: info = ' '

end module fhi
! 
!     ----------------------------------------------------------
subroutine read_fhi(iunps)
  !     ----------------------------------------------------------
  ! 
  use fhi
  implicit none
  integer, parameter    :: Nl=7  ! max number of l-components
  integer :: iunps
  real(8) :: r, rhoc, drhoc, d2rhoc
  !
  
  integer               :: l, i, idum, mesh

  ! Start reading file

  read(iunps,'(a)') info
  read(info,*,iostat=i) Zval, l
  if ( i /= 0 ) then
     write (6,'("read_fhi: assuming abinit format")')
     read(iunps,*) Zatom, Zion, pspdat
     read(iunps,*) pspcod, pspxc, lmax_,lloc_, mmax, r2well
     if (pspcod /= 6) then
        write (6,'("read_fhi: unknown PP type ",i1,"...stopping")') pspcod
        stop
     end if
     read(iunps,*) rchrg, fchrg, qchrg
     !
     read(iunps,*)
     read(iunps,*)
     read(iunps,*)
     !
     read(iunps,*) Zval, l
     if (abs(Zion-Zval) > 1.0d-8) then
        write (6,'("read_fhi: Zval/Zion mismatch...stopping")')
        stop
     end if
     if (l-1 /= lmax_) then
        write (6,'("read_fhi: lmax mismatch...stopping")')
        stop
     end if
  end if
  lmax_ = l - 1

  if (lmax_+1 > Nl) then
     write (6,'("read_fhi: too many l-components...stopping")')
     stop
  end if

  do i=1,10
     read(iunps,*)     ! skipping 11 lines 
  end do

  allocate( comp(0:lmax_) )

  do l=0,lmax_
     comp(l)%lcomp = l
     read(iunps,*) comp(l)%nmesh, comp(l)%amesh
     if (mmax > 0 .and. mmax /= comp(l)%nmesh) then
        write (6,'("read_fhi: mismatched number of grid points...stopping")')
        stop
     end if
     if ( l > 0) then
        if (comp(l)%nmesh /= comp(0)%nmesh .or.   &
            comp(l)%amesh /= comp(0)%amesh )      then
           write(6,'("read_fhi: different radial grids not allowed...stopping")')
           stop
        end if
     end if
     mesh = comp(l)%nmesh
     allocate( comp(l)%wfc(mesh),            &      ! wave-functions
               comp(l)%pot(mesh),            &      ! potentials
               comp(l)%grid(mesh)            )      ! real space radial grid
     ! read the above quantities
     do i=1,mesh
        read(iunps,*) idum, comp(l)%grid(i),   &
                            comp(l)%wfc(i),    &
                            comp(l)%pot(i)       
     end do
  end do

  nlcc_ =.false.
  allocate(rho_atc_(comp(0)%nmesh))
  mesh = comp(0)%nmesh
  do i=1,mesh
     read(iunps,*,end=10, err=20) r, rho_atc_(i), drhoc, d2rhoc
     if ( abs( r - comp(0)%grid(i) ) > 1.d-6 ) then
        write(6,'("read_fhi: radial grid for core charge? stopping")')
        stop
     end if
  end do
  nlcc_ = .true.
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential with NLCC successfully read'
  !     ----------------------------------------------------------
  return
10 continue
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential without NLCC successfully read'
  !     ----------------------------------------------------------
  return
  !
20 write(6,'("read_fhi: error reading core charge")')
  stop
  !
100  write(6,'("read_fhi: error reading pseudopotential file")')
  stop

end subroutine read_fhi

!     ----------------------------------------------------------
subroutine convert_fhi
  !     ----------------------------------------------------------
  !
  use fhi
  use upf
  use funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  implicit none
  real(8), parameter :: rmax = 10.0
  real(8), allocatable :: aux(:)
  real(8) :: vll
  character (len=20):: dft  
  character (len=2), external:: atom_name
  integer :: lloc, kkbeta
  integer :: l, i, ir, iv
  !
  if (nint(Zatom) > 0) then
     psd = atom_name(nint(Zatom))
  else
     print '("Atom name > ",$)'
     read (5,'(a)') psd
  end if
  if ( lloc_ < 0 ) then
     print '("l local (max: ",i1,") > ",$)', lmax_
     read (5,*) lloc
  else
     lloc = lloc_
  end if
  if (pspxc == 7) then
     dft = 'PW'
  else
     if (pspxc > 0) then
        print '("DFT read from abinit file: ",i1)', pspxc
     end if
     print '("DFT > ",$)'
     read (5,'(a)') dft
  end if
  write(generated, '("Generated using Fritz-Haber code")')
  write(date_author,'("Author: unknown    Generation date: as well")')
  if (trim(info) /= ' ') then
     comment = trim(info)
  else
     comment = 'Info: automatically converted from FHI format'
  end if
  ! reasonable assumption
  rel = 1
  rcloc = 0.0
  nwfs  = lmax_+1
  allocate( els(nwfs), oc(nwfs), epseu(nwfs))
  allocate(lchi(nwfs), nns(nwfs) )
  allocate(rcut (nwfs), rcutus (nwfs))
  do i=1, nwfs
     print '("Wavefunction # ",i1,": label, occupancy > ",$)', i
     read (5,*) els(i), oc(i)
     nns (i)  = 0
     lchi(i)  = i-1
     rcut(i)  = 0.0
     rcutus(i)= 0.0
     epseu(i) = 0.0
  end do

  pseudotype = 'NC'
  nlcc = nlcc_
  zp   = Zval
  etotps = 0.0
  ecutrho=0.0
  ecutwfc=0.0
  if ( lmax_ == lloc) then
     lmax = lmax_-1
  else
     lmax = lmax_
  end if
  nbeta= lmax_
  mesh = comp(0)%nmesh
  ntwfc= nwfs
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  do i=1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  end do
  call set_dft_from_name(dft)
  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  allocate(rab(mesh))
  allocate(  r(mesh))
  r = comp(0)%grid
  rab = r * log( comp(0)%amesh )

  if (nlcc) then
     allocate (rho_atc(mesh))
     rho_atc(:) = rho_atc_(:) / (4.d0*3.141592653589793d0)
  end if

  allocate (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0 = 2.d0*comp(lloc)%pot

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
              ! FHI potentials are in Hartree 
              betar(ir,iv) = 2.d0 * comp(l)%wfc(ir) * &
                   ( comp(l)%pot(ir) - comp(lloc)%pot(ir) )
              aux(ir) = comp(l)%wfc(ir) * betar(ir,iv)
           end do
           call simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        end if
     enddo

  end if

  allocate (rho_at(mesh))
  rho_at = 0.d0
  do i=1,nwfs
     l=lchi(i)
     rho_at = rho_at + ocw(i) * comp(l)%wfc ** 2
  end do
  
  allocate (chi(mesh,ntwfc))
  do i=1,ntwfc
     chi(:,i) = comp(i-1)%wfc(:)
  end do
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  return
end subroutine convert_fhi
