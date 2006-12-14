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
program oldcp2upf  
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in the old CP90 format
  !     (without core correction) to unified pseudopotential format
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open (unit = 1, file = filein, status = 'old', form = 'formatted')
  call read_oldcp(1)
  close (1)

  ! convert variables read from old CP90 format into those needed
  ! by the upf format - add missing quantities

  call convert_oldcp

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

stop
20 call errore ('oldcp2upf', 'Reading pseudo file name ', 1)

end program oldcp2upf

module oldcp
  !
  ! All variables read from old CP90 file format
  !
  real(8) :: amesh, z, zv
  integer :: exfact, lloc, nbeta_, mesh_
  real(8) :: wrc1, rc1, wrc2, rc2, rcl(3,3), al(3,3), bl(3,3)
  real(8), allocatable :: r_(:), vnl(:,:), chi_(:,:)
  !
  !------------------------------

end module oldcp
! 
!     ----------------------------------------------------------
subroutine read_oldcp(iunps)
  !     ----------------------------------------------------------
  ! 
  use oldcp
  implicit none
  integer :: iunps
  !
  real(8), external :: erf
  integer :: i, l, j, jj
  !
  read(iunps,*, end=10, err=10) z, zv, nbeta_, lloc, exfact
  if (z < 1 .or. z > 100 .or. zv < 1 .or. zv > 25 ) &
       call errore ('read_oldcp','wrong potential read',1)
  read(iunps,*, end=10, err=10) wrc1, rc1, wrc2, rc2
  read(iunps,*, end=10, err=10) ( ( rcl(i,l), al(i,l), &
                     bl(i,l), i = 1, 3), l = 1, 3)
  read(iunps,*, end=10, err=10) mesh_, amesh
  allocate(r_(mesh_))
  allocate (chi_(mesh_,nbeta_))
  do l = 1, nbeta_
     if (l > 1) read(iunps,*, end=10, err=10) mesh_, amesh
     do j = 1, mesh_
        read(iunps,*, end=10, err=10) jj, r_(j), chi_(j,l)
     end do
  end do
  !
  !  convert analytic to numeric form
  !
  allocate (vnl(mesh_,0:nbeta_))
  do l=0,nbeta_
     !
     !  DO NOT USE f90 ARRAY SYNTAX: erf IS NOT AN INTRINSIC FUNCTION!!!
     !
     do j=1, mesh_
        vnl(j,l)= - (wrc1*erf(sqrt(rc1)*r_(j)) + &
                     wrc2*erf(sqrt(rc2)*r_(j)) ) * zv/r_(j)
     end do
     !
     do i=1,3
        vnl(:,l)= vnl(:,l)+ (al(i,l+1)+ bl(i,l+1)*r_(:)**2) * &
             exp(-rcl(i,l+1)*r_(:)**2)
     end do
  end do

  return
10 call errore('read_oldcp','error in reading file',1)

end subroutine read_oldcp

!     ----------------------------------------------------------
subroutine convert_oldcp
  !     ----------------------------------------------------------
  !
  use oldcp
  use upf
  implicit none
  real(8), parameter :: rmax = 10.0d0
  real(8), allocatable :: aux(:)
  real(8) :: vll
  character (len=20):: dft  
  character (len=2), external :: atom_name
  integer :: kkbeta
  integer :: l, i, ir, iv
  !
  write(generated, '("Generated using unknown code")')
  write(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from old CP90 format'
  ! reasonable assumption
  if (z > 18) then
     rel = 1
  else
     rel = 0
  end if
  rcloc = 0.0d0
  nwfs  = nbeta_
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
  psd   = atom_name (nint(z))
  pseudotype = 'NC'
  nlcc = .false.
  zp = nint(zv)
  etotps =0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  lmax  = nbeta_ - 1
  nbeta = nbeta_
  mesh  = mesh_
  ntwfc = nwfs
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  do i=1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  end do
  !
  if ( exfact.eq.0) then
     iexch=1; icorr=1; igcx=0; igcc=0 ! Perdew-Zunger
  else if ( exfact.eq.1) then
     iexch=1; icorr=3; igcx=1; igcc=3 ! Becke-Lee-Yang-Parr
  else if ( exfact.eq.2) then
     iexch=1; icorr=1; igcx=1; igcc=0 ! Becke88 exchange
  else if (exfact.eq.-5.or.exfact.eq.3) then
     iexch=1; icorr=1; igcx=1; igcc=1 ! Becke88-Perdew 86
  else if (exfact.eq.-6.or.exfact.eq.4) then
     iexch=1; icorr=4; igcx=2; igcc=2 ! Perdew-Wang 91
  else if (exfact.eq. 5) then
     iexch=1; icorr=4; igcx=3; igcc=4 ! Perdew-Becke-Erkerhof
  else
     call errore('convert','Wrong xc in pseudopotential',1)
  end if

  allocate(rab(mesh))
  allocate(  r(mesh))
  r = r_
  rab = r * log( amesh )
  !
  !  convert analytic to numeric form
  !
  !
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
end subroutine convert_oldcp
