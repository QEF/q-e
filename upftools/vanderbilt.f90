
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
module Vanderbilt
  !
  ! All variables read from Vanderbilt's file format
  ! 
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  integer :: nvalps, nang, nbeta_, kkbeta, nchi, ifpcor, keyps, &
       mesh_, iver(3), idmy(3), nnlz, ifqopt, nqf_, irel,  npf, &
       nlc, lloc
  real(8) ::  z_, zp_, exfact, etot, eloc, rcloc_, rpcor, &
       qtryc, ptryc, rinner1_
  real(8), allocatable::  wwnlps(:), eeps(:), rinner_(:), rc(:), &
       beta(:,:), ddd0(:,:), ddd(:,:), qqq_(:,:), eee(:), rho_atc_(:), &
       r_(:), rab_(:), rho_at_(:), qfunc_(:,:,:), vloc(:), vloc_(:), &
       wf(:,:), qfcoef_(:,:,:,:)
  integer, allocatable :: lll_(:), nnlzps(:), iptype(:)
  Character(len=20):: title
end module Vanderbilt
! 
!     ----------------------------------------------------------
subroutine read_uspp(iunit)
  !     ----------------------------------------------------------
  !
  use Vanderbilt
  implicit none
  integer :: iunit
  !
  integer :: i, j, k, lp
  real(8) :: rinner1
  !
  !
  read (iunit) (iver(i),i=1,3),(idmy(i),i=1,3)
  read (iunit) title, z_, zp_, exfact, nvalps, mesh_, etot

  allocate(nnlzps(nvalps), wwnlps(nvalps), eeps(nvalps))
  read (iunit) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)

  read (iunit) keyps, ifpcor, rinner1

  if ( iver(1) .eq. 1 ) then
     nang = nvalps
     nqf_ = 3
     nlc = 5
  elseif ( iver(1) .eq. 2 ) then
     nang = nvalps
     nqf_ = 3
     nlc = 2 * nvalps - 1
  else if ( iver(1) .ge. 3 ) then
     read (iunit) nang, lloc, eloc, ifqopt, nqf_, qtryc
     nlc = 2 * nang - 1
  endif

  allocate(rinner_(2*nang-1))
  rinner_(1) = rinner1
  rinner1_ = rinner1
  if (10*iver(1)+iver(2).ge.51) &
       read (iunit) (rinner_(i),i=1,nang*2-1)

  if ( iver(1) .ge. 4 ) then
     read (iunit) irel
  else
     irel = 0
  end if

  allocate(rc(nang))
  read (iunit) (rc(i),i=1,nang)

  read (iunit) nbeta_,kkbeta
  !
  allocate(beta(kkbeta,nbeta_))
  allocate(qfunc_(kkbeta,nbeta_,nbeta_))
  allocate(ddd0(nbeta_,nbeta_))
  allocate(ddd (nbeta_,nbeta_))
  allocate(qqq_(nbeta_,nbeta_))
  allocate(lll_(nbeta_))
  allocate(eee(nbeta_))
  allocate(qfcoef_(nqf_,nlc,nbeta_,nbeta_))
  !
  do j=1,nbeta_
     read (iunit) lll_(j),eee(j),(beta(i,j),i=1,kkbeta)
     do k=j,nbeta_
        read (iunit) ddd0(j,k),ddd(j,k),qqq_(j,k), &
             (qfunc_(i,j,k),i=1,kkbeta), &
             ((qfcoef_(i,lp,j,k),i=1,nqf_),lp=1,2*nang-1)
     end do
  end do
  !
  allocate(iptype(nbeta_)) 
  if (10*iver(1)+iver(2).ge.72) &
       read (iunit) (iptype(j),j=1,nbeta_),npf,ptryc
  !
  allocate(vloc_(mesh_))
  read (iunit) rcloc_,(vloc_(i),i=1,mesh_)
  !
  allocate(rho_atc_(mesh_))
  if (ifpcor.gt.0) then
     read (iunit) rpcor
     read (iunit) (rho_atc_(i),i=1,mesh_)
  end if
  !
  allocate(rho_at_(mesh_), vloc(mesh_))
  read (iunit) (vloc(i),i=1,mesh_)
  read (iunit) (rho_at_(i),i=1,mesh_)

  allocate(r_(mesh_), rab_(mesh_))
  read (iunit) (r_(i),i=1,mesh_)
  read (iunit) (rab_(i),i=1,mesh_)
  if (iver(1) .ge. 6) then
     nchi = nvalps
     if (iver(1) .ge. 7)  read (iunit) nchi
     allocate(wf(mesh_,nchi))
     read (iunit) ((wf(i,j), i=1,mesh_),j=1,nchi)
  end if
  !
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
end subroutine read_uspp
!     ----------------------------------------------------------
!     ----------------------------------------------------------
subroutine read_vdb(iunit)
  !     ----------------------------------------------------------
  !
  use Vanderbilt
  implicit none
  integer :: iunit
  !
  integer :: i, j, k, lp
  real(8) :: rinner1
  !
  !
  read(iunit, *) (iver(i),i=1,3),(idmy(i),i=1,3)
  read(iunit,'(a20,3f15.9)' ) title, z_, zp_, exfact
  read(iunit, *) nvalps, mesh_, etot

  allocate(nnlzps(nvalps), wwnlps(nvalps), eeps(nvalps))
  do i = 1,nvalps
     read(iunit, *) nnlzps(i), wwnlps(i), eeps(i)
  end do

  read(iunit, *)  keyps, ifpcor, rinner1

  if ( iver(1) .eq. 1 ) then
     nang = nvalps
     nqf_ = 3
     nlc = 5
  elseif ( iver(1) .eq. 2 ) then
     nang = nvalps
     nqf_ = 3
     nlc = 2 * nvalps - 1
  else if ( iver(1) .ge. 3 ) then
     read(iunit, *) nang, lloc, eloc, ifqopt, nqf_, qtryc
     nlc = 2 * nang - 1
  endif

  allocate(rinner_(2*nang-1))
  rinner_(1) = rinner1
  if (10*iver(1)+iver(2).ge.51) &
       read (iunit, *) (rinner_(i),i=1,nang*2-1)
  if ( iver(1) .ge. 4 ) then
     read (iunit, *) irel
  else
     irel = 0
  end if

  allocate(rc(nang))
  read(iunit, *) ( rc(i), i=1,nang)

  read (iunit,* ) nbeta_, kkbeta

  allocate(beta(kkbeta,nbeta_))
  allocate(qfunc_(kkbeta,nbeta_,nbeta_))
  allocate(ddd0(nbeta_,nbeta_))
  allocate(ddd (nbeta_,nbeta_))
  allocate(qqq_(nbeta_,nbeta_))
  allocate(lll_(nbeta_))
  allocate(eee (nbeta_))
  allocate(qfcoef_(nqf_,nlc,nbeta_,nbeta_))

  do j=1,nbeta_
     read ( iunit, *) lll_(j) 
     read ( iunit, *) eee(j), ( beta(i,j), i=1,kkbeta )
     do k=j,nbeta_
        read( iunit, *) ddd0(j,k), ddd(j,k), qqq_(j,k), &
             (qfunc_(i,j,k),i=1,kkbeta),&
             ((qfcoef_(i,lp,j,k),i=1,nqf_),lp=1,2*nang-1)
     enddo
  enddo

  allocate(iptype(nbeta_)) 
  if (10*iver(1)+iver(2).ge.72) then
     read ( iunit, * ) (iptype(i), i=1,nbeta_)
     read ( iunit, * )  npf, ptryc
  end if

  allocate(vloc_(mesh_))
  read(iunit, *) rcloc_, ( vloc_(i), i=1,mesh_)
  
  allocate(rho_atc_(mesh_))
  if ( ifpcor.gt.0 ) then 
     read(iunit, *) rpcor
     read(iunit, *) ( rho_atc_(i), i=1,mesh_)
  endif

  allocate(rho_at_(mesh_), vloc(mesh_))
  read(iunit, *)  (vloc(i), i=1,mesh_)
  read(iunit, *)  (rho_at_(i), i=1,mesh_)

  allocate(r_(mesh_),rab_(mesh_))
  read(iunit, *)  (r_(i), i=1,mesh_)
  read(iunit, *)  (rab_(i),i=1,mesh_)

  if (iver(1) .ge. 6) then
     nchi = nvalps
     if (iver(1) .ge. 7)  read (iunit, *) nchi
     allocate(wf(mesh_,nchi))
     read (iunit, *) ((wf(i,j), i=1,mesh_),j=1,nchi)
  end if

  return
end subroutine read_vdb

subroutine convert_uspp
  !     ----------------------------------------------------------
  !
  use Vanderbilt
  use constants, only : fpi
  use upf
  implicit none
  integer i
  character(len=1), dimension(0:3) :: convel=(/'S','P','D','F'/)

  write(generated, '("Generated using Vanderbilt code, version ",3i3)') iver
  write(date_author,'("Author: unknown    Generation date:",3i5)') idmy
  write(comment,'("Automatically converted from original format")') 
  if (irel == 0) then
    rel = 0
  else if (irel == 1) then
    rel = 2
  else if (irel == 2) then
    rel = 1
  end if
  rcloc = rcloc_
  nwfs = nvalps
  allocate( els(nwfs), oc(nwfs), epseu(nwfs))
  allocate(lchi(nwfs), nns(nwfs) )
  allocate(rcut (nwfs), rcutus (nwfs))
  do i=1, nwfs
     nns (i)  = nnlzps(i)/100
     lchi(i)  = mod (nnlzps(i)/10,10)
     rcut(i)  = rinner1_
     rcutus(i)= rc(lchi(i)+1) 
     oc (i) = wwnlps(i)
     write(els(i),'(i1,a1)') nns(i), convel(lchi(i))
     epseu(i) = eeps(i)
  end do
  deallocate (nnlzps, rc, wwnlps, eeps)

  psd = title
  if (keyps.le.2) then
     pseudotype = 'NC'
  else
     pseudotype = 'US'
  end if
  nlcc = ifpcor.gt.0
  zp = zp_
  etotps = etot
  ecutrho=0.0d0
  ecutwfc=0.0d0
  lmax = nang - 1
  mesh = mesh_
  nbeta = nbeta_
  if (nvalps .ne. nchi) then
     print *, 'WARNING: verify info on atomic wavefunctions'
  end if
  ntwfc = nchi
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  do i=1, min(ntwfc,nwfs)
     elsw(i) = els(i)
     ocw(i)  = oc (i)
     lchiw(i)=lchi(i)
  end do
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
     write (6,'("convert: wrong xc in pseudopotential ",f12.6)') exfact
     stop
  end if

  allocate (r(mesh), rab(mesh))
  r  =  r_
  rab=rab_
  deallocate (r_, rab_)

  allocate (rho_atc(mesh))
  ! Vanderbilt rho_core(r) =  4pi*r^2*rho_core(r) UPF
  rho_atc (1) = 0.d0
  rho_atc (2:mesh) = rho_atc_(2:mesh) / fpi / r(2:mesh)**2
  deallocate (rho_atc_)

  allocate (vloc0(mesh))
  vloc0(2:mesh) = vloc_(2:mesh)/r(2:mesh)
  vloc0(1) = vloc0(2)
  deallocate (vloc_)

  allocate(ikk2(nbeta), lll(nbeta))
  ikk2 = kkbeta
  lll  = lll_
  deallocate (lll_)
  allocate(betar(kkbeta,nbeta))
  betar = beta
  deallocate (beta)

  allocate(dion(nbeta,nbeta))
  dion = ddd0
  deallocate (ddd0)

  allocate(qqq(nbeta,nbeta))
  qqq = qqq_
  deallocate (qqq_)

  allocate(qfunc(mesh,nbeta,nbeta))
  qfunc(1:kkbeta,:,:) = qfunc_(1:kkbeta,:,:)
  qfunc(kkbeta+1:mesh,:,:) = 0.d0
  deallocate (qfunc_)

  nqf = nqf_
  nqlc= nlc
  allocate(rinner(nqlc))
  rinner = rinner_
  deallocate(rinner_)
  allocate(qfcoef(nqf,nqlc,nbeta,nbeta))
  qfcoef = qfcoef_
  deallocate (qfcoef_)

  allocate (rho_at(mesh))
  rho_at = rho_at_
  deallocate (rho_at_)

  allocate (chi(mesh,ntwfc))
  chi = wf
  deallocate (wf)

  return
end subroutine convert_uspp

