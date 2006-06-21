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
program ncpp2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in PWSCF format
  !     (norm-conserving) to unified pseudopotential format

  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open(unit=1,file=filein,status='old',form='formatted')
  call read_ncpp(1)
  close (unit=1)

  ! convert variables read from NCPP format into those needed
  ! by the upf format - add missing quantities

  call convert_ncpp

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in US format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)
  stop
20 call errore ('ncpp2upf', 'Reading pseudo file name ', 1)

end program ncpp2upf

module ncpp
  !
  ! All variables read from NCPP file format
  ! 
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  character(len=20) :: dft_
  character(len=2)  :: psd_
  real(8) :: zp_
  integer nlc, nnl, lmax_, lloc, nchi
  logical :: numeric, bhstype, nlcc_
  real(8) :: alpc(2), cc(2), alps(3,0:3), aps(6,0:3)
  real(8) :: a_nlcc, b_nlcc, alpha_nlcc

  real(8) :: zmesh, xmin, dx
  real(8), allocatable::  r_(:), rab_(:)
  integer :: mesh_

  real(8), allocatable::  vnl(:,:), rho_atc_(:), rho_at_(:)
  integer, allocatable:: lchi_(:)
  real(8), allocatable:: chi_(:,:),  oc_(:)

end module ncpp
! 
!     ----------------------------------------------------------
subroutine read_ncpp(iunps)
  !     ----------------------------------------------------------
  ! 
  use ncpp
  use upf , only : els
  implicit none
  integer :: iunps
  !
  character(len=1), dimension(0:3) :: convel=(/'S','P','D','F'/)
  character(len=2) :: label
  real(8), parameter:: pi=3.141592653589793d0
  real (8) :: x, erf
  integer :: l, i, ir, nb, n
  character (len=255) line
  external erf

  read(iunps, '(a)', end=300, err=300 ) dft_
  if (dft_(1:2).eq.'**') dft_ = 'PZ'

  read (iunps, *, err=300) psd_, zp_, lmax_, nlc, nnl, nlcc_, &
       lloc, bhstype
  if ( nlc.gt.2 .or. nnl.gt.3) &
       call errore( 'read_ncpp','Wrong nlc or nnl',1 )
  if ( nlc* nnl .lt. 0 ) &
       call errore( 'read_ncpp','nlc*nnl < 0 ? ',1 )
  if ( zp_.le.0d0 ) &
      call errore( 'read_ncpp','Wrong zp ',1 )
  if ( lmax_.gt.3.or.lmax_.lt.0 ) &
      call errore( 'read_ncpp','Wrong lmax ',1 )

  if (lloc.eq.-1000) lloc=lmax_
  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric = nlc.le.0 .and. nnl.le.0

  if (.not.numeric) then
     !
     !   read pseudopotentials in analytic form
     !
     read(iunps, *, err=300) &
          ( alpc(i), i=1, 2 ), ( cc(i), i=1,2 )
     if ( abs(cc(1)+cc(2)-1.d0).gt.1.0d-6) &
          call errore ('read_ncpp','wrong pseudopotential coefficients',1)
     do l = 0, lmax_
        read (iunps, *, err=300) &
             ( alps(i,l),i=1,3 ), (aps(i,l),i=1,6)
     enddo

     if (nlcc_) then
        read(iunps, *, err=300) &
             a_nlcc, b_nlcc, alpha_nlcc
        if (alpha_nlcc.le.0.d0) &
             call errore('read_ncpp','nlcc but alpha=0',1)
     end if

     if (bhstype) call bachel(alps,aps,1,lmax_)
  end if

  read(iunps, *, err=300) zmesh, xmin, dx, mesh_, nchi

  if ( mesh_.le.0) call errore( 'read_ncpp', 'mesh too small', 1)
  if ( (nchi.lt.lmax_   .and. lloc.eq.lmax_).or.          &
       (nchi.lt.lmax_+1 .and. lloc.ne.lmax_)     )        &
       call errore( 'read_ncpp', 'wrong no. of wfcts', 1 )
  !
  !    compute the radial mesh
  !
  allocate(  r_(mesh_))
  allocate(rab_(mesh_))

  do ir = 1, mesh_
     x = xmin + DBLE(ir-1) * dx
     r_  (ir) = exp(x) / zmesh
     rab_(ir) = dx * r_(ir)
  end do

  allocate(vnl(mesh_,0:lmax_))
  if (numeric) then
     !
     !  read pseudopotentials in numeric form
     !
     do l = 0, lmax_
        read(iunps, '(a)', err=300)
        read(iunps, *, err=300)  (vnl(ir,l),ir=1,mesh_)
     enddo

     allocate(rho_atc_(mesh_))
     if(nlcc_) then
        read(iunps, *, err=300) ( rho_atc_(ir), ir=1,mesh_ )
     endif

  else
     !
     !  convert analytic to numeric form
     !
     do l=0,lmax_
     !
     !  DO NOT USE f90 ARRAY SYNTAX: erf IS NOT AN INTRINSIC FUNCTION!!!
     !
        do ir=1,mesh_
           vnl(ir,l)= - ( cc(1)*erf(sqrt(alpc(1))*r_(ir)) +     &
                          cc(2)*erf(sqrt(alpc(2))*r_(ir)) ) * zp_/r_(ir)
        end do

        do n=1,nnl
           vnl(:,l)= vnl(:,l)+ (aps(n,l)+ aps(n+3,l)*r_(:)**2 )* &
                exp(-alps(n,l)*r_(:)**2)
        end do
        !
        ! convert to Rydberg
        !
        vnl(:,l) = vnl(:,l)*2.0
     end do

     allocate(rho_atc_(mesh_))
     if (nlcc_) then
        rho_atc_(:) =(a_nlcc+b_nlcc*(r_(:)**2))*exp(-alpha_nlcc*r_(:)**2)
        where(abs(rho_atc_) < 1.0e-15)
           rho_atc_ = 0
        end where
     end if
  endif
  !
  ! subtract the local part
  !
   do l = 0, lmax_
     if ( l.ne.lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
  enddo
  !
  ! read pseudowavefunctions
  !
  allocate(lchi_(nchi), els(nchi))
  allocate(oc_(nchi))
  allocate(chi_(mesh_,nchi))
  do nb = 1, nchi
     ! read wavefunction label and store for later
     read(iunps, '(a)', err=300) line
     read(iunps, *, err=300) lchi_( nb), oc_( nb )
     !
     !     Test lchi and occupation numbers
     !
     if ( nb.le.lmax_.and.lchi_(nb)+1.ne.nb) &
          call errore('read_ncpp','order of wavefunctions',nb)
     if (lchi_(nb).gt.lmax_ .or. lchi_(nb).lt.0) &
          call errore('read_ncpp','wrong lchi',nb)
     if ( oc_(nb).lt.0.d0 .or.            &
          oc_(nb).gt.2.d0*(2*lchi_(nb)+1)) &
             call errore('read_ncpp','wrong oc',nb)
     !
     ! parse and check wavefunction label
     read(line,'(14x,a2)', err=222, end=222) label
     if (label(2:2).ne.convel(lchi_(nb))) goto 222
     do l = 0, lmax_
        if (label(2:2).eq.convel(l)) then
           els(nb) = label(1:2)
           goto 223
        endif
     end do
222  continue
     els(nb)   = '*'//convel(lchi_(nb)) 
223  continue
     !
     ! finally read the wavefunction
     read(iunps, *, err=300) (chi_(ir,nb),ir=1,mesh_)
  enddo
  !
  !    compute the atomic charges
  !
  allocate(rho_at_(mesh_))
  rho_at_(:)=0.d0
  do nb = 1, nchi
     if( oc_(nb).ne.0.d0) &
          rho_at_(:) = rho_at_(:) + oc_(nb)*chi_(:,nb)**2
  end do
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  return

300 call errore('read_ncpp','pseudo file is empty or wrong',1)

end subroutine read_ncpp

!     ----------------------------------------------------------
subroutine convert_ncpp
  !     ----------------------------------------------------------
  use ncpp
  use upf
  use funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  implicit none
  real(8), parameter :: rmax = 10.0
  real(8), allocatable :: aux(:)
  real(8) :: vll
  integer :: kkbeta, l, iv, ir, i

  write(generated, '("Generated using ld1 code (maybe, or maybe not)")')
  write(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from PWSCF format'
  ! reasonable assumption
  if (zmesh > 18) then
     rel = 1
  else
     rel = 0
  end if
  rcloc = 0.0
  nwfs  = nchi 
  allocate( oc(nwfs), epseu(nwfs))
  allocate(lchi(nwfs), nns(nwfs) )
  allocate(rcut (nwfs), rcutus (nwfs))
  do i=1, nwfs
     nns (i)  = 0
     lchi(i)  = lchi_(i)
     rcut(i)  = 0.0
     rcutus(i)= 0.0
     oc (i)   = oc_(i)
     epseu(i) = 0.0
  end do
  deallocate (lchi_, oc_)

  psd = psd_
  pseudotype = 'NC'
  nlcc = nlcc_
  zp = zp_
  etotps = 0.0
  ecutrho=0.0
  ecutwfc=0.0
  if ( lmax_ == lloc) then
     lmax = lmax_-1
  else
     lmax = lmax_
  end if
  nbeta= lmax_
  mesh = mesh_
  ntwfc= nchi
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  do i=1, nchi
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  end do
  call set_dft_from_name(dft_)
  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  allocate(rab(mesh))
  allocate(  r(mesh))
  rab = rab_
  r = r_

  allocate (rho_atc(mesh))
  rho_atc = rho_atc_
  deallocate (rho_atc_)

  allocate (vloc0(mesh))
  vloc0(:) = vnl(:,lloc)

  if (nbeta > 0) then

     allocate(ikk2(nbeta), lll(nbeta))
     kkbeta=mesh
     do ir = 1,mesh
        if ( r(ir) > rmax ) then
           kkbeta=ir
           exit
        end if
     end do

! make sure kkbeta is odd as required for simpson
     if(mod(kkbeta,2) == 0) kkbeta=kkbeta-1 
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
     do i=1,nchi
        l=lchi(i)
        if (l.ne.lloc) then
           iv=iv+1
           lll(iv)=l
           do ir=1,kkbeta
              betar(ir,iv)=chi_(ir,i)*vnl(ir,l)
              aux(ir) = chi_(ir,i)**2*vnl(ir,l)
           end do
           call simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        end if
        if(iv >= nbeta) exit  ! skip additional pseudo wfns
     enddo
     deallocate (vnl, aux)
!
!   redetermine ikk2
!
     do iv=1,nbeta
        ikk2(iv)=kkbeta
        do ir = kkbeta,1,-1
           if ( abs(betar(ir,iv)) > 1.d-12 ) then
              ikk2(iv)=ir
              exit
           end if
        end do
     end do
  end if
  allocate (rho_at(mesh))
  rho_at = rho_at_
  deallocate (rho_at_)
  
  allocate (chi(mesh,ntwfc))
  chi = chi_
  deallocate (chi_)

  return
end subroutine convert_ncpp
