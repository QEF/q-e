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
program rrkj2upf  
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in "rrkj3" format
  !     (Rabe-Rappe-Kaxiras-Joannopoulos with 3 Bessel functions)
  !     to unified pseudopotential format
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open (unit = 1, file = filein, status = 'old', form = 'formatted')
  call read_rrkj(1)
  close (1)

  ! convert variables read from rrkj3 format into those needed
  ! by the upf format - add missing quantities

  call convert_rrkj

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

stop
20 write (6,'("rrkj2upf: error reading pseudopotential file name")')
   stop

end program rrkj2upf

module rrkj3
  !
  ! All variables read from RRKJ3 file format
  ! 
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  character(len=75):: titleps
  character (len=2), allocatable :: els_(:)
  integer :: pseudotype_, iexch_, icorr_, igcx_, igcc_, mesh_, &
       nwfs_, nbeta_, lmax_
  logical :: rel_, nlcc_
  real (8) :: zp_, etotps_, xmin, rmax, zmesh, dx, rcloc_
  integer, allocatable:: lchi_(:), nns_(:), ikk2_(:)
  real (8), allocatable :: rcut_(:), rcutus_(:), oc_(:), &
       beta(:,:), dion_(:,:), qqq_(:,:), ddd(:,:), qfunc_(:,:,:), &
       rho_atc_(:), rho_at_(:), chi_(:,:), vloc_(:)
end module rrkj3
! 
!     ----------------------------------------------------------
subroutine read_rrkj(iunps)
  !     ----------------------------------------------------------
  ! 
  use rrkj3
  implicit none
  integer :: iunps
  integer :: nb, mb, n, ir, ios

  !---  > Start the header reading
  read (iunps, '(a75)', err = 100) titleps  
  read (iunps, *, err = 100)  pseudotype_
  read (iunps, *, err = 100) rel_, nlcc_
  read (iunps, *, err=100) iexch_, icorr_, igcx_, igcc_
  read (iunps, '(2e17.11,i5)') zp_, etotps_, lmax_
  read (iunps, '(4e17.11,i5)', err=100) xmin, rmax, zmesh, dx, mesh_
  read (iunps, *, err=100) nwfs_, nbeta_

  allocate(rcut_(nwfs_), rcutus_(nwfs_))
  read (iunps, *, err=100) (rcut_(nb), nb=1,nwfs_)
  read (iunps, *, err=100) (rcutus_(nb), nb=1,nwfs_)

  allocate(els_(nwfs_), nns_(nwfs_), lchi_(nwfs_), oc_(nwfs_))
  do nb = 1, nwfs_ 
     read (iunps, '(a2,2i3,f6.2)', err = 100) els_(nb), &
          nns_(nb), lchi_(nb) , oc_(nb)
  enddo

  allocate(ikk2_(nbeta_))
  allocate(beta( mesh_,nbeta_))
  allocate(dion_(nbeta_,nbeta_))
  allocate(ddd (nbeta_,nbeta_))
  allocate(qqq_(nbeta_,nbeta_))
  allocate(qfunc_(mesh_,nbeta_,nbeta_))

  do nb = 1, nbeta_
     read (iunps, *, err = 100) ikk2_(nb)   
     read (iunps, *, err = 100) (beta (ir, nb) , ir = 1,ikk2_(nb) )
     do ir = ikk2_(nb) + 1, mesh_
        beta (ir, nb) = 0.d0  
     enddo
     do mb = 1, nb  
        read (iunps, *, err = 100) dion_(nb, mb)
        dion_(mb, nb) = dion_(nb, mb)  
        if (pseudotype_.eq.3) then  
           read (iunps, *, err = 100) qqq_(nb, mb)
           qqq_(mb, nb) = qqq_(nb, mb)  
           read (iunps, *, err = 100) (qfunc_(n,nb, mb), n = 1, mesh_)
           do n = 1, mesh_   
              qfunc_(n, mb, nb) = qfunc_(n, nb, mb)  
           enddo
        else 
           qqq_(nb, mb) = 0.d0  
           qqq_(mb, nb) = 0.d0  
           do n = 1, mesh_   
              qfunc_(n, nb, mb) = 0.d0  
              qfunc_(n, mb, nb) = 0.d0  
           enddo
        endif
     enddo
  enddo
  !
  !     read the local potential
  !
  allocate(vloc_(mesh_))
  read (iunps, *, err = 100) rcloc_, (vloc_(ir ) , ir = 1, mesh_ )
  !
  !     read the atomic charge
  !
  allocate(rho_at_(mesh_))
  read (iunps, *, err=100) (rho_at_(ir), ir=1,mesh_)
  !
  !     if present read the core charge
  !
  allocate(rho_atc_(mesh_))
  if (nlcc_) then  
     read (iunps, *, err=100) (rho_atc_(ir), ir=1, mesh_)
  endif
  !
  !     read the pseudo wavefunctions of the atom
  !
  allocate(chi_(mesh_,nwfs_))
  read (iunps, *, err=100) ( (chi_(ir,nb), ir = 1,mesh_) , nb = 1, nwfs_)
  !
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
  return
100 write (6,'("read_rrkj: error reading pseudopotential file")')
    stop

end subroutine read_rrkj

subroutine convert_rrkj
  !     ----------------------------------------------------------
  !
  use rrkj3
  use upf
  implicit none
  integer i, n
  real(8), parameter:: pi=3.141592653589793d0
  real(8) :: x


  write(generated, '("Generated using Andrea Dal Corso code (rrkj3)")')
  write(date_author,'("Author: Andrea Dal Corso   Generation date: unknown")')
  comment = 'Info:'//titleps
  if (rel_) then
     rel = 1
  else
     rel = 0
  end if
  rcloc = rcloc_
  nwfs = nwfs_
  allocate( els(nwfs), oc(nwfs), epseu(nwfs))
  allocate(lchi(nwfs), nns(nwfs) )
  allocate(rcut (nwfs), rcutus (nwfs))
  do i=1, nwfs
     nns (i)  = nns_(i)
     lchi(i)  = lchi_(i)
     rcut(i)  = rcut_(i)
     rcutus(i)= rcutus_(i)
     oc (i)   = oc_(i)
     els(i)   = els_(i)
     epseu(i) = 0.0
  end do
  deallocate (els_, oc_, rcutus_, rcut_, nns_)

  psd  = titleps (7:8)  
  if (pseudotype_.eq.3) then
     pseudotype = 'US'
  else
     pseudotype = 'NC'
  endif
  nlcc = nlcc_
  zp = zp_
  etotps = etotps_
  ecutrho=0.0
  ecutwfc=0.0
  lmax = lmax_
  mesh = mesh_
  nbeta = nbeta_
  ntwfc = 0
  do i=1, nwfs
     if (oc(i) .gt. 1.0e-12) ntwfc = ntwfc + 1
  end do
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  n = 0
  do i=1, nwfs
     if (oc(i) .gt. 1.0e-12) then
        n = n + 1
        elsw(n) = els(i)
        ocw (n) = oc (i)
        lchiw(n)=lchi(i)
     end if
  end do
  iexch = iexch_
  icorr = icorr_
  igcx  = igcx_
  igcc  = igcc_

  allocate(rab(mesh))
  allocate(  r(mesh))
  ! define logarithmic mesh
  do i = 1, mesh
     x = xmin + DBLE(i-1) * dx
     r  (i) = exp(x) / zmesh
     rab(i) = dx * r(i)
  end do

  allocate (rho_atc(mesh))
  ! rrkj rho_core(r) =  4pi*r^2*rho_core(r) UPF
  rho_atc (:) = rho_atc_(:) / 4.0 / pi / r(:)**2
  deallocate (rho_atc_)

  allocate (vloc0(mesh))
  vloc0 = vloc_
  deallocate (vloc_)

  allocate(ikk2(nbeta), lll(nbeta))
  ikk2 = ikk2_
  lll  = lchi_
  deallocate (ikk2_, lchi_)
!  kkbeta  = 0  
!  do nb=1,nbeta
!     kkbeta  = max (kkbeta , ikk2(nb) )  
!  end do
  allocate(betar(mesh,nbeta))
  betar = 0.0
  do i=1, nbeta
     betar(1:ikk2(i),i) = beta(1:ikk2(i),i)
  end do
  deallocate (beta)

  allocate(dion(nbeta,nbeta))
  dion = dion_
  deallocate (dion_)

  allocate(qqq(nbeta,nbeta))
  qqq = qqq_
  deallocate (qqq_)

  allocate(qfunc(mesh,nbeta,nbeta))
  qfunc = qfunc_

  nqf = 0
  nqlc= 0

  allocate (rho_at(mesh))
  rho_at = rho_at_
  deallocate (rho_at_)

  allocate (chi(mesh,ntwfc))
  n = 0
  do i=1, nwfs
     if (oc(i) .gt. 1.0e-12) then
        n = n + 1
        chi(:,n) = chi_(:,i)
     end if
  end do
  deallocate (chi_)

  return
end subroutine convert_rrkj

