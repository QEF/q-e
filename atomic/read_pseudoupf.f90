!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
!---------------------------------------------------------------------
subroutine read_pseudoupf 
  !---------------------------------------------------------------------
  !
  !   read "is"-th pseudopotential in the Unified Pseudopotential Format
  !   from unit "iunps" - convert and copy to internal PWscf variables
  !   return error code in "ierr" (success: ierr=0)
  !
  ! PWSCF modules
  !
  use constants, only : fpi
  use kinds, only : dp
  use ld1inc, only : file_pseudo, zval, nlcc, pseudotype, etots, lmax, &
                     zed, nbeta, betas, lls, jjs, ikk, els, rcut, rcutus, &
                     lloc, vpsloc, grid, nwfs, bmat, qq, qvan, qvanl, rhoc, &
                     rhos, phis, which_augfun
  use funct, only: set_dft_from_name
  !
  use pseudo_types
  use read_upf_module
  !
  implicit none
  !
  integer :: iunps, ierr, ibeta, jbeta, kbeta, l, l1, l2
  !
  !     Local variables
  !
  integer :: nb, ios
  TYPE (pseudo_upf) :: upf
  !
  !
  iunps=2
  open(unit=iunps,file=file_pseudo,status='old',form='formatted', &
       err=100, iostat=ios)
100   call errore('read_pseudoupf','open error on file '//file_pseudo,ios)

  call read_pseudo_upf(iunps, upf, ierr)
  !
  if (ierr .ne. 0) &
     call errore('read_pseudoupf','reading pseudo upf',abs(ierr))
  !
  zval  = upf%zp
  nlcc = upf%nlcc
  call set_dft_from_name (upf%dft)

  if (upf%typ.eq.'NC') then
     pseudotype=2
  else
     pseudotype=3
  endif
  etots=upf%etotps
  lmax = upf%lmax
  grid%mesh = upf%mesh
  grid%r  (1:grid%mesh) = upf%r  (1:upf%mesh)
  grid%rab(1:grid%mesh) = upf%rab(1:upf%mesh)
  grid%r2 (1:grid%mesh) = grid%r(1:grid%mesh)**2
  grid%sqr(1:grid%mesh) = sqrt(grid%r(1:grid%mesh))
  if (.not.upf%has_so) then
     if (grid%r(1) > 0.0_dp) then
        !
        ! r(i+1) = exp(xmin)/zmesh * exp(i*dx)
        !
        grid%dx=log(grid%r(2)/grid%r(1))
        grid%rmax=grid%r(grid%mesh)
        grid%xmin=log(zed*grid%r(1))
        grid%zmesh=zed
     else
        !
        ! r(i+1) = exp(xmin)/zmesh * ( exp(i*dx) - 1 )
        !
        grid%dx=log( (grid%r(3)-grid%r(2)) / grid%r(2) )
        grid%rmax=grid%r(grid%mesh)
        grid%zmesh=zed
        grid%xmin=log(zed*grid%r(2)**2/(grid%r(3)-2.0_dp*grid%r(2)))
     end if
  else
     grid%dx=upf%dx
     grid%xmin=upf%xmin
     grid%zmesh=upf%zmesh
     grid%rmax=exp(grid%xmin+(grid%mesh-1)*grid%dx)/grid%zmesh
  endif
  if (abs(exp(grid%xmin+(grid%mesh-1)*grid%dx)/zed-grid%rmax).gt.1.e-6_dp) &
       call errore('read_pseudoup','mesh not supported',1)

  nwfs = upf%nwfc

  nbeta= upf%nbeta
  lls(1:nbeta)=upf%lll(1:nbeta)

  if (upf%has_so) then
     jjs(1:nbeta)=upf%jjj(1:nbeta)
  else
     jjs=0.0_dp
  endif
  !
  !
  do nb=1,nbeta
     ikk(nb)=upf%kbeta(nb)
     els(nb)=upf%els_beta(nb)
     rcut(nb)=upf%rcut(nb)
     rcutus(nb)=upf%rcutus(nb)
  end do
  betas(1:grid%mesh, 1:nbeta) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  bmat(1:nbeta, 1:nbeta) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  if (pseudotype.eq.3) then
     qq(1:nbeta,1:nbeta) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
     do ibeta=1,nbeta
        do jbeta=ibeta,nbeta
           kbeta = jbeta * (jbeta-1) / 2 + ibeta
           if (upf%q_with_l) then
              l1=upf%lll(ibeta)
              l2=upf%lll(jbeta)
              do l=abs(l1-l2), l1+l2
                 qvanl(1:upf%mesh,ibeta,jbeta,l)=upf%qfuncl(1:upf%mesh,kbeta,l)
                 if (ibeta /= jbeta) qvanl (1:grid%mesh, jbeta, ibeta, l)= &
                                    upf%qfuncl(1:upf%mesh,kbeta,l)
              enddo
              qvan(1:upf%mesh,ibeta,jbeta)=upf%qfuncl(1:upf%mesh,kbeta,0)
              if (ibeta /= jbeta) qvan(1:grid%mesh, jbeta, ibeta)= &
                                    upf%qfuncl(1:upf%mesh,kbeta,0)
              which_augfun='BESSEL'
           else
              qvan (1:grid%mesh, ibeta, jbeta) = upf%qfunc(1:upf%mesh,kbeta)
              if (ibeta /= jbeta) qvan (1:grid%mesh, jbeta, ibeta)= &
                                       upf%qfunc(1:upf%mesh,kbeta)
              which_augfun='AE'
           endif
        enddo
     enddo
  else
     qq=0.0_dp
     qvan=0.0_dp
  endif
  !
  !
  if (upf%nlcc) then
     rhoc(1:grid%mesh) = upf%rho_atc(1:upf%mesh)*fpi*grid%r2(1:upf%mesh)
  else
     rhoc(:) = 0.0_dp
  end if
  rhos=0.0_dp
  rhos (1:grid%mesh,1) = upf%rho_at (1:upf%mesh)
  phis(1:grid%mesh,1:nwfs)=upf%chi(1:grid%mesh,1:nwfs)
  !!! TEMP
  lloc = -1
  vpsloc(1:grid%mesh) = upf%vloc(1:upf%mesh)
  !!!

  CALL deallocate_pseudo_upf( upf )

end subroutine read_pseudoupf
