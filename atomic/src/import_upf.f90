!
! Copyright (C) 2004-2008 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine import_upf ( )
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
  use radial_grids, only : ndmx, radial_grid_type, allocate_radial_grid, &
                           nullify_radial_grid, deallocate_radial_grid
  use ld1inc, only : file_pseudo, zval, nlcc, pseudotype, etots, lmax, lsave_wfc,&
                     zed, nbeta, betas, lls, jjs, ikk, els, rcut, rcutus, &
                     lloc, vpsloc, grid, nwfs, bmat, qq, qvan, qvanl, rhoc, &
                     rhos, phis, which_augfun, lpaw, rmatch_augfun, pawsetup, psipaw
  use funct, only: set_dft_from_name
  !
  use pseudo_types
  use paw_type
  use upf_module
  !
  implicit none
  !
  integer :: iunps, ierr, ibeta, jbeta, kbeta, l, l1, l2
  !
  !     Local variables
  !
  integer :: nb, ios
  TYPE (pseudo_upf) :: upf
  TYPE (radial_grid_type), TARGET :: rgrid
  !
  CALL nullify_pseudo_upf( upf )
  CALL nullify_radial_grid( rgrid )
  !
  ! upf%grid => rgrid  is be associated in read_upf
  !
  iunps=2
  open(unit=iunps,file=file_pseudo,status='old',form='formatted', &
       err=100, iostat=ios)
100   call errore('import_upf','open error on file '//file_pseudo,ios)
  call read_upf(upf, rgrid, ierr, unit=iunps)
  !
  if (ierr>0) &
     call errore('import_upf','reading pseudo upf',abs(ierr))
  !
  zval  = upf%zp
  nlcc = upf%nlcc
  call set_dft_from_name (upf%dft)

  if (upf%typ.eq.'NC'.OR.upf%typ.eq.'SL') then
     pseudotype=2
  else
     pseudotype=3
  endif
  lpaw = upf%tpawp

  etots=upf%etotps
  lmax = upf%lmax
  grid%mesh = upf%mesh
  call allocate_radial_grid(grid, grid%mesh)
  grid%r  (1:grid%mesh) = upf%r  (1:upf%mesh)
  grid%rab(1:grid%mesh) = upf%rab(1:upf%mesh)
  grid%r2 (1:grid%mesh) = grid%r(1:grid%mesh)**2
  grid%sqr(1:grid%mesh) = sqrt(grid%r(1:grid%mesh))
  if (.not.upf%has_so) then
     if (grid%r(1) > 0.0_dp) then
        !
        ! r(i+1) = exp(xmin)/zmesh * exp(i*dx)
        !
        grid%dx=log(grid%r(grid%mesh)/grid%r(1))/(grid%mesh-1)
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
           if (upf%q_with_l .or. lpaw) then
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
              which_augfun='PSQ'
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
  if (upf%nlcc) then
     rhoc(1:grid%mesh) = upf%rho_atc(1:upf%mesh)*fpi*grid%r2(1:upf%mesh)
  else
     rhoc(:) = 0.0_dp
  end if
  rhos=0.0_dp
  rhos (1:grid%mesh,1) = upf%rho_at (1:upf%mesh)

  !phis(1:grid%mesh,1:nwfs)=upf%chi(1:grid%mesh,1:nwfs)
  if(upf%has_wfc) then
      lsave_wfc = .true.
      phis(1:grid%mesh,1:nbeta)=upf%pswfc(1:grid%mesh,1:nbeta)
      psipaw(1:grid%mesh,1:nbeta)=upf%aewfc(1:grid%mesh,1:nbeta)
  endif
  !!! TEMP
  lloc = -1
  vpsloc(1:grid%mesh) = upf%vloc(1:upf%mesh)
  !!!
  ! paw:
  if (lpaw) then
    which_augfun = upf%paw%augshape
    rmatch_augfun = upf%paw%raug
    call allocate_pseudo_paw( pawsetup, grid%mesh, nbeta, lmax )
    CALL nullify_radial_grid( pawsetup%grid )
    call allocate_radial_grid(pawsetup%grid,grid%mesh)
    call set_pawsetup( pawsetup, upf )
  endif

  CALL deallocate_pseudo_upf( upf )
  CALL deallocate_radial_grid( rgrid )


end subroutine import_upf

SUBROUTINE set_pawsetup(pawset_, upf_)
USE kinds, ONLY : DP
USE constants, ONLY : fpi
USE paw_type, ONLY : paw_t
USE pseudo_types, ONLY: pseudo_upf
USE ld1_parameters,   ONLY: nwfsx
USE radial_grids, ONLY : radial_grid_copy
USE atomic_paw,    ONLY : compute_nonlocal_coeff_ion
IMPLICIT NONE
TYPE(paw_t), INTENT(INOUT) :: pawset_
TYPE(pseudo_upf), INTENT(IN) :: upf_
REAL(DP), ALLOCATABLE :: ddd_(:,:)
INTEGER :: mesh, nbeta,ih,jh,ijh

   nbeta=upf_%nbeta
   mesh=upf_%mesh
   pawset_%augfun=0.0_DP
   pawset_%augmom=0.0_DP
   pawset_%enl(:) = 0.0_DP
   if (upf_%has_so) then
      pawset_%jj(1:nbeta) = upf_%jjj(1:nbeta)
      pawset_%rel=2
   else
      pawset_%jj(:) = 0.0_DP
      pawset_%rel=1
   endif
   pawset_%l(1:nbeta) = upf_%lll(1:nbeta)
   pawset_%ikk(1:nbeta) = upf_%kbeta(1:nbeta)
   pawset_%oc(1:nbeta) = upf_%paw%oc(1:nbeta)
   pawset_%aewfc(1:mesh,1:nbeta) = upf_%aewfc(1:mesh,1:nbeta)
   pawset_%pswfc(1:mesh,1:nbeta) = upf_%pswfc(1:mesh,1:nbeta)
   IF (upf_%has_so) &
   pawset_%aewfc_rel(1:mesh,1:nbeta) = upf_%paw%aewfc_rel(1:mesh,1:nbeta)
   pawset_%proj(1:mesh,1:nbeta) = upf_%beta(1:mesh,1:nbeta)

   DO ih = 1,nbeta
   DO jh = ih,nbeta
      ijh = jh * (jh-1) / 2 + ih
      pawset_%augfun(1:mesh,ih,jh,0:upf_%paw%lmax_aug) = &
                           upf_%qfuncl(1:mesh,ijh,0:upf_%paw%lmax_aug)
      IF ( ih /= jh ) &
      pawset_%augfun(1:mesh,jh,ih,0:upf_%paw%lmax_aug) = &
                           upf_%qfuncl(1:mesh,ijh,0:upf_%paw%lmax_aug)
   ENDDO
   ENDDO

   pawset_%augmom(1:nbeta,1:nbeta,0:upf_%paw%lmax_aug) = & 
                  upf_%paw%augmom(1:nbeta,1:nbeta,0:upf_%paw%lmax_aug)
   pawset_%aeccharge(1:mesh) = upf_%paw%ae_rho_atc(1:mesh)*fpi*upf_%grid%r2(1:mesh)
   pawset_%psccharge(1:mesh) = upf_%rho_atc(1:mesh)*fpi*upf_%grid%r2(1:mesh)
   pawset_%pscharge(1:mesh) = upf_%rho_at(1:mesh)
   pawset_%aeloc(1:mesh) = upf_%paw%ae_vloc(1:mesh)
   pawset_%psloc(1:mesh) = upf_%vloc(1:mesh)
   pawset_%kdiff(1:nbeta,1:nbeta) = 0.0_DP
   pawset_%dion (1:nbeta,1:nbeta) = upf_%dion(1:nbeta,1:nbeta)
   pawset_%symbol=upf_%psd
   pawset_%zval=upf_%zp
   pawset_%z=upf_%zmesh
   pawset_%nlcc=upf_%nlcc
   pawset_%nwfc=upf_%nbeta
   pawset_%irc=upf_%kkbeta
   pawset_%lmax=upf_%lmax
   pawset_%rmatch_augfun=upf_%paw%raug
   CALL radial_grid_copy(upf_%grid, pawset_%grid)
!
!  The kinetic energy must be recalculated
!
   ALLOCATE(ddd_(nwfsx,nwfsx))

   CALL compute_nonlocal_coeff_ion(ddd_, pawset_)
   pawset_%kdiff(1:nbeta,1:nbeta) = upf_%dion(1:nbeta,1:nbeta)- &
                                    ddd_(1:nbeta,1:nbeta)

   DEALLOCATE(ddd_)

END SUBROUTINE set_pawsetup
