!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine init_cond (nregion, flag)
!
!  This subroutine computes and allocates the z mesh and
!  the local potential for the left/right leads or for the
!  scattering region
!
!  input:
!    nregion  -  number of regions to divide the unit cell
!    flag     -  'l'/'s'/'r'/'t' if the unit cell containes
!                the left lead/scat. reg./right lead/all of them
!
  USE io_global,  ONLY : stdout
  USE uspp_param, ONLY : upf, nbetam
  USE atom,       ONLY: rgrid
  USE ions_base,  ONLY : atm, nat, ityp, ntyp => nsp, tau
  USE cell_base,  ONLY : at, bg, omega, alat
  USE ener,       ONLY : ef
  USE gvecw,      ONLY : ecutwfc
  USE fft_base,   ONLY : dfftp, dffts
  USE noncollin_module,       ONLY : noncolin, npol
  USE cond

  implicit none

  character(len=1) :: flag
  integer :: nregion, nrztot, iz, naux, k, mmax, &
             nt, ib, ir, nrz1, info, na
  real(DP), parameter :: epsbeta=1.d-4, eps=1.d-8
  real(DP) :: zlen, dz1, dz2, bd1, bd2, bd3, bd4, bd5, &
                   bmax
  real(DP), allocatable :: ztot(:), rsph(:,:), dwid(:), &
              nrzreg(:)
  complex(DP), allocatable :: vppottot(:,:,:,:)

  nrx = dffts%nr1
  nry = dffts%nr2
  nrztot = dffts%nr3
!  if(nrztot/2*2.eq.nrztot) nrztot = nrztot+1
  zlen = at(3,3)
  dz1 = zlen/nrztot
  sarea = abs(at(1,1)*at(2,2)-at(2,1)*at(1,2))*alat**2
  if (abs(ecut2d).le.eps) ecut2d = ecutwfc

  allocate ( ztot(nrztot+1) )
  allocate ( rsph(nbetam, ntyp) )
  allocate ( dwid(5) )
  allocate ( nrzreg(4) )

  bd1 = 0.d0
  bd3 = zlen
  bd4 = zlen
  bd5 = zlen
  if(flag.eq.'l') then
    bd2 = bdl
  elseif(flag.eq.'s') then
    bd2 = bds
  elseif(flag.eq.'r') then
    bd2 = bdr
  elseif(flag.eq.'t') then
    bd2 = bdl
    if(nregion.gt.1) then
      bd3 = bds
    endif
    if(nregion.gt.2) then
      bd4 = bdr
    endif
  endif
  if(bd2.le.1.d-6) bd2 = zlen
  if(bd3.le.1.d-6) bd3 = zlen
  if(bd4.le.1.d-6) bd4 = zlen

  dwid(1) = bd1
  dwid(2) = bd2
  dwid(3) = bd3
  dwid(4) = bd4
  dwid(5) = bd5

  nrz1 = 0
  do iz = 2, 5
     naux=(dwid(iz)+dz1*0.5d0)/dz1-nrz1
     nrzreg(iz-1) = naux
     if (naux.gt.0) then
       dz2=(dwid(iz)-dwid(iz-1))/naux
       do k=1, naux
         ztot(nrz1+k)=dwid(iz-1)+(k-1)*dz2
       enddo
       nrz1 = nrz1 + naux
     endif
  enddo
  if(nrz1.ne.nrztot) CALL errore ('in init_cond','wrong nrztot',info)
  ztot(nrztot+1) = zlen

  allocate (vppottot(nrztot, nrx*nry, npol, npol))

  call poten(vppottot,nrztot,ztot)


!
! to determine radii of nonlocal spheres
!
  mmax = 0
  do nt=1, ntyp
    do ib=1, upf(nt)%nbeta
      mmax = max(mmax, upf(nt)%lll(ib))
      bmax=0.d0
      do ir=2, rgrid(nt)%mesh
         bmax=max(bmax, upf(nT)%beta(ir,ib)/rgrid(nt)%r(ir))
      enddo
      ir=rgrid(nt)%mesh
      do while (abs(upf(nt)%beta(ir,ib)/rgrid(nt)%r(ir)).le.epsbeta*bmax)
        ir=ir-1
      enddo
      rsph(ib,nt)=rgrid(nt)%r(ir)/alat
    enddo
  enddo

  if (mmax.gt.3) call errore ('allocate','for l>3 -orbitals  &
        &                  the method is not implemented',1)
!
! We set up the radii of US spheres to be the same (to avoid
! the problem with the spheres crossing or not the boundaries)
!
  do nt=1, ntyp
   if (upf(nt)%tvanp) then
      bmax=0.d0
      do ib=1, upf(nt)%nbeta
         bmax=max(bmax, rsph(ib,nt))
      enddo
      do ib=1, upf(nt)%nbeta
         rsph(ib,nt)=bmax
      enddo
   endif
  enddo

!
! Move all atoms into the unit cell
!
  do na=1, nat
    if(tau(3,na).gt.zlen) tau(3,na)=tau(3,na)-zlen
    if(tau(3,na).le.0) tau(3,na)=tau(3,na)+zlen
  enddo

!----------------
! Some output

  write(stdout,*)
  if(flag.eq.'l') then
    write(stdout,'(''===== INPUT FILE containing the left lead ====='')')
  elseif(flag.eq.'s') then
    write(stdout,'(''===== INPUT FILE containing the scat. region ====='')')
  elseif(flag.eq.'r') then
    write(stdout,'(''===== INPUT FILE containing the right lead ====='')')
  elseif(flag.eq.'t') then
    write(stdout,'(''===== INPUT FILE containing all the regions ====='')')
  endif

  write(stdout,'(/,5x,''GEOMETRY:'')')
  write (stdout, 100) alat, omega, sarea, zlen, nat, ntyp
100 format (/,5x,                                                     &
     &   'lattice parameter (alat)  = ',f12.4,'  a.u.',/,5x,          &
     &   'the volume                = ',f12.4,' (a.u.)^3',/,5x,       &
     &   'the cross section         = ',f12.4,' (a.u.)^2',/,5x,       &
     &   'l of the unit cell        = ',f12.4,' (alat)',/,5x,         &
     &   'number of atoms/cell      = ',i12,/,5x,                     &
     &   'number of atomic types    = ',i12,/,5x)

  write(stdout,'(5x,"crystal axes: (cart. coord. in units of alat)",/,&
      &     3(15x,"a(",i1,") = (",3f8.4," )  ",/ ) )')                &
      &     ( na, (at(nt,na), nt=1,3), na=1,3)

  write(stdout,'(/,3x,"Cartesian axes")')
  write(stdout, '(/,5x,"site n.     atom        ",  &
      &           "          positions (alat units)")')
  write(stdout, '(7x,i4,8x,a6," tau(",i3,")=(",3f8.4,"  )")')   &
      &         ( na,atm(ityp(na)),na,                          &
      &         ( tau(nt,na),nt=1,3),na=1,nat )
  write (stdout, 300) dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x,    &
                 dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x
300   format (/,5x,                                         &
        &      'nr1s                      = ',i12,/,5x,     &
        &      'nr2s                      = ',i12,/,5x,     &
        &      'nr3s                      = ',i12,/,5x,     &
        &      'nr1sx                     = ',i12,/,5x,     &
        &      'nr2sx                     = ',i12,/,5x,     &
        &      'nr3sx                     = ',i12,/,5x,     &
        &      'nr1                       = ',i12,/,5x,     &
        &      'nr2                       = ',i12,/,5x,     &
        &      'nr3                       = ',i12,/,5x,     &
        &      'nr1x                      = ',i12,/,5x,     &
        &      'nr2x                      = ',i12,/,5x,     &
        &      'nr3x                      = ',i12,/,5x)

  write(stdout,*) '_______________________________'
  write(stdout,*) ' Radii of nonlocal spheres: '
  write(stdout, '(/,5x,"type       ibeta     ang. mom.",  &
      &           "          radius (alat units)")')
  write(stdout, '(7x,a6,3x,i3,7x,i3,14x,f12.4)')                     &
      &        ( ( atm(nt), ib, upf(nt)%lll(ib), rsph(ib,nt),        &
      &         ib=1,upf(nt)%nbeta ), nt=1,ntyp)

!-----------------------------

  if(flag.eq.'l') then
    nrzl = nrzreg(1)
    allocate( vppotl(nrzl, nrx * nry, npol, npol) )
    allocate( zl(nrzl+1) )
    call potz_split(vppottot,ztot,vppotl,zl,nrztot,nrzl,nrx*nry,npol,0)
    call init_orbitals(zlen,bd1,bd2,zl,nrzl,rsph,1)
    efl = ef
  elseif(flag.eq.'s') then
    nrzs = nrzreg(1)
    allocate( vppots(nrzs, nrx * nry, npol, npol) )
    allocate( zs(nrzs+1) )
    call potz_split(vppottot,ztot,vppots,zs,nrztot,nrzs,nrx*nry,npol,0)
    call init_orbitals(zlen,bd1,bd2,zs,nrzs,rsph,2)
    efs = ef
  elseif(flag.eq.'r') then
    nrzr = nrzreg(1)
    allocate( vppotr(nrzr, nrx * nry, npol, npol) )
    allocate( zr(nrzr+1) )
    call potz_split(vppottot,ztot,vppotr,zr,nrztot,nrzr,nrx*nry,npol,0)
    call init_orbitals(zlen,bd1,bd2,zr,nrzr,rsph,3)
    efr = ef
  elseif(flag.eq.'t') then
    nrzl = nrzreg(1)
    allocate( vppotl(nrzl, nrx * nry, npol, npol) )
    allocate( zl(nrzl+1) )
    call potz_split(vppottot,ztot,vppotl,zl,nrztot,nrzl,nrx*nry,npol,0)
    call init_orbitals(zlen,bd1,bd2,zl,nrzl,rsph,1)
    if(nregion.gt.1) then
      nrzs = nrzreg(2)
      allocate( vppots(nrzs, nrx * nry, npol, npol) )
      allocate( zs(nrzs+1) )
      call potz_split(vppottot,ztot,vppots,zs,nrztot,nrzs,nrx*nry,  &
                      npol,nrzl)
      call init_orbitals(zlen,bd2,bd3,zs,nrzs,rsph,2)
    endif
    if(nregion.gt.2) then
      nrzr = nrzreg(3)
      allocate( vppotr(nrzr, nrx * nry, npol, npol) )
      allocate( zr(nrzr+1) )
      call potz_split(vppottot,ztot,vppotr,zr,nrztot,nrzr,nrx*nry,  &
                      npol,nrzl+nrzs)
      call init_orbitals(zlen,bd3,bd4,zr,nrzr,rsph,3)
    endif
    efl = ef
    efs = ef
    efr = ef
  endif

  deallocate(dwid)
  deallocate (ztot)
  deallocate (rsph)
  deallocate (nrzreg)
  deallocate (vppottot)

  return
end subroutine init_cond


subroutine potz_split(vppottot,ztot,vppot,z,nrztot,nrz,nrxy,npol,iz0)
!
! vppottot and ztot --> vppot and z
!
  use kinds, only : DP
  implicit none
  integer :: nrztot, nrz, nrxy, npol, iz0, iz, ixy, ipol1, ipol2
  real(DP) :: ztot(nrztot+1), z(nrz+1), zinit
  complex(DP) ::   vppottot (nrztot, nrxy, npol, npol), &
                        vppot (nrz, nrxy, npol, npol)
  do iz = 1, nrz
   do ixy = 1, nrxy
     do ipol1 = 1, npol
      do ipol2 = 1, npol
       vppot(iz,ixy,ipol1,ipol2) = vppottot(iz0+iz,ixy,ipol1,ipol2)
      enddo
     enddo
   enddo
  enddo

  zinit = ztot(iz0+1)
  do iz = 1, nrz+1
   z(iz) = ztot(iz0+iz) - zinit
  enddo

  return
end subroutine potz_split

