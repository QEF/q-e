!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine cond_out ()
  use io_global, only : stdout
  USE ions_base, ONLY: atm
  use lsda_mod, only: nspin
  USE noncollin_module, ONLY : noncolin, npol
  use spin_orb, only : lspinorb
  use cond

!---------------------------
! Some output 
!---------------------------

  implicit none

  integer :: iorb, ipol, k

  write(stdout,'(''----- General information -----'')')
  write(stdout,*)
  if(ikind.eq.0) then
    write(stdout,'(''----- Complex band structure calculation -----'')')
  elseif(ikind.eq.1) then
    write(stdout,'(''--- T calc. with identical leads (ikind=1) --- '')')
  elseif(ikind.eq.2) then
    write(stdout,'(''--- T calc. with different leads (ikind=2) --- '')')
  endif

  if(nspin.eq.2) then
    write(6,'(/,9x, ''LSDA calculations, spin index ='',i6)') iofspin
  endif

  if(nspin.eq.4) then
    write(6,'(/,9x, ''Noncollinear calculations'')')
    if(lspinorb)   &
    write(6,'(/,9x, ''Noncollinear calculations with spin-orbit'')')
  endif

  write (6, 300) nrx, nry, nz1
300   format (/,5x,                                         &
        &      'nrx                       = ',i12,/,5x,     &
        &      'nry                       = ',i12,/,5x,     &
        &      'nz1                       = ',i12,/,5x)
  write (6, 301) energy0, denergy, nenergy, ecut2d, ewind, epsproj
301   format (/,5x,                                         &
        &      'energy0               = ',1pe15.1,/,5x,     &
        &      'denergy               = ',1pe15.1,/,5x,     &
        &      'nenergy               = ',i10,/,5x,         &
        &      'ecut2d                = ',1pe15.1,/,5x,     &
        &      'ewind                 = ',1pe15.1,/,5x,     &
        &      'epsproj               = ',1pe15.1,/,5x)


  if(ikind.eq.1) then
    write(stdout,'(''----- Information about left/right lead -----'')')
  else
    write(stdout,'(''----- Information about left lead ----- '')')
  endif

  write (6, 200) nocrosl, noinsl, norbl, norbf, nrzl
200  format (/,5x,                                              &
       &      'nocros                   = ',i12,/,5x,          &
       &      'noins                    = ',i12,/,5x,          &
       &      'norb                     = ',i12,/,5x,          &
       &      'norbf                    = ',i12,/,5x,          &
       &      'nrz                      = ',i12,/,5x)
  write(6, '(6x,''iorb      type   ibeta   ang. mom.'',3x,  &
       &        ''m       position (a_0)'')')
  write(6,'(5x,i4,4x,i5,5x,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i4,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,tblml(1,iorb), tblml(2,iorb), tblml(3,iorb),&
       &      tblml(4,iorb), iorb,                                 &
       &    (taunewl(ipol,iorb),ipol=1,3), iorb=1, norbl )
  if(norbl.le.80) then
    write(6,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')
    do k=1, nrzl
       write(6,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,zl(k),zl(k+1),zl(k+1)-zl(k),(crosl(iorb,k),iorb=1,norbl)
    enddo
  endif

  if(ikind.eq.2) then
    write(stdout,'(''----- Information about right lead -----'')')
    write (6, 200) nocrosr, noinsr, norbr, norbf, nrzr
    write(6, '(6x,''iorb      type   ibeta   ang. mom.'',3x,  &
       &        ''m       position (a_0)'')')
    write(6,'(5x,i4,4x,i5,5x,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i4,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,tblmr(1,iorb), tblmr(2,iorb), tblmr(3,iorb),&
       &      tblmr(4,iorb), iorb,                                 &
       &    (taunewr(ipol,iorb),ipol=1,3), iorb=1, norbr )
    if(norbr.le.80) then
      write(6,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')
      do k=1, nrzr
         write(6,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,zr(k),zr(k+1),zr(k+1)-zr(k),(crosr(iorb,k),iorb=1,norbr)
      enddo
    endif
  endif

  if(ikind.gt.0) then
    write(6,'(''----- Information about scattering region -----'')')
    write (6, 201) noinss, norbs, norbf, nrzs
201  format (/,5x,                                              &
       &      'noins                    = ',i12,/,5x,          &
       &      'norb                     = ',i12,/,5x,          &
       &      'norbf                    = ',i12,/,5x,          &
       &      'nrz                      = ',i12,/,5x)
    write(6, '(6x,''iorb      type   ibeta   ang. mom.'',3x,  &
       &        ''m       position (a_0)'')')
    write(6,'(5x,i4,4x,i5,5x,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i4,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,tblms(1,iorb), tblms(2,iorb), tblms(3,iorb),&
       &      tblms(4,iorb), iorb,                                 &
       &    (taunews(ipol,iorb),ipol=1,3), iorb=1, norbs )
    if(norbs.le.80) then
      write(6,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')
      do k=1, nrzs
         write(6,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,zs(k),zs(k+1),zs(k+1)-zs(k),(cross(iorb,k),iorb=1,norbs)
      enddo
    endif
  endif 

  return
end subroutine cond_out


