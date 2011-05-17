!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE cond_out ()
  USE io_global, ONLY : stdout
  USE ions_base, ONLY: atm
  USE lsda_mod, ONLY: nspin
  USE noncollin_module, ONLY : noncolin, npol
  USE spin_orb, ONLY : lspinorb, domag
  USE cond

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
    write(stdout,'(''--- T calc. with different leads (i2) --- '')')
  endif

  if(nspin.eq.2) then
    write(stdout,'(/,9x, ''LSDA calculations, spin index ='',i6)') iofspin
  endif

  if(nspin.eq.4) then
    if(lspinorb) then
      IF (domag) THEN
         WRITE( stdout, '(5X, "Noncollinear calculation with spin-orbit",/)')
      ELSE
         WRITE( stdout, '(5X, "Non magnetic calculation with spin-orbit",/)')
      ENDIF
    else
      write(stdout,'(/,9x, ''Noncollinear calculations'')')
    endif
  endif

  write (stdout, 300) nrx, nry, nz1
300   format (/,5x,                                         &
        &      'nrx                       = ',i12,/,5x,     &
        &      'nry                       = ',i12,/,5x,     &
        &      'nz1                       = ',i12,/,5x)
  write (stdout, 301) energy0, denergy, nenergy, ecut2d, ewind, epsproj
301   format (/,5x,                                         &
        &      'energy0               = ',1pe15.1,/,5x,     &
        &      'denergy               = ',1pe15.1,/,5x,     &
        &      'nenergy               = ',i10,/,5x,         &
        &      'ecut2d                = ',1pe15.1,/,5x,     &
        &      'ewind                 = ',1pe15.1,/,5x,     &
        &      'epsproj               = ',1pe15.1,/,5x)

!
!   Information about the k points
!
 WRITE( stdout, '(/5x,"number of k_|| points=",i5)') nkpts

 WRITE( stdout, '(23x,"cryst. coord. ")')
 DO k = 1, nkpts
    WRITE( stdout, '(8x,"k(",i5,") = (",2f12.7,"), wk =",f12.7)') k, &
         (xyk (ipol, k) , ipol = 1, 2) , wkpt (k)
 ENDDO


 IF (start_k.GT.1 .OR. last_k.LT.nkpts)   &
   WRITE(stdout,'(5x,"WARNING: computing from k(",i5,") to k(",i5,")"/)') &
   start_k, last_k

  if(ikind.eq.1) then
    write(stdout,'(''----- Information about left/right lead -----'')')
  else
    write(stdout,'(''----- Information about left lead ----- '')')
  endif

  write (stdout, 200) nocrosl, noinsl, norbl, norbf, nrzl
200  format (/,5x,                                              &
       &      'nocros                   = ',i12,/,5x,          &
       &      'noins                    = ',i12,/,5x,          &
       &      'norb                     = ',i12,/,5x,          &
       &      'norbf                    = ',i12,/,5x,          &
       &      'nrz                      = ',i12,/,5x)
  write(stdout, '(6x,''iorb      type   ibeta   ang. mom.'',3x,  &
       &        ''m       position (alat)'')')
  write(stdout,'(5x,i4,4x,i5,5x,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i4,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,tblml(1,iorb), tblml(2,iorb), tblml(3,iorb),&
       &      tblml(4,iorb), iorb,                                 &
       &    (taunewl(ipol,iorb),ipol=1,3), iorb=1, norbl )
  if(norbl.le.80) then
    write(stdout,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')
    do k=1, nrzl
       write(stdout,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,zl(k),zl(k+1),zl(k+1)-zl(k),(crosl(iorb,k),iorb=1,norbl)
    enddo
  endif

  if(ikind.eq.2) then
    write(stdout,'(''----- Information about right lead -----'')')
    write (stdout, 200) nocrosr, noinsr, norbr, norbf, nrzr
    write(stdout, '(6x,''iorb      type   ibeta   ang. mom.'',3x,  &
       &        ''m       position (alat)'')')
    write(stdout,'(5x,i4,4x,i5,5x,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i4,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,tblmr(1,iorb), tblmr(2,iorb), tblmr(3,iorb),&
       &      tblmr(4,iorb), iorb,                                 &
       &    (taunewr(ipol,iorb),ipol=1,3), iorb=1, norbr )
    if(norbr.le.80) then
      write(stdout,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')
      do k=1, nrzr
         write(stdout,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,zr(k),zr(k+1),zr(k+1)-zr(k),(crosr(iorb,k),iorb=1,norbr)
      enddo
    endif
  endif

  if(ikind.gt.0) then
    write(stdout,'(''----- Information about scattering region -----'')')
    write (stdout, 201) noinss, norbs, norbf, nrzs
201  format (/,5x,                                              &
       &      'noins                    = ',i12,/,5x,          &
       &      'norb                     = ',i12,/,5x,          &
       &      'norbf                    = ',i12,/,5x,          &
       &      'nrz                      = ',i12,/,5x)
    write(stdout, '(6x,''iorb      type   ibeta   ang. mom.'',3x,  &
       &        ''m       position (alat)'')')
    write(stdout,'(5x,i4,4x,i5,5x,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i4,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,tblms(1,iorb), tblms(2,iorb), tblms(3,iorb),&
       &      tblms(4,iorb), iorb,                                 &
       &    (taunews(ipol,iorb),ipol=1,3), iorb=1, norbs )
    if(norbs.le.80) then
      write(stdout,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')
      do k=1, nrzs
         write(stdout,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,zs(k),zs(k+1),zs(k+1)-zs(k),(cross(iorb,k),iorb=1,norbs)
      enddo
    endif
  endif

  return
end subroutine cond_out


