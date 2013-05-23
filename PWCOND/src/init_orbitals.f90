!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Mar. 2005 : In each region, the orbitals are ordered according to the
!             z coordinate of the atomic positions. (ADC)
!
subroutine init_orbitals (zlen, bd1, bd2, z, nrz, rsph, lsr)
!
! Calculates and allocates some variables describing the nonlocal
! orbitals
!
! input:
!   zlen     -  the length of the unit cell in the z direction
!   bd1, bd2 -  two boundaries of the region under interest
!   z(nrz)   -  mesh in the z direction
!   rsph     -  radii of the spheres
!   lsr      -  1/2/3 if the region is left lead/scat. reg./right lead
!

  use cond
  use lsda_mod, only: nspin
  use noncollin_module, only : noncolin
  use spin_orb, only: lspinorb
  use ions_base,  only : atm, nat, ityp, ntyp=>nsp, tau
  use uspp_param, only : upf, nbetam
  use uspp,       only : deeq, deeq_nc, qq, qq_so
  use atom,       only : rgrid

  implicit none

  integer :: noins, lnocros, rnocros, nocros, norb, na, nt, ih, ih1,&
             ioins, ilocros, irocros, orbin, orbfin, ib, lsr, nrz,  &
             m, k, ipol, iorb, iorb1, ind, is, ips
  integer, allocatable :: orbind(:,:), tblm(:,:), cros(:,:), natih(:,:)
  real(DP), parameter :: eps=1.d-8
  real(DP) :: ledge, redge, ledgel, redgel, ledger, redger, &
                   bd1, bd2, zlen, z(nrz+1), rsph(nbetam, ntyp)
  real(DP), allocatable :: taunew(:,:), zpseu(:,:,:)

  complex(DP), allocatable :: zpseu_nc(:,:,:,:)

  allocate ( orbind(nat,nbetam) )
  orbind = -1

!---------------------
! Calculate number of crossing and inside-lying orbitals
!
  noins = 0
  lnocros = 0
  rnocros = 0
  do na = 1, nat
     nt = ityp(na)
     do ib = 1, upf(nt)%nbeta
        ledge = tau(3,na)-rsph(ib,nt)
        ledgel = ledge-zlen
        ledger = ledge+zlen
        redge = tau(3,na)+rsph(ib,nt)
        redgel = redge-zlen
        redger = redge+zlen
        if (ledge.le.bd1.and.redge.gt.bd2) &
            call errore ('init_orbitals','Too big atomic spheres',1)
        if (ledge.gt.bd1.and.redge.le.bd2) then
           noins = noins+2*upf(nt)%lll(ib)+1
           orbind(na,ib) = 0

        elseif(ledge.le.bd1.and.redge.gt.bd1) then
           lnocros = lnocros+2*upf(nt)%lll(ib)+1
           orbind(na,ib) = 1
           if(ledger.le.bd2.and.redger.gt.bd2) then
             rnocros = rnocros+2*upf(nt)%lll(ib)+1
             orbind(na,ib) = 2
           endif

        elseif(ledger.le.bd2.and.redger.gt.bd2) then
           rnocros = rnocros+2*upf(nt)%lll(ib)+1
           orbind(na,ib) = 3


        elseif(ledge.le.bd2.and.redge.gt.bd2) then
           rnocros = rnocros+2*upf(nt)%lll(ib)+1
           orbind(na,ib) = 4
           if(ledgel.le.bd1.and.redgel.gt.bd1) then
             lnocros = lnocros+2*upf(nt)%lll(ib)+1
             orbind(na,ib) = 5
           endif

        elseif(ledgel.le.bd1.and.redgel.gt.bd1) then
           lnocros = lnocros+2*upf(nt)%lll(ib)+1
           orbind(na,ib) = 6

        endif
     enddo
  enddo
  norb = noins + lnocros + rnocros
  nocros = (lnocros + rnocros)/2
!------------------------------------

!-----------------------------
! Formation of some orbital arrays
!
  IF (norb>0) THEN
     allocate( taunew(4,norb) )
     allocate( tblm(4,norb) )
     allocate( natih(2,norb) )
     allocate( cros(norb, nrz) )
     if (noncolin) then
        allocate(zpseu_nc(2, norb, norb, nspin))
     else
        allocate( zpseu(2, norb, norb) )
     endif
  ENDIF

  ilocros = 0
  ioins =  lnocros
  irocros = ioins + noins

  do na = 1, nat
    nt = ityp(na)
    ih = 0
    do ib = 1, upf(nt)%nbeta
      do m = 1,2*upf(nt)%lll(ib) + 1
        ih = ih+1
        if(orbind(na,ib).eq.0) then
          ioins = ioins+1
          natih(1,ioins)=na
          natih(2,ioins)=ih
          tblm(1,ioins) = nt
          tblm(2,ioins) = ib
          tblm(3,ioins) = upf(nt)%lll(ib)
          tblm(4,ioins) = m
          do ipol = 1, 3
            taunew(ipol,ioins)=tau(ipol,na)
          enddo
          taunew(4,ioins) = rsph(ib,nt)
        endif
        if(orbind(na,ib).eq.1.or.orbind(na,ib).eq.2) then
          ilocros = ilocros + 1
          natih(1,ilocros)=na
          natih(2,ilocros)=ih
          tblm(1,ilocros) = nt
          tblm(2,ilocros) = ib
          tblm(3,ilocros) = upf(nt)%lll(ib)
          tblm(4,ilocros) = m
          do ipol = 1, 3
            taunew(ipol,ilocros)=tau(ipol,na)
          enddo
          taunew(4,ilocros) = rsph(ib,nt)
        endif
        if(orbind(na,ib).eq.2.or.orbind(na,ib).eq.3) then
          irocros = irocros + 1
          natih(1,irocros)=na
          natih(2,irocros)=ih
          tblm(1,irocros) = nt
          tblm(2,irocros) = ib
          tblm(3,irocros) = upf(nt)%lll(ib)
          tblm(4,irocros) = m
          do ipol = 1, 2
            taunew(ipol,irocros)=tau(ipol,na)
          enddo
          taunew(3,irocros) = tau(3,na) + zlen
          taunew(4,irocros) = rsph(ib,nt)
        endif
        if(orbind(na,ib).eq.4.or.orbind(na,ib).eq.5) then
          irocros = irocros + 1
          natih(1,irocros)=na
          natih(2,irocros)=ih
          tblm(1,irocros) = nt
          tblm(2,irocros) = ib
          tblm(3,irocros) = upf(nt)%lll(ib)
          tblm(4,irocros) = m
          do ipol = 1, 3
            taunew(ipol,irocros)=tau(ipol,na)
          enddo
          taunew(4,irocros) = rsph(ib,nt)
        endif
        if(orbind(na,ib).eq.5.or.orbind(na,ib).eq.6) then
          ilocros = ilocros + 1
          natih(1,ilocros)=na
          natih(2,ilocros)=ih
          tblm(1,ilocros) = nt
          tblm(2,ilocros) = ib
          tblm(3,ilocros) = upf(nt)%lll(ib)
          tblm(4,ilocros) = m
          do ipol = 1, 2
            taunew(ipol,ilocros)=tau(ipol,na)
          enddo
          taunew(3,ilocros) = tau(3,na) - zlen
          taunew(4,ilocros) = rsph(ib,nt)
        endif
      enddo
    enddo
  enddo

!
!  order orbital in order of increasing taunew
!
  do iorb=1,lnocros
     do iorb1=iorb+1,lnocros
        if (taunew(3,iorb1).lt.taunew(3,iorb)-1.d-8) then
           do ind=iorb,iorb1-1
              call exchange(natih(1,ind),tblm(1,ind),taunew(1,ind), &
                           natih(1,iorb1),tblm(1,iorb1),taunew(1,iorb1) )
           enddo
        endif
     enddo
  enddo

  do iorb=lnocros+1,lnocros+noins
     do iorb1=iorb+1,lnocros+noins
        if (taunew(3,iorb1).lt.taunew(3,iorb)-1.d-8) then
           do ind=iorb,iorb1-1
              call exchange(natih(1,ind),tblm(1,ind),taunew(1,ind), &
                         natih(1,iorb1),tblm(1,iorb1),taunew(1,iorb1) )
           enddo
        endif
     enddo
  enddo

  do iorb=lnocros+noins+1,lnocros+noins+rnocros
     do iorb1=iorb+1,lnocros+noins+rnocros
        if (taunew(3,iorb1).lt.taunew(3,iorb)-1.d-8)  then
           do ind=iorb,iorb1-1
              call exchange(natih(1,ind),tblm(1,ind),taunew(1,ind), &
                            natih(1,iorb1),tblm(1,iorb1),taunew(1,iorb1) )
           enddo
        endif
     enddo
  enddo

  do iorb = 1, norb
    taunew(3,iorb) = taunew(3,iorb) - bd1
  enddo
!--------------------------

!-------------------------
! to form the array containing the information does the orbital
! cross the given slab or not.
!
  do iorb=1, norb
    ledge = taunew(3,iorb)-taunew(4,iorb)
    redge = taunew(3,iorb)+taunew(4,iorb)
    do k=1, nrz
      if (ledge.gt.z(k+1).or.redge.lt.z(k)) then
         cros(iorb,k)=0
      else
         cros(iorb,k)=1
      endif
    enddo
  enddo
!----------------------------

!----------------------------
!    To form zpseu
!
  IF (norb>0) THEN
     if (noncolin) then
        zpseu_nc=(0.d0,0.d0)
     else
        zpseu = 0.d0
     endif
  ENDIF

  orbin = 1
  orbfin = lnocros+noins
  do k = 1, 2
   do iorb = orbin, orbfin
     nt = tblm(1,iorb)
     ib = tblm(2,iorb)
     if(upf(nt)%tvanp.or.lspinorb) then
       na = natih(1,iorb)
       ih = natih(2,iorb)
       do iorb1 = orbin, orbfin
         if (na.eq.natih(1,iorb1)) then
           ih1 = natih(2,iorb1)
           if (noncolin) then
             do is=1, nspin
               if(lspinorb) then
                zpseu_nc(1,iorb,iorb1,is)=deeq_nc(ih,ih1,na,is)
                zpseu_nc(2,iorb,iorb1,is)=qq_so(ih,ih1,is,nt)
               else
                zpseu_nc(1,iorb,iorb1,is)=deeq_nc(ih,ih1,na,is)
                zpseu_nc(2,iorb,iorb1,is)=qq(ih,ih1,nt)
               endif
             enddo
           else
             zpseu(1,iorb,iorb1)=deeq(ih,ih1,na,iofspin)
             zpseu(2,iorb,iorb1) = qq(ih,ih1,nt)
           endif
         endif
       enddo
     else
       do iorb1=orbin,orbfin
          if (natih(1,iorb)==natih(1,iorb1)) then
             na=natih(1,iorb1)
             ih=natih(2,iorb)
             ih1 = natih(2,iorb1)
             if (noncolin) then
                zpseu_nc(1,iorb,iorb1,1)=deeq(ih,ih1,na,1)
                zpseu_nc(1,iorb,iorb1,4)=deeq(ih,ih1,na,4)
             else
                zpseu(1,iorb,iorb1)=deeq(ih,ih1,na,1)
             end if
          end if
       end do
     endif
   enddo
   orbin = lnocros+noins+1
   orbfin = norb
  enddo
!--------------------------

!--------------------------
! Allocation
!
  if(lsr.eq.1) then
    norbl = norb
    nocrosl = nocros
    noinsl = noins
    if(ikind.eq.1) then
      norbr = norb
      nocrosr = nocros
      noinsr = noins
    endif
    IF (norbl>0) THEN
       allocate( taunewl(4,norbl) )
       allocate( tblml(4,norbl) )
       allocate( crosl(norbl, nrzl) )
       if (noncolin) then
          allocate(zpseul_nc(2, norbl, norbl, nspin))
       else
          allocate( zpseul(2, norbl, norbl) )
       endif
       taunewl = taunew
       tblml = tblm
       crosl = cros
       if (noncolin) then
          zpseul_nc = zpseu_nc
       else
          zpseul = zpseu
       endif
       do ips=1, ntyp
          rl(1:rgrid(ips)%mesh,ips) = rgrid(ips)%r(1:rgrid(ips)%mesh)
          rabl(1:rgrid(ips)%mesh,ips) = rgrid(ips)%rab(1:rgrid(ips)%mesh)
          betarl(1:rgrid(ips)%mesh,1:upf(ips)%nbeta,ips) = &
             upf(ips)%beta(1:rgrid(ips)%mesh,1:upf(ips)%nbeta)
       end do
    ENDIF
    norbf = norbl
  elseif(lsr.eq.2) then
    norbs = norb
    noinss = noins
    IF (norbs>0) THEN
       allocate( taunews(4,norbs) )
       allocate( tblms(4,norbs) )
       allocate( cross(norbs, nrzs) )
       if (noncolin) then
          allocate(zpseus_nc(2, norbs, norbs, nspin))
       else
          allocate( zpseus(2, norbs, norbs) )
       endif
       taunews = taunew
       tblms = tblm
       cross = cros
       if (noncolin) then
          zpseus_nc = zpseu_nc
       else
          zpseus = zpseu
       endif
       do ips=1, ntyp
          rs(1:rgrid(ips)%mesh,ips) = rgrid(ips)%r(1:rgrid(ips)%mesh)
          rabs(1:rgrid(ips)%mesh,ips) = rgrid(ips)%rab(1:rgrid(ips)%mesh)
          betars(1:rgrid(ips)%mesh,1:upf(ips)%nbeta,ips) = &
                  upf(ips)%beta(1:rgrid(ips)%mesh,1:upf(ips)%nbeta)
       end do
    ENDIF
    norbf = max(norbf,norbs)
  elseif(lsr.eq.3) then
    norbr = norb
    nocrosr = nocros
    noinsr = noins
    IF (norbr>0) THEN
       allocate( taunewr(4,norbr) )
       allocate( tblmr(4,norbr) )
       allocate( crosr(norbr, nrzr) )
       if (noncolin) then
          allocate(zpseur_nc(2, norbr, norbr, nspin))
       else
          allocate( zpseur(2, norbr, norbr) )
       endif
       taunewr = taunew
       tblmr = tblm
       crosr = cros
       if (noncolin) then
          zpseur_nc = zpseu_nc
       else
          zpseur = zpseu
       endif
       do ips=1,ntyp
          rr(1:rgrid(ips)%mesh,ips) = rgrid(ips)%r(1:rgrid(ips)%mesh)
          rabr(1:rgrid(ips)%mesh,ips) = rgrid(ips)%rab(1:rgrid(ips)%mesh)
          betarr(1:rgrid(ips)%mesh,1:upf(ips)%nbeta,ips) = &
            upf(ips)%beta(1:rgrid(ips)%mesh,1:upf(ips)%nbeta)
       enddo
    ENDIF
    norbf = max(norbf,norbr)
  endif
!---------------------------

!-- if LDA+U  
  call plus_u_setup (natih, lsr)
!--

  deallocate (orbind)
  if (norb>0) THEN
     deallocate (taunew)
     deallocate (tblm)
     deallocate (natih)
     deallocate (cros)
     if (noncolin) then
        deallocate (zpseu_nc)
     else
        deallocate (zpseu)
     endif
  endif
  return
end subroutine init_orbitals

subroutine exchange(natih1,tblm1,taunew1,natih2,tblm2,taunew2)

use kinds, only : dp
implicit none

integer :: natih1(2),natih2(2),tblm1(4),tblm2(4)
real(DP) ::taunew1(4),taunew2(4), rdum

integer :: i, idum

do i=1,2
   idum=natih1(i)
   natih1(i)=natih2(i)
   natih2(i)=idum
enddo
do i=1,4
   idum=tblm1(i)
   tblm1(i)=tblm2(i)
   tblm2(i)=idum
enddo
do i=1,4
   rdum=taunew1(i)
   taunew1(i)=taunew2(i)
   taunew2(i)=rdum
enddo
return
end subroutine exchange

