
!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine init_cond
!
! This subroutine sets up some variables of PWCOND
!
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, atm, tau 
  USE io_global,  ONLY : stdout
  USE pwcom
  USE noncollin_module, ONLY : noncolin
  USE uspp_param,    ONLY : dion, nbeta, lll, tvanp
  USE cond
  implicit none
  integer :: nt, ib, ir, is, na, iorb, iorb1, m, ih, ih1, ioins, &
             iocros, k, ipol, apol 
  real(kind=DP) :: ledge, redge, ra, zorb1, zorb2
  real(kind=DP), parameter :: eps=1.d-8

!
! We order the atomic orbitals
!

!
! Left tip
!
  iocros=0
  ioins=nocrosl
  do na=1, nat
     nt=ityp(na)
     ih=0
     do ib=1, nbeta(nt)
        ledge=tau(3,na)-rsph(ib,nt)
        redge=tau(3,na)+rsph(ib,nt)
        if (tau(3,na).gt.bdl1+eps.and.tau(3,na).le.bdl2+eps) then
          if (ledge.gt.bdl1.and.redge.le.bdl2) then
           do m=1,2*lll(ib,nt)+1
              ioins=ioins+1
              ih=ih+1 
              natih(ioins,1)=na 
              natih(ioins,2)=ih 
              itnew(ioins)=nt
              nbnew(ioins)=ib
              ls(ioins)=lll(ib,nt)
              mnew(ioins)=m
              do ipol=1,3
                 taunew(ipol,ioins)=tau(ipol,na)
              enddo
           enddo
          else
           do m=1,2*lll(ib,nt)+1
              iocros=iocros+1           
              ih=ih+1                    
              natih(iocros,1)=na 
              natih(iocros,2)=ih 
              itnew(iocros)=nt
              nbnew(iocros)=ib
              ls(iocros)=lll(ib,nt)
              mnew(iocros)=m
              natih(iocros+noinsl+nocrosl,1)=na 
              natih(iocros+noinsl+nocrosl,2)=ih 
              itnew(iocros+noinsl+nocrosl)=nt
              nbnew(iocros+noinsl+nocrosl)=ib
              ls(iocros+noinsl+nocrosl)=lll(ib,nt)
              mnew(iocros+noinsl+nocrosl)=m
              do ipol=1,2
               taunew(ipol,iocros)=tau(ipol,na)
               taunew(ipol,iocros+noinsl+nocrosl)=tau(ipol,na)
              enddo
              if (ledge.le.bdl1) then
               taunew(3,iocros)=tau(3,na)
               taunew(3,iocros+noinsl+nocrosl)=tau(3,na)+bdl2-bdl1
              else
               taunew(3,iocros)=tau(3,na)-(bdl2-bdl1)
               taunew(3,iocros+noinsl+nocrosl)=tau(3,na)                
              endif
           enddo
          endif
        endif
     enddo
  enddo

  if (ikind.ne.0) then
!
! If the scattering problem
!
    ioins=2*nocrosl+noinsl
    do na=1, nat
      nt=ityp(na)
      ih=0
      do ib=1, nbeta(nt)
        ledge=tau(3,na)-rsph(ib,nt)
        redge=tau(3,na)+rsph(ib,nt)
        if (ledge.gt.bdl2.and.redge.le.bdr1) then 
          do m=1,2*lll(ib,nt)+1
            ioins=ioins+1
            ih=ih+1 
            natih(ioins,1)=na 
            natih(ioins,2)=ih 
            itnew(ioins)=nt
            nbnew(ioins)=ib
            ls(ioins)=lll(ib,nt)
            mnew(ioins)=m
            do ipol=1,3
              taunew(ipol,ioins)=tau(ipol,na)
            enddo
          enddo
        endif
      enddo
    enddo
    iocros=2*nocrosl+noinsl+noinss
    ioins=iocros+nocrosr

    if (ikind.eq.1) then
!
!   If tips are identical
!       
      do ib=1, nocrosr
        iocros=iocros+1
        natih(iocros,1)=natih(ib,1)
        natih(iocros,2)=natih(ib,2)
        itnew(iocros)=itnew(ib)
        nbnew(iocros)=nbnew(ib)
        ls(iocros)=ls(ib)
        mnew(iocros)=mnew(ib)
        do ipol=1,2
          taunew(ipol,iocros)=taunew(ipol,ib)
        enddo
        taunew(3,iocros)=taunew(3,ib)+(bdr1-bdl1)
      enddo                                   

    else
!
!   Two tips are different
!
      do na=1, nat
        nt=ityp(na)
        ih=0
        do ib=1, nbeta(nt)
          ledge=tau(3,na)-rsph(ib,nt)
          redge=tau(3,na)+rsph(ib,nt)
          if (ledge.gt.bdr1.and.redge.le.bdr2) then
            do m=1,2*lll(ib,nt)+1
              ioins=ioins+1
              ih=ih+1
              natih(ioins,1)=na
              natih(ioins,2)=ih
              itnew(ioins)=nt
              nbnew(ioins)=ib
              ls(ioins)=lll(ib,nt)
              mnew(ioins)=m
              do ipol=1,3
                taunew(ipol,ioins)=tau(ipol,na)
              enddo
            enddo
          endif
          if (ledge.le.bdr1.and.redge.gt.bdr1) then
            do m=1,2*lll(ib,nt)+1
              iocros=iocros+1
              ih=ih+1
              natih(iocros,1)=na
              natih(iocros,2)=ih
              itnew(iocros)=nt
              nbnew(iocros)=ib
              ls(iocros)=lll(ib,nt)
              mnew(iocros)=m
              natih(iocros+noinsr+nocrosr,1)=na
              natih(iocros+noinsr+nocrosr,2)=ih
              itnew(iocros+noinsr+nocrosr)=nt
              nbnew(iocros+noinsr+nocrosr)=ib
              ls(iocros+noinsr+nocrosr)=lll(ib,nt)
              mnew(iocros+noinsr+nocrosr)=m
              do ipol=1,2
               taunew(ipol,iocros)=tau(ipol,na)
               taunew(ipol,iocros+noinsr+nocrosr)=tau(ipol,na)
              enddo
              taunew(3,iocros)=tau(3,na)
              taunew(3,iocros+noinsr+nocrosr)=tau(3,na)+ &
                                              (bdr2-bdr1)
            enddo
          endif
        enddo
      enddo                                     

    endif
  endif

!
!    To form zpseu
!
  zpseu=0.d0
  if (noncolin) zpseu_nc=0.d0

  iocros=nocrosr
  if (ikind.eq.0) iocros=0

  do iorb=nocrosl+1, norb-iocros
     ib=nbnew(iorb)
     nt=itnew(iorb) 
     if(tvanp(nt).or.lspinorb) then
      na=natih(iorb,1)
      ih=natih(iorb,2)
      do iorb1=nocrosl+1, norb-iocros
        if (na.eq.natih(iorb1,1)) then
          ih1=natih(iorb1,2)               
          do is=1, nspin
            if (noncolin) then
               zpseu_nc(iorb,iorb1,is)=deeq_nc(ih,ih1,na,is) 
            else
               zpseu(iorb,iorb1,is)=deeq(ih,ih1,na,is) 
            endif
          enddo
        endif 
      enddo
     else
      do is=1, nspin
       if (noncolin) then
          zpseu_nc(iorb,iorb,is)=dion(ib,ib,nt) 
       else
          zpseu(iorb,iorb,is)=dion(ib,ib,nt) 
       endif
      enddo
     endif
  enddo 

  do iorb=1, nocrosl
    do iorb1=1, nocrosl
      do is=1, nspin
        if (noncolin) then
           zpseu_nc(iorb,iorb1,is)= &
             zpseu_nc(nocrosl+noinsl+iorb,nocrosl+noinsl+iorb1,is)
        else
           zpseu(iorb,iorb1,is)= &
             zpseu(nocrosl+noinsl+iorb,nocrosl+noinsl+iorb1,is)
        endif
      enddo
    enddo
  enddo             

  iocros=0  
  if (ikind.eq.2) iocros=2*nocrosl+noinsl+noinss

  if (ikind.ne.0) then
    do iorb=1, nocrosr
      do iorb1=1, nocrosr
        do is=1, nspin
          if (noncolin) then
             zpseu_nc(norb-nocrosr+iorb,norb-nocrosr+iorb1,is)= &
                 zpseu_nc(iocros+iorb,iocros+iorb1,is)  
          else
             zpseu(norb-nocrosr+iorb,norb-nocrosr+iorb1,is)= &
                 zpseu(iocros+iorb,iocros+iorb1,is)  
          endif
        enddo         
      enddo
    enddo 
  endif
          
!
! to form the array containing the information does the orbital
! cross the given slab or not. 
!
  do iorb=1, norb
    ra=rsph(nbnew(iorb),itnew(iorb))
    zorb1=taunew(3,iorb)-ra
    zorb2=taunew(3,iorb)+ra    
    do k=1, nrz
      if (zorb1.gt.(z(k+1)-eps).or.zorb2.lt.(z(k)+eps)) then
         cross(iorb,k)=0
      else
         cross(iorb,k)=1
      endif
    enddo
  enddo

!
! Some output by PWCOND
!
  WRITE( stdout,'(9x, ''PWCOND calculation is now starting...'')') 
  if(ikind.eq.0) then
     WRITE( stdout,'(/,2x,''CBS calculation (ikind=0)'')')
     WRITE( stdout,'(/,2x,''left  boundary is '',f9.4)') bdl1
     WRITE( stdout,'(2x,''right boundary is '',f9.4)') bdl2
  elseif(ikind.eq.1) then
     WRITE( stdout,'(/,2x,''T calculation with identical tips (ikind=1)'')')
     WRITE( stdout,'(/,2x,''left  boundary of the tip is    '',f9.4)') bdl1
     WRITE( stdout,'(2x,''right boundary of the tip is    '',f9.4)') bdl2
     WRITE( stdout,'(2x,''the end of the scatt. region is '',f9.4)') bdr1
  elseif(ikind.eq.2) then 
     WRITE( stdout,'(/,2x,''T calculation with different tips (ikind=2)'')')
     WRITE( stdout,'(/,2x,''left  boundary of the left  tip is '',f9.4)') bdl1
     WRITE( stdout,'(2x,''right boundary of the left  tip is '',f9.4)') bdl2 
     WRITE( stdout,'(2x,''left  boundary of the right tip is '',f9.4)') bdr1
     WRITE( stdout,'(2x,''right boundary of the right tip is '',f9.4)') bdr2
  endif

  WRITE( stdout,'(/,5x,''GEOMETRY:'')') 
  WRITE( stdout, 100) alat, omega, sarea, zl, nat, ntyp
100 format (/,5x,                                                     &
     &   'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x,          &
     &   'the volume                = ',f12.4,' (a.u.)^3',/,5x,       &
     &   'the cross section         = ',f12.4,' (a.u.)^2',/,5x,       &
     &   'l of the unit cell        = ',f12.4,' (a_0)',/,5x,          &
     &   'number of atoms/cell      = ',i12,/,5x,                     &
     &   'number of atomic types    = ',i12,/,5x)

  WRITE( stdout,'(5x,''crystal axes: (cart. coord. in units of a_0)'',/,    &
      &     3(15x,''a('',i1,'') = ('',3f8.4,'' )  '',/ ) )')          &
      &     ( apol, (at(ipol,apol), ipol=1,3), apol=1,3)      

  WRITE( stdout,'(/,3x,''Cartesian axes'')')
  WRITE( stdout, '(/,5x,''site n.     atom        '',  &
      &           ''          positions (a_0 units)'')')
  WRITE( stdout, '(7x,i3,8x,a6,'' tau('',i2,'')=('',3f8.4,''  )'')')  &
      &         ( na,atm(ityp(na)),na,                          &
      &         ( tau(ipol,na),ipol=1,3),na=1,nat )          
  WRITE( stdout, 300) nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,     &
                 nr1, nr2, nr3, nrx1, nrx2, nrx3,           &
                 nrx, nry, nrz, nz1
300   format (/,5x,                                         &
        &      'nr1s                      = ',i12,/,5x,     &
        &      'nr2s                      = ',i12,/,5x,     &
        &      'nr3s                      = ',i12,/,5x,     &
        &      'nrx1s                     = ',i12,/,5x,     &
        &      'nrx2s                     = ',i12,/,5x,     &
        &      'nrx3s                     = ',i12,/,5x,     &
        &      'nr1                       = ',i12,/,5x,     &
        &      'nr2                       = ',i12,/,5x,     &
        &      'nr3                       = ',i12,/,5x,     &
        &      'nrx1                      = ',i12,/,5x,     &
        &      'nrx2                      = ',i12,/,5x,     &
        &      'nrx3                      = ',i12,/,5x,     &
        &      'nrx                       = ',i12,/,5x,     &
        &      'nry                       = ',i12,/,5x,     &
        &      'nrz                       = ',i12,/,5x,     &
        &      'nz1                       = ',i12,/,5x)

  WRITE( stdout,*) '_______________________________'
  WRITE( stdout,*) ' Radii of nonlocal spheres: '
  WRITE( stdout, '(/,5x,''type       ibeta     ang. mom.'',  &
      &           ''          radius (a_0 units)'')')
  WRITE( stdout, '(7x,a6,3x,i3,7x,i3,14x,f12.4)')                     &
      &        ( ( atm(nt), ib, lll(ib,nt), rsph(ib,nt),        &
      &         ib=1,nbeta(nt) ), nt=1,ntyp)

  WRITE( stdout, 200) norb, norbf, nocrosl, noinsl, noinss, nocrosr, &
                 noinsr
200  format (/,5x,                                              &
       &      'norb                      = ',i12,/,5x,          &
       &      'norbf                     = ',i12,/,5x,          &
       &      'nocrosl                   = ',i12,/,5x,          &
       &      'noinsl                    = ',i12,/,5x,          &
       &      'noinss                    = ',i12,/,5x,          &
       &      'nocrosr                   = ',i12,/,5x,          &
       &      'noinsr                    = ',i12,/,5x)

  WRITE( stdout, '(5x,''iorb  type   ibeta   ang. mom.'',3x,  & 
       &        ''m       position (a_0)'')')
  WRITE( stdout,'(5x,i3,4x,a6,i3,6x,i3,6x,i3,''   taunew('',          &
       &    i3,'')=('',3f8.4,'')'')')                           &
       &    ( iorb,atm(itnew(iorb)), nbnew(iorb), ls(iorb),     & 
       &      mnew(iorb), iorb,                                 &     
       &    (taunew(ipol,iorb),ipol=1,3), iorb=1, norb )          

  WRITE( stdout, 301) ef*rytoev, energy0, denergy, nenergy, ecut2d, &
                      ewind, epsproj 
301   format (/,5x,                                         &
        &      'Fermi energy          = ',1pe15.5,/,5x,     &
        &      'energy0               = ',1pe15.1,/,5x,     &
        &      'denergy               = ',1pe15.1,/,5x,     & 
        &      'nenergy               = ',i10,/,5x,         &
        &      'ecut2d                = ',1pe15.1,/,5x,     &
        &      'ewind                 = ',1pe15.1,/,5x,     &
        &      'epsproj               = ',1pe15.1,/,5x) 

  if(nspin.ne.1) then
    WRITE( stdout,'(/,9x, ''Calculations are done for: '')')
    WRITE( stdout,'(2x, ''spin index = '', i6)') iofspin 
  endif

  if(norb.le.80) then
    WRITE( stdout,'(4x,''k slab'',3x,'' z(k)  z(k+1)'',             &
          &   5x,''crossing(iorb=1,norb)'')')                     
    do k=1, nrz
       WRITE( stdout,'(2x,i3,2x,3f7.4,3x,80i1)')                   &
         k,z(k),z(k+1),z(k+1)-z(k),(cross(iorb,k),iorb=1,norb)
    enddo  
  endif

  return
end subroutine init_cond
