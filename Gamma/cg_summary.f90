!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cg_summary
  !-----------------------------------------------------------------------
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%           summarize input data          %%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  use pwcom
  use cgcom
  !
  implicit none
  integer :: nu, mu, i,l, na, nt
  !
  write (6,'(/5x,a75)') title
  write (6,9010) crystal,alat,omega,nat,ecutwfc,gcutm,tr2_ph
  !
  write (6,9020) (i,celldm(i),i=1,6)
  write (6,9030) ngm,nr1,nr2,nr3,nks
  write (6,9040)
  write (6,9050) (na,atm(ityp(na)),amass(ityp(na))/amconv,          &
       & (tau(i,na),i=1,3),na=1,nat)
  do nt = 1,ntyp
     write (6,9060) nlc(nt), nnl(nt)
     write (6,9070) nt,psd(nt),zp(nt),lmax(nt),lloc(nt)
     write (6,9080)
     write (6,'(/5x,"core")')
     write (6,9090) (alpc(i,nt),i=1,2)
     write (6,9100) (cc(i,nt),i=1,2)
     do l = 0,lmax(nt)
        write (6,'(/5x,"l = ",i2)') l
        write (6,9090) (alps(i,l,nt),i=1,3)
        write (6,9100) (aps(i,l,nt),i=1,3)
        write (6,9110) (aps(i,l,nt),i=4,6)
     end do
  end do
  write (6,9115)
  do nt = 1,ntyp
     write (6,9116) atm(nt),zv(nt),psd(nt)
  end do
  write (6,                                                         &
       &'(//5x,"atomic displacements are normalized to unity"/)')
  if (nmodes.lt.3*nat) then
     write (6,                                                     &
          &    '(5x,"phonon polarizations are as follows:"/)')
     do nu = 1,nmodes
        write (6,'(" mode # ",i3)') nu
        write (6,'(3(" (",f6.3,2f7.3,") "))')                   &
             &     ( u(mu,nu), mu = 1,3*nat)
     end do
  end if
  !
  return
  !
9010 format (//5x,'crystal is ',a20                                    &
       &        //5x,'lattice parameter     = ',f12.4                     &
       &         /5x,'unit-cell volume      = ',f12.4                     &
       &         /5x,'number of atoms /cell = ',i12                       &
       &         /5x,'kinetic-energy cutoff = ',f12.4                     &
       &         /5x,'g-space cutoff (gcutm)= ',f12.4                     &
       &         /5x,'convergence threshold = ',1pe12.1/)
  !
9020 format(/ 2 ( 3x,3(2x,'celldm(',i1,')=',f11.7) / ))
9030 format (5x,'ngm =',i6,'  nr1 =',i5,'  nr2 =',i5,'  nr3 =',i5,     &
       &          '  nks  =',i5)
9040 format (/5x,'site no     atom    mass',27x,'tau')
9050 format (7x,i2,9x,a2,3x,f8.4,9x,3f11.7)
9060 format (/15x,'atomic pseudopotential parameters',                 &
       &       ':  nlc =',i4,' nnl =',i4/)
9070 format (/5x,'atom',i2,' is ',a2,'   zval =',f5.1,'   lmax=',i2,   &
       &   '   lloc=',i2)
9080 format (/14x,'i=',7x,'1',10x,'2',10x,'3')
9090 format (5x,'alpha =',4x,3g11.5)
9100 format (5x,'a(i)  =',4x,3g11.5)
9110 format (5x,'a(i+3)=',4x,3g11.5)
9115 format (/5x,'atomic species       valence     pseudopotential')
9116 format (5x,a6,9x,f10.2,8x,5 (a2,'(',f5.2,')'))
  !
end subroutine cg_summary
