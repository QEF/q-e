!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program q2r
  !
#include "machine.h"
  !
  implicit none
  !
  integer, parameter ::  nax=16, nrx1=8, nrx2=8, nrx3=8
  real(kind=8), parameter :: eps=1.0d-5
  integer ::  nr1, nr2, nr3, nr(3)
  !
  character(len=20) :: crystal
  character(len=80) :: title
  character(len=80) :: filin,filj,filf,fild
  character(len=3)  :: atm(nax)
  !
  logical :: lq,lrigid,zasr, lrigid_save 
  integer :: m1, m2, m3, l1, l2, l3, i, j, j1, j2, na1, na2, ipol
  integer :: nat, nq, ntyp, iq, icar, nfile, nqtot, ifile, nqs
  integer :: na, nt
  !
  integer :: m(3), nc(nrx1,nrx2,nrx3),ibrav,ityp(nax)
  !
  real(kind=8) :: celldm(6), at(3,3), bg(3,3), tau(3,nax)
  real(kind=8) :: q(3,48),omega, xq, amass(nax), resi,sum
  real(kind=8) :: epsil(3,3),zeu(3,3,nax)
  !
  complex(kind=8) :: phiq(3,3,nax,nax,48)
  complex(kind=8) :: phid(nrx1,nrx2,nrx3,3,3,nax,nax)
  !
  namelist / input / nr1,nr2,nr3,fild, zasr
  !
  nr1=0
  nr2=0
  nr3=0
  !
  read (5,input)
  !
  ! check input
  !
  if (nr1 > nrx1) call errore ('q2r',' nr1 too big, increase nrx1',nrx1)
  if (nr2 > nrx2) call errore ('q2r',' nr2 too big, increase nrx2',nrx2)
  if (nr3 > nrx3) call errore ('q2r',' nr3 too big, increase nrx3',nrx3)
  if (nr1 < 1) call errore ('q2r',' nr1 wrong or missing',1)
  if (nr2 < 1) call errore ('q2r',' nr2 wrong or missing',1)
  if (nr3 < 1) call errore ('q2r',' nr3 wrong or missing',1)
  !
  !
  ! copy nrX -> nr(X)
  !
  nr(1) = nr1
  nr(2) = nr2
  nr(3) = nr3
  !
  ! D matrix (analytical part)
  !
  !
  nqtot = 0
  !
  do l1=1,nr1
     do l2=1,nr2
        do l3=1,nr3
           nc(l1,l2,l3)=0
        end do
     end do
  end do
  !
  ! Reciprocal space dyn.mat. read from file
  !
  read (5,*) nfile
  do ifile=1,nfile
     read(5,'(a)') filin
     write (6,*) ' reading dyn.mat. from file ',filin
     open(unit=1,file=filin,status='old',form='formatted')
     call read_file(nqs,q,phiq,nax,epsil,zeu,lrigid,  &
                    ntyp,nat,ibrav,celldm,atm,amass,ityp,tau)
     if (ifile.eq.1) then
        lrigid_save=lrigid
        call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3))
        call volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
        call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
        if (lrigid.and.zasr) then
           do i=1,3
              do j=1,3
                 sum=0.d0
                 do na=1,nat
                    sum=sum+zeu(i,j,na)
                 end do
                 do na=1,nat
                    zeu(i,j,na)=zeu(i,j,na)-sum/nat
                 end do
              end do
           end do
        end if
     end if
     if (lrigid.and..not.lrigid_save) call errore('main',            &
          &          'in this case Gamma must be the first file ',1)
     !
     write (6,*) ' nqs= ',nqs
     close(unit=1)
     do nq = 1,nqs
        write(6,'(a,3f12.8)') ' q= ',(q(i,nq),i=1,3)
        lq = .true.
        do ipol=1,3
           xq = 0.0
           do icar=1,3
              xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
           end do
           lq = lq .and. (abs(nint(xq) - xq) .lt. eps)
           iq = nint(xq)
           !
           m(ipol)= mod(iq,nr(ipol)) + 1
           if (m(ipol) .lt. 1) m(ipol) = m(ipol) + nr(ipol) 
        end do
        if (.not.lq) call errore('init','q not allowed',1)
        if(nc(m(1),m(2),m(3)).eq.0) then
           nc(m(1),m(2),m(3))=1
           if (lrigid) call rgd_blk (nax,nat,phiq(1,1,1,1,nq),q(1,nq), &
                                     tau,epsil,zeu,bg,omega,-1.d0)
           call trasl(phid,phiq,nq,nrx1,nrx2,nrx3,nat,m(1),m(2),m(3),nax)
           nqtot=nqtot+1
        else
           write (*,'(3i4)') (m(i),i=1,3)
           call errore('init',' nc already filled: wrong q grid or wrong nr',1)
        end if
     end do
  end do
  !
  ! Check grid dimension
  !
  if (nqtot .eq. nr1*nr2*nr3) then
     write (6,'(/5x,a,i4)') ' q-space grid ok, #points = ',nqtot
  else
     call errore('init',' missing q-point(s)!',1)
  end if
  !
  ! dyn.mat. FFT
  !
  do j1=1,3
     do j2=1,3
        do na1=1,nat
           do na2=1,nat
              call tolerant_cft3(phid(1,1,1,j1,j2,na1,na2), &
                   nr1,nr2,nr3,nrx1,nrx2,nrx3,1)
              call DSCAL(2*nrx1*nrx2*nrx3,1.d0/(nr1*nr2*nr3),       &
                         phid(1,1,1,j1,j2,na1,na2),1)
           end do
        end do
     end do
  end do
  !
  ! Real space force constants written to file (analytical part)
  !
  resi = 0
  open(unit=2,file=fild,status='unknown',form='formatted')
  write(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
  do nt = 1,ntyp
     write(2,*) nt," '",atm(nt),"' ",amass(nt)
  end do
  do na=1,nat
     write(2,'(2i5,3f15.7)') na,ityp(na),(tau(j,na),j=1,3)
  end do
  write (2,*) lrigid
  if (lrigid) then
     write(2,'(3f15.7)') ((epsil(i,j),j=1,3),i=1,3)
     do na=1,nat
        write(2,'(i5)') na
        write(2,'(3f15.7)') ((zeu(i,j,na),j=1,3),i=1,3)
     end do
  end if
  write (2,'(4i4)') nr1, nr2, nr3 
  do j1=1,3
     do j2=1,3
        do na1=1,nat
           do na2=1,nat
              do m1=1,nr1
                 do m2=1,nr2
                    do m3=1,nr3
                       resi = resi + dabs(dimag(phid(m1,m2,m3,j1,j2,na1,na2)))
                    end do
                 end do
              end do
              write (2,'(4i4)') j1,j2,na1,na2
              write (2,'(3i4,2x,1pe18.11)')   &
                   (((m1,m2,m3,real(phid(m1,m2,m3,j1,j2,na1,na2)), &
                   m1=1,nr1),m2=1,nr2),m3=1,nr3)
           end do
        end do
     end do
  end do
  close(2)
  write (6,"(/5x,' fft-check: imaginary sum = ',e12.7)") resi
  ! 
end program q2r
!
!-----------------------------------------------------------------------
subroutine read_file(nqs,xq,phi,nax,epsil,zeu,lrigid,             &
     &           ntyp,nat,ibrav,celldm,atm,amass,ityp,tau)
  !-----------------------------------------------------------------------
  !
  implicit none
  !
  ! I/O variables
  logical :: lrigid
  integer :: nqs, nax, ntyp, nat, ibrav, ityp(nax)
  real(kind=8) :: epsil(3,3),zeu(3,3,nax)
  real(kind=8) :: xq(3,48), celldm(6), amass(nax), tau(3,nax)
  complex(kind=8) :: phi(3,3,nax,nax,48)
  character(len=3) atm(nax)
  ! local variables
  integer :: ntyp1,nat1,ibrav1,ityp1
  integer :: i, j, na, nb, nt
  real(kind=8) :: tau1(3), amass1, celldm1(6),q2
  real(kind=8) :: phir(3),phii(3)
  complex(kind=8) dcmplx
  character(len=75) :: line
  character(len=3)  :: atm1
  logical :: first
  data first/.true./
  save first
  !
  read(1,*) 
  read(1,*) 
  if (first) then
     !
     ! read cell information from file
     !
     read(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     if (nat.gt.nax) call errore('read_f','nax too small',nat)
     if (ntyp.gt.nat) call errore('read_f','ntyp.gt.nat!!',ntyp)
     do nt = 1,ntyp
        read(1,*) i,atm(nt),amass(nt)
        if (i.ne.nt) call errore('read_f','wrong data read',nt)
     end do
     do na=1,nat
        read(1,*) i,ityp(na),(tau(j,na),j=1,3)
        if (i.ne.na) call errore('read_f','wrong data read',na)
     end do
     !
     first=.false.
     lrigid=.false.
     !
  else
     !
     ! check cell information with previous one
     !
     read(1,*) ntyp1,nat1,ibrav1,(celldm1(i),i=1,6)
     if (ntyp1.ne.ntyp) call errore('read_f','wrong ntyp',1)
     if (nat1.ne.nat) call errore('read_f','wrong nat',1)
     if (ibrav1.ne.ibrav) call errore('read_f','wrong ibrav',1)
     do i=1,6
        if(celldm1(i).ne.celldm(i)) call errore('read_f','wrong celldm',i)
     end do
     do nt = 1,ntyp
        read(1,*) i,atm1,amass1
        if (i.ne.nt) call errore('read_f','wrong data read',nt)
        if (atm1.ne.atm(nt)) call errore('read_f','wrong atm',nt)
        if (amass1.ne.amass(nt)) call errore('read_f','wrong amass',nt)
     end do
     do na=1,nat
        read(1,*) i,ityp1,(tau1(j),j=1,3)
        if (i.ne.na) call errore('read_f','wrong data read',na)
        if (ityp1.ne.ityp(na)) call errore('read_f','wrong ityp',na)
        if (tau1(1).ne.tau(1,na)) call errore('read_f','wrong tau1',na)
        if (tau1(2).ne.tau(2,na)) call errore('read_f','wrong tau2',na)
        if (tau1(3).ne.tau(3,na)) call errore('read_f','wrong tau3',na)
     end do
  end if
  !
  !
  nqs = 0
100 continue
  read(1,*)
  read(1,'(a)') line
  if (line(6:14).ne.'Dynamical') then
     if (nqs.eq.0) call errore('read',' stop with nqs=0 !!',1)
     q2 = xq(1,nqs)**2 + xq(2,nqs)**2 + xq(3,nqs)**2
     if (q2.ne.0.d0) return
     do while (line(6:15).ne.'Dielectric') 
        read(1,'(a)',err=200, end=200) line
     end do
     lrigid=.true.
     read(1,*) ((epsil(i,j),j=1,3),i=1,3)
     read(1,*)
     read(1,*)
     read(1,*)
     write (*,*) 'macroscopic fields =',lrigid
     write (*,'(3f10.5)') ((epsil(i,j),j=1,3),i=1,3)
     do na=1,nat
        read(1,*)
        read(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
        write (*,*) ' na= ', na
        write (*,'(3f10.5)') ((zeu(i,j,na),j=1,3),i=1,3)
     end do
     return
200  write (*,*) ' Dielectric Tensor not found'
     lrigid=.false.     
     return
  end if
  !
  nqs = nqs + 1
  read(1,*) 
  read(1,'(a)') line
  read(line(11:75),*) (xq(i,nqs),i=1,3)
  read(1,*) 
  !
  do na=1,nat
     do nb=1,nat
        read(1,*) i,j
        if (i.ne.na) call errore('read_f','wrong na read',na)
        if (j.ne.nb) call errore('read_f','wrong nb read',nb)
        do i=1,3
           read (1,*) (phir(j),phii(j),j=1,3)
           do j = 1,3
              phi(i,j,na,nb,nqs) = dcmplx(phir(j),phii(j))
           end do
        end do
     end do
  end do
  !
  go to 100
  !
end subroutine read_file
!
!---------------------------------------------------------------------
subroutine trasl (phi,phiq,nq,nrx1,nrx2,nrx3,nat,m1,m2,m3,nax)
  !---------------------------------------------------------------------
  !
  implicit none
  integer:: j1,j2, m1, m2, m3, nrx1, nrx2, nrx3, na1, na2, nat, nax, nq
  !
  complex(kind=8) :: phi(nrx1,nrx2,nrx3,3,3,nax,nax)
  complex(kind=8) :: phiq(3,3,nax,nax,48)
  !
  do j1=1,3
     do j2=1,3
        do na1=1,nat
           do na2=1,nat
              phi(m1,m2,m3,j1,j2,na1,na2) = &
                   0.5 * (      phiq(j1,j2,na1,na2,nq) +  &
                          conjg(phiq(j2,j1,na2,na1,nq)))
           end do
        end do
     end do
  end do
  !
  return 
end subroutine trasl

# if defined __AIX || defined __FFTW || defined __SGI
#  define __FFT_MODULE_DRV
# endif

!-------------------------------------------------------------------
subroutine tolerant_cft3(f,nr1,nr2,nr3,nrx1,nrx2,nrx3,iflg)
  !-------------------------------------------------------------------
  !
  !  cft3 called for vectors with arbitrary maximal dimensions
  !

#if defined __FFT_MODULE_DRV
  use fft_scalar, only : cfft3d
#endif

  implicit none
  integer :: nr1,nr2,nr3,nrx1,nrx2,nrx3,iflg
  complex(kind=8) :: f(nrx1,nrx2,nrx3)
#if defined __FFT_MODULE_DRV
  complex(kind=8) :: ftmp(nrx1*nrx2*nrx3)
#endif
  integer :: i0,i1,i2,i3
  !
  i0=0
  do i3=1,nr3
     do i2=1,nr2
        do i1=1,nr1
           i0 = i0 + 1
#if defined __FFT_MODULE_DRV
           ftmp(i0) = f(i1,i2,i3)
#else
           f(i0,1,1) = f(i1,i2,i3)
#endif
        enddo
     enddo
  enddo
  !
  if (nr1.ne.1 .or. nr2.ne.1 .or. nr3.ne.1)                         &
#if defined __FFT_MODULE_DRV
       &    call cfft3d(ftmp,nr1,nr2,nr3,nr1,nr2,nr3,1)
#else
       &    call cft_3(f,nr1,nr2,nr3,nr1,nr2,nr3,1,iflg)
#endif
  !
  i0=nr1*nr2*nr3
  do i3=nr3,1,-1
     do i2=nr2,1,-1
        do i1=nr1,1,-1
#if defined __FFT_MODULE_DRV
           f(i1,i2,i3) = ftmp(i0)
#else
           f(i1,i2,i3) = f(i0,1,1)
#endif
           i0 = i0 - 1
        enddo
     enddo
  enddo

  !
  return
end subroutine tolerant_cft3
