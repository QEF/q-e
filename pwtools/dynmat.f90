!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
      program dynmat
!--------------------------------------------------------------------
!
!  read a dynamical matrix at q=0 , diagonalise it, calculate IR
!  Input data (namelist "input")
!
!  flmat  character  input file containing the dynamical matrix
!                    (default: flmat='dynmat')
!  q(3)      real    calculate LO modes (add nonanalytic terms) along
!                    the direction q (default: q=(0,0,0) )
!  amass(nt) real    mass for atom type nt, a.u.
!                    (default: amass is read from file flmat)
!  asr    logical    impose Acoustic Sum Rule (default:asr=.true.)
!  flout  character  output file containing frequencies and modes
!                    (default: flle='dynout')
!  flmol  character  as above, in a format suitable for 'molden'
!                    (default: flmol='moldout') 
!
      implicit none
      integer nax
      parameter (nax=30)
      character(len=50) flmat, flout, flmol
      character(len=3) atm(nax)
      logical asr, lread, gamma
      complex(kind=8) dyn(3,3,nax,nax), z(3*nax,3*nax)
      real(kind=8) tau(3,nax), amass(nax), amass_(nax), zstar(3,3,nax),&
           eps0(3,3), a0, omega, amconv, q(3), q_(3), w2(3*nax)
      integer ityp(nax), itau(nax), nat, na, nt, ntyp, nu, iout
      namelist /input/ amass,asr,flmat,flout,flmol,q
!
!
      asr  =.true.
      flout=' '
      flmat='dynmat'
      flout='dynout'
      flmol='moldout'
      amass(:)=0.0
      q(1)=0.0
      q(2)=0.0
      q(3)=0.0
!
      read (5,input)
!
      inquire(file=flmat,exist=lread)
      if (lread) then
         write(6,'(/5x,a,a)') 'Reading Dynamical Matrix from file ',&
              & flmat
      else
         write(6,'(/5x,a)') 'file not found', flmat
         stop
      end if
!
      call readmat (flmat,asr,nax,nat,ntyp,ityp,atm,a0,omega, amass_&
           &,tau,zstar,eps0,dyn,q_)
!
      gamma = abs(q_(1)**2+q_(2)**2+q_(3)**2).lt.1.0e-8
      amconv = 1.66042e-24/9.1095e-28*0.5
      do nt=1, ntyp
         if (amass(nt).gt.0.0) then
            amass(nt)=amass(nt)*amconv
         else 
            amass(nt)=amass_(nt)
         end if
      end do
!
      if (gamma) then
         do na=1,nat
            itau(na)=na
         end do
         call nonanal (nax,nat,dyn,q,itau,nax,eps0,zstar,omega)
      end if
!
      call dyndiag(nax,nat,amass,ityp,dyn,w2,z)
!
      if (flout.eq.' ') then
         iout=6
      else
         iout=4
         open (unit=iout,file=flout,status='unknown',form='formatted')
      end if
      call writemodes(nax,nat,q_,w2,z,iout)
      if(iout .ne. 6) close(unit=iout)
!
      call writemolden(nax,nat,atm,a0,tau,ityp,w2,z,flmol)
!
      if (gamma) call writeIR (nax, nat, w2, z, zstar)
!
      stop
      end
!
!-----------------------------------------------------------------------
      subroutine readmat (flmat,asr,nax,nat,ntyp,ityp,atm,a0,           &
     &                    omega,amass,tau,zstar,eps0,dynr,q)
!-----------------------------------------------------------------------
!
      implicit none
      character(len=*) flmat
      integer nax, nat, ntyp, ityp(nax)
      character(len=3) atm(ntyp)
      real(kind=8) amass(ntyp), tau(3,nax), a0, omega
      real(kind=8) dynr(2,3,3,nax,nax), eps0(3,3), zstar(3,3,nax),q(3)
      logical asr
!
      character(len=80) line
      real(kind=8)  at(3,3), celldm(6), sum
      integer ibrav, nt, na, nb, naa, nbb, i, j
      logical qfinito
!
!
      open (unit=1,file=flmat,status='old',form='formatted')
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,*) ntyp,nat,ibrav,celldm
      if (nat.gt.nax) stop ' too many atoms '
      a0=celldm(1)
      call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3))
      call volume(a0,at(1,1),at(1,2),at(1,3),omega)
      do nt=1,ntyp
         read(1,*) i,atm(nt),amass(nt)
      end do
      do na=1,nat
         read(1,*) i,ityp(na), (tau(j,na),j=1,3)
      end do
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,'(a)') line
      read(line(11:80),*) (q(i),i=1,3)
      qfinito=q(1).ne.0.0 .or. q(2).ne.0.0 .or. q(3).ne.0.0 
      read(1,'(a)') line
      do na = 1,nat
         do nb = 1,nat
            read (1,*) naa, nbb
            if (na.ne.naa .or. nb.ne.nbb) then
               print *, 'na, nb, naa, nbb =', na,nb,naa,nbb
               stop ' sgrunt!'
            end if
            read (1,*) ((dynr(1,i,j,na,nb),                             &
     &           dynr(2,i,j,na,nb), j=1,3), i=1,3)
         end do
      end do
!
      if (.not.qfinito) then
         read(1,*)
         read(1,'(a)') line
         if (line(1:23).ne.'     Dielectric Tensor:') then
            do na=1,nat
               do j=1,3
                  do i=1,3
                     zstar(i,j,na)=0.0
                  end do 
               end do
            end do
            do j=1,3
               do i=1,3
                  eps0(i,j)=0.0
               end do
               eps0(j,j)=1.0
            end do
         else
            read(1,*)
            read(1,*) ((eps0(i,j), j=1,3), i=1,3)
            read(1,*)
            read(1,*)
            read(1,*)
            do na = 1,nat
               read(1,*)
               read(1,*) ((zstar(i,j,na), j=1,3),i=1,3)
            end do
         end if
      end if
      if (asr) then
!
! ASR on effective charges
!
            do i=1,3
               do j=1,3
                  sum=0.0
                  do na=1,nat
                     sum = sum + zstar(i,j,na)
                  end do
                  do na=1,nat
                     zstar(i,j,na) = zstar(i,j,na) - sum/nat
                  end do
               end do
            end do
!
! ASR on dynamical matrix
!
            do i=1,3
               do j=1,3
                  do na=1,nat
                     sum=0.0
                     do nb=1,nat
                        if (na.ne.nb)                                   &
     &                       sum=sum+dynr(1,i,j,na,nb)
                     end do
                     dynr(1,i,j,na,na) = -sum
                  end do
               end do
            end do
      end if
!
      close(unit=1)
!
      return
      end
!
!-----------------------------------------------------------------------
subroutine writeIR (nax, nat, w2, z, zstar)
  !-----------------------------------------------------------------------
  !
  !   write IR cross sections
  !   on input: z = eigendisplacements
  !
 implicit none
 ! input
 integer nax, nat
 real(kind=8) w2(3*nat), zstar(3,3,nat)
 complex(kind=8) z(3*nax,3*nat)
 ! local
 integer na, nu, ipol, jpol
 real(kind=8) :: infrared(3*nat)
 real(kind=8) :: polar(3), rydcm1, cm1thz, freq, irmax
 !
 !  conversion factors RYD=>THZ, RYD=>1/CM e 1/CM=>THZ
 !
 rydcm1 = 13.6058*8065.5
 cm1thz = 241.796/8065.5
 !
 irmax=0.d0
 do nu = 1,3*nat
    do ipol=1,3
       polar(ipol)=0.0
    end do
    do na=1,nat
       do ipol=1,3
          do jpol=1,3
             polar(ipol) = polar(ipol) +  &
                  zstar(ipol,jpol,na)*z((na-1)*3+jpol,nu)
          end do
       end do
    end do
    infrared(nu) = sqrt(polar(1)**2+polar(2)**2+polar(3)**2)
    irmax = max(irmax,infrared(nu))
 end do
 !
 write (6,'(/5x,''Max IR cross section: '',e12.4/)') irmax
 !
 do nu = 1,3*nat
    freq = sqrt(abs(w2(nu)))*rydcm1
    if (w2(nu).lt.0.0) freq = -freq
    write (6,9010) nu, freq, freq*cm1thz, infrared(nu)/irmax
 end do
 !
 return
 !
9010 format(5x,'omega(',i2,') =',f10.2,' [cm-1] = ',f12.4, ' [THz]',&
            5x,' IR = ',f12.4)
 !
end subroutine writeIR
