!
! Copyright (C) 2001-2004 PWSCF group
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
!  and Raman cross sections (if Z* and Raman tensor available) 
!
!  Input data (namelist "input")
!
!  fildyn  character input file containing the dynamical matrix
!                    (default: fildyn='matdyn')
!  q(3)      real    calculate LO modes (add nonanalytic terms) along
!                    the direction q (default: q=(0,0,0) )
!  amass(nt) real    mass for atom type nt, a.u.
!                    (default: amass is read from file fildyn)
!  asr     logical   impose Acoustic Sum Rule (default:asr=.true.)
!  filout  character output file containing frequencies and modes
!                    (default: filout='dynmat.out')
!  filmol  character as above, in a format suitable for 'molden'
!                    (default: filmol='molden.out') 
!
      implicit none
      integer, parameter :: nax=30
      character(len=50):: fildyn, filout, filmol
      character(len=3) :: atm(nax)
      logical :: asr, lread, gamma
      complex(kind=8) :: dyn(3,3,nax,nax), z(3*nax,3*nax)
      real(kind=8) :: tau(3,nax), amass(nax), amass_(nax), zstar(3,3,nax),&
           eps0(3,3), a0, omega, amconv, q(3), q_(3), w2(3*nax)
      real(kind=8) :: dchi_dtau(3,3,3,nax)
      integer :: ityp(nax), itau(nax), nat, na, nt, ntyp, nu, iout
      namelist /input/ amass, asr, fildyn, filout, filmol, q
!
!
      asr  =.true.
      fildyn='matdyn'
      filout='dynmat.out'
      filmol='molden.out'
      amass(:)=0.0
      q(1)=0.0
      q(2)=0.0
      q(3)=0.0
!
      read (5,input)
!
      inquire(file=fildyn,exist=lread)
      if (lread) then
         write(6,'(/5x,a,a)') 'Reading Dynamical Matrix from file ',&
              & fildyn
      else
         write(6,'(/5x,a)') 'file not found', fildyn
         stop
      end if
!
      call readmat (fildyn,asr,nax,nat,ntyp,ityp,atm,a0,omega, amass_&
           &,tau,zstar,eps0,dyn,dchi_dtau,q_)
!
      gamma = abs(q_(1)**2+q_(2)**2+q_(3)**2).lt.1.0e-8
      amconv = 1.66042e-24/9.1095e-28*0.5
      do nt=1, ntyp
         if (amass(nt) > 0.0) then
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
      if (filout.eq.' ') then
         iout=6
      else
         iout=4
         open (unit=iout,file=filout,status='unknown',form='formatted')
      end if
      call writemodes(nax,nat,q_,w2,z,iout)
      if(iout .ne. 6) close(unit=iout)
!
      call writemolden(nax,nat,atm,a0,tau,ityp,w2,z,filmol)
!
      if (gamma) call RamanIR &
           (nax, nat, omega, w2, z, zstar, eps0, dchi_dtau)
!
      stop
      end
!
!-----------------------------------------------------------------------
      subroutine readmat (fildyn,asr,nax,nat,ntyp,ityp,atm,a0,           &
           omega,amass,tau,zstar,eps0,dynr,dchi_dtau,q)
!-----------------------------------------------------------------------
!
      implicit none
      character(len=*) fildyn
      integer nax, nat, ntyp, ityp(nax)
      character(len=3) atm(ntyp)
      real(kind=8) amass(ntyp), tau(3,nax), a0, omega
      real(kind=8) dynr(2,3,3,nax,nax), eps0(3,3), zstar(3,3,nax), &
           dchi_dtau(3,3,3,nax), q(3)
      logical asr
!
      character(len=80) line
      real(kind=8)  at(3,3), celldm(6), sum
      integer ibrav, nt, na, nb, naa, nbb, i, j, k
      logical qfinito, noraman
!
!
      noraman=.true.
      open (unit=1,file=fildyn,status='old',form='formatted')
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,*) ntyp,nat,ibrav,celldm
      if (nat.gt.nax) stop ' too many atoms '
      a0=celldm(1)
      call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
      at = at / a0 !  bring at in units of alat
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
 20         read(1,'(a)',end=10,err=10) line
            if (line(1:17) == '     Raman tensor') go to 25
            go to 20
 25         read(1,*,end=10,err=10)
            do na = 1,nat
               do i = 1, 3
                  read(1,*,end=10,err=10)
                  read(1,*,end=10,err=10) &
                       ((dchi_dtau(k,j,i,na), j=1,3), k=1,3)
               end do
            end do
            write(6,'(/5x,a)') 'Raman cross sections read'
            noraman=.false.
10          continue
         end if
      end if
      if (noraman) dchi_dtau=0.d0
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
subroutine RamanIR (nax, nat, omega, w2, z, zstar, eps0, dchi_dtau)
  !-----------------------------------------------------------------------
  !
  !   write IR and Raman cross sections
  !   on input: z = eigendisplacements (normalized as <z|M|z>)
  !             zstar = effective charges (units of e)
  !             dchi_dtau = derivatives of chi wrt atomic displacement
  !                         (units: A^2)
 implicit none
 ! input
 integer nax, nat
 real(kind=8) omega, w2(3*nat), zstar(3,3,nat), eps0(3,3), &
      dchi_dtau(3,3,3,nat), chi(3,3)
 complex(kind=8) z(3*nax,3*nat)
 ! local
 integer na, nu, ipol, jpol, lpol
 logical noraman
 real(kind=8), pointer :: infrared(:), raman(:,:,:)
 real(kind=8):: polar(3), rydcm1, cm1thz, freq, r1fac, irfac
 real(kind=8):: cmfac, alpha, beta2
 !
 !  conversion factors Ry => THz, Ry=>cm^(-1) e cm^(-1)=>THz
 !
 rydcm1 = 13.6058*8065.5
 cm1thz = 241.796/8065.5
 !
 !   conversion factor from (Ry au for mass)^(-1) to amu(-1)
 !
 r1fac = 911.444
 !
 !   conversion factor for IR cross sections from Ry au to (Debye/A)^2/amu
 !
 irfac = 10514.0155
 !
 write (6,'(/5x,"Polarizability (A^3 units)")')
 !
 !  correction to molecular polarizabilities from Clausius-Mossotti formula
 !  (for anisoptropic systems the 
 !
 cmfac = 3.d0 / ( 2.d0 + (eps0(1,1) + eps0(2,2) + eps0(3,3))/3.d0 )
 !
 write (6,'(/5x,"Multiply by",f9.6," for Clausius-Mossotti correction")') cmfac
 do jpol=1,3
    do ipol=1,3
       if (ipol == jpol) then
          chi(ipol,jpol) = (eps0(ipol,jpol)-1.d0) 
       else
          chi(ipol,jpol) = eps0(ipol,jpol)
       end if
    end do
 end do
 do ipol=1,3
    write (6,'(5x,3f12.6)') (chi(ipol,jpol)*0.529177**3*omega/4.0/3.1415926, &
         jpol=1,3)
 end do
 !
 allocate(infrared (3*nat))
 allocate(raman(3,3,3*nat))
 !
 noraman=.true.
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
    !
    infrared(nu) = 2.d0*(polar(1)**2+polar(2)**2+polar(3)**2)*irfac
    !
    do ipol=1,3
       do jpol=1,3
          raman(ipol,jpol,nu)=0.0
          do na=1,nat
             do lpol=1,3
                raman(ipol,jpol,nu) = raman(ipol,jpol,nu) + &
                     dchi_dtau(ipol,jpol,lpol,na) * z((na-1)*3+lpol,nu) 
             end do
          end do
          noraman=noraman .and. abs(raman(ipol,jpol,nu)).lt.1.d-12
       end do
    end do
    !   Raman cross sections are in units of bohr^4/(Ry mass unit)
 end do
 !
 write (6,'(/5x,"IR cross sections are in (D/A)^2/amu units")')
 if (noraman) then
    write (6,'(/"#  mode   [cm-1]     [THz]       IR")')
 else
    write (6,'(5x,"Raman cross sections are in A^4/amu units")')
    write (6,'(/5x,"Multiply Raman by",f9.6," for Clausius-Mossotti", &
         & " correction")') cmfac**2
    write (6,'(/"#  mode   [cm-1]     [THz]      IR       Raman     depol")')
 end if
 !
 do nu = 1,3*nat
    !
    freq = sqrt(abs(w2(nu)))*rydcm1
    if (w2(nu).lt.0.0) freq = -freq
    !
    ! alpha, beta2: see PRB 54, 7830 (1996) and refs quoted therein
    !
    if (noraman) then
       write (6,'(i5,f10.2,f12.4,2f10.4)') &
         nu, freq, freq*cm1thz, infrared(nu)
    else
       alpha = (raman(1,1,nu) + raman(2,2,nu) + raman(3,3,nu))/3.d0
       beta2 = ( (raman(1,1,nu) - raman(2,2,nu))**2 + &
                 (raman(1,1,nu) - raman(3,3,nu))**2 + &
                 (raman(2,2,nu) - raman(3,3,nu))**2 + 6.d0 * &
          (raman(1,2,nu)**2 + raman(1,3,nu)**2 + raman(2,3,nu)**2) )/2.d0
       write (6,'(i5,f10.2,f12.4,3f10.4)') &
            nu, freq, freq*cm1thz, infrared(nu), &
            (45.d0*alpha**2 + 7.0d0*beta2)*r1fac, &
             3.d0*beta2/(45.d0*alpha**2 + 4.0d0*beta2)
    end if
 end do
 !
 return
 !
end subroutine RamanIR
