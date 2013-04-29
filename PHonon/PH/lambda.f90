!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Last edition: September 5, 2008
! Edition author:  Eyvaz Isaev
! Department of Theoretical Physics, Moscow State Institute of Steel and Alloys, Russia
! Department of Physics, Chemistry and Biophysics (IFM), Linkoping University, Sweden
! Materials Theory Group, Institute of Physics and Materials Science,  Uppsala University, Sweden
! Eyvaz.Isaev@fysik.uu.se, isaev@ifm.liu.se, eyvaz_isaev@yahoo.com

!
program elph

  ! read files 'filelph' produced by phonon (one for each q-point)
  ! sum over q-points to produce the electron-phonon coefficients:
  ! lambda (the one of BCS superconductivity) and alpha^2*F(omega)
  ! T_c using Allen-Dynes formula

  implicit none
  integer, parameter:: npk=200, nsigx=50, nmodex=100, nex=200
  integer :: nks, ios, iuelph, ngauss, ngauss1, ngaussq, nsig, nmodes
  integer :: ik, ng, mu, nu, i
  real(kind=8) :: q(3,npk), wk(npk), degauss(nsigx), w2(nmodex), &
       dosef(nsigx), ef(nsigx), lambdaq(nmodex,nsigx),  &
       lambda(nsigx), alpha2F(nex,nsigx), logavg
  real(kind=8) qread(3), dosef1, ef1, degauss1, gammaq, lambda2, &
       degaussq, emax, deltae, e, omega, sum
  character(len=80) :: filelph
  real(kind=8), external :: w0gauss

  ! INPUT from standard input:
  !    emax  degaussq  ngaussq
  !    nks
  !    q(1,1)    q(2,1)    q(3,1)    wk(1)
  !      ...       ...       ...      ...
  !    q(1,nks)  q(2,nks)  q(3,nks)  wk(nks)
  !    filelph(1)
  !     ...
  !    filelph(nks)
  !
  ! emax (THz)    : alpha2F is plotted from 0 to "emax" in "nex" steps
  ! degaussq (THz): gaussian smearing for sum over q
  !                 NB: not the same used in phonon !
  ! ngaussq       : 0 for simple gaussian, 1 for Methfessel-Paxton etc.
  ! nks           : number of q-points used in the sum
  ! q, wk         : q-points and weights
  ! filelph       : output files from phonon, one for each q-point
  !                 May contain "nsig" calculations done with different
  !                 broadenings for the sum over k - all of them are used
  !
  ! OUTPUT in xmgr-readable format: files 'lambda.dat' and 'alpha2F.dat'
  !

  real*8 mustar, omegalog(20), Tc, x

  read(5,*) emax, degaussq, ngaussq
  deltae=emax/(nex-1)
  read(5,*) nks
  if (nks.gt.npk) call errore('lambda','too many q-points',npk)
  sum=0.d0
  do ik=1,nks
     read(5,*) q(1,ik), q(2,ik), q(3,ik), wk(ik)
     sum = sum + wk(ik)
  end do
  do ik=1,nks
     wk(ik)=wk(ik)/sum
  end do

  iuelph=4
  do ik=1,nks
     read(5,'(a)') filelph
     call remove_comments_from_string(filelph)
     open(unit=iuelph,file=filelph,status='old',iostat=ios)
     read (iuelph,*) qread(1),qread(2),qread(3), nsig, nmodes
!     if ( (qread(1)-q(1,ik))**2 + &
!          (qread(2)-q(2,ik))**2 + &
!          (qread(3)-q(3,ik))**2 .gt. 1.d-6) &
!          call errore('lambda','inconsistent q read',ik)
     if (nsig.le.0.or.nsig.gt.nsigx) &
          call errore('lambda','wrong/too many gauss.broad.',nsigx)
     if (nmodes.le.0.or.nmodes.gt.nmodex) &
          call errore('lambda','wrong # or too many modes',nmodex)
     if (ik.eq.1) then
        do ng=1,nsig
           lambda(ng)=0.d0
           do i=1,nex
              alpha2F(i,ng)=0.d0
           end do
        end do
     end if
     ! read omega^2
     read(iuelph,*) (w2(nu),nu=1,nmodes)
     ! read data
     do ng=1,nsig
        read (iuelph,9000) degauss1, ngauss1
        if (ik.eq.1) then
           degauss(ng)=degauss1
           ngauss =ngauss1
        else
           if (degauss(ng).ne.degauss1.or.ngauss.ne.ngauss1) &
              call errore('lambda','inconsistent gauss.broad. read',ik)
        end if
        read (iuelph,9005) dosef1, ef1
        if (ik.eq.1) then
           dosef(ng)=dosef1
           ef(ng)=ef1
        else
           if (dosef(ng).ne.dosef1.or.ef(ng).ne.ef1) &
              call errore('lambda','inconsistent DOS(Ef) read',ik)
        end if
        do mu=1,nmodes
           read (iuelph,9010) nu, lambdaq(mu,ng), gammaq
           if (nu.ne.mu) call errore('lambda','wrong mode read',mu)
           ! sum over k-points
           lambda(ng) = lambda(ng) + wk(ik)*lambdaq(mu,ng)
           do i=1,nex
              e=(i-1)*deltae
              ! 1 Ry = 3289.828 THz
              omega=sqrt(w2(mu))*3289.828
              alpha2F(i,ng) = alpha2F(i,ng) + &
                   wk(ik) * lambdaq(mu,ng) * omega * 0.5d0 * &
                   w0gauss((e-omega)/degaussq,ngaussq)/degaussq
           end do
        end do
     end do
     close(unit=iuelph)

  end do

  open(unit=iuelph,file='lambda.dat',status='unknown',form='formatted')
  write(iuelph,9014)
  do ng=1,nsig
     ! lambda2 is used as a check
     ! logavg is the logarithmic average of omega used in McMillan's formula(?)
     lambda2=0.d0
     logavg =0.d0
     do i=2,nex
        e=(i-1)*deltae
        lambda2=lambda2 + alpha2F(i,ng)/e
        logavg =logavg + alpha2F(i,ng)*log(e)/e
     end do
     lambda2=lambda2*2.d0*deltae
     logavg =logavg*2.d0 *deltae
     ! 1 THz = 48 K
     logavg=exp(logavg/lambda2)*47.9924d0
     omegalog(ng)=logavg
     write(6,9015) lambda(ng), lambda2, logavg,dosef(ng),degauss(ng)
     write(iuelph,9016) &
          degauss(ng), lambda(ng), lambda2, logavg,dosef(ng)
  end do
  close(unit=iuelph)

 read(5,*) mustar

  write(6,'("lambda", 8x, "omega_log", 10x, "T_c")')
  do i =1, nsig
        x=lambda(i)
         Tc = omegalog(i)/1.2*exp(-1.04*(1+x)/(x-mustar*(1+0.62*x)))
  write(6,'(f10.5,5x,f9.3,10x,f9.3)')  lambda(i), omegalog(i),  Tc
  enddo

  open(unit=iuelph,file='alpha2F.dat',status='unknown', &
       form='formatted')
  write(iuelph,9020) (degauss(ng),ng=1,nsig)
  do i=1,nex
     e=(i-1)*deltae
     write(iuelph,9025) e,(alpha2F(i,ng),ng=1,nsig)
  end do
  close(unit=iuelph)

  stop
9000 format(26x,f7.3,12x,i4)
9005 format(10x,f10.6,32x,f10.6)
9010 format(12x,i5,2x,f8.4,9x,f8.2)
9014 format('# degauss   lambda    int alpha2F  <log w>     N(Ef)')
9015 format(5x,'lambda =',f9.6,' (',f10.6,')  <log w>=',f9.3,'K  ', &
            'N(Ef)=',f9.6,' at degauss=',f5.3)
9016 format(f7.3,2f12.6,f10.3,2f12.6)
9020 format('# E(THz)',10(f10.3))
9025 format(f8.4,10(f10.5))

end program elph

