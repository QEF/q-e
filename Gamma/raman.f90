!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program cg_raman
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom
  use io_files
  use cgcom
  use mp,             ONLY : mp_end
  use global_version
  use para

  implicit none

  real(kind=DP), allocatable :: deps_dtau(:,:,:,:), dynout(:,:)
  real(kind=DP), allocatable :: w2(:)
  character(len=9) :: cdate, ctime, code = 'RAMAN'
  logical :: exst
  integer :: i
  external date_and_tim
  !
  call init_clocks(.true.)
  call start_clock('raman')
  call startup( nd_nmbr, code, version_number )
  !
  gamma_only = .true.
  call cg_readin
  !
  call cg_setup
  !
  !  calculate eps0, zstar, dynamical matrix for the unperturbed system
  !
  allocate  ( dynout(3*nat,3*nat))    
  allocate  ( zstar( 3, 3, nat))    
  allocate  ( w2( 3*nat))    
  !
  call cg_eps0dyn(w2,dynout)
  !
  if (raman) then
     if (first.eq.1.and.last.eq.nmodes) then
        !
        !  calculate d eps0/d tau with finite differences
        !
        allocate ( deps_dtau( 3, 3, 3, nat))    
        call cg_deps(deps_dtau)
        !
        !  calculate nonresonant raman intensities for all modes
        !
        if (trans) call raman_cs(dynout,deps_dtau)
     else
        !
        !  calculate nonresonant raman intensities for selected modes
        !
        call raman_cs2(w2,dynout)
     end if
  end if
  !
  call stop_clock('raman')
  call print_clock(' ')
  !
  !  close everything and stop
  !
  if (epsil) then
     iubar=1
     call seqopn (iubar,'filbar1','unformatted',exst)
     close(unit=iubar,status='delete')
     call seqopn (iubar,'filbar2','unformatted',exst)
     close(unit=iubar,status='delete')
     call seqopn (iubar,'filbar3','unformatted',exst)
     close(unit=iubar,status='delete')
     iudwf=10
     call seqopn (iudwf,'fildwx1','unformatted',exst)
     close(unit=iudwf,status='delete')
     call seqopn (iudwf,'fildwx2','unformatted',exst)
     close(unit=iudwf,status='delete')
     call seqopn (iudwf,'fildwx3','unformatted',exst)
     close(unit=iudwf,status='delete')
  end if
#ifdef __PARA
  !
  !  parallel case: delete only once !
  !
  if (me.eq.1) then
#endif
!!!      open (unit=iunres,file='restartph',status='unknown')
!!!      close(unit=iunres,status='delete')
!!!      open (unit=iunres,file='restart_e',status='unknown')
!!!      close(unit=iunres,status='delete')
!!!      open (unit=iunres,file='restart_d',status='unknown')
!!!      close(unit=iunres,status='delete')
#ifdef __PARA
  end if
  ! call mpi_finalize(i)
#endif
  call mp_end()
  stop
  !
9000 format (/5x,'Program ',a12,' starts ...',/5x,                     &
       &            'Today is ',a9,' at ',a9)
end program cg_raman
!
!-----------------------------------------------------------------------
subroutine cg_deps(deps_dtau)
  !-----------------------------------------------------------------------
  !
  !  calculate d eps0/d tau with finite differences
  !
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE io_files,      ONLY : iunres
  use pwcom
  use cgcom
#ifdef __PARA
  use para
#endif
  implicit none
  real(kind=DP) :: deps_dtau(3,3,3,nat)
  !
  real(kind=DP) :: delta4(4), coeff4(4), delta2(2), coeff2(2), &
       delta, coeff
  integer iudyn, nd, na, ipol, nd_, na_, ipol_, jpol, kpol
  data delta2/-1.d0, 1.d0/, coeff2/-0.5d0, 0.5d0/
  data delta4/-2.d0, -1.d0, 1.d0, 2.d0/
  data coeff4/ 0.08333333333333d0,-0.66666666666666d0,              &
       &       0.66666666666667d0,-0.08333333333337d0 /
  !
  call start_clock('cg_deps')
  !
  !  Read partial results (if any)
  !
  open (unit=iunres,file='restart_d',form='formatted',status='unknown')
  read(iunres,*,err=1,end=1) na_,ipol_,nd_
  read(iunres,*,err=1,end=1) deps_dtau
  close(unit=iunres)
  if (na_.le.na) then
     WRITE( stdout,'(5x,"Restarting from atom ",i2,",  pol ",i1,      &
          &        ", nd=",i1)') na_,ipol_,nd_
  else
     WRITE( stdout,'(5x,"Reading saved data")')
  end if
  go to 2
1 close(unit=iunres)
  na_  =1
  ipol_=1
  nd_  =1
  deps_dtau(:,:,:,:) = 0.d0
  WRITE( stdout,'(5x,"Starting over from the beginning")')
2 continue
  !
  do na=na_,nat
     do ipol=1,3
        if (na.eq.na_.and.ipol.lt.ipol_) go to 11
        do nd=1,nderiv
           !
           !  Skip results from previous run (if any)
           !
           if (na.eq.na_.and.ipol.eq.ipol_.and.nd.lt.nd_) go to 12
           ! choose type of derivative formula (2- or 4-point)
           if (nderiv.eq.2) then
              delta=delta2(nd)
              coeff=coeff2(nd)
           else
              delta=delta4(nd)
              coeff=coeff4(nd)
           end if
           !
           ! Displaced atomic positions (remember that tau are in units of a0)
           !
           tau(ipol,na) =  tau(ipol,na) + delta*deltatau/alat
           !
           !
           call cg_neweps
           !
           tau(ipol,na) =  tau(ipol,na) - delta*deltatau/alat
           !
           do kpol=1,3
              do jpol=1,3
                 deps_dtau(kpol,jpol,ipol,na) =  &
                      deps_dtau(kpol,jpol,ipol,na) + &
                      epsilon0(kpol,jpol)*coeff/deltatau
              end do
           end do
           !
           !  Save partial results
           !
#ifdef __PARA
           !
           !  parallel case: write only once !
           !
           if (me.eq.1) then
#endif
              open (unit=iunres,file='restart_d',form='formatted',     &
                   status='unknown')
              if (nd.ne.nderiv) then
                 write(iunres,*) na,ipol,nd+1
              else if(ipol.ne.3) then
                 write(iunres,*) na,ipol+1,1
              else
                 write(iunres,*) na+1,1,1
              end if
              write(iunres,*) deps_dtau
              close(unit=iunres)
#ifdef __PARA
           endif
#endif
12         continue
        end do
11      continue
     end do
  end do
  !
  iudyn = 20
  open(unit=iudyn,file=fildyn,form='formatted',status='old',position='append')
  WRITE( stdout,'(/5x, "Raman tensors (atomic)"/)')
  write (iudyn,'(/5x,"Raman: D eps_{alpha,beta}/D tau_{s,gamma}"/)')
  do na=1,nat
     do ipol=1,3
        WRITE( stdout,'(/5x,"D eps(i,j)",5x,3e14.6     &
             &    /5x,"----------  =  ",3e14.6   &
             &    /5x,"D tau(",i2,")_",i1,4x,3e14.6)')             &
             &  (( deps_dtau(kpol,jpol,ipol,na), jpol=1,3), kpol=1,2),&
             &     na,ipol, (deps_dtau(3,jpol,ipol,na), jpol=1,3)
        write (iudyn,'("atom # ",i4,"   pol. ",i2)') na, ipol
        write (iudyn,'(3e24.12)') &
             ( (deps_dtau(kpol,jpol,ipol,na), jpol=1,3), kpol=1,3)
     end do
  end do
  WRITE( stdout,*)
  close (unit=iudyn)
  !
  return
end subroutine cg_deps
!
!-----------------------------------------------------------------------
subroutine cg_eps0dyn(w2,dynout)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE io_files,      ONLY : iunres
  use pwcom
  use cgcom
#ifdef __PARA
  use para
#endif
  implicit none
  real(kind=DP) :: w2(3*nat), dynout(3*nat,3*nat)
  !
  integer :: na, i,j,  nt, iudyn, mode_done
  !
  !   calculate linear response to macroscopic fields
  !
  if (epsil) then
     !
     !   verify if already calculated
     !
     open (unit=iunres,file='restart_e',form='formatted', status='unknown')
     read(iunres,*,end=1,err=1) epsilon0
     read(iunres,*,end=1,err=1) zstar
     close(unit=iunres)
     go to 2
     !
1    close(unit=iunres)
     call macro
     call solve_e
     !
     !   calculate the dielectric tensor and effective charges
     !
     call dielec(.true.)
     call generate_effective_charges  (nat,nsym,s,irt,at,bg, &
          n_diff_sites,equiv_atoms,has_equivalent,zstar)
     !
     !   save on file results
     !
#ifdef __PARA
     if (me.eq.1) then
#endif
     open (unit=iunres,file='restart_e',form='formatted',status='unknown')
     write(iunres,*) epsilon0
     write(iunres,*) zstar
     close(unit=iunres)
#ifdef __PARA
     end if
#endif
     !
     WRITE( stdout,'(/5x,"estimated dielectric constants =",3f10.3,  &
       &        /37x,3f10.3/37x,3f10.3)') ((epsilon0(i,j),j=1,3),i=1,3)
     WRITE( stdout,'(/5x,"z*(",i2,")",3f10.3,/11x,3f10.3/11x,3f10.3)') &
             (na, ((zstar(i,j,na),j=1,3),i=1,3), na=1,nat)
  end if
  !
2 continue
  !
  !   calculate linear response to lattice distorsions
  !
  if (trans) then
     !
     !   verify if already calculated
     !
     open (unit=iunres,file='restartph',form='formatted', status='unknown')
     read (iunres,*,err=10,end=10) mode_done
     if (mode_done.eq.nmodes+1) then
        read(iunres,*,end=10,err=10) dynout
        read(iunres,*,end=10,err=10) w2
        close(unit=iunres)
        go to 20
     end if
     !
10   close(unit=iunres)
     call solve_ph
     !
     !   get the complete dynamical matrix from the irreducible part
     !
     if (nmodes.eq.3*nat) call generate_dynamical_matrix &
          (nat,nsym,s,irt,at,bg, n_diff_sites,equiv_atoms,has_equivalent,dyn)
     !
     !   impose asr on the dynamical matrix
     !
     if (asr) call set_asr(nat,nasr,dyn)
     !
     ! diagonalize the dynamical matrix
     !
     call dyndiar(dyn,3*nat,nmodes,u,nat,ityp,amass,w2,dynout)
     !
     !  find new equilibrium positions (in the harmonic approximation)
     !         if (lforce)
     !     &        call equilib(nat,tau,force,nmodes,w2,dyn,3*nat,dtau,alat)
  end if
#ifdef __PARA
  if (me.ne.1) go to 20
#endif
  if (trans) call writedyn
  open (unit=iunres,file='restartph',form='formatted',status='unknown')
  write(iunres,*) nmodes+1
  write(iunres,*) dynout
  write(iunres,*) w2
  close(unit=iunres)
20 continue
  !
  return
  !
end subroutine cg_eps0dyn
!
!-----------------------------------------------------------------------
subroutine cg_neweps
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE io_global,  ONLY :  stdout
  use pwcom
  use cgcom
  !
  integer :: i, j
  real(kind=DP) :: rhotot, dmxc
  !
  !  recalculate self-consistent potential etc
  !
  call newscf
  !
  !  new derivative of the xc potential
  !
  dmuxc(:) = 0.d0
  do i = 1,nrxx
     rhotot = rho(i,current_spin)+rho_core(i)
     if ( rhotot.gt. 1.d-30 ) dmuxc(i)= dmxc( rhotot)
     if ( rhotot.lt.-1.d-30 ) dmuxc(i)=-dmxc(-rhotot)
  end do
  !
  !  re-initialize data needed for gradient corrections
  !
  call cg_setupdgc
  !
  !   calculate linear response to macroscopic fields
  !
  call macro
  !
  call solve_e
  !
  call dielec(.false.)
  !
  WRITE( stdout,'(/5x,"displaced atomic positions :")')
  WRITE( stdout,'(5x,3f12.6)') ((tau(i,j),i=1,3),j=1,nat)
  !
  WRITE( stdout,'(/5x,"estimated dielectric constants =",3f10.3,       &
       &        /37x,3f10.3/37x,3f10.3)') ((epsilon0(i,j),j=1,3),i=1,3)
  WRITE( stdout,*)
  !
end subroutine cg_neweps
!
!-----------------------------------------------------------------------
subroutine newscf
  !-----------------------------------------------------------------------
  !
  use pwcom
  use funct
  USE io_files,      ONLY : iunigk, iunwfc, input_drho, output_drho
  USE control_flags, ONLY : restart, reduce_io, lscf, istep, iprint, &
                            order, david, max_cg_iter, isolve, tr2,  &
                            ethr, mixing_beta, nmix, niter
  !
  implicit none
  integer :: iter
  !
  call start_clock('PWSCF')
  !
  !  set all kind of stuff needed by self-consistent (re-)calculation
  !
  dft='Same as Before'
  restart  =.false.
  reduce_io=.true.
  lscf=.true.
  lda_plus_u=.false.
  doublegrid=.false.
  lmovecell=.false.
  qcutz=0.0
  istep=1
  iprint=10000
  order=1
  input_drho=' '
  output_drho=' '
  !
  !  since we use only Gamma we don't need symmetries
  !
  nsym=1
  !
  ! these must be tuned for fast convergence
  !
  david = 4
  nbndx = max (nbndx, david*nbnd)
  max_cg_iter=20
  isolve=0
  tr2 =1.d-8
  ethr=1.d-8
  mixing_beta=0.7
  nmix=4
  niter=50
  !
  call openfil
  !
  call hinit1
  call electrons
  !
  close(unit=iunwfc, status='keep')
  close(unit=iunigk, status='delete')
  !
  call stop_clock('PWSCF')
  !
  return
end subroutine newscf
!
!-----------------------------------------------------------------------
subroutine raman_cs(dynout,deps_dtau)
  !-----------------------------------------------------------------------
  !
  !  calculate Raman cross section
  !
#include "machine.h"
  USE io_global,  ONLY : stdout
  use pwcom
  use cgcom
  !
  real(kind=DP) :: dynout(3*nat,3*nat), deps_dtau(3,3,3,nat)
  !
  integer :: nu, na, ipol, jpol, lpol
  real(kind=DP), allocatable :: raman_activity(:,:,:)
  !
  !
  allocate  ( raman_activity( 3, 3, nmodes))    
  WRITE( stdout,'(/5x, "Raman tensor for mode nu : dX_{alpha,beta}/d nu"/)')
  do nu=1,nmodes
     !
     do jpol=1,3
        do ipol=1,3
           raman_activity(ipol,jpol,nu) = 0.0
           do na=1,nat
              do lpol=1,3
                 raman_activity(ipol,jpol,nu) = raman_activity(ipol,jpol,nu) +&
                      deps_dtau(ipol,jpol,lpol,na) * dynout((na-1)*3+lpol,nu)
              end do
           end do
        end do
     end do
     WRITE( stdout,'(i5,3x,3e14.6,2(/8x,3e14.6))') &
             nu,( ( raman_activity(ipol,jpol,nu),jpol=1,3), ipol=1,3)
  end do
  deallocate(raman_activity)
  !
  return
end subroutine raman_cs
!
!-----------------------------------------------------------------------
subroutine raman_cs2(w2,dynout)
  !-----------------------------------------------------------------------
  !
  !  calculate d eps0/d u  (u=phonon mode) with finite differences
  !
#include "machine.h"
  USE io_global,  ONLY :  stdout
  USE io_files,      ONLY : iunres
  use pwcom
  use cgcom
#ifdef __PARA
  use para
#endif
  implicit none
  real(kind=DP) :: dynout(3*nat,3*nat), w2(3*nat)
  !
  real(kind=DP), allocatable :: raman_activity(:,:,:), infrared(:)
  real(kind=DP) :: delta4(4), coeff4(4), delta2(2), coeff2(2), &
       delta, norm, coeff
  integer iudyn, nd, nu, nd_, nu_, na, ipol, jpol
  data delta2/-1.d0, 1.d0/, coeff2/-0.5d0, 0.5d0/
  data delta4/-2.d0, -1.d0, 1.d0, 2.d0/
  data coeff4/ 0.08333333333333d0,-0.66666666666666d0,              &
       &       0.66666666666667d0,-0.08333333333337d0 /
  real(kind=8):: polar(3), rydcm1, cm1thz, freq, r1fac, r2fac, irfac
  real(kind=8):: alpha, beta2
  !
  call start_clock('raman_cs2')
  !
  !  Read partial results (if any)
  !
  allocate ( raman_activity( 3, 3, last-first+1))    
  open (unit=iunres,file='restart_d',form='formatted',status='unknown')
  read(iunres,*,err=1,end=1) nu_,nd_
  read(iunres,*,err=1,end=1) raman_activity
  close(unit=iunres)
  if (nu_.le.nu) then
     WRITE( stdout,'(5x,"Restarting from mode ",i3,", nd=",i1)') &
          nu_,nd_
  else
     WRITE( stdout,'(5x,"Reading saved data")')
  end if
  go to 2
1 close(unit=iunres)
  nu_=1
  nd_=1
  raman_activity(:,:,:) = 0.d0
  WRITE( stdout,'(5x,"Starting over from the beginning")')
2 continue
  !
  do nu=first,last
     if (nu.lt.nu_) go to 11
     !
     ! eigendisplacements are normalized as <u|M|U>=1
     ! we want to add a normalized eigendisplacement instead: <u|u>=1
     !
     norm = 0
     do na=1,nat
        do ipol=1,3
           norm = norm + dynout(3*(na-1)+ipol,nu)**2
        end do
     end do
     norm = sqrt(norm)
     !
     do nd=1,nderiv
        !
        !  Skip results from previous run (if any)
        !
        if (nu.eq.nu_.and.nd.lt.nd_) go to 12
        ! choose type of derivative formula (2- or 4-point)
        if (nderiv.eq.2) then
           delta=delta2(nd)
           coeff=coeff2(nd)
        else
           delta=delta4(nd)
           coeff=coeff4(nd)
        end if
        !
        ! Displaced atomic positions (remember that tau are in units of a0)
        !
        do na=1,nat
           do ipol=1,3
              tau(ipol,na) =  tau(ipol,na) + delta * deltatau * &
                   dynout(3*(na-1)+ipol,nu) / norm / alat
           end do
        end do
        !
        call cg_neweps
        !
        ! reset atomic positions to equilibrium value
        !
        do na=1,nat
           do ipol=1,3
              tau(ipol,na) =  tau(ipol,na) - delta * deltatau * &
                   dynout(3*(na-1)+ipol,nu) / norm / alat
           end do
        end do
        !
        ! calculate derivative, multiply by norm in order to have
        ! raman cross section for a mode normalized as <u|M|u>=1
        !
        do ipol=1,3
           do jpol=1,3
              raman_activity(ipol,jpol,nu-first+1) =  &
                   raman_activity(ipol,jpol,nu-first+1) + &
                   epsilon0(ipol,jpol)*coeff/deltatau * norm
           end do
        end do
        !
        !  Save partial results
        !
#ifdef __PARA
        !
        !  parallel case: write only once !
        !
        if (me.eq.1) then
#endif
           open (unit=iunres,file='restart_d',form='formatted',     &
                status='unknown')
           if (nd.ne.nderiv) then
              write(iunres,*) nu,nd+1
           else
              write(iunres,*) nu+1,1
           end if
           write(iunres,*) raman_activity
           close(unit=iunres)
#ifdef __PARA
        endif
#endif
12      continue
     end do
11   continue
  end do
  !
  do nu=first,last
     WRITE( stdout,'(i5,3x,3e14.6,2(/8x,3e14.6))') &
          nu,( ( raman_activity(ipol,jpol,nu-first+1),jpol=1,3), ipol=1,3)
  end do
  !
  !  conversion factors RYD=>THZ, RYD=>1/CM e 1/CM=>THZ
  !
  rydcm1 = 13.6058*8065.5
  cm1thz = 241.796/8065.5
  !
  !   conversion factor for IR cross sections from
  !   (Ry atomic units * e^2)  to  (Debye/A)^2/amu
  !   1 Ry mass unit = 2 * mass of one electron = 2 amu
  !   1 e = 4.80324x10^(-10) esu = 4.80324 Debye/A
  !     (1 Debye = 10^(-18) esu*cm = 0.2081928 e*A)
  !
  irfac = 4.80324**2/2.d0
  !
  !   conversion factor for Raman cross sections from Ry atomic units
  !   to A^4/amu
  !
  r1fac = 0.529177**4/2.d0
  !
  !  derivatives of epsilon are translated into derivatives of molecular
  !  polarizabilities by assuming a Clausius-Mossotti behavior
  !  (for anisotropic system epsilon is replaced by its trace)
  !
  r2fac = 3.d0*omega/(4.d0*3.1415926) *  &
       3.d0 / ( 2.d0 + (epsilon0(1,1) + epsilon0(2,2) + epsilon0(3,3))/3.d0 )
  !
  allocate (infrared(3*nat))    
  !
  do nu = 1,3*nat
     do ipol=1,3
        polar(ipol)=0.0
     end do
     do na=1,nat
        do ipol=1,3
           do jpol=1,3
              polar(ipol) = polar(ipol) +  &
                   zstar(ipol,jpol,na)*dynout((na-1)*3+jpol,nu)
           end do
        end do
     end do
     !
     ! infrared is in units of e^2/(Ry mass unit)
     !
     infrared(nu) = polar(1)**2+polar(2)**2+polar(3)**2
     !
  end do
  !
  WRITE( stdout,'(/5x,"IR cross sections are in (D/A)^2/amu units")')
  WRITE( stdout,'(/5x,"(multiplied by 1000)")')
  WRITE( stdout,'(5x,"Raman cross sections are in A^4/amu units")')
  WRITE( stdout,'(/"#  mode   [cm-1]     [THz]       IR      Raman")')
  !
  do nu = 1,3*nat
     !
     freq = sqrt(abs(w2(nu)))*rydcm1
     if (w2(nu).lt.0.0) freq = -freq
     !
     ! alpha, beta2: see PRB 54, 7830 (1996) and refs quoted therein
     !
     if( nu >= first .and. nu<= last ) then
        nu_ = nu-first+1
        alpha = (raman_activity(1,1,nu_) + &
                 raman_activity(2,2,nu_) + &
                 raman_activity(3,3,nu_))/3.d0
        beta2 = ( (raman_activity(1,1,nu_) - raman_activity(2,2,nu_))**2 + &
                  (raman_activity(1,1,nu_) - raman_activity(3,3,nu_))**2 + &
                  (raman_activity(2,2,nu_) - raman_activity(3,3,nu_))**2 + &
                  6.d0 * (raman_activity(1,2,nu_)**2 + &
                          raman_activity(1,3,nu_)**2 + &
                          raman_activity(2,3,nu_)**2) )/2.d0
     else
        alpha = 0
        beta2 = 0
     end if
     WRITE( stdout,'(i5,f10.2,f12.4,2f10.4)') &
          nu, freq, freq*cm1thz, infrared(nu)*irfac*1000, &
          (45.d0*alpha**2 + 7.0d0*beta2)*r1fac*r2fac
  end do
  !
  deallocate (infrared)
  deallocate (raman_activity)
  return
end subroutine raman_cs2

