!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
program dipole
  !-----------------------------------------------------------------------
  !
  !      DESCRIPTION of the INPUT: see file INPUT_CHDENS in pwdocs/
  !      (only for variables: nfile filepp weight x0)
  ! 
  !      Calculation of the dipole moment:
  !          to be used for an isolated molecule in a box;
  !          the molecule must be at the center of the box.
  !          The code computes the dipole on the Wigner-Seitz cell of
  !          the Bravais lattice. The 3d box must contain this cell 
  !          otherwise meaningless numbers are obtained
  !
  !      Calculation of Makov-Payne correction for charged supercells:
  !          - not thoroughly tested
  !          - works only for clusters embedded within a cubic supercell
  !          - the cluster (and the plotting box) MUST be CENTERED 
  !            around (0,0,0), otherwise meaningless results are printed
  !          - always check that the printed total charge is the right one
  !          - for impurities in bulk crystals the correction should work
  !            as well, but the Madelung constant of the considered lattice
  !            must be used and the correction has to be divided by the
  !             crystal dielectric constant.
  !      Ref.: G. Makov and M.C. Payne, PRB 51, 4014 (1995)
  !      (note that Eq. 15 has the wrong sign for the quadrupole term)
  !      Contributed by Giovanni Cantele
  !
#include "f_defs.h"
  USE io_global,  ONLY : stdout
  USE parameters, ONLY : ntypx
  USE constants,  ONLY :  pi, fpi
  USE cell_base
  USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp, tau, zv
  USE char
  USE lsda_mod,   ONLY: nspin
  USE gvect
  USE gsmooth
  USE scf, ONLY: rho
  USE wavefunctions_module,  ONLY: psic
  USE io_files, ONLY: nd_nmbr

  implicit none
  integer, parameter :: nfilemax = 7
  ! maximum number of files with charge

  integer :: nfile, plot_num

  real(kind=DP) :: x0 (3), weight (nfilemax), dipol(0:3), quadrupol

  character (len=256) :: filename (nfilemax)

  real(kind=DP) :: celldms (6), gcutmsa, duals, ecuts, zvs(ntypx), ats(3,3)
  real(kind=DP), allocatable :: taus (:,:), rhor(:)
  integer :: ibravs, nrx1sa, nrx2sa, nrx3sa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  integer, allocatable :: ityps (:)
  integer :: ios, ipol, ifile, na, ir, i, j
  character (len=3) :: atms(ntypx)
  character (len=256) :: filepp(nfilemax)

  namelist /input/  nfile, filepp, weight, x0
  !
  call start_postproc (nd_nmbr)
  !
  !   set the DEFAULT values
  !
  nfile         = 1
  filepp(1)   = ' '
  weight(1)     = 1.0d0
  x0(:)         = 0.d0
  !
  !    read and check input data
  !
  CALL input_from_file ( )
  !
  ! reading the namelist input
  !
  read (5, input, err = 200, iostat = ios)
200 call errore ('dipole', 'reading input namelist', abs (ios) )

  ! check for number of files
  if (nfile < 1 .or. nfile > nfilemax) &
       call errore ('dipole ', 'nfile is wrong ', 1)
  !
  ! Read the header and allocate objects
  !

  call read_io_header(filepp (1), title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num )
  !
  allocate(tau (3, nat))
  allocate(ityp(nat))
  allocate(rhor(nrx1*nrx2*nrx3))
  !
  alat = celldm (1)
  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2
  doublegrid = dual.gt.4.0d0
  if (doublegrid) then
     gcutms = 4.d0 * ecutwfc / tpiba2
  else
     gcutms = gcutm
  endif

  nspin = 1
  if (ibrav > 0) then
    call latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
    at = at / alat !  bring at in units of alat
  end if

  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  call volume (alat, at(1,1), at(1,2), at(1,3), omega)

  call set_fft_dim ( )

  allocate (rho (nrx1*nrx2*nrx3,1))
  !
  ! Read first file
  !
  call plot_io (filepp (1), title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
                plot_num, atm, ityp, zv, tau, rho(1,1), -1)
  !
  rhor (:) = weight (1) * rho (:,1)
  !
  ! Read following files (if any), verify consistency
  ! Note that only rho is read; all other quantities are discarded
  !
  do ifile = 2, nfile
     allocate  (taus( 3 , nat))    
     allocate  (ityps( nat))    
     !
     call plot_io (filepp (ifile), title, nrx1sa, nrx2sa, nrx3sa, &
          nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
          duals, ecuts, plot_num, atms, ityps, zvs, taus, rho(1,1), - 1)
     !
     deallocate (ityps)
     deallocate (taus)
     !
     if (nats.gt.nat) call errore ('dipole', 'wrong file order? ', 1)
     if (nrx1.ne.nrx1sa.or.nrx2.ne.nrx2sa) call &
          errore ('dipole', 'incompatible nrx1 or nrx2', 1)
     if (nr1.ne.nr1sa.or.nr2.ne.nr2sa.or.nr3.ne.nr3sa) call &
          errore ('dipole', 'incompatible nr1 or nr2 or nr3', 1)
     if (ibravs.ne.ibrav) call errore ('dipole', 'incompatible ibrav', 1)
     if (gcutmsa.ne.gcutm.or.duals.ne.dual.or.ecuts.ne.ecutwfc ) &
          call errore ('dipole', 'incompatible gcutm or dual or ecut', 1)
     do i = 1, 6
        if (abs( celldm (i)-celldms (i) ) .gt. 1.0e-7 ) call errore &
             ('dipole', 'incompatible celldm', 1)
     enddo
     !
     rhor (:) = rhor (:) + weight (ifile) * rho (:,1)
  enddo
  !
  call dipole_fast (celldm (1), at, nat, omega, tau, &
       nrx1, nrx2, nrx3, nr1, nr2, nr3, rhor, &
       x0, dipol, quadrupol)

  call write_dipol(dipol,quadrupol,tau,nat,alat,zv,ntyp,ityp, &
       ibrav)
 
  deallocate(rhor)
  deallocate(tau)
  deallocate(ityp)
  call stop_pp
end program dipole
!
!-----------------------------------------------------------------------
subroutine dipole_fast (alat, at, nat, omega, tau, nrx1, nrx2, nrx3, &
     nr1, nr2, nr3, rhor, x0, dipol, quadrupol)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  integer :: nat, nrx1, nrx2, nrx3, nr1, nr2, nr3
  real(kind=DP) :: alat, omega, tau (3, nat), at (3, 3), &
       rhor(nrx1,nrx2,nrx3), x0 (3), dipol(0:3), quadrupol

  real(kind=DP) :: r(3)
  integer :: i,j,k, i1,j1,k1, ipol

  dipol=0.d0
  quadrupol=0.d0
  do k1 = -nr3/2, nr3/2
     k = k1 + 1
     if (k < 1) k = k + nr3
     do j1 = -nr2/2, nr2/2
        j = j1 + 1
        if (j < 1) j = j + nr2
        do i1 = -nr1/2, nr1/2  
           i = i1 + 1
           if (i < 1) i = i + nr1
           dipol(0) = dipol(0) + rhor (i, j, k) 
           do ipol=1,3
!!! x0 ???
              r(ipol) = i1*at(ipol,1)/nr1 + &
                        j1*at(ipol,2)/nr2 + &
                        k1*at(ipol,3)/nr3
           enddo
           do ipol=1,3
              dipol(ipol)=dipol(ipol) + r(ipol)*rhor(i,j,k)
              quadrupol = quadrupol + rhor (i, j, k)*r(ipol)**2
           enddo
        enddo
     enddo
  enddo

  dipol(0) = dipol(0) * omega / (nr1*nr2*nr3)
  do ipol=1,3
     dipol(ipol)=dipol(ipol) * omega / (nr1*nr2*nr3) * alat
  enddo
  quadrupol = quadrupol * omega / (nr1*nr2*nr3) * alat**2 

  return
end subroutine dipole_fast
!
!------------------------------------------------------------
subroutine write_dipol(dipol_el,quadrupol_el,tau,nat,alat,zv, &
     ntyp,ityp,ibrav)
  !-----------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE kinds, only : dp
  USE constants,  ONLY :  pi, rytoev
  implicit none

  integer :: nat, ntyp, ityp(nat), ibrav
  real(kind=dp) :: dipol_el(0:3), quadrupol_el, tau(3,nat), zv(ntyp), alat

  real(kind=dp) :: debye, dipol_ion(3), quadrupol_ion, dipol(3), quadrupol
  real(kind=DP) :: rhotot, zvtot, corr1, corr2, qq, AA, BB

  integer :: na, ipol

  !  Note that the definition of the Madelung constant used here
  !  differs from the "traditional" one found in the literature. See
  !  Lento, Mozos, Nieminen, J. Phys.: Condens. Matter 14 (2002), 2637-2645

  real(kind=DP), parameter:: Madelung(3) = (/ 2.8373, 2.8883, 2.885/)

  !
  !   compute ion dipole moments
  !
  dipol_ion=0.d0
  quadrupol_ion=0.d0
  zvtot=0.d0
  do na=1,nat
     zvtot = zvtot+zv(ityp(na))
     do ipol=1,3
        dipol_ion(ipol)=dipol_ion(ipol)+zv(ityp(na))* tau(ipol,na)*alat
        quadrupol_ion  =quadrupol_ion  +zv(ityp(na))*(tau(ipol,na)*alat)**2
     enddo
  enddo
  !
  !   compute ionic+electronic total charge, dipole and quadrupole moments
  !
  qq = -dipol_el(0) + zvtot
  dipol(:) = -dipol_el(1:3) + dipol_ion(:)
  quadrupol = -quadrupol_el + quadrupol_ion

  !
  !  Makov-Payne correction, PRB 51, 43014 (1995)
  !  Note that Eq. 15 has the wrong sign for the quadrupole term 
  !
  corr1 = -Madelung(ibrav)/alat * qq**2
  AA = quadrupol
  BB = dipol(1)**2 + dipol(2)**2 + dipol(3)**2
  corr2 = (4.d0/3.d0*pi) * (qq*AA-BB) / alat**3
  if (abs(qq) < 1.d-3) then
      corr2 = 0.0d0
  endif

  !
  !  Charge inside the Wigner-Seitz cell
  !
  WRITE( stdout, '(/4x," Charge density inside the Wigner-Seitz cell:",&
       &3f14.8," el.")')  dipol_el(0)
  !
  !  print the electron dipole moment calculated by the plotting 3d routines
  !  A positive dipole goes from the - charge to the + charge.
  !
  WRITE( stdout, '(/4x,"Electrons dipole moments",3f14.8," a.u.")')  &
       (-dipol_el(ipol),ipol=1,3)
  !
  ! print the ionic and total dipole moment
  !
  WRITE( stdout, '(4x,"     Ions dipole moments",3f14.8," a.u.")') &
       (dipol_ion(ipol),ipol=1,3)
  WRITE( stdout,'(4x,"    Total dipole moments",3f14.8," a.u.")') &
       (-dipol(ipol),ipol=1,3)
  !
  !   Print the same information in Debye
  !
  debye=2.54176d0

  WRITE( stdout,'(/4x,"Electrons dipole moments",3f14.8," Debye")') &
       (-dipol_el(ipol)*debye,ipol=1,3)
  WRITE( stdout,'(4x,"     Ions dipole moments",3f14.8," Debye")') &
       (dipol_ion(ipol)*debye,ipol=1,3)
  WRITE( stdout,'(4x,"    Total dipole moments",3f14.8," Debye")') &
       (-dipol(ipol)*debye,ipol=1,3)
  !
  ! print the Makov-Payne correction
  !
  if ( ibrav < 1 .or. ibrav > 3 ) then
     WRITE(stdout,'(/4x,"Makov-Payne correction only for cubic lattices")')
     return
  endif

  WRITE( stdout, '(//8x,"*********    MAKOV-PAYNE CORRECTION    ********")')
  WRITE(stdout,'(/4x,"Makov-Payne correction with Madelung constant = ",f8.4)')&
       Madelung(ibrav)
  WRITE( stdout, '(/,"Warning: results are meaningless if the cluster is &
       & not centered within the 3d box.")')
  WRITE( stdout, '(/4x,"Electron charge: ",f15.8," el.")') -dipol_el(0)
  WRITE( stdout, '(4x,"   Ionic charge: ",f15.8," el.")') zvtot
  WRITE( stdout, '(4x,"   Total charge: ",f15.8," el.")') qq
  !
  ! print the electronic,ionic and total quadrupol moment
  !
  WRITE( stdout, '(/4x,"Electrons quadrupole moment",f20.8," a.u.")')  &
       -quadrupol_el
  WRITE( stdout, '(4x,"     Ions quadrupole moment",f20.8," a.u.")') &
       quadrupol_ion
  WRITE( stdout,'(4x,"    Total quadrupole moment",f20.8," a.u.")') &
       quadrupol
  !
  !  print the results
  !
  WRITE( stdout,'(/4x,"Makov-Payne correction ",f14.8," Ry = ",f6.3, &
       & " eV (1st order, 1/a0)")') -corr1, (-corr1)*rytoev
  WRITE( stdout,'(4x,"                       ",f14.8," Ry = ",f6.3, &
       & " eV (2nd order, 1/a0^3)")') -corr2, (-corr2)*rytoev
     WRITE( stdout,'(4x,"                       ",f14.8," Ry = ",f6.3, &
          " eV (total)")') -corr1-corr2, (-corr1-corr2)*rytoev

  return

end subroutine write_dipol
