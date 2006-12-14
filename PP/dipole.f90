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
  !      For isolated systems, e.g. molecules in a supercell
  !      The molecule should be centered around the origin. If not,
  !      use variable x0 to re-center the molecule: the results
  !      should not depend on the exact choice of the origin
  !      if the cell is sufficiently large
  !      Calculation of Makov-Payne correction for charged supercells
  !      as in : G. Makov and M.C. Payne, PRB 51, 4014 (1995)
  !      (note that Eq. 15 has the wrong sign for the quadrupole term)
  !      Contributed by Giovanni Cantele and Paolo Cazzato
  !
  !      DESCRIPTION of the INPUT:
  !
  !-&input           Namelist &input; contains
  !
  !      filepp      file containing the 3D charge (produced by pp.x)
  !                  (REQUIRED)
  !      x0(3)       center the box around point x0 (in alat units)
  !                  ( default: (0,0,0) )
  !
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

  real(DP) :: x0 (3), dipol(0:3), quadrupol
  real(DP), allocatable :: rhor(:)
  integer :: ios, plot_num
  character (len=256) :: filepp

  namelist /input/  filepp, x0
  !
  call start_postproc (nd_nmbr)
  !
  !   set the DEFAULT values
  !
  filepp   = ' '
  x0(:)    = 0.d0
  !
  !    read and check input data
  !
  CALL input_from_file ( )
  !
  ! reading the namelist input
  !
  read (5, input, err = 200, iostat = ios)
200 call errore ('dipole', 'reading input namelist', abs (ios) )

  !
  ! Read the header and allocate objects
  !
  call read_io_header(filepp, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num )
  !
  allocate(tau (3, nat))
  allocate(ityp(nat))
  allocate(rhor(nrx1*nrx2*nrx3))
  !
  call latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
  alat = celldm (1) ! define alat
  at = at / alat    ! bring at in units of alat

  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2
  doublegrid = dual.gt.4.0d0
  if (doublegrid) then
     gcutms = 4.d0 * ecutwfc / tpiba2
  else
     gcutms = gcutm
  endif
  !
  nspin = 1
  !
  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  call volume (alat, at(1,1), at(1,2), at(1,3), omega)
  !
  call set_fft_dim ( )
  !
  allocate (rho (nrx1*nrx2*nrx3,1))
  !
  ! Read file
  !
  call plot_io (filepp, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
                plot_num, atm, ityp, zv, tau, rho(1,1), -1)
  !
  rhor (:) = rho (:,1)
  !
  call dipole_fast (celldm (1), at, bg, omega, nat, tau, &
       nrx1, nrx2, nrx3, nr1, nr2, nr3, rhor, &
       x0, dipol, quadrupol)
  !
  call write_dipol(x0,dipol,quadrupol,tau,nat,alat,zv,ntyp,ityp, &
       ibrav)
  !
  deallocate(rhor)
  deallocate(tau)
  deallocate(ityp)
  !
  call stop_pp
  !
end program dipole
!
!-----------------------------------------------------------------------
subroutine dipole_fast (alat, at, bg, omega, nat, tau, nrx1, nrx2, nrx3, &
     nr1, nr2, nr3, rhor, x0, dipol, quadrupol)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  integer :: nat, nrx1, nrx2, nrx3, nr1, nr2, nr3
  real(DP) :: alat, omega, tau (3, nat), at (3, 3), bg (3, 3), &
       rhor(nrx1,nrx2,nrx3), x0 (3), dipol(0:3), quadrupol

  real(DP) :: r(3), w1, w2, w3
  integer :: i,j,k, i0,j0,k0, i1,j1,k1, nr1m,nr2m,nr3m, ipol

  dipol=0.d0
  quadrupol=0.d0

  ! find FFT grid point closer to the desired center, X0
  ! add 1 because r=0 corresponds to (1,1,1)

  i0 = mod( nint( (x0(1)*bg(1,1)+x0(2)*bg(2,1)+x0(3)*bg(3,1))*nr1), nr1) + 1
  j0 = mod( nint( (x0(1)*bg(1,2)+x0(2)*bg(2,2)+x0(3)*bg(3,2))*nr2), nr2) + 1
  k0 = mod( nint( (x0(1)*bg(1,3)+x0(2)*bg(2,3)+x0(3)*bg(3,3))*nr3), nr3) + 1
  !
  ! redefine X0 so that what is used in the ionic term is consistent with
  ! what is used for the dipole term
  !
  x0(1)=(i0-1)*at(1,1)/ nr1 +(j0-1)*at(1,2)/ nr2 +(k0-1)*at(1,3)/ nr3
  x0(2)=(i0-1)*at(2,1)/ nr1 +(j0-1)*at(2,2)/ nr2 +(k0-1)*at(2,3)/ nr3
  x0(3)=(i0-1)*at(3,1)/ nr1 +(j0-1)*at(3,2)/ nr2 +(k0-1)*at(3,3)/ nr3
  !
  WRITE( stdout,'(5x,"Origin translated by: ",3f8.4," (alat units)")') x0

  nr1m = nr1/2
  nr2m = nr2/2
  nr3m = nr3/2
  do k1 = -nr3m, nr3m
     k = mod (k1 + k0, nr3)
     if (k < 1) k = k + nr3
     !
     ! correct weight at boundaries : if nr is even, we are summing both
     ! points at -nr/2 and nr/2, and they are the same because of periodicity
     !
     IF ( 2*nr3m == nr3 .AND. ABS(k1) == nr3m ) THEN
        w3 = 0.5d0
     else
        w3 = 1.d0
     end if
     do j1 = -nr2m, nr2m
        j = mod (j1 + j0, nr2)
        if (j < 1) j = j + nr2
        !
        ! correct weight at boundaries
        !
        IF ( 2*nr2m == nr2 .AND. ABS(j1) == nr2m ) THEN
           w2 = 0.5d0 * w3
        else
           w2 = 1.d0 * w3
        end if
        do i1 = -nr1m, nr1m
           i = mod(i1 + i0, nr1)
           if (i < 1) i = i + nr1
           !
           ! correct weight at boundaries
           !
           IF ( 2*nr1m == nr1 .AND. ABS(i1) == nr1m ) THEN
              w1 = 0.5d0 * w2
           else
              w1 = 1.d0 * w2
           end if
           !
           ! dipol(0) = charge density
           !
           dipol(0) = dipol(0) + w1*rhor (i, j, k) 
           !
           ! r-vectors are centered around r=0
           ! it is the FFT cell that is centered around x0
           !
           do ipol=1,3
              r(ipol) = i1*at(ipol,1)/nr1 + &
                        j1*at(ipol,2)/nr2 + &
                        k1*at(ipol,3)/nr3
           enddo
           !
           do ipol=1,3
              dipol(ipol)=dipol(ipol) + r(ipol)*w1*rhor(i,j,k)
              quadrupol = quadrupol + r(ipol)**2*w1*rhor (i, j, k)
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
subroutine write_dipol(x0,dipol_el,quadrupol_el,tau,nat,alat,zv, &
     ntyp,ityp,ibrav)
  !-----------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE kinds, only : dp
  USE constants,  ONLY :  pi, rytoev
  implicit none

  integer :: nat, ntyp, ityp(nat), ibrav
  real(DP) :: dipol_el(0:3), quadrupol_el, tau(3,nat), zv(ntyp), alat
  real(DP) :: x0(3)
  real(DP) :: debye, dipol_ion(3), quadrupol_ion, dipol(3), quadrupol
  real(DP) :: rhotot, zvtot, corr1, corr2, qq, AA, BB

  integer :: na, ipol

  !  Note that the definition of the Madelung constant used here
  !  differs from the "traditional" one found in the literature. See
  !  Lento, Mozos, Nieminen, J. Phys.: Condens. Matter 14 (2002), 2637-2645

  real(DP), parameter:: Madelung(3) = (/ 2.8373d0, 2.8883d0, 2.885d0/)

  !
  !   compute ion dipole moments
  !
  dipol_ion=0.d0
  quadrupol_ion=0.d0
  zvtot=0.d0
  do na=1,nat
     zvtot = zvtot+zv(ityp(na))
     do ipol=1,3
        dipol_ion(ipol)=dipol_ion(ipol)+zv(ityp(na))* (tau(ipol,na)-x0(ipol))*alat
        quadrupol_ion  =quadrupol_ion  +zv(ityp(na))*((tau(ipol,na)-x0(ipol))*alat)**2
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
  WRITE( stdout, '(/5x,"Charge density inside the Wigner-Seitz cell:",&
       &3f14.8," el.")')  dipol_el(0)
  !
  debye=2.54176d0
  !
  !  A positive dipole goes from the - charge to the + charge.
  !
  WRITE( stdout, '( 5x,"Dipole moment:")')
  WRITE( stdout, '( 5x,"Elect",3f8.4," au,  ", 3f8.4," Debye")') &
       (-dipol_el(ipol),ipol=1,3), (-dipol_el(ipol)*debye, ipol=1,3)
  WRITE( stdout, '( 5x,"Ionic",3f8.4," au,  ", 3f8.4," Debye")') &
       ( dipol_ion(ipol),ipol=1,3), (dipol_ion(ipol)*debye,ipol=1,3)
  WRITE( stdout, '( 5x,"Total",3f8.4," au,  ", 3f8.4," Debye")') &
       ( dipol(ipol),ipol=1,3),     ( dipol(ipol)*debye,   ipol=1,3)
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
  WRITE( stdout, '(/4x,"Electron charge: ",f15.8," el.")') -dipol_el(0)
  WRITE( stdout, '( 4x,"   Ionic charge: ",f15.8," el.")') zvtot
  WRITE( stdout, '( 4x,"   Total charge: ",f15.8," el.")') qq
  !
  ! print the electronic,ionic and total quadrupol moment
  !
  WRITE( stdout, '(/4x,"Electrons quadrupole moment",f20.8," a.u.")')  &
       -quadrupol_el
  WRITE( stdout, '( 4x,"     Ions quadrupole moment",f20.8," a.u.")') &
       quadrupol_ion
  WRITE( stdout, '( 4x,"    Total quadrupole moment",f20.8," a.u.")') &
       quadrupol
  !
  !  print the results
  !
  WRITE( stdout,'(/4x,"Makov-Payne correction ",f14.8," Ry = ",f6.3, &
     &  " eV (1st order, 1/a0)")') -corr1, (-corr1)*rytoev
  WRITE( stdout,'(4x,"                       ",f14.8," Ry = ",f6.3, &
     &  " eV (2nd order, 1/a0^3)")') -corr2, (-corr2)*rytoev
  WRITE( stdout,'(4x,"                       ",f14.8," Ry = ",f6.3, &
     &  " eV (total)")') -corr1-corr2, (-corr1-corr2)*rytoev

  return

end subroutine write_dipol
