!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#if defined __NEW_PUNCH

!-----------------------------------------------------------------------
subroutine read_file
  !-----------------------------------------------------------------------
  !
  !     This routine allocates space for all quantities already computed
  !     in the pwscf program and reads them from the data file.
  !
  !
#include "machine.h"
  USE kinds, ONLY: dp
  USE basis, ONLY: nat, ntyp, tau, ityp, natomwfc
  USE cell_base, ONLY: tpiba2, bg
  USE force_mod, ONLY: force
  USE klist, ONLY: nkstot, nks, xk, wk
  USE lsda_mod, ONLY: lsda, nspin, current_spin, isk
  USE wvfct, ONLY: nbnd, nbndx, et, wg
  USE symme, ONLY: irt 
  USE ktetra, ONLY: tetra, ntetra 
  USE extfield, ONLY:  forcefield, tefield
  USE cellmd, ONLY: cell_factor, lmovecell
  USE gvect, ONLY: gg, ecutwfc, ngm, g, nr1, nr2, nr3, eigts1, eigts2, eigts3
  USE gsmooth, ONLY: ngms, nls, nrx1s, nr1s, nr2s, nr3s
  USE scf, ONLY: rho, vr
  USE vlocal, ONLY: strf
  use io_files, only: tmp_dir, prefix, iunpun
  use restart_module, only: readfile_new
#ifdef __PARA
  use para
#endif
  implicit none
  !
  integer, parameter :: nax =1000 ! an unlikely large number of atoms
  integer :: i, ik, ibnd, ios, ierr
  !
  real(kind=DP), allocatable :: et_g(:,:), wg_g(:,:)
  real(kind=DP) :: rdum(1,1)
  integer :: kunittmp
  !
  ! choose the fortran unit to attach to the file
  !
  iunpun = 4
  !
  !  a value of zero cause the parameter to be read from the ".save" file
  !
  kunittmp = 0

  !  here we read the variables that dimension the system
  !  in parallel execution, only root proc reads the file
  !  and then broadcasts the values to all other procs
  !
  call readfile_new( 'dim', iunpun, rdum, rdum, kunittmp, 0, 0, ierr )
  IF( ierr /= 0 ) THEN
    call errore ('read_file', 'problem reading file '// &
      &      trim(tmp_dir)//trim(prefix)//'.save', ierr)
  END IF
  !
#ifdef __PARA
  kunit = kunittmp
#endif
  !
  !  allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  if ( nat <= 0 .or. nat > nax ) &
       call errore ('read_file', 'wrong number of atoms', 1)
  !
  allocate( et_g(nbnd,  nkstot), wg_g(nbnd,  nkstot) )

  allocate(tau (3, nat) )
  allocate(ityp (nat) )
  allocate(force (3, nat) )
  if (tefield) allocate(forcefield (3, nat) )
  allocate (irt( 48, nat))    
  allocate (tetra(4, MAX(ntetra,1)))    
  !
  !     here we read all the variables defining the system
  !     in parallel execution, only root proc read the file
  !     and then broadcast the values to all ather procs
  !
  call readfile_new( 'nowave', iunpun, et_g, wg_g, kunittmp, 0, 0, ierr )
  IF( ierr /= 0 ) THEN
    call errore ('read_file', 'problem reading file '// &
      &      trim(tmp_dir)//trim(prefix)//'.save', ierr)
  END IF
  !
  !
#ifdef __PARA
  kunit = kunittmp
  ! parallel execution: distribute across pools k-points and
  ! related variables (not a smart implementation)
  nks = nkstot
  ! nks and nkstot are redefined by the following routine
  call divide_et_impera (xk, wk, isk, lsda, nkstot, nks)
#endif
  !
  !  check whether LSDA
  !
  if (lsda) then
     nspin = 2
  else
     nspin = 1
     current_spin = 1
  endif
  cell_factor = 1.d0
  lmovecell = .false.
  !
  !   allocate memory for G- and R-space fft arrays
  !
  call allocate_fft
  call ggen
  !
  !    allocate the potential
  !
  call allocate_locpot
  call allocate_nlpot
  !
  !    allocate wavefunctions and related quantities (including et and wg)
  !
  ! TEMP: dimension natomwfc required in projwave by variable becp
  !       this is not a good reason - eventually nbndx must be = nbnd
  !
  nbndx = max(nbnd,natomwfc)
  call allocate_wfc
  !
  et = et_g
  wg = wg_g
  !
  deallocate( et_g, wg_g )
  !
#ifdef __PARA
  call poolscatter (nbnd , nkstot, et, nks, et)
  call poolscatter (nbnd , nkstot, wg, nks, wg)
#endif
  !
  ! read the charge density
  !
  call io_pot ( - 1, trim(prefix)//'.rho', rho, nspin)
  !
  ! read the potential
  !
  call io_pot ( - 1, trim(prefix)//'.pot', vr, nspin)
  !
  ! re-calculate the local part of the pseudopotential vltot
  ! and the core correction charge (if any) - This is done here
  ! for compatibility with the previous version of read_file
  !
  call init_vloc
  call struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
       nr3, strf, eigts1, eigts2, eigts3)
  call setlocal
  call set_rhoc
  !
  return
end subroutine read_file

#else

!-----------------------------------------------------------------------
subroutine read_file
  !-----------------------------------------------------------------------
  !
  !     This routine allocates space for all quantities already computed
  !     in the pwscf program and reads them from the data file.
  !
  !
#include "machine.h"
  USE kinds, ONLY: dp
  USE basis, ONLY: nat, ntyp, tau, ityp, natomwfc
  USE cell_base, ONLY: tpiba2, bg
  USE force_mod, ONLY: lforce, force
  USE klist, ONLY: nkstot, nks, xk, wk
  USE lsda_mod, ONLY: lsda, nspin, current_spin, isk
  USE wvfct, ONLY: nbnd, nbndx, et, wg
  USE symme, ONLY: irt 
  USE ktetra, ONLY: ltetra, tetra, ntetra
  USE extfield, ONLY:  forcefield, tefield
  USE cellmd, ONLY: cell_factor, lmovecell
  USE gvect, ONLY: gg, ecutwfc, ngm, g, nr1, nr2, nr3, eigts1, eigts2, eigts3
  USE gsmooth, ONLY: ngms, nls, nrx1s, nr1s, nr2s, nr3s
  USE scf, ONLY: rho, vr
  USE vlocal, ONLY: strf
  use io_files, only: tmp_dir, prefix, iunpun
#ifdef __PARA
  use para
#endif
  implicit none
  !
  integer, parameter :: nax =1000 ! an unlikely large number of atoms
  integer :: i, ik, ibnd, ios
  !
  iunpun = 4
  open (unit = iunpun, file = trim(tmp_dir)//trim(prefix)//'.pun', &
       form = 'unformatted', status = 'old', iostat = ios)
  call errore ('read_file', 'problem reading file '// &
       &      trim(tmp_dir)//trim(prefix)//'.pun', ios)
  !
  !     here we read all the variables describing the system
  !     in parallel execution, all processors read the same file
  !
  call saveall (iunpun, - 1)
  !
  !  allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  if (nat.le.0.or.nat.gt.nax) &
       call errore ('read_file', 'wrong number of atoms', 1)
  allocate(tau (3, nat) )
  allocate(ityp (nat) )
  allocate(force (3, nat) )
  if (tefield) allocate(forcefield (3, nat) )
  allocate (irt( 48, nat))    
  if (ltetra) allocate (tetra(4, ntetra))    
  !
  read (iunpun, err = 100, iostat = ios) tau
  read (iunpun, err = 101, iostat = ios) ityp
  read (iunpun, err = 102, iostat = ios) irt
  if (lforce) read (iunpun, err = 103, iostat = ios) force
  if (ltetra) read (iunpun, err = 104, iostat = ios) tetra
  !
  read (iunpun, err = 105, iostat = ios) ( (xk (i, ik), i = 1, 3), &
       ik = 1, nkstot)
  read (iunpun, err = 105, iostat = ios) (  wk (ik), ik = 1, nkstot)
  read (iunpun, err = 105, iostat = ios) ( isk (ik), ik = 1, nkstot)
#ifdef __PARA
  read (iunpun) kunit
  ! parallel execution: distribute across pools k-points and
  ! related variables (not a smart implementation)
  nks = nkstot
  ! nks and nkstot are redefined by the following routine
  call divide_et_impera (xk, wk, isk, lsda, nkstot, nks)
#endif
  !
  !  check whether LSDA
  !
  if (lsda) then
     nspin = 2
  else
     nspin = 1
     current_spin = 1
  endif
  cell_factor = 1.d0
  lmovecell = .false.
  !
  !   allocate memory for G- and R-space fft arrays
  !
  call allocate_fft
  call ggen
  !
  !    allocate the potential
  !
  call allocate_locpot
  call allocate_nlpot
  !
  !    allocate wavefunctions and related quantities (including et and wg)
  !
  ! TEMP: dimension natomwfc required in projwave by variable becp
  !       this is not a good reason - eventually nbndx must be = nbnd
  !
  nbndx = max(nbnd,natomwfc)
  call allocate_wfc
  !
  read (iunpun, err = 106, iostat = ios) ( (et (ibnd, ik), ibnd = 1, &
       nbnd), ik = 1, nkstot)
  read (iunpun, err = 106, iostat = ios) ( (wg (ibnd, ik), ibnd = 1, &
       nbnd), ik = 1, nkstot)
  close (iunpun)
  !
#ifdef __PARA
  call poolscatter (nbnd , nkstot, et, nks, et)
  call poolscatter (nbnd , nkstot, wg, nks, wg)
#endif
  !
  ! read the charge density
  !
  call io_pot ( - 1, trim(prefix)//'.rho', rho, nspin)
  !
  ! read the potential
  !
  call io_pot ( - 1, trim(prefix)//'.pot', vr, nspin)
  !
  ! re-calculate the local part of the pseudopotential vltot
  ! and the core correction charge (if any) - This is done here
  ! for compatibility with the previous version of read_file
  !
  call init_vloc
  call struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
       nr3, strf, eigts1, eigts2, eigts3)
  call setlocal
  call set_rhoc
  !
  return
100 call errore ('read_file', 'reading tau', abs (ios) )
101 call errore ('read_file', 'reading ityp', abs (ios) )
102 call errore ('read_file', 'reading irt', abs (ios) )
103 call errore ('read_file', 'reading forces', abs (ios) )
104 call errore ('read_file', 'reading tetrahedra', abs (ios) )
105 call errore ('read_file', 'reading k-points', abs (ios) )
106 call errore ('read_file', 'reading eigenvalues', abs (ios) )
end subroutine read_file

#endif
