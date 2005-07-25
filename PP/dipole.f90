!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program dipole
  !-----------------------------------------------------------------------
  !
  !      DESCRIPTION of the INPUT: see file INPUT_CHDENS in pwdocs/
  !      (only for variables: nfile filepp weight e1,e2,e3 nx,ny,nz x0)
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

  integer :: ounit, ios, ipol, nfile, ifile, nx, ny, nz, &
       na, ir, i, j, ig, plot_num

  real(kind=DP) :: e1(3), e2(3), e3(3), x0 (3), radius, m1, m2, m3, &
       weight (nfilemax)

  character (len=256) :: filepol, filename (nfilemax)

  real(kind=DP) :: celldms (6), gcutmsa, duals, ecuts, zvs(ntypx), ats(3,3)
  real(kind=DP), allocatable :: taus (:,:), rhor(:)
  integer :: ibravs, nrx1sa, nrx2sa, nrx3sa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  integer, allocatable :: ityps (:)
  character (len=3) :: atms(ntypx)
  character (len=256) :: filepp(nfilemax)
  real(kind=DP) :: rhodum, dipol(0:3), quadrupol, rhotot
  complex(kind=DP), allocatable:: rhog (:)
  ! rho or polarization in G space
  logical :: fast3d

  namelist /input/  &
       nfile, filepp, weight, e1, e2, e3, nx, ny, nz, x0
  !
  call start_postproc (nd_nmbr)
  !
  !   set the DEFAULT values
  !
  nfile         = 1
  filepp(1)   = ' '
  weight(1)     = 1.0d0
  radius        = 1.0d0
  e1(:)         = 0.d0
  e2(:)         = 0.d0
  e3(:)         = 0.d0
  x0(:)         = 0.d0
  nx            = 0
  ny            = 0
  nz            = 0
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
  if (nfile.le.0.or.nfile.gt.nfilemax) &
       call errore ('dipole ', 'nfile is wrong ', 1)

  ! 3D plot : check variables

  if ( e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3) > 1d-6 .or. &
       e1(1)*e3(1) + e1(2)*e3(2) + e1(3)*e3(3) > 1d-6 .or. &
       e2(1)*e3(1) + e2(2)*e3(2) + e2(3)*e3(3) > 1d-6 )    &
       call errore ('dipole', 'e1, e2, e3 are not orthogonal', 1)

  !
  ! Read the header and allocate objects
  !

  call read_io_header(filepp (1), title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num )
  !
  ! ... see comment above
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
  if (ibrav.gt.0) then
    call latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
    at = at / alat !  bring at in units of alat
  end if

  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  call volume (alat, at(1,1), at(1,2), at(1,3), omega)

  call set_fft_dim ( )

  call allocate_fft ( )
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
  ounit = 6

  !
  !    At this point we start the calculations, first we normalize the 
  !    vectors defining the plotting region. 
  !    If these vectors have 0 length, replace them with crystal axis
  !

  m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  if (abs(m1) < 1.d-6) then
     e1 (:) = at(:,1)
     m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  end if
  e1 (:) = e1 (:) / m1
  !
  m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  if (abs(m2) < 1.d-6) then
     e2 (:) = at(:,2)
     m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  end if
  e2 (:) = e2 (:) / m2
  !
  m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  if (abs(m3) < 1.d-6) then
     e3 (:) = at(:,3)
     m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  end if
  e3 (:) = e3 (:) / m3
  !
  !    and rebuild G-vectors in reciprocal space
  !
  call ggen
  !
  !    here we compute the fourier component of the quantity to plot
  !
  psic(:) = DCMPLX (rhor(:), 0.d0)
  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  !    we store the fourier components in the array rhog
  !

  allocate (rhog( ngm))    
  do ig = 1, ngm
     rhog (ig) = psic (nl (ig) )
  enddo

  ! are vectors defining the plotting region aligned along xyz ?

  fast3d = ( e1(2) == 0.d0  .and.  e1(3) == 0.d0) .and. &
           ( e2(1) == 0.d0  .and.  e2(3) == 0.d0) .and. &
           ( e3(1) == 0.d0  .and.  e3(2) == 0.d0) 

  ! are crystal axis aligned along xyz ?

  fast3d = fast3d .and. &
       ( at(2,1) == 0.d0  .and.  at(3,1) == 0.d0) .and. &
       ( at(1,2) == 0.d0  .and.  at(3,2) == 0.d0) .and. &
       ( at(1,3) == 0.d0  .and.  at(2,3) == 0.d0) 

  !
  if (fast3d) then
     
     call dipole_fast (celldm (1), at, nat, tau, &
          nrx1, nrx2, nrx3, nr1, nr2, nr3, rhor, &
          bg, m1, m2, m3, x0, e1, e2, e3, dipol(0), quadrupol)
  else
     if (nx<=0 .or. ny <=0 .or. nz <=0) &
          call errore("dipole","nx,ny,nz, required",1)

     call dipole_3d (celldm (1), at, nat, tau, ngm, g, rhog, &
          nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, dipol(0), quadrupol )
  end if

  call write_dipol(dipol(0),quadrupol,tau,nat,alat,zv,ntyp,ityp, &
       ibrav, rhotot)
 
  deallocate(rhor)
  deallocate(rhog)
  deallocate(tau)
  deallocate(ityp)
  call stop_pp
end program dipole
!
!-----------------------------------------------------------------------
subroutine dipole_3d (alat, at, nat, tau, ngm, g, rhog, &
     nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, dipol, quadrupol )
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  use constants, only:  pi 
  implicit none
  integer :: nat, ngm, nx, ny, nz
  ! number of atoms, number of G vectors, number of points along x, y, z

  real(kind=DP) :: alat, tau(3,nat), at(3,3), g(3,ngm), x0(3), &
                   e1(3), e2(3), e3(3), m1, m2, m3, dipol(0:3), quadrupol
  ! lattice parameter
  ! atomic positions
  ! lattice vectors
  ! G-vectors
  ! origin
  ! vectors e1,e2,e3 defining the parallelepiped
  ! moduli of e1,e2,e3
  ! electronic dipole & quadrupole moments

  complex(kind=DP) :: rhog (ngm)
  ! rho or polarization in G space
  integer :: i, j, k, ig

  real(kind=DP) :: deltax, deltay, deltaz
  ! steps along e1, e2, e3
  complex(kind=DP), allocatable :: eigx (:), eigy (:), eigz (:)
  real(kind=DP), allocatable :: carica (:,:,:), rws(:,:)
  real(kind=dp) :: wsweight, r(3), omega, fact
  integer :: ipol, na, nrwsx, nrws

  allocate (eigx(  nx))    
  allocate (eigy(  ny))    
  allocate (eigz(  nz))    
  allocate (carica( nx , ny , nz))    

  deltax = m1 / nx 
  deltay = m2 / ny 
  deltaz = m3 / nz 

  carica = 0.d0
  do ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=exp(iG*e2), eigz=exp(iG*e3)
     ! These factors are calculated and stored in order to save CPU time
     !
     do i = 1, nx
        eigx (i) = exp( (0.d0,1.d0) * 2.d0 * pi * ( (i-1) * deltax * &
             (e1(1)*g(1,ig)+e1(2)*g(2,ig)+e1(3)*g(3,ig)) + &
             ( x0(1)*g(1,ig)+ x0(2)*g(2,ig)+ x0(3)*g(3,ig)) ) )
     enddo
     do j = 1, ny
        eigy (j) = exp( (0.d0,1.d0) * 2.d0 * pi * (j-1) * deltay * &
             (e2(1)*g(1,ig)+e2(2)*g(2,ig)+e2(3)*g(3,ig)) )
     enddo
     do k = 1, nz
        eigz (k) = exp( (0.d0,1.d0) * 2.d0 * pi * (k-1) * deltaz * &
             (e3(1)*g(1,ig)+e3(2)*g(2,ig)+e3(3)*g(3,ig)) )
     enddo
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              carica (i, j, k) = carica (i, j, k) + &
                   DREAL (rhog (ig) * eigz (k) * eigy (j) * eigx (i) )
           enddo
        enddo
     enddo

  enddo
  !
  !    compute the dipole of the charge
  !
  call volume(alat,e1(1),e2(1),e3(1),omega)

  nrwsx=125
  allocate(rws(0:3,nrwsx))
  call wsinit(rws,nrwsx,nrws,at) 

  fact=0.d0
  dipol=0.d0
  quadrupol=0.d0
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           do ipol=1,3
              r(ipol)=x0(ipol)+(i-1)*e1(ipol)*deltax + &
                               (j-1)*e2(ipol)*deltay + &
                               (k-1)*e3(ipol)*deltaz
           enddo
           fact=wsweight(r,rws,nrws)
           dipol(0) = dipol(0) + fact*carica (i, j, k) 
           do ipol=1,3
              dipol(ipol)=dipol(ipol)+fact*r(ipol)*carica(i,j,k)
              quadrupol=quadrupol + fact*carica (i, j, k)*r(ipol)**2
           enddo
        enddo
     enddo
  enddo

  dipol(0) = dipol(0)  * omega * deltax * deltay * deltaz
  do ipol=1,3
     dipol(ipol)=dipol(ipol)  * omega * deltax * deltay * deltaz * alat
  enddo
  quadrupol = quadrupol * omega * deltax * deltay * deltaz * alat**2 

  deallocate (carica)
  deallocate (rws)
  deallocate (eigz)
  deallocate (eigy)
  deallocate (eigx)
  return
end subroutine dipole_3d
!
!-----------------------------------------------------------------------
subroutine dipole_fast (alat, at, nat, tau, nrx1, nrx2, nrx3, nr1, nr2, nr3,&
     rho, bg, m1, m2, m3, x0, e1, e2, e3, dipol, quadrupol)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  integer :: nat, nrx1, nrx2, nrx3, nr1, nr2, nr3

  real(kind=DP) :: alat, tau (3, nat), at (3, 3), rho(nrx1,nrx2,nrx3), &
       bg (3, 3), e1(3), e2(3), e3(3), x0 (3), m1, m2, m3, dipol(0:3), &
       quadrupol

  integer :: nx, ny, nz, nx0, ny0, nz0, nx1, ny1, nz1, i, j, k, i1, j1, k1
  real(kind=DP), allocatable :: carica (:,:,:), rws(:,:)
  real(kind=DP) :: deltax, deltay, deltaz, debye
  real(kind=dp) :: wsweight, r(3), omega, fact
  integer :: ipol, na, nrwsx, nrws

  ! find FFT grid point closer to X0 (origin of the parallelepiped)
  ! (add 1 because r=0 correspond to n=1)

  nx0 = nint ( (x0(1)*bg(1,1) + x0(2)*bg(2,1) + x0(3)*bg(3,1) )*nr1) + 1
  ny0 = nint ( (x0(1)*bg(1,2) + x0(2)*bg(2,2) + x0(3)*bg(3,2) )*nr2) + 1
  nz0 = nint ( (x0(1)*bg(1,3) + x0(2)*bg(2,3) + x0(3)*bg(3,3) )*nr3) + 1
  !
  if ( e1(2) .ne. 0.d0  .or.  e1(3) .ne. 0.d0 .or. &
       e2(1) .ne. 0.d0  .or.  e2(3) .ne. 0.d0 .or. &
       e3(1) .ne. 0.d0  .or.  e3(2) .ne. 0.d0 )   &
       call errore ('dipole_fast','need vectors along x,y,z',1)

  ! find FFT grid points closer to X0 + e1, X0 + e2, X0 + e3
  ! (the opposite vertex of the parallelepiped)

  nx1 = nint ( ((x0(1)+m1)*bg(1,1)+x0(2)*bg(2,1)+x0(3)*bg(3,1) )*nr1)
  ny1 = nint ( (x0(1)*bg(1,2)+(x0(2)+m2)*bg(2,2)+x0(3)*bg(3,2) )*nr2)
  nz1 = nint ( (x0(1)*bg(1,3)+x0(2)*bg(2,3)+(x0(3)+m3)*bg(3,3) )*nr3)

  nx = nx1 - nx0 + 1
  ny = ny1 - ny0 + 1
  nz = nz1 - nz0 + 1

  allocate ( carica(nx, ny, nz) )    

  carica = 0.d0
  do k = nz0, nz1
     k1 = mod(k, nr3)
     if (k1.le.0) k1 = k1 + nr3
     do j = ny0, ny1
        j1 = mod(j, nr2)
        if (j1.le.0) j1 = j1 + nr2
        do i = nx0, nx1
           i1 = mod(i, nr1)
           if (i1.le.0) i1 = i1 + nr1
           carica (i-nx0+1, j-ny0+1, k-nz0+1) = rho(i1, j1, k1)
        enddo
     enddo
  enddo
  !
  ! recalculate m1, m2, m3 (the sides of the parallelepiped divided by alat)
  ! consistent with the FFT grid
  !
  WRITE( stdout,'(5x,"Requested parallelepiped sides : ",3f8.4)') m1, m2,m3
  m1 = nx * sqrt (at(1, 1) **2 + at(2, 1) **2 + at(3, 1) **2) / nr1
  m2 = ny * sqrt (at(1, 2) **2 + at(2, 2) **2 + at(3, 2) **2) / nr2
  m3 = nz * sqrt (at(1, 3) **2 + at(2, 3) **2 + at(3, 3) **2) / nr3
  WRITE( stdout,'(5x,"Redefined parallelepiped sides : ",3f8.4)') m1, m2,m3
  !
  ! recalculate x0 (the origin of the parallelepiped)
  ! consistent with the FFT grid
  !
  WRITE( stdout,'(5x,"Requested parallelepiped origin: ",3f8.4)') x0
  x0(1)=(nx0-1)*at(1,1)/ nr1 +(ny0-1)*at(1,2)/ nr2 +(nz0-1)*at(1,3)/ nr3
  x0(2)=(nx0-1)*at(2,1)/ nr1 +(ny0-1)*at(2,2)/ nr2 +(nz0-1)*at(2,3)/ nr3
  x0(3)=(nx0-1)*at(3,1)/ nr1 +(ny0-1)*at(3,2)/ nr2 +(nz0-1)*at(3,3)/ nr3
  WRITE( stdout,'(5x,"Redefined parallelepiped origin: ",3f8.4)') x0

  deltax = m1/nx 
  deltay = m2/ny 
  deltaz = m3/nz 
  !
  !    Here we check the value of the resulting charge
  !    and compute the dipole 
  !
  call volume(alat,at(1,1),at(1,2),at(1,3),omega)

  nrwsx=125
  allocate(rws(0:3,nrwsx))
  call wsinit(rws,nrwsx,nrws,at) 
  fact=0.d0
  dipol=0.d0
  quadrupol=0.d0
  do k = 1, nz  
     do j = 1, ny  
        do i = 1, nx  
           do ipol=1,3
              r(ipol)=x0(ipol)+(i-1)*e1(ipol)*deltax + &
                               (j-1)*e2(ipol)*deltay + &
                               (k-1)*e3(ipol)*deltaz
           enddo
           fact=wsweight(r,rws,nrws)
           dipol(0) = dipol(0) + fact*carica (i, j, k) 
           do ipol=1,3
              dipol(ipol)=dipol(ipol)+fact*r(ipol)*carica(i,j,k)
              quadrupol=quadrupol + fact*carica (i, j, k)*r(ipol)**2
           enddo
        enddo
     enddo
  enddo

  dipol(0) = dipol(0) * omega * deltax * deltay * deltaz 
  do ipol=1,3
     dipol(ipol)=dipol(ipol) * omega * deltax * deltay * deltaz * alat 
  enddo
  quadrupol = quadrupol * omega * deltax * deltay * deltaz * alat**2 

  if (omega > m1*m2*m3*alat**3) &
     WRITE( stdout,*) 'Warning: the box is too small to calculate dipole'
  !
  deallocate (carica)
  deallocate (rws)
  return
end subroutine dipole_fast
!
!------------------------------------------------------------
subroutine write_dipol(dipol_el,quadrupol_el,tau,nat,alat,zv, &
     ntyp,ityp,ibrav, rhotot)
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
  qq = -rhotot+zvtot
  dipol = -dipol_el(1:3)+dipol_ion
  quadrupol = -quadrupol_el+quadrupol_ion

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
  WRITE( stdout, '(/4x,"Electron charge: ",f15.8," el.")') -rhotot
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
