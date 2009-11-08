! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

!----------------------------------------------------------------------- 
PROGRAM wannier_composition
!----------------------------------------------------------------------- 
! 
! This program plots charge density of selected wannier function in
! IBM Data Explorer format
  
  use io_global, only: stdout, ionode, ionode_id
  use kinds,         ONLY : DP 
  USE io_files,      ONLY : prefix, tmp_dir, trimcheck
  use wannier_new,   ONLY : nwan, plot_wan_num, plot_wan_spin
  USE mp,            ONLY : mp_bcast
  USE io_global,     ONLY : ionode, stdout
  USE mp_global,     ONLY : mp_startup
  USE control_flags, ONLY : use_task_groups, ortho_para
  USE environment,   ONLY : environment_start

  implicit none
  CHARACTER(len=256) :: outdir 
  integer :: ios,nc(3),n0(3)
  namelist /inputpp/ outdir, prefix, nwan, plot_wan_num, plot_wan_spin, nc, n0
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( use_task_groups, ortho_para )
#endif
  CALL environment_start ( 'WANNIER_PLOT' )
 
  ios = 0
  !
  IF ( ionode ) THEN
     !
     !   set default values for variables in namelist
     !
     CALL get_env( 'ESPRESSO_TMPDIR', outdir )
     IF ( TRIM( outdir ) == ' ' ) outdir = './'
     prefix ='pwscf'
     nwan = 0
     plot_wan_spin=1
 
     nc(1) = 3
     nc(2) = 3
     nc(3) = 3
     n0(1) = -1
     n0(2) = -1
     n0(3) = -1
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat=ios )
     !
     tmp_dir = trimcheck (outdir)
  END IF
  !
  CALL mp_bcast( ios, ionode_id )
  IF ( ios /= 0 ) CALL errore('wannier_ham','reading inputpp namelist',ABS(ios))
  call read_file 
  call openfil_pp

  call wannier_init(.true.)

  !debug
  write(stdout,'(5x,"Calling plot_wannier for wannier",i3)') plot_wan_num
  !end of debug
  call plot_wannier(nc,n0)
  !debug
  write(stdout,'(5x,"Calling plot_atoms")')
  !end of debug
  call plot_atoms()
 
  call stop_pp 
 
  call wannier_clean()
  
END PROGRAM wannier_composition

SUBROUTINE plot_wannier(nc,n0)

  use io_global, only: stdout, ionode, ionode_id
  use io_files
  use kinds, only: DP 
  use wannier_new, only: nwan,plot_wan_num,plot_wan_spin
  use klist, only: nks, xk, wk
  use lsda_mod, only: isk, current_spin, lsda, nspin
  use wvfct, only: nbnd, npwx, igk, npw, g2kin
  use constants,  ONLY : rytoev , tpi
  use buffers
  USE symme,           ONLY : nsym
  USE ldaU,             ONLY : swfcatom
  use gvect
  use gsmooth
  use cell_base
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau, atm, zv
  USE vlocal,    ONLY : strf


  implicit none
  integer, intent(in) :: nc(3), n0(3)
  integer :: i,j, k, ik, n, ir, ios, n1, n2, n3,i1,j1,k1
  COMPLEX(DP) :: phase
  COMPLEX(DP), allocatable :: wan_func(:,:), pp_ort(:,:), psic(:), psic3(:,:,:), psic3_0(:,:,:), psic_sum(:,:,:,:), paux(:,:)
  real(DP), allocatable :: rho(:,:,:,:), raux(:)
  real(DP) :: r(3)
  
  IF (nsym.GT.1) THEN
     call errore('wannier_cmptn','k-points set is in the irreducible brillouin zone - not implemented',1)
  END IF

  allocate(wan_func(npwx,nwan))
  allocate(psic(nrxxs))
  allocate(psic3(nrx1s,nrx2s,nrx3s))
  allocate(psic3_0(nrx1s,nrx2s,nrx3s))
  allocate(psic_sum(nc(1)*nrx1s,nc(2)*nrx2s,nc(3)*nrx3s,nspin))
  allocate(rho(nc(1)*nrx1s,nc(2)*nrx2s,nc(3)*nrx3s,nspin))

  call init_us_1 
  call init_at_1 

  CALL struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
       strf, eigts1, eigts2, eigts3)

  current_spin = 1
  wan_func = ZERO
  psic3 = ZERO
  psic3_0 = ZERO
  psic_sum = ZERO

  do ik = 1, nks
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     if (lsda) current_spin  = isk(ik)
     
     wan_func = ZERO
     call get_buffer( wan_func, nwordwf, iunwf, ik)
     
     psic(1:nrxxs) = ZERO
     rho = ZERO
     do j = 1, npw
        psic (nls (igk (j) ) ) = wan_func (j, plot_wan_num)
     end do
     
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     
     do k=1, nrx3s
        do j=1,nrx2s
           do i=1,nrx1s
              n = i + (j-1)*nrx1s + (k-1)*nrx2s*nrx1s
              psic3_0(i,j,k) = psic(n)
           end do
        end do
     end do
     
     do k=1, (nrx3s-1)*nc(3)
        do j=1, (nrx2s-1)*nc(2)
           do i=1, (nrx1s-1)*nc(1)
              ! r = n0(1)*at(1,:)+n0(2)*at(2,:)+n0(3)*at(3,:)
              ! r = r + DBLE(i-1)*at(1,:)/DBLE(nrx1s-1)+DBLE(j-1)*at(2,:)/DBLE(nrx2s-1)+DBLE(k-1)*at(3,:)/DBLE(nrx3s-1)
              r = n0(1)*at(:,1)+n0(2)*at(:,2)+n0(3)*at(:,3)
              r = r + DBLE(i-1)*at(:,1)/DBLE(nrx1s-1) + &
                      DBLE(j-1)*at(:,2)/DBLE(nrx2s-1) + &
                      DBLE(k-1)*at(:,3)/DBLE(nrx3s-1)
              phase = cos(tpi*(xk(1,ik)*r(1)+xk(2,ik)*r(2)+xk(3,ik)*r(3))) + &
          (0.d0,1.d0)*sin(tpi*(xk(1,ik)*r(1)+xk(2,ik)*r(2)+xk(3,ik)*r(3)))
              
              i1 = i - FLOOR(DBLE(i-0.01)/DBLE(nrx1s-1))*(nrx1s-1)
              j1 = j - FLOOR(DBLE(j-0.01)/DBLE(nrx2s-1))*(nrx2s-1)
              k1 = k - FLOOR(DBLE(k-0.01)/DBLE(nrx3s-1))*(nrx3s-1)
              psic_sum(i,j,k,current_spin) = psic_sum(i,j,k,current_spin)+ &
                   CMPLX(wk(ik),0.d0,KIND=DP)*psic3_0(i1,j1,k1)*phase
           end do
        end do
     end do
     
  end do !ik

  rho = 0.d0

  do n=1, nspin
     do i=1, nrx1s*nc(1)
        do j=1, nrx2s*nc(2)
           do k=1,nrx3s*nc(3)
              rho(i,j,k,n) = dreal(psic_sum(i,j,k,n))**2+aimag(psic_sum(i,j,k,n))**2
           end do
        end do
     end do
  end do
 
  open (10, file='wannier.plot.dx', err = 100, iostat = ios)
100 call errore ('plot_wannier', 'Opening out file', abs (ios) )

  ! I want to write .dx file for dataexplorer
  write(10,'(a36,3i6)') 'object 1 class gridpositions counts ', nrx3s*nc(3), nrx2s*nc(2), nrx1s*nc(1)
  write(10,*) 'origin', n0(1)*at(:,1)+n0(2)*at(:,2)+n0(3)*at(:,3)
  !  write(10,'(a5, 3f9.5)') 'delta', (at(3,i)/(1.d0*(nrx3s-1)),i=1,3)
  !  write(10,'(a5, 3f9.5)') 'delta', (at(2,i)/(1.d0*(nrx2s-1)),i=1,3)
  !  write(10,'(a5, 3f9.5)') 'delta', (at(1,i)/(1.d0*(nrx1s-1)),i=1,3)
  write(10,'(a5, 3f9.5)') 'delta', (at(i,1)/(1.d0*(nrx3s-1)),i=1,3)
  write(10,'(a5, 3f9.5)') 'delta', (at(i,2)/(1.d0*(nrx2s-1)),i=1,3)
  write(10,'(a5, 3f9.5)') 'delta', (at(i,3)/(1.d0*(nrx1s-1)),i=1,3)
  write(10,'(a38,3i6)') 'object 2 class gridconnections counts ', nrx3s*nc(3), nrx2s*nc(2), nrx1s*nc(1)
  write(10,*) 'attribute "element type" string "cubes"'
  write(10,*) 'attribute "ref" string "positions"'
  write(10,'(a44,i10,a13)') 'object 3 class array type float rank 0 items', nrx3s*nc(3)*nrx2s*nc(2)*nrx1s*nc(1), 'data follows'
  
  do i=1, nrx3s*nc(3)
     do j=1,nrx2s*nc(2)
        do k=1,nrx1s*nc(1)
           write(10,'(f13.7)') rho(k,j,i,plot_wan_spin)
           ! write(10,'(f13.7)') aimag(psic_sum(k,j,i,plot_wan_spin))
        end do
     end do
  end do

  write(10,'(a34)') 'attribute "dep" string "positions"'
  write(10,*) 'object "regular positions regular connections" class field'
  write(10,*) 'component "positions" value 1'
  write(10,*) 'component "connections" value 2'
  write(10,*) 'component "data" value 3'
  write(10,*) 'end'
  
  close(10)
  
  deallocate(wan_func)
  deallocate(psic)
  deallocate(psic3)
  deallocate(psic3_0)
  deallocate(psic_sum)
  deallocate(rho)


END SUBROUTINE plot_wannier


SUBROUTINE plot_atoms
  use io_global, only: stdout
  use kinds, only: DP 
  use ions_base, only: tau, nat, ityp, zv
  implicit none
  integer :: i,na, ios
  
  open (20, file='atoms.plot.dx', err = 200, iostat = ios)
200 call errore ('plot_wannier', 'Opening out atoms file', abs (ios) )
  
  write(20,*) 'object 1 class array type float rank 1 shape 3 items', nat,' data follows'
  do na = 1, nat
     write(20,'(3f9.5)') (tau(i,na),i=1,3)
  enddo
  write(20,*) 'object 2 class array type float rank 0 items', nat,' data follows'
  do na = 1, nat
     write(20,*) zv(ityp(na))
  enddo
  write(20,*) 'attribute "dep" string "positions"'
  write(20,*) 'object "irregular positions" class field'
  write(20,*) 'component "positions" value 1'
  write(20,*) 'component "data" value 2'
  write(20,*) 'end'
  close(20)
END SUBROUTINE plot_atoms
