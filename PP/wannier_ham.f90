! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

!----------------------------------------------------------------------- 
PROGRAM wannier_ham
!----------------------------------------------------------------------- 
! 
! This program generates Hamiltonian matrix on Wannier-functions basis
  
  use io_global, only: stdout, ionode, ionode_id
  use kinds, only: DP 
  USE io_files,   ONLY : prefix, tmp_dir, trimcheck
  use wannier_new, only: nwan, use_energy_int
  use ktetra
  USE mp,         ONLY : mp_bcast
  USE read_cards_module, ONLY : read_cards
  USE mp_global,     ONLY : mp_startup
  USE environment,   ONLY : environment_start

  implicit none
  CHARACTER(len=256) :: outdir 
  integer :: ios
  logical :: plot_bands, u_matrix
  real(DP) :: U,J
  namelist /inputpp/ outdir, prefix, nwan, plot_bands, use_energy_int, u_matrix
  namelist /Umatrix/ U,J
 
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'WANNIER_HAM')
  !
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
     plot_bands = .FALSE.
     u_matrix=.FALSE.
     !
     U=0.d0
     J=0.d0
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat=ios )
     IF(u_matrix) READ (5, Umatrix, iostat=ios )
     !
     tmp_dir = trimcheck (outdir)
     
     CALL read_cards('WANNIER_AC')
     
  END IF
  !
  CALL mp_bcast( ios, ionode_id )
  IF ( ios /= 0 ) CALL errore('wannier_ham','reading inputpp namelist',ABS(ios))
  call read_file 
  call openfil_pp

  call wannier_init(.FALSE.)

  call new_hamiltonian(plot_bands)
  
  if(u_matrix) call wannier_u_matrix(U,J)

  call stop_pp 
  
  call wannier_clean()
  
END PROGRAM wannier_ham

SUBROUTINE new_hamiltonian(plot_bands)

  use io_global, only: stdout, ionode, ionode_id
  use io_files
  use kinds, only: DP 
  use wannier_new, only: nwan, pp, wannier_occ, wannier_energy,wan_in
  use klist, only: nks, xk, wk
  use lsda_mod, only: isk, current_spin, lsda, nspin
  use wvfct, only: nbnd, npwx, igk, npw, g2kin, et
  use gvect
  use cell_base, only: tpiba2
  use constants,  ONLY : rytoev , tpi
  use buffers
  USE symm_base,  ONLY : nsym

  implicit none
  logical :: plot_bands
  integer :: i,j,k,ik, n, ios, i1, i2, outfile, n_from, n_to
  complex(DP) :: wan_func(npwx,nwan), ham(nwan,nwan,nspin), v(nwan,nwan)
  complex(DP), allocatable :: hamk(:,:,:), hamh(:,:,:)
  real(DP), allocatable :: ek(:,:)
  real(DP) :: e(nwan), x, hoping(3)

  ! HMLT file unit
  outfile = 114

  allocate(ek(nwan,nks))
  allocate(hamk(nwan,nwan,nks))
  allocate(hamh(nwan,nwan,nspin))
  hamk = ZERO
  hamh = ZERO

  hoping(1) = 0.
  hoping(2) = 0.
  hoping(3) = 0.
  ek(:,:) = 0.d0
  
  IF (nsym.GT.1) THEN
     write(stdout,'(/5x,a103/)') &
       'WARNING: k-points set is in the irreducible brillouin zone.',&
       ' Wannier energies and occupations are wrong!'
  END IF

  current_spin = 1
  call init_us_1 
  call init_at_1 
  
  ! Generating igk for orthoatwfc()
  REWIND( iunigk )
  DO ik = 1, nks
     CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     IF ( nks > 1 ) WRITE( iunigk ) igk
  END DO
  !
  
  CALL orthoatwfc()

  wan_func = ZERO
  pp = ZERO
  ham = ZERO

  do ik = 1, nks
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     if (lsda) current_spin  = isk(ik)
     call wannier_proj(ik,wan_func)
     
     pp = ZERO
     call get_buffer( pp, nwordwpp, iunwpp, ik)

     hamk(:,:,ik) = ZERO

     do i=1, nwan
        do j=1,nwan
           n_from = INT (wan_in(i,current_spin)%bands_from )
           n_to   = INT (wan_in(i,current_spin)%bands_to )
           do n = n_from, n_to
              ! On-site hamiltonian
              ham(i,j,current_spin) = ham(i,j,current_spin) + &
                 pp(i,n)*CMPLX(et(n,ik),0.d0,KIND=DP)*CONJG(pp(j,n))*wk(ik)
              ! Hoping integrals
              hamh(i,j,current_spin) = hamh(i,j,current_spin) + &
                 pp(i,n)*CMPLX(et(n,ik),0.d0,KIND=DP)*CONJG(pp(j,n))*wk(ik)*&
                 cdexp( (0.d0,1.d0)*tpi* (xk(1,ik)*hoping(1) + &
                     xk(2,ik)*hoping(2) + xk(3,ik)*hoping(3)) )
              ! Current k-point hamiltonian
              hamk(i,j,ik) = hamk(i,j,ik) + pp(i,n)*CONJG(pp(j,n))* &
                             CMPLX(et(n,ik),0.d0,KIND=DP)
              !Overlap mtrx in current k-point (for debug purposes)
           end do
        end do
     end do

     if (plot_bands) call cdiagh(nwan,hamk(:,:,ik),nwan,ek(:,ik),v)
     
     !Hermicity check
     do i=1,nwan
        do j=1,nwan
           if(abs(hamk(i,j,ik)-CONJG(hamk(j,i,ik))).ge.1.d-8) then
              write(stdout,'(5x,"Wrong elements", 2i3," in",i4," k-point")') i,j,ik
              call errore ('wannier_ham', 'Hamiltonian is not hermitian', ik)
           end if
        end do
     end do
  end do !ik
  
  !Compute wannier parameters
  call wannier_occupancies(wannier_occ)
  call wannier_enrg(wannier_energy)
  
  !output computed
  do j=1, nspin
     write(stdout,'(/5x,a4,i2,a)') 'Spin', j,':'
     do i=1, nwan
        write(stdout,'(7x,a8,i3)') 'Wannier#',i
        write(stdout,'(9x,a11,f5.3)') 'occupation:',wannier_occ(i,i,j)
        write(stdout,'(9x,a7,f7.3,a3)') 'energy:',wannier_energy(i,j)*rytoev,' eV'
     end do
     write(stdout,'(7x,a26/)')'Wannier occupation matrix:'
     do i=1,nwan
        write(stdout,'(7x,50f7.3)') (wannier_occ(i,k,j),k=1,nwan)
     end do
  end do
  !end of output
  
  ! write HMLT file
  open (outfile, file = 'hamilt', status = 'unknown', form = 'formatted', err = 300, iostat = ios)
300 call errore ('HMLT', 'Opening hamilt', abs (ios) )
  
  call wannier_hamiltonian_JK(nwan,hamk,outfile)
  
  close(outfile)

  if(nspin.eq.1) then
     ham = 5.d-1*ham
     hamh = 5.d-1*hamh
  end if

  do i=1, nspin
     write(stdout,*) ' '

     call cdiagh(nwan,ham(:,:,i),nwan,e,v)
     write(stdout,'(5x,a39)') 'Projected Hamiltonian eigenvalues (eV):'
     write(stdout,'(6x,a5,i1,4x,50f9.4)') 'spin', i, (e(j)*rytoev,j=1,nwan)
     write(stdout,*) ' '
     
! hopings integrals
     if(ANY(hoping.ne.0.d0)) then
        write(stdout,'(5x,a44,3f6.2,a5)') 'Hopings from the atom in origin to direction', (hoping(j),j=1,3), 'are:'
        do j=1,nwan
           write(stdout,'(5x,20f9.5)') (dreal(hamh(j,n,i))*rytoev, n=1, nwan)
        end do
        write(stdout,*) ' '
     end if

! additional check: hamiltonian should be hermitian
     if(SUM(dimag(hamh)).ge.1d-9) then
        write(stdout,*) 'ATTENTION! Hamiltonian is NOT hermitian'
        write(stdout,*) 'Imaginary part is:' 
        do j=1,nwan
           write(stdout,'(20f9.5)') (dimag(hamh(j,n,i))*rytoev, n=1, nwan)
        end do
        write(stdout,*) '---' 
     end if
  end do

  if(plot_bands) call plot_wannier_bands(ek)

  deallocate(ek)
  deallocate(hamk)
  deallocate(hamh)

END SUBROUTINE new_hamiltonian

SUBROUTINE plot_wannier_bands(ek)
! This routine produces three files wannier_bands.dat, original_bands.dat
! and wannier_bands.plot to visual check how generated Wannier-Hamiltonian
! reproduses original bands structure. To check just type 'gnuplot wannier_bands.plot'
! in your terminal window. Of course one can use another ploting software for that purpose,
! for example 'xmgrace original_bands.dat wannier_bands.dat'

  USE constants, ONLY: rytoev
  use io_global, only: stdout, ionode, ionode_id
  use io_files
  use kinds, only: DP 
  use klist, only: nks, xk
  use lsda_mod, only: nspin
  use wvfct, only: nbnd, et
  use wannier_new, only: nwan
  use ener, only: ef

  implicit none
  REAL(DP), intent(in) :: ek(nwan,nks)
  
  INTEGER :: i,j,k,ik,ios
  REAL(DP) :: x, emax, emin
 
  open (unit = 113, file = 'wannier_bands.dat', status = 'unknown', form = 'formatted', err = 400, iostat = ios)
  open (unit = 114, file = 'original_bands.dat', status = 'unknown', form = 'formatted', err = 401, iostat = ios)
  open (unit = 115, file = 'wannier_bands.plot', status = 'unknown', form = 'formatted', err = 402, iostat = ios)
400 call errore ('plot_wannier_bands', 'wannier_bands.dat', abs (ios) )
401 call errore ('plot_wannier_bands', 'original_bands.dat', abs (ios) )
402 call errore ('plot_wannier_bands', 'wannier_bands.plot', abs (ios) )

  emax = ek(1,1)
  emin = ek(1,1)

  do i=1, nwan
     x = 0.d0
     do ik=1, nks/nspin
        ! find limits for pretty plotting
        if (emax.lt.ek(i,ik)*rytoev) emax = ek(i,ik)*rytoev
        if (emin.gt.ek(i,ik)*rytoev) emin = ek(i,ik)*rytoev
        !
        write(113,'(2f15.9)') x, ek(i,ik)*rytoev
        if (ik.ne.nks) then 
           x = x + SQRT((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
        end if
     end do
     write(113, '(2a)') '  '
  end do
  do i=1, nbnd
     x = 0.d0
     do ik=1, nks/nspin
        write(114,'(2f15.9)') x, et(i,ik)*rytoev
        if (ik.ne.nks) then 
           x = x + SQRT((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
        end if
     end do
     write(114, '(2a)') '  '
  end do
  
  if (nspin.eq.2) then
     do i=1, nwan
        x = 0.d0
        do ik=nks/2+1, nks
           ! find limits for pretty plotting
           if (emax.lt.ek(i,ik)*rytoev) emax = ek(i,ik)*rytoev
           if (emin.gt.ek(i,ik)*rytoev) emin = ek(i,ik)*rytoev
           !
           write(113,'(2f15.9)') x, ek(i,ik)*rytoev
           if (ik.ne.nks) then 
              x = x + SQRT((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
           end if
        end do
        write(113, '(2a)') '  '
     end do
     do i=1, nbnd
        x = 0.d0
        do ik=nks/2+1, nks
           write(114,'(2f15.9)') x, et(i,ik)*rytoev
           if (ik.ne.nks) then 
              x = x + SQRT((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
           end if
        end do
        write(114, '(2a)') '  '
     end do
  end if
  
  write(115,*)'reset'
  write(115,*)'set term post eps'
  write(115,*)'set output "wannier_bands.eps"'
  write(115,*)'unset xtics'
  write(115,'(a12,f7.3,a,f7.3,a)')'set yrange [',emin-1.5,':',emax+1.5,']'
  write(115,*)'set style line 1 lt 1 lc rgb "black" lw 2'
  write(115,*)'set style line 2 lt 2 lc rgb "red" lw 2'
  write(115,*)'set style line 3 lt 1 lc rgb "green" lw 1'
  write(115,*)'set ylabel "Energy (eV)"'
  write(115,*)'plot \\'
  write(115,*)'"original_bands.dat" title "LDA bands" with lines linestyle 1,\\'
  write(115,*)'"wannier_bands.dat" title "Wannier bands" with lines linestyle 2,\\'
  write(115,'(f7.3,a44)') ef*rytoev,'title "Fermi energy" with lines linestyle 3'

  close(113)
  close(114)
  close(115)
  
END SUBROUTINE plot_wannier_bands
