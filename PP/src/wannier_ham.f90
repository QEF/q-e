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

  USE io_global, ONLY: stdout, ionode, ionode_id
  USE kinds, ONLY: DP
  USE io_files,   ONLY : prefix, tmp_dir
  USE wannier_new, ONLY: nwan, use_energy_int
  USE mp,         ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE read_cards_module, ONLY : read_cards
  USE mp_global,     ONLY : mp_startup
  USE environment,   ONLY : environment_start, environment_end

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256) :: outdir, form
  INTEGER :: ios
  LOGICAL :: plot_bands
  NAMELIST /inputpp/ outdir, prefix, nwan, plot_bands, use_energy_int, form

  ! initialise environment
  !
#if defined(__MPI)
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
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix ='pwscf'
     nwan = 0
     plot_bands = .false.
     form = 'default'
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat=ios )
     !
     tmp_dir = trimcheck (outdir)

     CALL read_cards('WANNIER_AC')

  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios /= 0 ) CALL errore('wannier_ham','reading inputpp namelist',abs(ios))
  CALL read_file
  CALL openfil_pp

  CALL wannier_init(.false.)

  CALL new_hamiltonian(form, plot_bands)

  CALL environment_end ( 'WANNIER_HAM')

  CALL stop_pp

  CALL wannier_clean()

END PROGRAM wannier_ham

SUBROUTINE new_hamiltonian(form, plot_bands)

  USE io_global, ONLY: stdout, ionode, ionode_id
  USE io_files
  USE kinds, ONLY: DP
  USE wannier_new, ONLY: nwan, pp, wannier_occ, wannier_energy, wan_in
  USE klist, ONLY: nks, xk, wk
  USE lsda_mod, ONLY: isk, current_spin, lsda, nspin
  USE wvfct, ONLY: nbnd, npwx, et
  USE gvect
  USE constants,  ONLY : rytoev , tpi
  USE buffers
  USE symm_base,  ONLY : nsym

  IMPLICIT NONE
  LOGICAL :: plot_bands
  CHARACTER(len=256), INTENT(IN) :: form
  INTEGER :: i,j,k,ik, n, ios, i1, i2, n_from, n_to, seconds
  COMPLEX(DP) :: wan_func(npwx,nwan), ham(nwan,nwan,nspin), v(nwan,nwan)
  COMPLEX(DP), ALLOCATABLE :: hamk(:,:,:), hamh(:,:,:)
  real(DP), ALLOCATABLE :: ek(:,:)
  real(DP) :: e(nwan), x, hoping(3), nelec
  REAL(DP), EXTERNAL :: cclock
  CHARACTER(20) :: fmt

  ALLOCATE(ek(nwan,nks))
  ALLOCATE(hamk(nwan,nwan,nks))
  ALLOCATE(hamh(nwan,nwan,nspin))
  hamk = ZERO
  hamh = ZERO

  hoping(1) = 0.
  hoping(2) = 0.
  hoping(3) = 0.
  ek(:,:) = 0.d0
  write(fmt,*) nwan

  IF (nsym>1) THEN
     WRITE(stdout,'(/5x,a103/)') &
       'WARNING: k-points set is in the irreducible brillouin zone.',&
       ' Wannier energies and occupations are wrong!'
  ENDIF

  current_spin = 1
  CALL init_us_1
  CALL init_at_1

  CALL orthoatwfc( .true. )

  wan_func = ZERO
  pp = ZERO
  ham = ZERO

  DO ik = 1, nks
     write(stdout,*) '       Computing k-point', ik
     IF (lsda) current_spin  = isk(ik)
     CALL wannier_proj(ik,wan_func)

     pp = ZERO
     CALL get_buffer( pp, nwordwpp, iunwpp, ik)

     hamk(:,:,ik) = ZERO

     DO i=1, nwan
        DO j=1,nwan
           n_from = int (wan_in(i,current_spin)%bands_from )
           n_to   = int (wan_in(i,current_spin)%bands_to )
           DO n = n_from, n_to
              ! On-site hamiltonian
              ham(i,j,current_spin) = ham(i,j,current_spin) + &
                 pp(i,n)*cmplx(et(n,ik),0.d0,kind=DP)*conjg(pp(j,n))*wk(ik)
              ! Hoping integrals
              hamh(i,j,current_spin) = hamh(i,j,current_spin) + &
                 pp(i,n)*cmplx(et(n,ik),0.d0,kind=DP)*conjg(pp(j,n))*wk(ik)*&
                 exp( (0.d0,1.d0)*tpi* (xk(1,ik)*hoping(1) + &
                     xk(2,ik)*hoping(2) + xk(3,ik)*hoping(3)) )
              ! Current k-point hamiltonian
              hamk(i,j,ik) = hamk(i,j,ik) + pp(i,n)*conjg(pp(j,n))* &
                             cmplx(et(n,ik),0.d0,kind=DP)
           ENDDO
        ENDDO
     ENDDO

     IF (plot_bands) CALL cdiagh(nwan,hamk(:,:,ik),nwan,ek(:,ik),v)

     !Hermicity check
     DO i=1,nwan
        DO j=1,nwan
           IF(abs(hamk(i,j,ik)-conjg(hamk(j,i,ik)))>=1.d-8) THEN
              WRITE(stdout,'(5x,"Wrong elements", 2i3," in",i4," k-point")') i,j,ik
              CALL errore ('wannier_ham', 'Hamiltonian is not hermitian', ik)
           ENDIF
        ENDDO
     ENDDO
  ENDDO !ik

  !Compute wannier parameters
  CALL wannier_occupancies(wannier_occ)
  CALL wannier_enrg(wannier_energy)

  !output computed
  DO j=1, nspin
     WRITE(stdout,'(/5x,a4,i2,a)') 'Spin', j,':'
     DO i=1, nwan
        WRITE(stdout,'(7x,a8,i3)') 'Wannier#',i
        WRITE(stdout,'(9x,a11,f5.3)') 'occupation:',wannier_occ(i,i,j)
        WRITE(stdout,'(9x,a7,f7.3,a3)') 'energy:',wannier_energy(i,j)*rytoev,' eV'
     ENDDO
     WRITE(stdout,'(7x,a26/)')'Wannier occupation matrix:'
     DO i=1,nwan
        WRITE(stdout,'(7x,'// ADJUSTL(fmt) //'f7.3)') (wannier_occ(i,k,j),k=1,nwan)
     ENDDO
  ENDDO
  !end of output

  ! write HMLT file
  IF (form == 'amulet') THEN

    seconds = cclock()

    CALL write_hamiltonian_amulet(nwan,hamk,seconds,114)

    nelec = 0.d0
    DO i=1,nwan
      nelec = nelec + SUM(wannier_occ(i,i,:))
    END DO

    IF (nspin == 1) nelec = nelec*2.d0

    CALL write_systemdata_amulet(seconds,nelec,118)

  ELSE

    CALL write_hamiltonian_default(nwan,hamk,114)
  
  END IF

  IF(nspin==1) THEN
     ham = 5.d-1*ham
     hamh = 5.d-1*hamh
  ENDIF

  DO i=1, nspin
     WRITE(stdout,*) ' '

     CALL cdiagh(nwan,ham(:,:,i),nwan,e,v)
     WRITE(stdout,'(5x,a39)') 'Projected Hamiltonian eigenvalues (eV):'
     WRITE(stdout,'(6x,a5,i1,4x,'// ADJUSTL(fmt) //'f9.4)') 'spin', i, (e(j)*rytoev,j=1,nwan)
     WRITE(stdout,*) ' '

! hopings integrals
     IF(any(hoping/=0.d0)) THEN
        WRITE(stdout,'(5x,a44,3f6.2,a5)') 'Hopings from the atom in origin to direction', (hoping(j),j=1,3), 'are:'
        DO j=1,nwan
           WRITE(stdout,'(5x,'// ADJUSTL(fmt) //'f9.5)') (dreal(hamh(j,n,i))*rytoev, n=1, nwan)
        ENDDO
        WRITE(stdout,*) ' '
     ENDIF

! additional check: hamiltonian should be hermitian
     IF(sum(dimag(hamh))>=1d-9) THEN
        WRITE(stdout,*) 'ATTENTION! Hamiltonian is NOT hermitian'
        WRITE(stdout,*) 'Imaginary part is:'
        DO j=1,nwan
           WRITE(stdout,'('// ADJUSTL(fmt) //'f9.5)') (dimag(hamh(j,n,i))*rytoev, n=1, nwan)
        ENDDO
        WRITE(stdout,*) '---'
     ENDIF
  ENDDO

  IF(plot_bands) CALL plot_wannier_bands(ek)

  DEALLOCATE(ek)
  DEALLOCATE(hamk)
  DEALLOCATE(hamh)

END SUBROUTINE new_hamiltonian

SUBROUTINE plot_wannier_bands(ek)
! This routine produces three files wannier_bands.dat, original_bands.dat
! and wannier_bands.plot to visual check how generated Wannier-Hamiltonian
! reproduses original bands structure. To check just type 'gnuplot wannier_bands.plot'
! in your terminal window. Of course one can use another ploting software for that purpose,
! for example 'xmgrace original_bands.dat wannier_bands.dat'

  USE constants, ONLY: rytoev
  USE io_global, ONLY: stdout, ionode, ionode_id
  USE io_files
  USE kinds, ONLY: DP
  USE klist, ONLY: nks, xk
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: nbnd, et
  USE wannier_new, ONLY: nwan
  USE ener, ONLY: ef

  IMPLICIT NONE
  REAL(DP), INTENT(in) :: ek(nwan,nks)

  INTEGER :: i,j,k,ik,ios
  REAL(DP) :: x, emax, emin
  CHARACTER(len=1) :: backslash

  OPEN (unit = 113, file = 'wannier_bands.dat', status = 'unknown', form = 'formatted', err = 400, iostat = ios)
  OPEN (unit = 114, file = 'original_bands.dat', status = 'unknown', form = 'formatted', err = 401, iostat = ios)
  OPEN (unit = 115, file = 'wannier_bands.plot', status = 'unknown', form = 'formatted', err = 402, iostat = ios)
400 CALL errore ('plot_wannier_bands', 'wannier_bands.dat', abs (ios) )
401 CALL errore ('plot_wannier_bands', 'original_bands.dat', abs (ios) )
402 CALL errore ('plot_wannier_bands', 'wannier_bands.plot', abs (ios) )

  emax = ek(1,1)
  emin = ek(1,1)

  DO i=1, nwan
     x = 0.d0
     DO ik=1, nks/nspin
        ! find limits for pretty plotting
        IF (emax<ek(i,ik)*rytoev) emax = ek(i,ik)*rytoev
        IF (emin>ek(i,ik)*rytoev) emin = ek(i,ik)*rytoev
        !
        WRITE(113,'(2f15.9)') x, ek(i,ik)*rytoev
        IF (ik/=nks) THEN
           x = x + sqrt((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
        ENDIF
     ENDDO
     WRITE(113, '(2a)') '  '
  ENDDO
  DO i=1, nbnd
     x = 0.d0
     DO ik=1, nks/nspin
        WRITE(114,'(2f15.9)') x, et(i,ik)*rytoev
        IF (ik/=nks) THEN
           x = x + sqrt((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
        ENDIF
     ENDDO
     WRITE(114, '(2a)') '  '
  ENDDO

  IF (nspin==2) THEN
     DO i=1, nwan
        x = 0.d0
        DO ik=nks/2+1, nks
           ! find limits for pretty plotting
           IF (emax<ek(i,ik)*rytoev) emax = ek(i,ik)*rytoev
           IF (emin>ek(i,ik)*rytoev) emin = ek(i,ik)*rytoev
           !
           WRITE(113,'(2f15.9)') x, ek(i,ik)*rytoev
           IF (ik/=nks) THEN
              x = x + sqrt((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
           ENDIF
        ENDDO
        WRITE(113, '(2a)') '  '
     ENDDO
     DO i=1, nbnd
        x = 0.d0
        DO ik=nks/2+1, nks
           WRITE(114,'(2f15.9)') x, et(i,ik)*rytoev
           IF (ik/=nks) THEN
              x = x + sqrt((xk(1,ik)-xk(1,ik+1))**2+(xk(2,ik)-xk(2,ik+1))**2+(xk(3,ik)-xk(3,ik+1))**2)
           ENDIF
        ENDDO
        WRITE(114, '(2a)') '  '
     ENDDO
  ENDIF

  WRITE(115,*)'reset'
  WRITE(115,*)'set terminal postscript eps enhanced color'
  WRITE(115,*)'set output "wannier_bands.eps"'
  WRITE(115,*)'unset xtics'
  WRITE(115,'(a12,f7.3,a,f7.3,a)')'set yrange [',emin-1.5,':',emax+1.5,']'
  WRITE(115,*)'set style line 1 lt 1 lc rgb "black" lw 2'
  WRITE(115,*)'set style line 2 lt 2 lc rgb "red" lw 2'
  WRITE(115,*)'set style line 3 lt 1 lc rgb "green" lw 1'
  WRITE(115,*)'set ylabel "Energy (eV)"'
  backslash = char(92)
  WRITE(115,*)'plot '//backslash
  WRITE(115,*)'"original_bands.dat" title "LDA bands" with lines linestyle 1,'&
            & //backslash
  WRITE(115,*)'"wannier_bands.dat" title "Wannier bands" with lines linestyle 2,'&
            & //backslash
  WRITE(115,'(f7.3,a44)') ef*rytoev,'title "Fermi energy" with lines linestyle 3'

  CLOSE(113)
  CLOSE(114)
  CLOSE(115)

END SUBROUTINE plot_wannier_bands
