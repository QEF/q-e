
subroutine dielectric(sh,input,kpnts,energy)
  !------------------------------------------------------------------------
  !
  ! This subroutine diagonalizes the k-dependent Hamiltonian for every k in interp_grid
  ! and calculates the IP dielectric function
  !
  USE simple_ip_objects
  USE input_simple_ip
  USE tetra_ip, ONLY : tetrahedra1, weights_delta1
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : rytoev, pi
  USE mp_world,  ONLY : world_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(shirley) :: sh
  TYPE(input_options_simple_ip) :: input
  TYPE(energies) :: energy
  TYPE(kpoints) :: kpnts
  TYPE(eigen) :: eig
  REAL(kind=DP) :: q(3)
  INTEGER :: ik, idir, ii , counting , iun
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: tmp , tmp2
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: commut_interp
  REAL(kind=DP), EXTERNAL :: efermig , efermit
  REAL(kind=DP), EXTERNAL ::  w0gauss
  REAL(kind=DP), EXTERNAL ::  wgauss
  INTEGER ::  is , iw , iband1, iband2 , num_energy_dos
  REAL(kind=DP) :: nelec, ef, arg, tpiovera, alpha, diff_occ, delta_e, w
  REAL(kind=DP), DIMENSION(6) :: plasma_freq
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: wk
  INTEGER, DIMENSION(:), ALLOCATABLE :: isk
  REAL(kind=DP), DIMENSION(:,:), ALLOCATABLE :: energy_sum
  REAL(kind=DP), DIMENSION(:,:), ALLOCATABLE :: focc, focc_tmp, eps_im, eps_re, eps_tot_re, eps_tot_im, eels
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: wgrid, wgrid_dos, dos , jdos , eps_avg_re, eps_avg_im 
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: refractive_index , extinct_coeff , reflectivity
  REAL(kind=DP) :: lorentzian, matrix_element2 , norm_epsilon
  REAL(kind=DP) :: delta_e_thr = 1.0E-6   ! in eV
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE  :: matrix_element
  INTEGER :: num_tetra
  INTEGER, allocatable :: tetra(:,:)
  REAL(kind=DP), allocatable :: weights(:) 
  CHARACTER(len=256), DIMENSION(6) :: headline
  CHARACTER(256) :: str
  INTEGER, EXTERNAL :: find_free_unit
  !
  call start_clock('dielectric')
  !
  call initialize_eigen(eig)
  !
  allocate(tmp(sh%ntot_e,sh%num_bands))
  allocate(tmp2(sh%num_bands,sh%num_bands))
  if (input%nonlocal_commutator .and. sh%nonlocal_commutator) then
    allocate(commut_interp(sh%ntot_e,sh%ntot_e,3))
  endif
  allocate(focc(sh%num_bands,kpnts%nk),focc_tmp(sh%num_bands,energy%nk_loc))
  allocate(wgrid(input%nw),jdos(input%nw),eps_im(input%nw,3),eps_re(input%nw,3),eels(input%nw,3))
  allocate(eps_tot_im(input%nw,3),eps_tot_re(input%nw,3))
  allocate(eps_avg_im(input%nw),eps_avg_re(input%nw),refractive_index(input%nw),extinct_coeff(input%nw),reflectivity(input%nw))
  allocate(matrix_element(sh%num_bands,sh%num_bands))
  !
  !!!!!!! Interpolation of the bands
  counting = 0  !DEBUG
  tpiovera = 2.d0*pi/sh%alat
  write(stdout,*) ' '
  write(stdout,*) 'Computing the band structure...'
  write(stdout,*) ' '
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ik,q) FIRSTPRIVATE(eig) 
  do ik=1,energy%nk_loc
     !
     q(1:3) = kpnts%qk(1:3,ik + energy%ik_first - 1)
     call diagonalization(q,sh,input,eig,ik,kpnts)
     energy%energy(1:sh%num_bands,ik) = eig%energy(1:sh%num_bands) ! Save bands in energy%energy
     !
     counting = counting + 1  !DEBUG
     write(stdout,*)  '*****************************************************************'  !DEBUG
     write(stdout,*) 'k-point:', counting
     write(stdout,*) 'k-point coordinates:', q  !DEBUG
     write(stdout,*) 'Interpolated bands:' , energy%energy(1:sh%num_bands,ik)*rytoev !DEBUG
     write(stdout,*) ' '
  enddo
  !$OMP END PARALLEL DO
  !!!!!!! End interpolation of the bands
  !
  allocate(energy_sum(sh%num_bands,kpnts%nk))
  energy_sum = 0.d0
  ! Gather all the bands divided in the various processors by subgroups of k-points
  do ik =1, energy%nk_loc
    energy_sum(1:sh%num_bands,ik+energy%ik_first-1) = energy%energy(1:sh%num_bands, ik)
  enddo
  call mp_sum(energy_sum,world_comm)
  !
  !!!!!!! Fermi energy calculation
  if (input%fermi_energy == -1) then
    write(stdout,*) ' '
    write(stdout,*) 'Computing the Fermi energy...'
    nelec = sh%nelec
    allocate(wk(1:kpnts%nk),isk(1:kpnts%nk))  ! kpnts%nk = total number of k-points in interp_grid
    wk(1:kpnts%nk) = 2.d0/(sh%npol*dble(kpnts%nk))   ! uniform weights
    isk(1:kpnts%nk) = 1
    is = 0
    !
    if (input%tetrahedron_method) then
        ! Tetrahedra
        num_tetra = 6 * kpnts%nk
        allocate(tetra(4,num_tetra))
        allocate(weights(kpnts%nk))
        call tetrahedra1(kpnts%nkgrid(1), kpnts%nkgrid(2), kpnts%nkgrid(3), num_tetra, tetra)
        ef = efermit (energy_sum, sh%num_bands, kpnts%nk, nelec, sh%nspin, num_tetra, tetra, is, isk)
        write(stdout,'(a,f10.5,a)') ' Fermi energy (tetrahedra) = ', ef*rytoev , ' eV'
    else
        ! Broadening methods
        ef = efermig (energy_sum, sh%num_bands, kpnts%nk, nelec, wk, input%fermi_degauss, input%fermi_ngauss, is, isk)
        write(stdout,'(a,f10.5,a)') ' Fermi energy (broadening method) = ', ef*rytoev , ' eV'
    endif
    deallocate(wk,isk)
  !
  else
  !
    write(stdout,*) 'Reading the Fermi energy from input...'
    ef = input%fermi_energy
    write(stdout,*) 'Fermi energy = ', ef*rytoev
  endif
  !!!!!!! End Fermi energy calculation
  !
  ! Occupation of the states (with Fermi-Dirac distribution)
  focc = 0.0d0
  do ik=1,energy%nk_loc
    do ii=1,sh%num_bands
      arg = (ef - energy%energy(ii,ik))/input%elec_temp
      focc_tmp(ii,ik) = wgauss(arg,-99)
    enddo
  enddo
  do ik=1,energy%nk_loc
    focc(1:sh%num_bands,ik+energy%ik_first-1) = focc_tmp(1:sh%num_bands,ik)
  enddo
  call mp_sum(focc,world_comm) 
  !
  !!! Calculate DOS  (units: states/eV)
  headline(6) = "       Energy grid [eV]               DOS"
  num_energy_dos = int( ( maxval(energy_sum(sh%num_bands,:)) - minval(energy_sum(1,:)) )&
                    & /input%delta_energy_dos + 0.5) + 1
  allocate(wgrid_dos(1:num_energy_dos),dos(1:num_energy_dos))
  wgrid_dos = 0.0d0
  do iw = 1, num_energy_dos
      wgrid_dos(iw) = minval(energy_sum(1,:)) + (iw-1) * input%delta_energy_dos ! energy grid for DOS
  enddo
  ! 
  dos = 0.d0
  if (input%tetrahedron_method) then
    ! Calculate DOS (Tetrahedron method)
    write(stdout,*) 'Computing the DOS (Tetrahedron method)...'
    do iw = 1, num_energy_dos
      w = wgrid_dos(iw)
      do ii=1, sh%num_bands
        call weights_delta1(w, energy_sum(ii,:), num_tetra, tetra, kpnts%nk, weights)
        do ik=1,kpnts%nk
          dos(iw) = dos(iw) + weights(ik)
        enddo
      enddo
    enddo
    dos(1:num_energy_dos) = 2.d0 * dos(1:num_energy_dos) / sh%npol
    if (ionode) then
      write(stdout,"(/,1x, 'Writing output on file...' )")
      call writetofile(input%prefix,"dos_tetrahedra",headline(6),num_energy_dos,wgrid_dos,1,dos/rytoev)
    endif
  else
    ! Calculate DOS (Broadening methods)
    write(stdout,*) ' '
    write(stdout,*) 'Computing the DOS (Broadening method)...'
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ik,ii,iw,w,arg) &
    !$OMP FIRSTPRIVATE(eig) REDUCTION(+:dos)
    do ik=1,energy%nk_loc
        do ii=1, sh%num_bands
           do iw = 1, num_energy_dos
             w = wgrid_dos(iw)
             !! Approximate Dirac-delta with the derivative of the Fermi-Dirac function
             arg = (energy%energy(ii, ik) - w)/input%inter_broadening
             dos(iw) = dos(iw) + w0gauss(arg,input%drude_ngauss)/input%inter_broadening
           enddo
         enddo
    enddo
    !$OMP END PARALLEL DO
    call mp_sum(dos,world_comm)
    dos(1:num_energy_dos) = 2.d0*dos(1:num_energy_dos) / (sh%npol*dble(kpnts%nk))
    if (ionode) then
        write(stdout,"(/,1x, 'Writing output on file...' )")
        call writetofile(input%prefix,"dos_broadening",headline(6),num_energy_dos,wgrid_dos,1,dos/rytoev)
    endif
  endif
  deallocate(wgrid_dos,dos,energy_sum)
  !!! End calculation of DOS
  !
  ! Calculate energy grid for dielectric function (in Ry)
  alpha = (input%wmax - input%wmin) / REAL(input%nw-1, KIND=DP)
  wgrid = 0.d0
  do iw = 1, input%nw
      wgrid(iw) = input%wmin + (iw-1) * alpha
  enddo
  !
  ! Calculate dielectric function
  jdos = 0.d0
  plasma_freq = 0.d0
  eps_im = 0.d0
  eps_re = 0.d0
  counting = 0
  write(stdout,*) ' '
  write(stdout,*) 'Computing the IP optical properties...'
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP PRIVATE(ik,q,idir,tmp,tmp2,ii,matrix_element,iband1,iband2,diff_occ,delta_e,matrix_element2,iw,w,arg) &
  !$OMP FIRSTPRIVATE(eig,commut_interp) REDUCTION(+:eps_im, eps_re, jdos, plasma_freq)
  do ik=1,energy%nk_loc
     !
     counting = counting + 1
     write(stdout,*) 'k-point:', counting
     q(1:3) = kpnts%qk(1:3,ik + energy%ik_first - 1)
     call diagonalization(q,sh,input,eig,ik,kpnts)
     energy%energy(1:sh%num_bands,ik) = eig%energy(1:sh%num_bands) ! Save bands in energy%energy
     !
     !call start_clock('optic_elements')
     do idir=1,3
        !
        ! Calculate: 1/2 \sum_ij d_i* d_j K1_ij (local matrix: <f_nk|p|f_n'k> for every n,n')
        call ZGEMM('N','N',sh%ntot_e,sh%num_bands,sh%ntot_e,(1.d0,0.d0), &
        & sh%h1(1,1,idir),sh%ntot_e,eig%wave_func,sh%ntot_e,(0.d0,0.d0),tmp,sh%ntot_e)
        call ZGEMM('C','N',sh%num_bands,sh%num_bands,sh%ntot_e,(1.d0,0.d0), eig%wave_func,sh%ntot_e,&
        & tmp,sh%ntot_e,(0.d0,0.d0),tmp2,sh%num_bands)
        do ii=1,sh%num_bands
            energy%energy_der(idir,ii,ik) = dble(tmp2(ii,ii))/2.d0 ! We need to divide by 2 (from the definition of K1_ij)
        enddo
        !
        ! Local part: (k + G) --> -i[r,H(k)] = (hbar/m)*(p + hbar*k)
        ! (The factor 2 comes from the QE units: m=1/2, hbar=1)
        ! Save <nk|p|nk> in the variable energy_der
        energy%energy_der(idir,1:sh%num_bands,ik) = 2.0*( tpiovera*q(idir) + energy%energy_der(idir,1:sh%num_bands,ik) )
        !
        matrix_element = 0.d0
        if (.not. input%nonlocal_commutator .or. .not. sh%nonlocal_commutator) then
           ! Dielectric function calculation (local part: <nk|p|n'k>)
           do iband1=1, sh%num_bands
            do iband2=1, sh%num_bands
             !
             if (iband1 == iband2) cycle
             !
             delta_e = energy%energy(iband2,ik) - energy%energy(iband1,ik)  ! Note: we define dE = E_i - E_j  while df = f_j - f_i
             if (abs(delta_e) <= delta_e_thr/rytoev ) cycle    ! We do not consider transitions with very small dE
             diff_occ = focc(iband1,ik + energy%ik_first - 1) - focc(iband2,ik + energy%ik_first - 1)
             !
             matrix_element(iband1,iband2) = tmp2(iband1,iband2)
             !matrix_element2 = |<nk|p|n'k>|^2 / (E_nk - E_n'k)^2
             matrix_element2 = dble( matrix_element(iband1,iband2) * conjg(matrix_element(iband1,iband2)) ) / delta_e**2
             do iw = 1, input%nw
                   w = wgrid(iw)
                   arg = (delta_e - w)/input%inter_broadening
                   eps_re(iw,idir) = eps_re(iw,idir) + diff_occ*matrix_element2* ( &
                                    & (w - delta_e)/( (w-delta_e)**2 + input%inter_broadening**2 ) )
                   eps_im(iw,idir) = eps_im(iw,idir) + diff_occ*matrix_element2* &
                                    & lorentzian(arg*input%inter_broadening, input%inter_broadening)
                   ! NOTE: with the formula above both resonant and anti-resonant terms are included
                   !
                   if (idir == 1) then
                     ! JDOS calculation (equal to eps_im with matrix_element=1)
                     jdos(iw) = jdos(iw) + diff_occ*w0gauss(arg,input%drude_ngauss)/input%inter_broadening
                     !jdos(iw) = jdos(iw) + diff_occ * lorentzian(arg*input%inter_broadening, input%inter_broadening) / delta_e**2
                   endif
                   !
             enddo
             !
            enddo
           enddo
           ! End dielectric function calculation (local part)
           !
        ! Non-local matrix: <f_nk|[r,V_nl]|f_n'k>
        elseif (input%nonlocal_commutator .and. sh%nonlocal_commutator) then
             !
             matrix_element(1:sh%num_bands,1:sh%num_bands) = tmp2(1:sh%num_bands,1:sh%num_bands) ! store local contribution
             !
             if (input%nonlocal_interpolation) then
               call trilinear_parallel_commut(kpnts,ik,sh,commut_interp)
               call ZGEMM('N','N',sh%ntot_e,sh%num_bands,sh%ntot_e,(1.d0,0.d0), &
               & commut_interp(1,1,idir),sh%ntot_e,eig%wave_func,sh%ntot_e,(0.d0,0.d0),tmp,sh%ntot_e)
               call ZGEMM('C','N',sh%num_bands,sh%num_bands,sh%ntot_e,(1.d0,0.d0), eig%wave_func,sh%ntot_e,&
               & tmp,sh%ntot_e,(0.d0,0.d0),tmp2,sh%num_bands)
             else
               call ZGEMM('N','N',sh%ntot_e,sh%num_bands,sh%ntot_e,(1.d0,0.d0), &
               & sh%commut(1,1,idir,ik),sh%ntot_e,eig%wave_func,sh%ntot_e,(0.d0,0.d0),tmp,sh%ntot_e)
               call ZGEMM('C','N',sh%num_bands,sh%num_bands,sh%ntot_e,(1.d0,0.d0), eig%wave_func,sh%ntot_e,&
               & tmp,sh%ntot_e,(0.d0,0.d0),tmp2,sh%num_bands)
             endif
             !
             ! Add the non-local contribution: energy_der = <nk|p|nk> + <nk|[r,V_nl]|nk>
             do ii=1,sh%num_bands
                 energy%energy_der(idir,ii,ik) = energy%energy_der(idir,ii,ik) + dble(tmp2(ii,ii))
             enddo
             !
             ! Dielectric function calculation (non-local part)
             do iband1=1, sh%num_bands
              do iband2=1, sh%num_bands
                 !
                 if (iband1 == iband2) cycle
                 !
                 delta_e = energy%energy(iband2,ik) - energy%energy(iband1,ik)
                 if (abs(delta_e) <= delta_e_thr/rytoev ) cycle    ! We do not consider transitions with very small dE
                 diff_occ = focc(iband1,ik + energy%ik_first - 1) - focc(iband2,ik + energy%ik_first - 1)
                 !
                 ! matrix_element =  <nk|p|n'k> + <nk|[r,V_nl]|n'k>
                 matrix_element(iband1,iband2) = matrix_element(iband1,iband2) + tmp2(iband1,iband2)
                 matrix_element2 = dble( matrix_element(iband1,iband2) * conjg(matrix_element(iband1,iband2)) )  / delta_e**2
                 !
                 do iw = 1, input%nw
                   w = wgrid(iw)
                   arg = (delta_e - w)/input%inter_broadening
                   eps_re(iw,idir) = eps_re(iw,idir) + diff_occ*matrix_element2* ( &
                                    & (w - delta_e)/( (w-delta_e)**2 + input%inter_broadening**2 ) )
                   eps_im(iw,idir) = eps_im(iw,idir) + diff_occ*matrix_element2* &
                                    & lorentzian(arg*input%inter_broadening, input%inter_broadening)
                   ! NOTE: with the formula above both resonant and anti-resonant terms are included
                   !
                   if (idir == 1) then
                     ! JDOS calculation (eps_im with matrix_element=1)
                     jdos(iw) = jdos(iw) + diff_occ*w0gauss(arg,input%drude_ngauss)/input%inter_broadening
                     !jdos(iw) = jdos(iw) + diff_occ * lorentzian(arg*input%inter_broadening, input%inter_broadening) / delta_e**2
                   endif
                   !
                 enddo
                 !
               enddo
             enddo
             ! End dielectric function calculation (non-local part)
             !
        endif
        !
     enddo
     !
     ! Calculate the Drude plasma frequency
     do ii=1,sh%num_bands
            ! w0gauss = 1.0d0 / (2.0d0 + exp ( - x) + exp ( + x) ) for Fermi-Dirac smearing (ngauss=-99)
            ! dirac_delta = -df/dE = w0gauss/degauss
            ! 1 = xx
            ! 2 = yy
            ! 3 = zz
            ! 4 = xy
            ! 5 = xz
            ! 6 = yz
            arg = (ef - energy%energy(ii,ik))/input%drude_degauss
            plasma_freq(1) = plasma_freq(1) + energy%energy_der(1,ii,ik)**2*w0gauss(arg, &
            & input%drude_ngauss)/input%drude_degauss
            plasma_freq(2) = plasma_freq(2) + energy%energy_der(2,ii,ik)**2*w0gauss(arg, &
            & input%drude_ngauss)/input%drude_degauss
            plasma_freq(3) = plasma_freq(3) + energy%energy_der(3,ii,ik)**2*w0gauss(arg, &
            & input%drude_ngauss)/input%drude_degauss
            plasma_freq(4) = plasma_freq(4) + energy%energy_der(1,ii,ik)*energy%energy_der(2,ii,&
            & ik)*w0gauss(arg,input%drude_ngauss)/input%drude_degauss
            plasma_freq(5) = plasma_freq(5) + energy%energy_der(1,ii,ik)*energy%energy_der(3,ii,&
            & ik)*w0gauss(arg,input%drude_ngauss)/input%drude_degauss
            plasma_freq(6) = plasma_freq(6) + energy%energy_der(2,ii,ik)*energy%energy_der(3,ii,&
            & ik)*w0gauss(arg,input%drude_ngauss)/input%drude_degauss
     enddo
     ! End calculation of Drude plasma frequency
     !call stop_clock('optic_elements')
     !
  enddo
  !$OMP END PARALLEL DO
  !
  call mp_sum(plasma_freq,world_comm)
  call mp_sum(jdos,world_comm)
  call mp_sum(eps_im,world_comm)
  call mp_sum(eps_re,world_comm)
  !
  ! Drude plasma frequency in eV
  plasma_freq(1:6) = sqrt( 16.d0 * pi * rytoev**2 * abs(plasma_freq(1:6)) / ( sh%npol * dble(kpnts%nk) * sh%omega ) )
  !
  write(stdout,*) ' '
  write(stdout,'(a,f10.5,a)') ' Drude plasma frequency (xx) = ', plasma_freq(1) , ' eV'
  write(stdout,'(a,f10.5,a)') ' Drude plasma frequency (yy) = ', plasma_freq(2) , ' eV'
  write(stdout,'(a,f10.5,a)') ' Drude plasma frequency (zz) = ', plasma_freq(3) , ' eV'
  write(stdout,'(a,f10.5,a)') ' Drude plasma frequency (xy) = ', plasma_freq(4) , ' eV'
  write(stdout,'(a,f10.5,a)') ' Drude plasma frequency (xz) = ', plasma_freq(5) , ' eV'
  write(stdout,'(a,f10.5,a)') ' Drude plasma frequency (yz) = ', plasma_freq(6) , ' eV'
  !
  ! JDOS/w^2
  ! If the matrix elements are ~ 1 --> eps_2(w) ~ JDOS(w)/w^2
  jdos(1:input%nw) = 16.d0 * pi**2 * jdos(1:input%nw) / (sh%npol * dble(kpnts%nk) * sh%omega ) / rytoev
  jdos(1:input%nw) = jdos(1:input%nw) / wgrid(1:input%nw)**2
  !
  ! Dielectric function (interband)
  eps_im(1:input%nw,1:3) = 16.d0 * pi**2 * eps_im(1:input%nw,1:3) / (sh%npol * dble(kpnts%nk) * sh%omega )
  eps_re(1:input%nw,1:3) = 1.d0 - 16.d0 * pi * eps_re(1:input%nw,1:3) / (sh%npol * dble(kpnts%nk) * sh%omega )
  !
  ! Total dielectric function (interband + Drude)
  eps_tot_im = 0.d0
  eps_tot_re = 0.d0
  do idir=1,3
    do iw = 1, input%nw
      w = wgrid(iw)
      eps_tot_im(iw,idir) = eps_im(iw,idir) + (plasma_freq(idir)**2/rytoev**2) &
                        & * input%intra_broadening / (w*(w**2 + input%intra_broadening**2))
      eps_tot_re(iw,idir) = eps_re(iw,idir) - (plasma_freq(idir)**2/rytoev**2) &
                        & / (w**2 + input%intra_broadening**2)
    enddo
  enddo
  !
  ! EELS
  eels = 0.d0
  do idir=1,3
    do iw = 1, input%nw
      eels(iw,idir) = eps_tot_im(iw,idir) / (eps_tot_im(iw,idir)**2 + eps_tot_re(iw,idir)**2)
    enddo
  enddo
  !
  ! Average total dielectric function (average over x, y, z)
  eps_avg_im = 0.d0
  eps_avg_re = 0.d0
  do iw = 1, input%nw
      eps_avg_im(iw) = ( eps_tot_im(iw,1) + eps_tot_im(iw,2) + eps_tot_im(iw,3)  ) / 3.d0
      eps_avg_re(iw) = ( eps_tot_re(iw,1) + eps_tot_re(iw,2) + eps_tot_re(iw,3)  ) / 3.d0
  enddo
  !
  ! Refractive index (from the average dielectric function)
  refractive_index = 0.d0
  extinct_coeff = 0.d0
  do iw = 1, input%nw
      norm_epsilon = sqrt(eps_avg_re(iw)**2 + eps_avg_im(iw)**2)
      refractive_index(iw) = sqrt( ( eps_avg_re(iw) + norm_epsilon  ) / 2.d0 )
      extinct_coeff(iw) = sqrt( ( -eps_avg_re(iw) + norm_epsilon  ) / 2.d0 )
  enddo
  !
  ! Reflectivity
  reflectivity = 0.d0
  do iw = 1, input%nw
      reflectivity(iw) = ( (refractive_index(iw) - 1.d0)**2 + extinct_coeff(iw)**2 ) / &
                        & ( (refractive_index(iw) + 1.d0)**2 + extinct_coeff(iw)**2 )
  enddo
  !
  ! Write outputs to file
  headline(1) = "       Energy grid [eV]               Re(eps)_x                Re(eps)_y                Re(eps)_z"
  headline(2) = "       Energy grid [eV]               Im(eps)_x                Im(eps)_y                Im(eps)_z"
  headline(3) = "       Energy grid [eV]               JDOS/w^2"
  headline(4) = "       Energy grid [eV]               EELS_x                   EELS_y                   EELS_z"
  headline(5) = "       Energy grid [eV]               Re(eps)                  Im(eps)               refractive_index &
         &         extinct_coeff             reflectivity"
  !
  if (ionode) then
    write(stdout,"(/,1x, 'Writing output on file...' )")
    call writetofile(input%prefix,"eps_inter_re",headline(1),input%nw,wgrid,3,eps_re)
    call writetofile(input%prefix,"eps_inter_im",headline(2),input%nw,wgrid,3,eps_im)
    call writetofile(input%prefix,"jdos",headline(3),input%nw,wgrid,1,jdos)
    call writetofile(input%prefix,"eels",headline(4),input%nw,wgrid,3,eels)
    call writetofile(input%prefix,"eps_tot_re",headline(1),input%nw,wgrid,3,eps_tot_re)
    call writetofile(input%prefix,"eps_tot_im",headline(2),input%nw,wgrid,3,eps_tot_im)
  endif
  !
  if (ionode) then
    str = TRIM(input%prefix) // "." // TRIM("optical_constants") // ".dat"
    iun=find_free_unit()
    open(iun,FILE=TRIM(str))
    !
    write(iun,"(a)") "# "// TRIM(headline(5))
    write(iun,"(a)") "#"
    !
    do iw = 1, input%nw
       write(iun,"(10f25.6)") wgrid(iw)*rytoev, eps_avg_re(iw), eps_avg_im(iw), refractive_index(iw), &
                              & extinct_coeff(iw), reflectivity(iw)
    enddo
    !
    close(iun)
    !
    write(stdout,*) 'File ', TRIM(str), ' written'
  endif
  !
  if (input%tetrahedron_method) then
    deallocate(tetra,weights)
  endif
  deallocate(tmp,tmp2)
  if (input%nonlocal_commutator .and. sh%nonlocal_commutator) then
    deallocate(commut_interp)
  endif
  deallocate(focc,focc_tmp,wgrid,jdos,eps_im,eps_re,eps_tot_re,eps_tot_im)
  deallocate(eels,eps_avg_im,eps_avg_re,refractive_index,extinct_coeff,reflectivity)
  deallocate(matrix_element)
  call deallocate_eigen(eig)
  !
  write(stdout,*) ' '
  write(stdout,*) '*************************************'
  write(stdout,*) 'Optical properties computed and saved'
  write(stdout,*) '*************************************'
  !
  ! Write in output the input parameters
  write(stdout,*) '                     '
  write(stdout,*) '                     '
  write(stdout,*) 'INPUT PARAMETERS:'
  write(stdout,*)           'interpolation k-grid = ', input%interp_grid(1:3)
  if (.not. input%nonlocal_commutator .or. .not. sh%nonlocal_commutator) then
    write(stdout,*)           'nonlocal_commutator =    False'
  elseif (input%nonlocal_commutator .and. sh%nonlocal_commutator) then
    write(stdout,*)           'nonlocal_commutator =    True'
  endif
  if (input%tetrahedron_method) then
    write(stdout,*)           'tetrahedron_method =    True'
  endif
  if (input%nonlocal_interpolation) then
    write(stdout,*)           'nonlocal_interpolation =    True'
  endif
  write(stdout,'(a,f10.5)') ' fermi_degauss [eV] = ', input%fermi_degauss*rytoev
  write(stdout,*)           'fermi_ngauss = ',  input%fermi_ngauss
  write(stdout,'(a,f10.5)') ' drude_degauss [eV] = ', input%drude_degauss*rytoev
  write(stdout,'(a,f10.5)') ' elec_temp [eV] = ', input%elec_temp*rytoev
  write(stdout,'(a,f10.5)') ' inter broadening [eV] = ', input%inter_broadening*rytoev
  write(stdout,'(a,f10.5)') ' intra broadening [eV] = ', input%intra_broadening*rytoev
  write(stdout,'(a,f10.5)') ' s_bands [a.u.] = ', sh%s_bands
  write(stdout,*) ' '
  !
  !call stop_clock('dielectric')
  !
end subroutine dielectric

subroutine writetofile(prefix,name,headline,nw,wgrid,ncolumn,var)
  !------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE constants, ONLY : rytoev
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*),   INTENT(IN) :: prefix
  CHARACTER(LEN=*),   INTENT(IN) :: name
  CHARACTER(LEN=*),   INTENT(IN) :: headline
  INTEGER,            INTENT(IN) :: nw, ncolumn
  REAL(kind=DP),           INTENT(IN) :: wgrid(nw)
  REAL(kind=DP),           INTENT(IN) :: var(nw,ncolumn)
  !
  CHARACTER(256) :: str
  INTEGER        :: iw, iun
  INTEGER, EXTERNAL :: find_free_unit


  str = TRIM(prefix) // "." // TRIM(name) // ".dat"
  iun=find_free_unit()
  open(iun,FILE=TRIM(str))
  !
  write(iun,"(a)") "# "// TRIM(headline)
  write(iun,"(a)") "#"
  !
  do iw = 1, nw
     write(iun,"(10f25.6)") wgrid(iw)*rytoev, var(iw,1:ncolumn)
  enddo
  !
  close(iun)
  !
  WRITE(STDOUT,*) 'File ', TRIM(str), ' written'
end subroutine writetofile




function lorentzian (x, smr)
  !-----------------------------------------------------------------------
  !
  !     this function computes the Lorentzian function at the point x with
  !     smearing smr.
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE io_global, ONLY : stdout
  implicit none

  REAL(kind=DP) :: lorentzian
  REAL(kind=DP) :: x, smr
  ! output: the value of the function   (lorentzian)
  ! input: the argument of the function (x)
  ! input: the smearing of the function (smr)

  lorentzian = smr / (  pi*(  x**2 + smr**2  ) )
  return
end function lorentzian



