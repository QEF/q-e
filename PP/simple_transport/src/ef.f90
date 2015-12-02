      program efermi
      !---------------------------------------------------------------------------------
      ! Written by Burak Himmetoglu (burakhmmtgl@gmail.com)
      ! Uses some parts of the PW distribution
      !
      ! This is a simple code that reads the electron concentration and band
      ! structure (from a2Fsave file) and then computes the Fermi level. 
      ! The system must be insulating.
      !  
      ! Description of the input card:
      !
      ! fil_a2F            : prefix.a2Fsave file
      ! fil_info           : file containing direct and reciprocal lattice vectors
      ! T                  : Temperature in K 
      ! vol                : Volume in au^3
      ! cbm, vbm           : Conduction band min. and Valence band max. (in eV)
      ! doping             : Electron concentration (in cm^-3). For hole doping, use negative value
      ! ndop               : Number of doping levels considered 
      ! nthreads           : Number of threads for OpenMP parallelization
      ! lscissors          : if .true. conduction bands are shifted by a constant energy up
      ! shift              : value of the shift applied to conduction bands
      ! cbm_i              : the initial conduction band for which the shift is applied ( all bands with index >= cbm_i are shifted up)

!$    use omp_lib

      implicit none
      !
      integer :: i,j,iq,ik,imu,itemp,nu,nqtot,nsig,nat,nk1fit,nk2fit,   &
     &           nk3fit, nkfit, nksfit_real, nbnd, nksfit, npk, nsym,   &
     &           s(3,3,48),ns,nrot,ibnd,ndop,cbm_i
      ! 
      ! OMP 
      integer :: TID, nthreads
      double precision :: t0

      !
      double precision :: at(3,3), bg(3,3),vol,T, fd, nelec1,           &
     &                    cbm, vbm, ef_mid, doping(10), nelec,          &
     &                    Ef_IBZ(10), shift
      ! 
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                                 wk(:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      logical :: lscissors
      !
      character*20 :: fil_kp, fil_a2F, fil_info
      !
      double precision, PARAMETER :: Rytocm1 = 109737.57990d0,          &
     &                               RytoGHz = 3.289828D6,              &
     &                               RytoTHz = RytoGHz/1000.d0,         &
     &                               RytoeV = 13.60569253d0,            &
     &                               tpi = 6.283185307d0,               &
     &                               convsig = 2.89d5,                  &
     &                               KtoRy = 1.d0/38.681648/300/RytoeV, &
     &                               autocm = 5.2917721092d-9,          &
     &                               convRH = 9.2522431d-13
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 

      namelist /input/ fil_info, fil_a2F, vol, cbm, vbm, T,             &
     &                 doping, ndop, nthreads, lscissors, shift, cbm_i

      read(5,input)

      t0 = 0.0

      !$ call omp_set_num_threads(nthreads)

      !$ t0 = omp_get_wtime()
      npk = 40000
      !
      ! Convert temperature from K to Ryd
      T = T * KtoRy 
      !
      ! Read a2F
      !
      open(11,file=fil_a2F,status='unknown')
      !
      read(11,*) nbnd, nksfit
      !
      allocate(etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
      !
      ! Read band structure and k-point data
      read(11,*) etfit
      read(11,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
      read(11,*) wkfit
      read(11,*) nk1fit, nk2fit, nk3fit
      read(11,* ) nsym
      do ns=1,nsym
         read(11,*)  ((s(i,j,ns),j=1,3),i=1,3)
      enddo
      !
      close(11)
      !
      ! Regular grid size
      nkfit=nk1fit*nk2fit*nk3fit
      ! Read info file on k-points (and lattice)
      !
      open(11,file=fil_info,status='unknown')
      !
      read(11,*)
      read(11,*) ((at(i,j),i=1,3),j=1,3)
      !
      read(11,*)
      read(11,*)
      !
      read(11,*) ((bg(i,j),i=1,3),j=1,3)
      ! 
      close(11)
      !
      allocate (eqkfit(nkfit), sfit(nkfit),wk(nkfit))
      !         eqkfit : pointers to band energies in uniform grid
      !         sfit   : pointers to symmetries
      ! 
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      ! 
      ! eqkfit(nk) : maps IBZ to full grid. The full grid is in crystal coords
      ! 
      ! If lscissors is true, then shift band energies (move conduction bands higher in energy)
      if (lscissors) then
        etfit(cbm_i:nbnd,:) = etfit(cbm_i:nbnd,:) + shift/RytoeV
        cbm = cbm + shift
      end if
      !
      ! Determine the number of electrons (neutral, insulating system)
      ef_mid = (cbm + vbm)/2.0/RytoeV ! Fermi level at the middle of the gap
      !
      nelec1 = 0.d0
      !
      ! Irreducable BZ
      nelec1 = sumk (etfit,nbnd,nksfit,wkfit,T,1,ef_mid) 
      !
      do i=1,ndop
         ! Determine number of electrons from doping levels (given in cm-3)
         nelec = doping(i) * vol * autocm**3 + nelec1
         ! 
         ! Compute efermi using the bisection method
         Ef_IBZ(i) = fermi_en (etfit,wkfit,nbnd,nksfit,nelec,T,1) 
         ! 
      end do
      !
      ! Memory clean 
      !
      deallocate (xkfit,etfit,wkfit,wk,eqkfit,sfit)
      !
      !$t0 = omp_get_wtime() - t0
      !
      if (t0 >0) write(6,"(A,e14.6)") 'Walltime= ', t0
      ! 
      do i=1,ndop
         write(6,"(A,e14.6,f14.6)") 'doping, Efermi', doping(i), Ef_IBZ(i) * RytoeV
      end do
      !
      !
      contains
      ! -----------------------------------------------
      ! FUNCTIONS
      ! -----------------------------------------------

      double precision function sumk (et,nbnd,nks,wk,degauss,ngauss,ee)

!$    use omp_lib

      implicit none

      integer, intent(in) :: nbnd, nks, ngauss
      double precision, intent(in) :: ee,degauss
      ! wk: weight of kpt
      ! ee: E
      ! degauss: temperature/degauss
      double precision, intent(in) :: et(nbnd,nks), wk(nks)

      double precision :: fd, x, arg, a, hp, hd, sqrtpm1

      integer :: i, j, ik, n, ni

      sqrtpm1 = 1.0d0/1.77245385090551602729d0

      sumk = 0.d0
   
      n = 1 ! First order Gauss-Hermite

      if (ngauss .eq. 1) then ! FD
         !$omp parallel do      &
         !$omp default(shared)  &
         !$omp private(i,ik,fd) &
         !$omp reduction(+ : sumk)
         do i=1,nbnd
            do ik=1,nks
               fd = 1.d0 / ( 1.d0+exp( (et(i,ik)-ee)/degauss ) )
               sumk = sumk + wk(ik) * fd
            end do
         end do
         !$omp end parallel do
      else if (ngauss .eq. 0) then ! Gaussian (for zero temperature)
         !$omp parallel do      &
         !$omp default(shared)  &
         !$omp private(i,ik,fd,x) &
         !$omp reduction(+ : sumk)
         do i=1,nbnd
            do ik=1,nks
               x = (ee-et(i,ik))/degauss/dsqrt(2.d0)
               fd = 0.5d0 * (1.d0 + derf( x ) )
               sumk = sumk + wk(ik) * fd 
            end do
         end do
         !$omp end parallel do
      else if ( ngauss .eq. 2 ) then ! MP
         !$omp parallel do          &
         !$omp default(shared)      &
         !$omp private(i,j,ik,fd,x) &
         !$omp reduction(+ : sumk)
         do i=1,nbnd
            do ik=1,nks
               ! 
               x = (ee-et(i,ik))/degauss/dsqrt(2.d0)
               fd = 0.5d0 * (1.d0 + derf(x) )
               hd = 0.d0
               arg = min (200.d0, x**2)
               hp = exp ( - arg)
               ni = 0
               a = sqrtpm1
               !
              do j = 1, n
                 hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
                 ni = ni + 1
                 a = - a / (DBLE (j) * 4.0d0)
                 fd = fd - a * hd
                 hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
                 ni = ni + 1
              end do
              ! 
              sumk = sumk + wk(ik) * fd
            end do
         end do
         !$omp end parallel do

      end if

      end function sumk 
      !      
      double precision function fermi_en (eig,wk,nbnd,nktot,nelec,temperature,ngauss)

      implicit none

      integer, intent(in) :: nbnd, nktot, ngauss
      double precision, intent(in) :: eig(nbnd,nktot),wk(nktot), nelec, &
     &                                temperature

      ! Local variables
      integer :: i, ik
      double precision :: Elw, Eup, Ef, sumkup, sumklw, sumkmid
      integer, PARAMETER :: maxiter = 300
      double precision, PARAMETER :: eps = 1d-10
 
      Elw = eig(1,1)
      Eup = eig(nbnd,1)
      do ik = 2, nktot
         Elw = min ( Elw, eig (1, ik) ) 
         Eup = max ( Eup, eig (nbnd, ik) ) 
      end do
      Eup = Eup + 2*temperature
      Elw = Elw - 2*temperature

      sumkup = sumk(eig, nbnd, nktot, wk, temperature, ngauss, Eup)
      sumklw = sumk(eig, nbnd, nktot, wk, temperature, ngauss, Elw)
      !
      if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps ) then
         fermi_en = -99
         return
      end if 
      !
      do 100 i=1, maxiter
         Ef = 0.5d0 * (Eup + Elw)
         sumkmid = sumk(eig, nbnd, nktot, wk, temperature, ngauss, Ef)
         if ( abs (sumkmid-nelec) < eps ) then
            fermi_en = Ef
            exit
         else if ( (sumkmid-nelec) < -eps) then
            Elw = Ef
         else
            Eup = Ef
         end if
100   continue

      end function fermi_en
      !
      end program efermi
