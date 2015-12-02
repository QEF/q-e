      program dos
      !
      ! Written by Burak Himmetoglu (burakhmmtgl@gmail.com)
      ! Uses some parts of the PW distribution
      !
      ! Computation of DOS using the variable smearing technique. 
      ! See: J. R. Yates, X. Wang, D. Vanderbilt, and I. Souza, PRB 75, 195121 (2007)
      ! 
      ! Description of the input card:
      !
      ! fil_a2F    : prefix.a2Fsave file
      ! fil_info   : file containing direct and reciprocal lattice vectors
      ! T          : Temperature in K 
      ! Nen        : Number of energy values.
      ! emin, emax : Min and Max energies.
      ! lsoc       : if .true. then the band structure is non-collinear 
      ! lfix       : if .true. DOS with constant smearing is also computed.
      ! degauss    : if lsoc is .true. this is the constant smearing value in eV. 
      ! cbm, vbm   : Conduction band min. and Valence band max. (in eV)
      ! aa         : Parameter determining the k-dependent smearing deg_nk : deg_nk = aa * delta_k * (de_nk/delta_k). 
      !              Larger aa will lead to smoother DOS. Needs to be tested for various grids.
      ! nthreads   : Number of threads for OpenMP parallelization.
      ! lscissors          : if .true. conduction bands are shifted by a constant energy up
      ! shift              : value of the shift applied to conduction bands
      ! cbm_i              : the initial conduction band for which the shift is applied ( all bands with index >= cbm_i are shifted up)

!$    USE omp_lib
      USE smearing_mod

      implicit none
      !
      !
      integer :: i,j,iq,ik,ie,nk1fit,nk2fit,nk3fit,nkfit,               &
     &           nbnd, nksfit, npk, nsym,                               &
     &           s(3,3,48),ns,nrot,ibnd,io,Nen, ind_k, cbm_i
      ! 
      ! OpenMP
      integer :: TID, nthreads
      double precision :: t0

      !
      double precision :: at(3,3), bg(3,3),deg,degauss,wk,T,            &
     &                    cbm, vbm, ef_mid, doping, nelec, shift,       &
     &                    en(1000), de(1000), ne(1000),emin, emax, aa
      ! 
      !
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                                 vk(:,:,:), dfk(:,:,:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      logical :: lsoc, lfix, lscissors
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
                                     degmin = 1e-3                     
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 

      namelist /input/ fil_info, fil_a2F, Nen, aa, T, degauss,          &
     &                 cbm, vbm, emin, emax, nthreads, lsoc, lfix,      &
     &                 cbm_i, lscissors, shift

      read(5,input)

      !Set number of threads
      !$ call omp_set_num_threads(nthreads)

      t0 = 0.0

      npk = 40000
      !
      ! Energy grid
      do i=1,Nen
         en(i) = (emax-emin)/(Nen-1) * (i-1) + emin
      end do
      ! 
      ! Convert to Ryd
      en = en / RytoeV
      T  = T * KtoRy
      ! 
      ! Read a2Fsave
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
      ! 
      ! Read info file on k-points (and lattice)
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
      allocate(eqkfit(nkfit),sfit(nkfit))
      allocate(dfk(nbnd,nkfit,3),vk(nbnd,nkfit,3))
      !         eqkfit : pointers to band energies in uniform grid
      !         sfit   : pointers to symmetries
      ! 
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      ! 
      ! eqkfit(nk) : maps IBZ to full grid. The full grid is in crystal coords
      !
      ! Deallocate unnecessary variables
      deallocate(wkfit,sfit,xkfit)
      !
      ! Weights for regular grid (include lsoc)
      if ( lsoc .eqv. .true.) then 
         wk = 1.0/nkfit
      else
         wk = 2.0/nkfit
      end if
      !
      ! If lscissors is true, then shift band energies (move conduction bands higher in energy)
      if (lscissors) then
        etfit(cbm_i:nbnd,:) = etfit(cbm_i:nbnd,:) + shift/RytoeV
        cbm = cbm + shift
      end if
      !
      !Start timing here
      !$ t0 = omp_get_wtime()
      !
      ! Call band velocities and forward derivatives
      call vband_ibz( nk1fit,nk2fit,nk3fit,nbnd,nksfit,etfit,eqkfit,at, vk, dfk)
      !
      !Check number of electrons
      nelec = 0.d0
      ef_mid = (cbm+vbm)/2.d0/RytoeV ! Fermi level set to midgap
      !
      !$omp parallel do default(shared) &
      !$omp collapse(2) &
      !$omp private(ik,ibnd,ind_k,deg) &
      !$omp reduction(+: nelec)
      do ibnd=1,nbnd
         !
         do ik=1,nkfit
            !
            ind_k = eqkfit(ik)
            deg = sig0(nk1fit,nk2fit,nk3fit,dfk(ibnd,ik,:),aa)
            nelec = nelec + wk * w1gauss((ef_mid-etfit(ibnd,ind_k))/deg,0)
         end do
      end do
      !$omp end parallel do
      !
      !Compute density of states and number of electrons at a given energy (en(i))
      de = 0.d0
      ne = 0.d0
      do ie=1,Nen 
         !
         !$omp parallel do default(shared) &
         !$omp collapse(2) &
         !$omp private(ik,ibnd,deg,ind_k) &
         !$omp reduction(+: de, ne)  
         do ik=1,nkfit
            !
            do ibnd=1,nbnd
               !
               ind_k = eqkfit(ik)
               deg = sig0(nk1fit,nk2fit,nk3fit,dfk(ibnd,ik,:),aa)
               !
               de(ie) = de(ie) + wk * w0gauss((etfit(ibnd,ind_k)-en(ie))/deg)/deg
               ne(ie) = ne(ie) + wk * w1gauss((etfit(ibnd,ind_k)-en(ie))/deg,0)
               !
            end do ! ibnd
         end do ! ik
         !$omp end parallel do
      end do ! ie
      !
      ! End timing
      !$ t0 = omp_get_wtime() - t0
      !
      open(11,file='dos.out',status='unknown')
      if (t0 >0) then
         write(11,"(A,I5)") "#Number of threads = ", nthreads
         write(11,"(A,e14.6)") "#Integration time(s) = ", t0
      end if
      !
      do ie=1,Nen
         write(11,"(3f14.6)") en(ie)*RytoeV, de(ie)/RytoeV, ne(ie)
      end do
      close(11)
      !
      if (lfix .eqv. .true.) then
         !Compute density of states and number of electrons at a given energy (en(i)) for 
         !fixed smearing of degauss
         de = 0.d0
         ne = 0.d0
         deg = degauss / RytoeV
         do ie=1,Nen
            !
            !$omp parallel do default(shared) &
            !$omp collapse(2) &
            !$omp private(ik,ibnd,ind_k) &
            !$omp reduction(+: de,ne)  
            do ik=1,nkfit
               !
               do ibnd=1,nbnd
                  !
                  ind_k = eqkfit(ik)
                  de(ie) = de(ie) + wk * w0gauss((etfit(ibnd,ind_k)-en(ie))/deg)/deg
                  ne(ie) = ne(ie) + wk * w1gauss((etfit(ibnd,ind_k)-en(ie))/deg,0)
                  !
               end do ! ibnd
            end do ! ik
            !$omp end parallel do
         end do ! ie
         !
         open(11,file='dos.fix.out',status='unknown')
         do ie=1,Nen
            write(11,"(3f14.6)") en(ie)*RytoeV, de(ie)/RytoeV, ne(ie)
         end do
         close(11)
      end if ! lfix
      !
      ! Free memory  
      deallocate(eqkfit,etfit,vk,dfk)
      ! 
      end program dos
