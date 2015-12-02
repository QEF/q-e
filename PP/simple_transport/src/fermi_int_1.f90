      program fermi_int_1
      !
      ! Written by Burak Himmetoglu (burakhmmtgl@gmail.com)
      ! Uses some parts of the PW distribution
      !
      ! This program computes the transport integrals at a given Fermi level using  
      ! constant scattering rate, and results in conductivity and Seebeck coefficients.
      !
      ! Integrals that are computed:
      !
      ! I0 = \sum_{n,k} (-df_nk/de) 
      ! I1(i,j) = \sum_{n,k} tau (-df_nk/de) v_nk,i v_nk,j
      ! I2(i,j) = \sum_{n,k} tau (-df_nk/de) (e_nk - e_F) v_nk,i v_nk,j
      !
      ! The difference between fermi_int_0 is that the calculation proceeds in a reduced
      ! grid such that |e_nk-Ef| <= cut*kT. This allows faster integration times for a 
      ! given Fermi level (Ef).
      !
      ! Description of the input card:
      !
      ! fil_a2F            : prefix.a2Fsave file
      ! fil_info           : file containing direct and reciprocal lattice vectors
      ! T                  : Temperature in K 
      ! alat               : lattice parameter in a.u. (Bohr)
      ! vol                : volume of lattice in a.u.^3
      ! cut                : Reduction parameter for integrals (i.e. Integration reduced to the region satisfying |e_nk-Ef| <= cut*kT)
      ! invtau             : inverse tau in eV
      ! phband_i, phband_f : starting and ending indices of bands of interest. The bands must lie between efmin & efmax
      ! lsoc               : if .true. then the band structure is non-collinear
      ! nthreads           : Number of threads for OpenMP parallelization
      ! lscissors          : if .true. conduction bands are shifted by a constant energy up
      ! shift              : value of the shift applied to conduction bands
      ! cbm_i              : the initial conduction band for which the shift is applied ( all bands with index >= cbm_i are shifted up)

!$    USE omp_lib        

      implicit none

      integer :: i,j,k,nu,ik,ikk,nk1fit,nk2fit,nk3fit,nkfit,            &
     &           nbnd, nksfit, npk, nsym, Nmu, imu,                     &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f,          &
     &           nphband, n, nn, jbnd, ibnd_ph, ind_k, cbm_i
      !
      double precision :: wk, at(3,3), bg(3,3), efermi, alat,           &
     &                    T, wo(3), al(3), invtau,aa,cut,deg,           &
     &                    invT, fd, dfd, fac, vol, shift
      ! 
      logical :: lsoc, lscissors
      !
      double precision, allocatable :: xkfit(:,:),etfit(:,:),wkfit(:),  &
     &                                 dfk(:,:,:),vk(:,:,:)
      !
      double precision :: I0,I1(3,3),I2(3,3),sig(3,3),Se(3,3),inv_I1(3,3)
      !
      integer, allocatable :: eqkfit(:), sfit(:), iflag(:,:),nkeff(:)
      !
      ! OMP variables
      double precision :: t0
      !
      integer :: nthreads
      !
      character*20 :: fil_a2F, fil_info
      !
      double precision, PARAMETER :: Rytocm1 = 109737.57990d0,          &
     &                               RytoGHz = 3.289828D6,              &
     &                               RytoTHz = RytoGHz/1000.d0,         &
     &                               RytoeV = 13.60569253d0,            &
     &                               tpi = 6.283185307d0,               &
     &                               convsig = 2.89d5,                  &
     &                               KtoRy = 1.d0/38.681648/300/RytoeV
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 
      !
      namelist /input/ fil_info, fil_a2F, T, phband_i, phband_f,cut,    &
     &                 efermi, invtau, alat, vol, nthreads, lsoc,       &
     &                 cbm_i, shift, lscissors

      read(5,input)

      !Set number of threads
      !$ call omp_set_num_threads(nthreads)

      t0 = 0.0

      npk = 40000
      !
      ! Convert from eV to Ryd
      T = T * KtoRy
      invtau = invtau / RytoeV
      efermi = efermi / RytoeV
      !
      ! Total number of bands of interest (usually the number of relevant conduction/valence bands)
      nphband = phband_f - phband_i + 1
      !
      ! Read the a2Fsave file
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
      allocate (eqkfit(nkfit),sfit(nkfit))
      allocate (iflag(nphband,nkfit),nkeff(nphband))
      allocate (dfk(nphband,nkfit,3),vk(nphband,nkfit,3))
      !
      ! Find the map between IBZ and the full-grid
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      !
      ! Deallocate unnecessary variables
      deallocate(wkfit,sfit,xkfit)
      !
      ! Weights for regular grid
      if ( lsoc .eqv. .true.) then
         wk = 1.0/nkfit
      else
         wk = 2.0/nkfit
      end if
      !
      ! If lscissors is true, then shift band energies (move conduction bands higher in energy)
      if (lscissors) then
        etfit(cbm_i:nbnd,:) = etfit(cbm_i:nbnd,:) + shift/RytoeV
        !cbm = cbm + shift
      end if
      !
      ! Call band velocities and forward derivatives
      call vband_ibz( nk1fit,nk2fit,nk3fit,nphband,nksfit,etfit(phband_i:phband_f,:),eqkfit,at, vk, dfk)
      !
      ! Include the 2pi/a factor
      vk = vk / tpi * alat
      !
      ! Determine the reduced grid for each band
      call reducegrid(nkfit,nksfit,nphband,etfit(phband_i:phband_f,:),eqkfit,efermi, &
     &                cut,T, nkeff,iflag)
      ! 
      ! Initiate Fermi integrals
      I0 = 0.0
      I1 = 0.0
      I2 = 0.0
      invT = 1.0/T 
      !
      !$ t0 = omp_get_wtime()  
      ! Loop over bands and kpoints
      !$omp parallel default(shared) &
      !$omp private(ibnd,ibnd_ph,ik,ikk,ind_k,i,j,fd,dfd,fac)
      do ibnd=phband_i,phband_f 
         !
         ibnd_ph = ibnd - phband_i + 1
         ! 
         !$omp do reduction(+: I0, I1, I2)
         do ik=1,nkeff(ibnd_ph)
            !
            ikk = iflag(ibnd_ph,ik)  ! ikk is in full-grid (just reduced)
            ind_k = eqkfit(ikk)      ! ind_k is in IBZ
            !
            ! Fermi factors
            fac = etfit(ibnd,ind_k) - efermi
            fd = 1.0/( exp(fac*invT) + 1.0 )
            dfd = invT * fd * (1.0 - fd)
            !
            ! Compute I0 (related to Thomas-Fermi screening)
            I0 = I0 + wk * dfd
            !
            ! Compute I1, I2
            do i=1,3
               do j=1,3
                  I1(i,j) = I1(i,j) + wk * dfd / invtau * vk(ibnd_ph,ikk,i) &
     &                               * vk(ibnd_ph,ikk,j)  
                  ! 
                  I2(i,j) = I2(i,j) + wk * dfd * fac / invtau * vk(ibnd_ph,ikk,i) &
     &                               * vk(ibnd_ph,ikk,j)
                  !
               end do ! j
            end do ! i
            !
         end do ! ik
         !$omp end do
         !
      end do ! ibnd
      !$omp end parallel
      !
      !Total integration time
      !$ t0 = omp_get_wtime() - t0 
      !
      ! Conductivity (convert to units of Ohm^-1 cm^-1) 
      sig = I1 * convsig / vol
      !
      !Compute Seebeck
      Se = 0.0
      call inv33(I1, inv_I1)
      do i=1,3
         do j=1,3
            do k=1,3
               Se(i,j) = Se(i,j) + I2(i,k)*inv_I1(k,j)
            end do
         end do
      end do 
      ! Units (convert to units of V/K)
      Se = Se * (-RytoeV * KtoRy / T)
      !
      ! Write D(Ef), Conductivity and Seebeck into file
      open(11,file='sig_1.out',status='unknown')
      open(12,file='Se_1.out',status='unknown')
      open(13,file='Def_1.out',status='unknown')
      write(11,"(10e14.6)") efermi * RytoeV, ((sig(i,j),i=1,3),j=1,3)
      write(12,"(10e14.6)") efermi * RytoeV, ((Se(i,j),i=1,3),j=1,3)
      write(13,"(2e14.6)") efermi * RytoeV, I0 / RytoeV
      close(11)
      close(12)
      close(13) 
      !
      !
      open(11,file='report_1',status='unknown')
      if ( t0 > 0 ) then
         write(11,"(A,I2)") "Number of threads = ", nthreads
         write(11,"(A,e14.6)") "Integration time(s) =", t0
      end if
      !
      write(11,"(A,I5)") "Number of kpoints (regular grid) = ", nkfit
      do ibnd=1,nphband
         write(11,"(A,2I5)") "band, reduced grid size", ibnd, nkeff(ibnd) 
      end do
      close(11)
      ! Free memory
      deallocate(etfit,eqkfit,dfk,vk,nkeff,iflag)

      contains
         !
         ! Determinant of a 3x3 matrix
         double precision function det33 (a)
         
         implicit none
         double precision, intent(in) :: a(3,3) 
        
         ! Determinant of a
         det33 = a(1,1)*( a(2,2)*a(3,3)-a(2,3)*a(3,2) ) -               &
     &           a(1,2)*( a(1,1)*a(3,3)-a(1,3)*a(3,1) ) +               &
     &           a(1,3)*( a(2,1)*a(3,2)-a(2,2)*a(3,1) ) 

         end function det33
         !
         !
         ! Inverse of a 3x3 matrix (analytical)
         subroutine inv33(a, inv_a)

         implicit none

         double precision, intent(in) :: a(3,3) 
         double precision, intent(out) :: inv_a(3,3)

         double precision :: deta_inv
 
         deta_inv = 1.0 / det33(a)

         inv_a(1,1) = deta_inv * ( a(2,2)*a(3,3)-a(2,3)*a(3,2) )
         inv_a(1,2) = deta_inv * ( a(1,3)*a(3,2)-a(1,2)*a(3,3) )       
         inv_a(1,3) = deta_inv * ( a(1,2)*a(2,3)-a(1,3)*a(2,2) )       

         inv_a(2,1) = deta_inv * ( a(2,3)*a(3,1)-a(2,1)*a(3,3) )
         inv_a(2,2) = deta_inv * ( a(1,1)*a(3,3)-a(1,3)*a(3,1) )       
         inv_a(2,3) = deta_inv * ( a(1,3)*a(2,1)-a(1,1)*a(2,3) )       

         inv_a(3,1) = deta_inv * ( a(2,1)*a(3,2)-a(2,2)*a(3,1) )
         inv_a(3,2) = deta_inv * ( a(1,2)*a(3,1)-a(1,1)*a(3,2) )       
         inv_a(3,3) = deta_inv * ( a(1,1)*a(2,2)-a(1,2)*a(2,1) )       

         end subroutine inv33


      end program fermi_int_1
     
