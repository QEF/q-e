subroutine build_spectrum(ampl,en,ipol)
!this subroutine builds up the absorption spectrum
!and prints it on file

USE exciton
USE io_global,   ONLY : stdout,ionode
USE bse_wannier, ONLY : spectra_e_min,spectra_e_max,n_eig,spectra_nstep,spectra_broad
USE constants,   ONLY : RYTOEV, PI
USE cell_base,   ONLY : omega
USE io_files,    ONLY : tmp_dir,prefix

implicit none

REAL(DP), INTENT(in)  :: ampl(n_eig), en(n_eig)

INTEGER, EXTERNAL :: find_free_unit

REAL(DP), ALLOCATABLE :: omega_g(:),absorption(:),broad_abs(:),excdos(:)
COMPLEX(DP), ALLOCATABLE :: den(:,:),cspectrum(:),campl(:)
REAL(DP)    :: eta,step,prefac,sumdos,norm
COMPLEX(DP) :: lambda_sum
INTEGER     :: i,j,iun,ipol 

LOGICAL     :: debug

call start_clock('build_spectrum')
debug=.true.
eta=0.001

if(debug) then
   if(ionode) then
     do i=1,n_eig
        write(stdout,*) '#',i,'E=',en(i),'A=',ampl(i)
     enddo
   endif
endif
 
allocate(omega_g(spectra_nstep))
allocate(absorption(spectra_nstep))
allocate(excdos(spectra_nstep))
allocate(broad_abs(spectra_nstep))
allocate(cspectrum(spectra_nstep))
allocate(den(spectra_nstep,n_eig))
allocate(campl(n_eig))

if (ipol==1) then
!convert energy range in Ry
   spectra_e_min=spectra_e_min/RYTOEV
   spectra_e_max=spectra_e_max/RYTOEV
endif

!build the omega grid
step=(spectra_e_max-spectra_e_min)/dble(spectra_nstep-1)

do i=0, spectra_nstep-1
   den(i+1,1:n_eig)=dcmplx(1.d0,0.d0)/(dcmplx(en(1:n_eig),0.d0)-dcmplx(spectra_e_min+dble(i)*step,eta))
   omega_g(i+1)=(spectra_e_min+dble(i)*step)*RYTOEV
enddo

!compute the absorption spectrum
campl(1:n_eig)=dcmplx(ampl(1:n_eig),0.d0)
!campl(1:n_eig)=dcmplx(1.d0,0.d0)

cspectrum(1:spectra_nstep)=(0.d0,0.d0)

call zgemm('N','N',spectra_nstep,1,n_eig,(1.d0,0.d0),den,spectra_nstep,campl,n_eig,(0.d0,0.d0),cspectrum,spectra_nstep)

prefac=8.d0*PI/omega
!prefac=1.d0
absorption(1:spectra_nstep)=prefac*aimag(cspectrum(1:spectra_nstep))

!add gaussian broadening
broad_abs(1:spectra_nstep)=0.d0
do i=1,spectra_nstep 
     norm=0.d0
     do j=1,spectra_nstep
        broad_abs(i)=broad_abs(i)+&
                     absorption(j)*exp(-((omega_g(i)-omega_g(j))**2)/(2.d0*spectra_broad**2))
        norm=norm+exp(-((omega_g(i)-omega_g(j))**2)/(2.d0*spectra_broad**2))
     enddo
     broad_abs(i)=broad_abs(i)/norm
enddo

!compute DOS using a lorentzian
if (ipol==1) then
   excdos(1:spectra_nstep)=0.d0
   do i=0,spectra_nstep-1
       do j=1,n_eig
          excdos(i+1)=excdos(i+1)+2.d0*eta/(PI*((en(j)-spectra_e_min-dble(i)*step)**2+eta**2))
       enddo
   enddo
   excdos(1:spectra_nstep)=excdos(1:spectra_nstep)/(2.d0*n_eig)
endif



write(*,*) 'Absorption' 
write(*,*) 'Energy(eV)   Eps2 Eps2(Nogaussbroad)'

!print absorption spectrum on file
if(ionode) then
   iun = find_free_unit()

   if (ipol==1) then
      open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps2x.dat', status='unknown', form='formatted')
   elseif (ipol==2) then
      open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps2y.dat', status='unknown', form='formatted')
   elseif (ipol==3) then
      open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps2z.dat', status='unknown', form='formatted')
   endif

   do i=1,spectra_nstep
      write(iun,*) omega_g(i),broad_abs(i), absorption (i)
      write(*,*) omega_g(i),broad_abs(i), absorption(i)
   enddo
   close(iun)
endif

!print excdos spectrum on file
if(ionode.and.(ipol==1)) then
   sumdos=0.d0
   iun = find_free_unit()
   open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.excdos.dat', status='unknown', form='formatted')
   do i=1,spectra_nstep
      write(iun,*) omega_g(i),excdos(i)
      sumdos=sumdos+excdos(i)
   enddo
   close(iun)
   write(*,*) 'sumdos=',sumdos/dble(spectra_nstep)
endif

deallocate(omega_g,absorption,broad_abs,cspectrum,den,campl,excdos)

call stop_clock('build_spectrum')
end subroutine

