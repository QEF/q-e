subroutine print_spectrum(aspectrum,im)
!this subroutine applies a gaussian broadening and prints the absorption
!spectrum on file

USE exciton
USE io_global,   ONLY : stdout,ionode
USE bse_wannier, ONLY : spectra_e_min,spectra_e_max,n_eig,spectra_nstep,spectra_broad,l_lanczos
USE constants,   ONLY : RYTOEV, PI
USE cell_base,   ONLY : omega
USE io_files,    ONLY : tmp_dir,prefix

implicit none

INTEGER, EXTERNAL :: find_free_unit

REAL(kind=DP), INTENT(inout)  :: aspectrum(spectra_nstep,3)
INTEGER           :: ipol

REAL(DP), ALLOCATABLE :: omega_g(:),broad_abs(:,:)
REAL(DP)    :: step,prefac,sumdos,norm
INTEGER     :: i,j,iun 

LOGICAL     :: debug, im

call start_clock('print_spectrum')
debug=.false.
 
allocate(omega_g(spectra_nstep))
allocate(broad_abs(spectra_nstep,3))

!build the omega grid (in eV)
step=(spectra_e_max-spectra_e_min)/dble(spectra_nstep-1)

do i=0, spectra_nstep-1
   omega_g(i+1)=(spectra_e_min+dble(i)*step)
enddo


prefac=4.d0*PI/omega
!prefac=1.d0
do ipol=1,3
   aspectrum(1:spectra_nstep,ipol)=prefac*aspectrum(1:spectra_nstep,ipol)

   broad_abs(1:spectra_nstep, ipol)=0.d0
   do i=1,spectra_nstep 
        norm=0.d0
        do j=1,spectra_nstep
           broad_abs(i,ipol)=broad_abs(i,ipol)+&
                    aspectrum(j,ipol)*exp(-((omega_g(i)-omega_g(j))**2)/(2.d0*spectra_broad**2))
           norm=norm+exp(-((omega_g(i)-omega_g(j))**2)/(2.d0*spectra_broad**2))
        enddo
        broad_abs(i, ipol)=broad_abs(i, ipol)/norm
   enddo
enddo

!print absorption aspectrum on file
if(im) then
   do ipol=1,3
      if(ionode) then
         iun = find_free_unit()

         if (ipol==1) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps2x.dat', status='unknown', form='formatted')
         elseif (ipol==2) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps2y.dat', status='unknown', form='formatted')
         elseif (ipol==3) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps2z.dat', status='unknown', form='formatted')
         endif
      
!         write(*,*) '# Energy(eV)   Eps2 Eps2(Nogaussbroad)'
         write(iun,*) '# Energy(eV)   Eps2 Eps2(Nogaussbroad)'
         do i=1,spectra_nstep
            write(iun,*) omega_g(i),broad_abs(i,ipol),aspectrum(i,ipol)
            !write(*,*) omega_g(i),broad_abs(i,ipol),aspectrum(i,ipol)
         enddo
         close(iun)
      endif
   enddo 
else
   do ipol=1,3
      if(ionode) then
         iun = find_free_unit()

         if (ipol==1) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps1x.dat', status='unknown', form='formatted')
         elseif (ipol==2) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps1y.dat', status='unknown', form='formatted')
         elseif (ipol==3) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'.eps1z.dat', status='unknown', form='formatted')
         endif
      
!         write(*,*) '# Energy(eV)   Eps1 Eps1(Nogaussbroad)'
         write(iun,*) '# Energy(eV)   Eps1 Eps1(Nogaussbroad)'
         do i=1,spectra_nstep
            write(iun,*) omega_g(i),broad_abs(i,ipol),aspectrum(i,ipol)
            !write(*,*) omega_g(i),broad_abs(i,ipol),aspectrum(i,ipol)
         enddo
         close(iun)
      endif
   enddo 
endif

deallocate(omega_g,broad_abs)
call stop_clock('print_spectrum')

end subroutine

