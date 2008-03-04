!
! Copyright (C) 2005  Eyvaz Isaev
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!       Program is designed to map the Fermi Surface using XCrySDen
!       See www.xcrysden.org 
!
!       Eyvaz Isaev, 2004-2008 
!       eyvaz_isaev@yahoo.com, isaev@ifm.liu.se, Eyvaz.Isaev@fysik.uu.se
!       
!       Theoretical Physics Department,
!       Moscow State Institute of Steel and Alloys  
!      
!       Department of Physics, Chemistry, and Biology (IFM), Linkoping University, Sweden,
!       
!       Condensed Matter Theory Group, Uppsala University, Sweden
!
! Description:
! The program reads output files for self-consistent and band structure
! calculations produced by PWscf. The first file is used to extract 
! reciprocal basis vectors and the Fermi level, as well as bands number.
! The output file Bands.bxsf is  written so that it can be used directly 
! in conjunction with XCrySDen to visualize the Fermi Surface.
!
! Spin-polarized calculations are allowed
!
!----------------------------------------------------------------------- 
      PROGRAM bands_FS
!----------------------------------------------------------------------- 
!       
      implicit real*8(a-h,o-z)
  
      parameter (max_kpoints=100000, max_bands=500)

      real,     allocatable :: e_up(:,:),e_down(:,:)
      real                  :: kx,ky,kz
      real                  :: x(3),y(3),z(3)
      integer               :: kpoints, nbands

      character*1 dummy1, dummy2
      character*100 line 
      character*13 dummy3
      character*24 nkpt
      character*33 n_bands
      character*32 Band_structure
      character*13 kpoint
      character*80 sysname
      character*22 Magnetic
      character*7  SPINUP
      logical      lsda

!
      nkpt='     number of k points='
      n_bands='nbnd'
      Band_structure='     Band Structure Calculation'
      kpoint='          k ='
!
        Magnetic='     Starting magnetic'
	SPINUP='SPIN UP'
	lsda=.false.

!
! Read input information
!
        open(12,file='input_FS')

        read(12,*) n_start, n_last
        read(12,*) E_fermi
        read(12,*) sysname
        read(12,*) na,nb, nc
        read(12,*) x(1),x(2),x(3)
        read(12,*) y(1),y(2),y(3)
        read(12,*) z(1),z(2),z(3)
      
        print*,'E_Fermi=', E_Fermi
        x0=0.
        y0=0.
        z0=0.

	close(12)
!	

        do while( .true. )
	read(5,'(a)',end=110) line
	if(line(1:22).eq.Magnetic) then 
        lsda=.true. 
	print*, line(1:22)
	goto 110
	endif
	enddo
110	continue
	
	print*, 'LSDA====', lsda

	rewind(5)

        do while( .true. )
         read(5,'(a)') line
         if(line(1:24).eq.nkpt) then 
            read(line,'(24x,i5)') n_kpoints
            print *,n_kpoints
            goto 101
         endif
       enddo
 101   if(n_kpoints.gt.max_kpoints) then
         stop 'Toooooooo many k-points'
       endif

        do while (.true.)
        read(5,'(a)') line
        if(line(22:25).eq.n_bands) then
        print*, line(22:25)
        backspace(5)
        read(5,'(29x,i6)') nbands
        print*,'nbands=',nbands
        goto 121
        endif
        enddo
121     if(nbands.gt.max_bands) then
	stop 'Tooooooooo many bands'
	endif

	print*, 'n_bands==', nbands
	print*, ' lsda==', lsda
   
	if(lsda.eqv..true.) then   ! begin for lsda calculations

	n_kpoints=n_kpoints/2

	print*, 'in lsda calculations==', n_kpoints

	allocate (e_up(n_kpoints,nbands)) 
	allocate (e_down(n_kpoints,nbands))
	
	do while (.true.) 
	read(5,'(a)') line  
	if(line(9:15).eq.SPINUP) then
        print*, line(9:15)
	goto 105
	endif 
	enddo
105	continue

	read(5,*) 

! Reading spin-up energies
!
	do k1=1,n_kpoints
        read(5,*) 
	read(5,'(13x,3f7.4)')  kx, ky, kz
	read(5,*) 
	read(5,*,end=99) (e_up(k1,j),j=1,nbands)
	enddo
99      continue


	read(5,*)
	read(5,*)
	read(5,*)

! Reading Spin-down bands

	do k1=1,n_kpoints
        read(5,*)
	read(5,*)
	read(5,*)
	read(5,*,end=99) (e_down(k1,j),j=1,nbands)
	enddo

         open(11,file='Bands_up.bxsf',form='formatted')
         
!     Write  header file here
         
         write(11, '(" BEGIN_INFO")')
         write(11, '("   #")') 
         write(11, '("   # this is a Band-XCRYSDEN-Structure-File")')
         write(11, '("   # aimed at Visualization of Fermi Surface")')
         write(11, '("   #")') 
         write(11, '("   # Case:   ",A)')     Sysname
         write(11, '("   #")') 
         write(11, '(" Fermi Energy:    ", f12.4)') E_Fermi 
         write(11, '(" END_INFO")') 
         
         write(11, '(" BEGIN_BLOCK_BANDGRID_3D")')
         write(11, '(" band_energies")')
         write(11, '(" BANDGRID_3D_BANDS")')
         write(11, '(I5)')  n_last-n_start+1
         write(11, '(3I5)') na+1, nb+1, nc+1
         write(11, '(3f10.6)') x0,  y0,  z0
         write(11, '(3f10.6)') x(1), x(2), x(3)
         write(11, '(3f10.6)') y(1), y(2), y(3)
         write(11, '(3f10.6)') z(1), z(2), z(3)
         
         do i=n_start, n_last
            write(11, '("BAND:", i4)') i
            write(11, '(6f10.4)') (e_up(j,i),j=1,n_kpoints)
         enddo
         
!     Write 2 last lines 
         write(11, '(" END_BANDGRID_3D")')
         write(11, '(" END_BLOCK_BANDGRID_3D")')
         
         close(11)


         open(11,file='Bands_down.bxsf',form='formatted')
         
!     Write  header file here
         
         write(11, '(" BEGIN_INFO")')
         write(11, '("   #")') 
         write(11, '("   # this is a Band-XCRYSDEN-Structure-File")')
         write(11, '("   # aimed at Visualization of Fermi Surface")')
         write(11, '("   #")') 
         write(11, '("   # Case:   ",A)')     Sysname
         write(11, '("   #")') 
         write(11, '(" Fermi Energy:    ", f12.4)') E_Fermi 
         write(11, '(" END_INFO")') 
         
         write(11, '(" BEGIN_BLOCK_BANDGRID_3D")')
         write(11, '(" band_energies")')
         write(11, '(" BANDGRID_3D_BANDS")')
         write(11, '(I5)')  n_last-n_start+1
         write(11, '(3I5)') na+1, nb+1, nc+1
         write(11, '(3f10.6)') x0,  y0,  z0
         write(11, '(3f10.6)') x(1), x(2), x(3)
         write(11, '(3f10.6)') y(1), y(2), y(3)
         write(11, '(3f10.6)') z(1), z(2), z(3)
         
         do i=n_start, n_last
            write(11, '("BAND:", i4)') i
            write(11, '(6f10.4)') (e_down(j,i),j=1,n_kpoints)
         enddo
         
!     Write 2 last lines 
         write(11, '(" END_BANDGRID_3D")')
         write(11, '(" END_BLOCK_BANDGRID_3D")')

         
         close(11)

	deallocate (e_up)
	deallocate (e_down)

	print*,'LSDA FINISHED!!!!'	 

!!! end for LSDA calculations

	else     ! end of lsda section
!	
	allocate (e_up(n_kpoints,nbands))
 
      do while( .true. )
         read(5,'(a)') line  
         if(line(1:13).eq.kpoint) then
            print *, line(1:13)
            backspace(5)
            goto 103
         endif 
      enddo
 103  continue

      backspace(5)
       
      do k1=1,n_kpoints
         read(5,*)
         read(5,'(13X,3f7.4)')  kx, ky, kz
         read(5,*)
         read(5,*,end=98) (e_up(k1,j),j=1,nbands)
      enddo
       
 98   continue       

         open(11,file='Bands.bxsf',form='formatted')
         
!     Write  header file here
         
         write(11, '(" BEGIN_INFO")')
         write(11, '("   #")') 
         write(11, '("   # this is a Band-XCRYSDEN-Structure-File")')
         write(11, '("   # aimed at Visualization of Fermi Surface")')
         write(11, '("   #")') 
         write(11, '("   # Case:   ",A)')     Sysname
         write(11, '("   #")') 
         write(11, '(" Fermi Energy:    ", f12.4)') E_Fermi 
         write(11, '(" END_INFO")') 
         
         write(11, '(" BEGIN_BLOCK_BANDGRID_3D")')
         write(11, '(" band_energies")')
         write(11, '(" BANDGRID_3D_BANDS")')
         write(11, '(I5)')  n_last-n_start+1
         write(11, '(3I5)') na+1, nb+1, nc+1
         write(11, '(3f10.6)') x0,  y0,  z0
         write(11, '(3f10.6)') x(1), x(2), x(3)
         write(11, '(3f10.6)') y(1), y(2), y(3)
         write(11, '(3f10.6)') z(1), z(2), z(3)
         
         do i=n_start, n_last
            write(11, '("BAND:", i4)') i
            write(11, '(6f10.4)') (e_up(j,i),j=1,n_kpoints)
         enddo
         
!     Write 2 last lines 
         write(11, '(" END_BANDGRID_3D")')
         write(11, '(" END_BLOCK_BANDGRID_3D")')
         
         close(11)
         
	deallocate (e_up)
 
	endif
	
         stop 
         end               
       

