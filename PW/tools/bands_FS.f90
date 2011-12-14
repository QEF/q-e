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
!       Eyvaz Isaev, 2004-2009
!       eyvaz_isaev@yahoo.com, isaev@ifm.liu.se
!       
!       Theoretical Physics Department,
!       Moscow State Institute of Steel and Alloys, Russia  
!      
!       Department of Physics, Chemistry, and Biology (IFM), Linkoping University, Sweden,
!       
!       Division of  Materials Theory, Institute of Physics and Materials Sciene, 
!       Uppsala University, Sweden
!
! Description:
! The program reads output files for band structure calculations produced by PWscf. 
! Input_FS file contains reciprocal basis vectors, the Fermi level, grids numbers, and 
! System name extracted from self-consistent output file (See Input_FS).
! The output file(s) Bands_FS.bxsf (non spin-polarized) or Bands_FS_up.bxsf and Bands_FS_down.bxsf 
! (spin-polarized)  is (are)  written so that it can be used directly 
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
      real                  :: x(3),y(3),z(3)
      real,     allocatable :: valence(:)
      integer               :: n_kpoints, nbands
      integer               :: KS_number

      character*100 line 
      character*24 nkpt
      character*33 n_bands
      character*38 Band_structure
      character*13 kpoint
      character*80 sysname
      character*22 Magnetic
      character*9  blank
      character*16 KS_states

      logical      lsda

!
      nkpt='     number of k points='
      n_bands='nbnd'
      Band_structure='     End of band structure calculation'
      kpoint='          k ='
      blank='         '
      KS_states='Kohn-Sham states'
!
        Magnetic='     Starting magnetic'
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

        do while( .true. )
        read(5,'(a)',end=110) line
        if(line(16:31).eq.KS_states) then 
        goto 110
        endif
        enddo
110        continue

        Backspace(5)
        read(5,'(36x,I9)') KS_number
        print*, 'KS_number==', KS_number  
        
        if(n_last.gt.KS_number) then
        write(6,'("n_last > number of Kohn-Sham states")')
        write(6,'("Wrong input: you have specifed more bands than number of Kohn-Sham states")')
        stop
        endif
        
        print*, 'LSDA====', lsda

        rewind(5)

        do while( .true. )
        read(5,'(a)',end=111) line
        if(line(1:22).eq.Magnetic) then 
        lsda=.true. 
        goto 111
        endif
        enddo
111        continue
        
        print*, 'LSDA====', lsda

        rewind(5)

        do while( .true. )
         read(5,'(a)') line
         if(line(1:24).eq.nkpt) then 
            backspace(5)
            read(line,'(24x,i6)') n_kpoints
            goto 101
         endif
       enddo
 101   if(n_kpoints.gt.max_kpoints) then
         stop 'Toooooooo many k-points'
       endif
       
!     End of band structure calculation

        do while( .true. )
         read(5,'(a)',end=102) line
         if(line(1:38).eq.Band_Structure) then 
            goto 102
         endif
       enddo
 102   continue

        print*, '  lsda==', lsda
   
! Find bands number, nbands
!                

        read(5,*) 
        read(5,*) 
        read(5,*) 

        if(lsda.eqv..true.) then

        read(5,*) 
        read(5,*) 
        read(5,*) 

        endif        

        nlines=0
3        read(5,'(a)',end=4) line
        if(line(1:11).ne.blank) then
        nlines=nlines+1
        goto 3
        
        else
        
        goto 4
        endif
4        continue
        
        print*,'nlines==', nlines

        do k=1,nlines+1
        backspace(5)
        enddo
        
        nbands=0
        do k=1,nlines
        read(5,'(a)') line                
            do j=1,8
!            
! 9 is due to output format for e(n,k): 2X, 8f9.4            
!
            if(line((3+9*(j-1)):(3+9*j)).ne.blank) then
            nbands=nbands+1
            endif
            enddo
            
        enddo
        
        print*, 'nbands==', nbands

        if(lsda.eqv..true.) then   ! begin for lsda calculations

        n_kpoints=n_kpoints/2

        print*, 'kpoints=', n_kpoints

        allocate (e_up(n_kpoints,nbands)) 
        allocate (e_down(n_kpoints,nbands))

! back nlines+1 positions (number of eigenvalues lines plus one blank line)
!
        do k=1,nlines+1
        backspace(5)
        enddo
!
! back 3 positions for k-points
!
        backspace(5)
        backspace(5)
        backspace(5)

! Now ready to start
!
        read(5,*) 
!
! Reading spin-up energies
!        
        do k1=1,n_kpoints

        read(5,*) 
        read(5,*) 
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
        read(5,*,end=96) (e_down(k1,j),j=1,nbands)
        enddo
96      continue
         open(11,file='Bands_FS_up.bxsf',form='formatted')
         
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


         open(11,file='Bands_FS_down.bxsf',form='formatted')
         
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

        print*,'SPIN-POLARIZED CASE: FINISHED!!!!'         

!!! end for LSDA calculations

        else     ! end of lsda section
!        
        allocate (e_up(n_kpoints,nbands))
 
! back nlines+1 positions (number of eigenvalues lines plus one blank line)
!
        print*, 'nlines==', nlines

        do k=1,nlines+1
        backspace(5)
        enddo

!
! back 3 positions for k-points
!

        backspace(5)
        backspace(5)
        backspace(5)

        print*, 'n_kpoints===', n_kpoints        
               
      do k1=1,n_kpoints

        read(5,*) 
        read(5,*) 
        read(5,*) 

         read(5,*,end=98) (e_up(k1,j),j=1,nbands)
!         read(5,'(2x,8f9.4)',end=98) (e_up(k1,j),j=1,nbands)
      enddo
       
 98   continue       

         open(11,file='Bands_FS.bxsf',form='formatted')
         
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

        print*,'NON-SPIN-POLARIZED CASE: FINISHED!!!!'         
 
        endif
        
         stop 
      END PROGRAM bands_FS

