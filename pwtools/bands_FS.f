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
!       Eyvaz Isaev, 2004-2005 
!       eyvaz_isaev@yahoo.com, e.isaev@misis.ru, Eyvaz.Isaev@fysik.uu.se
!       Theoretical Physics Department,
!       Moscow State Institute of Steel and Alloys  
!
! Description:
! The program reads output files for self-consistent and band structure
! calculations produced by PWscf. The first file is used to extract 
! reciprocal basis vectors and the Fermi level, as well as bands number.
! The output file Bands.bxsf is  written so that it can be used directly 
! in conjunction with XCrySDen to visualize the Fermi Surface.
!
!----------------------------------------------------------------------- 
      PROGRAM bands_FS
!----------------------------------------------------------------------- 
!       
      implicit real*8(a-h,o-z)
  
      parameter (max_kpoints=100000, max_bands=500)

      real*8 e1(max_kpoints,max_bands) 
      real*8 kx(max_kpoints), ky(max_kpoints), kz(max_kpoints), 
     &     qk(max_kpoints)

      real*8 x(3),y(3),z(3)
      character*1 dummy1, dummy2
      character*100 line 
      character*13 dummy3
      character*24 nkpt
      character*4 n_bands
      character*32 Band_structure
      character*13 kpoint
      character*80 sysname
      character*5 Calc_type
!
      nkpt='     number of k points='         
      n_bands='nbnd'
      Band_structure='     Band Structure Calculation'
      kpoint='          k ='
!     
      open(9,file='Bands.out')               

      do while( .true. )
         read(9,'(a)') line
         if(line(1:24).eq.nkpt) then 
            print*, line(1:24)
            backspace(9)
            read(9,'(24x,i5)') n_kpoints
            print*,n_kpoints
            goto 100
         endif
      enddo
 100  continue

      do while( .true. )
         read(9,'(a)') line  
         if(line(22:25).eq.n_bands) then
            print*, line(22:25)
            backspace(9)
            read(9,'(29x,i6)') nbands
            print*,'nbands=',nbands
            goto 101
         endif 
      enddo
 101  continue

      if(n_kpoints.gt.max_kpoints) then
         stop 'Toooooooo many k-points'
      endif

      if(nbands.gt.max_bands) then
         stop 'Toooooooo many bands'
      endif
      
      do while( .true. )
         read(9,'(a)') line  
         if(line(1:31).eq.Band_structure) then
            print*, line(1:31)
            goto 102
         endif 
      enddo
 102  continue
       
      do while( .true. )
         read(9,'(a)') line  
         if(line(1:13).eq.kpoint) then
            print*, line(1:13)
            backspace(9)
            goto 103
         endif 
      enddo
 103  continue

      backspace(9)
       
      do k1=1,n_kpoints
         read(9,*)
         read(9,'(13X,3f7.4)')  kx(k1), ky(k1), kz(k1)
         read(9,*)
         read(9,*,end=99) (e1(k1,j),j=1,nbands)
         if(k1.eq.1) then 
            qk(k1)=0.
         else
            qk(k1)=qk(k1-1)+ 
     &           sqrt((kx(k1)-kx(k1-1))**2 +
     &           (ky(k1)-ky(k1-1))**2 + 
     &           (kz(k1)-kz(k1-1))**2)
        endif       
        
      enddo
       
 99   continue
       
      close(9)

      open(12,file='input_FS')
      read(12,*) n_start, n_last
      read(12,*) E_fermi
      read(12,*) sysname
      read(12,*) Calc_Type
      read(12,*) na,nb, nc
      read(12,*) x(1),x(2),x(3)
      read(12,*) y(1),y(2),y(3)
      read(12,*) z(1),z(2),z(3)
      
      print*,'E_Fermi=', E_Fermi
      x0=0.
      y0=0.
      z0=0.

      if(Calc_Type.eq.'Bands') then
         open(10,file='Bands_structure.out',form='formatted')
         
         do i=1, nbands
            do j=1, n_kpoints
               write(10, '(2f10.4)') qk(j), e1(j,i)
            enddo
            write(10,*)
         enddo
         close(10)
         
         stop
         
      else if(Calc_type.eq.'FS') then
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
            write(11, '(6f10.4)') (e1(j,i),j=1,n_kpoints)
         enddo
         
!     Write 2 last lines 
         write(11, '(" END_BANDGRID_3D")')
         write(11, '(" END_BLOCK_BANDGRID_3D")')
         
         close(11)
         
         stop 
      endif
       
      end               
       
