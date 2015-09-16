!
! Copyright (C) 2005  Eyvaz Isaev
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
! Used for k-points generation for the Fermi Surface construction
! Eyvaz Isaev, 2005
! eyvaz_isaev@yahoo.com, e.isaev@misis.ru
! Theoretical Physics Department 
! Moscow State Institute of Steel and Alloys 
!----------------------------------------------------------------------- 
      PROGRAM kvecs_FS
!----------------------------------------------------------------------- 

      implicit real*8(a-h,o-z)
      dimension x(3),y(3),z(3), rijk(100,100,100,3)
      character*80 sysname
!
      read(5,*) x(1),x(2),x(3)
      read(5,*) y(1),y(2),y(3)
      read(5,*) z(1),z(2),z(3)
      read(5,*) na,nb,nc
      read(5,*) sysname
!
      fna=dble(na)
      fnb=dble(nb)
      fnc=dble(nc)
      jj=0
      DO I=0,na
         I1=i+1
         DO J=0,nb
            j1=j+1
            DO K=0,nc
               K1=k+1
               Rijk(I1,j1,k1,1)=I*X(1)/fna + J*Y(1)/fnb + K*Z(1)/fnc
               Rijk(I1,j1,k1,2)=I*X(2)/fna + J*Y(2)/fnb + K*Z(2)/fnc
               Rijk(I1,j1,k1,3)=I*X(3)/fna + J*Y(3)/fnb + K*Z(3)/fnc
! 
               jj=jj+1
            END DO
         END DO
      END DO
!     
! 3    format('i1,j1,k1=',3i4,'  Rijk=',3f9.3)
!     
      print *,'jj=',jj
!     
      wk=1.0
      open(9,file='kvecs_'//sysname)
      write(9,'(i6)') jj
!     
      DO I=1,na+1
         DO J=1,nb+1
            DO K=1,nc+1
               write(9,'(3f12.6,f6.2)') rijk(i,j,k,1),rijk(i,j,k,2),rijk(i,j,k,3), wk
            END DO
         END DO
      END DO

      close (9)
      stop
      end
