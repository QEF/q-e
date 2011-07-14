! For GWW
!
! Author: P. Umari
! Modified by G. Stenuit
!
subroutine distance_wannier
!
!this subroutine calculates the distance in a.u.
!between couples of wanniers and write on file
!
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : find_free_unit, prefix
  USE wannier_gw
  USE io_global,            ONLY : ionode, stdout
  USE fft_base,      ONLY : dfftp
  USE cell_base, ONLY: at, alat

  implicit none

  INTEGER :: ii,jj, i
  INTEGER :: iun, n(3)
  INTEGER :: rspacel(3)
  INTEGER, EXTERNAL :: ndistance
  REAL(kind=DP) :: d(3), dist

  if(ionode) then

     rspacel(1)=dfftp%nr1
     rspacel(2)=dfftp%nr2
     rspacel(3)=dfftp%nr3


     iun = find_free_unit()
     open( unit= iun, file=trim(prefix)//'.w_distance', status='unknown',form='unformatted')
     write(iun) nbnd_normal
     do ii=1,nbnd_normal
        do jj=ii,nbnd_normal
           do i=1,3
              n(i)=ndistance(w_centers(i,ii),w_centers(i,jj),rspacel(i))
              d(i)=dble(n(i))*at(i,i)*alat
           enddo
           dist=dsqrt(d(1))**2.d0+d(2)**2.d0+d(3)**2.d0
           write(iun) dist
        enddo
     enddo
     close(iun)
  endif

  return

end subroutine distance_wannier
