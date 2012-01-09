! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE wannier_u_matrix(U,hJ)

  USE io_global, ONLY: stdout, ionode, ionode_id
  USE io_files
  USE kinds, ONLY: DP
  USE wannier_new, ONLY: nwan, pp, wannier_occ, wannier_energy,wan_in

  IMPLICIT NONE
  INTEGER i,j, k,c, iwan, l
  real(DP) :: U, hJ, u2(10,10)

  INTEGER :: atoms(10)
  real(DP) :: rotm(10,10), unew(10,10), tmp

  WRITE(stdout,'(5x,a34)') 'Generation of interaction matrix U'
  WRITE(stdout,'(5x,a29)') '(works only for nspin=1 case)'
  WRITE(stdout,*)
  u2 = 0.d0
  CALL mk_u(2,5,U,hJ,u2)

  !rotation from TB-LMTO basis to our new

  rotm = 0.d0
  c = 0
  DO iwan=1, nwan
     DO j=1,wan_in(iwan,1)%ning
        IF(wan_in(iwan,1)%ing(j)%l==2) THEN
           c = c+1
           SELECT CASE(wan_in(iwan,1)%ing(j)%m)
           CASE(1)
              rotm(c,3) = wan_in(iwan,1)%ing(j)%c
           CASE(2)
              rotm(c,4) = wan_in(iwan,1)%ing(j)%c
           CASE(3)
              rotm(c,2) = wan_in(iwan,1)%ing(j)%c
           CASE(4)
              rotm(c,5) = wan_in(iwan,1)%ing(j)%c
           CASE(5)
              rotm(c,1) = wan_in(iwan,1)%ing(j)%c
           END SELECT
        ENDIF
     ENDDO
  ENDDO

  IF(c>5) CALL errore('Too many interactiong atoms - cant construct U matrix',c)

  DO i=1,5
     DO j=1,5
        rotm(i+5,j+5) = rotm(i,j)
     ENDDO
  ENDDO

  DO i = 1,10
     DO j = 1, 10
        tmp = 0.d0
        DO k=1,10
           DO l=1,10
              tmp=tmp+rotm(i,k)*u2(k,l)*rotm(j,l)
           ENDDO
        ENDDO
        unew(i,j)=tmp
     ENDDO
  ENDDO

!output
  DO i=1,c
     WRITE(stdout,'(5x,10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
  ENDDO
  DO i=6,5+c
     WRITE(stdout,'(5x,10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
  ENDDO
  WRITE(stdout,*)

  OPEN(70,file='umatrix',status='unknown',form='formatted')
  DO i=1,c
     WRITE(70,'(10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
  ENDDO
  DO i=6,5+c
     WRITE(70,'(10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
  ENDDO
  WRITE(70,*)
  CLOSE(70)

END SUBROUTINE wannier_u_matrix
