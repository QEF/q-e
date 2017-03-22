!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE run_dist ( exit_status )
!----------------------------------------------------------------------
  !
  ! Find distances, nearest neighbors, angles, taking into account periodicity
  ! Requires as input: lattice vectors, types and positions of atoms
  ! Must be run on a signle process only. Output in file "dist.out"
  !
  USE kinds,     ONLY : dp
  USE constants, ONLY : pi, bohr_radius_angs
  USE cell_base, ONLY : at, bg, alat
  USE ions_base, ONLY : atm, nat, ityp, tau, nsp
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(out) :: exit_status
  INTEGER, PARAMETER:: ounit=4, ndistx=9999, nn=4
  real(dp), PARAMETER:: dmin=0.01_dp, dmax=3.0_dp
  INTEGER :: ndist, nsp1, nsp2, na, nb, n, nd, nn1, nn2, nn3, i, nbad
  INTEGER,  ALLOCATABLE :: atom1(:), atom2(:), idx(:)
  CHARACTER(len=3 ) :: atm1, atm2
  CHARACTER(len=80) :: filename, line
  CHARACTER(len=1)  :: other_cell(ndistx)
  real(dp), ALLOCATABLE :: d(:)
  real(dp) :: dr(3), dd, dn1, dn2, dn3, scalef, arg
  real(dp) :: angolo(nn*(nn-1)/2), drv(3), drn(3,nn), temp, rtemp(3)
  !
  exit_status=0
  OPEN(unit=ounit,file='dist.out',status='unknown',form='formatted')
  WRITE(stdout,'(/,5x,"Output written to file dist.out")')
  !
  scalef=bohr_radius_angs*alat
  !
  ALLOCATE (d(ndistx), atom1(ndistx), atom2(ndistx), idx(ndistx))
  ndist=0
  nbad =0
  DO na=1,nat
     DO nb=na+1,nat
        dr(:) = (tau(1,na)-tau(1,nb))*bg(1,:) + &
                (tau(2,na)-tau(2,nb))*bg(2,:) + &
                (tau(3,na)-tau(3,nb))*bg(3,:)
        DO nn1=-1,1
           dn1=dr(1)-nn1
           DO nn2=-1,1
              dn2=dr(2)-nn2
              DO nn3=-1,1
                 dn3=dr(3)-nn3
                 dd = scalef * sqrt( &
                         ( dn1*at(1,1)+dn2*at(1,2)+dn3*at(1,3) )**2 + &
                         ( dn1*at(2,1)+dn2*at(2,2)+dn3*at(2,3) )**2 + &
                         ( dn1*at(3,1)+dn2*at(3,2)+dn3*at(3,3) )**2 )
                 IF (dd < dmin) THEN
                    nbad=nbad+1
                    IF (nn1==0 .and. nn2==0 .and. nn3==0) THEN
                       WRITE(ounit,60) na,nb
                    ELSE
                       WRITE(ounit,61) na,nb
                    ENDIF
                 ELSEIF (dd < dmax) THEN
                    ndist=ndist+1
                    IF (ndist > ndistx) THEN
                       WRITE(stdout,62) ndistx, dmax
                       GOTO 20
                    ENDIF
                    atom1(ndist)=na
                    atom2(ndist)=nb
                    d(ndist)= dd
                    IF (nn1==0 .and. nn2==0 .and. nn3==0) THEN
                       other_cell(ndist)=' '
                    ELSE
                       other_cell(ndist)='*'
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
20 CONTINUE
  !
  IF (nbad > 0) THEN
     WRITE(stdout,59) nbad
     exit_status=2
     CLOSE(unit=ounit,status='keep')
     RETURN
  ENDIF
  !
  idx(1)=0.0
  IF (ndist>0) CALL hpsort(ndist,d,idx)
  !
  WRITE(ounit,50) dmax
  !
  DO nd=1,ndist
     na=atom1(idx(nd))
     nb=atom2(idx(nd))
     atm1=trim(atm(ityp(na)))
     atm2=trim(atm(ityp(nb)))
     WRITE(ounit,200) na,nb,adjustr(atm1),atm2,d(nd), other_cell(idx(nd))
  ENDDO
  !
  WRITE(ounit,70) nn
  !
  ! look for nearest neighbors
  !
  DO na=1,nat
     !
     ! ndist keeps tracks of how many neighbors have been found
     !
     ndist=0
     DO nd=1,nn
        d(nd)=100000.0
        drn(:,nd)=0.0
     ENDDO
     DO nb=1,nat
        dr(:)=(tau(1,na)-tau(1,nb))*bg(1,:) + &
              (tau(2,na)-tau(2,nb))*bg(2,:) + &
              (tau(3,na)-tau(3,nb))*bg(3,:)
        DO nn1=-1,1
           dn1=dr(1)-nn1
           DO nn2=-1,1
              dn2=dr(2)-nn2
              DO nn3=-1,1
                 dn3=dr(3)-nn3
                 dd = scalef* sqrt( &
                        ( dn1*at(1,1)+dn2*at(1,2)+dn3*at(1,3) )**2 + &
                        ( dn1*at(2,1)+dn2*at(2,2)+dn3*at(2,3) )**2 + &
                        ( dn1*at(3,1)+dn2*at(3,2)+dn3*at(3,3) )**2 )
                 drv(:) = tau(:,na)-tau(:,nb) - &
                          (nn1*at(:,1)+nn2*at(:,2)+nn3*at(:,3))
!
! the "first" neighbor is the atom itself
!
                 IF (dd>0.01) THEN
! straight insertion: look for first nn neighbors
                    DO nd=1,nn
                       IF (dd<d(nd)) THEN
! swap d(nd) with dd
                          temp = d(nd)
                          d(nd)= dd
                          dd   = temp
! do the same for delta r
                          rtemp(:)  = drn(:,nd)
                          drn(:,nd) = drv(:)
                          drv(:)    = rtemp(:)
!
                          ndist=min(ndist+1,nn)
                       ENDIF
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!
     IF (ndist/=nn) THEN
        CALL infomsg ('dist','internal error 1')
        exit_status=1
     ENDIF
!
! calculate angles with nearest neighbors
!
      nd=0
      DO nn1=1,nn
         DO nn2=nn1+1,nn
            nd=nd+1
            arg = scalef**2 * ( drn(1,nn1)*drn(1,nn2) + &
                                drn(2,nn1)*drn(2,nn2) + &
                                drn(3,nn1)*drn(3,nn2) ) / d(nn1) / d(nn2)
            IF(abs(arg)>1.d0) arg = sign(1.d0, arg)
            angolo(nd) = 360/(2*pi) * acos ( arg )
         ENDDO
      ENDDO
      IF (nd/=nn*(nn-1)/2) THEN
        CALL infomsg ('dist','internal error 2')
        exit_status=1
     ENDIF
!
! dd is the distance from the origin
!
      dd = sqrt(tau(1,na)**2 + tau(2,na)**2 + tau(3,na)**2)*scalef
      WRITE(ounit,250) na, atm(ityp(na)), (d(nn1),nn1=1,nn)
      WRITE(ounit,300) dd, (angolo(nn1),nn1=1,nn*(nn-1)/2)
   ENDDO
   DEALLOCATE (d, atom1, atom2, idx )
   !
   CLOSE(unit=ounit,status='keep')
   RETURN
   !
  50  FORMAT('Distances between atoms, up to dmax=',f6.2,' A   (* = with lattice translation)',/,'  #1  #2   bond       d')
  59  FORMAT(/,80('*'),/,' Fatal error: ',i3,' overlapping atoms (see file dist.out)',/,80('*'))
  60  FORMAT(' Atom',i4,' and',i4,' overlap')
  61  FORMAT(' Atom',i4,' and',i4,' overlap (with lattice translation)')
  62  FORMAT(/,80('*'),/,' Serious warning: more than ',i4,' distances smaller than',f6.2,' A found',/,80('*'))
  70  FORMAT(/,'Nearest neighbors for each atom (up to ',i1,')',/)
 200  FORMAT(2i4,a4,'-',a3,f10.5,' A ',a1)
 250  FORMAT(' atom ',i3,' species ',a3,': neighbors at ',4f8.3,' A')
 300  FORMAT(9x,'d(center):',f6.3,' A  angles  :',6f8.1)
!
END SUBROUTINE run_dist
