!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
PROGRAM dist
!----------------------------------------------------------------------
  !
  ! find distances, nearest neighbors, angles 
  ! taking into account periodicity
  !
  ! Same input file as in pw.x or cp.x
  !
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE mp_world,          ONLY : nproc
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  INTEGER :: exit_status = 0
  !
  CALL mp_startup ( )
  CALL environment_start ( 'DIST' )
  IF ( nproc > 1 ) CALL errore ('dist','run on a single process!',1)
  !
  CALL read_input_file ('PW', input_file_ )
  ! ... convert to internal variables
  CALL iosys()
  !
  CALL run_dist  ( )
  !
  ! CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
  STOP
  !
END PROGRAM dist
!----------------------------------------------------------------------
SUBROUTINE run_dist ( ) 
!----------------------------------------------------------------------
  !
  USE kinds,     ONLY : dp
  USE constants, ONLY : pi, bohr_radius_angs
  USE cell_base, ONLY : at, bg, alat
  USE ions_base, ONLY : atm, nat, ityp, tau, nsp
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  integer, parameter:: ounit=4, ndistx=9999, nn=4
  real(dp), parameter:: dmin=0.01_dp, dmax=3.0_dp
  integer :: nsp1, nsp2, na, nb, n, nd, nn1, nn2, nn3, i, nbad
  integer :: atom1(ndistx), atom2(ndistx), idx(ndistx), ndist
  character(len=3 ) :: atm1, atm2
  character(len=80) :: filename, line
  character(len=1)  :: other_cell(ndistx)
  real(dp) :: d(ndistx)
  real(dp) :: dr(3), dd, dn1, dn2, dn3, scalef, arg
  real(dp) :: angolo(nn*(nn-1)/2), drv(3), drn(3,nn), temp, rtemp(3)
  !
  open(unit=ounit,file='dist.out',status='unknown',form='formatted')
  write(stdout,'(/,5x,"Output written to file dist.out")')
  !
  scalef=bohr_radius_angs*alat
  !
  ndist=0
  nbad =0
  do na=1,nat
     do nb=na+1,nat
        dr(:) = (tau(1,na)-tau(1,nb))*bg(1,:) + &
                (tau(2,na)-tau(2,nb))*bg(2,:) + &
                (tau(3,na)-tau(3,nb))*bg(3,:)
        do nn1=-1,1
           dn1=dr(1)-nn1
           do nn2=-1,1
              dn2=dr(2)-nn2
              do nn3=-1,1
                 dn3=dr(3)-nn3
                 dd = scalef * sqrt( &
                         ( dn1*at(1,1)+dn2*at(1,2)+dn3*at(1,3) )**2 + &
                         ( dn1*at(2,1)+dn2*at(2,2)+dn3*at(2,3) )**2 + &
                         ( dn1*at(3,1)+dn2*at(3,2)+dn3*at(3,3) )**2 )
                 if (dd < dmin) then
                    nbad=nbad+1
                    if (nn1==0 .and. nn2==0 .and. nn3==0) then
                       write(ounit,60) na,nb
                    else
                       write(ounit,61) na,nb
                    end if
                 else if (dd < dmax) then
                    ndist=ndist+1
                    if (ndist > ndistx) then
                       write(stdout,62) ndistx, dmax
                       go to 20 
                    end if
                    atom1(ndist)=na
                    atom2(ndist)=nb
                    d(ndist)= dd
                    if (nn1==0 .and. nn2==0 .and. nn3==0) then
                       other_cell(ndist)=' '
                    else
                       other_cell(ndist)='*'
                    end if
                 end if
              end do
           end do
        end do
     end do
  end do
20 continue
  !
  if (nbad > 0) then
     write(stdout,59) nbad
     close(unit=ounit,status='keep')
     return
  end if
  !
  idx(1)=0.0
  if (ndist.gt.0) call hpsort(ndist,d,idx)
  !
  write(ounit,50) dmax
  !
  do nd=1,ndist
     na=atom1(idx(nd))
     nb=atom2(idx(nd))
     atm1=trim(atm(ityp(na)))
     atm2=trim(atm(ityp(nb)))
     write(ounit,200) na,nb,adjustr(atm1),atm2,d(nd), other_cell(idx(nd))
  end do
  !
  write(ounit,70) nn
  !
  ! look for nearest neighbors
  !
  do na=1,nat
     !
     ! ndist keeps tracks of how many neighbors have been found
     !
     ndist=0
     do nd=1,nn
        d(nd)=100000.0
        drn(:,nd)=0.0
     end do
     do nb=1,nat
        dr(:)=(tau(1,na)-tau(1,nb))*bg(1,:) + &               
              (tau(2,na)-tau(2,nb))*bg(2,:) + &              
              (tau(3,na)-tau(3,nb))*bg(3,:)
        do nn1=-1,1
           dn1=dr(1)-nn1
           do nn2=-1,1
              dn2=dr(2)-nn2
              do nn3=-1,1
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
                 if (dd.gt.0.01) then
! straight insertion: look for first nn neighbors
                    do nd=1,nn
                       if (dd.lt.d(nd)) then
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
                       end if
                    end do
                 end if
              end do
           end do
        end do
     end do
!
     if (ndist.ne.nn) call errore ('dist','internal error',1)
!
! calculate angles with nearest neighbors
!
      nd=0
      do nn1=1,nn
         do nn2=nn1+1,nn
            nd=nd+1
            arg = scalef**2 * ( drn(1,nn1)*drn(1,nn2) + &
                                drn(2,nn1)*drn(2,nn2) + &
                                drn(3,nn1)*drn(3,nn2) ) / d(nn1) / d(nn2)
            if(abs(arg)>1.d0) arg = sign(1.d0, arg)
            angolo(nd) = 360/(2*pi) * acos ( arg )
         end do
      end do
      if (nd.ne.nn*(nn-1)/2) call errore('dist','internal err.',2)
!
! dd is the distance from the origin
!
      dd = sqrt(tau(1,na)**2 + tau(2,na)**2 + tau(3,na)**2)*scalef
      write(ounit,250) na, atm(ityp(na)), (d(nn1),nn1=1,nn)
      write(ounit,300) dd, (angolo(nn1),nn1=1,nn*(nn-1)/2)
   end do
   !
   close(unit=ounit,status='keep')
   return
   !
  50  format('Distances between atoms, up to dmax=',f6.2,' A   (* = with lattice translation)',/,'  #1  #2   bond       d')
  59  format(/,80('*'),/,' Fatal error: ',i3,' overlapping atoms (see file dist.out)',/,80('*'))
  60  format(' Atom',i4,' and',i4,' overlap')
  61  format(' Atom',i4,' and',i4,' overlap (with lattice translation)')
  62  format(/,80('*'),/,' Serious warning: more than ',i4,' distances smaller than',f6.2,' A found',/,80('*'))
  70  format(/,'Nearest neighbors for each atom (up to ',i1,')',/)
 200  format(2i4,a4,'-',a3,f10.5,' A ',a1)
 250  format(' atom ',i3,' species ',a3,': neighbors at ',4f8.3,' A')
 300  format(9x,'d(center):',f6.3,' A  angles  :',6f8.1)
!
end subroutine run_dist
