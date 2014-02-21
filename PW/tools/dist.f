!
! Copyright (C) 2003-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
      program dist
!----------------------------------------------------------------------
!
! find distances, nearest neighbors, angles 
! taking into account periodicity
!
! Input file:
!   first line contains  nat, ibrav                if ibrav = 0
!                        nat, ibrav, celldm(1-6)   if ibrav/= 0
!   if (ibrav.eq.0) :
!       card CELL_PARAMETERS
!       a(1,1) a(2,1) a(3,1)
!       a(1,2) a(2,2) a(3,2)
!       a(1,3) a(2,3) a(3,3)
!       [ alat units if celldm(1) was specified, a.u. if celldm(1)=0 ]
!   end if
!   card ATOMIC_POSITIONS
!   cards with atomic positions
! See documentation of pw.x input for more details
!
      implicit none
      integer nspx, ndistx, nax, nnx
      parameter (nspx=6, ndistx=1000, nax=200, nnx=4)
      integer nat, ibrav, nsp, ityp(nax)
      integer nsp1, nsp2, nat0, na, nb, n, nd, nn, nn1, nn2, nn3
      integer iout, i
      integer atom1(ndistx), atom2(ndistx), idx(ndistx), ndist
      character*3 atm(nspx), atm1, atm2
      character*80 filename, line
      character*1  capital, other_cell(ndistx)
      real*8 at(3,3), bg(3,3), celldm(6), omega, d(ndistx)
      real*8 tau(3,nax), dr(3), dd, dn1, dn2, dn3, dmin, dmax, scalef
      real*8 angolo(nnx*(nnx-1)/2), drv(3), drn(3,nnx), temp, rtemp(3)
      real*8 fact, pi, arg
      parameter (fact=0.529177d0, pi=3.141592653589793d0)
      logical crys, matches
      external capital, matches
!
!
      call get_file ( filename ) 
      open(unit=1,file=filename,form='formatted',status='old')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin data reading
      read(1,'(a)') line
      read(line,*) nat, ibrav
      if (nat.le.0.or.nat.gt.nax) call errore ('dist','wrong nat',1)
      if (ibrav.lt.0.or.ibrav.gt.14) call errore('dist','wrong ibrav',1)
      if (ibrav.eq.0) then
! read cell unit vectors from input
         read(1,'(a)') line
         do i=1,len(line)
            line(i:i) = capital(line(i:i))
         end do
         if (matches('CELL_PARAMETERS',line)) then
            READ (1,*) ( ( at(i,n), i=1,3 ), n=1,3 )
         else
            call errore ('dist','expecting card CELL_PARAMETERS',1)
         end if
!
         if (celldm (1).eq.0.d0) then
! set the lattice parameter from primitive vectors at
            if ( matches('ANGSTROM', line) ) then
! input at in angstrom
               celldm (1) = sqrt(at(1,1)**2+at(1,2)**2+at(1,3)**2)/fact
            else
! input at in atomic units (default)
               celldm (1) = sqrt(at(1,1)**2+at(1,2)**2+at(1,3)**2)
            end if
            do i=1,3
               at(i,1) = at(i,1)/celldm(1)
               at(i,2) = at(i,2)/celldm(1)
               at(i,3) = at(i,3)/celldm(1)
            end do
         end if
      else
         read(line,*) nat, ibrav, celldm(:)
! generate cell unit vectors 
         call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
         do i=1,3
            at(i,1) = at(i,1)/celldm(1)
            at(i,2) = at(i,2)/celldm(1)
            at(i,3) = at(i,3)/celldm(1)
         end do
      end if
      call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
!
      read(1,'(a)') line
      do i=1,len(line)
         line(i:i) = capital(line(i:i))
      end do
      if (matches('ATOMIC_POSITIONS',line)) then
         if ( matches('ALAT', line) ) then
            scalef = 1.d0
            crys = .false.
         else if ( matches('BOHR', line) ) then
            scalef = celldm(1)
            crys = .false.
         else if ( matches('CRYSTAL', line) ) then
            scalef = 1.d0
            crys = .true.
         else if ( matches('ANGSTROM', line) ) then
            scalef = celldm(1) * fact
            crys = .false.
         else
            scalef = 1.d0
            crys = .false.
         end if
      else
         read(line,*) scalef
      end if
      if (scalef.le.0.d0 .or. scalef.gt.1000.d0) 
     &     call errore ('dist','wrong scalef',1)

      nsp = 0
      do na=1,nat
         read(1,*) atm1, (tau(i,na),i=1,3)
         do n = 1, nsp
            if ( atm1 .eq. atm(n) ) then
               ityp(na) = n
               go to 10
            end if
         end do
         nsp = nsp + 1
         if (nsp .gt. nspx) call errore ('dist','too many types',1)
         atm(nsp) = atm1
         ityp(na) = nsp
 10      continue
      enddo
      if (crys) then
         call cryst_to_cart ( nat , tau, at, 1)
      else
         do na=1,nat
            do i=1,3
               tau(i,na)=tau(i,na)/scalef
            end do
         end do
      end if
      close(unit=1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End data reading
      print '(''output file (CR=output) > '',$)'
      read(5,'(a)',end=20,err=20) filename
      if (filename.eq.' ') then
         iout=6
      else
         iout=1
         open(unit=1,file=filename,form='formatted')
      end if
!
      scalef=fact*celldm(1)
 30   continue
      if (nsp.gt.1) then
         do n = 1, nsp
            print '(''species # '',i1,'' : '', a3)', n, atm(n)
         end do
         print '(''indices of species 1 and 2 > '',$)'
         read(5,*,end=20,err=20)  nsp1, nsp2
         if (nsp1 .le. 0 .or. nsp2 .le. 0 .or.            
     &      nsp1 .gt. nspx .or. nsp2 .gt. nspx) then
            print '(''  wrong indices (type control-D to exit)'')'
            go to 30
         end if
      else
         nsp1 = nsp
         nsp2 = nsp
      end if
      atm1=atm(nsp1)
      atm2=atm(nsp2)
      print '(''min and max distance  > '',$)'
      read(5,*,end=20,err=20) dmin, dmax
      if (dmin.ge.dmax.or.dmax.le.0.0) go to 20
!
!!!      dmin=max(dmin,0.01)
      ndist=0
      do na=1,nat
         if (atm(ityp(na)) .eq. atm1) then

            if (atm1.eq.atm2) then
               nat0=na+1
            else
               nat0=1
            end if
            do nb=nat0,nat
               if (atm(ityp(nb)) .eq. atm2) then

                  do i=1,3
                     dr(i) = (tau(1,na)-tau(1,nb))*bg(1,i) +   
     &                       (tau(2,na)-tau(2,nb))*bg(2,i) +  
     &                       (tau(3,na)-tau(3,nb))*bg(3,i)
                  enddo
                  do nn1=-2,2
                     dn1=dr(1)-nn1
                     do nn2=-2,2
                        dn2=dr(2)-nn2
                        do nn3=-2,2
                           dn3=dr(3)-nn3
                           dd = scalef* sqrt(                           
     &                       ( dn1*at(1,1)+dn2*at(1,2)+dn3*at(1,3) )**2+
     &                       ( dn1*at(2,1)+dn2*at(2,2)+dn3*at(2,3) )**2+
     &                       ( dn1*at(3,1)+dn2*at(3,2)+dn3*at(3,3) )**2)
                           if(dd.ge.dmin.and.dd.le.dmax) then
                              ndist=ndist+1
                              if (ndist.gt.ndistx) 
     &                           call errore ('dist','wrong ndist',1)
                              atom1(ndist)=na
                              atom2(ndist)=nb
                              d(ndist)= dd
                              if (nn1.eq.0.and.nn2.eq.0.and.nn3.eq.0)   
     &                             then
                                 other_cell(ndist)=' '
                              else
                                 other_cell(ndist)='*'
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end if
      end do
!
      idx(1)=0.0
      if (ndist.gt.0) call hpsort(ndist,d,idx)
!
      if (iout.eq.1)                                                    
     &     write(iout,100) atm1, atm2, dmin,dmax
      do nd=1,ndist
         write(iout,200) atom1(idx(nd)), atom2(idx(nd)),            
     &        other_cell(idx(nd)), d(nd)
      end do
!
      go to 30
!
 20   nn = nnx
      print '(/''number of neighbors (max '',i1,'') > '', $)', nnx
      read(5,*,end=21,err=21) n
      nn = n
      if (nn.gt.nnx) call errore ('dist','too many neighbors',1)
 21   continue
!
! look for nearest neighbors
!
      do na=1,nat
!!!         if (atm(ityp(na)) .eq. atm1) then
!
! ndist (.le.nnx) keep tracks of how many neighbors have been found
!
            ndist=0
            do nd=1,nn
               d(nd)=100000.0
               do i=1,3
                  drn(i,nd)=0.0
               end do
            end do
            do nb=1,nat
               do i=1,3
                  dr(i)=(tau(1,na)-tau(1,nb))*bg(1,i) +                 
     &                  (tau(2,na)-tau(2,nb))*bg(2,i) +                
     &                  (tau(3,na)-tau(3,nb))*bg(3,i)
               end do
               do nn1=-1,1
                  dn1=dr(1)-nn1
                  do nn2=-1,1
                     dn2=dr(2)-nn2
                     do nn3=-1,1
                        dn3=dr(3)-nn3
                        dd = scalef* sqrt(                              
     &                      ( dn1*at(1,1)+dn2*at(1,2)+dn3*at(1,3) )**2 +
     &                      ( dn1*at(2,1)+dn2*at(2,2)+dn3*at(2,3) )**2 +
     &                      ( dn1*at(3,1)+dn2*at(3,2)+dn3*at(3,3) )**2 )
                        do i=1,3
                           drv(i) = tau(i,na)-tau(i,nb) -               
     &                          (nn1*at(i,1)+nn2*at(i,2)+nn3*at(i,3))
                        end do
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
                                 do i=1,3
                                    rtemp(i)  = drn(i,nd)
                                 end do
                                 do i=1,3
                                    drn(i,nd) = drv(i)
                                 end do
                                 do i=1,3
                                    drv(i)    = rtemp(i)
                                 end do
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
                  arg = scalef**2 *
     &                 ( drn(1,nn1)*drn(1,nn2) +
     &                   drn(2,nn1)*drn(2,nn2) +
     &                   drn(3,nn1)*drn(3,nn2) ) / d(nn1) / d(nn2)
                  if(abs(arg)>1.d0) arg = sign(1.d0, arg)
                  angolo(nd) = 360/(2*pi) * acos ( arg )
               end do
            end do
            if (nd.ne.nn*(nn-1)/2) call errore('dist','internal err.',2)
!
! dd is the distance from the origin
!
            dd = sqrt(tau(1,na)**2 + tau(2,na)**2 + tau(3,na)**2)*scalef
            write(iout,250) atm(ityp(na)), na, (d(nn1),nn1=1,nn)
            write(iout,300) dd, (angolo(nn1),nn1=1,nn*(nn-1)/2)
!!!         end if
      end do
!
      stop
!
 100  format(' species: ',a3,' - ',a3,3x,f6.2,' < D <',f6.2)
 200  format(' atoms:', 2i6,a1, ' distance =',f10.5,' A')
 250  format(a3,i3,':  neighbors at ',4f8.3,' A')
 300  format(9x,'d(center):',f6.3,' A  angles  :',6f8.1)
!
!
      end
