!----------------------------------------------------------------------
      program dist
!----------------------------------------------------------------------
!
! find distances, nearest neighbors, angles 
! taking into account periodicity
!
! Input file:
!   first line contains  nat, ibrav, celldm(1-6)
!   followed by card ATOMIC_POSITIONS
!   followed by cards with atomic positions
! See documentation of pw.x input for more details
!
      implicit none
      integer nspx, ndistx, nax, nnx
      parameter (nspx=6, ndistx=1000, nax=200, nnx=4)
      integer nat, ibrav, nsp, ityp(nax)
      integer nsp1, nsp2, nat0, na, nb, n, nd, nn, nn1, nn2, nn3
      integer iout, i, iargc
      integer atom1(ndistx), atom2(ndistx), index(ndistx), ndist
      character*3 atm(nspx), atm1, atm2
      character*80 filename, line
      character*1  capital, other_cell(ndistx)
      real*8 at(3,3), bg(3,3), celldm(6), omega, d(ndistx)
      real rr(3,nax), dr(3), dd, dn1, dn2, dn3, dmin, dmax, scale
      real angolo(nnx*(nnx-1)/2), drv(3), drn(3,nnx), temp, rtemp(3)
      real fact /0.529177/, pi/3.1415926/
      logical matches
      external capital, matches, iargc
!
!
      i=iargc()
      if (i.eq.0) then
         print '(''  input file > '',$)'
         read(5,'(a)',end=20,err=20)  filename
      else if(i.eq.1) then
         call getarg(1,filename)
      else
         print '(''   usage: dist  [input file] '')'
      end if
      open(unit=1,file=filename,form='formatted',status='old',iostat=i)
      if (i.ne.0) then
         print '(''   file not found '')'
         stop
      end if
!
      read(1,*) nat, ibrav, celldm
      if (nat.le.0.or.nat.gt.nax) stop ' nat!'
      if (ibrav.le.0.or.ibrav.gt.14) stop ' ibrav!'
      call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
      do i=1,3
         at(i,1) = at(i,1)/celldm(1)
         at(i,2) = at(i,2)/celldm(1)
         at(i,3) = at(i,3)/celldm(1)
      end do
      call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
!
      read(1,'(a)') line
      do i=1,len(line)
         line(i:i) = capital(line(i:i))
      end do
      if (matches('ATOMIC_POSITIONS',line)) then
         if ( matches('ALAT', line) ) then
            scale = 1.0
         else if ( matches('BOHR', line) ) then
            scale = celldm(1)
         else if ( matches('CRYSTAL', line) ) then
            stop 'option crystal not implemented'
         else if ( matches('ANGSTROM', line) ) then
            scale = celldm(1) * fact
         else
            scale = 1.0
         end if
      else
         read(line,*) scale
      end if

      if (scale.le.0.0) stop ' scale ! '
      fact=fact*celldm(1)
      nsp = 0
      do na=1,nat
         read(1,*) atm1, (rr(i,na),i=1,3)
         do n = 1, nsp
            if ( atm1 .eq. atm(n) ) then
               ityp(na) = n
               go to 10
            end if
         end do
         nsp = nsp + 1
         if (nsp .gt. nspx)  stop ' too many atom types!'
         atm(nsp) = atm1
         ityp(na) = nsp
 10      continue
         do i=1,3
            rr(i,na)=rr(i,na)/scale
         end do
      enddo
      close(unit=1)
!
      print '(''  output file (CR=output) > '',$)'
      read(5,'(a)',end=20,err=20) filename
      if (filename.eq.' ') then
         iout=6
      else
         iout=1
         open(unit=1,file=filename,form='formatted')
      end if
!
 30   continue
      if (nsp.gt.1) then
         do n = 1, nsp
            print '(''  species # '',i1,'' : '', a3)', n, atm(n)
         end do
         print '(''  indeces of species 1 and 2 > '',$)'
         read(5,*,end=20,err=20)  nsp1, nsp2
         if (nsp1 .le. 0 .or. nsp2 .le. 0 .or.                          &
     &       nsp1 .gt. nspx .or. nsp2 .gt. nspx) go to 30
      else
         nsp1 = nsp
         nsp2 = nsp
      end if
      atm1=atm(nsp1)
      atm2=atm(nsp2)
      print '(''  min and max distance  > '',$)'
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
                     dr(i) = (rr(1,na)-rr(1,nb))*bg(1,i) +              &
     &                       (rr(2,na)-rr(2,nb))*bg(2,i) +              &
     &                       (rr(3,na)-rr(3,nb))*bg(3,i)
                  enddo
                  do nn1=-2,2
                     dn1=dr(1)-nn1
                     do nn2=-2,2
                        dn2=dr(2)-nn2
                        do nn3=-2,2
                           dn3=dr(3)-nn3
                           dd = fact * sqrt(                            &
     &                                (dn1*at(1,1) +                    &
     &                                 dn2*at(1,2) +                    &
     &                                 dn3*at(1,3) )**2 +               &
     &                                (dn1*at(2,1) +                    &
     &                                 dn2*at(2,2) +                    &
     &                                 dn3*at(2,3) )**2 +               &
     &                                (dn1*at(3,1) +                    &
     &                                 dn2*at(3,2) +                    &
     &                                 dn3*at(3,3) )**2 )    
                           if(dd.ge.dmin.and.dd.le.dmax) then
                              ndist=ndist+1
                              if (ndist.gt.ndistx) stop ' ndist !'
                              atom1(ndist)=na
                              atom2(ndist)=nb
                              d(ndist)= dd
                              if (nn1.eq.0.and.nn2.eq.0.and.nn3.eq.0)   &
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
      index(1)=0.0
      if (ndist.gt.0) call hpsort(ndist,d,index)
!
      if (iout.eq.1)                                                    &
     &     write(iout,100) atm1, atm2, dmin,dmax
      do nd=1,ndist
         write(iout,200) atom1(index(nd)), atom2(index(nd)),            &
     &        other_cell(index(nd)), d(nd)
      end do
!
      go to 30
!
 20   print '(''  number of neighbors (max 4) for species '',           &
     &     a3,'' > '', $)', atm1
      read(5,*) nn
      if (nn.gt.nnx) stop ' too many neighbors' 
!
! look for nearest neighbors
!
      do na=1,nat
         if (atm(ityp(na)) .eq. atm1) then
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
                  dr(i)=(rr(1,na)-rr(1,nb))*bg(1,i) +                   &
     &                  (rr(2,na)-rr(2,nb))*bg(2,i) +                   &
     &                  (rr(3,na)-rr(3,nb))*bg(3,i)
               end do
               do nn1=-1,1
                  dn1=dr(1)-nn1
                  do nn2=-1,1
                     dn2=dr(2)-nn2
                     do nn3=-1,1
                        dn3=dr(3)-nn3
                        dd = fact * sqrt(                               &
     &                       (dn1*at(1,1) +                             &
     &                        dn2*at(1,2) +                             &
     &                        dn3*at(1,3) )**2 +                        &
     &                       (dn1*at(2,1) +                             &
     &                        dn2*at(2,2) +                             &
     &                        dn3*at(2,3) )**2 +                        &
     &                       (dn1*at(3,1) +                             &
     &                        dn2*at(3,2) +                             &
     &                        dn3*at(3,3) )**2 )       
                        do i=1,3
                           drv(i) = rr(i,na)-rr(i,nb) -                 &
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
            if (ndist.ne.nn) stop ' what?'
!
! calculate angles with nearest neighbors
!
            nd=0
            do nn1=1,nn
               do nn2=nn1+1,nn
                  nd=nd+1
                  angolo(nd) = 360/(2*pi) * acos ( fact**2 *            &
     &                 ( drn(1,nn1)*drn(1,nn2) +                        &
     &                   drn(2,nn1)*drn(2,nn2) +                        &
     &                   drn(3,nn1)*drn(3,nn2) ) / d(nn1) / d(nn2) )
               end do
            end do
            if (nd.ne.nn*(nn-1)/2) stop ' what??'
!
! dd is the distance from the origin
!
            dd = sqrt( rr(1,na)**2 + rr(2,na)**2 + rr(3,na)**2 ) * fact
            write(iout,250) atm(ityp(na)), na, (d(nn1),nn1=1,nn)
            write(iout,300) dd, (angolo(nn1),nn1=1,nn*(nn-1)/2)
         end if
      end do
!
      stop
!
 100  format(2x, ' species: ',a3,' - ',a3,3x,f6.2,' < D <',f6.2)
 200  format(2x, ' atoms:', 2i6,a1, ' distance =',f10.5,' A')
 250  format(2x,  a3, i3, '  d(neigh):',8f8.3)
 300  format(2x, ' d(O) :',f6.3,' A     angles  :',8f8.1)
!
!
      end
!
! NB: "capital" and "matches" are in Module/parser.f90
!     and in upftools/nclib.f90 as well !!!
!-----------------------------------------------------------------------
      function matches (string1, string2)  
!---------------------------------------------------------------------
!
      implicit none  
      character*(*) string1, string2  
      logical matches
      integer len1, len2, l  

      len1 = len(string1)  
      len2 = len(string2)  
      do l = 1, len2 - len1 + 1  
         if (string1 (1:len1) .eq.string2 (l:l + len1 - 1) ) then  
            matches = .true.  
            return  
         endif
      enddo
      matches = .false.  
      return  
      end
!
!-----------------------------------------------------------------------
      function capital (character)  
!---------------------------------------------------------------------
!
!   converts character to capital if lowercase
!   copy character to output in all other cases
!
      implicit none  
      character*1 capital, character
!
      character*26 minuscole/'abcdefghijklmnopqrstuvwxyz'/,             &
     &             maiuscole/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      integer  i
!
      do i=1,26
         if (character.eq.minuscole(i:i)) then
            capital=maiuscole(i:i)
            return
         end if
      end do
      capital = character  
!
      return  
      end

