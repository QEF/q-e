!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program special_points
  !-----======================--------------------------------------------
  !
  !     calculates special points for any structure,
  !     the default definition for the mesh is a shift of 1/(2n_i)
  !     where the length of b_i is equal to 1
  !_______________________________________________________________________
  !
  use kinds,     only: dp
  use cell_base, only: at, bg
  use symm_base, only: set_sym_bl, s, nrot
  implicit none
  integer, parameter :: nptx=20000
  character(len=30) :: filout
  character(len=1)  :: answer
  real(dp) ::  celldm(6), xk(3,nptx), xkw(nptx), omega
  integer  ::  k(3,nptx), kw(nptx), ieq(nptx), i,j,l, n1,n2,n3
  integer  ::  ibrav, nmax(3), nshift(3), nstart(3),n,n6,nf,nk,nptot
  logical  ::  aflag, sflag
  !
  write(*,1)
1 format(/,5x,'***************************************************',/,&
           5x,'*                                                 *',/,&
           5x,'*       Welcome to the special points world!      *',/,&
           5x,'*________________________________________________ *',/,&
           5x,'*    1 = cubic p (sc )      8 = orthor p (so )    *',/,&
           5x,'*    2 = cubic f (fcc)      9 = orthor base-cent. *',/,&
           5x,'*    3 = cubic i (bcc)     10 = orthor face-cent. *',/,&
           5x,'*    4 = hex & trig p      11 = orthor body-cent. *',/,&
           5x,'*    5 = trigonal   r      12 = monoclinic  p     *',/,&
           5x,'*    6 = tetrag p (st )    13 = monocl base-cent. *',/,&
           5x,'*    7 = tetrag i (bct)    14 = triclinic   p     *',/,&
           5x,'***************************************************',/ )
  !
  !.....default values
  ! 
  celldm(1)=1.d0
  do i=1,3
     nshift(i)=0
  enddo
  !
  write(*,'(5x,a)', advance="no") 'bravais lattice  >> '
  read(*,*) ibrav
  !   
  write(*,'(5x,a)',advance="no") 'filout [mesh_k]  >> '
  read(*,'(a)') filout
  if (filout.eq.' ') filout='mesh_k'
  open(unit=1,file=filout,status='unknown')
  open(unit=2,file='info',status='unknown')
  !
  if(ibrav.eq.4 .or. ibrav.gt.5) then
     write(*,'(5x,a)',advance="no") 'enter celldm(3)  >> '
     read(*,*) celldm(3)
  end if
  if(ibrav.ge.8) then
     write(*,'(5x,a)',advance="no") 'enter celldm(2)  >> '
     read(*,*) celldm(2)
  end if
  if(ibrav.eq.5 .or. ibrav.ge.12) then
     write(*,'(5x,a)',advance="no") 'enter celldm(4)  >> '
     read(*,*) celldm(4)
  end if
  if(ibrav.eq.14) then
     write(*,'(5x,a)')   'enter celldm(5)  >> cos(ac)'
     write(*,'(5x,a)',advance="no") 'enter celldm(5)  >> '
     read(*,*) celldm(5)
     write(*,'(5x,a)')   'enter celldm(6)  >> cos(ab)'
     write(*,'(5x,a)',advance="no") 'enter celldm(6)  >> '
     read(*,*) celldm(6)
  end if
  !
  write(*,'(5x,a)',advance="no") 'mesh: n1 n2 n3   >> '
  read(*,*) nmax
  nptot=nmax(1)*nmax(2)*nmax(3)
  if(nptot.gt.nptx) then
     write(*,'(5x,i6)') nptx
     call errore('kpoints','nptx too small for this mesh',1)
  endif
  write(*,'(5x,a)',advance="no") 'mesh: k1 k2 k3 (0 no shift, 1 shifted) >> '
  read(*,*) nshift(1), nshift(2), nshift(3)
  !
  write(*,'(5x,a)',advance="no") 'write all k? [f] >> '
  read(*,'(a1)') answer
  aflag= answer.eq.'t'.or.answer.eq.'T' .or.                        &
         answer.eq.'y'.or.answer.eq.'Y' .or.                        &
         answer.eq.'1'
  !
  call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  !
  ! normalize at to celldm(1) ( a0 for cubic lattices )
  !
  do i = 1, 3
     at( i, 1 ) = at( i, 1 ) / celldm( 1 )
     at( i, 2 ) = at( i, 2 ) / celldm( 1 )
     at( i, 3 ) = at( i, 3 ) / celldm( 1 )
  enddo
  !
  call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  ! 
  write(2,'(2x,''crystal axis:  ''/3(2x,''('',3f7.4,'')  ''/) )') &
                 ((at(i,j), i=1,3), j=1,3)
  write(2,'(2x,''reciprocal axis:  ''/3(2x,''('',3f7.4,'')  ''/) )') &
                 ((bg(i,j), i=1,3), j=1,3)
  write(2,*)' Omega (in a^3 units) = ',omega
  !
  !.......................................................................
  !
  call set_sym_bl ( )
  !
  write(2,'(//,1x,i3,2x,a19)') nrot,'symmetry operations' 
  do n6=0,(nrot-1)/6
     nf=min(nrot-6*n6,6)
     write(2,'(1x)')
     do i=1,3
        write(2,'(6(3i3,2x))') ((s(i,j,n6*6+n), j=1,3), n=1,nf)
     end do
  end do
  !
  sflag=.false.
  do i=1,3 
     ! shifted grid
     if(nshift(i).eq.1) then
        nshift(i)=2
        nmax(i)=nshift(i)*nmax(i)
        nstart(i)=1
        sflag=.true.
     else
        ! unshifted grid
        nstart(i)=0
        nshift(i)=1
     end if
  enddo
  !
  n=0
  do n3=nstart(3),nmax(3)-1,nshift(3)
     do n2=nstart(2),nmax(2)-1,nshift(2)
        do n1=nstart(1),nmax(1)-1,nshift(1)
           n=n+1
           k(1,n)=n1
           k(2,n)=n2
           k(3,n)=n3
           kw(n)=1
           ieq(n)=0
           call check(n,k,kw,ieq,s,nrot,nmax) 
        enddo
     enddo
  enddo
  !
  nk=0
  write(2,'(/)')
  do j=1,n
     if(kw(j).gt.0.or.aflag) then
        nk=nk+1
        xkw(nk)=kw(j)
        do l=1,3
           xk(l,nk)=0.d0
           do i=1,3
              xk(l,nk)=xk(l,nk)+k(i,j)*bg(l,i)/nmax(i)
           enddo
        end do
        write(2,2) j,k(1,j),k(2,j),k(3,j),kw(j),ieq(j)
2       format('  k(',i3,')=( ',i2,' ',i2,' ',i2,' ) --- weight=',  &
                      i3,' |folds in point #',i3)
     endif
  enddo
  !
  write(*,'(/5x,a)',advance="no") '# of k-points   == '
  write(*,'(i5,a5,i5)') nk,'  of ',n
  write(*,'(2x)')
  !
  write(1,'(i5)') nk
  do j=1,nk
     if(aflag.and.kw(j).eq.0) then
        write(1,'(i5,1x,3f11.7,f7.2,i4)') j,(xk(l,j),l=1,3),xkw(j),ieq(j)
     else
        write(1,'(i5,1x,3f11.7,f7.2)') j,(xk(l,j),l=1,3),xkw(j)
     end if
  end do
  !
  if(.not.sflag.and.kw(1).ne.1) then
     write(*,'(5x,a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,'(5x,a)') '!the considered mesh has not the correct symmetry!!'
     write(*,'(5x,a/)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  endif
  !
  close(unit=1)
  close(unit=2)
  !
end program special_points
!
!-----------------------------------------------------------------------
subroutine check(n,k,kw,ieq,s,nrot,nmax)
  !-----------------------------------------------------------------------
  !
  integer k(3,n),kw(n), s(3,3,nrot),kr(3),ieq(n),nmax(3)
  logical flag
  !
  irot=1
  flag=.true.
  do while(irot.le.nrot.and.flag)
     kr(1)=0
     kr(2)=0
     kr(3)=0
     call ruotaijk ( s(1,1,irot),k(1,n),k(2,n),k(3,n),kr(1),kr(2),kr(3) )
     do j=1,3
        do while(kr(j).ge.nmax(j)) 
           kr(j)=kr(j)-nmax(j)
        enddo
        do while(kr(j).le.-1) 
           kr(j)=kr(j)+nmax(j)
        enddo
     enddo
     np=1
     do while(flag.and.np.le.n-1)
        if( kr(1).eq.k(1,np) .and. &
            kr(2).eq.k(2,np) .and. &
            kr(3).eq.k(3,np) ) then
           kw(n)=0
           naux =np
           do while(kw(naux).eq.0)
              naux=ieq(naux)
           enddo
           ieq(n)=naux
           kw(naux)=kw(naux)+1
           flag=.false. 
        endif
        np=np+1
     enddo
     irot=irot+1
  enddo
  !
  return
end subroutine check
!
!-----------------------------------------------------------------------
subroutine ruotaijk(s,i,j,k,ri,rj,rk)
  !-----------------------------------------------------------------------
  !
  implicit real*8 (a-h, o-z)
  integer  s(3,3),i,j,k,ri,rj,rk
  !
  ri=s(1,1)*i+s(1,2)*j+s(1,3)*k
  rj=s(2,1)*i+s(2,2)*j+s(2,3)*k
  rk=s(3,1)*i+s(3,2)*j+s(3,3)*k
  !
  return
end subroutine ruotaijk
