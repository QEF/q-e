! dftd3 program for computing the dispersion energy and forces from cart
! and atomic numbers as described in
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
!
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456
! (for BJ-damping)
!
! Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germany
!
! Repackaging of the original code without any change in the functionality:
!
! Copyright (C) 2016, Bálint Aradi
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

module dftd3_extras
  use dftd3_common
  use dftd3_core
  implicit none


contains

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine printoptions
    write(*,*) 'dftd3 <coord filename> [-options]'
    write(*,*) 'options:'
    write(*,*) '-func <functional name in TM style>'
    write(*,*) '-grad'
    write(*,*) '-anal (pair analysis)'
    write(*,*) ' file <fragemt> with atom numbers'
    write(*,*) ' is read for a fragement based '
    write(*,*) ' analysis (one fragment per line,'
    write(*,*) ' atom ranges (e.g. 1-14 17-20) are allowed)'
    write(*,*) '-noprint'
    write(*,*) '-pbc (periodic boundaries; reads VASP-format)'
    write(*,*) '-abc (compute E(3))'
    write(*,*) '-cnthr (neglect threshold in Bohr for CN, default=40)'
    write(*,*) '-cutoff (neglect threshold in Bohr for E_disp, &
        & default=95)'
    write(*,*) '-old (DFT-D2)'
    write(*,*) '-zero (DFT-D3 original zero-damping)'
    write(*,*) '-bj (DFT-D3 with Becke-Johnson finite-damping)'
    write(*,*) '-zerom (revised DFT-D3 original zero-damping)'
    write(*,*) '-bjm (revised DFT-D3 with Becke-Johnson damping)'
    write(*,*) '-tz (use special parameters for TZ-type calculations)'
    write(*,*) 'variable parameters can be read from <current-director&
        &y>/.dftd3par.local'
    write(*,*) ' or '
    write(*,*) 'variable parameters read from ~/.dftd3par.<hostname>'
    write(*,*) 'if -func is used, -zero or -bj or -old is required!"'
    stop
  end subroutine printoptions


  subroutine rdpar(dtmp,version,s6,s18,rs6,rs18,alp)
    real(wp) s6,s18,rs6,rs18,alp
    integer version
    character*(*) dtmp
    character*80 ftmp,homedir
    logical ex
    real(wp) xx(10)
    integer nn
    ! initialize
    s6 =0
    s18=0
    rs6=0
    rs18=0
    alp =0
    ! read parameter file from current directory
    inquire(file='.dftd3par.local',exist=ex)
    if (ex)then
      write(dtmp,'(a)')'.dftd3par.local'
      open(unit=43,file=dtmp)
      read(43,'(a)',end=10)ftmp
      call readl(ftmp,xx,nn)
      if (nn.eq.6) then
        s6 =xx(1)
        rs6 =xx(2)
        s18 =xx(3)
        rs18=xx(4)
        alp =xx(5)
        version=idint(xx(6))
      end if
10    close(43)
      return
    end if
    ! read parameter file from home directory
    call system('hostname > .tmpx')
    open(unit=43,file='.tmpx')
    read(43,'(a)')ftmp
    close(43,status='delete')
    call get_environment_variable("HOME", homedir)
    write (*,*) trim(homedir)
    write(dtmp,'(a,''/.dftd3par.'',a)')trim(homedir),trim(ftmp)
    inquire(file=dtmp,exist=ex)
    if (ex)then
      open(unit=42,file=dtmp)
      read(42,'(a)',end=9)ftmp
      call readl(ftmp,xx,nn)
      if (nn.eq.6) then
        s6 =xx(1)
        rs6 =xx(2)
        s18 =xx(3)
        rs18=xx(4)
        alp =xx(5)
        version=idint(xx(6))
      end if
9     close(42)
    end if

  end subroutine rdpar


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! analyse all pairs
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine adisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
      & rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal, &
      & autoang,rthr,cn_thr,s6,s18,etot)
    integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp) xyz(3,*),r0ab(max_elem,max_elem),r2r4(*),s6
    real(wp) rs6,rs8,rs10,alp6,alp8,alp10,autokcal,etot,s18,autoang
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3),rcov(max_elem)

    integer iat,jat,i,j,k,nbin
    real(wp) R0,r,r2,r6,r8,tmp,alp,dx,dy,dz,c6,c8,c10
    real(wp) damp6,damp8,damp10,r42,rr,check,rthr,cn_thr,rvdw
    real(wp) cn(n),i6,e6,e8,e10,edisp
    real(wp) dist(0:15),li(0:15,2)
    real(wp) xx(500),eg(10000)
    integer grplist(500,20)
    integer grpn(20),at(n)
    integer ngrp,dash
    integer iiii, jjjj, iii, jjj, ii, jj, ni, nj
    integer iout(500)
    logical ex
    character*80 atmp

    real(wp),dimension(:,:), allocatable :: ed
    allocate(ed(n,n))


    ! distance bins
    li(0,1)=0
    li(0,2)=1.5
    li(1,1)=1.5
    li(1,2)=2
    li(2,1)=2
    li(2,2)=2.3333333333
    li(3,1)=2.3333333333
    li(3,2)=2.6666666666
    li(4,1)=2.6666666666
    li(4,2)=3.0
    li(5,1)=3.0
    li(5,2)=3.3333333333
    li(6,1)=3.3333333333
    li(6,2)=3.6666666666
    li(7,1)=3.6666666666
    li(7,2)=4.0
    li(8,1)=4.0
    li(8,2)=4.5
    li(9,1)=4.5
    li(9,2)=5.0
    li(10,1)=5.0
    li(10,2)=5.5
    li(11,1)=5.5
    li(11,2)=6.0
    li(12,1)=6.0
    li(12,2)=7.0
    li(13,1)=7.0
    li(13,2)=8.0
    li(14,1)=8.0
    li(14,2)=9.0
    li(15,1)=9.0
    li(15,2)=10.0
    nbin=15

    call ncoord(n,rcov,iz,xyz,cn,cn_thr)

    write(*,*)
    write(*,*)'analysis of pair-wise terms (in kcal/mol)'
    write(*,'(''pair'',2x,''atoms'',9x,''C6'',14x,''C8'',12x, &
        &''E6'',7x,''E8'',7x,''Edisp'')')
    e8=0
    ed=0
    dist=0
    check=0
    do iat=1,n-1
      do jat=iat+1,n

        dx=xyz(1,iat)-xyz(1,jat)
        dy=xyz(2,iat)-xyz(2,jat)
        dz=xyz(3,iat)-xyz(3,jat)
        r2=(dx*dx+dy*dy+dz*dz)
        !THR
        if (r2.gt.rthr) cycle
        r =sqrt(r2)
        R0=r0ab(iz(jat),iz(iat))
        rr=R0/r
        r6=r2**3

        if(version.eq.3)then
          ! DFT-D3 zero-damp
          tmp=rs6*rr
          damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
          tmp=rs8*rr
          damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
        else
          ! DFT-D3M zero-damp
          tmp=(r/(rs6*R0))+rs8*R0
          damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
          tmp=(r/R0)+rs8*R0
          damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
        endif

        if (version.eq.2)then
          c6=c6ab(iz(jat),iz(iat),1,1,1)
          damp6=1.d0/(1.d0+exp(-alp6*(r/(rs6*R0)-1.0d0)))
          e6 =s6*autokcal*c6*damp6/r6
          e8=0.0
        else
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat), &
              & cn(iat),cn(jat),c6)
        end if

        if((version.eq.3).or.(version.eq.5))then
          e6 =s6*autokcal*c6*damp6/r6
          r8 =r6*r2
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          c8 =3.0d0*c6*r42
          e8 =s18*autokcal*c8*damp8/r8
        end if

        if((version.eq.4).or.(version.eq.6))then
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          c8 =3.0d0*c6*r42
          ! use BJ radius
          R0=sqrt(c8/c6)
          rvdw=rs6*R0+rs8
          e6 =s6*autokcal*c6/(r6+rvdw**6)
          r8 =r6*r2
          e8 =s18*autokcal*c8/(r8+rvdw**8)
        end if

        edisp=-(e6+e8)
        ed(iat,jat)=edisp
        ed(jat,iat)=edisp

        write(*,'(2i4,2x,2i3,2D16.6,2F9.4,F10.5)') &
            & iat,jat,iz(iat),iz(jat),c6,c8, &
            & -e6,-e8,edisp

        check=check+edisp
        rr=r*autoang
        do i=0,nbin
          if (rr.gt.li(i,1).and.rr.le.li(i,2)) dist(i)=dist(i)+edisp
        end do
      end do
    end do

    write(*,'(/''distance range (Angstroem) analysis'')')
    write(*,'( ''writing histogram data to <histo.dat>'')')
    open(unit=11,file='histo.dat')
    do i=0,nbin
      write(*,'(''R(low,high), Edisp, %tot :'',2f4.1,F12.5,F8.2)') &
          & li(i,1),li(i,2),dist(i),100.*dist(i)/etot
      write(11,*)(li(i,1)+li(i,2))*0.5,dist(i)
    end do
    close(11)

    write(*,*) 'checksum (Edisp) ',check
    if (abs(check-etot).gt.1.d-3)stop'something is weired in adisp'

    inquire(file='fragment',exist=ex)
    if (.not.ex) return
    write(*,'(/''fragment based analysis'')')
    write(*,'( ''reading file <fragment> ...'')')
    open(unit=55,file='fragment')
    i=0
    at=0
111 read(55,'(a)',end=222) atmp
    call readfrag(atmp,iout,j)
    if (j.gt.0)then
      i=i+1
      grpn(i)=j
      do k=1,j
        grplist(k,i)=iout(k)
        at(grplist(k,i))=at(grplist(k,i))+1
      end do
    end if
    goto 111
222 continue
    ngrp=i
    k=0
    do i=1,n
      if (at(i).gt.1) stop 'something is weird in file <fragment>'
      if (at(i).eq.0)then
        k=k+1
        grplist(k,ngrp+1)=i
      end if
    end do
    if (k.gt.0) then
      ngrp=ngrp+1
      grpn(ngrp)=k
    end if
    ! Implemented display of atom ranges instead of whole list of atoms
    write(*,*)'group # atoms '
    dash=0
    do i=1,ngrp
      write(*,'(i4,3x,i4)',advance='no')i,grplist(1,i)
      do j=2,grpn(i)
        if (grplist(j,i).eq.(grplist(j-1,i)+1)) then
          if (dash.eq.0)then
            write(*,'(A1)',advance='no')'-'
            dash=1
          end if
        else
          if (dash.eq.1)then
            write(*,'(i4)',advance='no') grplist(j-1,i)
            dash=0
          end if
          write(*,'(i4)',advance='no') grplist(j,i)
        end if
      end do
      if (dash.eq.1)then
        write(*,'(i4)',advance='no') grplist(j-1,i)
        dash=0
      end if
      write(*,*)''
    end do

    ! old display list code
    ! write(*,*)'group # atoms '
    ! do i=1,ngrp
    ! write(*,'(i4,3x,100i3)')i,(grplist(j,i),j=1,grpn(i))
    ! end do

    eg=0
    iii=0
    do i=1,ngrp
      ni=grpn(i)
      iii=iii+1
      jjj=0
      do j=1,ngrp
        nj=grpn(j)
        jjj=jjj+1
        do ii=1,ni
          iiii=grplist(ii,i)
          do jj=1,nj
            jjjj=grplist(jj,j)
            if (jjjj.lt.iiii)cycle
            eg(lin(iii,jjj))=eg(lin(iii,jjj))+ed(iiii,jjjj)
          end do
        end do
      end do
    end do

    ! call prmat(6,eg,ngrp,0,'intra- + inter-group dispersion energies')
    write(*,*)' group i j Edisp'
    k=0
    check=0
    do i=1,ngrp
      do j=1,i
        k=k+1
        check=check+eg(k)
        write(*,'(5x,i4,'' --'',i4,F8.2)')i,j,eg(k)
      end do
    end do
    write(*,*) 'checksum (Edisp) ',check

  end subroutine adisp


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! load C6 coefficients from file
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine loadc6(fname,maxc,max_elem,c6ab,maxci)
    integer maxc,max_elem,maxci(max_elem)
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    character*(*) fname
    character*1 atmp
    character*80 btmp

    real(wp) x,y,f,cn1,cn2,cmax,xx(10)
    integer iat,jat,i,n,l,j,k,il,iadr,jadr,nn

    c6ab=-1
    maxci=0

    ! read file
    open(unit=1,file=fname)
    read(1,'(a)')btmp
10  read(1,*,end=11) y,iat,jat,cn1,cn2
    call limit(iat,jat,iadr,jadr)
    maxci(iat)=max(maxci(iat),iadr)
    maxci(jat)=max(maxci(jat),jadr)
    c6ab(iat,jat,iadr,jadr,1)=y
    c6ab(iat,jat,iadr,jadr,2)=cn1
    c6ab(iat,jat,iadr,jadr,3)=cn2

    c6ab(jat,iat,jadr,iadr,1)=y
    c6ab(jat,iat,jadr,iadr,2)=cn2
    c6ab(jat,iat,jadr,iadr,3)=cn1
    ! end if
    goto 10
11  continue
    close(1)

  end subroutine loadc6



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! load DFT-D2 data
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine loadoldpar(autoang,max_elem,maxc,c6ab,r0ab,c6)
    integer max_elem,maxc
    real(wp) r0ab(max_elem,max_elem)
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp) autoang

    real(wp) c6(86),r0(86)
    integer i,j

    ! the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799
    ! refer to the following values multiplied by 1.1 (rs6 in this code)
    ! H, He
    r0(1:86) = (/ 0.91d0,0.92d0, &
        & 0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, &
        & 1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0, &
        & 1.35d0,1.34d0, &
        & 1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
        & 1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
        & 1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0, &
        & 1.48d0,1.46d0, &
        & 1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
        & 1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
        & 1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0, &
        & 1.638d0,1.602d0,1.564d0,1.594d0,1.594d0,1.594d0,1.594d0, &
        & 1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0, &
        & 1.594d0,1.594d0,1.594d0, &
        & 1.625d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0, &
        & 1.611d0, &
        & 1.598d0,1.805d0,1.767d0,1.725d0,1.823d0,1.810d0,1.749d0/)
    ! Li-Ne
    ! Na-Ar
    ! K, Ca old
    ! Sc-Zn
    ! Ga-Kr
    ! Rb, Sr
    ! Y-Cd
    ! In, Sn, Sb, Te, I, Xe
    ! Cs,Ba,La,Ce-Lu
    ! Hf, Ta-Au
    ! Hg,Tl,Pb,Bi,Po,At,Rn

    c6(1:86) = (/0.14d0,0.08d0, &
        & 1.61d0,1.61d0,3.13d0,1.75d0,1.23d0,0.70d0,0.75d0,0.63d0, &
        & 5.71d0,5.71d0,10.79d0,9.23d0,7.84d0,5.57d0,5.07d0,4.61d0, &
        & 10.8d0,10.8d0,10.8d0,10.8d0,10.8d0, &
        & 10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,16.99d0, &
        & 17.10d0,16.37d0,12.64d0,12.47d0,12.01d0,24.67d0,24.67d0, &
        & 24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0, &
        & 24.67d0,24.67d0,24.67d0,37.32d0,38.71d0,38.44d0,31.74d0, &
        & 31.50d0,29.99d0,315.275d0,226.994d0,176.252d0, &
        & 140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0, &
        & 140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0, &
        & 105.112d0, &
        & 81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0, &
        & 57.364d0,57.254d0,63.162d0,63.540d0,55.283d0,57.171d0,56.64d0 /)

    c6ab = -1
    do i=1,86
      do j=1,i
        r0ab(i,j)=(r0(i)+r0(j))/autoang
        r0ab(j,i)=(r0(i)+r0(j))/autoang
        c6ab(i,j,1,1,1)=dsqrt(c6(i)*c6(j))
        c6ab(j,i,1,1,1)=dsqrt(c6(i)*c6(j))
      end do
    end do

  end subroutine loadoldpar

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! read atomic data (radii, r4/r2)
  ! not used
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine rdatpar(fname,max_elem,val)
    integer max_elem
    real(wp) val(max_elem)
    character*(*) fname

    integer i
    real(wp) dum1

    val = 0

    open(unit=142,file=fname)
502 read(142,*,end=602) i,dum1
    if (i.gt.0)then
      if (i.gt.max_elem)call stoprun('wrong cardinal number (rdatpar)')
      val(i)=dum1
    end if
    goto 502
602 close(142)

    do i=1,max_elem
      if (val(i).lt.1.d-6)then
        write(*,*)'atom ',i
        call stoprun( 'atomic value missing' )
      end if
    end do

  end subroutine rdatpar

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! read radii
  ! not used
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine rdab(fname,autoang,max_elem,ab)
    real(wp) autoang
    integer max_elem
    real(wp) ab(max_elem,max_elem)
    character*(*) fname

    integer i,j
    real(wp) dum1

    ab = 0

    open(unit=142,file=fname)
502 read(142,*,end=602) dum1,i,j
    if (i.gt.0.and.j.gt.0)then
      if (i.gt.max_elem) call stoprun( 'wrong cardinal number (rdab)')
      if (j.gt.max_elem) call stoprun( 'wrong cardinal number (rdab)')
      ab(i,j)=dum1/autoang
      ab(j,i)=dum1/autoang
    end if
    goto 502
602 close(142)

  end subroutine rdab

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! read coordinates in au or Ang. if its a xmol file
  ! redone by S.E. to avoid some input errors. Looks for $coord, ang, bohr
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine rdcoord(fname,n,xyz,iat,fix,fdum)
    real(wp) xyz(3,*)
    integer iat(*),n
    character*(*) fname
    !fix:array of fixed coordinates, fdum: whether
    logical fix(n),fdum

    real(wp) floats(3),f
    character*80 line
    character*80 strings(3)
    integer j,ich,cs,cf,ncheck

    f=0.0d0
    ich=142
    open(unit=ich,file=fname)
    ncheck=0
    rewind(ich)
    do
      read(ich,'(a)',end=200)line
      if (line.ne."") exit
    end do

    call readline(line,floats,strings,cs,cf)
    if (cf.eq.1.and.floats(1).gt.0) then
      f=1.0d0 / autoang
      read(ich,'(A)',end=200)line
    else if (index(line,'$coord').ne.0) then
      f=1.0d0
    else if (index(line,'ang').ne.0) then
      f=1.0d0 / autoang
    else if (index(line,'bohr').ne.0) then
      f=1.0d0
    end if
    if (f.lt.1.0d0) then
      call stoprun('Coordinate format not recognized!')
    end if
    do
      read(ich,'(a)',end=200)line
      if (index(line,'$redu').ne.0) exit
      if (index(line,'$user').ne.0) exit
      if (index(line,'$end' ).ne.0) exit
      call readline(line,floats,strings,cs,cf)
      if (cf.ne.3) cycle
      call elem(strings(1),j)
      if (j.eq.0) then
        fdum=.true.
        !ignores dummies and unknown elements
        cycle
      end if
      ncheck=ncheck+1
      xyz(1,ncheck)=floats(1)*f
      xyz(2,ncheck)=floats(2)*f
      xyz(3,ncheck)=floats(3)*f
      iat(ncheck)=j
      !fixes coordinate i
      if (strings(2).ne.'')fix(ncheck)=.true.
      ! write(*,321)floats(1:3),strings(1:3),j,fix(ncheck) !debug print

    end do

200 continue

    if (n.ne.ncheck) then
      write(*,*)n,'/=',ncheck
      call stoprun('error reading coord file')
    end if
    close(ich)

    ! 321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
  end subroutine rdcoord


  subroutine rdatomnumber(fname,n)
    integer n
    character*(*) fname

    real(wp) floats(3),f
    character*80 line
    character*80 strings(3)
    integer j,ich,cs,cf

    f=0.0d0
    ich=53
    open(unit=ich,file=fname)
    n=0
300 read(ich,'(a)',end=200)line
    if (line.eq."") goto 300
    call readline(line,floats,strings,cs,cf)
    if (cf.eq.1.and.floats(1).gt.0.and.cs.eq.0) then
      f=1.0d0 / autoang
      ! write(*,*)floats(1)
      n=int(floats(1))
      close(ich)
      return
    else if (index(line,'$coord').ne.0) then
      f=1.0d0
    else if (index(line,'ang').ne.0) then
      f=1.0d0 / autoang
    else if (index(line,'bohr').ne.0) then
      f=1.0d0
    end if
    if (f.lt.1.0d0) then
      call stoprun('Coordinate format not recognized!')
    end if
    do
      read(ich,'(a)',end=200)line
      if (index(line,'$redu').ne.0) exit
      if (index(line,'$user').ne.0) exit
      if (index(line,'$end' ).ne.0) exit
      call readline(line,floats,strings,cs,cf)
      if (cf.ne.3) exit
      call elem(strings(1),j)
      if (j.eq.0) cycle


      n=n+1
    end do

200 continue

    close(ich)

    ! 321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
  end subroutine rdatomnumber




  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine outg(nat,g,fname)
    integer nat,i
    real(wp) g(3,nat)
    character*(*) fname

    open(unit=142,file=fname)

    ! write(*,*)'Gradient : ', fname
    ! write(*,*)
    do i=1,nat
      write(142,'(3E22.14)')g(1:3,i)
      ! write(*,'(3D22.14)')g(1:3,i)
    end do

    close(142)

  end subroutine outg





  !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! add edisp in file energy
  ! and g to file gradient
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine wregrad(nat,xyz,iat,edisp,g)
    integer nat,iat(*)
    real(wp) edisp
    real(wp) g (3,*)
    real(wp) xyz(3,*)

    integer i,j,nn,nl
    character*128 a,a1
    character*20 fname
    real(wp) xx(10),gsum,x,y,z
    real(wp), dimension(:,:), allocatable :: gr
    logical ex

    allocate(gr(3,nat))

    fname='dftd3_gradient'
    inquire(file='gradient',exist=ex)
    if (.not.ex) then
      write(*,*) 'no gradient file found to add Gdisp!'
      write(*,*) 'hence written to file dftd3_gradient'
      call outg(nat,g,fname)
      return
    end if
    ! write file gradient
    j=0
    open(unit=42,file='gradient')
201 read(42,'(a)',end=301)a1
    j=j+1
    if (index(a1,'cycle').ne.0)nl=j
    goto 201
301 continue

    if (nl.lt.2)then
      write(*,*) 'illegal gradient file to add Gdisp!'
      write(*,*) 'hence written to file dftd3_gradient'
      call outg(nat,g,fname)
      return
    end if

    rewind 42
    do i=1,nl
      read(42,'(a)')a1
    end do
    do i=1,nat
      read(42,*)x,y,z
      if (abs(x-xyz(1,i)).gt.1.d-8) &
          & call stoprun( 'gradient read error x')
      if (abs(y-xyz(2,i)).gt.1.d-8) &
          & call stoprun( 'gradient read error y')
      if (abs(z-xyz(3,i)).gt.1.d-8) &
          & call stoprun( 'gradient read error z')
      xyz(1,i)=x
      xyz(2,i)=y
      xyz(3,i)=z
    end do
    gsum=0
    do i=1,nat
      read(42,*)gr(1,i),gr(2,i),gr(3,i)
      g(1:3,i)=g(1:3,i)+gr(1:3,i)
    end do
    gsum=sqrt(sum(g(1:3,1:nat)**2))

    rewind 42
    open(unit=43,file='gradient.tmp')
    j=0
401 read(42,'(a)',end=501)a1
    j=j+1
    if (j.lt.nl)then
      write(43,'(a)')trim(a1)
    else
      call readl(a1,xx,nn)
      j=idint(xx(1))
      write(43,'('' cycle = '',i6,4x,''SCF energy ='',F18.11,3x, &
          & ''|dE/dxyz| ='',F10.6)')j,xx(2)+edisp,gsum
      do i=1,nat
        write(43,'(3(F20.14,2x),4x,a2)')xyz(1,i),xyz(2,i),xyz(3,i), &
            & esym(iat(i))
      end do
      do i=1,nat
        write(43,'(3D22.13)')g(1,i),g(2,i),g(3,i)
      end do
      a='$end'
      write(43,'(a)')trim(a)
      goto 501
    end if
    goto 401
501 continue
    close(42)
    close(43)

    call system('mv gradient.tmp gradient')

    ! write file energy
    j=1
    open(unit=42,file='energy')
    open(unit=43,file='energy.tmp')
50  read(42,'(a)',end=100)a
    call readl(a,xx,nn)
    if (nn.gt.3)nl=j
    j=j+1
    goto 50
100 continue

    rewind 42
    j=0
20  read(42,'(a)',end=200)a
    j=j+1
    if (j.lt.nl)then
      write(43,'(a)')trim(a)
      call readl(a,xx,nn)
    else
      call readl(a,xx,nn)
      xx(2)=xx(2)+edisp
      write(43,'(i6,4F20.11)')idint(xx(1)),xx(2:nn)
      a='$end'
      write(43,'(a)')trim(a)
      goto 200
    end if
    goto 20
200 continue
    close(42)
    close(43)

    call system('mv energy.tmp energy')

  end subroutine wregrad


  ! *****************************************************************

  function ESYM(I)
    INTEGER :: I
    
    CHARACTER*2 ESYM
    CHARACTER*2 ELEMNT(94)
    DATA ELEMNT/'h ','he', &
        & 'li','be','b ','c ','n ','o ','f ','ne', &
        & 'na','mg','al','si','p ','s ','cl','ar', &
        & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
        & 'zn','ga','ge','as','se','br','kr', &
        & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
        & 'cd','in','sn','sb','te','i ','xe', &
        & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
        & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
        & 'au','hg','tl','pb','bi','po','at','rn', &
        & 'fr','ra','ac','th','pa','u ','np','pu'/
    ESYM=ELEMNT(I)
    RETURN
  end function ESYM


  ! *****************************************************************
  ! Reads a given line
  ! *****************************************************************

  subroutine READL(A1,X,N)
    CHARACTER*(*) A1
    REAL(WP), DIMENSION(*) :: X
    INTEGER :: N

    INTEGER :: I, IS, IE, IB
    I=0
    IS=1
10  I=I+1
    X(I)=READAA(A1,IS,IB,IE)
    if (IB.GT.0 .AND. IE.GT.0) then
      IS=IE
      GOTO 10
    end if
    N=I-1
    RETURN
  end subroutine READL

  ! *****************************************************************
  ! this seems to be part of an old MOPAC code, extracts strings and n
  ! *****************************************************************

  function READAA(A,ISTART,IEND,IEND2)
    CHARACTER*(*) A
    INTEGER ISTART, IEND, IEND2
    REAL(WP) READAA

    INTEGER :: NINE, IZERO, MINUS, IdoT, ND, NE, IBL, IDIG, NL
    INTEGER :: I, II, J, N, M
    REAL(WP) :: C1, C2, ONE, X
    
    NINE=ICHAR('9')
    IZERO=ICHAR('0')
    MINUS=ICHAR('-')
    IdoT=ICHAR('.')
    ND=ICHAR('D')
    NE=ICHAR('E')
    IBL=ICHAR(' ')
    IEND=0
    IEND2=0
    IDIG=0
    C1=0
    C2=0
    ONE=1.D0
    X = 1.D0
    NL=LEN(A)
    do J=ISTART,NL-1
      N=ICHAR(A(J:J))
      M=ICHAR(A(J+1:J+1))
      if (N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IdoT) then
        GOTO 20
      end if
      if (N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO &
          & .OR. M.EQ.IdoT)) then
        GOTO 20
      end if
    end do
    READAA=0.D0
    RETURN
20  CONTINUE
    IEND=J
    do I=J,NL
      N=ICHAR(A(I:I))
      if (N.LE.NINE.AND.N.GE.IZERO) then
        IDIG=IDIG+1
        if (IDIG.GT.10) GOTO 60
        C1=C1*10+N-IZERO
      ELSEif (N.EQ.MINUS.AND.I.EQ.J) then
        ONE=-1.D0
      ELSEif (N.EQ.IdoT) then
        GOTO 40
      ELSE
        GOTO 60
      end if
    end do
40  CONTINUE
    IDIG=0
    do II=I+1,NL
      N=ICHAR(A(II:II))
      if (N.LE.NINE.AND.N.GE.IZERO) then
        IDIG=IDIG+1
        if (IDIG.GT.10) GOTO 60
        C2=C2*10+N-IZERO
        X = X /10
      ELSEif (N.EQ.MINUS.AND.II.EQ.I) then
        X=-X
      ELSE
        GOTO 60
      end if
    end do
    !
    ! PUT THE PIECES TOGETHER
    !
60  CONTINUE
    READAA= ONE * ( C1 + C2 * X)
    do J=IEND,NL
      N=ICHAR(A(J:J))
      IEND2=J
      if (N.EQ.IBL)RETURN
      if (N.EQ.ND .OR. N.EQ.NE)GOTO 57
    end do
    RETURN

57  C1=0.0D0
    ONE=1.0D0
    do I=J+1,NL
      N=ICHAR(A(I:I))
      IEND2=I
      if (N.EQ.IBL)GOTO 70
      if (N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
      if (N.EQ.MINUS)ONE=-1.0D0
    end do
61  CONTINUE
70  READAA=READAA*10**(ONE*C1)
    RETURN
  end function READAA

  subroutine prmat(iuout,r,n,m,head)
    integer :: iuout, n, m
    real(wp) r
    character*(*) head
    dimension r(*)

    integer :: nkpb, ibl, ir, j1, k1s, kd, j2, i, k, k1, k2, kk, j, i1, i2, ij
    ! subroutine prints matrix r,which is supposed
    ! to have dimension n,m when m is nonzero and
    ! ((n+1)*n)/2 when m is zero

    write(iuout,1001) head
    nkpb=10
    if (m)10,10,80
    !
10  continue
    ibl=n/nkpb
    ir=n-ibl*nkpb
    j1=1
    k1s=1
    kd=0
    if (ibl.eq.0) go to 50
    j2=nkpb
    do i=1,ibl
      write(iuout,1002)(j,j=j1,j2)
      k1=k1s
      k2=k1
      kk=0
      do j=j1,j2
        write(iuout,1003)j,(r(k),k=k1,k2)
        kk=kk+1
        k1=k1+kd+kk
        k2=k1+kk
      end do
      j1=j1+nkpb
      if (j1.gt.n) return
      j2=j2+nkpb
      k2=k1-1
      k1=k2+1
      k2=k1+(nkpb-1)
      k1s=k2+1
      kk=kd+nkpb
      do j=j1,n
        write(iuout,1003)j,(r(k),k=k1,k2)
        kk=kk+1
        k1=k1+kk
        k2=k2+kk
      end do
      kd=kd+nkpb
    end do
50  if (ir.eq.0) go to 70
    k1=k1s
    j2=j1+ir-1
    kk=0
    k2=k1
    write(iuout,1002)(j,j=j1,j2)
    write(iuout,1003)
    do j=j1,j2
      write(iuout,1003)j,(r(k),k=k1,k2)
      kk=kk+1
      k1=k1+kd+kk
      k2=k1+kk
    end do
70  return
80  ibl=m/nkpb
    ir=m-ibl*nkpb
    i2=0
    k2=0
    if (ibl.eq.0) go to 100
    do i=1,ibl
      i1=(i-1)*n*nkpb+1
      i2=i1+(nkpb-1)*n
      k1=k2+1
      k2=k1+(nkpb-1)
      write(iuout,1002)(k,k=k1,k2)
      do j=1,n
        write(iuout,1003)j,(r(ij),ij=i1,i2,n)
        i1=i1+1
90      i2=i1+(nkpb-1)*n
      end do
    end do
100 if (ir.eq.0) go to 120
    i1=ibl*n*nkpb+1
    i2=i1+(ir-1)*n
    k1=k2+1
    k2=m
    write(iuout,1002)(k,k=k1,k2)
    write(iuout,1003)
    do j=1,n
      write(iuout,1003)j,(r(ij),ij=i1,i2,n)
      i1=i1+1
      i2=i1+(ir-1)*n
    end do

120 write(iuout,1003)
    return
1001 format(/,a)
1002 format(/,' ',4x,10(3x,i4,3x),/)
1003 format(' ',i4,10f10.3)
  end subroutine prmat

  !cccccccccccccccccccccccccccccccccccccccccccccc
  ! readfrag: will read a list of numbers c
  ! from character. c
  ! line: string containing the numbers c
  ! iout: array which returns integer c
  ! numbers c
  ! n: number of returned integers c
  ! S.E., 17.02.2011 c
  !cccccccccccccccccccccccccccccccccccccccccccccc

  subroutine readfrag(line,iout,n)
    character*80 line
    character*12 str1,str2
    integer iout(500)
    integer*4 n,i,j,k,sta,sto
    character*11 nums

    ! write(*,*) 'In readfrag:'
    ! write(*,*) 'Line reads: ',line
    n=0
    iout=0
    str1=''
    str2=''
    nums='0123456789-'
    do i=1,80
      ! Check for allowed characters (part of nums) and add to number-string s
      ! If char NOT allowed, its a separator and str1 will be processed

      if (index(nums,line(i:i)).ne.0) then
        str1=trim(str1)//line(i:i)
      else

        if (str1.ne.'') then

          ! If str1 is simple number, cast to integer and put in iout. increase n.
          if (index(str1,'-').eq.0) then
            n=n+1
            read(str1,'(I4)') iout(n)
            str1=''

          end if
          ! If str1 contains '-', its treated as a number range.
          if (index(str1,'-').ne.0) then

            do k=1,12
              ! Determine beginning number
              if (str1(k:k).ne.'-') then
                str2=trim(str2)//str1(k:k)
              end if
              ! '-' Marks end of beginning number. cast to integer and store in sta
              if (str1(k:k).eq.'-') then
                read(str2,'(I4)') sta
                str2=''
              end if

            end do
            ! Get the rest, store in sto
            read(str2,'(I4)') sto
            str2=''
            ! Write all numbers between sta and sto to iout and increase n
            do k=sta,sto
              n=n+1
              iout(n)=k
            end do
            str1=''
          end if
        end if
      end if

    end do

  end subroutine readfrag

  ! Input Geometry sanity check via CNs, not used
  subroutine checkcn(n,iz,cn,c6ab,max_elem,maxc)

    integer iz(*),i,n
    logical check
    real(wp) cn(*),maxcn
    integer max_elem,maxc
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)

    check=.false.
    do i=1,n
      if (iz(i).lt.10) then
        if (iz(i).ne.2) then
          maxcn=maxval(c6ab(iz(i),1,1:5,1,2))
          if (cn(i).gt.maxcn*2.0) then
            check=.true.
          end if
        end if
      end if
    end do
    if (check) then
      write(*,*)'----------------------------------------------------'
      write(*,*)'!!!! SOME CN ARE VERY LARGE. CHECK COORDINATES !!!!'
      write(*,*)'----------------------------------------------------'
    end if
  end subroutine checkcn

  ! Input Geometry sanity check (to avoid au/Angtstrom mixups) S.E. 16
  subroutine checkrcov(n,iz,rcov,xyz)
    logical check
    integer iz(*),n,i,j
    real(wp) rcov(94),dist,dx,dy,dz,thr,xyz(3,*),r
    check=.false.
    do i=1,n-1
      do j=i+1,n
        dx=xyz(1,i)-xyz(1,j)
        dy=xyz(2,i)-xyz(2,j)
        dz=xyz(3,i)-xyz(3,j)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        thr=0.6*(rcov(iz(i))+rcov(iz(j)))
        if (r.lt.thr) then
          check=.true.
        end if
      end do
    end do
    if (check) then
      write(*,*)'--------------------------------------------------'
      write(*,*)'!! SOME DISTANCES VERY SHORT. CHECK COORDINATES !!'
      write(*,*)'--------------------------------------------------'
    end if
  end subroutine checkrcov




  !reads a line cuts the at blanks and tabstops and returns all floats and
  subroutine readline(line,floats,strings,cs,cf)
    real(wp) floats(3)
    character*80 line
    character*80 strings(3)

    real(wp) num
    character*80 stmp,str
    character*1 digit
    integer i,ty,cs,cf

    stmp=''
    cs=1
    cf=1
    strings(:)=''
    do i=1,len(trim(line))
      digit=line(i:i)
      if (digit.ne.' '.and.digit.ne.char(9)) then
        stmp=trim(stmp)//trim(digit)
      elseif (stmp.ne.'')then
        call checktype(stmp,num,str,ty)
        if (ty.eq.0) then
          floats(cf)=num
          cf=cf+1
        elseif (ty.eq.1) then
          strings(cs)=str
          cs=cs+1
        else
          write(*,*)'Problem in checktype, must abort'
          exit
        end if
        stmp=''
      end if
      if (i.eq.len(trim(line))) then
        call checktype(stmp,num,str,ty)
        if (ty.eq.0) then
          floats(cf)=num
          cf=cf+1
        elseif (ty.eq.1) then
          strings(cs)=str
          cs=cs+1
        else
          write(*,*)'Problem in checktype, must abort'
          exit
        end if
        stmp=''
      end if
    end do
    cs=cs-1
    cf=cf-1
  end subroutine readline


  !this checks the type of the string and returns it cast to real or as st
  subroutine checktype(field,num,str,ty)
    character*80 field,str,tstr
    real(wp) num
    integer ty

    integer i,e
    logical is_num

    ty=99
    str=''
    is_num=.false.
    read(field,'(F20.10)',IOSTAT=e)num
    if (e.eq.0)is_num=.true.
    if (field.eq.'q'.or.field.eq.'Q')is_num=.false.
    if (field.eq.'e'.or.field.eq.'E')is_num=.false.
    if (field.eq.'d'.or.field.eq.'D')is_num=.false.
    if (is_num)then
      if (index(field,'.').ne.0) then
        read(field,'(F20.10)')num
        ty=0
      else
        str=trim(field)//'.0'
        read(str,'(F20.10)')num
        str=''
        ty=0
      end if
    else
      str=field
      ty=1
    end if

  end subroutine checktype




  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! B E G I N O F P B C P A R T
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! read coordinates in Angst and converts them to au
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcrdcoord(fname,lattice,n,xyz,iat,autoang)
    real(wp) :: xyz(3,*)
    real(wp), INTENT(OUT) ::lattice(3,3)
    integer, INTENT(out) :: iat(*)
    integer, INTENT(in) :: n
    character*(*), INTENT(IN) :: fname
    logical :: selective=.FALSE.
    logical :: cartesian=.TRUE.
    real(wp), INTENT(IN) ::autoang

    real(wp) xx(10),scalar
    character*200 line
    character*80 args(90),args2(90)

    integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck


    lattice=0

    ich=142
    open(unit=ich,file=fname)
    rewind(ich)
    ncheck=0
    ntype=0
    read(ich,'(a)',end=200)line
    call parse(line,' ',args,ntype)
    read(ich,'(a)',end=200)line
    call readl(line,xx,nn)
    scalar=xx(1)/autoang
    ! write(*,'(F8.6)')scalar
    do i=1,3
      read(ich,'(a)',end=200)line
      call readl(line,xx,nn)
      if (nn < 3) call stoprun( 'Error reading unit cell vectors' )
      lattice(1,i)=xx(1)*scalar
      lattice(2,i)=xx(2)*scalar
      lattice(3,i)=xx(3)*scalar
    end do
    read(ich,'(a)',end=200)line
    line=adjustl(line)
    call readl(line,xx,nn)
    if (nn.eq.0) then
      call parse(line,' ',args,ntype)
      read(ich,'(a)',end=200)line
      line=adjustl(line)
      call readl(line,xx,nn)
    end if
    ! call elem(args(1),i_dummy2)
    ! if (i_dummy2<1 .OR. i_dummy2>94) then
    ! args=args2
    ! end if
    if (nn.NE.ntype ) then
      call stoprun( 'Error reading number of atomtypes')
    end if
    ncheck=0
    do i=1,nn
      i_dummy1=INT(xx(i))
      call elem(args(i),i_dummy2)
      if (i_dummy2<1 .OR. i_dummy2>94) &
          & call stoprun( 'Error: unknown element.')
      do j=1,i_dummy1
        ncheck=ncheck+1
        iat(ncheck)=i_dummy2
      end do
    end do
    if (n.ne.ncheck) call stoprun('Error reading Number of Atoms')

    read(ich,'(a)',end=200)line
    line=adjustl(line)
    if (line(:1).EQ.'s' .OR. line(:1).EQ.'S') then
      selective=.TRUE.
      read(ich,'(a)',end=200)line
      line=adjustl(line)
    end if

    ! write(*,*)line(:1)
    cartesian=(line(:1).EQ.'c' .OR. line(:1).EQ.'C' .OR.&
        &line(:1).EQ.'k' .OR. line(:1).EQ.'K')
    do i=1,n
      read(ich,'(a)',end=200)line
      call readl(line,xx,nn)
      if (nn.NE.3) call stoprun( 'Error reading coordinates.')

      if (cartesian) then
        xyz(1,i)=xx(1)*scalar
        xyz(2,i)=xx(2)*scalar
        xyz(3,i)=xx(3)*scalar
      ELSE
        xyz(1,i)=lattice(1,1)*xx(1)+lattice(1,2)*&
            & xx(2)+lattice(1,3)*xx(3)
        xyz(2,i)=lattice(2,1)*xx(1)+lattice(2,2)*xx(2)+lattice(2,3)*&
            & xx(3)
        xyz(3,i)=lattice(3,1)*xx(1)+lattice(3,2)*xx(2)+lattice(3,3)*&
            & xx(3)
      end if

      ! write(*,'(3F20.10,1X,I3)')xyz(:,i),iat(i) !debug printout

    end do


200 continue

    close(ich)
  end subroutine pbcrdcoord



  subroutine pbcrdatomnumber(fname,n)
    integer, INTENT(out) :: n
    character*(*), INTENT(IN) :: fname
    logical :: selective=.FALSE.
    logical :: cartesian=.TRUE.

    real(wp) xx(10),scalar,fdum
    character*80 line,args(90),args2(90)

    integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2

    ich=142
    open(unit=ich,file=fname)
    n=0
    ntype=0
    read(ich,'(a)',end=200)line
    call parse(line,' ',args,ntype)
    read(ich,'(a)',end=200)line
    call readl(line,xx,nn)
    ! write(*,'(F8.6)')scalar
    do i=1,3
      read(ich,'(a)',end=200)line
      call readl(line,xx,nn)
      if (nn < 3) call stoprun( 'Error reading unit cell vectors' )
    end do
    read(ich,'(a)',end=200)line
    line=adjustl(line)
    call readl(line,xx,nn)
    if (nn.eq.0) then
      call parse(line,' ',args,ntype)
      read(ich,'(a)',end=200)line
      line=adjustl(line)
      call readl(line,xx,nn)
    end if
    ! call elem(args(1),i_dummy2)
    ! if (i_dummy2<1 .OR. i_dummy2>94) then
    ! args=args2
    ! end if
    if (nn.NE.ntype ) then
      ! if (nn.NE.ntype2) then
      call stoprun( 'Error reading number of atomtypes')
      ! ELSE
      ! ntype=ntype2
      ! end if
    end if
    n=0
    do i=1,nn
      i_dummy1=INT(xx(i))
      n=n+i_dummy1
    end do

200 continue

    close(ich)
  end subroutine pbcrdatomnumber





  ! Input Geometry sanity check for pbc (to avoid au/Angtstrom mixups)
  subroutine pbccheckrcov(n,iz,rcov,xyz,lat)
    logical check
    integer iz(*),n,i,j,taux,tauy,tauz
    real(wp) rcov(94),dist,dx,dy,dz,thr,xyz(3,*),r,lat(3,3),tau(3)
    check=.false.
    do i=1,n-1
      do j=i+1,n
        thr=0.5*(rcov(iz(i))+rcov(iz(j)))
        do taux=-1,1
          do tauy=-1,1
            do tauz=-1,1
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              dx=xyz(1,i)-xyz(1,j)+tau(1)
              dy=xyz(2,i)-xyz(2,j)+tau(2)
              dz=xyz(3,i)-xyz(3,j)+tau(3)

              r=sqrt(dx*dx+dy*dy+dz*dz)
              if (r.lt.thr) then
                check=.true.
                ! write(*,*)'short distance',i,'(',iz(i),') and ',j,'(',iz(j),'):'
                ! write(*,*)r,' < ',thr
                ! write(*,*)
              end if
            end do
          end do
        end do
      end do
    end do
    if (check) then
      write(*,*)'--------------------------------------------------'
      write(*,*)'!! SOME DISTANCES VERY SHORT. CHECK COORDINATES !!'
      write(*,*)'--------------------------------------------------'
    end if
  end subroutine pbccheckrcov



  subroutine pbcwregrad(nat,g,g_lat)
    integer nat,i
    real(wp) g(3,nat)
    real(wp) g_lat(3,3)

    open(unit=142,file='dftd3_gradient')

    ! write(*,*)'Gradient:' !Jonas
    ! write(*,*) !Jonas
    do i=1,nat
      write(142,'(3E22.14)')g(1:3,i)
      ! write(*,'(3D22.14)')g(1:3,i) !Jonas
    end do

    close(142)

    open(unit=143,file='dftd3_cellgradient')

    ! write(*,*)'Gradient:' !Jonas
    ! write(*,*) !Jonas
    do i=1,3
      write(143,'(3E22.14)')g_lat(1:3,i)
      ! write(*,'(3D22.14)')g(1:3,i) !Jonas
    end do

    close(143)
  end subroutine pbcwregrad

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! analyse all pairs
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcadisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
      & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal,&
      & autoang,rthr,rep_v,cn_thr,rep_cn,s6,s18,etot,lat)
    integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp) xyz(3,*),r0ab(max_elem,max_elem),r2r4(*),s6
    real(wp) rs6,rs8,rs10,alp6,alp8,alp10,autokcal,etot,s18,autoang
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3),rcov(max_elem)
    real(wp) lat(3,3)
    integer rep_v(3),rep_cn(3)

    integer iat,jat,i,j,k,nbin,taux,tauy,tauz
    real(wp) R0,r,r2,r6,r8,tmp,alp,dx,dy,dz,c6,c8,c10
    real(wp) damp6,damp8,damp10,r42,rr,check,rthr,cn_thr,rvdw
    real(wp) cn(n),i6,e6,e8,e10,edisp
    real(wp),allocatable :: dist(:),li(:,:)
    real(wp) xx(500),eg(10000)
    integer grplist(500,20)
    integer grpn(20),at(n)
    integer ngrp,dash
    integer iiii, jjjj, iii, jjj, ii, jj, ni, nj
    integer iout(500)
    logical ex
    character*80 atmp
    real(wp) tau(3)

    real(wp),dimension(:,:), allocatable :: ed
    allocate(ed(n,n))


    ! distance bins
    nbin=17
    allocate(dist(0:nbin))
    allocate(li(0:nbin,2))

    li(0,1)=0
    li(0,2)=1.5
    li(1,1)=1.5
    li(1,2)=2
    li(2,1)=2
    li(2,2)=2.3333333333
    li(3,1)=2.3333333333
    li(3,2)=2.6666666666
    li(4,1)=2.6666666666
    li(4,2)=3.0
    li(5,1)=3.0
    li(5,2)=3.3333333333
    li(6,1)=3.3333333333
    li(6,2)=3.6666666666
    li(7,1)=3.6666666666
    li(7,2)=4.0
    li(8,1)=4.0
    li(8,2)=4.5
    li(9,1)=4.5
    li(9,2)=5.0
    li(10,1)=5.0
    li(10,2)=5.5
    li(11,1)=5.5
    li(11,2)=6.0
    li(12,1)=6.0
    li(12,2)=7.0
    li(13,1)=7.0
    li(13,2)=8.0
    li(14,1)=8.0
    li(14,2)=9.0
    li(15,1)=9.0
    li(15,2)=10.0
    li(16,1)=10.0
    li(16,2)=20.0
    li(17,1)=20.0
    li(17,2)=dsqrt(rthr)*autoang


    call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,cn_thr)

    write(*,*)
    write(*,*)'analysis of pair-wise terms (in kcal/mol)'
    write(*,'(''pair'',2x,''atoms'',9x,''C6'',14x,''C8'',12x,&
        &''E6'',7x,''E8'',7x,''Edisp'')')
    e8=0
    ed=0
    dist=0
    check=0
    do iat=1,n
      do jat=iat,n

        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
              dx=xyz(1,iat)-xyz(1,jat)+tau(1)
              dy=xyz(2,iat)-xyz(2,jat)+tau(2)
              dz=xyz(3,iat)-xyz(3,jat)+tau(3)
              r2=(dx*dx+dy*dy+dz*dz)
              !THR
              if (r2.gt.rthr.or.r2.lt.0.5) cycle
              r =sqrt(r2)
              R0=r0ab(iz(jat),iz(iat))
              rr=R0/r
              r6=r2**3

              if(version.eq.3)then
                tmp=rs6*rr
                damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
                tmp=rs8*rr
                damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
              else
                tmp=(r/(R0*rs6)+R0*rs8)**(-alp6)
                damp6 =1.d0/( 1.d0+6.d0*tmp )
                tmp=(r/(R0)+R0*rs8)**(-alp8)
                damp8 =1.d0/( 1.d0+6.d0*tmp )
              endif

              if (version.eq.2)then
                c6=c6ab(iz(jat),iz(iat),1,1,1)
                damp6=1.d0/(1.d0+exp(-alp6*(r/(rs6*R0)-1.0d0)))
                if (iat.eq.jat) then
                  e6 =s6*autokcal*c6*damp6/r6
                else
                  e6 =s6*autokcal*c6*damp6/r6
                end if
                e8=0.0d0
              else
                call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
                    & cn(iat),cn(jat),c6)
              end if

              if((version.eq.3).or.(version.eq.5))then
                r8 =r6*r2
                r42=r2r4(iz(iat))*r2r4(iz(jat))
                c8 =3.0d0*c6*r42
                if (iat.eq.jat) then
                  e6 =s6*autokcal*c6*damp6/r6*0.5
                  e8 =s18*autokcal*c8*damp8/r8*0.5
                else
                  e6 =s6*autokcal*c6*damp6/r6
                  e8 =s18*autokcal*c8*damp8/r8
                end if
              end if

              if((version.eq.4).or.(version.eq.6))then
                r42=r2r4(iz(iat))*r2r4(iz(jat))
                c8 =3.0d0*c6*r42
                ! use BJ radius
                R0=dsqrt(c8/c6)
                rvdw=rs6*R0+rs8
                r8 =r6*r2
                if (iat.eq.jat) then
                  e6 =s6*autokcal*c6/(r6+rvdw**6)*0.5
                  e8 =s18*autokcal*c8/(r8+rvdw**8)*0.5
                else
                  e6 =s6*autokcal*c6/(r6+rvdw**6)
                  e8 =s18*autokcal*c8/(r8+rvdw**8)
                end if
              end if

              edisp=-(e6+e8)
              ed(iat,jat)=edisp
              ed(jat,iat)=edisp

              ! write(*,'(2i4,2x,2i3,2D16.6,2F9.4,F10.5)')
              ! . iat,jat,iz(iat),iz(jat),c6,c8,
              ! . -e6,-e8,edisp

              check=check+edisp
              rr=r*autoang
              do i=0,nbin
                if (rr.gt.li(i,1).and.rr.le.li(i,2)) dist(i)=dist(i)+edisp
              end do
            end do
          end do
        end do
      end do
    end do

    write(*,'(/''distance range (Angstroem) analysis'')')
    write(*,'( ''writing histogram data to <histo.dat>'')')
    open(unit=11,file='histo.dat')
    do i=0,nbin
      write(*,'(''R(low,high), Edisp, %tot :'',2f5.1,F12.5,F8.2)')&
          & li(i,1),li(i,2),dist(i),100.*dist(i)/etot
      write(11,*)(li(i,1)+li(i,2))*0.5,dist(i)
    end do
    close(11)

    write(*,*) 'checksum (Edisp) ',check
    if (abs(check-etot).gt.1.d-3)stop'something is weired in adisp'

    deallocate(dist,li)
    return








    inquire(file='fragment',exist=ex)
    if (ex) return
    write(*,'(/''fragment based analysis'')')
    write(*,'( ''reading file <fragment> ...'')')
    open(unit=55,file='fragment')
    i=0
    at=0
111 read(55,'(a)',end=222) atmp
    call readfrag(atmp,iout,j)
    if (j.gt.0)then
      i=i+1
      grpn(i)=j
      do k=1,j
        grplist(k,i)=iout(k)
        at(grplist(k,i))=at(grplist(k,i))+1
      end do
    end if
    goto 111
222 continue
    ngrp=i
    k=0
    do i=1,n
      if (at(i).gt.1) stop 'something is weird in file <fragment>'
      if (at(i).eq.0)then
        k=k+1
        grplist(k,ngrp+1)=i
      end if
    end do
    if (k.gt.0) then
      ngrp=ngrp+1
      grpn(ngrp)=k
    end if
    ! Implemented display of atom ranges instead of whole list of atoms
    write(*,*)'group # atoms '
    dash=0
    do i=1,ngrp
      write(*,'(i4,3x,i4)',advance='no')i,grplist(1,i)
      do j=2,grpn(i)
        if (grplist(j,i).eq.(grplist(j-1,i)+1)) then
          if (dash.eq.0)then
            write(*,'(A1)',advance='no')'-'
            dash=1
          end if
        else
          if (dash.eq.1)then
            write(*,'(i4)',advance='no') grplist(j-1,i)
            dash=0
          end if
          write(*,'(i4)',advance='no') grplist(j,i)
        end if
      end do
      if (dash.eq.1)then
        write(*,'(i4)',advance='no') grplist(j-1,i)
        dash=0
      end if
      write(*,*)''
    end do

    ! old display list code
    ! write(*,*)'group # atoms '
    ! do i=1,ngrp
    ! write(*,'(i4,3x,100i3)')i,(grplist(j,i),j=1,grpn(i))
    ! end do

    eg=0
    iii=0
    do i=1,ngrp
      ni=grpn(i)
      iii=iii+1
      jjj=0
      do j=1,ngrp
        nj=grpn(j)
        jjj=jjj+1
        do ii=1,ni
          iiii=grplist(ii,i)
          do jj=1,nj
            jjjj=grplist(jj,j)
            if (jjjj.lt.iiii)cycle
            eg(lin(iii,jjj))=eg(lin(iii,jjj))+ed(iiii,jjjj)
          end do
        end do
      end do
    end do

    ! call prmat(6,eg,ngrp,0,'intra- + inter-group dispersion energies')
    write(*,*)' group i j Edisp'
    k=0
    check=0
    do i=1,ngrp
      do j=1,i
        k=k+1
        check=check+eg(k)
        write(*,'(5x,i4,'' --'',i4,F8.2)')i,j,eg(k)
      end do
    end do
    write(*,*) 'checksum (Edisp) ',check

    deallocate(dist,li)
  end subroutine pbcadisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  REAL(WP) function volume(lat)
    REAL(WP), INTENT(in) ::lat(3,3)
    REAL(WP) zwerg

    zwerg=lat(1,1)*lat(2,2)*lat(3,3)+lat(1,2)*lat(2,3)*lat(3,1)+&
        & lat(1,3)*lat(2,1)*lat(3,2)-lat(1,3)*lat(2,2)*lat(3,1)-&
        & lat(1,2)*lat(2,1)*lat(3,3)-lat(1,1)*lat(2,3)*lat(3,2)
    volume=abs(zwerg)
  end function volume




  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! string pars procedures
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine parse(str,delims,args,nargs)
    character(len=*),intent(inout) :: str
    character(len=*),intent(in) :: delims
    character(len=len_trim(str)) :: strsav
    character(len=*),dimension(:),intent(inout) :: args
    integer, intent(out) :: nargs

    integer :: na, i, lenstr, k

    strsav=str
    call compact(str)
    na=size(args)
    do i=1,na
      args(i)=' '
    end do
    nargs=0
    lenstr=len_trim(str)
    if (lenstr==0) return
    k=0

    do
      if (len_trim(str) == 0) exit
      nargs=nargs+1
      call split(str,delims,args(nargs))
      call removebksl(args(nargs))
    end do
    str=strsav

  end subroutine parse

  !**********************************************************************

  subroutine compact(str)

    ! Converts multiple spaces and tabs to single spaces; deletes control ch
    ! removes initial spaces.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str)):: outstr

    integer :: lenstr, isp, k, i, ich

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    isp=0
    k=0

    do i=1,lenstr
      ch=str(i:i)
      ich=iachar(ch)

      select case(ich)

      case(9,32)
        if (isp==0) then
          k=k+1
          outstr(k:k)=' '
        end if
        isp=1

      case(33:)
        k=k+1
        outstr(k:k)=ch
        isp=0

      end select

    end do

    str=adjustl(outstr)

  end subroutine compact

  !**********************************************************************

  subroutine removesp(str)


    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    integer :: lenstr, k, i, ich

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0

    do i=1,lenstr
      ch=str(i:i)
      ich=iachar(ch)
      select case(ich)
      case(0:32)
        cycle
      case(33:)
        k=k+1
        outstr(k:k)=ch
      end select
    end do

    str=adjustl(outstr)

  end subroutine removesp



  subroutine split(str,delims,before,sep)

    ! Routine finds the first instance of a character from 'delims' in the
    ! the string 'str'. The characters before the found delimiter are
    ! output in 'before'. The characters after the found delimiter are
    ! output in 'str'. The optional output character 'sep' contains the
    ! found delimiter. A delimiter in 'str' is treated like an ordinary
    ! character if it is preceded by a backslash (\). If the backslash
    ! character is desired in 'str', then precede it with another backslash.

    character(len=*),intent(inout) :: str,before
    character(len=*),intent(in) :: delims
    character,optional :: sep
    logical :: pres
    character :: ch,cha

    integer :: lenstr, k, ibsl, i, ipos, iposa

    pres=present(sep)
    str=adjustl(str)
    call compact(str)
    lenstr=len_trim(str)
    if (lenstr == 0) return
    k=0
    ibsl=0
    before=' '
    do i=1,lenstr
      ch=str(i:i)
      if (ibsl == 1) then
        k=k+1
        before(k:k)=ch
        ibsl=0
        cycle
      end if
      if (ch == '\') then
        k=k+1
        before(k:k)=ch
        ibsl=1
        cycle
      end if
      ipos=index(delims,ch)
      if (ipos == 0) then
        k=k+1
        before(k:k)=ch
        cycle
      end if
      if (ch /= ' ') then
        str=str(i+1:)
        if (pres) sep=ch
        exit
      end if
      cha=str(i+1:i+1)
      iposa=index(delims,cha)
      if (iposa > 0) then
        str=str(i+2:)
        if (pres) sep=cha
        exit
      else
        str=str(i+1:)
        if (pres) sep=ch
        exit
      end if
    end do
    if (i >= lenstr) str=''
    str=adjustl(str)
    return

  end subroutine split

  !**********************************************************************

  subroutine removebksl(str)

    ! Removes backslash (\) characters. Double backslashes (\\) are replaced
    ! by a single backslash.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    integer :: lenstr, k, ibsl, i

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0
    ibsl=0

    do i=1,lenstr
      ch=str(i:i)
      if (ibsl == 1) then
        k=k+1
        outstr(k:k)=ch
        ibsl=0
        cycle
      end if
      if (ch == '\') then
        ibsl=1
        cycle
      end if
      k=k+1
      outstr(k:k)=ch
    end do

    str=adjustl(outstr)

  end subroutine removebksl


  ! Uncomment this, if you work with the nagfor compiler
  !subroutine system(command)
  !  character(*), intent(in) :: command
  !
  !  call execute_command_line(command)
  !
  !end subroutine system


end module dftd3_extras
