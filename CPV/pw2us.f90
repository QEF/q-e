!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

module ncpp
  character(len=50) dft
  integer lmax, lloc
  logical nlcc
  real(kind=8) zp
  real(kind=8), allocatable::  vnl(:,:), rho_atc(:), rho_at(:)
  integer, allocatable:: lchi(:)
  character(len=2)  psd
end module ncpp

module bhs
  integer nlc, nnl
  logical bhstype
  real(kind=8)  alpc(2), cc(2), alps(3,0:3), aps(6,0:3)
  real(kind=8)  a_nlcc, b_nlcc, alpha_nlcc
end module bhs

module grid
  real(kind=8)  zmesh, xmin, dx
  real(kind=8), allocatable:: r(:), rab(:)
  integer mesh
end module grid

module wavefunctions
  integer nchi
  real(kind=8), allocatable:: chi(:,:),  oc(:)
end module wavefunctions

module ultrasoft
  integer nvales, nang, nbeta, kkbeta, ifpcor, keyps
  real(kind=8) z, exfact
  real(kind=8) ,allocatable:: betar(:), dion(:,:), ddd(:,:), &
       qqq(:,:), qfunc(:,:,:), vloc0(:)
  integer nnlz, ifqopt, nqf, irel, iptype, npf
  real(kind=8) etotpseu, wwnl, ee, rinner, eloc, dummy, rc, rcloc, eee
  character(len=20) title
end module ultrasoft


program pw2us

  implicit none
  character(len=32) file_pw, file_us
  logical pwformat, bhsformat

5 print '(''Input PP file in PW format > '',$)'
  read '(a)', file_pw
  pwformat=file_pw.ne.' '
  if (.not. pwformat) then
     print '(''Input PP file in BHS format > '',$)'
     read '(a)', file_pw
     bhsformat=file_pw.ne.' '
     if (.not.bhsformat) stop
  else
     inquire (file=file_pw,exist=pwformat)
     if (.not. pwformat) go to 5
  end if

  open(unit=1,file=file_pw,status='old',form='formatted')
  if (pwformat) then
     call read_ncpp
  else
     call read_bhs
  end if

  close (unit=1)

  file_us=trim(file_pw)//'.van'
  print '(''Output PP file in US format :  '',a)', file_us

  open(unit=2,file=file_us,status='unknown',form='formatted')
  call write_us
  close (unit=2)

end program pw2us

subroutine read_ncpp

  use ncpp
  use bhs
  use grid
  use wavefunctions

  implicit none
  real(kind=8) x, erf
  integer nb, ios, i, l, n, ir
  logical numeric
  external erf
  
  read( 1, '(a)', end=300, err=300, iostat=ios ) dft
  if (dft(1:2).eq.'**') dft='PZ'

  read ( 1, *, err=300, iostat=ios ) psd, zp, lmax, nlc, nnl, nlcc, &
       lloc, bhstype
  if ( nlc.gt.2 .or. nnl.gt.3) &
       call error( 'read_ncpp','Wrong nlc or nnl',1 )
  if ( nlc* nnl .lt. 0 ) &
       call error( 'read_ncpp','nlc*nnl < 0 ? ',1 )
  if ( zp.le.0d0 ) &
      call error( 'read_ncpp','Wrong zp ',1 )
  if ( lmax.gt.3.or.lmax.lt.0 ) &
      call error( 'read_ncpp','Wrong lmax ',1 )

  if (lloc.eq.-1000) lloc=lmax
  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric = nlc.le.0 .and. nnl.le.0

  if (.not.numeric) then
     !
     !   read here pseudopotentials in analytic form
     !
     read( 1, *, err=300, iostat=ios ) &
          ( alpc(i), i=1, 2 ), ( cc(i), i=1,2 )
     if ( abs(cc(1)+cc(2)-1.d0).gt.1.0d-6) &
          call error ('read_ncpp','wrong pseudopotential coefficients',1)
     do l = 0, lmax
        read ( 1, *, err=300, iostat=ios ) &
             ( alps(i,l),i=1,3 ), (aps(i,l),i=1,6)
     enddo

     if (nlcc) then
        read( 1, *, err=300, iostat=ios ) &
             a_nlcc, b_nlcc, alpha_nlcc
        if (alpha_nlcc.le.0.d0) &
             call error('read_ncpp','nlcc but alpha=0',1)
     end if

     if (bhstype) call bachel(alps,aps,1,lmax)

  end if
  read( 1, *, err=300, iostat=ios ) zmesh, xmin, dx, mesh, nchi
  if ( mesh.le.0) call error( 'read_ncpp', 'mesh too small', 1)
  if ( (nchi.lt.lmax   .and. lloc.eq.lmax).or.          &
       (nchi.lt.lmax+1 .and. lloc.ne.lmax)     )        &
       call error( 'read_ncpp', 'wrong no. of wfcts', 1 )
  !
  !    compute the radial mesh
  !
  allocate(  r(mesh))
  allocate(rab(mesh))
  !  r  (0) = 0.d0
  !  rab(0) = 0.d0
  do ir = 1, mesh
     x = xmin + float(ir-1) * dx
     r  (ir) = exp(x) / zmesh
     rab(ir) = dx * r(ir)
  end do

  allocate(vnl(mesh,0:lmax))

  if (numeric) then
     !
     !  Here pseudopotentials in numeric form are read
     !
     do l = 0, lmax
        read( 1, '(a)', err=300, iostat=ios )
        read( 1, *, err=300, iostat=ios )  (vnl(ir,l),ir=1,mesh)
        !        vnl(0,l) = (r(2)*vnl(1,l)-r(1)*vnl(2,l))/(r(2)-r(1))
     enddo

     if(nlcc) then
        allocate(rho_atc(mesh))
        read( 1, *, err=300, iostat=ios ) ( rho_atc(ir), ir=1,mesh )
     endif

  else
     !
     !  convert analytic to numeric form
     !
     do l=0,lmax
        !        vnl(0,l) = - ( cc(1)*sqrt(alpc(1))+cc(2)*sqrt(alpc(2)))*zp
        do ir=1,mesh
           vnl(ir,l)= - ( cc(1)*erf(sqrt(alpc(1))*r(ir)) +     &
                          cc(2)*erf(sqrt(alpc(2))*r(ir)) ) * zp/r(ir)
        end do

        do n=1,nnl
           vnl(:,l)= vnl(:,l)+ (aps(n,l)+ aps(n+3,l)*r(:)**2 )* &
                exp(-alps(n,l)*r(:)**2)
        end do
        !
        ! convert to Rydberg
        !
        vnl(:,l) = vnl(:,l)*2.0
     end do

     if (nlcc) then
        allocate(rho_atc(mesh))
        rho_atc(:) =(a_nlcc+b_nlcc*r(:)**2)*exp(-alpha_nlcc*r(:)**2)
     end if
  endif
  !
  ! subract the local part
  !
   do l = 0, lmax
     if ( l.ne.lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
  enddo
  !
  ! Here pseudowavefunctions in numeric form are read
  !
  allocate(lchi(nchi))
  allocate(oc(nchi))
  allocate(chi(mesh,nchi))
  do nb = 1, nchi
     read( 1, '(a)', err=300, iostat=ios )
     read( 1, *, err=300, iostat=ios ) lchi( nb), oc( nb )
     !
     !     Test lchi and occupation numbers
     !
     if ( nb.le.lmax.and.lchi(nb)+1.ne.nb) &
          call error('read_ncpp','order of wavefunctions',nb)
     if (lchi(nb).gt.lmax .or. lchi(nb).lt.0) &
          call error('read_ncpp','wrong lchi',nb)
     if ( oc(nb).lt.0.d0 .or.            &
          oc(nb).gt.2.d0*(2*lchi(nb)+1)) &
             call error('read_ncpp','wrong oc',nb)
     read( 1, *, err=300, iostat=ios ) & 
          (chi(ir,nb),ir=1,mesh)
  enddo
  !
  !    compute the atomic charges
  !
  allocate(rho_at(mesh))
  rho_at(:)=0.d0
  do nb = 1, nchi
     if( oc(nb).ne.0.d0) &
          rho_at(:) = rho_at(:) + oc(nb)*chi(:,nb)**2
  end do
  return

300 call error('read_ncpp','pseudo file is empty or wrong',1)

end subroutine read_ncpp

subroutine read_bhs

  use ncpp
  use bhs
  use grid
  use wavefunctions

  implicit none
  real(kind=8) erf
  integer nb, ios, is, i, n, l, ir
  external erf

  read ( 1, *, err=300, iostat=ios ) zmesh, zp, lmax, lloc
  if ( zmesh.le.0.or.zmesh.gt.120.0 ) &
      call error( 'read_bhs','Wrong z ',1 )
  if ( zp.le.0.or.zp.gt.25.0 ) &
      call error( 'read_bhs','Wrong zp ',1 )
  if ( lloc.gt.3.or.lloc.lt.0 ) &
      call error( 'read_bhs','Wrong lloc ',1 )

  read( 1, *, err=300, iostat=ios ) &
          ( cc(i), alpc(i), i=1,2 )
  if ( abs(cc(1)+cc(2)-1.d0).gt.1.0d-6) &
       call error ('read_bhs','wrong pseudopotential coefficients',1)
  do l = 0, lmax
     read ( 1, *, err=300, iostat=ios ) &
          ( alps(i,l), aps(i,l), aps(i+3,l), i=1,3)
  enddo

  read(1, *, err=300, iostat=ios) mesh, dx

  allocate(  r(mesh))
  allocate(rab(mesh))
  !
  ! Here pseudowavefunctions in numeric form are read
  !
  nchi=lmax
  allocate(chi(mesh,nchi))
  allocate(lchi(nchi))
  allocate(oc(nchi))
  do nb = 1, nchi
     lchi(nb)=nb
     oc(nb)=0.0
     if (nb.ne.1) read(1, *, err=300, iostat=ios) mesh, dx
     do ir=1,mesh
        read( 1, *, err=300, iostat=ios ) is, r(ir), chi(ir,nb)
        if (is.ne.ir) &
             call error('read_bhs','error 1 at line',ir)
     end do
  enddo

  allocate(rho_at(mesh))
  rho_at(:)=0.d0
  !
  !    compute the radial mesh derivative
  !
  do ir = 1, mesh
     rab(ir) = log(dx) * r(ir)
  end do
  !
  allocate(vnl(mesh,0:lmax))
  !
  !  convert analytic to numeric form
  !
  do l=0,lmax
     !        vnl(0,l) = - ( cc(1)*sqrt(alpc(1))+cc(2)*sqrt(alpc(2)))*zp
     do ir=1,mesh
        vnl(ir,l)= - ( cc(1)*erf(sqrt(alpc(1))*r(ir)) +     &
                       cc(2)*erf(sqrt(alpc(2))*r(ir)) ) * zp/r(ir)
     end do
     do n=1,3
        vnl(:,l)= vnl(:,l)+ (aps(n,l)+ aps(n+3,l)*r(:)**2 )* &
             exp(-alps(n,l)*r(:)**2)
     end do
     !
     ! convert to Rydberg
     !
     vnl(:,l) = vnl(:,l)*2.0
  end do
  !
  ! subract the local part
  !
  do l = 0, lmax
     if ( l.ne.lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
  enddo

  return

300 call error('read_bhs','pseudo file is empty or wrong',1)

end subroutine read_bhs

subroutine write_us

  use ncpp
  use grid
  use wavefunctions
  use ultrasoft

  implicit none
  real(kind=8), parameter:: pi=3.141592653589793d0
  integer iv, jv, l, lp, ir, convert_dft
  external convert_dft

  write(2,'(6i5)') 7,3,2,0,0,0
  ! iver, idmy
 
 
  title=psd//'                  '
  print '(''assumed title: '',a)',psd

  z=zmesh
  exfact=convert_dft(dft)
  write(2, '(a20,3f15.9)' ) title, z, zp, exfact

  nvales=nchi
  etotpseu=0
  write(2, '(2i5,1pe19.11)' ) nvales,mesh,etotpseu

  ee  =0.0
  do iv=1,nchi
     nnlz=10*lchi(iv)
     wwnl=oc(iv)
     write(2, '(i5,2f15.9)') nnlz, wwnl, ee
  end do

  keyps=0
  if (nlcc) then
     ifpcor=1
  else
     ifpcor=0
  end if
  rinner=0.0
  write(2, '(2i5,f15.9)' ) keyps, ifpcor, rinner

! note that in the Vanderbilt program l runs from 1 to lmax+1
  if (lloc.eq.lmax) then
     nang=lmax
  else
     nang=lmax+1
  end if
  eloc=0.0
  ifqopt=0
  nqf=0
  dummy=0.0
  write(2, '(2i5,f9.5,2i5,f9.5)' ) nang, lloc, eloc, ifqopt, nqf, dummy

  write(2, * ) (rinner, lp=1,2*nang-1)

  if (z.ge.19) then
     irel=1
  else
     irel=0
  end if
  write(2,'(i5)') irel

  rc=0.0
  write(2,'(1p4e19.11)') ( rc, l=1,nang )

  nbeta=lmax
! the definition of nbeta takes into account the absence of l=lloc
! from the projectors
  kkbeta=mesh
  write (2,'(2i5)' ) nbeta, kkbeta

  allocate(betar(kkbeta))
  allocate(qfunc(kkbeta,nbeta,nbeta))
  allocate(dion(nbeta,nbeta))
  allocate(ddd (nbeta,nbeta))
  allocate(qqq (nbeta,nbeta))
  iv=0
  do l=0,lmax
     if (l.ne.lloc) then
        iv=iv+1
        write(2, '(i5)') l
           
        eee=0.0
        betar(:)=chi(:,l+1)*vnl(:,l)
        write(2, '(1p4e19.11)') eee, ( betar(ir), ir=1,kkbeta )
           
        do jv=iv,nbeta
           if (jv.eq.iv) then
              betar(:)=chi(:,l+1)*betar(:)
              call simpson(mesh,betar,rab,dion(iv,jv))
              dion(iv,jv)=1.0/dion(iv,jv)
           else
              dion(iv,jv)= 0.0
           end if
           ddd(iv,jv)= 0.0
           qqq(iv,jv)= 0.0
           qfunc(:,iv,jv)=0.0
           write(2, '(1p4e19.11)' ) dion(iv,jv), &
                ddd(iv,jv), qqq(iv,jv), &
                (qfunc(ir,iv,jv),ir=1,kkbeta)
        enddo
     end if
  enddo
 
  iptype=0
  write(2, '(6i5)' ) (iptype, iv=1,nbeta)
     
  npf=0
  write(2, '(i5,f15.9)') npf, dummy
     
  rcloc=0.0
  allocate(vloc0(mesh))
  vloc0(:)=r(:)*vnl(:,lloc)
  write(2, '(1p4e19.11)' ) rcloc, ( vloc0(ir), ir=1,mesh )

  if ( ifpcor.eq.1 ) then 
     write(2, '(1p4e19.11)') dummy
     ! Vanderbilt rho_atc(r) =  4pi*r^2*rho_atc(r) PWSCF
     write(2, '(1p4e19.11)') ( rho_atc(ir)*r(ir)**2*4.d0*pi, ir=1,mesh )
  endif

  vloc0(:)=0.0
  write(2, '(1p4e19.11)')  (vloc0(ir), ir=1,mesh)

  write(2, '(1p4e19.11)')  (rho_at(ir), ir=1,mesh)
 
  write(2, '(1p4e19.11)')  (r(ir),ir=1,mesh)
  write(2, '(1p4e19.11)')  (rab(ir),ir=1,mesh)

  write(2, '(i5)') nvales

  write(2, '(1p4e19.11)') ((chi(ir,iv),ir=1,mesh),iv=1,nvales)

  return
end subroutine write_us


subroutine simpson( mesh, func, rab, intg )

  implicit none
  integer mesh
  real(kind=8)  func(mesh), rab(mesh), intg

  real(kind=8) c(4)
  integer I

  if ( mesh .lt. 8 ) call error('simpson','few mesh points',8)

  c(1) = 109.0 / 48.d0
  c(2) = -5.d0 / 48.d0
  c(3) = 63.d0 / 48.d0
  c(4) = 49.d0 / 48.d0 

  intg = ( func(1)*rab(1) + func(mesh  )*rab(mesh  ) )*c(1) &
       + ( func(2)*rab(2) + func(mesh-1)*rab(mesh-1) )*c(2) &
       + ( func(3)*rab(3) + func(mesh-2)*rab(mesh-2) )*c(3) &
       + ( func(4)*rab(4) + func(mesh-3)*rab(mesh-3) )*c(4)
  do i=5,mesh-4
     intg = intg + func(i)*rab(i)
  end do

  return
end subroutine simpson

subroutine error(a,b,n)
  character(len=*) a,b

  write(6,'(//'' program '',a,'':'',a,''.'',8x,i8,8x,''stop'')') a,b,n
  stop
end subroutine error

!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
real(kind=8) function erf (x)  
  !---------------------------------------------------------------------
  !
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  !
  implicit none  
  real(kind=8) :: x, x2, p1 (4), q1 (4), erfc  
  external erfc  
  data p1 / 2.42667955230532d2, 2.19792616182942d1, &
       6.99638348861914d0, - 3.56098437018154d-2 /
  data q1 / 2.15058875869861d2, 9.11649054045149d1, &
       1.50827976304078d1, 1.00000000000000d0 /
  !
  if (abs (x) .gt.6.d0) then  
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1 with 16-byte words
     !
     erf = sign (1.d0, x)  
  else  
     if (abs (x) .le.0.47d0) then  
        x2 = x**2  
        erf = x * (p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 ( &
             4) ) ) ) / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 ( &
             4) ) ) )
     else  
        erf = 1.d0 - erfc (x)  
     endif
  endif
  !
  return  
end function erf
!
!---------------------------------------------------------------------
real(kind=8) function erfc (x)  
  !---------------------------------------------------------------------
  !
  !     erfc(x) = 1-erf(x)  - See comments in erf
  !
  implicit none  
  real(kind=8) :: x, ax, x2, xm2, erf, p2 (8), q2 (8), p3 (5), q3 (5), &
       pim1
  external erf  
  data p2 / 3.00459261020162d2, 4.51918953711873d2, &
       3.39320816734344d2, 1.52989285046940d2, 4.31622272220567d1, &
       7.21175825088309d0, 5.64195517478974d-1, - 1.36864857382717d-7 /
  data q2 / 3.00459260956983d2, 7.90950925327898d2, &
       9.31354094850610d2, 6.38980264465631d2, 2.77585444743988d2, &
       7.70001529352295d1, 1.27827273196294d1, 1.00000000000000d0 /
  data p3 / - 2.99610707703542d-3, - 4.94730910623251d-2, - &
       2.26956593539687d-1, - 2.78661308609648d-1, - 2.23192459734185d-2 &
       /
  data q3 / 1.06209230528468d-2, 1.91308926107830d-1, &
       1.05167510706793d0, 1.98733201817135d0, 1.00000000000000d0 /

  data pim1 / 0.564189583547756d0 /  
  !        ( pim1= sqrt(1/pi) )
  ax = abs (x)  
  if (ax.gt.26.d0) then  
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     erfc = 0.d0  
  elseif (ax.gt.4.d0) then  
     x2 = x**2  
     xm2 = (1.d0 / ax) **2  
     erfc = (1.d0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  elseif (ax.gt.0.47d0) then  
     x2 = x**2  
     erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  else  
     erfc = 1.d0 - erf (ax)  
  endif
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  if (x.lt.0.d0) erfc = 2.d0 - erfc  
  !
  return  
end function erfc
!---------------------------------------------------------------------
real(kind=8) function freq (x)  
  !---------------------------------------------------------------------
  !
  !     freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
  !             - See comments in erf
  !
  real(kind=8) :: x, c, erf, erfc  
  external erf  

  data c / 0.707106781186548d0 /  
  !        ( c= sqrt(1/2) )
  freq = 0.5d0 * erfc ( - x * c)  
  !
  return  
end function freq

subroutine bachel(alps,aps,npseu,lmax)

  implicit none

  integer  npseu, lmax(npseu)
  real(kind=8)   alps(3,0:3,npseu), aps(6,0:3,npseu)

  integer  np, lmx, l, i, j, k, ia, ka, nik
  real(kind=8), parameter:: pi=3.141592653589793d0
  real(kind=8)    s(6,6), alpl, alpi, ail

  do np=1,npseu
     lmx=lmax(np)
     do l=0,lmx
        do k=1,6
           ka= mod(k-1,3)+1
           alpl= alps(ka,l,np)
           do i=1,k
              ia= mod(i-1,3)+1
              alpi= alps(ia,l,np)
              ail=alpi+alpl
              s(i,k)= sqrt(pi/ail)/4.d0/ail
              nik=int((k-1)/3)+int((i-1)/3)+1
              do j=2, nik
                 s(i,k)= s(i,k)/2.d0/ail*(2*j-1)
              end do
           end do
        end do

        do i=1,6
           do j=i,6
              do k=1,i-1
                 s(i,j)=s(i,j)-s(k,i)*s(k,j)
              end do
              if(i.eq.j) then
                 s(i,i)=sqrt(s(i,i))
              else
                 s(i,j)=s(i,j)/s(i,i)
              end if
           end do
        end do

        aps(6,l,np)=-aps(6,l,np)/s(6,6)
        do i=5,1,-1
           aps(i,l,np)=-aps(i,l,np)
           do k=i+1,6
              aps(i,l,np)=aps(i,l,np)-aps(k,l,np)*s(i,k)
           end do
           aps(i,l,np)=aps(i,l,np)/s(i,i)
        end do
     end do
  end do

  return
end subroutine bachel

integer function convert_dft (dft)
  character(len=*) dft
  integer i, index, ichar, len
  character(len=1) char
!
! convert to lowercase
!
  do i=1,len(trim(dft))
     index = ichar(dft(i:i))
     if (index.ge.65 .and. index.le.90) dft(i:i) = char(index+32)
  end do
!
  if (trim(dft).eq.'pz'.or.trim(dft).eq."'pz'") then
     convert_dft=0
  else if (trim(dft).eq.'blyp'.or.trim(dft).eq."'blyp'") then
     convert_dft=1
  else if (trim(dft).eq.'b88'.or.trim(dft).eq."'b88'") then
     convert_dft=2
  else if (trim(dft).eq.'bp'.or.trim(dft).eq."'bp'") then
     convert_dft=3
  else if (trim(dft).eq.'pw91'.or.trim(dft).eq."'pw91'") then
     convert_dft=4
  else if (trim(dft).eq.'pbe'.or.trim(dft).eq."'pbe'") then
     convert_dft=5
  else
     print *, '***dft unknown: ',dft,'***'
     print *, '***exfact is set to bogus value -9, change by hand!***'
     convert_dft=-9
  end if
!
  return
end function convert_dft
