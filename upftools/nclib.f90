!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine simpson(mesh,func,rab,asum)
  !-----------------------------------------------------------------------
!
  !     simpson's rule integrator for function stored on the
  !     radial logarithmic mesh
  !
  
  implicit none
  
  integer :: i, mesh
  real(kind=8) ::  rab(mesh), func(mesh), f1, f2, f3, r12, asum
      
      !     routine assumes that mesh is an odd number so run check
      !     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
      !       write(*,*) '***error in subroutine radlg'
!       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
!       stop
!     endif

  asum = 0.0d0
  r12 = 1.0d0 / 12.0d0
  f3  = func(1) * rab(1) * r12
  
  do i = 2,mesh-1,2
     f1 = f3
     f2 = func(i) * rab(i) * r12
     f3 = func(i+1) * rab(i+1) * r12
     asum = asum + 4.0d0*f1 + 16.0d0*f2 + 4.0d0*f3
  enddo
  
  return
end subroutine simpson

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
       6.99638348861914d0, -3.56098437018154d-2 /
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
       7.21175825088309d0, 5.64195517478974d-1, -1.36864857382717d-7 /
  data q2 / 3.00459260956983d2, 7.90950925327898d2, &
       9.31354094850610d2, 6.38980264465631d2, 2.77585444743988d2, &
       7.70001529352295d1, 1.27827273196294d1, 1.00000000000000d0 /
  data p3 / -2.99610707703542d-3, -4.94730910623251d-2, &
       -2.26956593539687d-1, -2.78661308609648d-1, -2.23192459734185d-2 &
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
!
!-----------------------------------------------------------------------
subroutine which_dft (dft, iexch, icorr, igcx, igcc)  
  !-----------------------------------------------------------------------
  !
  implicit none  
  ! input
  character (len=*) :: dft  
  ! output
  integer :: iexch, icorr, igcx, igcc  
  ! data
  integer :: nxc, ncc, ngcx, ngcc  
  parameter (nxc = 1, ncc = 9, ngcx = 3, ngcc = 4)  
  character (len=3) :: exc, corr
  character (len=4) :: gradx, gradc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)
  ! local
  integer :: len, l, i, notset  
  character (len=1) :: capital
  character (len=50):: dftout * 50  
  logical :: matches  
  external capital, matches  
  data notset / -1 /  
  data exc / 'NOX', 'SLA' /  
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
       'OBW', 'GL' /
  data gradx / 'NOGX', 'B88', 'GGX', 'PBE' /  


  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBE' /  
  ! convert to uppercase
  len = len_trim(dft)  
  dftout = ' '  
  do l = 1, len  
     dftout (l:l) = capital (dft (l:l) )  
  enddo
  !  exchange

  iexch = notset  
  do i = 0, nxc  
     if (matches (exc (i), dftout) ) call set_dft_value (iexch, i)  
  enddo
  !  correlation

  icorr = notset  
  do i = 0, ncc  
     if (matches (corr (i), dftout) ) call set_dft_value (icorr, i)  
  enddo
  !  gradient correction, exchange

  igcx = notset  
  do i = 0, ngcx  
     if (matches (gradx (i), dftout) ) call set_dft_value (igcx, i)  
  enddo
  !  gradient correction, correlation

  igcc = notset  
  do i = 0, ngcc  
     if (matches (gradc (i), dftout) ) call set_dft_value (igcc, i)  
  enddo
  ! special case : BLYP => B88 for gradient correction on exchange


  if (matches ('BLYP', dftout) ) call set_dft_value (igcx, 1)  
  ! special case : PBE
  if (matches ('PBE', dftout) ) then  
     call set_dft_value (iexch, 1)  
     call set_dft_value (icorr, 4)  
  endif
  ! special case : BP = B88 + P86
  if (matches ('BP', dftout) ) then  
     call set_dft_value (igcx, 1)  
     call set_dft_value (igcc, 1)  
  endif
  ! special case : PW91 = GGX + GGC
  if (matches ('PW91', dftout) ) then  
     call set_dft_value (igcx, 2)  
     call set_dft_value (igcc, 2)  
  endif
  ! Default value: Slater exchange


  if (iexch.eq.notset) call set_dft_value (iexch, 1)  
  ! Default value: Perdew-Zunger correlation

  if (icorr.eq.notset) call set_dft_value (icorr, 1)  
  ! Default value: no gradient correction on exchange

  if (igcx.eq.notset) call set_dft_value (igcx, 0)  
  ! Default value: no gradient correction on correlation

  if (igcc.eq.notset) call set_dft_value (igcc, 0)  
  !

  dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
       &//gradc (igcc)
  !cc      write (6,'(a)') dftout
  return  
end subroutine which_dft
!
!-----------------------------------------------------------------------
subroutine set_dft_value (m, i)  
  !-----------------------------------------------------------------------
  !
  implicit none  
  ! input / output
  integer :: m, i  
  ! local
  integer :: notset  

  parameter (notset = -1)  

  if (m.ne.notset.and.m.ne.i) call errore ('set_dft_value', &
       'two conflicting matching values', 1)

  m = i  
  return  

end subroutine set_dft_value
!-----------------------------------------------------------------------
logical function matches (string1, string2)  
  !-----------------------------------------------------------------------
  !
  implicit none  
  character (len=*) :: string1, string2  
  integer :: len1, len2, l  


  len1 = len_trim(string1)  
  len2 = len_trim(string2)  
  do l = 1, len2 - len1 + 1  
     if (string1 (1:len1) .eq.string2 (l:l + len1 - 1) ) then  
        matches = .true.  
        return  
     endif

  enddo

  matches = .false.  
  return  
end function matches
!
!-----------------------------------------------------------------------
function capital (character)  
  !-----------------------------------------------------------------------
  !
  !   converts character to capital if lowercase
  !   copy character to output in all other cases
  !
  implicit none  
  character (len=1) :: capital, character
  !
  character(len=26) :: minuscole='abcdefghijklmnopqrstuvwxyz', &
                       maiuscole='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer :: i
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
end function capital
!
! ------------------------------------------------------------------
function atom_name(atomic_number)
! ------------------------------------------------------------------
  !
  integer :: atomic_number
  character(len=2) :: atom_name
  !
  character(len=2), dimension (94) :: elements
  data elements/' H',                              'He', &
                'Li','Be',' B',' C',' N',' O',' F','Ne', &  
                'Na','Mg','Al','Si',' P',' S','Cl','Ar', &
                ' K','Ca','Sc','Ti',' V','Cr','Mn',      &
                          'Fe','Co','Ni','Cu','Zn',      &
                          'Ga','Ge','As','Se','Br','Kr', &
                'Rb','Sr',' Y','Zr','Nb','Mo','Tc',      &
                          'Ru','Rh','Pd','Ag','Cd',      &
                          'In','Sn','Sb','Te',' I','Xe', &
                'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', &
                               'Tb','Dy','Ho','Er','Tm','Yb','Lu', &
                          'Hf','Ta',' W','Re','Os', &
                               'Ir','Pt','Au','Hg', &
                          'Tl','Pb','Bi','Po','At','Rn', &
                'Fr','Ra','Ac','Th','Pa',' U','Np','Pu' /
  !
  if (atomic_number.lt.1.or.atomic_number.gt.94) then
     call errore('atom_name','invalid atomic number', &
          1000+atomic_number)
  else
     atom_name=elements(atomic_number)
  end if
  return

end function atom_name
