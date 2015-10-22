module coef_CG

  use kinds, only : dp
contains
!    ARTURO QUIRANTES SIERRA
!    Department of Applied Physics, Faculty of Sciences
!    University of Granada, 18071 Granada (SPAIN)
!    http://www.ugr.es/local/aquiran/codigos.htm
!    aquiran@ugr.es
!
!    Last update: 20 May 2.003    
!
!    Subroutine NED
!    to calculate Clebsch-Gordan coefficients
!
!    You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine
!    Input:
!        AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices)
!    Output:
!        CG=C-G(AJ,BJ,CJ,AM,BM,CM)
!
    DOUBLE PRECISION function ClebschG(AJ,BJ,CJ,AM,BM,CM)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION Q(100,100)
    DOUBLE PRECISION CG
    DOUBLE PRECISION AJ,BJ,CJ,AM,BM,CM
    INTEGER ZZ
    ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
    DO 2 I=1,ZZ
        Q(I,1)=1.D0
        Q(I,I)=1.D0
2    CONTINUE
    DO 3 I=2,ZZ-1
    DO 3 K=2,I
        Q(I+1,K)=Q(I,K-1)+Q(I,K)
3    CONTINUE
    CG=0.D0
    JA=AJ+AM+1.01D0
    MA=AJ-AM+1.01D0
    JB=BJ+BM+1.01D0
    MB=BJ-BM+1.01D0
    JC=CJ+CM+1.01D0
    MC=CJ-CM+1.01D0
    LA=BJ+CJ-AJ+1.01D0
    LB=CJ+AJ-BJ+1.01D0
    LC=AJ+BJ-CJ+1.01D0
    LT=AJ+BJ+CJ+1.01D0
    D=DABS(AM+BM-CM)-0.01D0
    IF (D) 10,10,20
10    LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
    IF (LD) 20,20,30
30    JA2=AJ+AJ+AM+AM
    JB2=BJ+BJ+BM+BM
    JC2=CJ+CJ-CM-CM
    I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
    IF (I2) 20,40,20
40    FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
    FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
    FN=FN/Q(JA+MA-1,JA)
    FN=FN/Q(JB+MB-1,JB)
    FN=FN/Q(JC+MC-1,JC)
    K0=MAX(0,LC-JA,LC-MB)+1
    K1=MIN(LC,MA,JB)
    X=0.D0
    DO 50 K=K0,K1
        X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
50    CONTINUE
    IP=K1+LB+JC
    P=1-2*(IP-IP/2*2)
    CG=P*X*DSQRT(FN)
!    What we´ve calculated is a Wigner 3-j coefficient
!    Next, we´ll turn it into a Clebsch-Gordan coefficient
    CG=CG*DSQRT(2*CJ+1)*(-1)**IDNINT(AJ-BJ-CM)
20    CONTINUE
    ClebschG = CG
    RETURN
    END function ClebschG

end module coef_CG


module coef_gaunt

 
  use kinds, only : dp
  USE coef_CG
 
contains
 
! Clebsch-gordan coefficient  eq. 3.18, Rose

      real(kind=DP) function cgc(l1,l2,l3,m1,m2)

      implicit real(kind=dp) (a-h,o-z)

      integer, dimension(5):: nd, num

      cgc = 0.d0

      m3 = m1 + m2

! Arguments of factorials
      num( 1 ) = l3 + l1 - l2
      num( 2 ) = l3 - l1 + l2
      num( 3 ) = l1 + l2 - l3
      num( 4 ) = l3 + m3
      num( 5 ) = l3 - m3
      nd( 1 ) = l1 + l2 + l3 + 1
      nd( 2 ) = l1 - m1
      nd( 3 ) = l1 + m1
      nd( 4 ) = l2 - m2
      nd( 5 ) = l2 + m2

! Check triangle and projection conditions
      do i = 1, 5
        if( num(i) < 0 .or. nd(i) < 0 ) return
      end do

      nff = 0
      nmin = max( 0, l2 + m3 - l1 )
      ff = 1.0d0

! Two sets of factorial products
      n = 5
      do nfac = 1, 2
        n1 = n - 1
! Arrange arguments in descending order
        do i = 1, n1
          inum = i
          id = i
          i1 = i + 1
          do j = i1, n
            if( num( j ) > num( inum ) ) inum = j
            if( nd( j ) > nd( id ) ) id = j
          end do
          ntemp = num( i )
          num( i ) = num( inum )
          num( inum ) = ntemp
          ntemp = nd( i )
          nd( i ) = nd( id )
          nd( id ) = ntemp
        end do
! Compute factorial ratios
        do i = 1, n
          if( num( i ) < nd( i ) ) then
            jm = nd( i )
            if( jm == 1 ) cycle
            j0 = num( i ) + 1
            if( num( i ) == 0 ) j0 = 2
            do j = j0, jm
              if( abs ( ff ) < 1.0e-20 ) then
                ff = ff * 1.e20
                nff = nff - 2
              endif
              ff = ff / j
            end do
          elseif( num( i ) > nd( i ) ) then
            jm = num( i )
            if( jm == 1 ) cycle
            j0 = nd( i ) + 1
            if( nd( i ) == 0 ) j0 = 2
            do j = j0, jm
              if( abs ( ff ) > 1.0e20 ) then
                ff = ff / 1.0e20
                nff = nff + 2
              endif
              ff = j * ff
            end do
          endif
        end do

        if ( nfac == 2 ) exit

        nff = nff / 2
        ff = sqrt( ( 2 * l3 + 1 ) * ff )

! Second set of factorial arguments
        num( 1 ) = l2 + l3 + m1 - nmin
        num( 2 ) = l1 - m1 + nmin
        num( 3 ) = 0
        nd( 1 ) = nmin
        if( nmin == 0 ) nd( 1 ) = l1 - l2 - m3
        nd( 2 ) = l3 - l1 + l2 - nmin
        nd( 3 ) = l3 + m3 - nmin
        n = 3
      end do

      if( mod( nmin + l2 + m2, 2 ) /= 0 ) ff = - ff
      ff = ff * 1.0e10**nff
      cgcp = ff
      nmax = min( l3 - l1 + l2, l3 + m3 )

      do nu = nmin+1, nmax
        ff= - ff * ( (l1-m1+nu) * (l3-l1+l2-nu+1) * (l3+m3-nu+1) ) &
           / dble( nu * (nu+l1-l2-m3) * (l2+l3+m1-nu+1) )
        cgcp = cgcp + ff
      end do

      cgc = cgcp

   return
  end function cgc

  real(kind=DP) function gaunt(l1,m1,l2,m2,l3,m3)

! adapted from FDMNES (Y. Joly)
! Calculates the Gaunt coef for cplx harmonics: Y1*Y2Y3

  implicit none
! external cgc
  integer l1,m1,l2,m2,l3,m3
  real(kind=DP) fac,dl1,dm1,dl2,dm2,dl3,dm3
  real(kind=DP) :: pi = 3.14159265
  if( m1 /= m2 + m3 .or. mod(l1+l2+l3,2) == 1 &
      .or. l1 < abs(l2-l3) .or. l1 > abs(l2+l3) ) then

    gaunt = 0.d0

  else

    fac = sqrt( (2*l3 + 1.d0) * (2*l2 + 1.d0)    &
            / ( (2*l1 + 1.d0) * 4 * pi ) )
    dl1 = dble( l1 )
    dl2 = dble( l2 )
    dl3 = dble( l3 )
    dm2 = dble( m2 )
    dm3 = dble( m3 )

!    gaunt = fac * ClebschG(dl3,dl2,dl1,dm3,dm2,dm3-dm2) &
!                * ClebschG(dl3,dl2,dl1,0.d0,0.d0,0.d0)
    gaunt = fac * cgc(l3,l2,l1,m3,m2) * cgc(l3,l2,l1,0,0)
  endif

  return
  end function gaunt

  real(kind=DP) function gaunt4Y(l1,m1,l2,m2,l3,m3,l4,m4)

! calculates \int d\Omega Y(l1,m1)* Y(l2,m2) Y(l3,m3)* Y(l4,m4) 
! uses : Y(L1)* Y(L2) = \sum_L  gaunt(L2,L1,L) Y(L)  where L = (l,m)
! hence: gaunt4Y(L1,L2,L3,L4) = \sum_L   G3(L2,L1,L) x G3(L4,L3,L)
! derivation in OB PhD thesis (appendix)

  implicit none
! external gaunt
  integer, intent(in):: l1, m1, l2, m2, l3, m3, l4, m4
  integer:: l, lmin, lmax, m


  gaunt4Y = 0.d0

  lmin = max( abs(l2-l1), abs(l3-l4) )
  lmax = min( l2+l1, l3+l4 )

  do l = lmin,lmax,2
    do m = -l,l
      gaunt4Y = gaunt4Y + gaunt(l2,m2,l1,m1,l,m)     &
                        * gaunt(l4,m4,l3,m3,l,m)
    end do
  end do

  return
  end function gaunt4Y



end module coef_gaunt


  ! initializes seuilK
  ! Compile par Gwyn Williams,
  ! http://xray.uu.se/hypertext/EBindEnergies.html

module edge_energy

  use kinds, ONLY : dp
  implicit none

contains
  real(dp) function getE(symbol,edge)
 
   character(len=2) :: symbol, edge
  
    ! a table with all the energies of K1, L1, L2, L3  edges of the atoms in eV
  
    type seuil
       character(len=2) :: name
       real(dp) :: seuil_en
    end type seuil
  
    integer, parameter :: Size_tab = 92
  
    type(seuil), parameter :: seuilK1_tab(Size_tab) = (/ &
        seuil('H',    13.6), seuil('He',    24.6),seuil('Li',    54.7),&
        seuil('Be',   111.5),seuil('B',   188.0),&!5
        seuil('C',   284.2), seuil('N',   409.9),seuil('O',   543.1),&
        seuil('F',   696.7), seuil('Ne',   870.2),&! 10
        seuil('Na',  1070.8),seuil('Mg',  1303.0),seuil('Al',  1559.0),&
        seuil('Si',  1839.0),seuil('P',  2145.5),&! 15
        seuil('S',  2472.0), seuil('Cl',  2822.4),seuil('Ar',  3205.9),&
        seuil('K',  3608.4), seuil('Ca',  4038.5),&! 20
        seuil('Sc',  4492.0),seuil('Ti',  4966.0),seuil('V',  5465.0),&
        seuil('Cr',  5989.0),seuil('Mn',  6539.0),&! 25
        seuil('Fe',  7112.0),seuil('Co',  7709.0),seuil('Ni',  8333.0),&
        seuil('Cu',  8979.0),seuil('Zn',  9659.0),&! 30
        seuil('Ga', 10367.0),seuil('Ge', 11103.0),seuil('As', 11867.0),&
        seuil('Se', 12658.0),seuil('Br', 13474.0),&! 35
        seuil('Kr', 14326.0),seuil('Rb', 15200.0),seuil('Sr', 16105.0),&
        seuil('Y', 17038.0), seuil('Zr', 17998.0),&! 40
        seuil('Nb', 18986.0),seuil('Mo', 20000.0),seuil('Tc', 21044.0),&
        seuil('Ru', 22117.0),seuil('Rh', 23220.0),&! 45
        seuil('Pd', 24350.0),seuil('Ag', 25514.0),seuil('Cd', 26711.0),&
        seuil('In', 27940.0),seuil('Sn', 29200.0),&! 50
        seuil('Sb', 30491.0),seuil('Te', 31814.0),seuil('I', 33169.0),&
        seuil('Xe', 34561.0),seuil('Cs', 35985.0),&! 55
        seuil('Ba', 37441.0),seuil('La', 38925.0),seuil('Ce', 40443.0),&
        seuil('Pr', 41991.0),seuil('Nd', 43569.0),&! 60
        seuil('Pm', 45184.0),seuil('Sm', 46834.0),seuil('Eu', 48519.0),&
        seuil('Gd', 50239.0),seuil('Tb', 51996.0),&! 65
        seuil('Dy', 53789.0),seuil('Ho', 55618.0),seuil('Er', 57486.0),&
        seuil('Tm', 59390.0),seuil('Yb', 61332.0),&! 70
        seuil('Lu', 63314.0),seuil('Hf', 65351.0),seuil('Ta', 67416.0),&
        seuil('W', 69525.0), seuil('Re', 71676.0),&! 75
        seuil('Os', 73871.0),seuil('Ir', 76111.0),seuil('Pt', 78395.0),&
        seuil('Au', 80725.0),seuil('Hg', 83102.0),&! 80
        seuil('Tl', 85530.0),seuil('Pb', 88005.0),seuil('Bi', 90526.0),&
        seuil('Po', 93105.0),seuil('At', 95730.0),&! 85
        seuil('Rn', 98404.0),seuil('Fr',101137.0),seuil('Ra',103922.0),&
        seuil('Ac',106755.0),seuil('Th',109651.0),&! 90
        seuil('Pa',112601.0),seuil('U',115606.0) /)
  
  
    type(seuil), parameter :: seuilL1_tab(Size_tab) = (/ &
        seuil('H',    0.),& 
        seuil('He',   0.),&
        seuil('Li',   0.),&
        seuil('Be',   0.),&
        seuil('B',    0.),&!5
        seuil('C',    0.),& 
        seuil('N',  37.3),&
        seuil('O',  41.6),&
        seuil('F',    0.),&
        seuil('Ne', 48.5),&! 10
        seuil('Na', 63.5),&
        seuil('Mg', 88.6),&
        seuil('Al',117.8),&
        seuil('Si',149.7),&
        seuil('P', 189.0),&! 15
        seuil('S', 230.9),& 
        seuil('Cl',270.0),&
        seuil('Ar',326.3),&
        seuil('K', 378.6),& 
        seuil('Ca',438.4),&! 20
        seuil('Sc',498.0),&
        seuil('Ti',560.9),&
        seuil('V', 626.7),&
        seuil('Cr',696.0),&
        seuil('Mn',769.1),&! 25
        seuil('Fe',844.6),&
        seuil('Co',925.1),&
        seuil('Ni',1008.6),&
        seuil('Cu',1096.7),&
        seuil('Zn',1196.2),&! 30
        seuil('Ga', 1299.0),&
        seuil('Ge', 1414.6),&
        seuil('As', 1527.0),&
        seuil('Se', 1652.0),&
        seuil('Br', 1782.0),&! 35
        seuil('Kr', 1921.0),&
        seuil('Rb', 2065.0),&
        seuil('Sr', 2216.0),&
        seuil('Y',  2373.0),&
        seuil('Zr', 2532.0),&! 40
        seuil('Nb', 2698.0),&
        seuil('Mo', 2866.0),&
        seuil('Tc', 3043.0),&
        seuil('Ru', 3224.0),&
        seuil('Rh', 3412.0),&! 45
        seuil('Pd', 3604.0),&
        seuil('Ag', 3806.0),&
        seuil('Cd', 4018.0),&
        seuil('In', 4238.0),&
        seuil('Sn', 4465.0),&! 50
        seuil('Sb', 4698.0),&
        seuil('Te', 4939.0),&
        seuil('I',  5188.0),&
        seuil('Xe', 5453.0),&
        seuil('Cs', 5714.0),&! 55
        seuil('Ba', 5989.0),&
        seuil('La', 6266.0),&
        seuil('Ce', 6548.0),&
        seuil('Pr', 6835.0),&
        seuil('Nd', 7126.0),&! 60
        seuil('Pm', 7428.0),&
        seuil('Sm', 7737.0),&
        seuil('Eu', 8052.0),&
        seuil('Gd', 8376.0),&
        seuil('Tb', 8708.0),&! 65
        seuil('Dy', 9046.0),&
        seuil('Ho', 9394.0),&
        seuil('Er', 9751.0),&
        seuil('Tm', 10116.0),&
        seuil('Yb', 10486.0),&! 70
        seuil('Lu', 10870.0),&
        seuil('Hf', 11271.0),&
        seuil('Ta', 11682.0),&
        seuil('W',  12100.0),&
        seuil('Re', 12527.0),&! 75
        seuil('Os', 12968.0),&
        seuil('Ir', 13419.0),&
        seuil('Pt', 13880.0),&
        seuil('Au', 14353.0),&
        seuil('Hg', 14839.0),&! 80
        seuil('Tl', 15347.0),&
        seuil('Pb', 15861.0),&
        seuil('Bi', 16388.0),&
        seuil('Po', 16939.0),&
        seuil('At', 17493.0),&! 85
        seuil('Rn', 18049.0),&
        seuil('Fr', 18639.0),&
        seuil('Ra', 19237.0),&
        seuil('Ac', 19840.0),&
        seuil('Th', 20472.0),&! 90
        seuil('Pa', 21105.0),&
        seuil('U' , 21757.0) /)
  
    type(seuil), parameter :: seuilL2_tab(Size_tab) = (/ &
        seuil('H',  0.0),& 
        seuil('He', 0.0),&
        seuil('Li', 0.0),&
        seuil('Be', 0.0),&
        seuil('B',  0.0),&!5
        seuil('C',  0.0),& 
        seuil('N',  0.0),&
        seuil('O',  0.0),&
        seuil('F', 48.5),&
        seuil('Ne',63.5),&! 10
        seuil('Na',88.6),&
        seuil('Mg',17.8),&
        seuil('Al',49.7),&
        seuil('Si',89.0),&
        seuil('P', 136.0),&! 15
        seuil('S', 163.6),& 
        seuil('Cl',202.0),&
        seuil('Ar',250.6),&
        seuil('K', 297.3),&
        seuil('Ca',349.7),&! 20
        seuil('Sc',403.6),&
        seuil('Ti',460.2),&
        seuil('V', 519.8),&
        seuil('Cr',583.8),&
        seuil('Mn',649.9),&! 25
        seuil('Fe',719.9),&
        seuil('Co',793.2),&
        seuil('Ni',870.0),&
        seuil('Cu',952.3),&
        seuil('Zn', 1044.9),&! 30
        seuil('Ga', 1143.2),&
        seuil('Ge', 1248.1),&
        seuil('As', 1359.1),&
        seuil('Se', 1474.3),&
        seuil('Br', 1596.0),&! 35
        seuil('Kr', 1730.9),&
        seuil('Rb', 1864.0),&
        seuil('Sr', 2007.0),&
        seuil('Y',  2156.0),&
        seuil('Zr', 2307.0),&! 40
        seuil('Nb', 2465.0),&
        seuil('Mo', 2625.0),&
        seuil('Tc', 2793.0),&
        seuil('Ru', 2967.0),&
        seuil('Rh', 3146.0),&! 45
        seuil('Pd', 3330.0),&
        seuil('Ag', 3524.0),&
        seuil('Cd', 3727.0),&
        seuil('In', 3938.0),&
        seuil('Sn', 4156.0),&! 50
        seuil('Sb', 4380.0),&
        seuil('Te', 4612.0),&
        seuil('I',  4852.0),&
        seuil('Xe', 5107.0),&
        seuil('Cs', 5359.0),&! 55
        seuil('Ba', 5624.0),&
        seuil('La', 5891.0),&
        seuil('Ce', 6164.0),&
        seuil('Pr', 6440.0),&
        seuil('Nd', 6722.0),&! 60
        seuil('Pm', 7013.0),&
        seuil('Sm', 7312.0),&
        seuil('Eu', 7617.0),&
        seuil('Gd', 7930.0),&
        seuil('Tb', 8252.0),&! 65
        seuil('Dy', 8581.0),&
        seuil('Ho', 8918.0),&
        seuil('Er', 9264.0),&
        seuil('Tm', 9617.0),&
        seuil('Yb', 9978.0),&! 70
        seuil('Lu', 10349.0),&
        seuil('Hf', 10739.0),&
        seuil('Ta', 11136.0),&
        seuil('W',  11544.0),&
        seuil('Re', 11959.0),&! 75
        seuil('Os', 12385.0),&
        seuil('Ir', 12824.0),&
        seuil('Pt', 13273.0),&
        seuil('Au', 13734.0),&
        seuil('Hg', 14209.0),&! 80
        seuil('Tl', 14698.0),&
        seuil('Pb', 15200.0),&
        seuil('Bi', 15711.0),&
        seuil('Po', 16244.0),&
        seuil('At', 16785.0),&! 85
        seuil('Rn', 17337.0),&
        seuil('Fr', 17907.0),&
        seuil('Ra', 18484.0),&
        seuil('Ac', 19083.0),&
        seuil('Th', 19693.0),&! 90
        seuil('Pa', 20314.0),&
        seuil('U' , 20948.0) /)
  
    type(seuil), parameter :: seuilL3_tab(Size_tab) = (/ &
        seuil('H', 0.0),& 
        seuil('He',0.0),&
        seuil('Li',0.0),&
        seuil('Be',0.0),&
        seuil('B', 0.0),&!5
        seuil('C', 0.0),& 
        seuil('N', 0.0),&
        seuil('O', 0.0),&
        seuil('F', 0.0),&
        seuil('Ne',21.6),&! 10
        seuil('Na',30.5),&
        seuil('Mg',49.2),&
        seuil('Al',72.5),&
        seuil('Si',99.2),&
        seuil('P', 135.0),&! 15
        seuil('S', 162.5),& 
        seuil('Cl',200.0),&
        seuil('Ar',248.4),&
        seuil('K', 294.6),&
        seuil('Ca',346.2),&! 20
        seuil('Sc',398.7),&
        seuil('Ti',453.8),&
        seuil('V', 512.1),&
        seuil('Cr',574.1),&
        seuil('Mn',638.7),&! 25
        seuil('Fe',706.8),&
        seuil('Co',778.1),&
        seuil('Ni',852.7),&
        seuil('Cu',932.7),&
        seuil('Zn', 1021.8),&! 30
        seuil('Ga', 1116.4),&
        seuil('Ge', 1217.0),&
        seuil('As', 1323.6),&
        seuil('Se', 1433.9),&
        seuil('Br', 1550.0),&! 35
        seuil('Kr', 1678.4),&
        seuil('Rb', 1804.0),&
        seuil('Sr', 1940.0),&
        seuil('Y',  2080.0),&
        seuil('Zr', 2223.0),&! 40
        seuil('Nb', 2371.0),&
        seuil('Mo', 2520.0),&
        seuil('Tc', 2677.0),&
        seuil('Ru', 2838.0),&
        seuil('Rh', 3004.0),&! 45
        seuil('Pd', 3173.0),&
        seuil('Ag', 3351.0),&
        seuil('Cd', 3538.0),&
        seuil('In', 3730.0),&
        seuil('Sn', 3929.0),&! 50
        seuil('Sb', 4132.0),&
        seuil('Te', 4341.0),&
        seuil('I',  4557.0),&
        seuil('Xe', 4786.0),&
        seuil('Cs', 5012.0),&! 55
        seuil('Ba', 5247.0),&
        seuil('La', 5483.0),&
        seuil('Ce', 5723.0),&
        seuil('Pr', 5964.0),&
        seuil('Nd', 6208.0),&! 60
        seuil('Pm', 6459.0),&
        seuil('Sm', 6716.0),&
        seuil('Eu', 6977.0),&
        seuil('Gd', 7243.0),&
        seuil('Tb', 7510.0),&! 65
        seuil('Dy', 7790.0),&
        seuil('Ho', 8071.0),&
        seuil('Er', 8358.0),&
        seuil('Tm', 8648.0),&
        seuil('Yb', 8944.0),&! 70
        seuil('Lu', 9244.0),&
        seuil('Hf', 9561.0),&
        seuil('Ta', 9881.0),&
        seuil('W',  10207.0),&
        seuil('Re', 10535.0),&! 75
        seuil('Os', 10871.0),&
        seuil('Ir', 11215.0),&
        seuil('Pt', 11564.0),&
        seuil('Au', 11919.0),&
        seuil('Hg', 12284.0),&! 80
        seuil('Tl', 12658.0),&
        seuil('Pb', 13035.0),&
        seuil('Bi', 13419.0),&
        seuil('Po', 13814.0),&
        seuil('At', 14214.0),&! 85
        seuil('Rn', 14619.0),&
        seuil('Fr', 15031.0),&
        seuil('Ra', 15444.0),&
        seuil('Ac', 15871.0),&
        seuil('Th', 16300.0),&! 90
        seuil('Pa', 16733.0),&
        seuil('U' , 17166.0) /)
  
   character (len=2) :: sym1
    integer :: i
    
    sym1 = trim(adjustl(symbol))
  
    do i = 1, Size_tab
       if (sym1.eq.trim(adjustl(seuilK1_tab(i)%name)) .and. edge.eq.'K1' ) then
          getE = seuilK1_tab(i)%seuil_en
          return
       else if (sym1.eq.trim(adjustl(seuilL1_tab(i)%name)) .and. edge.eq.'L1' ) then
          getE = seuilL1_tab(i)%seuil_en
          return 
       else if (sym1.eq.trim(adjustl(seuilL2_tab(i)%name)) .and. edge.eq.'L2' ) then
          getE = seuilL2_tab(i)%seuil_en
          return 
       else if (sym1.eq.trim(adjustl(seuilL3_tab(i)%name)) .and. edge.eq.'L3' ) then
          getE = seuilL3_tab(i)%seuil_en
          return 
       end if
    end do
  
    if( edge /= 'K1' .and. edge /= 'L1' .and. edge /= 'L2' .and. edge /= 'L3' ) then 
      write(*,'(3A)') 'For the',edge,'the tables are missing !!'
    else
      write(*,'(3A)') 'Could not find the element ', trim(symbol), &
         ' in the table of edge energies!'
    end if
    stop
  
  end function getE

end module edge_energy
