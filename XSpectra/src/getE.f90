  ! initializes seuilK
  ! Compile par Gwyn Williams,
  ! http://xray.uu.se/hypertext/EBindEnergies.html

function getE(symbol,edge)
  use kinds, ONLY : dp
  implicit none
  real(dp) :: getE 
  character(len=2) :: symbol, edge

  ! a table with all the energies of K1, L1, L2, L3  edges of the atoms in eV

  type seuil
     character(len=2) :: name
     real(dp) :: seuil_en
  end type seuil

  integer, parameter :: Size_tab = 92

  type(seuil), parameter :: seuilK1_tab(Size_tab) = (/ &
      seuilK1('H',    13.6), seuilK1('He',    24.6),seuilK1('Li',    54.7),&
      seuilK1('Be',   111.5),seuilK1('B',   188.0),&!5
      seuilK1('C',   284.2), seuilK1('N',   409.9),seuilK1('O',   543.1),&
      seuilK1('F',   696.7), seuilK1('Ne',   870.2),&! 10
      seuilK1('Na',  1070.8),seuilK1('Mg',  1303.0),seuilK1('Al',  1559.0),&
      seuilK1('Si',  1839.0),seuilK1('P',  2145.5),&! 15
      seuilK1('S',  2472.0), seuilK1('Cl',  2822.4),seuilK1('Ar',  3205.9),&
      seuilK1('K',  3608.4), seuilK1('Ca',  4038.5),&! 20
      seuilK1('Sc',  4492.0),seuilK1('Ti',  4966.0),seuilK1('V',  5465.0),&
      seuilK1('Cr',  5989.0),seuilK1('Mn',  6539.0),&! 25
      seuilK1('Fe',  7112.0),seuilK1('Co',  7709.0),seuilK1('Ni',  8333.0),&
      seuilK1('Cu',  8979.0),seuilK1('Zn',  9659.0),&! 30
      seuilK1('Ga', 10367.0),seuilK1('Ge', 11103.0),seuilK1('As', 11867.0),&
      seuilK1('Se', 12658.0),seuilK1('Br', 13474.0),&! 35
      seuilK1('Kr', 14326.0),seuilK1('Rb', 15200.0),seuilK1('Sr', 16105.0),&
      seuilK1('Y', 17038.0), seuilK1('Zr', 17998.0),&! 40
      seuilK1('Nb', 18986.0),seuilK1('Mo', 20000.0),seuilK1('Tc', 21044.0),&
      seuilK1('Ru', 22117.0),seuilK1('Rh', 23220.0),&! 45
      seuilK1('Pd', 24350.0),seuilK1('Ag', 25514.0),seuilK1('Cd', 26711.0),&
      seuilK1('In', 27940.0),seuilK1('Sn', 29200.0),&! 50
      seuilK1('Sb', 30491.0),seuilK1('Te', 31814.0),seuilK1('I', 33169.0),&
      seuilK1('Xe', 34561.0),seuilK1('Cs', 35985.0),&! 55
      seuilK1('Ba', 37441.0),seuilK1('La', 38925.0),seuilK1('Ce', 40443.0),&
      seuilK1('Pr', 41991.0),seuilK1('Nd', 43569.0),&! 60
      seuilK1('Pm', 45184.0),seuilK1('Sm', 46834.0),seuilK1('Eu', 48519.0),&
      seuilK1('Gd', 50239.0),seuilK1('Tb', 51996.0),&! 65
      seuilK1('Dy', 53789.0),seuilK1('Ho', 55618.0),seuilK1('Er', 57486.0),&
      seuilK1('Tm', 59390.0),seuilK1('Yb', 61332.0),&! 70
      seuilK1('Lu', 63314.0),seuilK1('Hf', 65351.0),seuilK1('Ta', 67416.0),&
      seuilK1('W', 69525.0), seuilK1('Re', 71676.0),&! 75
      seuilK1('Os', 73871.0),seuilK1('Ir', 76111.0),seuilK1('Pt', 78395.0),&
      seuilK1('Au', 80725.0),seuilK1('Hg', 83102.0),&! 80
      seuilK1('Tl', 85530.0),seuilK1('Pb', 88005.0),seuilK1('Bi', 90526.0),&
      seuilK1('Po', 93105.0),seuilK1('At', 95730.0),&! 85
      seuilK1('Rn', 98404.0),seuilK1('Fr',101137.0),seuilK1('Ra',103922.0),&
      seuilK1('Ac',106755.0),seuilK1('Th',109651.0),&! 90
      seuilK1('Pa',112601.0),seuilK1('U',115606.0) /)


  type(seuil), parameter :: seuilL1_tab(Size_tab) = (/ &
      seuilL1('H',    0.),& 
      seuilL1('He',   0.),&
      seuilL1('Li',   0.),&
      seuilL1('Be',   0.),&
      seuilL1('B',    0.),&!5
      seuilL1('C',    0.),& 
      seuilL1('N',  37.3),&
      seuilL1('O',  41.6),&
      seuilL1('F',    0.),&
      seuilL1('Ne', 48.5),&! 10
      seuilL1('Na', 63.5),&
      seuilL1('Mg', 88.6),&
      seuilL1('Al',117.8),&
      seuilL1('Si',149.7),&
      seuilL1('P', 189.0),&! 15
      seuilL1('S', 230.9),& 
      seuilL1('Cl',270.0),&
      seuilL1('Ar',326.3),&
      seuilL1('K', 378.6),& 
      seuilL1('Ca',438.4),&! 20
      seuilL1('Sc',498.0),&
      seuilL1('Ti',560.9),&
      seuilL1('V', 626.7),&
      seuilL1('Cr',696.0),&
      seuilL1('Mn',769.1),&! 25
      seuilL1('Fe',844.6),&
      seuilL1('Co',925.1),&
      seuilL1('Ni',1008.6),&
      seuilL1('Cu',1096.7),&
      seuilL1('Zn',1196.2),&! 30
      seuilL1('Ga', 1299.0),&
      seuilL1('Ge', 1414.6),&
      seuilL1('As', 1527.0),&
      seuilL1('Se', 1652.0),&
      seuilL1('Br', 1782.0),&! 35
      seuilL1('Kr', 1921.0),&
      seuilL1('Rb', 2065.0),&
      seuilL1('Sr', 2216.0),&
      seuilL1('Y',  2373.0),&
      seuilL1('Zr', 2532.0),&! 40
      seuilL1('Nb', 2698.0),&
      seuilL1('Mo', 2866.0),&
      seuilL1('Tc', 3043.0),&
      seuilL1('Ru', 3224.0),&
      seuilL1('Rh', 3412.0),&! 45
      seuilL1('Pd', 3604.0),&
      seuilL1('Ag', 3806.0),&
      seuilL1('Cd', 4018.0),&
      seuilL1('In', 4238.0),&
      seuilL1('Sn', 4465.0),&! 50
      seuilL1('Sb', 4698.0),&
      seuilL1('Te', 4939.0),&
      seuilL1('I',  5188.0),&
      seuilL1('Xe', 5453.0),&
      seuilL1('Cs', 5714.0),&! 55
      seuilL1('Ba', 5989.0),&
      seuilL1('La', 6266.0),&
      seuilL1('Ce', 6548.0),&
      seuilL1('Pr', 6835.0),&
      seuilL1('Nd', 7126.0),&! 60
      seuilL1('Pm', 7428.0),&
      seuilL1('Sm', 7737.0),&
      seuilL1('Eu', 8052.0),&
      seuilL1('Gd', 8376.0),&
      seuilL1('Tb', 8708.0),&! 65
      seuilL1('Dy', 9046.0),&
      seuilL1('Ho', 9394.0),&
      seuilL1('Er', 9751.0),&
      seuilL1('Tm', 10116.0),&
      seuilL1('Yb', 10486.0),&! 70
      seuilL1('Lu', 10870.0),&
      seuilL1('Hf', 11271.0),&
      seuilL1('Ta', 11682.0),&
      seuilL1('W',  12100.0),&
      seuilL1('Re', 12527.0),&! 75
      seuilL1('Os', 12968.0),&
      seuilL1('Ir', 13419.0),&
      seuilL1('Pt', 13880.0),&
      seuilL1('Au', 14353.0),&
      seuilL1('Hg', 14839.0),&! 80
      seuilL1('Tl', 15347.0),&
      seuilL1('Pb', 15861.0),&
      seuilL1('Bi', 16388.0),&
      seuilL1('Po', 16939.0),&
      seuilL1('At', 17493.0),&! 85
      seuilL1('Rn', 18049.0),&
      seuilL1('Fr', 18639.0),&
      seuilL1('Ra', 19237.0),&
      seuilL1('Ac', 19840.0),&
      seuilL1('Th', 20472.0),&! 90
      seuilL1('Pa', 21105.0),&
      seuilL1('U' , 21757.0) /)

  type(seuil), parameter :: seuilL2_tab(Size_tab) = (/ &
      seuilL2('H',  0.0),& 
      seuilL2('He', 0.0),&
      seuilL2('Li', 0.0),&
      seuilL2('Be', 0.0),&
      seuilL2('B',  0.0),&!5
      seuilL2('C',  0.0),& 
      seuilL2('N',  0.0),&
      seuilL2('O',  0.0),&
      seuilL2('F', 48.5),&
      seuilL2('Ne',63.5),&! 10
      seuilL2('Na',88.6),&
      seuilL2('Mg',17.8),&
      seuilL2('Al',49.7),&
      seuilL2('Si',89.0),&
      seuilL2('P', 136.0),&! 15
      seuilL2('S', 163.6),& 
      seuilL2('Cl',202.0),&
      seuilL2('Ar',250.6),&
      seuilL2('K', 297.3),&
      seuilL2('Ca',349.7),&! 20
      seuilL2('Sc',403.6),&
      seuilL2('Ti',460.2),&
      seuilL2('V', 519.8),&
      seuilL2('Cr',583.8),&
      seuilL2('Mn',649.9),&! 25
      seuilL2('Fe',719.9),&
      seuilL2('Co',793.2),&
      seuilL2('Ni',870.0),&
      seuilL2('Cu',952.3),&
      seuilL2('Zn', 1044.9),&! 30
      seuilL2('Ga', 1143.2),&
      seuilL2('Ge', 1248.1),&
      seuilL2('As', 1359.1),&
      seuilL2('Se', 1474.3),&
      seuilL2('Br', 1596.0),&! 35
      seuilL2('Kr', 1730.9),&
      seuilL2('Rb', 1864.0),&
      seuilL2('Sr', 2007.0),&
      seuilL2('Y',  2156.0),&
      seuilL2('Zr', 2307.0),&! 40
      seuilL2('Nb', 2465.0),&
      seuilL2('Mo', 2625.0),&
      seuilL2('Tc', 2793.0),&
      seuilL2('Ru', 2967.0),&
      seuilL2('Rh', 3146.0),&! 45
      seuilL2('Pd', 3330.0),&
      seuilL2('Ag', 3524.0),&
      seuilL2('Cd', 3727.0),&
      seuilL2('In', 3938.0),&
      seuilL2('Sn', 4156.0),&! 50
      seuilL2('Sb', 4380.0),&
      seuilL2('Te', 4612.0),&
      seuilL2('I',  4852.0),&
      seuilL2('Xe', 5107.0),&
      seuilL2('Cs', 5359.0),&! 55
      seuilL2('Ba', 5624.0),&
      seuilL2('La', 5891.0),&
      seuilL2('Ce', 6164.0),&
      seuilL2('Pr', 6440.0),&
      seuilL2('Nd', 6722.0),&! 60
      seuilL2('Pm', 7013.0),&
      seuilL2('Sm', 7312.0),&
      seuilL2('Eu', 7617.0),&
      seuilL2('Gd', 7930.0),&
      seuilL2('Tb', 8252.0),&! 65
      seuilL2('Dy', 8581.0),&
      seuilL2('Ho', 8918.0),&
      seuilL2('Er', 9264.0),&
      seuilL2('Tm', 9617.0),&
      seuilL2('Yb', 9978.0),&! 70
      seuilL2('Lu', 10349.0),&
      seuilL2('Hf', 10739.0),&
      seuilL2('Ta', 11136.0),&
      seuilL2('W',  11544.0),&
      seuilL2('Re', 11959.0),&! 75
      seuilL2('Os', 12385.0),&
      seuilL2('Ir', 12824.0),&
      seuilL2('Pt', 13273.0),&
      seuilL2('Au', 13734.0),&
      seuilL2('Hg', 14209.0),&! 80
      seuilL2('Tl', 14698.0),&
      seuilL2('Pb', 15200.0),&
      seuilL2('Bi', 15711.0),&
      seuilL2('Po', 16244.0),&
      seuilL2('At', 16785.0),&! 85
      seuilL2('Rn', 17337.0),&
      seuilL2('Fr', 17907.0),&
      seuilL2('Ra', 18484.0),&
      seuilL2('Ac', 19083.0),&
      seuilL2('Th', 19693.0),&! 90
      seuilL2('Pa', 20314.0),&
      seuilL2('U' , 20948.0) /)

  type(seuil), parameter :: seuilL3_tab(Size_tab) = (/ &
      seuilL3('H', 0.0),& 
      seuilL3('He',0.0),&
      seuilL3('Li',0.0),&
      seuilL3('Be',0.0),&
      seuilL3('B', 0.0),&!5
      seuilL3('C', 0.0),& 
      seuilL3('N', 0.0),&
      seuilL3('O', 0.0),&
      seuilL3('F', 0.0),&
      seuilL3('Ne',21.6),&! 10
      seuilL3('Na',30.5),&
      seuilL3('Mg',49.2),&
      seuilL3('Al',72.5),&
      seuilL3('Si',99.2),&
      seuilL3('P', 135.0),&! 15
      seuilL3('S', 162.5),& 
      seuilL3('Cl',200.0),&
      seuilL3('Ar',248.4),&
      seuilL3('K', 294.6),&
      seuilL3('Ca',346.2),&! 20
      seuilL3('Sc',398.7),&
      seuilL3('Ti',453.8),&
      seuilL3('V', 512.1),&
      seuilL3('Cr',574.1),&
      seuilL3('Mn',638.7),&! 25
      seuilL3('Fe',706.8),&
      seuilL3('Co',778.1),&
      seuilL3('Ni',852.7),&
      seuilL3('Cu',932.7),&
      seuilL3('Zn', 1021.8),&! 30
      seuilL3('Ga', 1116.4),&
      seuilL3('Ge', 1217.0),&
      seuilL3('As', 1323.6),&
      seuilL3('Se', 1433.9),&
      seuilL3('Br', 1550.0),&! 35
      seuilL3('Kr', 1678.4),&
      seuilL3('Rb', 1804.0),&
      seuilL3('Sr', 1940.0),&
      seuilL3('Y',  2080.0),&
      seuilL3('Zr', 2223.0),&! 40
      seuilL3('Nb', 2371.0),&
      seuilL3('Mo', 2520.0),&
      seuilL3('Tc', 2677.0),&
      seuilL3('Ru', 2838.0),&
      seuilL3('Rh', 3004.0),&! 45
      seuilL3('Pd', 3173.0),&
      seuilL3('Ag', 3351.0),&
      seuilL3('Cd', 3538.0),&
      seuilL3('In', 3730.0),&
      seuilL3('Sn', 3929.0),&! 50
      seuilL3('Sb', 4132.0),&
      seuilL3('Te', 4341.0),&
      seuilL3('I',  4557.0),&
      seuilL3('Xe', 4786.0),&
      seuilL3('Cs', 5012.0),&! 55
      seuilL3('Ba', 5247.0),&
      seuilL3('La', 5483.0),&
      seuilL3('Ce', 5723.0),&
      seuilL3('Pr', 5964.0),&
      seuilL3('Nd', 6208.0),&! 60
      seuilL3('Pm', 6459.0),&
      seuilL3('Sm', 6716.0),&
      seuilL3('Eu', 6977.0),&
      seuilL3('Gd', 7243.0),&
      seuilL3('Tb', 7510.0),&! 65
      seuilL3('Dy', 7790.0),&
      seuilL3('Ho', 8071.0),&
      seuilL3('Er', 8358.0),&
      seuilL3('Tm', 8648.0),&
      seuilL3('Yb', 8944.0),&! 70
      seuilL3('Lu', 9244.0),&
      seuilL3('Hf', 9561.0),&
      seuilL3('Ta', 9881.0),&
      seuilL3('W',  10207.0),&
      seuilL3('Re', 10535.0),&! 75
      seuilL3('Os', 10871.0),&
      seuilL3('Ir', 11215.0),&
      seuilL3('Pt', 11564.0),&
      seuilL3('Au', 11919.0),&
      seuilL3('Hg', 12284.0),&! 80
      seuilL3('Tl', 12658.0),&
      seuilL3('Pb', 13035.0),&
      seuilL3('Bi', 13419.0),&
      seuilL3('Po', 13814.0),&
      seuilL3('At', 14214.0),&! 85
      seuilL3('Rn', 14619.0),&
      seuilL3('Fr', 15031.0),&
      seuilL3('Ra', 15444.0),&
      seuilL3('Ac', 15871.0),&
      seuilL3('Th', 16300.0),&! 90
      seuilL3('Pa', 16733.0),&
      seuilL3('U' , 17166.0) /)

  integer :: i

  do i = 1, Size_tab
     if (symbol.eq.seuilK1_tab(i)%name .and. edge.eq.'K1' ) then
        getE = seuilK1_tab(i)%seuil_en
        return
     else if (symbol.eq.seuilL1_tab(i)%name .and. edge.eq.'L1' )
        getE = seuilL1_tab(i)%seuil_en
        return 
     else if (symbol.eq.seuilL2_tab(i)%name .and. edge.eq.'L2' )
        getE = seuilL2_tab(i)%seuil_en
        return 
     else if (symbol.eq.seuilL3_tab(i)%name .and. edge.eq.'L3' )
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
