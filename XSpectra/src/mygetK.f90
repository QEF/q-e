  ! initializes seuilK
  ! Compile par Gwyn Williams,
  ! http://xray.uu.se/hypertext/EBindEnergies.html
function mygetK(sym)
  use kinds, ONLY : dp
  implicit none
  real(dp) :: mygetK
  character(len=2), intent(in) :: sym
  ! seuilK: a table with all the energies of K edges of the atoms in eV
  type seuilK
     character(len=2) :: name
     real(dp) :: Kenergy
  end type seuilK

  integer, parameter :: Size_tab = 92

  type(seuilK), parameter :: seuilK_tab(Size_tab) = (/ &
      seuilK('H',    13.6), seuilK('He',    24.6),seuilK('Li',    54.7),&
      seuilK('Be',   111.5),seuilK('B',   188.0),&!5
      seuilK('C',   284.2), seuilK('N',   409.9),seuilK('O',   543.1),&
      seuilK('F',   696.7), seuilK('Ne',   870.2),&! 10
      seuilK('Na',  1070.8),seuilK('Mg',  1303.0),seuilK('Al',  1559.0),&
      seuilK('Si',  1839.0),seuilK('P',  2145.5),&! 15
      seuilK('S',  2472.0), seuilK('Cl',  2822.4),seuilK('Ar',  3205.9),&
      seuilK('K',  3608.4), seuilK('Ca',  4038.5),&! 20
      seuilK('Sc',  4492.0),seuilK('Ti',  4966.0),seuilK('V',  5465.0),&
      seuilK('Cr',  5989.0),seuilK('Mn',  6539.0),&! 25
      seuilK('Fe',  7112.0),seuilK('Co',  7709.0),seuilK('Ni',  8333.0),&
      seuilK('Cu',  8979.0),seuilK('Zn',  9659.0),&! 30
      seuilK('Ga', 10367.0),seuilK('Ge', 11103.0),seuilK('As', 11867.0),&
      seuilK('Se', 12658.0),seuilK('Br', 13474.0),&! 35
      seuilK('Kr', 14326.0),seuilK('Rb', 15200.0),seuilK('Sr', 16105.0),&
      seuilK('Y', 17038.0), seuilK('Zr', 17998.0),&! 40
      seuilK('Nb', 18986.0),seuilK('Mo', 20000.0),seuilK('Tc', 21044.0),&
      seuilK('Ru', 22117.0),seuilK('Rh', 23220.0),&! 45
      seuilK('Pd', 24350.0),seuilK('Ag', 25514.0),seuilK('Cd', 26711.0),&
      seuilK('In', 27940.0),seuilK('Sn', 29200.0),&! 50
      seuilK('Sb', 30491.0),seuilK('Te', 31814.0),seuilK('I', 33169.0),&
      seuilK('Xe', 34561.0),seuilK('Cs', 35985.0),&! 55
      seuilK('Ba', 37441.0),seuilK('La', 38925.0),seuilK('Ce', 40443.0),&
      seuilK('Pr', 41991.0),seuilK('Nd', 43569.0),&! 60
      seuilK('Pm', 45184.0),seuilK('Sm', 46834.0),seuilK('Eu', 48519.0),&
      seuilK('Gd', 50239.0),seuilK('Tb', 51996.0),&! 65
      seuilK('Dy', 53789.0),seuilK('Ho', 55618.0),seuilK('Er', 57486.0),&
      seuilK('Tm', 59390.0),seuilK('Yb', 61332.0),&! 70
      seuilK('Lu', 63314.0),seuilK('Hf', 65351.0),seuilK('Ta', 67416.0),&
      seuilK('W', 69525.0), seuilK('Re', 71676.0),&! 75
      seuilK('Os', 73871.0),seuilK('Ir', 76111.0),seuilK('Pt', 78395.0),&
      seuilK('Au', 80725.0),seuilK('Hg', 83102.0),&! 80
      seuilK('Tl', 85530.0),seuilK('Pb', 88005.0),seuilK('Bi', 90526.0),&
      seuilK('Po', 93105.0),seuilK('At', 95730.0),&! 85
      seuilK('Rn', 98404.0),seuilK('Fr',101137.0),seuilK('Ra',103922.0),&
      seuilK('Ac',106755.0),seuilK('Th',109651.0),&! 90
      seuilK('Pa',112601.0),seuilK('U',115606.0) /)

  character(len=2) :: sym1
  integer :: i

  do i = 1, Size_tab
!-------------------------------------------------------------------------------
! paranoid comparison, should prevent any problem with spaces
!-----------------------------------------------------------------
     sym1 = trim(adjustl(seuilK_tab(i)%name))
     if ( trim(adjustl(sym)) == sym1 ) then
        mygetK = seuilK_tab(i)%Kenergy
        return
     end if
  end do

  write(*,'(A)') 'Could not find element >'//trim(adjustl(sym))// &
       '< in the table of K edge energies!'
  stop

end function mygetK
