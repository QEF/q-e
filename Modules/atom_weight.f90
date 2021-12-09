!
! Copyright (C) 2004-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ------------------------------------------------------------------
function atom_weight(atomic_number)
  ! ------------------------------------------------------------------
  !! Provides the standard atomic weights, given the atomic number.
  !
  USE kinds, ONLY : dp
  implicit none
  integer :: atomic_number
  !! input: atomic number
  real(DP) :: atom_weight
  !! output: atomic weight
  
  real(DP) :: weights(109)
  data weights/ 1.00794_DP,                           4.00260_DP, &
                6.941_DP,9.01218_DP,10.811_DP,12.0107_DP,14.00674_DP, &
                15.9994_DP,18.99840_DP,20.1797_DP, &
                22.98977_DP,24.3050_DP,26.98154_DP,28.0855_DP,30.97376_DP, &
                32.066_DP,35.4527_DP,39.948_DP, &
                39.0983_DP,40.078_DP,44.95591_DP,47.867_DP,50.9415_DP, &
                51.9961_DP,54.93805_DP, 55.845_DP,      &
                58.93320_DP,58.6934_DP,63.546_DP,65.39_DP,      &
                69.723_DP,72.61_DP,74.92160_DP,78.96_DP,79.904_DP,83.80_DP, &
                85.4678_DP,87.62_DP,88.90585_DP,91.224_DP,92.90638_DP, &
                95.94_DP,98._DP,      &
                101.07_DP,102.90550_DP,106.42_DP,107.8682_DP,112.411_DP,    &
                114.818_DP,118.710_DP,121.760_DP,127.60_DP,126.90447_DP, &
                131.29_DP, &
                132.90545_DP,137.327_DP,138.9055_DP,140.116_DP,140.90765_DP, &
                144.24_DP,145._DP,150.36_DP,151.964_DP,157.25_DP, &
                158.92534_DP,162.50_DP,164.93032_DP,167.26_DP,   &
                168.93421_DP,173.04_DP,174.967_DP, &
                178.49_DP,180.9479_DP,183.84_DP,186.207_DP,190.23_DP, &
                192.217_DP,195.078_DP,196.96655_DP,200.59_DP,           &
                204.3833_DP,207.2_DP,208.98038_DP,209._DP,210._DP,222._DP, &
                223._DP,226._DP,227._DP,232.0381_DP,231.03588_DP, &
                238.0289_DP,237._DP,244._DP, &
                243._DP,247._DP,247._DP,251._DP,252._DP,257._DP,  &
                258._DP,259._DP,262._DP,261._DP,262._DP,266._DP,264._DP, &
                277._DP,268._DP  /

  if (atomic_number < 1 .or. atomic_number > 109) then
     call errore('atom_weight','invalid atomic number',abs(atomic_number))
  else
     atom_weight=weights(atomic_number)
  end if
  return

end function atom_weight
!
