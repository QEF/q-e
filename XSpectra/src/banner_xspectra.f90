!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE banner_xspectra()

  USE io_global,   ONLY : stdout

  IMPLICIT NONE 
  WRITE(stdout,'(/,5x,a)')&
  '-------------------------------------------------------------------------'
  WRITE (stdout,'(20x,2a1,2x,4a1,17x,a1)') "_","_","_","_","_","_","_"
  WRITE (stdout,&
  '(20x,a1,x,2a1,x,a1,x,3a1,x,2a1,3x,3a1,2x,4a1,x,2a1,x,a,x,2a1,x,2a1,x,a1)')&
  "\\","\\","/","/","_","\\","_","_","_","_","_","_","_","_","_","|","|","_","_",&
  "_","_","_","_","_"
  WRITE (stdout,&
  '(21x,a1,2x,3(2a1,x),3(a1,x),2a1,2(x,3a1),x,4a1,x,2a1,x,a1)')& 
  "\\","/","\\","\\","|","'","_","\\","/","_","\\","/","_","_","|","_","_","|",&
  "\'","_","_","/","_","\`","|"
  WRITE (stdout,&
  '(21x,a1,2x,3a1,x,a1,x,3a1,x,a1,2x,3a1,x,4a1,x,3a1,2(x,a1),x,3a1,x,a1)') &
  "/","\\","_","\\","\\","|","_",")","|","_","_","/","(","_","_","|","|","_","|",&
  "|","|","(","_","|","|"
  WRITE (stdout,&
  '(20x,9a1,x,4a1,x,16a1,2x,6a1)') &
  "/","_","/","\\","_","\\","_","_","/",".","_","_","/","\\","_","_","_","|","\\",&
  "_","_","_","|","\\","_","_","|","_","|","\\","_","_",",","_","|"
  WRITE (stdout,'(28x,3a1)') "|","_","|"
  
  WRITE (stdout, '(/,5x,a)')&
  'In publications arising from the use of XSpectra, please cite:'
  WRITE (stdout, '(6x,a)')&
  '- O. Bunau and M. Calandra,'
  WRITE (stdout, '(8x,a)')&
  'Phys. Rev. B 87, 205105 (2013)'
  WRITE (stdout, '(6x,a)')&
  '- Ch. Gougoussis, M. Calandra, A. P. Seitsonen, F. Mauri,'
  WRITE (stdout, '(8x,a)')&
  'Phys. Rev. B 80, 075102 (2009)'
  WRITE (stdout, '(6x,a)')&
  '- M. Taillefumier, D. Cabaret, A. M. Flank, and F. Mauri,'
  WRITE (stdout, '(8x,a)')&
  'Phys. Rev. B 66, 195107 (2002)' 

  
end SUBROUTINE banner_xspectra 
