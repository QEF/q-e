  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------
  SUBROUTINE dasmio ( mat, nsize, lrec, iun, nrec, iop )  
  !--------------------------------------------------------------
  !
  ! A simple wrapper to the davcio routine to read/write square 
  ! matrices instead of vectors (Direc Access Square Matrix I/O)
  ! by  FG
  !
  ! iop =-1 : read
  ! iop = 1 : write
  !--------------------------------------------------------------
  USE kinds, ONLY : DP
  implicit none
  integer :: nsize, lrec, iun, nrec, iop, i
  complex(kind=DP):: mat(nsize,nsize), aux ( nsize*nsize ) 
  !
  IF ( iop .eq. -1 ) then
     !
     !  read matrix
     !
     CALL davcio ( aux, lrec, iun, nrec, -1 )
     DO i = 1, nsize * nsize
        mat ( 1 + (i-1)/nsize , i - nsize*((i-1)/nsize) ) = aux (i)
     ENDDO
     !
  ELSEif ( iop .eq. 1 ) then
     !
     !  write matrix
     !
     DO i = 1, nsize * nsize
        aux (i) = mat ( 1 + (i-1)/nsize , i - nsize*((i-1)/nsize) ) 
     ENDDO
    CALL davcio ( aux, lrec, iun, nrec, +1 )
    !
 ELSE
    !
    CALL errore ('dasmio','iop not permitted',1)
    !
 ENDIF
 !
END SUBROUTINE dasmio
!
