
SUBROUTINE enFdos
  
  USE kinds,        ONLY : DP
  USE parameters,   ONLY : npk
  USE cell_base,    ONLY : at, bg
  USE control_flags,ONLY : lphonon
  USE klist,        ONLY : nks, nkstot, xk, wk
  USE ions_base,    ONLY : nat
  USE wvfct,        ONLY : et, nbnd
  USE ktetra,       ONLY : nk1, nk2, nk3 
  USE lsda_mod,     ONLY : isk
  USE symme,        ONLY : s, nsym, minus_q, irt
  implicit none
  !
  INTEGER :: iuna2Fsave  = 40, i, j, ik, ns, na
  logical  ::  exst
  !
  !
  CALL seqopn( iuna2Fsave, 'a2Fsave', 'FORMATTED', exst )
  !===========================================
  !
  WRITE( iuna2Fsave, * ) nbnd, nkstot
  WRITE( iuna2Fsave, * ) et
  WRITE( iuna2Fsave, * ) ((xk(i,ik), i=1,3), ik=1,nkstot)
  WRITE( iuna2Fsave, * ) wk(1:nkstot)
  WRITE( iuna2Fsave, * ) nk1, nk2, nk3
  !
  WRITE( iuna2Fsave, * ) nsym
  do ns=1,nsym
     WRITE( iuna2Fsave, * )  ((s(i,j,ns),j=1,3),i=1,3) 
  enddo
  WRITE( iuna2Fsave, * )  ((irt(ns,na),ns=1,nsym),na=1,nat)
  !
  CLOSE( UNIT = iuna2Fsave, STATUS = 'KEEP' )
  !  
  RETURN
END SUBROUTINE enFdos
