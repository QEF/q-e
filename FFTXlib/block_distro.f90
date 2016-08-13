!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=
   MODULE block_distro
!=----------------------------------------------------------------------=

!  ... Added by Eric Pascolo


   IMPLICIT NONE
   PRIVATE
   SAVE

   PUBLIC :: map_blocks, find_max,write_matrix_strange_idx

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=


SUBROUTINE map_blocks(maps,row_w,col_w,ub,lb,myid,num_of_core)

  IMPLICIT NONE
  INTEGER, INTENT(in) :: ub(3),lb(3),row_w,col_w,myid,num_of_core
  INTEGER, INTENT(inout) :: maps( lb(1): ub(1), lb(2):ub(2) )
  INTEGER :: i,nprow,npcol
  INTEGER :: row_el,col_el,block_num
  INTEGER, ALLOCATABLE :: proc_idx_matrix(:,:)
  
! GENERAZIONE DIMENSIONI BLOCCHI
  CALL set_block(nprow,npcol,num_of_core)
  
! GENERAZIONE INFO BLOCCHI
  CALL get_info_block(row_w,col_w,nprow,npcol,row_el,col_el,block_num)

 
!  WRITE(6,*) '  PARAMETRI DISTRIBUZIONE G scalapack like'
!  WRITE(6,*) '  -----------------------------------------'
!  WRITE(6,*) '  Number of CORE', num_of_core
!  WRITE(6,*) '  Number of ROW', ub(1)-lb(1)
!  WRITE(6,*) '  Number of ROW wave', row_w
!  WRITE(6,*) '  Number of COL', ub(2)-lb(2)
!  WRITE(6,*) '  Number of COL wave', col_w
!  WRITE(6,*) '  Number of core/ROW', nprow
!  WRITE(6,*) '  Number of core/COL', npcol
!  WRITE(6,*) '  SIZE ROW block', row_el
!  WRITE(6,*) '  SIZE COL block', col_el
!  WRITE(6,*) '  NUMBER OF BLOCK', block_num
!  WRITE(6,*) '  -----------------------------------------'

  ALLOCATE(proc_idx_matrix(nprow,npcol))
  
  proc_idx_matrix = 0
  maps = 0
        
! GENERAZIONE MATRICE PER DISTRIBUZIONE PROCESSORI 
  CALL set_matrix_processor(nprow,npcol,proc_idx_matrix)
  
! SCRITTURA SU FILE MATRICE PER DISTRIBUZIONE PROCESSORI
!  IF(myid .eq. 0) CALL write_matrix(proc_idx_matrix,nprow,npcol,"./mappe/pattern_proc.matrix")
  
! GENERAZIONE MATRICE MAPPA 
  CALL distro_matrix_processor(ub,lb,nprow,npcol,row_el,col_el,proc_idx_matrix,maps)
  
!  IF(myid .eq. 0) CALL write_matrix_strange_idx(maps,ub,lb,"./mappe/mappa_core_scalike.matrix")
    
  DEALLOCATE(proc_idx_matrix)
    
END SUBROUTINE map_blocks

!=----------------------------------------------------------------------=

SUBROUTINE set_block(npr,npc,num_of_core)
! GENERAZIONE DIMENSIONI BLOCCHI
  IMPLICIT NONE
  INTEGER, INTENT(in) :: num_of_core
  INTEGER, INTENT(inout) :: npr,npc
  INTEGER :: i,sqrtnp

  sqrtnp = INT( SQRT( REAL( num_of_core ) + 0.1 ) )
  
  DO i = 1, sqrtnp + 1
    IF( MOD( num_of_core, i ) == 0 ) npr = i
  ENDDO

  npc = num_of_core / npr

END SUBROUTINE set_block

!=----------------------------------------------------------------------=

SUBROUTINE get_info_block(row,col,npr,npc,elpr,elpc,epb)
! CALCOLO INFO BLOCCHI
  IMPLICIT NONE
  INTEGER, INTENT(in) :: row,col
  INTEGER, INTENT(in) :: npr,npc
  INTEGER, INTENT(out) :: elpr,elpc,epb
  INTEGER :: i,sqrtnp
  
  elpr = (row/npr);
  elpc = (col/npc);
  epb =  ((row*col)/(2*elpr*elpc));
  
END SUBROUTINE get_info_block

!=----------------------------------------------------------------------=

SUBROUTINE set_matrix_processor(npr,npc,pmatrix)
! GENERAZIONE MATRICE PROCESSORI
  IMPLICIT NONE
  INTEGER, INTENT(in) :: npr,npc
  INTEGER, INTENT(out) :: pmatrix(:,:)
  INTEGER :: i,j,p,row,col
  
  p = 1
  DO i=1,npc
    DO j=1,npr
        pmatrix(j,i) = p
        p= p+1
    ENDDO
  ENDDO    
  
  
END SUBROUTINE set_matrix_processor

!=----------------------------------------------------------------------=

SUBROUTINE distro_matrix_processor(ub,lb,npr,npc,elpr,elpc,pmatrix,proc_distro_matrix)

! DISTRO

  IMPLICIT NONE
  INTEGER, INTENT(in) :: npr,npc,elpr,elpc,ub(3),lb(3)
  INTEGER, INTENT(in) :: pmatrix(npr,npc)
  INTEGER, INTENT(inout) :: proc_distro_matrix(lb(1): ub(1), lb(2):ub(2))
  INTEGER :: i,j,i0,j0,i0b,j0b,ib,jb,jb1,row,col
  proc_distro_matrix = 0
  ib = 0
  jb = 0
  row = ub(1)
  col = ub(2)
  
!   RIEMPIMENTO CENTRO RIGHT
  DO i=elpr,row,elpr
    
    ib = MOD(ib,npr)+1;
    i0=i - elpr;
    jb = 0
    DO j=elpc,col,elpc
           
      j0=j - elpc;
      jb = MOD(jb,npc)+1;
      
      proc_distro_matrix(i0:i,j0:j) = pmatrix(ib,jb); 
      j0b = j
    ENDDO
    i0b = i
  ENDDO

     
  ib = 0
  
!   RIEMPIMENTO BORDO RIGHT LATERALE
  DO i=0,row
  
      proc_distro_matrix(i,j0b:col) = proc_distro_matrix(i,j0b-1) 
        
  ENDDO
    
!   RIEMPIMENTO BORDO RIGHT SOPRA
  DO j=0,col
      proc_distro_matrix(i0b:row,j) =  proc_distro_matrix(i0b-1,j);
  ENDDO
  
  !RIEMPIMENTO CENTRO LEFT
  ib = 0
  row = ub(1)
  col = lb(2)

  DO i=elpr,row,elpr
    
    ib = MOD(ib,npr)+1;
    jb1 = 0
    DO j=-elpc-1,col,-elpc
       
      i0=i - elpr;
      j0=j + elpc;
      jb = npc - MOD(jb1,npc);
      proc_distro_matrix(i0:i,j:j0) = pmatrix(ib,jb);
      jb1 = jb1+1
      j0b = j

    ENDDO
    
    i0b = i
    
  ENDDO
  
  ib = 2
  
  !RIEMPIMENTO BORDO LEFT LATERALE
  DO i=0,row
      proc_distro_matrix(i,col:j0b) = proc_distro_matrix(i,j0b+1) 
  ENDDO
    
  !RIEMPIMENTO BORDO LEFT SOTTO
  DO j=-1,col,-1
      proc_distro_matrix(i0b:ub(1),j) =  proc_distro_matrix(i0b-1,j);
  ENDDO
  
  
END SUBROUTINE distro_matrix_processor

!=----------------------------------------------------------------------=

SUBROUTINE find_max(ub,lb,mtw,row,col,my,ncore)

  INTEGER,INTENT(OUT) :: row,col
  INTEGER,INTENT(IN)  :: ub(:),lb(:),mtw(lb(1):ub(1),lb(2):ub(2))
  INTEGER,INTENT(IN)  :: my,ncore
  INTEGER :: i,j,max_proc,vector(0:ub(2)),ncr,ncc
  
  CALL set_block(ncr,ncc,ncore)
  vector = 0
  max_proc = max(ncr,ncc)
  
  WRITE(6,*) 'max proc',max_proc,'ncr',ncr,'ncc',ncc
  WRITE(6,*) 'ub1',ub(1),'ub2',ub(2)
  
  DO j=0,ub(2)
  
    DO i=0,ub(1)
    
    IF (mtw(i,j).ne. 0)THEN
      vector(j) = i
    ELSE
      EXIT
    ENDIF
    
    ENDDO
        
  ENDDO
  
  
  IF(my .eq. 1) WRITE(*,*) 'vector',vector(:)
  
  
  DO i=0,ub(2)
  
     IF( abs(i-vector(i)) .le. 1 .AND. i>max_proc .AND. vector(i)>max_proc ) THEN
     row = vector(i)
     col = i
     EXIT
     
     ENDIF
  
  ENDDO
  
 
  
END SUBROUTINE

!=----------------------------------------------------------------------=

SUBROUTINE write_matrix(mtw,righe,colonne,percorso)

    IMPLICIT NONE
    INTEGER :: i, j, righe,colonne
    character(LEN=*), INTENT(in) :: percorso
    INTEGER, INTENT(in) :: mtw(righe,colonne)
    
    OPEN(unit=115, file=percorso)
    DO i=1,righe
       WRITE(115,'(1000I6)') mtw(i,1:colonne)
    END DO
    CLOSE(115)

END SUBROUTINE write_matrix

!=----------------------------------------------------------------------=

SUBROUTINE write_matrix_strange_idx(mtw,ub,lb,percorso)

  IMPLICIT NONE
  INTEGER :: i, j, ub(3),lb(3)
  character(LEN=*), INTENT(in) :: percorso
  INTEGER, INTENT(in) :: mtw(lb(1):ub(1),lb(2):ub(2))
    
  OPEN(unit=116, file=percorso)

    DO i=ub(1),lb(1),-1
       WRITE(116,'(1000I6)') mtw(i,lb(2):ub(2))
    END DO
    
  CLOSE(116)

END SUBROUTINE write_matrix_strange_idx

!=----------------------------------------------------------------------=
   END MODULE block_distro
!=----------------------------------------------------------------------=
