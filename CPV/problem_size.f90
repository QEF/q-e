!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  MODULE problem_size

     IMPLICIT NONE
     SAVE
     PRIVATE
     PUBLIC :: cpsizes

  CONTAINS

      SUBROUTINE cpsizes() 
                                                                        
      USE kinds
      use ions_base,          only: nat, nsp
      use electrons_base,     only: nx => nbnd, nspin
      use gvecw,              only: ngwx
      use gvect,              only: ngmx
      use smallbox_gvec,              only: ngb
      use uspp_param,         only: nhm
      use uspp,               only: nkb
      USE io_global,          ONLY: ionode
      USE io_global,          ONLY: stdout
      USE fft_base,           ONLY: dfftp, dffts
 
      implicit none 
                                                                        
      integer nr1x, nr2x, nr3x, nr1_l, nr2_l, nr3_l
      integer nbyte 
      integer nbyte_alloc 
      integer itmp 

      nr1_l = dfftp%nr1x
      nr2_l = dfftp%nr2x
      nr3_l = dfftp%npl

      nr1x  = dfftp%nr1x
      nr2x  = dfftp%nr2x
      nr3x  = dfftp%nr3x

      nbyte         = 0 
      nbyte_alloc   = 0 
                                                                        
                                                                        
! ... Atoms type
      nbyte = nbyte + 8* 3 * 14 * nat

! ... GVEC
      nbyte = nbyte + 8 * 10 * ngb
      nbyte = nbyte + 8 * 13 * ngmx

! ... Pseudo
      nbyte = nbyte + 8 * 5 * nkb * nx * nspin

! ... C0 CM CP
      nbyte = nbyte + 3 * 16 * ngwx * nx * nspin

! ... ei1 ei2 ei3, eigr, sfac
      nbyte = nbyte + 3 * 16 * MAX( nr1x, nr2x, nr3x ) * nat
      nbyte = nbyte + 16 * ngwx * nat
      nbyte = nbyte + 16 * ngmx * nsp

! ... rhoe and vpot ( nr1_l, nr2_l, nr3_l, nspin )
      nbyte = nbyte + ( 8 + 16 ) * NR1_L * NR2_L * NR3_L * nspin 

! ... TEMPORARY ALLOCATED MEMORY
                                                                        
! ... ortho                                                             
      itmp = 8 * 8 * NX * NX
      if(itmp.gt.nbyte_alloc) nbyte_alloc = itmp 

! ... pvofrho & pstress                                                 
      itmp = 8 * ( NR1_L * NR2_L * NR3_L * 8 + &
     &       nat * NX * nhm * 6 + 6 * ngmx + 6*ngwx +                   &
     &       ngwx*nhm*nsp + 2*ngwx*nat )
      if(itmp.gt.nbyte_alloc) nbyte_alloc = itmp 

      IF(ionode) THEN
        WRITE( stdout,10) nbyte + nbyte_alloc 
      END IF
                                                                        
   10 FORMAT(//,3X,'Estimated Sizes of the problem',/                   &
     &         ,3X,'------------------------------',/                   &
     &         ,3X,'dimension of the problem (byte/pe) : ',I12)
                                                                        
      RETURN 
      END SUBROUTINE cpsizes

  END MODULE problem_size
