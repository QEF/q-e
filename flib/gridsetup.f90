!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
!-----------------------------------------------------------------------
!

      subroutine GRIDSETUP(NPROC,NPROW,NPCOL)

!
! This subroutine factorizes the number of processors (NPROC)
! into NPROW and NPCOL,  that are the sizes of the 2D processors mesh.
!
! Written by Carlo Cavazzoni
!

      IMPLICIT NONE

      integer nproc,nprow,npcol
      integer sqrtnp,i

      if(nproc.lt.2) then
        npcol = 1
        nprow = 1
      else
        sqrtnp = int( sqrt( dble(nproc) ) + 1 )
        do i=1,sqrtnp
          if(mod(nproc,i).eq.0) nprow = i
        end do
        npcol = nproc/nprow
      endif

      return
      end

